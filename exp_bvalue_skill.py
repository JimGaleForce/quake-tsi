"""P-2.2 / K-405 -- MAGNITUDE AS THE RESPONDER, tested as forecast skill.

Every statistic this program owns is an OCCURRENCE statistic: it asks when and where,
never how big. K-405 says that is a blind spot with a physical reason behind it -- a
threshold system whose rate is unmoved can still have its SIZE DISTRIBUTION moved. This
arm asks the smallest honest version of that question, and it asks it as forecast skill
rather than as a correlation:

    Does a LOCAL b-value, estimated only from the past, forecast the size of the next
    earthquake better than one global b?

The unit is bits/event, walk-forward, on the exploration split only (S-3: the gate is
incremental bits, never AUC, never a correlation coefficient).

---------------------------------------------------------------------------------------
THE ARTIFACT THAT COULD FAKE THIS, NAMED BEFORE THE RUN (PLAYBOOK rule 7)
---------------------------------------------------------------------------------------
SPATIAL AND TEMPORAL VARIATION IN COMPLETENESS MASQUERADING AS VARIATION IN b. This is
the dominant failure mode of the entire b-value literature and it is the reason
Gulia & Wiemer sits at CONTESTED in this program's ledger (HYPOTHESIS_LEDGER.md
L3195-3205). The mechanism is exact and it is not subtle: if the true Mc in a sub-region
is higher than the Mc our estimator assigns, the surviving magnitude sample is depleted
at its low end, the Aki-Utsu mean magnitude rises, and b-hat falls. A map of b(x) built
this way is a map of the seismic NETWORK, not of the crust. Because network density is
spatially smooth and persistent in time, a local estimator built from the past will
"forecast" it, and the forecast will earn real bits. Those bits are an instrument
reading, not a fact about the Earth.

Three declared controls, all reported whether they help or hurt:

  (a) THE FLOOR-RAISE. Re-run with every floor (local and global, estimator and
      likelihood) raised by +0.3. If the local-b gain collapses when the sample is
      pushed further above completeness, the gain was completeness.
  (b) MAGNITUDE ROUNDING. The likelihood must be bin-aware or the comparison is not a
      comparison of probabilities at all. The native rounding of each catalogue is
      MEASURED (not assumed) and the likelihood is re-run at that native bin width as a
      control on the 0.1-mag primary.
  (c) THE DAY/NIGHT DETECTION ARTIFACT BELOW M2.5. SoCal completeness is diurnal at
      small magnitudes. The gain restricted to scored events with M >= 2.5 -- far above
      any plausible Mc anywhere in either network -- is reported separately.

Two further guards that are part of the design rather than controls:

  * S-4 (n-bias). b-hat's variance and its small-sample bias both scale with the number
    of events used. A local estimator whose n varies with local rate would manufacture
    exactly the structure we are looking for. n is therefore FIXED AT 250 BY
    CONSTRUCTION -- k-nearest, not a fixed radius -- so the estimator's sampling
    properties are identical everywhere. (The count above the local Mc still varies;
    the shuffle null absorbs that, since it runs the identical estimator.)
  * PLAYBOOK rule 1 (sequence, not event). Clustering is KEPT here on purpose -- the
    magnitudes of aftershocks are the target, not the nuisance -- so the events are
    emphatically not independent. Every CI is a block bootstrap over sequences.

---------------------------------------------------------------------------------------
THE NULL (PLAYBOOK rules 2/5)
---------------------------------------------------------------------------------------
Magnitudes are shuffled among events sharing a (0.5 deg cell, calendar year) group. That
destroys any real b structure while preserving each cell-epoch's magnitude distribution,
i.e. preserving completeness. The identical pipeline -- identical neighbour sets,
identical estimator, identical scoring rule -- must then return ~0 bits/event. 100
shuffles give the null distribution and the p of the real gain.

---------------------------------------------------------------------------------------
SENSITIVITY BEFORE THE NULL IS QUOTED (PLAYBOOK rule 5)
---------------------------------------------------------------------------------------
A known smooth b(x) field is injected into a magnitude-shuffled catalogue and the
identical pipeline is run. The minimum delta-b detectable at 80% power is reported. A
null quoted without this number is a statement about the instrument, not the Earth.

---------------------------------------------------------------------------------------
FROZEN BEFORE THE RUN
---------------------------------------------------------------------------------------
  catalogues        QTM_12dev.txt (full, clustering kept), SCSN_original_catalog.txt
  exploration       first 70% of each catalogue's own time span; last 30% NEVER read
  training          first 60% of the exploration span (global b, all tuning)
  scoring           remaining 40% of the exploration span, walk-forward
  local sample      k = 250 nearest prior events, strictly t < t_i, <= 30 km horizontal,
                    within the prior 730 days; fewer than 250 -> no local estimate
  local Mc          maximum curvature on 0.1-mag bins, + 0.2
  local b           Aki-Utsu on the subset >= local Mc, bin-width corrected
  global Mc, b      the same rule on the whole training window
  scored events     M_i >= local Mc AND M_i >= global Mc
  likelihood        bin-aware truncated Gutenberg-Richter, floor = local Mc for BOTH
                    models (same support), upper truncation M = 8.0
  primary statistic mean (log2 p_local - log2 p_global), bits/event
  CI                95%, block bootstrap over sequences, 2000 resamples
  SUCCESS RULE      gain > 0.02 bits/event with the CI excluding 0, on BOTH catalogues

---------------------------------------------------------------------------------------
ADDED AFTER THE FIRST RUN, AND LABELLED AS SUCH
---------------------------------------------------------------------------------------
Nothing below changes the primary, the null or controls (a)/(b)/(c); with the fixed
seeds those reproduce bit-identically. All of it is reported as SECONDARY.

  * A second power definition in the sensitivity arm. The frozen definition (gain > 0.02
    bits AND above the shuffle null's 95th percentile) returned ZERO power at every
    injected delta up to 1.0 on QTM. That is itself the finding -- the frozen 0.02-bit
    bar is far above what any plausible b(x) contrast can pay -- but it leaves rule 5
    without a usable number, so a matched reference (the 95% point of the delta = 0
    injection replicates) is reported alongside it.
  * MERTON PRIOR ART, folded in on instruction. Herrmann & Marzocchi (2021) report that
    the QTM catalogue has the same Lilliefors completeness as SCSN, near M 3.24: the
    magnitudes are NOT exponential down to a maximum-curvature floor. So this arm
    measures a Lilliefors-style exponential goodness-of-fit p on each training window,
    scans the floor upward until the exponential survives, reports the primary gain
    restricted to M >= 3.0, and states which floor the headline is defensible on.
    van der Elst (2021) b-positive -- the estimator the field now prefers under
    incompleteness -- is run as a secondary local estimator with everything else
    identical.

Exploration split only. HOLDOUT IS NOT TOUCHED.
"""

from __future__ import annotations

import datetime as _dt
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "xue_lu_zenodo"
OUT_JSON = HERE / "results_bvalue_skill.json"

# ----------------------------------------------------------------- FROZEN CONSTANTS
EXPLORE_FRAC = 0.70           # exploration = first 70% of the catalogue span
TRAIN_FRAC_OF_EXPLORE = 0.60  # training = first 60% of the exploration span
K_LOCAL = 250                 # fixed n for the local estimator (S-4)
RADIUS_KM = 30.0              # horizontal radius for the local sample
LOOKBACK_DAYS = 730.0         # temporal lookback for the local sample
BIN_W = 0.1                   # analysis magnitude bin width
MC_OFFSET = 0.2               # maximum curvature + 0.2
MMAX_TRUNC = 8.0              # upper truncation of the GR likelihood
MIN_NSUB = 30                 # minimum events above local Mc for a usable local b
B_BOUNDS = (0.2, 3.0)         # a local b outside this is treated as unusable
FLOOR_RAISE = 0.3             # control (a)
HIGH_M_CUT = 2.5              # control (c)
HIGH_M_CUT_2 = 3.0            # Merton follow-up: above Herrmann & Marzocchi's Mc=3.24-ish
LILLIE_N = 20000              # fixed-n subsample for the exponentiality test
LILLIE_SIM = 200              # Monte-Carlo replicates for the Lilliefors critical value
LILLIE_SCAN = 20              # floors scanned upward, in 0.2-mag steps
BPOS_DMC = 0.1                # b-positive difference threshold (van der Elst 2021)
NULL_CELL_DEG = 0.5           # shuffle group: 0.5 deg cell x calendar year
N_NULL = 100
N_BOOT = 2000
SEQ_LINK_KM = 50.0            # sequence blocks for the bootstrap
SEQ_LINK_DAYS = 30.0
COARSE_BLOCK_DEG = 0.5        # robustness block bootstrap
COARSE_BLOCK_DAYS = 180.0
SCORING_CAP = 120000          # fixed-seed subsample of the scoring window
SEED = 20260902
SUCCESS_BITS = 0.02
INJECT_FLOOR_BELOW = 0.5      # injection floor sits this far BELOW the global Mc, so
#                               that every local sample above its own Mc is a pure GR
#                               and delta = 0 therefore has to return ~0 bits
INJECT_DELTAS = [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.45, 0.60, 0.80, 1.00]
INJECT_REPS = 15
POWER_TARGET = 0.80
MAP_CELL_DEG = 0.2            # descriptive b(x) map
DEPTH_BIN_KM = 2.0            # descriptive b vs depth
DESC_MIN_N = 250
DESC_BOOT = 200

LN10 = np.log(10.0)
LN2 = np.log(2.0)

RAW_COLS = ["yr", "mo", "dy", "hr", "mi", "sec", "eid", "lat", "lon", "depth", "mag"]
DEC_COLS = ["yr", "mo", "dy", "hr", "mi", "sec", "lat", "lon", "depth", "mag", "eid"]


# ------------------------------------------------------------------------ loading
def load_catalog(fname):
    """Auto-detect column order exactly as xue_lu_crosstest.load_catalog / exp_h_etas."""
    probe = pd.read_csv(DATA / fname, sep=r"\s+", header=None, nrows=1000)
    is_raw = probe[6].abs().max() > 90
    cols = RAW_COLS if is_raw else DEC_COLS
    df = pd.read_csv(DATA / fname, sep=r"\s+", header=None, names=cols)
    assert df.lat.abs().max() <= 90 and df.lon.abs().max() <= 180, \
        "column detection failed for %s" % fname
    ts = pd.to_datetime(dict(year=df.yr, month=df.mo, day=df.dy,
                             hour=df.hr, minute=df.mi, second=0), utc=True) \
        + pd.to_timedelta(df["sec"].astype(float), unit="s")
    df["t"] = (ts - pd.Timestamp(0, tz="UTC")).dt.total_seconds().to_numpy() / 86400.0
    df["year"] = ts.dt.year.to_numpy()
    df = df.sort_values("t", kind="mergesort").reset_index(drop=True)
    assert len(df) > 100000, "expected a large catalogue, got %d rows" % len(df)
    return df, ("raw" if is_raw else "declustered")


def measure_rounding(mag):
    """MEASURE the native magnitude bin width rather than assume it (control b)."""
    for w in (0.2, 0.1, 0.05, 0.02, 0.01, 0.001):
        if np.mean(np.abs(mag / w - np.round(mag / w)) < 1e-6) > 0.999:
            return w
    return None


def project_km(lat, lon):
    lat0 = float(np.median(lat))
    lon0 = float(np.median(lon))
    x = (lon - lon0) * 111.320 * np.cos(np.radians(lat0))
    y = (lat - lat0) * 110.574
    return x, y, lat0, lon0


# ---------------------------------------------------------------- b-value estimator
def maxc_mc(mag, binw=BIN_W, offset=MC_OFFSET):
    """Maximum-curvature Mc: the modal bin of the non-cumulative GR histogram, + offset."""
    q = np.round(np.asarray(mag) / binw).astype(np.int64)
    qmin = int(q.min())
    counts = np.bincount(q - qmin)
    return (int(np.argmax(counts)) + qmin) * binw + offset


def aki_b(mag, mc, binw=BIN_W):
    """Aki-Utsu maximum-likelihood b with the Utsu/Bender bin-width correction."""
    mag = np.asarray(mag, dtype=float)
    sub = mag[mag >= mc - 1e-9]
    if sub.size < MIN_NSUB:
        return float("nan"), int(sub.size)
    denom = sub.mean() - (mc - binw / 2.0)
    if denom <= 0:
        return float("nan"), int(sub.size)
    return 1.0 / (LN10 * denom), int(sub.size)


def gr_bin_logprob(mbin, mc, b, binw=BIN_W, mmax=MMAX_TRUNC):
    """log_e P(magnitude falls in the bin centred on mbin) under a truncated GR.

    Support: bins from the one centred on `mc` up to the one centred on `mmax`, so the
    two competing models (local b, global b) see EXACTLY the same support.
    """
    beta = np.asarray(b, dtype=float) * LN10
    lo = np.asarray(mc, dtype=float) - binw / 2.0
    a = np.asarray(mbin, dtype=float) - binw / 2.0 - lo
    top = np.asarray(mmax, dtype=float) + binw / 2.0 - lo
    # exp(-beta*a) - exp(-beta*(a+binw)) = exp(-beta*a) * (1 - exp(-beta*binw))
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        log_num = -beta * a + np.log(-np.expm1(-beta * binw))
        log_den = np.log(-np.expm1(-beta * top))
        return log_num - log_den


# -------------------- MERTON FOLLOW-UP: is the magnitude distribution even exponential?
def _lilliefors_D(x):
    """Kolmogorov-Smirnov distance to an exponential whose rate is fitted to x itself."""
    x = np.sort(x)
    n = x.size
    lam = 1.0 / x.mean()
    F = -np.expm1(-lam * x)
    i = np.arange(1, n + 1)
    return float(max(np.max(i / n - F), np.max(F - (i - 1) / n)))


def lilliefors_exp_p(mag, floor, binw=BIN_W, n_sub=LILLIE_N, n_sim=LILLIE_SIM,
                     seed=SEED):
    """Lilliefors-style test that magnitudes above `floor` are exponential (i.e. GR).

    Herrmann & Marzocchi (2021) apply exactly this test and report that BOTH the QTM and
    the SCSN SoCal catalogues are only Lilliefors-complete near M 3.24 -- far above any
    maximum-curvature floor. If that holds on our own training windows, then the local
    and global "b" in this arm are both fits of an exponential to a distribution that is
    not exponential. The relative log-score comparison remains fair; the b-value
    INTERPRETATION does not survive.

    Discreteness is broken by uniform dithering inside the magnitude bin, as in H&M.
    Critical values come from Monte Carlo, because the rate is estimated from the sample.
    """
    rng = np.random.default_rng(seed)
    sub = np.asarray(mag, dtype=float)
    sub = sub[sub >= floor - 1e-9]
    if sub.size < 500:
        return {"floor": float(floor), "n": int(sub.size), "p": None,
                "note": "too few events"}
    n_used = int(min(n_sub, sub.size))
    sub = rng.choice(sub, size=n_used, replace=False)
    x = sub + rng.uniform(-binw / 2.0, binw / 2.0, n_used) - (floor - binw / 2.0)
    x = x[x > 0]
    D = _lilliefors_D(x)
    lam = 1.0 / x.mean()
    sim = rng.exponential(1.0 / lam, size=(n_sim, x.size))
    Ds = np.array([_lilliefors_D(r) for r in sim])
    return {"floor": float(floor), "n": int(x.size), "D": D,
            "p": float((np.sum(Ds >= D) + 1) / (n_sim + 1)),
            "b_implied": float(lam / LN10)}


def lilliefors_completeness_scan(mag, start_floor, binw=BIN_W, steps=LILLIE_SCAN,
                                 seed=SEED):
    """Raise the floor in 0.2-mag steps until the exponential is no longer rejected.
    This is Herrmann & Marzocchi's Mc, computed on OUR training window."""
    out = []
    mc_lil = None
    for k in range(steps):
        f = start_floor + 0.2 * k
        r = lilliefors_exp_p(mag, f, binw=binw, seed=seed + k)
        out.append(r)
        if r["p"] is not None and r["p"] > 0.05 and mc_lil is None:
            mc_lil = f
            break
        if r.get("n", 0) < 500:
            break
    return {"scan": out, "lilliefors_mc": mc_lil}


def b_positive(q, dmc_bins=1, binw=BIN_W):
    """van der Elst (2021) b-positive on a 1-D sequence of magnitude bin indices in
    TIME order: use only the positive differences between successive events, which are
    exponentially distributed with the same b and are immune to a changing Mc."""
    d = np.diff(np.asarray(q, dtype=np.int64))
    d = d[d >= dmc_bins]
    if d.size < MIN_NSUB:
        return float("nan"), int(d.size)
    denom = d.mean() * binw - (dmc_bins * binw - binw / 2.0)
    if denom <= 0:
        return float("nan"), int(d.size)
    return 1.0 / (LN10 * denom), int(d.size)


def b_positive_rows(qn_time_ordered, dmc_bins=1, binw=BIN_W):
    """Vectorised b-positive for every row of an [n, k] array already in time order."""
    d = np.diff(qn_time_ordered.astype(np.int64), axis=1)
    m = d >= dmc_bins
    ns = m.sum(axis=1)
    s = np.where(m, d, 0).sum(axis=1) * binw
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = s / np.maximum(ns, 1) - (dmc_bins * binw - binw / 2.0)
        b = np.where(ns > 0, 1.0 / (LN10 * denom), np.nan)
    ok = (ns >= MIN_NSUB) & np.isfinite(b) & (b > B_BOUNDS[0]) & (b < B_BOUNDS[1])
    return b, ns, ok


# -------------------------------------------------------- vectorised local estimator
def local_mc_b(qn, binw=BIN_W, offset_bins=2, raise_bins=0):
    """Local Mc (max curvature + 0.2) and Aki-Utsu b for every row of an [n, k] array
    of magnitude BIN INDICES (integer, units of `binw`).

    `raise_bins` implements control (a): every local floor raised by that many bins.
    Returns (mc_bin, b, n_sub, ok).
    """
    n = qn.shape[0]
    qmin = int(qn.min())
    nb = int(qn.max()) - qmin + 1
    mc_bin = np.empty(n, dtype=np.int64)
    b = np.empty(n, dtype=np.float64)
    n_sub = np.empty(n, dtype=np.int32)
    step = int(max(2000, min(20000, 2.0e7 // max(nb, 1))))
    for a in range(0, n, step):
        blk = qn[a:a + step].astype(np.int64)
        z = blk - qmin
        nr = z.shape[0]
        flat = (np.arange(nr, dtype=np.int64)[:, None] * nb + z).ravel()
        hist = np.bincount(flat, minlength=nr * nb).reshape(nr, nb)
        modal = np.argmax(hist, axis=1)
        mcb = modal + qmin + offset_bins + raise_bins
        mask = blk >= mcb[:, None]
        ns = mask.sum(axis=1)
        s = np.where(mask, blk, 0).sum(axis=1) * binw
        with np.errstate(divide="ignore", invalid="ignore"):
            mean = s / np.maximum(ns, 1)
            denom = mean - (mcb * binw - binw / 2.0)
            bb = np.where(ns > 0, 1.0 / (LN10 * denom), np.nan)
        mc_bin[a:a + step] = mcb
        b[a:a + step] = bb
        n_sub[a:a + step] = ns
    ok = (n_sub >= MIN_NSUB) & np.isfinite(b) & (b > B_BOUNDS[0]) & (b < B_BOUNDS[1])
    return mc_bin, b, n_sub, ok


# ------------------------------------------------------------------ neighbour search
def neighbour_sets(pool, t, x, y, chunk=1500, kq=600, verbose=True):
    """[n_pool, K_LOCAL] int32 of the K_LOCAL nearest strictly-prior events within
    RADIUS_KM and LOOKBACK_DAYS, plus a validity mask.

    Exact by construction: the KD-tree covers only events strictly before the chunk, the
    in-chunk priors are handled by brute force, and any row where the tree query may have
    truncated at k is redone by exact brute force.
    """
    n = pool.size
    out = np.full((n, K_LOCAL), -1, dtype=np.int32)
    ok = np.zeros(n, dtype=bool)
    t0 = time.time()
    n_redo = 0
    for a in range(0, n, chunk):
        idx = pool[a:a + chunk]
        nr = idx.size
        i0 = int(idx[0])
        lo = int(np.searchsorted(t, t[i0] - LOOKBACK_DAYS))
        pts = np.c_[x[idx], y[idx]]
        if i0 > lo:
            kk = min(kq, i0 - lo)
            tree = cKDTree(np.c_[x[lo:i0], y[lo:i0]])
            dd, ii = tree.query(pts, k=kk, distance_upper_bound=RADIUS_KM, workers=4)
            dd = np.atleast_2d(dd)
            ii = np.atleast_2d(ii)
            finite = np.isfinite(dd)
            gi = lo + np.where(finite, np.minimum(ii, i0 - lo - 1), 0)
            good = finite & (t[gi] >= (t[idx] - LOOKBACK_DAYS)[:, None])
            dd_old = np.where(good, dd, np.inf)
            gi_old = np.where(good, gi, -1)
            full = (finite.sum(axis=1) >= kk) & (kk < i0 - lo)
        else:
            dd_old = np.full((nr, 1), np.inf)
            gi_old = np.full((nr, 1), -1, dtype=np.int64)
            full = np.zeros(nr, dtype=bool)
        i1 = int(idx[-1])
        rec = np.arange(i0, i1 + 1)
        if rec.size:
            dxr = x[idx][:, None] - x[rec][None, :]
            dyr = y[idx][:, None] - y[rec][None, :]
            dr = np.sqrt(dxr * dxr + dyr * dyr)
            tt = t[idx][:, None]
            okr = (t[rec][None, :] < tt) & (t[rec][None, :] >= tt - LOOKBACK_DAYS) \
                & (dr <= RADIUS_KM)
            dd_rec = np.where(okr, dr, np.inf)
            gi_rec = np.where(okr, rec[None, :], -1)
        else:
            dd_rec = np.full((nr, 1), np.inf)
            gi_rec = np.full((nr, 1), -1, dtype=np.int64)
        D = np.concatenate([dd_old, dd_rec], axis=1)
        G = np.concatenate([gi_old, gi_rec], axis=1)
        cnt = np.isfinite(D).sum(axis=1)
        redo = np.where(full & (cnt < K_LOCAL))[0]
        for r in redo:
            n_redo += 1
            i = int(idx[r])
            lo_i = int(np.searchsorted(t, t[i] - LOOKBACK_DAYS))
            if i - lo_i < K_LOCAL:
                cnt[r] = 0
                continue
            dxx = x[lo_i:i] - x[i]
            dyy = y[lo_i:i] - y[i]
            drr = np.sqrt(dxx * dxx + dyy * dyy)
            sel = np.where(drr <= RADIUS_KM)[0]
            if sel.size < K_LOCAL:
                cnt[r] = sel.size
                continue
            pick = sel[np.argsort(drr[sel], kind="stable")[:K_LOCAL]] + lo_i
            out[a + r] = pick
            ok[a + r] = True
            cnt[r] = 0
        take = (cnt >= K_LOCAL)
        if take.any():
            sub = np.where(take)[0]
            part = np.argpartition(D[sub], K_LOCAL - 1, axis=1)[:, :K_LOCAL]
            out[a + sub] = np.take_along_axis(G[sub], part, axis=1).astype(np.int32)
            ok[a + sub] = True
        if verbose and (a // chunk) % 25 == 0:
            print("    neighbours %7d/%7d  ok=%.3f  redo=%d  %.0fs"
                  % (min(a + nr, n), n, ok[:a + nr].mean(), n_redo, time.time() - t0),
                  flush=True)
    assert ok.sum() > 0, "expected many events with a local estimate, got 0"
    assert (out[ok] >= 0).all(), "negative neighbour index survived"
    assert (out[ok] < len(t)).all(), "out-of-range neighbour index"
    assert (out[ok] < pool[ok][:, None]).all(), "a non-prior neighbour leaked in"
    assert (t[out[ok]] >= (t[pool[ok]] - LOOKBACK_DAYS)[:, None]).all(), \
        "a neighbour outside the 730 d lookback leaked in"
    return out, ok, n_redo


# --------------------------------------------------------------- sequence blocking
def sequence_ids(t, lat, lon, mag, link_km=SEQ_LINK_KM, link_days=SEQ_LINK_DAYS):
    """Root-of-the-triggering-chain sequence id, the same rule as
    exp_nearcritical.classify / exp_fluid_driven._sequence_ids at a SoCal-appropriate
    radius. Applied to the scored events, which is the population the bootstrap
    resamples."""
    n = t.size
    assert np.all(np.diff(t) >= 0), "sequence_ids needs time-sorted input"
    seq = np.arange(n)
    lo_all = np.searchsorted(t, t - link_days, side="left")
    r2 = link_km * link_km
    coslat = np.cos(np.radians(lat))
    for ii in range(n):
        lo = int(lo_all[ii])
        if lo >= ii:
            continue
        cand = slice(lo, ii)
        m_ok = mag[cand] >= mag[ii]
        if not m_ok.any():
            continue
        dx = (lon[cand] - lon[ii]) * 111.320 * coslat[ii]
        dy = (lat[cand] - lat[ii]) * 110.574
        hit = m_ok & (dx * dx + dy * dy <= r2)
        if hit.any():
            seq[ii] = seq[lo + int(np.flatnonzero(hit)[-1])]
    return seq


def block_bootstrap_mean(vals, blocks, n_boot=N_BOOT, seed=SEED):
    """95% CI for the mean of `vals` by resampling BLOCKS with replacement."""
    vals = np.asarray(vals, dtype=float)
    uniq, inv = np.unique(blocks, return_inverse=True)
    nb = uniq.size
    sums = np.bincount(inv, weights=vals, minlength=nb)
    cnts = np.bincount(inv, minlength=nb).astype(float)
    rng = np.random.default_rng(seed)
    lo_l, hi_l = [], []
    step = 200
    means = np.empty(n_boot)
    for a in range(0, n_boot, step):
        nn = min(step, n_boot - a)
        draws = rng.integers(0, nb, size=(nn, nb))
        means[a:a + nn] = sums[draws].sum(axis=1) / np.maximum(cnts[draws].sum(axis=1), 1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5)), int(nb)


# --------------------------------------------------------------------- scoring pass
def score_pass(qmag_all, nbr, nbr_ok, pool, q_train, raise_bins=0, binw=BIN_W,
               mmax=MMAX_TRUNC):
    """One complete walk-forward pass of the frozen pipeline over a magnitude field.

    `qmag_all` is the whole catalogue's magnitudes as integer bin indices (units binw).
    Returns (gain_bits per scored event, boolean mask of scored pool rows, diagnostics).
    """
    # --- global Mc and b from the training window only
    counts = np.bincount(q_train - int(q_train.min()))
    g_mc_bin = int(np.argmax(counts)) + int(q_train.min()) + int(round(MC_OFFSET / binw)) \
        + raise_bins
    sub = q_train[q_train >= g_mc_bin]
    assert sub.size >= MIN_NSUB, "expected a usable global sample, got %d" % sub.size
    g_den = sub.mean() * binw - (g_mc_bin * binw - binw / 2.0)
    if g_den <= 0:
        return None
    b_glob = 1.0 / (LN10 * g_den)

    # --- local Mc and b on the frozen neighbour sets
    qn = qmag_all[nbr[nbr_ok]]
    mc_bin, b_loc, n_sub, ok_loc = local_mc_b(qn, binw=binw,
                                              offset_bins=int(round(MC_OFFSET / binw)),
                                              raise_bins=raise_bins)
    rows = np.where(nbr_ok)[0]
    qi = qmag_all[pool[rows]]
    scored = ok_loc & (qi >= mc_bin) & (qi >= g_mc_bin)
    if scored.sum() == 0:
        return None
    mb = qi[scored] * binw
    mc = mc_bin[scored] * binw
    lp_loc = gr_bin_logprob(mb, mc, b_loc[scored], binw=binw, mmax=mmax)
    lp_glo = gr_bin_logprob(mb, mc, b_glob, binw=binw, mmax=mmax)
    gain = (lp_loc - lp_glo) / LN2
    assert np.isfinite(gain).all(), "non-finite bits produced"
    est = {"mc_bin": mc_bin[scored], "b_loc": b_loc[scored], "b_glob": float(b_glob),
           "g_mc_bin": int(g_mc_bin)}
    diag = {"b_global": float(b_glob), "mc_global": float(g_mc_bin * binw),
            "n_local_ok": int(ok_loc.sum()), "n_scored": int(scored.sum()),
            "b_local_median": float(np.median(b_loc[scored])),
            "b_local_iqr": [float(np.percentile(b_loc[scored], 25)),
                            float(np.percentile(b_loc[scored], 75))],
            "mc_local_median": float(np.median(mc)),
            "n_sub_median": float(np.median(n_sub[scored]))}
    return gain, rows[scored], diag, est


# ------------------------------------------------------------------ the null shuffle
def shuffle_groups(lat, lon, year, cell_deg=NULL_CELL_DEG):
    gi = (np.floor(lat / cell_deg).astype(np.int64) * 100000
          + np.floor(lon / cell_deg).astype(np.int64) * 10
          + 0)
    key = gi * 10000 + year.astype(np.int64)
    order = np.argsort(key, kind="stable")
    _, starts = np.unique(key[order], return_index=True)
    return order, starts


def shuffled_mag(q, order, starts, rng):
    """Permute magnitudes within each (0.5 deg cell, calendar year) group."""
    out = q.copy()
    vals = q[order]
    ends = np.append(starts[1:], order.size)
    perm = np.empty_like(vals)
    for s, e in zip(starts, ends):
        if e - s > 1:
            perm[s:e] = rng.permutation(vals[s:e])
        else:
            perm[s:e] = vals[s:e]
    out[order] = perm
    return out


# ------------------------------------------------------- sensitivity: b(x) injection
def inject_b_field(q_shuf, lon, floor_bin, delta, rng, binw=BIN_W, mmax=MMAX_TRUNC):
    """Overwrite every magnitude at or above `floor_bin` with a fresh draw from a
    truncated GR whose b varies smoothly with longitude, b(x) = 1.0 +/- delta.

    `floor_bin` deliberately sits INJECT_FLOOR_BELOW under the global Mc so that every
    local sample, above whatever local Mc the estimator picks, is a pure GR with a
    single b. If the floor were placed at the global Mc, local samples whose Mc fell
    below it would straddle a kink -- and a local b would beat a global b at delta = 0
    purely by fitting that kink. Events below the floor keep their shuffled magnitude,
    so the maximum-curvature estimator still reads a realistic completeness roll-off."""
    out = q_shuf.copy()
    g_mc_bin = floor_bin
    sel = np.where(q_shuf >= floor_bin)[0]
    if sel.size == 0:
        return out
    u = 2.0 * (lon[sel] - lon.min()) / (lon.max() - lon.min()) - 1.0
    b = 1.0 + delta * u
    beta = b * LN10
    J = int(round((mmax - g_mc_bin * binw) / binw))
    z = rng.random(sel.size)
    span = -np.expm1(-beta * binw * (J + 1))
    j = np.floor(-np.log1p(-z * span) / (beta * binw)).astype(np.int64)
    j = np.clip(j, 0, J)
    out[sel] = g_mc_bin + j
    return out


# ------------------------------------------------------------- descriptive estimates
def _fixed_n_b(q, blocks, rng, n_fixed=DESC_MIN_N, n_boot=DESC_BOOT, binw=BIN_W):
    """Fixed-n=250 Mc/b point estimate plus a block-bootstrap CI (blocks resampled with
    replacement until 250 events are accumulated). DESCRIPTIVE ONLY."""
    if q.size < n_fixed:
        return None
    pick = rng.choice(q.size, size=n_fixed, replace=False)
    mc = maxc_mc(q[pick] * binw, binw=binw)
    bb, ns = aki_b(q[pick] * binw, mc, binw=binw)
    if not np.isfinite(bb):
        return None
    uniq, inv = np.unique(blocks, return_inverse=True)
    srt = np.argsort(inv, kind="stable")
    bounds = np.searchsorted(inv[srt], np.arange(uniq.size + 1))
    groups = [srt[bounds[g]:bounds[g + 1]] for g in range(uniq.size)]
    boots = []
    for _ in range(n_boot):
        acc, tot, guard = [], 0, 0
        while tot < n_fixed and guard < 5000:
            g = groups[rng.integers(0, len(groups))]
            acc.append(g)
            tot += g.size
            guard += 1
        if tot < n_fixed:
            continue
        sel = np.concatenate(acc)[:n_fixed]
        mcb = maxc_mc(q[sel] * binw, binw=binw)
        v, _ = aki_b(q[sel] * binw, mcb, binw=binw)
        if np.isfinite(v):
            boots.append(v)
    if len(boots) < 20:
        return None
    return {"b": float(bb), "mc": float(mc), "n_above_mc": int(ns),
            "ci": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]}


def descriptive_maps(df_tr, q_tr, seed=SEED):
    """Training-window-only b(x) map (0.2 deg cells) and b vs depth (2 km bins).
    DESCRIPTIVE. No claim is attached to either; they exist so Jim can look."""
    rng = np.random.default_rng(seed)
    lat = df_tr["lat"].to_numpy()
    lon = df_tr["lon"].to_numpy()
    dep = df_tr["depth"].to_numpy()
    t = df_tr["t"].to_numpy()
    tblock = np.floor(t / COARSE_BLOCK_DAYS).astype(np.int64)
    cell_lat = np.floor(lat / MAP_CELL_DEG).astype(np.int64)
    cell_lon = np.floor(lon / MAP_CELL_DEG).astype(np.int64)
    pairs = np.stack([cell_lat, cell_lon], axis=1)
    uniq, inv = np.unique(pairs, axis=0, return_inverse=True)
    inv = np.asarray(inv).ravel()
    srt = np.argsort(inv, kind="stable")
    bnd = np.searchsorted(inv[srt], np.arange(len(uniq) + 1))
    bmap = []
    for g in range(len(uniq)):
        m = srt[bnd[g]:bnd[g + 1]]
        if m.size < DESC_MIN_N:
            continue
        r = _fixed_n_b(q_tr[m], tblock[m], rng)
        if r is None:
            continue
        r.update({"lat": float(uniq[g, 0] * MAP_CELL_DEG + MAP_CELL_DEG / 2),
                  "lon": float(uniq[g, 1] * MAP_CELL_DEG + MAP_CELL_DEG / 2),
                  "n_cell": int(m.size)})
        bmap.append(r)
    cblock = (np.floor(lat / COARSE_BLOCK_DEG).astype(np.int64) * 1000000
              + np.floor(lon / COARSE_BLOCK_DEG).astype(np.int64) * 1000 + tblock)
    dprof = []
    dbin = np.floor(dep / DEPTH_BIN_KM).astype(np.int64)
    dsrt = np.argsort(dbin, kind="stable")
    for k in np.unique(dbin):
        if k < 0 or k * DEPTH_BIN_KM > 40:
            continue
        lo_i, hi_i = np.searchsorted(dbin[dsrt], [k, k + 1])
        m = dsrt[lo_i:hi_i]
        if m.size < DESC_MIN_N:
            continue
        r = _fixed_n_b(q_tr[m], cblock[m], rng)
        if r is None:
            continue
        r.update({"depth_km_lo": float(k * DEPTH_BIN_KM),
                  "depth_km_hi": float((k + 1) * DEPTH_BIN_KM), "n_bin": int(m.size)})
        dprof.append(r)
    return bmap, dprof


# ------------------------------------------------------------------------ one arm
def run_catalogue(fname, verbose=True):
    t_start = time.time()
    df, order = load_catalog(fname)
    mag = df["mag"].to_numpy()
    native_binw = measure_rounding(mag)
    t = df["t"].to_numpy()
    lat = df["lat"].to_numpy()
    lon = df["lon"].to_numpy()
    year = df["year"].to_numpy()

    span = t.max() - t.min()
    cut_expl = t.min() + EXPLORE_FRAC * span
    expl = t < cut_expl
    n_holdout = int((~expl).sum())
    dfe = df[expl].reset_index(drop=True)
    te = dfe["t"].to_numpy()
    cut_train = te.min() + TRAIN_FRAC_OF_EXPLORE * (te.max() - te.min())
    is_train = te < cut_train
    is_score = ~is_train
    assert is_train.sum() > 10000 and is_score.sum() > 10000, \
        "expected populated train/score windows, got %d/%d" % (is_train.sum(),
                                                               is_score.sum())

    q_all = np.round(dfe["mag"].to_numpy() / BIN_W).astype(np.int32)
    q_train = q_all[is_train]
    x, y, lat0, lon0 = project_km(dfe["lat"].to_numpy(), dfe["lon"].to_numpy())
    te = dfe["t"].to_numpy()

    score_idx = np.flatnonzero(is_score)
    rng0 = np.random.default_rng(SEED)
    subsampled = False
    if score_idx.size > SCORING_CAP:
        score_idx = np.sort(rng0.choice(score_idx, size=SCORING_CAP, replace=False))
        subsampled = True
    if verbose:
        print("  %s [%s cols] n=%d  explore=%d (holdout %d reserved)"
              % (fname, order, len(df), int(expl.sum()), n_holdout), flush=True)
        print("  train n=%d  score-window n=%d  scoring pool n=%d%s  native binw=%s"
              % (int(is_train.sum()), int(is_score.sum()), score_idx.size,
                 " (SUBSAMPLED)" if subsampled else "", native_binw), flush=True)

    nbr, nbr_ok, n_redo = neighbour_sets(score_idx, te, x, y, verbose=verbose)
    assert nbr_ok.sum() > 1000, "expected many local estimates, got %d" % nbr_ok.sum()
    if verbose:
        print("  local estimate available for %d / %d (%.1f%%), %d exact redos"
              % (nbr_ok.sum(), nbr_ok.size, 100.0 * nbr_ok.mean(), n_redo), flush=True)

    # ---------------- PRIMARY
    prim = score_pass(q_all, nbr, nbr_ok, score_idx, q_train)
    assert prim is not None, "expected a primary pass, got None"
    gain, rows_scored, diag, est = prim
    assert gain.size > 1000, "expected ~thousands scored, got %d" % gain.size
    ev = score_idx[rows_scored]
    seq = sequence_ids(te[ev], dfe["lat"].to_numpy()[ev], dfe["lon"].to_numpy()[ev],
                       dfe["mag"].to_numpy()[ev])
    lo, hi, n_seq = block_bootstrap_mean(gain, seq)
    coarse = (np.floor(dfe["lat"].to_numpy()[ev] / COARSE_BLOCK_DEG).astype(np.int64)
              * 1000000
              + np.floor(dfe["lon"].to_numpy()[ev] / COARSE_BLOCK_DEG).astype(np.int64)
              * 1000
              + np.floor(te[ev] / COARSE_BLOCK_DAYS).astype(np.int64))
    clo, chi, n_cb = block_bootstrap_mean(gain, coarse)
    primary = {"gain_bits_per_event": float(gain.mean()),
               "ci95_sequence_block": [lo, hi], "n_sequences": n_seq,
               "ci95_coarse_spacetime_block": [clo, chi], "n_coarse_blocks": n_cb,
               "n_scored": int(gain.size), **diag}
    if verbose:
        print("  PRIMARY gain = %+.4f bits/event  CI [%+.4f, %+.4f]  n=%d seq=%d"
              % (gain.mean(), lo, hi, gain.size, n_seq), flush=True)

    # ---------------- control (c): M >= 2.5, and the Merton-requested M >= 3.0
    def _subset(cut):
        mhi = dfe["mag"].to_numpy()[ev] >= cut
        if mhi.sum() < 50:
            return {"n_scored": int(mhi.sum()), "note": "too few events to quote"}
        hlo, hhi, hseq = block_bootstrap_mean(gain[mhi], seq[mhi])
        return {"gain_bits_per_event": float(gain[mhi].mean()),
                "ci95_sequence_block": [hlo, hhi], "n_scored": int(mhi.sum()),
                "n_sequences": hseq}

    high_m = _subset(HIGH_M_CUT)
    high_m3 = _subset(HIGH_M_CUT_2)
    if verbose:
        print("  M>=2.5 subset: %s" % json.dumps(high_m), flush=True)
        print("  M>=3.0 subset: %s" % json.dumps(high_m3), flush=True)

    # ---------------- MERTON FOLLOW-UP 1: is the magnitude distribution exponential?
    mag_tr = dfe["mag"].to_numpy()[is_train]
    lil_bw = native_binw if native_binw else BIN_W   # dither at the MEASURED resolution
    lil_at_floor = lilliefors_exp_p(mag_tr, diag["mc_global"], binw=lil_bw)
    lil_scan = lilliefors_completeness_scan(mag_tr, diag["mc_global"], binw=lil_bw)
    if verbose:
        print("  Lilliefors exponentiality at the global floor M%.1f: p=%s (n=%s); "
              "first non-rejected floor = %s"
              % (diag["mc_global"], lil_at_floor.get("p"), lil_at_floor.get("n"),
                 lil_scan["lilliefors_mc"]), flush=True)

    # ---------------- MERTON FOLLOW-UP 2: b-positive (van der Elst 2021) as the local
    # estimator, everything else identical. b-positive uses only the POSITIVE magnitude
    # differences between successive events, which are immune to a changing Mc.
    rows_all = np.where(nbr_ok)[0]
    qi_all = q_all[score_idx[rows_all]]
    mc_bin_all, _b_all, _ns_all, ok_loc_all = local_mc_b(q_all[nbr[nbr_ok]])
    g_mc_bin = int(round(diag["mc_global"] / BIN_W))
    # sorting the neighbour indices puts the local sample back into TIME order, which is
    # what b-positive needs
    bpos, ns_bp, ok_bp = b_positive_rows(q_all[np.sort(nbr[nbr_ok], axis=1)])
    b_glob_pos, n_gp = b_positive(q_train)
    bp = {"estimator": "van der Elst (2021) b-positive, dMc = %.1f" % BPOS_DMC,
          "b_global_positive": None if not np.isfinite(b_glob_pos) else float(b_glob_pos),
          "n_global_positive_diffs": int(n_gp)}
    sc_bp = ok_loc_all & ok_bp & (qi_all >= mc_bin_all) & (qi_all >= g_mc_bin)
    if np.isfinite(b_glob_pos) and sc_bp.sum() > 500:
        mb = qi_all[sc_bp] * BIN_W
        mcv = mc_bin_all[sc_bp] * BIN_W
        gbp = (gr_bin_logprob(mb, mcv, bpos[sc_bp])
               - gr_bin_logprob(mb, mcv, b_glob_pos)) / LN2
        assert np.isfinite(gbp).all(), "non-finite bits in the b-positive arm"
        evb = score_idx[rows_all[sc_bp]]
        seqb = sequence_ids(te[evb], dfe["lat"].to_numpy()[evb],
                            dfe["lon"].to_numpy()[evb], dfe["mag"].to_numpy()[evb])
        lb, hb, sb = block_bootstrap_mean(gbp, seqb)
        bp.update({"gain_bits_per_event": float(gbp.mean()),
                   "ci95_sequence_block": [lb, hb], "n_scored": int(sc_bp.sum()),
                   "n_sequences": sb,
                   "b_local_positive_median": float(np.median(bpos[sc_bp])),
                   "n_positive_diffs_median": float(np.median(ns_bp[sc_bp]))})
    if verbose:
        print("  b-positive secondary: gain = %s"
              % bp.get("gain_bits_per_event", "n/a"), flush=True)

    # ---------------- control (a): every floor raised by +0.3
    rb = int(round(FLOOR_RAISE / BIN_W))
    pr = score_pass(q_all, nbr, nbr_ok, score_idx, q_train, raise_bins=rb)
    if pr is None:
        raised = {"note": "no scored events after the floor raise"}
    else:
        g2, rows2, d2, _e2 = pr
        ev2 = score_idx[rows2]
        seq2 = sequence_ids(te[ev2], dfe["lat"].to_numpy()[ev2],
                            dfe["lon"].to_numpy()[ev2], dfe["mag"].to_numpy()[ev2])
        l2, h2, s2 = block_bootstrap_mean(g2, seq2)
        raised = {"gain_bits_per_event": float(g2.mean()),
                  "ci95_sequence_block": [l2, h2], "n_scored": int(g2.size),
                  "n_sequences": s2, **d2}
    if verbose:
        print("  +0.3 floor control: gain = %s"
              % raised.get("gain_bits_per_event", "n/a"), flush=True)

    # ---------------- control (b): the likelihood at the MEASURED native bin width
    # The ESTIMATOR stays on 0.1-mag bins (max curvature over 0.01-mag bins on 250
    # events is pure noise); only the likelihood's bin width changes, which is the
    # thing the control is actually about.
    native = {"native_binw": native_binw,
              "note": "same 0.1-mag Mc/b estimator, likelihood bins at the measured "
                      "native magnitude resolution"}
    if native_binw is not None and native_binw < BIN_W:
        m_nat = dfe["mag"].to_numpy()[ev]
        mc_f = est["mc_bin"] * BIN_W
        keep = m_nat >= mc_f - 1e-9
        if keep.sum() > 100:
            lpl = gr_bin_logprob(m_nat[keep], mc_f[keep], est["b_loc"][keep],
                                 binw=native_binw)
            lpg = gr_bin_logprob(m_nat[keep], mc_f[keep], est["b_glob"],
                                 binw=native_binw)
            gnat = (lpl - lpg) / LN2
            assert np.isfinite(gnat).all(), "non-finite bits in the native-bin control"
            ln, hn, sn = block_bootstrap_mean(gnat, seq[keep])
            native.update({"gain_bits_per_event": float(gnat.mean()),
                           "ci95_sequence_block": [ln, hn],
                           "n_scored": int(keep.sum()), "n_sequences": sn,
                           "n_dropped_below_floor_by_rounding": int((~keep).sum())})
    if verbose:
        print("  native-binwidth control: gain = %s"
              % native.get("gain_bits_per_event", "n/a"), flush=True)

    # ---------------- the null
    ordr, starts = shuffle_groups(dfe["lat"].to_numpy(), dfe["lon"].to_numpy(), year[expl])
    rng = np.random.default_rng(SEED + 1)
    null_gains, null_n = [], []
    t_null = time.time()
    for s in range(N_NULL):
        qs = shuffled_mag(q_all, ordr, starts, rng)
        p = score_pass(qs, nbr, nbr_ok, score_idx, qs[is_train])
        if p is None:
            continue
        null_gains.append(float(p[0].mean()))
        null_n.append(int(p[0].size))
        if verbose and (s + 1) % 20 == 0:
            print("    null %3d/%d  mean=%+.4f  %.0fs"
                  % (s + 1, N_NULL, np.mean(null_gains), time.time() - t_null),
                  flush=True)
    null_gains = np.asarray(null_gains)
    assert null_gains.size >= N_NULL // 2, \
        "expected ~%d null replicates, got %d" % (N_NULL, null_gains.size)
    p_real = float((np.sum(null_gains >= gain.mean()) + 1) / (null_gains.size + 1))
    null_block = {"n_shuffles": int(null_gains.size),
                  "mean": float(null_gains.mean()), "sd": float(null_gains.std(ddof=1)),
                  "pct": [float(np.percentile(null_gains, q)) for q in (2.5, 50, 95, 97.5)],
                  "median_n_scored": float(np.median(null_n)),
                  "p_real_gain_one_sided": p_real,
                  "distribution": [round(float(v), 6) for v in null_gains]}
    if verbose:
        print("  NULL mean=%+.5f sd=%.5f  p(real) = %.4f"
              % (null_gains.mean(), null_gains.std(ddof=1), p_real), flush=True)

    # ---------------- sensitivity: injected b(x)
    null95 = float(np.percentile(null_gains, 95))
    rngi = np.random.default_rng(SEED + 2)
    floor_bin = int(round((diag["mc_global"] - INJECT_FLOOR_BELOW) / BIN_W))
    lon_e = dfe["lon"].to_numpy()
    power = []
    all_gains = {}
    t_sens = time.time()
    for delta in INJECT_DELTAS:
        gains = []
        for _ in range(INJECT_REPS):
            qs = shuffled_mag(q_all, ordr, starts, rngi)
            qi = inject_b_field(qs, lon_e, floor_bin, delta, rngi)
            p = score_pass(qi, nbr, nbr_ok, score_idx, qi[is_train])
            if p is None:
                continue
            gains.append(float(p[0].mean()))
        all_gains[delta] = np.asarray(gains)
        if verbose:
            print("    inject delta=%.2f  mean gain=%+.4f  %.0fs"
                  % (delta, np.mean(gains) if gains else float("nan"),
                     time.time() - t_sens), flush=True)
    # Two power definitions, both reported.
    #   frozen   -- the pre-registered success rule: gain > 0.02 bits AND above the
    #               shuffle null's 95th percentile.
    #   matched  -- the weaker, purely statistical question "can this pipeline SEE the
    #               injected field at all": gain above the 95% point of the delta = 0
    #               injection replicates, which is the correct matched reference for the
    #               injection path. Added AFTER the frozen definition returned zero power
    #               at every delta; it changes nothing upstream of this block.
    base = all_gains[INJECT_DELTAS[0]]
    thr_matched = float(base.mean() + 1.645 * base.std(ddof=1)) if base.size > 2 else None
    for delta in INJECT_DELTAS:
        g = all_gains[delta]
        pw_f = float(np.mean((g > SUCCESS_BITS) & (g > null95))) if g.size else 0.0
        pw_m = (float(np.mean(g > thr_matched)) if (g.size and thr_matched is not None)
                else None)
        power.append({"delta_b": delta, "power_frozen_rule": pw_f,
                      "power_matched_reference": pw_m,
                      "mean_gain_bits": float(g.mean()) if g.size else None,
                      "sd_gain_bits": float(g.std(ddof=1)) if g.size > 1 else None,
                      "n_reps": int(g.size),
                      "gains": [round(float(v), 6) for v in g]})

    def _cross(key):
        for a, bq in zip(power[:-1], power[1:]):
            if a[key] is None or bq[key] is None:
                continue
            if a[key] < POWER_TARGET <= bq[key]:
                f = (POWER_TARGET - a[key]) / max(bq[key] - a[key], 1e-9)
                return float(a["delta_b"] + f * (bq["delta_b"] - a["delta_b"]))
        if power and power[0][key] is not None and power[0][key] >= POWER_TARGET:
            return float(power[0]["delta_b"])
        return None

    mdd = _cross("power_frozen_rule")
    mdd_matched = _cross("power_matched_reference")
    if verbose:
        print("    min delta_b at %.0f%% power: frozen rule = %s, matched ref = %s"
              % (100 * POWER_TARGET, mdd, mdd_matched), flush=True)

    # ---------------- descriptive (training window only)
    df_tr = dfe[is_train].reset_index(drop=True)
    bmap, dprof = descriptive_maps(df_tr, q_train)
    if verbose:
        print("  descriptive: %d map cells, %d depth bins  (%.0fs total)"
              % (len(bmap), len(dprof), time.time() - t_start), flush=True)

    passed = (primary["gain_bits_per_event"] > SUCCESS_BITS
              and primary["ci95_sequence_block"][0] > 0.0)
    return {
        "catalogue": fname, "column_order": order, "n_rows": int(len(df)),
        "n_exploration": int(expl.sum()), "n_holdout_reserved_unread": n_holdout,
        "span_days": float(span),
        "n_train": int(is_train.sum()), "n_score_window": int(is_score.sum()),
        "scoring_pool": int(score_idx.size), "scoring_pool_subsampled": subsampled,
        "n_with_local_estimate": int(nbr_ok.sum()),
        "neighbour_exact_redos": int(n_redo),
        "projection_origin": {"lat0": lat0, "lon0": lon0},
        "primary": primary,
        "control_a_floor_plus_0.3": raised,
        "control_b_native_binwidth": native,
        "control_c_M_ge_2.5": high_m,
        "merton_M_ge_3.0": high_m3,
        "merton_lilliefors_exponentiality": {
            "at_global_floor": lil_at_floor,
            "floor_scan": lil_scan,
            "dither_binw": lil_bw,
            "reference": "Herrmann & Marzocchi (2021): QTM and SCSN SoCal are only "
                         "Lilliefors-complete near M 3.24, i.e. magnitudes are NOT "
                         "exponential down to a maximum-curvature floor.",
        },
        "merton_b_positive_secondary": bp,
        "null": null_block,
        "sensitivity": {"power_curve": power,
                        "min_delta_b_at_80pct_power_frozen_rule": mdd,
                        "min_delta_b_at_80pct_power_matched_reference": mdd_matched,
                        "matched_reference_threshold_bits": thr_matched,
                        "null_95pct_bits": null95,
                        "injection_floor_mag": float(floor_bin * BIN_W),
                        "note": "b(x) = 1.0 +/- delta linear in longitude, injected "
                                "above (global Mc - 0.5) into a magnitude-shuffled "
                                "catalogue; identical pipeline. delta=0 is the "
                                "false-positive check."},
        "descriptive_b_map_0.2deg_fixed_n250": bmap,
        "descriptive_b_vs_depth_2km_fixed_n250": dprof,
        "defensible_floor_statement": (
            "The relative log-score comparison (local b vs global b, identical support) "
            "is fair at any floor. The b-VALUE interpretation is defensible only above "
            "the Lilliefors floor measured here (%s); at the maximum-curvature floor "
            "(M %.1f) the magnitude distribution of this catalogue is not exponential, "
            "so 'b' is a two-parameter summary of a non-GR distribution, not a "
            "Gutenberg-Richter slope."
            % (lil_scan["lilliefors_mc"], diag["mc_global"])),
        "success_rule_met": bool(passed),
        "runtime_sec": round(time.time() - t_start, 1),
    }


def main():
    t0 = time.time()
    out = {
        "experiment": "exp_bvalue_skill",
        "hypothesis": "K-405 / PLAYBOOK P-2.2 -- magnitude is the responder: does a "
                      "LOCAL b estimated only from the past forecast the size "
                      "distribution of the next events better than one global b?",
        "state_class": "first-run, exploration split",
        "run_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "frozen_constants": {
            "EXPLORE_FRAC": EXPLORE_FRAC,
            "TRAIN_FRAC_OF_EXPLORE": TRAIN_FRAC_OF_EXPLORE,
            "K_LOCAL": K_LOCAL, "RADIUS_KM": RADIUS_KM,
            "LOOKBACK_DAYS": LOOKBACK_DAYS, "BIN_W": BIN_W, "MC_OFFSET": MC_OFFSET,
            "MMAX_TRUNC": MMAX_TRUNC, "MIN_NSUB": MIN_NSUB, "B_BOUNDS": list(B_BOUNDS),
            "FLOOR_RAISE": FLOOR_RAISE, "HIGH_M_CUT": HIGH_M_CUT,
            "NULL_CELL_DEG": NULL_CELL_DEG, "N_NULL": N_NULL, "N_BOOT": N_BOOT,
            "SEQ_LINK_KM": SEQ_LINK_KM, "SEQ_LINK_DAYS": SEQ_LINK_DAYS,
            "COARSE_BLOCK_DEG": COARSE_BLOCK_DEG, "COARSE_BLOCK_DAYS": COARSE_BLOCK_DAYS,
            "SCORING_CAP": SCORING_CAP, "SEED": SEED, "SUCCESS_BITS": SUCCESS_BITS,
            "INJECT_DELTAS": INJECT_DELTAS, "INJECT_REPS": INJECT_REPS,
            "POWER_TARGET": POWER_TARGET,
            "INJECT_FLOOR_BELOW": INJECT_FLOOR_BELOW, "MAP_CELL_DEG": MAP_CELL_DEG,
            "DEPTH_BIN_KM": DEPTH_BIN_KM, "DESC_MIN_N": DESC_MIN_N,
            "DESC_BOOT": DESC_BOOT,
            "HIGH_M_CUT_2": HIGH_M_CUT_2, "LILLIE_N": LILLIE_N,
            "LILLIE_SIM": LILLIE_SIM, "BPOS_DMC": BPOS_DMC,
        },
        "success_rule": "gain > 0.02 bits/event with the 95%% sequence-block CI "
                        "excluding 0, on BOTH catalogues",
        "named_observer_artifact": "spatial/temporal Mc variation masquerading as b "
                                   "variation; controls (a) +0.3 floor raise, "
                                   "(b) measured native magnitude binning, "
                                   "(c) M>=2.5 subset above the diurnal-detection band",
        "catalogues": {},
    }
    for fname in ("QTM_12dev.txt", "SCSN_original_catalog.txt"):
        print("=== %s ===" % fname, flush=True)
        out["catalogues"][fname] = run_catalogue(fname)
    both = all(v["success_rule_met"] for v in out["catalogues"].values())
    out["verdict"] = {
        "success_rule_met_on_both": bool(both),
        "per_catalogue": {k: {"gain": v["primary"]["gain_bits_per_event"],
                              "ci95": v["primary"]["ci95_sequence_block"],
                              "null_p": v["null"]["p_real_gain_one_sided"],
                              "met": v["success_rule_met"]}
                          for k, v in out["catalogues"].items()},
    }
    out["total_runtime_sec"] = round(time.time() - t0, 1)
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print("wrote %s  (%.0fs)" % (OUT_JSON.name, out["total_runtime_sec"]), flush=True)


if __name__ == "__main__":
    main()
