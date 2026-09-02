"""K-038 ARM A -- does a continuous "dynamic-stress weather" covariate add FORECAST SKILL
over a per-region temporal ETAS, out of sample, on the 14 K-034 fluid/geothermal regions?

K-034 licensed the pipeline as a detector of remote dynamic triggering (certified floor
~34 kPa, suggestive ~10 kPa). This arm converts that positive control into a forecasting
claim: incremental bits/event of [ETAS x exp(beta*g(D))] over [ETAS], walk-forward, on the
EXPLORATION split only.

================================================================================
FROZEN BEFORE THE RUN (PLAYBOOK rule 7 requires the artifacts be named here first)
================================================================================

AMPLITUDE AXIS -- reused VERBATIM from K034_SEALED_LITERATURE.md, not refit:
  van der Elst & Brodsky (2010) eq. (6): log10 A20[um] = M - 1.66 log10(D[deg]) - 2,
  V = 2*pi*A20/T, T = 20 s, sigma = G*V/c, G = 30 GPa, c = 3.5 km/s.
  (K-034 used Ms; a global ComCat catalogue has no Ms column, so ComCat `mag` is used --
  this is exactly K-034's carried "Mw-substituted axis". DECLARED DEVIATION D1.)

COVARIATE (frozen):  D_r(t) = sum_k sigma(M_k, r_k) * exp(-(t - t_arr_k)/tau), t >= t_arr_k
  t_arr_k = t_k + r_k / 3.5 km/s   (surface-wave arrival, K-034's c)
  LINK      g(D) = log1p(D / D0),  D0 = 10 kPa   (K-034's suggestive floor)
  TAU FAMILY (declared, size 2, max-statistic paid over it): tau in {1 d, 5 d}

TRIGGERS: ComCat worldwide M >= 5.5, 1985-01-01 .. 2023-01-01, one FDSN page per year.
  EXCLUSIONS (per target region): trigger closer than 300 km to the target BOX;
  trigger closer than 2 rupture lengths (Wells & Coppersmith 1994 subsurface rupture
  length, all slip types: log10 L[km] = -2.44 + 0.59 M -- DECLARED DEVIATION D2, this
  coefficient pair is from memory and is UNVERIFIED against the paper in session; it sets
  only an exclusion gate, never the statistic, and it reproduces K-034's four hand-entered
  lengths as 74/56/56/166 km against their 70/48/50/340 km).
  GLOBAL EXCLUSION: any trigger inside ANY of the 14 K-034 target boxes is dropped.

MODELS (per region, M0 = 1.5, the catalogue floor):
  A: lambda_A(t) = mu_{epoch(t)} + K * sum_{t_i<t} 10^(alpha (M_i - M0)) (t - t_i + c)^-p
  B: lambda_B(t) = lambda_A(t) * exp(beta * g(D_r(t)))
  Epoch mu is piecewise-constant on 5-year epochs from 1985-01-01, in BOTH models, IF the
  frozen completeness rule fires: FROZEN RULE -- if in any region max-over-epochs minus
  min-over-epochs of the maximum-curvature Mc is >= 0.3, epoch mu is used everywhere.
  Scoring-window epochs inherit the LAST TRAINED epoch's mu (they are not estimable from
  training; both models inherit identically, so the increment is unaffected to first order).

SPLIT: exploration = first 70% of each region's own time span (engine/splits.py rule).
  training = first 60% of exploration; scoring = remaining 40% of exploration, walk-forward
  with theta and beta frozen from training. THE LAST 30% OF THE SPAN IS NEVER READ.

STATISTIC (frozen, primary): pooled incremental bits/event of B over A on the scoring
  window, maximised over the tau family. Secondary: per region; by K-034 class (A/B vs C).
  SUCCESS (Popper, K-038 verdict): >= 0.01 bits/event OOS pooled AND p < 0.05 under a
  circular-shift null (whole trigger catalogue shifted by 200 random offsets >= 1 yr,
  wrapping in the trigger span; beta refit each time; max-statistic over the tau family).
  CI: block bootstrap over sequences (25 km / 7 d local declustering, per exp_fluid_driven).

================================================================================
ARTIFACTS THAT COULD FAKE THIS RESULT -- NAMED BEFORE THE RUN (PLAYBOOK rule 7)
================================================================================
A1. CODA BLINDNESS. Large teleseismic surface waves saturate local records, so SMALL local
    events are MISSED for minutes-to-hours right after the arrival. This biases AGAINST a
    positive beta, so it cannot manufacture skill -- but it can hide it, and it can invert
    the sign at short tau. READING REPORTED: counts of M<2.0 vs M>=2.0 in the 6 h after the
    20 largest predicted amplitudes, against the matched 6 h windows one day earlier.
A2. IN-BOX AFTERSHOCK CONTAMINATION. A trigger inside (or beside) the target box makes its
    own aftershocks look like "remote triggering". HANDLED BY CONSTRUCTION: the 300 km box
    gate, the 2-rupture-length gate, and the global K-034-box exclusion.
A3. COMPLETENESS DRIFT 1985->2022. Network growth lowers Mc and raises the apparent rate;
    if that drift happened to correlate with global M>=5.5 activity it would masquerade as
    skill. HANDLED: 5-year epoch mu in both models (frozen rule above), Mc measured and
    reported per epoch per region.
A4. DAY/NIGHT DETECTION ARTIFACT. Cultural noise makes Mc diurnal. IRRELEVANT HERE and
    stated so: the covariate decays on tau = 1-5 d and trigger origin times are uniformly
    distributed over the solar day, so a diurnal detection modulation is orthogonal to
    D_r(t) at day scales. No correction applied; no correction needed.
A5. (added) The ETAS triggering term can ABSORB dynamic triggering -- a remotely triggered
    burst has its own aftershocks, which model A explains. This biases AGAINST B. Noted.

================================================================================
ADDED AFTER THE FIRST RUN -- FLAGGED, EXPLORATORY, SETS NO FLAG (K-034 D3 precedent)
================================================================================
X1. TYPE-I ERROR AT beta = 0 (`false_positive_rate_beta0`). Added on prior art supplied by
    Merton: Hardebeck, DeSalvio, Fan & Barbour (2025), JGR, doi:10.1029/2025JB031566, show
    standard remote-triggering tests run 3.5-8.5% false positives on NO-triggering synthetic
    catalogues because the null does not preserve clustering. Measured here by simulating
    no-triggering catalogues from each region's fitted ETAS intensity path and running the
    identical pipeline, circular-shift null and tau family included.
X2. SINGLE POOLED beta (`pooled_beta_exploratory`). The pre-registered model spends one beta
    per region; with 14 weakly-informative training windows the out-of-sample penalty for a
    noise-fitted parameter is large. This variant shares one beta across all 14 regions. It
    runs on the coarse Monte-Carlo grid and 40 shifts, is EXPLORATORY, and sets no flag.
NOTE on the Monte-Carlo grid: the primary statistic uses 0.25-day quadrature cells; the
    Monte-Carlo arms (X1, X2, sensitivity) use 1-day cells for tractability, which slightly
    smears the tau = 1 d covariate. Declared, not corrected.
PRIOR ART, recorded so nothing here reads as a new phenomenon: continuous-field framings of
    dynamic triggering are already published (DeSalvio & Fan on QTM; Miyazawa; Brodsky; Guo).
    The only contribution here is the incremental-bits-per-event forecast metric on an
    out-of-sample split.

Run: python -u exp_dyntrig_skill.py            (downloads triggers on first use)
"""
from __future__ import annotations

import io
import json
import os
import ssl
import time
import hashlib
import urllib.request
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
DK034 = os.path.join(HERE, "data", "k034")
DGLOB = os.path.join(HERE, "data", "global_m55")
OUT = os.path.join(HERE, "results_dyntrig_skill.json")

# ----------------------------- FROZEN CONSTANTS -----------------------------
G_SHEAR = 30e9          # Pa                        (K034_SEALED_LITERATURE.md)
C_PHASE = 3500.0        # m/s
T_SW = 20.0             # s
D0_PA = 10e3            # Pa, link scale = K-034 suggestive floor
TAUS = (1.0, 5.0)       # days -- DECLARED FAMILY OF 2
EXCL_BOX_KM = 300.0
RUPT_A, RUPT_B = -2.44, 0.59      # Wells & Coppersmith 1994 RLD, all slip types
TRIG_MINMAG = 5.5
TRIG_T0, TRIG_T1 = "1985-01-01", "2023-01-01"
M0 = 1.5
EXPLORE_FRAC = 0.70
TRAIN_FRAC_OF_EXPLORE = 0.60
EPOCH_YEARS = 5.0
EPOCH_ORIGIN = "1985-01-01"
MC_DRIFT_TRIGGER = 0.30
GRID_DT_D = 0.25        # covariate quadrature cell, days
N_SHIFT = 200
SHIFT_MIN_YR = 1.0
N_BOOT = 2000
# Type-I-error arm (Hardebeck, DeSalvio, Fan & Barbour 2025, JGR 10.1029/2025JB031566:
# standard remote-triggering tests run 3.5-8.5% false positives on no-triggering synthetics
# because the null does not preserve clustering). Scaled to the wall-clock budget.
N_FPR_SIM = 60
N_FPR_SHIFT = 40
FPR_CELL_DAYS = 1.0
FPR_BINS = 60
DECLUSTER_KM_LOCAL = 25.0
DECLUSTER_DAYS_LOCAL = 7.0
TRUNC_W_FIT_D = 1000.0
BETA_BOUNDS = (-5.0, 5.0)
SEED = 20260902
LN2 = np.log(2.0)
LN10 = np.log(10.0)
N_PROC = 8

CLASS = {'long_valley': 'A', 'coso': 'A', 'geysers': 'A', 'yellowstone': 'A',
         'salton_brawley': 'A', 'lassen': 'B', 'mono_west_nv_mina': 'B',
         'little_skull_mtn': 'B', 'cedar_city_ut': 'B', 'smith_valley_nv': 'C',
         'parkfield': 'C', 'mendocino': 'C', 'wasatch_slc': 'C', 'san_jacinto': 'C'}

LANDERS = dict(t0="1992-06-28T11:57:34Z", lat=34.200, lon=-116.436, M=7.3)
POSCTRL_CELLS = ("long_valley", "yellowstone", "geysers")


# ----------------------------- geometry / amplitude -----------------------------
def gcdist_km(la1, lo1, la2, lo2):
    p = np.pi / 180.0
    la1, lo1, la2, lo2 = (np.asarray(x, float) for x in (la1, lo1, la2, lo2))
    a = (np.sin((la2 - la1) * p / 2) ** 2 +
         np.cos(la1 * p) * np.cos(la2 * p) * np.sin((lo2 - lo1) * p / 2) ** 2)
    return 2 * 6371.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def dist_to_box_km(lat, lon, box):
    """Great-circle distance to the nearest point of box [lat0,lat1,lon0,lon1]; 0 inside."""
    la0, la1, lo0, lo1 = box
    return gcdist_km(lat, lon, np.clip(lat, la0, la1), np.clip(lon, lo0, lo1))


def in_box(lat, lon, box):
    la0, la1, lo0, lo1 = box
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    return (lat >= la0) & (lat <= la1) & (lon >= lo0) & (lon <= lo1)


def sigma_primary_pa(M, r_km):
    """vdE&B (2010) eq. (6) far-field surface-wave peak dynamic stress. FROZEN, not refit."""
    Ddeg = np.asarray(r_km, float) / 111.195
    A20_um = 10.0 ** (np.asarray(M, float) - 1.66 * np.log10(Ddeg) - 2.0)
    V = 2 * np.pi * (A20_um * 1e-6) / T_SW
    return G_SHEAR * V / C_PHASE


def rupture_len_km(M):
    return 10.0 ** (RUPT_A + RUPT_B * np.asarray(M, float))


def glink(D_pa):
    return np.log1p(np.asarray(D_pa, float) / D0_PA)


# ----------------------------- covariate construction -----------------------------
def decay_series_ref(eval_t, trig_t, trig_amp, tau):
    """Reference O(N*M) implementation of D(t) = sum_k amp_k exp(-(t-trig_t_k)/tau),
    over trig_t_k <= t. Used by the tests to certify the fast path."""
    eval_t = np.asarray(eval_t, float)
    trig_t = np.asarray(trig_t, float)
    trig_amp = np.asarray(trig_amp, float)
    out = np.zeros(eval_t.size)
    for i, t in enumerate(eval_t):
        m = trig_t <= t
        if m.any():
            out[i] = float(np.sum(trig_amp[m] * np.exp(-(t - trig_t[m]) / tau)))
    return out


def decay_series(eval_t, trig_t, trig_amp, tau):
    """Fast exact D(t). eval_t and trig_t must be sorted ascending.

    Blocked cumulative-sum form: within a block of width 20*tau the exponentials are
    representable, and a scalar `carry` transports the decayed state across blocks."""
    eval_t = np.asarray(eval_t, float)
    trig_t = np.asarray(trig_t, float)
    trig_amp = np.asarray(trig_amp, float)
    out = np.zeros(eval_t.size)
    if eval_t.size == 0 or trig_t.size == 0:
        return out
    blk = 20.0 * tau
    t_lo = float(min(eval_t[0], trig_t[0])) - 1e-9
    t_hi = float(max(eval_t[-1], trig_t[-1])) + blk
    edges = np.arange(t_lo, t_hi + blk, blk)
    idx_e = np.clip(np.searchsorted(edges, eval_t, side="right") - 1, 0, len(edges) - 2)
    carry = 0.0
    ti = 0
    for b in range(len(edges) - 1):
        e0, e1 = float(edges[b]), float(edges[b + 1])
        sel = np.nonzero(idx_e == b)[0]
        tj = int(np.searchsorted(trig_t, e1, side="left"))
        tt = trig_t[ti:tj]
        aa = trig_amp[ti:tj]
        if sel.size:
            te = eval_t[sel]
            dec = np.exp(-(te - e0) / tau)
            base = carry * dec
            if tt.size:
                cw = np.cumsum(aa * np.exp((tt - e0) / tau))
                k = np.searchsorted(tt, te, side="right")
                add = np.where(k > 0, cw[np.clip(k - 1, 0, cw.size - 1)], 0.0)
                base = base + add * dec
            out[sel] = base
        carry = carry * np.exp(-(e1 - e0) / tau)
        if tt.size:
            carry += float(np.sum(aa * np.exp(-(e1 - tt) / tau)))
        ti = tj
    return out


# ----------------------------- ETAS with epoch mu -----------------------------
def _pair_sums(t, w, m, c, p, lo, hi, W, want_grad, chunk=512):
    n_t = hi - lo
    S = np.zeros(n_t)
    Sa = np.zeros(n_t) if want_grad else None
    Sc = np.zeros(n_t) if want_grad else None
    Sp = np.zeros(n_t) if want_grad else None
    finite_W = np.isfinite(W)
    for x0 in range(lo, hi, chunk):
        x1 = min(x0 + chunk, hi)
        hlo = int(np.searchsorted(t, t[x0] - W, side="left")) if finite_W else 0
        if hlo >= x1:
            continue
        th = t[hlo:x1]
        dt = t[x0:x1, None] - th[None, :]
        valid = dt > 0
        if finite_W:
            valid &= dt <= W
        dtc = np.where(valid, dt + c, 1.0)
        ldtc = np.log(dtc)
        q = np.exp(-p * ldtc) * valid * w[hlo:x1][None, :]
        o = slice(x0 - lo, x1 - lo)
        S[o] = q.sum(axis=1)
        if want_grad:
            Sa[o] = (q * m[hlo:x1][None, :]).sum(axis=1)
            Sc[o] = -p * (q / dtc).sum(axis=1)
            Sp[o] = -(q * ldtc).sum(axis=1)
    return S, Sa, Sc, Sp


def _G_terms(t_src, T0, T1, c, p, W, want_grad):
    lo = np.maximum(T0, t_src)
    hi = np.minimum(T1, t_src + W) if np.isfinite(W) else np.full_like(t_src, T1)
    live = hi > lo
    A = np.where(live, hi - t_src + c, 1.0)
    B = np.where(live, lo - t_src + c, 1.0)
    lA, lB = np.log(A), np.log(B)
    s = 1.0 - p
    small = abs(s) < 1e-6
    if small:
        G = (lA - lB) + s * (lA ** 2 - lB ** 2) / 2.0 + s * s * (lA ** 3 - lB ** 3) / 6.0
    else:
        As, Bs = np.exp(s * lA), np.exp(s * lB)
        G = (As - Bs) / s
    G = G * live
    if not want_grad:
        return G, None, None
    dGdc = (np.exp(-p * lA) - np.exp(-p * lB)) * live
    if small:
        dGds = (lA ** 2 - lB ** 2) / 2.0 + s * (lA ** 3 - lB ** 3) / 3.0
    else:
        As, Bs = np.exp(s * lA), np.exp(s * lB)
        dGds = ((lA * As - lB * Bs) * s - (As - Bs)) / (s * s)
    return G, dGdc, -dGds * live


def etas_ll(theta, t_all, m_all, lo, hi, T0, T1, W, ep_tgt, ep_len, want_grad=True, chunk=512):
    """theta = [mu_0..mu_{E-1}, K, alpha, c, p]. ep_tgt = epoch index per TARGET event
    (length hi-lo); ep_len = days of [T0,T1) falling in each epoch."""
    E = len(ep_len)
    mu = np.asarray(theta[:E], float)
    K, alpha, c, p = (float(v) for v in theta[E:])
    w = np.exp(alpha * LN10 * m_all)
    S, Sa, Sc, Sp = _pair_sums(t_all, w, m_all, c, p, lo, hi, W, want_grad, chunk)
    lam = mu[ep_tgt] + K * S
    if not np.all(np.isfinite(lam)) or np.any(lam <= 0):
        return -1e18, (np.zeros(E + 4) if want_grad else None)
    src = t_all <= T1
    G, dGdc, dGdp = _G_terms(t_all[src], T0, T1, c, p, W, want_grad)
    wG = w[src] * G
    Lam = float(mu @ ep_len) + K * float(wG.sum())
    LL = float(np.log(lam).sum() - Lam)
    if not want_grad:
        return LL, None
    inv = 1.0 / lam
    g_mu = np.bincount(ep_tgt, weights=inv, minlength=E) - ep_len
    g_K = float((S * inv).sum() - wG.sum())
    g_a = K * float((Sa * inv).sum()) - K * float((m_all[src] * wG).sum())
    g_c = K * float((Sc * inv).sum()) - K * float((w[src] * dGdc).sum())
    g_p = K * float((Sp * inv).sum()) - K * float((w[src] * dGdp).sum())
    return LL, np.concatenate([g_mu, [g_K, g_a * LN10, g_c, g_p]])


def _obj(x, *args):
    theta = np.exp(x)
    LL, g = etas_ll(theta, *args)
    if not np.isfinite(LL):
        return 1e18, np.zeros_like(x)
    return -LL, -g * theta


# ----------------------------- beta likelihood -----------------------------
def beta_profile(gj, Acell, gcell, beta):
    """log L(beta) - log L_A, up to the terms that do not involve beta.

    Model B multiplies lambda_A by exp(beta*g). On the scoring/training window:
      LL_B - LL_A = beta*sum_j g_j - [ sum_k A_k exp(beta g_k) - sum_k A_k ]
    where A_k = integral of lambda_A over quadrature cell k (computed analytically, so
    lambda_A's Omori spikes are exact) and g_k is the smooth covariate at the cell midpoint.
    """
    beta = float(beta)
    return float(beta * gj.sum() - (np.sum(Acell * np.exp(beta * gcell)) - Acell.sum()))


def fit_beta(gj, Acell, gcell):
    r = minimize_scalar(lambda b: -beta_profile(gj, Acell, gcell, b),
                        bounds=BETA_BOUNDS, method="bounded",
                        options={"xatol": 1e-6})
    return float(r.x), float(-r.fun)


# ----------------------------- misc -----------------------------
def max_curvature_mc(mags, binw=0.1):
    mags = np.asarray(mags, float)
    if mags.size < 200:
        return None
    lo = np.floor(mags.min() * 10) / 10 - binw / 2
    edges = np.arange(lo, mags.max() + binw, binw)
    h, _ = np.histogram(mags, bins=edges)
    i = int(np.argmax(h))
    return float(round((edges[i] + edges[i + 1]) / 2, 2))


def sequence_ids(t, lat, lon, mag, days=DECLUSTER_DAYS_LOCAL, km=DECLUSTER_KM_LOCAL):
    """Cluster IDs by the program's local rule (exp_fluid_driven._sequence_ids, local radii)."""
    n = t.size
    seq = np.arange(n)
    recent = []
    for ii in np.argsort(t):
        while recent and t[ii] - t[recent[0]] > days:
            recent.pop(0)
        for jj in reversed(recent):
            if mag[jj] >= mag[ii] and gcdist_km(lat[jj], lon[jj], lat[ii], lon[ii]) <= km:
                seq[ii] = seq[jj]
                break
        recent.append(ii)
    return seq


# ----------------------------- download -----------------------------
def download_triggers(verbose=True):
    os.makedirs(DGLOB, exist_ok=True)
    csv_p = os.path.join(DGLOB, "global_m55.csv")
    log_p = os.path.join(DGLOB, "download_log.jsonl")
    man_p = os.path.join(DGLOB, "manifest.json")
    if os.path.exists(csv_p) and os.path.exists(man_p):
        if verbose:
            print("[dl] cached " + csv_p, flush=True)
        return csv_p
    try:
        import certifi
        ctx = ssl.create_default_context(cafile=certifi.where())
    except Exception:                                    # noqa: BLE001
        ctx = ssl.create_default_context()
    url = "https://earthquake.usgs.gov/fdsnws/event/1/query"
    frames, logrecs = [], []
    for y in range(1985, 2023):
        t0, t1 = "%d-01-01" % y, "%d-01-01" % (y + 1)
        q = (url + "?format=csv&starttime=" + t0 + "&endtime=" + t1 +
             "&minmagnitude=" + str(TRIG_MINMAG) + "&orderby=time-asc&limit=20000")
        body, status, err = None, None, None
        for attempt in range(1, 6):
            try:
                with urllib.request.urlopen(q, context=ctx, timeout=300) as r:
                    status = int(r.status)
                    body = r.read().decode("utf-8")
                break
            except Exception as e:                       # noqa: BLE001
                err = "%s: %s" % (type(e).__name__, e)
                time.sleep(4 * attempt)
        if body is None:
            raise RuntimeError("FDSN page %s->%s FAILED after 5 attempts: %s" % (t0, t1, err))
        if status != 200:
            raise RuntimeError("FDSN page %s->%s HTTP %s" % (t0, t1, status))
        df = pd.read_csv(io.StringIO(body)) if body.strip() else pd.DataFrame()
        n = len(df)
        # FAILURE-FIRST: worldwide M>=5.5 runs ~100-400 events/yr; 0 rows is never legitimate.
        if n == 0:
            raise RuntimeError("FDSN page %s->%s returned 0 rows; expected >0" % (t0, t1))
        if n >= 20000:
            raise RuntimeError("FDSN page %s->%s hit the 20000 cap (n=%d)" % (t0, t1, n))
        logrecs.append(dict(page=t0 + "->" + t1, url=q, http_status=status, rows=int(n),
                            ts_utc=pd.Timestamp.utcnow().isoformat()))
        frames.append(df)
        if verbose:
            print("[dl] %s -> %s  HTTP %s  rows=%d" % (t0, t1, status, n), flush=True)
    cat = pd.concat(frames, ignore_index=True)
    n_pre = len(cat)
    cat = cat.drop_duplicates(subset=["id"]).sort_values("time").reset_index(drop=True)
    if len(cat) < 5000:
        raise RuntimeError("global M>=5.5 1985-2023: expected ~8000-14000 rows, got %d" % len(cat))
    cat.to_csv(csv_p, index=False)
    with open(log_p, "w", encoding="utf-8") as fh:
        for r in logrecs:
            fh.write(json.dumps(r, sort_keys=True) + "\n")
    h = hashlib.sha256(open(csv_p, "rb").read()).hexdigest()
    json.dump(dict(source="USGS ComCat FDSN event query", minmag=TRIG_MINMAG,
                   t0=TRIG_T0, t1=TRIG_T1, pages=len(logrecs),
                   rows_before_dedup=int(n_pre), rows=int(len(cat)), sha256=h,
                   first=str(cat["time"].iloc[0]), last=str(cat["time"].iloc[-1]),
                   downloaded_utc=pd.Timestamp.utcnow().isoformat()),
              open(man_p, "w"), indent=1)
    print("[dl] TOTAL rows=%d (pre-dedup %d) sha256=%s" % (len(cat), n_pre, h[:16]), flush=True)
    return csv_p


# ----------------------------- region assembly -----------------------------
DAY = 86400.0


def load_region(name):
    df = pd.read_csv(os.path.join(DK034, name + ".csv"),
                     usecols=["time", "latitude", "longitude", "mag"])
    t = pd.to_datetime(df["time"], utc=True, format="mixed")
    df = df.assign(t=(t - pd.Timestamp(0, tz="UTC")).dt.total_seconds().to_numpy() / DAY)
    df = df.dropna(subset=["t", "mag", "latitude", "longitude"]).sort_values("t")
    df = df[df["mag"] >= M0 - 1e-9].reset_index(drop=True)
    if len(df) < 200:
        raise RuntimeError("region %s: expected >=200 events, got %d" % (name, len(df)))
    return df


def epoch_edges_days(t_lo, t_hi):
    org = (pd.Timestamp(EPOCH_ORIGIN, tz="UTC") - pd.Timestamp(0, tz="UTC")).total_seconds() / DAY
    step = EPOCH_YEARS * 365.25
    k0 = int(np.floor((t_lo - org) / step))
    k1 = int(np.floor((t_hi - org) / step))
    return org + step * np.arange(k0, k1 + 2)


def epoch_lengths(edges, T0, T1):
    lo = np.clip(edges[:-1], T0, T1)
    hi = np.clip(edges[1:], T0, T1)
    return np.maximum(hi - lo, 0.0)


def build_trigger_set(trig, box, clat, clon):
    """Apply the frozen exclusions and return (t_arrival_days, amplitude_Pa) sorted."""
    d_box = dist_to_box_km(trig["latitude"].to_numpy(), trig["longitude"].to_numpy(), box)
    r_cen = gcdist_km(trig["latitude"].to_numpy(), trig["longitude"].to_numpy(), clat, clon)
    twoL = 2.0 * rupture_len_km(trig["mag"].to_numpy())
    keep = (d_box >= EXCL_BOX_KM) & (d_box >= twoL) & (~trig["in_any_k034_box"].to_numpy())
    sub = trig[keep]
    r = r_cen[keep]
    amp = sigma_primary_pa(sub["mag"].to_numpy(), r)
    t_arr = sub["t"].to_numpy() + (r / 3.5) / DAY          # r/3.5 km/s in days
    o = np.argsort(t_arr)
    return t_arr[o], amp[o], r[o], sub["mag"].to_numpy()[o], int(keep.sum()), int(len(trig))


def cell_integrals(edges, t_src, w_src, c, p, mu_vec, ep_edges, chunk=64):
    """A_k = integral of lambda_A over cell k, EXACT for the Omori part.

    F(t) = K*sum_i w_i H(t - t_i), H(s) = ((s+c)^(1-p) - c^(1-p))/(1-p), H(s<=0)=0.
    A_k = mu_{epoch(k)} * dt_k + F(edge_{k+1}) - F(edge_k).   (K folded into w_src)
    """
    s = 1.0 - p
    cs = c ** s if abs(s) > 1e-12 else np.log(c)
    F = np.zeros(edges.size)
    for x0 in range(0, edges.size, chunk):
        x1 = min(x0 + chunk, edges.size)
        d = edges[x0:x1, None] - t_src[None, :]
        pos = d > 0
        dc = np.where(pos, d + c, c)
        if abs(s) > 1e-12:
            H = (dc ** s - cs) / s
        else:
            H = np.log(dc) - cs
        F[x0:x1] = (np.where(pos, H, 0.0) * w_src[None, :]).sum(axis=1)
    dt = np.diff(edges)
    mid = 0.5 * (edges[:-1] + edges[1:])
    ei = np.clip(np.searchsorted(ep_edges, mid, side="right") - 1, 0, len(mu_vec) - 1)
    A = mu_vec[ei] * dt + np.diff(F)
    return np.maximum(A, 1e-300)


def make_edges(T0, T1, dt=GRID_DT_D):
    e = np.arange(T0, T1, dt)
    if e[-1] < T1:
        e = np.append(e, T1)
    return e


STARTS = [(0.50, 0.020, 1.00, 0.010, 1.10),
          (0.20, 0.005, 1.50, 0.050, 1.20),
          (0.80, 0.050, 0.80, 0.005, 1.05)]


def fit_region(args):
    name, box, trig_rec, offsets, use_epoch_mu = args
    t_run = time.time()
    trig = pd.DataFrame(trig_rec)
    df = load_region(name)
    t_all_full = df["t"].to_numpy()
    T_first, T_last = float(t_all_full[0]), float(t_all_full[-1])
    span = T_last - T_first
    T_expl = T_first + EXPLORE_FRAC * span
    T_train = T_first + EXPLORE_FRAC * TRAIN_FRAC_OF_EXPLORE * span
    # HOLDOUT BLINDNESS: nothing beyond T_expl is ever touched again.
    df = df[df["t"] < T_expl].reset_index(drop=True)
    t = df["t"].to_numpy()
    m = df["mag"].to_numpy() - M0
    lat = df["latitude"].to_numpy()
    lon = df["longitude"].to_numpy()
    T_fit0 = T_first + 365.0
    if not (T_fit0 < T_train < T_expl):
        raise RuntimeError("region %s: degenerate split" % name)

    # ---- Mc per 5-yr epoch, exploration only ----
    ep_all = epoch_edges_days(T_fit0, T_expl)
    mc_epoch = {}
    for i in range(len(ep_all) - 1):
        sel = (t >= ep_all[i]) & (t < ep_all[i + 1])
        if sel.sum() >= 200:
            mc_epoch["%d" % i] = dict(n=int(sel.sum()),
                                      mc=max_curvature_mc(df["mag"].to_numpy()[sel]))
    mcs = [v["mc"] for v in mc_epoch.values() if v["mc"] is not None]
    mc_range = (max(mcs) - min(mcs)) if len(mcs) >= 2 else 0.0

    # ---- epochs used by the model ----
    if use_epoch_mu:
        ep_edges = epoch_edges_days(T_fit0, T_expl)
    else:
        ep_edges = np.array([T_fit0 - 1.0, T_expl + 1.0])
    E_all = len(ep_edges) - 1
    ep_len_fit = epoch_lengths(ep_edges, T_fit0, T_train)
    live = ep_len_fit > 0
    n_live = int(live.sum())
    ep_idx_map = np.cumsum(live) - 1          # global epoch -> fitted-parameter index
    ep_idx_map = np.where(live, ep_idx_map, n_live - 1)   # later epochs inherit the last fitted mu

    fit_lo = int(np.searchsorted(t, T_fit0, "left"))
    fit_hi = int(np.searchsorted(t, T_train, "left"))
    n_fit = fit_hi - fit_lo
    if n_fit < 100:
        raise RuntimeError("region %s: only %d training events" % (name, n_fit))
    ep_tgt_fit = ep_idx_map[np.clip(np.searchsorted(ep_edges, t[fit_lo:fit_hi], "right") - 1,
                                    0, E_all - 1)]
    ep_len_fit_live = ep_len_fit[live]
    rate = n_fit / (T_train - T_fit0)

    args_tr = (t, m, fit_lo, fit_hi, T_fit0, T_train, TRUNC_W_FIT_D,
               ep_tgt_fit, ep_len_fit_live, True, 512)
    lb = np.log(np.concatenate([np.full(n_live, 1e-6), [1e-6, 0.5, 1e-6, 0.8]]))
    ub = np.log(np.concatenate([np.full(n_live, 1e3), [1e2, 2.5, 1e1, 2.0]]))
    bounds = list(zip(lb, ub))
    best = None
    for (mf, K, al, c, p) in STARTS:
        x0 = np.log(np.clip(np.concatenate([np.full(n_live, mf * rate), [K, al, c, p]]),
                            np.exp(lb) * 1.001, np.exp(ub) * 0.999))
        r = minimize(_obj, x0, args=args_tr, jac=True, method="L-BFGS-B", bounds=bounds,
                     options={"maxiter": 200, "maxfun": 400, "ftol": 1e-12, "gtol": 1e-8})
        if best is None or -r.fun > best[0]:
            best = (float(-r.fun), np.exp(r.x))
    LL_A_train_trunc, theta = best
    mu_live = theta[:n_live]
    K, alpha, c, p = (float(v) for v in theta[n_live:])
    mu_vec = mu_live[ep_idx_map]              # one mu per GLOBAL epoch (carry-forward)

    # ---- quadrature cells + exact lambda_A integrals ----
    w_all = K * np.exp(alpha * LN10 * m)
    e_tr = make_edges(T_fit0, T_train)
    e_sc = make_edges(T_train, T_expl)
    src_tr = t < T_train
    A_tr = cell_integrals(e_tr, t[src_tr], w_all[src_tr], c, p, mu_vec, ep_edges)
    A_sc = cell_integrals(e_sc, t, w_all, c, p, mu_vec, ep_edges)
    mid_tr = 0.5 * (e_tr[:-1] + e_tr[1:])
    mid_sc = 0.5 * (e_sc[:-1] + e_sc[1:])

    sc_lo = int(np.searchsorted(t, T_train, "left"))
    t_ev_tr = t[fit_lo:fit_hi]
    t_ev_sc = t[sc_lo:]
    n_sc = t_ev_sc.size
    if n_sc < 50:
        raise RuntimeError("region %s: only %d scoring events" % (name, n_sc))

    # ---- model-A skill vs Poisson on the scoring window (sanity, not the statistic) ----
    S_sc, _, _, _ = _pair_sums(t, np.exp(alpha * LN10 * m), m, c, p, sc_lo, len(t),
                               np.inf, False, 256)
    ep_sc = np.clip(np.searchsorted(ep_edges, t_ev_sc, "right") - 1, 0, E_all - 1)
    lam_sc = mu_vec[ep_sc] + K * S_sc
    LL_A_sc = float(np.log(lam_sc).sum() - A_sc.sum())
    lam0 = n_fit / (T_train - T_fit0)
    LL_pois = n_sc * np.log(lam0) - lam0 * (T_expl - T_train)
    bits_A_vs_pois = (LL_A_sc - LL_pois) / n_sc / LN2

    # ---- triggers for this region ----
    clat, clon = (box[0] + box[1]) / 2.0, (box[2] + box[3]) / 2.0
    t_arr, amp, r_km, tmag, n_keep, n_tot = build_trigger_set(trig, box, clat, clon)
    if n_keep < 1000:
        raise RuntimeError("region %s: only %d triggers survive exclusion (expected >1000)"
                           % (name, n_keep))

    Tg0, Tg1 = float(trig["t"].min()), float(trig["t"].max())
    Gspan = Tg1 - Tg0

    out = dict(region=name, cls=CLASS[name], n_events_explore=int(len(t)),
               n_train=int(n_fit), n_score=int(n_sc),
               t_first=str(pd.Timestamp(T_first * DAY, unit="s", tz="UTC")),
               t_explore_end=str(pd.Timestamp(T_expl * DAY, unit="s", tz="UTC")),
               t_train_end=str(pd.Timestamp(T_train * DAY, unit="s", tz="UTC")),
               span_days=float(span), mc_by_epoch=mc_epoch, mc_range=float(mc_range),
               n_epochs_fitted=int(n_live),
               params=dict(mu=[float(v) for v in mu_live], K=K, alpha=alpha, c=c, p=p),
               LL_A_train_truncW=float(LL_A_train_trunc),
               bits_A_vs_poisson_score=float(bits_A_vs_pois),
               n_triggers_kept=int(n_keep), n_triggers_total=int(n_tot),
               trigger_max_amp_kPa=float(np.max(amp) / 1e3),
               dist_km_min=float(np.min(r_km)), taus={})

    per_tau = {}
    for tau in TAUS:
        g_ev_tr = glink(decay_series(t_ev_tr, t_arr, amp, tau))
        g_ev_sc = glink(decay_series(t_ev_sc, t_arr, amp, tau))
        g_c_tr = glink(decay_series(mid_tr, t_arr, amp, tau))
        g_c_sc = glink(decay_series(mid_sc, t_arr, amp, tau))
        beta, dLL_tr = fit_beta(g_ev_tr, A_tr, g_c_tr)
        dLL_sc = beta_profile(g_ev_sc, A_sc, g_c_sc, beta)
        dInt = float(np.sum(A_sc * (np.exp(beta * g_c_sc) - 1.0)))
        contrib = (beta * g_ev_sc - dInt / n_sc) / LN2       # per-event bits, for the bootstrap
        # ---- circular-shift null (beta REFIT on training each time) ----
        nullbits = np.empty(len(offsets))
        for i, off in enumerate(offsets):
            ta = Tg0 + np.mod(t_arr - Tg0 + off, Gspan)
            o = np.argsort(ta)
            ta, am = ta[o], amp[o]
            gj = glink(decay_series(t_ev_tr, ta, am, tau))
            gc = glink(decay_series(mid_tr, ta, am, tau))
            b2, _ = fit_beta(gj, A_tr, gc)
            gjs = glink(decay_series(t_ev_sc, ta, am, tau))
            gcs = glink(decay_series(mid_sc, ta, am, tau))
            nullbits[i] = beta_profile(gjs, A_sc, gcs, b2) / n_sc / LN2
        per_tau[str(tau)] = dict(
            beta=float(beta), dLL_train=float(dLL_tr), dLL_score=float(dLL_sc),
            bits_per_event=float(dLL_sc / n_sc / LN2),
            bits_train=float(dLL_tr / n_fit / LN2),
            mean_g_score=float(g_ev_sc.mean()), max_g_score=float(g_ev_sc.max()),
            null_bits=nullbits.tolist(), contrib=contrib.tolist(),
            dLL_score_total=float(dLL_sc))
    out["taus"] = per_tau
    # ---- payload for the sensitivity / type-I-error arms ----
    # Coarse 1-day quadrature cells (the covariate is smooth at tau >= 1 d) plus, for each
    # tau, the covariate field for the OBSERVED trigger catalogue (row 0) and for the first
    # N_FPR_SHIFT circular shifts (rows 1..S). Cells are binned on g so the beta fits inside
    # the Monte-Carlo loops are O(FPR_BINS) instead of O(cells).
    def _binmap(A, g, nb=FPR_BINS):
        gmax = float(g.max())
        if gmax <= 0:
            idx = np.zeros(g.size, np.int16)
            return idx, np.array([A.sum()]), np.array([0.0])
        edges = np.linspace(0.0, gmax * (1 + 1e-9), nb + 1)
        idx = np.clip(np.searchsorted(edges, g, "right") - 1, 0, nb - 1).astype(np.int16)
        Ab = np.bincount(idx, weights=A, minlength=nb)
        Gb = np.bincount(idx, weights=A * g, minlength=nb)
        Gb = np.where(Ab > 0, Gb / np.maximum(Ab, 1e-300), 0.0)
        return idx, Ab, Gb

    ce_tr = make_edges(T_fit0, T_train, FPR_CELL_DAYS)
    ce_sc = make_edges(T_train, T_expl, FPR_CELL_DAYS)
    cA_tr = cell_integrals(ce_tr, t[src_tr], w_all[src_tr], c, p, mu_vec, ep_edges)
    cA_sc = cell_integrals(ce_sc, t, w_all, c, p, mu_vec, ep_edges)
    cm_tr = 0.5 * (ce_tr[:-1] + ce_tr[1:])
    cm_sc = 0.5 * (ce_sc[:-1] + ce_sc[1:])
    n_tr_real = np.histogram(t_ev_tr, bins=ce_tr)[0].astype(float)
    n_sc_real = np.histogram(t_ev_sc, bins=ce_sc)[0].astype(float)
    if abs(n_tr_real.sum() - n_fit) > 0.5 or abs(n_sc_real.sum() - n_sc) > 0.5:
        raise RuntimeError("region %s: coarse-cell binning lost events (%d/%d, %d/%d)"
                           % (name, n_tr_real.sum(), n_fit, n_sc_real.sum(), n_sc))
    S_fpr = min(N_FPR_SHIFT, len(offsets))
    sens_pl = {}
    for tau in TAUS:
        rows = []
        for si in range(S_fpr + 1):
            if si == 0:
                ta, am = t_arr, amp
            else:
                tashift = Tg0 + np.mod(t_arr - Tg0 + offsets[si - 1], Gspan)
                o = np.argsort(tashift)
                ta, am = tashift[o], amp[o]
            rows.append((glink(decay_series(cm_tr, ta, am, tau)),
                         glink(decay_series(cm_sc, ta, am, tau))))
        it, At, Gt = zip(*[_binmap(cA_tr, r[0]) for r in rows])
        isc, Asc, Gsc = zip(*[_binmap(cA_sc, r[1]) for r in rows])
        sens_pl[str(tau)] = dict(
            A_tr_cells=cA_tr, A_sc_cells=cA_sc,
            n_tr_real=n_tr_real, n_sc_real=n_sc_real,
            g_tr_cells_obs=rows[0][0].astype(np.float32),
            g_sc_cells_obs=rows[0][1].astype(np.float32),
            idx_tr=[x for x in it], idx_sc=[x for x in isc],
            A_tr_bins=[x for x in At], g_tr_bins=[x for x in Gt],
            A_sc_bins=[x for x in Asc], g_sc_bins=[x for x in Gsc],
            n_shift=S_fpr)
    out["_sens"] = sens_pl

    # ---- sequences for the block bootstrap (scoring window only) ----
    seq = sequence_ids(t_ev_sc, lat[sc_lo:], lon[sc_lo:], df["mag"].to_numpy()[sc_lo:])
    out["seq"] = seq.tolist()
    out["n_sequences_score"] = int(len(np.unique(seq)))

    # ---- A1 coda-blindness reading ----
    win = 0.25
    sel_expl = (t_arr >= T_fit0) & (t_arr < T_expl - 1.0)
    ta_e, am_e = t_arr[sel_expl], amp[sel_expl]
    top = np.argsort(am_e)[::-1][:20]
    mags = df["mag"].to_numpy()
    small_a = small_b = big_a = big_b = 0
    for k in top:
        t0 = ta_e[k]
        for (lo_, hi_, is_after) in ((t0, t0 + win, True), (t0 - 1.0, t0 - 1.0 + win, False)):
            s = (t >= lo_) & (t < hi_)
            ns = int((mags[s] < 2.0).sum())
            nb = int((mags[s] >= 2.0).sum())
            if is_after:
                small_a += ns
                big_a += nb
            else:
                small_b += ns
                big_b += nb
    out["coda_blindness"] = dict(
        window_hours=6.0, n_top_amplitudes=int(len(top)),
        top_amp_kPa=[float(x / 1e3) for x in am_e[top]],
        n_M_lt_2_after=small_a, n_M_ge_2_after=big_a,
        n_M_lt_2_matched_1d_before=small_b, n_M_ge_2_matched_1d_before=big_b)

    # ---- positive control: Landers ----
    if name in POSCTRL_CELLS:
        tl = (pd.Timestamp(LANDERS["t0"]) - pd.Timestamp(0, tz="UTC")).total_seconds() / DAY
        rL = float(gcdist_km(LANDERS["lat"], LANDERS["lon"], clat, clon))
        tarrL = tl + (rL / 3.5) / DAY
        pc = dict(dist_km=rL, sigma_kPa=float(sigma_primary_pa(LANDERS["M"], rL) / 1e3),
                  in_training_window=bool(T_fit0 <= tarrL < T_train))
        nbg = int(((t >= tarrL - 90.0) & (t < tarrL)).sum())
        for wname, wd in (("0-1d", 1.0), ("0-5d", 5.0)):
            npost = int(((t >= tarrL) & (t < tarrL + wd)).sum())
            expct = nbg / 90.0 * wd
            pc[wname] = dict(n_post=npost, n_bg_90d=nbg, expected=float(expct),
                             rate_ratio=(float(npost / expct) if expct > 0 else None))
        out["positive_control_landers"] = pc

    out["runtime_s"] = round(time.time() - t_run, 1)
    print("[region] %-20s n_sc=%6d  bits(t1)=%+.4f  bits(t5)=%+.4f  beta=%.3f/%.3f  %.0fs"
          % (name, n_sc, per_tau["1.0"]["bits_per_event"], per_tau["5.0"]["bits_per_event"],
             per_tau["1.0"]["beta"], per_tau["5.0"]["beta"], out["runtime_s"]), flush=True)
    return out


REGION_BOX = {}


def main():
    t_start = time.time()
    rng = np.random.default_rng(SEED)
    man = json.load(open(os.path.join(DK034, "manifest.json")))
    global REGION_BOX
    REGION_BOX = {k: v["box"] for k, v in man["cells"].items()}
    assert set(REGION_BOX) == set(CLASS), "manifest cells != declared 14 regions"

    csv_p = download_triggers()
    trig = pd.read_csv(csv_p, usecols=["time", "latitude", "longitude", "mag", "id"])
    tt = pd.to_datetime(trig["time"], utc=True, format="mixed")
    trig = trig.assign(t=(tt - pd.Timestamp(0, tz="UTC")).dt.total_seconds().to_numpy() / DAY)
    trig = trig.dropna(subset=["t", "latitude", "longitude", "mag"]).sort_values("t")
    trig = trig[trig["mag"] >= TRIG_MINMAG - 1e-9].reset_index(drop=True)
    if len(trig) < 5000:
        raise RuntimeError("trigger catalogue: expected >5000 rows, got %d" % len(trig))
    inbox = np.zeros(len(trig), bool)
    for b in REGION_BOX.values():
        inbox |= in_box(trig["latitude"].to_numpy(), trig["longitude"].to_numpy(), b)
    trig["in_any_k034_box"] = inbox
    print("[trig] %d M>=%.1f events, %d inside a K-034 box (dropped globally)"
          % (len(trig), TRIG_MINMAG, int(inbox.sum())), flush=True)

    # ---- FROZEN completeness rule, decided before any skill number exists ----
    mc_probe = {}
    for name in CLASS:
        df = load_region(name)
        ta = df["t"].to_numpy()
        T0f, T1f = float(ta[0]), float(ta[-1])
        Te = T0f + EXPLORE_FRAC * (T1f - T0f)
        d = df[df["t"] < Te]
        ee = epoch_edges_days(T0f, Te)
        vals = []
        for i in range(len(ee) - 1):
            s = (d["t"].to_numpy() >= ee[i]) & (d["t"].to_numpy() < ee[i + 1])
            if s.sum() >= 200:
                v = max_curvature_mc(d["mag"].to_numpy()[s])
                if v is not None:
                    vals.append(v)
        rates = []
        for i in range(len(ee) - 1):
            lo_ = max(ee[i], T0f); hi_ = min(ee[i + 1], Te)
            if hi_ - lo_ > 180:
                s_ = (d["t"].to_numpy() >= lo_) & (d["t"].to_numpy() < hi_)
                rates.append(float(int(s_.sum()) / (hi_ - lo_)))
        mc_probe[name] = dict(mc_by_epoch=vals,
                              mc_range=float(max(vals) - min(vals)) if len(vals) >= 2 else 0.0,
                              rate_per_day_by_epoch=rates,
                              rate_ratio_max_min=(float(max(rates) / min(rates))
                                                  if rates and min(rates) > 0 else None))
    max_range = max(v["mc_range"] for v in mc_probe.values())
    use_epoch_mu = bool(max_range >= MC_DRIFT_TRIGGER)
    print("[mc] max across-epoch Mc range = %.2f -> epoch mu %s"
          % (max_range, "ON" if use_epoch_mu else "OFF"), flush=True)

    Tg0, Tg1 = float(trig["t"].min()), float(trig["t"].max())
    Gspan = Tg1 - Tg0
    offsets = rng.uniform(SHIFT_MIN_YR * 365.25, Gspan - SHIFT_MIN_YR * 365.25, N_SHIFT)

    trec = {k: trig[k].to_numpy() for k in
            ["t", "latitude", "longitude", "mag", "in_any_k034_box"]}
    jobs = [(n, REGION_BOX[n], trec, offsets, use_epoch_mu) for n in sorted(CLASS)]
    results = []
    with ProcessPoolExecutor(max_workers=N_PROC) as ex:
        for r in ex.map(fit_region, jobs):
            results.append(r)
    if len(results) != 14:
        raise RuntimeError("expected 14 region results, got %d" % len(results))

    # ---------------- pooled statistic ----------------
    N_tot = sum(r["n_score"] for r in results)
    pooled = {}
    for tau in TAUS:
        k = str(tau)
        dll = sum(r["taus"][k]["dLL_score_total"] for r in results)
        pooled[k] = dict(bits_per_event=float(dll / N_tot / LN2), dLL=float(dll))
    tau_best = max(pooled, key=lambda k: pooled[k]["bits_per_event"])
    primary = pooled[tau_best]["bits_per_event"]

    # ---------------- circular-shift null, max-statistic over the tau family ----------------
    null_pool = {}
    for tau in TAUS:
        k = str(tau)
        nb = np.zeros(N_SHIFT)
        for r in results:
            nb += np.asarray(r["taus"][k]["null_bits"]) * r["n_score"]
        null_pool[k] = nb / N_tot
    null_max = np.maximum.reduce([null_pool[str(t_)] for t_ in TAUS])
    obs_max = max(pooled[str(t_)]["bits_per_event"] for t_ in TAUS)
    p_null = float((1 + int((null_max >= obs_max).sum())) / (1 + N_SHIFT))
    p_per_tau = {str(t_): float((1 + int((null_pool[str(t_)] >= pooled[str(t_)]["bits_per_event"]).sum()))
                                / (1 + N_SHIFT)) for t_ in TAUS}

    # ---------------- block bootstrap over sequences ----------------
    contrib, blocks, wts = [], [], []
    off = 0
    for r in results:
        c_ = np.asarray(r["taus"][tau_best]["contrib"])
        s_ = np.asarray(r["seq"]) + off
        off = int(s_.max()) + 1 if s_.size else off
        contrib.append(c_)
        blocks.append(s_)
    contrib = np.concatenate(contrib)
    blocks = np.concatenate(blocks)
    uniq, inv = np.unique(blocks, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    inv_s = inv[order]
    c_s = contrib[order]
    starts = np.searchsorted(inv_s, np.arange(len(uniq)), "left")
    ends = np.searchsorted(inv_s, np.arange(len(uniq)), "right")
    csum = np.concatenate([[0.0], np.cumsum(c_s)])
    blk_sum = csum[ends] - csum[starts]
    blk_n = (ends - starts).astype(float)
    boot = np.empty(N_BOOT)
    for b in range(N_BOOT):
        idx = rng.integers(0, len(uniq), len(uniq))
        boot[b] = blk_sum[idx].sum() / blk_n[idx].sum()
    ci = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]

    # ---------------- by class ----------------
    by_class = {}
    for cl in ("AB", "C"):
        sel = [r for r in results if (r["cls"] in ("A", "B")) == (cl == "AB")]
        n = sum(r["n_score"] for r in sel)
        d = sum(r["taus"][tau_best]["dLL_score_total"] for r in sel)
        by_class[cl] = dict(regions=[r["region"] for r in sel], n_score=int(n),
                            bits_per_event=float(d / n / LN2))

    # ---------------- sensitivity (PLAYBOOK rule 5) ----------------
    sens = sensitivity(results, null_max, rng)
    fpr = false_positive_rate(results, rng)
    pooled_beta = pooled_beta_arm(results)
    for r in results:
        r.pop("_sens", None)

    # ---------------- assemble ----------------
    coda = dict(n_M_lt_2_after=sum(r["coda_blindness"]["n_M_lt_2_after"] for r in results),
                n_M_ge_2_after=sum(r["coda_blindness"]["n_M_ge_2_after"] for r in results),
                n_M_lt_2_matched_1d_before=sum(
                    r["coda_blindness"]["n_M_lt_2_matched_1d_before"] for r in results),
                n_M_ge_2_matched_1d_before=sum(
                    r["coda_blindness"]["n_M_ge_2_matched_1d_before"] for r in results))
    passed = bool(primary >= 0.01 and p_null < 0.05)

    res = dict(
        experiment="K-038 ARM A -- dynamic-stress weather as a forecast covariate",
        state_class="first-run, exploration split",
        run_utc=pd.Timestamp.utcnow().isoformat(),
        frozen_constants=dict(
            amplitude_axis="van der Elst & Brodsky (2010) eq.(6), verbatim from "
                           "K034_SEALED_LITERATURE.md; G=30 GPa, c=3.5 km/s, T=20 s",
            G_shear_Pa=G_SHEAR, c_phase_m_s=C_PHASE, T_surface_wave_s=T_SW,
            link="g(D) = log1p(D/D0)", D0_Pa=D0_PA, tau_family_days=list(TAUS),
            trigger_minmag=TRIG_MINMAG, trigger_span=[TRIG_T0, TRIG_T1],
            exclusion_box_km=EXCL_BOX_KM, exclusion_rupture_lengths=2.0,
            rupture_length_relation="log10 L[km] = %.2f + %.2f M (Wells & Coppersmith 1994 "
                                    "RLD, all slip types) -- UNVERIFIED in session" % (RUPT_A, RUPT_B),
            k034_box_exclusion=True, M0=M0,
            explore_frac=EXPLORE_FRAC, train_frac_of_explore=TRAIN_FRAC_OF_EXPLORE,
            epoch_years=EPOCH_YEARS, epoch_origin=EPOCH_ORIGIN,
            mc_drift_trigger=MC_DRIFT_TRIGGER, grid_dt_days=GRID_DT_D,
            n_shifts=N_SHIFT, shift_min_years=SHIFT_MIN_YR, n_bootstrap=N_BOOT,
            decluster_km=DECLUSTER_KM_LOCAL, decluster_days=DECLUSTER_DAYS_LOCAL,
            seed=SEED,
            success_rule=">= 0.01 bits/event pooled OOS AND p < 0.05 under the "
                         "circular-shift null, max-statistic over the tau family"),
        completeness=dict(max_across_epoch_mc_range=float(max_range),
                          epoch_mu_used=use_epoch_mu, per_region=mc_probe),
        triggers=dict(n_rows=int(len(trig)), n_inside_k034_boxes=int(inbox.sum()),
                      span=[str(pd.Timestamp(Tg0 * DAY, unit="s", tz="UTC")),
                            str(pd.Timestamp(Tg1 * DAY, unit="s", tz="UTC"))]),
        headline=dict(
            primary_bits_per_event=float(primary), tau_selected_days=float(tau_best),
            pooled_by_tau={k: v["bits_per_event"] for k, v in pooled.items()},
            n_scored_events=int(N_tot),
            n_sequences=int(len(uniq)),
            ci95_block_bootstrap_over_sequences=ci,
            p_circular_shift_max_over_tau=p_null,
            p_circular_shift_per_tau=p_per_tau,
            null_bits_mean=float(null_max.mean()),
            null_bits_p95=float(np.percentile(null_max, 95)),
            PASS=passed),
        by_k034_class=by_class,
        per_region=[{k: v for k, v in r.items()
                     if k not in ("seq",)} for r in results],
        coda_blindness_pooled=coda,
        positive_control=[r["positive_control_landers"] | {"region": r["region"]}
                          for r in results if "positive_control_landers" in r],
        sensitivity=sens,
        false_positive_rate_beta0=fpr,
        pooled_beta_exploratory=pooled_beta,
        prior_art=dict(
            note=("NOT a new phenomenon. Continuous-field framings of dynamic triggering are "
                  "already published; the only thing this arm contributes is the "
                  "incremental-bits-per-event forecast metric on an out-of-sample split."),
            continuous_field_framings=["DeSalvio & Fan (QTM)", "Miyazawa", "Brodsky", "Guo"],
            type_I_error_reference="Hardebeck, DeSalvio, Fan & Barbour 2025, JGR, "
                                   "doi:10.1029/2025JB031566"),
        artifacts=dict(
            A1_coda_blindness="reported in coda_blindness_pooled; biases AGAINST a positive beta",
            A2_in_box_aftershock_contamination="handled by the 300 km box gate, the "
                                               "2-rupture-length gate and the K-034-box exclusion",
            A3_completeness_drift="handled by 5-yr epoch mu in BOTH models (see completeness)",
            A4_day_night="irrelevant at tau = 1-5 d; trigger origin times are uniform over the "
                         "solar day, so a diurnal detection modulation is orthogonal to D_r(t)",
            A5_etas_absorption="ETAS explains the aftershocks of remotely triggered events, "
                               "which biases AGAINST model B"),
        runtime_minutes=round((time.time() - t_start) / 60.0, 2))

    # strip the bulky null/contrib vectors from the per-region dump
    for r in res["per_region"]:
        for k in r["taus"]:
            r["taus"][k] = {kk: vv for kk, vv in r["taus"][k].items()
                            if kk not in ("null_bits", "contrib")}
        r.pop("coda_blindness", None)

    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=2)

    print("=" * 74)
    print("K-038 ARM A  |  pooled OOS incremental bits/event = %+.5f (tau=%s d)"
          % (primary, tau_best))
    print("  95%% block-bootstrap CI over %d sequences: [%+.5f, %+.5f]" % (len(uniq), ci[0], ci[1]))
    print("  circular-shift null p (max over tau family) = %.4f" % p_null)
    print("  by class: A+B %+.5f   C %+.5f"
          % (by_class["AB"]["bits_per_event"], by_class["C"]["bits_per_event"]))
    print("  min detectable beta at 80%% power = %s (least favourable tau %s d)"
          % (sens.get("min_beta_80pct_power"), sens.get("least_favourable_tau")))
    print("  false-positive rate at beta=0 over %d synthetic catalogues = %.3f (nominal %.3f)"
          % (fpr["n_sim"], fpr["false_positive_rate"], fpr["nominal_alpha"]))
    print("  EXPLORATORY single pooled beta: %+.5f bits/event, p = %.3f"
          % (pooled_beta["bits_max_over_tau"],
             pooled_beta["p_circular_shift_max_over_tau"]))
    print("  SUCCESS RULE (>=0.01 bits AND p<0.05): %s" % ("PASS" if passed else "FAIL"))
    print("  runtime %.1f min -> %s" % (res["runtime_minutes"], os.path.basename(OUT)))
    print("=" * 74)
    return res


def _fit_beta_bins(sum_g, Ab, gb):
    r = minimize_scalar(lambda b: -(b * sum_g - float(np.sum(Ab * (np.exp(b * gb) - 1.0)))),
                        bounds=BETA_BOUNDS, method="bounded", options={"xatol": 1e-6})
    return float(r.x)


def _score_sim(pl, tn, si, n_tr, n_sc):
    """Identical pipeline on a simulated catalogue: fit beta on the training cells under
    covariate field `si` (0 = observed triggers, si>0 = circular shift si), then score on
    the scoring cells. Returns (delta log-likelihood, n_scored, beta_hat)."""
    q = pl[tn]
    gb_t = q["g_tr_bins"][si]
    gb_s = q["g_sc_bins"][si]
    cnt_t = np.bincount(q["idx_tr"][si], weights=n_tr, minlength=gb_t.size)
    cnt_s = np.bincount(q["idx_sc"][si], weights=n_sc, minlength=gb_s.size)
    b = _fit_beta_bins(float(cnt_t @ gb_t), q["A_tr_bins"][si], gb_t)
    dll = (b * float(cnt_s @ gb_s)
           - float(np.sum(q["A_sc_bins"][si] * (np.exp(b * gb_s) - 1.0))))
    return dll, float(n_sc.sum()), b


def sensitivity(results, null_max, rng, n_sim=200,
                betas=(0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0)):
    """PLAYBOOK rule 5 -- injection recovery, measured at the LEAST FAVOURABLE tau.

    Generator: per region the FITTED ETAS intensity path supplies the rate on 1-day cells,
    and events are drawn as an inhomogeneous Poisson process with intensity
    lambda_A(t) * exp(beta_true * g(D_r(t))). The downstream pipeline is identical (beta fit
    on the training cells, frozen, scored on the scoring cells, pooled over regions).
    Threshold = the 95th percentile of the REAL circular-shift null of the pooled statistic.
    DECLARED LIMITATION: the generator has no ETAS feedback -- injected events do not spawn
    their own aftershocks, part of which model A would absorb -- so this power estimate is
    mildly OPTIMISTIC.
    """
    payload = [r["_sens"] for r in results]
    crit = float(np.percentile(null_max, 95))
    out = dict(threshold_bits=crit, n_sim=n_sim, betas=list(betas),
               generator="inhomogeneous Poisson on the fitted ETAS intensity path x exp(beta g)",
               by_tau={})
    for tn in [str(t_) for t_ in TAUS]:
        pw, meanbits = [], []
        for bt in betas:
            hits, acc = 0, []
            for _ in range(n_sim):
                dll_tot, n_tot = 0.0, 0.0
                for pl in payload:
                    q = pl[tn]
                    n_tr = rng.poisson(q["A_tr_cells"] * np.exp(bt * q["g_tr_cells_obs"]))
                    n_sc = rng.poisson(q["A_sc_cells"] * np.exp(bt * q["g_sc_cells_obs"]))
                    if n_tr.sum() < 10 or n_sc.sum() < 10:
                        continue
                    d, n, _ = _score_sim(pl, tn, 0, n_tr, n_sc)
                    dll_tot += d
                    n_tot += n
                bits = dll_tot / n_tot / LN2 if n_tot else 0.0
                acc.append(bits)
                if bits > crit:
                    hits += 1
            pw.append(hits / n_sim)
            meanbits.append(float(np.mean(acc)))
        idx = next((i for i, v in enumerate(pw) if v >= 0.80), None)
        out["by_tau"][tn] = dict(power=pw, mean_pooled_bits=meanbits,
                                 min_beta_80pct=(betas[idx] if idx is not None else None))
        print("[sens] tau=" + tn + " power=" + repr(pw), flush=True)
    vals = [v["min_beta_80pct"] for v in out["by_tau"].values()]
    out["least_favourable_tau"] = max(
        out["by_tau"], key=lambda k: (out["by_tau"][k]["min_beta_80pct"] is None,
                                      out["by_tau"][k]["min_beta_80pct"] or 0))
    out["min_beta_80pct_power"] = (None if any(v is None for v in vals) else max(vals))
    return out


def pooled_beta_arm(results):
    """EXPLORATORY, added AFTER the primary run and flagged as such (the K-034 D3 precedent).

    The pre-registered model fits ONE beta PER REGION. With 14 regions that is 14 parameters
    estimated on training windows that carry very little information about beta, and the
    out-of-sample penalty for a noise-fitted parameter is real: the circular-shift null of the
    primary statistic is itself centred well above zero. This arm asks the same question with a
    SINGLE beta shared by all 14 regions. It runs on the 1-day Monte-Carlo cells and the 40
    Monte-Carlo shifts (not the 0.25-day cells and 200 shifts of the primary), so it is coarser
    as well as exploratory. It sets no flag.
    """
    payload = [r["_sens"] for r in results]
    tns = [str(t_) for t_ in TAUS]
    S = int(payload[0][tns[0]]["n_shift"])

    def score(tn, si, use_real=True):
        Sg_tr = 0.0
        parts_tr, parts_sc, Sg_sc, ntot = [], [], 0.0, 0.0
        for pl in payload:
            q = pl[tn]
            gt, gs = q["g_tr_bins"][si], q["g_sc_bins"][si]
            ct = np.bincount(q["idx_tr"][si], weights=q["n_tr_real"], minlength=gt.size)
            cs = np.bincount(q["idx_sc"][si], weights=q["n_sc_real"], minlength=gs.size)
            Sg_tr += float(ct @ gt)
            Sg_sc += float(cs @ gs)
            ntot += float(q["n_sc_real"].sum())
            parts_tr.append((q["A_tr_bins"][si], gt))
            parts_sc.append((q["A_sc_bins"][si], gs))

        def negll(b):
            return -(b * Sg_tr - sum(float(np.sum(A * (np.exp(b * g) - 1.0)))
                                     for A, g in parts_tr))
        r = minimize_scalar(negll, bounds=BETA_BOUNDS, method="bounded",
                            options={"xatol": 1e-7})
        b = float(r.x)
        dll = b * Sg_sc - sum(float(np.sum(A * (np.exp(b * g) - 1.0))) for A, g in parts_sc)
        return b, dll / ntot / LN2

    obs, nul, betas = {}, {}, {}
    for tn in tns:
        b0, v0 = score(tn, 0)
        obs[tn], betas[tn] = v0, b0
        nul[tn] = np.array([score(tn, si)[1] for si in range(1, S + 1)])
    o = max(obs.values())
    nm = np.maximum.reduce([nul[tn] for tn in tns])
    return dict(exploratory=True, n_shifts=S, cell_days=FPR_CELL_DAYS,
                beta_by_tau=betas, bits_by_tau=obs,
                bits_max_over_tau=float(o),
                p_circular_shift_max_over_tau=float((1 + int((nm >= o).sum())) / (1 + S)),
                null_mean=float(nm.mean()), null_p95=float(np.percentile(nm, 95)))


def false_positive_rate(results, rng, n_sim=N_FPR_SIM):
    """Type-I error of the WHOLE pipeline on no-triggering synthetic catalogues.

    Prior art (via Merton): Hardebeck, DeSalvio, Fan & Barbour (2025), JGR Solid Earth,
    doi:10.1029/2025JB031566 -- standard remote-triggering tests return 3.5-8.5% false
    positives on synthetics with NO triggering, because the null does not preserve
    clustering. This arm measures the equivalent number for THIS pipeline.

    Per simulated catalogue: draw events per region as Poisson on the fitted ETAS intensity
    path with beta_true = 0 (so the clustering of the real fitted intensity IS preserved),
    then run the identical procedure -- fit beta on training, score out of sample, compare
    the pooled statistic to the SAME circular-shift null (trigger catalogue shifted, the
    simulated catalogue held fixed), max-statistic over the tau family.
    SCALED TO BUDGET: N_FPR_SIM catalogues x N_FPR_SHIFT shifts, so the achievable nominal
    alpha is 2/(1+N_FPR_SHIFT), not exactly 0.05.
    """
    payload = [r["_sens"] for r in results]
    tn0 = str(TAUS[0])
    S = int(payload[0][tn0]["n_shift"])
    alpha_nom = 2.0 / (1.0 + S)
    ps = []
    for k in range(n_sim):
        draws = [(rng.poisson(pl[tn0]["A_tr_cells"]), rng.poisson(pl[tn0]["A_sc_cells"]))
                 for pl in payload]
        obs, nul = {}, {}
        for tn in [str(t_) for t_ in TAUS]:
            arr = np.zeros(S + 1)
            ntot = 0.0
            for pl, (n_tr, n_sc) in zip(payload, draws):
                if n_tr.sum() < 10 or n_sc.sum() < 10:
                    continue
                for si in range(S + 1):
                    d, n, _ = _score_sim(pl, tn, si, n_tr, n_sc)
                    arr[si] += d
                ntot += float(n_sc.sum())
            arr = arr / ntot / LN2
            obs[tn], nul[tn] = arr[0], arr[1:]
        o = max(obs.values())
        nm = np.maximum.reduce([nul[tn] for tn in [str(t_) for t_ in TAUS]])
        ps.append(float((1 + int((nm >= o).sum())) / (1 + S)))
        if (k + 1) % 10 == 0:
            print("[fpr] %d/%d sims, running FPR = %.3f"
                  % (k + 1, n_sim, float((np.array(ps) < 0.05).mean())), flush=True)
    ps = np.array(ps)
    return dict(n_sim=int(n_sim), n_shifts=int(S),
                nominal_alpha=float(alpha_nom),
                false_positive_rate=float((ps < 0.05).mean()),
                fpr_at_nominal_alpha=float((ps <= alpha_nom + 1e-12).mean()),
                p_values=ps.tolist(),
                literature=("Hardebeck, DeSalvio, Fan & Barbour 2025, JGR, "
                            "doi:10.1029/2025JB031566 report 3.5-8.5 percent false positives "
                            "for standard remote-triggering tests on no-triggering synthetics"),
                note="scaled to the wall-clock budget")


if __name__ == "__main__":
    main()

