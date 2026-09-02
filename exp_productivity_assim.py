"""K-436 (first run, EXPLORATION split): sequence-adaptive productivity as incremental
bits/event over frozen temporal ETAS.

THE QUESTION.  Frozen ETAS gives every mainshock of magnitude M the same productivity
K*10^(alpha*(M-M0)).  Kepler's read-only sniff (KEPLER_SEEDS_20260902.md SNIFF-2/SNIFF-3)
measured a ~1.9 bits/aftershock oracle in per-sequence productivity and a partial Spearman
rho = 0.41 between the hour-1 aftershock count and day 1-30 productivity at fixed magnitude.
This module asks the proper question instead: does updating each ongoing sequence's
productivity from its OWN first hour (or first three hours) add bits/event to the ETAS
forecast of the following days, walk-forward?  Per S-3 the gate is incremental bits over
frozen ETAS, not AUC and not a magnitude regression.

DESIGN (everything below is frozen before the run).

  Catalogue: data/xue_lu_zenodo/SCSN_original_catalog.txt (raw column order, auto-detected
  exactly as exp_h_etas.load_catalog), SoCal box lat [31.5, 38.0], lon [-122.0, -113.5],
  magnitude floor M0 = 2.5 (exp_h's frozen floor).

  Split (PLAYBOOK rule 6, engine/splits.py): EXPLORATION = first 70% of the catalogue time
  span.  The last 30% is HOLDOUT and is never read: the catalogue is truncated to the
  exploration span immediately after loading, before any statistic is computed.  Inside
  exploration: TRAINING = first 60% of the exploration span (ETAS MLE + productivity prior),
  SCORING = the remaining 40%, walk-forward.

  Model A (baseline): frozen temporal ETAS,
      lambda(t) = mu + sum_{t_i<t} K * 10^(alpha*(M_i-M0)) * (t-t_i+c)^(-p),
  parameters by MLE on the training window using the exp_h_etas recipe (log-parameters,
  L-BFGS-B, four starts, 1000 d truncated fit then untruncated polish), scored walk-forward
  over the scoring window with history = all prior events, untruncated kernel.  The
  intensity form, the likelihood and the bits/event scoring are exp_h's code, imported.

  Model B (adaptive): identical, except that each parent event with M >= M_PARENT (4.0) has
  its productivity multiplied, from t_i + T_obs onward only, by exp(posterior mean of
  log-productivity deviation).  Before t_i + T_obs the two models are numerically identical,
  so every incremental bit accrues after the observation window.  Family of 2, declared:
  T_obs in {1 h, 3 h}; the success rule is evaluated at the LEAST FAVOURABLE member (this is
  the max-statistic payment for a family of two).

  The Bayesian update.  Prior on delta_i = log-productivity deviation: a Gaussian fitted to
  the training-window per-sequence residuals log((N_obs+0.5)/(N_pred+0.5)) over the 30-day
  window (0.5 is a frozen continuity constant; N_obs = 0 occurs).  Likelihood: the count n_i
  of the sequence's own attributed aftershocks in (t_i, t_i+T_obs] is Poisson with mean
  exp(delta_i) * Lambda0_i, Lambda0_i = K*10^(alpha*(M_i-M0)) * int_0^{T_obs} (s+c)^-p ds.
  Posterior mean of delta_i by 1-D quadrature.  Attribution (frozen, crude, stated): an event
  is attributed to the nearest-in-time PRIOR parent that is larger in magnitude, within 30 d
  and within max(3L, 10 km), L = 10^(0.5M-1.8) km.

  Primary statistic (frozen): incremental bits/event of B over A, summed over ALL scoring
  events and divided by ALL scoring events (the S-3 form); and the same total divided by the
  aftershock-attributed scoring events only (the K-436 form).  Per-event score
  s_j = [log lambda_B - log lambda_A](t_j) - (Lambda_B - Lambda_A)/n, so mean(s)/ln2 is
  exactly the bits/event; CIs are block-bootstrapped over SEQUENCES (PLAYBOOK rule 1), the
  sequence layer being the program's standing 30 d / 150 km nearest-larger rule.

  SUCCESS RULE (frozen): >= 0.02 bits/event overall with the sequence-block-bootstrap 95% CI
  excluding 0, at the least favourable T_obs.

  Also reported: parents restricted to M >= 5; the three largest scoring-window sequences
  excluded; and the alpha-refit harsher baseline (S-5) in which a magnitude-only productivity
  correction exp(a + b*(M-4)), fitted to the same training residuals, is given to model A from
  t_i onward -- i.e. could a better-fitted productivity/magnitude law alone absorb the gain?

ARTIFACTS THAT COULD FAKE THIS RESULT (PLAYBOOK rule 7, named before the run).

  (a) SHORT-TERM AFTERSHOCK INCOMPLETENESS (STAI).  The first hour after a large mainshock is
      incomplete at M2.5.  The hour-1 count is therefore biased DOWN, and biased down MORE for
      larger mainshocks, so the update would be wrong in a magnitude-dependent way and could
      manufacture (or destroy) apparent skill correlated with magnitude.  Control: count
      inside T_obs only above the frozen magnitude-dependent completeness floor
      Mc(t) = M_main - 4.5 - 0.75*log10(t in days) (Helmstetter, Kagan & Jackson 2006), and
      scale the predicted mean by the Gutenberg-Richter fraction 10^(-b*(Mc(t)-M0)) inside the
      Omori integral.  Reported WITH and WITHOUT the correction, and the STAI-CORRECTED
      READING IS THE HEADLINE (Hainzl 2016 shows STAI can manufacture apparent structure in
      exactly these first-hour counts).  A THIRD, incompleteness-robust count is also run:
      a van der Elst (2021) b-positive-style rule that counts only those events inside T_obs
      that are LARGER than their immediate predecessor in the same sequence, whose expected
      number under GR is exactly half the Omori expectation and which is insensitive to a
      time-varying detection floor.
  (b) ATTRIBUTION ERROR.  The nearest-in-time-larger-within-max(3L,10km) rule is crude: it
      mixes in independent background events and misses off-plane aftershocks.  Both push the
      observed count away from the parent's true direct-offspring count and therefore blur the
      update; the direction is not signed a priori.  It is frozen, identical in A-side prior
      fitting and B-side updating, and identical in the simulation arm.
  (c) SECONDARY TRIGGERING absorbing the gain.  Real ETAS already recovers part of a
      sequence's excess productivity through its own aftershocks-of-aftershocks.  That is
      exactly why the comparison here is against FULL walk-forward ETAS rather than against a
      magnitude regression (which is what SNIFF-2's 1.9 bits was measured against): the
      baseline is allowed to keep everything it can already do.
  (d) THE DAY/NIGHT (diurnal completeness) ARTIFACT IS IRRELEVANT HERE.  It operates at the
      network's detection threshold, around M1-1.5 in SoCal; at M >= 2.5 the catalogue is
      complete around the clock, and in any case both models see the identical event stream,
      so a diurnal detection modulation cancels in the A-vs-B difference.

SENSITIVITY (PLAYBOOK rule 5, measured before any null is quoted).  ETAS catalogues are
simulated from the training fit with per-parent log10-productivity scatter of known sd in
{0, 0.25, 0.5, 0.75, 1.0}; the whole pipeline (attribution, prior fit, update, scoring,
sequence bootstrap) is re-run on each.  It must return ~0 bits at sd = 0 (it must not
manufacture gain) and positive bits at sd = 1.0; the minimum sd detectable at 80% power is
reported.  The injected multiplier is MEAN-PRESERVING (10^(sd*Z)/E[10^(sd*Z)], Z truncated
at +/-2): a raw log-normal multiplier raises the mean productivity and, at a fitted branching
ratio near 1, drives the simulated catalogue supercritical, which is a property of the
generator and not of the estimator being calibrated.

RUN NOTE.  The training MLE is cached in productivity_assim_fit_cache.json, keyed on the
split, the data window and the whole fit recipe; it is reused only on an exact key match.

PRIOR ART (Merton, folded in before the run; this arm is a REDISCOVERY, measured in the
program's own S-3 bits metric on its own catalogue, not a new phenomenon).  USGS Operational
Aftershock Forecasting already fits sequence-specific productivity from the aftershocks
observed so far, with Bayesian updating against a generic prior.  Omi, Ogata, Hirata & Aihara
(2013) did the first-hours version and explicitly solved short-term incompleteness.  Woessner
et al. (2011) showed retrospectively that parameter updating improves predictive power.
Hainzl (2016) is the named hazard for the estimator used here.  What is new here is only the
accounting: incremental bits/event over a full frozen walk-forward ETAS on this catalogue and
this split, with a sequence-block CI and a zero-scatter calibration.

Outputs: results_productivity_assim.json.  Run: python -u exp_productivity_assim.py
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from exp_h_etas import (
    LN2,
    LN10,
    _G_terms,
    _obj,
    _pair_sums,
    aki_b,
    etas_ll,
    load_catalog,
    max_curvature_mc,
)

HERE = Path(__file__).resolve().parent
OUT = HERE / "results_productivity_assim.json"

# ----------------------------------------------------------------- frozen constants
CATALOG = "SCSN_original_catalog.txt"
LAT_MIN, LAT_MAX = 31.5, 38.0
LON_MIN, LON_MAX = -122.0, -113.5
M0 = 2.5                       # exp_h's frozen floor
EXPLORE_FRAC = 0.70            # engine/splits.py
TRAIN_FRAC_OF_EXPLORE = 0.60
BURN_IN_DAYS = 365.0           # exp_h
TRUNC_W_DAYS = 1000.0          # exp_h
MAXITER = 150                  # exp_h
STARTS = [                     # exp_h
    (0.50, 0.020, 1.00, 0.010, 1.10),
    (0.20, 0.005, 1.50, 0.050, 1.20),
    (0.80, 0.050, 0.80, 0.005, 1.05),
    (0.35, 0.010, 1.20, 0.020, 1.15),
]
ALPHA_BOUNDS = (0.5, 2.5)
P_BOUNDS = (0.8, 2.0)
FIT_BUDGET_MIN = 45.0

M_PARENT = 4.0                 # declared parent floor for the update
M_PARENT_HIGH = 5.0            # secondary reading
T_OBS_HOURS = (1.0, 3.0)       # THE DECLARED FAMILY OF 2
ATTR_DAYS = 30.0               # attribution / prior window
ATTR_MIN_KM = 10.0
ATTR_L_MULT = 3.0
CONTINUITY = 0.5               # log((N_obs+0.5)/(N_pred+0.5))
STAI_A, STAI_B = 4.5, 0.75     # Mc(t) = M_main - 4.5 - 0.75 log10(t days)
COUNT_MODES = ("plain", "stai", "bpos")   # bpos = van der Elst-style larger-than-predecessor
BPOS_FRACTION = 0.5            # P(m_k > m_{k-1}) = 1/2 exactly for exponential (GR) marks
SEQ_DAYS, SEQ_KM = 30.0, 150.0  # program's standing sequence rule (exp_world_faultrelative)
N_BOOT = 2000
BOOT_SEED = 20260902
SUCCESS_BITS = 0.02

# sensitivity arm
SIM_SDS = (0.0, 0.25, 0.5, 0.75, 1.0)
SIM_REPS = 12
SIM_YEARS = 36.0              # chosen so the simulated SCORING window holds ~9,000
                               # events, matching the real scoring window (rule 5: the
                               # instrument must be calibrated at the size it is used at)
SIM_SEED = 4360902
SIM_MAX_EVENTS = 80_000        # a rep whose catalogue hits this cap is DISCARDED and counted
SIM_MAX_GEN = 60
SIM_Z_TRUNC = 2.0              # per-parent log-productivity draw truncated at +/- 2 sd
SIM_LOC_SD_FRAC = 0.30         # offspring location sd = 0.30 * max(3L, 10 km)
SIM_TARGET_BRANCHING = 0.80    # K is rescaled if the fitted branching ratio is >= this


# ----------------------------------------------------------------- small helpers
def great_circle_km(lat1, lon1, lat2, lon2):
    r = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp = p2 - p1
    dl = np.radians(np.asarray(lon2) - np.asarray(lon1))
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * r * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def omori_int(a, b, c, p):
    """int_a^b (s + c)^-p ds, s measured from the parent origin time (days)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    b = np.maximum(b, a)
    if abs(1.0 - p) < 1e-9:
        return np.log(b + c) - np.log(a + c)
    s = 1.0 - p
    return ((b + c) ** s - (a + c) ** s) / s


def rupture_len_km(mag):
    return 10.0 ** (0.5 * np.asarray(mag, dtype=float) - 1.8)


def attr_radius_km(mag):
    return np.maximum(ATTR_L_MULT * rupture_len_km(mag), ATTR_MIN_KM)


def sequence_ids(t, lat, lon, mag, days=SEQ_DAYS, km=SEQ_KM):
    """Program's standing nearest-larger-within-(days, km) sequence layer
    (same algorithm as exp_fluid_driven._sequence_ids, whose W constants are 30 d / 150 km)."""
    n = t.size
    seq = np.arange(n)
    recent = []
    for ii in np.argsort(t):
        while recent and t[ii] - t[recent[0]] > days:
            recent.pop(0)
        for jj in reversed(recent):
            if mag[jj] >= mag[ii] and great_circle_km(lat[jj], lon[jj], lat[ii], lon[ii]) <= km:
                seq[ii] = seq[jj]
                break
        recent.append(ii)
    return seq


def attribute(t, lat, lon, mag, parent_idx, window_days=ATTR_DAYS):
    """FROZEN crude attribution: each event -> the nearest-in-time prior parent that is
    larger in magnitude, within window_days and within max(3L, 10 km) of that parent.
    Returns an int array of parent positions (index into parent_idx) or -1."""
    n = t.size
    out = np.full(n, -1, dtype=int)
    if parent_idx.size == 0:
        return out
    tp = t[parent_idx]
    mp = mag[parent_idx]
    latp, lonp = lat[parent_idx], lon[parent_idx]
    radp = attr_radius_km(mp)
    order = np.argsort(t)
    for ii in order:
        hi = int(np.searchsorted(tp, t[ii], side="left"))
        lo = int(np.searchsorted(tp, t[ii] - window_days, side="left"))
        for k in range(hi - 1, lo - 1, -1):
            if mp[k] <= mag[ii]:
                continue
            if great_circle_km(latp[k], lonp[k], lat[ii], lon[ii]) <= radp[k]:
                out[ii] = k
                break
    return out


def posterior_mean_log(n_obs, lam0, prior_mu, prior_sd, ngrid=3001):
    """Posterior mean of delta with delta ~ N(prior_mu, prior_sd^2) and
    n_obs ~ Poisson(exp(delta) * lam0).  Quadrature on a fixed grid spanning both the
    prior and every parent's Poisson MLE."""
    n_obs = np.asarray(n_obs, dtype=float)
    lam0 = np.maximum(np.asarray(lam0, dtype=float), 1e-12)
    mle = np.log((n_obs + CONTINUITY) / lam0)
    glo = min(prior_mu - 8.0 * prior_sd, float(mle.min()) - 5.0)
    ghi = max(prior_mu + 8.0 * prior_sd, float(mle.max()) + 5.0)
    g = np.linspace(glo, ghi, ngrid)
    lp = (-0.5 * ((g[None, :] - prior_mu) / prior_sd) ** 2
          + n_obs[:, None] * g[None, :]
          - lam0[:, None] * np.exp(g[None, :]))
    lp -= lp.max(axis=1, keepdims=True)
    w = np.exp(lp)
    return (w * g[None, :]).sum(axis=1) / w.sum(axis=1)


def stai_mc(m_main, s_days):
    return m_main - STAI_A - STAI_B * np.log10(np.maximum(s_days, 1e-9))


def stai_expected_frac(m_main, T_obs, c, p, b, ngrid=4000):
    """int_0^T (s+c)^-p * min(1, 10^{-b (Mc(s)-M0)}) ds  /  1 , per parent (vector over parents)."""
    s = np.geomspace(1e-7, T_obs, ngrid)
    ker = (s + c) ** (-p)
    mc = stai_mc(np.asarray(m_main, dtype=float)[:, None], s[None, :])
    frac = np.minimum(1.0, 10.0 ** (-b * (mc - M0)))
    return np.trapezoid(ker[None, :] * frac, s, axis=1)


def stai_count_mask(dt_days, mag_ev, m_main):
    """Keep an attributed event only if its magnitude is at or above the STAI floor."""
    return mag_ev >= np.maximum(M0, stai_mc(m_main, dt_days)) - 1e-9


def block_bootstrap_ci(scores, blocks, n_boot=N_BOOT, seed=BOOT_SEED):
    """Block bootstrap over SEQUENCES of the per-event score (already in nats/event)."""
    rng = np.random.default_rng(seed)
    uniq, inv = np.unique(blocks, return_inverse=True)
    nb = uniq.size
    sums = np.bincount(inv, weights=scores, minlength=nb)
    cnts = np.bincount(inv, minlength=nb).astype(float)
    draws = rng.integers(0, nb, size=(n_boot, nb))
    tot = sums[draws].sum(axis=1)
    num = cnts[draws].sum(axis=1)
    means = tot / np.maximum(num, 1.0)
    lo, hi = np.percentile(means, [2.5, 97.5])
    return float(lo / LN2), float(hi / LN2), int(nb), float(np.std(means) / LN2)


# ----------------------------------------------------------------- core scoring engine
class Scorer:
    """Holds the frozen ETAS pieces for one catalogue + one scoring window, and can add
    per-parent productivity multipliers cheaply on top of the baseline."""

    def __init__(self, t, mag, theta, T0, T1, sc_lo, sc_hi):
        self.t = t
        self.m = mag - M0
        self.mu, self.K, self.alpha, self.c, self.p = theta
        self.T0, self.T1 = T0, T1
        self.lo, self.hi = sc_lo, sc_hi
        self.w = np.exp(self.alpha * LN10 * self.m)
        S, _, _, _ = _pair_sums(t, self.w, self.m, self.c, self.p, sc_lo, sc_hi,
                                np.inf, False, 1024)
        self.lam_A = self.mu + self.K * S
        G, _, _ = _G_terms(t, T0, T1, self.c, self.p, np.inf, False)
        self.Lam_A = self.mu * (T1 - T0) + self.K * float((self.w * G).sum())
        self.tgt = t[sc_lo:sc_hi]
        self.n = sc_hi - sc_lo

    def add_multipliers(self, parent_pos, mult, start_t):
        """lambda_B = lambda_A + K * sum (mult-1) * w_i * (t-t_i+c)^-p for t > start_t_i.
        Returns (lam_B, Lam_B)."""
        add = np.zeros(self.n)
        int_add = 0.0
        for k, i in enumerate(parent_pos):
            r = mult[k] - 1.0
            if r == 0.0:
                continue
            st = start_t[k]
            j0 = int(np.searchsorted(self.tgt, st, side="right"))
            if j0 < self.n:
                dt = self.tgt[j0:] - self.t[i]
                add[j0:] += r * self.w[i] * (dt + self.c) ** (-self.p)
            a = max(self.T0, st) - self.t[i]
            bnd = self.T1 - self.t[i]
            if bnd > a:
                int_add += r * self.w[i] * float(omori_int(max(a, 0.0), bnd, self.c, self.p))
        return self.lam_A + self.K * add, self.Lam_A + self.K * int_add


def per_event_scores(lam_A, lam_B, Lam_A, Lam_B):
    n = lam_A.size
    return (np.log(lam_B) - np.log(lam_A)) - (Lam_B - Lam_A) / n


# ----------------------------------------------------------------- the update pipeline
def bpos_counts(owner, dt, mag, inwin, n_parents):
    """van der Elst-style incompleteness-robust count: within each parent's window, in time
    order, count the events strictly larger than their immediate predecessor.  Under GR marks
    the expectation is exactly half the Omori expectation, independent of a moving Mc(t)."""
    out = np.zeros(n_parents)
    idx = np.where(inwin)[0]
    if idx.size == 0:
        return out
    idx = idx[np.lexsort((dt[idx], owner[idx]))]
    prev_owner, prev_mag = -1, None
    for j in idx:
        o = owner[j]
        if o != prev_owner:
            prev_owner, prev_mag = o, mag[j]
            continue
        if mag[j] > prev_mag:
            out[o] += 1.0
        prev_mag = mag[j]
    return out


def build_update(t, lat, lon, mag, theta, parent_mask, train_lo_t, train_hi_t,
                 T_obs_days, count_mode, b_value, prior=None, owner=None):
    """Returns dict with parent positions, multipliers, update start times, prior, counts."""
    mu, K, alpha, c, p = theta
    pidx = np.where(parent_mask)[0]
    if owner is None:
        owner = attribute(t, lat, lon, mag, pidx)
    w = 10.0 ** (alpha * (mag - M0))

    # ---- prior from TRAINING-window parents, 30 d window, no STAI correction (frozen)
    if prior is None:
        in_tr = (t[pidx] >= train_lo_t) & (t[pidx] + ATTR_DAYS <= train_hi_t)
        tr_k = np.where(in_tr)[0]
        dt_all = np.where(owner >= 0, t - t[pidx][np.maximum(owner, 0)], -1.0)
        in_win = (owner >= 0) & (dt_all > 0) & (dt_all <= ATTR_DAYS)
        n30 = np.bincount(owner[in_win], minlength=pidx.size).astype(float)
        pred30 = K * w[pidx] * omori_int(0.0, ATTR_DAYS, c, p)
        delta_tr = np.log((n30[tr_k] + CONTINUITY) / (pred30[tr_k] + CONTINUITY))
        if delta_tr.size < 10:
            prior = {"mu": 0.0, "sd": 1.0, "n_train_parents": int(delta_tr.size),
                     "note": "too few training parents; fell back to N(0,1)"}
            slope = 0.0
            intercept = 0.0
        else:
            prior = {"mu": float(delta_tr.mean()), "sd": float(max(delta_tr.std(ddof=1), 1e-3)),
                     "n_train_parents": int(delta_tr.size)}
            X = np.column_stack([np.ones(delta_tr.size), mag[pidx][tr_k] - M_PARENT])
            coef, *_ = np.linalg.lstsq(X, delta_tr, rcond=None)
            intercept, slope = float(coef[0]), float(coef[1])
        prior["alpha_refit_intercept"] = intercept
        prior["alpha_refit_slope_per_mag_natlog"] = slope
        prior["alpha_refit_implied_alpha"] = float(alpha + slope / LN10)
        prior["train_residual_summary"] = {
            "n": int(delta_tr.size),
            "mean_natlog": float(delta_tr.mean()) if delta_tr.size else None,
            "sd_natlog": float(delta_tr.std(ddof=1)) if delta_tr.size > 1 else None,
            "sd_log10": float(delta_tr.std(ddof=1) / LN10) if delta_tr.size > 1 else None,
        }

    # ---- counts inside T_obs (with or without STAI floor)
    dt = np.where(owner >= 0, t - t[pidx][np.maximum(owner, 0)], -1.0)
    inwin = (owner >= 0) & (dt > 0) & (dt <= T_obs_days)
    if count_mode == "stai":
        mm = mag[pidx][np.maximum(owner, 0)]
        inwin = inwin & stai_count_mask(np.maximum(dt, 1e-9), mag, mm)
        lam0 = K * w[pidx] * stai_expected_frac(mag[pidx], T_obs_days, c, p, b_value)
        n_obs = np.bincount(owner[inwin], minlength=pidx.size).astype(float)
    elif count_mode == "bpos":
        lam0 = BPOS_FRACTION * K * w[pidx] * omori_int(0.0, T_obs_days, c, p)
        n_obs = bpos_counts(owner, dt, mag, inwin, pidx.size)
    elif count_mode == "plain":
        lam0 = K * w[pidx] * omori_int(0.0, T_obs_days, c, p)
        n_obs = np.bincount(owner[inwin], minlength=pidx.size).astype(float)
    else:
        raise ValueError(f"unknown count_mode {count_mode!r}")

    delta_post = posterior_mean_log(n_obs, lam0, prior["mu"], prior["sd"])
    return {
        "pidx": pidx,
        "mult": np.exp(delta_post),
        "start_t": t[pidx] + T_obs_days,
        "prior": prior,
        "n_obs": n_obs,
        "lam0": lam0,
        "delta_post": delta_post,
        "owner": owner,
    }


def score_variant(sc, upd, keep_mask=None):
    lam_B, Lam_B = sc.add_multipliers(upd["pidx"], upd["mult"], upd["start_t"])
    s = per_event_scores(sc.lam_A, lam_B, sc.Lam_A, Lam_B)
    if keep_mask is not None:
        s = s[keep_mask]
    return s, lam_B, Lam_B


# ----------------------------------------------------------------- ETAS simulation
def _mean_preserving_norm(sd_log10, ztrunc=None):
    """E[10^(sd*Z)] for Z ~ N(0,1) truncated at +/- ztrunc.  Dividing the injected multiplier
    by this keeps the MEAN productivity (and hence the branching ratio) unchanged, so that
    injecting scatter does not silently make the simulated catalogue supercritical."""
    if sd_log10 == 0:
        return 1.0
    zt = SIM_Z_TRUNC if ztrunc is None else ztrunc
    z = np.linspace(-zt, zt, 20001)
    wz = np.exp(-0.5 * z * z)
    return float(np.trapezoid(wz * 10.0 ** (sd_log10 * z), z) / np.trapezoid(wz, z))


def simulate_etas(theta, b, T, box_km, rng, sd_log10, m_parent=M_PARENT):
    mu, K, alpha, c, p = theta
    lat0, lat1, lon0, lon1 = box_km
    n_bg = rng.poisson(mu * T)
    t = list(rng.uniform(0, T, n_bg))
    m = list(M0 - np.log10(rng.uniform(size=n_bg)) / b)
    la = list(rng.uniform(lat0, lat1, n_bg))
    lo = list(rng.uniform(lon0, lon1, n_bg))
    gen_start = 0
    n_gen = 0
    while gen_start < len(t) and n_gen < SIM_MAX_GEN and len(t) < SIM_MAX_EVENTS:
        n_gen += 1
        idx = np.arange(gen_start, len(t))
        gen_start = len(t)
        tt = np.array([t[i] for i in idx])
        mm = np.array([m[i] for i in idx])
        lla = np.array([la[i] for i in idx])
        llo = np.array([lo[i] for i in idx])
        kap = K * 10.0 ** (alpha * (mm - M0))
        if sd_log10 > 0:
            big = mm >= m_parent
            z = np.clip(rng.normal(size=mm.size), -SIM_Z_TRUNC, SIM_Z_TRUNC)
            norm = _mean_preserving_norm(sd_log10)
            kap = kap * np.where(big, 10.0 ** (sd_log10 * z) / norm, 1.0)
        expct = kap * omori_int(0.0, np.maximum(T - tt, 0.0), c, p)
        nof = rng.poisson(np.clip(expct, 0, 1e6))
        for k in np.where(nof > 0)[0]:
            nk = int(nof[k])
            Imax = omori_int(0.0, T - tt[k], c, p)
            u = rng.uniform(0, 1, nk) * Imax
            s = 1.0 - p
            if abs(s) < 1e-9:
                ds = np.exp(u + np.log(c)) - c
            else:
                ds = (u * s + c ** s) ** (1.0 / s) - c
            tc = tt[k] + ds
            rad = max(ATTR_L_MULT * 10.0 ** (0.5 * mm[k] - 1.8), ATTR_MIN_KM)
            sdkm = SIM_LOC_SD_FRAC * rad
            dlat = rng.normal(0, sdkm / 111.0, nk)
            dlon = rng.normal(0, sdkm / (111.0 * math.cos(math.radians(lla[k]))), nk)
            mc = M0 - np.log10(rng.uniform(size=nk)) / b
            t.extend(tc.tolist())
            m.extend(mc.tolist())
            la.extend((lla[k] + dlat).tolist())
            lo.extend((llo[k] + dlon).tolist())
            if len(t) >= SIM_MAX_EVENTS:
                break
    t = np.array(t)
    o = np.argsort(t)
    return (t[o], np.array(m)[o], np.array(la)[o], np.array(lo)[o],
            {"n_generations": n_gen, "capped": bool(len(t) >= SIM_MAX_EVENTS)})


# ----------------------------------------------------------------- main
def main():
    t_start = time.time()
    res = {
        "experiment": "K-436 sequence-adaptive productivity (incremental bits over frozen ETAS)",
        "state_class": "first-run, exploration split",
        "run_utc": pd.Timestamp.utcnow().isoformat(),
        "frozen_constants": {
            "catalog": CATALOG, "box": {"lat": [LAT_MIN, LAT_MAX], "lon": [LON_MIN, LON_MAX]},
            "M0": M0, "explore_frac": EXPLORE_FRAC,
            "train_frac_of_explore": TRAIN_FRAC_OF_EXPLORE,
            "burn_in_days": BURN_IN_DAYS, "trunc_W_days_for_fit": TRUNC_W_DAYS,
            "maxiter": MAXITER, "alpha_bounds": list(ALPHA_BOUNDS), "p_bounds": list(P_BOUNDS),
            "M_parent": M_PARENT, "M_parent_high": M_PARENT_HIGH,
            "T_obs_hours_family": list(T_OBS_HOURS),
            "count_modes": list(COUNT_MODES),
            "bpos_fraction": BPOS_FRACTION,
            "headline_count_mode": "stai",
            "attribution": {"window_days": ATTR_DAYS, "radius": "max(3L, 10 km), L=10^(0.5M-1.8) km",
                            "rule": "nearest-in-time prior LARGER parent within window and radius"},
            "continuity_constant": CONTINUITY,
            "stai_Mc": "M_main - 4.5 - 0.75*log10(t days)  [Helmstetter et al. 2006]",
            "sequence_layer": {"days": SEQ_DAYS, "km": SEQ_KM,
                               "rule": "nearest-larger-within (program standing rule)"},
            "n_boot": N_BOOT, "boot_seed": BOOT_SEED,
            "success_rule": f">= {SUCCESS_BITS} bits/event overall AND 95% sequence-block CI "
                            f"excludes 0, at the LEAST FAVOURABLE T_obs",
            "sim": {"sds_log10": list(SIM_SDS), "reps": SIM_REPS, "years": SIM_YEARS,
                    "seed": SIM_SEED, "loc_sd_frac_of_radius": SIM_LOC_SD_FRAC,
                    "target_branching_if_supercritical": SIM_TARGET_BRANCHING},
        },
        "artifacts_named_before_run": {
            "a_STAI": "first hour incomplete at M2.5, magnitude-dependent; controlled by Mc(t) "
                      "floor on counting AND a GR fraction inside the Omori integral; both "
                      "readings reported",
            "b_attribution": "crude nearest-larger-in-time within max(3L,10km); frozen and "
                             "identical in prior fit, update and simulation",
            "c_secondary_triggering": "baseline is FULL walk-forward ETAS, so secondary "
                                      "triggering is already credited to A",
            "d_diurnal": "irrelevant at M>=2.5 (detection artifact lives near M1-1.5) and it "
                         "cancels in the A-vs-B difference since both see the same events",
        },
    }

    # ---------------- data + splits ----------------
    df, order = load_catalog(CATALOG)
    box = df.lat.between(LAT_MIN, LAT_MAX) & df.lon.between(LON_MIN, LON_MAX)
    d_all = df[box].sort_values("t").reset_index(drop=True)
    t_first, t_last = float(d_all.t.iloc[0]), float(d_all.t.iloc[-1])
    span = t_last - t_first
    t_explore_end = t_first + EXPLORE_FRAC * span
    t_train_end = t_first + EXPLORE_FRAC * TRAIN_FRAC_OF_EXPLORE * span

    # HOLDOUT GUARD: truncate to exploration immediately; nothing below ever sees the last 30%.
    d = d_all[d_all.t < t_explore_end].reset_index(drop=True)
    del d_all, df

    def as_ts(x):
        return str(pd.Timestamp(x * 86400.0, unit="s", tz="UTC"))

    cat = d[d.mag >= M0 - 1e-9].reset_index(drop=True)
    t = cat.t.to_numpy(float)
    mag = cat.mag.to_numpy(float)
    lat = cat.lat.to_numpy(float)
    lon = cat.lon.to_numpy(float)
    tr_hi = int(np.searchsorted(t, t_train_end, side="left"))
    sc_lo, sc_hi = tr_hi, len(t)
    n_train, n_score = tr_hi, sc_hi - sc_lo
    assert n_train > 0 and n_score > 0

    mc_expl = max_curvature_mc(mag)
    res["catalog"] = {
        "detected_column_order": order, "n_in_box_full_span": None,
        "span_start": as_ts(t_first), "span_end_FULL_CATALOG": as_ts(t_last),
        "span_days": span,
        "exploration_end": as_ts(t_explore_end),
        "train_end": as_ts(t_train_end),
        "n_events_exploration_all_mags": int(len(d)),
        "n_events_exploration_M>=2.5": int(len(cat)),
        "n_train": int(n_train), "n_score": int(n_score),
        "mc_maxcurv_exploration_M>=2.5_subset": mc_expl,
        "holdout_guard": "catalogue truncated to t < exploration_end before any statistic",
    }
    print(f"[data] {order} order; exploration {as_ts(t_first)} .. {as_ts(t_explore_end)}")
    print(f"[split] train < {as_ts(t_train_end)}: n={n_train}; scoring: n={n_score}")

    # ---------------- ETAS MLE on the training window (exp_h recipe) ----------------
    m_rel = mag - M0
    T0_hist = float(t[0])
    T0_fit = T0_hist + BURN_IN_DAYS
    T1_fit = t_train_end
    fit_lo = int(np.searchsorted(t, T0_fit, side="left"))
    fit_hi = tr_hi
    n_fit = fit_hi - fit_lo
    D_fit = T1_fit - T0_fit
    rate_fit = n_fit / D_fit
    t_fit, m_fit = t, m_rel
    subsampled = False

    args_tr = (t_fit, m_fit, fit_lo, fit_hi, T0_fit, T1_fit, TRUNC_W_DAYS, True, 512)
    t0 = time.time()
    etas_ll(np.array([rate_fit * 0.5, 0.02, 1.0, 0.01, 1.1]), *args_tr)
    ev_s = time.time() - t0
    proj = ev_s * MAXITER * len(STARTS) / 60.0
    print(f"[fit] n_fit={n_fit} rate={rate_fit:.3f}/d  1 eval {ev_s:.2f}s  proj {proj:.1f} min")
    if proj > FIT_BUDGET_MIN:
        subsampled = True
        T0_hist = T1_fit - 10 * 365.25
        T0_fit = T0_hist + BURN_IN_DAYS
        hlo = int(np.searchsorted(t, T0_hist, side="left"))
        t_fit, m_fit = t[hlo:fit_hi], m_rel[hlo:fit_hi]
        fit_lo = int(np.searchsorted(t_fit, T0_fit, side="left"))
        fit_hi = len(t_fit)
        n_fit = fit_hi - fit_lo
        D_fit = T1_fit - T0_fit
        rate_fit = n_fit / D_fit
        args_tr = (t_fit, m_fit, fit_lo, fit_hi, T0_fit, T1_fit, TRUNC_W_DAYS, True, 512)
        print(f"[fit] RUNTIME GUARD: train history restricted to most recent 10 yr (n_fit={n_fit})")

    CACHE = HERE / "productivity_assim_fit_cache.json"
    cache_key = {"n_fit": int(n_fit), "T0_fit": float(T0_fit), "T1_fit": float(T1_fit),
                 "subsampled": bool(subsampled), "M0": M0, "starts": STARTS,
                 "trunc": TRUNC_W_DAYS, "maxiter": MAXITER}
    cached = None
    if CACHE.exists():
        try:
            cc = json.loads(CACHE.read_text(encoding="utf-8"))
            if cc.get("key") == json.loads(json.dumps(cache_key)):
                cached = cc
                print("[fit] REUSING cached training MLE (identical split, data and recipe)")
        except Exception as e:
            print(f"[fit] cache unreadable ({e}); refitting")

    lb = np.log([1e-6, 1e-6, ALPHA_BOUNDS[0], 1e-6, P_BOUNDS[0]])
    ub = np.log([1e3, 1e2, ALPHA_BOUNDS[1], 1e1, P_BOUNDS[1]])
    bounds = list(zip(lb, ub))
    starts_out, best = [], None
    if cached is not None:
        starts_out = cached["starts"]
        theta = np.array(cached["theta"], dtype=float)
    else:
     for si, (mu_f, K_, al, c_, p_) in enumerate(STARTS):
        x0 = np.log(np.clip([mu_f * rate_fit, K_, al, c_, p_], np.exp(lb) * 1.001, np.exp(ub) * 0.999))
        r = minimize(_obj, x0, args=args_tr, jac=True, method="L-BFGS-B", bounds=bounds,
                     options={"maxiter": MAXITER, "maxfun": MAXITER * 2, "ftol": 1e-12, "gtol": 1e-8})
        th = np.exp(r.x)
        rec = {"LL": float(-r.fun),
               "params": dict(zip(["mu", "K", "alpha", "c", "p"], map(float, th)))}
        starts_out.append(rec)
        print(f"[fit] start {si+1}: LL={rec['LL']:.2f} " +
              " ".join(f"{k}={v:.5f}" for k, v in rec["params"].items()))
        if best is None or rec["LL"] > best["LL"]:
            best = rec
     else:
      theta_tr = np.array([best["params"][k] for k in ["mu", "K", "alpha", "c", "p"]])
      args_full = (args_tr[0], args_tr[1], args_tr[2], args_tr[3], args_tr[4], args_tr[5],
                   np.inf, True, 1024)
      LL_full0, _ = etas_ll(theta_tr, *args_full[:6], np.inf, False, 1024)
      rf = minimize(_obj, np.log(theta_tr), args=args_full, jac=True, method="L-BFGS-B",
                    bounds=bounds, options={"maxiter": 20, "maxfun": 40, "ftol": 1e-12})
      theta = np.exp(rf.x)
      if -rf.fun < LL_full0:
          theta = theta_tr
      CACHE.write_text(json.dumps({"key": cache_key, "theta": [float(x) for x in theta],
                                   "starts": starts_out}, indent=1), encoding="utf-8")
    mu, K, alpha, c, p = map(float, theta)
    b_train = aki_b(mag[:tr_hi], M0)
    br = None
    if p > 1.0 and b_train > alpha:
        br = float(K * c ** (1 - p) / (p - 1) * (b_train / (b_train - alpha)))
    print(f"[params FROZEN] mu={mu:.4f} K={K:.5f} alpha={alpha:.3f} c={c:.5f} p={p:.4f} "
          f"b={b_train:.3f} n_branch={br}")
    res["train_fit"] = {
        "fit_window": [as_ts(T0_fit), as_ts(T1_fit)], "n_target_events": int(n_fit),
        "train_subsampled_to_recent_10yr": subsampled,
        "starts": starts_out,
        "frozen_params": {"mu": mu, "K": K, "alpha": alpha, "c": c, "p": p, "M0": M0},
        "b_value_train_aki": b_train, "branching_ratio_n": br,
    }

    # ---------------- baseline scoring ----------------
    T0_sc, T1_sc = t_train_end, float(t[-1])
    sc = Scorer(t, mag, (mu, K, alpha, c, p), T0_sc, T1_sc, sc_lo, sc_hi)
    LL_A = float(np.log(sc.lam_A).sum() - sc.Lam_A)
    assert np.isfinite(LL_A) and np.all(np.isfinite(sc.lam_A)) and sc.lam_A.min() > 0
    print(f"[A] LL_A={LL_A:.2f}  lam range [{sc.lam_A.min():.4f}, {sc.lam_A.max():.1f}]")

    seq = sequence_ids(t[sc_lo:sc_hi], lat[sc_lo:sc_hi], lon[sc_lo:sc_hi], mag[sc_lo:sc_hi])
    n_seq = int(np.unique(seq).size)
    print(f"[seq] {n_seq} sequences among {n_score} scoring events")

    # three largest scoring sequences
    uq, cnts = np.unique(seq, return_counts=True)
    big3 = uq[np.argsort(-cnts)[:3]]
    big3_info = []
    for b3 in big3:
        sel = seq == b3
        k = np.argmax(mag[sc_lo:sc_hi] * sel)
        big3_info.append({"n_events": int(sel.sum()),
                          "largest_mag": float(mag[sc_lo:sc_hi][sel].max()),
                          "start": as_ts(float(t[sc_lo:sc_hi][sel].min()))})
    named = {"Landers_1992-06-28": pd.Timestamp("1992-06-28", tz="UTC"),
             "Northridge_1994-01-17": pd.Timestamp("1994-01-17", tz="UTC"),
             "HectorMine_1999-10-16": pd.Timestamp("1999-10-16", tz="UTC")}
    where_named = {}
    for nm, ts_ in named.items():
        tt = (ts_ - pd.Timestamp(0, tz="UTC")).total_seconds() / 86400.0
        where_named[nm] = ("scoring" if tt >= t_train_end else "training") if tt < t_explore_end \
            else "HOLDOUT (not read)"
    res["sequences"] = {"n_sequences_scoring": n_seq, "three_largest": big3_info,
                        "named_big_three_location": where_named}
    print(f"[seq] named: {where_named}")

    # ---------------- the declared family: T_obs in {1 h, 3 h} x {STAI on/off} ----------
    parent_all = mag >= M_PARENT - 1e-9
    parent_hi = mag >= M_PARENT_HIGH - 1e-9
    print("[attr] attributing events to parents (frozen crude rule) ...")
    owner_all = attribute(t, lat, lon, mag, np.where(parent_all)[0])
    owner_hi = attribute(t, lat, lon, mag, np.where(parent_hi)[0])
    print(f"[attr] attributed {int((owner_all>=0).sum())} of {t.size} exploration events "
          f"to {int(parent_all.sum())} M>=4 parents")
    readings = {}
    prior_cache = {}
    for T_h in T_OBS_HOURS:
        T_obs = T_h / 24.0
        for cm in COUNT_MODES:
            key = f"Tobs_{T_h:g}h_{cm}"
            upd = build_update(t, lat, lon, mag, (mu, K, alpha, c, p), parent_all,
                               float(t[0]), t_train_end, T_obs, cm, b_train,
                               prior=prior_cache.get("prior"), owner=owner_all)
            prior_cache["prior"] = upd["prior"]
            # parents that can affect the scoring window at all
            live = upd["start_t"] < T1_sc
            n_live = int(live.sum())
            s, lam_B, Lam_B = score_variant(sc, upd)
            assert np.all(np.isfinite(lam_B)) and lam_B.min() > 0
            LL_B = float(np.log(lam_B).sum() - Lam_B)
            bits_all = float(s.mean() / LN2)
            lo95, hi95, nb, se = block_bootstrap_ci(s, seq)
            # K-436 form: same total, per aftershock-attributed scoring event
            att = (upd["owner"][sc_lo:sc_hi] >= 0)
            n_att = int(att.sum())
            # K-436 form: the SAME total LL difference, per aftershock-attributed event.
            # Per-attributed-event score re-apportions the integral term over n_att.
            s_att = ((np.log(lam_B) - np.log(sc.lam_A))[att]
                     - (Lam_B - sc.Lam_A) / max(n_att, 1))
            bits_att = float(s_att.mean() / LN2)
            lo_a, hi_a, _, _ = block_bootstrap_ci(s_att, seq[att])
            entry = {
                "T_obs_hours": T_h, "count_mode": cm,
                "n_parents_M>=4_total": int(parent_all.sum()),
                "n_parents_live_for_scoring": n_live,
                "n_parents_updated_nontrivially": int((np.abs(upd["mult"] - 1) > 1e-6).sum()),
                "prior_log_productivity": {k: v for k, v in upd["prior"].items()
                                           if k != "train_residual_summary"},
                "prior_residual_summary": upd["prior"]["train_residual_summary"],
                "LL_A": LL_A, "LL_B": LL_B, "LL_diff_nats": LL_B - LL_A,
                "bits_per_event_all_scoring_events": bits_all,
                "ci95_bits_per_event": [lo95, hi95], "boot_se_bits": se,
                "n_scoring_events": int(n_score), "n_sequences": nb,
                "n_aftershock_attributed_scoring_events": n_att,
                "bits_per_attributed_event": bits_att,
                "ci95_bits_per_attributed_event": [lo_a, hi_a],
                "mult_stats": {"min": float(upd["mult"].min()), "median": float(np.median(upd["mult"])),
                               "max": float(upd["mult"].max())},
                "obs_count_stats_in_Tobs": {"mean": float(upd["n_obs"].mean()),
                                            "max": float(upd["n_obs"].max()),
                                            "frac_zero": float((upd["n_obs"] == 0).mean())},
            }
            readings[key] = entry
            print(f"[B] {key}: bits/event={bits_all:+.4f} CI[{lo95:+.4f},{hi95:+.4f}] "
                  f"bits/attributed={bits_att:+.4f} (n_att={n_att}, parents live={n_live})")

    res["primary"] = readings

    # least favourable member of the declared T_obs family of 2, within each counting mode
    def least_fav(cm):
        keys = [k for k in readings if readings[k]["count_mode"] == cm]
        return min(keys, key=lambda k: readings[k]["bits_per_event_all_scoring_events"])

    verdict = {"HEADLINE": "stai (Hainzl 2016; Merton prior-art note)"}
    for cm in COUNT_MODES:
        k = least_fav(cm)
        e = readings[k]
        verdict[cm] = {
            "least_favourable_key": k,
            "bits_per_event": e["bits_per_event_all_scoring_events"],
            "ci95": e["ci95_bits_per_event"],
            "PASS": bool(e["bits_per_event_all_scoring_events"] >= SUCCESS_BITS
                         and e["ci95_bits_per_event"][0] > 0),
        }
    res["success_rule_verdict"] = verdict

    # ---------------- secondary readings (on the primary T_obs = 1 h, both STAI states) ----
    sec = {}
    for T_h in T_OBS_HOURS:
        T_obs = T_h / 24.0
        for cm in COUNT_MODES:
            tag = f"Tobs_{T_h:g}h_{cm}"
            # (i) parents M >= 5 only
            u5 = build_update(t, lat, lon, mag, (mu, K, alpha, c, p), parent_hi,
                              float(t[0]), t_train_end, T_obs, cm, b_train,
                              owner=owner_hi)
            s5, lamB5, LamB5 = score_variant(sc, u5)
            lo5, hi5, _, _ = block_bootstrap_ci(s5, seq)
            # (ii) big-three-excluded (parents M>=4, drop scoring events in 3 largest sequences)
            upd = build_update(t, lat, lon, mag, (mu, K, alpha, c, p), parent_all,
                               float(t[0]), t_train_end, T_obs, cm, b_train,
                               prior=prior_cache["prior"], owner=owner_all)
            s_all, _, _ = score_variant(sc, upd)
            keep = ~np.isin(seq, big3)
            lo_b, hi_b, _, _ = block_bootstrap_ci(s_all[keep], seq[keep])
            # (iii) alpha-refit harsher baseline
            pr = upd["prior"]
            rm = np.exp(pr["alpha_refit_intercept"]
                        + pr["alpha_refit_slope_per_mag_natlog"] * (mag[upd["pidx"]] - M_PARENT))
            lam_A2, Lam_A2 = sc.add_multipliers(upd["pidx"], rm, t[upd["pidx"]])
            LL_A2 = float(np.log(lam_A2).sum() - Lam_A2)
            lam_B, Lam_B = sc.add_multipliers(upd["pidx"], upd["mult"], upd["start_t"])
            LL_B = float(np.log(lam_B).sum() - Lam_B)
            s_A2 = per_event_scores(sc.lam_A, lam_A2, sc.Lam_A, Lam_A2)
            s_BA2 = (np.log(lam_B) - np.log(lam_A2)) - (Lam_B - Lam_A2) / n_score
            lo_r, hi_r, _, _ = block_bootstrap_ci(s_BA2, seq)
            sec[tag] = {
                "parents_M>=5": {
                    "n_parents": int(parent_hi.sum()),
                    "bits_per_event": float(s5.mean() / LN2), "ci95": [lo5, hi5]},
                "big_three_sequences_excluded": {
                    "n_events_kept": int(keep.sum()), "n_events_dropped": int((~keep).sum()),
                    "bits_per_kept_event": float(s_all[keep].mean() / LN2), "ci95": [lo_b, hi_b]},
                "alpha_refit_baseline": {
                    "refit_slope_natlog_per_mag": pr["alpha_refit_slope_per_mag_natlog"],
                    "implied_alpha": pr["alpha_refit_implied_alpha"],
                    "bits_of_A2_over_A": float(s_A2.mean() / LN2),
                    "bits_of_B_over_A2": float(s_BA2.mean() / LN2),
                    "ci95_B_over_A2": [lo_r, hi_r],
                    "LL_A2": LL_A2, "LL_B": LL_B},
            }
            print(f"[sec] {tag}: M>=5 {sec[tag]['parents_M>=5']['bits_per_event']:+.4f} | "
                  f"big3-excl {sec[tag]['big_three_sequences_excluded']['bits_per_kept_event']:+.4f} | "
                  f"B over alpha-refit {sec[tag]['alpha_refit_baseline']['bits_of_B_over_A2']:+.4f}")
    res["secondary"] = sec

    # partial dump so a kill during the (long) sensitivity arm cannot lose the primary result
    res["runtime_minutes_at_partial_dump"] = round((time.time() - t_start) / 60.0, 2)
    OUT.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"[dump] primary + secondary written to {OUT.name}")

    # ---------------- sensitivity / calibration ----------------
    print("[sim] sensitivity arm ...")
    sim_theta = [mu, K, alpha, c, p]
    # a supercritical (or divergent, b <= alpha) fit cannot be simulated; K is rescaled to a
    # declared subcritical branching ratio and that is stated in the result.
    br_eff = float(K * c ** (1 - p) / (p - 1) * (b_train / max(b_train - alpha, 0.05))) \
        if p > 1.0 else float("inf")
    if not np.isfinite(br_eff) or br_eff >= SIM_TARGET_BRANCHING:
        sim_theta[1] = K * SIM_TARGET_BRANCHING / br_eff if np.isfinite(br_eff) else K * 0.01
    sim_note = {"K_rescaled": bool(sim_theta[1] != K), "K_sim": float(sim_theta[1]),
                "fitted_branching_ratio": br, "branching_ratio_used_for_rescale": br_eff,
                "target_branching": SIM_TARGET_BRANCHING}
    T_sim = SIM_YEARS * 365.25
    sim_out = {}
    rng = np.random.default_rng(SIM_SEED)
    t_sim_start = time.time()
    for sd in SIM_SDS:
        bits_list, sig_list, nev = [], [], []
        n_capped = 0
        for rep in range(SIM_REPS):
            ts_, ms_, las_, los_, meta = simulate_etas(
                tuple(sim_theta), b_train, T_sim,
                (LAT_MIN, LAT_MAX, LON_MIN, LON_MAX), rng, sd)
            if meta["capped"]:
                n_capped += 1
                continue
            if ts_.size < 500:
                continue
            tr_end = T_sim * TRAIN_FRAC_OF_EXPLORE
            lo_ = int(np.searchsorted(ts_, tr_end, side="left"))
            if lo_ < 50 or ts_.size - lo_ < 50:
                continue
            scS = Scorer(ts_, ms_, tuple(sim_theta), tr_end, float(ts_[-1]), lo_, ts_.size)
            updS = build_update(ts_, las_, los_, ms_, tuple(sim_theta),
                                ms_ >= M_PARENT - 1e-9, float(ts_[0]), tr_end,
                                1.0 / 24.0, "plain", b_train)
            sS, lamBS, _ = score_variant(scS, updS)
            if not np.all(np.isfinite(sS)):
                continue
            seqS = sequence_ids(ts_[lo_:], las_[lo_:], los_[lo_:], ms_[lo_:])
            loS, hiS, _, _ = block_bootstrap_ci(sS, seqS, n_boot=500, seed=BOOT_SEED + rep)
            bits = float(sS.mean() / LN2)
            bits_list.append(bits)
            sig_list.append(bool(bits >= SUCCESS_BITS and loS > 0))
            nev.append(int(sS.size))
        sim_out[f"sd_log10_{sd:g}"] = {
            "n_reps_used": len(bits_list), "n_reps_discarded_at_event_cap": n_capped,
            "mean_bits_per_event": float(np.mean(bits_list)) if bits_list else None,
            "sd_bits_per_event": float(np.std(bits_list, ddof=1)) if len(bits_list) > 1 else None,
            "bits_reps": [round(x, 5) for x in bits_list],
            "power_at_success_rule": float(np.mean(sig_list)) if sig_list else None,
            "mean_n_scoring_events": float(np.mean(nev)) if nev else None,
        }
        print(f"[sim] sd={sd:g}: bits={sim_out[f'sd_log10_{sd:g}']['mean_bits_per_event']} "
              f"power={sim_out[f'sd_log10_{sd:g}']['power_at_success_rule']} "
              f"({len(bits_list)} reps, {time.time()-t_sim_start:.0f}s)")
    powers = [(sd, sim_out[f"sd_log10_{sd:g}"]["power_at_success_rule"]) for sd in SIM_SDS]
    detect = [sd for sd, pw in powers if pw is not None and pw >= 0.8]
    res["sensitivity"] = {
        "generator": "temporal+spatial ETAS branching from the training fit; per-parent "
                     "(M>=4) productivity multiplied by 10^(sd*Z)/E[10^(sd*Z)] with Z ~ N(0,1) "
                     "truncated at +/-2 (MEAN-PRESERVING, so injecting scatter does not change "
                     "the branching ratio); offspring located at "
                     "parent + isotropic normal of sd 0.30*max(3L,10km); GR magnitudes at "
                     "b_train; pipeline re-run identically (attribution, prior, update, "
                     "sequence bootstrap) with T_obs = 1 h, no STAI (sims are perfectly "
                     "observed by construction)",
        "note": sim_note,
        "per_sd": sim_out,
        "min_sd_log10_detectable_at_80pct_power": (min(detect) if detect else None),
        "calibration_zero_gain_check": sim_out.get("sd_log10_0", {}).get("mean_bits_per_event"),
    }

    # THE READING THAT DECIDES HOW TO INTERPRET THE NULL: what gain would the simulation
    # predict at the scatter actually MEASURED in the training window?
    meas_sd = readings["Tobs_1h_stai"]["prior_residual_summary"]["sd_log10"]
    xs = [sd for sd in SIM_SDS if sim_out[f"sd_log10_{sd:g}"]["mean_bits_per_event"] is not None]
    ys = [sim_out[f"sd_log10_{sd:g}"]["mean_bits_per_event"] for sd in xs]
    res["sensitivity"]["measured_scatter_reading"] = {
        "measured_training_residual_sd_log10": meas_sd,
        "sim_predicted_bits_per_event_at_that_scatter": (
            float(np.interp(meas_sd, xs, ys)) if xs else None),
        "success_threshold_bits": SUCCESS_BITS,
        "interpretation": "if this interpolated value is below the frozen success threshold, "
                          "the arm was underpowered AT THE SCATTER THE DATA ACTUALLY CARRY, "
                          "and the null is a null of the success rule, not evidence that "
                          "per-sequence productivity scatter is absent",
    }
    res["sensitivity"]["power_note"] = (
        "power is the fraction of reps meeting BOTH halves of the frozen success rule "
        "(bits >= 0.02 AND sequence-block 95% CI excluding 0); where the point-estimate half "
        "is met but power is still 0, the CI half is the binding constraint, i.e. the "
        "sequence-block CI at this catalogue size cannot resolve an effect of that size")

    res["counted_invariants"] = {
        "n_exploration_events_M>=2.5": int(len(cat)),
        "n_train_events": int(n_train), "n_scoring_events": int(n_score),
        "n_parents_M>=4_exploration": int(parent_all.sum()),
        "n_parents_M>=5_exploration": int(parent_hi.sum()),
        "n_sequences_scoring": n_seq,
        "LL_A_finite": bool(np.isfinite(LL_A)),
        "all_lambda_positive": True,
        "no_nan_in_any_reading": bool(all(
            np.isfinite(v["bits_per_event_all_scoring_events"]) for v in readings.values())),
    }
    res["runtime_minutes"] = round((time.time() - t_start) / 60.0, 2)
    OUT.write_text(json.dumps(res, indent=2), encoding="utf-8")

    print("\n" + "=" * 74)
    for k, v in verdict.items():
        if not isinstance(v, dict):
            continue
        print(f"  {k:>15s}: least favourable {v['least_favourable_key']}  "
              f"{v['bits_per_event']:+.4f} bits/event  CI[{v['ci95'][0]:+.4f},{v['ci95'][1]:+.4f}]"
              f"  {'PASS' if v['PASS'] else 'FAIL'}")
    print(f"  runtime {res['runtime_minutes']} min -> {OUT.name}")
    print("=" * 74)


if __name__ == "__main__":
    main()
