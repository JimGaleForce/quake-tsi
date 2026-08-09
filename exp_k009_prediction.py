"""K-009 PRE-SCORING PREDICTION REGISTER (Popper ruling, round 2, binding).

Turns K-009 from a bare null test into a two-generator discrimination.  This script
writes results_k009_prediction.json and computes NOTHING about the residual field.

  Generator A -- W-003 (the null unifier, current champion).
      "P2 -- THE RESIDUALS ARE WHITE. K-009 must return spatial correlation length
       ~= 0 and temporal correlation time ~= 0 in the ETAS residual field, against
       the ETAS-sim null."
      => predicted residual correlation TIME  = 0 weeks (inside the sim-null envelope)
      => predicted residual correlation LENGTH = 0 km   (inside the sim-null envelope)

  Generator B -- W-001 / W-002 (rate-and-state with a latent local state).
      "W-001/W-002 both predict K-009 residuals are RED, with a correlation time
       equal to the local t_a -- months to years."
      t_a = A*sigma / stressing-rate is, in Dieterich (1994), exactly the time at which
      an aftershock sequence's rate returns to the background rate.  That is directly
      measurable, so generator B's prediction is quantitative and zero-free-parameter.

INPUT DISCIPLINE.  Everything here is computed on the TRAIN window only
(< 2010-01-01), the same window used to freeze EXP-H's temporal parameters and this
experiment's background field.  Nothing from the 2010-2018 residual window enters, and
no residual correlation statistic of any kind is computed in this file.

Two independent estimators of t_a are registered:
  (1) EMPIRICAL (primary).  Stack aftershock sequences of train-window mainshocks and
      find where the stacked rate crosses back through the pre-mainshock background
      rate.  This is Dieterich's t_a by definition.
  (2) ETAS-IMPLIED (secondary).  From the frozen EXP-H parameters alone, the time at
      which one mainshock's Omori triggering rate falls to the local background rate:
      K*10^(alpha (M-M0)) (t+c)^(-p) = mu_local  =>  t = (K*10^(a dM)/mu_local)^(1/p) - c.
      Note this one is magnitude-dependent, which the rate-and-state t_a is not; it is
      reported as a consistency check, not as t_a.

Run: python -u exp_k009_prediction.py
"""

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from exp_k009_residual_whiteness import (load_catalog, to_km, Grid, adaptive_background,
                                         CATALOG, LAT_MIN, LAT_MAX, LON_MIN, LON_MAX,
                                         T_TEST0, DT_DAYS, MIN_TRAIN_EVENTS_PER_CELL)

HERE = Path(__file__).resolve().parent
OUT = HERE / "results_k009_prediction.json"
EXPH = HERE / "results_exp_h.json"

MS_MMIN = 4.5           # mainshock selection floor
MS_RADIUS_KM = 25.0     # aftershock collection radius
MS_ISOLATION_DAYS = 365.0
MS_ISOLATION_KM = 100.0
MAX_T_DAYS = 4000.0
PRE_WIN_DAYS = (100.0, 3000.0)   # background estimated from -3000 .. -100 d


def main():
    t0 = time.time()
    exph = json.loads(EXPH.read_text())
    fp = exph["train_fit"]["frozen_params"]
    mu0, K, alpha, c, p, M0 = fp["mu"], fp["K"], fp["alpha"], fp["c"], fp["p"], fp["M0"]

    df, order = load_catalog(CATALOG)
    box = df.lat.between(LAT_MIN, LAT_MAX) & df.lon.between(LON_MIN, LON_MAX)
    d = df[box].sort_values("t").reset_index(drop=True)
    cat = d[d.mag >= M0 - 1e-9].reset_index(drop=True)
    t_split = (T_TEST0 - pd.Timestamp(0, tz="UTC")).total_seconds() / 86400.0
    tr = cat[cat.t < t_split].reset_index(drop=True)
    print(f"[data] {order}; M>={M0} train(<2010) = {len(tr)} events")

    t = tr.t.to_numpy(); m = tr.mag.to_numpy()
    x, y = to_km(tr.lat.to_numpy(), tr.lon.to_numpy())

    # ---------------- mainshock selection (isolated, train window, room for 4000 d) ----
    cand = np.nonzero((m >= MS_MMIN) & (t > t[0] + PRE_WIN_DAYS[1] + 30) &
                      (t < t_split - MAX_T_DAYS))[0]
    ms = []
    for i in cand:
        w = (np.abs(t - t[i]) <= MS_ISOLATION_DAYS) & (t != t[i])
        near = w & (np.sqrt((x - x[i]) ** 2 + (y - y[i]) ** 2) <= MS_ISOLATION_KM)
        if np.any(m[near] > m[i]):
            continue
        ms.append(i)
    ms = np.array(ms, dtype=int)
    print(f"[mainshocks] {len(ms)} isolated M>={MS_MMIN} train mainshocks "
          f"(no larger event within {MS_ISOLATION_KM} km / +-{MS_ISOLATION_DAYS} d)")

    # ---------------- stacked aftershock rate vs background ----------------------------
    edges = np.logspace(np.log10(0.01), np.log10(MAX_T_DAYS), 40)
    counts = np.zeros(len(edges) - 1)
    n_ms_used = 0
    bg_rate_sum = 0.0
    per_ms = []
    for i in ms:
        r = np.sqrt((x - x[i]) ** 2 + (y - y[i]) ** 2)
        near = r <= MS_RADIUS_KM
        dt = t - t[i]
        pre = near & (dt <= -PRE_WIN_DAYS[0]) & (dt >= -PRE_WIN_DAYS[1])
        bg = pre.sum() / (PRE_WIN_DAYS[1] - PRE_WIN_DAYS[0])     # events/day in the zone
        post = near & (dt > 0) & (dt <= MAX_T_DAYS)
        if post.sum() < 20:
            continue
        h, _ = np.histogram(dt[post], bins=edges)
        counts += h
        bg_rate_sum += bg
        n_ms_used += 1
        # per-mainshock t_a: first log-bin where the rate falls to the local background
        rate_i = h / np.diff(edges)
        ctr_i = np.sqrt(edges[:-1] * edges[1:])
        if bg > 0:
            above = np.nonzero((rate_i > bg) & (ctr_i > 1.0))[0]
            if len(above):
                # first bin after the last above-background bin, past 1 day
                k = above[-1]
                per_ms.append(float(ctr_i[min(k + 1, len(ctr_i) - 1)]))
    rate = counts / np.diff(edges) / max(n_ms_used, 1)
    bg_mean = bg_rate_sum / max(n_ms_used, 1)
    ctr = np.sqrt(edges[:-1] * edges[1:])
    below = np.nonzero((rate <= bg_mean) & (ctr > 1.0))[0]
    if len(below):
        j = below[0]
        # log-linear interpolation of the crossing between bins j-1 and j
        if j > 0 and rate[j - 1] > bg_mean > 0 and rate[j] > 0:
            f = (np.log(rate[j - 1]) - np.log(bg_mean)) / (np.log(rate[j - 1]) - np.log(rate[j]))
            t_a_emp = float(np.exp(np.log(ctr[j - 1]) + f * (np.log(ctr[j]) - np.log(ctr[j - 1]))))
        else:
            t_a_emp = float(ctr[j])
    else:
        t_a_emp = float("nan")
    print(f"[t_a empirical] stacked over {n_ms_used} mainshocks; zone background "
          f"{bg_mean:.4f} ev/d; stacked rate crosses background at t_a = {t_a_emp:.1f} d "
          f"= {t_a_emp/30.44:.2f} months = {t_a_emp/DT_DAYS:.1f} weeks")
    pm = np.array(per_ms, dtype=float)
    if len(pm):
        print(f"[t_a empirical] per-mainshock spread: n={len(pm)} median={np.median(pm):.0f} d "
              f"IQR {np.percentile(pm,25):.0f}-{np.percentile(pm,75):.0f} d")

    # ---------------- ETAS-implied return-to-background time ---------------------------
    train_lat, train_lon = tr.lat.to_numpy(), tr.lon.to_numpy()
    grid = Grid(train_lat, train_lon)
    p_bg = adaptive_background(train_lat, train_lon, grid.clat, grid.clon, k=4)
    # local background rate inside a MS_RADIUS_KM zone (cells within the radius, typical cell)
    cell_area = grid.area_km2
    zone_cells = max(1, int(np.pi * MS_RADIUS_KM ** 2 / cell_area))
    mu_zone_typ = mu0 * np.median(np.sort(p_bg)[-int(0.2 * grid.n_cells):]) * zone_cells
    etas_ta = {}
    for M in (4.5, 5.5, 6.5):
        A = K * 10.0 ** (alpha * (M - M0))
        etas_ta[f"M{M}"] = float((A / mu_zone_typ) ** (1.0 / p) - c)
    print(f"[t_a ETAS-implied] typical active-zone background {mu_zone_typ:.4f} ev/d; "
          + "  ".join(f"M{k[1:]}: {v:.0f} d" for k, v in etas_ta.items()))

    # ---------------- registered predictions -------------------------------------------
    T_pred_weeks = t_a_emp / DT_DAYS
    pred = {
        "experiment": "K-009 prediction register (pre-scoring)",
        "ruling": ("Popper round-2 binding ruling: K-009 must be a two-generator discrimination. "
                   "This file is written BEFORE the spec-window residual correlation statistics "
                   "are scored, and is to be publicly committed before scoring proceeds."),
        "state_class": "prediction-only (no residual statistic computed in this file)",
        "run_utc": pd.Timestamp.utcnow().isoformat(),
        "input_discipline": "train window only (< 2010-01-01); the 2010-2018 residual window is untouched",
        "frozen_params_source": "results_exp_h.json :: train_fit.frozen_params",
        "frozen_params": fp,

        "generator_A_W003": {
            "quote": ("W-003-P2: THE RESIDUALS ARE WHITE. K-009 must return spatial correlation "
                      "length ~= 0 and temporal correlation time ~= 0 in the ETAS residual field, "
                      "against the ETAS-sim null."),
            "predicted_correlation_TIME_weeks": 0.0,
            "predicted_correlation_LENGTH_km": 0.0,
            "operational": ("measured T and L both inside the ETAS-sim null 2.5-97.5 envelope, and "
                            "lag-1 weekly ACF excess < 0.05 over the null 97.5th percentile")},

        "generator_B_W001_W002": {
            "quote": ("W-001/W-002 both predict K-009 residuals are RED, with a correlation time "
                      "equal to the local t_a -- months to years. t_a = A*sigma / stressing-rate is "
                      "the time at which an aftershock sequence returns to background (Dieterich 1994)."),
            "predicted_correlation_TIME_days": t_a_emp,
            "predicted_correlation_TIME_weeks": T_pred_weeks,
            "predicted_correlation_TIME_months": t_a_emp / 30.44,
            "acceptance_band_factor": 2.0,
            "operational": ("measured T within a factor of 2 of the predicted t_a, i.e. "
                            f"{T_pred_weeks/2:.1f} - {T_pred_weeks*2:.1f} weeks "
                            f"({t_a_emp/2:.0f} - {t_a_emp*2:.0f} days), AND above the ETAS-sim "
                            "null 97.5th percentile"),
            "estimator_uncertainty": {
                "stacked_crossing_days": t_a_emp,
                "per_mainshock_median_days": float(np.median(pm)) if len(pm) else None,
                "per_mainshock_IQR_days": [float(np.percentile(pm, 25)), float(np.percentile(pm, 75))] if len(pm) else None,
                "note": ("the two empirical estimators disagree by roughly an order of magnitude: "
                         "the stacked crossing is robust but averages over heterogeneous zones, "
                         "while the per-mainshock 'last bin above background' is biased long by "
                         "single late events in wide log bins. Honest range for t_a in this "
                         "catalogue is ~200 d to ~5 yr, i.e. exactly W-001/W-002's 'months to "
                         "years'. The factor-2 acceptance band below is anchored on the stacked "
                         "estimator (the primary); a measurement falling between 27.6 and 250 "
                         "weeks is inside the full estimator range and must be reported as "
                         "'consistent with t_a but not sharply discriminating'."),
                "full_estimator_range_weeks": [t_a_emp / DT_DAYS,
                                               float(np.median(pm) / DT_DAYS) if len(pm) else None]},
            "predicted_correlation_LENGTH_km": None,
            "length_note": ("t_a fixes a TIME, not a LENGTH. W-001/W-002 do not supply a "
                            "zero-parameter length prediction; the K-009 hypothesis text's own "
                            "prior expectation (Moran's I excess at 20-60 km) is registered here "
                            "as a weaker, non-derived expectation and must not be scored as if it "
                            "were derived."),
            "length_prior_expectation_km": [20.0, 60.0]},

        "t_a_empirical": {
            "method": ("stacked aftershock rate of isolated train-window mainshocks vs the "
                       "pre-mainshock background rate in the same zone; t_a = crossing time"),
            "mainshock_Mmin": MS_MMIN, "radius_km": MS_RADIUS_KM,
            "isolation": {"days": MS_ISOLATION_DAYS, "km": MS_ISOLATION_KM},
            "background_window_days_before": list(PRE_WIN_DAYS),
            "n_mainshocks_selected": int(len(ms)), "n_mainshocks_stacked": int(n_ms_used),
            "zone_background_rate_per_day": float(bg_mean),
            "t_a_days": t_a_emp, "t_a_weeks": T_pred_weeks, "t_a_months": t_a_emp / 30.44,
            "per_mainshock_t_a_days": {
                "n": int(len(pm)),
                "median": float(np.median(pm)) if len(pm) else None,
                "p25": float(np.percentile(pm, 25)) if len(pm) else None,
                "p75": float(np.percentile(pm, 75)) if len(pm) else None},
            "stacked_curve": {"bin_center_days": ctr.tolist(),
                              "rate_per_day_per_mainshock": rate.tolist(),
                              "background_per_day": float(bg_mean)}},

        "t_a_etas_implied": {
            "method": ("time at which one mainshock's frozen-ETAS Omori triggering rate equals the "
                       "typical active-zone background rate: (K 10^(alpha dM)/mu_zone)^(1/p) - c"),
            "mu_zone_per_day": float(mu_zone_typ), "zone_cells": int(zone_cells),
            "t_days_by_magnitude": etas_ta,
            "caveat": ("magnitude-dependent, so it is NOT the rate-and-state t_a; registered as a "
                       "consistency check on the order of magnitude only")},

        "scoring_plan": {
            "outcome_W003_wins": "T and L inside the ETAS-sim null envelope",
            "outcome_W001_W002_wins": ("T above the null 97.5th percentile AND within a factor of 2 "
                                       f"of {T_pred_weeks:.1f} weeks"),
            "outcome_neither": ("T above the null envelope but outside the factor-2 band around "
                                "the predicted t_a -- residuals are red but not at the "
                                "rate-and-state timescale; both generators lose and the number "
                                "is reported as a bare measurement"),
            "note": ("Popper's own K-009 success rule still gates any 'there is weather' claim: "
                     "lag-1 ACF excess >= 0.05 over the sim-null 97.5th percentile, stable across "
                     "the kernel swap, and surviving the rho_sta partial.")},

        "ORDERING_DISCLOSURE": {
            "status": "PREDICTION POST-REGISTERED (ordering violation, self-reported)",
            "what_happened": (
                "Two smoke-test runs of exp_k009_residual_whiteness.py were executed BEFORE this "
                "ruling arrived, on a REDUCED 2-year sub-window (2010-01-01 .. 2012-01-01, not the "
                "spec's 2010-2018) with 2 simulations, to debug the pipeline. Those runs printed "
                "measured residual correlation statistics, which the executing agent has therefore "
                "seen. The spec-window (2010-2018) scoring has NOT been run and results_k009.json "
                "does not exist."),
            "what_was_seen": (
                "2010-2012 sub-window, adaptive-k4 background: lag-1 pooled ACF approx +0.16 to "
                "+0.18, exponential-fit correlation time approx 9.8-10.2 weeks, correlation length "
                "approx 27-29 km; 2-sim null lag-1 ACF approx -0.015 (post-bugfix run). These are "
                "debug-run numbers on a quarter-length window with a 2-sample null and are not "
                "results."),
            "why_it_matters": (
                "The predicted t_a in this file is derived only from train-window (<2010) "
                "aftershock stacking and the frozen EXP-H parameters -- no residual statistic "
                "enters its computation -- but the executing agent was not blind when writing it. "
                "The prediction is therefore POST-registered, not pre-registered, and the result "
                "must carry that caveat."),
            "not_backdated": True},
    }
    OUT.write_text(json.dumps(pred, indent=2, default=float))
    print(f"\n[prediction registered] -> {OUT.name}  ({time.time()-t0:.0f}s)")
    print("=" * 78)
    print(f"  W-003        predicts  T = 0 weeks, L = 0 km (white)")
    print(f"  W-001/W-002  predicts  T = t_a = {t_a_emp:.0f} d = {t_a_emp/30.44:.2f} months "
          f"= {T_pred_weeks:.1f} weeks  (accept within factor 2: "
          f"{T_pred_weeks/2:.1f}-{T_pred_weeks*2:.1f} wk)")
    print(f"  ORDERING: prediction POST-registered (smoke-run statistics already seen) -- disclosed")
    print("=" * 78)


if __name__ == "__main__":
    main()
