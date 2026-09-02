"""POWER OF THE EXTENDED LEARNED ARM. What the block-3 null is worth. Priced 0.

STANDING RULE 5: a statistic must be shown to RESPOND before its null may be quoted.
`exp_learned_ext.py` reports dAUC(C vs A) for a non-tidal forcing block against a
within-stratum permutation null. That number is worth nothing until this module has
measured how big an effect would have to be before the design saw it. K-118 is the
precedent: a statistic can pass through ZERO sensitivity at an unlucky nuisance phase,
which makes its null void rather than merely weak.

THE INJECTION, identical in form to `exp_learned_sensitivity.py` so the two numbers are
comparable. Within each stratum -- the same K+1 candidate times the real design chose
among -- re-select which member is the case with probability proportional to
(1 + eps * s * z_driver), where z_driver is the driver standardised across that stratum's
OWN candidates and s = +/-1. eps = 0.10 is a ten percent preference for high-driver
moments over the event's own fortnight. Then the ENTIRE pipeline runs unchanged: fit model
A, fit model C, take dAUC, compare against the permutation null measured by the main
module.

The injection lives only WITHIN the matched set, so it cannot borrow power from site,
season, rate or clustering -- exactly the quantities the matching was built to remove.
This is the smallest honest version of the hypothesis.

LEAST FAVOURABLE PHASE, which rule 5 requires and which the committed tidal sensitivity
did not do. Block-3 drivers are not sinusoids, so "phase" here is not an angle: the
nuisance degrees of freedom that can help or hurt are the SIGN of the response and the
CHOICE of driver. Both are swept, and the headline number is the MINIMUM over the sweep --
the smallest eps detected in at least 80 percent of trials AT THE WORST cell -- not the
best cell and not the average.

DRIVERS. `lod_d1` is the one-day change in length of day: the mechanical channel, and the
only part of the LOD signal this design can see at all (the matched fortnight differences
away everything slower). `kp` is the placebo channel, quantised to thirds, whose within-
stratum spread is frequently zero -- strata with no spread receive no injected preference
at all, which is itself a real and reportable limit on what the placebo arm can detect.

WHAT THE ANSWER DECIDES. If a few-percent modulation is detected reliably at the WORST
cell, the block-3 null forecloses that covariate class down to a few percent, and the
P-2.4 survivor rate in the main module is a calibrated false-positive rate rather than a
curiosity. If even the largest injected effect is invisible, the extended arm is an
underpowered instrument, its null must not be quoted as a bound, and that must be said
plainly rather than filed as a null.

STATE CLASS: first-run, exploration split. The holdout is not touched.
"""

from __future__ import annotations

import datetime as _dt
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import engine_ext_forcing as F
import exp_highn as HN
import exp_learned as L
import exp_learned_ext as X
import exp_world_harmonics as W

OUT_JSON = HERE / "results_learned_ext_sensitivity.json"
MAIN_JSON = HERE / "results_learned_ext.json"
ZEN = HERE / "data" / "xue_lu_zenodo"

EPS = (0.05, 0.10, 0.20, 0.40)
SIGNS = (+1, -1)
TRIALS = 10                      # 80% power is exactly 8/10 trials clearing threshold
DRIVERS = ("lod_d1", "kp")


def main():
    t_start = time.time()
    if not MAIN_JSON.exists():
        raise SystemExit("run exp_learned_ext.py first: %s is missing" % MAIN_JSON.name)
    main_out = json.loads(MAIN_JSON.read_text(encoding="utf-8"))
    c_null = main_out["contrasts"]["C_vs_A"]
    null_mean, null_sd = c_null["null_mean"], c_null["null_sd"]
    thr = null_mean + 2.0 * null_sd
    observed = c_null["observed_dAUC"]
    print("block-3 null (C vs A): mean %+.6f sd %.6f  -> 2sd threshold %+.6f "
          "(observed on real data %+.6f)" % (null_mean, null_sd, thr, observed),
          flush=True)

    rng = np.random.default_rng(L.RNG_SEED)
    t_raw, la_raw, lo_raw, dp_raw, mg_raw = HN.load_zenodo(ZEN / "QTM_decluster_m0.1.txt")
    t, la, lo, dp, mg, _nh = HN.split(t_raw, la_raw, lo_raw, dp_raw, mg_raw)
    n, K = t.size, L.CONTROLS_PER_CASE
    if n != X.EXPECT["n_cases_explore"]:
        raise SystemExit("counted invariant failed: %d cases, expected %d"
                         % (n, X.EXPECT["n_cases_explore"]))
    print("QTM declustered exploration: %d events, %d candidates per stratum"
          % (n, K + 1), flush=True)

    # the SAME strata as the main module: same seed, same draw order
    off = rng.uniform(-L.NULL_HALF_WINDOW_DAYS, L.NULL_HALF_WINDOW_DAYS, size=(n, K))
    t_all = np.concatenate([t, (t[:, None] + off).ravel()])
    la_all = np.concatenate([la, np.repeat(la, K)])
    lo_all = np.concatenate([lo, np.repeat(lo, K)])
    stratum = np.concatenate([np.arange(n), np.repeat(np.arange(n), K)])
    member = np.empty((n, K + 1), dtype=np.int64)
    member[:, 0] = np.arange(n)
    member[:, 1:] = n + np.arange(n)[:, None] * K + np.arange(K)[None, :]

    cols = L.build_features(t_all, la_all, lo_all)
    series, load_audit = F.load_forcing()
    # EPOCH FIX (supervisor, 2026-09-02): HN.load_zenodo now returns days since
    # 1970-01-01Z, so the explicit epoch here is the Unix epoch, not SPAN_START.
    jd_true = F.jd_from_days_since(t_all, _dt.datetime(1970, 1, 1,
                                                       tzinfo=_dt.timezone.utc))
    b3, fill = F.event_block(jd_true, series)
    allcols = dict(cols)
    allcols.update(b3)

    Xa = np.column_stack([allcols[c] for c in L.NUISANCE])
    Xc = np.column_stack([allcols[c] for c in L.NUISANCE + list(F.BLOCK3)])

    order = np.argsort(t)
    cut = order[int(L.TRAIN_FRAC * n)]
    tr_s = np.zeros(n, dtype=bool)
    tr_s[t < t[cut]] = True
    tr = tr_s[stratum]
    te = ~tr
    print("  train rows %d   test rows %d" % (tr.sum(), te.sum()), flush=True)

    out = {
        "arm": "power of the extended (block-3) learned arm",
        "state_class": "first-run, exploration split",
        "priced_tests": 0, "holdout_touched": False,
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "contrast_under_test": "C vs A (nuisance + block3 vs nuisance)",
        "observed_dauc_real_data": observed,
        "null_mean": null_mean, "null_sd": null_sd,
        "detection_threshold_dauc_2sd": thr,
        "n_permutations_behind_null": c_null["n_permutations"],
        "trials_per_cell": TRIALS, "eps_grid": list(EPS), "signs": list(SIGNS),
        "drivers": list(DRIVERS),
        "block3_fill_report": fill,
        "forcing_rows": {k: v["n_rows"] for k, v in load_audit["sources"].items()},
        "power": {},
    }

    for drv in DRIVERS:
        z, usable = X.stratum_z(allcols[drv], member)
        out["power"][drv] = {"n_usable_strata": int(usable.sum()),
                             "n_strata": int(n),
                             "usable_fraction": float(usable.mean()),
                             "cells": {}}
        print("\ndriver %s: %d/%d strata have non-degenerate within-stratum spread"
              % (drv, usable.sum(), n), flush=True)
        for sgn in SIGNS:
            for eps in EPS:
                w = np.clip(1.0 + eps * sgn * z, 1e-9, None)
                w /= w.sum(axis=1, keepdims=True)
                cw = np.cumsum(w, axis=1)
                daucs = []
                for tri in range(TRIALS):
                    r2 = np.random.default_rng(9000 + tri)
                    pick = (cw < r2.random((n, 1))).sum(axis=1).clip(0, K)
                    yv = np.zeros(t_all.size)
                    yv[member[np.arange(n), pick]] = 1.0
                    a, _ = L.fit_score(Xa[tr], yv[tr], Xa[te], yv[te], tri)
                    c, _ = L.fit_score(Xc[tr], yv[tr], Xc[te], yv[te], tri)
                    daucs.append(c - a)
                daucs = np.asarray(daucs)
                det = float(np.mean(daucs > thr))
                key = "sign=%+d,eps=%.2f" % (sgn, eps)
                out["power"][drv]["cells"][key] = {
                    "sign": int(sgn), "eps": float(eps),
                    "mean_dauc": float(daucs.mean()),
                    "sd_dauc": float(daucs.std(ddof=1)),
                    "fraction_clearing_2sd": det}
                print("  %-8s %-16s  dAUC %+.5f +/- %.5f   clears 2sd in %3.0f%% of trials"
                      % (drv, key, daucs.mean(), daucs.std(ddof=1), 100 * det),
                      flush=True)

    # LEAST FAVOURABLE CELL: for each eps, the WORST power over driver and sign. The
    # reported minimum detectable effect is the smallest eps whose WORST cell still
    # reaches 80 percent, which is what rule 5's "minimum over nuisance phase" means here.
    worst = {}
    for eps in EPS:
        vals = []
        for drv in DRIVERS:
            for sgn in SIGNS:
                vals.append((out["power"][drv]["cells"]["sign=%+d,eps=%.2f" % (sgn, eps)]
                             ["fraction_clearing_2sd"], drv, sgn))
        v = min(vals)
        worst["eps=%.2f" % eps] = {"worst_power": v[0], "worst_driver": v[1],
                                   "worst_sign": int(v[2])}
    out["least_favourable_cell_by_eps"] = worst
    mde_worst = min((eps for eps in EPS
                     if worst["eps=%.2f" % eps]["worst_power"] >= 0.8), default=None)
    mde_best = min((eps for eps in EPS for drv in DRIVERS for sgn in SIGNS
                    if out["power"][drv]["cells"]["sign=%+d,eps=%.2f" % (sgn, eps)]
                    ["fraction_clearing_2sd"] >= 0.8), default=None)
    out["min_detectable_eps_at_least_favourable_phase"] = mde_worst
    out["min_detectable_eps_at_most_favourable_cell"] = mde_best

    if mde_worst is not None:
        out["verdict"] = (
            "DESIGN IS SENSITIVE ON BLOCK 3. At the LEAST FAVOURABLE cell (driver and sign "
            "swept), the smallest injected within-stratum modulation detected in at least "
            "80%% of trials is eps = %.2f. The observed dAUC(C vs A) on real data is "
            "%+.5f against a null of %+.5f +/- %.5f, so the block-3 null forecloses that "
            "covariate class down to roughly a %.0f%% within-fortnight modulation -- for "
            "every interaction a gradient-boosted model can represent, not one declared "
            "projection at a time. It forecloses NOTHING at periods longer than the "
            "matched fortnight, which is where the decadal LOD claim lives."
            % (mde_worst, observed, null_mean, null_sd, 100 * mde_worst))
    else:
        out["verdict"] = (
            "UNDERPOWERED AT THE LEAST FAVOURABLE PHASE; THE BLOCK-3 NULL MUST NOT BE "
            "QUOTED AS A BOUND WITHOUT THIS CAVEAT. No eps on the grid is detected in at "
            "least 80%% of trials at its worst driver/sign cell (best cell reaches 80%% at "
            "eps = %s). Recorded as a machinery result, not as evidence about the Earth. "
            "The likely cause is the same one the tidal arm has: inside a matched stratum "
            "a case and its four controls share a site and a fortnight, so achievable "
            "discrimination is intrinsically tiny, and for `kp` it is worse still because "
            "the index is quantised and many strata have zero within-stratum spread."
            % (("%.2f" % mde_best) if mde_best is not None else "none"))

    out["runtime_seconds"] = time.time() - t_start
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print("\n" + "=" * 78)
    print(out["verdict"])
    print("\nruntime %.1f s   wrote %s" % (out["runtime_seconds"], OUT_JSON.name))


if __name__ == "__main__":
    main()
