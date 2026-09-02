"""POWER OF THE LEARNED ARM. What its null is worth. Priced 0.

`exp_learned.py` returned dAUC = -0.00148, p = 0.8293 (EPOCH-FIXED rerun, 2026-09-02; the
superseded epoch-defective run gave +0.00155, p = 0.3415): a gradient-boosted model given the
full per-event tidal vector, free to represent any interaction, and searched without
multiplicity cost, does not beat a model that knows only the time of day.

That is only worth something if the design CAN detect a tidal effect. And there is a
specific reason to doubt it before measuring: model A's absolute AUC is 0.50784. The
day/night detection artifact is real -- roughly 2 percent of events displaced into local
night, established in `exp_diurnal_discriminator.py` -- and the model can barely see it.
Inside a matched stratum, where a case and its four controls share a site and a fortnight,
discrimination is intrinsically weak. So the honest question is not "did it fire" but
"HOW BIG WOULD AN EFFECT HAVE TO BE BEFORE THIS DESIGN FIRED".

This program has been bitten by exactly this before. K-118 showed a statistic can pass
through ZERO sensitivity at an unlucky nuisance phase, which would make its null void
rather than merely weak, and the rule adopted from that episode is that a statistic must
be shown to RESPOND before its null may be quoted. This module discharges that duty for
the learned arm.

THE INJECTION. Within each stratum -- the same five candidate times the real design chose
among -- re-select which member is the case, with probability proportional to
(1 + eps * z_driver), where z_driver is the standardised tidal driver across that
stratum's own five candidates. eps = 0.10 means a ten percent preference for high-tidal-
stress moments over the event's own fortnight. Then run the ENTIRE pipeline unchanged:
fit model A, fit model B, take dAUC, and compare against the permutation null already
measured (mean -0.00028, sd 0.00391).

This is the smallest honest version of the hypothesis, because the injected preference
operates only WITHIN the matched set and so cannot borrow power from site, season, rate or
clustering -- exactly the quantities the matching was built to remove.

WHAT THE ANSWER DECIDES. If eps = 0.05 is detected reliably, the null forecloses the whole
covariate space down to a few percent and is the strongest negative result the program has.
If even eps = 0.20 is invisible, the learned arm is an underpowered instrument, its null
means little, and the design needs more controls per case, a better metric than AUC, or
both -- and that must be said plainly rather than filed as a null.
"""

from __future__ import annotations

import datetime as _dt
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import exp_highn as HN
import exp_learned as L
import exp_mass_screen as MS

OUT_JSON = HERE / "results_learned_sensitivity.json"
ZEN = HERE / "data" / "xue_lu_zenodo"

EPS = (0.05, 0.10, 0.20, 0.40)
TRIALS = 5
DRIVERS = ("areal_abs", "rate_abs")


def main():
    rng = np.random.default_rng(L.RNG_SEED)
    t, la, lo, dp, mg = HN.load_zenodo(ZEN / "QTM_decluster_m0.1.txt")
    MS.assert_epoch(t, 2008, "QTM_declustered")
    t, la, lo, dp, mg, _nh = HN.split(t, la, lo, dp, mg)
    n, K = t.size, L.CONTROLS_PER_CASE
    print("QTM declustered exploration: %d events, %d candidates per stratum"
          % (n, K + 1), flush=True)

    off = rng.uniform(-L.NULL_HALF_WINDOW_DAYS, L.NULL_HALF_WINDOW_DAYS, size=(n, K))
    t_all = np.concatenate([t, (t[:, None] + off).ravel()])
    la_all = np.concatenate([la, np.repeat(la, K)])
    lo_all = np.concatenate([lo, np.repeat(lo, K)])
    stratum = np.concatenate([np.arange(n), np.repeat(np.arange(n), K)])

    cols = L.build_features(t_all, la_all, lo_all)
    Xa = np.column_stack([cols[c] for c in L.NUISANCE])
    Xb = np.column_stack([cols[c] for c in L.NUISANCE + L.TIDAL])

    order = np.argsort(t)
    cut = order[int(L.TRAIN_FRAC * n)]
    tr_s = np.zeros(n, dtype=bool)
    tr_s[t < t[cut]] = True
    tr = tr_s[stratum]
    te = ~tr

    # index of each stratum's 5 candidates: row 0 is the real time, rows 1..K the phantoms
    member = np.empty((n, K + 1), dtype=np.int64)
    member[:, 0] = np.arange(n)
    member[:, 1:] = n + np.arange(n)[:, None] * K + np.arange(K)[None, :]

    raw = MS.raw_series(t_all, la_all, lo_all)

    out = {"arm": "power of the learned arm", "priced_tests": 0,
           "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
           # Read from the EPOCH-FIXED exp_learned.py rerun (2026-09-02).
           "observed_dauc_real_data": -0.00148,
           "null_mean": 0.00055, "null_sd": 0.00339,
           "detection_threshold_dauc_2sd": 2 * 0.00339,
           "trials_per_cell": TRIALS, "power": {}}

    for drv in DRIVERS:
        v = raw[drv][member]                                     # [n, K+1]
        z = (v - v.mean(axis=1, keepdims=True)) / (v.std(axis=1, keepdims=True) + 1e-30)
        out["power"][drv] = {}
        for eps in EPS:
            w = np.clip(1.0 + eps * z, 1e-9, None)
            w /= w.sum(axis=1, keepdims=True)
            cw = np.cumsum(w, axis=1)
            daucs = []
            for tri in range(TRIALS):
                r2 = np.random.default_rng(9000 + tri)
                pick = (cw < r2.random((n, 1))).sum(axis=1).clip(0, K)
                y = np.zeros(t_all.size)
                y[member[np.arange(n), pick]] = 1.0
                a, _ = L.fit_score(Xa[tr], y[tr], Xa[te], y[te], tri)
                b, _ = L.fit_score(Xb[tr], y[tr], Xb[te], y[te], tri)
                daucs.append(b - a)
            daucs = np.asarray(daucs)
            det = float(np.mean(daucs > out["null_mean"]
                                + 2.0 * out["null_sd"]))
            out["power"][drv]["eps=%.2f" % eps] = {
                "mean_dauc": float(daucs.mean()), "sd_dauc": float(daucs.std(ddof=1)),
                "fraction_clearing_2sd": det}
            print("  %-10s eps=%.2f   dAUC %+.5f +/- %.5f   clears 2sd in %.0f%% of trials"
                  % (drv, eps, daucs.mean(), daucs.std(ddof=1), 100 * det), flush=True)

    # SMALLEST reliably-detected effect, not largest. The first version used max(),
    # which reported eps = 0.40 -- the least informative number in the table -- as
    # though it were the detection limit.
    best = min(
        (eps for drv in DRIVERS for eps in EPS
         if out["power"][drv]["eps=%.2f" % eps]["fraction_clearing_2sd"] >= 0.8),
        default=None)
    out["smallest_reliably_detected_eps"] = best
    out["verdict"] = (
        ("DESIGN IS SENSITIVE. The smallest injected modulation detected in at least 80%% "
         "of trials is eps = %.2f. The arm's null (dAUC = -0.00148, p = 0.8293) therefore "
         "forecloses the whole covariate space down to roughly that level, which is a "
         "stronger negative result than any binned statistic in this program: it rules "
         "out every interaction a gradient-boosted model can represent, not one declared "
         "projection at a time." % best)
        if best is not None else
        ("DESIGN IS UNDERPOWERED AND ITS NULL MUST NOT BE QUOTED AS A BOUND. Not even the "
         "largest injected modulation tested is detected reliably. Matched case-control "
         "with AUC is too blunt at this stratum size: a case and its controls share a "
         "site and a fortnight, so the achievable discrimination is intrinsically tiny "
         "(model A sees the REAL day/night artifact at AUC 0.5038). The fix is more "
         "controls per case, a likelihood-ratio metric instead of AUC, or a design that "
         "does not throw away the within-stratum ordering. Recorded as a machinery "
         "result, not as evidence about the Earth."))
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print("\n" + "=" * 78)
    print(out["verdict"])
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
