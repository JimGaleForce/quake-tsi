"""CALIBRATED POOLED-z MAX-STATISTIC. The test the world scan never ran. Priced 0.

Fable audit finding 5, and it is a genuine POWER GAP in a committed result rather than a
bug.

THE GAP. Every world-scan verdict is `P(max |z| over 13 regions x 16 statistics)`. That
is the right correction for an effect concentrated in ONE region-statistic cell. It is
the WRONG one for a weak effect present in EVERY region in the SAME statistic, which is
exactly what a real global mechanism would look like:

    a true effect at z = 1.0 in each of 13 regions
      -> pooled Stouffer z = 3.6         (strongly significant)
      -> cell max |z|      = ~2.5        (invisible; the null 95th is ~3.5-3.7)

So the declared test is blind to the most physically plausible alternative. The scan
computed Stouffer pooled z and Cochran's Q, printed them, and never gave them a
calibrated family p. This module supplies exactly that and nothing else.

THE FIX. The per-region null columns are already standardised per statistic and drawn
INDEPENDENTLY across regions at each replicate index -- and independence is precisely
what a pooled sum requires, so the pooled z can be recomputed inside each replicate:

    z_pooled[stat, rep] = sum_r w_r z_r[stat, rep] / sqrt(sum_r w_r^2),  w_r = sqrt(n_r)

and the family correction is `P(max over STATISTICS of |z_pooled|)`. That corrects over
16 statistics instead of 208 cells, which is the correct denominator for this
alternative and is a far lower bar.

WHAT THIS IS NOT. It is not a new statistic and not a new declaration: the pooled z was
already computed and reported in both committed artifacts. What was missing was its
null. So this is priced 0 -- the same quantity, finally calibrated.

AND IT CANNOT MANUFACTURE A SURVIVOR QUIETLY. Both the cell-max and the pooled results
are reported side by side, for both the declustered PRIMARY and the full SECONDARY arm.
The full arm is known to be contaminated by aftershock dependence -- its pooled bearing
numbers sit near +2.3 with Cochran Q above 0.05, which is precisely the shape this test
is built to detect and precisely the shape dependence also produces. If the pooled test
fires in the full arm and not the declustered one, that is dependence again, not a
finding.
"""

from __future__ import annotations

import datetime as _dt
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import exp_world_harmonics as W

OUT_JSON = HERE / "results_world_pooled_calibration.json"


def run_arm(declustered):
    rng = np.random.default_rng(W.RNG_SEED)
    regions = {}
    for path in sorted(glob.glob(str(HERE / "data" / "comcat_world" / "*.csv"))):
        name = os.path.basename(path)[:-4]
        t, la, lo, _h, _f = W.load_region(path, declustered=declustered)
        r = W.run_region(name, t, la, lo, rng, verbose=False)
        if r is not None:
            regions[name] = r
    names = sorted(regions)
    w = np.array([np.sqrt(regions[r]["n"]) for r in names])
    denom = np.sqrt((w * w).sum())

    obs, null = {}, {}
    for k in W.STATS:
        zo = np.array([regions[r]["per_statistic"][k]["z"] for r in names])
        obs[k] = float((w * zo).sum() / denom)
        zn = np.stack([regions[r]["_null_cols"][k] for r in names])   # [R_regions, rep]
        null[k] = (w[:, None] * zn).sum(axis=0) / denom
    for r in regions:
        regions[r].pop("_null_cols", None)

    nullmat = np.stack([np.abs(null[k]) for k in W.STATS])
    gw = nullmat.max(axis=0)
    obs_max = max(abs(obs[k]) for k in W.STATS)
    p = float((np.sum(gw >= obs_max) + 1) / (W.N_REPLICATES + 1))
    top = max(W.STATS, key=lambda k: abs(obs[k]))
    per = {k: {"pooled_z_observed": obs[k],
               "pooled_z_null_sd": float(null[k].std(ddof=1)),
               "exact_p_this_statistic":
                   float((np.sum(np.abs(null[k]) >= abs(obs[k])) + 1)
                         / (W.N_REPLICATES + 1))}
           for k in W.STATS}
    return {
        "declustered": declustered,
        "n_regions": len(names), "regions": names,
        "n_events": int(sum(regions[r]["n"] for r in names)),
        "weights": "sqrt(n) per region",
        "per_statistic": per,
        "max_statistic_over_STATISTICS": {
            "observed_max_abs_pooled_z": float(obs_max),
            "null_max_p95": float(np.quantile(gw, 0.95)),
            "p": p,
            "where": top,
            "n_corrected_over": len(W.STATS),
        },
    }


def main():
    out = {
        "arm": "calibrated pooled-z max-statistic (Fable audit finding 5)",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 0,
        "why_priced_zero": ("the pooled z was already computed and reported in both "
                            "committed artifacts; what was missing was its null. Same "
                            "quantity, finally calibrated. No statistic is added."),
        "the_gap": ("the declared verdict was P(max |z| over 208 region-statistic "
                    "CELLS), which is blind to a weak effect present in every region "
                    "in the SAME statistic: z = 1.0 in each of 13 regions pools to 3.6 "
                    "but has a cell max of only ~2.5 against a null 95th of ~3.5."),
    }
    for label, dec in (("declustered_PRIMARY", True), ("full_SECONDARY", False)):
        print("running %s ..." % label, flush=True)
        out[label] = run_arm(dec)
        m = out[label]["max_statistic_over_STATISTICS"]
        print("  pooled max |z| = %.3f at %s ; null 95th %.3f ; p = %.4f"
              % (m["observed_max_abs_pooled_z"], m["where"], m["null_max_p95"],
                 m["p"]), flush=True)

    dp = out["declustered_PRIMARY"]["max_statistic_over_STATISTICS"]["p"]
    fp = out["full_SECONDARY"]["max_statistic_over_STATISTICS"]["p"]
    out["verdict"] = {
        "declustered_primary_p": dp, "full_secondary_p": fp,
        "reading": ("The PRIMARY arm governs. If the full arm fires and the "
                    "declustered one does not, that is aftershock dependence again -- "
                    "the same thing that killed the cell-max survivor -- and not a "
                    "finding."),
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print("\n" + "=" * 74)
    print("POOLED CALIBRATION: declustered PRIMARY p = %.4f ; full SECONDARY p = %.4f"
          % (dp, fp))
    print("  (cell-max on the same declustered arm gave p = 0.9876)")
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
