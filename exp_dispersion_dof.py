"""IS THE QUIET REAL? The dispersion gate's own d.o.f., measured rather than assumed. Priced 0.

Jim asked whether zero survivors across the whole program is suspiciously quiet. It is
exactly the right question, and this program already had the measurement that answers it
and did not follow it up: `exp_dispersion_gate.py` found the declustered per-region z's
have VARIANCE 0.719 against an expected 1.000 and called it p = 6e-5, then filed it as
"worth chasing separately".

If that number is real it is the most important fact in the program, because it does not
mean "no signal". It means OUR NULLS ARE TOO WIDE -- every null standard deviation
inflated by about 18 percent, every z shrunk by the same, every test under-powered and
every published bound looser than the data supports. Manufacturing quiet is a far more
serious failure than finding nothing.

BUT THE GATE'S p IS NOT TRUSTWORTHY AS COMPUTED, and the defect is in the supervisor's own
code. The chi-square there treats the 280 region-statistic z's as INDEPENDENT. They are
obviously not: R_1 through R_4 of one angle are four views of one series, the level flag
and the quadrant flag share the same areal strain, and every region is scored against the
same ephemeris over overlapping epochs. Correlation does not bias the variance ESTIMATE,
but it collapses the degrees of freedom, and the entire p depends on the d.o.f.

For z ~ N(0, R) with R the correlation matrix over cells, the sum of squares has

    E[sum z^2] = p        and       Var[sum z^2] = 2 * trace(R^2)

so the effective d.o.f. is p^2 / trace(R^2) -- equal to p only when R is the identity.
This module MEASURES R from the null ensemble itself (the null columns are draws of the
cell z's under H0, so their sample correlation IS an estimate of R), and redoes the test
with the honest d.o.f.

THE ANSWER DECIDES WHERE TO LOOK NEXT. If the under-dispersion survives, stop hunting
physics and go fix the null construction, because everything measured so far is measured
with a mis-calibrated instrument. If it does not survive, the quiet is ordinary and the
gate's own alarming p was an artifact of assuming independence -- which is worth knowing
with equal force, since that p is currently in a committed artifact.
"""

from __future__ import annotations

import datetime as _dt
import glob
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import exp_world_harmonics as W

OUT_JSON = HERE / "results_dispersion_dof.json"


def collect(declustered):
    rng = np.random.default_rng(W.RNG_SEED)
    zs, cols, labels = [], [], []
    for path in sorted(glob.glob(str(HERE / "data" / "comcat_world" / "*.csv"))):
        name = os.path.basename(path)[:-4]
        t, la, lo, _h, _f = W.load_region(path, declustered=declustered)
        r = W.run_region(name, t, la, lo, rng, verbose=False)
        if r is None:
            continue
        for k in W.STATS:
            z = r["per_statistic"][k]["z"]
            if not np.isfinite(z):
                continue
            zs.append(float(z))
            cols.append(r["_null_cols"][k])
            labels.append("%s/%s" % (name, k))
    return np.asarray(zs), np.stack(cols), labels


def analyse(zs, cols, label):
    p = zs.size
    S = float((zs * zs).sum())

    # R estimated from the null ensemble: each column is a draw of that cell's z under H0
    C = np.corrcoef(cols)
    C = np.nan_to_num(C, nan=0.0)
    trR2 = float((C * C).sum())
    dof_eff = p * p / trR2

    var_obs = float(zs.var(ddof=1))
    # honest test: S ~ mean p, sd sqrt(2 trR2)
    sd_S = np.sqrt(2.0 * trR2)
    zscore = (S - p) / sd_S
    # eigen spectrum -> a second, independent read on redundancy
    ev = np.linalg.eigvalsh(C)
    ev = np.clip(ev, 0.0, None)
    eff_rank = float((ev.sum() ** 2) / (ev * ev).sum())

    return {
        "arm": label, "n_cells": int(p),
        "sum_of_squares": S,
        "variance_of_z": var_obs,
        "naive_dof_assumed_by_the_gate": int(p),
        "trace_R_squared": trR2,
        "effective_dof": dof_eff,
        "participation_ratio_effective_rank": eff_rank,
        "z_of_sum_of_squares": zscore,
        "two_sided_p_honest": float(math.erfc(abs(zscore) / math.sqrt(2.0))),
        "reading": ("UNDER-dispersed beyond chance" if zscore < -1.96 else
                    ("OVER-dispersed beyond chance" if zscore > 1.96 else
                     "consistent with a correctly calibrated null once the "
                     "correlation between cells is accounted for")),
    }


def main():
    out = {"arm": "dispersion gate d.o.f., measured not assumed",
           "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
           "priced_tests": 0,
           "why_priced_zero": ("a re-analysis of an already-computed quantity with the "
                               "correct degrees of freedom; no new statistic, no new "
                               "catalogue read"),
           "by_arm": {}}
    for label, dec in (("declustered_PRIMARY", True), ("full_SECONDARY", False)):
        print("collecting %s ..." % label, flush=True)
        zs, cols, _lab = collect(dec)
        rec = analyse(zs, cols, label)
        out["by_arm"][label] = rec
        print("  cells %d   var(z) %.3f   naive dof %d   EFFECTIVE dof %.1f"
              % (rec["n_cells"], rec["variance_of_z"],
                 rec["naive_dof_assumed_by_the_gate"], rec["effective_dof"]), flush=True)
        print("  sum z^2 = %.1f   z = %+.2f   honest p = %.4f   -> %s"
              % (rec["sum_of_squares"], rec["z_of_sum_of_squares"],
                 rec["two_sided_p_honest"], rec["reading"]), flush=True)

    d = out["by_arm"]["declustered_PRIMARY"]
    out["verdict"] = {
        "gate_reported": ("variance 0.719 with p_under = 6e-5, computed as if the 280 "
                          "region-statistic z's were independent"),
        "measured_effective_dof": d["effective_dof"],
        "honest_p": d["two_sided_p_honest"],
        "consequence": d["reading"],
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
