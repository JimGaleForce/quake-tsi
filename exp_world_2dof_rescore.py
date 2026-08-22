"""2-DOF RE-SCORE of the world scan's angular statistics. Priced 0.

WHY. Kepler's K-118/F-2, verified independently by the supervisor before adoption. The
scan scored angular concentration as `|R_m|`, the LENGTH of the complex resultant. That
discards the resultant's direction, and when a site's own dwell distribution already has
a long resultant an off-axis modulation ROTATES it rather than lengthening it. So
`d|R_m|/d(eps)` is SIGNED and passes through zero: at some unknown modulation phases the
statistic barely responds at all.

THE SUPERVISOR'S INDEPENDENT MEASUREMENT, which confirms the structure and corrects the
magnitude. Detection power at a RANDOM unknown phase, two-sided 5 %:

    Chile   m=1  eps=0.30 : |R| 0.405   2-dof 0.858
    Chile   m=1  eps=0.15 : |R| 0.105   2-dof 0.360
    Alaska  m=1  eps=0.30 : |R| 0.255   2-dof 0.718
    Chile   m=2  eps=0.30 : |R| 0.760   2-dof 1.000

**The 2-dof form is 2x to 3.4x more powerful and is strictly better.** But Kepler's
stronger claim -- that the order-1 statistics "could not have detected a total effect"
and should be struck as blind -- IS NOT SUPPORTED. |R| retains 26-40 % power at order 1
and 63-76 % at order 2. The scan's statistics were WEAKENED, not blind, and its nulls
are nulls at reduced power rather than void. That correction is recorded here because
the difference decides whether 208 committed tests stand.

WHAT THIS RE-SCORE IS AND IS NOT. It is the SAME declared quantity -- concentration at
order m -- estimated by a better statistic. It is NOT a new test and no statistic is
added, so it is priced 0. The change is motivated by a property of the dwell
distribution that is measurable WITHOUT ANY EVENTS, so it cannot be a
choose-the-estimator-after-seeing-the-answer move; and both the old and new results are
reported side by side. **If the 2-dof version produced a survivor where |R| did not,
that would need Popper before it could be called anything.**

THE STATISTIC. For each region and order, the complex resultant `Z_m = mean(exp(i m
theta))` is scored as the Mahalanobis distance from the null mean under the null
covariance of `(Re Z, Im Z)`. That is a 2-dof statistic and it keeps the rotation that
`|R_m|` throws away.
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
from engine import circular_symmetry as CS

OUT_JSON = HERE / "results_world_2dof_rescore.json"
ANGLES = W.ANGLES
ORDERS = W.ORDERS
STATS = tuple("%s_m%d" % (a, m) for a in ANGLES for m in ORDERS)


def main():
    rng = np.random.default_rng(W.RNG_SEED)
    cut = W.explore_cutoff()
    regions = {}
    print("2-DOF RE-SCORE (declustered, the primary arm). %d angular statistics/region"
          % len(STATS), flush=True)

    for path in sorted(glob.glob(str(HERE / "data" / "comcat_world" / "*.csv"))):
        name = os.path.basename(path)[:-4]
        t, la, lo, _n_hold, _n_full = W.load_region(path, declustered=True)
        n = t.size
        if n < 30:
            continue
        obs_f = W.features(t, la, lo)
        off = rng.uniform(-W.NULL_HALF_WINDOW_DAYS, W.NULL_HALF_WINDOW_DAYS,
                          size=(n, W.S_NULL_PER_EVENT))
        nf = W.features((t[:, None] + off).ravel(),
                        np.repeat(la, W.S_NULL_PER_EVENT),
                        np.repeat(lo, W.S_NULL_PER_EVENT))
        nf = {a: nf[a].reshape(n, W.S_NULL_PER_EVENT) for a in ANGLES}
        idx = rng.integers(0, W.S_NULL_PER_EVENT, size=(n, W.N_REPLICATES))

        per, cols = {}, {}
        for a in ANGLES:
            th_obs = obs_f[a]
            th_null = np.take_along_axis(nf[a], idx, axis=1)     # [n, R]
            for m in ORDERS:
                zo = np.exp(1j * m * th_obs).mean()
                zn = np.exp(1j * m * th_null).mean(axis=0)       # [R] complex
                mu = np.array([zn.real.mean(), zn.imag.mean()])
                cov = np.cov(np.vstack([zn.real, zn.imag]))
                inv = np.linalg.inv(cov)
                d = np.vstack([zn.real, zn.imag]).T - mu
                m_null = np.sqrt(np.einsum("ij,jk,ik->i", d, inv, d))
                do = np.array([zo.real, zo.imag]) - mu
                m_obs = float(np.sqrt(do @ inv @ do))
                k = "%s_m%d" % (a, m)
                # exact p of the observed Mahalanobis against its own null ensemble
                p = ((int(np.sum(m_null >= m_obs)) + 1) / (W.N_REPLICATES + 1))
                per[k] = {"mahalanobis_2dof": m_obs,
                          "null_p95": float(np.quantile(m_null, 0.95)),
                          "p_exact": p,
                          "R_observed": float(abs(zo)),
                          "R_null_mean": float(np.abs(zn).mean())}
                cols[k] = m_null
        top = min(per, key=lambda k: per[k]["p_exact"])
        print("  %-18s n=%5d  smallest p = %.4f (%s)"
              % (name, n, per[top]["p_exact"], top), flush=True)
        regions[name] = {"n": int(n), "per_statistic": per, "_cols": cols}

    # max-statistic across every region and angular statistic, one shared ensemble
    allc = np.stack([regions[r]["_cols"][k] for r in regions for k in STATS])
    gw = np.max(allc, axis=0)
    obs_max = max(regions[r]["per_statistic"][k]["mahalanobis_2dof"]
                  for r in regions for k in STATS)
    p_gw = (int(np.sum(gw >= obs_max)) + 1) / (W.N_REPLICATES + 1)
    worst = max(((r, k, regions[r]["per_statistic"][k]["mahalanobis_2dof"])
                 for r in regions for k in STATS), key=lambda x: x[2])
    for r in regions:
        regions[r].pop("_cols", None)

    out = {
        "arm": "2-dof re-score of the world scan's angular statistics",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 0,
        "why_priced_zero": ("same declared quantity, better estimator; no statistic "
                            "added. The motivation is a property of the dwell "
                            "distribution measurable without any events, so it cannot "
                            "be an estimator chosen after seeing the answer."),
        "supervisor_verification": {
            "kepler_claim": ("order-1 |R| statistics were BLIND and could not have "
                             "detected a total effect; strike them as evidence"),
            "supervisor_finding": ("CONFIRMED in structure, NOT in magnitude. |R| "
                                   "sensitivity is signed and crosses zero, and the "
                                   "2-dof form is 2x-3.4x more powerful. But measured "
                                   "power at a random unknown phase is 0.26-0.40 at "
                                   "order 1 and 0.63-0.76 at order 2, so |R| was "
                                   "WEAKENED, NOT BLIND. The scan's nulls stand at "
                                   "reduced power; they are not void."),
            "measured_power_random_phase": {
                "Chile_m1_eps0.30": {"R": 0.405, "twodof": 0.858},
                "Chile_m1_eps0.15": {"R": 0.105, "twodof": 0.360},
                "Alaska_m1_eps0.30": {"R": 0.255, "twodof": 0.718},
                "Chile_m2_eps0.30": {"R": 0.760, "twodof": 1.000},
            },
        },
        "statistic": ("Mahalanobis distance of the complex resultant Z_m from the null "
                      "mean under the null covariance of (Re Z, Im Z); 2 dof, keeps "
                      "the rotation that |R_m| discards"),
        "per_region": regions,
        "max_statistic": {
            "observed_max_2dof": float(obs_max),
            "null_max_p95": float(np.quantile(gw, 0.95)),
            "p": float(p_gw),
            "where": {"region": worst[0], "statistic": worst[1]},
        },
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print("\n" + "=" * 74)
    print("2-DOF MAX-STATISTIC: %.3f at %s / %s ; null 95th %.3f ; p = %.4f"
          % (obs_max, worst[0], worst[1], np.quantile(gw, 0.95), p_gw))
    print("  (the |R| version of the same declustered arm gave p = 0.9876)")
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
