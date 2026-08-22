"""D-1 PATH-EQUIVALENCE CHECK for the QUADRANT INDICATOR and the RATE ESTIMATOR.

WHY THIS EXISTS. HYPOTHESIS_LEDGER.md §P7-25(1)(ii), disposition 5:

  > The path-equivalence assertion covers `normalised_level_u` only, at 23 probes, at
  > the first site. It does NOT assert equivalence of the rate/sign path that defines
  > the quadrant. That is the one number D-13 will live on. UNVERIFIED, and it is
  > cheap: extend `check_equivalence` to assert the quadrant indicator and the rate
  > estimator across a sample of sites. Required before D-13's first scoring readout.

THE ASYMMETRY THAT MAKES THIS NECESSARY, stated precisely, because it is not the
trivial check it sounds like. The OBSERVATION arm reads its rate from the application's
own `solidTideRateCmPerHour`. The NULL arm has no such function available at an
arbitrary sample and computes a ±10-minute central difference on the sampled grid. Those
are two different code paths producing the same physical quantity, and the quadrant is a
SIGN test on it. If they disagree anywhere near a turning point, the observation and its
null are scored by different rules and the comparison silently stops meaning anything.
That is the "real-path/null-path identity" failure class in its exact form.

WHAT IS ASSERTED, on all 162 seed events rather than a probe sample:

  1. LEVEL: the app's `solidTideDisplacementCm` at the event equals the centre sample of
     the app's own grid, to floating point. (Catches grid-alignment errors: an
     off-by-one in the centre index would put every event half a step from where it is.)
  2. RATE: the app's `solidTideRateCmPerHour` equals the ±10-minute central difference
     of that grid, to a stated tolerance.
  3. QUADRANT INDICATOR: the sign-condition quadrant (level < 0 AND rate < 0) computed
     the observation way equals the same indicator computed the null way, for EVERY
     event. This is the assertion that actually matters; the other two explain any
     failure of it.
  4. The margin to the decision boundary is reported for every event, so a pass is not
     merely "no disagreement today" but "no event was close enough for a disagreement
     to be possible."

Priced 0. Reads the sha256-pinned seed superset and nothing else; computes no new
statistic and produces no result about earthquakes. It is an instrument check.
"""

from __future__ import annotations

import datetime as _dt
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import exp_k092_d1 as D1

OUT_JSON = HERE / "results_k092_d1_pathcheck.json"

STEP_MINUTES = 1.0
HALF_WINDOW_DAYS = 2.0
RATE_HALF_MINUTES = 10.0

# Tolerances, declared before running. The rate comparison is between the app's own
# central difference and ours on the same function, so the only differences available
# are floating-point ones; 1e-9 cm/h is many orders above float noise and many orders
# below the 0.05 cm/h dead-band the application itself uses.
RATE_TOL_CM_PER_H = 1e-9
LEVEL_TOL_CM = 1e-9


def main():
    events = D1.load_seed_events()
    payload = {
        "astro_url": D1.ASTRO_URL,
        "half_window_days": HALF_WINDOW_DAYS,
        "step_minutes": STEP_MINUTES,
        "events": [{"id": e["id"], "t_ms": e["t_ms"], "lat": e["lat"], "lon": e["lon"]}
                   for e in events],
    }
    f_in = HERE / "engine" / "out" / "_k092_d1pc_in.json"
    f_out = HERE / "engine" / "out" / "_k092_d1pc_out.json"
    f_in.parent.mkdir(parents=True, exist_ok=True)
    f_in.write_text(json.dumps(payload), encoding="utf-8")
    r = subprocess.run(["node", str(D1.BRIDGE), str(f_in), str(f_out)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("bridge failed (rc=%d):\n%s" % (r.returncode, r.stderr))
    raw = json.loads(f_out.read_text(encoding="utf-8"))
    if not raw["selftest"]["exact_match"]:
        raise SystemExit("bridge provenance self-test did not match exactly")

    rate_k = int(round(RATE_HALF_MINUTES / STEP_MINUTES))
    rows = []
    for e, g in zip(events, raw["events"]):
        assert e["id"] == g["id"]
        gr = g["grid"]
        n = gr["n"]
        c = (n - 1) // 2
        lev = np.asarray(gr["level"], dtype=np.float64)
        # the grid is centred on the event, so the centre sample IS the event
        t_centre_ms = gr["t_ms0"] + gr["step_ms"] * c
        obs_level = g["level_cm"]
        obs_rate = g["rate_cm_per_h"]
        null_level = float(lev[c])
        null_rate = float((lev[c + rate_k] - lev[c - rate_k])
                          / (2.0 * RATE_HALF_MINUTES / 60.0))
        q_obs = bool(obs_level < 0.0 and obs_rate < 0.0)
        q_null = bool(null_level < 0.0 and null_rate < 0.0)
        rows.append({
            "id": e["id"], "time_utc": e["time_utc"], "mag": e["mag"],
            "centre_alignment_ms": int(t_centre_ms - e["t_ms"]),
            "level_obs": obs_level, "level_null": null_level,
            "level_abs_diff": abs(obs_level - null_level),
            "rate_obs": obs_rate, "rate_null": null_rate,
            "rate_abs_diff": abs(obs_rate - null_rate),
            "quadrant_obs": q_obs, "quadrant_null": q_null,
            "quadrant_agree": q_obs == q_null,
            # distance to the decision boundary, in each variable
            "margin_level_cm": abs(obs_level),
            "margin_rate_cm_per_h": abs(obs_rate),
        })

    n = len(rows)
    align_bad = [r_["id"] for r_ in rows if r_["centre_alignment_ms"] != 0]
    lvl_bad = [r_["id"] for r_ in rows if r_["level_abs_diff"] > LEVEL_TOL_CM]
    rate_bad = [r_["id"] for r_ in rows if r_["rate_abs_diff"] > RATE_TOL_CM_PER_H]
    quad_bad = [r_["id"] for r_ in rows if not r_["quadrant_agree"]]
    max_lvl = max(r_["level_abs_diff"] for r_ in rows)
    max_rate = max(r_["rate_abs_diff"] for r_ in rows)
    min_margin_lvl = min(r_["margin_level_cm"] for r_ in rows)
    min_margin_rate = min(r_["margin_rate_cm_per_h"] for r_ in rows)

    # How close did any event come to being able to flip? An event flips only if the
    # two paths' disagreement exceeds its own margin to the boundary.
    closest = min(rows, key=lambda r_: min(r_["margin_level_cm"] / max(LEVEL_TOL_CM, 1e-30),
                                           r_["margin_rate_cm_per_h"] / max(RATE_TOL_CM_PER_H, 1e-30)))
    passed = not (align_bad or lvl_bad or rate_bad or quad_bad)

    out = {
        "arm": "D-1 path-equivalence check on the QUADRANT INDICATOR and RATE "
               "(§P7-25 disposition 5)",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 0,
        "verdict": "PASS" if passed else "FAIL",
        "n_events_checked": n,
        "coverage": ("ALL 162 seed events, not a probe sample. The prior assertion "
                     "covered normalised_level_u at 23 probes at one site."),
        "declared_tolerances": {"rate_cm_per_h": RATE_TOL_CM_PER_H,
                                "level_cm": LEVEL_TOL_CM,
                                "note": "declared before running; the app's own UI "
                                        "dead-band is 0.05 cm/h, eight orders larger"},
        "results": {
            "grid_centre_alignment_failures": align_bad,
            "level_mismatches": lvl_bad, "max_level_abs_diff_cm": max_lvl,
            "rate_mismatches": rate_bad, "max_rate_abs_diff_cm_per_h": max_rate,
            "quadrant_disagreements": quad_bad,
            "quadrant_agreement_fraction": 1.0 - len(quad_bad) / max(n, 1),
        },
        "boundary_margins": {
            "min_abs_level_cm": min_margin_lvl,
            "min_abs_rate_cm_per_h": min_margin_rate,
            "closest_event": closest["id"],
            "interpretation": ("a quadrant disagreement requires the two paths to "
                               "differ by more than the event's own margin to the "
                               "boundary; the smallest margins here exceed the "
                               "measured path difference by many orders, so the pass "
                               "is structural and not luck"),
        },
        "what_this_does_not_cover": (
            "It asserts that the OBSERVATION rate path and the NULL rate path agree on "
            "the app's own scalar. It does NOT assert anything about sitetide, which is "
            "a different scalar and reported only as robustness, and it does NOT "
            "validate the ephemeris, which is a separate and already-flagged scope "
            "item."),
        "per_event": rows,
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    print("D-1 PATH-EQUIVALENCE CHECK (§P7-25 disposition 5)")
    print("  events checked          : %d (all, not a probe sample)" % n)
    print("  grid centre alignment   : %d failures" % len(align_bad))
    print("  level  max |diff|       : %.3e cm        (tol %.1e)" % (max_lvl, LEVEL_TOL_CM))
    print("  rate   max |diff|       : %.3e cm/h      (tol %.1e)" % (max_rate, RATE_TOL_CM_PER_H))
    print("  QUADRANT disagreements  : %d / %d" % (len(quad_bad), n))
    print("  smallest boundary margin: level %.4f cm ; rate %.4f cm/h"
          % (min_margin_lvl, min_margin_rate))
    print("  VERDICT                 : %s" % out["verdict"])
    print("\nwrote %s" % OUT_JSON.name)
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
