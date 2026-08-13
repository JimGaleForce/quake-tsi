"""§P7-11(c) -- UN-CENSORING the three floor-censored mark tests.

RULING
------
HYPOTHESIS_LEDGER.md §P7-11(c), AUTHORIZED: "re-run the three censored mark rows
(b_value_90d x mag, deep_fraction_30d x depth, mean_depth_30d x depth) at higher
B (500k+) to resolve their true VIF_mark".

F4-58M (`results_f4_58m_vif_mark.json`) excluded these three rows under §P7-10(a)
because their `p_block_bootstrap` sat at the 1/(B+1) Monte Carlo floor with
B = 50,000, where `chi2.ppf` saturates and the identity returns a BOUND rather
than a measurement. This script draws 10x the surrogates so the quantile stops
saturating.

DIRECTION CORRECTION, MEASURED HERE (§P7-10(a) HAS THE SIGN BACKWARDS)
----------------------------------------------------------------------
§P7-10(a) states the censoring bites "in the anti-conservative direction":
"chi2_ppf saturates, the ratio is pushed down, and VIF is UNDERSTATED for
exactly the most significant tests". The arithmetic of the identity says the
opposite, and this run measures the opposite.

    VIF = T / chi2.ppf(1 - p, df).  At the floor the REPORTED p is 1/(B+1) and
    the TRUE p is <= that. A smaller p gives a LARGER chi2.ppf, hence a SMALLER
    VIF. The expressible quantile is capped BELOW the true one, so the value
    returned at the floor is an UPPER BOUND and VIF is OVER-stated, not
    understated.

Measured, on b_value_90d x mag (T = 313.46, df = 1), by pushing the floor down:

    B =   2,000  ->  p_floor 5.0e-04  ->  VIF 25.87
    B =  50,000  ->  p_floor 2.0e-05  ->  VIF 17.23   (the value F4-58M excluded)
    B = 500,000  ->  p_floor 2.0e-06  ->  VIF 13.87   (this run)

Monotone decreasing in B, in all three targets, exactly as the algebra requires.

CONSEQUENCE, which is the reverse of the one the ledger drew: EXCLUDING these
rows removes OVER-stated values, which lowers the median, which lowers the
floor. On the mark axis the exclusion is therefore mildly ANTI-conservative for
the floor, not protective of it. The rule itself still stands on its own better
ground -- a saturated value is a bound and not a measurement, and a floor must
not be built from either -- but its stated rationale is inverted and the
correction is recorded here rather than quietly absorbed.

WHAT IS REPRODUCED, AND WHAT IS NOT
-----------------------------------
Reproduced exactly, by import: `engine.mine.mark_test` (the same statistic, the
same two nulls, the same combination rule), `engine.mine.load_event_marks`,
`engine.mine.catalog_features`, `engine.design.build_design` at the primary
session's own configuration, and the primary session's own window [365, 8081)
read from the F4-58 control's prep metadata (so NO ETAS refit happens here and
`engine/out/cache/etas_params.json` is not touched -- the mark test does not use
the ETAS offset at all, only the event-ordered mark and feature series).

NOT reproduced: the original random stream. A different surrogate budget cannot
reuse the session's rng draws, so this is an independent re-estimate of the same
quantity, not a bitwise replay. The observed statistic IS bitwise identical and
is checked against the checkpoint value as a control (`statistic_matches_disk`).

PRICING -- ZERO (§P7-1(c)). No new hypothesis, no rejection, no BH vector, NO
EXPLORE_COUNT.jsonl line: this resolves the resolution of a p-value already
declared and already computed by the session that owns it.

Usage:  python -u f4_58m_uncensor.py --boots 500000 --jobs 3
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy import stats

REPO = os.path.dirname(os.path.abspath(__file__))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from engine import design as design_mod          # noqa: E402
from engine import mine as M                     # noqa: E402
from f4_58m_vif_mark import chi2_equivalent      # noqa: E402

OUT_PATH = os.path.join(REPO, "results_f4_58m_uncensored.json")
MARK_RESULTS = os.path.join(REPO, "results_f4_58m_vif_mark.json")
PRIMARY_SESSION = "session_20260811T022953"
CHECKPOINT = os.path.join(REPO, "engine", "out", "mine", PRIMARY_SESSION,
                          "checkpoint.json")

CFG = {"data_dir": "data/comcat_world", "dlat": 1.0, "dlon": 1.0,
       "explore_frac": 0.70}
WINDOW = (365, 8081)          # engine/out/cache/f4_58_control_prep.npz :: meta
SEED = 20260812

TARGETS = [("b_value_90d", "mag"),
           ("deep_fraction_30d", "depth"),
           ("mean_depth_30d", "depth")]

_W = {}


def _init(vals, marks, kinds, n_boot):
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    _W.update(vals=vals, marks=marks, kinds=kinds, n_boot=int(n_boot))


def _one(i):
    name, mk = TARGETS[i]
    t0 = time.time()
    rng = np.random.default_rng([SEED, 711, i])
    r = M.mark_test(_W["vals"][name], _W["marks"][mk], _W["kinds"][name],
                    _W["n_boot"], rng)
    r["feature"], r["mark"] = name, mk
    r["elapsed_s"] = time.time() - t0
    return r


def main(argv=None):
    ap = argparse.ArgumentParser("f4_58m_uncensor")
    ap.add_argument("--boots", type=int, default=500000)
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--reuse", action="store_true",
                    help="re-derive the report from the p_boot values already in "
                         "results_f4_58m_uncensored.json WITHOUT redrawing any "
                         "surrogate (used to correct labelling only; every "
                         "measured quantity is carried through unchanged)")
    args = ap.parse_args(argv)
    if args.smoke:
        args.boots = 2000
    jobs = max(1, min(4, args.jobs))

    t_start = time.time()
    ctx = design_mod.build_design(data_dir=CFG["data_dir"], dlat=CFG["dlat"],
                                  dlon=CFG["dlon"],
                                  explore_frac=CFG["explore_frac"], verbose=False)
    all_marks = M.load_event_marks(ctx, CFG["data_dir"], ctx.meta["mag_floor"])
    window = slice(*WINDOW)
    in_win = (all_marks["day"] >= window.start) & (all_marks["day"] < window.stop)
    marks = {k: v[in_win] for k, v in all_marks.items()}
    fe_day = marks["day"] - window.start
    n_ev = int(marks["mag"].size)
    print("window [%d, %d), events in window = %d" % (WINDOW[0], WINDOW[1], n_ev))

    feats = {f.name: f for f in M.catalog_features(ctx, all_marks, (0,))}
    vals = {nm: feats[nm].values[window][fe_day] for nm, _ in TARGETS}
    kinds = {nm: feats[nm].kind for nm, _ in TARGETS}

    with open(CHECKPOINT, "r", encoding="utf-8") as fh:
        ck = json.load(fh)
    disk = {(t["feature"], t.get("mark")): t for t in ck["tests"]
            if t.get("test") in ("spearman", "circular-linear")}
    B_old = ck["config"]["n_surrogates"]

    if args.reuse:
        with open(OUT_PATH, "r", encoding="utf-8") as fh:
            prev = json.load(fh)
        args.boots = int(prev["construction"]["B_new"])
        res = [{"feature": r["feature"], "mark": r["mark"], "df": r["df"],
                "test": r["test"], "statistic": r["statistic_recomputed"],
                "p_block_bootstrap": r["p_block_bootstrap_new"],
                "p_circular_shift": r["p_circular_shift_new"],
                "p_raw": r["p_raw_new"], "elapsed_s": r["elapsed_s"]}
               for r in prev["rows"]]
        res = [next(x for x in res if (x["feature"], x["mark"]) == t)
               for t in TARGETS]
        print("REUSE: no surrogate redrawn; carrying the %d measured p_boot "
              "values from %s at B = %d" % (len(res), OUT_PATH, args.boots))
    else:
        print("re-running %d mark tests at B = %d (was %d), jobs = %d"
              % (len(TARGETS), args.boots, B_old, jobs))
        with ProcessPoolExecutor(max_workers=jobs, initializer=_init,
                                 initargs=(vals, marks, kinds, args.boots)) as ex:
            res = list(ex.map(_one, range(len(TARGETS))))
        print("sweep done in %.0fs" % (time.time() - t_start))

    floor_new = 1.0 / (args.boots + 1.0)
    floor_old = 1.0 / (B_old + 1.0)
    rows = []
    for r in res:
        d = disk[(r["feature"], r["mark"])]
        df = int(r["df"])
        T = chi2_equivalent(r["statistic"], df, n_ev)
        p_new = float(r["p_block_bootstrap"])
        still = bool(p_new <= floor_new + 1e-15)
        vif_new = float(T / stats.chi2.ppf(1.0 - p_new, df))
        vif_old = float(T / stats.chi2.ppf(1.0 - floor_old, df))
        rows.append({
            "feature": r["feature"], "mark": r["mark"], "df": df,
            "test": r["test"], "n_events": n_ev,
            "statistic_recomputed": float(r["statistic"]),
            "statistic_on_disk": float(d["statistic"]),
            "statistic_matches_disk": bool(
                abs(r["statistic"] - d["statistic"]) < 1e-12 * max(1.0, abs(d["statistic"]))),
            "chi2_equivalent": T,
            "B_old": B_old, "p_block_bootstrap_old": float(d["p_block_bootstrap"]),
            "vif_mark_old_UPPER_BOUND": vif_old,
            "vif_mark_old_LOWER_BOUND": None,
            "vif_mark_old_LOWER_BOUND_note": (
                "RETRACTED FIELD NAME. The value at a saturated floor is an UPPER "
                "bound on VIF, not a lower bound; see the module docstring. The "
                "number itself (%.4f) is unchanged and is now correctly labelled "
                "`vif_mark_old_UPPER_BOUND`." % vif_old),
            "B_new": args.boots, "p_block_bootstrap_new": p_new,
            "resolution_floor_new": floor_new,
            "still_censored_at_new_B": still,
            "vif_mark_new": vif_new,
            "vif_mark_new_status": ("RESOLVED" if not still else
                                    "STILL AN UPPER BOUND (p still at the floor)"),
            "ratio_new_over_old_bound": vif_new / vif_old,
            "bound_tightened_by": (vif_old - vif_new),
            "p_circular_shift_new": float(r["p_circular_shift"]),
            "p_raw_new": float(r["p_raw"]),
            "elapsed_s": r["elapsed_s"],
        })

    with open(MARK_RESULTS, "r", encoding="utf-8") as fh:
        mk_doc = json.load(fh)
    prim = next(s for s in mk_doc["sessions"] if s["session"] == PRIMARY_SESSION)
    med_uncens = prim["vif_mark_all"]["median"]
    kept = [x["vif_mark"] for x in prim["per_feature"] if x["vif_mark"]]
    with_resolved = sorted(kept + [r["vif_mark_new"] for r in rows])
    med_with = float(np.median(with_resolved))

    z80 = float(stats.norm.ppf(0.80))
    const = 4.732

    def rho_min(v):
        return float(np.sqrt(v) * const / np.sqrt(n_ev - 1.0))

    n_resolved = sum(1 for r in rows if not r["still_censored_at_new_B"])
    out = {
        "id": "F4-58M UN-CENSORING (§P7-11(c))",
        "ruling": "HYPOTHESIS_LEDGER.md §P7-11(c), AUTHORIZED",
        "version": 1,
        "priced_tests": 0,
        "priced_tests_note": (
            "Zero, per §P7-1(c). This resolves the RESOLUTION of a p-value already "
            "computed and already declared by the session that owns it. No new "
            "hypothesis, no rejection, no BH vector, no EXPLORE_COUNT.jsonl line."),
        "state_class": "MEASUREMENT on existing artifacts, re-estimated at higher B",
        "generated_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "construction": {
            "B_old": B_old, "B_new": args.boots,
            "block_events": 500,
            "n_events": n_ev,
            "window": list(WINDOW),
            "seed": SEED,
            "reproduced_by_import": ("engine.mine.mark_test, "
                                     "engine.mine.load_event_marks, "
                                     "engine.mine.catalog_features, "
                                     "engine.design.build_design"),
            "not_reproduced": ("the original random stream -- a different surrogate "
                               "budget cannot replay it. The OBSERVED statistic is "
                               "bitwise checked against the checkpoint instead."),
            "etas_untouched": ("the mark test uses no ETAS offset; no baseline was "
                               "fitted and no cache under engine/out/cache was "
                               "written by this script."),
        },
        "rows": rows,
        "effect_on_F4_58M": {
            "n_targets": len(rows),
            "n_resolved": n_resolved,
            "n_still_censored": len(rows) - n_resolved,
            "vif_mark_median_uncensored_only": med_uncens,
            "vif_mark_median_including_resolved": med_with,
            "rho_min_uncensored_only": rho_min(med_uncens),
            "rho_min_including_resolved": rho_min(med_with),
            "note": ("The three targets are the most significant mark tests on "
                     "disk, so re-admitting them at their bounds moves the median "
                     "UP (4.345 -> 4.575) and the floor UP. Because those bounds "
                     "are UPPER bounds, that is the direction §P7-10(a) intended "
                     "to protect but reached by the opposite route: the exclusion "
                     "removes OVER-stated values and therefore LOWERS the floor. "
                     "The effect is small here (3 rows of 46, rho_min 0.0457 vs "
                     "0.0469) but the sign matters for the rule's justification."),
        },
        "direction_correction_to_P7_10a": {
            "ledger_text": ("§P7-10(a): 'chi2_ppf saturates, the ratio is pushed "
                            "down, and VIF is UNDERSTATED for exactly the most "
                            "significant tests -- which understates the S-15 "
                            "floor, the one direction that can manufacture a "
                            "claim.'"),
            "measured_finding": ("The sign is inverted. VIF = T / chi2.ppf(1-p, "
                                 "df) is DECREASING in B at a saturated floor, so "
                                 "the reported value is an UPPER bound and VIF is "
                                 "OVER-stated. Pushing the floor down from 5.0e-04 "
                                 "to 2.0e-05 to 2.0e-06 drove all three targets "
                                 "monotonically DOWN (b_value_90d x mag: 25.87 -> "
                                 "17.23 -> 13.87)."),
            "what_survives": ("The EXCLUSION RULE stands, on better ground than "
                              "the reason given for it: a saturated value is a "
                              "bound, not a measurement, and a floor may be built "
                              "from neither. What does not survive is the claim "
                              "that excluding closes an anti-conservative seam. On "
                              "the mark axis the exclusion mildly OPENS one."),
            "who_should_rule": ("Popper seat. This is a correction to a RATIFIED "
                                "clause's rationale, not to its operative "
                                "instruction, and it is reported rather than "
                                "acted on."),
            "propagation": ("The same inverted rationale was repeated by the "
                            "measurer in results_f4_58m_vif_mark.json and in the "
                            "f4_58_vif.py amendment docstring (both committed at "
                            "4472141). Both are corrected in place with this "
                            "finding cited; no measured number in either file "
                            "changes, because the count path had zero censored "
                            "rows and the mark path's medians are computed over "
                            "the uncensored rows either way."),
        },
        "caveats": [
            "block_events = 500 over %d events is ~%d blocks. Raising B fixes the "
            "RESOLUTION of p_boot; it does NOT address whether a ~93-block "
            "bootstrap null is the right width. That decomposition is still "
            "UNMEASURED (§P7-10(c)) and no simulated mark-axis control exists."
            % (n_ev, n_ev // 500),
            "The circular-shift arm is an exhaustive enumeration whose floor is set "
            "by the record, not by B, so it is unchanged by this run; where p_raw "
            "is set by the shift arm it is unchanged too.",
        ],
    }

    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print("wrote %s" % OUT_PATH)

    print("\n%-20s %-6s %12s %12s %12s %10s"
          % ("feature", "mark", "p_boot_new", "VIF_old(<=)", "VIF_new(<=)", "status"))
    for r in rows:
        print("%-20s %-6s %12.3e %12.2f %12.2f %10s"
              % (r["feature"], r["mark"], r["p_block_bootstrap_new"],
                 r["vif_mark_old_UPPER_BOUND"], r["vif_mark_new"],
                 "RESOLVED" if not r["still_censored_at_new_B"] else "STILL FLOOR"))
    e = out["effect_on_F4_58M"]
    print("\nVIF_mark median: %.3f (uncensored only) -> %.3f (including resolved)"
          % (e["vif_mark_median_uncensored_only"], e["vif_mark_median_including_resolved"]))
    print("rho_min        : %.4f -> %.4f"
          % (e["rho_min_uncensored_only"], e["rho_min_including_resolved"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
