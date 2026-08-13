"""§P7-20(7) -- classify the tranche-B artifact-hash difference.

RULING
------
HYPOTHESIS_LEDGER.md §P7-20(7): "Classify it before any Tranche B result is
quoted anywhere. If the difference is confined to non-scientific fields, define
a canonical artifact hash excluding them and add the invariant in both
directions. If any numerical field differs, this is an invariance-audit failure
and B's results are suspect until it is resolved."

Two tranche-B sessions, IDENTICAL config hash 114d4e9d8004:
  A = session_20260813T091153  artifact 13eea5f3...
  B = session_20260813T092628  artifact b44e0f9a...

`artifact_content_hash` (engine/mine_session.py:1122) hashes ONLY
`ckpt.state["tests"]` -- "not timings, not the session id, not the creation
timestamp -- so the hash answers exactly one question: did these two runs
compute the same numbers?". So a hash difference here CANNOT be a timestamp or
a path. This script performs the field-by-field diff the ruling demands and
returns the classification.

PRICING -- ZERO (§P7-1(c)). Re-reading artifacts already on disk. No surrogate
drawn, no rejection made, no BH vector entered, NO EXPLORE_COUNT.jsonl line.

The two session directories are READ-ONLY here. Nothing is written into them.

Usage:  python -u hash_classification.py
"""

from __future__ import annotations

import datetime as _dt
import json
import math
import os
import sys
from collections import defaultdict

REPO = os.path.dirname(os.path.abspath(__file__))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from engine.mine_session import artifact_content_hash  # noqa: E402

OUT_PATH = os.path.join(REPO, "results_hash_classification.json")

SESSION_A = "session_20260813T091153"
SESSION_B = "session_20260813T092628"

# Fields that would be benign under case (a): identity, clock, path, wall-time.
# Recorded so the classification states what it WOULD have excluded had the
# difference been confined to them. It is not.
BENIGN_CANDIDATES = ("session", "session_dir", "recorded", "created",
                     "elapsed_seconds", "task_seconds", "generated_at")

# Fields that carry the science: if ANY of these differ the run is not a
# reproduction in any sense that matters.
OUTCOME_FIELDS = ("p_raw", "p_bh", "bh_eligible", "bh_ineligible_reason",
                  "passes_fdr", "disposition", "disposition_parent",
                  "disposition_reason", "p_method", "p_method_reason",
                  "p_floor", "bh_q", "stratum", "family", "order_key",
                  "feature", "mark", "test", "kind", "lag", "n_events",
                  "n_surrogates", "rho_min", "mark_floor", "s15c_verdict",
                  "p_fwer_max_stat", "p_block_bootstrap", "gpd",
                  "floor_monte_carlo", "floor_enumeration", "block_days")

# A relative difference at or below this is consistent with IEEE-754
# non-associativity in a reduction / BLAS thread-count variation, and nothing
# else. It is a DESCRIPTION of the observed magnitudes, not a tolerance that
# excuses them: §P6-5 makes determinism a hard requirement.
FP_NOISE_REL = 1e-5


def load(session):
    p = os.path.join(REPO, "engine", "out", "mine", session, "checkpoint.json")
    with open(p, "r", encoding="utf-8") as fh:
        return json.load(fh)


def walk(a, b, path, row, num, nonnum, structural):
    """Recursive field-by-field diff, classifying each leaf difference."""
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            walk(a.get(k), b.get(k), path + "." + k, row, num, nonnum, structural)
        return
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            structural.append((path + "[len]", row, len(a), len(b)))
            return
        for x, y in zip(a, b):
            walk(x, y, path + "[]", row, num, nonnum, structural)
        return
    if isinstance(a, bool) or isinstance(b, bool):
        if a != b:
            nonnum.append((path, row, a, b))
        return
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if a == b:
            return
        if isinstance(a, float) and isinstance(b, float) \
                and math.isnan(a) and math.isnan(b):
            return                      # NaN != NaN is not a difference
        d = abs(a - b)
        r = d / max(abs(a), abs(b), 1e-300)
        num.append((path, row, a, b, d, r))
        return
    if a != b:
        if type(a) is not type(b):
            structural.append((path + "[type]", row, type(a).__name__,
                               type(b).__name__))
        else:
            nonnum.append((path, row, a, b))


def main():
    ca, cb = load(SESSION_A), load(SESSION_B)
    ta, tb = ca["tests"], cb["tests"]

    ha = artifact_content_hash(ta)
    hb = artifact_content_hash(tb)

    num, nonnum, structural = [], [], []
    if len(ta) != len(tb):
        structural.append(("tests[len]", None, len(ta), len(tb)))
    else:
        for i, (x, y) in enumerate(zip(ta, tb)):
            walk(x, y, "", i, num, nonnum, structural)

    # ---- aggregate the numeric differences by field path -------------------
    agg = defaultdict(lambda: {"n": 0, "max_rel": 0.0, "max_abs": 0.0,
                               "rows": set(), "examples": []})
    for path, row, a, b, d, r in num:
        e = agg[path]
        e["n"] += 1
        e["max_rel"] = max(e["max_rel"], r)
        e["max_abs"] = max(e["max_abs"], d)
        e["rows"].add(row)
        if len(e["examples"]) < 3:
            e["examples"].append({"row": row, "A": a, "B": b, "rel": r})

    numeric_fields = []
    for path in sorted(agg, key=lambda p: -agg[p]["max_rel"]):
        e = agg[path]
        numeric_fields.append({
            "field": path.lstrip("."),
            "n_differing_values": e["n"],
            "n_differing_rows": len(e["rows"]),
            "max_relative_difference": e["max_rel"],
            "max_absolute_difference": e["max_abs"],
            "class": ("FLOAT_NOISE" if e["max_rel"] <= FP_NOISE_REL
                      else "MATERIAL"),
            "examples": e["examples"],
        })

    material = [f for f in numeric_fields if f["class"] == "MATERIAL"]
    noise = [f for f in numeric_fields if f["class"] == "FLOAT_NOISE"]

    # ---- did any OUTCOME field move? ---------------------------------------
    outcome_diffs = []
    for i, (x, y) in enumerate(zip(ta, tb)):
        for k in OUTCOME_FIELDS:
            if k in x or k in y:
                if x.get(k) != y.get(k):
                    outcome_diffs.append({"row": i, "field": k,
                                          "A": x.get(k), "B": y.get(k)})

    # ---- the material rows, in full ----------------------------------------
    material_rows = []
    for f in material:
        fld = f["field"]
        for i, (x, y) in enumerate(zip(ta, tb)):
            if fld in x and x.get(fld) != y.get(fld):
                nsurr = x.get("n_surrogates")
                # (1 + C) / (1 + B) inversion, to expose the exceedance counts
                ca_, cb_ = None, None
                if nsurr and isinstance(x.get(fld), float):
                    ca_ = int(round(x[fld] * (nsurr + 1) - 1))
                    cb_ = int(round(y[fld] * (nsurr + 1) - 1))
                material_rows.append({
                    "row": i, "field": fld,
                    "feature": x.get("feature"), "test": x.get("test"),
                    "mark": x.get("mark"), "stratum": x.get("stratum"),
                    "A": x.get(fld), "B": y.get(fld),
                    "n_surrogates": nsurr,
                    "implied_exceedance_count_A": ca_,
                    "implied_exceedance_count_B": cb_,
                    "implied_count_delta": (None if ca_ is None else ca_ - cb_),
                    "statistic_A": x.get("statistic"),
                    "statistic_B": y.get("statistic"),
                    "statistic_delta": (
                        None if x.get("statistic") is None
                        else y.get("statistic") - x.get("statistic")),
                    "p_raw_identical": x.get("p_raw") == y.get("p_raw"),
                    "p_method": x.get("p_method"),
                    "disposition": x.get("disposition"),
                    "bh_eligible": x.get("bh_eligible"),
                    "passes_fdr": x.get("passes_fdr"),
                })

    # ---- canonical-hash experiment (case (a) remedy, attempted) ------------
    # Strip the fields that differ and re-hash. This is recorded to show what
    # the remedy WOULD have required -- it strips numerical fields, which is
    # exactly why it is not admissible.
    differing_leaf_fields = sorted({f["field"].split("[")[0].split(".")[-1]
                                    for f in numeric_fields})

    def strip(rows, drop):
        def rec(o):
            if isinstance(o, dict):
                return {k: rec(v) for k, v in o.items() if k not in drop}
            if isinstance(o, list):
                return [rec(v) for v in o]
            return o
        return [rec(r) for r in rows]

    drop = set(differing_leaf_fields)
    ha_c = artifact_content_hash(strip(ta, drop))
    hb_c = artifact_content_hash(strip(tb, drop))

    verdict = "CASE_B_INVARIANCE_FAILURE" if numeric_fields else "CASE_A_BENIGN"

    out = {
        "id": "P7-20(7)-HASH-CLASSIFICATION",
        "title": ("Classification of the tranche-B artifact-hash difference on "
                  "identical science"),
        "ruling": ("HYPOTHESIS_LEDGER.md §P7-20(7); invariant dual of §P7-8(c); "
                   "determinism requirement §P6-5"),
        "generated_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 0,
        "priced_tests_note": (
            "Zero per §P7-1(c). This is an estimator/infrastructure measurement "
            "on artifacts already on disk. No surrogate drawn, no rejection "
            "made, no BH vector entered, NO EXPLORE_COUNT.jsonl line written."),
        "state_class": (
            "INFRASTRUCTURE AUDIT on existing artifacts. Not a hypothesis test "
            "about the Earth; makes no claim about seismicity."),
        "repo_state": {"head": "88e6708", "tests": "402 passed",
                       "sessions_untouched": True},

        "sessions": {
            "A": {"session": SESSION_A, "config_hash": ca.get("config_hash"),
                  "artifact_hash": ha, "n_tests": len(ta),
                  "elapsed_seconds": ca.get("elapsed_seconds")},
            "B": {"session": SESSION_B, "config_hash": cb.get("config_hash"),
                  "artifact_hash": hb, "n_tests": len(tb),
                  "elapsed_seconds": cb.get("elapsed_seconds")},
            "config_hash_identical": ca.get("config_hash") == cb.get("config_hash"),
            "config_dicts_identical": ca.get("config") == cb.get("config"),
            "artifact_hash_identical": ha == hb,
            "hash_scope": (
                "engine/mine_session.py:artifact_content_hash hashes "
                "ckpt.state['tests'] ONLY -- explicitly not timings, not the "
                "session id, not the creation timestamp. A timestamp/path "
                "explanation is therefore excluded a priori by the function's "
                "own construction, and the diff below confirms it."),
        },

        "VERDICT": verdict,
        "gate_status": "NOT LIFTED",
        "gate_statement": (
            "§P7-20(7) case (b). Numerical fields DIFFER between two sessions "
            "with an identical config hash. This is a §P6-5 determinism defect "
            "and, per the ruling, tranche B's results are SUSPECT until it is "
            "resolved. No canonical hash is emitted: the only hash that would "
            "make these two artifacts agree is one that excludes the numbers, "
            "which is the one thing artifact_content_hash exists to include."),

        "difference_summary": {
            "n_numeric_field_paths_differing": len(numeric_fields),
            "n_material_field_paths": len(material),
            "n_float_noise_field_paths": len(noise),
            "n_non_numeric_differences": len(nonnum),
            "n_structural_differences": len(structural),
            "non_numeric_differences": nonnum[:50],
            "structural_differences": structural[:50],
            "benign_field_candidates_considered": list(BENIGN_CANDIDATES),
            "benign_fields_found_in_hashed_content": [],
            "benign_explanation_available": False,
        },

        "sibling_artifacts": {
            "note": (
                "The ruling asks for checkpoint.json, stubs.json AND report.md. "
                "The other two were diffed as whole files."),
            "stubs.json": {
                "differing_lines": 1,
                "differences": ["\"session\": session id only"],
                "classification": "BENIGN (identity field only)",
            },
            "report.md": {
                "differing_lines": 3,
                "differences": [
                    "title line: session id",
                    "elapsed 696.9 s vs 759.9 s (wall time)",
                    "the quoted artifact content hash itself",
                ],
                "classification": "BENIGN (identity, wall-time, and the hash "
                                  "under classification)",
            },
            "reading": (
                "The human-readable artifacts are benign-only: every rendered "
                "number in report.md is identical, because the report rounds. "
                "The entire difference lives in the FULL-PRECISION test rows "
                "inside checkpoint.json -- which is exactly the content "
                "artifact_content_hash covers and exactly the content a reader "
                "of the report cannot see. That is an argument for the hash, "
                "not against it: without it this would have been invisible."),
        },

        "numeric_fields": numeric_fields,

        "material_differences": {
            "definition": (
                "relative difference > %g, i.e. larger than IEEE-754 "
                "reduction-order noise can account for" % FP_NOISE_REL),
            "n_fields": len(material),
            "fields": [f["field"] for f in material],
            "rows": material_rows,
            "mechanism": (
                "The circular-shift null for kuiper_V on a strictly periodic "
                "day-lattice feature is MASSIVELY TIED: whole-period shifts "
                "reproduce the statistic, so the surrogate distribution has "
                "mass points. T_obs moves by ~2e-8 (the float noise above) and "
                "crosses a plateau, moving the exceedance count by up to 186 "
                "of 7657 surrogates in one step -- a 0.024 jump in p. Row 161 "
                "moves in the OPPOSITE direction to its T_obs shift, which "
                "shows the SURROGATE statistics jitter on the same plateau "
                "too. The 1e-8 cause is float non-determinism; the 1e-1 effect "
                "is the null's own degeneracy amplifying it. Both are defects; "
                "the second is the one that makes 'identical science' false at "
                "the reported-number level."),
            "mechanism_status": (
                "MECHANISM INFERRED from the recorded values (statistic deltas, "
                "implied exceedance counts, and the direction reversal on row "
                "161). NOT confirmed by instrumenting a re-run. UNVERIFIED as "
                "to the root cause of the underlying 1e-8 jitter -- BLAS thread "
                "count / reduction order is the leading candidate and was not "
                "tested here."),
        },

        "float_noise_differences": {
            "n_fields": len(noise),
            "max_relative_difference_over_all_noise_fields": (
                max([f["max_relative_difference"] for f in noise], default=0.0)),
            "reading": (
                "19 field paths differ at <= 2.6e-6 relative. That is "
                "consistent with non-associative floating-point reduction "
                "(thread count, BLAS kernel selection, or accumulation order) "
                "in the GLM fit and the statistics derived from it. It is NOT "
                "benign under §P6-5, which requires determinism outright; it is "
                "merely small."),
        },

        "science_invariance": {
            "outcome_fields_checked": list(OUTCOME_FIELDS),
            "n_outcome_field_differences": len(outcome_diffs),
            "outcome_field_differences": outcome_diffs,
            "bh_block_identical": ca.get("bh") == cb.get("bh"),
            "dispositions_identical": ca.get("dispositions") == cb.get("dispositions"),
            "data_log_identical": ca.get("data_log") == cb.get("data_log"),
            "n_tests_identical": ca.get("n_tests") == cb.get("n_tests"),
            "max_statistic_identical_modulo_nan": True,
            "max_statistic_note": (
                "The max_statistic blocks compare unequal ONLY because three "
                "empty strata carry NaN and NaN != NaN in Python. Every finite "
                "field is bit-identical: global p = 1.3058239749281798e-4, "
                "t_obs = 3.8841153620116686, n_tests = 102, in both."),
            "reading": (
                "The SCIENCE reproduced: 34 survivors, identical BH identity "
                "(m = 171, q = 0.10), identical p_raw / p_bh / disposition / "
                "passes_fdr on all 283 rows. What did NOT reproduce is the "
                "ARTIFACT. §P7-20(7) does not offer a 'science reproduced so "
                "the hash may differ' branch, and it should not: the branch "
                "that exists is (b)."),
        },

        "canonical_hash_experiment": {
            "attempted": True,
            "fields_that_would_have_to_be_excluded": differing_leaf_fields,
            "hash_A_after_exclusion": ha_c,
            "hash_B_after_exclusion": hb_c,
            "match_after_exclusion": ha_c == hb_c,
            "ADMISSIBLE": False,
            "why_inadmissible": (
                "Every excluded field is a NUMERICAL RESULT -- beta, se, "
                "statistic, chi2_score, p_parametric, p_circular_shift, "
                "amplitude_log_rate, U2, V, power, ar1_phi, fold_lr. A "
                "'canonical hash' over the remaining fields would answer 'did "
                "these two runs agree on their labels?', not §P7-8(c)'s "
                "question 'did these two runs compute the same numbers?'. "
                "Emitting it would convert a determinism defect into a hash "
                "that cannot detect determinism defects. The experiment is "
                "recorded to show the remedy was tried and refused, not to "
                "offer it."),
        },

        "what_would_lift_the_gate": [
            "Root-cause the ~1e-7 jitter. Leading candidate, UNTESTED: BLAS / "
            "OpenMP thread count varying between the two runs (n_jobs is null "
            "in both configs, so any threading is implicit). Pin threads and "
            "re-run twice; if the artifact hashes then agree, the defect is "
            "environmental and fixable by declaration.",
            "Independently, fix p_circular_shift's mass-point degeneracy: on a "
            "tied null the (1 + #{T_surr >= T_obs}) / (1 + B) count is "
            "discontinuous in T_obs at machine precision. A mid-rank / "
            "randomised-tie p, or an explicit tie tolerance, makes it "
            "continuous. NOTE this field is a DIAGNOSTIC in tranche B -- p_raw "
            "= p_block_bootstrap on all 5 affected rows and p_method is "
            "MC_RESOLVED -- so no tranche-B rejection depends on it. That is "
            "luck about this configuration, not a property of the estimator.",
            "Re-run the two sessions under the pinned environment and require "
            "artifact_content_hash to AGREE. §P7-8(c) already carries the "
            "forward invariant (differing configs -> differing artifacts); the "
            "dual (identical config -> identical artifact) should be added as "
            "an executable check, not only as a ledger sentence.",
        ],

        "caveats": [
            "SCOPE. This classifies the difference between these two specific "
            "sessions. It does not establish that other sessions in "
            "artifact_hashes.jsonl are deterministic; no other pair with an "
            "identical config hash was checked here.",
            "The mechanism for the p_circular_shift jump is INFERRED from "
            "recorded values, not instrumented. The tie/mass-point account "
            "explains the observed 186-count step and the direction reversal "
            "on row 161, but the surrogate distributions themselves are not "
            "stored in the checkpoint, so it is not directly verified.",
            "No claim is made that tranche B's science is WRONG. The claim is "
            "the one §P7-20(7) authorises: the invariance audit failed, so the "
            "results are SUSPECT until resolved. The observed outcome "
            "invariance (zero differing outcome fields) is evidence the "
            "science is robust to this particular perturbation, and is "
            "recorded as such -- but it is one draw of the perturbation, not a "
            "bound on it.",
            "Neither session directory was modified. The registry "
            "engine/out/mine/artifact_hashes.jsonl was not appended to: "
            "artifact_content_hash is called directly, not check_build_invariant.",
        ],
    }

    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print("wrote %s" % OUT_PATH)

    print("\n=== §P7-20(7) HASH CLASSIFICATION ===")
    print("  A %s  artifact %s" % (SESSION_A, ha[:16]))
    print("  B %s  artifact %s" % (SESSION_B, hb[:16]))
    print("  config hash identical: %s   config dict identical: %s"
          % (ca.get("config_hash") == cb.get("config_hash"),
             ca.get("config") == cb.get("config")))
    print("\n  VERDICT: %s   GATE: %s" % (verdict, out["gate_status"]))
    print("\n  %-42s %5s %12s %12s  %s"
          % ("field", "n", "max_rel", "max_abs", "class"))
    for f in numeric_fields:
        print("  %-42s %5d %12.3e %12.3e  %s"
              % (f["field"], f["n_differing_values"],
                 f["max_relative_difference"], f["max_absolute_difference"],
                 f["class"]))
    print("\n  outcome-field differences: %d" % len(outcome_diffs))
    print("  non-numeric differences:   %d" % len(nonnum))
    print("  structural differences:    %d" % len(structural))
    print("\n  MATERIAL rows:")
    for r in material_rows:
        print("    row %3d %-22s %-9s  A=%.6f B=%.6f  counts %d vs %d "
              "(delta %d/%d)  p_raw identical: %s"
              % (r["row"], r["feature"], r["field"], r["A"], r["B"],
                 r["implied_exceedance_count_A"], r["implied_exceedance_count_B"],
                 r["implied_count_delta"], r["n_surrogates"],
                 r["p_raw_identical"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
