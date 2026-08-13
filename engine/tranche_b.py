"""TRANCHE B -- THE DECLARATION, PREPARED AND NOT RUN.

§P7-14(d): *"AUTHORIZED NOW: the §P7-2(b) precondition builds... GATED: B's RUN
begins when F9-19/G-M1 arm (ii) reports."*  §P7-15(b) then confirms the run gate as
read: *"B's run is gated only on (i) its own §P7-2(b) precondition builds and (ii)
per-statistic G-M1 clearance for the NEW statistics and aggregations B introduces."*

So this module PREPARES the declaration -- the count enumeration, the strata file,
the frozen config -- and REFUSES to execute it. `main()` prints the declaration;
`run()` raises `TrancheBRunGated`. The refusal is in code rather than in a comment
because a driver that could be run by typing one flag is a driver that will be.

WHAT A "PREPARED DECLARATION" HAS TO GET RIGHT
----------------------------------------------
Everything that is frozen BEFORE a run is seen (S-9, §P6-3 rules 3 and 5): the test
count as an exact integer, the stratum partition and its `m_s`/`q_s`, the budget
identity `sum_s m_s q_s = m q`, and the config hash that binds them. All four are
produced here and none of them may be edited after B's first result is looked at --
that is what makes the difference between a declaration and a description.

THE COUNT RECONCILIATION, AND THE DISCREPANCY IT SURFACES
---------------------------------------------------------
§P7-10(c) states the declared integer as *"Kepler's ~1,000 (17 second-moment + 34
omnibus + ~161 mark + 68 already-declared §P5-5(3) ladder + 32 linear ladder + 17
two-stage + 31 bilinear + declared overhead)"*, and MINING_CATALOG's Tranche B
section gives the same list. **Those seven named items sum to 360, not to ~1,000.**
`enumerate_declared()` computes the sum from the item list rather than quoting the
headline, and reports the gap. **This module does not resolve it.** Choosing a
denominator is an adjudication with direct consequences for every BH threshold in the
tranche, and §P7-10(c) reserves the integer to the Popper seat ("the exact integer
frozen in the config hash before the run"). The build's job is to make the
discrepancy impossible to miss, which is what `COUNT_DISCREPANCY_NOTE` is for.

THE 68 ARE RE-OCCUPIED, NOT NEW
-------------------------------
§P7-3: *"68 of it is the §P5-5(3) tranche-3 ladder, already declared in 2026-08-11.
Already-declared is not free -- it still occupies its 68 slots in this session's BH
denominator -- but it may not be presented as new scope, and the report must say
which 68 they are."*  Handled exactly as Tranche A handled its 550 (`tranche_a.py`):
carried in the denominator, excluded from the new-scope figure, and labelled.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import types

from . import (circstat, clocks, floors, marks_ext, mine as M,
               mine_session as ms, observer, s15c, splits, stability,
               strata as strata_mod)

TRANCHE_B_ID = "TRANCHE-B-STATISTICS"
STRATA_FILE = os.path.join("engine", "configs", "strata_tranche_b.json")
RESULTS_JSON = "results_tranche_b_declaration.json"

# §P7-10(c) / MINING_CATALOG Tranche B. Each item carries its own arithmetic, its
# catalog source, and whether it is NEW scope or a RE-OCCUPIED already-declared slot.
# `built` records whether this build implements the arm -- an honest scope flag, not
# a status field: a declared stratum with no implementation would enter B's
# denominator and produce zero rows, and §P6-3 rule 3 counts those as non-rejections
# against their declared m_s. That is legal and it is not free, so it is printed.
COMPOSITION = (
    {"key": "second_moment", "catalog": "F9-01", "stratum": "tb_second_moment",
     "test_kind": "moment2", "n": 17, "new": True, "built": True,
     "arithmetic": "17 cyclic features x 1 new statistic",
     "note": "the second circular moment, as the 2-df score on the doubled angle"},
    {"key": "omnibus", "catalog": "F9-04", "stratum": "tb_omnibus",
     "test_kind": "omnibus", "n": 34, "new": True, "built": True,
     "arithmetic": "17 cyclic features x 2 statistics (Kuiper V, Watson U^2)",
     "note": "shape sensitivity, NOT band coverage (§P7-3(2))"},
    {"key": "mark_axis", "catalog": "F9-10", "stratum": "tb_mark_axis",
     "test_kind": "markx", "n": 161, "new": True, "built": True,
     "arithmetic": "7 marks x 23 features",
     "note": "its OWN §P6-3 stratum, reallocation justified by measurement "
             "(§P7-11(c)); sub-daily arm gated on F7-01/02/03 (§P7-3(3))"},
    {"key": "ladder_reoccupied", "catalog": "F9-05", "stratum": "tb_ladder",
     "test_kind": "ladder", "n": 68, "new": False, "built": False,
     "arithmetic": "the §P5-5(3) tranche-3 ladder, DECLARED 2026-08-11",
     "note": "RE-OCCUPIED, not new scope (§P7-3). Occupies its 68 slots in the BH "
             "denominator and may not be presented as new."},
    {"key": "linear_ladder", "catalog": "F2-18..F2-25", "stratum": "tb_linear_ladder",
     "test_kind": "ladder_linear", "n": 32, "new": True, "built": False,
     "arithmetic": "8 linear cyclic features x 4 new rungs",
     "note": "not a new estimator; inherits the existing ladder discipline "
             "(§P7-3(4)), so it carries no S-17 recovery demand of its own"},
    {"key": "two_stage", "catalog": "F10-14", "stratum": "tb_two_stage",
     "test_kind": "two_stage", "n": 17, "new": True, "built": False,
     "arithmetic": "1 model x 17 cyclic features",
     "note": "K-088's object; inherits §P5-3's two mandated null repairs and may "
             "NOT be reported as a fresh design (§P7-3(5))"},
    {"key": "bilinear", "catalog": "F10-08", "stratum": "tb_bilinear",
     "test_kind": "bilinear", "n": 31, "new": True, "built": False,
     "arithmetic": "the full bilinear perigee-syzygy interaction x 31 lags",
     "note": "priced as a transform per §P7-5(5)"},
    {"key": "observer_controls", "catalog": "F7-01/02/03", "stratum": "tb_controls",
     "test_kind": "glm", "n": None, "new": True, "built": True,
     "arithmetic": "the count-path observer control features x 1 lag "
                   "(computed from engine/observer.py, not retyped)",
     "note": "PRICED, on §P7-2(a)'s own reasoning that an unpriced control is an "
             "unaudited channel -- which sits in tension with MINING_CATALOG "
             "F7-01's 'Price: 31 as a covariate; 0 as a control'. Flagged, priced "
             "the conservative way, not resolved here."},
)

POPPER_HEADLINE_COUNT = 1000
COUNT_DISCREPANCY_NOTE = (
    "§P7-10(c) declares Tranche B at ~1,000 and itemises it as '17 second-moment + "
    "34 omnibus + ~161 mark + 68 already-declared ladder + 32 linear ladder + 17 "
    "two-stage + 31 bilinear + declared overhead'. THOSE SEVEN ITEMS SUM TO 360. "
    "The residual 'declared overhead' would therefore have to be ~640 -- 178% of the "
    "named scope -- which is not an overhead, it is a majority. MINING_CATALOG's own "
    "Tranche B section gives the same seven numbers and the same ~1,000 headline, so "
    "the gap originates in the catalog and was carried into the ruling verbatim. "
    "THIS BUILD DOES NOT RESOLVE IT: the declared integer is reserved to the Popper "
    "seat (§P7-10(c): 'the exact integer frozen in the config hash before the run'), "
    "and the choice moves every BH threshold in the tranche. Reported, flagged, and "
    "left open.")

RUN_GATE_NOTE = (
    "§P7-14(d) + §P7-15(b): Tranche B's BUILD is authorized and its RUN is gated on "
    "(i) the §P7-2(b) precondition builds and (ii) PER-STATISTIC G-M1 clearance for "
    "the new statistics and aggregations B introduces. `engine/recovery_b.py` "
    "produces the evidence for (ii) on SIMULATED catalogues; it does not grant the "
    "clearance, which is the Popper seat's. This module refuses to run either way.")


class TrancheBRunGated(RuntimeError):
    """Someone tried to execute a declaration whose run is gated."""


# ------------------------------------------------------------- enumeration ----
def observer_control_count():
    """How many COUNT-PATH observer control tests B declares. Derived, not retyped."""
    import numpy as np
    n_days = 400
    rng = np.random.default_rng(0)
    d = np.sort(rng.uniform(0, n_days, size=4000))
    marks = {"day_float": d, "mag": 4.5 + rng.exponential(0.4, d.size),
             "day": np.floor(d).astype(np.int64)}
    feats = observer.observer_features(marks, n_days)
    return len(observer.count_path_features(feats))


def enumerate_declared():
    """The count reconciliation. Sums the ITEMS; never quotes the headline."""
    items = []
    for it in COMPOSITION:
        n = it["n"]
        if n is None:
            n = observer_control_count()
        items.append(dict(it, n=int(n)))
    new = sum(i["n"] for i in items if i["new"])
    reoccupied = sum(i["n"] for i in items if not i["new"])
    total = new + reoccupied
    built = sum(i["n"] for i in items if i["built"])
    return {
        "items": items,
        "n_new_scope": new,
        "n_reoccupied_already_declared": reoccupied,
        "bh_denominator_m": total,
        "n_in_arms_this_build_implements": built,
        "n_in_arms_declared_but_not_built": total - built,
        "popper_headline": POPPER_HEADLINE_COUNT,
        "gap_vs_headline": POPPER_HEADLINE_COUNT - total,
        "agrees_with_headline": abs(POPPER_HEADLINE_COUNT - total) <= 50,
        "discrepancy_note": COUNT_DISCREPANCY_NOTE,
        "reoccupied_identification": (
            "the 68 are F9-05, the §P5-5(3) tranche-3 harmonic ladder declared "
            "2026-08-11. §P7-3 requires the report to say WHICH 68 they are; this "
            "is that sentence."),
        "not_built_flag": (
            "arms declared but NOT implemented in this build: %s. A declared "
            "stratum with no rows is legal -- §P6-3 rule 3 counts errored/absent "
            "tests as NON-REJECTIONS against their declared m_s -- but it is not "
            "free, and B must either build them or re-declare without them before "
            "the run."
            % ", ".join("%s (%s, %d)" % (i["key"], i["catalog"], i["n"])
                        for i in items if not i["built"])),
    }


# ------------------------------------------------------------ strata file -----
def strata_document(q=M.FDR_Q):
    """The declared §P6-3 partition for Tranche B. Flat q_s: no reallocation, no prior.

    §P6-3 rule 2: any DEPARTURE from the flat budget is a prior that must carry its
    justification in the file itself, before the run. B declares none -- every
    stratum gets the declared q -- so there is nothing to justify and the budget
    identity `sum_s m_s q_s = m q` holds by construction. §P7-11(c) reallocated the
    MARK arm's q_s in Tranche A on measured grounds; B does not inherit that
    reallocation, because B's mark arm is a different arm (7 marks, event-time) and
    the measurement that justified the old allocation is not a measurement of this
    one.
    """
    enum = enumerate_declared()
    strata = []
    for it in enum["items"]:
        strata.append({
            "name": it["stratum"],
            "feature_family": None,
            "test_kind": it["test_kind"],
            "region": None,
            "m_s": int(it["n"]),
            "q_s": float(q),
            "note": "",
        })
    return {
        "note": ("TRANCHE B (HYPOTHESIS_LEDGER.md §P7-3, §P7-10(c), §P7-14(d), "
                 "§P7-15(b)), declared BEFORE the run and frozen in the config hash "
                 "by its own sha256 (§P6-3(3)+(5)). BH denominator m = %d = %d new "
                 "scope + %d re-occupied already-declared slots (the F9-05 "
                 "§P5-5(3) tranche-3 ladder, 2026-08-11). q_s is FLAT at the "
                 "declared q in every stratum: no budget is reallocated, so there "
                 "is no prior to justify under §P6-3(2). The mark axis sits in its "
                 "OWN stratum (§P7-11(c)) rather than in the v1 `mark` stratum, "
                 "because 7 event-time marks are not the 2 day-binned marks that "
                 "allocation was measured on. COUNT DISCREPANCY, CARRIED NOT "
                 "RESOLVED: %s"
                 % (enum["bh_denominator_m"], enum["n_new_scope"],
                    enum["n_reoccupied_already_declared"], COUNT_DISCREPANCY_NOTE)),
        "q": float(q),
        "catch_all": "tb_second_moment",
        "strata": strata,
    }


def write_strata_file(path=STRATA_FILE, q=M.FDR_Q):
    """Write the declared partition and verify the budget identity it must satisfy."""
    doc = strata_document(q)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    raw = json.dumps(doc, indent=1).encode("utf-8")
    with open(path, "wb") as fh:
        fh.write(raw)
    part = strata_mod.load_partition(path, q=q)
    strata_mod.assert_budget_identity(part["strata"], q)
    return {"path": path.replace("\\", "/"),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "m": int(part["m"]), "q": float(q),
            "n_strata": len(part["strata"])}


# ---------------------------------------------------------------- config ------
def build_args(seed=20260813, jobs=8, subdaily=False, strata=STRATA_FILE):
    """The frozen Tranche B config inputs. One declared value each (S-9)."""
    return types.SimpleNamespace(
        mag_target=4.5, dlat=1.0, dlon=1.0, explore_frac=0.70,
        data_dir="data/comcat_world", no_download=True,
        seed=int(seed), tranche1=False, controls=False, strata=strata,
        ladder=False, gpd=False, regsum=False, regions=False, jobs=int(jobs),
        tranche_b=True, tb_second_moment=True, tb_omnibus=True, tb_mark_axis=True,
        tb_subdaily=bool(subdaily),
        tb_declared=enumerate_declared()["bh_denominator_m"])


def frozen_config(seed=20260813, subdaily=False, strata=STRATA_FILE):
    """The config that would be hashed if B ran. Built, hashed, and NOT executed."""
    cfg = ms.build_config(build_args(seed=seed, subdaily=subdaily, strata=strata),
                          ms.OVERNIGHT)
    enum = enumerate_declared()
    assert cfg["tranche_b"]["enabled"] is True
    assert cfg["strata"]["m"] == enum["bh_denominator_m"], (
        "strata file declares m = %s, enumeration says %s"
        % (cfg["strata"]["m"], enum["bh_denominator_m"]))
    return cfg


def run(*_a, **_k):
    """REFUSED. B's run is gated; this exists so the refusal is executable."""
    raise TrancheBRunGated(RUN_GATE_NOTE)


# ----------------------------------------------------------- the declaration --
def declaration(seed=20260813, subdaily=False, write=True):
    """Everything that must be frozen before B runs, assembled in one artifact."""
    enum = enumerate_declared()
    sf = write_strata_file() if write else {"path": STRATA_FILE, "sha256": None,
                                            "m": enum["bh_denominator_m"]}
    cfg = frozen_config(seed=seed, subdaily=subdaily) if write else None
    sweep = s15c.sweep()
    return {
        "id": TRANCHE_B_ID,
        "title": "EQ-24 v2 Tranche B -- the statistic tranche, DECLARATION PREPARED "
                 "AND NOT RUN",
        "banner": M.GENERATOR_NOT_EVIDENCE,
        "standing_warning_eq24": M.MIGNAN_BROCCARDO,
        "prepared_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "run_status": "GATED -- NOT RUN",
        "run_gate": RUN_GATE_NOTE,
        "count_reconciliation": enum,
        "strata_file": sf,
        "config": cfg,
        "config_hash": (splits.config_hash(cfg) if cfg else None),
        "floors": {
            "count_path_alpha": floors.ALPHA_TRANCHE_B,
            "count_path_vif": floors.MEASURED_VIF_DF2_PHASE,
            "count_path_A_min_at_46585": floors.a_min(
                floors.MEASURED_VIF_DF2_PHASE, floors.ALPHA_TRANCHE_B, 46585),
            "mark_axis_formula": "rho_min = sqrt(VIF_mark) * 4.732 / sqrt(n - 1)",
            "mark_axis_vif_fallback": floors.VIF_MARK_FALLBACK,
            "mark_axis_rho_min_at_46585": floors.rho_min(46585),
            "s15_headline_clauses": [
                "(i) the floor formula, its VIF and the tranche's own alpha are "
                "printed in the headline, not a footnote",
                "(ii) the fraction of declared tests UNMEASURABLE for amplitude at "
                "the declared effect of interest is PRINTED, not omitted -- on the "
                "count path at a ~15% floor this is essentially all of them",
                "(iii) DETECTION results are adjudicated against surrogate-"
                "calibrated nulls and are NOT gated by the amplitude floor, which "
                "is why the tranche survives its own floor",
                "(iv) a null from this tranche bounds non-sinusoidal and "
                "second-moment structure at ~15%, NOT at the ~5.6% this programme "
                "is accustomed to quoting",
                "(v) the mark-axis floor is stated separately, in correlation "
                "units, with its own VIF_mark",
                "PER-FEATURE AMPLITUDE QUOTES ON THE COUNT PATH ARE UNRESOLVED AND "
                "ARE NOT PRINTED (§P7-10(c)).",
            ],
        },
        "s15c_sweep": {k: v for k, v in sweep.items() if k != "rows"},
        "s15c_unmeasurable": sweep["unmeasurable"],
        "observer_controls": {
            "required_for_subdaily": list(observer.REQUIRED_FOR_SUBDAILY),
            "rule": observer.SUBDAILY_GATE_RULE,
            "banner": observer.OBSERVER_CONTROL_BANNER,
            "subdaily_arm_requested": bool(subdaily),
        },
        "f8_15_random_clock": {
            "rule_id": clocks.F8_15_RULE_ID, "modes": list(clocks.F8_15_MODES),
            "mandatory_note": clocks.MANDATORY_NOTE,
            "status": ("BUILT and FROZEN now per §P7-4, whose consumer is Tranche C. "
                       "Tranche B declares no clock scan, so no F8-15 control is "
                       "owed by B -- `clocks.assert_random_clock_control` is what "
                       "will refuse a C scan that forgets one."),
        },
        "mark_axis": {
            "marks": list(marks_ext.MARK_NAMES),
            "definitions": dict(marks_ext.MARK_DEFINITIONS),
            "subdaily_note": marks_ext.SUBDAILY_NOTE,
        },
        "resampling_stability_condition": {
            "rule": "§P7-14(b), MANDATORY on every BH survivor and the "
                    "max-statistic detection",
            "null_baseline": dict(stability.NULL_BASELINE),
            "thresholds": dict(stability.THRESHOLDS),
            "label_not_kill": stability.LABEL_NOT_KILL,
            "lower_bound": stability.LOWER_BOUND_SENTENCE,
            "r4_not_sufficient": stability.R4_NOT_SUFFICIENT,
        },
        "per_statistic_g_m1": {
            "F9-01": "two-lobed positive control AND the sinusoid negative control",
            "F9-04": "narrow-arc positive and day-lattice negative, band by band",
            "F9-10": "G-M1 arm (i) re-run on the MARK PATH specifically",
            "F10-14": "inherits §P5-3's null repairs; not re-designed",
            "harness": "engine/recovery_b.py -- SIMULATION ONLY",
            "adjudication": "the Popper seat's; this build supplies evidence only",
        },
        "what_none_of_this_licenses": (
            "any run. This is a DECLARATION prepared before the fact, on no data at "
            "all. It licenses nothing, bounds nothing, and may not be entered for "
            "or against any ledger entry."),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=20260813)
    ap.add_argument("--subdaily", action="store_true")
    ap.add_argument("--json", default=RESULTS_JSON)
    a = ap.parse_args(argv)
    d = declaration(seed=a.seed, subdaily=a.subdaily)
    e = d["count_reconciliation"]
    print("=" * 78)
    print("TRANCHE B -- DECLARATION PREPARED, RUN GATED")
    print("=" * 78)
    print("| item | catalog | tests | scope | built |")
    print("| --- | --- | ---: | --- | --- |")
    for it in e["items"]:
        print("| %s | %s | %d | %s | %s |"
              % (it["key"], it["catalog"], it["n"],
                 "NEW" if it["new"] else "re-occupied",
                 "yes" if it["built"] else "NO"))
    print("")
    print("new scope                : %d" % e["n_new_scope"])
    print("re-occupied (F9-05, 68)  : %d" % e["n_reoccupied_already_declared"])
    print("BH denominator m         : %d" % e["bh_denominator_m"])
    print("Popper headline          : ~%d   -> gap %d"
          % (e["popper_headline"], e["gap_vs_headline"]))
    print("strata file              : %s (sha256 %s, m = %s)"
          % (d["strata_file"]["path"], str(d["strata_file"]["sha256"])[:16],
             d["strata_file"]["m"]))
    print("config hash              : %s" % d["config_hash"])
    print("count-path floor A_min   : %.4f at N = 46,585, alpha = %.1e"
          % (d["floors"]["count_path_A_min_at_46585"],
             d["floors"]["count_path_alpha"]))
    print("mark-axis floor rho_min  : %.4f at n = 46,585, VIF_mark = %.3f"
          % (d["floors"]["mark_axis_rho_min_at_46585"],
             d["floors"]["mark_axis_vif_fallback"]))
    print("S-15(c) unmeasurable     : %d entries (cut at period > %.0f d)"
          % (d["s15c_sweep"]["n_unmeasurable_by_window"],
             d["s15c_sweep"]["cut_period_days"]))
    print("")
    print("RUN STATUS: %s" % d["run_status"])
    print(RUN_GATE_NOTE)
    print("")
    print("COUNT DISCREPANCY: " + COUNT_DISCREPANCY_NOTE)
    with open(a.json, "w", encoding="utf-8") as fh:
        json.dump(d, fh, indent=1, default=float)
    print("\nwrote %s" % a.json)
    return d


if __name__ == "__main__":
    main()
