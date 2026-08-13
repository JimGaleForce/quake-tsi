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

THE COUNT, RECONCILED EXACTLY AT §P7-16 -- FOUR CLASSES, NOT ONE NUMBER
-----------------------------------------------------------------------
§P7-10(c) declared B at ~1,000 and itemised it as seven arms that sum to **360**.
The build reported that gap rather than choosing a denominator; §P7-16 resolved it,
and the resolution is a taxonomy rather than an arithmetic correction. What used to
be one undifferentiated "declared count" is four things with four different
accounting consequences:

  **PRICED -- 189.** 17 second-moment (F9-01) + 34 omnibus (F9-04) + 138 mark axis
  (F9-10). Enters the BH denominator; can be rejected. This is `m`.

  **DEFERRED -- 148.** The four arms this build does not implement (F9-05 68,
  F2-18..F2-25 32, F10-14 17, F10-08 31). They get **no strata here** and return as
  **ONE B-2 declaration**. A declared stratum with no implementation consumes budget
  and produces nothing, taxing the arms that were built for scope that was never
  going to run -- and §P7-16 ruled over-declaration a **defect, not a virtue**
  (a 5.3x needless threshold tax, and a resolvability count about a fiction).

  **UNPRICED -- 9.** The F7 observer controls, under §P7-16's new general rule:
  *a control is PRICED if it can produce a survivor mistakable for a finding about
  the world, and UNPRICED if it can only calibrate a reference or condemn our own
  instrument.* They are declared at `m_s = 0` -- present so a control row has a
  stratum to route to, outside the denominator so it can neither be rejected nor
  consume budget. (The build's F9-20-style instinct to price them was right for
  SURVIVABLE controls and wrong for these; the distinction is now the rule.)

  **DE-DUPLICATED -- 23.** `log_moment` is rank-identical to `mag` under both
  declared mark statistics, proved numerically. **De-duplicate, do not keep**: a
  surviving duplicate reads as phantom replication.

189 + 148 + 9 + 23 = 369. **The ~1,000 headline is WITHDRAWN** -- an unsummed
aggregate published beside its own itemisation -- and may not be quoted.

THE ASSERTION THAT MAKES THIS CHECKABLE
---------------------------------------
§P7-16 mandated `sum_s m_s == m` beside the budget identity
(`strata.assert_partition_total`). `assert_budget_identity` derives `m` from the
strata and so is satisfied by any consistent table, including one whose total
disagrees with the integer the tranche was declared at -- which is precisely the
failure that produced the ~1,000. The partition file now carries `"m": 189` and
refuses to load if the table does not add up to it.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import types

from . import (circstat, clocks, dispositions as disp, floors, marks_ext, mine as M,
               mine_session as ms, observer, s15c, splits, stability,
               strata as strata_mod)

# THE MARK AXIS'S FEATURE COUNT, and the §P7-18 correction to it.
#
# MINING_CATALOG F9-10 prices "7 marks x 23+ features". Three of those 23 are
# DOWNLOAD features (`mine.DOWNLOADS`), and Tranche B's frozen config declares
# `no_download=True` -- so the catalog's 23 was never buildable under B's own
# config, and 18 priced mark slots would have entered the denominator and never
# executed. §P7-18 lowered the axis to the 20 features B actually builds and
# DEFERRED the three BY NAME.
MARK_AXIS_FEATURES_CATALOG = 23
MARK_AXIS_FEATURES = 20

# Deferred to B-2 BY NAME (§P7-18 requires the names, not a count).
DEFERRED_FEATURES_B2 = (
    {"feature": "Ap_geomagnetic", "family": 3, "source": "GFZ daily Ap index",
     "why_absent": "DOWNLOAD feature; B declares no_download=True",
     "s15c": None},
    {"feature": "F107_solar_flux", "family": 3,
     "source": "Penticton F10.7 cm solar radio flux",
     "why_absent": "DOWNLOAD feature; B declares no_download=True",
     "s15c": ("AND it is UNMEASURABLE-BY-WINDOW on its own band anyway: F10.7's "
              "11 y solar cycle gives 1.92 cycles in the 7,716 d window, below "
              "S-15(c)'s 3. This is the feature whose z = 32 corpse the clause "
              "retroactively explains -- deferring it costs the tranche nothing it "
              "could have identified.")},
    {"feature": "length_of_day", "family": 3, "source": "IERS finals.all LOD",
     "why_absent": "DOWNLOAD feature; B declares no_download=True",
     "s15c": None},
)

DEFERRAL_REASON_CLASS = "config-determined-not-result-determined"
RATCHET_NOTE = (
    "THE RATCHET (§P7-18). The declared integer FIXES AT CONFIG-HASH FREEZE. This "
    "re-issue from 189 to 171 is legitimate precisely because its reason is "
    "CONFIG-DETERMINED, NOT RESULT-DETERMINED: the three features were dropped "
    "because B's own frozen `no_download=True` cannot build them, a fact knowable -- "
    "and known -- before any datum was touched, and B has not been run. A change "
    "with the same arithmetic but a result-determined reason would be S-9's forking "
    "path wearing a correction's clothes. The distinction is the whole of the "
    "ratchet, and it is recorded on the artifact rather than in a commit message so "
    "it travels with the number.")

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
     "test_kind": "moment2", "n": 17, "class": "PRICED", "built": True,
     "arithmetic": "17 cyclic features x 1 new statistic",
     "note": "the second circular moment, as the 2-df score on the doubled angle"},
    {"key": "omnibus", "catalog": "F9-04", "stratum": "tb_omnibus",
     "test_kind": "omnibus", "n": 34, "class": "PRICED", "built": True,
     "arithmetic": "17 cyclic features x 2 statistics (Kuiper V, Watson U^2)",
     "note": "shape sensitivity, NOT band coverage (§P7-3(2))"},
    {"key": "mark_axis", "catalog": "F9-10", "stratum": "tb_mark_axis",
     "test_kind": "markx", "n": 120, "class": "PRICED", "built": True,
     "arithmetic": "6 SCORED marks x 20 BUILDABLE features (7 declared marks less "
                   "the de-duplicated `log_moment` (§P7-16); 23 catalog features "
                   "less the 3 DOWNLOAD features B's own config cannot build, "
                   "deferred to B-2 by name (§P7-18))",
     "note": "its OWN §P6-3 stratum, reallocation justified by measurement "
             "(§P7-11(c)); sub-daily arm gated on F7-01/02/03 (§P7-3(3))"},

    # ---- DEFERRED: returns later as ONE B-2 declaration, NOT as strata here ----
    {"key": "ladder_reoccupied", "catalog": "F9-05", "stratum": None,
     "test_kind": "ladder", "n": 68, "class": "DEFERRED", "built": False,
     "arithmetic": "the §P5-5(3) tranche-3 ladder, DECLARED 2026-08-11",
     "note": "already-declared; deferred to B-2 with the other unbuilt arms"},
    {"key": "linear_ladder", "catalog": "F2-18..F2-25", "stratum": None,
     "test_kind": "ladder_linear", "n": 32, "class": "DEFERRED", "built": False,
     "arithmetic": "8 linear cyclic features x 4 new rungs",
     "note": "not a new estimator; inherits the existing ladder discipline "
             "(§P7-3(4)), so it carries no S-17 recovery demand of its own"},
    {"key": "two_stage", "catalog": "F10-14", "stratum": None,
     "test_kind": "two_stage", "n": 17, "class": "DEFERRED", "built": False,
     "arithmetic": "1 model x 17 cyclic features",
     "note": "K-088's object; inherits §P5-3's two mandated null repairs and may "
             "NOT be reported as a fresh design (§P7-3(5))"},
    {"key": "bilinear", "catalog": "F10-08", "stratum": None,
     "test_kind": "bilinear", "n": 31, "class": "DEFERRED", "built": False,
     "arithmetic": "the full bilinear perigee-syzygy interaction x 31 lags",
     "note": "priced as a transform per §P7-5(5)"},

    # ---- UNPRICED: F7 controls, under §P7-16's new general rule ----------------
    {"key": "observer_controls", "catalog": "F7-01/02/03", "stratum": "tb_controls",
     "test_kind": "glm", "n": None, "class": "UNPRICED", "built": True,
     "arithmetic": "the count-path observer control features x 1 lag "
                   "(computed from engine/observer.py, not retyped)",
     "note": "UNPRICED under §P7-16's general rule (below). m_s = 0: declared so "
             "the rows have a stratum to route to, outside the priced denominator "
             "so they cannot be rejected and cannot consume budget."},
)

DEFERRED_LABEL = "B-2"
DEFERRED_NOTE = (
    "§P7-16: the four arms this build does not implement -- F9-05 (68, "
    "already-declared), F2-18..F2-25 (32), F10-14 (17), F10-08 (31), 148 tests in "
    "total -- are DEFERRED and return later as ONE B-2 declaration. They are NOT "
    "declared as empty strata here: a declared stratum with no implementation "
    "consumes budget through its m_s and produces nothing, which taxes the arms "
    "that did get built for scope that was never going to run. Deferring them is "
    "the honest denominator; re-declaring them together is the honest scope.")

UNPRICED_RULE = (
    "§P7-16, the general rule, stated in full because it generalises past this "
    "tranche: A CONTROL IS UNPRICED IF IT CAN ONLY CALIBRATE A REFERENCE OR CONDEMN "
    "OUR OWN INSTRUMENT. Multiplicity is owed on rejections we might MAKE, not on "
    "checks that can only unmake them -- the same principle §P7-5(3) applied to "
    "LORO. The F9-20 battery is priced because its controls are SURVIVABLE: a "
    "control there can 'win' and its survivor count is read as a false-positive "
    "rate, so it must face the same threshold as the real arm or F10-25's ratio is "
    "meaningless. The F7 observer controls are not that: a hit on `obs_mc_drift_365d` "
    "does not become a finding, it condemns a co-moving science row. Nothing about "
    "them can promote anything, so they owe nothing. "
    "The build's own instinct -- price them, on §P7-2(a)'s 'an unpriced control is "
    "an unaudited channel' -- was RIGHT FOR SURVIVABLE CONTROLS and wrong here; the "
    "distinction is now the rule rather than a case."
)

POPPER_RULED_DENOMINATOR = 171
POPPER_RULED_DENOMINATOR_SUPERSEDED = 189      # §P7-16, superseded by §P7-18
POPPER_HEADLINE_COUNT = 1000
HEADLINE_RESOLUTION_NOTE = (
    "§P7-10(c) declared Tranche B at ~1,000 and itemised it as '17 second-moment + "
    "34 omnibus + ~161 mark + 68 already-declared ladder + 32 linear ladder + 17 "
    "two-stage + 31 bilinear + declared overhead'. Those seven items sum to 360, so "
    "the residual 'declared overhead' would have had to be ~640 -- 178% of the named "
    "scope. §P7-16 RESOLVED it: the overhead was never scope. What is left is 189 "
    "priced, 148 deferred to a single B-2 declaration, 9 unpriced controls and 23 "
    "de-duplicated. The ~1,000 headline is SUPERSEDED and may not be quoted.")

RUN_GATE_NOTE = (
    "§P7-14(d) + §P7-15(b): Tranche B's BUILD is authorized and its RUN is gated on "
    "(i) the §P7-2(b) precondition builds and (ii) PER-STATISTIC G-M1 clearance for "
    "the new statistics and aggregations B introduces. As of §P7-16 the run is ALSO "
    "gated on (iii) the recovery-versus-amplitude curve §P7-15(a) assigned at price "
    "0 -- IN FLIGHT -- and (iv) the Popper seat's ruling on the arm (i) anomaly (80% "
    "recovery at >= 2x the floor is what the formula predicts at 1x, so either the "
    "plants span a range or the global floor formula is ~2x anti-conservative). "
    "`engine/recovery_b.py` produces the evidence for (ii) on SIMULATED catalogues; "
    "it does not grant the clearance. This module refuses to run either way.")


class TrancheBRunGated(RuntimeError):
    """Someone tried to execute a declaration whose run is gated."""


# ------------------------------------------------------------- enumeration ----
# What B's frozen config actually builds, derived from the generators rather than
# retyped. `no_download=True` and `tranche1=False`, so: 17 ephemeris cyclic features
# at lag 0, 3 family-4 catalogue features across the full lag grid, no F9-20 controls.
def session_feature_inventory():
    """The feature names a Tranche B session would carry. No data required."""
    import datetime as _d
    import numpy as _np
    eph = M.ephemeris_features(_d.datetime(2000, 1, 1), 400)
    cyclic = [f.name for f in eph if f.kind == "phase" or f.period_hint]
    other = ["b_value_90d", "deep_fraction_30d", "mean_depth_30d"]
    rng = _np.random.default_rng(0)
    d = _np.sort(rng.uniform(0, 400, size=4000))
    marks = {"day_float": d, "mag": 4.5 + rng.exponential(0.4, d.size),
             "day": _np.floor(d).astype(_np.int64)}
    ctl = [f.name for f in
           observer.count_path_features(observer.observer_features(marks, 400))]
    return {"cyclic": cyclic, "other_science": other, "controls": ctl}


def row_enumeration(mark_axis_features=None):
    """§P7-17: every row B's session would execute, against the four dispositions.

    Two views, because they answer two different questions and conflating them is
    how a shortfall hides:

      DECLARED VIEW   the mark axis at its DECLARED feature count
                      (MARK_AXIS_FEATURES) -- this is what the partition prices, and
                      what `assert_declared_row_count` checks m against.
      EXECUTED VIEW   the mark axis at the features B's FROZEN CONFIG actually
                      builds. `no_download=True`, so the download-derived features
                      F9-10's price assumed are not there.
    """
    inv = session_feature_inventory()
    n_exec_features = len(inv["cyclic"]) + len(inv["other_science"])
    v1_marks = ("mag", "depth")
    declared = disp.enumerate_session_rows(
        inv["cyclic"], inv["other_science"], inv["controls"],
        n_period_peaks=int(ms.OVERNIGHT["n_peaks"]),
        other_lags=len(ms.OVERNIGHT["lags"]),
        mark_axis_features=(mark_axis_features or MARK_AXIS_FEATURES),
        scored_marks=marks_ext.SCORED_MARK_NAMES)
    executed = disp.enumerate_session_rows(
        inv["cyclic"], inv["other_science"], inv["controls"],
        n_period_peaks=int(ms.OVERNIGHT["n_peaks"]),
        other_lags=len(ms.OVERNIGHT["lags"]),
        mark_axis_features=n_exec_features,
        scored_marks=marks_ext.SCORED_MARK_NAMES)
    shortfall = (declared["totals"][disp.DECLARED]
                 - executed["totals"][disp.DECLARED])
    supp = suppression_map_declaration(inv, v1_marks)
    return {
        "feature_inventory": inv,
        "declared_view": declared,
        "executed_view": executed,
        "mark_axis_features_declared": (mark_axis_features or MARK_AXIS_FEATURES),
        "mark_axis_features_executed": n_exec_features,
        "declared_minus_executed": int(shortfall),
        "shortfall_flag": (None if shortfall == 0 else
                           ("F9-10 is priced at %d features (MINING_CATALOG's "
                            "'7 marks x 23+ features'), but B's frozen config sets "
                            "`no_download=True` and therefore builds %d. %d priced "
                            "mark slots would never execute. That is LEGAL -- §P6-3 "
                            "rule 3 counts an unexecuted test as a NON-REJECTION "
                            "against its declared m_s -- and it is not free: those "
                            "slots tax the rows that do run. It is REPORTED, not "
                            "resolved: lowering the mark axis to %d moves the "
                            "integer §P7-16 just fixed, and the integer is the "
                            "Popper seat's."
                            % ((mark_axis_features or MARK_AXIS_FEATURES),
                               n_exec_features, shortfall,
                               len(marks_ext.SCORED_MARK_NAMES) * n_exec_features))),
        "n_genuinely_new": declared["n_genuinely_new_beyond_the_declared_statistics"],
        "genuinely_new_note": declared["genuinely_new_note"],
        "suppression_map": supp,
    }


def suppression_map_declaration(inv=None, v1_marks=("mag", "depth")):
    """§P7-18's suppression map, PROVED on representative rows before the run.

    The map is produced by the same `dispositions.suppression_map` the session will
    run, on one row pair per (feature, mark) the session will emit -- so the
    declaration's map and the session's map cannot disagree about what was
    suppressed or why.
    """
    inv = inv or session_feature_inventory()
    feats = list(inv["cyclic"]) + list(inv["other_science"])
    rows = []
    for f in feats:
        for mk in v1_marks:
            rows.append({"test": "spearman", "feature": f, "mark": mk,
                         "lag": None, "p_raw": None})
            rows.append({"test": "spearman", "feature": f, "mark": mk,
                         "lag": None, "mark_axis": "F9-10", "subdaily": False,
                         "p_raw": None})
    kept, suppressed, mp = disp.suppression_map(rows, subdaily=False)
    # and the counterfactual that proves the suppression is CONDITIONAL, not blanket
    rows_sd = []
    for f in feats[:1]:
        for mk in v1_marks[:1]:
            rows_sd.append({"test": "spearman", "feature": f, "mark": mk,
                            "lag": None, "p_raw": None})
            rows_sd.append({"test": "spearman", "feature": f, "mark": mk,
                            "lag": None, "mark_axis": "F9-10", "subdaily": True,
                            "p_raw": None})
    _k2, s2, mp2 = disp.suppression_map(rows_sd, subdaily=True)
    return {
        "n_suppressed": mp["n_suppressed"],
        "active": mp["active"],
        "rule": mp["rule"],
        "shape_fields": mp["shape_fields"],
        "n_features": len(feats), "v1_marks": list(v1_marks),
        "entries_sample": mp["entries"][:4],
        "entries_all_proved_same_shape": all(
            "all %d shape coordinates equal" % len(mp["shape_fields"]) in e["proof"]
            for e in mp["entries"]),
        "counterfactual_subdaily_on": {
            "n_suppressed": len(s2),
            "verdict": ("the suppression LIFTS BY ITSELF: with the sub-daily arm on "
                        "the F9-10 row reads the feature at event times, the shapes "
                        "differ, and the v1 row is no longer a duplicate of "
                        "anything. %d suppressed, as it must be." % len(s2)),
        },
    }


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
    """The §P7-16 reconciliation: 189 priced + 148 deferred + 9 unpriced + 23 removed.

    Three classes and a removal, each with a different accounting consequence, all
    of which used to be one undifferentiated "declared count":

      PRICED    enters the BH denominator and can be rejected.        -> m = 189
      DEFERRED  does not enter anything; returns as ONE B-2
                declaration when it is built.                         -> 148
      UNPRICED  declared with m_s = 0 so its rows have a stratum to
                route to; can neither be rejected nor consume budget. -> 9
      REMOVED   de-duplicated, so it is not a test at all.            -> 23
    """
    items = []
    for it in COMPOSITION:
        n = it["n"]
        if n is None:
            n = observer_control_count()
        items.append(dict(it, n=int(n)))
    priced = sum(i["n"] for i in items if i["class"] == "PRICED")
    deferred = sum(i["n"] for i in items if i["class"] == "DEFERRED")
    unpriced = sum(i["n"] for i in items if i["class"] == "UNPRICED")
    # §P7-16 de-duplicated `log_moment` against the CATALOG's 23-feature axis, and
    # 23 is the number its 369 reconciliation is built from -- so that is the figure
    # the reconciliation keeps. §P7-18 then lowered the axis to 20, on which the same
    # de-duplication removes 20; the other 3 left with their features to B-2. Both
    # are reported, because quoting either alone makes one of the two rulings not
    # add up.
    removed = len(marks_ext.DEDUPLICATED_MARKS) * MARK_AXIS_FEATURES_CATALOG
    removed_now = len(marks_ext.DEDUPLICATED_MARKS) * MARK_AXIS_FEATURES
    return {
        "items": items,
        "bh_denominator_m": priced,
        "n_priced": priced,
        "n_deferred_to_B2": deferred,
        "n_unpriced_controls": unpriced,
        "n_deduplicated_removed": removed,
        "n_deduplicated_removed_on_current_axis": removed_now,
        "n_deferred_features_B2": len(DEFERRED_FEATURES_B2),
        "deferred_features_B2": [dict(f) for f in DEFERRED_FEATURES_B2],
        "mark_axis_features_catalog": MARK_AXIS_FEATURES_CATALOG,
        "mark_axis_features_declared": MARK_AXIS_FEATURES,
        "superseded_denominator": POPPER_RULED_DENOMINATOR_SUPERSEDED,
        "change_reason_class": DEFERRAL_REASON_CLASS,
        "ratchet": RATCHET_NOTE,
        "ruled_denominator": POPPER_RULED_DENOMINATOR,
        "agrees_with_ruling": priced == POPPER_RULED_DENOMINATOR,
        "reconciliation": (
            "%d priced + %d deferred (B-2) + %d unpriced + %d de-duplicated = %d, "
            "against the ~1,000 headline §P7-10(c) carried from the catalog. "
            "RECONCILED EXACTLY at §P7-16; the residual is the catalog's own "
            "'declared overhead', which was never scope. §P7-18 then lowered the "
            "mark axis from 6 x 23 to 6 x 20 -- the 3 DOWNLOAD features B's own "
            "frozen config cannot build, deferred to B-2 BY NAME -- so m moved "
            "189 -> 171, %s."
            % (priced, deferred, unpriced, removed,
               priced + deferred + unpriced + removed, DEFERRAL_REASON_CLASS)),
        "deferred_note": DEFERRED_NOTE,
        "unpriced_rule": UNPRICED_RULE,
        "deduplication_rule": marks_ext.DEDUPLICATION_RULE,
        "built_check": (
            "every PRICED arm is implemented in this build: %s"
            % all(i["built"] for i in items if i["class"] == "PRICED")),
    }


# ------------------------------------------------------------ strata file -----
def strata_document(q=M.FDR_Q):
    """The declared §P6-3 partition for Tranche B, at the §P7-16 ruled integer.

    THREE THINGS THIS FILE DOES THAT THE FIRST DRAFT DID NOT.

    1. **It declares its own total.** `"m": 189` sits beside the strata table and
       `strata.assert_partition_total` refuses to load the file if the table does
       not add up to it. The first draft's headline and its itemisation disagreed by
       631 and nothing in the code could see it; now nothing can carry that.
    2. **It declares no strata for the DEFERRED arms.** A stratum with a declared
       `m_s` and no implementation consumes budget and produces nothing, taxing the
       arms that were built. The four unbuilt arms return as ONE B-2 declaration.
    3. **The F7 controls sit at `m_s = 0`.** They are DECLARED -- `stratum_of` must
       have somewhere to route a control row, and an undeclared row raises -- but
       they are outside the priced denominator, so they can neither be rejected
       (`_bh_within` returns no rejections at `m_s <= 0`) nor consume budget.

    §P6-3 rule 2: `q_s` is FLAT at the declared `q` in every priced stratum, so no
    budget is reallocated and there is no prior to justify. §P7-11(c) reallocated the
    MARK arm's `q_s` in Tranche A on measured grounds; B does not inherit that -- B's
    mark arm is a different arm (6 event-time marks, not 2 day-binned ones) and the
    measurement that justified the old allocation is not a measurement of this one.
    """
    enum = enumerate_declared()
    strata = []
    for it in enum["items"]:
        if it["class"] == "DEFERRED":
            continue
        unpriced = it["class"] == "UNPRICED"
        strata.append({
            "name": it["stratum"],
            "feature_family": (observer.OBSERVER_FAMILY if unpriced else None),
            "test_kind": it["test_kind"],
            "region": None,
            "m_s": (0 if unpriced else int(it["n"])),
            "q_s": float(q),
            "note": (UNPRICED_RULE if unpriced else ""),
        })
    return {
        "note": ("TRANCHE B (HYPOTHESIS_LEDGER.md §P7-3, §P7-10(c), §P7-14(d), "
                 "§P7-15(b)), RE-ISSUED at the §P7-18 ruled integer (189 -> 171, "
                 "change reason CONFIG-DETERMINED-NOT-RESULT-DETERMINED), declared "
                 "BEFORE the run and frozen in the config hash by its own sha256 "
                 "(§P6-3(3)+(5)). PRICED BH denominator m = %d = 17 second-moment "
                 "(F9-01) + 34 omnibus (F9-04) + 120 mark axis (F9-10, 6 scored "
                 "marks x 20 BUILDABLE features; the 3 DOWNLOAD features are "
                 "deferred to B-2 BY NAME). RECONCILIATION against the ~1,000 "
                 "headline: "
                 "%s DEFERRED: %s UNPRICED: %s DE-DUPLICATED: %s SUPPRESSED: "
                 "No budget is reallocated, so there is no prior to justify under "
                 "§P6-3(2). The mark axis sits in its OWN stratum (§P7-11(c)) "
                 "rather than in the v1 `mark` stratum, because 6 event-time marks "
                 "are not the 2 day-binned marks that allocation was measured on."
                 "%s RATCHET: %s DEFERRED FEATURES (B-2, by name): %s"
                 % (enum["bh_denominator_m"], enum["reconciliation"],
                    DEFERRED_NOTE, UNPRICED_RULE,
                    marks_ext.DEDUPLICATION_RULE,
                    disp.SUPPRESSION_RULE, RATCHET_NOTE,
                    ", ".join(f["feature"] for f in DEFERRED_FEATURES_B2))),
        "q": float(q),
        "m": int(enum["bh_denominator_m"]),
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
    part = strata_mod.load_partition(path, q=q)          # asserts sum_s m_s == m
    strata_mod.assert_budget_identity(part["strata"], q)  # asserts sum_s m_s q_s = m q
    return {"path": path.replace("\\", "/"),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "m": int(part["m"]), "m_declared": part["m_declared"],
            "declared_total_check": part["declared_total_check"],
            "q": float(q), "n_strata": len(part["strata"])}


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
    assert cfg["strata"]["m_declared"] == POPPER_RULED_DENOMINATOR, (
        "the partition's declared total is %s, the §P7-16 ruled integer is %s"
        % (cfg["strata"]["m_declared"], POPPER_RULED_DENOMINATOR))
    # §P7-17's third identity: the number of DECLARED rows must equal m. Asserted
    # HERE, on the declaration, because that is when the integer can still move.
    strata_mod.assert_declared_row_count(
        row_enumeration()["declared_view"]["totals"][disp.DECLARED],
        enum["bh_denominator_m"])
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
            "declared_marks": list(marks_ext.MARK_NAMES),
            "scored_marks": list(marks_ext.SCORED_MARK_NAMES),
            "deduplicated_marks": list(marks_ext.DEDUPLICATED_MARKS),
            "deduplication_rule": marks_ext.DEDUPLICATION_RULE,
            "n_scored_tests": len(marks_ext.SCORED_MARK_NAMES) * MARK_AXIS_FEATURES,
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
        "row_enumeration_P7_17": row_enumeration(),
        "ratchet_P7_18": {
            "denominator_now": POPPER_RULED_DENOMINATOR,
            "denominator_superseded": POPPER_RULED_DENOMINATOR_SUPERSEDED,
            "change_reason_class": DEFERRAL_REASON_CLASS,
            "rule": RATCHET_NOTE,
            "fixes_at": "config-hash freeze",
            "deferred_features_B2_by_name": [dict(f) for f in DEFERRED_FEATURES_B2],
        },
        "suppression_map_P7_18": row_enumeration()["suppression_map"],
        "open_items_before_the_run": [
            ("THE RUN GATE IS NOT DISCHARGED: the recovery-versus-amplitude curve "
             "is IN FLIGHT and the arm (i) anomaly is unruled (§P7-15(a), §P7-16)."),
        ],
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
    print("TRANCHE B -- DECLARATION RE-ISSUED AT THE §P7-16 RULED INTEGER, RUN GATED")
    print("=" * 78)
    print("| item | catalog | tests | class | stratum | m_s | built |")
    print("| --- | --- | ---: | --- | --- | ---: | --- |")
    for it in e["items"]:
        ms_ = ("-" if it["class"] == "DEFERRED"
               else "0" if it["class"] == "UNPRICED" else str(it["n"]))
        print("| %s | %s | %d | %s | %s | %s | %s |"
              % (it["key"], it["catalog"], it["n"], it["class"],
                 it["stratum"] or "(none -- B-2)", ms_,
                 "yes" if it["built"] else "NO"))
    print("")
    print("PRICED BH denominator m  : %d   (ruled %d -- %s)"
          % (e["bh_denominator_m"], e["ruled_denominator"],
             "AGREES" if e["agrees_with_ruling"] else "DISAGREES"))
    print("DEFERRED to one B-2 decl : %d" % e["n_deferred_to_B2"])
    print("UNPRICED (F7 controls)   : %d  (m_s = 0, outside the denominator)"
          % e["n_unpriced_controls"])
    print("DE-DUPLICATED, removed   : %d  (log_moment x %d features)"
          % (e["n_deduplicated_removed"], MARK_AXIS_FEATURES))
    print("reconciliation           : %s" % e["reconciliation"])
    print("")
    print("strata file              : %s (sha256 %s)"
          % (d["strata_file"]["path"], str(d["strata_file"]["sha256"])[:16]))
    print("  sum_s m_s == m         : %s" % (d["strata_file"]["declared_total_check"],))
    print("config hash              : %s" % d["config_hash"])
    print("count-path floor A_min   : %.4f at N = 46,585, alpha = %.1e"
          % (d["floors"]["count_path_A_min_at_46585"],
             d["floors"]["count_path_alpha"]))
    print("mark-axis floor rho_min  : %.4f at n = 46,585, VIF_mark = %.3f"
          % (d["floors"]["mark_axis_rho_min_at_46585"],
             d["floors"]["mark_axis_vif_fallback"]))
    print("S-15(c) unmeasurable     : %d entries (cut at period > %.0f d; "
          "period grid UNCHANGED at %.0f d, band reported not clamped)"
          % (d["s15c_sweep"]["n_unmeasurable_by_window"],
             d["s15c_sweep"]["cut_period_days"], ms.PERIOD_MAX))
    print("")
    _re = d["row_enumeration_P7_17"]
    print("")
    print("§P7-17 ROW ENUMERATION -- every row B's session would execute")
    print("| rows | test | disposition | parent / why |")
    print("| ---: | --- | --- | --- |")
    for r in _re["declared_view"]["rows"]:
        print("| %d | `%s` | %s | %s |"
              % (r["n"], r["test"], r["disposition"],
                 r["parent"] or r["detail"]))
    t = _re["declared_view"]["totals"]
    print("")
    print("DECLARED (= m)           : %d   (assertion count(DECLARED) == m: %s)"
          % (t[disp.DECLARED],
             "PASS" if t[disp.DECLARED] == e["bh_denominator_m"] else "FAIL"))
    print("COMPONENT-OF             : %d   (never standalone; attached to parents)"
          % t[disp.COMPONENT])
    print("REPLICATION-OF-DECLARED  : %d   (not re-priced; p-invariance wired)"
          % t[disp.REPLICATION])
    print("UNPRICED-CONTROL         : %d   (§P7-16, m_s = 0)"
          % t[disp.UNPRICED_CONTROL])
    print("TOTAL rows executed      : %d" % t["TOTAL"])
    print("GENUINELY NEW rows       : %d   -> the integer does NOT move"
          % _re["n_genuinely_new"])
    print("SUPPRESSED (mapped)      : %d   (v1 mark rows, proved same-shape)"
          % t[disp.SUPPRESSED])
    if _re["shortfall_flag"]:
        print("")
        print("SHORTFALL FLAG: " + _re["shortfall_flag"])
    else:
        print("shortfall                : 0 -- every priced slot is buildable "
              "under B's own frozen config")
    print("")
    print("§P7-18 RATCHET: m 189 -> 171, reason class %s" % DEFERRAL_REASON_CLASS)
    print("  deferred to B-2 BY NAME:")
    for f in DEFERRED_FEATURES_B2:
        print("    %-18s (family %d, %s)%s"
              % (f["feature"], f["family"], f["why_absent"],
                 "  [+ S-15(c) UNMEASURABLE-BY-WINDOW]" if f["s15c"] else ""))
    _sm = _re["suppression_map"]
    print("  suppression map: %d v1 mark rows suppressed, all proved same-shape "
          "on %s; active = %s; with the sub-daily arm ON it lifts (%d suppressed)"
          % (_sm["n_suppressed"], "/".join(_sm["shape_fields"]), _sm["active"],
             _sm["counterfactual_subdaily_on"]["n_suppressed"]))
    print("")
    print("RUN STATUS: %s" % d["run_status"])
    print(RUN_GATE_NOTE)
    with open(a.json, "w", encoding="utf-8") as fh:
        json.dump(d, fh, indent=1, default=float)
    print("\nwrote %s" % a.json)
    return d


if __name__ == "__main__":
    main()
