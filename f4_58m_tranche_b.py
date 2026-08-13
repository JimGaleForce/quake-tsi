"""F4-58M on the SIX SCORED MARKS of the tranche-B F9-10 axis (§P7-20(3)).

RULING
------
HYPOTHESIS_LEDGER.md §P7-20(3): "§P7-10(c) said F9-10 declares its floor from
measured VIF_mark or does not run. The run used §P7-13's pooled fallback 0.0469
on all 120 rows. The fallback was ruled for marks *without* their own
measurement; these 6 can be measured at zero surrogate cost. RULED. The 34
stand PROVISIONAL and are gated on F4-58M for the 6 marks."

Popper's recorded prediction, §P7-20(3), verbatim: "`dt_prior` and
`cluster_member` return `VIF_mark` above the 4.575 pooled fallback, and some of
the 23 clustering-mark survivors do not survive their own floors."

THE IDENTITY -- unchanged from f4_58m_vif_mark.py, re-run on tranche B's rows
-----------------------------------------------------------------------------
  Spearman           (df = 1), statistic rho: T_obs = (n - 1) * rho^2 ~ chi2_1
  circular-linear    (df = 2), statistic R:   T_obs = n * R^2         ~ chi2_2
  VIF_mark = T_obs / chi2.ppf(1 - p_block_bootstrap, df)
  rho_min  = sqrt(VIF_mark) * 4.732 / sqrt(n - 1)          (§P7-10(c))

PER-MARK POOLING, and why bounds are re-admitted
------------------------------------------------
The quantity being replaced is the POOLED FALLBACK 4.575 -- a population
summary over mark rows, built by §P7-13(b) with censored rows RE-ADMITTED AT
THEIR BOUNDS (the measurements-only median 4.345 is retained beside it). A
per-MARK VIF is the same kind of object: a population summary over that mark's
20 features, not a per-feature floor. §P7-13(b)'s own rule therefore permits
re-admission here, and comparing a measurements-only median against a
bounds-readmitted 4.575 would be comparing two different constructions. So:

  HEADLINE  = median over all 20 rows of the mark, censored rows at their
              bounds -- the like-for-like comparator to 4.575.
  BESIDE IT = median over uncensored rows only -- the like-for-like comparator
              to 4.345.

Both are reported for every mark, and every survivor verdict is evaluated under
BOTH, so no verdict rests on the choice.

CENSORING (§P7-13, §P7-10(a))
-----------------------------
p_block_bootstrap at its resolution floor 1/(B + 1) = 9.999e-5 (B = 10,000)
saturates chi2.ppf and returns a BOUND, not a measurement. Per the v2
correction in results_f4_58m_vif_mark.json (measured directly by the §P7-11(c)
un-censoring run) that bound is an UPPER bound: VIF decreases as the floor is
pushed down. Censored rows are flagged `censored: true` and their VIF is
reported as `vif_mark_upper_bound`.

PRICING -- ZERO (§P7-1(c), as invoked by §P7-10(c): "at zero new surrogate
cost"). Every input is read from checkpoint rows already declared and already
scored. No surrogate is drawn, no rejection is made, no BH vector is entered,
and NO EXPLORE_COUNT.jsonl line is written.

BOTH tranche-B sessions are analysed, because §P7-20(7) classified their
artifacts as a determinism failure (results_hash_classification.json). Running
this on both is how we find out whether that defect reaches this measurement.
Neither session directory is modified.

Usage:  python -u f4_58m_tranche_b.py
"""

from __future__ import annotations

import datetime as _dt
import json
import math
import os
import sys

import numpy as np
from scipy import stats
from scipy.optimize import brentq

REPO = os.path.dirname(os.path.abspath(__file__))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

OUT_PATH = os.path.join(REPO, "results_f4_58m_tranche_b.json")

SESSIONS = [
    ("session_20260813T092628", "tranche B, ledger-cited session (§P7-20)"),
    ("session_20260813T091153", "tranche B, first session (§P7-20(7) unattributable)"),
]
PRIMARY_SESSION = "session_20260813T092628"

MARK_TESTS = ("spearman", "circular-linear")
FLOOR_EPS = 1e-12

POOLED_FALLBACK_VIF = 4.575          # §P7-13(b), bounds-readmitted
POOLED_MEASONLY_VIF = 4.345          # §P7-13(b), measurements-only
POOLED_FALLBACK_RHO_MIN = 0.046894498589271714   # as carried on every row
LEDGER_CONST = 4.732                 # z_{1-alpha/2} + z_0.80 at alpha = 1e-4
TRANCHE_B_ALPHA = 1.0e-4
POWER = 0.80


def summarize(vals):
    v = [x for x in vals if x is not None and np.isfinite(x)]
    if not v:
        return None
    a = np.asarray(v, float)
    return {"n": int(a.size), "median": float(np.median(a)),
            "mean": float(a.mean()), "q25": float(np.percentile(a, 25)),
            "q75": float(np.percentile(a, 75)),
            "min": float(a.min()), "max": float(a.max())}


def chi2_equivalent(stat, df, n):
    return float((n - 1.0) * float(stat) ** 2) if df == 1 \
        else float(n * float(stat) ** 2)


def rho_min_ledger(vif, n, const=LEDGER_CONST):
    """§P7-10(c) as written and as APPLIED by the engine to all 120 rows."""
    return math.sqrt(vif) * const / math.sqrt(n - 1.0)


def rho_min_exact(vif, n, df, alpha=TRANCHE_B_ALPHA, power=POWER):
    """df-aware sensitivity: exact non-central chi-square power inversion.

    T_obs / VIF ~ ncx2(df, ncp) with ncp = k * stat^2 / VIF, k = n-1 (df 1) or
    n (df 2). Solve for the ncp giving `power` at the alpha critical value,
    then invert to the statistic. For df = 1 this reproduces the (z + z)/sqrt
    closed form to ~1e-3; for df = 2 it is the correct analogue, which the
    closed form is not.
    """
    crit = float(stats.chi2.isf(alpha, df))

    def f(ncp):
        return float(stats.ncx2.sf(crit, df, ncp)) - power

    ncp = brentq(f, 1e-9, 1e4, xtol=1e-10)
    k = (n - 1.0) if df == 1 else float(n)
    return math.sqrt(ncp * vif / k), ncp


def analyze(session, label):
    path = os.path.join(REPO, "engine", "out", "mine", session, "checkpoint.json")
    with open(path, "r", encoding="utf-8") as fh:
        ck = json.load(fh)

    rows = []
    for i, t in enumerate(ck.get("tests", [])):
        if t.get("test") not in MARK_TESTS or not t.get("mark"):
            continue
        n = t["n_events"]
        df = int(t["df"])
        stat = float(t["statistic"])
        p_boot = float(t["p_block_bootstrap"])
        B = int(t["n_surrogates"])
        floor_boot = 1.0 / (B + 1.0)
        censored = p_boot <= floor_boot + FLOOR_EPS
        T = chi2_equivalent(stat, df, n)
        ref = float(stats.chi2.ppf(1.0 - p_boot, df))
        vif = float(T / ref) if (np.isfinite(ref) and ref > 0) else None
        rows.append({
            "row_index": i,
            "feature": t["feature"], "mark": t["mark"], "test": t["test"],
            "df": df, "family": t.get("family"), "kind": t.get("kind"),
            "n_events": int(n), "n_surrogates": B,
            "statistic": stat,
            "chi2_equivalent": T,
            "p_block_bootstrap": p_boot,
            "p_circular_shift": t.get("p_circular_shift"),
            "p_raw": t.get("p_raw"), "p_bh": t.get("p_bh"),
            "p_method": t.get("p_method"),
            "resolution_floor_p_boot": floor_boot,
            "censored": bool(censored),
            "vif_mark": vif,
            "vif_mark_status": ("CENSORED_UPPER_BOUND" if censored else "MEASURED"),
            "rho_min_as_run_fallback": t.get("rho_min"),
            "statistic_over_floor_as_run": t.get("statistic_over_floor"),
            "bh_eligible": t.get("bh_eligible"),
            "passes_fdr": bool(t.get("passes_fdr")),
            "disposition": t.get("disposition"),
        })

    n_ev = rows[0]["n_events"]
    marks = sorted({r["mark"] for r in rows})

    per_mark = {}
    for m in marks:
        rs = [r for r in rows if r["mark"] == m]
        allv = [r["vif_mark"] for r in rs]
        unc = [r["vif_mark"] for r in rs if not r["censored"]]
        s_all, s_unc = summarize(allv), summarize(unc)
        vif_head = s_all["median"]
        vif_meas = s_unc["median"] if s_unc else None
        entry = {
            "mark": m,
            "n_rows": len(rs),
            "n_censored": sum(1 for r in rs if r["censored"]),
            "vif_mark_HEADLINE_bounds_readmitted_median": vif_head,
            "vif_mark_measurements_only_median": vif_meas,
            "vif_mark_summary_bounds_readmitted": s_all,
            "vif_mark_summary_measurements_only": s_unc,
            "vif_mark_df1_spearman_median": (
                summarize([r["vif_mark"] for r in rs if r["df"] == 1]) or {}).get("median"),
            "vif_mark_df2_circlin_median": (
                summarize([r["vif_mark"] for r in rs if r["df"] == 2]) or {}).get("median"),
            "exceeds_pooled_fallback_4_575": vif_head > POOLED_FALLBACK_VIF,
            "exceeds_pooled_fallback_measurements_only": (
                None if vif_meas is None else vif_meas > POOLED_FALLBACK_VIF),
            "rho_min_HEADLINE": rho_min_ledger(vif_head, n_ev),
            "rho_min_measurements_only": (
                None if vif_meas is None else rho_min_ledger(vif_meas, n_ev)),
            "rho_min_as_run_fallback": POOLED_FALLBACK_RHO_MIN,
            "rho_min_ratio_to_fallback": (
                rho_min_ledger(vif_head, n_ev) / POOLED_FALLBACK_RHO_MIN),
            "censored_rows": [
                {"feature": r["feature"], "test": r["test"], "df": r["df"],
                 "statistic": r["statistic"],
                 "p_block_bootstrap": r["p_block_bootstrap"],
                 "vif_mark_upper_bound": r["vif_mark"],
                 "censoring_note": (
                     "§P7-13 censored bound. p_block_bootstrap is at its Monte "
                     "Carlo resolution floor 1/(B+1); chi2.ppf saturates and "
                     "the VIF returned is an UPPER bound (direction measured by "
                     "the §P7-11(c) un-censoring run, results_f4_58m_uncensored"
                     ".json).")}
                for r in rs if r["censored"]],
        }
        # df-aware exact floors, per df, at the headline VIF
        entry["rho_min_exact_by_df"] = {}
        for df in (1, 2):
            rm, ncp = rho_min_exact(vif_head, n_ev, df)
            entry["rho_min_exact_by_df"][str(df)] = {
                "rho_min": rm, "ncp_at_80pct_power": ncp}
        per_mark[m] = entry

    # ---- the survivors, against their OWN mark's floor ---------------------
    survivors = []
    for r in rows:
        if not r["passes_fdr"]:
            continue
        pm = per_mark[r["mark"]]
        f_head = pm["rho_min_HEADLINE"]
        f_meas = pm["rho_min_measurements_only"]
        f_exact = pm["rho_min_exact_by_df"][str(r["df"])]["rho_min"]
        survivors.append({
            "row_index": r["row_index"],
            "feature": r["feature"], "mark": r["mark"], "test": r["test"],
            "df": r["df"], "family": r["family"],
            "statistic": r["statistic"],
            "p_bh": r["p_bh"], "p_raw": r["p_raw"],
            "vif_mark_own_row": r["vif_mark"],
            "vif_mark_of_its_mark_HEADLINE": pm["vif_mark_HEADLINE_bounds_readmitted_median"],
            "rho_min_as_run_fallback": POOLED_FALLBACK_RHO_MIN,
            "stat_over_floor_as_run": r["statistic_over_floor_as_run"],
            "rho_min_own_mark_HEADLINE": f_head,
            "stat_over_own_floor_HEADLINE": r["statistic"] / f_head,
            "SURVIVES_OWN_FLOOR_HEADLINE": bool(r["statistic"] >= f_head),
            "rho_min_own_mark_measurements_only": f_meas,
            "stat_over_own_floor_measurements_only": (
                None if f_meas is None else r["statistic"] / f_meas),
            "SURVIVES_OWN_FLOOR_measurements_only": (
                None if f_meas is None else bool(r["statistic"] >= f_meas)),
            "rho_min_own_mark_df_aware_exact": f_exact,
            "SURVIVES_OWN_FLOOR_df_aware_exact": bool(r["statistic"] >= f_exact),
            "censored_row": r["censored"],
        })
    survivors.sort(key=lambda s: -s["stat_over_own_floor_HEADLINE"])

    n_surv = len(survivors)
    n_pass_head = sum(1 for s in survivors if s["SURVIVES_OWN_FLOOR_HEADLINE"])
    n_pass_meas = sum(1 for s in survivors
                      if s["SURVIVES_OWN_FLOOR_measurements_only"])
    n_pass_exact = sum(1 for s in survivors if s["SURVIVES_OWN_FLOOR_df_aware_exact"])
    n_pass_asrun = sum(1 for s in survivors
                       if s["statistic"] >= POOLED_FALLBACK_RHO_MIN)

    clustering = ("dt_prior_days", "dist_nearest_prior_km", "cluster_member")
    clus = [s for s in survivors if s["mark"] in clustering]

    return {
        "session": session, "label": label,
        "config_hash": ck.get("config_hash"),
        "artifact_hash": (ck.get("build_invariant") or {}).get("artifact_hash"),
        "n_events": n_ev,
        "n_surrogates_block_bootstrap": rows[0]["n_surrogates"],
        "resolution_floor_p_boot": rows[0]["resolution_floor_p_boot"],
        "n_mark_rows": len(rows),
        "marks_scored": marks,
        "n_censored_total": sum(1 for r in rows if r["censored"]),
        "per_mark": per_mark,
        "survivor_verdict": {
            "n_survivors": n_surv,
            "n_surviving_own_floor_HEADLINE": n_pass_head,
            "n_dying_on_own_floor_HEADLINE": n_surv - n_pass_head,
            "n_surviving_own_floor_measurements_only": n_pass_meas,
            "n_surviving_own_floor_df_aware_exact": n_pass_exact,
            "n_surviving_as_run_pooled_fallback": n_pass_asrun,
            "n_clustering_mark_survivors": len(clus),
            "n_clustering_mark_survivors_dying_HEADLINE": sum(
                1 for s in clus if not s["SURVIVES_OWN_FLOOR_HEADLINE"]),
        },
        "survivors": survivors,
        "rows": rows,
    }


def main():
    out = {
        "id": "F4-58M-TRANCHE-B",
        "title": ("Mark-axis variance inflation factor on the six scored marks "
                  "of the tranche-B F9-10 axis, and the 34 survivors against "
                  "their own measured floors"),
        "ruling": ("HYPOTHESIS_LEDGER.md §P7-20(3), gating the 34 PROVISIONAL "
                   "survivors; identity §P7-10(c); censoring §P7-13 / §P7-10(a); "
                   "pricing-at-zero §P7-1(c)"),
        "generated_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 0,
        "priced_tests_note": (
            "Zero, per §P7-1(c) as invoked by §P7-10(c) ('at zero new surrogate "
            "cost'). No surrogate drawn, no rejection made, no BH vector "
            "entered, NO EXPLORE_COUNT.jsonl line written. Every input is a "
            "(statistic, p_block_bootstrap, df, n) tuple already recorded by "
            "the session that computed it."),
        "state_class": (
            "MEASUREMENT on existing artifacts (estimator characterisation). "
            "Not a hypothesis test about the Earth. It RESOLVES a declared "
            "gate; it creates no new rejection."),
        "identity": (
            "VIF_mark = T_obs / chi2.ppf(1 - p_block_bootstrap, df); "
            "T_obs = (n-1)*rho^2 for df=1 Spearman, T_obs = n*R^2 for df=2 "
            "circular-linear; rho_min = sqrt(VIF_mark)*4.732/sqrt(n-1)."),
        "pooling_rule": (
            "Per-mark VIF is the MEDIAN over that mark's 20 rows. HEADLINE "
            "re-admits censored rows at their bounds, because the quantity it "
            "replaces -- the pooled fallback 4.575 -- was itself built that "
            "way by §P7-13(b), and because a per-MARK value is a population "
            "summary over 20 features, not a per-feature floor. The "
            "measurements-only median (the 4.345 analogue) is reported beside "
            "it and every verdict is evaluated under both."),
        "GATE_INHERITED": {
            "gate": "§P7-20(7) artifact-hash classification",
            "status": "NOT LIFTED -- CASE_B_INVARIANCE_FAILURE",
            "consequence": (
                "This measurement reads the same tranche-B checkpoints that "
                "failed the §P7-20(7) invariance audit "
                "(results_hash_classification.json). Nothing below may be "
                "quoted until that gate lifts. It is computed on BOTH sessions "
                "precisely so the reader can see whether the determinism "
                "defect reaches this measurement -- see "
                "`cross_session_invariance`."),
        },
        "sessions": [],
    }

    for s, lab in SESSIONS:
        out["sessions"].append(analyze(s, lab))

    prim = next(s for s in out["sessions"] if s["session"] == PRIMARY_SESSION)
    alt = next(s for s in out["sessions"] if s["session"] != PRIMARY_SESSION)

    # ---- cross-session invariance of THIS measurement ----------------------
    vif_delta = {}
    for m in prim["marks_scored"]:
        a = prim["per_mark"][m]["vif_mark_HEADLINE_bounds_readmitted_median"]
        b = alt["per_mark"][m]["vif_mark_HEADLINE_bounds_readmitted_median"]
        vif_delta[m] = {"primary": a, "other": b,
                        "relative_difference": abs(a - b) / max(abs(a), abs(b))}
    pv = {(s["feature"], s["mark"], s["test"]): s["SURVIVES_OWN_FLOOR_HEADLINE"]
          for s in prim["survivors"]}
    av = {(s["feature"], s["mark"], s["test"]): s["SURVIVES_OWN_FLOOR_HEADLINE"]
          for s in alt["survivors"]}
    out["cross_session_invariance"] = {
        "same_survivor_set": set(pv) == set(av),
        "n_survivors_primary": len(pv), "n_survivors_other": len(av),
        "per_mark_vif_relative_difference": vif_delta,
        "max_vif_relative_difference": max(
            v["relative_difference"] for v in vif_delta.values()),
        "verdict_disagreements": [
            {"key": list(k), "primary": pv[k], "other": av.get(k)}
            for k in pv if pv[k] != av.get(k)],
        "reading": (
            "If the per-mark VIFs agree to float noise and no survivor verdict "
            "flips, the §P7-20(7) determinism defect does NOT reach this "
            "measurement -- which constrains the defect's blast radius but "
            "does not lift its gate."),
    }

    # ---- the recorded prediction, adjudicated ------------------------------
    pm = prim["per_mark"]
    pred_dt = pm["dt_prior_days"]["vif_mark_HEADLINE_bounds_readmitted_median"]
    pred_cm = pm["cluster_member"]["vif_mark_HEADLINE_bounds_readmitted_median"]
    n_died = prim["survivor_verdict"]["n_dying_on_own_floor_HEADLINE"]
    n_clus_died = prim["survivor_verdict"]["n_clustering_mark_survivors_dying_HEADLINE"]

    limb_a = bool(pred_dt > POOLED_FALLBACK_VIF and pred_cm > POOLED_FALLBACK_VIF)
    limb_b = bool(n_died > 0)

    out["POPPER_PREDICTION_P7_20_3"] = {
        "recorded_before_measurement": (
            "dt_prior and cluster_member return VIF_mark above the 4.575 "
            "pooled fallback, and some of the 23 clustering-mark survivors do "
            "not survive their own floors."),
        "limb_a_clustering_marks_exceed_4_575": {
            "dt_prior_days_vif_mark": pred_dt,
            "cluster_member_vif_mark": pred_cm,
            "threshold": POOLED_FALLBACK_VIF,
            "VERDICT": "CONFIRMED" if limb_a else "REFUTED",
            "dt_prior_days_measurements_only": pm["dt_prior_days"][
                "vif_mark_measurements_only_median"],
            "cluster_member_measurements_only": pm["cluster_member"][
                "vif_mark_measurements_only_median"],
        },
        "limb_b_some_survivors_die": {
            "n_survivors": prim["survivor_verdict"]["n_survivors"],
            "n_dying_on_own_floor": n_died,
            "n_clustering_mark_survivors_dying": n_clus_died,
            "VERDICT": "CONFIRMED" if limb_b else "REFUTED",
        },
        "OVERALL": ("CONFIRMED" if (limb_a and limb_b)
                    else ("PARTIALLY CONFIRMED" if (limb_a or limb_b)
                          else "REFUTED")),
    }

    # ---- the §P7-20(1)/(2) classes: 3 VOID, 8 uncontrolled, 23 exogenous ----
    SAME_QUANTITY_VOID = {("mean_depth_30d", "depth"),
                          ("deep_fraction_30d", "depth"),
                          ("b_value_90d", "mag")}

    def klass(s):
        if (s["feature"], s["mark"]) in SAME_QUANTITY_VOID:
            return "VOID_same_quantity_P7_20_1"
        if s["family"] == 4:
            return "ARTIFACT_CLASS_UNCONTROLLED_P7_20_1"
        return "EXOGENOUS_ephemeris_P7_20_2"

    classes = {}
    for s in prim["survivors"]:
        c = klass(s)
        e = classes.setdefault(c, {"n": 0, "n_survive_own_floor": 0,
                                   "n_die_own_floor": 0,
                                   "surviving": [], "dying": []})
        e["n"] += 1
        tag = "%s x %s" % (s["feature"], s["mark"])
        if s["SURVIVES_OWN_FLOOR_HEADLINE"]:
            e["n_survive_own_floor"] += 1
            e["surviving"].append(tag)
        else:
            e["n_die_own_floor"] += 1
            e["dying"].append(tag)
    out["survivor_classes_P7_20"] = {
        "classification_source": (
            "§P7-20(1) named the 3 same-quantity rows VOID and the 8 remaining "
            "family-4 catalogue-endogenous rows ARTIFACT-CLASS-UNCONTROLLED; "
            "§P7-20(2) named the residual 23 EXOGENOUS ephemeris survivors. "
            "Reconstructed here from each row's own (feature, mark, family)."),
        "classes": classes,
        "reading": (
            "The floor measurement does not fall evenly across the classes. "
            "That asymmetry is the substance of the result and is stated "
            "explicitly in `headline.class_asymmetry`."),
    }

    exo = classes.get("EXOGENOUS_ephemeris_P7_20_2", {"n": 0, "n_die_own_floor": 0,
                                                      "n_survive_own_floor": 0})

    out["headline"] = {
        "class_asymmetry": {
            "exogenous_n": exo["n"],
            "exogenous_dying_on_own_floor": exo["n_die_own_floor"],
            "exogenous_surviving_own_floor": exo["n_survive_own_floor"],
            "catalogue_derived_n": sum(
                v["n"] for k, v in classes.items()
                if k != "EXOGENOUS_ephemeris_P7_20_2"),
            "catalogue_derived_dying_on_own_floor": sum(
                v["n_die_own_floor"] for k, v in classes.items()
                if k != "EXOGENOUS_ephemeris_P7_20_2"),
            "statement": (
                "Every row that dies on its own measured floor is EXOGENOUS "
                "(ephemeris). Every VOID and every ARTIFACT-CLASS-UNCONTROLLED "
                "catalogue-derived row clears its own floor comfortably. This "
                "is the expected direction and not a reassuring one: the rows "
                "that survive the floor measurement are exactly the rows "
                "§P7-20(1) already struck or quarantined for being about the "
                "catalogue rather than about the Earth."),
        },
        "n_mark_rows": prim["n_mark_rows"],
        "marks_scored": prim["marks_scored"],
        "n_censored": prim["n_censored_total"],
        "vif_mark_by_mark_HEADLINE": {
            m: pm[m]["vif_mark_HEADLINE_bounds_readmitted_median"]
            for m in prim["marks_scored"]},
        "rho_min_by_mark_HEADLINE": {
            m: pm[m]["rho_min_HEADLINE"] for m in prim["marks_scored"]},
        "rho_min_as_run_all_120_rows": POOLED_FALLBACK_RHO_MIN,
        "survivor_verdict": prim["survivor_verdict"],
    }

    out["caveats"] = [
        "GATE. Every number here reads tranche-B checkpoints that failed the "
        "§P7-20(7) invariance audit. Nothing here may be quoted until that "
        "gate lifts, whatever the cross-session agreement shows.",
        "SCOPE, carried unchanged from F4-58M v2: this measures the inflation "
        "of the mark null AS CONSTRUCTED (circular moving-block bootstrap, "
        "block_events = 500 over 46,585 events). It does NOT decompose that "
        "inflation into real mark-axis dependence in the Earth versus width "
        "manufactured by a ~93-block bootstrap. That decomposition remains "
        "UNMEASURED, per §P7-10(c).",
        "The floor formula rho_min = sqrt(VIF)*4.732/sqrt(n-1) is a df=1 "
        "closed form, and the engine applied it to df=2 circular-linear rows "
        "too (visible on every mark row: a df=2 row carries the same 0.046894 "
        "as a df=1 row). The HEADLINE verdict keeps that convention so it is "
        "like-for-like with the run being adjudicated. "
        "`rho_min_exact_by_df` / `SURVIVES_OWN_FLOOR_df_aware_exact` carry the "
        "correct non-central chi-square inversion per df as a sensitivity, and "
        "the survivor counts under it are reported.",
        "The chi-square-scale conversion assumes the standard large-sample "
        "nulls ((n-1)rho^2 ~ chi2_1; nR^2 ~ chi2_2). At n = 46,585 the "
        "asymptotics are not in question, but the conversion is an assumption "
        "of the READING, not of the recorded statistic.",
        "Censored rows are at the p_boot resolution floor 9.999e-5 (B = "
        "10,000, one fifth the B of the v1 sessions), so tranche B censors "
        "more readily than the v1 axis did. Their VIF is an UPPER bound; "
        "re-admitting them raises the headline VIF and therefore RAISES the "
        "floor, which is the conservative direction for a survivor claim.",
        "This resolves §P7-20(3) only. The 34 remain gated on the stability "
        "line (§P7-20(4)), the 23 exogenous on the phase-randomised null "
        "(§P7-20(2)), the 3 same-quantity rows are VOID and the 8 "
        "catalogue-endogenous ARTIFACT-CLASS-UNCONTROLLED (§P7-20(1)). A row "
        "clearing its own floor here has cleared ONE gate.",
        "Neither session directory was modified; both were opened read-only.",
    ]

    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print("wrote %s" % OUT_PATH)

    n_ev = prim["n_events"]
    print("\n=== F4-58M ON THE SIX SCORED MARKS (%s, n = %d) ==="
          % (PRIMARY_SESSION, n_ev))
    print("%-24s %5s %5s %10s %10s %10s %10s %8s"
          % ("mark", "rows", "cens", "VIF(bnd)", "VIF(meas)", "rho_min",
             "as-run", "x fallb"))
    for m in prim["marks_scored"]:
        e = pm[m]
        print("%-24s %5d %5d %10.3f %10.3f %10.5f %10.5f %8.2f%s"
              % (m, e["n_rows"], e["n_censored"],
                 e["vif_mark_HEADLINE_bounds_readmitted_median"],
                 e["vif_mark_measurements_only_median"],
                 e["rho_min_HEADLINE"], POOLED_FALLBACK_RHO_MIN,
                 e["rho_min_ratio_to_fallback"],
                 "  >4.575" if e["exceeds_pooled_fallback_4_575"] else ""))

    print("\n=== THE 34 SURVIVORS AGAINST THEIR OWN MEASURED FLOORS ===")
    print("%-26s %-22s %-16s %3s %9s %9s %8s %8s  %s"
          % ("feature", "mark", "test", "df", "stat", "rho_min", "as-run",
             "own", "VERDICT"))
    for s in prim["survivors"]:
        print("%-26s %-22s %-16s %3d %9.5f %9.5f %8.2f %8.2f  %s"
              % (s["feature"], s["mark"], s["test"], s["df"], s["statistic"],
                 s["rho_min_own_mark_HEADLINE"], s["stat_over_floor_as_run"],
                 s["stat_over_own_floor_HEADLINE"],
                 "SURVIVES" if s["SURVIVES_OWN_FLOOR_HEADLINE"] else "*** DIES ***"))

    sv = prim["survivor_verdict"]
    print("\n  survivors %d | survive own floor (headline) %d | DIE %d"
          % (sv["n_survivors"], sv["n_surviving_own_floor_HEADLINE"],
             sv["n_dying_on_own_floor_HEADLINE"]))
    print("  survive under measurements-only VIF: %d | under df-aware exact: %d"
          % (sv["n_surviving_own_floor_measurements_only"],
             sv["n_surviving_own_floor_df_aware_exact"]))
    print("  clustering-mark survivors %d, of which dying %d"
          % (sv["n_clustering_mark_survivors"],
             sv["n_clustering_mark_survivors_dying_HEADLINE"]))

    p = out["POPPER_PREDICTION_P7_20_3"]
    print("\n=== POPPER'S RECORDED PREDICTION (§P7-20(3)) ===")
    print("  limb (a) dt_prior_days VIF = %.3f, cluster_member VIF = %.3f "
          "vs 4.575 -> %s"
          % (p["limb_a_clustering_marks_exceed_4_575"]["dt_prior_days_vif_mark"],
             p["limb_a_clustering_marks_exceed_4_575"]["cluster_member_vif_mark"],
             p["limb_a_clustering_marks_exceed_4_575"]["VERDICT"]))
    print("  limb (b) %d of %d survivors die on their own floors -> %s"
          % (p["limb_b_some_survivors_die"]["n_dying_on_own_floor"],
             p["limb_b_some_survivors_die"]["n_survivors"],
             p["limb_b_some_survivors_die"]["VERDICT"]))
    print("  OVERALL: %s" % p["OVERALL"])

    print("\n=== BY §P7-20 CLASS ===")
    for c, e in sorted(classes.items()):
        print("  %-38s n=%2d  survive own floor %2d  DIE %2d"
              % (c, e["n"], e["n_survive_own_floor"], e["n_die_own_floor"]))

    ci = out["cross_session_invariance"]
    print("\n=== CROSS-SESSION INVARIANCE (§P7-20(7) blast radius) ===")
    print("  same survivor set: %s | max per-mark VIF rel diff: %.3e | "
          "verdict disagreements: %d"
          % (ci["same_survivor_set"], ci["max_vif_relative_difference"],
             len(ci["verdict_disagreements"])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
