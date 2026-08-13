"""F4-58M -- the MARK-AXIS variance inflation factor, from existing session output.

RULING
------
HYPOTHESIS_LEDGER.md §P7-10(c): "Precondition, priced at 0: F4-58M, the mark-axis
VIF. The same identity applied to the mark tests already on disk --
`VIF_mark = stat_obs / null_ppf(1 - p_block_bootstrap)` -- at zero new surrogate
cost, exactly as F4-58 was. It is needed because `block_events = 500` over
~46,585 events is ~93 blocks and the mark null is inflated by its own
construction, unmeasured."

PRICING -- ZERO, on §P7-1(c)'s own logic (cited, not re-argued): this is a
re-reading of statistics already computed and already declared by the sessions
below. It draws no surrogate, makes no rejection, and enters no BH vector. NO
EXPLORE_COUNT.jsonl line is written by this script.

THE IDENTITY, MADE EXPLICIT ON THE MARK AXIS
--------------------------------------------
The count-path identity is VIF = chi2_obs / chi2.ppf(1 - p, df). To apply the
SAME identity to a correlation statistic, the statistic is first put on the
chi-square scale, which is where a variance inflation is multiplicative:

  * Spearman (df = 1), statistic |rho|:  T_obs = (n - 1) * rho^2   ~ chi2_1
  * circular-linear (df = 2), statistic R: T_obs = n * R^2         ~ chi2_2

Both are the standard large-sample null distributions (Var(rho) = 1/(n-1) under
independence; Mardia's circular-linear R with n*R^2 -> chi2_2). At n = 46,585 the
difference between n and n - 1 is 2e-5 relative and is immaterial. Then

    VIF_mark = T_obs / chi2.ppf(1 - p_block_bootstrap, df)

which is EXACTLY the §P7-1(c) identity, and it is the quantity that makes the
§P7-10(c) floor formula exact rather than analogical, because

    Var_inflated(rho) = VIF_mark / (n - 1)
      =>  rho_min = sqrt(VIF_mark) * (z_alpha + z_0.80) / sqrt(n - 1).

CENSORING (§P7-10(a)) APPLIES HERE AND, UNLIKE ON THE COUNT PATH, IT BITES
--------------------------------------------------------------------------
Rows whose `p_block_bootstrap` sits at its Monte Carlo resolution floor
1/(B + 1) are EXCLUDED from every median and reported separately with their
count. On the mark axis these are the three most significant tests in each
session, which is precisely the anti-conservative case §P7-10(a) names.

Usage:  python -u f4_58m_vif_mark.py
"""

from __future__ import annotations

import datetime as _dt
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np
from scipy import stats

REPO = os.path.dirname(os.path.abspath(__file__))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from f4_58_vif import FLOOR_EPS, SESSIONS, summarize  # noqa: E402

OUT_PATH = os.path.join(REPO, "results_f4_58m_vif_mark.json")
PRIMARY_SESSION = "session_20260811T022953"

MARK_TESTS = ("spearman", "circular-linear")

# §P7-10(c): Tranche B alpha and the 80%-power constant.
TRANCHE_B_ALPHA = 1.0e-4
Z80 = float(stats.norm.ppf(0.80))                      # 0.8416
Z_ALPHA_B = float(stats.norm.isf(TRANCHE_B_ALPHA))     # 3.7190 one-sided
# The ledger's constant 4.732 is z + z_0.80 with z = 3.891, i.e. the TWO-SIDED
# z at alpha = 1e-4 (Phi^-1(1 - alpha/2) = 3.8906). Both are carried; the
# ledger's own value is the operative one.
Z_ALPHA_B_TWOSIDED = float(stats.norm.isf(TRANCHE_B_ALPHA / 2.0))
LEDGER_CONST = 4.732


def chi2_equivalent(stat, df, n):
    """Correlation statistic -> chi-square-scale statistic under the null."""
    if df == 1:
        return float((n - 1.0) * float(stat) ** 2)
    return float(n * float(stat) ** 2)


def analyze(session, label):
    path = os.path.join(REPO, "engine", "out", "mine", session, "checkpoint.json")
    with open(path, "r", encoding="utf-8") as fh:
        ck = json.load(fh)
    B = ck.get("config", {}).get("n_surrogates")
    floor_boot = (1.0 / (B + 1.0)) if B else None

    rows = []
    status = defaultdict(int)
    for t in ck.get("tests", []):
        if t.get("test") not in MARK_TESTS:
            continue
        n = t.get("n_events")
        stat = t.get("statistic")
        df = t.get("df")
        p_boot = t.get("p_block_bootstrap")
        if n is None or stat is None or df is None or p_boot is None:
            status["MISSING_INPUT"] += 1
            continue
        T = chi2_equivalent(stat, df, n)
        censored = bool(floor_boot is not None and p_boot <= floor_boot + FLOOR_EPS)
        if p_boot >= 1.0:
            vif, st = None, "P_ONE_UNRESOLVED"
        else:
            ref = float(stats.chi2.ppf(1.0 - p_boot, df))
            if not np.isfinite(ref) or ref <= 0:
                vif, st = None, "BAD_REF"
            else:
                vif, st = float(T / ref), "OK"
        p_raw = t.get("p_raw")
        vif_raw = None
        if p_raw is not None and 0.0 < p_raw < 1.0:
            r2 = float(stats.chi2.ppf(1.0 - p_raw, df))
            if np.isfinite(r2) and r2 > 0:
                vif_raw = float(T / r2)
        status[st] += 1
        rows.append({
            "feature": t.get("feature"),
            "family": t.get("family"),
            "kind": t.get("kind"),
            "mark": t.get("mark"),
            "test": t.get("test"),
            "df": int(df),
            "n_events": int(n),
            "statistic": float(stat),
            "effect": t.get("effect"),
            "chi2_equivalent": T,
            "p_block_bootstrap": float(p_boot),
            "p_circular_shift": t.get("p_circular_shift"),
            "p_raw": p_raw,
            "vif_mark": vif,
            "vif_mark_from_p_raw": vif_raw,
            "status": st,
            "censored": censored,
        })

    ok = [r for r in rows if r["status"] == "OK"]
    keep = [r for r in ok if not r["censored"]]
    cens = [r for r in ok if r["censored"]]

    def by(pred):
        return summarize([r["vif_mark"] for r in keep if pred(r)])

    per_feature = {}
    for r in keep:
        per_feature.setdefault((r["feature"], r["mark"]), []).append(r)
    feature_table = []
    for (feat, mark), rs in sorted(per_feature.items()):
        s = summarize([r["vif_mark"] for r in rs])
        feature_table.append({
            "feature": feat, "mark": mark, "test": rs[0]["test"],
            "df": rs[0]["df"], "n_events": rs[0]["n_events"],
            "statistic": rs[0]["statistic"],
            "p_block_bootstrap": rs[0]["p_block_bootstrap"],
            "vif_mark": None if s is None else s["median"],
        })
    feature_table.sort(key=lambda r: (-(r["vif_mark"] or 0.0)))

    return {
        "session": session,
        "label": label,
        "config_hash": ck.get("config_hash"),
        "grid": ck.get("config", {}).get("grid"),
        "n_surrogates_block_bootstrap": B,
        "resolution_floor_p_boot": floor_boot,
        "block_events": 500,
        "block_events_source": "engine/mine.py:mark_test default block_events=500",
        "n_mark_rows": len(rows),
        "status_counts": dict(status),
        "censoring_P7_10a": {
            "n_scored": len(ok),
            "n_censored_excluded": len(cens),
            "censored_rows": [
                {k: r[k] for k in ("feature", "mark", "test", "df", "statistic",
                                   "p_block_bootstrap", "p_raw", "chi2_equivalent",
                                   "vif_mark", "vif_mark_from_p_raw")}
                for r in cens
            ],
            "censored_vif_summary": summarize([r["vif_mark"] for r in cens]),
            "direction": (
                "EXCLUDED because chi2.ppf saturates at the floor and the VIF "
                "returned there is a BOUND, not a measurement. CORRECTED "
                "2026-08-12 by the §P7-11(c) un-censoring run: it is an UPPER "
                "bound, not a lower one. VIF = T / chi2.ppf(1 - p, df) DECREASES "
                "as the floor is pushed down, measured at B = 2,000 / 50,000 / "
                "500,000 on all three rows (b_value_90d x mag: 25.87 -> 17.23 -> "
                "13.87). See results_f4_58m_uncensored.json :: "
                "direction_correction_to_P7_10a."),
            "direction_RETRACTED_2026_08_12": (
                "Earlier text in this field read 'the VIF returned there is a "
                "LOWER BOUND ... including them would understate the median and so "
                "understate the floor (anti-conservative)'. That is inverted and "
                "is retracted. No number in this file changes: the medians are "
                "computed over the uncensored rows under either reading."),
        },
        "vif_mark_all": by(lambda r: True),
        "vif_mark_df1_spearman": by(lambda r: r["df"] == 1),
        "vif_mark_df2_circular_linear": by(lambda r: r["df"] == 2),
        "vif_mark_by_mark": {
            m: by(lambda r, m=m: r["mark"] == m)
            for m in sorted({r["mark"] for r in keep if r["mark"]})
        },
        "vif_mark_null_only_p_gt_0.10": by(lambda r: (r["p_block_bootstrap"] or 0) > 0.10),
        "per_feature": feature_table,
        "rows": rows,
    }


def rho_min(vif, n, const=LEDGER_CONST):
    return math.sqrt(vif) * const / math.sqrt(n - 1.0)


def main():
    out = {
        "id": "F4-58M",
        "title": "Mark-axis variance inflation factor (VIF_mark) from existing session output",
        "ruling": ("HYPOTHESIS_LEDGER.md §P7-10(c) -- zero-priced precondition for "
                   "Tranche B; censoring rule §P7-10(a); pricing-at-zero §P7-1(c)"),
        "version": 2,
        "version_note": (
            "v2, 2026-08-12: no measured value changed. The DIRECTION of the "
            "§P7-10(a) censoring bias is corrected from 'understated / lower "
            "bound' to 'over-stated / upper bound', on the direct evidence of the "
            "§P7-11(c) un-censoring run (results_f4_58m_uncensored.json). v1's "
            "medians, floors and per-feature table are unchanged because they were "
            "always computed over the UNCENSORED rows."),
        "priced_tests": 0,
        "priced_tests_note": (
            "Zero, per §P7-1(c) as invoked by §P7-10(c) ('at zero new surrogate cost, "
            "exactly as F4-58 was'). No surrogate is drawn, no rejection is made, no "
            "BH vector is entered, and NO EXPLORE_COUNT.jsonl line is written. Every "
            "input is read from checkpoint.json rows already declared by the sessions "
            "that computed them."),
        "generated_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "identity": ("VIF_mark = T_obs / chi2.ppf(1 - p_block_bootstrap, df), with "
                     "T_obs = (n-1)*rho^2 for df=1 Spearman and T_obs = n*R^2 for "
                     "df=2 circular-linear -- the §P7-1(c) identity with the "
                     "correlation statistic put on its own chi-square scale."),
        "state_class": ("MEASUREMENT on existing artifacts (estimator characterisation, "
                        "not a hypothesis test about the Earth)."),
        "sessions": [],
    }
    for s, lab in SESSIONS:
        try:
            out["sessions"].append(analyze(s, lab))
        except FileNotFoundError:
            print("missing session: %s" % s)

    prim = next(s for s in out["sessions"] if s["session"] == PRIMARY_SESSION)
    n_ev = prim["rows"][0]["n_events"]

    med_all = prim["vif_mark_all"]["median"]
    med_df1 = prim["vif_mark_df1_spearman"]["median"]
    med_df2 = prim["vif_mark_df2_circular_linear"]["median"]

    # ---- the operative mark floor, §P7-10(c) --------------------------------
    ladder = {}
    for lab, v in (("VIF_mark = 1 (no inflation)", 1.0),
                   ("VIF_mark = 4", 4.0), ("VIF_mark = 9", 9.0),
                   ("VIF_mark = 24 (count-path scalar, for comparison only)", 24.0),
                   ("VIF_mark MEASURED, all mark tests (median)", med_all),
                   ("VIF_mark MEASURED, df=1 Spearman (median)", med_df1),
                   ("VIF_mark MEASURED, df=2 circular-linear (median)", med_df2),
                   ("VIF_mark MEASURED, 75th pct all mark tests (conservative)",
                    prim["vif_mark_all"]["q75"])):
        ladder[lab] = {"vif_mark": v, "rho_min": rho_min(v, n_ev)}

    per_feature_floor = [
        {"feature": f["feature"], "mark": f["mark"], "df": f["df"],
         "vif_mark": f["vif_mark"], "rho_min": rho_min(f["vif_mark"], f["n_events"]),
         "observed_statistic": f["statistic"],
         "observed_over_floor": (f["statistic"] / rho_min(f["vif_mark"], f["n_events"]))}
        for f in prim["per_feature"] if f["vif_mark"]
    ]

    out["operative_mark_floor_P7_10c"] = {
        "formula": "rho_min = sqrt(VIF_mark) * (z_alpha + z_0.80) / sqrt(n - 1)",
        "ledger_constant": LEDGER_CONST,
        "ledger_constant_decomposition": (
            "z_alpha + z_0.80 = %.4f + %.4f = %.4f, z_alpha two-sided at "
            "alpha = 1.0e-4" % (Z_ALPHA_B_TWOSIDED, Z80, Z_ALPHA_B_TWOSIDED + Z80)),
        "one_sided_alternative_constant": Z_ALPHA_B + Z80,
        "alpha": TRANCHE_B_ALPHA,
        "n_events": n_ev,
        "n_events_source": ("checkpoint.json mark rows carry n_events = 46585 "
                            "directly; cross-checks against n_surrogates = 46386 "
                            "= n - 2*min_shift + 1 with min_shift = 100"),
        "ladder": ladder,
        "per_feature": per_feature_floor,
        "declaration_rule": ("§P7-10(c): 'F9-10 declares its floor from the measured "
                            "VIF_mark before it runs, or it does not run.'"),
    }

    out["headline"] = {
        "n_mark_rows_primary": prim["n_mark_rows"],
        "n_censored_excluded_primary": prim["censoring_P7_10a"]["n_censored_excluded"],
        "vif_mark_median_all": med_all,
        "vif_mark_median_df1": med_df1,
        "vif_mark_median_df2": med_df2,
        "vif_mark_iqr_all": [prim["vif_mark_all"]["q25"], prim["vif_mark_all"]["q75"]],
        "rho_min_at_measured_median": rho_min(med_all, n_ev),
        "rho_min_at_vif_1": rho_min(1.0, n_ev),
        "count_path_vif_df2_for_contrast": 24.08,
    }

    out["caveats"] = [
        "SCOPE. This measures the inflation of the MARK null as it is currently "
        "constructed (circular moving-block bootstrap over the time-ordered event "
        "sequence, block_events = 500 over 46,585 events = ~93 blocks). It does NOT "
        "decompose that inflation into (a) real mark-axis dependence in the Earth "
        "and (b) width manufactured by a ~93-block bootstrap. §P7-10(c) names that "
        "decomposition as UNMEASURED, and it stays unmeasured here: F4-58M is the "
        "count-path F4-58 item (i) analogue, not the item (ii) simulated control.",
        "The chi-square-scale conversion uses the standard large-sample nulls "
        "((n-1)*rho^2 ~ chi2_1; n*R^2 ~ chi2_2). At n = 46,585 the asymptotics are "
        "not in question, but the conversion is an assumption of the READING, not "
        "of the recorded statistic.",
        "VIF_mark is computed against p_block_bootstrap, per §P7-10(c)'s literal "
        "formula. `vif_mark_from_p_raw` is carried per row for reference; p_raw is "
        "the conservative max(circular-shift, block-bootstrap) and where the shift "
        "arm binds it gives a LARGER VIF, so the reported medians are the smaller "
        "(more conservative for the floor's denominator, less conservative for the "
        "floor itself) of the two readings.",
        "Mark rows carry no lag axis (lag = null in every row on disk), so no lag "
        "reduction is involved and none is claimed.",
        "The three censored rows are the same three (feature, mark) pairs in all "
        "four sessions -- b_value_90d x mag, deep_fraction_30d x depth, "
        "mean_depth_30d x depth. Their VIF is an UPPER BOUND (17.2, 28.2, 38.7 in "
        "the primary session at B = 50,000; tightened to 13.9, 22.7, 31.1 at "
        "B = 500,000 by the §P7-11(c) run, still at the floor) and is reported, "
        "not medianed. CORRECTION 2026-08-12: an earlier version of this caveat "
        "called these LOWER bounds and cited them as evidence that the §P7-10(a) "
        "seam is anti-conservative. The sign is inverted, measured directly by "
        "pushing B up; excluding these rows removes OVER-stated values and so "
        "mildly LOWERS the floor. The exclusion rule still stands -- a bound is "
        "not a measurement -- but not for the reason given.",
    ]

    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print("wrote %s" % OUT_PATH)

    for s in out["sessions"]:
        print("\n=== %s" % s["session"])
        print("    mark rows %d, scored %d, CENSORED-EXCLUDED %d (floor p_boot = %.3e)"
              % (s["n_mark_rows"], s["censoring_P7_10a"]["n_scored"],
                 s["censoring_P7_10a"]["n_censored_excluded"],
                 s["resolution_floor_p_boot"]))
        print("    VIF_mark all:", s["vif_mark_all"])
        print("    VIF_mark df1:", s["vif_mark_df1_spearman"])
        print("    VIF_mark df2:", s["vif_mark_df2_circular_linear"])
        for c in s["censoring_P7_10a"]["censored_rows"]:
            print("      CENSORED %-18s %-6s p_boot=%.3e  VIF<=%.2f"
                  % (c["feature"], c["mark"], c["p_block_bootstrap"], c["vif_mark"]))

    print("\n=== OPERATIVE MARK FLOOR (alpha = 1e-4, n = %d) ===" % n_ev)
    for k, v in ladder.items():
        print("  %-56s VIF=%8.3f  rho_min = %.4f" % (k, v["vif_mark"], v["rho_min"]))

    print("\n=== per-feature (primary), sorted by VIF_mark ===")
    print("%-24s %-6s %3s %10s %10s %10s %8s"
          % ("feature", "mark", "df", "stat", "VIF_mark", "rho_min", "obs/floor"))
    for f in per_feature_floor:
        print("%-24s %-6s %3d %10.5f %10.3f %10.5f %8.2f"
              % (f["feature"], f["mark"], f["df"], f["observed_statistic"],
                 f["vif_mark"], f["rho_min"], f["observed_over_floor"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
