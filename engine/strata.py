"""Stratified (weighted) BH and the partition-invariant max-statistic -- §P6-3.

THE CORRECTION §P6-3 MAKES TO THE STATED MOTIVE, up front, because it is the whole
point of this module: running BH at a flat q independently inside each of S families
controls the AVERAGE OVER FAMILIES of FDR (Benjamini & Bogomolov 2014), NOT the
overall FDR across all tests. "So dead families do not tax live ones" is not what
that buys. Overall control requires the weighted-BH constraint

    sum_s ( m_s * q_s )  ==  m * q          (m = sum_s m_s, q = 0.10)

which IS weighted BH with pre-specified weights and DOES control global FDR at q
exactly. This module asserts that identity and REFUSES TO RUN when it fails
(`assert_budget_identity`), which is rule 1.

The five rules and where they live:

  RULE 1  the global accounting identity, enforced in code, refusal on violation,
          and the (stratum, m_s, q_s, threshold_s) table printed in the report.
                                    -> `assert_budget_identity`, `stratified_bh`,
                                       `budget_table_lines`
  RULE 2  the honest trade, stated in the report: a dead family STILL consumes
          budget through its m_s; the only way a live family gets more is if another
          gets less, and that reallocation is a PRIOR, frozen in the config hash
          with its justification written into the partition's own note field BEFORE
          the run. Reallocating after seeing which strata look alive is S-9's
          forking path with extra steps and is forbidden.
                                    -> `HONEST_TRADE`, `load_partition` (requires
                                       a non-empty `note` on every reallocation)
  RULE 3  strata frozen in the config hash before the run, defined as
          (feature_family x test_kind [x region]), m_s a declared integer per
          stratum; tests that ERROR OUT count as NON-REJECTIONS AGAINST THEIR
          DECLARED m_s -- the denominator is the declaration, not the execution.
                                    -> `load_partition` (hash), `stratum_of`,
                                       `stratified_bh` (m_s is the denominator)
  RULE 4  the anti-repartition guarantee is S-8, not BH: the sim-calibrated
          max-statistic p per stratum AND globally, the global one INVARIANT to how
          the family is partitioned, and no stratum may be reported without the
          global figure adjacent.
                                    -> `max_statistic_p`, `max_statistic_report`
  RULE 5  re-partitioning is a NEW SESSION: new config hash, new EXPLORE_COUNT
          line, and the prior partition's report is neither deleted nor amended.
                                    -> the partition file's sha256 is hash-affecting
                                       (`load_partition` -> cfg["strata"]), so a
                                       changed partition cannot resume an old
                                       session directory.

DEFAULT IS UNSTRATIFIED. Flat BH over one implicit stratum remains the default and
is bit-identical to the pre-§P6-3 path. Stratification is opt-in via `--strata
<partition.json>` and the partition file's CONTENT is hash-affecting.
"""

import hashlib
import json
import math
import os

import numpy as np

UNSTRATIFIED = "UNSTRATIFIED"

HONEST_TRADE = (
    "THE HONEST TRADE (§P6-3(2)). Under the identity sum_s m_s q_s = m q a DEAD "
    "family still consumes budget through its m_s. The only way a live family gets "
    "a larger q_s is if another family gets a smaller one, and that reallocation is "
    "a PRIOR: it is frozen in the config hash with its justification -- physics or a "
    "prior ledger result -- written into the partition's note field BEFORE the run. "
    "Reallocating after seeing which strata look alive is S-9's forking path with "
    "extra steps, and is forbidden.")

MAXSTAT_NOTE = (
    "THE ANTI-REPARTITION GUARANTEE IS S-8, NOT BH (§P6-3(4)). The global "
    "max-statistic p is computed over the WHOLE declared family from surrogates "
    "carried on a COMMON replicate index, so it is invariant to how that family is "
    "partitioned -- which is exactly what makes stratification a discipline rather "
    "than a knob. Any cross-stratum statement is printed adjacent to it.")


# ------------------------------------------------------------- partition loading
def _canon(v):
    return None if v is None else str(v)


def stratum_key(feature_family, test_kind, region=None):
    """(feature_family x test_kind [x region]) -- the declared partition axes, §P6-3(3).

    `region` is RESERVED and defaulted to None now, so that adding per-region strata
    per §P6-4 later does not re-key the strata that already exist.
    """
    parts = [_canon(feature_family), _canon(test_kind)]
    if region is not None:
        parts.append(_canon(region))
    return "|".join("*" if p is None else p for p in parts)


def stratum_of(row, partition):
    """The stratum a result row belongs to. Falls back to the declared catch-all."""
    if partition is None:
        return UNSTRATIFIED
    key = stratum_key(row.get("family"), _test_kind(row), row.get("region"))
    if key in partition["by_key"]:
        return partition["by_key"][key]
    key2 = stratum_key(row.get("family"), _test_kind(row), None)
    if key2 in partition["by_key"]:
        return partition["by_key"][key2]
    if partition.get("catch_all"):
        return partition["catch_all"]
    raise KeyError(
        f"no declared stratum for (family={row.get('family')!r}, "
        f"test_kind={_test_kind(row)!r}, region={row.get('region')!r}) and the "
        f"partition declares no catch_all. Every test must be priced in a declared "
        f"stratum before the run (§P6-3(3)). Refusing to run.")


def _test_kind(row):
    """The COARSE test kind: the three kinds this engine runs."""
    t = row.get("test")
    if t == "glm_poisson_offset_etas":
        return "glm"
    if t == "lomb_scargle_peak":
        return "period"
    if t in ("spearman", "circular-linear") or row.get("mark") is not None:
        return "mark"
    return str(t)


def load_partition(path, q=0.10):
    """Read, validate and HASH a declared partition file. Content is hash-affecting.

    Schema (JSON):
      {"q": 0.10,
       "catch_all": "other",            # optional; a name that must also appear below
       "strata": [
         {"name": "eph-glm", "feature_family": 1, "test_kind": "glm",
          "region": null, "m_s": 240, "q_s": 0.12,
          "note": "why this stratum gets more/less than flat -- REQUIRED whenever
                   q_s != q, per §P6-3(2)"},
         ...]}

    The returned dict is what goes into the config (and therefore into the config
    hash): the file's sha256, its declared strata, and nothing derived from data.
    """
    with open(path, "rb") as fh:
        raw = fh.read()
    sha = hashlib.sha256(raw).hexdigest()
    doc = json.loads(raw.decode("utf-8"))
    strata = doc.get("strata") or []
    if not strata:
        raise ValueError(f"partition file {path} declares no strata")
    q_decl = float(doc.get("q", q))
    seen, by_key, out = set(), {}, []
    for s in strata:
        name = str(s["name"])
        if name in seen:
            raise ValueError(f"partition file {path}: duplicate stratum {name!r}")
        seen.add(name)
        m_s = int(s["m_s"])
        q_s = float(s["q_s"])
        if m_s < 0:
            raise ValueError(f"stratum {name!r}: m_s must be >= 0")
        if not (0.0 < q_s <= 1.0):
            raise ValueError(f"stratum {name!r}: q_s must be in (0, 1]")
        note = str(s.get("note", "")).strip()
        # RULE 2: any DEPARTURE from the flat budget is a prior and must carry its
        # justification in the file itself, before the run.
        if abs(q_s - q_decl) > 1e-12 and not note:
            raise ValueError(
                f"stratum {name!r} sets q_s = {q_s} != q = {q_decl} with an empty "
                f"note. §P6-3(2): a reallocation is a PRIOR and its justification "
                f"(physics, prior ledger result) must be written into the config "
                f"BEFORE the run. Refusing to run.")
        rec = {"name": name, "feature_family": s.get("feature_family"),
               "test_kind": s.get("test_kind"), "region": s.get("region"),
               "m_s": m_s, "q_s": q_s, "note": note}
        out.append(rec)
        by_key[stratum_key(rec["feature_family"], rec["test_kind"],
                           rec["region"])] = name
    catch_all = doc.get("catch_all")
    if catch_all is not None and str(catch_all) not in seen:
        raise ValueError(f"catch_all {catch_all!r} is not a declared stratum")
    return {"path": os.path.basename(path), "sha256": sha, "q": q_decl,
            "strata": out, "by_key": by_key,
            "catch_all": None if catch_all is None else str(catch_all),
            "m": int(sum(r["m_s"] for r in out)),
            "note": str(doc.get("note", ""))}


# ------------------------------------------------------- rule 1: budget identity
class BudgetIdentityError(ValueError):
    """Raised when sum_s m_s q_s != m q. The engine refuses to run."""


def assert_budget_identity(strata, q, rtol=1e-9):
    """RULE 1. sum_s (m_s * q_s) == m * q, or the engine refuses to run.

    This is not a sanity check on a config file, it is the thing that makes the
    procedure control global FDR at all. Per-stratum BH at a FLAT q -- the obvious
    implementation, and the one the motive as originally stated describes -- fails
    this identity whenever S > 1 with unequal m_s, and controls only the average
    over families. So the assertion is the difference between weighted BH and a
    procedure with no global guarantee, and it refuses rather than warns.
    """
    m = int(sum(int(s["m_s"]) for s in strata))
    lhs = float(sum(int(s["m_s"]) * float(s["q_s"]) for s in strata))
    rhs = float(m) * float(q)
    if not math.isclose(lhs, rhs, rel_tol=rtol, abs_tol=1e-12):
        raise BudgetIdentityError(
            f"§P6-3(1) VIOLATED: sum_s m_s q_s = {lhs:.12g} but m q = {rhs:.12g} "
            f"(m = {m}, q = {q}). Weighted BH controls global FDR at q ONLY under "
            f"this identity; without it the procedure controls the average over "
            f"families and nothing global. Refusing to run.\n"
            + "\n".join(f"    {s['name']}: m_s = {s['m_s']}, q_s = {s['q_s']}"
                        for s in strata))
    return {"m": m, "q": float(q), "sum_ms_qs": lhs, "m_q": rhs, "ok": True}


def flat_partition(m, q=0.10, name=UNSTRATIFIED):
    """The default: one stratum holding the whole declared family. Satisfies rule 1."""
    return [{"name": name, "feature_family": None, "test_kind": None,
             "region": None, "m_s": int(m), "q_s": float(q), "note": ""}]


# ------------------------------------------------------------ rule 1/3: the BH --
def _bh_within(pvals, m_s, q_s):
    """BH step-up inside one stratum with the DECLARED m_s as denominator.

    RULE 3, the load-bearing detail: `m_s` is the DECLARATION, not `len(pvals)`. A
    test that errored out, or that was never executed, counts as a NON-REJECTION
    against its declared m_s -- it does not shrink the denominator, which would
    quietly make every surviving test easier to reject.
    """
    p = np.asarray(pvals, dtype=np.float64)
    n = p.size
    m_s = int(m_s)
    if n == 0 or m_s <= 0:
        return (np.array([]), np.array([], dtype=bool), 0.0)
    if n > m_s:
        raise ValueError(
            f"stratum received {n} executed tests but declares m_s = {m_s}. The "
            f"declared count must be an upper bound on what runs in it; a "
            f"denominator smaller than the family is exactly the multiplicity "
            f"understatement this engine refuses. Refusing to run.")
    order = np.argsort(p, kind="stable")
    ps = p[order]
    adj = ps * m_s / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    out = np.empty(n)
    out[order] = adj
    passed = out <= q_s
    # the largest p that is rejected -- the stratum's operating threshold
    thresh = float(ps[passed[order]].max()) if bool(passed.any()) else 0.0
    return out, passed, thresh


def stratified_bh(pvals, stratum_names, strata, q, eligible=None):
    """Weighted BH: BH inside each declared stratum at its own q_s, identity asserted.

    `eligible` is an optional boolean mask. An INELIGIBLE test (§P6-2(1): a test
    whose p is UNRESOLVED) still counts against its declared m_s -- same accounting
    as an errored test -- but can never be rejected.
    """
    ident = assert_budget_identity(strata, q)
    p = np.asarray(pvals, dtype=np.float64)
    n = p.size
    names = list(stratum_names)
    if len(names) != n:
        raise ValueError("pvals and stratum_names differ in length")
    elig = (np.ones(n, dtype=bool) if eligible is None
            else np.asarray(eligible, dtype=bool))
    by_name = {s["name"]: s for s in strata}
    unknown = sorted(set(names) - set(by_name))
    if unknown:
        raise KeyError(f"results carry undeclared strata {unknown}")
    qv = np.ones(n)
    passed = np.zeros(n, dtype=bool)
    table = []
    for s in strata:
        idx = np.array([i for i, nm in enumerate(names) if nm == s["name"]],
                       dtype=int)
        run_idx = np.array([i for i in idx if elig[i]], dtype=int)
        qq, pp, thr = _bh_within(p[run_idx], s["m_s"], s["q_s"]) if run_idx.size \
            else (np.array([]), np.array([], dtype=bool), 0.0)
        if run_idx.size:
            qv[run_idx] = qq
            passed[run_idx] = pp
        table.append({
            "stratum": s["name"], "m_s": int(s["m_s"]), "q_s": float(s["q_s"]),
            "threshold_s": float(thr),
            "bh_line_smallest": float(s["q_s"] / max(int(s["m_s"]), 1)),
            "n_executed": int(idx.size),
            "n_eligible": int(run_idx.size),
            "n_not_rejected_by_accounting": int(s["m_s"] - run_idx.size),
            "n_pass": int(pp.sum()) if run_idx.size else 0,
            "expected_false_discoveries": float(s["m_s"] * s["q_s"]),
            "note": s.get("note", ""),
        })
    return qv, passed, {"table": table, "identity": ident}


def budget_table_lines(table, identity):
    """RULE 1's printed table. Markdown, used verbatim by report and console."""
    lines = ["| stratum | m_s | q_s | BH line (smallest) | threshold_s | executed |"
             " ineligible/errored | survivors | E[false] |",
             "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for r in table:
        thr = r["threshold_s"]
        thr_s = f"{thr:.3g}" if thr > 0 else "none (no rejection)"
        lines.append(
            f"| `{r['stratum']}` | {r['m_s']} | {r['q_s']:.4g} | "
            f"{r['bh_line_smallest']:.3g} | {thr_s} | "
            f"{r['n_executed']} | {r['n_not_rejected_by_accounting']} | "
            f"{r['n_pass']} | {r['expected_false_discoveries']:.2f} |")
    lines.append("")
    lines.append(f"Identity check (§P6-3(1)): sum_s m_s q_s = "
                 f"{identity['sum_ms_qs']:.6g}, m q = {identity['m_q']:.6g} "
                 f"(m = {identity['m']}, q = {identity['q']}). ASSERTED IN CODE; "
                 f"the engine refuses to run on violation.")
    return lines


# --------------------------------------------------- rule 4: the max-statistic --
def _null_ranks(col_null, obs):
    """-log10 of the empirical p of each replicate WITHIN its own test's null.

    Standardisation is not cosmetic here: the max is taken ACROSS tests, and raw
    score statistics from a 1-df GLM, a rank correlation and a Lomb-Scargle peak are
    not on a common scale, so a max over the raw numbers would be a max over units.
    Ranking each column inside its own null puts every test on the p-value scale,
    which is the only scale on which "the most extreme test in the family" means
    anything. This is the min-p / max-statistic construction.
    """
    x = np.asarray(col_null, dtype=np.float64).ravel()
    n = x.size
    # p of a value v within this column's null, leave-one-out style:
    #   p(v) = (1 + #{x >= v}) / (1 + n)
    order = np.argsort(-x, kind="stable")
    ranks = np.empty(n, dtype=np.float64)
    srt = x[order]
    # #{x >= srt[i]} via searchsorted on the descending-sorted array
    ge = np.searchsorted(-srt, -srt, side="right")
    ranks[order] = ge
    p_null = (1.0 + ranks) / (1.0 + n)
    ge_obs = float((x >= float(obs)).sum())
    p_obs = (1.0 + ge_obs) / (1.0 + n)
    return -np.log10(p_null), -np.log10(p_obs)


def max_statistic_p(null_matrix, observed, columns=None):
    """S-8's sim-calibrated max-statistic p over a set of tests sharing a replicate index.

    `null_matrix` is (n_replicates, n_tests): entry [j, i] is test i's statistic on
    the SAME surrogate replicate j. That shared index is what makes this a joint
    null and therefore exact under arbitrary dependence between tests -- no PRDS
    assumption, no Bonferroni, no correction factor.

    `columns` restricts to a subset (one stratum). PARTITION INVARIANCE: called with
    columns=None the answer depends on the set of tests and the surrogates and on
    NOTHING about how they were grouped, which is §P6-3(4) exactly.
    """
    A = np.asarray(null_matrix, dtype=np.float64)
    obs = np.asarray(observed, dtype=np.float64).ravel()
    if A.ndim != 2 or A.shape[1] != obs.size:
        raise ValueError(f"null_matrix {A.shape} incompatible with {obs.size} observed")
    cols = np.arange(A.shape[1]) if columns is None else np.asarray(columns, dtype=int)
    if cols.size == 0:
        return {"p": float("nan"), "n_tests": 0, "n_replicates": int(A.shape[0]),
                "t_obs": float("nan"), "floor": float("nan")}
    T = np.empty((A.shape[0], cols.size))
    t_obs = np.empty(cols.size)
    for k, i in enumerate(cols):
        T[:, k], t_obs[k] = _null_ranks(A[:, i], obs[i])
    tn = T.max(axis=1)
    to = float(t_obs.max())
    n = tn.size
    ge = int((tn >= to).sum())
    return {"p": float((1.0 + ge) / (1.0 + n)), "n_tests": int(cols.size),
            "n_replicates": int(n), "t_obs": to,
            "floor": float(1.0 / (n + 1.0)),
            "statistic": "max over tests of -log10(empirical p within own null)"}


def max_statistic_report(null_matrix, observed, stratum_names, strata=None):
    """RULE 4. Per-stratum AND global, with the global one carried on every row.

    The returned rows deliberately each carry `global_p`: §P6-3(4) says NO STRATUM
    MAY BE REPORTED ALONE, so the data structure makes the global figure impossible
    to drop when a caller renders one row.
    """
    names = list(stratum_names)
    glob = max_statistic_p(null_matrix, observed)
    rows = []
    order = ([s["name"] for s in strata] if strata else
             sorted(set(names), key=names.index))
    for nm in order:
        cols = [i for i, x in enumerate(names) if x == nm]
        r = max_statistic_p(null_matrix, observed, columns=cols) if cols else \
            {"p": float("nan"), "n_tests": 0,
             "n_replicates": int(np.asarray(null_matrix).shape[0]),
             "t_obs": float("nan"), "floor": float("nan")}
        r["stratum"] = nm
        r["global_p"] = glob["p"]
        r["global_n_tests"] = glob["n_tests"]
        rows.append(r)
    return {"global": glob, "per_stratum": rows, "note": MAXSTAT_NOTE}
