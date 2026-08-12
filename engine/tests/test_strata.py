"""§P6-3, rule by rule, plus the default-is-unstratified pin.

  rule 1  sum_s m_s q_s == m q, ASSERTED, and the engine REFUSES on violation
  rule 2  the honest trade: a reallocation without a written justification is refused
  rule 3  strata frozen in the config hash; errored/ineligible tests count as
          non-rejections against the DECLARED m_s
  rule 4  the max-statistic per stratum AND globally, the global one
          PARTITION-INVARIANT, and never a stratum row without it
  rule 5  re-partitioning is a new session: a changed partition file is a changed
          config hash
"""

import json
import types

import numpy as np
import pytest

from engine import mine as M
from engine import mine_session as ms
from engine import splits
from engine import strata as S


def _partition_doc(q=0.10):
    return {
        "q": q,
        "note": "phase-2a test partition",
        "strata": [
            {"name": "eph-glm", "feature_family": 1, "test_kind": "glm",
             "region": None, "m_s": 100, "q_s": 0.15,
             "note": "ephemeris GLM: prior ledger weight (K-089-R tranche 1)"},
            {"name": "cat-glm", "feature_family": 2, "test_kind": "glm",
             "region": None, "m_s": 100, "q_s": 0.05,
             "note": "catalogue-derived GLM: dead in v1, down-weighted a priori"},
        ],
    }


def _write(tmp_path, doc, name="partition.json"):
    p = tmp_path / name
    p.write_text(json.dumps(doc), encoding="utf-8")
    return str(p)


# ------------------------------------------------------------------ rule 1 -----
def test_rule1_identity_holds_for_a_valid_partition(tmp_path):
    part = S.load_partition(_write(tmp_path, _partition_doc()))
    ident = S.assert_budget_identity(part["strata"], 0.10)
    assert ident["sum_ms_qs"] == pytest.approx(ident["m_q"])
    assert ident["m"] == 200


def test_rule1_engine_refuses_to_run_when_the_identity_fails(tmp_path):
    doc = _partition_doc()
    doc["strata"][1]["q_s"] = 0.10          # 100*0.15 + 100*0.10 != 200*0.10
    part = S.load_partition(_write(tmp_path, doc))
    with pytest.raises(S.BudgetIdentityError) as e:
        S.assert_budget_identity(part["strata"], 0.10)
    assert "Refusing to run" in str(e.value)
    assert "average over families" in str(e.value)


def test_rule1_refusal_reaches_build_config(tmp_path):
    """The refusal is not a library nicety: it stops a RUN from being configured."""
    doc = _partition_doc()
    doc["strata"][0]["q_s"] = 0.30
    path = _write(tmp_path, doc)
    args = types.SimpleNamespace(
        mag_target=4.5, dlat=1.0, dlon=1.0, explore_frac=0.7,
        data_dir="data", no_download=True, seed=1, tranche1=False, ladder=False,
        gpd=False, strata=path)
    with pytest.raises(S.BudgetIdentityError):
        ms.build_config(args, ms.QUICK)


def test_rule1_flat_partition_satisfies_the_identity():
    S.assert_budget_identity(S.flat_partition(137, 0.10), 0.10)


def test_rule1_budget_table_prints_every_declared_column(tmp_path):
    part = S.load_partition(_write(tmp_path, _partition_doc()))
    p = np.full(20, 0.5)
    names = ["eph-glm"] * 10 + ["cat-glm"] * 10
    _q, _pa, meta = S.stratified_bh(p, names, part["strata"], 0.10)
    lines = S.budget_table_lines(meta["table"], meta["identity"])
    body = "\n".join(lines)
    for col in ("m_s", "q_s", "threshold_s"):
        assert col in body
    assert "eph-glm" in body and "cat-glm" in body
    assert "sum_s m_s q_s" in body


# ------------------------------------------------------------------ rule 2 -----
def test_rule2_a_reallocation_without_a_written_justification_is_refused(tmp_path):
    doc = _partition_doc()
    doc["strata"][0]["note"] = ""
    with pytest.raises(ValueError) as e:
        S.load_partition(_write(tmp_path, doc))
    assert "PRIOR" in str(e.value) and "BEFORE the run" in str(e.value)


def test_rule2_a_flat_stratum_needs_no_note(tmp_path):
    doc = {"q": 0.10, "strata": [
        {"name": "a", "feature_family": 1, "test_kind": "glm", "m_s": 50,
         "q_s": 0.10},
        {"name": "b", "feature_family": 2, "test_kind": "glm", "m_s": 50,
         "q_s": 0.10}]}
    part = S.load_partition(_write(tmp_path, doc))
    S.assert_budget_identity(part["strata"], 0.10)


def test_rule2_the_honest_trade_text_is_carried_for_the_report():
    assert "DEAD family" in S.HONEST_TRADE
    assert "forbidden" in S.HONEST_TRADE


# ------------------------------------------------------------------ rule 3 -----
def test_rule3_partition_content_is_hash_affecting(tmp_path):
    base = types.SimpleNamespace(
        mag_target=4.5, dlat=1.0, dlon=1.0, explore_frac=0.7, data_dir="data",
        no_download=True, seed=1, tranche1=False, ladder=False, gpd=False,
        strata=None)
    h_flat = splits.config_hash(ms.build_config(base, ms.QUICK))
    a = _write(tmp_path, _partition_doc(), "a.json")
    base.strata = a
    h_a = splits.config_hash(ms.build_config(base, ms.QUICK))
    doc_b = _partition_doc()
    doc_b["strata"][0]["note"] += " (revised wording)"
    b = _write(tmp_path, doc_b, "b.json")
    base.strata = b
    h_b = splits.config_hash(ms.build_config(base, ms.QUICK))
    assert h_flat != h_a != h_b and h_a != h_b


def test_rule3_errored_tests_count_as_non_rejections_against_declared_m_s():
    """The denominator is the DECLARATION, not the execution.

    Ten tests declared, only four executed. Their q-values must be the ones BH
    computes at m_s = 10, i.e. exactly what padding the missing six with p = 1 would
    give -- not the smaller denominators that would make every survivor cheaper.
    """
    p_exec = np.array([0.001, 0.02, 0.04, 0.5])
    strata = [{"name": "s", "m_s": 10, "q_s": 0.10}]
    q_short, pass_short, _ = S.stratified_bh(p_exec, ["s"] * 4, strata, 0.10)
    padded = np.concatenate([p_exec, np.ones(6)])
    q_pad, _pass_pad, _ = S.stratified_bh(padded, ["s"] * 10, strata, 0.10)
    assert np.allclose(q_short, q_pad[:4])
    # and the executed-vs-declared gap is reported, not hidden
    _q, _pa, meta = S.stratified_bh(p_exec, ["s"] * 4, strata, 0.10)
    assert meta["table"][0]["n_not_rejected_by_accounting"] == 6


def test_rule3_ineligible_tests_also_count_against_m_s():
    p = np.array([1e-6, 0.02, 0.04, 0.5])
    elig = np.array([False, True, True, True])
    strata = [{"name": "s", "m_s": 4, "q_s": 0.10}]
    q, passed, meta = S.stratified_bh(p, ["s"] * 4, strata, 0.10, eligible=elig)
    assert passed[0] == False           # noqa: E712 -- ineligible can never reject
    assert meta["table"][0]["n_eligible"] == 3
    assert meta["table"][0]["n_not_rejected_by_accounting"] == 1


def test_rule3_a_stratum_smaller_than_its_executed_count_is_refused():
    with pytest.raises(ValueError):
        S.stratified_bh(np.full(5, 0.5), ["s"] * 5,
                        [{"name": "s", "m_s": 3, "q_s": 0.10}], 0.10)


def test_rule3_undeclared_stratum_is_refused():
    with pytest.raises(KeyError):
        S.stratified_bh(np.array([0.1]), ["ghost"],
                        [{"name": "s", "m_s": 1, "q_s": 0.10}], 0.10)


def test_rule3_stratum_of_uses_family_x_kind_and_reserves_region(tmp_path):
    part = S.load_partition(_write(tmp_path, _partition_doc()))
    assert S.stratum_of({"family": 1, "test": "glm_poisson_offset_etas"},
                        part) == "eph-glm"
    assert S.stratum_of({"family": 2, "test": "glm_poisson_offset_etas"},
                        part) == "cat-glm"
    # region is a declared axis that simply is not used yet
    assert S.stratum_key(1, "glm") != S.stratum_key(1, "glm", "socal")
    with pytest.raises(KeyError):
        S.stratum_of({"family": 9, "test": "lomb_scargle_peak"}, part)


# ------------------------------------------------------------------ rule 4 -----
def _maxstat_fixture(seed=0, n_rep=500, n_tests=6):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n_rep, n_tests)) ** 2
    obs = np.array([float(np.quantile(A[:, i], 0.99)) for i in range(n_tests)])
    return A, obs


def test_rule4_global_max_statistic_is_partition_invariant():
    A, obs = _maxstat_fixture()
    part_a = ["x", "x", "x", "y", "y", "y"]
    part_b = ["p", "q", "p", "q", "p", "q"]
    ra = S.max_statistic_report(A, obs, part_a)
    rb = S.max_statistic_report(A, obs, part_b)
    assert ra["global"]["p"] == rb["global"]["p"]
    assert ra["global"]["n_tests"] == rb["global"]["n_tests"] == 6


def test_rule4_no_stratum_row_can_be_reported_without_the_global_number():
    A, obs = _maxstat_fixture()
    rep = S.max_statistic_report(A, obs, ["x", "x", "x", "y", "y", "y"])
    assert rep["per_stratum"]
    for r in rep["per_stratum"]:
        assert r["global_p"] == rep["global"]["p"]
        assert "n_tests" in r


def test_rule4_max_statistic_is_a_valid_p_and_respects_its_floor():
    A, obs = _maxstat_fixture(seed=5)
    r = S.max_statistic_p(A, obs)
    assert r["floor"] == pytest.approx(1.0 / (A.shape[0] + 1.0))
    assert r["floor"] <= r["p"] <= 1.0


def test_rule4_max_statistic_is_uniform_under_the_null():
    """A joint null drawn from the same machinery must give a ~uniform p."""
    rng = np.random.default_rng(17)
    ps = []
    for _ in range(120):
        A = rng.standard_normal((300, 5)) ** 2
        obs = rng.standard_normal(5) ** 2
        ps.append(S.max_statistic_p(A, obs)["p"])
    ps = np.array(ps)
    assert 0.03 < float((ps <= 0.10).mean()) < 0.25


def test_rule4_a_planted_common_signal_is_detected():
    rng = np.random.default_rng(19)
    A = rng.standard_normal((800, 4)) ** 2
    obs = np.full(4, float(A.max()) * 1.5)
    assert S.max_statistic_p(A, obs)["p"] == pytest.approx(1.0 / 801.0)


def test_rule4_engine_max_statistic_matrix_shares_one_replicate_index():
    """The columns must be the SAME surrogate worlds, or the max is meaningless."""
    W = _mini_payload()
    tests = [{"test": "glm_poisson_offset_etas", "feature": f.name, "lag": 0}
             for f in W["feats"]]
    # a mark test and a period peak are present and must be EXCLUDED, not folded in
    tests.append({"test": "lomb_scargle_peak", "feature": "period_29d"})
    tests.append({"test": "spearman", "feature": W["feats"][0].name,
                  "mark": "mag"})
    A, obs, idx = ms.max_statistic_matrix(W["feats"], W["window"], W["counts"],
                                          W["offset"], tests)
    assert A.shape[1] == len(idx) == len(tests) - 2
    assert obs.size == A.shape[1]
    # every column is evaluated on the SAME admissible shift set, so a row of A is
    # one surrogate world seen by every test at once
    assert A.shape[0] == W["counts"].size - 2 * 30 + 1


def _mini_payload(n_days=900, n_feats=3):
    rng = np.random.default_rng(4)
    offset = 0.9 * (1 + 0.2 * np.sin(2 * np.pi * np.arange(n_days) / 365.25))
    counts = rng.poisson(offset).astype(float)
    t0 = __import__("datetime").datetime(2000, 1, 1)
    feats = M.ephemeris_features(t0, n_days)[:n_feats]
    return {"feats": feats, "window": slice(0, n_days), "counts": counts,
            "offset": offset}


# ------------------------------------------------------------------ rule 5 -----
def test_rule5_repartitioning_produces_a_different_config_hash(tmp_path):
    base = types.SimpleNamespace(
        mag_target=4.5, dlat=1.0, dlon=1.0, explore_frac=0.7, data_dir="data",
        no_download=True, seed=1, tranche1=False, ladder=False, gpd=False,
        strata=_write(tmp_path, _partition_doc(), "p1.json"))
    h1 = splits.config_hash(ms.build_config(base, ms.QUICK))
    doc2 = _partition_doc()
    doc2["strata"][0]["m_s"] = 120
    doc2["strata"][1]["m_s"] = 80
    doc2["strata"][0]["q_s"] = 0.125
    doc2["strata"][1]["q_s"] = 0.0625
    base.strata = _write(tmp_path, doc2, "p2.json")
    h2 = splits.config_hash(ms.build_config(base, ms.QUICK))
    assert h1 != h2


# --------------------------------------------- the default is UNSTRATIFIED -----
def test_default_flat_path_reproduces_plain_benjamini_hochberg():
    """The pin that keeps stratification opt-in: with one stratum holding the whole
    declared family, `stratified_bh` IS `M.benjamini_hochberg`."""
    rng = np.random.default_rng(23)
    p = rng.uniform(0.0, 1.0, 200) ** 3
    q_ref, pass_ref = M.benjamini_hochberg(p, M.FDR_Q)
    q_new, pass_new, meta = S.stratified_bh(
        p, [S.UNSTRATIFIED] * p.size, S.flat_partition(p.size, M.FDR_Q), M.FDR_Q)
    assert np.allclose(q_ref, q_new)
    assert np.array_equal(pass_ref, pass_new)
    assert meta["identity"]["ok"] is True


def test_default_config_has_no_strata_key():
    args = types.SimpleNamespace(
        mag_target=4.5, dlat=1.0, dlon=1.0, explore_frac=0.7, data_dir="data",
        no_download=True, seed=1, tranche1=False, ladder=False, gpd=False,
        strata=None)
    cfg = ms.build_config(args, ms.QUICK)
    assert "strata" not in cfg and "gpd" not in cfg
