"""§P6-2, rule by rule. Every numbered rule that can be unit-tested, is.

  rule 1(a) AD gate rejects a non-GPD tail      test_rule1a_*
  rule 1(b) xi-stability gate                   test_rule1b_*
  rule 2    CI-upper, never the point estimate  test_rule2_*
  rule 3    one-decade floor                    test_rule3_*
  rule 4    all-shifts prohibition (structural) test_rule4_*
  rule 5    label on every row + census         test_rule5_*
  rule 6    candidate, not stub; N >= 10/p      test_rule6_*
  rule 7    the harness exists and is honest    engine/tests/test_gpd_audit.py
"""

import json
import math
import os

import numpy as np
import pytest

from engine import gpd_tail as G
from engine import mine_session as ms

FAST = dict(n_boot=60, ad_boot=60)


def _exp_surrogates(n=2000, seed=0):
    """Exponential(1) statistics: an exactly-GPD tail (xi = 0), the easy case."""
    return np.random.default_rng(seed).standard_exponential(n)


# ---------------------------------------------------------------- rule 1(a) ----
def test_rule1a_ad_gate_rejects_a_non_gpd_tail():
    """A tail that is DISCRETE and lumpy is not GPD, and the gate must say so.

    The exceedances here take a handful of distinct values with big gaps. No GPD
    fits that, and the Anderson-Darling statistic against the fitted GPD is huge,
    so the parametric bootstrap must return p_AD ~ 0 and the estimator must refuse.
    """
    rng = np.random.default_rng(3)
    bulk = rng.uniform(0.0, 1.0, 1800)
    lumps = np.repeat(np.array([2.0, 2.001, 5.0, 5.001, 9.0]), 40)
    S = np.concatenate([bulk, lumps])
    out = G.gpd_tail_p(S, 12.0, 2000, "block_bootstrap",
                       np.random.default_rng(4), **FAST)
    assert out["p_method"] == G.P_UNRESOLVED
    assert "Anderson-Darling" in out["reason"]
    assert out["ad"]["p_ad"] < G.GPD_AD_ALPHA
    assert out["ad"]["pass"] is False


def test_rule1a_ad_gate_passes_a_genuine_gpd_tail():
    S = _exp_surrogates(seed=7)
    u = float(np.quantile(S, 0.90))
    exc = S[S > u] - u
    xi, beta = G.fit_gpd(exc)
    ad = G.ad_gate(exc, xi, beta, np.random.default_rng(8), n_boot=100)
    assert ad["pass"] is True and ad["p_ad"] >= G.GPD_AD_ALPHA
    # the method must be the parametric bootstrap, not a table lookup
    assert "parametric" in ad["method"] and "refit" in ad["method"]


def test_rule1a_ad_statistic_grows_when_the_fit_is_wrong():
    """A2 is only a gate if it is actually sensitive to misfit."""
    exc = np.random.default_rng(9).standard_exponential(400)
    xi, beta = G.fit_gpd(exc)
    good = G.ad_statistic(exc, xi, beta)
    bad = G.ad_statistic(exc, xi, beta * 4.0)
    assert bad > good * 3


# ---------------------------------------------------------------- rule 1(b) ----
def test_rule1b_xi_stability_passes_on_a_single_gpd():
    S = _exp_surrogates(n=4000, seed=11)
    xs = G.xi_stability(S, np.random.default_rng(12), n_boot=80)
    assert xs["pass"] is True
    assert [r["q"] for r in xs["rows"]] == list(G.GPD_STABILITY_QS)


def test_rule1b_xi_stability_fails_when_the_shape_walks_with_the_threshold():
    """Bulk light, extreme tail heavy: xi at the top 5% cannot agree with the top 20%."""
    rng = np.random.default_rng(13)
    bulk = rng.uniform(0.0, 2.0, 3600)
    heavy = 2.0 + 30.0 * (rng.pareto(0.45, 400) + 1.0)
    S = np.concatenate([bulk, heavy])
    xs = G.xi_stability(S, np.random.default_rng(14), n_boot=80)
    assert xs["pass"] is False


def test_rule1b_gate_blocks_extrapolation_end_to_end():
    rng = np.random.default_rng(15)
    S = np.concatenate([rng.uniform(0.0, 2.0, 3600),
                        2.0 + 30.0 * (rng.pareto(0.45, 400) + 1.0)])
    out = G.gpd_tail_p(S, float(S.max()) * 5, 4000, "block_bootstrap",
                       np.random.default_rng(16), **FAST)
    assert out["p_method"] == G.P_UNRESOLVED
    assert out["reason"] is not None


# ------------------------------------------------------------------ rule 2 -----
def test_rule2_reported_number_is_the_ci_upper_not_the_point_estimate():
    S = _exp_surrogates(n=2000, seed=21)
    s_obs = float(np.quantile(S, 0.90)) + 6.0
    out = G.gpd_tail_p(S, s_obs, 2000, "block_bootstrap",
                       np.random.default_rng(22), **FAST)
    assert out["p_method"] == G.P_GPD
    assert out["p_report"] == out["p_ci_upper"]
    assert out["p_ci_lower"] <= out["p_point"] <= out["p_ci_upper"]
    # the conservative direction: what is reported is LARGER than the point estimate
    assert out["p_report"] > out["p_point"]


def test_rule2_price_row_enters_bh_at_the_ci_upper():
    S = _exp_surrogates(n=2000, seed=23)
    s_obs = float(np.quantile(S, 0.90)) + 6.0
    pr = G.price_row(1.0 / 2001, 2000, s_obs, S, "block_bootstrap",
                     np.random.default_rng(24), enabled=True, **FAST)
    assert pr["p_method"] == G.P_GPD
    assert pr["p_bh"] == pr["gpd"]["p_ci_upper"]
    assert pr["p_bh"] != pr["gpd"]["p_point"]


# ------------------------------------------------------------------ rule 3 -----
def test_rule3_floor_is_one_decade_below_the_mc_floor():
    assert G.gpd_floor(2000) == pytest.approx(0.1 / 2001.0)
    assert G.gpd_floor(9999) == pytest.approx(0.1 / 10000.0)


def test_rule3_below_the_floor_is_unresolved_never_very_significant():
    S = _exp_surrogates(n=2000, seed=31)
    # 25 e-folds past the threshold: the extrapolation is many decades out
    s_obs = float(np.quantile(S, 0.90)) + 25.0
    out = G.gpd_tail_p(S, s_obs, 2000, "block_bootstrap",
                       np.random.default_rng(32), **FAST)
    assert out["p_method"] == G.P_UNRESOLVED
    assert "one-decade" in out["reason"]
    assert out["p_report"] == pytest.approx(1.0 / 2001.0)
    # and the row is INELIGIBLE for BH rejection (§P6-2(1))
    pr = G.price_row(1.0 / 2001, 2000, s_obs, S, "block_bootstrap",
                     np.random.default_rng(33), enabled=True, **FAST)
    assert pr["p_method"] == G.P_UNRESOLVED and pr["bh_eligible"] is False


# ------------------------------------------------------------------ rule 4 -----
@pytest.mark.parametrize("nt", ["all_shifts", "circular_shift", "circular-shift",
                                "enumeration"])
def test_rule4_all_shifts_null_is_structurally_forbidden(nt):
    with pytest.raises(G.ForbiddenNullError) as e:
        G.assert_null_permitted(nt)
    assert "FORBIDDEN" in str(e.value)


def test_rule4_gpd_tail_p_raises_on_the_forbidden_null():
    S = _exp_surrogates(seed=41)
    with pytest.raises(G.ForbiddenNullError):
        G.gpd_tail_p(S, 10.0, 2000, "all_shifts", np.random.default_rng(42))


def test_rule4_undeclared_and_missing_null_types_are_refused_too():
    with pytest.raises(G.ForbiddenNullError):
        G.assert_null_permitted("something_new")
    with pytest.raises(G.ForbiddenNullError):
        G.assert_null_permitted(None)


@pytest.mark.parametrize("nt", ["block_bootstrap", "ar1", "permutation",
                                "etas_sim"])
def test_rule4_permitted_substrates_are_accepted(nt):
    assert G.assert_null_permitted(nt) == nt


# ------------------------------------------------------------------ rule 5 -----
def test_rule5_census_always_has_the_three_numbers():
    c = G.census([{"p_method": G.P_MC_RESOLVED}, {"p_method": G.P_GPD},
                  {"p_method": G.P_GPD}])
    assert set(G.P_METHODS) <= set(c)
    assert c[G.P_GPD] == 2 and c[G.P_UNRESOLVED] == 0
    line = G.census_line(c, 3)
    for m in G.P_METHODS:
        assert m in line


def test_rule5_label_is_emitted_even_when_extrapolation_is_off():
    pr = G.price_row(1.0 / 201, 200, 5.0, None, "block_bootstrap", None,
                     enabled=False)
    assert pr["p_method"] == G.P_UNRESOLVED       # censored, and honestly so
    assert pr["bh_eligible"] is True              # but not disqualified
    pr2 = G.price_row(0.3, 200, 1.0, None, "block_bootstrap", None, enabled=False)
    assert pr2["p_method"] == G.P_MC_RESOLVED


def test_rule5_a_binding_enumeration_floor_is_labelled_unresolved_but_eligible():
    """§P6-2(4) forbids extrapolating off the shift null, so a p censored THERE is
    honestly UNRESOLVED -- and nothing failed a gate, so it stays BH-eligible."""
    # the case is real whenever there are FEWER admissible shifts than surrogates:
    # 7,641 shifts (floor 1.31e-4) against N_max = 10,000 (floor 1.00e-4), so the
    # shift null is the binding component AND is sitting on its own floor
    pr = G.price_row(1 / 10001.0, 10000, 3.0, None, "block_bootstrap", None,
                     other_p=1 / 7642.0, other_floor=1 / 7642.0, enabled=False)
    assert pr["p_method"] == G.P_UNRESOLVED
    assert pr["bh_eligible"] is True
    assert "record length" in pr["p_method_reason"].lower()
    assert pr["p_bh"] == pytest.approx(1 / 7642.0)
    # a shift p ABOVE its floor is resolved, as before
    pr2 = G.price_row(1 / 10001.0, 10000, 3.0, None, "block_bootstrap", None,
                      other_p=0.2, other_floor=1 / 7642.0, enabled=False)
    assert pr2["p_method"] == G.P_MC_RESOLVED


def test_rule5_label_propagates_to_stubs_json(tmp_path):
    tests = [
        {"feature": "f_mc", "family": 1, "kind": "linear", "lag": 0,
         "test": "glm_poisson_offset_etas", "beta": [0.2], "se": [0.01],
         "amplitude_log_rate": 0.2, "pct_rate_modulation": 22.0,
         "p_raw": 1e-4, "p_bh": 1e-4, "bh_q": 0.01, "passes_fdr": True,
         "p_method": G.P_MC_RESOLVED, "stratum": "S1", "n_surrogates": 2000},
        {"feature": "f_gpd", "family": 1, "kind": "linear", "lag": 3,
         "test": "glm_poisson_offset_etas", "beta": [0.3], "se": [0.01],
         "amplitude_log_rate": 0.3, "pct_rate_modulation": 35.0,
         "p_raw": 5e-4, "p_bh": 8e-5, "bh_q": 0.02, "passes_fdr": True,
         "p_method": G.P_GPD, "stratum": "S1", "n_surrogates": 2000},
        {"feature": "f_dead", "family": 2, "kind": "linear", "lag": 0,
         "test": "glm_poisson_offset_etas", "beta": [0.0], "se": [0.01],
         "amplitude_log_rate": 0.0, "pct_rate_modulation": 0.0,
         "p_raw": 0.5, "p_bh": 0.5, "bh_q": 0.9, "passes_fdr": False,
         "p_method": G.P_UNRESOLVED, "stratum": "S2", "n_surrogates": 2000},
    ]
    path = ms.write_stubs(str(tmp_path), {"preset": "test"}, tests)
    doc = json.load(open(path, encoding="utf-8"))
    # rule 5: every emitted stub row carries the label
    assert all("p_method" in s for s in doc["stubs"])
    assert doc["p_method_census"][G.P_GPD] == 1
    assert G.P_MC_RESOLVED in doc["p_method_census_line"]


# ------------------------------------------------------------------ rule 6 -----
def test_rule6_brute_force_n_is_at_least_ten_over_p():
    assert G.brute_force_n(1e-4) == 100000
    assert G.brute_force_n(0.003) == math.ceil(10 / 0.003)
    assert G.brute_force_n(1e-4) * 1e-4 >= 10.0
    with pytest.raises(ValueError):
        G.brute_force_n(0.0)


def test_rule6_a_gpd_survivor_is_a_candidate_and_never_a_stub(tmp_path):
    tests = [
        {"feature": "f_gpd", "family": 1, "kind": "linear", "lag": 3,
         "test": "glm_poisson_offset_etas", "beta": [0.3], "se": [0.01],
         "amplitude_log_rate": 0.3, "pct_rate_modulation": 35.0,
         "p_raw": 5e-4, "p_bh": 8e-5, "bh_q": 0.02, "passes_fdr": True,
         "p_method": G.P_GPD, "stratum": "S1", "n_surrogates": 2000},
    ]
    doc = json.load(open(ms.write_stubs(str(tmp_path), {"preset": "t"}, tests),
                         encoding="utf-8"))
    assert doc["n_stubs"] == 0
    assert doc["n_gpd_candidates"] == 1
    cand = doc["gpd_candidates_requiring_brute_force"][0]
    assert cand["status"] == G.CANDIDATE_LABEL == "CANDIDATE-REQUIRES-BRUTE-FORCE"


def test_rule6_a_confirmed_candidate_becomes_a_stub(tmp_path):
    """Relabelled to MC_RESOLVED by the confirmation, it is a stub again -- and the
    brute-force p is the number carried."""
    tests = [
        {"feature": "f_gpd", "family": 1, "kind": "linear", "lag": 3,
         "test": "glm_poisson_offset_etas", "beta": [0.3], "se": [0.01],
         "amplitude_log_rate": 0.3, "pct_rate_modulation": 35.0,
         "p_raw": 5e-4, "p_bh": 8e-5, "bh_q": 0.02, "passes_fdr": True,
         "p_method": G.P_MC_RESOLVED, "p_brute_force": 6.1e-5,
         "stratum": "S1", "n_surrogates": 2000},
    ]
    doc = json.load(open(ms.write_stubs(str(tmp_path), {"preset": "t"}, tests),
                         encoding="utf-8"))
    assert doc["n_stubs"] == 1 and doc["n_gpd_candidates"] == 0
    assert doc["stubs"][0]["p_brute_force"] == 6.1e-5
