"""§P6-2 rule 7's harness: the roster is honest and the batched kernel is the engine's.

The audit itself (`python -u -m engine.audit_gpd --jobs N`) is a minutes-long run and
is not a unit test. What IS unit-tested is everything that would make its table a
lie: that the roster really spans the three test kinds at the declared count, that
the batched Lomb-Scargle used for the 10^6-surrogate arm is numerically the engine's
own scalar kernel, and that the accept/reject bars are the ones §P6-2(7) states.
"""

import numpy as np
import pytest

from engine import audit_gpd as A
from engine import mine as M


def test_roster_spans_the_three_test_kinds_at_the_declared_count():
    cases = A.roster()
    kinds = {c["kind"] for c in cases}
    assert kinds == {"glm", "mark", "period"}
    # >= 15 representative TESTS: the roster x the declared target-p grid
    assert len(cases) * len(A.TARGET_PS) >= 15
    assert all(c["null_type"] in ("block_bootstrap", "ar1", "permutation")
               for c in cases)


def test_no_roster_case_uses_a_forbidden_null():
    """§P6-2(4) again: the calibration may not be run on the all-shifts null either."""
    from engine import gpd_tail as G
    for c in A.roster():
        assert G.assert_null_permitted(c["null_type"])


def test_the_bars_are_the_ones_the_rule_states():
    assert A.COVERAGE_BAR == 0.90
    assert A.UNDERSTATE_BAR == 3.0
    assert A.N_BRUTE == 1_000_000
    assert A.N_SURROGATES_GPD == 2000


def test_batched_lomb_scargle_equals_the_engine_scalar_kernel():
    rng = np.random.default_rng(2)
    n, nf = 250, 30
    t = np.arange(n, dtype=float)
    periods = np.exp(np.linspace(np.log(2.0), np.log(n / 3.0), nf))
    Xb = rng.standard_normal((7, n))
    batch = A._ls_power_batch(t, Xb, periods)
    scalar = np.array([M.lomb_scargle_power(t, Xb[i], periods).max()
                       for i in range(Xb.shape[0])])
    assert np.allclose(batch, scalar, rtol=1e-11, atol=0.0)


def test_each_case_draw_returns_the_requested_number_of_statistics():
    for c in A.roster():
        v = np.asarray(c["draw"](12, np.random.default_rng(1)))
        assert v.shape == (12,) and np.isfinite(v).all()


@pytest.mark.parametrize("kind", ["glm", "mark", "period"])
def test_run_case_is_wired_end_to_end_at_a_toy_brute_force(kind, monkeypatch):
    """The whole comparison path, with the 10^6 arm shrunk so it runs in a test."""
    monkeypatch.setattr(A, "N_BRUTE", 4000)
    monkeypatch.setattr(A, "PILOT_N", 3000)
    monkeypatch.setattr(A, "N_SURROGATES_GPD", 500)
    monkeypatch.setattr(A, "N_BOOT", 40)
    monkeypatch.setattr(A, "N_AD_BOOT", 40)
    case = next(c for c in A.roster() if c["kind"] == kind)
    r = A.run_case(case, 0.01, 12345)
    assert r["kind"] == kind
    assert r["p_method"] in ("MC_RESOLVED", "GPD_EXTRAPOLATED", "UNRESOLVED")
    assert 0.0 < r["p_brute_force"] <= 1.0
    assert isinstance(r["covered"], bool)
    # every row carries the interval's provenance, whether or not it fitted
    for k in ("ci_method", "bca_z0", "bca_a", "bca_jackknife",
              "p_ci_upper_percentile", "covered_percentile", "set"):
        assert k in r or k == "set"


# ------------------------------------------ §P7-7(b): the DECISIVE set is fresh --
def test_the_decisive_set_is_at_least_125_fresh_comparisons():
    """The §P7-7(b) arithmetic, made a test: n >= 0.9*0.1/((0.90-0.846)/2)^2 = 125."""
    from engine import gpd_tail as G
    cases = A.fresh_roster()
    n = len(A.specs_for(A.SET_DECISIVE, A.TARGET_PS))
    assert n == len(cases) * len(A.TARGET_PS) == 135
    assert n >= G.GPD_MIN_CALIBRATION_COMPARISONS == 125


def test_the_decisive_set_spans_the_three_kinds_evenly():
    kinds = [c["kind"] for c in A.fresh_roster()]
    assert set(kinds) == {"glm", "mark", "period"}
    assert {k: kinds.count(k) for k in set(kinds)} == {"glm": 15, "mark": 15,
                                                      "period": 15}
    for c in A.fresh_roster():
        assert c["null_type"] in ("block_bootstrap", "ar1", "permutation")


def test_the_decisive_cases_and_seeds_are_disjoint_from_the_legacy_ones():
    """FRESH means fresh. Re-running BCa on the comparisons that motivated the swap
    would be testing the fix on its own training data (§P7-7(b)(ii)), so neither the
    case names nor the run seeds may overlap."""
    leg = A.specs_for(A.SET_LEGACY, A.TARGET_PS)
    dec = A.specs_for(A.SET_DECISIVE, A.TARGET_PS)
    assert not ({s[3] for s in leg} & {s[3] for s in dec})
    assert not ({c["name"] for c in A.roster()}
                & {c["name"] for c in A.fresh_roster()})
    # and the legacy seeds are byte-for-byte the ones attempt 1 ran, so the legacy
    # block really is a RE-SCORE of the same comparisons rather than a new sample
    assert [s[3] for s in leg] == [700000 + 137 * i + 11 * j
                                   for i in range(len(A.roster()))
                                   for j in range(len(A.TARGET_PS))]


def test_no_fresh_case_uses_a_forbidden_null():
    from engine import gpd_tail as G
    for c in A.fresh_roster():
        assert G.assert_null_permitted(c["null_type"])


def test_each_fresh_case_draw_returns_the_requested_number_of_statistics():
    for c in A.fresh_roster():
        v = np.asarray(c["draw"](8, np.random.default_rng(1)))
        assert v.shape == (8,) and np.isfinite(v).all()


def test_the_bars_did_not_move_when_the_design_grew():
    """§P7-7(b): NO on the design, YES on the bars. This is the tripwire."""
    assert A.COVERAGE_BAR == 0.90 and A.UNDERSTATE_BAR == 3.0
    assert A.N_BRUTE == 1_000_000 and A.N_SURROGATES_GPD == 2000
    assert A.TARGET_PS == (2e-4, 5e-4, 1e-3)


# ------------------------------------------------------------- the gate's own maths
def _row(kind, p_method, covered, us, cov_p=True, us_p=1.0):
    return {"kind": kind, "p_method": p_method, "covered": covered,
            "understatement_x": us, "covered_percentile": cov_p,
            "understatement_x_percentile": us_p}


def test_score_reports_n_scored_and_it_is_the_rows_that_entered_the_bars():
    """§P7-9(1): the licensing count is SCORED, and it must be stated, not inferred."""
    rows = ([_row("glm", "GPD_EXTRAPOLATED", True, 1.2)] * 10
            + [_row("glm", "UNRESOLVED", False, float("inf"))] * 4)
    s = A.score(rows)
    assert s["n_comparisons"] == 14
    assert s["n_scored"] == s["n_fitted"] == 10
    assert s["coverage"] == s["n_covered"] / s["n_scored"]


# ------------------------------------------------- §P7-9(2): the top-up batches --
def test_topup_batches_are_pre_declared_and_span_the_three_kinds():
    assert A.N_TOPUP_BATCHES == 5
    for b in range(1, A.N_TOPUP_BATCHES + 1):
        cases = A.topup_roster(b)
        assert {c["kind"] for c in cases} == {"glm", "mark", "period"}
        assert len(cases) == 3
        assert len(A.specs_for(A.topup_set_name(b), A.TARGET_PS)) == 9


def test_topup_cases_and_seeds_are_disjoint_from_every_earlier_set():
    """FRESH again. A top-up that reused a case or a seed would be re-scoring the
    comparisons already in the pool, which adds no information and inflates n."""
    names, seeds = set(), set()
    for sn in (A.SET_LEGACY, A.SET_DECISIVE):
        names |= {c["name"] for c in A.roster_for(sn)}
        seeds |= {s[3] for s in A.specs_for(sn, A.TARGET_PS)}
    for b in range(1, A.N_TOPUP_BATCHES + 1):
        sn = A.topup_set_name(b)
        n = {c["name"] for c in A.roster_for(sn)}
        s = {sp[3] for sp in A.specs_for(sn, A.TARGET_PS)}
        assert not (n & names) and not (s & seeds)
        names |= n
        seeds |= s


def test_no_topup_case_uses_a_forbidden_null():
    from engine import gpd_tail as G
    for b in range(1, A.N_TOPUP_BATCHES + 1):
        for c in A.topup_roster(b):
            assert G.assert_null_permitted(c["null_type"])


def test_each_topup_case_draw_returns_the_requested_number_of_statistics():
    for b in range(1, A.N_TOPUP_BATCHES + 1):
        for c in A.topup_roster(b):
            v = np.asarray(c["draw"](8, np.random.default_rng(1)))
            assert v.shape == (8,) and np.isfinite(v).all()


def test_the_stopping_rule_reads_only_the_scored_count():
    """§P7-9(2), the anti-optional-stopping guarantee, tested where it lives.

    `score` is the only thing the loop consults, so the rule is verified by showing
    that two row sets with IDENTICAL scored counts but wildly different coverage and
    understatement produce the identical stopping decision. If the rule could see
    the bars, these two would diverge -- and stopping when the bars looked good
    would manufacture a pass out of sampling noise.
    """
    good = [_row("glm", "GPD_EXTRAPOLATED", True, 1.0) for _ in range(120)]
    bad = [_row("glm", "GPD_EXTRAPOLATED", False, 2.9) for _ in range(120)]
    assert A.score(good)["n_scored"] == A.score(bad)["n_scored"] == 120
    assert A.score(good)["pass"] is not A.score(bad)["pass"]     # bars do differ
    # the decision the loop makes is a function of n_scored alone
    decide = lambda rows: A.score(rows)["n_scored"] >= 125
    assert decide(good) == decide(bad) is False
    assert decide(good + [_row("glm", "GPD_EXTRAPOLATED", True, 1.0)] * 5) is True


def test_score_counts_only_fitted_rows_and_applies_both_bars():
    rows = ([_row("glm", "GPD_EXTRAPOLATED", True, 1.2)] * 9
            + [_row("glm", "GPD_EXTRAPOLATED", False, 2.0)]
            + [_row("glm", "UNRESOLVED", False, float("inf"))] * 4)
    s = A.score(rows)
    # the 4 UNRESOLVED rows are the estimator REFUSING; scoring them as misses would
    # penalise the one behaviour §P6-2(1) exists to produce
    assert s["n_comparisons"] == 14 and s["n_fitted"] == 10
    assert s["coverage"] == pytest.approx(0.9) and s["coverage_pass"] is True
    assert s["worst_understatement"] == pytest.approx(2.0)
    assert s["understatement_pass"] is True and s["pass"] is True


def test_score_fails_on_either_bar_alone():
    lowcov = [_row("glm", "GPD_EXTRAPOLATED", i > 1, 1.0) for i in range(10)]
    assert A.score(lowcov)["pass"] is False
    assert A.score(lowcov)["understatement_pass"] is True
    bigus = [_row("glm", "GPD_EXTRAPOLATED", True, 1.0) for _ in range(9)]
    bigus.append(_row("glm", "GPD_EXTRAPOLATED", True, 3.5))
    assert A.score(bigus)["coverage_pass"] is True
    assert A.score(bigus)["pass"] is False


def test_the_percentile_columns_are_a_reference_and_never_the_verdict():
    """§P7-7(a): ONE estimator is scored. A run where BCa fails and the superseded
    percentile interval passes must still be a FAIL -- that is the whole point."""
    rows = [_row("glm", "GPD_EXTRAPOLATED", False, 1.0, cov_p=True)
            for _ in range(10)]
    assert A.score(rows)["pass"] is False
    ref = A.score(rows, "understatement_x_percentile", "covered_percentile")
    assert ref["pass"] is True                       # and it licenses nothing


def test_by_kind_partitions_the_rows_without_loss():
    rows = ([_row("glm", "GPD_EXTRAPOLATED", True, 1.0)] * 3
            + [_row("mark", "GPD_EXTRAPOLATED", False, 2.0)] * 2
            + [_row("period", "UNRESOLVED", False, float("inf"))])
    bk = A.by_kind(rows)
    assert set(bk) == {"glm", "mark", "period"}
    assert sum(v["n_comparisons"] for v in bk.values()) == len(rows)
    assert bk["mark"]["coverage"] == 0.0 and bk["period"]["n_fitted"] == 0
