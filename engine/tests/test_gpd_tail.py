"""§P6-2, rule by rule. Every numbered rule that can be unit-tested, is.

  rule 1(a) AD gate rejects a non-GPD tail      test_rule1a_*
  rule 1(b) xi-stability gate                   test_rule1b_*
  rule 2    CI-upper, never the point estimate  test_rule2_*
  rule 3    one-decade floor                    test_rule3_*
  rule 4    all-shifts prohibition (structural) test_rule4_*
  rule 5    label on every row + census         test_rule5_*
  rule 6    candidate, not stub; N >= 10/p      test_rule6_*
  rule 7    the harness exists and is honest    engine/tests/test_gpd_audit.py

and the §P7-7 BCa interval that rule 2 now reports:

  z0 / a arithmetic on a known synthetic          test_bca_z0_and_acceleration_*
  the degenerate cases, each one named            test_bca_degenerate_*
  BCa is THE estimator, percentile is a reference test_bca_replaces_*
  one recipe for point, bootstrap and jackknife   test_bca_jackknife_*
  §P7-7(b) the 125-comparison licensing floor     test_rule7_*
"""

import json
import math
import os

import numpy as np
import pytest
from scipy import stats

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


# --------------------------------------------------- §P7-7(a): the BCa interval --
def _bca_by_hand(theta_hat, boot, jack, level=0.95):
    """The BCa endpoints written out again, independently, from the formulas.

    Deliberately a second implementation rather than a call into the module: a test
    that reuses the code under test only proves the code equals itself.
    """
    boot = np.asarray(boot, float)
    B = boot.size
    frac = (np.count_nonzero(boot < theta_hat)
            + 0.5 * np.count_nonzero(boot == theta_hat)) / B
    z0 = stats.norm.ppf(frac)
    jack = np.asarray(jack, float)
    d = jack.mean() - jack
    a = float(np.sum(d ** 3) / (6.0 * np.sum(d * d) ** 1.5))
    z_lo = stats.norm.ppf((1 - level) / 2)
    z_hi = stats.norm.ppf(1 - (1 - level) / 2)
    al = float(stats.norm.cdf(z0 + (z0 + z_lo) / (1 - a * (z0 + z_lo))))
    ah = float(stats.norm.cdf(z0 + (z0 + z_hi) / (1 - a * (z0 + z_hi))))
    return float(z0), a, al, ah


def test_bca_z0_and_acceleration_match_the_formulas_on_a_known_synthetic():
    """A synthetic where every input is known by construction, ties included.

    boot = {1/400, ..., 400/400} and theta_hat = 0.25, so exactly 99 replicates lie
    strictly below and exactly one is an exact tie. The tie is SPLIT, giving
    frac = 99.5/400 = 0.24875 -- a pile-up at the point estimate must not be counted
    wholly on one side, which is the difference between a defensible z0 and one that
    drifts with the discreteness of the resample.
    """
    boot = np.arange(1, 4001, dtype=float) / 4000.0
    jack = np.array([1.0, 2.0, 3.0, 4.0, -6.0])      # deliberately skewed
    theta_hat = 0.25                                 # = 1000/4000, an exact tie
    assert np.count_nonzero(boot < theta_hat) == 999
    assert np.count_nonzero(boot == theta_hat) == 1
    z0, a, al, ah = _bca_by_hand(theta_hat, boot, jack)
    assert z0 == pytest.approx(stats.norm.ppf(999.5 / 4000.0))
    assert a > 0
    # both adjusted levels land INSIDE the bootstrap's own support, so this case
    # tests the formula and nothing but the formula
    assert 0.5 / 4000 < al < ah < 1 - 0.5 / 4000

    got = G.bca_interval(theta_hat, boot, jack, level=0.95)
    assert got["fallback"] is False and got["notes"] == []
    assert got["z0"] == pytest.approx(z0, rel=1e-12)
    assert got["a"] == pytest.approx(a, rel=1e-12)
    assert got["alpha_lo"] == pytest.approx(al, rel=1e-9)
    assert got["alpha_hi"] == pytest.approx(ah, rel=1e-9)
    assert got["lo"] == pytest.approx(float(np.percentile(boot, 100 * al)))
    assert got["hi"] == pytest.approx(float(np.percentile(boot, 100 * ah)))
    assert got["n_boot"] == 4000 and got["n_jack"] == 5


def test_bca_adjusted_levels_are_clamped_to_the_bootstrap_support():
    """The ECDF holds nothing past its own extremes, so an alpha outside
    [1/2B, 1-1/2B] is clamped -- and the clamp is announced, not silent."""
    boot = np.arange(1, 401, dtype=float) / 400.0
    jack = np.array([1.0, 2.0, 3.0, 4.0, 10.0])
    _z0, _a, al, _ah = _bca_by_hand(0.25, boot, jack)
    assert al < 0.5 / 400                            # the raw level is off the end
    got = G.bca_interval(0.25, boot, jack)
    assert got["alpha_lo"] == pytest.approx(0.5 / 400)
    assert any("clamped" in n for n in got["notes"])
    assert got["fallback"] is False                  # clamped is still BCa


def test_bca_acceleration_is_bounded_by_one_sixth_so_the_map_rarely_degenerates():
    """Cauchy-Schwarz on the jackknife influence values: |sum d^3| <= (sum d^2)^1.5,
    so |a| <= 1/6 for ANY input. That is why the degenerate-denominator branch is a
    guard rather than a routine outcome -- 1 - a(z0+z) can only go non-positive if
    |z0 + z| > 6, which the z0 clamp keeps out of reach at any realistic B."""
    rng = np.random.default_rng(99)
    for _ in range(200):
        jack = rng.standard_normal(rng.integers(3, 60)) * rng.uniform(0.1, 100)
        jack[rng.integers(0, jack.size)] *= rng.uniform(1, 500)   # force outliers
        got = G.bca_interval(0.5, np.arange(1, 401, dtype=float) / 400.0, jack)
        assert abs(got["a"]) <= 1.0 / 6.0 + 1e-12


def test_bca_equals_percentile_when_there_is_no_bias_and_no_skew():
    """The sanity anchor: z0 = 0 and a = 0 must reduce BCa to the percentile CI."""
    boot = np.arange(1, 401, dtype=float)
    theta_hat = 200.5                                # exactly the median: z0 = 0
    jack = np.array([1.0, 2.0, 3.0, 4.0, 5.0])       # symmetric: a = 0
    got = G.bca_interval(theta_hat, boot, jack)
    assert got["z0"] == pytest.approx(0.0, abs=1e-12)
    assert got["a"] == pytest.approx(0.0, abs=1e-12)
    assert got["lo"] == pytest.approx(got["percentile"][0])
    assert got["hi"] == pytest.approx(got["percentile"][1])
    assert got["fallback"] is False


def test_bca_degenerate_all_bootstrap_below_the_point_estimate_is_clamped():
    """z0 = Phi^-1(1) = +inf. Clamped to 1 - 1/(2B), flagged, and still a BCa row."""
    boot = np.arange(1, 401, dtype=float)
    jack = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    got = G.bca_interval(1e9, boot, jack)            # every replicate is below
    assert np.isfinite(got["z0"]) and got["z0"] > 2.0
    assert got["z0"] == pytest.approx(stats.norm.ppf(1 - 0.5 / 400))
    assert any("one side of the point estimate" in n for n in got["notes"])
    assert got["fallback"] is False                  # clamped, not abandoned
    # a huge positive z0 pushes BOTH ends up: the conservative direction for a p
    assert got["hi"] >= got["percentile"][1]


def test_bca_degenerate_all_bootstrap_above_the_point_estimate_is_clamped():
    boot = np.arange(1, 401, dtype=float)
    got = G.bca_interval(-1e9, boot, np.array([1.0, 2.0, 3.0, 4.0, 5.0]))
    assert got["z0"] == pytest.approx(stats.norm.ppf(0.5 / 400))
    assert any("one side of the point estimate" in n for n in got["notes"])
    assert got["lo"] <= got["percentile"][0]


def test_bca_degenerate_zero_acceleration_when_the_jackknife_has_no_spread():
    """sum d^2 = 0 is a division by zero, not a small number. a := 0, i.e. plain BC."""
    boot = np.arange(1, 401, dtype=float)
    got = G.bca_interval(150.0, boot, np.full(50, 7.0))
    assert got["a"] == 0.0
    assert any("identical" in n for n in got["notes"])
    assert got["fallback"] is False
    z0 = got["z0"]
    z_hi = stats.norm.ppf(0.975)
    assert got["alpha_hi"] == pytest.approx(float(stats.norm.cdf(2 * z0 + z_hi)))


def test_bca_degenerate_too_few_jackknife_replicates_sets_acceleration_to_zero():
    got = G.bca_interval(150.0, np.arange(1, 401, dtype=float), np.array([3.0, 4.0]))
    assert got["a"] == 0.0 and got["n_jack"] == 2
    assert any("(< 3)" in n for n in got["notes"])


def test_bca_degenerate_map_falls_back_to_percentile_and_says_so():
    """1 - a(z0 + z) <= 0 makes the BCa map undefined; that row must be FLAGGED.

    A fallback that looked like a BCa result would let a degenerate arithmetic case
    masquerade as the declared estimator, which is exactly the bookkeeping §P7-7(a)
    forbids. Since |a| <= 1/6 always, the branch needs |z0 + z| > 6, which is reached
    here by asking for an absurdly extreme LEVEL rather than by a fake acceleration.
    """
    boot = np.arange(1, 401, dtype=float)
    jack = np.concatenate([np.zeros(40), [1.0]])     # a ~ -0.16, near the bound
    level = float(1.0 - 2.0 * stats.norm.sf(7.0))    # z = +-7
    got = G.bca_interval(200.5, boot, jack, level=level)
    assert got["a"] < -0.15
    assert got["fallback"] is True
    assert (got["lo"], got["hi"]) == tuple(got["percentile"])
    assert any("degenerate" in n for n in got["notes"])


def test_bca_degenerate_too_few_bootstrap_replicates_yields_no_interval():
    got = G.bca_interval(1.0, np.arange(5, dtype=float), np.arange(5, dtype=float))
    assert got["lo"] is None and got["hi"] is None
    assert any("< 20" in n for n in got["notes"])


def test_bca_replaces_the_percentile_interval_everywhere_it_is_reported():
    """§P7-7(a): ONE estimator. BCa is what p_report, p_ci_* and BH see."""
    S = _exp_surrogates(n=2000, seed=21)
    s_obs = float(np.quantile(S, 0.90)) + 6.0
    out = G.gpd_tail_p(S, s_obs, 2000, "block_bootstrap",
                       np.random.default_rng(22), **FAST)
    assert out["p_method"] == G.P_GPD
    b = out["bca"]
    assert b["method"].startswith("BCa")
    assert out["p_ci_lower"] == b["lo"] and out["p_ci_upper"] == b["hi"]
    assert out["p_report"] == out["p_ci_upper"]
    # the percentile endpoints exist ONLY as a labelled reference, under their own
    # names, and they are NOT what is reported
    assert out["p_ci_upper_percentile"] == b["percentile"][1]
    assert out["p_ci_lower_percentile"] == b["percentile"][0]
    assert out["p_report"] != out["p_ci_upper_percentile"]
    c = G.compact(out)
    assert c["ci_method"] == "BCa"
    assert c["p_ci_upper"] == b["hi"]


def test_bca_there_is_no_switch_that_selects_the_percentile_interval():
    """The forking path is closed structurally: no argument chooses the estimator.

    §P7-7(a): what would make the swap inadmissible is running BCa and percentile
    and reporting whichever passed. So there must be no `ci_method=` / `use_bca=`
    knob anywhere on the public estimator surface -- not defaulted to BCa, ABSENT.
    """
    import inspect
    for fn in (G.gpd_tail_p, G.price_row, G.compact):
        names = set(inspect.signature(fn).parameters)
        assert not (names & {"ci_method", "use_bca", "bca", "interval",
                             "ci_type", "percentile"}), fn.__name__


def test_bca_jackknife_is_over_the_surrogate_set_and_re_derives_everything():
    """The estimand is the extrapolated p, so a jackknife replicate re-derives
    u, zeta, xi and beta -- it may NOT condition on the full-sample values."""
    S = _exp_surrogates(n=300, seed=51)
    s_obs = float(np.quantile(S, 0.90)) + 3.0
    vals, mode = G._jackknife_ps(S, s_obs)
    assert mode == "leave-one-out" and vals.size == 300
    # Deleting a point moves u and zeta, not only (xi, beta). Both checks matter:
    # every deletion moves u OFF the full-sample threshold, and deleting a point
    # from inside the tail moves it somewhere different from deleting one from the
    # bulk -- which is precisely the threshold uncertainty a conditional jackknife
    # would have thrown away.
    q = G.GPD_THRESHOLD_Q
    full = G._tail_p_from_sample(S, s_obs, q)
    u_bulk = G._tail_p_from_sample(np.delete(S, int(np.argmin(S))), s_obs, q)["u"]
    u_tail = G._tail_p_from_sample(np.delete(S, int(np.argmax(S))), s_obs, q)["u"]
    assert u_bulk != pytest.approx(full["u"])
    assert u_tail != pytest.approx(u_bulk)
    # zeta is re-derived from each replicate too, though at an EMPIRICAL quantile it
    # is pinned near 1 - q by construction and barely moves; the threshold
    # uncertainty this jackknife has to carry travels through u and the refit that
    # follows it, which is what the assertions above check.
    assert full["zeta"] == pytest.approx(1.0 - q, abs=0.01)
    # the replicates genuinely differ -- an estimator conditioned on the full-sample
    # u and (xi, beta) would produce a constant here and an acceleration of exactly 0
    assert vals.std() > 0
    assert len(set(np.round(vals, 15))) > 10


def test_bca_jackknife_switches_to_grouped_above_the_declared_cap():
    """Above the cap the exact jackknife's linear cost would dominate; the grouped
    delete-d jackknife is the declared substitute, and it names itself on the row."""
    S = np.random.default_rng(52).standard_exponential(G.GPD_JACKKNIFE_MAX_POINTS + 50)
    s_obs = float(np.quantile(S, 0.90)) + 3.0
    vals, mode = G._jackknife_ps(S, s_obs, groups=40)
    assert mode == "grouped-40" and vals.size <= 40 and vals.size >= 30
    assert G.GPD_JACKKNIFE_GROUPS > 0 and G.GPD_JACKKNIFE_MAX_POINTS >= 2000


def test_bca_tiny_exceedance_counts_cannot_produce_an_interval():
    """The min-exceedance gate fires before any of this, and the helper refuses too."""
    S = _exp_surrogates(n=100, seed=53)             # top 10% = 10 exceedances
    s_obs = float(S.max()) + 1.0
    out = G.gpd_tail_p(S, s_obs, 100, "block_bootstrap",
                       np.random.default_rng(54), **FAST)
    assert out["p_method"] == G.P_UNRESOLVED
    assert "exceedances" in out["reason"]
    assert out["p_ci_upper"] is None and out.get("bca") is None
    # and the shared recipe refuses on its own terms
    assert G._tail_p_from_sample(S[:5], s_obs, 0.90, min_exceedances=25) is None
    assert G._tail_p_from_sample(S, float(np.quantile(S, 0.5)), 0.90) is None


def test_bca_one_recipe_is_shared_by_point_bootstrap_and_jackknife():
    """The point estimate is literally the helper's output on the full sample."""
    S = _exp_surrogates(n=2000, seed=21)
    s_obs = float(np.quantile(S, 0.90)) + 6.0
    out = G.gpd_tail_p(S, s_obs, 2000, "block_bootstrap",
                       np.random.default_rng(22), **FAST)
    assert out["p_method"] == G.P_GPD
    direct = G._tail_p_from_sample(S, s_obs, G.GPD_THRESHOLD_Q)
    assert out["p_point"] == pytest.approx(direct["p"], rel=1e-12)
    assert out["u"] == pytest.approx(direct["u"], rel=1e-12)
    assert out["zeta"] == pytest.approx(direct["zeta"], rel=1e-12)


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


# ------------------------------------------------------------------ rule 7 -----
def _artifact(tmp_path, **kw):
    doc = {"rule": "§P6-2(7)", "ci_method": "BCa", "pass": True,
           "coverage": 0.94, "coverage_bar": 0.90, "worst_understatement": 1.8,
           "understatement_bar": 3.0, "n_comparisons": 144, "n_fitted": 131,
           "n_scored": 131}
    doc.update(kw)
    p = tmp_path / "cal.json"
    p.write_text(json.dumps(doc), encoding="utf-8")
    return str(p)


def test_rule7_licence_defaults_to_the_bca_artifact_not_the_superseded_one():
    """§P7-7(b): the legacy percentile artifact may not be the licence, ever."""
    assert G.CALIBRATION_PATH.replace("\\", "/").endswith("audit_gpd_bca.json")
    assert G.LEGACY_CALIBRATION_PATH.replace("\\", "/").endswith("audit_gpd.json")
    assert G.CALIBRATION_PATH != G.LEGACY_CALIBRATION_PATH


def test_rule7_missing_artifact_refuses_to_run(tmp_path):
    with pytest.raises(G.NotCalibratedError) as e:
        G.assert_calibrated(str(tmp_path / "nope.json"))
    assert "no calibration artifact" in str(e.value)


def test_rule7_a_reject_verdict_refuses_to_run(tmp_path):
    with pytest.raises(G.NotCalibratedError) as e:
        G.assert_calibrated(_artifact(tmp_path, **{"pass": False,
                                                   "coverage": 0.846}))
    assert "REJECT" in str(e.value)


def test_rule7_an_underpowered_set_cannot_license_even_when_it_passes(tmp_path):
    """§P7-7(b) made binding. 26 scored comparisons cannot resolve a 90% bar, so a
    PASS on 26 is not a pass -- and this is what keeps the frozen legacy subset from
    ever quietly becoming the licence."""
    path = _artifact(tmp_path, n_comparisons=27, n_fitted=26, n_scored=26,
                     coverage=1.0)
    with pytest.raises(G.NotCalibratedError) as e:
        G.assert_calibrated(path)
    msg = str(e.value)
    assert "§P7-9(1)" in msg and str(G.GPD_MIN_CALIBRATION_COMPARISONS) in msg
    assert G.GPD_MIN_CALIBRATION_COMPARISONS == 125


def test_rule7_the_125_is_counted_in_scored_comparisons_not_run_ones(tmp_path):
    """§P7-9(1). The exact case that withheld the first BCa licence: 135 RUN looks
    like a pass and 122 SCORED is three short. Gated-out rows entered neither bar,
    so counting them toward the design would count an absence as evidence."""
    path = _artifact(tmp_path, n_comparisons=135, n_fitted=122, n_scored=122)
    with pytest.raises(G.NotCalibratedError) as e:
        G.assert_calibrated(path)
    msg = str(e.value)
    assert "SCORED only 122" in msg and "135 run" in msg
    # and 125 scored out of the very same 135 run does license
    ok = G.assert_calibrated(_artifact(tmp_path, n_comparisons=135, n_fitted=125,
                                       n_scored=125))
    assert ok["n_scored"] == 125 and ok["n_comparisons"] == 135


def test_rule7_an_artifact_that_hides_its_scored_count_is_refused(tmp_path):
    """A run count alone cannot be checked against a design stated in scored
    comparisons, so the absence is refused rather than resolved favourably."""
    doc = {"rule": "§P6-2(7)", "ci_method": "BCa", "pass": True, "coverage": 0.94,
           "coverage_bar": 0.90, "worst_understatement": 1.8,
           "understatement_bar": 3.0, "n_comparisons": 500}
    p = tmp_path / "no_scored.json"
    p.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(G.NotCalibratedError) as e:
        G.assert_calibrated(str(p))
    assert "n_scored" in str(e.value)


def test_rule7_a_licence_is_issued_to_the_calibrated_estimator_only(tmp_path):
    with pytest.raises(G.NotCalibratedError) as e:
        G.assert_calibrated(_artifact(tmp_path, ci_method="percentile"))
    assert "not the declared BCa estimator" in str(e.value)


def test_rule7_a_passing_decisive_artifact_is_accepted_and_is_hashed(tmp_path):
    path = _artifact(tmp_path)
    cal = G.assert_calibrated(path)
    assert cal["pass"] is True and cal["n_comparisons"] == 144
    assert cal["n_scored"] == 131 >= G.GPD_MIN_CALIBRATION_COMPARISONS
    assert cal["ci_method"] == "BCa"
    # the summary carries a content hash: a re-calibration is a NEW licence
    assert len(cal["sha256"]) == 64
    _artifact(tmp_path, coverage=0.95)               # same path, new content
    assert G.load_calibration(path)["sha256"] != cal["sha256"]


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
