"""Acceptance tests for engine/gate_calibration.py (SP-7 v2).

The point of these is that the criterion must be DERIVED. So the tests check
properties and identities, not values: that the band brackets the mean, that its
achieved error rate never exceeds the declared one, that it widens with n and narrows
with the declared rate, that it reproduces SP-7 v1's own arithmetic when handed v1's
design, and that the two-sided verdict actually distinguishes the two failure modes.
A test that pinned a threshold integer would defeat the module's whole purpose.
"""

from __future__ import annotations

import math

import pytest

from engine import gate_calibration as GC


def test_binomial_tails_agree_with_scipy():
    sp = pytest.importorskip("scipy.stats")
    for n, p in ((900, 1 / 401.0), (900, 1 / 300.0), (60, 0.1), (1000, 0.02)):
        for k in (0, 1, 3, 5, 12):
            assert GC.binom_sf(k, n, p) == pytest.approx(sp.binom.sf(k - 1, n, p),
                                                         rel=1e-9, abs=1e-15)
            assert GC.binom_cdf(k, n, p) == pytest.approx(sp.binom.cdf(k, n, p),
                                                          rel=1e-9, abs=1e-15)


def test_pmf_sums_to_one():
    assert sum(GC.binom_pmf(k, 40, 0.07) for k in range(41)) == pytest.approx(1.0)


def test_effective_alpha_is_the_achievable_rate_not_the_declared_one():
    """floor(alpha*(B+1))/(B+1). The first version returned max(alpha, 1/(B+1)).

    For the SP-7 design (alpha = 1/300, B = 400) the achievable p-values are
    k/401, so `p <= 1/300` fires exactly when p = 1/401 and the TRUE false-positive
    rate is 1/401 = 0.0024938, not 1/300 = 0.0033333. Getting this wrong overstates
    the expected promotion count by 34 percent and the error propagates into the band.
    """
    assert GC.effective_alpha(1 / 300.0, 400) == pytest.approx(1 / 401.0)
    assert GC.effective_alpha(0.05, 999) == pytest.approx(50 / 1000.0)
    assert GC.effective_alpha(1e-6, 999999) == pytest.approx(1e-6)
    # finer than the resampling can resolve: the test cannot fire at any effect size
    assert GC.effective_alpha(1e-6, 400) == 0.0


# --------------------------------------------------------------- the band itself --
def test_band_brackets_the_mean():
    """A criterion that excluded its own expected value would be the v1 defect again."""
    for n, p in ((900, 1 / 401.0), (900, 1 / 300.0), (3000, 0.01), (100, 0.05)):
        b = GC.calibration_band(n, p, 0.01)
        assert b["k_fail_low"] < b["expected_count"] < b["k_fail_high"], (n, p, b)


def test_achieved_false_fail_never_exceeds_declared():
    for rate in (0.001, 0.01, 0.05, 0.10):
        b = GC.calibration_band(900, 1 / 401.0, rate)
        assert b["achieved_false_fail_total"] <= rate + 1e-12


def test_band_widens_as_the_declared_rate_tightens():
    tight = GC.calibration_band(900, 1 / 401.0, 0.001)
    loose = GC.calibration_band(900, 1 / 401.0, 0.10)
    assert tight["k_fail_high"] >= loose["k_fail_high"]
    assert tight["k_fail_low"] <= loose["k_fail_low"]


def test_no_threshold_is_hardcoded_the_band_tracks_the_design():
    """Double the trials and the band must move. A written-down cap could not."""
    small = GC.calibration_band(900, 1 / 401.0, 0.01)
    big = GC.calibration_band(9000, 1 / 401.0, 0.01)
    assert big["k_fail_high"] > small["k_fail_high"]
    assert big["expected_count"] == pytest.approx(10.0 * small["expected_count"])


def test_v1_cap_of_three_is_the_mean_not_a_bound():
    """The arithmetic that motivated the amendment, asserted rather than asserted-to.

    N * m * (q/m) = N * q = 3.0 at N = 30, q = 0.10, REGARDLESS of m. A cap set at the
    mean rejects a correctly calibrated searcher about a third of the time.
    """
    for m in (10, 30, 100, 1000):
        n_trials = 30 * m
        alpha = 0.10 / m
        assert n_trials * alpha == pytest.approx(3.0)
    # Poisson(3) upper tail beyond a cap of 3
    poisson_sf_3 = 1.0 - math.exp(-3.0) * (1 + 3 + 4.5 + 4.5)
    assert poisson_sf_3 == pytest.approx(0.3528, abs=5e-4)
    # and the derived band at that design must be strictly wider than the v1 cap
    b = GC.calibration_band(900, 0.10 / 30.0, 0.01)
    assert b["k_fail_high"] > 3


# ------------------------------------------------------------------ power terms --
def test_detectable_inflation_is_above_one_and_monotone_in_n():
    b1 = GC.calibration_band(900, 1 / 401.0, 0.01)
    b2 = GC.calibration_band(9000, 1 / 401.0, 0.01)
    r1 = GC.detectable_inflation(900, 1 / 401.0, b1["k_fail_high"])
    r2 = GC.detectable_inflation(9000, 1 / 401.0, b2["k_fail_high"])
    assert r1 > 1.0 and r2 > 1.0
    assert r2 < r1, "more null catalogues must detect a smaller miscalibration"


def test_detectable_inflation_actually_achieves_the_power():
    n, p = 900, 1 / 401.0
    b = GC.calibration_band(n, p, 0.01)
    r = GC.detectable_inflation(n, p, b["k_fail_high"], power=0.80)
    assert GC.binom_sf(b["k_fail_high"], n, r * p) >= 0.80 - 1e-6
    assert GC.binom_sf(b["k_fail_high"], n, 0.98 * r * p) < 0.80 + 1e-6


def test_required_catalogues_inverts_detectable_inflation():
    got = GC.required_catalogues(m=30, alpha_eff=1 / 401.0, r_min=3.0, power=0.80)
    assert got is not None
    assert got["detectable_inflation"] <= 3.0 + 1e-9
    smaller = GC.calibration_band(30 * (got["n_catalogues"] - 1), 1 / 401.0, 0.01)
    r_smaller = GC.detectable_inflation(30 * (got["n_catalogues"] - 1), 1 / 401.0,
                                        smaller["k_fail_high"], 0.80)
    assert r_smaller is None or r_smaller > 3.0, "must be the MINIMUM n"


# -------------------------------------------------------------------- verdicts --
def test_verdict_distinguishes_the_two_failure_modes():
    n, p = 900, 1 / 401.0
    b = GC.calibration_band(n, p, 0.01)
    assert GC.verdict(b["k_fail_high"], n, p)["verdict"] == "FAIL-HIGH"
    assert GC.verdict(b["k_fail_low"], n, p)["verdict"] == "FAIL-LOW"
    mid = int(round(n * p))
    assert GC.verdict(mid, n, p)["verdict"] == "PASS"


def test_zero_promotions_is_a_failure_not_a_triumph():
    """'We saw none' is evidence of a bug, not of virtue -- v1 could not say this."""
    v = GC.verdict(0, 9000, 1 / 401.0)
    assert v["verdict"] == "FAIL-LOW"
    assert "over-conservative" in v["interpretation"]


def test_vacuity_dominates_every_other_verdict():
    v = GC.verdict(5, 900, 1 / 401.0, vacuous=True, vacuity_reason="control never fired")
    assert v["verdict"] == "VACUOUS-FAIL" and not v["passed"]


def test_verdict_carries_its_own_power_number():
    v = GC.verdict(3, 900, 1 / 401.0)
    assert v["power"]["detectable_inflation_R"] > 1.0
    assert v["band"]["achieved_false_fail_total"] <= 0.01 + 1e-12


def test_upper_tail_is_correct_far_below_the_mode():
    """The recorded bug: k far below n*p underflowed and returned ~0 instead of ~1.

    This is the regime detectable_inflation probes at large R, so the wrong answer
    silently turned a detectable miscalibration into an undetectable one and made
    required_catalogues return None.
    """
    sp = pytest.importorskip("scipy.stats")
    for n, k, p in ((9000, 37, 50 / 401.0), (9000, 37, 100 / 401.0),
                    (900, 8, 100 / 401.0)):
        assert GC.binom_sf(k, n, p) == pytest.approx(sp.binom.sf(k - 1, n, p),
                                                     rel=1e-9, abs=1e-15)
        assert GC.binom_sf(k, n, p) > 0.99
        # the guarded fallback must agree with the exact path in this regime too
        assert GC._binom_sf_recurrence(k, n, p) > 0.99


def test_fallback_agrees_with_the_exact_path():
    for n, p in ((900, 1 / 401.0), (900, 3 / 401.0), (9000, 1 / 401.0)):
        for k in (1, 2, 5, 8, 20, 37):
            if k > n:
                continue
            assert GC._binom_sf_recurrence(k, n, p) == pytest.approx(
                GC.binom_sf(k, n, p), rel=1e-9, abs=1e-15)
