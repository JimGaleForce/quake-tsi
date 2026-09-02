"""Tests for exp_productivity_assim.py (K-436).

Three things are tested, and they are the three that can silently corrupt the result:
  1. the Bayesian productivity update on a synthetic Poisson count,
  2. the Omori-window integral (both the closed form and its agreement with exp_h's
     _G_terms, which is what the likelihood actually uses),
  3. the zero-gain calibration: on a tiny simulated catalogue with NO per-sequence
     productivity scatter, the pipeline must not manufacture bits.

Run only this file:  python -u -m pytest tests/test_productivity_assim.py -q
"""

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import exp_productivity_assim as PA  # noqa: E402
from exp_h_etas import _G_terms  # noqa: E402


# ------------------------------------------------------------------ 1. the update
def _brute_posterior_mean(n_obs, lam0, mu, sd):
    f = lambda x: np.exp(-0.5 * ((x - mu) / sd) ** 2 + n_obs * x - lam0 * np.exp(x))
    num = quad(lambda x: x * f(x), mu - 40 * sd - 20, mu + 40 * sd + 20, limit=400)[0]
    den = quad(f, mu - 40 * sd - 20, mu + 40 * sd + 20, limit=400)[0]
    return num / den


def test_posterior_matches_quadrature():
    n = np.array([0.0, 1.0, 5.0, 40.0])
    lam0 = np.array([0.4, 0.4, 2.0, 3.0])
    got = PA.posterior_mean_log(n, lam0, prior_mu=-0.1, prior_sd=0.9)
    want = np.array([_brute_posterior_mean(a, b, -0.1, 0.9) for a, b in zip(n, lam0)])
    assert np.allclose(got, want, atol=2e-3), (got, want)


def test_posterior_limits_and_monotonicity():
    lam0 = np.full(1, 2.0)
    # vague prior -> close to the Poisson MLE log(n/lam0)
    vague = PA.posterior_mean_log(np.array([20.0]), lam0, 0.0, 50.0)[0]
    assert abs(vague - np.log(20.0 / 2.0)) < 0.05
    # very tight prior -> stays at the prior mean
    tight = PA.posterior_mean_log(np.array([20.0]), lam0, 0.3, 1e-3)[0]
    assert abs(tight - 0.3) < 0.02
    # monotone increasing in the observed count
    ns = np.array([0.0, 1.0, 2.0, 5.0, 10.0, 30.0])
    post = PA.posterior_mean_log(ns, np.full(ns.size, 1.5), 0.0, 1.0)
    assert np.all(np.diff(post) > 0)
    # a zero count must pull the productivity DOWN relative to the prior mean
    assert PA.posterior_mean_log(np.array([0.0]), np.array([3.0]), 0.0, 1.0)[0] < 0.0


# ------------------------------------------------------------------ 2. the Omori integral
@pytest.mark.parametrize("p", [0.85, 1.0, 1.0 + 1e-10, 1.12, 1.6])
def test_omori_int_matches_numeric(p):
    c = 0.014
    for a, b in [(0.0, 1.0 / 24.0), (0.0, 3.0 / 24.0), (0.0, 30.0), (2.0, 40.0)]:
        want = quad(lambda s: (s + c) ** (-p), a, b, limit=200)[0]
        got = float(PA.omori_int(a, b, c, p))
        assert abs(got - want) <= 1e-6 * max(1.0, abs(want)), (p, a, b, got, want)


def test_omori_int_degenerate_and_vectorised():
    assert float(PA.omori_int(5.0, 5.0, 0.01, 1.2)) == 0.0
    v = PA.omori_int(np.zeros(3), np.array([1.0, 2.0, 3.0]), 0.02, 1.1)
    assert v.shape == (3,) and np.all(np.diff(v) > 0)


def test_omori_int_agrees_with_exp_h_G_terms():
    """_G_terms is what the ETAS likelihood integrates; the update must use the same kernel."""
    c, p = 0.0143, 1.118
    t_src = np.array([0.0, 5.0, 12.5])
    T0, T1 = 3.0, 20.0
    G, _, _ = _G_terms(t_src, T0, T1, c, p, np.inf, False)
    mine = np.array([float(PA.omori_int(max(T0 - ts, 0.0), T1 - ts, c, p)) for ts in t_src])
    assert np.allclose(G, mine, rtol=1e-10, atol=1e-12), (G, mine)


# ------------------------------------------------------------------ 3. STAI pieces
def test_stai_floor_and_fraction():
    # Helmstetter et al. 2006 form, checked by hand at one point
    assert abs(PA.stai_mc(7.0, 0.01) - (7.0 - 4.5 + 1.5)) < 1e-9
    # a big mainshock loses more of its first-hour offspring than a small one
    c, p, b = 0.014, 1.12, 1.0
    big = PA.stai_expected_frac([7.0], 1 / 24.0, c, p, b)[0]
    small = PA.stai_expected_frac([4.0], 1 / 24.0, c, p, b)[0]
    plain = float(PA.omori_int(0.0, 1 / 24.0, c, p))
    assert big < small <= plain * (1 + 1e-9)
    assert small > 0.5 * plain  # M4 is barely affected at M0 = 2.5


# ------------------------------------------------------------------ 4. zero-gain calibration
def test_zero_gain_on_tiny_simulated_catalogue():
    """sd = 0: the pipeline must not manufacture bits out of Poisson noise."""
    # SUBCRITICAL by construction: n = K c^(1-p)/(p-1) * b/(b-alpha) = 0.154
    theta = (0.6, 0.005, 0.6, 0.015, 1.15)
    b = 1.0
    rng = np.random.default_rng(11)
    bits = []
    for _ in range(4):
        t, m, la, lo, meta = PA.simulate_etas(
            theta, b, 900.0, (33.0, 35.0, -119.0, -117.0), rng, sd_log10=0.0)
        assert t.size > 300, t.size
        tr_end = 900.0 * 0.6
        lo_i = int(np.searchsorted(t, tr_end, side="left"))
        sc = PA.Scorer(t, m, theta, tr_end, float(t[-1]), lo_i, t.size)
        assert np.all(np.isfinite(sc.lam_A)) and sc.lam_A.min() > 0
        upd = PA.build_update(t, la, lo, m, theta, m >= PA.M_PARENT - 1e-9,
                              float(t[0]), tr_end, 1.0 / 24.0, "plain", b)
        s, lam_B, Lam_B = PA.score_variant(sc, upd)
        assert np.all(np.isfinite(lam_B)) and lam_B.min() > 0
        bits.append(float(s.mean() / PA.LN2))
    mean_bits = float(np.mean(bits))
    assert abs(mean_bits) < 0.05, (mean_bits, bits)


def test_multiplier_of_one_reproduces_baseline_exactly():
    theta = (0.5, 0.005, 0.6, 0.015, 1.15)
    rng = np.random.default_rng(3)
    t = np.sort(rng.uniform(0, 400, 800))
    m = 2.5 - np.log10(rng.uniform(size=800)) / 1.0
    sc = PA.Scorer(t, m, theta, 200.0, float(t[-1]), int(np.searchsorted(t, 200.0)), t.size)
    pidx = np.where(m >= 4.0)[0]
    lam_B, Lam_B = sc.add_multipliers(pidx, np.ones(pidx.size), t[pidx] + 1 / 24.0)
    assert np.allclose(lam_B, sc.lam_A, rtol=0, atol=0)
    assert Lam_B == sc.Lam_A


def test_attribution_rule_picks_nearest_in_time_larger_parent():
    # two parents 5 km apart in time order; the aftershock must go to the later, larger one
    t = np.array([0.0, 1.0, 1.5])
    mag = np.array([5.0, 5.5, 3.0])
    lat = np.array([34.0, 34.0, 34.01])
    lon = np.array([-118.0, -118.0, -118.0])
    owner = PA.attribute(t, lat, lon, mag, np.array([0, 1]))
    assert owner[2] == 1
    # move it far away -> unattributed
    lat2 = np.array([34.0, 34.0, 39.0])
    owner2 = PA.attribute(t, lat2, lon, mag, np.array([0, 1]))
    assert owner2[2] == -1


def test_bpos_count_is_larger_than_predecessor_within_parent():
    # parent 0 gets 4 aftershocks with magnitudes 3, 4, 3.5, 3.6 -> larger-than-predecessor: 2
    owner = np.array([-1, 0, 0, 0, 0])
    dt = np.array([-1.0, 0.001, 0.002, 0.003, 0.004])
    mag = np.array([5.0, 3.0, 4.0, 3.5, 3.6])
    inwin = np.array([False, True, True, True, True])
    out = PA.bpos_counts(owner, dt, mag, inwin, 1)
    assert out[0] == 2.0
    # under GR marks the expectation is half the Omori expectation (frozen constant)
    assert PA.BPOS_FRACTION == 0.5
