"""Unit tests for exp_dyntrig_skill (K-038 arm A).

Covers the three things that would silently corrupt the headline number:
  1. the covariate construction (the fast blocked decay against a brute-force reference,
     the frozen amplitude axis against K-034's own published readings, and the exclusions),
  2. the beta likelihood and its exact-integral quadrature identity,
  3. injection recovery of a known beta on a tiny synthetic case.

Run ONLY this file:  python -u -m pytest tests/test_dyntrig_skill.py -q
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import exp_dyntrig_skill as X  # noqa: E402


# ----------------------------------------------------------------- covariate
def test_decay_series_matches_reference():
    rng = np.random.default_rng(0)
    trig_t = np.sort(rng.uniform(0, 500, 60))
    amp = rng.uniform(1e2, 5e4, 60)
    ev = np.sort(rng.uniform(0, 500, 400))
    for tau in (1.0, 5.0, 0.3):
        fast = X.decay_series(ev, trig_t, amp, tau)
        ref = X.decay_series_ref(ev, trig_t, amp, tau)
        assert np.allclose(fast, ref, rtol=1e-9, atol=1e-9 * max(1.0, ref.max()))


def test_decay_series_single_trigger_is_exact_exponential():
    tau = 5.0
    ev = np.array([-1.0, 0.0, 1.0, 5.0, 20.0])
    out = X.decay_series(ev, np.array([0.0]), np.array([1000.0]), tau)
    assert out[0] == 0.0                                   # causal: nothing before arrival
    assert np.allclose(out[1:], 1000.0 * np.exp(-ev[1:] / tau))


def test_decay_series_is_causal_and_additive():
    t1, t2 = np.array([0.0]), np.array([10.0])
    a = np.array([100.0])
    ev = np.sort(np.array([-5.0, 0.5, 9.9, 10.1, 30.0]))
    both = X.decay_series(ev, np.array([0.0, 10.0]), np.array([100.0, 100.0]), 3.0)
    sep = (X.decay_series(ev, t1, a, 3.0) + X.decay_series(ev, t2, a, 3.0))
    assert np.allclose(both, sep)
    assert both[0] == 0.0


def test_amplitude_axis_reproduces_k034_readings():
    """K-034 reported denali->yellowstone 3110 km / Ms 8.5 at 33.8 kPa (certified floor),
    landers->cedar_city 487 km / Ms 7.3 at 46.3 kPa, landers->yellowstone 1253 km at 9.6 kPa."""
    assert X.sigma_primary_pa(8.5, 3110.0) / 1e3 == pytest.approx(33.8, rel=0.02)
    assert X.sigma_primary_pa(7.3, 487.0) / 1e3 == pytest.approx(46.3, rel=0.03)
    assert X.sigma_primary_pa(7.3, 1253.0) / 1e3 == pytest.approx(9.64, rel=0.02)


def test_amplitude_axis_scaling_is_the_frozen_form():
    # doubling distance must divide the amplitude by 2**1.66
    a1 = X.sigma_primary_pa(7.0, 500.0)
    a2 = X.sigma_primary_pa(7.0, 1000.0)
    assert a1 / a2 == pytest.approx(2.0 ** 1.66, rel=1e-9)
    # one magnitude unit is exactly a factor of 10 on this axis
    assert X.sigma_primary_pa(8.0, 500.0) / a1 == pytest.approx(10.0, rel=1e-9)


def test_link_is_frozen_log1p_at_10kPa():
    assert X.D0_PA == 10e3
    assert X.glink(0.0) == 0.0
    assert X.glink(10e3) == pytest.approx(np.log(2.0))


def test_rupture_length_brackets_k034_hand_entered_values():
    # K-034 used L = 70 / 48 / 50 km for Landers 7.3, Hector Mine 7.4(Ms), Ridgecrest 7.1
    assert 50.0 < X.rupture_len_km(7.3) < 110.0
    assert 30.0 < X.rupture_len_km(7.1) < 80.0
    assert X.rupture_len_km(8.0) > X.rupture_len_km(7.0)


def test_box_distance_and_membership():
    box = [37.4, 37.9, -119.1, -118.6]          # long_valley, from data/k034/manifest.json
    assert X.in_box(37.6, -118.8, box)
    assert not X.in_box(37.0, -118.8, box)
    assert X.dist_to_box_km(37.6, -118.8, box) == pytest.approx(0.0, abs=1e-9)
    d = X.dist_to_box_km(34.2, -116.436, box)   # Landers epicentre
    assert 300.0 < d < 600.0                     # survives the 300 km gate, as K-034 found


def test_exclusions_drop_close_and_in_box_triggers():
    import pandas as pd
    box = [37.4, 37.9, -119.1, -118.6]
    trig = pd.DataFrame(dict(
        t=[0.0, 1.0, 2.0, 3.0],
        latitude=[37.6, 37.45, 34.2, 63.517],
        longitude=[-118.8, -118.9, -116.436, -147.444],
        mag=[6.0, 6.0, 7.3, 7.9],
        in_any_k034_box=[True, False, False, False]))
    t_arr, amp, r, m, n_keep, n_tot = X.build_trigger_set(trig, box, 37.65, -118.85)
    assert n_tot == 4
    # row 0 dropped (inside a K-034 box), row 1 dropped (<300 km), rows 2 and 3 kept
    assert n_keep == 2
    assert sorted(np.round(m, 1)) == [7.3, 7.9]
    # arrival is later than origin, by r / 3.5 km/s
    assert np.all(t_arr > np.array([2.0, 3.0]) - 1e-12)


# ----------------------------------------------------------------- likelihood
def test_beta_profile_is_zero_at_beta_zero():
    A = np.array([1.0, 2.0, 3.0])
    g = np.array([0.0, 0.5, 1.0])
    gj = np.array([0.2, 0.4])
    assert X.beta_profile(gj, A, g, 0.0) == pytest.approx(0.0)


def test_beta_profile_matches_direct_poisson_loglik_difference():
    """On a piecewise-constant intensity the quadrature identity must be exact:
    LL_B - LL_A = beta*sum_j g_j - sum_k A_k (exp(beta g_k) - 1)."""
    rng = np.random.default_rng(3)
    K = 40
    dt = 0.25
    lamA = rng.uniform(0.1, 3.0, K)
    A = lamA * dt
    g = rng.uniform(0.0, 1.2, K)
    n = rng.poisson(lamA * dt * 3)
    gj = np.repeat(g, n)
    beta = 0.37
    LL_A = float((n * np.log(lamA)).sum() - (lamA * dt).sum())
    lamB = lamA * np.exp(beta * g)
    LL_B = float((n * np.log(lamB)).sum() - (lamB * dt).sum())
    assert X.beta_profile(gj, A, g, beta) == pytest.approx(LL_B - LL_A, rel=1e-12)


def test_beta_profile_is_concave_and_fit_finds_the_maximum():
    rng = np.random.default_rng(5)
    A = rng.uniform(0.1, 1.0, 200)
    g = rng.uniform(0.0, 1.0, 200)
    gj = rng.uniform(0.0, 1.0, 400)
    bhat, LLhat = X.fit_beta(gj, A, g)
    for d in (-0.05, -0.01, 0.01, 0.05):
        assert X.beta_profile(gj, A, g, bhat + d) <= LLhat + 1e-9


def test_cell_integrals_reproduce_a_pure_background():
    """With K = 0 the exact cell integral must be mu * dt in every cell."""
    edges = np.arange(0.0, 10.0 + 1e-12, 0.25)
    mu = np.array([0.7])
    ep = np.array([-1e6, 1e6])
    A = X.cell_integrals(edges, np.array([-5.0]), np.array([0.0]), 0.01, 1.1, mu, ep)
    assert np.allclose(A, 0.7 * 0.25)


def test_cell_integrals_match_numeric_quadrature_of_one_omori_burst():
    mu = np.array([0.0])
    ep = np.array([-1e6, 1e6])
    c, p = 0.02, 1.15
    t_src = np.array([1.0])
    w = np.array([2.5])                                  # already K * 10^(alpha dM)
    edges = np.arange(0.0, 6.0 + 1e-12, 0.25)
    A = X.cell_integrals(edges, t_src, w, c, p, mu, ep)
    for k in (3, 5, 10, 20):                             # cells after the burst
        lo, hi = edges[k], edges[k + 1]
        if hi <= t_src[0]:
            assert A[k] == pytest.approx(0.0, abs=1e-12)
            continue
        xs = np.linspace(max(lo, t_src[0]), hi, 20001)
        num = np.trapezoid(w[0] * (xs - t_src[0] + c) ** (-p), xs)
        assert A[k] == pytest.approx(num, rel=2e-4)


def test_etas_ll_gradient_matches_finite_differences():
    rng = np.random.default_rng(7)
    t = np.sort(rng.uniform(0.0, 300.0, 400))
    m = rng.exponential(0.5, 400)
    ep_edges = np.array([-1e6, 150.0, 1e6])
    T0, T1 = 20.0, 300.0
    lo = int(np.searchsorted(t, T0))
    hi = len(t)
    ep_tgt = (t[lo:hi] >= 150.0).astype(int)
    ep_len = X.epoch_lengths(ep_edges, T0, T1)
    theta = np.array([1.0, 1.3, 0.02, 1.0, 0.01, 1.1])
    args = (t, m, lo, hi, T0, T1, np.inf, ep_tgt, ep_len)
    LL, g = X.etas_ll(theta, *args, True, 256)
    for i in range(len(theta)):
        h = theta[i] * 1e-6
        tp = theta.copy(); tp[i] += h
        tm = theta.copy(); tm[i] -= h
        num = (X.etas_ll(tp, *args, False, 256)[0] - X.etas_ll(tm, *args, False, 256)[0]) / (2 * h)
        assert g[i] == pytest.approx(num, rel=2e-4, abs=1e-5)


# ----------------------------------------------------------------- injection recovery
def test_injection_recovery_on_a_tiny_synthetic_case():
    """Inject a known beta into a piecewise-constant intensity built from real-shaped
    dynamic-stress weather, then recover it with the same fit_beta used in the run."""
    rng = np.random.default_rng(11)
    dt = 0.25
    edges = np.arange(0.0, 4000.0 + 1e-12, dt)
    mid = 0.5 * (edges[:-1] + edges[1:])
    trig_t = np.sort(rng.uniform(0.0, 4000.0, 1200))
    amp = 10.0 ** rng.uniform(2.0, 4.6, 1200)            # 0.1 - 40 kPa, the observed range
    tau = 5.0
    g = X.glink(X.decay_series(mid, trig_t, amp, tau))
    assert g.max() > 0.3 and g.mean() > 0.0
    lamA = 4.0                                            # events/day
    A = np.full(mid.size, lamA * dt)
    for beta_true in (0.0, 0.25, -0.25):
        n = rng.poisson(A * np.exp(beta_true * g))
        gj = np.repeat(g, n)
        bhat, _ = X.fit_beta(gj, A, g)
        assert bhat == pytest.approx(beta_true, abs=0.06), (beta_true, bhat)


def test_injection_recovery_is_unbiased_over_replicates():
    rng = np.random.default_rng(13)
    dt = 0.25
    mid = np.arange(0.0, 2000.0, dt) + dt / 2
    trig_t = np.sort(rng.uniform(0.0, 2000.0, 600))
    amp = 10.0 ** rng.uniform(2.0, 4.6, 600)
    g = X.glink(X.decay_series(mid, trig_t, amp, 1.0))
    A = np.full(mid.size, 3.0 * dt)
    beta_true = 0.3
    hats = []
    for _ in range(30):
        n = rng.poisson(A * np.exp(beta_true * g))
        hats.append(X.fit_beta(np.repeat(g, n), A, g)[0])
    assert np.mean(hats) == pytest.approx(beta_true, abs=0.05)


def test_zero_injection_gives_zero_expected_bits():
    """The scoring statistic must be centred at ~0 when there is no signal."""
    rng = np.random.default_rng(17)
    dt = 0.25
    mid = np.arange(0.0, 3000.0, dt) + dt / 2
    trig_t = np.sort(rng.uniform(0.0, 3000.0, 900))
    amp = 10.0 ** rng.uniform(2.0, 4.6, 900)
    g = X.glink(X.decay_series(mid, trig_t, amp, 5.0))
    half = mid.size // 2
    A = np.full(mid.size, 2.0 * dt)
    bits = []
    for _ in range(40):
        n = rng.poisson(A)
        b, _ = X.fit_beta(np.repeat(g[:half], n[:half]), A[:half], g[:half])
        gj_sc = np.repeat(g[half:], n[half:])
        bits.append(X.beta_profile(gj_sc, A[half:], g[half:], b) / max(1, n[half:].sum()) / np.log(2))
    assert abs(float(np.mean(bits))) < 0.005


# ----------------------------------------------------------------- guards
def test_sequence_ids_group_a_burst_and_separate_distant_events():
    t = np.array([0.0, 0.1, 0.5, 100.0])
    lat = np.array([37.5, 37.51, 37.52, 37.5])
    lon = np.array([-118.8, -118.81, -118.82, -118.8])
    mag = np.array([4.0, 2.0, 2.0, 3.0])
    s = X.sequence_ids(t, lat, lon, mag)
    assert s[0] == s[1] == s[2]
    assert s[3] != s[0]


def test_frozen_constants_are_the_k034_values():
    assert X.G_SHEAR == 30e9
    assert X.C_PHASE == 3500.0
    assert X.T_SW == 20.0
    assert X.TAUS == (1.0, 5.0)
    assert X.EXCL_BOX_KM == 300.0
    assert X.TRIG_MINMAG == 5.5
    assert X.EXPLORE_FRAC == 0.70 and X.TRAIN_FRAC_OF_EXPLORE == 0.60
    assert len(X.CLASS) == 14
