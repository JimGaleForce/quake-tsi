"""Tests for exp_neural_tpp.py.

Run ONLY this file (other workers share the machine):
    python -u -m pytest tests/test_neural_tpp.py -q

What is covered, and why each one exists:
  * the GRU encoder is STRICTLY CAUSAL -- the single failure mode that would manufacture
    arbitrary bits out of nothing;
  * the ETAS <-> TPP likelihood conversion agrees with the closed form on a tiny synthetic
    (a homogeneous Poisson process, where log f(tau) = log mu - mu tau exactly), and the
    segment machinery reproduces the un-segmented interval integral;
  * the log-normal mixture density integrates to 1 over tau in (0, inf);
  * the torch ETAS port reproduces exp_h_etas.etas_ll in value and gradient;
  * build_segments tiles every interval exactly (no gaps, no overlap).
"""

import math
import os
import sys

import numpy as np
import pytest
import torch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import exp_h_etas as eh                      # noqa: E402
import exp_neural_tpp as X                   # noqa: E402

LN10 = math.log(10.0)
DEV = torch.device("cpu")


# --------------------------------------------------------------------- causality
def test_encoder_is_strictly_causal():
    """Perturbing event j must leave every output at position < j bit-identical."""
    torch.manual_seed(0)
    model = X.NeuralTPP(n_feat=4, hidden=16, layers=1, n_mix=8, dropout=0.0).double()
    # the head is deliberately zero-initialised (a constant marginal density at step 0), so
    # give it real weights before asking whether outputs react at all
    with torch.no_grad():
        model.head.weight.normal_(0.0, 0.5)
    x = torch.randn(1, 40, 4, dtype=torch.float64)
    with torch.no_grad():
        p0, _ = model(x)
        for j in (5, 17, 31):
            x2 = x.clone()
            x2[0, j] += 3.7
            p1, _ = model(x2)
            assert torch.equal(p0[0, :j], p1[0, :j]), f"leak: outputs before {j} changed"
            assert not torch.allclose(p0[0, j], p1[0, j]), f"position {j} did not react"


def test_full_pass_matches_single_shot_and_is_causal():
    """The chunked, state-carrying evaluation pass equals one unchunked pass."""
    torch.manual_seed(1)
    model = X.NeuralTPP(4, 16, 1, 8, 0.0).double()
    N = 300
    Xf = np.random.default_rng(0).normal(size=(N, 4))
    Y = np.random.default_rng(1).normal(size=N)
    a = X.nn_full_pass(model, Xf, Y, DEV, 8, chunk=1000)
    b = X.nn_full_pass(model, Xf, Y, DEV, 8, chunk=37)
    assert np.isnan(a[0]) and np.isnan(b[0])
    assert np.allclose(a[1:], b[1:], atol=1e-10)

    # and: logf at position i must not depend on anything at index > i
    Xf2 = Xf.copy()
    Xf2[200:] += 5.0
    c = X.nn_full_pass(model, Xf2, Y, DEV, 8, chunk=1000)
    assert np.allclose(a[1:200], c[1:200], atol=1e-12)


# --------------------------------------------------------------------- mixture density
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_mixture_density_integrates_to_one(seed):
    """int_0^inf f(tau) dtau = 1, evaluated by substituting tau = e^y."""
    torch.manual_seed(seed)
    K = 12
    params = torch.randn(1, 3 * K, dtype=torch.float64)
    y = torch.linspace(-30, 30, 400001, dtype=torch.float64)
    logf_tau = X.NeuralTPP.log_f_tau(params.expand(y.numel(), -1), y, K)
    # dtau = e^y dy, so int f(tau) dtau = int f(tau(y)) e^y dy
    integrand = torch.exp(logf_tau + y)
    total = torch.trapz(integrand, y).item()
    assert abs(total - 1.0) < 1e-6, f"mixture integrates to {total}"


def test_mixture_is_a_proper_lognormal_mixture():
    """With one component the density must equal the analytic log-normal."""
    K = 1
    mu, sd = 0.7, 1.3
    raw_sd = math.log(math.expm1(sd))
    params = torch.tensor([[0.0, mu, raw_sd]], dtype=torch.float64)
    tau = torch.tensor([0.05, 0.5, 2.0, 17.0], dtype=torch.float64)
    y = torch.log(tau)
    got = X.NeuralTPP.log_f_tau(params.expand(4, -1), y, K).numpy()
    sd_eff = torch.nn.functional.softplus(torch.tensor(raw_sd)).item() + 1e-4
    want = (-0.5 * ((y.numpy() - mu) / sd_eff) ** 2 - math.log(sd_eff)
            - 0.5 * math.log(2 * math.pi) - y.numpy())
    assert np.allclose(got, want, atol=1e-12)


# --------------------------------------------------------------------- ETAS <-> TPP
def test_etas_tpp_conversion_poisson_limit():
    """K = 0 makes ETAS a homogeneous Poisson process: log f(tau) = log mu - mu tau."""
    rng = np.random.default_rng(3)
    n = 400
    mu = 2.5
    t = np.cumsum(rng.exponential(1.0 / mu, size=n))
    m = rng.exponential(0.4, size=n)
    E = X.TorchETAS(t, m, DEV)
    theta = np.array([mu, 0.0, 1.0, 0.01, 1.1])
    lo, hi = 100, n
    lam = E.lam_at(t[lo:hi], theta, np.inf, src_hi=hi)
    Lam = E.seg_A(t[lo - 1:hi - 1], t[lo:hi], theta, np.inf, src_hi=hi)
    logf = np.log(lam) - Lam
    tau = np.diff(t[lo - 1:hi])
    want = math.log(mu) - mu * tau
    assert np.allclose(logf, want, atol=1e-10), np.max(np.abs(logf - want))


def test_interval_integral_matches_exp_h_G_terms():
    """seg_A must reproduce exp_h_etas._G_terms on the same interval and parameters."""
    rng = np.random.default_rng(4)
    n = 500
    t = np.sort(rng.uniform(0, 200, size=n))
    m = rng.exponential(0.5, size=n)
    theta = np.array([0.3, 0.02, 1.1, 0.02, 1.15])
    mu, K, alpha, c, p = theta
    E = X.TorchETAS(t, m, DEV)
    for W in (np.inf, 30.0):
        for (T0, T1) in [(50.0, 51.3), (10.0, 190.0), (120.0, 120.05)]:
            got = E.seg_A(np.array([T0]), np.array([T1]), theta, W)[0]
            G, _, _ = eh._G_terms(t[t <= T1], T0, T1, c, p, W, False)
            w = np.exp(alpha * LN10 * m[t <= T1])
            want = mu * (T1 - T0) + K * float((w * G).sum())
            assert abs(got - want) < 1e-9 * max(1.0, abs(want)), (W, T0, T1, got, want)


def test_segments_tile_intervals_exactly():
    """Splitting at solar-hour boundaries must lose nothing and overlap nothing."""
    rng = np.random.default_rng(5)
    t = np.sort(rng.uniform(1000.0, 1010.0, size=200))
    lo, hi, b, idx = X.build_segments(t, ref_lon=-116.5)
    assert np.all(hi > lo)
    assert b.min() >= 0 and b.max() < X.N_HOUR_BINS
    for i in range(len(t) - 1):
        s = idx == i
        if t[i + 1] <= t[i]:
            continue
        assert abs(lo[s].min() - t[i]) < 1e-12
        assert abs(hi[s].max() - t[i + 1]) < 1e-12
        assert abs((hi[s] - lo[s]).sum() - (t[i + 1] - t[i])) < 1e-12
        # contiguity
        o = np.argsort(lo[s])
        assert np.allclose(lo[s][o][1:], hi[s][o][:-1], atol=1e-12)


def test_segment_sum_reproduces_unsegmented_integral():
    """sum of segment integrals over an interval == the interval integral (g == 1 case)."""
    rng = np.random.default_rng(6)
    n = 300
    t = np.sort(rng.uniform(0, 60, size=n))
    m = rng.exponential(0.5, size=n)
    theta = np.array([0.4, 0.03, 1.0, 0.02, 1.2])
    E = X.TorchETAS(t, m, DEV)
    lo_idx, hi_idx = 100, n
    bounds = t[lo_idx - 1:hi_idx]
    slo, shi, sbin, sidx = X.build_segments(bounds, ref_lon=-116.5)
    A = E.seg_A(slo, shi, theta, np.inf, src_hi=hi_idx)
    Lam_seg = np.zeros(hi_idx - lo_idx)
    np.add.at(Lam_seg, sidx, A)
    Lam_dir = E.seg_A(t[lo_idx - 1:hi_idx - 1], t[lo_idx:hi_idx], theta, np.inf, src_hi=hi_idx)
    assert np.allclose(Lam_seg, Lam_dir, rtol=1e-10, atol=1e-12)


def test_torch_etas_matches_exp_h_value_and_gradient():
    rng = np.random.default_rng(7)
    n = 900
    t = np.sort(rng.uniform(0, 400, size=n))
    m = rng.exponential(0.5, size=n)
    out = X.verify_torch_matches_numpy(t, m, DEV, log=lambda *a: None)
    assert out["rel_LL"] < 1e-9
    assert out["max_rel_grad"] < 1e-6


# --------------------------------------------------------------------- split invariants
def test_split_has_no_future_leakage():
    """The frozen split must put every test event strictly after every training event."""
    t = np.arange(0.0, 1000.0, 0.01)          # 100k events, uniform
    T_first, T_last = t[0], t[-1]
    span = T_last - T_first
    T_expl = T_first + X.EXPLORE_FRAC * span
    e = T_expl - T_first
    T_tr = T_first + X.TRAIN_FRAC * e
    T_va = T_first + (X.TRAIN_FRAC + X.VAL_FRAC) * e
    tr_hi = int(np.searchsorted(t, T_tr, "left"))
    va_hi = int(np.searchsorted(t, T_va, "left"))
    te_hi = int(np.searchsorted(t, T_expl, "left"))
    assert 0 < tr_hi < va_hi < te_hi < len(t)
    assert t[va_hi] > t[:tr_hi].max()
    assert te_hi / len(t) == pytest.approx(X.EXPLORE_FRAC, abs=1e-3)
    # the holdout is never inside the exploration slice
    assert len(t) - te_hi > 0
