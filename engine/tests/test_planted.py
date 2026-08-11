"""End-to-end test on a synthetic catalogue with a PLANTED covariate effect.

Counted invariants (SPEC):
  * the fitted beta recovers the planted beta within tolerance;
  * bits/event > 0 on the planted covariate;
  * bits/event ~ 0 (and beta ~ 0 within a few SE) on a scrambled copy.
"""

from __future__ import annotations

import numpy as np
import pytest

from engine import baseline, covariates, score
from engine.tests._synth import make_ctx, make_clustered_catalog

BETA_TRUE = 0.75


def _planted(seed=11, n_side=12, n_days=1600):
    rng = np.random.default_rng(seed)
    n_cells = n_side * n_side
    # exogenous covariate: smooth AR(1)-ish field, standardised
    z = rng.normal(size=(n_cells, n_days))
    for t in range(1, n_days):
        z[:, t] = 0.9 * z[:, t - 1] + 0.44 * z[:, t]
    z = ((z - z.mean()) / z.std()).astype(np.float32)

    mu = 0.02 * (1.0 + 3.0 * rng.random(n_cells))
    lam = mu[:, None] * np.exp(BETA_TRUE * z)
    ctx, y, _ = make_ctx(n_side=n_side, n_days=n_days, seed=seed, lam=lam, rng=rng)
    return ctx, y, z, mu


def test_planted_beta_recovered_and_bits_positive():
    ctx, y, z, _mu = _planted()
    explore = slice(0, ctx.n_explore_days)
    base = baseline.ClimatologyV1().fit(ctx, y, explore)

    fit = score.fit_poisson(base.mu, z[:, explore], y[:, explore])
    assert fit["converged"]
    assert abs(fit["beta"] - BETA_TRUE) < 6 * fit["se_beta"], (
        f"beta {fit['beta']:.4f} +/- {fit['se_beta']:.4f} missed planted {BETA_TRUE}"
    )
    assert abs(fit["beta"] - BETA_TRUE) < 0.05

    # out-of-sample: fit on exploration, score on the held-back tail
    hold = slice(ctx.n_explore_days, ctx.n_days)
    ig = score.information_gain(base.mu, z[:, hold], y[:, hold], fit["a"], fit["beta"])
    assert ig["bits_per_event"] > 0.05, ig["bits_per_event"]

    mol = score.molchan(base.mu, z[:, hold], y[:, hold], fit["a"], fit["beta"])
    assert mol["molchan_skill"] > 0.0


def test_scrambled_covariate_has_no_skill():
    ctx, y, z, _mu = _planted()
    rng = np.random.default_rng(99)
    z_s = z.copy()
    rng.shuffle(z_s, axis=1)               # destroy the cell-day pairing in time
    rng.shuffle(z_s, axis=0)

    explore = slice(0, ctx.n_explore_days)
    hold = slice(ctx.n_explore_days, ctx.n_days)
    base = baseline.ClimatologyV1().fit(ctx, y, explore)
    fit = score.fit_poisson(base.mu, z_s[:, explore], y[:, explore])
    assert abs(fit["beta"]) < 4 * fit["se_beta"], (
        f"scrambled covariate got beta {fit['beta']:.4f} +/- {fit['se_beta']:.4f}"
    )
    ig = score.information_gain(base.mu, z_s[:, hold], y[:, hold], fit["a"], fit["beta"])
    assert abs(ig["bits_per_event"]) < 0.01, ig["bits_per_event"]


def test_recent_rate_finds_planted_clustering():
    """Gap #2, demonstrated on synthetic data: a self-exciting catalogue gives the
    recent_rate covariate real bits against a time-stationary climatology."""
    ctx, y = make_clustered_catalog()
    explore = slice(0, ctx.n_explore_days)
    hold = slice(ctx.n_explore_days, ctx.n_days)
    base = baseline.ClimatologyV1().fit(ctx, y, explore)
    z, stats = covariates.compute("recent_rate", ctx, {"radius_km": 150.0, "days": 30})
    burn = stats["burn_in"]
    fit = score.fit_poisson(base.mu, z[:, burn:ctx.n_explore_days],
                            y[:, burn:ctx.n_explore_days])
    assert fit["beta"] > 0, fit
    ig = score.information_gain(base.mu, z[:, hold], y[:, hold], fit["a"], fit["beta"])
    assert ig["bits_per_event"] > 0.02, ig["bits_per_event"]


def test_etas_baseline_slot_is_implemented():
    """v1 shipped this slot raising NotImplementedError; v1.1 fills it. The ETAS
    tests proper live in engine/tests/test_etas.py."""
    assert baseline.EtasBaseline is baseline.EtasV1
    assert baseline.get_baseline("etas") is baseline.EtasV1
    assert baseline.get_baseline("climatology") is baseline.ClimatologyV1
    assert not hasattr(baseline, "ETAS_MESSAGE")


def test_covariate_registry_shape_and_finiteness():
    ctx, _y, _lam = make_ctx(n_side=8, n_days=800, seed=5)
    # synthetic events are all M4.6 -> anisotropy needs a matching m_min
    for name, params in (("recent_rate", {}), ("quiescence", {}),
                         ("anisotropy", {"m_min": 4.5})):
        z, stats = covariates.compute(name, ctx, params)
        assert z.shape == (ctx.n_cells, ctx.n_days)
        assert np.isfinite(z).all()
        assert stats["burn_in"] > 0
