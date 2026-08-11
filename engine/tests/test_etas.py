"""Tests for the ETAS baseline (engine/baseline.py:EtasV1).

Three counted invariants, matching the three things that can silently be wrong in a
self-exciting baseline:

  * CAUSALITY -- lambda(c, t) may use only events strictly before day t. Scrambling
    everything from a cut day onward must not move lambda at or before the cut.
  * RECOVERY  -- on a synthetic catalogue that really is ETAS plus a planted
    multiplicative covariate effect, fitting ETAS as the offset must recover the
    planted beta, and the SAME covariate must lose most of its apparent skill when
    the baseline stops being a time-stationary climatology (the gate, in miniature).
  * CACHE     -- fitted parameters round-trip through etas_params.json and reproduce
    the identical intensity field; a stale key is rejected rather than reused.
"""

from __future__ import annotations

import json
import os

import numpy as np
import pytest

from engine import baseline, covariates, design, score
from engine.tests._synth import make_ctx, make_etas_catalog

FIXED = {"nu": 0.6, "K": 0.4, "alpha": 1.0, "c": 0.05, "p": 1.15, "sigma": 80.0}


def _shuffled_future_ctx(ctx, t_cut, seed=7):
    rng = np.random.default_rng(seed)
    cell = ctx.ev_cell.copy()
    day = ctx.ev_day.copy()
    mag = ctx.ev_mag.copy()
    fut = day >= t_cut
    assert fut.sum() > 50, "need a meaningful future block to scramble"
    cell[fut] = rng.permutation(cell[fut])
    mag[fut] = rng.permutation(mag[fut])
    day[fut] = rng.integers(t_cut, ctx.n_days, size=int(fut.sum()))
    return design._make_ctx(ctx.grid, cell, day, mag, dict(ctx.meta), verbose=False)


# ------------------------------------------------------------- causality ---
def test_etas_intensity_is_causal():
    """SPEC plugin contract, extended to the baseline: no peeking past day t."""
    ctx, y, _lam = make_ctx(n_side=8, n_days=900, seed=3)
    t_cut = 700
    mu = baseline.ClimatologyV1().fit(ctx, y, slice(0, ctx.n_explore_days)).mu

    lam0 = baseline.EtasKernel(ctx, mu, trunc_days=120).intensity(FIXED)
    ctx2 = _shuffled_future_ctx(ctx, t_cut)
    lam1 = baseline.EtasKernel(ctx2, mu, trunc_days=120).intensity(FIXED)

    assert np.array_equal(lam0[:, : t_cut + 1], lam1[:, : t_cut + 1]), (
        "ETAS intensity moved on days <= the cut after the future was scrambled "
        "-> the baseline is not causal"
    )
    # sanity: the scramble did change the future, so this test can actually fail
    assert not np.array_equal(lam0[:, t_cut + 5:], lam1[:, t_cut + 5:])


def test_etas_triggering_starts_the_day_after_the_event():
    """A single event contributes nothing on its own day and everything after it."""
    ctx, y, _ = make_ctx(n_side=4, n_days=40, seed=2)
    ctx2 = design._make_ctx(ctx.grid, np.array([0]), np.array([10]),
                            np.array([6.0]), dict(ctx.meta), verbose=False)
    mu = np.full(ctx2.n_cells, 1e-6)
    lam = baseline.EtasKernel(ctx2, mu, trunc_days=20).intensity(FIXED)
    trig = lam - FIXED["nu"] * mu[:, None]
    assert trig[:, :11].max() < 1e-9, "triggering leaked into the event's own day"
    assert trig[0, 11] > 0
    assert trig[0, 11] > trig[0, 12] > trig[0, 20]        # Omori decay
    assert trig[:, 31:].max() < 1e-9                      # truncation at 20 d


# ---------------------------------------------------- planted under ETAS ---
BETA_TRUE = 0.60


@pytest.fixture(scope="module")
def planted_etas():
    return make_etas_catalog(n_side=8, n_days=900, seed=4, beta_true=BETA_TRUE)


def _fit_etas(ctx, y, tmp_path, name="etas_params.json", **kw):
    b = baseline.EtasV1(cache_path=os.path.join(str(tmp_path), name), verbose=False,
                        polish=False, trunc_days=120, **kw)
    return b.fit(ctx, y, slice(0, ctx.n_explore_days))


def test_planted_effect_recovered_under_etas_baseline(planted_etas, tmp_path):
    ctx, y, z = planted_etas
    base = _fit_etas(ctx, y, tmp_path)
    sl = slice(base.burn_in_days, ctx.n_explore_days)

    fit = score.fit_poisson(base.rate(sl), z[:, sl], y[:, sl])
    assert fit["converged"]
    assert abs(fit["beta"] - BETA_TRUE) < 6 * fit["se_beta"], (
        f"beta {fit['beta']:.4f} +/- {fit['se_beta']:.4f} missed planted {BETA_TRUE}"
    )
    assert abs(fit["beta"] - BETA_TRUE) < 0.06
    ig = score.information_gain(base.rate(sl), z[:, sl], y[:, sl], fit["a"], fit["beta"])
    assert ig["bits_per_event"] > 0.02, ig["bits_per_event"]

    # a scrambled copy of the same covariate must buy nothing against the same offset
    rng = np.random.default_rng(5)
    zs = z.copy()
    rng.shuffle(zs, axis=1)
    rng.shuffle(zs, axis=0)
    fs = score.fit_poisson(base.rate(sl), zs[:, sl], y[:, sl])
    assert abs(fs["beta"]) < 5 * fs["se_beta"], fs


def test_etas_baseline_beats_climatology_on_a_clustered_catalogue(planted_etas,
                                                                  tmp_path):
    ctx, y, _z = planted_etas
    sl = slice(120, ctx.n_explore_days)
    base = _fit_etas(ctx, y, tmp_path, name="beats.json")
    clim = baseline.ClimatologyV1().fit(ctx, y, slice(0, ctx.n_explore_days))
    ll_clim, _ = score.baseline_ll(clim.mu, y[:, sl])
    ll_etas, _ = score.baseline_ll(base.rate(sl), y[:, sl])
    assert ll_etas > ll_clim, (ll_etas, ll_clim)
    assert base.fit_info["bits_per_event_vs_climatology"] > 0.05, base.fit_info


def test_recent_rate_skill_collapses_under_etas(planted_etas, tmp_path):
    """THE GATE, in miniature: aftershock clustering must stop paying once the
    baseline models it. Same catalogue, same covariate, two baselines."""
    ctx, y, _z = planted_etas
    sl = slice(120, ctx.n_explore_days)
    zr, _ = covariates.compute("recent_rate", ctx, {"radius_km": 150.0, "days": 30})

    clim = baseline.ClimatologyV1().fit(ctx, y, slice(0, ctx.n_explore_days))
    f_c = score.fit_poisson(clim.mu, zr[:, sl], y[:, sl])
    b_c = score.information_gain(clim.mu, zr[:, sl], y[:, sl],
                                 f_c["a"], f_c["beta"])["bits_per_event"]

    etas = _fit_etas(ctx, y, tmp_path, name="collapse.json")
    f_e = score.fit_poisson(etas.rate(sl), zr[:, sl], y[:, sl])
    b_e = score.information_gain(etas.rate(sl), zr[:, sl], y[:, sl],
                                 f_e["a"], f_e["beta"])["bits_per_event"]

    assert b_c > 0.015, f"climatology should hand out free skill here, got {b_c}"
    assert b_e < 0.35 * b_c, f"no collapse: climatology {b_c:.4f} -> etas {b_e:.4f}"


# ----------------------------------------------------------- param cache ---
def test_params_cache_round_trip(planted_etas, tmp_path):
    ctx, y, _z = planted_etas
    path = os.path.join(str(tmp_path), "rt.json")
    b1 = baseline.EtasV1(cache_path=path, verbose=False, polish=False, trunc_days=120)
    b1.fit(ctx, y, slice(0, ctx.n_explore_days))
    assert b1.fit_info["source"] == "fitted"
    assert os.path.exists(path)

    rec = json.load(open(path, encoding="utf-8"))
    assert set(rec["params"]) == set(baseline.PARAM_NAMES)
    assert rec["fit_window_days"][1] == ctx.n_explore_days
    assert rec["omori_trunc_days"] == 120
    assert np.isfinite(rec["ll_etas"]) and rec["ll_etas"] > rec["ll_climatology"]

    # second construction must LOAD, not refit, and reproduce lambda bit-for-bit
    b2 = baseline.EtasV1(cache_path=path, verbose=False, polish=False, trunc_days=120)
    b2.fit(ctx, y, slice(0, ctx.n_explore_days))
    assert b2.fit_info["source"] == "cache"
    assert b2.params == b1.params
    assert np.array_equal(b1.lam, b2.lam)

    # a cache written for a different design must be REJECTED, not silently reused
    rec2 = dict(rec)
    rec2["key"] = dict(rec["key"], design_cache_key="something-else")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(rec2, fh)
    b3 = baseline.EtasV1(cache_path=path, verbose=False, polish=False, trunc_days=120,
                         fit_if_missing=False)
    with pytest.raises(FileNotFoundError):
        b3.fit(ctx, y, slice(0, ctx.n_explore_days))


def test_baseline_caveat_switches_per_baseline():
    from engine import baseline_caveat, canonical_baseline
    assert "clustering NOT absorbed" in baseline_caveat("climatology")
    assert baseline_caveat("etas") == (
        "baseline=etas-v1 (isotropic kernel; anisotropy/mechanism NOT absorbed)")
    assert baseline_caveat("etas-v1") == baseline_caveat("etas")
    assert canonical_baseline("climatology") == "climatology-v1"
    with pytest.raises(KeyError):
        baseline_caveat("nope")


def test_score_accepts_1d_and_2d_offsets():
    """The intercept-refit convention must be identical for both offset shapes."""
    ctx, y, _ = make_ctx(n_side=6, n_days=400, seed=8)
    mu = baseline.ClimatologyV1().fit(ctx, y, slice(0, ctx.n_days)).mu
    z = np.zeros((ctx.n_cells, ctx.n_days), dtype=np.float32)
    z[: ctx.n_cells // 2] = 1.0
    mu2 = np.repeat(mu[:, None], ctx.n_days, axis=1)     # same rate, 2-D spelling
    f1 = score.fit_poisson(mu, z, y)
    f2 = score.fit_poisson(mu2, z, y)
    assert abs(f1["beta"] - f2["beta"]) < 1e-8
    assert abs(f1["ll"] - f2["ll"]) < 1e-6
    i1 = score.information_gain(mu, z, y, f1["a"], f1["beta"])
    i2 = score.information_gain(mu2, z, y, f2["a"], f2["beta"])
    assert abs(i1["bits_per_event"] - i2["bits_per_event"]) < 1e-9
    m1 = score.molchan(mu, z, y, f1["a"], f1["beta"])
    m2 = score.molchan(mu2, z, y, f2["a"], f2["beta"])
    assert abs(m1["molchan_skill"] - m2["molchan_skill"]) < 1e-9
    with pytest.raises(ValueError):
        score.fit_poisson(mu2[:, :10], z, y)              # window mismatch caught
