"""SPEC: plugins must be causal -- z[:, t] may use only events strictly before day t.

The check: take a synthetic catalogue, scramble everything from day T onward
(cells AND days reshuffled among themselves), recompute the covariate, and assert
z[:, :T+1] is bit-identical. A covariate that peeks at the future fails this.
"""

from __future__ import annotations

import numpy as np
import pytest

from engine import covariates, design
from engine.tests._synth import make_ctx

# synthetic events are all M4.6, so anisotropy is given a matching m_min
COVARIATES = [
    ("recent_rate", {}),
    ("quiescence", {}),
    ("anisotropy", {"m_min": 4.5}),
]


def _shuffled_future_ctx(ctx, t_cut, seed=7):
    rng = np.random.default_rng(seed)
    cell = ctx.ev_cell.copy()
    day = ctx.ev_day.copy()
    fut = day >= t_cut
    assert fut.sum() > 50, "need a meaningful future block to scramble"
    cell[fut] = rng.permutation(cell[fut])
    day[fut] = rng.integers(t_cut, ctx.n_days, size=int(fut.sum()))
    meta = dict(ctx.meta)
    return design._make_ctx(ctx.grid, cell, day, ctx.ev_mag, meta, verbose=False)


@pytest.mark.parametrize("name,params", COVARIATES, ids=[c[0] for c in COVARIATES])
def test_covariate_is_causal(name, params):
    ctx, _y, _lam = make_ctx(n_side=8, n_days=900, seed=3)
    # The cut sits AFTER the exploration window (630 days). `quiescence` normalises by
    # an exploration-period climatology, which is a training-time constant fitted on
    # the whole exploration window -- exactly like the baseline mu. That is a
    # fit-once quantity, not a per-day peek; day-level causality is what is tested.
    t_cut = 700
    assert t_cut > ctx.n_explore_days
    fn = covariates.get(name)

    def run(c):
        p = dict(fn.defaults)
        p.update(params)
        return np.asarray(fn(c, p))

    z0 = run(ctx)
    ctx2 = _shuffled_future_ctx(ctx, t_cut)
    z1 = run(ctx2)

    assert np.array_equal(np.asarray(z0)[:, : t_cut + 1], np.asarray(z1)[:, : t_cut + 1]), (
        f"{name} changed on days <= {t_cut} after the future was scrambled -> not causal"
    )
    # sanity: the scramble did change something later on, so the test can actually fail
    assert not np.array_equal(np.asarray(z0)[:, t_cut + 5:], np.asarray(z1)[:, t_cut + 5:])


def test_mine_catalog_features_are_causal():
    """EQ-24 mine mode, family 4: the same contract, the same shuffle.

    Family 1/2 (ephemeris) features are deterministic functions of t -- no
    rearrangement of the catalogue can move them -- and are documented as exempt.
    Family 3 (downloaded indices) are lagged one day by construction. Family 4 is
    catalogue-derived and therefore has to sign the contract here.
    """
    from engine import mine as M

    ctx, _y, _lam = make_ctx(n_side=8, n_days=900, seed=3)
    t_cut = 700

    def marks_of(c):
        rng = np.random.default_rng(0)
        return {"day": np.asarray(c.ev_day), "mag": np.asarray(c.ev_mag),
                "depth": rng.random(c.ev_day.size) * 200.0}

    m0 = marks_of(ctx)
    ctx2 = _shuffled_future_ctx(ctx, t_cut)
    # keep the same per-event depths attached to the same events after the shuffle
    m1 = dict(m0)
    m1["day"] = np.asarray(ctx2.ev_day)

    f0 = {f.name: f.values for f in M.catalog_features(ctx, m0, (0,))}
    f1 = {f.name: f.values for f in M.catalog_features(ctx2, m1, (0,))}
    assert set(f0) == {"b_value_90d", "deep_fraction_30d", "mean_depth_30d"}
    for name in f0:
        assert np.array_equal(f0[name][:t_cut + 1], f1[name][:t_cut + 1]), (
            f"mine feature {name} changed on days <= {t_cut} after the future was "
            f"scrambled -> not causal")
    assert any(not np.array_equal(f0[n][t_cut + 5:], f1[n][t_cut + 5:]) for n in f0), (
        "the scramble changed nothing at all -- this test cannot fail")


def test_mine_ephemeris_features_are_exempt_and_deterministic():
    """The exemption, asserted: ephemeris features do not read the catalogue."""
    import datetime as dt

    from engine import mine as M

    ctx, _y, _lam = make_ctx(n_side=8, n_days=900, seed=3)
    ctx2 = _shuffled_future_ctx(ctx, 700)
    a = M.ephemeris_features(dt.datetime(1995, 1, 1), ctx.n_days)
    b = M.ephemeris_features(dt.datetime(1995, 1, 1), ctx2.n_days)
    assert a and all(f.causality_exempt for f in a)
    for fa, fb in zip(a, b):
        assert np.array_equal(fa.values, fb.values)


def test_trailing_sum_excludes_today():
    c = np.zeros((1, 6), dtype=np.float32)
    c[0, 2] = 5.0
    s = design.EngineContext.trailing_sum(c, 3)
    assert s[0].tolist() == [0.0, 0.0, 0.0, 5.0, 5.0, 5.0]
    s1 = design.EngineContext.trailing_sum(c, 1)
    assert s1[0].tolist() == [0.0, 0.0, 0.0, 5.0, 0.0, 0.0]


def test_mine_lagged_cyclic_features_stay_causal_and_exempt():
    """K-089-R tranche 1: turning the lag axis on must not weaken the contract.

    Two things are asserted, because the tranche adds an axis, not a data source:
      1. the lagged designs are still deterministic in t -- scrambling the future of
         the catalogue changes nothing, so the family-1/2 exemption still holds;
      2. `Feature.design(window, lag)` reads STRICTLY EARLIER days: the lag-L design
         over [a, b) is the lag-0 design over [a-L, b-L), never anything at or after
         the day being predicted. A negative lag (peeking forward) is not on the grid.
    """
    import datetime as dt

    from engine import mine as M

    ctx, _y, _lam = make_ctx(n_side=8, n_days=900, seed=3)
    ctx2 = _shuffled_future_ctx(ctx, 700)
    kw = dict(lags=M.TRANCHE1_LAGS, lag_features="tranche1")
    a = M.ephemeris_features(dt.datetime(1995, 1, 1), ctx.n_days, **kw)
    b = M.ephemeris_features(dt.datetime(1995, 1, 1), ctx2.n_days, **kw)
    assert sum(1 for f in a if len(f.lags) > 1) == 13
    win = slice(400, 900)
    for fa, fb in zip(a, b):
        assert fa.causality_exempt and fa.family in (1, 2)
        for lag in fa.lags:
            assert lag >= 0, "a negative lag would read the future"
            da, db = fa.design(win, lag), fb.design(win, lag)
            assert np.array_equal(da, db), (
                f"{fa.name} lag {lag} moved when the catalogue's future was scrambled")
            shifted = fa.design(slice(win.start - lag, win.stop - lag), 0)
            assert np.array_equal(da, shifted), (
                f"{fa.name} lag {lag} is not a strict backward shift")
    # and the lag actually does something: at least one scanned feature moves
    moved = [f for f in a if len(f.lags) > 1
             and not np.array_equal(f.design(win, 0), f.design(win, 30))]
    assert len(moved) == 13, "a lag that changes no design cannot be a new test"
