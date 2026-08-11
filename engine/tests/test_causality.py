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


def test_trailing_sum_excludes_today():
    c = np.zeros((1, 6), dtype=np.float32)
    c[0, 2] = 5.0
    s = design.EngineContext.trailing_sum(c, 3)
    assert s[0].tolist() == [0.0, 0.0, 0.0, 5.0, 5.0, 5.0]
    s1 = design.EngineContext.trailing_sum(c, 1)
    assert s1[0].tolist() == [0.0, 0.0, 0.0, 5.0, 0.0, 0.0]
