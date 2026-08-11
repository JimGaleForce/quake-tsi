"""Synthetic catalogue helpers for the test suite (no CSVs, no network)."""

from __future__ import annotations

import numpy as np

from engine import design, grid as gridmod


def make_ctx(n_side=12, n_days=1500, explore_frac=0.7, seed=0, lam=None, rng=None):
    """Build an EngineContext over an n_side x n_side lattice of 1-deg cells.

    `lam` (n_cells, n_days) gives the Poisson intensity; default is a smooth
    inhomogeneous field. Returns (ctx, y_counts, lam).
    """
    rng = rng or np.random.default_rng(seed)
    lat0, lon0 = 0.0, 0.0
    ii, jj = np.meshgrid(np.arange(n_side), np.arange(n_side), indexing="ij")
    clat = lat0 + ii.ravel() + 0.5
    clon = lon0 + jj.ravel() + 0.5
    n_cells = clat.size

    if lam is None:
        base = 0.01 * (1.0 + 3.0 * rng.random(n_cells))
        lam = np.repeat(base[:, None], n_days, axis=1)

    y = rng.poisson(lam).astype(np.float32)

    cell_idx, day_idx = np.nonzero(y)
    rep = y[cell_idx, day_idx].astype(int)
    ev_cell = np.repeat(cell_idx, rep)
    ev_day = np.repeat(day_idx, rep)
    ev_lat = clat[ev_cell]
    ev_lon = clon[ev_cell]

    g = gridmod.Grid(1.0, 1.0)
    g.build_domain(np.concatenate([clat, ev_lat]), np.concatenate([clon, ev_lon]))
    ec = g.cell_index(ev_lat, ev_lon)
    assert (ec >= 0).all()

    meta = dict(n_days=n_days, n_explore_days=int(round(explore_frac * n_days)),
                mag_floor=4.5, dlat=1.0, dlon=1.0)
    ctx = design._make_ctx(g, ec, ev_day, np.full(ec.size, 4.6), meta, verbose=False)
    return ctx, ctx.day_counts(4.5), lam


def make_clustered_catalog(n_side=10, n_days=1200, seed=1, background=0.004,
                           n_aftershocks=6, decay_days=12.0):
    """A crude self-exciting (ETAS-lite) catalogue: each event spawns a decaying
    cluster of nearby aftershocks. Used to show that a covariate proxying recent
    nearby activity buys real bits against a time-stationary climatology."""
    rng = np.random.default_rng(seed)
    ii, jj = np.meshgrid(np.arange(n_side), np.arange(n_side), indexing="ij")
    clat = ii.ravel() + 0.5
    clon = jj.ravel() + 0.5
    n_cells = clat.size

    rate = background * (1.0 + 2.0 * rng.random(n_cells))
    y = rng.poisson(np.repeat(rate[:, None], n_days, axis=1)).astype(np.int64)

    parents = np.array(np.nonzero(y)).T
    for c, d in parents:
        k = rng.poisson(n_aftershocks)
        if k == 0:
            continue
        dt = np.ceil(rng.exponential(decay_days, k)).astype(int)
        dd = d + dt
        ok = dd < n_days
        cc = np.full(ok.sum(), c)
        # ~70% stay in the parent cell, rest move to an immediate neighbour
        move = rng.random(ok.sum()) > 0.7
        cc = np.where(move, np.clip(cc + rng.choice([-n_side, -1, 1, n_side], cc.size), 0,
                                    n_cells - 1), cc)
        np.add.at(y, (cc, dd[ok]), 1)

    cell_idx, day_idx = np.nonzero(y)
    rep = y[cell_idx, day_idx]
    ev_cell = np.repeat(cell_idx, rep)
    ev_day = np.repeat(day_idx, rep)

    g = gridmod.Grid(1.0, 1.0)
    g.build_domain(clat, clon)
    meta = dict(n_days=n_days, n_explore_days=int(0.7 * n_days), mag_floor=4.5,
                dlat=1.0, dlon=1.0)
    ctx = design._make_ctx(g, ev_cell, ev_day, np.full(ev_cell.size, 4.6), meta,
                           verbose=False)
    return ctx, ctx.day_counts(4.5)
