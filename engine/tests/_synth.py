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


def make_etas_catalog(n_side=8, n_days=900, seed=4, nu_mu=0.02, branching=0.55,
                      p=1.15, c0=0.05, sigma_km=80.0, trunc=120, beta_true=0.0,
                      z=None, z_phi=0.3, explore_frac=0.7):
    """Simulate a space-time ETAS catalogue on the engine's own cell-day lattice.

    lambda(c, t) = (mu(c) + K * sum_i (t - t_i + c0)^-p * G_sigma(d(c, cell_i)))
                   * exp(beta_true * z(c, t))

    The planted covariate multiplies the TOTAL intensity, which is exactly the model
    score.py fits on top of a baseline offset (lambda = offset * exp(a + beta z)), so
    a correct ETAS baseline + Poisson fit must recover beta_true. Returns
    (ctx, y_counts, z).
    """
    rng = np.random.default_rng(seed)
    ii, jj = np.meshgrid(np.arange(n_side), np.arange(n_side), indexing="ij")
    clat = ii.ravel() + 0.5
    clon = jj.ravel() + 0.5
    n_cells = clat.size

    g = gridmod.Grid(1.0, 1.0)
    g.build_domain(clat, clon)
    d = g.pair_distance_km().astype(np.float64)
    G = np.exp(-0.5 * (d / sigma_km) ** 2)
    G /= G.sum(axis=0, keepdims=True)

    if z is None:
        # Deliberately SHORT-memory in time (z_phi=0.3). A strongly autocorrelated
        # covariate is partly proxied by the ETAS triggering term itself -- the
        # baseline eats some of the planted effect and beta comes back attenuated,
        # which is the correct behaviour but makes a poor recovery test.
        z = rng.normal(size=(n_cells, n_days))
        for t in range(1, n_days):
            z[:, t] = z_phi * z[:, t - 1] + np.sqrt(1 - z_phi ** 2) * z[:, t]
        z = (z - z.mean()) / z.std()
    z = np.asarray(z, dtype=np.float64)

    mu = nu_mu * (1.0 + 2.0 * rng.random(n_cells))
    kern = (np.arange(1, trunc + 1) + c0) ** (-p)
    # `branching` is the expected number of direct offspring per event; it MUST be
    # < 1 or the simulation is supercritical and the intensity runs away.
    assert 0.0 <= branching < 1.0, "supercritical ETAS: branching must be < 1"
    K = branching / kern.sum()
    A = np.zeros((n_cells, n_days + trunc + 1))
    y = np.zeros((n_cells, n_days), dtype=np.int64)
    for t in range(n_days):
        lam = (mu + K * (G @ A[:, t])) * np.exp(beta_true * z[:, t])
        n = rng.poisson(lam)
        y[:, t] = n
        hit = np.nonzero(n)[0]
        if hit.size:
            A[hit, t + 1:t + 1 + trunc] += n[hit][:, None] * kern[None, :]

    cell_idx, day_idx = np.nonzero(y)
    rep = y[cell_idx, day_idx]
    ev_cell = np.repeat(cell_idx, rep)
    ev_day = np.repeat(day_idx, rep)
    meta = dict(n_days=n_days, n_explore_days=int(round(explore_frac * n_days)),
                mag_floor=4.5, dlat=1.0, dlon=1.0)
    ctx = design._make_ctx(g, ev_cell, ev_day, np.full(ev_cell.size, 4.6), meta,
                           verbose=False)
    return ctx, ctx.day_counts(4.5), z.astype(np.float32)


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
