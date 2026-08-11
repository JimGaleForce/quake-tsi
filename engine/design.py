"""Precompute + cache the per-cell/per-day design arrays.

Cache lives in engine/out/cache/design_<key>.npz. Key is a hash of
(data files + mtimes, grid size, mag floor, split fraction).

The heavy object is a (n_cells, n_days) day-count array; every covariate helper is
expressed as (neighbourhood matrix) @ (trailing-window cumsum of that array), which
is a couple of float32 GEMMs rather than a python loop over cells.
"""

from __future__ import annotations

import hashlib
import json
import os
import time

import numpy as np

from . import datasets, grid as gridmod

CACHE_DIR = os.path.join("engine", "out", "cache")


def _cache_key(data_dir, dlat, dlon, mag_floor, explore_frac):
    files = sorted(
        (os.path.basename(p), int(os.path.getmtime(p)), os.path.getsize(p))
        for p in __import__("glob").glob(os.path.join(data_dir, "*.csv"))
    )
    blob = json.dumps(
        {"files": files, "dlat": dlat, "dlon": dlon, "floor": mag_floor,
         "explore_frac": explore_frac, "v": 1},
        sort_keys=True,
    )
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


class EngineContext:
    """Everything a covariate plugin is allowed to touch."""

    def __init__(self, **kw):
        self.__dict__.update(kw)

    # ---- plugin-facing helpers -------------------------------------------
    def day_counts(self, m_min: float) -> np.ndarray:
        """(n_cells, n_days) float32 count of catalog events with mag >= m_min."""
        key = round(float(m_min), 3)
        c = self._daycount_cache.get(key)
        if c is None:
            sel = self.ev_mag >= m_min
            c = np.zeros((self.n_cells, self.n_days), dtype=np.float32)
            np.add.at(c, (self.ev_cell[sel], self.ev_day[sel]), 1.0)
            self._daycount_cache[key] = c
        return c

    @staticmethod
    def trailing_sum(counts: np.ndarray, days: int) -> np.ndarray:
        """S[:, t] = sum of counts over days [t-days, t-1]  (STRICTLY causal).

        z[:, t] therefore uses only events on days strictly before t.
        """
        days = int(days)
        if days < 1:
            raise ValueError("window must be >= 1 day")
        n_days = counts.shape[1]
        cs = np.zeros((counts.shape[0], n_days + 1), dtype=np.float64)
        np.cumsum(counts, axis=1, out=cs[:, 1:])
        t = np.arange(n_days)
        hi = cs[:, t]                                 # cumulative through day t-1
        lo = cs[:, np.maximum(t - days, 0)]
        return (hi - lo).astype(np.float32)

    def neighbour_matrix(self, radius_km: float, r_inner_km: float = 0.0) -> np.ndarray:
        """(n_cells, n_cells) float32 0/1 mask of cells within an annulus."""
        key = (round(float(r_inner_km), 3), round(float(radius_km), 3))
        w = self._nbr_cache.get(key)
        if w is None:
            d = self.dist_km
            w = ((d >= r_inner_km) & (d <= radius_km)).astype(np.float32)
            self._nbr_cache[key] = w
        return w

    def past_counts(self, m_min: float, radius_km: float, days: int,
                    r_inner_km: float = 0.0) -> np.ndarray:
        """(n_cells, n_days) count of M>=m_min events within radius, in past `days`.

        Strictly causal: day t sees only days < t.
        """
        s = self.trailing_sum(self.day_counts(m_min), days)
        w = self.neighbour_matrix(radius_km, r_inner_km)
        return (w @ s).astype(np.float32)

    def past_azimuth_resultant(self, m_min: float, r_inner_km: float,
                               r_outer_km: float, days: int):
        """Circular resultant of azimuths cell->source over an annulus.

        v1 approximation (declared caveat): the azimuth is taken from the target cell
        centre to the SOURCE CELL centre, not to the individual hypocentre. At 1 deg
        cells with a 100-600 km annulus that is a <~10 deg bearing error; it does NOT
        resolve fault structure. Returns (R, N) both (n_cells, n_days).
        """
        s = self.trailing_sum(self.day_counts(m_min), days)
        mask = self.neighbour_matrix(r_outer_km, r_inner_km)
        az = self.az_rad
        n = (mask @ s).astype(np.float32)
        cx = ((mask * np.cos(az)) @ s).astype(np.float32)
        cy = ((mask * np.sin(az)) @ s).astype(np.float32)
        r = np.sqrt(cx * cx + cy * cy) / np.maximum(n, 1.0)
        return r.astype(np.float32), n


def build_design(data_dir=datasets.DEFAULT_DATA_DIR, dlat=1.0, dlon=1.0,
                 mag_floor=datasets.CATALOG_MAG_FLOOR, explore_frac=0.70,
                 rebuild=False, verbose=True, cache_dir=CACHE_DIR):
    """Build (or load) the cached design, then return an EngineContext."""
    key = _cache_key(data_dir, dlat, dlon, mag_floor, explore_frac)
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"design_{key}.npz")

    if os.path.exists(path) and not rebuild:
        if verbose:
            print(f"design cache HIT  {path}")
        z = np.load(path, allow_pickle=True)
        meta = json.loads(str(z["meta"]))
        g = gridmod.Grid(meta["dlat"], meta["dlon"])
        g.cell_key = z["cell_key"]
        g.cell_lat = z["cell_lat"]
        g.cell_lon = z["cell_lon"]
        g.n_cells = int(g.cell_key.size)
        g.key_to_idx = None
        return _make_ctx(g, z["ev_cell"], z["ev_day"], z["ev_mag"], meta, verbose=verbose)

    t_start = time.time()
    if verbose:
        print(f"design cache MISS {path} -- building")
    cat, rep = datasets.load_catalog(data_dir, mag_floor)
    if verbose:
        for line in rep.lines():
            print(line)

    days_f, lat, lon, mag, t0 = datasets.catalog_arrays(cat)
    n_days = int(np.floor(days_f.max())) + 1
    day = gridmod.day_index(days_f, n_days)

    n_explore_days = int(round(explore_frac * n_days))
    if not (0 < n_explore_days < n_days):
        raise ValueError("explore_frac gives a degenerate split")

    g = gridmod.Grid(dlat, dlon)
    train = day < n_explore_days
    n_cells = g.build_domain(lat[train], lon[train])
    cell = g.cell_index(lat, lon)

    n_out = int((cell < 0).sum())
    n_out_holdout = int(((cell < 0) & ~train).sum())
    keep = cell >= 0
    if keep.sum() == 0:
        raise ValueError("bulk transform produced 0 in-domain events")

    ev_cell = cell[keep].astype(np.int32)
    ev_day = day[keep].astype(np.int32)
    ev_mag = mag[keep].astype(np.float32)

    # --- Invariants (SPEC "Invariants"), counted and loud ------------------
    distinct_cells = int(np.unique(ev_cell).size)
    assert distinct_cells > 100, f"only {distinct_cells} distinct occupied cells"
    assert np.isfinite(ev_mag).all() and np.isfinite(ev_day).all()
    assert ev_day.min() >= 0 and ev_day.max() < n_days

    meta = dict(
        dlat=dlat, dlon=dlon, mag_floor=mag_floor, explore_frac=explore_frac,
        n_days=n_days, n_explore_days=n_explore_days, t0=str(t0),
        n_events_total=int(rep["n_events"]), n_events_domain=int(keep.sum()),
        n_events_outside_domain=n_out, n_holdout_events_outside_domain=n_out_holdout,
        n_cells=n_cells, distinct_occupied_cells=distinct_cells,
        load_report={k: v for k, v in rep.items() if k != "per_file"},
        cache_key=key,
    )

    np.savez_compressed(
        path, cell_key=g.cell_key, cell_lat=g.cell_lat, cell_lon=g.cell_lon,
        ev_cell=ev_cell, ev_day=ev_day, ev_mag=ev_mag, meta=json.dumps(meta),
    )
    if verbose:
        print(f"design cache WROTE {path}  ({os.path.getsize(path)/1e6:.1f} MB, "
              f"{time.time()-t_start:.1f}s)")
    return _make_ctx(g, ev_cell, ev_day, ev_mag, meta, verbose=verbose)


def _make_ctx(g, ev_cell, ev_day, ev_mag, meta, verbose=True):
    ctx = EngineContext(
        grid=g,
        n_cells=int(g.n_cells),
        n_days=int(meta["n_days"]),
        n_explore_days=int(meta["n_explore_days"]),
        ev_cell=np.asarray(ev_cell, dtype=np.int64),
        ev_day=np.asarray(ev_day, dtype=np.int64),
        ev_mag=np.asarray(ev_mag, dtype=np.float64),
        meta=meta,
        _daycount_cache={},
        _nbr_cache={},
    )
    ctx.dist_km = g.pair_distance_km()
    ctx.az_rad = g.pair_azimuth_rad()

    counts = ctx.day_counts(meta["mag_floor"])
    assert np.isfinite(counts).all(), "design matrix has NaN/Inf"
    tot = float(counts.sum())
    assert abs(tot - ctx.ev_cell.size) < 0.5, f"count array {tot} != events {ctx.ev_cell.size}"

    if verbose:
        print("design:")
        print(f"  grid                 = {g.dlat} deg x {g.dlon} deg x 1 day")
        print(f"  n_cells (active)     = {ctx.n_cells}")
        print(f"  n_days               = {ctx.n_days}  (explore {ctx.n_explore_days}, "
              f"holdout {ctx.n_days - ctx.n_explore_days})")
        print(f"  space-time cells     = {ctx.n_cells * ctx.n_days:,}")
        print(f"  events in domain     = {ctx.ev_cell.size}")
        print(f"  events outside domain= {meta['n_events_outside_domain']} "
              f"(holdout-period share: {meta['n_holdout_events_outside_domain']})")
        print(f"  distinct occupied    = {meta['distinct_occupied_cells']}")
        print(f"  count-array sum      = {tot:.0f}  (== events in domain: OK)")
        print(f"  finite check         = OK (no NaN/Inf)")
    return ctx
