"""Grid spec and (lat, lon, t) -> cell/day indexing.

v1 grid: regular lat/lon cells (default 1 deg) x 1 day.

DOMAIN RESTRICTION (declared, v1): a global 1x1 deg grid is 64,800 cells x ~11,500
days = 7.5e8 space-time cells, which is not tractable and is ~98% ocean/craton with
zero rate. The forecast domain is therefore restricted to ACTIVE cells: cells that
contain at least one catalog event during the EXPLORATION period. Holdout events
that land outside the domain are counted and reported (they are, correctly, a cost
charged to the domain definition, not silently dropped).
"""

from __future__ import annotations

import numpy as np

EARTH_R_KM = 6371.0088


class Grid:
    def __init__(self, dlat: float = 1.0, dlon: float = 1.0):
        if dlat <= 0 or dlon <= 0:
            raise ValueError("cell size must be positive")
        self.dlat = float(dlat)
        self.dlon = float(dlon)
        # populated by build_domain
        self.cell_lat = None   # (n_cells,) lower-left lat index * dlat
        self.cell_lon = None
        self.key_to_idx = None
        self.n_cells = 0

    # --- integer cell keys -------------------------------------------------
    def ij(self, lat, lon):
        lat = np.asarray(lat, dtype="float64")
        lon = np.asarray(lon, dtype="float64")
        lon = ((lon + 180.0) % 360.0) - 180.0
        i = np.floor(lat / self.dlat).astype(np.int64)
        j = np.floor(lon / self.dlon).astype(np.int64)
        return i, j

    @staticmethod
    def _key(i, j):
        return (i.astype(np.int64) + 1000) * 100000 + (j.astype(np.int64) + 1000)

    def keys(self, lat, lon):
        i, j = self.ij(lat, lon)
        return self._key(i, j)

    # --- domain ------------------------------------------------------------
    def build_domain(self, lat, lon):
        """Define active cells from the given (training) events."""
        k = self.keys(lat, lon)
        uniq = np.unique(k)
        if uniq.size == 0:
            raise ValueError("bulk transform produced 0 cells (build_domain)")
        self.cell_key = uniq
        self.key_to_idx = {int(v): n for n, v in enumerate(uniq)}
        i = (uniq // 100000) - 1000
        j = (uniq % 100000) - 1000
        self.cell_lat = (i + 0.5) * self.dlat   # cell centers
        self.cell_lon = (j + 0.5) * self.dlon
        self.n_cells = int(uniq.size)
        return self.n_cells

    def cell_index(self, lat, lon):
        """Map events to cell index; -1 for events outside the active domain."""
        k = self.keys(lat, lon)
        idx = np.searchsorted(self.cell_key, k)
        idx_clipped = np.clip(idx, 0, self.cell_key.size - 1)
        hit = self.cell_key[idx_clipped] == k
        out = np.where(hit, idx_clipped, -1).astype(np.int64)
        return out

    # --- geometry ----------------------------------------------------------
    def cell_area_km2(self):
        h = EARTH_R_KM * np.deg2rad(self.dlat)
        w = EARTH_R_KM * np.deg2rad(self.dlon) * np.cos(np.deg2rad(self.cell_lat))
        return np.abs(h * w)

    def pair_distance_km(self):
        """(n_cells, n_cells) great-circle distance between cell centers, float32."""
        la = np.deg2rad(self.cell_lat)
        lo = np.deg2rad(self.cell_lon)
        x = np.cos(la) * np.cos(lo)
        y = np.cos(la) * np.sin(lo)
        z = np.sin(la)
        v = np.stack([x, y, z], axis=1)
        c = np.clip(v @ v.T, -1.0, 1.0)
        d = EARTH_R_KM * np.arccos(c)
        np.fill_diagonal(d, 0.0)   # arccos(~1) leaves ~1e-4 km of float noise otherwise
        return d.astype(np.float32)

    def pair_azimuth_rad(self):
        """(n_cells, n_cells) initial bearing FROM row cell TO column cell, float32."""
        la1 = np.deg2rad(self.cell_lat)[:, None]
        la2 = np.deg2rad(self.cell_lat)[None, :]
        dlo = np.deg2rad(self.cell_lon[None, :] - self.cell_lon[:, None])
        y = np.sin(dlo) * np.cos(la2)
        x = np.cos(la1) * np.sin(la2) - np.sin(la1) * np.cos(la2) * np.cos(dlo)
        return np.arctan2(y, x).astype(np.float32)


def day_index(days_float: np.ndarray, n_days: int) -> np.ndarray:
    d = np.floor(np.asarray(days_float, dtype="float64")).astype(np.int64)
    return np.clip(d, 0, n_days - 1)
