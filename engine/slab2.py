"""SLAB2 SUBDUCTION INTERFACE GEOMETRY: a fault plane for every event, for free.

WHY THIS EXISTS. Every angular statistic this program has run on real catalogues used
an ABSOLUTE, GEOGRAPHIC bearing, because a per-event fault plane was not available. The
world scan said so in its own declaration: "an effect organised relative to the FAULT
rather than to geography is DILUTED here." That is the single largest known weakness in
the scan, and the physics says fault-relative is exactly where an orientation effect
should live -- a fault plane does not care which way is north.

Kepler's K-101 named the way out: the usual constraint is "you need a focal mechanism
per event, so you are capped at GCMT", and that constraint is an artefact of insisting
the plane come from the EVENT. Take it from the SLAB MODEL instead. Slab2 (Hayes et al.
2018, USGS, doi 10.5066/F7PV6JNV) gives interface depth, strike and dip on a 0.05-degree
grid for 27 subduction zones worldwide, with ZERO free parameters and ZERO catalogue
input into the geometry.

WHAT THIS BUYS AND WHAT IT COSTS.

  BUYS: a declared strike and dip for any event within a modelled subduction zone,
  which turns every bearing statistic from geographic into fault-relative, and makes
  resolved Coulomb computable outside the one hand-declared Alaska geometry.

  COSTS, and these are real and must travel with any result:
    * The plane is the INTERFACE, so an event on a splay fault, an outer-rise normal
      fault, or an intraslab plane inherits a geometry that is not its own. The depth
      misfit `|event depth - interface depth|` is returned for exactly this reason, and
      a declared misfit cut is how an analysis restricts itself to plausible interface
      events.
    * Slab2 models SUBDUCTION ZONES ONLY. California, Iceland, Turkey, Iran and the
      Himalaya are wholly or partly outside it, so a fault-relative scan covers a
      SUBSET of the world scan's regions and the difference must be reported rather
      than quietly absorbed.
    * Assignment is nearest-node on the published grid. No interpolation across the
      slab edge, and NaN outside the modelled surface, so an event off the model
      returns nothing rather than an extrapolation.

THE FILES. The distribution ships both `.grd` (HDF5/netCDF-4) and `.xyz` (plain
lon,lat,value text). This module reads the `.xyz`, deliberately: it needs no HDF5
dependency at all, and a plain-text grid is auditable by eye. Grids are cached to
compressed `.npz` on first read because parsing 600k lines per field per slab is slow
and the parse is deterministic.

LONGITUDE CONVENTION, WHICH IS THE EASIEST THING TO GET WRONG HERE. Slab2 grids are on
0..360. ComCat is on -180..180. `_wrap360` is applied to every query and there is a
test that an event at -160 and the same event at 200 return identical geometry.
"""

from __future__ import annotations

import glob
import os
import re

import numpy as np

SLAB2_DIR = os.path.join("data", "slab2", "Slab2Distribute_Mar2018", "Slab2_TXT")
CACHE_DIR = os.path.join("data", "slab2", "_cache")

CITATION = ("Hayes, G.P. et al. (2018), Slab2 - A Comprehensive Subduction Zone "
            "Geometry Model, USGS data release, doi:10.5066/F7PV6JNV. Exogenous to "
            "the earthquake catalogue: zero free parameters, zero catalogue input "
            "into the geometry.")

SCOPE = ("Slab2 gives the SUBDUCTION INTERFACE. An event on a splay, outer-rise or "
         "intraslab plane inherits a geometry that is not its own, which is why "
         "depth_misfit_km is returned and why any analysis must declare a misfit cut. "
         "Non-subduction regions are NOT modelled and are returned as unassigned "
         "rather than approximated.")

_FIELDS = ("dep", "str", "dip")
_cache = {}


def _wrap360(lon):
    return np.mod(np.asarray(lon, dtype=np.float64), 360.0)


def available_slabs():
    pat = os.path.join(SLAB2_DIR, "*_slab2_dep_*.xyz")
    out = []
    for p in sorted(glob.glob(pat)):
        b = os.path.basename(p)
        if b.startswith("._"):
            continue
        m = re.match(r"([a-z]{3})_slab2_dep_", b)
        if m:
            out.append(m.group(1))
    return out


def _field_path(code, field):
    pat = os.path.join(SLAB2_DIR, "%s_slab2_%s_*.xyz" % (code, field))
    hits = [p for p in sorted(glob.glob(pat))
            if not os.path.basename(p).startswith("._")]
    if not hits:
        raise FileNotFoundError("no Slab2 %s grid for %r" % (field, code))
    return hits[0]


def _parse_xyz(path):
    """lon, lat, value -> (lon axis, lat axis, 2-D value array). NaN preserved."""
    raw = np.loadtxt(path, delimiter=",", dtype=np.float64)
    lon, lat, val = raw[:, 0], raw[:, 1], raw[:, 2]
    ulon = np.unique(lon)
    ulat = np.unique(lat)
    grid = np.full((ulat.size, ulon.size), np.nan, dtype=np.float32)
    i = np.searchsorted(ulat, lat)
    j = np.searchsorted(ulon, lon)
    grid[i, j] = val
    return ulon, ulat, grid


def load_slab(code):
    """Grids for one slab, cached in memory and on disk."""
    if code in _cache:
        return _cache[code]
    os.makedirs(CACHE_DIR, exist_ok=True)
    cpath = os.path.join(CACHE_DIR, "%s.npz" % code)
    if os.path.exists(cpath):
        z = np.load(cpath)
        rec = {"code": code, "lon": z["lon"], "lat": z["lat"],
               **{f: z[f] for f in _FIELDS}}
    else:
        rec = {"code": code}
        for f in _FIELDS:
            lon, lat, g = _parse_xyz(_field_path(code, f))
            rec["lon"], rec["lat"], rec[f] = lon, lat, g
        np.savez_compressed(cpath, lon=rec["lon"], lat=rec["lat"],
                            **{f: rec[f] for f in _FIELDS})
    _cache[code] = rec
    return rec


def _lookup(rec, lon360, lat):
    """Nearest-node lookup. Returns NaN outside the grid or on a NaN node."""
    lo, la = rec["lon"], rec["lat"]
    out = {f: np.full(np.shape(lon360), np.nan, dtype=np.float64) for f in _FIELDS}
    inside = ((lon360 >= lo[0]) & (lon360 <= lo[-1])
              & (lat >= la[0]) & (lat <= la[-1]))
    if not np.any(inside):
        return out
    j = np.clip(np.searchsorted(lo, lon360[inside]), 0, lo.size - 1)
    i = np.clip(np.searchsorted(la, lat[inside]), 0, la.size - 1)
    # searchsorted gives the insertion point; step back where that node is closer
    j = np.where((j > 0) & (np.abs(lo[np.maximum(j - 1, 0)] - lon360[inside])
                            < np.abs(lo[j] - lon360[inside])), j - 1, j)
    i = np.where((i > 0) & (np.abs(la[np.maximum(i - 1, 0)] - lat[inside])
                            < np.abs(la[i] - lat[inside])), i - 1, i)
    for f in _FIELDS:
        out[f][inside] = rec[f][i, j]
    return out


def assign(lat, lon, depth_km=None, codes=None):
    """Interface strike, dip and depth for each (lat, lon), plus the depth misfit.

    Where more than one slab covers a point, the one whose interface is CLOSEST IN
    DEPTH to the event is chosen when a depth is supplied, and otherwise the shallowest
    modelled interface. The rule is declared here rather than left to file order,
    because overlapping models at trench junctions are common and a first-match rule
    would make the answer depend on a directory listing.
    """
    lat = np.atleast_1d(np.asarray(lat, dtype=np.float64))
    lon360 = _wrap360(np.atleast_1d(lon))
    dep_ev = (np.atleast_1d(np.asarray(depth_km, dtype=np.float64))
              if depth_km is not None else None)
    n = lat.size
    best = {f: np.full(n, np.nan) for f in _FIELDS}
    best_code = np.full(n, "", dtype=object)
    best_score = np.full(n, np.inf)

    for code in (codes or available_slabs()):
        rec = load_slab(code)
        got = _lookup(rec, lon360, lat)
        ok = np.isfinite(got["dep"])
        if not np.any(ok):
            continue
        slab_dep = np.abs(got["dep"])          # Slab2 depths are negative downward
        score = (np.abs(slab_dep - dep_ev) if dep_ev is not None else slab_dep)
        take = ok & (score < best_score)
        if np.any(take):
            for f in _FIELDS:
                best[f][take] = got[f][take]
            best_score[take] = score[take]
            best_code[take] = code
    misfit = (np.abs(np.abs(best["dep"]) - dep_ev) if dep_ev is not None
              else np.full(n, np.nan))
    return {
        "strike_deg": best["str"],
        "dip_deg": best["dip"],
        "interface_depth_km": np.abs(best["dep"]),
        "depth_misfit_km": misfit,
        "slab_code": best_code,
        "assigned": np.isfinite(best["str"]),
        "citation": CITATION,
        "scope": SCOPE,
    }
