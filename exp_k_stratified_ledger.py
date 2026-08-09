"""EXP-K: style-stratified, fault-resolved stress ledger
(PATTERN_PROTOCOL.md, section EXP-K, frozen 2026-08-09).

Upgrades EXP-J's style-blind ledger (exp_j_stress_ledger.py).  Per 0.2-deg cell:

K1  full strain-rate tensor (exx, eyy, exy) recomputed from NGL MIDAS via
    strain_comparison.strain_grid_full (0.1 deg, nanostrain/yr) and aggregated 2x2 to 0.2 deg;
    styleness s = dilatation / (2 * max_shear); frozen classes at +/-0.25.
K2  CFM5.3 fault geometry: distance to the nearest trace point, that trace's strike
    (CFM segment header strike; local trace-derived strike reported alongside),
    fault-resolved shear loading rate, ON-FAULT if distance < 10 km.
K3  variable seismogenic thickness H (frozen coarse model) and the ledger chi under BOTH
    flat-11 km and variable H; sensitivity of the SILENT-LOADING list reported.
K4  strata = style x on/off-fault (6); Kruskal-Wallis of log10(chi) across styles;
    within-stratum persistence (Spearman log-chi train vs test).
K5  UNEXPLAINED-SILENT list = variable-H SILENT-LOADING cells minus SAF-creeping-proximal
    minus geothermal-proximal cells.  This is the deliverable hazard-candidate list.

Conventions reused verbatim from EXP-J: 0.2 deg cells aligned to -122.0+0.2k / 31.5+0.2k,
mu = 30 GPa, haversine-corrected cell area, Hanks-Kanamori moment from
SCSN_original_catalog.txt M >= 2.5, train < 2010 / test >= 2010, chi classes
SILENT-LOADING / COUPLED / OVERSHOOT with identical thresholds, +0.5 ev/yr rate regularizer.

Outputs results_exp_k.json and maps/exp_k_stratified.png.  Reads only; modifies nothing.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal, spearmanr

import tsi_map

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT_MAPS = HERE / "maps"

# ---------------------------------------------------------------- frozen params (EXP-J inherited)
BOX = {"lat": (31.5, 38.0), "lon": (-122.0, -113.5)}
CELL = 0.2
MC = 2.5
MAG_BIN = 0.1
MU = 30e9                        # Pa
H_FLAT = 11e3                    # m, EXP-J's flat seismogenic thickness
TRAIN_END = pd.Timestamp("2010-01-01", tz="UTC")
TRAIN_YEARS = 29.0
MIN_EVENTS_CHI = 20
CHI_SILENT = 0.01
CHI_COUPLED_HI = 1.0
RATE_REGULARIZER = 0.5           # events/yr, as EXP-J J2a
EARTH_R = 6371.0
SAF_CREEP = ((36.0, -120.6), (36.8, -121.6))
SEC_PER_YEAR = 365.25 * 86400.0

# ---------------------------------------------------------------- frozen params (EXP-K new)
STYLE_THRESH = 0.25              # |s| > 0.25 -> TRANSTENSIONAL / CONTRACTIONAL
ON_FAULT_KM = 10.0
H_GEOTHERMAL = 6e3               # m, within GEOTHERMAL_DEG of a tsi_map.GEOTHERMAL field
H_TRANSVERSE = 25e3              # m, inside the Transverse Ranges box
H_ELSEWHERE = 11e3               # m
GEOTHERMAL_DEG = 0.35            # degrees (protocol states "0.35 deg", not km)
TRANSVERSE_BOX = {"lat": (34.0, 34.8), "lon": (-119.8, -117.3)}
SAF_EXCLUDE_KM = 25.0
CFM_TRACES = DATA / "xue_lu_zenodo" / "CFM5.3_traces.lonLat"
CFM_GEOM = DATA / "xue_lu_zenodo" / "CFM5.3_traceslonLat_fault_geometry.txt"

STYLES = ["CONTRACTIONAL", "STRIKE-SLIP", "TRANSTENSIONAL"]
CLASSES = ["UNCLASSIFIED", "COUPLED", "OVERSHOOT", "SILENT-LOADING"]


# ---------------------------------------------------------------- helpers
def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1, la2, lo2 = map(np.radians, [lat1, lon1, lat2, lon2])
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * EARTH_R * np.arcsin(np.sqrt(a))


def load_catalog(fname):
    """Auto-detect column order (raw vs declustered) exactly like xue_lu_crosstest.load_catalog."""
    raw_cols = ["yr", "mo", "dy", "hr", "mi", "sec", "eid", "lat", "lon", "depth", "mag"]
    dec_cols = ["yr", "mo", "dy", "hr", "mi", "sec", "lat", "lon", "depth", "mag", "eid"]
    probe = pd.read_csv(DATA / fname, sep=r"\s+", header=None, nrows=1000)
    cols = raw_cols if probe[6].abs().max() > 90 else dec_cols
    df = pd.read_csv(DATA / fname, sep=r"\s+", header=None, names=cols)
    assert df.lat.abs().max() <= 90 and df.lon.abs().max() <= 180, f"column detection failed: {fname}"
    ts = pd.to_datetime(
        dict(year=df.yr, month=df.mo, day=df.dy, hour=df.hr, minute=df.mi, second=0), utc=True
    ) + pd.to_timedelta(df["sec"].astype(float), unit="s")
    return pd.DataFrame({"t": ts, "lat": df.lat, "lon": df.lon, "depth": df.depth, "mag": df.mag})


def dist_to_segment_km(lat, lon, seg):
    """Nearest-point distance to a great-circle segment, by dense sampling (>=1 km spacing)."""
    (a_lat, a_lon), (b_lat, b_lon) = seg
    n = max(200, int(haversine_km(a_lat, a_lon, b_lat, b_lon) * 5))
    f = np.linspace(0.0, 1.0, n)
    slat = a_lat + f * (b_lat - a_lat)
    slon = a_lon + f * (b_lon - a_lon)
    lat = np.asarray(lat)[:, None]
    lon = np.asarray(lon)[:, None]
    return haversine_km(lat, lon, slat[None, :], slon[None, :]).min(axis=1)


def circ_diff_deg_180(a, b):
    """Smallest angular difference between two strikes, modulo 180 deg (planes, not vectors)."""
    d = np.abs((np.asarray(a) - np.asarray(b)) % 180.0)
    return np.minimum(d, 180.0 - d)


# ---------------------------------------------------------------- K1: tensor + style
def load_tensor_grid():
    """Import strain_comparison by path (its module-level code needs no CLI) and call
    strain_grid_full(midas, ...) exactly as its main() does.  Returns 0.1 deg grids."""
    spec = importlib.util.spec_from_file_location("sc", HERE / "strain_comparison.py")
    sc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sc)
    midas = tsi_map.load_midas()
    lats, lons, exx, eyy, exy, n_sta = sc.strain_grid_full(
        midas, midas.ve * 1000, midas.vn * 1000, midas.se * 1000, midas.sn * 1000)
    return lats, lons, exx, eyy, exy, int(n_sta)


def aggregate_2x2(glats, glons, fields):
    """Aggregate 0.1 deg node fields to 0.2 deg cells (mean of the finite nodes in each cell),
    with cell edges aligned to the EXP-J convention (-122.0 + 0.2k / 31.5 + 0.2k)."""
    lat_edges = np.arange(BOX["lat"][0], BOX["lat"][1] - 1e-9, CELL)
    lon_edges = np.arange(BOX["lon"][0], BOX["lon"][1] - 1e-9, CELL)
    nlat, nlon = len(lat_edges), len(lon_edges)
    ilat = np.floor((glats - BOX["lat"][0]) / CELL + 1e-6).astype(int)
    ilon = np.floor((glons - BOX["lon"][0]) / CELL + 1e-6).astype(int)

    ok = np.ones_like(fields[0], dtype=bool)
    for f in fields:
        ok &= np.isfinite(f)
    sums = [np.zeros((nlat, nlon)) for _ in fields]
    cnt = np.zeros((nlat, nlon))
    for a in range(len(glats)):
        if not (0 <= ilat[a] < nlat):
            continue
        row_ok = ok[a] & (ilon >= 0) & (ilon < nlon)
        for s, f in zip(sums, fields):
            np.add.at(s[ilat[a]], ilon[row_ok], f[a][row_ok])
        np.add.at(cnt[ilat[a]], ilon[row_ok], 1.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = [np.where(cnt > 0, s / np.maximum(cnt, 1), np.nan) for s in sums]
    return lat_edges, lon_edges, out, cnt


def build_cells(lat_edges, lon_edges, exx_c, eyy_c, exy_c, cnt):
    rows = []
    for i in range(len(lat_edges)):
        for j in range(len(lon_edges)):
            if not np.isfinite(exx_c[i, j]):
                continue
            lat0, lon0 = lat_edges[i], lon_edges[j]
            latc, lonc = lat0 + CELL / 2, lon0 + CELL / 2
            dy_m = CELL * 111_320.0
            dx_m = CELL * 111_320.0 * np.cos(np.radians(latc))   # haversine-corrected area
            rows.append(dict(i=i, j=j, lat_c=float(latc), lon_c=float(lonc),
                             area_m2=float(dx_m * dy_m),
                             exx=float(exx_c[i, j]), eyy=float(eyy_c[i, j]),
                             exy=float(exy_c[i, j]), n_nodes=int(cnt[i, j])))
    led = pd.DataFrame(rows)
    led["dil"] = led.exx + led.eyy
    led["max_shear"] = np.hypot((led.exx - led.eyy) / 2.0, led.exy)
    with np.errstate(invalid="ignore", divide="ignore"):
        led["styleness"] = led.dil / (2.0 * led.max_shear)
    led["style"] = np.where(led.styleness > STYLE_THRESH, "TRANSTENSIONAL",
                            np.where(led.styleness < -STYLE_THRESH, "CONTRACTIONAL", "STRIKE-SLIP"))
    return led


# ---------------------------------------------------------------- K2: CFM fault geometry
def load_cfm():
    """Parse the CFM5.3 geometry file: segments headed by '> strike dip rake', then
    'lon lat elev' lines (INVENTORY.md and xue_lu_crosstest.load_fault_points agree on this).
    Fault NAMES come from the parallel CFM5.3_traces.lonLat file, whose segment headers are
    '> "name"' and whose point lines are identical.  Returns a dict of flat arrays."""
    seg_strike, seg_dip, seg_rake = [], [], []
    lon, lat, seg_id = [], [], []
    bad_headers = 0
    with open(CFM_GEOM, encoding="utf-8") as f:
        for line in f:
            p = line.split()
            if not p or p[0].startswith("#"):
                continue
            if p[0] == ">":
                try:
                    seg_strike.append(float(p[1])); seg_dip.append(float(p[2]))
                    seg_rake.append(float(p[3]))
                except (IndexError, ValueError):
                    seg_strike.append(np.nan); seg_dip.append(np.nan); seg_rake.append(np.nan)
                    bad_headers += 1
            elif seg_strike:
                try:
                    x, y = float(p[0]), float(p[1])
                except (IndexError, ValueError):
                    continue
                lon.append(x); lat.append(y); seg_id.append(len(seg_strike) - 1)

    names, tr_lon, tr_lat, tr_seg = [], [], [], []
    with open(CFM_TRACES, encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if s.startswith(">"):
                names.append(s[1:].strip().strip('"'))
            elif names:
                p = s.split()
                try:
                    tr_lon.append(float(p[0])); tr_lat.append(float(p[1]))
                except (IndexError, ValueError):
                    continue
                tr_seg.append(len(names) - 1)

    geom = dict(
        strike=np.array(seg_strike), dip=np.array(seg_dip), rake=np.array(seg_rake),
        lon=np.array(lon), lat=np.array(lat), seg=np.array(seg_id, dtype=int),
        names=names, n_segments=len(seg_strike), n_points=len(lon),
        bad_headers=bad_headers,
        traces_n_segments=len(names), traces_n_points=len(tr_lon),
    )
    # do the two files describe the same traces?  (they should: same tool, same order)
    same = (len(names) == len(seg_strike) and len(tr_lon) == len(lon)
            and np.allclose(np.array(tr_lon), geom["lon"], atol=1e-9)
            and np.allclose(np.array(tr_lat), geom["lat"], atol=1e-9))
    geom["traces_match_geometry"] = bool(same)
    geom["usable"] = bool(np.isfinite(geom["strike"]).all() and geom["n_points"] > 1000)

    # LOCAL strike at every trace point, from the neighbouring points of the same segment.
    # Azimuth measured clockwise from north: atan2(dEast, dNorth), reduced modulo 180 deg.
    n = geom["n_points"]
    idx = np.arange(n)
    prev = np.maximum(idx - 1, 0)
    nxt = np.minimum(idx + 1, n - 1)
    prev = np.where(geom["seg"][prev] == geom["seg"], prev, idx)
    nxt = np.where(geom["seg"][nxt] == geom["seg"], nxt, idx)
    lat_m = np.radians((geom["lat"][prev] + geom["lat"][nxt]) / 2.0)
    d_east = (geom["lon"][nxt] - geom["lon"][prev]) * 111.32 * np.cos(lat_m)
    d_north = (geom["lat"][nxt] - geom["lat"][prev]) * 110.57
    with np.errstate(invalid="ignore"):
        loc = np.degrees(np.arctan2(d_east, d_north)) % 180.0
    loc[(prev == nxt)] = np.nan          # isolated single-point segment
    geom["local_strike"] = loc
    return geom


def attach_fault_geometry(led, geom, chunk=200):
    """Nearest CFM trace point per cell centre (haversine, chunked)."""
    clat = led.lat_c.to_numpy()
    clon = led.lon_c.to_numpy()
    best_d = np.full(len(led), np.inf)
    best_k = np.zeros(len(led), dtype=int)
    for a in range(0, len(led), chunk):
        b = min(a + chunk, len(led))
        d = haversine_km(clat[a:b, None], clon[a:b, None],
                         geom["lat"][None, :], geom["lon"][None, :])
        k = d.argmin(axis=1)
        best_k[a:b] = k
        best_d[a:b] = d[np.arange(b - a), k]
    seg = geom["seg"][best_k]
    led = led.copy()
    led["d_fault_km"] = best_d
    led["fault_pt_idx"] = best_k
    led["fault_seg"] = seg
    led["fault_name"] = [geom["names"][s] if s < len(geom["names"]) else "" for s in seg]
    led["strike_cfm"] = geom["strike"][seg] % 180.0
    led["strike_local"] = geom["local_strike"][best_k]
    led["fault_dip"] = geom["dip"][seg]
    led["fault_rake"] = geom["rake"][seg]
    led["on_fault"] = led.d_fault_km < ON_FAULT_KM
    return led


def resolve_shear(led, strike_col):
    """Shear strain rate resolved onto a vertical plane of the given strike.

    Grid frame: x = EAST, y = NORTH; exx = d(vE)/dx, eyy = d(vN)/dy, exy = 0.5(d(vE)/dy+d(vN)/dx)
    (see tsi_map.strain_grid / strain_comparison.strain_grid_full).

    A fault STRIKE is an azimuth measured CLOCKWISE FROM NORTH.  The grid frame is a standard
    math frame (CCW from +x = east), so the in-plane angle is
        theta = 90 deg - strike_azimuth
    (0 deg strike = due north = +y = theta 90 deg; 90 deg strike = due east = +x = theta 0).

    With unit strike vector s = (cos t, sin t) and horizontal fault normal n = (-sin t, cos t),
    the resolved shear strain rate is
        e_ns = n . E . s = 0.5*(eyy - exx)*sin(2t) + exy*cos(2t)
    which is the frozen protocol formula.  We take |e_ns| (sense of slip is irrelevant to a
    loading-rate magnitude), so the +/-180 deg ambiguity of the strike is harmless.
    Bounded above by max_shear = sqrt(((exx-eyy)/2)^2 + exy^2) -- asserted in main().
    """
    t = np.radians(90.0 - led[strike_col].to_numpy(dtype=float))
    v = 0.5 * (led.eyy.to_numpy() - led.exx.to_numpy()) * np.sin(2 * t) \
        + led.exy.to_numpy() * np.cos(2 * t)
    return np.abs(v)


# ---------------------------------------------------------------- K3: variable H
def assign_H(led):
    """Frozen coarse model: 6 km geothermal-proximal, else 25 km in the Transverse Ranges box,
    else 11 km.  'Proximal' is 0.35 DEGREES (protocol wording), Euclidean in (lat, lon)."""
    names = list(tsi_map.GEOTHERMAL)
    glat = np.array([tsi_map.GEOTHERMAL[n][0] for n in names])
    glon = np.array([tsi_map.GEOTHERMAL[n][1] for n in names])
    dd = np.hypot(led.lat_c.to_numpy()[:, None] - glat[None, :],
                  led.lon_c.to_numpy()[:, None] - glon[None, :])
    led = led.copy()
    led["d_geothermal_deg"] = dd.min(axis=1)
    led["nearest_geothermal"] = [names[k] for k in dd.argmin(axis=1)]
    led["d_geothermal_km"] = haversine_km(
        led.lat_c.to_numpy()[:, None], led.lon_c.to_numpy()[:, None],
        glat[None, :], glon[None, :]).min(axis=1)
    geo = led.d_geothermal_deg < GEOTHERMAL_DEG
    tr = (led.lat_c.between(*TRANSVERSE_BOX["lat"]) & led.lon_c.between(*TRANSVERSE_BOX["lon"]))
    H = np.full(len(led), H_ELSEWHERE)
    H[tr.to_numpy()] = H_TRANSVERSE
    H[geo.to_numpy()] = H_GEOTHERMAL          # geothermal wins (listed first in the protocol)
    led["H_m"] = H
    led["H_zone"] = np.where(geo, "geothermal", np.where(tr, "transverse_ranges", "elsewhere"))
    led["d_saf_creep_km"] = dist_to_segment_km(led.lat_c.to_numpy(), led.lon_c.to_numpy(), SAF_CREEP)
    return led


# ---------------------------------------------------------------- seismicity (EXP-J conventions)
def add_seismicity(led, cat, nlon):
    c = cat[(cat.mag >= MC)
            & cat.lat.between(*BOX["lat"]) & cat.lon.between(*BOX["lon"])].copy()
    nlat_max = int(np.ceil((BOX["lat"][1] - BOX["lat"][0]) / CELL))
    ci = np.floor((c.lat.to_numpy() - BOX["lat"][0]) / CELL + 1e-9).astype(int)
    cj = np.floor((c.lon.to_numpy() - BOX["lon"][0]) / CELL + 1e-9).astype(int)
    keep = (ci >= 0) & (ci < nlat_max) & (cj >= 0) & (cj < nlon)
    c = c.loc[keep].copy()
    c["cell"] = ci[keep] * nlon + cj[keep]
    c["m0"] = 10.0 ** (1.5 * c.mag.to_numpy() + 9.05)          # Hanks-Kanamori, N*m
    c["is_train"] = c.t < TRAIN_END

    t_end = c.t.max()
    test_years = float((t_end - TRAIN_END).total_seconds() / SEC_PER_YEAR)

    led = led.copy()
    led["cell"] = led.i.to_numpy() * nlon + led.j.to_numpy()
    for tag, sub in (("train", c[c.is_train]), ("test", c[~c.is_train])):
        g = sub.groupby("cell")
        agg = pd.DataFrame({f"n_{tag}": g.size(), f"m0_{tag}": g.m0.sum(),
                            f"depth_med_{tag}": g.depth.median()})
        led = led.merge(agg, left_on="cell", right_index=True, how="left")
        led[f"n_{tag}"] = led[f"n_{tag}"].fillna(0).astype(int)
        led[f"m0_{tag}"] = led[f"m0_{tag}"].fillna(0.0)
    led["Mdot_seis_train"] = led.m0_train / TRAIN_YEARS
    led["Mdot_seis_test"] = led.m0_test / test_years
    led["rate_train"] = led.n_train / TRAIN_YEARS
    led["rate_test"] = led.n_test / test_years
    return led, test_years, t_end


def ledger_for_H(led, H_m, shear_col="max_shear"):
    """Mdot_geo = 2 mu H A edot (edot 1/yr = nanostrain/yr * 1e-9); chi and frozen classes."""
    d = led.copy()
    d["Mdot_geo"] = 2.0 * MU * np.asarray(H_m) * d.area_m2 * (d[shear_col].to_numpy() * 1e-9)
    with np.errstate(invalid="ignore", divide="ignore"):
        d["chi_train"] = d.Mdot_seis_train / d.Mdot_geo
        d["chi_test"] = d.Mdot_seis_test / d.Mdot_geo
    q75 = float(np.nanpercentile(d.Mdot_geo, 75))
    top = d.Mdot_geo >= q75
    meas = d.n_train >= MIN_EVENTS_CHI
    cls = np.full(len(d), "UNCLASSIFIED", dtype=object)
    reason = np.full(len(d), "", dtype=object)
    silent_low_chi = top & meas & (d.chi_train < CHI_SILENT)
    silent_lowcnt = top & ~meas
    coupled = meas & (d.chi_train >= CHI_SILENT) & (d.chi_train <= CHI_COUPLED_HI)
    overshoot = meas & (d.chi_train > CHI_COUPLED_HI)
    cls[silent_low_chi.to_numpy()] = "SILENT-LOADING"
    reason[silent_low_chi.to_numpy()] = "low_chi_top_quartile_loading"
    cls[silent_lowcnt.to_numpy()] = "SILENT-LOADING"
    reason[silent_lowcnt.to_numpy()] = "low_count_top_quartile_loading"
    cls[coupled.to_numpy()] = "COUPLED"; reason[coupled.to_numpy()] = "chi_in_[0.01,1]"
    cls[overshoot.to_numpy()] = "OVERSHOOT"; reason[overshoot.to_numpy()] = "chi_gt_1"
    rest = cls == "UNCLASSIFIED"
    reason[rest & (~meas).to_numpy()] = "low_count_not_top_quartile_loading"
    reason[rest & meas.to_numpy()] = "low_chi_not_top_quartile_loading"
    d["cls"] = cls
    d["reason"] = reason
    return d, q75


# ---------------------------------------------------------------- K4: strata
def stratum_stats(led):
    out = {}
    total = 0
    for style in STYLES:
        for onf in (True, False):
            sub = led[(led["style"] == style) & (led.on_fault == onf)]
            key = f"{style}|{'ON-FAULT' if onf else 'OFF-FAULT'}"
            total += len(sub)
            meas = sub[(sub.n_train >= MIN_EVENTS_CHI) & (sub.chi_train > 0)]
            both = sub[(sub.n_train >= MIN_EVENTS_CHI) & (sub.n_test >= MIN_EVENTS_CHI)
                       & (sub.chi_train > 0) & (sub.chi_test > 0)]
            pers = {"n": int(len(both)), "min_events_each_period": MIN_EVENTS_CHI}
            if len(both) >= 3:
                rho, p = spearmanr(np.log10(both.chi_train.to_numpy()),
                                   np.log10(both.chi_test.to_numpy()))
                pers.update({"spearman_rho": (None if not np.isfinite(rho) else float(rho)),
                             "p": (None if not np.isfinite(p) else float(p)),
                             "caveat": ("n < 10: reported for completeness, not interpretable"
                                        if len(both) < 10 else None)})
            else:
                pers["note"] = f"n={len(both)} < 3, Spearman not computed"
            ratio = (sub.rate_test + RATE_REGULARIZER) / (sub.rate_train + RATE_REGULARIZER)
            out[key] = {
                "style": style, "on_fault": bool(onf),
                "n_cells": int(len(sub)),
                "n_measurable_train": int(len(meas)),
                "median_chi_train_measurable": (float(np.median(meas.chi_train)) if len(meas) else None),
                "median_chi_train_all_positive": (float(np.median(sub.chi_train[sub.chi_train > 0]))
                                                  if (sub.chi_train > 0).any() else None),
                "median_Mdot_geo_Nm_per_yr": (float(np.median(sub.Mdot_geo)) if len(sub) else None),
                "median_resolved_shear_nstrain_yr": (float(np.median(sub.resolved_shear))
                                                     if len(sub) else None),
                "median_max_shear_nstrain_yr": (float(np.median(sub.max_shear)) if len(sub) else None),
                "median_n_train": (float(np.median(sub.n_train)) if len(sub) else None),
                "median_rate_ratio_test_over_train_reg0.5": (float(np.median(ratio))
                                                             if len(sub) else None),
                "class_counts": {k: int(v) for k, v in sub.cls.value_counts().items()},
                "persistence": pers,
            }
    return out, total


# ---------------------------------------------------------------- map
def draw_maps(led, lat_edges, lon_edges, unexplained_idx):
    OUT_MAPS.mkdir(exist_ok=True)
    nlat, nlon = len(lat_edges), len(lon_edges)
    ext = [lon_edges[0], lon_edges[-1] + CELL, lat_edges[0], lat_edges[-1] + CELL]

    def grid(values):
        g = np.full((nlat, nlon), np.nan)
        g[led.i.to_numpy(), led.j.to_numpy()] = np.asarray(values, dtype=float)
        return g

    fig, axes = plt.subplots(2, 2, figsize=(15.5, 13))

    # panel 1: style
    scode = grid([STYLES.index(s) for s in led["style"]])
    ax = axes[0, 0]
    cmap1 = matplotlib.colors.ListedColormap(["#b2182b", "#bdbdbd", "#2166ac"])
    im = ax.imshow(scode, origin="lower", extent=ext, aspect="auto", cmap=cmap1,
                   vmin=-0.5, vmax=2.5)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, ticks=range(3))
    cb.ax.set_yticklabels(STYLES, fontsize=8)
    ax.set_title("Style class  s = dilatation / (2 x max shear),  thresholds +/-0.25", fontsize=10)

    # panel 2: fault-resolved loading
    with np.errstate(invalid="ignore", divide="ignore"):
        ax = axes[0, 1]
        im = ax.imshow(grid(led.resolved_shear), origin="lower", extent=ext, aspect="auto",
                       cmap="viridis")
        fig.colorbar(im, ax=ax, shrink=0.85)
        ax.set_title("Fault-resolved shear loading rate onto nearest CFM strike "
                     "(nanostrain/yr)", fontsize=10)

        # panel 3: log chi, variable H
        ax = axes[1, 0]
        v = led.chi_train.to_numpy(dtype=float)
        im = ax.imshow(grid(np.log10(np.where(v > 0, v, np.nan))), origin="lower", extent=ext,
                       aspect="auto", cmap="coolwarm")
        fig.colorbar(im, ax=ax, shrink=0.85)
        ax.set_title("log10 chi = seismic release / loading (train), VARIABLE H", fontsize=10)

    # panel 4: classes + unexplained-silent highlight
    ax = axes[1, 1]
    code = grid([CLASSES.index(c) for c in led.cls])
    cmap = matplotlib.colors.ListedColormap(["#d9d9d9", "#4575b4", "#7b3294", "#d73027"])
    im = ax.imshow(code, origin="lower", extent=ext, aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, ticks=range(4))
    cb.ax.set_yticklabels(CLASSES, fontsize=8)
    u = led.loc[unexplained_idx]
    ax.plot(u.lon_c, u.lat_c, "o", ms=5, mfc="none", mec="lime", mew=1.4,
            label=f"unexplained silent (n={len(u)})")
    for name, (gla, glo) in tsi_map.GEOTHERMAL.items():
        if BOX["lat"][0] <= gla <= BOX["lat"][1] and BOX["lon"][0] <= glo <= BOX["lon"][1]:
            ax.plot(glo, gla, "k^", ms=8)
    ax.plot([SAF_CREEP[0][1], SAF_CREEP[1][1]], [SAF_CREEP[0][0], SAF_CREEP[1][0]],
            "k-", lw=2, label="SAF creeping segment")
    ax.legend(fontsize=8, loc="lower left")
    ax.set_title("Variable-H ledger classes; unexplained-silent cells circled", fontsize=10)

    for ax in axes.ravel():
        ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    fig.suptitle("EXP-K style-stratified fault-resolved ledger - SoCal 0.2 deg, train 1981-2009 "
                 f"(mu={MU/1e9:.0f} GPa, variable H 6/11/25 km, M>={MC})", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_MAPS / "exp_k_stratified.png", dpi=170)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    print("EXP-K style-stratified fault-resolved ledger")

    # ---------------- K1 tensor + style
    glats, glons, exx, eyy, exy, n_sta = load_tensor_grid()
    print(f"  MIDAS stations after QC: {n_sta}; 0.1 deg nodes with tensor: "
          f"{int(np.isfinite(exx).sum())}")
    lat_edges, lon_edges, (exx_c, eyy_c, exy_c), cnt = aggregate_2x2(
        glats, glons, [exx, eyy, exy])
    led = build_cells(lat_edges, lon_edges, exx_c, eyy_c, exy_c, cnt)
    n_cells = len(led)
    assert n_cells > 0, "no 0.2 deg cells with a finite tensor"
    print(f"  cells with finite tensor: {n_cells}")

    # cross-check the recomputed shear against EXP-J's stored 0.1 deg grid
    z = np.load(DATA / "socal_strain_grid.npz")
    ms_ref = z["max_shear_nstrain_yr"]
    ok = np.isfinite(ms_ref) & np.isfinite(exx)
    ms_new = np.hypot((exx - eyy) / 2.0, exy)
    shear_check = {"n_nodes": int(ok.sum()),
                   "max_abs_diff_nstrain_yr": float(np.max(np.abs(ms_new[ok] - ms_ref[ok]))),
                   "median_abs_diff_nstrain_yr": float(np.median(np.abs(ms_new[ok] - ms_ref[ok])))}
    print(f"  recomputed max_shear vs stored EXP-J grid: max|diff| = "
          f"{shear_check['max_abs_diff_nstrain_yr']:.3g} nanostrain/yr")

    style_counts = {k: int(v) for k, v in led["style"].value_counts().items()}
    print(f"  style counts: {style_counts}")

    # ---------------- K2 fault geometry
    geom = load_cfm()
    cfm_used = bool(geom["usable"])
    print(f"  CFM5.3: {geom['n_segments']} segments / {geom['n_points']} trace points; "
          f"header strike parsed for all = {bool(np.isfinite(geom['strike']).all())}; "
          f"traces file matches geometry file = {geom['traces_match_geometry']}")
    led = attach_fault_geometry(led, geom)
    strike_col = "strike_cfm" if cfm_used else "strike_local"
    if not cfm_used:
        print("  !! CFM geometry header strikes unusable - falling back to trace-derived strike")
    led["strike_used_deg"] = led[strike_col]
    led["resolved_shear"] = resolve_shear(led, strike_col)
    led["resolved_shear_local_strike"] = resolve_shear(led, "strike_local")
    slack = led.resolved_shear.to_numpy() - led.max_shear.to_numpy()
    assert np.nanmax(slack) < 1e-6, f"resolved shear exceeds max shear by {np.nanmax(slack)}"
    dstr = circ_diff_deg_180(led.strike_cfm.to_numpy(), led.strike_local.to_numpy())
    strike_agreement = {
        "n": int(np.isfinite(dstr).sum()),
        "median_abs_diff_deg_mod180": float(np.nanmedian(dstr)),
        "frac_within_20deg": float(np.nanmean(dstr[np.isfinite(dstr)] <= 20.0)),
    }
    print(f"  CFM header strike vs local trace strike: median |diff| "
          f"{strike_agreement['median_abs_diff_deg_mod180']:.1f} deg, "
          f"{100*strike_agreement['frac_within_20deg']:.0f}% within 20 deg")
    print(f"  ON-FAULT (<{ON_FAULT_KM:.0f} km) cells: {int(led.on_fault.sum())} / {n_cells}")

    # ---------------- K3 H model + seismicity
    led = assign_H(led)
    print("  H zones: " + str({k: int(v) for k, v in led.H_zone.value_counts().items()}))
    cat = load_catalog("xue_lu_zenodo/SCSN_original_catalog.txt")
    print(f"  catalog events: {len(cat)}  {cat.t.min().date()} -> {cat.t.max().date()}")
    led, test_years, t_end = add_seismicity(led, cat, len(lon_edges))

    flat, q75_flat = ledger_for_H(led, H_FLAT)
    var, q75_var = ledger_for_H(led, led.H_m.to_numpy())
    counts_flat = {k: int(v) for k, v in flat.cls.value_counts().items()}
    counts_var = {k: int(v) for k, v in var.cls.value_counts().items()}
    assert sum(counts_var.values()) == n_cells, "class counts do not sum to ledger cells"
    print(f"  classes flat-11 : {counts_flat}")
    print(f"  classes variable: {counts_var}")

    s_flat = set(flat.index[flat.cls == "SILENT-LOADING"])
    s_var = set(var.index[var.cls == "SILENT-LOADING"])
    jac = len(s_flat & s_var) / len(s_flat | s_var) if (s_flat | s_var) else None
    h_sens = {
        "flat_H_m": H_FLAT,
        "variable_H_model": {"geothermal_within_deg": GEOTHERMAL_DEG, "H_geothermal_m": H_GEOTHERMAL,
                             "transverse_ranges_box": TRANSVERSE_BOX, "H_transverse_m": H_TRANSVERSE,
                             "H_elsewhere_m": H_ELSEWHERE},
        "H_zone_counts": {k: int(v) for k, v in led.H_zone.value_counts().items()},
        "class_counts_flat": counts_flat, "class_counts_variable": counts_var,
        "loading_top_quartile_threshold_flat_Nm_per_yr": q75_flat,
        "loading_top_quartile_threshold_variable_Nm_per_yr": q75_var,
        "n_silent_flat": len(s_flat), "n_silent_variable": len(s_var),
        "n_silent_intersection": len(s_flat & s_var), "n_silent_union": len(s_flat | s_var),
        "jaccard_silent_flat_vs_variable": (None if jac is None else float(jac)),
        "added_by_variable_H": len(s_var - s_flat), "dropped_by_variable_H": len(s_flat - s_var),
    }
    print(f"  silent list: flat {len(s_flat)}, variable {len(s_var)}, "
          f"Jaccard {jac:.3f}" if jac is not None else "  silent list empty")

    # ---------------- K4 strata (primary ledger = variable H)
    strata, strata_total = stratum_stats(var)
    assert strata_total == n_cells, f"strata cells {strata_total} != ledger cells {n_cells}"
    for k, v in strata.items():
        print(f"    {k:28s} n={v['n_cells']:4d} meas={v['n_measurable_train']:3d} "
              f"med_chi={v['median_chi_train_measurable']}")

    meas = var[(var.n_train >= MIN_EVENTS_CHI) & (var.chi_train > 0)]
    groups = [np.log10(meas.chi_train[meas["style"] == s].to_numpy()) for s in STYLES]
    ns = [int(g.size) for g in groups]
    if min(ns) >= 2 and sum(ns) >= 6:
        H_stat, p_kw = kruskal(*[g for g in groups if g.size])
        kw = {"variable": "log10(chi_train), variable-H ledger",
              "min_train_events": MIN_EVENTS_CHI,
              "groups": {s: {"n": int(g.size),
                             "median_log10_chi": (float(np.median(g)) if g.size else None)}
                         for s, g in zip(STYLES, groups)},
              "H": float(H_stat), "p": float(p_kw), "significant_0.05": bool(p_kw < 0.05)}
        print(f"  Kruskal-Wallis log10(chi) across styles: n={ns} H={H_stat:.2f} p={p_kw:.3g}")
    else:
        kw = {"ERROR": f"insufficient cells per style group: {dict(zip(STYLES, ns))}"}
        print(f"  !! Kruskal-Wallis not run: {ns}")

    # ---------------- K5 unexplained-silent list
    silent = var.cls == "SILENT-LOADING"
    near_saf = var.d_saf_creep_km < SAF_EXCLUDE_KM
    near_geo = var.d_geothermal_deg < GEOTHERMAL_DEG
    unexp = silent & ~near_saf & ~near_geo
    assert set(var.index[unexp]) <= s_var, "unexplained-silent is not a subset of silent"
    u = var[unexp].sort_values("Mdot_geo", ascending=False)
    rows = [{
        "lat": float(r.lat_c), "lon": float(r.lon_c),
        "Mdot_geo_Nm_per_yr": float(r.Mdot_geo),
        "Mdot_seis_train_Nm_per_yr": float(r.Mdot_seis_train),
        "max_shear_nstrain_yr": float(r.max_shear),
        "resolved_shear_nstrain_yr": float(r.resolved_shear),
        "n_train": int(r.n_train), "n_test": int(r.n_test),
        "chi_train": float(r.chi_train),
        "chi_test": float(r.chi_test),
        "style": r.style, "styleness": float(r.styleness),
        "H_m": float(r.H_m), "H_zone": r.H_zone,
        "d_fault_km": float(r.d_fault_km),
        "nearest_fault": r.fault_name,
        "nearest_fault_strike_deg": (None if not np.isfinite(r.strike_cfm) else float(r.strike_cfm)),
        "nearest_fault_local_strike_deg": (None if not np.isfinite(r.strike_local)
                                           else float(r.strike_local)),
        "on_fault": bool(r.on_fault),
        "d_saf_creep_km": float(r.d_saf_creep_km),
        "d_geothermal_deg": float(r.d_geothermal_deg),
        "d_geothermal_km": float(r.d_geothermal_km),
        "reason": r.reason,
    } for r in u.itertuples()]
    print(f"  UNEXPLAINED SILENT: {len(rows)} cells "
          f"(silent {int(silent.sum())} - SAF-creep {int((silent & near_saf).sum())} "
          f"- geothermal {int((silent & ~near_saf & near_geo).sum())})")
    for r in rows[:15]:
        print(f"    {r['lat']:.1f} {r['lon']:.1f}  load={r['Mdot_geo_Nm_per_yr']:.3g}  "
              f"n_tr={r['n_train']:4d} "
              f"chi={r['chi_train']:.2e}  {r['style']:14s} d_f={r['d_fault_km']:5.1f} km "
              f"strike={r['nearest_fault_strike_deg']}  {'ON' if r['on_fault'] else 'OFF'}")

    # same list under the flat-11 ledger, for the H-sensitivity statement
    u_flat = set(flat.index[(flat.cls == "SILENT-LOADING")
                            & ~near_saf.to_numpy() & ~near_geo.to_numpy()])
    u_var = set(u.index)
    jac_u = (len(u_flat & u_var) / len(u_flat | u_var)) if (u_flat | u_var) else None
    h_sens["unexplained_silent_flat_n"] = len(u_flat)
    h_sens["unexplained_silent_variable_n"] = len(u_var)
    h_sens["jaccard_unexplained_flat_vs_variable"] = (None if jac_u is None else float(jac_u))

    # ---------------- output
    res = {
        "experiment": "EXP-K style-stratified fault-resolved ledger",
        "protocol": "PATTERN_PROTOCOL.md section EXP-K (frozen 2026-08-09)",
        "predecessor": "EXP-J (exp_j_stress_ledger.py / results_exp_j.json), style-blind",
        "cfm_geometry_used": cfm_used,
        "params": {
            "mu_Pa": MU, "Mc": MC, "cell_deg": CELL, "cell_alignment": "-122.0+0.2k / 31.5+0.2k",
            "train_end_utc": str(TRAIN_END), "train_years": TRAIN_YEARS,
            "test_years": test_years, "test_end": str(t_end),
            "min_events_chi": MIN_EVENTS_CHI, "chi_silent_lt": CHI_SILENT,
            "chi_coupled_hi": CHI_COUPLED_HI, "rate_regularizer_ev_per_yr": RATE_REGULARIZER,
            "style_threshold": STYLE_THRESH, "on_fault_km": ON_FAULT_KM,
            "saf_creeping_segment": SAF_CREEP, "saf_exclude_km": SAF_EXCLUDE_KM,
            "geothermal_exclude_deg": GEOTHERMAL_DEG,
            "chi_definition": "chi = Mdot_seis / Mdot_geo with Mdot_geo = 2*mu*H*A*max_shear "
                              "(EXP-J definition kept so chi is comparable); the fault-RESOLVED "
                              "shear is reported per cell and used for stratification/reporting, "
                              "not substituted into chi",
        },
        "k1_tensor_and_style": {
            "midas_stations_after_qc": n_sta,
            "n_cells": n_cells,
            "shear_crosscheck_vs_exp_j_grid": shear_check,
            "style_counts": style_counts,
            "styleness_quantiles": {q: float(np.nanpercentile(led.styleness, q))
                                    for q in (5, 25, 50, 75, 95)},
        },
        "k2_fault_geometry": {
            "cfm_geometry_used": cfm_used,
            "geometry_file": CFM_GEOM.name,
            "geometry_file_format": "GMT multi-segment; segment headers are '> strike dip rake' "
                                    "(3 numeric fields, all 557 headers parsed), followed by "
                                    "'lon lat elevation_m' point lines",
            "traces_file": CFM_TRACES.name,
            "traces_file_format": "same GMT multi-segment traces with '> \"fault_name\"' headers; "
                                  "point lines identical to the geometry file",
            "n_segments": geom["n_segments"], "n_trace_points": geom["n_points"],
            "n_named_segments": geom["traces_n_segments"],
            "traces_match_geometry_points": geom["traces_match_geometry"],
            "bad_headers": geom["bad_headers"],
            "strike_used": ("CFM segment-header strike (mod 180)" if cfm_used
                            else "FALLBACK: local strike from consecutive trace points"),
            "cfm_vs_local_strike_agreement": strike_agreement,
            "rotation": "theta = 90deg - strike_azimuth (grid x=east, y=north); "
                        "|e_ns| = |0.5*(eyy-exx)*sin(2*theta) + exy*cos(2*theta)|",
            "n_on_fault": int(led.on_fault.sum()), "n_off_fault": int((~led.on_fault).sum()),
            "median_d_fault_km": float(np.median(led.d_fault_km)),
            "median_resolved_over_max_shear": float(np.median(
                led.resolved_shear / led.max_shear.replace(0, np.nan))),
        },
        "k3_H_sensitivity": h_sens,
        "k4_strata": {
            "definition": "style (3) x on/off-fault (2); primary ledger = variable H",
            "n_cells_total_check": strata_total,
            "strata": strata,
            "kruskal_wallis_log_chi_across_styles": kw,
        },
        "k5_unexplained_silent": {
            "definition": "variable-H SILENT-LOADING cells, minus cells within "
                          f"{SAF_EXCLUDE_KM:.0f} km of the SAF creeping segment, minus cells "
                          f"within {GEOTHERMAL_DEG} deg of any tsi_map.GEOTHERMAL field",
            "n_silent_variable_H": int(silent.sum()),
            "n_excluded_saf_creep": int((silent & near_saf).sum()),
            "n_excluded_geothermal": int((silent & ~near_saf & near_geo).sum()),
            "n_unexplained": len(rows),
            "unexplained_reason_counts": {k: int(v) for k, v in
                                          var[unexp].reason.value_counts().items()},
            "unexplained_style_counts": {k: int(v) for k, v in
                                         var[unexp]["style"].value_counts().items()},
            "n_unexplained_on_fault": int(var[unexp].on_fault.sum()),
            "sorted_by": "Mdot_geo descending",
            "cells": rows,
        },
    }
    (HERE / "results_exp_k.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    draw_maps(var, lat_edges, lon_edges, list(u.index))
    print("-> results_exp_k.json, maps/exp_k_stratified.png")


if __name__ == "__main__":
    main()
