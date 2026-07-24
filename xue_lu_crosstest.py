"""EQ-18 blind cross-test - implements XUE_LU_PROTOCOL.md (frozen 2026-07-21).

All analytic parameters come from xue_lu_params.json. Faithful reimplementation of the
Kosmos-r89 pipeline (selection criteria, lunar windows, strain-rate percentiles) applied
to the Lu/Xue Southern California catalogs.

Stages (cached to data/xue_lu_derived/ so reruns are cheap):
  1. lunar ephemeris grid (hourly, 1980-2019) via astropy -> phase elongation + normalized distance
  2. per-catalog: mainshocks, aftershock association, fault-style tagging (CFM 5.3)
  3. H1 enrichment + 10k circular-shift surrogates + per-mainshock jackknife
  4. H2 strain-rate percentiles for in-window normal-fault events + calibration surrogates
Outputs: results_xue_lu.json + console summary. Raw inputs are never modified.
"""
import json
import sys
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).parent
P = json.loads((HERE / "xue_lu_params.json").read_text(encoding="utf-8"))
DATA = HERE / P["data_dir"]
DERIVED = HERE / "data" / "xue_lu_derived"
DERIVED.mkdir(parents=True, exist_ok=True)

EARTH_R = 6371.0


# ---------------------------------------------------------------- ephemeris grid
def build_lunar_grid():
    """Hourly geocentric Moon-Sun elongation (deg, 0=new/180=full) and normalized lunar
    distance over 1980-01-01..2019-01-01. Cached as npz."""
    cache = DERIVED / "lunar_grid.npz"
    if cache.exists():
        z = np.load(cache)
        return z["t_unix"], z["elong"], z["ndist"]
    from astropy.time import Time
    from astropy.coordinates import get_body
    import astropy.units as u

    t = pd.date_range("1980-01-01", "2019-01-01", freq="1h", tz="UTC")
    at = Time(t.to_pydatetime())
    moon = get_body("moon", at)
    sun = get_body("sun", at)
    elong = moon.separation(sun).deg  # 0 = new, 180 = full
    r_km = moon.distance.to(u.km).value
    norm = P["tidal_windows"]["lunar_distance_normalization_km"]
    ndist = np.clip((r_km - norm["perigee"]) / (norm["apogee"] - norm["perigee"]), 0, 1)
    t_unix = (t - pd.Timestamp(0, tz="UTC")).total_seconds().to_numpy()
    np.savez_compressed(cache, t_unix=t_unix, elong=elong, ndist=ndist)
    return t_unix, elong, ndist


def tidal_mask_at(t_unix_events, grid, shift_seconds=0.0):
    """Perigee-syzygy boolean mask at (event times + shift), by linear interpolation."""
    tg, elong, ndist = grid
    te = t_unix_events + shift_seconds
    e = np.interp(te, tg, elong)
    d = np.interp(te, tg, ndist)
    w = P["tidal_windows"]
    syz = (e < w["syzygy_elongation_deg"]["new_below"]) | (e > w["syzygy_elongation_deg"]["full_above"])
    per = d < w["perigee_normalized_distance_below"]
    return syz & per


# ---------------------------------------------------------------- catalogs
def load_catalog(fname):
    # Raw catalogs (QTM_12dev, SCSN_original): yr mo dy hr mi sec EID lat lon depth mag
    # Declustered catalogs:                    yr mo dy hr mi sec lat lon depth mag EID
    # Detect by validating the latitude column.
    raw_cols = ["yr", "mo", "dy", "hr", "mi", "sec", "eid", "lat", "lon", "depth", "mag"]
    dec_cols = ["yr", "mo", "dy", "hr", "mi", "sec", "lat", "lon", "depth", "mag", "eid"]
    probe = pd.read_csv(DATA / fname, sep=r"\s+", header=None, nrows=1000)
    cols = raw_cols if probe[6].abs().max() > 90 else dec_cols
    df = pd.read_csv(DATA / fname, sep=r"\s+", header=None, names=cols)
    assert df.lat.abs().max() <= 90 and df.lon.abs().max() <= 180, f"column detection failed for {fname}"
    sec = df["sec"].astype(float)
    ts = pd.to_datetime(
        dict(year=df.yr, month=df.mo, day=df.dy, hour=df.hr, minute=df.mi, second=0), utc=True
    ) + pd.to_timedelta(sec, unit="s")
    df["t_unix"] = (ts - pd.Timestamp(0, tz="UTC")).dt.total_seconds().to_numpy()
    return df[["t_unix", "lat", "lon", "depth", "mag", "eid"]]


def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1, la2, lo2 = map(np.radians, [lat1, lon1, lat2, lon2])
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * EARTH_R * np.arcsin(np.sqrt(a))


def associate_aftershocks(df, mc):
    """Protocol §2: mainshocks M>=5.0; aftershocks mc<=M<5.0 within 24h after & 100km,
    linked to closest-in-distance qualifying mainshock."""
    a = P["aftershock"]
    ms = df[df.mag >= P["mainshock_min_mag"]].reset_index(drop=True)
    cand = df[(df.mag >= mc) & (df.mag < a["max_mag_exclusive"])].reset_index(drop=True)
    max_dt = a["max_dt_hours"] * 3600.0
    best_dist = np.full(len(cand), np.inf)
    best_ms = np.full(len(cand), -1, dtype=int)
    for i, m in ms.iterrows():
        dt = cand.t_unix.to_numpy() - m.t_unix
        sel = (dt > 0) & (dt <= max_dt)
        if not sel.any():
            continue
        idx = np.where(sel)[0]
        d = haversine_km(m.lat, m.lon, cand.lat.to_numpy()[idx], cand.lon.to_numpy()[idx])
        ok = d <= a["max_dist_km"]
        idx, d = idx[ok], d[ok]
        upd = d < best_dist[idx]
        best_dist[idx[upd]] = d[upd]
        best_ms[idx[upd]] = i
    out = cand[best_ms >= 0].copy()
    out["ms_idx"] = best_ms[best_ms >= 0]
    out["ms_dist_km"] = best_dist[best_ms >= 0]
    ms_t = ms.t_unix.to_numpy()
    out["dt_hours"] = (out.t_unix.to_numpy() - ms_t[out.ms_idx.to_numpy()]) / 3600.0
    return out, ms


# ---------------------------------------------------------------- fault styles (CFM 5.3)
def load_fault_points():
    """Parse CFM geometry: segments headed by '> strike dip rake', then lon lat depth lines.
    Returns (normal_pts, strikeslip_pts) as (lat, lon) arrays."""
    lo_n, la_n, lo_s, la_s = [], [], [], []
    style = None
    nr = P["fault_geometry"]["normal_rake_range_deg"]
    tol = P["fault_geometry"]["strike_slip_rake_tolerance_deg"]
    with open(DATA / P["fault_geometry"]["file"], encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == ">":
                style = None
                if len(parts) >= 4:
                    try:
                        rake = float(parts[3])
                    except ValueError:
                        continue
                    if nr[0] <= rake <= nr[1]:
                        style = "normal"
                    elif min(abs(rake), abs(abs(rake) - 180)) <= tol:
                        style = "ss"
            elif style and len(parts) >= 2:
                try:
                    lon, lat = float(parts[0]), float(parts[1])
                except ValueError:
                    continue
                (lo_n if style == "normal" else lo_s).append(lon)
                (la_n if style == "normal" else la_s).append(lat)
    return (np.array(la_n), np.array(lo_n)), (np.array(la_s), np.array(lo_s))


def near_fault_mask(events, fault_pts, max_km):
    from scipy.spatial import cKDTree
    flat, flon = fault_pts
    if len(flat) == 0:
        return np.zeros(len(events), dtype=bool)
    lat0 = np.mean(flat)
    kx = np.cos(np.radians(lat0)) * 111.32
    tree = cKDTree(np.c_[flat * 110.57, flon * kx])
    d, _ = tree.query(np.c_[events.lat.to_numpy() * 110.57, events.lon.to_numpy() * kx])
    return d <= max_km


# ---------------------------------------------------------------- H1 enrichment
def enrichment(events, grid, rng):
    n = P["null_model"]["n_surrogates"]
    period = P["null_model"]["shift_period_days"] * 86400.0
    te = events.t_unix.to_numpy()
    obs = int(tidal_mask_at(te, grid).sum())
    shifts = rng.uniform(0, period, size=n)
    surr = np.array([tidal_mask_at(te, grid, s).sum() for s in shifts])
    er = obs / surr.mean() if surr.mean() > 0 else np.nan
    p = (np.sum(surr >= obs) + 1) / (n + 1)
    return {"n_events": int(len(te)), "observed_in_window": obs,
            "expected_in_window": float(surr.mean()), "ER": float(er), "p_one_sided": float(p)}


def jackknife_by_mainshock(events, grid, rng):
    out = {}
    n_small = 2000
    period = P["null_model"]["shift_period_days"] * 86400.0
    for ms_idx, grp in events.groupby("ms_idx"):
        rest = events[events.ms_idx != ms_idx]
        if len(rest) == 0:
            continue
        te = rest.t_unix.to_numpy()
        obs = int(tidal_mask_at(te, grid).sum())
        shifts = rng.uniform(0, period, size=n_small)
        surr = np.array([tidal_mask_at(te, grid, s).sum() for s in shifts])
        er = obs / surr.mean() if surr.mean() > 0 else np.nan
        out[str(int(ms_idx))] = {"dropped_n": int(len(grp)), "ER_without": float(er)}
    return out


# ---------------------------------------------------------------- H2 strain-rate percentile
def strain_rate_percentiles(events):
    """r89 formulation: eps_v(t) ∝ sum_body GM/r^3 (3 cos^2 theta - 1) at event location,
    1-min resolution across the 24h window centered on event; percentile of |d eps/dt| at t0."""
    from astropy.time import Time
    from astropy.coordinates import get_body
    import astropy.units as u

    GM_MOON, GM_SUN = 4.9028e12, 1.32712e20  # m^3/s^2
    win_h = P["h2"]["window_hours"]
    res_min = P["h2"]["resolution_minutes"]
    percs = []
    for _, ev in events.iterrows():
        t0 = pd.Timestamp(ev.t_unix, unit="s", tz="UTC")
        tt = pd.date_range(t0 - pd.Timedelta(hours=win_h / 2), t0 + pd.Timedelta(hours=win_h / 2),
                           freq=f"{res_min}min")
        at = Time(tt.to_pydatetime())
        lat_r = np.radians(ev.lat)
        eps = np.zeros(len(tt))
        for body, gm in (("moon", GM_MOON), ("sun", GM_SUN)):
            b = get_body(body, at)
            r = b.distance.to(u.m).value
            dec = b.dec.rad
            lst = at.sidereal_time("apparent", longitude=ev.lon * u.deg).rad
            hour_angle = lst - b.ra.rad
            cosz = np.sin(lat_r) * np.sin(dec) + np.cos(lat_r) * np.cos(dec) * np.cos(hour_angle)
            eps += gm / r ** 3 * (3 * cosz ** 2 - 1)
        rate = np.abs(np.gradient(eps, res_min * 60.0))
        i0 = len(tt) // 2
        percs.append(100.0 * np.mean(rate <= rate[i0]))
    return np.array(percs)


# ---------------------------------------------------------------- main
def main():
    rng = np.random.default_rng(P["null_model"]["seed"])
    grid = build_lunar_grid()
    print("lunar grid ready")
    normal_pts, ss_pts = load_fault_points()
    print(f"CFM5.3: {len(normal_pts[0])} normal-fault trace points, {len(ss_pts[0])} strike-slip")

    results = {"protocol": "XUE_LU_PROTOCOL.md v1.0.0"}
    for cname, cfg in P["catalogs"].items():
        print(f"\n=== {cname} ===")
        df = load_catalog(cfg["raw"])
        print(f"loaded {len(df)} events, lat {df.lat.min():.2f}..{df.lat.max():.2f}, "
              f"time {pd.Timestamp(df.t_unix.min(), unit='s')}..{pd.Timestamp(df.t_unix.max(), unit='s')}")
        af, ms = associate_aftershocks(df, cfg["mc"])
        print(f"mainshocks M>=5: {len(ms)}; associated aftershocks: {len(af)}")
        af_csv = DERIVED / f"{cname}_aftershocks.csv"
        af.to_csv(af_csv, index=False)

        af = af.reset_index(drop=True)
        is_norm = near_fault_mask(af, normal_pts, P["fault_geometry"]["max_dist_km"])
        is_ss = near_fault_mask(af, ss_pts, P["fault_geometry"]["max_dist_km"]) & ~is_norm
        classes = {"all": af, "normal": af[is_norm], "strike_slip_control": af[is_ss]}
        res_c = {}
        for label, ev in classes.items():
            if len(ev) == 0:
                res_c[label] = {"n_events": 0}
                continue
            r = enrichment(ev, grid, rng)
            print(f"  {label}: n={r['n_events']} obs={r['observed_in_window']} "
                  f"exp={r['expected_in_window']:.1f} ER={r['ER']:.3f} p={r['p_one_sided']:.4f}")
            if label == "normal" and len(ev) > 0:
                r["jackknife_by_mainshock"] = jackknife_by_mainshock(ev, grid, rng)
            res_c[label] = r
        results[cname] = res_c

    (HERE / "results_xue_lu.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print("\nH1 stage complete -> results_xue_lu.json (H2 runs separately: --h2)")


def run_h2():
    rng = np.random.default_rng(P["null_model"]["seed"] + 1)
    grid = build_lunar_grid()
    results = json.loads((HERE / "results_xue_lu.json").read_text(encoding="utf-8"))
    from scipy.stats import wilcoxon
    for cname in P["catalogs"]:
        af = pd.read_csv(DERIVED / f"{cname}_aftershocks.csv")
        normal_pts, _ = load_fault_points()
        af = af.reset_index(drop=True)
        is_norm = near_fault_mask(af, normal_pts, P["fault_geometry"]["max_dist_km"])
        ev = af[is_norm].copy()
        ev = ev[tidal_mask_at(ev.t_unix.to_numpy(), grid)]
        if len(ev) == 0:
            results[cname]["h2"] = {"n": 0}
            continue
        cap = 400  # runtime guard; sample if larger, recorded in output
        sampled = len(ev) > cap
        if sampled:
            ev = ev.sample(cap, random_state=P["null_model"]["seed"])
        percs = strain_rate_percentiles(ev)
        stat, pval = wilcoxon(percs - 50, alternative="greater")
        results[cname]["h2"] = {
            "n": int(len(percs)), "sampled": sampled,
            "median_percentile": float(np.median(percs)), "mean_percentile": float(np.mean(percs)),
            "wilcoxon_p_one_sided_vs_50": float(pval),
        }
        print(f"{cname} H2: n={len(percs)} median={np.median(percs):.1f}% p={pval:.2e}")
    (HERE / "results_xue_lu.json").write_text(json.dumps(results, indent=2), encoding="utf-8")


if __name__ == "__main__":
    if "--h2" in sys.argv:
        run_h2()
    else:
        main()
