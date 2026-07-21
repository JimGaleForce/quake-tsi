"""EQ-22 — Long Valley TSI prediction test, per LONG_VALLEY_PROTOCOL.md (frozen).

This script implements: catalog load, Mc estimate, GK-window declustering to a background set,
the astropy tidal-potential volumetric strain series at the caldera, and the phase-modulation
test (bias-corrected, 3,000 synthetics). The volumetric test is the protocol's SECONDARY endpoint
(always reported); the PRIMARY endpoint reuses whichever statistic validates on Coso and is run
separately once results_coso_fm.json exists.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.signal import find_peaks

HERE = Path(__file__).parent
RAW = HERE / "data" / "raw" / "long_valley"
DERIVED = HERE / "data"
LV = dict(lat=37.68, lon=-118.90)
ERA = ("1985-01-01", "2026-06-30")
DT = 600.0
N_BINS = 36
N_SYNTH = 3000
SEED = 20260724

GK_T = lambda m: 10 ** (0.032 * m + 2.7389) * 86400.0   # s
GK_L = lambda m: 10 ** (0.1238 * m + 0.983)             # km


def load_catalog():
    frames = []
    for f in sorted(RAW.glob("*.txt")):
        if f.stat().st_size == 0:
            continue
        df = pd.read_csv(f, sep="|")
        df.columns = [c.strip().lstrip("#").strip().lower() for c in df.columns]
        df = df.rename(columns={"latitude": "lat", "longitude": "lon", "magnitude": "mag",
                                "depth/km": "depth"})
        df["t"] = pd.to_datetime(df.time, utc=True, format="mixed")
        frames.append(df[["t", "lat", "lon", "depth", "mag"]])
    out = pd.concat(frames, ignore_index=True).dropna()
    out = out[(out.t >= pd.Timestamp(ERA[0], tz="UTC")) & (out.t <= pd.Timestamp(ERA[1], tz="UTC"))]
    out["t_unix"] = (out.t - pd.Timestamp(0, tz="UTC")).dt.total_seconds()
    return out.sort_values("t_unix").reset_index(drop=True)


def estimate_mc(mags):
    """Max-curvature + 0.2 (standard correction)."""
    bins = np.arange(0.95, mags.max() + 0.1, 0.1)
    h, _ = np.histogram(mags, bins=bins)
    return round(float(bins[np.argmax(h)] + 0.05 + 0.2), 1)


def nnd_background(df, b=1.0, df_frac=1.6, log_eta0=-5.0, lookback=4000):
    """Zaliapin & Ben-Zion (2013) nearest-neighbor declustering (protocol Amendment 1).
    eta_ij = dt_years * r_km^df * 10^(-b*m_i) over preceding events; background = eta > eta0."""
    la, lo, tm, mg = (df[c].to_numpy() for c in ["lat", "lon", "t_unix", "mag"])
    kx = 111.32 * np.cos(np.radians(LV["lat"]))
    x, y = lo * kx, la * 110.57
    n = len(df)
    eta = np.full(n, np.inf)
    yr = 365.25 * 86400.0
    for j in range(1, n):
        i0 = max(0, j - lookback)
        dt = (tm[j] - tm[i0:j]) / yr
        ok = dt > 0
        if not ok.any():
            continue
        r = np.hypot(x[i0:j] - x[j], y[i0:j] - y[j])
        r = np.maximum(r, 0.1)
        e = dt * (r ** df_frac) * (10.0 ** (-b * mg[i0:j]))
        e[~ok] = np.inf
        eta[j] = e.min()
    return df[(eta > 10.0 ** log_eta0)].reset_index(drop=True)


def gk_background(df):
    """Remove every event inside the GK space-time window of any larger event."""
    la, lo, tm, mg = (df[c].to_numpy() for c in ["lat", "lon", "t_unix", "mag"])
    n = len(df)
    removed = np.zeros(n, dtype=bool)
    order = np.argsort(-mg)
    kx = 111.32 * np.cos(np.radians(LV["lat"]))
    x, y = lo * kx, la * 110.57
    for i in order:
        if removed[i]:
            continue
        Tw, Lw = GK_T(mg[i]), GK_L(mg[i])
        cand = np.flatnonzero((np.abs(tm - tm[i]) <= Tw) & ~removed & (mg < mg[i]))
        if len(cand) == 0:
            continue
        d = np.hypot(x[cand] - x[i], y[cand] - y[i])
        removed[cand[d <= Lw]] = True
    return df[~removed].reset_index(drop=True)


def tidal_vol_series():
    """Volumetric tidal strain proxy at the caldera: sum GM/r^3 (3cos^2 z - 1), 600-s sampling."""
    cache = DERIVED / "lv_tidal_vol.npz"
    if cache.exists():
        z = np.load(cache)
        return z["t_unix"], z["eps"]
    from astropy.time import Time
    from astropy.coordinates import get_body
    import astropy.units as u
    GM = {"moon": 4.9028e12, "sun": 1.32712e20}
    lat_r = np.radians(LV["lat"])
    ts_all, eps_all = [], []
    years = pd.date_range(ERA[0], ERA[1], freq="YS", tz="UTC")
    edges = list(years) + [pd.Timestamp(ERA[1], tz="UTC")]
    for a, b in zip(edges[:-1], edges[1:]):
        tt = pd.date_range(a, b, freq=f"{int(DT)}s", inclusive="left", tz="UTC")
        if len(tt) == 0:
            continue
        at = Time(tt.to_pydatetime())
        eps = np.zeros(len(tt))
        for body, gm in GM.items():
            bo = get_body(body, at)
            r = bo.distance.to(u.m).value
            lst = at.sidereal_time("apparent", longitude=LV["lon"] * u.deg).rad
            ha = lst - bo.ra.rad
            cosz = (np.sin(lat_r) * np.sin(bo.dec.rad)
                    + np.cos(lat_r) * np.cos(bo.dec.rad) * np.cos(ha))
            eps += gm / r ** 3 * (3 * cosz ** 2 - 1)
        ts_all.append((tt - pd.Timestamp(0, tz="UTC")).total_seconds().to_numpy())
        eps_all.append(eps)
        print(f"  ephemeris {a.year} done")
    t_unix = np.concatenate(ts_all)
    eps = np.concatenate(eps_all).astype(np.float32)
    np.savez_compressed(cache, t_unix=t_unix, eps=eps)
    return t_unix, eps


def phase_series(v):
    pk, _ = find_peaks(v)
    tr, _ = find_peaks(-v)
    phase = np.full(len(v), np.nan, dtype=np.float32)
    for i in range(len(tr) - 1):
        a, b = tr[i], tr[i + 1]
        mids = pk[(pk > a) & (pk < b)]
        if len(mids) == 0:
            continue
        m = mids[np.argmax(v[mids])]
        phase[a:m] = -180.0 + 180.0 * (np.arange(a, m) - a) / (m - a)
        phase[m:b] = 180.0 * (np.arange(m, b) - m) / (b - m)
    return phase


def fit_pm(hist):
    d = hist / hist.sum() * N_BINS
    th = np.radians(np.linspace(-180, 180, N_BINS + 1)[:-1] + 180 / N_BINS)
    A = np.c_[np.ones_like(th), np.cos(th), np.sin(th)]
    coef, *_ = np.linalg.lstsq(A, d, rcond=None)
    return np.hypot(coef[1], coef[2]) / coef[0], np.degrees(np.arctan2(coef[2], coef[1]))


def modulation_test(ev_t, grid_t, phase, rng):
    idx = np.clip(np.searchsorted(grid_t, ev_t), 0, len(phase) - 1)
    ph = phase[idx]
    ph = ph[~np.isnan(ph)]
    bins = np.linspace(-180, 180, N_BINS + 1)
    valid = np.flatnonzero(~np.isnan(phase))
    p0h, _ = np.histogram(phase[rng.choice(valid, 300000)], bins=bins)
    p0h = p0h / p0h.sum() * N_BINS
    h, _ = np.histogram(ph, bins=bins)
    hn = (h / h.sum() * N_BINS) / p0h
    obs, phi = fit_pm(hn)
    syn = np.empty(N_SYNTH)
    for k in range(N_SYNTH):
        hh, _ = np.histogram(phase[rng.choice(valid, len(ph))], bins=bins)
        syn[k], _ = fit_pm((hh / hh.sum() * N_BINS) / p0h)
    return {"n": int(len(ph)), "Pm_over_P0": float(obs), "phi_deg": float(phi),
            "synth_p97.5": float(np.percentile(syn, 97.5)),
            "significant_95": bool(obs > np.percentile(syn, 97.5)),
            "empirical_p": float((np.sum(syn >= obs) + 1) / (N_SYNTH + 1))}


def main():
    rng = np.random.default_rng(SEED)
    df = load_catalog()
    print(f"Long Valley catalog: {len(df)} events {df.t.min().date()}..{df.t.max().date()}")
    mc = estimate_mc(df.mag.to_numpy())
    floor = max(1.0, mc)
    df = df[df.mag >= floor].reset_index(drop=True)
    print(f"Mc estimate {mc}; floor {floor}; events above floor: {len(df)}")
    bg = nnd_background(df)
    print(f"background after Zaliapin-Ben-Zion NND declustering (Amendment 1): {len(bg)}")
    gt, eps = tidal_vol_series()
    phase = phase_series(eps)
    res = {"mc": mc, "n_above_floor": int(len(df)), "n_background": int(len(bg))}
    res["background_volumetric"] = modulation_test(bg.t_unix.to_numpy(), gt, phase, rng)
    r = res["background_volumetric"]
    print(f"SECONDARY (volumetric phase, background set): n={r['n']} Pm/P0={r['Pm_over_P0']:.3f} "
          f"thr={r['synth_p97.5']:.3f} p={r['empirical_p']:.4f} sig={r['significant_95']}")
    res["all_events_volumetric"] = modulation_test(df.t_unix.to_numpy(), gt, phase, rng)
    r = res["all_events_volumetric"]
    print(f"exploratory (all events): n={r['n']} Pm/P0={r['Pm_over_P0']:.3f} "
          f"thr={r['synth_p97.5']:.3f} p={r['empirical_p']:.4f} sig={r['significant_95']}")
    (HERE / "results_long_valley.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    print("-> results_long_valley.json")


if __name__ == "__main__":
    main()
