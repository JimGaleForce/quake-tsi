"""EQ-1/EQ-18 step C — three-region extension (Apennines, Taupo, Iceland) + Parkfield control.

Mainshock-level perigee-syzygy timing test (the artifact-free version of the question) plus the
exploratory class-level enrichment with the circular-shift null. Fault-style assignment is NOT
attempted here (no CFM equivalent loaded); the Apennines belt is treated as predominantly
normal-faulting, Taupo as extensional — labeled exploratory, distinct from the frozen EQ-1
confirmatory run. Also lists the famous mainshocks and whether each fell in-window.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import binomtest

HERE = Path(__file__).parent
RAW = HERE / "data" / "raw"
DERIVED = HERE / "data"
P = json.loads((HERE / "protocol_params.json").read_text(encoding="utf-8"))

W = {"new_below": 30, "full_above": 150, "ndist_below": 0.25,
     "perigee_km": 356500, "apogee_km": 406700}
SYNODIC_S = 29.53059 * 86400
N_SURR = 10000
SEED = 20260723

FAMOUS = {
    "apennines": [("L'Aquila 2009", "2009-04-06"), ("Amatrice 2016", "2016-08-24"),
                  ("Visso 2016", "2016-10-26"), ("Norcia 2016", "2016-10-30")],
    "taupo": [("Kawerau 1987 (Edgecumbe)", "1987-03-02")],
    "iceland": [("South Iceland 2000a", "2000-06-17"), ("South Iceland 2000b", "2000-06-21"),
                ("South Iceland 2008", "2008-05-29")],
    "parkfield_control": [("Parkfield 2004", "2004-09-28")],
}


def build_grid():
    cache = DERIVED / "lunar_grid_1980_2027.npz"
    if cache.exists():
        z = np.load(cache)
        return z["t_unix"], z["elong"], z["ndist"]
    from astropy.time import Time
    from astropy.coordinates import get_body
    import astropy.units as u
    t = pd.date_range("1980-01-01", "2027-01-01", freq="1h", tz="UTC")
    at = Time(t.to_pydatetime())
    moon = get_body("moon", at)
    sun = get_body("sun", at)
    elong = moon.separation(sun).deg
    r = moon.distance.to(u.km).value
    nd = np.clip((r - W["perigee_km"]) / (W["apogee_km"] - W["perigee_km"]), 0, 1)
    tu = (t - pd.Timestamp(0, tz="UTC")).total_seconds().to_numpy()
    np.savez_compressed(cache, t_unix=tu, elong=elong, ndist=nd)
    return tu, elong, nd


def mask_at(te, grid, shift=0.0):
    tg, e, d = grid
    ee = np.interp(te + shift, tg, e)
    dd = np.interp(te + shift, tg, d)
    return ((ee < W["new_below"]) | (ee > W["full_above"])) & (dd < W["ndist_below"])


def load_region(key):
    frames = []
    d = RAW / key
    for f in sorted(d.iterdir()):
        if f.stat().st_size == 0:
            continue
        if f.suffix == ".csv":  # USGS
            try:
                df = pd.read_csv(f, usecols=["time", "latitude", "longitude", "depth", "mag"])
            except Exception:
                continue
            df = df.rename(columns={"latitude": "lat", "longitude": "lon"})
            df["t"] = pd.to_datetime(df.time, utc=True, format="mixed")
        else:  # FDSN text (INGV / GeoNet), pipe-delimited with # header
            try:
                df = pd.read_csv(f, sep="|", comment=None)
            except Exception:
                continue
            df.columns = [c.strip().lstrip("#").strip().lower() for c in df.columns]
            need = {"time", "latitude", "longitude", "magnitude"}
            if not need.issubset(df.columns):
                continue
            df = df.rename(columns={"latitude": "lat", "longitude": "lon", "magnitude": "mag",
                                    "depth/km": "depth"})
            df["t"] = pd.to_datetime(df.time, utc=True, format="mixed")
        frames.append(df[["t", "lat", "lon", "mag"]])
    out = pd.concat(frames, ignore_index=True).dropna(subset=["t", "lat", "lon", "mag"])
    out["t_unix"] = (out.t - pd.Timestamp(0, tz="UTC")).dt.total_seconds()
    return out.sort_values("t_unix").reset_index(drop=True)


def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1, la2, lo2 = map(np.radians, [lat1, lon1, lat2, lon2])
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def main():
    grid = build_grid()
    _, e_g, d_g = grid
    fwin = float((((e_g < 30) | (e_g > 150)) & (d_g < 0.25)).mean())
    print(f"window fraction of time: {fwin:.4f}")
    rng = np.random.default_rng(SEED)
    results = {"window_fraction": fwin}

    for key in ["apennines", "taupo", "iceland", "parkfield_control"]:
        df = load_region(key)
        ms = df[df.mag >= 5.0]
        te = ms.t_unix.to_numpy()
        k = int(mask_at(te, grid).sum())
        n = len(ms)
        p = binomtest(k, n, fwin, alternative="greater").pvalue if n else np.nan
        print(f"\n=== {key}: events {len(df)}, mainshocks M>=5 {n}, in-window {k} "
              f"({(k/n if n else 0):.1%} vs {fwin:.1%}), binom p={p:.3f}")

        famous = []
        for name, day in FAMOUS.get(key, []):
            t0, t1 = pd.Timestamp(day, tz="UTC"), pd.Timestamp(day, tz="UTC") + pd.Timedelta(days=1)
            sel = ms[(ms.t >= t0) & (ms.t < t1)]
            if len(sel):
                big = sel.loc[sel.mag.idxmax()]
                inw = bool(mask_at(np.array([big.t_unix]), grid)[0])
                famous.append({"name": name, "mag": float(big.mag), "in_window": inw})
                print(f"  {name}: M{big.mag:.1f} in-window={inw}")

        # exploratory class-level enrichment (aftershocks M<5, <24h, <=100km of an M>=5 mainshock)
        cand = df[df.mag < 5.0]
        best = np.full(len(cand), False)
        for _, m in ms.iterrows():
            dt = cand.t_unix.to_numpy() - m.t_unix
            sel = (dt > 0) & (dt <= 86400)
            if not sel.any():
                continue
            idx = np.where(sel)[0]
            d = haversine_km(m.lat, m.lon, cand.lat.to_numpy()[idx], cand.lon.to_numpy()[idx])
            best[idx[d <= 100]] = True
        af = cand[best]
        te_a = af.t_unix.to_numpy()
        obs = int(mask_at(te_a, grid).sum())
        shifts = rng.uniform(0, SYNODIC_S, size=N_SURR)
        surr = np.array([mask_at(te_a, grid, s).sum() for s in shifts])
        er = obs / surr.mean() if surr.mean() > 0 else np.nan
        pa = (np.sum(surr >= obs) + 1) / (N_SURR + 1)
        print(f"  aftershock class: n={len(af)} obs={obs} exp={surr.mean():.1f} "
              f"ER={er:.3f} p={pa:.4f}")
        results[key] = {"n_events": int(len(df)), "n_mainshocks": n, "ms_in_window": k,
                        "ms_binom_p": float(p), "famous": famous,
                        "aftershock_class": {"n": int(len(af)), "obs": obs,
                                             "exp": float(surr.mean()),
                                             "ER": float(er), "p": float(pa)}}

    (HERE / "results_three_region.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print("\n-> results_three_region.json")


if __name__ == "__main__":
    main()
