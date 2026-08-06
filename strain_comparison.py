"""EQ-20 follow-up - compare our SoCal strain-rate field against the Kreemer & Young (2022)
velocity compilation (Burgmann's 2026-07-29 suggestion).

Kreemer & Young's gridded MELD/Haines-Holt model is not publicly deposited (their Dataverse
release doi:10.7910/DVN/BICMWB contains the GPS velocity table only), so the quantitative
comparison runs our exact estimator (Gaussian-weighted plane fit, sigma=40 km, 0.1 deg grid -
tsi_map.py) on BOTH velocity sets: NGL MIDAS IGS14 (our original) and the K&Y compilation
(3,049 stations: NGL + campaign + other networks, already in mm/yr). Identical method on both
isolates velocity-data sensitivity; spatial-pattern agreement with their published Figure 2 and
the eps_min >= 47 nanostrain/yr SAF-system contour (their Fig. 8 criterion) checks the
realization against the paper. Invariants follow their definitions: second invariant
sqrt(e1^2+e2^2); dilatation e1+e2; eps_min = max(|e1|,|e2|).
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
import importlib.util

HERE = Path(__file__).parent
spec = importlib.util.spec_from_file_location("tm", HERE / "tsi_map.py")
tm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tm)

BOX, GRID_STEP, SIGMA_KM, MIN_STA = tm.BOX, tm.GRID_STEP, tm.SIGMA_KM, tm.MIN_STA
MAX_SU = tm.MAX_SU
KY_FILE = HERE / "data" / "kreemer_young" / "GPS_table.txt"
HIGH_STRAIN = 47.0  # nanostrain/yr, Kreemer & Young Fig. 8 contour on eps_min


def load_ky():
    df = pd.read_csv(KY_FILE, sep=r"\s+")
    df = df.rename(columns={"sdve": "se", "sdvn": "sn"})
    for c in ["lat", "lon", "ve", "vn", "se", "sn"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["lat", "lon", "ve", "vn", "se", "sn"]).reset_index(drop=True)
    assert df.lat.abs().max() <= 90 and df.lon.abs().max() <= 180
    return df  # ve/vn/se/sn in mm/yr


def strain_grid_full(df, ve_mm, vn_mm, se_mm, sn_mm):
    """tsi_map.strain_grid with the full gradient tensor kept (units mm/yr in, nanostrain/yr out)."""
    d = df[(df.lat.between(*BOX["lat"])) & (df.lon.between(*BOX["lon"]))
           & (se_mm < MAX_SU) & (sn_mm < MAX_SU)]
    ve = ve_mm[d.index].to_numpy()
    vn = vn_mm[d.index].to_numpy()
    lat0 = np.mean(BOX["lat"])
    kx = 111.32 * np.cos(np.radians(lat0))
    ky = 110.57
    x = d.lon.to_numpy() * kx
    y = d.lat.to_numpy() * ky
    lats = np.arange(BOX["lat"][0], BOX["lat"][1] + 1e-9, GRID_STEP)
    lons = np.arange(BOX["lon"][0], BOX["lon"][1] + 1e-9, GRID_STEP)
    exx = np.full((len(lats), len(lons)), np.nan)
    eyy = np.full_like(exx, np.nan)
    exy = np.full_like(exx, np.nan)
    for i, la in enumerate(lats):
        for j, lo in enumerate(lons):
            gx, gy = lo * kx, la * ky
            d2 = (x - gx) ** 2 + (y - gy) ** 2
            w = np.exp(-d2 / (2 * SIGMA_KM ** 2))
            sel = w > 0.01
            if w[sel].sum() < MIN_STA * 0.2 or sel.sum() < MIN_STA:
                continue
            X = np.c_[np.ones(sel.sum()), x[sel] - gx, y[sel] - gy]
            A = X * w[sel][:, None]
            try:
                ce, *_ = np.linalg.lstsq(A, ve[sel] * w[sel], rcond=None)
                cn, *_ = np.linalg.lstsq(A, vn[sel] * w[sel], rcond=None)
            except np.linalg.LinAlgError:
                continue
            # mm/yr per km = 1e-6/yr; ×1000 -> nanostrain/yr
            exx[i, j] = ce[1] * 1000
            eyy[i, j] = cn[2] * 1000
            exy[i, j] = 0.5 * (ce[2] + cn[1]) * 1000
    return lats, lons, exx, eyy, exy, len(d)


def invariants(exx, eyy, exy):
    mean = (exx + eyy) / 2
    r = np.sqrt(((exx - eyy) / 2) ** 2 + exy ** 2)
    e1, e2 = mean + r, mean - r
    return {
        "dil": exx + eyy,
        "max_shear": r,
        "inv2": np.sqrt(e1 ** 2 + e2 ** 2),
        "eps_min": np.maximum(np.abs(e1), np.abs(e2)),  # K&Y's eps_min (their naming)
    }


def compare(a, b):
    ok = ~np.isnan(a) & ~np.isnan(b)
    from scipy.stats import pearsonr, spearmanr
    return {
        "n_nodes": int(ok.sum()),
        "pearson_r": float(pearsonr(a[ok], b[ok])[0]),
        "spearman_rho": float(spearmanr(a[ok], b[ok])[0]),
        "median_abs_diff": float(np.nanmedian(np.abs(a[ok] - b[ok]))),
        "median_a": float(np.nanmedian(a[ok])),
        "median_b": float(np.nanmedian(b[ok])),
    }


def main():
    midas = tm.load_midas()
    lats, lons, *mid_t, n_mid = strain_grid_full(
        midas, midas.ve * 1000, midas.vn * 1000, midas.se * 1000, midas.sn * 1000)
    ky = load_ky()
    _, _, *ky_t, n_ky = strain_grid_full(ky, ky.ve, ky.vn, ky.se, ky.sn)
    print(f"stations after QC: MIDAS {n_mid}, K&Y {n_ky}")
    M, K = invariants(*mid_t), invariants(*ky_t)

    out = {"stations": {"midas": n_mid, "kreemer_young": n_ky},
           "fields": {k: compare(M[k], K[k]) for k in M}}

    # K&Y Fig. 8-style high-strain mask agreement (eps_min >= 47 nanostrain/yr)
    a = M["eps_min"] >= HIGH_STRAIN
    b = K["eps_min"] >= HIGH_STRAIN
    ok = ~np.isnan(M["eps_min"]) & ~np.isnan(K["eps_min"])
    inter, union = (a & b & ok).sum(), ((a | b) & ok).sum()
    out["high_strain_mask_47"] = {"midas_frac": float(a[ok].mean()), "ky_frac": float(b[ok].mean()),
                                  "jaccard": float(inter / union) if union else None}

    # site benchmarks
    sites = {}
    for name, (gla, glo) in tm.GEOTHERMAL.items():
        i = int(round((gla - lats[0]) / GRID_STEP))
        j = int(round((glo - lons[0]) / GRID_STEP))
        if 0 <= i < len(lats) and 0 <= j < len(lons):
            sites[name] = {src: {k: (None if np.isnan(F[k][i, j]) else round(float(F[k][i, j]), 2))
                                 for k in ("dil", "max_shear")}
                           for src, F in [("midas", M), ("ky", K)]}
    out["sites"] = sites

    (HERE / "results_strain_comparison.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out["fields"], indent=2))
    print(json.dumps(out["high_strain_mask_47"], indent=2))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = [("Max shear (nanostrain/yr)", "max_shear", "viridis", (0, 300)),
            ("Dilatation (nanostrain/yr)", "dil", "RdBu_r", (-100, 100))]
    fig, axes = plt.subplots(2, 3, figsize=(21, 12), constrained_layout=True)
    for r, (title, key, cmap, clim) in enumerate(rows):
        diff = M[key] - K[key]
        panels = [(f"NGL MIDAS - {title}", M[key], cmap, clim),
                  (f"Kreemer & Young velocities - {title}", K[key], cmap, clim),
                  ("MIDAS minus K&Y", diff, "RdBu_r", (-50, 50))]
        for c, (t, z, cm, cl) in enumerate(panels):
            ax = axes[r, c]
            im = ax.pcolormesh(lons, lats, z, cmap=cm, vmin=cl[0], vmax=cl[1], shading="auto")
            for name, (gla, glo) in tm.GEOTHERMAL.items():
                if BOX["lat"][0] <= gla <= BOX["lat"][1] and BOX["lon"][0] <= glo <= BOX["lon"][1]:
                    ax.plot(glo, gla, "r^", ms=8, mec="k")
            ax.set_title(t, fontsize=11)
            fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle("Same estimator (Gaussian sigma=40 km plane fit, 0.1 deg), two velocity sets - "
                 "our NGL MIDAS field vs the Kreemer & Young (2022) compilation", fontsize=13)
    fig.savefig(tm.OUT / "strain_comparison_kreemer_young.png", dpi=160)
    print("-> maps/strain_comparison_kreemer_young.png; results_strain_comparison.json")


if __name__ == "__main__":
    main()
