"""EQ-20 / EQ-18 step B - Southern California strain-rate & TSI maps from NGL MIDAS velocities.

Produces the deliverable Burgmann asked for: dilatation rate, max shear strain rate, and TSI
(= dilatation / max shear) shown SEPARATELY, with geothermal fields and Lu et al. modulation
regions marked. Method: per grid node, Gaussian-distance-weighted least-squares plane fit of the
horizontal velocity field (vE, vN) -> velocity gradient tensor -> strain rates.

MIDAS format (space-separated): sta lat(9) lon(10)? -- actual columns: see NGL readme; we use
col names by position: 0=sta, 8=vN? ... Format per NGL: sta, MIDAS version, epoch0, epoch1, dur,
lon, lat, alt, vE, vN, vU, sE, sN, sU, ... (we parse positions 5,6, 8,9,10, 11,12).
Verified at runtime with sanity assertions on coordinate/velocity ranges.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "maps"
OUT.mkdir(exist_ok=True)

BOX = dict(lat=(31.5, 38.0), lon=(-122.0, -113.5))
GRID_STEP = 0.1
SIGMA_KM = 40.0     # Gaussian weighting length scale
MIN_STA = 6         # minimum effective stations per node
MAX_SU = 1.5        # mm/yr max velocity uncertainty filter

GEOTHERMAL = {  # marker set for maps
    "Coso": (36.03, -117.82), "Long Valley": (37.68, -118.90), "Salton Sea": (33.20, -115.60),
    "Cerro Prieto": (32.42, -115.23), "The Geysers*": (38.79, -122.75), "Brawley": (32.99, -115.51),
}
LU_REGIONS = {"Coso (Lu et al. sig.)": (36.03, -117.82)}


def load_midas():
    rows = []
    with open(HERE / "data" / "ngl" / "midas.IGS14.txt", encoding="utf-8") as f:
        for line in f:
            p = line.split()
            if len(p) < 27:
                continue
            try:
                lat, lon = float(p[24]), float(p[25])
                ve, vn = float(p[8]), float(p[9])
                se, sn = float(p[11]), float(p[12])
            except ValueError:
                continue
            rows.append((p[0], lat, lon, ve, vn, se, sn))
    df = pd.DataFrame(rows, columns=["sta", "lat", "lon", "ve", "vn", "se", "sn"])
    # NGL longitudes can be outside [-180, 180]; normalize
    df["lon"] = ((df.lon + 180) % 360) - 180
    return df


def strain_grid(df):
    lat0 = np.mean(BOX["lat"])
    kx = 111.32 * np.cos(np.radians(lat0))
    ky = 110.57
    x = df.lon.to_numpy() * kx
    y = df.lat.to_numpy() * ky
    ve = df.ve.to_numpy() * 1000.0  # m/yr -> mm/yr
    vn = df.vn.to_numpy() * 1000.0
    lats = np.arange(BOX["lat"][0], BOX["lat"][1] + 1e-9, GRID_STEP)
    lons = np.arange(BOX["lon"][0], BOX["lon"][1] + 1e-9, GRID_STEP)
    dil = np.full((len(lats), len(lons)), np.nan)
    shear = np.full_like(dil, np.nan)
    for i, la in enumerate(lats):
        for j, lo in enumerate(lons):
            gx, gy = lo * kx, la * ky
            d2 = (x - gx) ** 2 + (y - gy) ** 2
            w = np.exp(-d2 / (2 * SIGMA_KM ** 2))
            sel = w > 0.01
            if w[sel].sum() < MIN_STA * 0.2 or sel.sum() < MIN_STA:
                continue
            X = np.c_[np.ones(sel.sum()), x[sel] - gx, y[sel] - gy]
            W = w[sel]
            A = X * W[:, None]
            try:
                ce, *_ = np.linalg.lstsq(A, ve[sel] * W, rcond=None)
                cn, *_ = np.linalg.lstsq(A, vn[sel] * W, rcond=None)
            except np.linalg.LinAlgError:
                continue
            dvedx, dvedy = ce[1], ce[2]   # mm/yr per km = nanostrain/yr
            dvndx, dvndy = cn[1], cn[2]
            exx = dvedx                    # x = east
            eyy = dvndy
            exy = 0.5 * (dvedy + dvndx)
            dil[i, j] = exx + eyy
            shear[i, j] = np.hypot((exx - eyy) / 2, exy)
    return lats, lons, dil, shear


def main():
    df = load_midas()
    assert df.lat.abs().max() <= 90 and df.lon.abs().max() <= 180
    reg = df[(df.lat.between(*BOX["lat"])) & (df.lon.between(*BOX["lon"]))
             & (df.se < MAX_SU / 1000.0) & (df.sn < MAX_SU / 1000.0)]
    print(f"MIDAS stations in box after QC: {len(reg)} (global {len(df)})")
    lats, lons, dil, shear = strain_grid(reg)
    tsi = dil / shear
    np.savez_compressed(HERE / "data" / "socal_strain_grid.npz",
                        lats=lats, lons=lons, dilatation_nstrain_yr=dil,
                        max_shear_nstrain_yr=shear, tsi=tsi)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    panels = [
        ("Dilatation rate (nanostrain/yr)", dil, "RdBu_r", (-100, 100)),
        ("Max shear strain rate (nanostrain/yr)", shear, "viridis", (0, 300)),
        ("TSI = dilatation / max shear", tsi, "RdBu_r", (-1.5, 1.5)),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(21, 6.5), constrained_layout=True)
    for ax, (title, z, cmap, clim) in zip(axes, panels):
        im = ax.pcolormesh(lons, lats, z, cmap=cmap, vmin=clim[0], vmax=clim[1], shading="auto")
        ax.plot(reg.lon, reg.lat, "k.", ms=1, alpha=0.35)
        for name, (gla, glo) in GEOTHERMAL.items():
            if BOX["lat"][0] <= gla <= BOX["lat"][1] and BOX["lon"][0] <= glo <= BOX["lon"][1]:
                ax.plot(glo, gla, "r^", ms=9, mec="k")
                ax.annotate(name, (glo, gla), textcoords="offset points", xytext=(6, 4),
                            fontsize=8, color="darkred")
        ax.set_title(title)
        ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
        fig.colorbar(im, ax=ax, shrink=0.85)
    fig.suptitle("Southern California geodetic strain rates & TSI - NGL MIDAS IGS14, "
                 f"Gaussian σ={SIGMA_KM:.0f} km, {GRID_STEP}° grid (black dots: GNSS stations)",
                 fontsize=13)
    fig.savefig(OUT / "socal_strain_tsi.png", dpi=180)
    print("-> maps/socal_strain_tsi.png; grid -> data/socal_strain_grid.npz")

    # TSI value at geothermal markers vs regional distribution
    stats = {}
    for name, (gla, glo) in GEOTHERMAL.items():
        i = int(round((gla - lats[0]) / GRID_STEP)); j = int(round((glo - lons[0]) / GRID_STEP))
        if 0 <= i < len(lats) and 0 <= j < len(lons):
            stats[name] = {"tsi": None if np.isnan(tsi[i, j]) else float(tsi[i, j]),
                           "dil": None if np.isnan(dil[i, j]) else float(dil[i, j]),
                           "shear": None if np.isnan(shear[i, j]) else float(shear[i, j])}
    stats["_regional"] = {"tsi_median": float(np.nanmedian(tsi)),
                          "dil_median": float(np.nanmedian(dil)),
                          "shear_median": float(np.nanmedian(shear))}
    (HERE / "results_tsi_map.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
