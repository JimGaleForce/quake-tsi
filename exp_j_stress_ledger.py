"""EXP-J: the forward stress ledger (PATTERN_PROTOCOL.md, section EXP-J, frozen 2026-08-09).

J1  per-0.2deg-cell geodetic loading rate vs seismic moment release (train period) -> chi,
    frozen classes SILENT-LOADING / COUPLED / OVERSHOOT.
J2  single test unblinding, two rival frozen hypotheses, both scored:
    (a) CATCH-UP   Mann-Whitney (one-sided, silent > coupled) on test/train event-rate ratios
    (b) PERSISTENCE Spearman log10(chi_train) vs log10(chi_test), cells >=20 events in both
J3  exploratory correlates of the negative space.

Outputs results_exp_j.json and maps/exp_j_ledger.png.  Reads only; modifies nothing.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

import tsi_map

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT_MAPS = HERE / "maps"

# ---------------------------------------------------------------- frozen params
BOX = {"lat": (31.5, 38.0), "lon": (-122.0, -113.5)}
CELL = 0.2
BIN_A = 0.4                      # exp_a tidal bins
MC = 2.5                         # completeness for the ledger
MAG_BIN = 0.1                    # -> Aki half-bin correction 0.05
MU = 30e9                        # Pa
H_SEIS = 11e3                    # m
TRAIN_END = pd.Timestamp("2010-01-01", tz="UTC")
TRAIN_YEARS = 29.0               # protocol-fixed train duration
MIN_EVENTS_CHI = 20
CHI_SILENT = 0.01
CHI_COUPLED_HI = 1.0
RATE_REGULARIZER = 0.5           # events/yr added to numerator AND denominator in J2a
EARTH_R = 6371.0
SAF_CREEP = ((36.0, -120.6), (36.8, -121.6))   # approximate creeping-section segment
SEC_PER_YEAR = 365.25 * 86400.0


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


def aki_b(mags):
    """Aki (1965) MLE b-value with half-bin correction; None if degenerate."""
    m = np.asarray(mags, dtype=float)
    if m.size < 20:
        return None
    denom = m.mean() - (MC - MAG_BIN / 2.0)
    if denom <= 0:
        return None
    return float(np.log10(np.e) / denom)


def med(x):
    x = np.asarray([v for v in x if v is not None and np.isfinite(v)], dtype=float)
    return float(np.median(x)) if x.size else None


# ---------------------------------------------------------------- J1: the ledger
def build_ledger():
    z = np.load(DATA / "socal_strain_grid.npz")
    glats, glons = z["lats"], z["lons"]
    shear = z["max_shear_nstrain_yr"]                       # nanostrain/yr

    # 0.2 deg cell edges aligned to the box corner
    lat_edges = np.arange(BOX["lat"][0], BOX["lat"][1] - 1e-9, CELL)
    lon_edges = np.arange(BOX["lon"][0], BOX["lon"][1] - 1e-9, CELL)

    # map each 0.1 deg node to its cell index (2x2 nodes per cell)
    ilat = np.floor((glats - BOX["lat"][0]) / CELL + 1e-6).astype(int)
    ilon = np.floor((glons - BOX["lon"][0]) / CELL + 1e-6).astype(int)

    nlat, nlon = len(lat_edges), len(lon_edges)
    ssum = np.zeros((nlat, nlon))
    scnt = np.zeros((nlat, nlon))
    ok = np.isfinite(shear)
    for a in range(len(glats)):
        if not (0 <= ilat[a] < nlat):
            continue
        row_ok = ok[a] & (ilon >= 0) & (ilon < nlon)
        np.add.at(ssum[ilat[a]], ilon[row_ok], shear[a][row_ok])
        np.add.at(scnt[ilat[a]], ilon[row_ok], 1.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        shear_cell = np.where(scnt > 0, ssum / np.maximum(scnt, 1), np.nan)

    rows = []
    for i in range(nlat):
        for j in range(nlon):
            if not np.isfinite(shear_cell[i, j]):
                continue
            lat0, lon0 = lat_edges[i], lon_edges[j]
            latc, lonc = lat0 + CELL / 2, lon0 + CELL / 2
            dy_m = CELL * 111_320.0
            dx_m = CELL * 111_320.0 * np.cos(np.radians(latc))
            rows.append(dict(i=i, j=j, lat0=float(lat0), lon0=float(lon0),
                             lat_c=float(latc), lon_c=float(lonc),
                             area_m2=float(dx_m * dy_m),
                             shear_nstrain_yr=float(shear_cell[i, j]),
                             n_nodes=int(scnt[i, j])))
    led = pd.DataFrame(rows)
    # Mdot_geo = 2 mu H A edot   (edot in 1/yr = nanostrain/yr * 1e-9)
    led["Mdot_geo"] = 2.0 * MU * H_SEIS * led.area_m2 * (led.shear_nstrain_yr * 1e-9)
    return led, lat_edges, lon_edges, shear_cell


def add_seismicity(led, cat, lat_edges, lon_edges):
    """Vectorized single pass: bin M>=2.5 events into cells for train and test."""
    nlat, nlon = len(lat_edges), len(lon_edges)
    c = cat[(cat.mag >= MC)
            & cat.lat.between(*BOX["lat"]) & cat.lon.between(*BOX["lon"])].copy()
    ci = np.floor((c.lat.to_numpy() - BOX["lat"][0]) / CELL + 1e-9).astype(int)
    cj = np.floor((c.lon.to_numpy() - BOX["lon"][0]) / CELL + 1e-9).astype(int)
    keep = (ci >= 0) & (ci < nlat) & (cj >= 0) & (cj < nlon)
    c = c.loc[keep].copy()
    c["cell"] = ci[keep] * nlon + cj[keep]
    c["m0"] = 10.0 ** (1.5 * c.mag.to_numpy() + 9.05)       # Hanks-Kanamori, N*m
    c["is_train"] = c.t < TRAIN_END

    t_end = c.t.max()
    test_years = float((t_end - TRAIN_END).total_seconds() / SEC_PER_YEAR)

    led = led.copy()
    led["cell"] = led.i.to_numpy() * nlon + led.j.to_numpy()

    for tag, sub in (("train", c[c.is_train]), ("test", c[~c.is_train])):
        g = sub.groupby("cell")
        agg = pd.DataFrame({f"n_{tag}": g.size(), f"m0_{tag}": g.m0.sum(),
                            f"depth_med_{tag}": g.depth.median(),
                            f"magsum_{tag}": g.mag.sum()})
        led = led.merge(agg, left_on="cell", right_index=True, how="left")
        led[f"n_{tag}"] = led[f"n_{tag}"].fillna(0).astype(int)
        led[f"m0_{tag}"] = led[f"m0_{tag}"].fillna(0.0)

    led["Mdot_seis_train"] = led.m0_train / TRAIN_YEARS
    led["Mdot_seis_test"] = led.m0_test / test_years
    led["chi_train"] = led.Mdot_seis_train / led.Mdot_geo
    led["chi_test"] = led.Mdot_seis_test / led.Mdot_geo
    led["rate_train"] = led.n_train / TRAIN_YEARS
    led["rate_test"] = led.n_test / test_years
    # Aki b-value per cell from train M>=2.5 mean magnitude
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_mag = led.magsum_train / led.n_train.replace(0, np.nan)
    led["b_train"] = np.where(led.n_train >= MIN_EVENTS_CHI,
                              np.log10(np.e) / (mean_mag - (MC - MAG_BIN / 2.0)), np.nan)
    led.loc[~np.isfinite(led.b_train), "b_train"] = np.nan
    return led, test_years, t_end


def classify(led):
    """Frozen classes.  SILENT-LOADING requires top-quartile loading and either chi<0.01 or
    an unmeasurable chi (n_train < 20) -- the latter recorded via `reason`."""
    q75 = float(np.nanpercentile(led.Mdot_geo, 75))
    top = led.Mdot_geo >= q75
    meas = led.n_train >= MIN_EVENTS_CHI

    cls = np.full(len(led), "UNCLASSIFIED", dtype=object)
    reason = np.full(len(led), "", dtype=object)

    silent_low_chi = top & meas & (led.chi_train < CHI_SILENT)
    silent_lowcnt = top & ~meas
    coupled = meas & (led.chi_train >= CHI_SILENT) & (led.chi_train <= CHI_COUPLED_HI)
    overshoot = meas & (led.chi_train > CHI_COUPLED_HI)

    cls[silent_low_chi.to_numpy()] = "SILENT-LOADING"; reason[silent_low_chi.to_numpy()] = "low_chi_top_quartile_loading"
    cls[silent_lowcnt.to_numpy()] = "SILENT-LOADING"; reason[silent_lowcnt.to_numpy()] = "low_count_top_quartile_loading"
    cls[coupled.to_numpy()] = "COUPLED"; reason[coupled.to_numpy()] = "chi_in_[0.01,1]"
    cls[overshoot.to_numpy()] = "OVERSHOOT"; reason[overshoot.to_numpy()] = "chi_gt_1"
    rest = cls == "UNCLASSIFIED"
    reason[rest & (~meas).to_numpy()] = "low_count_not_top_quartile_loading"
    reason[rest & meas.to_numpy()] = "low_chi_not_top_quartile_loading"

    led = led.copy()
    led["cls"] = cls
    led["reason"] = reason
    return led, q75


# ---------------------------------------------------------------- J3 correlates
def add_correlates(led):
    led = led.copy()
    led["d_saf_creep_km"] = dist_to_segment_km(led.lat_c.to_numpy(), led.lon_c.to_numpy(), SAF_CREEP)
    names = list(tsi_map.GEOTHERMAL)
    glat = np.array([tsi_map.GEOTHERMAL[n][0] for n in names])
    glon = np.array([tsi_map.GEOTHERMAL[n][1] for n in names])
    d = haversine_km(led.lat_c.to_numpy()[:, None], led.lon_c.to_numpy()[:, None],
                     glat[None, :], glon[None, :])
    led["d_geothermal_km"] = d.min(axis=1)
    led["nearest_geothermal"] = [names[k] for k in d.argmin(axis=1)]

    bins = pd.read_csv(HERE / "exp_a_train_bins.csv")
    key = {(round(r.bin_lat0, 3), round(r.bin_lon0, 3)): float(r.a_b) for r in bins.itertuples()}
    blat0 = BOX["lat"][0] + np.floor((led.lat_c - BOX["lat"][0]) / BIN_A) * BIN_A
    blon0 = BOX["lon"][0] + np.floor((led.lon_c - BOX["lon"][0]) / BIN_A) * BIN_A
    led["a_b"] = [key.get((round(a, 3), round(o, 3)), np.nan) for a, o in zip(blat0, blon0)]
    return led


def rank_block(led, mask_silent, cols):
    """Exploratory: Spearman vs log10(chi_train) and Mann-Whitney silent vs other."""
    out = {}
    lchi = np.log10(led.chi_train.where(led.chi_train > 0))
    for c in cols:
        v = led[c].to_numpy(dtype=float)
        entry = {}
        ok = np.isfinite(v) & np.isfinite(lchi.to_numpy()) & (led.n_train >= MIN_EVENTS_CHI).to_numpy()
        if ok.sum() >= 10:
            rho, p = spearmanr(v[ok], lchi.to_numpy()[ok])
            entry["spearman_vs_log_chi"] = {"rho": float(rho), "p": float(p), "n": int(ok.sum())}
        else:
            entry["spearman_vs_log_chi"] = {"note": f"insufficient n ({int(ok.sum())})"}
        a = v[mask_silent & np.isfinite(v)]
        b = v[(~mask_silent) & np.isfinite(v)]
        if a.size and b.size:
            u, p = mannwhitneyu(a, b, alternative="two-sided")
            entry["mannwhitney_silent_vs_other"] = {
                "n_silent": int(a.size), "n_other": int(b.size),
                "median_silent": float(np.median(a)), "median_other": float(np.median(b)),
                "U": float(u), "p_two_sided": float(p)}
        else:
            entry["mannwhitney_silent_vs_other"] = {
                "note": f"empty group (silent={a.size}, other={b.size})"}
        out[c] = entry
    return out


# ---------------------------------------------------------------- map
def draw_maps(led, lat_edges, lon_edges):
    OUT_MAPS.mkdir(exist_ok=True)
    nlat, nlon = len(lat_edges), len(lon_edges)
    ext = [lon_edges[0], lon_edges[-1] + CELL, lat_edges[0], lat_edges[-1] + CELL]

    def grid(col, transform=None):
        g = np.full((nlat, nlon), np.nan)
        v = led[col].to_numpy(dtype=float)
        if transform is not None:
            v = transform(v)
        g[led.i.to_numpy(), led.j.to_numpy()] = v
        return g

    with np.errstate(invalid="ignore", divide="ignore"):
        panels = [
            ("log10 loading rate  Mdot_geo (N m/yr)", grid("Mdot_geo", np.log10), "viridis"),
            ("log10 seismic release  Mdot_seis train (N m/yr)",
             grid("Mdot_seis_train", lambda x: np.log10(np.where(x > 0, x, np.nan))), "magma"),
            ("log10 coupling  chi = release / loading (train)",
             grid("chi_train", lambda x: np.log10(np.where(x > 0, x, np.nan))), "coolwarm"),
        ]

    classes = ["UNCLASSIFIED", "COUPLED", "OVERSHOOT", "SILENT-LOADING"]
    code = np.full((nlat, nlon), np.nan)
    for k, cname in enumerate(classes):
        m = (led.cls == cname).to_numpy()
        code[led.i.to_numpy()[m], led.j.to_numpy()[m]] = k

    fig, axes = plt.subplots(2, 2, figsize=(15, 13))
    for ax, (title, g, cmap) in zip(axes.ravel()[:3], panels):
        im = ax.imshow(g, origin="lower", extent=ext, aspect="auto", cmap=cmap)
        ax.set_title(title, fontsize=10)
        fig.colorbar(im, ax=ax, shrink=0.85)
    ax = axes.ravel()[3]
    cmap = matplotlib.colors.ListedColormap(["#d9d9d9", "#4575b4", "#7b3294", "#d73027"])
    im = ax.imshow(code, origin="lower", extent=ext, aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, ticks=range(4))
    cb.ax.set_yticklabels(classes, fontsize=8)
    ax.set_title("Frozen ledger classes (red triangles: geothermal fields)", fontsize=10)
    for name, (gla, glo) in tsi_map.GEOTHERMAL.items():
        if BOX["lat"][0] <= gla <= BOX["lat"][1] and BOX["lon"][0] <= glo <= BOX["lon"][1]:
            ax.plot(glo, gla, "r^", ms=9, mec="k")
            ax.annotate(name, (glo, gla), textcoords="offset points", xytext=(6, 4),
                        fontsize=8, color="darkred")
    for ax in axes.ravel():
        ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    fig.suptitle("EXP-J forward stress ledger - SoCal 0.2 deg cells, train 1981-2009 "
                 f"(mu={MU/1e9:.0f} GPa, H={H_SEIS/1e3:.0f} km, M>={MC})", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_MAPS / "exp_j_ledger.png", dpi=170)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    print("EXP-J forward stress ledger")
    led, lat_edges, lon_edges, _ = build_ledger()
    n_strain_cells = len(led)
    assert n_strain_cells > 0, "no 0.2 deg cells with finite strain"
    print(f"  cells with finite strain: {n_strain_cells} "
          f"({len(lat_edges)}x{len(lon_edges)} candidate grid)")

    cat = load_catalog("xue_lu_zenodo/SCSN_original_catalog.txt")
    print(f"  catalog events: {len(cat)}  {cat.t.min().date()} -> {cat.t.max().date()}")
    led, test_years, t_end = add_seismicity(led, cat, lat_edges, lon_edges)
    led, q75 = classify(led)
    led = add_correlates(led)

    counts = led.cls.value_counts().to_dict()
    assert sum(counts.values()) == len(led), "class counts do not sum to ledger cells"
    print(f"  train events in cells: {int(led.n_train.sum())}, "
          f"test: {int(led.n_test.sum())} over {test_years:.2f} yr")
    for k, v in sorted(counts.items()):
        print(f"    {k:16s} {v}")

    per_class = {}
    for cname, sub in led.groupby("cls"):
        per_class[cname] = {
            "n": int(len(sub)),
            "median_Mdot_geo_Nm_per_yr": float(np.median(sub.Mdot_geo)),
            "median_Mdot_seis_train_Nm_per_yr": float(np.median(sub.Mdot_seis_train)),
            "median_chi_train": (float(np.median(sub.chi_train[sub.chi_train > 0]))
                                 if (sub.chi_train > 0).any() else None),
            "median_n_train": float(np.median(sub.n_train)),
            "reasons": {k: int(v) for k, v in sub.reason.value_counts().items()},
        }

    silent = led.cls == "SILENT-LOADING"
    coupled = led.cls == "COUPLED"
    top10 = (led[silent].sort_values("Mdot_geo", ascending=False).head(10)
             [["lat_c", "lon_c", "Mdot_geo", "Mdot_seis_train", "chi_train",
               "n_train", "n_test", "reason"]])
    j1 = {
        "n_cells": int(n_strain_cells),
        "cell_deg": CELL,
        "loading_top_quartile_threshold_Nm_per_yr": q75,
        "train_years": TRAIN_YEARS,
        "test_years": test_years,
        "test_end": str(t_end),
        "class_counts": {k: int(v) for k, v in counts.items()},
        "per_class": per_class,
        "top10_silent_loading_by_loading": [
            {"lat": float(r.lat_c), "lon": float(r.lon_c),
             "Mdot_geo_Nm_per_yr": float(r.Mdot_geo),
             "Mdot_seis_train_Nm_per_yr": float(r.Mdot_seis_train),
             "chi_train": float(r.chi_train), "n_train": int(r.n_train),
             "n_test": int(r.n_test), "reason": r.reason}
            for r in top10.itertuples()],
    }

    # ---------------- J2a CATCH-UP
    R = RATE_REGULARIZER
    ratio = (led.rate_test + R) / (led.rate_train + R)
    a = ratio[silent].to_numpy()
    b = ratio[coupled].to_numpy()
    if a.size == 0 or b.size == 0:
        catch_up = {"ERROR": "EMPTY CLASS - test not run",
                    "n_silent": int(a.size), "n_coupled": int(b.size)}
        print("  !! J2a: EMPTY CLASS, Mann-Whitney not run")
    else:
        u, p = mannwhitneyu(a, b, alternative="greater")
        catch_up = {
            "hypothesis": "train SILENT-LOADING cells have higher test/train event-rate ratios "
                          "than COUPLED cells (one-sided, silent > coupled)",
            "rate_regularizer_events_per_yr": R,
            "n_silent": int(a.size), "n_coupled": int(b.size),
            "median_ratio_silent": float(np.median(a)),
            "median_ratio_coupled": float(np.median(b)),
            "median_rate_train_silent": float(np.median(led.rate_train[silent])),
            "median_rate_test_silent": float(np.median(led.rate_test[silent])),
            "median_rate_train_coupled": float(np.median(led.rate_train[coupled])),
            "median_rate_test_coupled": float(np.median(led.rate_test[coupled])),
            "U": float(u), "p_one_sided": float(p),
            "supported": bool(p < 0.05),
        }
        print(f"  J2a CATCH-UP: n_silent={a.size} n_coupled={b.size} "
              f"median ratio {np.median(a):.3f} vs {np.median(b):.3f}  p={p:.3g}")

    # ---------------- J2b PERSISTENCE
    both = (led.n_train >= MIN_EVENTS_CHI) & (led.n_test >= MIN_EVENTS_CHI) \
        & (led.chi_train > 0) & (led.chi_test > 0)
    nb = int(both.sum())
    if nb < 10:
        persistence = {"ERROR": f"only {nb} cells measurable in both periods - test not run"}
        print(f"  !! J2b: only {nb} cells measurable in both periods")
    else:
        x = np.log10(led.chi_train[both].to_numpy())
        y = np.log10(led.chi_test[both].to_numpy())
        rho, p = spearmanr(x, y)
        persistence = {
            "hypothesis": "cell character persists: log10(chi_train) correlates with log10(chi_test)",
            "min_events_each_period": MIN_EVENTS_CHI,
            "n": nb, "spearman_rho": float(rho), "p": float(p),
            "supported": bool(p < 0.05 and rho > 0),
        }
        print(f"  J2b PERSISTENCE: n={nb} rho={rho:.3f} p={p:.3g}")

    cu_ok = bool(catch_up.get("supported", False))
    pe_ok = bool(persistence.get("supported", False))
    verdict = (
        "both supported: persistence of coupling AND elevated catch-up in silent-loading cells"
        if cu_ok and pe_ok else
        "PERSISTENCE wins (prior confirmed); no catch-up signal" if pe_ok else
        "CATCH-UP wins (prior overturned); coupling does not persist" if cu_ok else
        "neither hypothesis supported")

    # ---------------- J3
    led["depth_med_train"] = led.depth_med_train.astype(float)
    cols = ["d_saf_creep_km", "d_geothermal_km", "depth_med_train", "b_train", "a_b"]
    j3 = {"label": "EXPLORATORY - not pre-registered as confirmatory; no multiplicity control",
          "n_cells_with_a_b": int(np.isfinite(led.a_b).sum()),
          "correlates": rank_block(led, silent.to_numpy(), cols)}

    res = {
        "experiment": "EXP-J forward stress ledger",
        "protocol": "PATTERN_PROTOCOL.md section EXP-J (frozen 2026-08-09)",
        "params": {"mu_Pa": MU, "H_m": H_SEIS, "Mc": MC, "cell_deg": CELL,
                   "train_end_utc": str(TRAIN_END), "min_events_chi": MIN_EVENTS_CHI,
                   "chi_silent_lt": CHI_SILENT, "chi_coupled_hi": CHI_COUPLED_HI,
                   "rate_regularizer": RATE_REGULARIZER,
                   "saf_creeping_segment": SAF_CREEP},
        "j1": j1,
        "j2": {"catch_up": catch_up, "persistence": persistence, "verdict": verdict},
        "j3": j3,
    }
    (HERE / "results_exp_j.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    draw_maps(led, lat_edges, lon_edges)
    print(f"  VERDICT: {verdict}")
    print("-> results_exp_j.json, maps/exp_j_ledger.png")


if __name__ == "__main__":
    main()
