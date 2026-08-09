"""EXP-C (exploratory): time-varying susceptibility as a stress gauge.

Frozen spec: OVERNIGHT_PREDICTION_PROTOCOL.md, EXP-C section.

For the 2 bins (0.4 deg x 0.4 deg, grid aligned to -122.0+0.4k / 31.5+0.4k) with the highest
total FM-matched event counts over the full 1981-2018 SCSN declustered M>=1.5 catalog:
sliding 3-year windows stepped 1 year (window start 1981..2015). In each window with
>= 60 FM-matched events, compute the per-event tau phase (per-event FM, orientation groups
rounded to 15 deg, reusing coso_fm_test.py's validated machinery) and a 36-bin sinusoid fit
-> Pm/P0(t), phi(t). Also record next-calendar-year event count (all declustered events in
the bin, not just FM-matched) and whether any M>=4 occurred in the bin in that next year.

This is EXPLORATORY / DESCRIPTIVE ONLY per the frozen protocol - no significance testing,
multiple-comparison burden acknowledged. Existing machinery (coso_fm_test.py,
coso_positive_control.py) is imported via importlib and NOT modified.
"""
import json
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

HERE = Path(__file__).parent

spec = importlib.util.spec_from_file_location("fmtest", HERE / "coso_fm_test.py")
fmtest = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fmtest)
cc = fmtest.cc  # coso_positive_control, loaded by coso_fm_test.py

# ---- config (frozen protocol values) ----
LAT0, LON0, BINSIZE = 31.5, -122.0, 0.4
SOCAL_BOX = dict(lat=(31.5, 38.0), lon=(-122.0, -113.5))
WINDOW_YEARS = 3
STEP_YEARS = 1
MIN_EVENTS_PER_WINDOW = 60
N_BINS = 36
CATALOG_FILE = "SCSN_decluster_m1.5.txt"
WINDOW_STARTS = list(range(1981, 2016))  # 1981..2015 inclusive, per spec

OUT_JSON = HERE / "results_exp_c.json"
OUT_PNG = HERE / "maps" / "exp_c_susceptibility_drift.png"


def bin_index(lat, lon):
    ilat = np.floor((lat - LAT0) / BINSIZE).astype(int)
    ilon = np.floor((lon - LON0) / BINSIZE).astype(int)
    return ilat, ilon


def bin_center(ilat, ilon):
    return (LAT0 + BINSIZE * ilat + BINSIZE / 2.0, LON0 + BINSIZE * ilon + BINSIZE / 2.0)


def fit36(hist):
    """Unweighted 36-bin sinusoid fit, identical in form to coso_fm_test.run's inner `fit`."""
    d = hist / hist.sum() * N_BINS
    th = np.radians(np.linspace(-180, 180, N_BINS + 1)[:-1] + 180 / N_BINS)
    A = np.c_[np.ones_like(th), np.cos(th), np.sin(th)]
    coef, *_ = np.linalg.lstsq(A, d, rcond=None)
    pm_p0 = np.hypot(coef[1], coef[2]) / coef[0]
    phi = np.degrees(np.arctan2(coef[2], coef[1]))
    return float(pm_p0), float(phi)


def main():
    # ---- load catalog + FM matches ----
    cat = cc.load_declustered(CATALOG_FILE)
    cat = cat[cat.lat.between(*SOCAL_BOX["lat"]) & cat.lon.between(*SOCAL_BOX["lon"])].copy()
    ilat, ilon = bin_index(cat.lat.to_numpy(), cat.lon.to_numpy())
    cat["ilat"], cat["ilon"] = ilat, ilon
    cat["year"] = cat.yr.astype(int)
    cat["time"] = pd.to_datetime(dict(year=cat.yr, month=cat.mo, day=cat.dy, hour=cat.hr,
                                       minute=cat.mi, second=0), utc=True) + \
                  pd.to_timedelta(cat.sec.astype(float), unit="s")

    fm = fmtest.load_fm()
    m = cat.merge(fm[["eid", "strike", "dip", "rake"]], on="eid", how="inner").copy()
    print(f"Catalog (SoCal box) events: {len(cat)}, FM-matched: {len(m)}")

    # ---- find the two bins with highest total FM-matched counts, full 1981-2018 period ----
    counts = m.groupby(["ilat", "ilon"]).size().sort_values(ascending=False)
    top2 = counts.index[:2].tolist()
    print("Top FM-matched bins (ilat, ilon, n, lat_center, lon_center):")
    for (ila, ilo) in top2:
        clat, clon = bin_center(ila, ilo)
        print(f"  ({ila},{ilo}) n={counts[(ila, ilo)]}  center=({clat:.2f}, {clon:.2f})")

    # ---- precompute upsampled base tidal-stress series (shared across all bins/windows) ----
    SXX, SYY, SXY, SZZ = fmtest.base_series()
    nfine = len(SXX)

    results = {"protocol": "OVERNIGHT_PREDICTION_PROTOCOL.md#EXP-C", "exploratory": True,
               "window_years": WINDOW_YEARS, "step_years": STEP_YEARS,
               "min_events_per_window": MIN_EVENTS_PER_WINDOW, "bins": []}

    fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True)

    for panel_i, (ila, ilo) in enumerate(top2):
        clat, clon = bin_center(ila, ilo)
        bin_fm = m[(m.ilat == ila) & (m.ilon == ilo)].copy()
        bin_all = cat[(cat.ilat == ila) & (cat.ilon == ilo)].copy()
        n_total_fm = len(bin_fm)
        print(f"\n--- Bin ({ila},{ilo}) center=({clat:.2f},{clon:.2f}): "
              f"{n_total_fm} FM-matched events, {len(bin_all)} total declustered events ---")

        # per-event tau phase, orientation groups rounded to 15 deg (matches coso_fm_test.run)
        bin_fm["strike_r"] = (np.round(bin_fm.strike / fmtest.ROUND_DEG) * fmtest.ROUND_DEG).astype(int)
        bin_fm["dip_r"] = (np.round(bin_fm.dip / fmtest.ROUND_DEG) * fmtest.ROUND_DEG).astype(int)
        bin_fm["rake_r"] = (np.round(bin_fm.rake / fmtest.ROUND_DEG) * fmtest.ROUND_DEG).astype(int)
        t_unix = (bin_fm["time"] - pd.Timestamp(0, tz="UTC")).dt.total_seconds().to_numpy()
        idx_ev = np.clip(((t_unix - cc.T0.timestamp()) / cc.DT).astype(np.int64), 0, nfine - 1)

        ph_all = np.full(len(bin_fm), np.nan)
        groups = list(bin_fm.groupby(["strike_r", "dip_r", "rake_r"]).groups.items())
        print(f"  orientation groups: {len(groups)}")
        for (st, di, ra), gidx in groups:
            cn, ctc = fmtest.combo_coeffs(st, di, ra)
            c = ctc  # tau (shear stress) component, per protocol
            series = c[0] * SXX + c[1] * SYY + c[2] * SXY + c[3] * SZZ
            phase = cc.phase_series(series)
            gpos = bin_fm.index.get_indexer(gidx)
            ph_all[gpos] = phase[idx_ev[gpos]]
        bin_fm["phase"] = ph_all
        bin_fm["year"] = bin_fm["time"].dt.year

        series_out = []
        for ws in WINDOW_STARTS:
            we = ws + WINDOW_YEARS  # window = [ws, we)
            w_start_ts = pd.Timestamp(f"{ws}-01-01", tz="UTC")
            w_end_ts = pd.Timestamp(f"{we}-01-01", tz="UTC")
            wmask = (bin_fm["time"] >= w_start_ts) & (bin_fm["time"] < w_end_ts)
            wevents = bin_fm[wmask & bin_fm.phase.notna()]
            n_fm = int(len(wevents))

            ny_start = pd.Timestamp(f"{we}-01-01", tz="UTC")
            ny_end = pd.Timestamp(f"{we + 1}-01-01", tz="UTC")
            nymask = (bin_all["time"] >= ny_start) & (bin_all["time"] < ny_end)
            next_year_count = int(nymask.sum())
            next_year_m4 = bool((bin_all.loc[nymask, "mag"] >= 4.0).any())

            entry = {"window_start": ws, "n_fm": n_fm,
                     "next_year_count": next_year_count, "next_year_m4": next_year_m4}
            if n_fm >= MIN_EVENTS_PER_WINDOW:
                hist, _ = np.histogram(wevents.phase.to_numpy(), bins=np.linspace(-180, 180, N_BINS + 1))
                pm_p0, phi = fit36(hist)
                entry["pm_p0"] = pm_p0
                entry["phi"] = phi
            else:
                entry["pm_p0"] = None
                entry["phi"] = None
            series_out.append(entry)

        valid = [e for e in series_out if e["pm_p0"] is not None]
        if len(valid) >= 3:
            pm_arr = [e["pm_p0"] for e in valid]
            ny_arr = [e["next_year_count"] for e in valid]
            rho, pval = spearmanr(pm_arr, ny_arr)
        else:
            rho, pval = float("nan"), float("nan")

        print(f"  windows with n_fm>={MIN_EVENTS_PER_WINDOW}: {len(valid)} / {len(series_out)}")
        print(f"  Spearman(Pm/P0, next_year_count) = {rho:.3f} (p={pval:.3f}, n={len(valid)}) "
              f"[EXPLORATORY - not a significance claim]")

        results["bins"].append({
            "ilat": int(ila), "ilon": int(ilo),
            "lat_center": clat, "lon_center": clon,
            "n_total_fm_matched": n_total_fm,
            "n_total_declustered": int(len(bin_all)),
            "spearman_rho_pmp0_vs_next_year_count": None if np.isnan(rho) else float(rho),
            "spearman_p_exploratory": None if np.isnan(pval) else float(pval),
            "n_windows_valid": len(valid),
            "series": series_out,
        })

        # ---- plot ----
        ax = axes[panel_i]
        ws_arr = np.array([e["window_start"] for e in series_out])
        pm_plot = np.array([e["pm_p0"] if e["pm_p0"] is not None else np.nan for e in series_out])
        ny_plot = np.array([e["next_year_count"] for e in series_out])
        m4_mask = np.array([e["next_year_m4"] for e in series_out])

        ax2 = ax.twinx()
        ax2.bar(ws_arr + WINDOW_YEARS, ny_plot, width=0.8, color="tab:gray", alpha=0.4,
                label="next-year event count (all events)")
        ax.plot(ws_arr, pm_plot, "o-", color="tab:blue", label="Pm/P0(t) (window start)")
        for ws_m4, ny_m4 in zip(ws_arr[m4_mask], ny_plot[m4_mask]):
            ax2.axvline(ws_m4 + WINDOW_YEARS, color="tab:red", alpha=0.5, linestyle="--", linewidth=1)
        ax.set_ylabel("Pm/P0 (3-yr window)", color="tab:blue")
        ax2.set_ylabel("next-year event count / M>=4 (dashed red)", color="tab:gray")
        ax.set_title(f"Bin ({clat:.2f}N, {clon:.2f}E)  n_FM_total={n_total_fm}  "
                      f"rho={rho:.2f} (exploratory)")
        ax.set_xlabel("window start year")
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc="upper left", fontsize=8)

    fig.suptitle("EXP-C (exploratory): time-varying tidal-phase susceptibility vs next-year "
                  "seismicity, two highest FM-count bins")
    fig.tight_layout()
    OUT_PNG.parent.mkdir(exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150)
    print(f"\n-> {OUT_PNG}")

    OUT_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"-> {OUT_JSON}")

    print("\n=== SUMMARY ===")
    for b in results["bins"]:
        print(f"Bin ({b['lat_center']:.2f}, {b['lon_center']:.2f}): "
              f"n_fm_total={b['n_total_fm_matched']}, "
              f"valid windows={b['n_windows_valid']}, "
              f"Spearman rho={b['spearman_rho_pmp0_vs_next_year_count']}, "
              f"p(exploratory)={b['spearman_p_exploratory']}")


if __name__ == "__main__":
    main()
