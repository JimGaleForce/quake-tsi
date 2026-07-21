"""EQ-18/EQ-22 — full SoCal 0.2-degree bin scan of tidal phase modulation (Lu Fig 4 layout).

Purpose (dual):
 1. Sensitivity diagnostic: if ANY bin reaches significance with our implementation, the pipeline
    can detect real modulations and the Coso/Long Valley nulls become interpretable.
 2. Label map: per-bin Pm/P0 + significance is the training label set for the supervised
    susceptibility approach; immediately correlated here against the geodetic strain grid.

Method: per-event focal-mechanism sigma_n phases (YSH/SCEDC mechanisms, preferred plane),
orientation groups rounded to 15 deg; per-group phase series from the Lu/Xue SPOTL-derived
stress tensor (upsample x5 -> 1200 s); per-bin significance vs 1,000 orientation-preserving
synthetic catalogs (reservoir sampling per group). Bins: 0.2 deg, min 100 FM events.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.signal import find_peaks
from scipy.interpolate import CubicSpline
import importlib.util

HERE = Path(__file__).parent
spec = importlib.util.spec_from_file_location("cc", HERE / "coso_positive_control.py")
cc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cc)
spec2 = importlib.util.spec_from_file_location("fmt", HERE / "coso_fm_test.py")
fmt = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(fmt)

UP = 5
DT = 6000.0 / UP
N_BINS = 36
N_SYNTH = 1000
RESERVOIR = 2000
BIN_DEG = 0.2
MIN_EVENTS = 100
SEED = 20260725


def base_series_up5():
    E, NU = 75e9, 0.25
    G = E / (2 * (1 + NU))
    LAM = E * NU / ((1 + NU) * (1 - 2 * NU))
    D = cc.DATA
    exx = np.loadtxt(D / "Tidal_N_0.txt") * 1e-9
    eyy = np.loadtxt(D / "Tidal_N_90.txt") * 1e-9
    exy = np.loadtxt(D / "Tidal_S_0.txt") * 1e-9
    evol = np.loadtxt(D / "Tidal_Vol.txt") * 1e-9
    ezz = evol - exx - eyy
    sxx = (LAM + 2 * G) * exx + LAM * eyy + LAM * ezz
    syy = LAM * exx + (LAM + 2 * G) * eyy + LAM * ezz
    szz = LAM * exx + LAM * eyy + (LAM + 2 * G) * ezz
    sxy = 2 * G * exy
    n = len(sxx)
    fine = np.arange(0, n - 1, 1.0 / UP)
    up = lambda v: CubicSpline(np.arange(n), v)(fine).astype(np.float32)
    return up(sxx), up(syy), up(sxy), up(szz)


def fit_pm_bins(hist_rows):
    """Vectorized Pm/P0 fit for an array of histograms (rows)."""
    d = hist_rows / hist_rows.sum(axis=1, keepdims=True) * N_BINS
    th = np.radians(np.linspace(-180, 180, N_BINS + 1)[:-1] + 180 / N_BINS)
    A = np.c_[np.ones_like(th), np.cos(th), np.sin(th)]
    coef = np.linalg.lstsq(A, d.T, rcond=None)[0]  # 3 x rows
    return np.hypot(coef[1], coef[2]) / coef[0]


def main():
    rng = np.random.default_rng(SEED)
    fm = fmt.load_fm()
    SXX, SYY, SXY, SZZ = base_series_up5()
    nfine = len(SXX)
    print(f"FM catalog: {len(fm)}; series samples: {nfine}")

    results = {}
    for cat_name, cat_file in [("QTM", "QTM_decluster_m0.1.txt"), ("SCSN", "SCSN_decluster_m1.5.txt")]:
        cat = cc.load_declustered(cat_file)
        m = cat.merge(fm[["eid", "strike", "dip", "rake"]], on="eid", how="inner").reset_index(drop=True)
        print(f"{cat_name}: {len(cat)} declustered, {len(m)} FM-matched")
        for c in ["strike", "dip", "rake"]:
            m[f"{c}_r"] = (np.round(m[c] / 15) * 15).astype(int)
        idx_ev = np.clip(((m.t_unix.to_numpy() - cc.T0.timestamp()) / DT).astype(np.int64), 0, nfine - 1)

        ev_phase = np.full(len(m), np.nan, dtype=np.float32)
        groups = list(m.groupby(["strike_r", "dip_r", "rake_r"]).groups.items())
        print(f"  orientation groups: {len(groups)}")
        reservoirs = np.full((len(groups), RESERVOIR), np.nan, dtype=np.float32)
        ev_group = np.full(len(m), -1, dtype=np.int32)
        for gi, ((st, di, ra), gidx) in enumerate(groups):
            cn, _ = fmt.combo_coeffs(st, di, ra)
            series = cn[0] * SXX + cn[1] * SYY + cn[2] * SXY + cn[3] * SZZ
            phase = cc.phase_series(series)
            gpos = m.index.get_indexer(gidx)
            ev_phase[gpos] = phase[idx_ev[gpos]]
            ev_group[gpos] = gi
            valid = np.flatnonzero(~np.isnan(phase))
            reservoirs[gi] = phase[rng.choice(valid, RESERVOIR)]
            if gi % 200 == 0:
                print(f"  group {gi}/{len(groups)}")

        ok = ~np.isnan(ev_phase)
        m, ev_phase, ev_group = m[ok].reset_index(drop=True), ev_phase[ok], ev_group[ok]
        m["blat"] = np.floor(m.lat / BIN_DEG) * BIN_DEG
        m["blon"] = np.floor(m.lon / BIN_DEG) * BIN_DEG

        rows = []
        edges = np.linspace(-180, 180, N_BINS + 1)
        for (bl, bo), bidx in m.groupby(["blat", "blon"]).groups.items():
            gpos = m.index.get_indexer(bidx)
            if len(gpos) < MIN_EVENTS:
                continue
            obs_h, _ = np.histogram(ev_phase[gpos], bins=edges)
            obs = fit_pm_bins(obs_h[None, :])[0]
            gr = ev_group[gpos]
            draws = reservoirs[gr[:, None], rng.integers(0, RESERVOIR, size=(len(gpos), N_SYNTH))]
            binned = np.clip(((draws + 180) / (360 / N_BINS)).astype(np.int64), 0, N_BINS - 1)
            syn_h = np.zeros((N_SYNTH, N_BINS))
            for b in range(N_BINS):
                syn_h[:, b] = (binned == b).sum(axis=0)
            syn = fit_pm_bins(syn_h)
            thr = float(np.percentile(syn, 97.5))
            rows.append({"lat": float(bl + BIN_DEG / 2), "lon": float(bo + BIN_DEG / 2),
                         "n": int(len(gpos)), "Pm_over_P0": float(obs), "thr97.5": thr,
                         "significant": bool(obs > thr),
                         "p": float((np.sum(syn >= obs) + 1) / (N_SYNTH + 1))})
        bins_df = pd.DataFrame(rows)
        nsig = int(bins_df.significant.sum()) if len(bins_df) else 0
        print(f"  {cat_name}: {len(bins_df)} bins with n>={MIN_EVENTS}; significant at 95%: {nsig}")
        if len(bins_df):
            print(bins_df.sort_values("p").head(8).to_string(index=False))
        bins_df.to_csv(HERE / f"binscan_{cat_name}.csv", index=False)
        results[cat_name] = {"n_bins": len(bins_df), "n_significant": nsig}

    # correlate with strain grid
    z = np.load(HERE / "data" / "socal_strain_grid.npz")
    glats, glons = z["lats"], z["lons"]
    tsi_g, dil_g, sh_g = z["tsi"], z["dilatation_nstrain_yr"], z["max_shear_nstrain_yr"]
    from scipy.stats import spearmanr
    for cat_name in results:
        bins_df = pd.read_csv(HERE / f"binscan_{cat_name}.csv")
        if len(bins_df) < 5:
            continue
        ii = np.clip(np.round((bins_df.lat - glats[0]) / 0.1).astype(int), 0, len(glats) - 1)
        jj = np.clip(np.round((bins_df.lon - glons[0]) / 0.1).astype(int), 0, len(glons) - 1)
        for name, g in [("TSI", tsi_g), ("dilatation", dil_g), ("max_shear", sh_g)]:
            v = g[ii, jj]
            okv = ~np.isnan(v)
            if okv.sum() > 4:
                rho, p = spearmanr(v[okv], bins_df.Pm_over_P0[okv])
                print(f"{cat_name}: Spearman({name}, Pm/P0) rho={rho:.3f} p={p:.3f} (n={okv.sum()})")
                results[cat_name][f"spearman_{name}"] = {"rho": float(rho), "p": float(p)}

    (HERE / "results_bin_scan.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print("-> results_bin_scan.json, binscan_QTM.csv, binscan_SCSN.csv")


if __name__ == "__main__":
    main()
