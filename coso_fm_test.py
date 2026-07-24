"""EQ-18 step A (part 3) - Coso positive control with PER-EVENT focal mechanisms (Lu et al. §2.5).

Each declustered event matched (by SCSN evid) to the Yang-Hauksson-Shearer / SCEDC focal
mechanism catalog gets the tidal stress tensor resolved onto its own mechanism (preferred plane).
Phases are assigned per event from that plane's stress series. Following Lu et al.: no background
bias correction is possible (heterogeneous orientations), so significance is evaluated against
3,000 synthetic catalogs that KEEP the orientation set and randomize event times.

sigma_n and tau are linear combinations of the base stress series (sxx, syy, sxy, szz), so we
precompute the four upsampled base series once and combine per orientation group (angles rounded
to 15 deg) - one find_peaks pass per group.
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

FM_DIR = HERE / "data" / "scedc_fm"
COSO = cc.COSO_BOX
N_SYNTH = 3000
N_BINS = 36
ROUND_DEG = 15
SEED = 20260722


def load_fm():
    rows = []
    files = [FM_DIR / "YSH_2010.hash"] + sorted(FM_DIR.glob("sc*.focmec"))
    for f in files:
        for line in open(f, encoding="utf-8", errors="replace"):
            p = line.split()
            if len(p) < 14:
                continue
            try:
                evid = int(p[6])
                lat, lon = float(p[7]), float(p[8])
                strike, dip, rake = float(p[11]), float(p[12]), float(p[13])
            except ValueError:
                continue
            rows.append((evid, lat, lon, strike, dip, rake))
    fm = pd.DataFrame(rows, columns=["eid", "lat", "lon", "strike", "dip", "rake"])
    return fm.drop_duplicates("eid")


def base_series():
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
    fine = np.arange(0, n - 1, 1.0 / cc.UPSAMPLE)
    up = lambda v: CubicSpline(np.arange(n), v)(fine).astype(np.float32)
    return up(sxx), up(syy), up(sxy), up(szz)


def combo_coeffs(strike, dip, rake):
    sr, s90, dr, d90, rr = map(np.radians, [strike, strike + 90, dip, 90 - dip, rake])
    n = np.array([np.cos(s90) * np.cos(d90), np.sin(s90) * np.cos(d90), -np.sin(d90)])
    s = np.array([np.cos(sr), np.sin(sr), 0.0])
    sd = np.array([np.cos(s90) * np.cos(dr), np.sin(s90) * np.cos(dr), np.sin(dr)])
    u = s * np.cos(rr) + sd * np.sin(rr)
    cn = (n[0] ** 2, n[1] ** 2, 2 * n[0] * n[1], n[2] ** 2)
    ct = (n[0] * u[0], n[1] * u[1], n[1] * u[0] + n[0] * u[1], n[2] * u[2])
    return cn, ct


def run(component_key):
    rng = np.random.default_rng(SEED)
    fm = load_fm()
    df = cc.load_declustered("QTM_decluster_m0.1.txt")
    df2 = cc.load_declustered("SCSN_decluster_m1.5.txt")
    results = {}
    SXX, SYY, SXY, SZZ = base_series()
    nfine = len(SXX)

    for cat_name, cat in [("QTM", df), ("SCSN", df2)]:
        coso = cat[(cat.lat.between(*COSO["lat"])) & (cat.lon.between(*COSO["lon"]))]
        m = coso.merge(fm[["eid", "strike", "dip", "rake"]], on="eid", how="inner")
        print(f"{cat_name}: Coso events {len(coso)}, with focal mechanisms {len(m)}")
        if len(m) < 50:
            results[cat_name] = {"n": int(len(m)), "note": "too few FM matches"}
            continue
        m = m.copy()
        for c in ["strike", "dip", "rake"]:
            m[f"{c}_r"] = (np.round(m[c] / ROUND_DEG) * ROUND_DEG).astype(int)
        idx_ev = np.clip(((m.t_unix.to_numpy() - cc.T0.timestamp()) / cc.DT).astype(np.int64),
                         0, nfine - 1)

        obs_hist = np.zeros(N_BINS)
        syn_hist = np.zeros((N_SYNTH, N_BINS))
        n_used = 0
        groups = list(m.groupby(["strike_r", "dip_r", "rake_r"]).groups.items())
        print(f"  orientation groups: {len(groups)}")
        for (st, di, ra), gidx in groups:
            cn, ct = combo_coeffs(st, di, ra)
            c = cn if component_key == "sigma_n" else ct
            series = c[0] * SXX + c[1] * SYY + c[2] * SXY + c[3] * SZZ
            phase = cc.phase_series(series)
            gpos = m.index.get_indexer(gidx)
            ph = phase[idx_ev[gpos]]
            ok = ~np.isnan(ph)
            n_used += ok.sum()
            h, _ = np.histogram(ph[ok], bins=np.linspace(-180, 180, N_BINS + 1))
            obs_hist += h
            valid = np.flatnonzero(~np.isnan(phase))
            ridx = rng.choice(valid, size=(N_SYNTH, int(ok.sum())))
            rph = phase[ridx]
            binned = np.clip(((rph + 180) / (360 / N_BINS)).astype(np.int64), 0, N_BINS - 1)
            for b in range(N_BINS):
                syn_hist[:, b] += (binned == b).sum(axis=1)

        def fit(h):
            d = h / h.sum() * N_BINS
            th = np.radians(np.linspace(-180, 180, N_BINS + 1)[:-1] + 180 / N_BINS)
            A = np.c_[np.ones_like(th), np.cos(th), np.sin(th)]
            coef, *_ = np.linalg.lstsq(A, d, rcond=None)
            return np.hypot(coef[1], coef[2]) / coef[0], np.degrees(np.arctan2(coef[2], coef[1]))

        obs, phi = fit(obs_hist)
        syn = np.array([fit(syn_hist[k])[0] for k in range(N_SYNTH)])
        r = {"n": int(n_used), "Pm_over_P0": float(obs), "phi_deg": float(phi),
             "synth_p97.5": float(np.percentile(syn, 97.5)),
             "synth_p99.5": float(np.percentile(syn, 99.5)),
             "significant_95": bool(obs > np.percentile(syn, 97.5)),
             "empirical_p": float((np.sum(syn >= obs) + 1) / (N_SYNTH + 1))}
        results[cat_name] = r
        print(f"  {cat_name} {component_key}: n={r['n']} Pm/P0={r['Pm_over_P0']:.3f} "
              f"thr97.5={r['synth_p97.5']:.3f} p={r['empirical_p']:.4f} phi={r['phi_deg']:.0f}"
              f"{' <-- SIGNIFICANT' if r['significant_95'] else ''}")
    return results


if __name__ == "__main__":
    out = {}
    for comp in ["sigma_n", "tau"]:
        print(f"=== component: {comp} ===")
        out[comp] = run(comp)
    (HERE / "results_coso_fm.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("-> results_coso_fm.json")
