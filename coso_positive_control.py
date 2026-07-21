"""EQ-18 step A — pipeline positive control: reproduce Lu et al. (2025) tidal phase modulation
at Coso using THEIR declustered catalogs, THEIR tidal strain series, and THEIR method
(phase histogram, inherent-background-bias correction, synthetic-catalog significance).

Expected from their paper: whole-SoCal catalog NOT significant; Coso (geothermal) significant.
If we reproduce both, our pipeline demonstrably detects real signals -> the EQ-18 null is credible.

Phase convention (their Fig 2): for each event, locate the surrounding trough-peak-trough cycle of
the stress series; phase = -180 deg at preceding trough, 0 at peak (max encouraging stress),
+180 at following trough, piecewise linear in time.
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.signal import find_peaks

HERE = Path(__file__).parent
DATA = HERE / "data" / "xue_lu_zenodo"
DERIVED = HERE / "data" / "xue_lu_derived"

T0 = pd.Timestamp("1981-01-01", tz="UTC")
DT_RAW = 6000.0   # s, native sampling of the Zenodo series
UPSAMPLE = 10     # cubic-spline upsampling factor -> 600 s effective sampling
DT = DT_RAW / UPSAMPLE
E_YOUNG, NU = 75e9, 0.25
G = E_YOUNG / (2 * (1 + NU))
LAM = E_YOUNG * NU / ((1 + NU) * (1 - 2 * NU))

COSO_BOX = dict(lat=(35.60, 36.25), lon=(-118.05, -117.60))
N_BINS = 36
N_SYNTH = 3000
SEED = 20260721


def load_stress():
    exx = np.loadtxt(DATA / "Tidal_N_0.txt") * 1e-9
    eyy = np.loadtxt(DATA / "Tidal_N_90.txt") * 1e-9
    exy = np.loadtxt(DATA / "Tidal_S_0.txt") * 1e-9
    evol = np.loadtxt(DATA / "Tidal_Vol.txt") * 1e-9
    ezz = evol - exx - eyy
    sxx = (LAM + 2 * G) * exx + LAM * eyy + LAM * ezz
    syy = LAM * exx + (LAM + 2 * G) * eyy + LAM * ezz
    szz = LAM * exx + LAM * eyy + (LAM + 2 * G) * ezz
    sxy = 2 * G * exy
    svol = (sxx + syy + szz) / 3.0  # mean stress (positive = dilatational, encourages failure)
    from scipy.interpolate import CubicSpline
    n = len(svol)
    fine = np.arange(0, n - 1, 1.0 / UPSAMPLE)
    out = {}
    for k, v in {"vol": svol, "sxx": sxx}.items():
        out[k] = CubicSpline(np.arange(n), v)(fine)
    return out


def phase_series(stress):
    """Assign a tidal phase to every sample index by trough-peak-trough interpolation."""
    pk, _ = find_peaks(stress)
    tr, _ = find_peaks(-stress)
    phase = np.full(len(stress), np.nan)
    ti = 0
    for i in range(len(tr) - 1):
        a, b = tr[i], tr[i + 1]
        mids = pk[(pk > a) & (pk < b)]
        if len(mids) == 0:
            continue
        m = mids[np.argmax(stress[mids])]
        phase[a:m] = -180.0 + 180.0 * (np.arange(a, m) - a) / (m - a)
        phase[m:b] = 180.0 * (np.arange(m, b) - m) / (b - m)
    return phase


def event_phases(t_unix, phase):
    idx = np.clip(((t_unix - T0.timestamp()) / DT).astype(int), 0, len(phase) - 1)
    ph = phase[idx]
    return ph[~np.isnan(ph)]


def pm_over_p0(phases, p0_hist):
    """Least-squares fit P(theta) = P0 + Pm cos(theta - phi) on the bias-corrected histogram."""
    bins = np.linspace(-180, 180, N_BINS + 1)
    h, _ = np.histogram(phases, bins=bins)
    h = h / h.sum() * N_BINS  # density normalized so uniform == 1
    dens = h / np.where(p0_hist > 0, p0_hist, np.nan)  # correct inherent background bias
    ok = ~np.isnan(dens)
    th = np.radians((bins[:-1] + bins[1:]) / 2)[ok]
    A = np.c_[np.ones_like(th), np.cos(th), np.sin(th)]
    coef, *_ = np.linalg.lstsq(A, dens[ok], rcond=None)
    p0f, a, b = coef
    pm = np.hypot(a, b)
    phi = np.degrees(np.arctan2(b, a))
    return pm / p0f, phi


def background_hist(phase, rng, n=200000):
    """P0(theta): phase distribution of uniform random times over the loading period."""
    valid = np.flatnonzero(~np.isnan(phase))
    idx = rng.choice(valid, size=n)
    bins = np.linspace(-180, 180, N_BINS + 1)
    h, _ = np.histogram(phase[idx], bins=bins)
    return h / h.sum() * N_BINS


def significance(phases, phase, p0_hist, rng):
    obs, phi = pm_over_p0(phases, p0_hist)
    valid = np.flatnonzero(~np.isnan(phase))
    synth = np.empty(N_SYNTH)
    for i in range(N_SYNTH):
        idx = rng.choice(valid, size=len(phases))
        synth[i], _ = pm_over_p0(phase[idx], p0_hist)
    p975 = np.percentile(synth, 97.5)
    p995 = np.percentile(synth, 99.5)
    return {"n": int(len(phases)), "Pm_over_P0": float(obs), "phi_deg": float(phi),
            "synth_p97.5": float(p975), "synth_p99.5": float(p995),
            "significant_95": bool(obs > p975), "significant_99": bool(obs > p995),
            "empirical_p": float((np.sum(synth >= obs) + 1) / (N_SYNTH + 1))}


def load_declustered(fname):
    cols = ["yr", "mo", "dy", "hr", "mi", "sec", "lat", "lon", "depth", "mag", "eid"]
    df = pd.read_csv(DATA / fname, sep=r"\s+", header=None, names=cols)
    ts = pd.to_datetime(dict(year=df.yr, month=df.mo, day=df.dy, hour=df.hr, minute=df.mi,
                             second=0), utc=True) + pd.to_timedelta(df.sec.astype(float), unit="s")
    df["t_unix"] = (ts - pd.Timestamp(0, tz="UTC")).dt.total_seconds().to_numpy()
    return df


def main():
    rng = np.random.default_rng(SEED)
    stress = load_stress()
    results = {}
    for comp in ["vol", "sxx"]:
        print(f"--- stress component: {comp} ---")
        phase = phase_series(stress[comp])
        p0h = background_hist(phase, rng)
        for cat in ["QTM_decluster_m0.1.txt", "SCSN_decluster_m1.5.txt"]:
            df = load_declustered(cat)
            whole = event_phases(df.t_unix.to_numpy(), phase)
            coso = df[(df.lat.between(*COSO_BOX["lat"])) & (df.lon.between(*COSO_BOX["lon"]))]
            coso_ph = event_phases(coso.t_unix.to_numpy(), phase)
            r_whole = significance(whole, phase, p0h, rng)
            r_coso = significance(coso_ph, phase, p0h, rng)
            key = f"{cat.split('_')[0]}_{comp}"
            results[key] = {"whole_catalog": r_whole, "coso": r_coso}
            print(f"{key}: whole n={r_whole['n']} Pm/P0={r_whole['Pm_over_P0']:.3f} "
                  f"(97.5%={r_whole['synth_p97.5']:.3f}, sig={r_whole['significant_95']}) | "
                  f"coso n={r_coso['n']} Pm/P0={r_coso['Pm_over_P0']:.3f} "
                  f"(97.5%={r_coso['synth_p97.5']:.3f}, sig={r_coso['significant_95']}, "
                  f"p={r_coso['empirical_p']:.4f}, phi={r_coso['phi_deg']:.0f} deg)")
    (HERE / "results_coso_control.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print("-> results_coso_control.json")


if __name__ == "__main__":
    main()
