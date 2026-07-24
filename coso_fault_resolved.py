"""EQ-18 step A (part 2) - fault-resolved positive control at Coso.

Projects the Lu/Xue tidal stress tensor onto representative CFM 5.3 Coso fault geometries
(their own rotation code, vectorized) and runs the phase-modulation test on sigma_n, shear
(rake direction), and Coulomb (mu grid, per Lu et al. Eq. 2) for Coso-box events.
"""
import json
import numpy as np
from pathlib import Path
import importlib.util

HERE = Path(__file__).parent
spec = importlib.util.spec_from_file_location("cc", HERE / "coso_positive_control.py")
cc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cc)

GEOMETRIES = {  # representative CFM5.3 geometries in the Coso box
    "normal_41_26_-110": (41, 26, -110),
    "normal_6_51_-110": (6, 51, -110),
    "normal_178_81_-110": (178, 81, -110),
    "ss_352_87_170": (352, 87, 170),
}
MUS = [0.0, 0.2, 0.4, 0.6]


def fault_stress(strike, dip, rake):
    """Vectorized version of their tidal_strain_rotate, from the raw strain files."""
    E, NU = 75e9, 0.25
    G = E / (2 * (1 + NU))
    LAM = E * NU / ((1 + NU) * (1 - 2 * NU))
    DATA = cc.DATA
    exx = np.loadtxt(DATA / "Tidal_N_0.txt") * 1e-9
    eyy = np.loadtxt(DATA / "Tidal_N_90.txt") * 1e-9
    exy = np.loadtxt(DATA / "Tidal_S_0.txt") * 1e-9
    evol = np.loadtxt(DATA / "Tidal_Vol.txt") * 1e-9
    ezz = evol - exx - eyy
    sxx = (LAM + 2 * G) * exx + LAM * eyy + LAM * ezz
    syy = LAM * exx + (LAM + 2 * G) * eyy + LAM * ezz
    szz = LAM * exx + LAM * eyy + (LAM + 2 * G) * ezz
    sxy = 2 * G * exy

    sr, s90, dr, d90, rr = map(np.radians, [strike, strike + 90, dip, 90 - dip, rake])
    n = np.array([np.cos(s90) * np.cos(d90), np.sin(s90) * np.cos(d90), -np.sin(d90)])
    s = np.array([np.cos(sr), np.sin(sr), 0])
    sd = np.array([np.cos(s90) * np.cos(dr), np.sin(s90) * np.cos(dr), np.sin(dr)])
    # T_i = n_j sigma_ji ; with sigma zz only on diagonal 3rd axis
    Tx = n[0] * sxx + n[1] * sxy
    Ty = n[0] * sxy + n[1] * syy
    Tz = n[2] * szz
    sigma_n = Tx * n[0] + Ty * n[1] + Tz * n[2]
    svx, svy, svz = Tx - sigma_n * n[0], Ty - sigma_n * n[1], Tz - sigma_n * n[2]
    tau = (svx * s[0] + svy * s[1]) * np.cos(rr) + (svx * sd[0] + svy * sd[1] + svz * sd[2]) * np.sin(rr)
    return sigma_n, tau


def upsample(v):
    from scipy.interpolate import CubicSpline
    n = len(v)
    return CubicSpline(np.arange(n), v)(np.arange(0, n - 1, 1.0 / cc.UPSAMPLE))


def main():
    rng = np.random.default_rng(cc.SEED + 5)
    df = cc.load_declustered("QTM_decluster_m0.1.txt")
    coso = df[(df.lat.between(*cc.COSO_BOX["lat"])) & (df.lon.between(*cc.COSO_BOX["lon"]))]
    print(f"Coso QTM declustered events: {len(coso)}")
    results = {}
    for gname, (st, di, ra) in GEOMETRIES.items():
        sn_raw, tau_raw = fault_stress(st, di, ra)
        series = {"sigma_n": upsample(sn_raw), "tau": upsample(tau_raw)}
        for mu in MUS:
            series[f"cfs_mu{mu}"] = upsample(tau_raw + mu * sn_raw)
        for comp, sarr in series.items():
            phase = cc.phase_series(sarr)
            p0h = cc.background_hist(phase, rng)
            ph = cc.event_phases(coso.t_unix.to_numpy(), phase)
            r = cc.significance(ph, phase, p0h, rng)
            key = f"{gname}:{comp}"
            results[key] = r
            flag = " <-- SIGNIFICANT" if r["significant_95"] else ""
            print(f"{key}: Pm/P0={r['Pm_over_P0']:.3f} thr={r['synth_p97.5']:.3f} "
                  f"p={r['empirical_p']:.4f} phi={r['phi_deg']:.0f}{flag}")
    (HERE / "results_coso_fault_resolved.json").write_text(json.dumps(results, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
