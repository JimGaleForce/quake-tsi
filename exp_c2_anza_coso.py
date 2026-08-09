"""EXP-C2 (exploratory SUPPLEMENT, outside the frozen protocol, clearly labeled as such):
phase-stationarity diagnostic for the two physically interesting regions - the EXP-A
train-selected Anza bin and Weifan's Coso Fig 4c box (which does not align to the 0.4-deg
protocol grid, hence exact bounds). Question: does the preferred phase phi hold steady across
sliding windows (a stationary property a static map could exploit), or wander (episodic
susceptibility -> argues for a tracking instrument)? Descriptive only.
"""
import json
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).parent
spec = importlib.util.spec_from_file_location("fmtest", HERE / "coso_fm_test.py")
fmtest = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fmtest)
cc = fmtest.cc

REGIONS = {
    "Anza (EXP-A selected)": dict(lat=(33.5, 33.9), lon=(-116.4, -116.0)),
    "Coso Fig4c (Weifan)": dict(lat=(36.2, 36.6), lon=(-118.0, -117.6)),
}
WIN, STEP, MIN_N, NB = 5, 1, 40, 36  # 5-yr windows for power; min 40 FM events

cat = cc.load_declustered("SCSN_decluster_m1.5.txt")
cat["time"] = pd.to_datetime(dict(year=cat.yr, month=cat.mo, day=cat.dy, hour=cat.hr,
                                  minute=cat.mi, second=0), utc=True) + \
              pd.to_timedelta(cat.sec.astype(float), unit="s")
fm = fmtest.load_fm()
m = cat.merge(fm[["eid", "strike", "dip", "rake"]], on="eid", how="inner").copy()
SXX, SYY, SXY, SZZ = fmtest.base_series()
nfine = len(SXX)


def fit36(hist):
    d = hist / hist.sum() * NB
    th = np.radians(np.linspace(-180, 180, NB + 1)[:-1] + 180 / NB)
    A = np.c_[np.ones_like(th), np.cos(th), np.sin(th)]
    coef, *_ = np.linalg.lstsq(A, d, rcond=None)
    return float(np.hypot(coef[1], coef[2]) / coef[0]), float(np.degrees(np.arctan2(coef[2], coef[1])))


out = {"exploratory_supplement": True, "not_in_frozen_protocol": True,
       "window_years": WIN, "min_fm_per_window": MIN_N, "regions": {}}
for name, box in REGIONS.items():
    r = m[m.lat.between(*box["lat"]) & m.lon.between(*box["lon"])].copy()
    for c in ["strike", "dip", "rake"]:
        r[f"{c}_r"] = (np.round(r[c] / fmtest.ROUND_DEG) * fmtest.ROUND_DEG).astype(int)
    t_unix = (r["time"] - pd.Timestamp(0, tz="UTC")).dt.total_seconds().to_numpy()
    idx_ev = np.clip(((t_unix - cc.T0.timestamp()) / cc.DT).astype(np.int64), 0, nfine - 1)
    ph = np.full(len(r), np.nan)
    for (st, di, ra), gidx in r.groupby(["strike_r", "dip_r", "rake_r"]).groups.items():
        cn, ct = fmtest.combo_coeffs(st, di, ra)
        series = ct[0] * SXX + ct[1] * SYY + ct[2] * SXY + ct[3] * SZZ
        phase = cc.phase_series(series)
        gpos = r.index.get_indexer(gidx)
        ph[gpos] = phase[idx_ev[gpos]]
    r["phase"] = ph
    rows = []
    for ws in range(1981, 2019 - WIN + 1):
        w = r[(r.time >= pd.Timestamp(f"{ws}-01-01", tz="UTC"))
              & (r.time < pd.Timestamp(f"{ws + WIN}-01-01", tz="UTC")) & r.phase.notna()]
        if len(w) < MIN_N:
            rows.append({"window_start": ws, "n": int(len(w)), "pm_p0": None, "phi": None})
            continue
        hist, _ = np.histogram(w.phase, bins=np.linspace(-180, 180, NB + 1))
        a, phi = fit36(hist)
        rows.append({"window_start": ws, "n": int(len(w)), "pm_p0": round(a, 3), "phi": round(phi, 1)})
    valid = [x for x in rows if x["pm_p0"] is not None]
    phis = np.radians([x["phi"] for x in valid])
    R = float(np.hypot(np.mean(np.cos(phis)), np.mean(np.sin(phis)))) if valid else None
    out["regions"][name] = {"n_fm_total": int(len(r)), "windows": rows,
                            "n_valid_windows": len(valid),
                            "phase_concentration_R": None if R is None else round(R, 3)}
    print(f"{name}: FM events {len(r)}, valid windows {len(valid)}, "
          f"phase concentration R={R if R is None else round(R, 2)} (1=stationary, 0=wandering)")
    for x in valid:
        print(f"  {x['window_start']}-{x['window_start'] + WIN}: n={x['n']} "
              f"Pm/P0={x['pm_p0']} phi={x['phi']}")

(HERE / "results_exp_c2.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print("-> results_exp_c2.json")
