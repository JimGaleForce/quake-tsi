"""EXP-E (exploratory SUPPLEMENT, outside the frozen protocol; design chosen after seeing
EXP-C2's phase-migration result - stated plainly): two checks of the "episodic susceptibility,
tracking instrument" hypothesis at Anza (lat 33.5-33.9, lon -116.4 to -116.0).

E1: cross-catalog replication - late-period (2012-2017) Anza phase in the INDEPENDENT QTM
    catalog (different detection pipeline) vs SCSN's.

E2: walk-forward adaptive phase prediction (causally clean: only past data predicts future).
    For each FM-matched Anza event from 1987 on with >= 40 events in its trailing 5-yr
    window, phi_pred = trailing-window sinusoid preferred phase (from unshifted data - the
    predictions are identical under the null). Statistic S_adapt = mean cos(phi_i - phi_pred).
    Null: per-target-year circular time-shifts (each rep draws one offset per calendar year,
    uniform 2-2000 d, applied to event times before phase lookup; predictions unchanged),
    1,000 reps, one-sided p. Static contrast: whole-train (<2010) phi scored on 2010+.
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

BOX = dict(lat=(33.5, 33.9), lon=(-116.4, -116.0))
NB, MIN_TRAIL, TRAIL_YR = 36, 40, 5
SEED = 20260812
N_REP = 1000

SXX, SYY, SXY, SZZ = fmtest.base_series()
nfine = len(SXX)
SPAN_S = nfine * cc.DT


def fit_phi(phases):
    hist, _ = np.histogram(phases, bins=np.linspace(-180, 180, NB + 1))
    d = hist / hist.sum() * NB
    th = np.radians(np.linspace(-180, 180, NB + 1)[:-1] + 180 / NB)
    A = np.c_[np.ones_like(th), np.cos(th), np.sin(th)]
    coef, *_ = np.linalg.lstsq(A, d, rcond=None)
    return float(np.degrees(np.arctan2(coef[2], coef[1]))), float(np.hypot(coef[1], coef[2]) / coef[0])


def load_region(catfile):
    cat = cc.load_declustered(catfile)
    cat["time"] = pd.to_datetime(dict(year=cat.yr, month=cat.mo, day=cat.dy, hour=cat.hr,
                                      minute=cat.mi, second=0), utc=True) + \
                  pd.to_timedelta(cat.sec.astype(float), unit="s")
    fm = fmtest.load_fm()
    r = cat[cat.lat.between(*BOX["lat"]) & cat.lon.between(*BOX["lon"])]
    r = r.merge(fm[["eid", "strike", "dip", "rake"]], on="eid", how="inner").copy()
    for c in ["strike", "dip", "rake"]:
        r[f"{c}_r"] = (np.round(r[c] / fmtest.ROUND_DEG) * fmtest.ROUND_DEG).astype(int)
    r["t_unix"] = (r["time"] - pd.Timestamp(0, tz="UTC")).dt.total_seconds()
    return r.reset_index(drop=True)


def to_idx(t_seconds):
    rel = np.mod(np.asarray(t_seconds, dtype=float) - cc.T0.timestamp(), SPAN_S - 86400.0)
    return np.clip((rel / cc.DT).astype(np.int64), 0, nfine - 1)


def phases_obs_and_shifted(r, offsets_by_rep):
    """Observed tau phase per event, plus a (N_REP x n) matrix of per-year-shifted phases.
    One phase_series computation per orientation group (the expensive part), all lookups
    gathered from it."""
    n = len(r)
    yrs = r.time.dt.year.to_numpy()
    t = r.t_unix.to_numpy().astype(float)
    idx_obs = to_idx(t)
    ph_obs = np.full(n, np.nan)
    ph_shift = np.full((len(offsets_by_rep), n), np.nan) if offsets_by_rep is not None else None
    for key, gidx in r.groupby(["strike_r", "dip_r", "rake_r"]).groups.items():
        cn, ct = fmtest.combo_coeffs(*key)
        series = ct[0] * SXX + ct[1] * SYY + ct[2] * SXY + ct[3] * SZZ
        phase = cc.phase_series(series)
        gpos = r.index.get_indexer(gidx)
        ph_obs[gpos] = phase[idx_obs[gpos]]
        if offsets_by_rep is not None:
            tg, yg = t[gpos], yrs[gpos]
            for k, offs in enumerate(offsets_by_rep):
                ph_shift[k, gpos] = phase[to_idx(tg + np.array([offs[int(y)] for y in yg]))]
    return ph_obs, ph_shift


out = {"exploratory_supplement": True, "post_hoc_design": True, "region": BOX,
       "trail_years": TRAIL_YR, "min_trail_events": MIN_TRAIL, "n_null_reps": N_REP}

# ---- load, with per-rep year offsets prepared up front for SCSN ----
scsn = load_region("SCSN_decluster_m1.5.txt")
rng = np.random.default_rng(SEED)
uy = np.unique(scsn.time.dt.year.to_numpy())
offsets_by_rep = [{int(y): float(rng.uniform(2, 2000) * 86400) for y in uy} for _ in range(N_REP)]
print(f"SCSN Anza FM events: {len(scsn)}; computing phases + {N_REP} shifted variants...")
ph_obs, ph_shift = phases_obs_and_shifted(scsn, offsets_by_rep)
scsn["phase"] = ph_obs

qtm = load_region("QTM_decluster_m0.1.txt")
qtm["phase"], _ = phases_obs_and_shifted(qtm, None)
print(f"QTM Anza FM events: {len(qtm)}")

# ---- E1: cross-catalog late-period phase ----
e1 = {}
for name, df in [("SCSN", scsn), ("QTM_all", qtm), ("QTM_m1.5", qtm[qtm.mag >= 1.5])]:
    late = df[(df.time >= pd.Timestamp("2012-01-01", tz="UTC"))
              & (df.time < pd.Timestamp("2018-01-01", tz="UTC")) & df.phase.notna()]
    if len(late) >= 30:
        phi, a = fit_phi(late.phase.to_numpy())
        e1[name] = {"n": int(len(late)), "phi": round(phi, 1), "pm_p0": round(a, 3)}
    else:
        e1[name] = {"n": int(len(late)), "phi": None, "pm_p0": None}
    print(f"E1 {name} 2012-2017 Anza: {e1[name]}")
out["E1_late_period_phase"] = e1

# ---- E2: walk-forward adaptive prediction (SCSN) ----
s_order = np.argsort(scsn.t_unix.to_numpy())
t_arr = scsn.t_unix.to_numpy()[s_order]
ph_arr = ph_obs[s_order]
ph_shift_o = ph_shift[:, s_order]
yrs_o = scsn.time.dt.year.to_numpy()[s_order]
ok = ~np.isnan(ph_arr)

phi_pred = np.full(len(t_arr), np.nan)
for i in range(len(t_arr)):
    if yrs_o[i] < 1987 or not ok[i]:
        continue
    sel = ok & (t_arr >= t_arr[i] - TRAIL_YR * 365.25 * 86400) & (t_arr < t_arr[i])
    if sel.sum() >= MIN_TRAIL:
        phi_pred[i], _ = fit_phi(ph_arr[sel])
use = ~np.isnan(phi_pred) & ok
n_used = int(use.sum())
S_obs = float(np.mean(np.cos(np.radians(ph_arr[use] - phi_pred[use]))))
null_S = np.nanmean(np.cos(np.radians(ph_shift_o[:, use] - phi_pred[use][None, :])), axis=1)
p_adapt = float((np.sum(null_S >= S_obs) + 1) / (N_REP + 1))

train_mask = scsn.time < pd.Timestamp("2010-01-01", tz="UTC")
phi_static, _ = fit_phi(ph_obs[train_mask.to_numpy() & ~np.isnan(ph_obs)])
test_mask = (~train_mask).to_numpy() & ~np.isnan(ph_obs)
S_static = float(np.mean(np.cos(np.radians(ph_obs[test_mask] - phi_static))))

out["E2_walk_forward"] = {
    "S_adaptive": round(S_obs, 4), "n_events": n_used,
    "p_vs_per_year_shift_null": p_adapt,
    "null_mean": round(float(np.mean(null_S)), 4),
    "null_p95": round(float(np.percentile(null_S, 95)), 4),
    "static_train_phi": round(phi_static, 1), "S_static_on_test": round(S_static, 4),
    "n_test_static": int(test_mask.sum()),
}
print(f"E2: S_adapt={S_obs:.4f} over n={n_used}, p={p_adapt:.4f} "
      f"(null mean {np.mean(null_S):.4f}, p95 {np.percentile(null_S, 95):.4f}); "
      f"static phi_train={phi_static:.1f} -> S_test={S_static:.4f} (n={int(test_mask.sum())})")

(HERE / "results_exp_e.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print("-> results_exp_e.json")
