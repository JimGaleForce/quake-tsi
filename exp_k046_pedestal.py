"""K-046 -- THE ACF FLOOR: is K-009's surviving excess a static per-cell map error?

FROZEN SPEC = Kepler's K-046 in HYPOTHESIS_LEDGER.md ("PROPOSED (Kepler) -- Round 2",
entry "K-046 -- THE ACF FLOOR").  His claim and his frozen success rule:

  Claim. The pooled residual ACF in results_k009.json is not a decaying correlation
  function. It is A*exp(-k/tau) + C with A = 0.0654, tau = 7.16 weeks, C = 0.0382,
  and dBIC = +42.2 against the single-exponential form.  "A wrong map does not decay."

  Statistic and success rule. Primary: the demeaned lag-1 excess rho'(1) - null p97.5.
  "If rho'(1) excess < 0.01 while the raw excess is 0.0935, the surviving content of
  K-009 is 'our background map is wrong by 3.8% of residual variance in a spatially
  organised way', not 'there is weather', and the assimilation thread closes on a
  static finding."  Secondary: dBIC preference for the +constant form over pure
  exponential, real vs sims.

  Expected effect if he is right: rho'(1) excess falls to ~0.00-0.015; the +constant
  model wins by dBIC > 20 on the real field and loses on every sim; m_i correlates
  with distance-to-fault in the Bray direction at |rho| > 0.2.
  If he is wrong: the demeaned ACF keeps a decaying excess with a measurable tau and
  K-009 becomes much stronger.

Arms implemented here (his test list):
  1. Decisive arm  -- pooled ACF on the cell-demeaned field, real + all 20 sims,
                      identical code path; plus the El-Mayor-excluded field.
  2. Decomposition -- Var(r) = Var_between-cell(time-means) + Var_within-cell;
                      confirm the between-cell fraction equals C.
  3. Shape scoring -- BIC for (i) pure exponential, (ii) exponential + constant,
                      (iii) two exponentials, on real and on each sim, raw and excl-2010.
  4. Localise      -- per-cell time-mean residual m_i regressed on distance to CFM5.3
                      fault traces (Bray et al. 2014 on-fault-under / off-fault-over).
  5. tau_slow      -- two-exponential fit on the sequence-excluded field with a
                      cell-bootstrap CI; only that number is admissible against the
                      t_a band (13.8-55.3 wk).
  Null             -- the identical pipeline on the 20 ETAS-sim catalogues, demeaned
                      identically, PLUS a map-error injection positive control: perturb
                      mu(x) by a frozen static multiplicative field of known amplitude,
                      simulate, and confirm the pipeline recovers C = the injected
                      between-cell variance fraction.

REPRODUCIBILITY: the grid, background fields and spatial kernel fit are deterministic
and identical to exp_k009_residual_whiteness.py, and the 20 null catalogues are
regenerated from the same seed with the same call sequence, so they are the SAME 20
catalogues K-009 scored.

Run: python -u exp_k046_pedestal.py
"""

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, minimize
from scipy.stats import spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import exp_k009_residual_whiteness as K9
from exp_k009_residual_whiteness import (load_catalog, to_km, Grid, adaptive_background,
                                         SpatialOps, expected_counts, residual_field,
                                         pooled_temporal_acf, simulate_st_etas, leading_eof,
                                         CATALOG, LAT_MIN, LAT_MAX, LON_MIN, LON_MAX,
                                         T_TEST0, T_TEST1, DT_DAYS, MAX_LAG_WEEKS,
                                         RNG_SEED, N_SIMS_TARGET, SIM_EVENT_CAP_FACTOR)

HERE = Path(__file__).resolve().parent
OUT = HERE / "results_k046.json"
FIG = HERE / "maps" / "k046_pedestal.png"
EXPH = HERE / "results_exp_h.json"
K009 = HERE / "results_k009.json"

KEPLER_A, KEPLER_TAU, KEPLER_C = 0.0654, 7.16, 0.0382
N_INJECT = 3            # map-error injection positive controls
INJECT_SD = [0.30, 0.60]
N_BOOT = 200
RUNTIME_BUDGET_MIN = 90.0


# ------------------------------------------------------------------ shape models
def m_exp(k, A, tau):
    return A * np.exp(-k / tau)


def m_exp_c(k, A, tau, C):
    return A * np.exp(-k / tau) + C


def m_exp2(k, A1, t1, A2, t2):
    return A1 * np.exp(-k / t1) + A2 * np.exp(-k / t2)


def _fit(model, k, y, p0s, bounds):
    best = None
    for p0 in p0s:
        try:
            popt, _ = curve_fit(model, k, y, p0=p0, bounds=bounds, maxfev=40000)
        except Exception:
            continue
        rss = float(np.sum((y - model(k, *popt)) ** 2))
        if best is None or rss < best[1]:
            best = (popt, rss)
    return best


def bic(rss, n, p):
    return float(n * np.log(max(rss, 1e-300) / n) + p * np.log(n))


def shape_scoring(acf, max_lag=MAX_LAG_WEEKS):
    """BIC for pure-exponential vs exponential+constant vs two-exponential."""
    k = np.arange(1, max_lag + 1, dtype=float)
    y = np.asarray(acf[1:max_lag + 1], dtype=float)
    n = len(k)
    out = {}
    f1 = _fit(m_exp, k, y, [(0.1, 5.0), (0.05, 20.0), (0.02, 2.0)],
              ([-1.0, 0.3], [2.0, 500.0]))
    f2 = _fit(m_exp_c, k, y, [(0.07, 7.0, 0.04), (0.1, 3.0, 0.02), (0.05, 20.0, 0.01)],
              ([-1.0, 0.3, -1.0], [2.0, 500.0, 1.0]))
    f3 = _fit(m_exp2, k, y, [(0.07, 3.0, 0.04, 60.0), (0.1, 1.0, 0.03, 30.0),
                             (0.05, 8.0, 0.02, 200.0)],
              ([-1.0, 0.3, -1.0, 0.3], [2.0, 200.0, 2.0, 3000.0]))
    for name, f, p in [("exp", f1, 2), ("exp_plus_const", f2, 3), ("two_exp", f3, 4)]:
        if f is None:
            out[name] = {"FAILED": True}
            continue
        popt, rss = f
        out[name] = {"params": [float(v) for v in popt], "rss": rss, "bic": bic(rss, n, p), "k": p}
    have = {kk: v for kk, v in out.items() if "bic" in v}
    if have:
        bmin = min(v["bic"] for v in have.values())
        for kk, v in have.items():
            v["dBIC_vs_best"] = v["bic"] - bmin
        out["best_model"] = min(have, key=lambda kk: have[kk]["bic"])
        if "exp" in have and "exp_plus_const" in have:
            # positive => +constant preferred
            out["dBIC_exp_minus_expconst"] = have["exp"]["bic"] - have["exp_plus_const"]["bic"]
    return out


def variance_split(R):
    """Var(r) = Var_between-cell(time-means) + Var_within-cell (about each cell's mean)."""
    m = R.mean(axis=1)
    gm = R.mean()
    v_tot = float(np.mean((R - gm) ** 2))
    v_between = float(np.mean((m - gm) ** 2))
    v_within = float(np.mean((R - m[:, None]) ** 2))
    return {"var_total": v_tot, "var_between_cell": v_between, "var_within_cell": v_within,
            "between_fraction": float(v_between / v_tot) if v_tot > 0 else float("nan"),
            "identity_check_between_plus_within_over_total":
                float((v_between + v_within) / v_tot) if v_tot > 0 else float("nan")}


def demean_cells(R):
    return R - R.mean(axis=1, keepdims=True)


def both_stats(R, sops, max_lag=MAX_LAG_WEEKS):
    """Raw and cell-demeaned pooled ACF + Moran's I through the identical code path."""
    Rd = demean_cells(R)
    a_raw = pooled_temporal_acf(R, max_lag)
    a_dm = pooled_temporal_acf(Rd, max_lag)
    mI_raw, _ = sops.moran(R)
    mI_dm, _ = sops.moran(Rd)
    cor_raw = sops.correlogram(R)
    cor_dm = sops.correlogram(Rd)
    return {"acf_raw": a_raw, "acf_demeaned": a_dm,
            "acf1_raw": float(a_raw[1]), "acf1_demeaned": float(a_dm[1]),
            "moran_raw": mI_raw, "moran_demeaned": mI_dm,
            "correlogram_raw": cor_raw, "correlogram_demeaned": cor_dm,
            "variance_split": variance_split(R)}


# ------------------------------------------------------------------ main
def main():
    t_start = time.time()
    res = {"experiment": "K-046",
           "spec": ("HYPOTHESIS_LEDGER.md :: 'PROPOSED (Kepler) -- Round 2' :: "
                    "'K-046 -- THE ACF FLOOR'"),
           "relation_to_K009": ("post-hoc re-analysis of the K-009 residual field with one "
                                "additional control (per-cell temporal demeaning) that the K-009 "
                                "spec, code and adjudications all omitted"),
           "state_class": "first-run (post-hoc re-analysis; NOT pre-registered)",
           "run_utc": pd.Timestamp.now("UTC").isoformat()}

    exph = json.loads(EXPH.read_text())
    fp = exph["train_fit"]["frozen_params"]
    b_val = exph["train_fit"]["b_value_train_aki"]
    M0 = fp["M0"]

    # ---------- rebuild the K-009 configuration (deterministic) ----------
    df, order = load_catalog(CATALOG)
    box = df.lat.between(LAT_MIN, LAT_MAX) & df.lon.between(LON_MIN, LON_MAX)
    dall = df[box].sort_values("t").reset_index(drop=True)
    cat = dall[dall.mag >= M0 - 1e-9].reset_index(drop=True)
    t0w = (T_TEST0 - pd.Timestamp(0, tz="UTC")).total_seconds() / 86400.0
    t1w = (T_TEST1 - pd.Timestamp(0, tz="UTC")).total_seconds() / 86400.0
    week_edges = np.arange(t0w, t1w + 1e-9, DT_DAYS)
    n_weeks = len(week_edges) - 1
    train = cat[cat.t < t0w]
    test = cat[(cat.t >= week_edges[0]) & (cat.t < week_edges[-1])]
    grid = Grid(train.lat.to_numpy(), train.lon.to_numpy())
    p_k4 = adaptive_background(train.lat.to_numpy(), train.lon.to_numpy(),
                               grid.clat, grid.clon, k=4)
    sops = SpatialOps(grid)
    k9 = json.loads(K009.read_text())
    sk = k9["spatial_kernel_fit"]
    PARS = dict(mu=fp["mu"], K=fp["K"], alpha=fp["alpha"], c=fp["c"], p=fp["p"], M0=M0,
                d=sk["d_km"], gamma=sk["gamma"], q=sk["q"])
    print(f"[setup] {grid.n_cells} cells x {n_weeks} weeks; spatial kernel reused from "
          f"results_k009.json (d={PARS['d']:.3f} gamma={PARS['gamma']:.3f} q={PARS['q']:.3f})")

    src = cat[cat.t < week_edges[-1]]
    obs = grid.bin_events(test.lat.to_numpy(), test.lon.to_numpy(), test.t.to_numpy(), week_edges)
    E_real, F_real = expected_counts(grid, src.t.to_numpy(), src.mag.to_numpy(),
                                     src.lat.to_numpy(), src.lon.to_numpy(),
                                     week_edges, PARS, p_k4)
    R_real = residual_field(obs, E_real)
    # sanity: must reproduce the K-009 numbers exactly
    a1_check = float(pooled_temporal_acf(R_real, 1)[1])
    print(f"[check] reproduced raw acf1 = {a1_check:.6f}  (results_k009.json: "
          f"{k9['real']['adaptive_k4']['acf1']:.6f})")
    res["reproduction_check"] = {"acf1_here": a1_check,
                                 "acf1_in_results_k009": k9["real"]["adaptive_k4"]["acf1"],
                                 "match": bool(abs(a1_check - k9["real"]["adaptive_k4"]["acf1"]) < 1e-9)}

    # ---------- ARM 1 + 2: real field, raw vs demeaned ----------
    st_real = both_stats(R_real, sops)
    vs = st_real["variance_split"]
    print(f"[real] acf1 raw={st_real['acf1_raw']:+.4f}  demeaned={st_real['acf1_demeaned']:+.4f}")
    print(f"[real] Moran raw={st_real['moran_raw']:+.5f}  demeaned={st_real['moran_demeaned']:+.5f}")
    print(f"[real] variance split: between-cell fraction = {vs['between_fraction']:.4f} "
          f"(Kepler predicts = C = {KEPLER_C})")

    # El-Mayor-excluded arm
    R_ex = R_real[:, 52:]
    st_ex = both_stats(R_ex, sops)
    print(f"[excl-2010] acf1 raw={st_ex['acf1_raw']:+.4f}  demeaned={st_ex['acf1_demeaned']:+.4f}; "
          f"between-cell fraction={st_ex['variance_split']['between_fraction']:.4f}")

    # ---------- ARM 3: shape scoring on the real curves ----------
    shp_real_raw = shape_scoring(st_real["acf_raw"])
    shp_real_dm = shape_scoring(st_real["acf_demeaned"])
    shp_ex_raw = shape_scoring(st_ex["acf_raw"])
    shp_ex_dm = shape_scoring(st_ex["acf_demeaned"])
    print(f"[shape real raw] best={shp_real_raw['best_model']}  "
          f"dBIC(exp - exp+const)={shp_real_raw.get('dBIC_exp_minus_expconst', float('nan')):+.1f}  "
          f"exp+const params A,tau,C = "
          f"{[round(v,4) for v in shp_real_raw['exp_plus_const']['params']]}")
    print(f"[shape real demeaned] best={shp_real_dm['best_model']}  "
          f"dBIC(exp - exp+const)={shp_real_dm.get('dBIC_exp_minus_expconst', float('nan')):+.1f}")

    # Kepler's stated decomposition vs our raw curve
    kk = np.arange(1, MAX_LAG_WEEKS + 1, dtype=float)
    yk = np.asarray(st_real["acf_raw"][1:MAX_LAG_WEEKS + 1])
    y_kep = m_exp_c(kk, KEPLER_A, KEPLER_TAU, KEPLER_C)
    rss_kep = float(np.sum((yk - y_kep) ** 2))
    ours = shp_real_raw["exp_plus_const"]["params"]
    res["kepler_decomposition_check"] = {
        "kepler_stated": {"A": KEPLER_A, "tau_weeks": KEPLER_TAU, "C": KEPLER_C,
                          "dBIC_claimed": 42.2},
        "our_refit": {"A": ours[0], "tau_weeks": ours[1], "C": ours[2]},
        "rss_kepler_params": rss_kep, "rss_our_refit": shp_real_raw["exp_plus_const"]["rss"],
        "max_abs_deviation_kepler_curve": float(np.max(np.abs(yk - y_kep))),
        "rms_deviation_kepler_curve": float(np.sqrt(np.mean((yk - y_kep) ** 2))),
        "dBIC_ours_exp_minus_expconst": shp_real_raw.get("dBIC_exp_minus_expconst"),
        "flat_tail_mean_lags20_52": float(np.mean(st_real["acf_raw"][20:53])),
        "flat_tail_sd_lags20_52": float(np.std(st_real["acf_raw"][20:53])),
        "flat_tail_slope_per_lag": float(np.polyfit(np.arange(20, 53),
                                                    st_real["acf_raw"][20:53], 1)[0])}
    print(f"[kepler] his A={KEPLER_A} tau={KEPLER_TAU} C={KEPLER_C}; our refit "
          f"A={ours[0]:.4f} tau={ours[1]:.2f} C={ours[2]:.4f}; "
          f"rms dev of his curve from ours = {res['kepler_decomposition_check']['rms_deviation_kepler_curve']:.5f}")
    print(f"[kepler] flat tail lags 20-52: mean={res['kepler_decomposition_check']['flat_tail_mean_lags20_52']:.4f} "
          f"+-{res['kepler_decomposition_check']['flat_tail_sd_lags20_52']:.4f}, "
          f"slope={res['kepler_decomposition_check']['flat_tail_slope_per_lag']:.2e}/lag")

    # ---------- NULL: the same 20 ETAS-sim catalogues, demeaned identically ----------
    hist = cat[cat.t < week_edges[0]]
    h_t, h_m = hist.t.to_numpy(), hist.mag.to_numpy()
    h_la, h_lo = hist.lat.to_numpy(), hist.lon.to_numpy()
    cap = int(SIM_EVENT_CAP_FACTOR * len(test))
    rng = np.random.default_rng(RNG_SEED)
    sims = []
    print(f"[null] regenerating the same {N_SIMS_TARGET} ETAS-sim catalogues (seed {RNG_SEED}) ...")
    for s in range(N_SIMS_TARGET):
        ts_ = time.time()
        st_, sm_, sla_, slo_, trunc = simulate_st_etas(
            rng, h_t, h_m, h_la, h_lo, PARS, p_k4, grid, week_edges[0], week_edges[-1], b_val, cap)
        inw = (st_ >= week_edges[0]) & (st_ < week_edges[-1])
        inbox = ((sla_ >= LAT_MIN) & (sla_ <= LAT_MAX) & (slo_ >= LON_MIN) & (slo_ <= LON_MAX))
        obs_s = grid.bin_events(sla_[inw & inbox], slo_[inw & inbox], st_[inw & inbox], week_edges)
        ss_t = np.concatenate([h_t, st_[inw]]); ss_m = np.concatenate([h_m, sm_[inw]])
        ss_la = np.concatenate([h_la, sla_[inw]]); ss_lo = np.concatenate([h_lo, slo_[inw]])
        o = np.argsort(ss_t)
        E_s, F_s = expected_counts(grid, ss_t[o], ss_m[o], ss_la[o], ss_lo[o],
                                   week_edges, PARS, p_k4)
        R_s = residual_field(obs_s, E_s)
        stt = both_stats(R_s, sops)
        stt["shape_raw"] = shape_scoring(stt["acf_raw"])
        stt["n_events_in_box"] = int((inw & inbox).sum())
        sims.append(stt)
        del F_s
        print(f"[null] sim {s+1}/{N_SIMS_TARGET}: acf1 raw={stt['acf1_raw']:+.4f} "
              f"demeaned={stt['acf1_demeaned']:+.4f}  Moran raw={stt['moran_raw']:+.5f} "
              f"dm={stt['moran_demeaned']:+.5f}  betweenfrac="
              f"{stt['variance_split']['between_fraction']:.4f} ({time.time()-ts_:.0f}s)")

    def env(vals):
        v = np.asarray([x for x in vals if np.isfinite(x)], dtype=float)
        return {"n": int(len(v)), "median": float(np.median(v)),
                "p2_5": float(np.percentile(v, 2.5)), "p97_5": float(np.percentile(v, 97.5)),
                "min": float(v.min()), "max": float(v.max())}

    NULL = {
        "acf1_raw": env([s["acf1_raw"] for s in sims]),
        "acf1_demeaned": env([s["acf1_demeaned"] for s in sims]),
        "moran_raw": env([s["moran_raw"] for s in sims]),
        "moran_demeaned": env([s["moran_demeaned"] for s in sims]),
        "between_fraction": env([s["variance_split"]["between_fraction"] for s in sims]),
        "dBIC_exp_minus_expconst": env([s["shape_raw"].get("dBIC_exp_minus_expconst", np.nan)
                                        for s in sims]),
    }
    n_sims_expconst_wins = int(sum(s["shape_raw"].get("best_model") == "exp_plus_const"
                                   for s in sims))

    # ---------- primary statistic ----------
    exc_raw = st_real["acf1_raw"] - NULL["acf1_raw"]["p97_5"]
    exc_dm = st_real["acf1_demeaned"] - NULL["acf1_demeaned"]["p97_5"]
    exc_dm_ex = st_ex["acf1_demeaned"] - NULL["acf1_demeaned"]["p97_5"]
    print(f"\n[PRIMARY] raw lag-1 excess = {exc_raw:+.4f}   "
          f"demeaned lag-1 excess = {exc_dm:+.4f}  (Kepler predicts < 0.015; "
          f"his 'closes the thread' threshold is < 0.01)")

    # ---------- ARM 5: two-exponential tau_slow on the sequence-excluded field ----------
    def tau_slow_from(acf):
        sh = shape_scoring(acf)
        if "two_exp" not in sh or "params" not in sh["two_exp"]:
            return None
        A1, t1, A2, t2 = sh["two_exp"]["params"]
        return float(max(t1, t2)), float(A2 if t2 >= t1 else A1)

    boot = []
    rngb = np.random.default_rng(12345)
    n_c = R_ex.shape[0]
    for _ in range(N_BOOT):
        idx = rngb.integers(0, n_c, n_c)
        a = pooled_temporal_acf(R_ex[idx], MAX_LAG_WEEKS)
        r = tau_slow_from(a)
        if r is not None and np.isfinite(r[0]):
            boot.append(r[0])
    ts_pt = tau_slow_from(st_ex["acf_raw"])
    tau_slow = {"point_weeks": ts_pt[0] if ts_pt else None,
                "slow_amplitude": ts_pt[1] if ts_pt else None,
                "boot_n": len(boot),
                "ci95_weeks": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]
                if len(boot) > 20 else None,
                "t_a_band_weeks": [13.8, 55.3],
                "note": ("Kepler: only this number is admissible against the t_a band; the 2.42-week "
                         "integral statistic is an integral over a curve with a pedestal and is not "
                         "an estimate of any timescale")}
    if tau_slow["ci95_weeks"]:
        lo, hi = tau_slow["ci95_weeks"]
        tau_slow["overlaps_t_a_band"] = bool(hi >= 13.8 and lo <= 55.3)
        print(f"[tau_slow] excl-2010 two-exp slow timescale = {ts_pt[0]:.1f} wk "
              f"(95% CI {lo:.1f}-{hi:.1f}); t_a band 13.8-55.3 -> "
              f"overlaps={tau_slow['overlaps_t_a_band']}")

    # ---------- ARM 4: localise the map error against CFM5.3 fault traces ----------
    m_i = R_real.mean(axis=1)
    tr = HERE / "data" / "xue_lu_zenodo" / "CFM5.3_traces.lonLat"
    fault = None
    if tr.exists():
        raw = np.genfromtxt(tr, comments="#", invalid_raise=False)
        raw = raw[np.isfinite(raw[:, 0]) & np.isfinite(raw[:, 1])]
        fl, fla = raw[:, 0], raw[:, 1]
        keep = ((fla > LAT_MIN - 0.5) & (fla < LAT_MAX + 0.5) &
                (fl > LON_MIN - 0.5) & (fl < LON_MAX + 0.5))
        fx, fy = to_km(fla[keep], fl[keep])
        from scipy.spatial import cKDTree
        dist, _ = cKDTree(np.column_stack([fx, fy])).query(
            np.column_stack([grid.cx, grid.cy]), k=1)
        rho_s, p_s = spearmanr(dist, m_i)
        onf = dist <= 5.0
        offf = dist > 15.0
        fault = {"n_trace_points": int(keep.sum()),
                 "spearman_m_i_vs_distance_to_fault": float(rho_s), "p_value": float(p_s),
                 "mean_m_i_on_fault_le5km": float(m_i[onf].mean()) if onf.any() else None,
                 "n_on_fault": int(onf.sum()),
                 "mean_m_i_off_fault_gt15km": float(m_i[offf].mean()) if offf.any() else None,
                 "n_off_fault": int(offf.sum()),
                 "bray_direction_expected": ("on-fault UNDER-predicted by a smoothed map "
                                             "(m_i > 0 near faults) and off-fault OVER-predicted "
                                             "(m_i < 0 away), i.e. NEGATIVE spearman with distance"),
                 "bray_signature_present": bool(rho_s < -0.2 and p_s < 0.05),
                 "kepler_threshold_abs_rho_gt_0.2": bool(abs(rho_s) > 0.2)}
        print(f"[fault] Spearman(m_i, distance to CFM5.3) = {rho_s:+.3f} (p={p_s:.2e}); "
              f"mean m_i on-fault(<=5km, n={int(onf.sum())}) = {m_i[onf].mean():+.4f}, "
              f"off-fault(>15km, n={int(offf.sum())}) = {m_i[offf].mean():+.4f}")
    res["map_error_localisation"] = fault

    # ---------- map-error injection POSITIVE CONTROL ----------
    print("[inject] map-error injection control: perturb mu(x) by a frozen static field, "
          "simulate, score with the UNperturbed map ...")
    inj = []
    d2 = ((grid.cx[:, None] - grid.cx[None, :]) ** 2 + (grid.cy[:, None] - grid.cy[None, :]) ** 2)
    Cs = np.exp(-np.sqrt(d2) / 35.0)
    ev, V = np.linalg.eigh(Cs)
    Lh = V * np.sqrt(np.clip(ev, 0, None))[None, :]
    rngi = np.random.default_rng(777)
    for sd in INJECT_SD:
        for rep in range(N_INJECT // len(INJECT_SD) + 1):
            if (time.time() - t_start) / 60.0 > RUNTIME_BUDGET_MIN * 0.88:
                print("[inject] RUNTIME GUARD: stopping injection controls")
                break
            g = Lh @ rngi.normal(size=grid.n_cells)
            g = g / max(g.std(), 1e-9) * sd
            pert = p_k4 * np.exp(g - 0.5 * sd ** 2)
            pert = pert / pert.sum()
            st_, sm_, sla_, slo_, _ = simulate_st_etas(
                rngi, h_t, h_m, h_la, h_lo, PARS, pert, grid,
                week_edges[0], week_edges[-1], b_val, cap)
            inw = (st_ >= week_edges[0]) & (st_ < week_edges[-1])
            inbox = ((sla_ >= LAT_MIN) & (sla_ <= LAT_MAX) &
                     (slo_ >= LON_MIN) & (slo_ <= LON_MAX))
            obs_s = grid.bin_events(sla_[inw & inbox], slo_[inw & inbox], st_[inw & inbox], week_edges)
            ss_t = np.concatenate([h_t, st_[inw]]); ss_m = np.concatenate([h_m, sm_[inw]])
            ss_la = np.concatenate([h_la, sla_[inw]]); ss_lo = np.concatenate([h_lo, slo_[inw]])
            o = np.argsort(ss_t)
            # scored with the ORIGINAL (unperturbed) map -- that is the map error
            E_s, _ = expected_counts(grid, ss_t[o], ss_m[o], ss_la[o], ss_lo[o],
                                     week_edges, PARS, p_k4)
            R_s = residual_field(obs_s, E_s)
            sti = both_stats(R_s, sops)
            shi = shape_scoring(sti["acf_raw"])
            inj.append({"inject_log_sd": sd,
                        "acf1_raw": sti["acf1_raw"], "acf1_demeaned": sti["acf1_demeaned"],
                        "between_fraction": sti["variance_split"]["between_fraction"],
                        "fitted_C": (shi["exp_plus_const"]["params"][2]
                                     if "params" in shi.get("exp_plus_const", {}) else None),
                        "best_model": shi.get("best_model"),
                        "moran_raw": sti["moran_raw"], "moran_demeaned": sti["moran_demeaned"]})
            print(f"[inject] sd={sd}: acf1 raw={sti['acf1_raw']:+.4f} dm={sti['acf1_demeaned']:+.4f}  "
                  f"between-frac={sti['variance_split']['between_fraction']:.4f}  "
                  f"fitted C={inj[-1]['fitted_C']}  best={shi.get('best_model')}")
            if len(inj) >= N_INJECT * len(INJECT_SD):
                break
    res["map_error_injection_control"] = {
        "design": ("mu(x) perturbed by a frozen static log-normal field (35 km correlation "
                   "length), catalogues simulated from the perturbed map, residuals scored "
                   "against the UNPERTURBED map -- a pure static map error with no dynamics"),
        "runs": inj,
        "PASSES": bool(inj and all(x["acf1_raw"] - x["acf1_demeaned"] > 0.01 for x in inj)),
        "reads": ("a pure static map error must produce a raw lag-1 excess that VANISHES on "
                  "demeaning, and a fitted C close to its between-cell variance fraction")}

    # ---------- VERDICT on Kepler's frozen rule ----------
    kepler_right = bool(exc_dm < 0.015)
    kepler_strong = bool(exc_dm < 0.01)
    expconst_wins_real = bool(shp_real_raw.get("best_model") == "exp_plus_const"
                              and shp_real_raw.get("dBIC_exp_minus_expconst", 0) > 20)
    if kepler_strong:
        verdict = ("K-046 CONFIRMED -- K-009's surviving excess is a STATIC PER-CELL MAP ERROR. "
                   "Reassign K-009 to 'background-map error measured'.")
    elif kepler_right:
        verdict = ("K-046 LARGELY CONFIRMED -- the demeaned excess collapses below 0.015 but not "
                   "below his 0.01 'closes the thread' threshold.")
    else:
        verdict = ("K-046 REFUTED -- a decaying excess survives per-cell demeaning; the pedestal "
                   "does not account for K-009, which is HARDENED.")
    res["verdict"] = {
        "quoted_rule": ("Primary: the demeaned lag-1 excess rho'(1) - null p97.5. If rho'(1) "
                        "excess < 0.01 while the raw excess is 0.0935, the surviving content of "
                        "K-009 is 'our background map is wrong by 3.8% of residual variance in a "
                        "spatially organised way', not 'there is weather'. Secondary: dBIC "
                        "preference for the +constant form over pure exponential, real vs sims."),
        "raw_lag1_excess": exc_raw, "demeaned_lag1_excess": exc_dm,
        "demeaned_lag1_excess_excl_2010": exc_dm_ex,
        "kepler_threshold_0.015_met": kepler_right,
        "kepler_threshold_0.010_met": kepler_strong,
        "expconst_wins_on_real_by_dBIC_gt20": expconst_wins_real,
        "n_sims_where_expconst_wins": n_sims_expconst_wins, "n_sims": len(sims),
        "VERDICT": verdict}

    res["real"] = {
        "raw": {"acf": list(map(float, st_real["acf_raw"])), "acf1": st_real["acf1_raw"],
                "moran": st_real["moran_raw"],
                "correlogram": [None if not np.isfinite(x) else float(x)
                                for x in st_real["correlogram_raw"]]},
        "demeaned": {"acf": list(map(float, st_real["acf_demeaned"])),
                     "acf1": st_real["acf1_demeaned"], "moran": st_real["moran_demeaned"],
                     "correlogram": [None if not np.isfinite(x) else float(x)
                                     for x in st_real["correlogram_demeaned"]]},
        "variance_split": vs,
        "shape_raw": shp_real_raw, "shape_demeaned": shp_real_dm,
        "per_cell_time_mean_stats": {
            "mean": float(m_i.mean()), "sd": float(m_i.std()),
            "min": float(m_i.min()), "max": float(m_i.max()),
            "frac_cells_negative": float((m_i < 0).mean())}}
    res["excl_2010"] = {
        "raw": {"acf1": st_ex["acf1_raw"], "moran": st_ex["moran_raw"],
                "acf": list(map(float, st_ex["acf_raw"]))},
        "demeaned": {"acf1": st_ex["acf1_demeaned"], "moran": st_ex["moran_demeaned"],
                     "acf": list(map(float, st_ex["acf_demeaned"]))},
        "variance_split": st_ex["variance_split"],
        "shape_raw": shp_ex_raw, "shape_demeaned": shp_ex_dm,
        "tau_slow": tau_slow}
    res["null"] = NULL
    res["null_per_sim"] = [{"acf1_raw": s["acf1_raw"], "acf1_demeaned": s["acf1_demeaned"],
                            "moran_raw": s["moran_raw"], "moran_demeaned": s["moran_demeaned"],
                            "between_fraction": s["variance_split"]["between_fraction"],
                            "best_model": s["shape_raw"].get("best_model"),
                            "dBIC_exp_minus_expconst": s["shape_raw"].get("dBIC_exp_minus_expconst"),
                            "n_events_in_box": s["n_events_in_box"]} for s in sims]
    res["null_acf_envelopes"] = {
        "lags": list(range(MAX_LAG_WEEKS + 1)),
        "raw_p97_5": np.percentile(np.array([s["acf_raw"] for s in sims]), 97.5, axis=0).tolist(),
        "demeaned_p97_5": np.percentile(np.array([s["acf_demeaned"] for s in sims]), 97.5, axis=0).tolist(),
        "demeaned_median": np.median(np.array([s["acf_demeaned"] for s in sims]), axis=0).tolist()}
    res["runtime_minutes"] = round((time.time() - t_start) / 60.0, 2)
    res["flags"] = {
        "post_hoc_not_pre_registered": ("K-046 is a post-hoc re-analysis of an already-scored "
                                        "result. Kepler himself flags this and invokes S-12(c): it "
                                        "is proposed as a frozen rule for the K-009R 2019+ re-run, "
                                        "not as a pre-registered test of this window."),
        "same_20_sims_as_K009": True, "n_sims": len(sims), "n_sims_spec_ideal": 500,
        "rho_sta_still_unavailable": "K-031 not run; his item 4(c) not testable",
        "ross_cochran_labels_absent": "his item 4(b) not on disk; not tested",
        "world_arm_unrun": "inherited from K-009",
        "bootstrap_for_tau_slow": f"cell-resampling bootstrap, n={N_BOOT}"}

    OUT.write_text(json.dumps(res, indent=2, default=float))

    # ---------- figure ----------
    FIG.parent.mkdir(exist_ok=True)
    lags = np.arange(MAX_LAG_WEEKS + 1)
    raw_p975 = np.percentile(np.array([s["acf_raw"] for s in sims]), 97.5, axis=0)
    dm_p975 = np.percentile(np.array([s["acf_demeaned"] for s in sims]), 97.5, axis=0)
    dm_p025 = np.percentile(np.array([s["acf_demeaned"] for s in sims]), 2.5, axis=0)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    a = ax[0, 0]
    a.plot(lags[1:], st_real["acf_raw"][1:], "-", lw=1.8, color="C0", label="real, RAW ACF")
    a.plot(lags[1:], m_exp_c(lags[1:].astype(float), *ours), "--", color="k", lw=1.2,
           label=f"our fit A={ours[0]:.3f} tau={ours[1]:.1f}wk C={ours[2]:.4f}")
    a.plot(lags[1:], y_kep, ":", color="crimson", lw=1.4,
           label=f"Kepler {KEPLER_A}exp(-k/{KEPLER_TAU})+{KEPLER_C}")
    a.axhline(ours[2], color="0.5", lw=0.8)
    a.fill_between(lags[1:], 0, raw_p975[1:], color="0.85", label="ETAS-sim null p97.5 (raw)")
    a.set_xlabel("lag (weeks)"); a.set_ylabel("pooled residual ACF")
    a.set_title(f"(i) RAW ACF has a flat pedestal\n"
                f"dBIC(exp - exp+const) = {shp_real_raw.get('dBIC_exp_minus_expconst', float('nan')):+.1f} "
                f"-> {shp_real_raw['best_model']}")
    a.legend(fontsize=7)

    a = ax[0, 1]
    a.fill_between(lags[1:], dm_p025[1:], dm_p975[1:], color="0.85",
                   label=f"ETAS-sim null 95% ({len(sims)} sims)")
    a.plot(lags[1:], st_real["acf_demeaned"][1:], "-", lw=1.8, color="C3",
           label="real, CELL-DEMEANED")
    a.plot(lags[1:], st_ex["acf_demeaned"][1:], "--", lw=1.2, color="C1",
           label="real, demeaned, excl-2010")
    a.axhline(0, color="k", lw=0.5)
    a.set_xlabel("lag (weeks)"); a.set_ylabel("pooled residual ACF (cell-demeaned)")
    a.set_title(f"(ii) THE DECISIVE ARM\ndemeaned lag-1 excess = {exc_dm:+.4f} "
                f"(raw {exc_raw:+.4f}); Kepler predicts < 0.015")
    a.legend(fontsize=7)

    a = ax[1, 0]
    labels = ["raw", "demeaned"]
    rv = [st_real["acf1_raw"], st_real["acf1_demeaned"]]
    nv = [NULL["acf1_raw"]["p97_5"], NULL["acf1_demeaned"]["p97_5"]]
    xpos = np.arange(2)
    a.bar(xpos - 0.18, rv, 0.36, label="real", color="C3")
    a.bar(xpos + 0.18, nv, 0.36, label="sim-null p97.5", color="0.6")
    for i, (r_, n_) in enumerate(zip(rv, nv)):
        a.text(i, max(r_, n_) + 0.003, f"excess {r_-n_:+.4f}", ha="center", fontsize=9)
    a.set_xticks(xpos); a.set_xticklabels(labels)
    a.axhline(0, color="k", lw=0.5)
    a.set_ylabel("lag-1 pooled ACF")
    a.set_title(f"(iii) lag-1 before/after per-cell demeaning\n"
                f"between-cell variance fraction = {vs['between_fraction']:.4f} "
                f"vs fitted C = {ours[2]:.4f}")
    a.legend(fontsize=8)

    a = ax[1, 1]
    sc = a.scatter(grid.clon, grid.clat, c=m_i, s=14, cmap="RdBu_r",
                   vmin=-np.percentile(np.abs(m_i), 98), vmax=np.percentile(np.abs(m_i), 98))
    plt.colorbar(sc, ax=a, label="per-cell time-mean residual m_i")
    a.set_xlabel("lon"); a.set_ylabel("lat")
    ttl = "(iv) the static map error m_i"
    if fault:
        ttl += f"\nSpearman vs distance-to-CFM5.3 = {fault['spearman_m_i_vs_distance_to_fault']:+.3f}"
    a.set_title(ttl)

    fig.suptitle(f"K-046 the ACF floor  |  {res['verdict']['VERDICT'][:96]}", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FIG, dpi=130)
    print(f"[fig] {FIG}")

    print("\n" + "=" * 78)
    print(f"K-046  |  raw lag-1 excess {exc_raw:+.4f}  ->  DEMEANED lag-1 excess {exc_dm:+.4f}")
    print(f"  between-cell variance fraction = {vs['between_fraction']:.4f}  "
          f"(fitted pedestal C = {ours[2]:.4f}; Kepler predicted {KEPLER_C})")
    print(f"  Moran's I  raw {st_real['moran_raw']:+.5f} -> demeaned {st_real['moran_demeaned']:+.5f} "
          f"(null demeaned p97.5 {NULL['moran_demeaned']['p97_5']:+.5f})")
    print(f"  shape: dBIC(exp - exp+const) real = "
          f"{shp_real_raw.get('dBIC_exp_minus_expconst', float('nan')):+.1f}; "
          f"exp+const wins on {n_sims_expconst_wins}/{len(sims)} sims")
    print(f"  VERDICT: {res['verdict']['VERDICT']}")
    print(f"  runtime {res['runtime_minutes']} min  ->  {OUT.name}, {FIG.name}")
    print("=" * 78)


if __name__ == "__main__":
    main()
