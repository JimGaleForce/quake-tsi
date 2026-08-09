"""EXP-H (confirmatory): temporal ETAS vs Poisson baselines.

Frozen protocol: PATTERN_PROTOCOL.md section "EXP-H".
  - Catalog: SCSN_original_catalog.txt (raw column order; auto-detected).
  - SoCal box lat [31.5, 38.0], lon [-122.0, -113.5].
  - Split: train < 2010-01-01 UTC, test 2010-01-01 -> catalog end.
  - Magnitude floor M >= 2.5 unless Mc is unstable at 2.5 over 1981-2018 (then 3.0).
  - lambda(t) = mu + sum_{t_i<t} K * 10^(alpha (M_i - M0)) * (t - t_i + c)^(-p)
  - MLE on train (L-BFGS-B, log-parameters, multiple starts).
  - TEST: single scoring, parameters frozen, walk-forward (history = all prior events).
  - Success (frozen): gain > 0.5 bits/event vs BOTH Poisson baselines.

Run: python -u exp_h_etas.py
"""

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "xue_lu_zenodo"
CATALOG = "SCSN_original_catalog.txt"
OUT = HERE / "results_exp_h.json"

LAT_MIN, LAT_MAX = 31.5, 38.0
LON_MIN, LON_MAX = -122.0, -113.5
SPLIT = pd.Timestamp("2010-01-01", tz="UTC")
BURN_IN_DAYS = 365.0          # train history-only lead-in before the fit target window
TRUNC_W_DAYS = 1000.0         # truncation window used for the (fast) train fit
RUNTIME_BUDGET_MIN = 90.0
LN2 = np.log(2.0)
LN10 = np.log(10.0)

# multi-start seeds: (mu_frac_of_rate, K, alpha, c, p)
STARTS = [
    (0.50, 0.020, 1.00, 0.010, 1.10),
    (0.20, 0.005, 1.50, 0.050, 1.20),
    (0.80, 0.050, 0.80, 0.005, 1.05),
    (0.35, 0.010, 1.20, 0.020, 1.15),
]
ALPHA_BOUNDS = (0.5, 2.5)
P_BOUNDS = (0.8, 2.0)
MAXITER = 150


# ------------------------------------------------------------------ data loading
def load_catalog(fname):
    """Auto-detect column order (raw vs declustered) exactly as xue_lu_crosstest.load_catalog."""
    raw_cols = ["yr", "mo", "dy", "hr", "mi", "sec", "eid", "lat", "lon", "depth", "mag"]
    dec_cols = ["yr", "mo", "dy", "hr", "mi", "sec", "lat", "lon", "depth", "mag", "eid"]
    probe = pd.read_csv(DATA / fname, sep=r"\s+", header=None, nrows=1000)
    is_raw = probe[6].abs().max() > 90
    cols = raw_cols if is_raw else dec_cols
    df = pd.read_csv(DATA / fname, sep=r"\s+", header=None, names=cols)
    assert df.lat.abs().max() <= 90 and df.lon.abs().max() <= 180, f"column detection failed for {fname}"
    assert df.lat.between(-90, 90).all() and df.lon.between(-180, 180).all()
    sec = df["sec"].astype(float)
    ts = pd.to_datetime(
        dict(year=df.yr, month=df.mo, day=df.dy, hour=df.hr, minute=df.mi, second=0), utc=True
    ) + pd.to_timedelta(sec, unit="s")
    df["ts"] = ts
    # pandas 3.0 safe: no astype(int64)/1e9
    df["t"] = (ts - pd.Timestamp(0, tz="UTC")).dt.total_seconds().to_numpy() / 86400.0
    return df[["ts", "t", "lat", "lon", "depth", "mag"]], ("raw" if is_raw else "declustered")


def max_curvature_mc(mags, binw=0.1):
    """Maximum-curvature Mc: magnitude of the modal bin of the non-cumulative GR histogram."""
    if len(mags) < 200:
        return None
    lo = np.floor(mags.min() * 10) / 10 - binw / 2
    edges = np.arange(lo, mags.max() + binw, binw)
    h, _ = np.histogram(mags, bins=edges)
    i = int(np.argmax(h))
    return float(round((edges[i] + edges[i + 1]) / 2, 2))


def aki_b(mags, mc, binw=0.01):
    m = mags[mags >= mc - 1e-9]
    return float(1.0 / (LN10 * (m.mean() - (mc - binw / 2.0))))


# ------------------------------------------------------------------ ETAS core
def _pair_sums(t, w, m, c, p, tgt_lo, tgt_hi, W, want_grad, chunk):
    """Sum over prior events i (t_i < t_j, and t_j - t_i <= W) for targets j in [tgt_lo, tgt_hi).

    Returns S (= sum w_i (dt+c)^-p) and, if want_grad, the companion sums needed for
    d/da, d/dc, d/dp of S (a = alpha*ln10).
    """
    n_t = tgt_hi - tgt_lo
    S = np.zeros(n_t)
    Sa = np.zeros(n_t) if want_grad else None
    Sc = np.zeros(n_t) if want_grad else None
    Sp = np.zeros(n_t) if want_grad else None
    finite_W = np.isfinite(W)
    for x0 in range(tgt_lo, tgt_hi, chunk):
        x1 = min(x0 + chunk, tgt_hi)
        hlo = int(np.searchsorted(t, t[x0] - W, side="left")) if finite_W else 0
        if hlo >= x1:
            continue
        th = t[hlo:x1]
        dt = t[x0:x1, None] - th[None, :]
        valid = dt > 0
        if finite_W:
            valid &= dt <= W
        dtc = np.where(valid, dt + c, 1.0)
        ldtc = np.log(dtc)
        q = np.exp(-p * ldtc)
        q *= valid
        q *= w[hlo:x1][None, :]
        o = slice(x0 - tgt_lo, x1 - tgt_lo)
        S[o] = q.sum(axis=1)
        if want_grad:
            Sa[o] = (q * m[hlo:x1][None, :]).sum(axis=1)
            Sc[o] = -p * (q / dtc).sum(axis=1)
            Sp[o] = -(q * ldtc).sum(axis=1)
    return S, Sa, Sc, Sp


def _G_terms(t_src, T0, T1, c, p, W, want_grad):
    """Per-source integral of the kernel over the target window, G_i, and dG/dc, dG/dp."""
    lo = np.maximum(T0, t_src)
    hi = np.minimum(T1, t_src + W) if np.isfinite(W) else np.full_like(t_src, T1)
    live = hi > lo
    A = np.where(live, hi - t_src + c, 1.0)
    B = np.where(live, lo - t_src + c, 1.0)
    lA, lB = np.log(A), np.log(B)
    s = 1.0 - p
    small = abs(s) < 1e-6
    if small:
        G = (lA - lB) + s * (lA**2 - lB**2) / 2.0 + s * s * (lA**3 - lB**3) / 6.0
    else:
        As, Bs = np.exp(s * lA), np.exp(s * lB)
        G = (As - Bs) / s
    G *= live
    if not want_grad:
        return G, None, None
    dGdc = (np.exp(-p * lA) - np.exp(-p * lB)) * live
    if small:
        dGds = (lA**2 - lB**2) / 2.0 + s * (lA**3 - lB**3) / 3.0
    else:
        As, Bs = np.exp(s * lA), np.exp(s * lB)
        dGds = ((lA * As - lB * Bs) * s - (As - Bs)) / (s * s)
    dGdp = -dGds * live
    return G, dGdc, dGdp


def etas_ll(theta, t_all, m_all, tgt_lo, tgt_hi, T0, T1, W, want_grad=True, chunk=512):
    """theta = (mu, K, alpha, c, p). Returns (LL, grad wrt theta) or (LL, None)."""
    mu, K, alpha, c, p = theta
    a = alpha * LN10
    w = np.exp(a * m_all)
    S, Sa, Sc, Sp = _pair_sums(t_all, w, m_all, c, p, tgt_lo, tgt_hi, W, want_grad, chunk)
    lam = mu + K * S
    if not np.all(np.isfinite(lam)) or np.any(lam <= 0):
        return -1e18, (np.zeros(5) if want_grad else None)
    src = t_all <= T1
    G, dGdc, dGdp = _G_terms(t_all[src], T0, T1, c, p, W, want_grad)
    wG = w[src] * G
    Lam = mu * (T1 - T0) + K * wG.sum()
    LL = float(np.log(lam).sum() - Lam)
    if not want_grad:
        return LL, None
    inv = 1.0 / lam
    g_mu = inv.sum() - (T1 - T0)
    g_K = (S * inv).sum() - wG.sum()
    g_a = K * (Sa * inv).sum() - K * (m_all[src] * wG).sum()
    g_c = K * (Sc * inv).sum() - K * (w[src] * dGdc).sum()
    g_p = K * (Sp * inv).sum() - K * (w[src] * dGdp).sum()
    return LL, np.array([g_mu, g_K, g_a * LN10, g_c, g_p])


def _obj(x, *args):
    theta = np.exp(x)
    LL, g = etas_ll(theta, *args)
    if not np.isfinite(LL):
        return 1e18, np.zeros(5)
    return -LL, -g * theta


# ------------------------------------------------------------------ main
def main():
    t_start = time.time()
    res = {"experiment": "EXP-H", "protocol": "PATTERN_PROTOCOL.md :: EXP-H",
           "state_class": "first-run", "run_utc": pd.Timestamp.utcnow().isoformat()}

    df, order = load_catalog(CATALOG)
    res["catalog"] = {"file": CATALOG, "detected_column_order": order, "n_rows": int(len(df))}
    box = df.lat.between(LAT_MIN, LAT_MAX) & df.lon.between(LON_MIN, LON_MAX)
    d = df[box].sort_values("t").reset_index(drop=True)
    res["catalog"].update(
        n_in_box=int(len(d)),
        box={"lat": [LAT_MIN, LAT_MAX], "lon": [LON_MIN, LON_MAX]},
        t_first=str(d.ts.iloc[0]), t_last=str(d.ts.iloc[-1]),
    )
    print(f"[data] {order} order, {len(df)} rows, {len(d)} in SoCal box, "
          f"{d.ts.iloc[0]} .. {d.ts.iloc[-1]}")

    # ---------- step 1: Mc check ----------
    mc_ev = {}
    decades = [(1981, 1990), (1990, 2000), (2000, 2010), (2010, 2019)]
    for y0, y1 in decades:
        s = d[(d.ts >= pd.Timestamp(f"{y0}-01-01", tz="UTC")) & (d.ts < pd.Timestamp(f"{y1}-01-01", tz="UTC"))]
        mc = max_curvature_mc(s.mag.to_numpy())
        mc_ev[f"{y0}-{y1}"] = {"n": int(len(s)), "mc_maxcurv": mc,
                               "mc_maxcurv_plus_0.2": None if mc is None else round(mc + 0.2, 2)}
    for lab, sub in [("train(<2010)", d[d.ts < SPLIT]), ("test(>=2010)", d[d.ts >= SPLIT])]:
        mc = max_curvature_mc(sub.mag.to_numpy())
        mc_ev[lab] = {"n": int(len(sub)), "mc_maxcurv": mc,
                      "mc_maxcurv_plus_0.2": None if mc is None else round(mc + 0.2, 2)}
    worst = max(v["mc_maxcurv_plus_0.2"] for v in mc_ev.values() if v["mc_maxcurv_plus_0.2"] is not None)
    M0 = 2.5 if worst <= 2.5 else 3.0
    res["mc_check"] = {"method": "maximum curvature (0.1-mag bins) + 0.2 safety, per decade and per split",
                       "per_period": mc_ev, "worst_case_mc_est": worst,
                       "chosen_floor_M0": M0,
                       "rule": "M0=2.5 if every period's (maxcurv+0.2) <= 2.5, else 3.0",
                       "stable_at_2.5": bool(worst <= 2.5)}
    print(f"[mc] worst (maxcurv+0.2) over periods = {worst}  -> floor M0 = {M0}")
    for k, v in mc_ev.items():
        print(f"      {k:>14s}  n={v['n']:>7d}  Mc~{v['mc_maxcurv']}  (+0.2 -> {v['mc_maxcurv_plus_0.2']})")

    cat = d[d.mag >= M0 - 1e-9].reset_index(drop=True)
    t_all = cat.t.to_numpy(float)
    m_all = (cat.mag.to_numpy(float) - M0)
    is_test = (cat.ts >= SPLIT).to_numpy()
    n_train_all, n_test = int((~is_test).sum()), int(is_test.sum())
    t_split = (SPLIT - pd.Timestamp(0, tz="UTC")).total_seconds() / 86400.0
    T_end = float(t_all[-1])
    print(f"[cat] M>={M0}: {len(cat)} events ({n_train_all} train / {n_test} test)")

    # ---------- step 2: train fit ----------
    T0_hist = float(t_all[0])
    T0_fit = T0_hist + BURN_IN_DAYS
    T1_fit = t_split
    fit_lo = int(np.searchsorted(t_all, T0_fit, side="left"))
    fit_hi = int(np.searchsorted(t_all, T1_fit, side="left"))
    n_fit = fit_hi - fit_lo
    D_fit = T1_fit - T0_fit
    rate_fit = n_fit / D_fit
    subsampled = False

    args_tr = (t_all, m_all, fit_lo, fit_hi, T0_fit, T1_fit, TRUNC_W_DAYS, True, 512)
    t0 = time.time()
    LL0, _ = etas_ll(np.array([rate_fit * 0.5, 0.02, 1.0, 0.01, 1.1]), *args_tr)
    ev_s = time.time() - t0
    proj_min = ev_s * MAXITER * len(STARTS) / 60.0
    print(f"[fit] window {pd.Timestamp(T0_fit*86400,unit='s',tz='UTC')} .. {SPLIT}  "
          f"n_fit={n_fit}  rate={rate_fit:.3f}/d  |  1 fun+grad = {ev_s:.2f}s  proj {proj_min:.1f} min")

    if proj_min > RUNTIME_BUDGET_MIN:
        subsampled = True
        T0_hist = T1_fit - 15 * 365.25
        T0_fit = T0_hist + BURN_IN_DAYS
        hist_lo = int(np.searchsorted(t_all, T0_hist, side="left"))
        # keep the full array but restrict targets AND sources: rebuild a train-only view
        t_tr = t_all[hist_lo:fit_hi]
        m_tr = m_all[hist_lo:fit_hi]
        fit_lo = int(np.searchsorted(t_tr, T0_fit, side="left"))
        fit_hi = len(t_tr)
        n_fit = fit_hi - fit_lo
        D_fit = T1_fit - T0_fit
        rate_fit = n_fit / D_fit
        args_tr = (t_tr, m_tr, fit_lo, fit_hi, T0_fit, T1_fit, TRUNC_W_DAYS, True, 512)
        t0 = time.time()
        etas_ll(np.array([rate_fit * 0.5, 0.02, 1.0, 0.01, 1.1]), *args_tr)
        ev_s = time.time() - t0
        proj_min = ev_s * MAXITER * len(STARTS) / 60.0
        print(f"[fit] RUNTIME GUARD TRIGGERED -> train fit restricted to most recent 15 yr "
              f"(n_fit={n_fit}, 1 eval {ev_s:.2f}s, proj {proj_min:.1f} min)")

    lb = np.log([1e-6, 1e-6, ALPHA_BOUNDS[0], 1e-6, P_BOUNDS[0]])
    ub = np.log([1e3, 1e2, ALPHA_BOUNDS[1], 1e1, P_BOUNDS[1]])
    bounds = list(zip(lb, ub))

    starts_out, best = [], None
    for si, (mu_f, K, al, c, p) in enumerate(STARTS):
        x0 = np.log(np.clip([mu_f * rate_fit, K, al, c, p], np.exp(lb) * 1.001, np.exp(ub) * 0.999))
        ts_ = time.time()
        r = minimize(_obj, x0, args=args_tr, jac=True, method="L-BFGS-B",
                     bounds=bounds, options={"maxiter": MAXITER, "maxfun": MAXITER * 2, "ftol": 1e-12, "gtol": 1e-8})
        th = np.exp(r.x)
        rec = {"start": {"mu": mu_f * rate_fit, "K": K, "alpha": al, "c": c, "p": p},
               "LL": float(-r.fun), "params": dict(zip(["mu", "K", "alpha", "c", "p"], map(float, th))),
               "nit": int(r.nit), "success": bool(r.success), "message": str(r.message),
               "seconds": round(time.time() - ts_, 1)}
        starts_out.append(rec)
        print(f"[fit] start {si+1}: LL={rec['LL']:.2f}  mu={th[0]:.4f} K={th[1]:.5f} "
              f"alpha={th[2]:.3f} c={th[3]:.5f} p={th[4]:.4f}  ({rec['nit']} it, {rec['seconds']}s)")
        if best is None or rec["LL"] > best["LL"]:
            best = rec

    theta_trunc = np.array([best["params"][k] for k in ["mu", "K", "alpha", "c", "p"]])
    print(f"[fit] best truncated-kernel LL = {best['LL']:.2f}")

    # truncation-error check + untruncated refinement (protocol model has no truncation)
    LL_tr_at_opt, _ = etas_ll(theta_trunc, *args_tr[:6], TRUNC_W_DAYS, False, 512)
    args_full = (args_tr[0], args_tr[1], args_tr[2], args_tr[3], args_tr[4], args_tr[5], np.inf, True, 512)
    t0 = time.time()
    LL_full_at_opt, _ = etas_ll(theta_trunc, *args_full[:6], np.inf, False, 1024)
    full_ev_s = time.time() - t0
    print(f"[trunc] LL(W=1000d)={LL_tr_at_opt:.2f}  LL(untruncated)={LL_full_at_opt:.2f}  "
          f"diff={LL_full_at_opt-LL_tr_at_opt:.2f}  (1 untrunc eval {full_ev_s:.1f}s)")

    refine_iters = int(max(5, min(25, (RUNTIME_BUDGET_MIN * 60 - (time.time() - t_start)) * 0.4 / max(full_ev_s * 3, 1))))
    print(f"[refine] untruncated L-BFGS-B polish from truncated optimum, maxiter={refine_iters}")
    tr = time.time()
    rf = minimize(_obj, np.log(theta_trunc), args=args_full, jac=True, method="L-BFGS-B",
                  bounds=bounds, options={"maxiter": refine_iters, "maxfun": refine_iters * 2, "ftol": 1e-12})
    theta = np.exp(rf.x)
    LL_refined = float(-rf.fun)
    print(f"[refine] LL={LL_refined:.2f}  mu={theta[0]:.4f} K={theta[1]:.5f} alpha={theta[2]:.3f} "
          f"c={theta[3]:.5f} p={theta[4]:.4f}  ({rf.nit} it, {time.time()-tr:.0f}s)")
    if LL_refined < LL_full_at_opt:   # never accept a worse point than the start
        theta, LL_refined = theta_trunc, LL_full_at_opt
        print("[refine] refinement did not improve; keeping truncated-fit parameters")

    mu, K, alpha, c, p = map(float, theta)
    b_train = aki_b(cat.mag.to_numpy()[~is_test], M0)
    br = None
    if p > 1.0 and b_train > alpha:
        br = float(K * c ** (1 - p) / (p - 1) * (b_train / (b_train - alpha)))

    res["train_fit"] = {
        "history_start": str(pd.Timestamp(T0_hist * 86400, unit="s", tz="UTC")),
        "fit_window_start": str(pd.Timestamp(T0_fit * 86400, unit="s", tz="UTC")),
        "fit_window_end": str(SPLIT), "burn_in_days": BURN_IN_DAYS,
        "n_target_events": int(n_fit), "fit_window_days": float(D_fit),
        "train_window_rate_per_day": float(rate_fit),
        "truncation_window_days_for_fit": TRUNC_W_DAYS,
        "truncation_error_LL_at_optimum": {"LL_truncated": LL_tr_at_opt, "LL_untruncated": LL_full_at_opt,
                                           "delta_LL": float(LL_full_at_opt - LL_tr_at_opt),
                                           "delta_LL_per_event": float((LL_full_at_opt - LL_tr_at_opt) / n_fit)},
        "starts": starts_out,
        "best_truncated": best["params"], "best_truncated_LL": best["LL"],
        "final_untruncated_refine": {"maxiter": refine_iters, "nit": int(rf.nit), "LL": LL_refined},
        "frozen_params": {"mu": mu, "K": K, "alpha": alpha, "c": c, "p": p, "M0": M0},
        "b_value_train_aki": b_train,
        "branching_ratio_n": br,
        "branching_ratio_formula": "n = K*c^(1-p)/(p-1) * b/(b-alpha)  [valid for p>1 and b>alpha]",
        "train_subsampled_to_recent_15yr": subsampled,
    }
    print(f"[params FROZEN] mu={mu:.4f}/d K={K:.5f} alpha={alpha:.3f} c={c:.5f} p={p:.4f} | "
          f"b={b_train:.3f} branching n={br}")

    # ---------- step 3: TEST (single scoring, walk-forward) ----------
    T0_te, T1_te = t_split, T_end
    D_te = T1_te - T0_te
    te_lo = int(np.searchsorted(t_all, T0_te, side="left"))
    te_hi = len(t_all)
    assert te_hi - te_lo == n_test

    print(f"[test] window {SPLIT} .. {pd.Timestamp(T1_te*86400,unit='s',tz='UTC')}  "
          f"({D_te:.1f} d, {n_test} events); computing walk-forward lambda (untruncated)...")
    t0 = time.time()
    w_all = np.exp(alpha * LN10 * m_all)
    S_te, _, _, _ = _pair_sums(t_all, w_all, m_all, c, p, te_lo, te_hi, np.inf, False, 1024)
    lam_te = mu + K * S_te
    G_te, _, _ = _G_terms(t_all, T0_te, T1_te, c, p, np.inf, False)
    Lam_te = mu * D_te + K * float((w_all * G_te).sum())
    LL_etas = float(np.log(lam_te).sum() - Lam_te)
    print(f"[test] sum log lambda = {np.log(lam_te).sum():.2f}, integral = {Lam_te:.2f}, "
          f"LL_ETAS = {LL_etas:.2f}  ({time.time()-t0:.0f}s)")

    rate_train_base = rate_fit
    rate_test_oracle = n_test / D_te
    baselines = {}
    for name, lam0 in [("poisson_train_rate", rate_train_base), ("poisson_test_rate_oracle", rate_test_oracle)]:
        LLp = n_test * np.log(lam0) - lam0 * D_te
        gain = (LL_etas - LLp) / n_test / LN2
        baselines[name] = {"rate_per_day": float(lam0), "LL": float(LLp),
                           "gain_bits_per_event": float(gain)}
        print(f"[test] vs {name}: rate={lam0:.4f}/d  LL={LLp:.2f}  gain={gain:+.3f} bits/event")

    # magnitude-band breakdown
    mags_te = cat.mag.to_numpy()[te_lo:te_hi]
    bands = {"2.5-3": (M0, 3.0), "3-4": (3.0, 4.0), ">=4": (4.0, 99.0)}
    if M0 >= 3.0:
        bands = {"3-4": (3.0, 4.0), ">=4": (4.0, 99.0)}
    band_out = {}
    for bn, (b0, b1) in bands.items():
        sel = (mags_te >= b0 - 1e-9) & (mags_te < b1)
        nb = int(sel.sum())
        e = {"n": nb}
        if nb:
            for name, lam0 in [("poisson_train_rate", rate_train_base),
                               ("poisson_test_rate_oracle", rate_test_oracle)]:
                logterm = float(np.mean(np.log(lam_te[sel]) - np.log(lam0)) / LN2)
                LLp = n_test * np.log(lam0) - lam0 * D_te
                integ_corr = float(-(Lam_te - lam0 * D_te) / n_test / LN2)
                e[name] = {"log_term_bits_per_event": logterm,
                           "with_uniform_integral_correction_bits": logterm + integ_corr}
        band_out[bn] = e
        print(f"[test] band {bn:>6s} n={nb:>5d}  " + "  ".join(
            f"{k.split('_')[1]}:{v['with_uniform_integral_correction_bits']:+.2f}b"
            for k, v in e.items() if isinstance(v, dict)))

    trig_frac_events = float(np.mean((lam_te - mu) / lam_te))
    trig_frac_integral = float((Lam_te - mu * D_te) / Lam_te)
    gains = [baselines[k]["gain_bits_per_event"] for k in baselines]
    passed = bool(min(gains) > 0.5)

    res["test"] = {
        "window_start": str(SPLIT), "window_end": str(pd.Timestamp(T1_te * 86400, unit="s", tz="UTC")),
        "window_days": float(D_te), "n_test_events": n_test,
        "scoring": "walk-forward, params frozen from train, history = all prior events (train + earlier test), untruncated kernel",
        "sum_log_lambda": float(np.log(lam_te).sum()), "integral_lambda": float(Lam_te),
        "LL_etas": LL_etas, "baselines": baselines,
        "magnitude_bands": band_out,
        "band_note": ("log_term = mean[log2 lambda_ETAS - log2 lambda_poisson] over band events; "
                      "the integral difference is a single global term, apportioned uniformly per test event "
                      "in with_uniform_integral_correction_bits"),
        "triggered_fraction_mean_over_test_events": trig_frac_events,
        "triggered_fraction_of_integrated_intensity": trig_frac_integral,
        "lambda_stats": {"min": float(lam_te.min()), "median": float(np.median(lam_te)),
                         "max": float(lam_te.max()), "mean": float(lam_te.mean())},
    }
    res["success_rule"] = {"rule": "gain > 0.5 bits/event vs BOTH baselines",
                           "gain_vs_train_rate": baselines["poisson_train_rate"]["gain_bits_per_event"],
                           "gain_vs_test_rate_oracle": baselines["poisson_test_rate_oracle"]["gain_bits_per_event"],
                           "PASS": passed}
    res["runtime_minutes"] = round((time.time() - t_start) / 60.0, 2)
    res["flags"] = {"train_subsampled_to_recent_15yr": subsampled,
                    "truncated_kernel_used_for_fit_only": True,
                    "test_scoring_untruncated": True}

    OUT.write_text(json.dumps(res, indent=2))
    print("\n" + "=" * 72)
    print(f"EXP-H  |  Mc floor M0={M0} (stable={worst<=2.5}, worst maxcurv+0.2={worst})")
    print(f"  params: mu={mu:.4f}/d  K={K:.5f}  alpha={alpha:.3f}  c={c:.5f} d  p={p:.4f}  "
          f"(b={b_train:.3f}, n={br if br is None else round(br,3)})")
    print(f"  test: {n_test} events over {D_te:.0f} d;  triggered fraction "
          f"{trig_frac_integral:.3f} (integral) / {trig_frac_events:.3f} (mean over events)")
    print(f"  gain vs Poisson(train rate) = {baselines['poisson_train_rate']['gain_bits_per_event']:+.3f} bits/event")
    print(f"  gain vs Poisson(test rate)  = {baselines['poisson_test_rate_oracle']['gain_bits_per_event']:+.3f} bits/event")
    print(f"  SUCCESS RULE (>0.5 bits vs BOTH): {'PASS' if passed else 'FAIL'}")
    print(f"  runtime {res['runtime_minutes']} min  ->  {OUT.name}")
    print("=" * 72)


if __name__ == "__main__":
    main()
