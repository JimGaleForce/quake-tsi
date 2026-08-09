"""EXP-I: patterns-within-patterns / meta-patterns (PATTERN_PROTOCOL.md, frozen 2026-08-09).

(i)   WITHIN  - sequence-shape recurrence: per-sequence Omori (K,c,p) + Aki b on train
                sequences (>=30 aftershocks); variance decomposition of (p,b) over
                0.5-degree neighborhoods (ICC-like ratio + 1000-shuffle permutation p);
                TEST: neighborhood-mean vs global-mean prediction of (p,b).
(ii)  OUTSIDE - meta-clustering of M>=5 mainshock times: CV of interevent times +
                Ripley-K-in-time at {30,100,365} d vs 1000 Poisson sims (train);
                hazard rule X* mined on train, frozen, scored on test.
(iii) DRIFT   - b-value (M>=2.5) and interevent CV (M>=3) in 5-yr train windows;
                monotone trend (|Spearman rho| > 0.8) frozen and checked on test.

Inputs are never modified. Output: results_exp_i.json + console summary.
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

HERE = Path(__file__).parent
DATA = HERE / "data" / "xue_lu_zenodo"
CATALOG = "SCSN_original_catalog.txt"

EARTH_R = 6371.0
BOX = dict(lat=(31.5, 38.0), lon=(-122.0, -113.5))
SPLIT = pd.Timestamp("2010-01-01", tz="UTC")
SPLIT_UNIX = (SPLIT - pd.Timestamp(0, tz="UTC")).total_seconds()
DAY = 86400.0

MS_MIN_MAG = 5.0
AFT_WINDOW_D = 100.0
AFT_RADIUS_KM = 50.0
MIN_AFT = 30
B_FLOOR = 2.5
B_BIN = 0.1
MIN_B_EVENTS = 30
NEIGH_DEG = 0.5
N_PERM = 1000
N_SIM = 1000
HAZARD_X_DAYS = [7.0, 30.0, 90.0, 365.0]
RIPLEY_LAGS = [30.0, 100.0, 365.0]

RNG = np.random.default_rng(20260809)


# ------------------------------------------------------------------ catalog
def load_catalog(fname):
    """Auto-detect column order (raw vs declustered) exactly like xue_lu_crosstest.py."""
    raw_cols = ["yr", "mo", "dy", "hr", "mi", "sec", "eid", "lat", "lon", "depth", "mag"]
    dec_cols = ["yr", "mo", "dy", "hr", "mi", "sec", "lat", "lon", "depth", "mag", "eid"]
    probe = pd.read_csv(DATA / fname, sep=r"\s+", header=None, nrows=1000)
    cols = raw_cols if probe[6].abs().max() > 90 else dec_cols
    df = pd.read_csv(DATA / fname, sep=r"\s+", header=None, names=cols)
    assert df.lat.abs().max() <= 90 and df.lon.abs().max() <= 180, "column detection failed"
    sec = df["sec"].astype(float)
    ts = pd.to_datetime(
        dict(year=df.yr, month=df.mo, day=df.dy, hour=df.hr, minute=df.mi, second=0), utc=True
    ) + pd.to_timedelta(sec, unit="s")
    df["t_unix"] = (ts - pd.Timestamp(0, tz="UTC")).dt.total_seconds().to_numpy()
    df = df[["t_unix", "lat", "lon", "depth", "mag"]].sort_values("t_unix").reset_index(drop=True)
    m = (
        df.lat.between(*BOX["lat"])
        & df.lon.between(*BOX["lon"])
    )
    return df[m].reset_index(drop=True)


def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1, la2, lo2 = map(np.radians, [lat1, lon1, lat2, lon2])
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * EARTH_R * np.arcsin(np.sqrt(a))


# ------------------------------------------------------------------ sequences (EXP-G definition)
def build_sequences(df):
    """Mainshocks M>=5.0; aftershocks = events within 100 d after and 50 km, magnitude below
    the mainshock, not already claimed by a larger earlier mainshock within its window."""
    t = df.t_unix.to_numpy()
    lat = df.lat.to_numpy()
    lon = df.lon.to_numpy()
    mag = df.mag.to_numpy()
    ms_idx = np.where(mag >= MS_MIN_MAG)[0]  # chronological

    owner = np.full(len(df), -1, dtype=int)       # index of owning mainshock
    owner_mag = np.full(len(df), -np.inf)         # magnitude of owning mainshock
    win = AFT_WINDOW_D * DAY

    for j in ms_idx:
        lo = np.searchsorted(t, t[j], side="right")
        hi = np.searchsorted(t, t[j] + win, side="right")
        if hi <= lo:
            continue
        cand = np.arange(lo, hi)
        cand = cand[mag[cand] < mag[j]]
        if cand.size == 0:
            continue
        d = haversine_km(lat[j], lon[j], lat[cand], lon[cand])
        cand = cand[d <= AFT_RADIUS_KM]
        if cand.size == 0:
            continue
        # claim unless already owned by a *larger* earlier mainshock
        take = cand[owner_mag[cand] <= mag[j]]
        owner[take] = j
        owner_mag[take] = mag[j]

    seqs = []
    for j in ms_idx:
        members = np.where(owner == j)[0]
        seqs.append(
            dict(
                ms_index=int(j),
                t_unix=float(t[j]),
                lat=float(lat[j]),
                lon=float(lon[j]),
                mag=float(mag[j]),
                n_aft=int(members.size),
                dt_days=(t[members] - t[j]) / DAY,
                aft_mag=mag[members],
            )
        )
    return seqs


# ------------------------------------------------------------------ Omori + b
def omori_mle(dt_days, T=AFT_WINDOW_D):
    """MLE of modified-Omori rate K/(t+c)^p on [0,T] given aftershock delays (days).
    K is profiled out analytically; (c,p) optimized in log space."""
    t = np.asarray(dt_days, float)
    t = t[(t > 0) & (t <= T)]
    n = t.size
    if n < 2:
        return None

    def integral(c, p):
        if abs(p - 1.0) < 1e-8:
            return np.log((T + c) / c)
        return ((T + c) ** (1 - p) - c ** (1 - p)) / (1 - p)

    def nll(theta):
        logc, p = theta
        c = np.exp(logc)
        I = integral(c, p)
        if not np.isfinite(I) or I <= 0:
            return 1e12
        # profiled K = n / I  ->  ll = n*log(K) - p*sum(log(t+c)) - K*I
        return -(n * np.log(n / I) - p * np.log(t + c).sum() - n)

    best = None
    for logc0 in (np.log(0.001), np.log(0.01), np.log(0.1), np.log(1.0)):
        for p0 in (0.8, 1.0, 1.3):
            r = optimize.minimize(
                nll, [logc0, p0], method="L-BFGS-B",
                bounds=[(np.log(1e-5), np.log(10.0)), (0.2, 3.0)],
            )
            if r.success and (best is None or r.fun < best.fun):
                best = r
    if best is None:
        return None
    c = float(np.exp(best.x[0]))
    p = float(best.x[1])
    K = float(n / integral(c, p))
    return dict(K=K, c=c, p=p, n_fit=int(n), nll=float(best.fun))


def aki_b(mags, floor=B_FLOOR, binw=B_BIN):
    m = np.asarray(mags, float)
    m = m[m >= floor - 1e-9]
    if m.size < MIN_B_EVENTS:
        return None, int(m.size)
    mc_eff = floor - binw / 2.0
    denom = m.mean() - mc_eff
    if denom <= 0:
        return None, int(m.size)
    return float(1.0 / (np.log(10) * denom)), int(m.size)


def fit_sequences(seqs):
    out = []
    for s in seqs:
        if s["n_aft"] < MIN_AFT:
            continue
        om = omori_mle(s["dt_days"])
        if om is None:
            continue
        b, n_b = aki_b(s["aft_mag"])
        out.append(
            dict(
                t_unix=s["t_unix"], lat=s["lat"], lon=s["lon"], mag=s["mag"],
                n_aft=s["n_aft"], K=om["K"], c=om["c"], p=om["p"],
                b=b, n_b=n_b,
            )
        )
    return out


# ------------------------------------------------------------------ (i)
def icc_ratio(values, labels):
    """Fraction of total sum-of-squares lying between groups (ICC-like)."""
    v = np.asarray(values, float)
    lab = np.asarray(labels)
    grand = v.mean()
    sst = ((v - grand) ** 2).sum()
    if sst <= 0:
        return np.nan
    ssb = 0.0
    for g in np.unique(lab):
        gi = v[lab == g]
        ssb += gi.size * (gi.mean() - grand) ** 2
    return float(ssb / sst)


def perm_p(values, labels, n_perm=N_PERM):
    obs = icc_ratio(values, labels)
    lab = np.asarray(labels)
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = icc_ratio(values, RNG.permutation(lab))
    return obs, float((1 + (null >= obs).sum()) / (n_perm + 1)), float(np.nanmean(null))


def cell_label(lat, lon, size=NEIGH_DEG):
    return f"{int(np.floor(lat / size))}_{int(np.floor(lon / size))}"


def part_i(train_fits, test_fits, log):
    res = {"n_train_sequences_fitted": len(train_fits), "n_test_sequences_fitted": len(test_fits)}

    # --- variance decomposition on train
    for name in ("p", "b"):
        vals, labs = [], []
        for f in train_fits:
            if f[name] is None:
                continue
            vals.append(f[name])
            labs.append(cell_label(f["lat"], f["lon"]))
        vals = np.array(vals, float)
        labs = np.array(labs)
        n_cells = len(np.unique(labs))
        n_multi = int(sum(1 for g in np.unique(labs) if (labs == g).sum() >= 2))
        if vals.size >= 3 and n_cells >= 2:
            obs, pval, nullmean = perm_p(vals, labs)
        else:
            obs, pval, nullmean = np.nan, np.nan, np.nan
        res[f"train_{name}"] = dict(
            n=int(vals.size), n_cells=n_cells, n_cells_with_2plus=n_multi,
            mean=float(np.mean(vals)) if vals.size else None,
            sd=float(np.std(vals, ddof=1)) if vals.size > 1 else None,
            icc_like_ssb_over_sst=None if np.isnan(obs) else obs,
            perm_null_mean=None if np.isnan(nullmean) else nullmean,
            perm_p=None if np.isnan(pval) else pval,
            n_permutations=N_PERM,
        )

    # --- prediction on test
    tr_lat = np.array([f["lat"] for f in train_fits])
    tr_lon = np.array([f["lon"] for f in train_fits])
    tr_p = np.array([f["p"] for f in train_fits])
    tr_b = np.array([np.nan if f["b"] is None else f["b"] for f in train_fits])
    glob_p = float(np.mean(tr_p))
    glob_b = float(np.nanmean(tr_b))
    sd_p = float(np.std(tr_p, ddof=1))
    sd_b = float(np.nanstd(tr_b, ddof=1))

    wins_p, wins_b, wins_joint = [], [], []
    detail = []
    n_nb_ok = 0
    for f in test_fits:
        d = np.hypot(tr_lat - f["lat"], tr_lon - f["lon"])
        sel = d <= NEIGH_DEG
        if sel.sum() < 2:
            continue
        n_nb_ok += 1
        nb_p = float(np.mean(tr_p[sel]))
        bsel = tr_b[sel]
        nb_b = float(np.nanmean(bsel)) if np.isfinite(bsel).sum() >= 2 else None
        rec = dict(t_unix=f["t_unix"], lat=f["lat"], lon=f["lon"], mag=f["mag"],
                   n_train_neighbors=int(sel.sum()), obs_p=f["p"], obs_b=f["b"],
                   nb_p=nb_p, nb_b=nb_b)
        e_nb_p = (f["p"] - nb_p) ** 2
        e_gl_p = (f["p"] - glob_p) ** 2
        rec["win_p"] = bool(e_nb_p < e_gl_p)
        wins_p.append(rec["win_p"])
        z_nb, z_gl = e_nb_p / sd_p ** 2, e_gl_p / sd_p ** 2
        if f["b"] is not None and nb_b is not None:
            e_nb_b = (f["b"] - nb_b) ** 2
            e_gl_b = (f["b"] - glob_b) ** 2
            rec["win_b"] = bool(e_nb_b < e_gl_b)
            wins_b.append(rec["win_b"])
            z_nb += e_nb_b / sd_b ** 2
            z_gl += e_gl_b / sd_b ** 2
        rec["win_joint"] = bool(z_nb < z_gl)
        wins_joint.append(rec["win_joint"])
        detail.append(rec)

    def binom(wins):
        k, n = int(sum(wins)), len(wins)
        if n == 0:
            return dict(n=0, wins=0, frac=None, binom_p_greater=None)
        bt = stats.binomtest(k, n, 0.5, alternative="greater")
        return dict(n=n, wins=k, frac=k / n, binom_p_greater=float(bt.pvalue))

    res["test_prediction"] = dict(
        n_test_predictable=n_nb_ok,
        global_train_mean_p=glob_p, global_train_mean_b=glob_b,
        train_sd_p=sd_p, train_sd_b=sd_b,
        p_param=binom(wins_p), b_param=binom(wins_b), joint=binom(wins_joint),
        per_sequence=detail,
    )

    j = res["test_prediction"]["joint"]
    if j["n"] == 0:
        verdict = ("FAIL/UNTESTABLE: no test sequence had >=2 train sequences within 0.5 deg "
                   "(n_test_predictable=0), so the frozen success criterion cannot be evaluated.")
    else:
        ok = (j["frac"] >= 0.60) and (j["binom_p_greater"] is not None and j["binom_p_greater"] < 0.05)
        verdict = (
            f"{'PASS' if ok else 'FAIL'}: neighborhood beat global on {j['wins']}/{j['n']} "
            f"({j['frac']:.1%}) predictable test sequences (joint standardized SE), binomial "
            f"p={j['binom_p_greater']:.3g}; criterion was >=60% and p<0.05."
        )
    res["verdict"] = verdict
    log(f"[i] {verdict}")
    return res


# ------------------------------------------------------------------ (ii)
def ripley_k_time(t_days, lags, span):
    n = t_days.size
    if n < 2:
        return {str(h): None for h in lags}
    ts = np.sort(t_days)
    out = {}
    for h in lags:
        lo = np.searchsorted(ts, ts - h, side="left")
        hi = np.searchsorted(ts, ts + h, side="right")
        cnt = (hi - lo - 1).sum()  # ordered pairs within h, self excluded
        out[str(h)] = float(span * cnt / (n * (n - 1)))
    return out


def part_ii(df, log, n_sim=N_SIM):
    ms = df[df.mag >= MS_MIN_MAG].reset_index(drop=True)
    tr = ms[ms.t_unix < SPLIT_UNIX]
    te = ms[ms.t_unix >= SPLIT_UNIX]
    cat_t0, cat_t1 = df.t_unix.min(), df.t_unix.max()
    tr_span = (SPLIT_UNIX - cat_t0) / DAY
    te_span = (cat_t1 - SPLIT_UNIX) / DAY

    tt = (tr.t_unix.to_numpy() - cat_t0) / DAY
    n_tr = tt.size
    dtt = np.diff(tt)
    cv = float(np.std(dtt, ddof=1) / np.mean(dtt)) if dtt.size > 1 else None

    k_obs = ripley_k_time(tt, RIPLEY_LAGS, tr_span)

    cv_null = np.empty(n_sim)
    k_null = {str(h): np.empty(n_sim) for h in RIPLEY_LAGS}
    for i in range(n_sim):
        s = np.sort(RNG.uniform(0, tr_span, n_tr))
        d = np.diff(s)
        cv_null[i] = np.std(d, ddof=1) / np.mean(d)
        kk = ripley_k_time(s, RIPLEY_LAGS, tr_span)
        for h in RIPLEY_LAGS:
            k_null[str(h)][i] = kk[str(h)]

    res = {
        "train": dict(
            n_mainshocks=int(n_tr), span_days=float(tr_span),
            rate_per_day=float(n_tr / tr_span),
            interevent_cv=cv,
            cv_poisson_null_mean=float(cv_null.mean()),
            cv_p_greater=float((1 + (cv_null >= cv).sum()) / (n_sim + 1)) if cv is not None else None,
            ripley_k=k_obs,
            ripley_k_poisson_expectation={str(h): 2.0 * h for h in RIPLEY_LAGS},
            ripley_k_null_mean={str(h): float(k_null[str(h)].mean()) for h in RIPLEY_LAGS},
            ripley_k_p_greater={
                str(h): float((1 + (k_null[str(h)] >= k_obs[str(h)]).sum()) / (n_sim + 1))
                for h in RIPLEY_LAGS
            },
            n_simulations=n_sim,
        )
    }

    # ---- mine hazard rule on train
    mining = []
    for X in HAZARD_X_DAYS:
        elig = tt[tt <= tr_span - X]  # full X-day window inside train
        if elig.size == 0:
            mining.append(dict(X_days=X, n_eligible=0))
            continue
        hit = np.array([((tt > a) & (tt <= a + X)).any() for a in elig])
        p_obs = float(hit.mean())
        rate = n_tr / tr_span
        p_pois = float(1 - np.exp(-rate * X))
        mining.append(dict(X_days=X, n_eligible=int(elig.size), n_hits=int(hit.sum()),
                           p_observed=p_obs, p_poisson=p_pois,
                           ratio=float(p_obs / p_pois) if p_pois > 0 else None))
    scored = [m for m in mining if m.get("ratio") is not None]
    x_star = max(scored, key=lambda m: m["ratio"])["X_days"] if scored else None
    res["train"]["hazard_mining"] = mining
    res["X_star_days"] = x_star

    # ---- score frozen X* on test
    tte = (te.t_unix.to_numpy() - SPLIT_UNIX) / DAY
    n_te = tte.size
    if x_star is None or n_te < 2:
        res["test"] = dict(n_mainshocks=int(n_te), note="insufficient test mainshocks")
        res["verdict"] = "FAIL/UNTESTABLE: could not score a frozen hazard rule on the test period."
        log("[ii] " + res["verdict"])
        return res

    elig = tte[tte <= te_span - x_star]
    hit = np.array([((tte > a) & (tte <= a + x_star)).any() for a in elig]) if elig.size else np.array([], bool)
    rate_te = n_te / te_span
    p0 = float(1 - np.exp(-rate_te * x_star))
    if elig.size:
        bt = stats.binomtest(int(hit.sum()), int(elig.size), p0, alternative="greater")
        pv = float(bt.pvalue)
    else:
        pv = None
    res["test"] = dict(
        n_mainshocks=int(n_te), span_days=float(te_span), rate_per_day=float(rate_te),
        X_star_days=x_star, n_eligible=int(elig.size), n_hits=int(hit.sum()) if elig.size else 0,
        p_observed=float(hit.mean()) if elig.size else None,
        poisson_base_rate=p0, binom_p_greater=pv,
    )
    ok = pv is not None and pv < 0.05
    res["verdict"] = (
        f"{'PASS' if ok else 'FAIL'}: with X*={x_star:g} d frozen from train, "
        f"P(next M>=5 within X* | M>=5) on test = {res['test']['p_observed']:.3f} "
        f"({res['test']['n_hits']}/{res['test']['n_eligible']}) vs Poisson base {p0:.3f}, "
        f"binomial p={pv:.3g}."
    )
    log("[ii] " + res["verdict"])
    return res


# ------------------------------------------------------------------ (iii)
def window_stats(df, t0, t1):
    w = df[(df.t_unix >= t0) & (df.t_unix < t1)]
    b, n_b = aki_b(w.mag.to_numpy())
    t3 = np.sort(w[w.mag >= 3.0].t_unix.to_numpy()) / DAY
    if t3.size > 2:
        d = np.diff(t3)
        cv = float(np.std(d, ddof=1) / np.mean(d))
    else:
        cv = None
    return dict(n_events=int(len(w)), b=b, n_b=n_b, cv_m3=cv, n_m3=int((w.mag >= 3.0).sum()))


def part_iii(df, log):
    t0 = df.t_unix.min()
    t_end = df.t_unix.max()
    W = 5 * 365.25 * DAY

    def edges(a, b):
        out, s = [], a
        while s < b - 1:
            e = min(s + W, b)
            out.append((s, e, bool((e - s) < W - 1)))  # third field: partial window
            s += W
        return out

    train_edges = edges(t0, SPLIT_UNIX)
    test_edges = edges(SPLIT_UNIX, t_end)

    def label(a, b):
        return (pd.Timestamp(a, unit="s", tz="UTC").strftime("%Y-%m-%d") + ".." +
                pd.Timestamp(b, unit="s", tz="UTC").strftime("%Y-%m-%d"))

    tr_w = [dict(window=label(a, b), partial=pt, **window_stats(df, a, b))
            for a, b, pt in train_edges]
    te_w = [dict(window=label(a, b), partial=pt, **window_stats(df, a, b))
            for a, b, pt in test_edges]

    res = dict(train_windows=tr_w, test_windows=te_w,
               n_train_windows=len(tr_w), n_test_windows=len(te_w))
    verdicts = []
    for key in ("b", "cv_m3"):
        # trend mined on FULL 5-yr train windows only; partial trailing windows are
        # reported descriptively and used only in the test continuation check.
        vals = [w[key] for w in tr_w if w[key] is not None and not w["partial"]]
        idx = np.arange(len(vals), dtype=float)
        if len(vals) >= 3:
            rho, pv = stats.spearmanr(idx, vals)
        else:
            rho, pv = np.nan, np.nan
        entry = dict(n_windows=len(vals), values=vals,
                     spearman_rho=None if np.isnan(rho) else float(rho),
                     spearman_p=None if np.isnan(pv) else float(pv),
                     monotone_trend_frozen=bool(np.isfinite(rho) and abs(rho) > 0.8))
        if entry["monotone_trend_frozen"]:
            sign = int(np.sign(rho))
            entry["frozen_sign"] = sign
            tvals = [w[key] for w in te_w if w[key] is not None]
            entry["test_values"] = tvals
            last_train = vals[-1]
            if len(tvals) >= 2:
                trho, _ = stats.spearmanr(np.arange(len(tvals), dtype=float), tvals)
                entry["test_spearman_rho"] = float(trho)
            else:
                entry["test_spearman_rho"] = None
            cont = [bool(np.sign(v - last_train) == sign) for v in tvals]
            entry["test_beyond_last_train_in_frozen_direction"] = cont
            entry["continuation_holds"] = bool(len(cont) > 0 and all(cont))
            verdicts.append(
                f"{key}: monotone train trend (rho={rho:.2f}) frozen sign {sign:+d}; "
                f"test windows continue in that direction: "
                f"{'YES' if entry['continuation_holds'] else 'NO'} ({cont})."
            )
        else:
            verdicts.append(
                f"{key}: no monotone train trend (rho="
                f"{'nan' if np.isnan(rho) else f'{rho:.2f}'}, |rho|<=0.8) - descriptive only."
            )
        res[key] = entry
    res["verdict"] = " ".join(verdicts)
    log("[iii] " + res["verdict"])
    return res


# ------------------------------------------------------------------ main
def main():
    t_start = time.time()
    lines = []

    def log(msg):
        print(msg, flush=True)
        lines.append(msg)

    df = load_catalog(CATALOG)
    log(f"catalog: {len(df)} events in SoCal box, "
        f"{pd.Timestamp(df.t_unix.min(), unit='s', tz='UTC').date()} .. "
        f"{pd.Timestamp(df.t_unix.max(), unit='s', tz='UTC').date()}")

    seqs = build_sequences(df)
    tr_seqs = [s for s in seqs if s["t_unix"] < SPLIT_UNIX]
    te_seqs = [s for s in seqs if s["t_unix"] >= SPLIT_UNIX]
    log(f"mainshocks M>=5.0: total {len(seqs)} (train {len(tr_seqs)}, test {len(te_seqs)})")
    log(f"sequences with >={MIN_AFT} aftershocks: train "
        f"{sum(1 for s in tr_seqs if s['n_aft'] >= MIN_AFT)}, test "
        f"{sum(1 for s in te_seqs if s['n_aft'] >= MIN_AFT)}")

    t_fit = time.time()
    train_fits = fit_sequences(tr_seqs)
    test_fits = fit_sequences(te_seqs)
    log(f"Omori/b fits: train {len(train_fits)}, test {len(test_fits)} "
        f"({time.time() - t_fit:.1f}s)")

    # runtime guard: project total simulation cost before running the heavy loops
    t0 = time.time()
    _ = ripley_k_time(np.sort(RNG.uniform(0, 10000, 100)), RIPLEY_LAGS, 10000.0)
    per_sim = max(time.time() - t0, 1e-6)
    projected_min = (per_sim * N_SIM + 0.0005 * N_PERM * 2) / 60.0
    n_sim = N_SIM
    guard = dict(projected_sim_minutes=float(projected_min), halved=False,
                 n_sim_used=N_SIM, n_perm_used=N_PERM)
    if projected_min > 90:
        n_sim = N_SIM // 2
        guard.update(halved=True, n_sim_used=n_sim)
    log(f"runtime guard: projected simulation time {projected_min:.2f} min, "
        f"halved={guard['halved']}")

    out = {
        "meta": dict(
            catalog=CATALOG, box=BOX, split_utc=str(SPLIT),
            mainshock_min_mag=MS_MIN_MAG, aftershock_window_days=AFT_WINDOW_D,
            aftershock_radius_km=AFT_RADIUS_KM, min_aftershocks=MIN_AFT,
            b_floor=B_FLOOR, b_bin=B_BIN, min_b_events=MIN_B_EVENTS,
            neighborhood_deg=NEIGH_DEG, seed=20260809, runtime_guard=guard,
            n_events_in_box=int(len(df)),
            n_mainshocks_total=len(seqs), n_mainshocks_train=len(tr_seqs),
            n_mainshocks_test=len(te_seqs),
        ),
        "i": part_i(train_fits, test_fits, log),
        "ii": part_ii(df, log, n_sim=n_sim),
        "iii": part_iii(df, log),
    }
    out["meta"]["runtime_seconds"] = float(time.time() - t_start)

    (HERE / "results_exp_i.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    log(f"wrote results_exp_i.json ({out['meta']['runtime_seconds']:.1f}s total)")


if __name__ == "__main__":
    main()
