"""EXP-M (confirmatory): world transfer test of tectonic-type ETAS parameter pools.

Frozen protocol: PATTERN_PROTOCOL.md section "EXP-M".
  - Data: ComCat global FDSN, M >= 4.5, 1995-01-01 -> now, 13 region boxes.
  - Burn-in 1995-2000 (history only); fit / walk-forward scoring 2000 -> now.
  - Temporal ETAS in EXP-H's form with M0 = 4.5:
        lambda(t) = mu + sum_{t_i<t} K * 10^(alpha (M_i - M0)) * (t - t_i + c)^(-p)
  - TRAIN (7): fit per region; pool by median (GLOBAL = all 7, TYPE = same type).
  - HOLDOUT (6): scored walk-forward only, bits/event vs that region's OWN
    scoring-period Poisson rate (local oracle), under
        (a) GLOBAL pool, (b) TYPE pool, (c) SoCal EXP-H frozen params rescaled to M0=4.5.
    Plus each holdout's own post-hoc fit (DESCRIPTIVE ONLY - the ceiling).
  - Success (frozen): TYPE beats GLOBAL in >= 5 of 6 holdouts (sign test); 4/6 suggestive.

Fitting machinery (_pair_sums / _G_terms / etas_ll / _obj) is taken from exp_h_etas.py
unchanged in form: truncated-kernel (1000 d) multi-start L-BFGS-B in log-parameters,
then an untruncated polish; scoring is untruncated unless the runtime guard fires.

Run: python -u exp_m_world_transfer.py
"""

import csv
import io
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
WORLD = DATA / "comcat_world"
LOGF = DATA / "comcat_world_log.txt"
OUT = HERE / "results_exp_m.json"
EXP_H = HERE / "results_exp_h.json"

FDSN = "https://earthquake.usgs.gov/fdsnws/event/1/query"
USER_AGENT = "quake-replication/exp_m_world_transfer (EXP-M; local research use)"

M0 = 4.5
MINMAG = 4.5
CATALOG_START = pd.Timestamp("1995-01-01", tz="UTC")
SCORE_START = pd.Timestamp("2000-01-01", tz="UTC")
TRUNC_W_DAYS = 1000.0          # truncation used for the (fast) fits, as in EXP-H
SCORE_TRUNC_DAYS = 500.0       # only used if the runtime guard fires
RUNTIME_BUDGET_MIN = 120.0
UNDERPOWERED_N = 300
LN2 = np.log(2.0)
LN10 = np.log(10.0)

ALPHA_BOUNDS = (0.3, 2.5)      # EXP-M protocol bounds
P_BOUNDS = (0.8, 2.0)
MAXITER = 150
REFINE_ITERS = 25

# multi-start seeds: (mu_frac_of_window_rate, K, alpha, c, p) - 4 starts (>= 3 required)
STARTS = [
    (0.50, 0.020, 1.00, 0.010, 1.10),
    (0.20, 0.005, 1.50, 0.050, 1.20),
    (0.80, 0.050, 0.80, 0.005, 1.05),
    (0.35, 0.010, 0.50, 0.020, 1.15),
]

PARAM_KEYS = ["mu", "K", "alpha", "c", "p"]

# ---------------------------------------------------------------- region boxes
# (name, role, type, lat_min, lat_max, lon_min, lon_max)  -- verbatim from PATTERN_PROTOCOL.md
REGIONS = [
    ("Japan",            "train",   "subduction",  30.0,  46.0,  129.0,  147.0),
    ("Chile",            "train",   "subduction", -46.0, -17.0,  -76.0,  -66.0),
    ("Indonesia",        "train",   "subduction", -11.0,   6.0,   95.0,  130.0),
    ("California",       "train",   "transform",   31.5,  42.0, -125.0, -113.0),
    ("Turkey",           "train",   "transform",   35.0,  42.0,   25.0,   45.0),
    ("Himalaya",         "train",   "collision",   25.0,  38.0,   70.0,   98.0),
    ("Iceland",          "train",   "rift",        62.0,  67.0,  -25.0,  -13.0),
    ("Alaska-Aleutians", "holdout", "subduction",  50.0,  62.0, -180.0, -140.0),
    ("Mexico",           "holdout", "subduction",  14.0,  20.0, -105.0,  -92.0),
    ("Philippines",      "holdout", "subduction",   5.0,  20.0,  120.0,  128.0),
    ("Caribbean",        "holdout", "transform",   17.0,  20.0,  -75.0,  -68.0),
    ("Iran",             "holdout", "collision",   26.0,  36.0,   44.0,   62.0),
    ("Greece-Aegean",    "holdout", "rift",        34.0,  41.0,   19.0,   29.0),
]


class RegionDownloadError(RuntimeError):
    """Raised when a region's download fails after retries -> region excluded LOUDLY."""


# ==================================================================== download
EXPECTED_HEAD = "time,latitude,longitude,depth,mag,magType"


def _fetch(url, timeout=240):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.getcode(), resp.read().decode("utf-8", errors="replace")


def fetch_window(reg, t0, t1, log_lines):
    """One [t0, t1) window for one region. Returns (DataFrame, n). Raises after 3 tries."""
    name, _role, _typ, la0, la1, lo0, lo1 = reg
    q = {
        "format": "csv",
        "starttime": t0.strftime("%Y-%m-%dT%H:%M:%S"),
        "endtime": t1.strftime("%Y-%m-%dT%H:%M:%S"),
        "minmagnitude": MINMAG,
        "minlatitude": la0, "maxlatitude": la1,
        "minlongitude": lo0, "maxlongitude": lo1,
        "orderby": "time-asc",
    }
    url = FDSN + "?" + urllib.parse.urlencode(q)
    last_err = None
    for attempt in (1, 2, 3):
        try:
            code, body = _fetch(url)
            if code != 200:
                raise RuntimeError(f"HTTP {code}")
            if body.startswith("time,") is False or not body.startswith(EXPECTED_HEAD):
                raise RuntimeError(f"unexpected header: {body[:120]!r}")
            df = pd.read_csv(io.StringIO(body), quoting=csv.QUOTE_MINIMAL, low_memory=False)
            need = {"time", "latitude", "longitude", "mag", "id", "type"}
            missing = need - set(df.columns)
            if missing:
                raise RuntimeError(f"missing columns {sorted(missing)}")
            n = len(df)
            if n >= 20000:
                raise RuntimeError(f"row count {n} at/over FDSN limit - window too wide")
            log_lines.append(f"{name}  {t0.date()} -> {t1.date()}  HTTP {code}  rows={n}  "
                             f"attempt={attempt}")
            return df, n
        except Exception as exc:                       # noqa: BLE001
            last_err = exc
            log_lines.append(f"{name}  {t0.date()} -> {t1.date()}  FAILED attempt={attempt}: {exc}")
            print(f"    [dl] {name} {t0.date()}->{t1.date()} FAILED attempt={attempt}: {exc}")
            if attempt < 3:
                time.sleep(5.0 * attempt)
    raise RegionDownloadError(f"{name}: window {t0.date()} -> {t1.date()} failed 3x ({last_err})")


def year_windows(t0, t1):
    out, a = [], t0
    while a < t1:
        b = min(a + pd.DateOffset(years=1), t1)
        out.append((a, b))
        a = b
    return out


def download_region(reg, now, log_lines):
    """Download/refresh one region, cached at data/comcat_world/<region>.csv."""
    name, _role, _typ, la0, la1, lo0, lo1 = reg
    WORLD.mkdir(parents=True, exist_ok=True)
    cache = WORLD / f"{name}.csv"
    start = CATALOG_START
    cached = None
    if cache.exists():
        cached = pd.read_csv(cache, low_memory=False)
        cached["time"] = pd.to_datetime(cached["time"], utc=True, format="mixed")
        if len(cached):
            start = max(cached["time"].max() - pd.Timedelta(days=30), CATALOG_START)
            log_lines.append(f"# {name}: cache hit {len(cached)} rows, last "
                             f"{cached['time'].max().isoformat()}; refetch from {start.date()}")
            print(f"  [cache] {name}: {len(cached)} rows, refetch from {start.date()}")
        else:
            cached = None

    wins = year_windows(start, now)
    frames, counts = [], []
    for a, b in wins:
        df, n = fetch_window(reg, a, b, log_lines)
        frames.append(df)
        counts.append(n)
        time.sleep(0.25)
    fresh = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["time", "id"])
    # counted invariant: per-window row counts must sum to merged rows BEFORE dedup
    assert sum(counts) == len(fresh), f"{name}: window sum {sum(counts)} != merged {len(fresh)}"
    log_lines.append(f"# {name}: invariant OK sum(window rows)={sum(counts)} == "
                     f"merged before dedup={len(fresh)} over {len(wins)} windows")

    if len(fresh):
        fresh["time"] = pd.to_datetime(fresh["time"], utc=True, format="mixed")
    if cached is not None and len(fresh):
        keep = [c for c in cached.columns if c in fresh.columns]
        merged = pd.concat([cached[keep], fresh[keep]], ignore_index=True)
    elif cached is not None:
        merged = cached
    else:
        merged = fresh

    n_pre = len(merged)
    merged = merged.drop_duplicates(subset="id", keep="last")
    n_dupes = n_pre - len(merged)
    merged = merged.sort_values("time", kind="mergesort").reset_index(drop=True)
    merged = merged[(merged.time >= CATALOG_START) & (merged.time <= now)].reset_index(drop=True)

    # ComCat's minmagnitude filter occasionally leaks a row whose PREFERRED magnitude is
    # below the floor (e.g. nc21364840, M3.38 ml, an automatic solution inside the
    # Indonesia box) and rows with a null magnitude. Drop them, but COUNT and LOG them -
    # never silently.
    n_sub_floor = int((merged.mag < MINMAG - 1e-9).sum()) if len(merged) else 0
    n_nan_mag = int(merged.mag.isna().sum()) if len(merged) else 0
    if n_sub_floor or n_nan_mag:
        bad_ids = merged.loc[merged.mag.isna() | (merged.mag < MINMAG - 1e-9), "id"].tolist()
        log_lines.append(f"# {name}: DROPPED {n_sub_floor} sub-floor + {n_nan_mag} null-mag "
                         f"rows leaked past the FDSN minmagnitude filter: {bad_ids}")
        print(f"  [clean] {name}: dropped {n_sub_floor} sub-floor / {n_nan_mag} null-mag rows "
              f"{bad_ids}")
        merged = merged[merged.mag >= MINMAG - 1e-9].reset_index(drop=True)

    assert merged["id"].is_unique, f"{name}: duplicate ids survived dedup"
    assert merged["time"].is_monotonic_increasing, f"{name}: times not monotonic"
    if len(merged):
        assert merged.mag.min() >= MINMAG - 1e-9, f"{name}: magnitude floor violated"
        assert merged.latitude.between(la0, la1).all(), f"{name}: latitude outside box"
        assert merged.longitude.between(lo0, lo1).all(), f"{name}: longitude outside box"

    merged.to_csv(cache, index=False)
    log_lines.append(f"# {name}: before dedup={n_pre}, dup ids removed={n_dupes}, "
                     f"final={len(merged)}, range "
                     f"{merged.time.min().isoformat() if len(merged) else 'NA'} -> "
                     f"{merged.time.max().isoformat() if len(merged) else 'NA'}")
    meta = {"windows": len(wins), "window_counts": counts, "rows_before_dedup": int(n_pre),
            "duplicate_ids_removed": int(n_dupes), "cache_file": str(cache),
            "n_events": int(len(merged))}
    return merged, meta


# ==================================================================== ETAS core
# (identical in form to exp_h_etas.py)
def _pair_sums(t, w, m, c, p, tgt_lo, tgt_hi, W, want_grad, chunk):
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
    return G, dGdc, -dGds * live


def etas_ll(theta, t_all, m_all, tgt_lo, tgt_hi, T0, T1, W, want_grad=True, chunk=512):
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


def aki_b(mags, mc, binw=0.1):
    m = np.asarray(mags, float)
    m = m[m >= mc - 1e-9]
    if len(m) < 20 or m.mean() <= (mc - binw / 2.0):
        return None
    return float(1.0 / (LN10 * (m.mean() - (mc - binw / 2.0))))


# ==================================================================== fit / score
def fit_region(t_all, m_all, T0_fit, T1_fit, label):
    """Multi-start truncated-kernel MLE then untruncated polish. Returns dict."""
    lo = int(np.searchsorted(t_all, T0_fit, side="left"))
    hi = int(np.searchsorted(t_all, T1_fit, side="right"))
    n_fit = hi - lo
    D = T1_fit - T0_fit
    rate = max(n_fit / D, 1e-8)
    args_tr = (t_all, m_all, lo, hi, T0_fit, T1_fit, TRUNC_W_DAYS, True, 512)
    lb = np.log([1e-8, 1e-8, ALPHA_BOUNDS[0], 1e-6, P_BOUNDS[0]])
    ub = np.log([1e3, 1e2, ALPHA_BOUNDS[1], 1e1, P_BOUNDS[1]])
    bounds = list(zip(lb, ub))

    starts_out, best = [], None
    for si, (mu_f, K, al, c, p) in enumerate(STARTS):
        x0 = np.log(np.clip([mu_f * rate, K, al, c, p], np.exp(lb) * 1.001, np.exp(ub) * 0.999))
        ts_ = time.time()
        r = minimize(_obj, x0, args=args_tr, jac=True, method="L-BFGS-B", bounds=bounds,
                     options={"maxiter": MAXITER, "maxfun": MAXITER * 2, "ftol": 1e-12,
                              "gtol": 1e-8})
        th = np.exp(r.x)
        rec = {"start": dict(zip(PARAM_KEYS, [mu_f * rate, K, al, c, p])),
               "LL": float(-r.fun), "params": dict(zip(PARAM_KEYS, map(float, th))),
               "nit": int(r.nit), "success": bool(r.success), "seconds": round(time.time() - ts_, 1)}
        starts_out.append(rec)
        print(f"    [fit:{label}] start {si+1}: LL={rec['LL']:.2f} mu={th[0]:.4f} K={th[1]:.5f} "
              f"alpha={th[2]:.3f} c={th[3]:.5f} p={th[4]:.4f} ({rec['seconds']}s)")
        if best is None or rec["LL"] > best["LL"]:
            best = rec

    theta_trunc = np.array([best["params"][k] for k in PARAM_KEYS])
    LL_tr, _ = etas_ll(theta_trunc, *args_tr[:6], TRUNC_W_DAYS, False, 1024)
    LL_full, _ = etas_ll(theta_trunc, *args_tr[:6], np.inf, False, 1024)
    args_full = (*args_tr[:6], np.inf, True, 1024)
    rf = minimize(_obj, np.log(theta_trunc), args=args_full, jac=True, method="L-BFGS-B",
                  bounds=bounds, options={"maxiter": REFINE_ITERS, "maxfun": REFINE_ITERS * 2,
                                          "ftol": 1e-12})
    theta = np.exp(rf.x)
    LL_ref = float(-rf.fun)
    kept_trunc = False
    if not np.isfinite(LL_ref) or LL_ref < LL_full:
        theta, LL_ref, kept_trunc = theta_trunc, LL_full, True
    print(f"    [fit:{label}] FINAL LL={LL_ref:.2f} mu={theta[0]:.4f} K={theta[1]:.5f} "
          f"alpha={theta[2]:.3f} c={theta[3]:.5f} p={theta[4]:.4f}"
          f"{'  (kept truncated optimum)' if kept_trunc else ''}")
    return {"params": dict(zip(PARAM_KEYS, map(float, theta))), "LL": LL_ref,
            "n_target_events": int(n_fit), "fit_window_days": float(D),
            "window_rate_per_day": float(rate),
            "truncated_fit_LL": LL_tr, "untruncated_LL_at_truncated_optimum": LL_full,
            "delta_LL_truncation": float(LL_full - LL_tr),
            "refine_improved": bool(not kept_trunc), "starts": starts_out}


def score_components(t_all, m_all, sc_lo, sc_hi, T0, T1, alpha, c, p, W):
    """S_j for scored targets and sum_i w_i G_i - everything that does not depend on (mu, K)."""
    w = np.exp(alpha * LN10 * m_all)
    S, _, _, _ = _pair_sums(t_all, w, m_all, c, p, sc_lo, sc_hi, W, False, 1024)
    src = t_all <= T1
    G, _, _ = _G_terms(t_all[src], T0, T1, c, p, W, False)
    return S, float((w[src] * G).sum())


def ll_from_components(S, wGsum, mu, K, T0, T1):
    lam = mu + K * S
    if not np.all(np.isfinite(lam)) or np.any(lam <= 0):
        return None, None
    return float(np.log(lam).sum() - (mu * (T1 - T0) + K * wGsum)), lam


# ==================================================================== main
def main():
    t_start = time.time()
    now = pd.Timestamp.utcnow().floor("h")
    res = {"experiment": "EXP-M", "protocol": "PATTERN_PROTOCOL.md :: EXP-M",
           "state_class": "first-run",
           "run_utc": pd.Timestamp.utcnow().isoformat(),
           "download_cutoff_utc": now.isoformat(),
           "settings": {"M0": M0, "catalog_start": str(CATALOG_START),
                        "burn_in": "1995-01-01 -> 2000-01-01 (history only)",
                        "score_window_start": str(SCORE_START),
                        "score_window_end": now.isoformat(),
                        "trunc_days_fit": TRUNC_W_DAYS,
                        "alpha_bounds": list(ALPHA_BOUNDS), "p_bounds": list(P_BOUNDS),
                        "n_starts": len(STARTS)}}

    # ---------------- step 1: download ----------------
    print("=" * 78)
    print(f"[EXP-M] downloading 13 regions, M>={MINMAG}, {CATALOG_START.date()} -> {now}")
    log_lines = [f"# EXP-M ComCat retrieval log  run_utc={pd.Timestamp.utcnow().isoformat()}",
                 f"# now={now.isoformat()}  M>={MINMAG}  start={CATALOG_START.isoformat()}"]
    cats, region_meta, excluded = {}, {}, {}
    for reg in REGIONS:
        name, role, typ = reg[0], reg[1], reg[2]
        print(f"  [region] {name} ({role}/{typ}) box lat[{reg[3]},{reg[4]}] lon[{reg[5]},{reg[6]}]")
        try:
            df, meta = download_region(reg, now, log_lines)
        except (RegionDownloadError, urllib.error.URLError, AssertionError) as exc:
            excluded[name] = {"role": role, "type": typ, "reason": str(exc)}
            log_lines.append(f"# {name}: EXCLUDED - {exc}")
            print(f"  !! EXCLUDED {name}: {exc}")
            continue
        meta.update(role=role, type=typ,
                    box={"lat": [reg[3], reg[4]], "lon": [reg[5], reg[6]]},
                    t_first=str(df.time.min()) if len(df) else None,
                    t_last=str(df.time.max()) if len(df) else None,
                    underpowered=bool(len(df) < UNDERPOWERED_N))
        cats[name] = df
        region_meta[name] = meta
        print(f"  [region] {name}: {len(df)} events"
              f"{'   *** UNDERPOWERED (<300) ***' if meta['underpowered'] else ''}")
    LOGF.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print(f"[log] {LOGF}")
    res["regions"] = region_meta
    res["excluded_regions"] = excluded
    res["download_log"] = str(LOGF)

    # numeric arrays per region
    arr = {}
    t_score0 = (SCORE_START - pd.Timestamp(0, tz="UTC")).total_seconds() / 86400.0
    t_end = (now - pd.Timestamp(0, tz="UTC")).total_seconds() / 86400.0
    for name, df in cats.items():
        d = df.sort_values("time").reset_index(drop=True)
        t = (d.time - pd.Timestamp(0, tz="UTC")).dt.total_seconds().to_numpy() / 86400.0
        m = d.mag.to_numpy(float) - M0
        arr[name] = (t, m, d.mag.to_numpy(float))
        region_meta[name]["n_burn_in_1995_2000"] = int((t < t_score0).sum())
        region_meta[name]["n_scoring_2000_now"] = int((t >= t_score0).sum())
        region_meta[name]["b_value_aki"] = aki_b(d.mag.to_numpy(float), M0)

    # ---------------- step 2: fit train regions ----------------
    print("\n" + "=" * 78)
    print("[EXP-M] fitting TRAIN regions (2000 -> now, 1995-2000 history-only burn-in)")
    fits = {}
    for reg in REGIONS:
        name, role, typ = reg[0], reg[1], reg[2]
        if role != "train" or name not in arr:
            continue
        t, m, mags = arr[name]
        if region_meta[name]["n_scoring_2000_now"] < 30:
            fits[name] = {"type": typ, "skipped": True,
                          "reason": f"only {region_meta[name]['n_scoring_2000_now']} target events"}
            print(f"  [fit] {name}: SKIPPED (too few target events)")
            continue
        print(f"  [fit] {name} ({typ}) n_total={len(t)} n_target={region_meta[name]['n_scoring_2000_now']}")
        f = fit_region(t, m, t_score0, t_end, name)
        f.update(type=typ, role="train", n_events_total=int(len(t)),
                 b_value_aki=region_meta[name]["b_value_aki"])
        fits[name] = f
    res["train_fits"] = fits

    # ---------------- step 3: pools ----------------
    good = {k: v for k, v in fits.items() if not v.get("skipped")}
    def _median_pool(names):
        return {k: float(np.median([good[n]["params"][k] for n in names])) for k in PARAM_KEYS}

    global_pool = _median_pool(list(good))
    type_names = {}
    for n, v in good.items():
        type_names.setdefault(v["type"], []).append(n)
    type_pools = {ty: {"params": _median_pool(ns), "n_train_regions": len(ns), "regions": ns}
                  for ty, ns in type_names.items()}

    # --- SoCal EXP-H frozen params, rescaled from M0=2.5 to M0=4.5 ---------------
    # Intensity:  lambda = mu + sum_i K * 10^(alpha (M_i - M0)) * (t-t_i+c)^-p
    # Changing the reference magnitude M0_old -> M0_new must leave each source's
    # productivity for a GIVEN magnitude M unchanged:
    #     K_old * 10^(alpha (M - M0_old))
    #   = K_old * 10^(alpha (M0_new - M0_old)) * 10^(alpha (M - M0_new))
    #   => K_new = K_old * 10^(alpha * (M0_new - M0_old))          [exact invariance]
    # NOTE: the task brief proposed K' = K * 10^(alpha*(M0_old - M0_new)); the algebra
    # above gives the OPPOSITE sign. We use the derived exact-invariance form as primary
    # and record the brief's form as `K_brief_sign_variant` for transparency.
    # mu is a background RATE of M >= M0 events, so it rescales by Gutenberg-Richter:
    #     mu_new = mu_old * 10^(-b (M0_new - M0_old)),  b = SoCal Aki b from EXP-H.
    # c, p, alpha are magnitude-reference invariant.
    with open(EXP_H) as fh:
        h = json.load(fh)
    hp = h["train_fit"]["frozen_params"]
    b_socal = float(h["train_fit"]["b_value_train_aki"])
    a_s, M0_old = float(hp["alpha"]), float(hp["M0"])
    K_new = float(hp["K"]) * 10.0 ** (a_s * (M0 - M0_old))
    K_brief = float(hp["K"]) * 10.0 ** (a_s * (M0_old - M0))
    mu_new = float(hp["mu"]) * 10.0 ** (-b_socal * (M0 - M0_old))
    socal = {"mu": mu_new, "K": K_new, "alpha": a_s, "c": float(hp["c"]), "p": float(hp["p"])}
    assert np.isfinite(K_new) and K_new > 0 and np.isfinite(mu_new) and mu_new > 0

    socal_block = {
        "source": "results_exp_h.json :: train_fit.frozen_params (SoCal, M0=2.5)",
        "original": dict(hp), "b_value_socal": b_socal,
        "K_rescale_derivation": ("K*10^(alpha(M-M0_old)) = K*10^(alpha(M0_new-M0_old))"
                                 "*10^(alpha(M-M0_new)) => K_new = K*10^(alpha*(M0_new-M0_old))"),
        "K_rescale_factor": float(10.0 ** (a_s * (M0 - M0_old))),
        "mu_rescale_derivation": "GR: rate(M>=M0_new) = rate(M>=M0_old)*10^(-b(M0_new-M0_old))",
        "mu_rescale_factor": float(10.0 ** (-b_socal * (M0 - M0_old))),
        "K_brief_sign_variant": K_brief,
        "note_on_brief": ("the task brief wrote K'=K*10^(alpha*(M0_old-M0_new)); the intensity "
                          "algebra gives the opposite sign. Exact-invariance form used."),
        "params_rescaled": socal,
        "caveat": ("SoCal params were estimated with M>=2.5 SOURCES; applied to an M>=4.5 "
                   "catalog the triggering sum omits all 2.5<=M<4.5 sources, so the transferred "
                   "intensity is a lower bound on what those params imply."),
    }
    res["pools"] = {"GLOBAL": {"params": global_pool, "n_train_regions": len(good),
                               "regions": list(good)},
                    "TYPE": type_pools, "SOCAL_EXP_H_rescaled": socal_block}
    print("\n[pools] GLOBAL  " + "  ".join(f"{k}={global_pool[k]:.5g}" for k in PARAM_KEYS))
    for ty, tp in type_pools.items():
        print(f"[pools] TYPE {ty:<11s} (n={tp['n_train_regions']}) " +
              "  ".join(f"{k}={tp['params'][k]:.5g}" for k in PARAM_KEYS))
    print("[pools] SoCal->M0=4.5  " + "  ".join(f"{k}={socal[k]:.5g}" for k in PARAM_KEYS) +
          f"   (K x{socal_block['K_rescale_factor']:.4g}, mu x{socal_block['mu_rescale_factor']:.4g})")

    # ---------------- step 4: runtime projection for scoring ----------------
    holdouts = [r[0] for r in REGIONS if r[1] == "holdout" and r[0] in arr]
    n_pairs = 0
    for name in holdouts:
        t, _, _ = arr[name]
        ns = region_meta[name]["n_scoring_2000_now"]
        n_pairs += ns * len(t)
    n_pass = 4                       # GLOBAL, TYPE, SoCal, own-fit ceiling
    proj_s = n_pairs * n_pass * 3e-8 * 6   # ~6 array ops per pair, ~30 ns each (rough)
    use_trunc = proj_s / 60.0 > RUNTIME_BUDGET_MIN
    W_score = SCORE_TRUNC_DAYS if use_trunc else np.inf
    print(f"\n[runtime] holdout scoring pairs ~{n_pairs:,} x {n_pass} passes; "
          f"projected ~{proj_s/60:.1f} min -> "
          f"{'500-day TRUNCATED scoring' if use_trunc else 'UNTRUNCATED scoring'}")
    res["runtime_guard"] = {"projected_scoring_minutes": round(proj_s / 60.0, 2),
                            "budget_minutes": RUNTIME_BUDGET_MIN,
                            "scoring_truncation_days": (SCORE_TRUNC_DAYS if use_trunc else None),
                            "scoring_untruncated": bool(not use_trunc)}

    # ---------------- step 5: score holdouts ----------------
    print("\n" + "=" * 78)
    print("[EXP-M] scoring HOLDOUT regions (walk-forward 2000 -> now)")
    scores, trunc_probe = {}, None
    for reg in REGIONS:
        name, role, typ = reg[0], reg[1], reg[2]
        if role != "holdout" or name not in arr:
            continue
        t, m, mags = arr[name]
        sc_lo = int(np.searchsorted(t, t_score0, side="left"))
        sc_hi = len(t)
        n_sc = sc_hi - sc_lo
        D = t_end - t_score0
        entry = {"type": typ, "n_events_total": int(len(t)),
                 "n_burn_in": int(sc_lo), "n_scored": int(n_sc),
                 "scoring_days": float(D), "b_value_aki": region_meta[name]["b_value_aki"],
                 "underpowered": region_meta[name]["underpowered"]}
        if n_sc < 10:
            entry["skipped"] = f"only {n_sc} scored events"
            scores[name] = entry
            print(f"  [score] {name}: SKIPPED ({n_sc} events)")
            continue
        rate_oracle = n_sc / D
        LL_pois = n_sc * np.log(rate_oracle) - rate_oracle * D
        entry["poisson_local_oracle"] = {"rate_per_day": float(rate_oracle), "LL": float(LL_pois)}

        sources = {"GLOBAL": global_pool,
                   "TYPE": type_pools[typ]["params"] if typ in type_pools else None,
                   "SOCAL_rescaled": socal}
        # descriptive-only own fit (the ceiling)
        own = fit_region(t, m, t_score0, t_end, f"{name}(own,POST-HOC)")
        entry["own_fit_POST_HOC_DESCRIPTIVE_ONLY"] = {
            "params": own["params"], "LL_in_sample": own["LL"],
            "label": "post-hoc in-sample fit on the holdout itself - NOT a forecast; ceiling only"}
        sources["OWN_FIT_ceiling"] = own["params"]

        entry["bits_per_event"] = {}
        for sname, par in sources.items():
            if par is None:
                entry["bits_per_event"][sname] = None
                continue
            S, wG = score_components(t, m, sc_lo, sc_hi, t_score0, t_end,
                                     par["alpha"], par["c"], par["p"], W_score)
            LL, lam = ll_from_components(S, wG, par["mu"], par["K"], t_score0, t_end)
            if LL is None:
                entry["bits_per_event"][sname] = None
                continue
            bits = (LL - LL_pois) / n_sc / LN2
            entry["bits_per_event"][sname] = {
                "LL": LL, "bits_per_event_vs_local_poisson": float(bits),
                "params": dict(par),
                "triggered_fraction_mean_over_events": float(np.mean((lam - par["mu"]) / lam))}
            print(f"  [score] {name:<17s} {sname:<16s} LL={LL:11.1f}  "
                  f"bits/event={bits:+.4f}")
            # one-off truncation-effect probe on a mid-size region under the TYPE pool
            if (trunc_probe is None and sname == "TYPE" and 300 <= n_sc <= 3000
                    and not use_trunc):
                S2, wG2 = score_components(t, m, sc_lo, sc_hi, t_score0, t_end,
                                           par["alpha"], par["c"], par["p"], SCORE_TRUNC_DAYS)
                LL2, _ = ll_from_components(S2, wG2, par["mu"], par["K"], t_score0, t_end)
                trunc_probe = {"region": name, "n_scored": int(n_sc), "LL_untruncated": LL,
                               "LL_trunc_500d": LL2, "delta_LL": float(LL2 - LL),
                               "delta_bits_per_event": float((LL2 - LL) / n_sc / LN2),
                               "note": "diagnostic only; scoring above is untruncated"}
                print(f"  [probe] 500-d truncation on {name}: dLL={LL2-LL:+.2f} "
                      f"({(LL2-LL)/n_sc/LN2:+.5f} bits/event)")
        scores[name] = entry
    res["holdout_scores"] = scores
    res["scoring_truncation_probe"] = trunc_probe

    # ---------------- step 6: sign test ----------------
    wins, ties, losses, usable = 0, 0, 0, []
    for name, e in scores.items():
        bp = e.get("bits_per_event") or {}
        g, ty = bp.get("GLOBAL"), bp.get("TYPE")
        if not g or not ty:
            continue
        usable.append(name)
        dg = ty["bits_per_event_vs_local_poisson"] - g["bits_per_event_vs_local_poisson"]
        if dg > 0:
            wins += 1
        elif dg < 0:
            losses += 1
        else:
            ties += 1
        e["TYPE_minus_GLOBAL_bits"] = float(dg)
    n_used = len(usable)
    verdict = ("SUCCESS" if wins >= 5 else "SUGGESTIVE" if wins == 4 else "FAIL")
    if n_used < 6:
        verdict += f" (on n={n_used} holdouts, not the frozen 6)"
    res["sign_test"] = {"rule": "TYPE pool beats GLOBAL pool in >=5 of 6 holdouts "
                                "(4/6 = suggestive only)",
                        "n_holdouts_scored": n_used, "holdouts_used": usable,
                        "TYPE_wins": wins, "GLOBAL_wins": losses, "ties": ties,
                        "verdict": verdict}

    res["flags"] = {
        "underpowered_regions": [n for n, m_ in region_meta.items() if m_["underpowered"]],
        "excluded_regions": list(excluded),
        "type_pools_with_n1": [ty for ty, tp in type_pools.items() if tp["n_train_regions"] == 1],
        "own_fit_is_post_hoc_descriptive": True,
        "socal_params_from_M2.5_sources_applied_to_M4.5_catalog": True,
        "scoring_untruncated": bool(not use_trunc),
        "mu_pooled_across_regions_is_a_known_weakness": (
            "mu is a per-region background RATE; medianing it across regions of very "
            "different seismicity is what the frozen protocol specifies, and it dominates "
            "the transferred LL. Recorded as-is."),
    }
    res["runtime_minutes"] = round((time.time() - t_start) / 60.0, 2)
    OUT.write_text(json.dumps(res, indent=2))

    # ---------------- console summary ----------------
    print("\n" + "=" * 78)
    print("EXP-M  |  per-region ETAS fits (M0=4.5, fit 2000->now, burn-in 1995-2000)")
    print(f"{'region':<18s}{'role':<9s}{'type':<12s}{'n_tot':>7s}{'n_fit':>7s}"
          f"{'mu':>9s}{'K':>9s}{'alpha':>7s}{'c':>9s}{'p':>7s}{'b':>6s}{'LL':>11s}")
    for reg in REGIONS:
        name, role, typ = reg[0], reg[1], reg[2]
        if name in excluded:
            print(f"{name:<18s}{role:<9s}{typ:<12s}  EXCLUDED: {excluded[name]['reason'][:40]}")
            continue
        if name not in region_meta:
            continue
        f = fits.get(name) or (scores.get(name, {}).get("own_fit_POST_HOC_DESCRIPTIVE_ONLY"))
        nt = region_meta[name]["n_events"]
        bv = region_meta[name]["b_value_aki"]
        if not f or f.get("skipped"):
            print(f"{name:<18s}{role:<9s}{typ:<12s}{nt:>7d}   (no fit)")
            continue
        pr = f["params"]
        ll = f.get("LL", f.get("LL_in_sample"))
        nf = f.get("n_target_events", region_meta[name]["n_scoring_2000_now"])
        tag = "" if role == "train" else "  [post-hoc]"
        print(f"{name:<18s}{role:<9s}{typ:<12s}{nt:>7d}{nf:>7d}{pr['mu']:>9.4f}{pr['K']:>9.4f}"
              f"{pr['alpha']:>7.3f}{pr['c']:>9.5f}{pr['p']:>7.3f}"
              f"{(bv if bv else float('nan')):>6.2f}{ll:>11.1f}{tag}")

    print("\nPOOLS (per-parameter medians)")
    print(f"{'pool':<22s}{'n':>3s}{'mu':>10s}{'K':>10s}{'alpha':>8s}{'c':>10s}{'p':>8s}")
    print(f"{'GLOBAL':<22s}{len(good):>3d}" + "".join(
        f"{global_pool[k]:>10.5f}" if k in ('mu', 'K', 'c') else f"{global_pool[k]:>8.3f}"
        for k in PARAM_KEYS))
    for ty, tp in sorted(type_pools.items()):
        p_ = tp["params"]
        print(f"{'TYPE:'+ty:<22s}{tp['n_train_regions']:>3d}" + "".join(
            f"{p_[k]:>10.5f}" if k in ('mu', 'K', 'c') else f"{p_[k]:>8.3f}" for k in PARAM_KEYS))
    print(f"{'SoCal(EXP-H)->4.5':<22s}{1:>3d}" + "".join(
        f"{socal[k]:>10.5f}" if k in ('mu', 'K', 'c') else f"{socal[k]:>8.3f}" for k in PARAM_KEYS))

    print("\nHOLDOUT bits/event vs OWN-period Poisson (local oracle)")
    print(f"{'region':<18s}{'type':<12s}{'n':>6s}{'GLOBAL':>10s}{'TYPE':>10s}"
          f"{'SoCal':>10s}{'own(ceil)':>11s}{'T-G':>9s}")
    for name, e in scores.items():
        bp = e.get("bits_per_event") or {}
        def g(k):
            v = bp.get(k)
            return f"{v['bits_per_event_vs_local_poisson']:>10.3f}" if v else f"{'--':>10s}"
        d = e.get("TYPE_minus_GLOBAL_bits")
        print(f"{name:<18s}{e['type']:<12s}{e.get('n_scored',0):>6d}{g('GLOBAL')}{g('TYPE')}"
              f"{g('SOCAL_rescaled')}{g('OWN_FIT_ceiling'):>11s}"
              f"{(f'{d:+.3f}' if d is not None else '--'):>9s}")

    st = res["sign_test"]
    print(f"\nSIGN TEST  TYPE > GLOBAL in {st['TYPE_wins']}/{st['n_holdouts_scored']} holdouts "
          f"(GLOBAL wins {st['GLOBAL_wins']}, ties {st['ties']})  ->  {st['verdict']}")
    if res["flags"]["underpowered_regions"]:
        print(f"FLAG underpowered (<{UNDERPOWERED_N} events): "
              f"{', '.join(res['flags']['underpowered_regions'])}")
    if excluded:
        print(f"FLAG EXCLUDED regions: {', '.join(excluded)}")
    print(f"FLAG type pools with n=1: {', '.join(res['flags']['type_pools_with_n1']) or 'none'}")
    print(f"runtime {res['runtime_minutes']} min  ->  {OUT.name}")
    print("=" * 78)


if __name__ == "__main__":
    main()
