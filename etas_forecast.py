"""EXP-L (application, not a hypothesis test): live 7-day forecast from frozen ETAS.

Frozen protocol: PATTERN_PROTOCOL.md section "EXP-L".

  - Parameters are TAKEN, NOT FIT: results_exp_h.json -> train_fit.frozen_params
    (mu, K, alpha, c, p, M0=2.5) and train_fit.b_value_train_aki (b = 1.0654...).
    Those were fit on SCSN 1982-2010 and scored once, walk-forward, on 2010-2018.
    NOTHING here is refit. That is the whole point of the exercise.

  - Intensity convention (identical to exp_h_etas.py, t in DAYS):
        lambda(t) = mu + sum_{t_i < t} K * 10^(alpha * (M_i - M0)) * (t - t_i + c)^(-p)

  - Catalog: USGS ComCat FDSN event service, format=csv, M >= 2.5,
    lat [31.5, 38.0], lon [-122.0, -113.5], 2010-01-01 -> now, paged in 1-year
    windows, every window's HTTP status and row count checked and logged.

  - Forecast horizon: 7 days from --now (default: current UTC).

Run:  python -u etas_forecast.py [--now 2026-08-09T12:00:00Z] [--refresh]
"""

import argparse
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

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
CACHE = DATA / "comcat_socal_m25.csv"
EXP_H = HERE / "results_exp_h.json"

# ---- frozen geometry / thresholds (same box as EXP-H) -------------------------
LAT_MIN, LAT_MAX = 31.5, 38.0
LON_MIN, LON_MAX = -122.0, -113.5
MINMAG = 2.5
CATALOG_START = "2010-01-01T00:00:00"

# ---- frozen forecast settings ------------------------------------------------
HORIZON_DAYS = 7.0
TRUNC_DAYS = 1000.0        # history truncation for lambda(now) and for seeding
N_SIMS = 1000
SEED = 20260813
MAG_CAP = 8.0              # GR upper truncation
MAX_EVENTS_PER_SIM = 200000
ALARM_WINDOW_DAYS = 7.0    # EXP-I(ii) alarm: any M >= 5 in the last 7 days
ALARM_MAG = 5.0

FDSN = "https://earthquake.usgs.gov/fdsnws/event/1/query"
USER_AGENT = "quake-replication/etas_forecast (EXP-L; contact: local research use)"
LN10 = np.log(10.0)


# ============================================================ frozen parameters
def load_frozen():
    with open(EXP_H) as fh:
        h = json.load(fh)
    fp = h["train_fit"]["frozen_params"]
    par = {
        "mu": float(fp["mu"]),
        "K": float(fp["K"]),
        "alpha": float(fp["alpha"]),
        "c": float(fp["c"]),
        "p": float(fp["p"]),
        "M0": float(fp["M0"]),
        "b": float(h["train_fit"]["b_value_train_aki"]),
    }
    assert par["p"] > 1.0, "kernel integral formula assumes p > 1"
    assert par["b"] > par["alpha"], "GR productivity expectation requires b > alpha"
    return par, h


# ==================================================================== download
def _fetch(url, timeout=180):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        code = resp.getcode()
        body = resp.read().decode("utf-8", errors="replace")
    return code, body


EXPECTED_HEAD = "time,latitude,longitude,depth,mag,magType"


def fetch_window(t0, t1, log_lines, attempt_note=""):
    """One [t0, t1) window. Returns DataFrame. Raises on unrecoverable failure."""
    q = {
        "format": "csv",
        "starttime": t0.strftime("%Y-%m-%dT%H:%M:%S"),
        "endtime": t1.strftime("%Y-%m-%dT%H:%M:%S"),
        "minmagnitude": MINMAG,
        "minlatitude": LAT_MIN,
        "maxlatitude": LAT_MAX,
        "minlongitude": LON_MIN,
        "maxlongitude": LON_MAX,
        "orderby": "time-asc",
    }
    url = FDSN + "?" + urllib.parse.urlencode(q)
    last_err = None
    for attempt in (1, 2):                      # one retry, per protocol
        try:
            code, body = _fetch(url)
            if code != 200:
                raise RuntimeError(f"HTTP {code}")
            if not body.startswith(EXPECTED_HEAD):
                raise RuntimeError(f"unexpected header: {body[:120]!r}")
            df = pd.read_csv(io.StringIO(body), quoting=csv.QUOTE_MINIMAL,
                             low_memory=False)
            need = {"time", "latitude", "longitude", "mag", "id", "type"}
            missing = need - set(df.columns)
            if missing:
                raise RuntimeError(f"missing columns {sorted(missing)}")
            n = len(df)
            if n >= 20000:
                raise RuntimeError(f"row count {n} at/over FDSN limit - window too wide")
            log_lines.append(
                f"{t0.date()} -> {t1.date()}  HTTP {code}  rows={n}  attempt={attempt}"
                f"{attempt_note}")
            print(f"  [dl] {t0.date()} -> {t1.date()}  HTTP {code}  rows={n}"
                  f"{'  (retry)' if attempt == 2 else ''}")
            return df, n
        except Exception as exc:                # noqa: BLE001 - retry once then abort
            last_err = exc
            log_lines.append(f"{t0.date()} -> {t1.date()}  FAILED attempt={attempt}: {exc}")
            print(f"  [dl] {t0.date()} -> {t1.date()}  FAILED attempt={attempt}: {exc}")
            if attempt == 1:
                time.sleep(5.0)
    raise SystemExit(
        f"ABORT: ComCat window {t0.date()} -> {t1.date()} failed twice ({last_err}). "
        "Refusing to continue with a truncated catalog. Re-run when the service "
        "is reachable; the cache on disk is unchanged.")


def year_windows(t0, t1):
    out, a = [], t0
    while a < t1:
        b = min(a + pd.DateOffset(years=1), t1)
        out.append((a, b))
        a = b
    return out


def download_catalog(now, refresh=False):
    """Merged catalog 2010-01-01 -> now, using data/comcat_socal_m25.csv as cache."""
    DATA.mkdir(exist_ok=True)
    log_lines = [f"# ComCat retrieval log  run_utc={pd.Timestamp.utcnow().isoformat()}",
                 f"# now={now.isoformat()}  box lat[{LAT_MIN},{LAT_MAX}] "
                 f"lon[{LON_MIN},{LON_MAX}]  M>={MINMAG}"]
    start = pd.Timestamp(CATALOG_START, tz="UTC")

    cached = None
    if CACHE.exists() and not refresh:
        cached = pd.read_csv(CACHE, low_memory=False)
        cached["time"] = pd.to_datetime(cached["time"], utc=True, format="mixed")
        if len(cached):
            refetch_from = cached["time"].max() - pd.Timedelta(days=30)
            refetch_from = max(refetch_from, start)
            log_lines.append(f"# cache hit: {len(cached)} rows, last event "
                             f"{cached['time'].max().isoformat()}; re-downloading from "
                             f"{refetch_from.isoformat()} (last event - 30 d)")
            print(f"[cache] {len(cached)} rows, last event {cached['time'].max()}; "
                  f"refetch from {refetch_from.date()}")
            start = refetch_from
        else:
            cached = None

    wins = year_windows(start, now)
    log_lines.append(f"# {len(wins)} windows")
    frames, counts = [], []
    for a, b in wins:
        df, n = fetch_window(a, b, log_lines)
        frames.append(df)
        counts.append(n)

    fresh = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    # counted invariant: window counts must equal merged rows before any dedup
    assert sum(counts) == len(fresh), \
        f"window count sum {sum(counts)} != merged rows {len(fresh)}"
    log_lines.append(f"# invariant OK: sum(window rows)={sum(counts)} == "
                     f"merged rows before dedup={len(fresh)}")

    if cached is not None:
        fresh["time"] = pd.to_datetime(fresh["time"], utc=True, format="mixed")
        keep = [c for c in cached.columns if c in fresh.columns]
        merged = pd.concat([cached[keep], fresh[keep]], ignore_index=True)
    else:
        merged = fresh
        merged["time"] = pd.to_datetime(merged["time"], utc=True, format="mixed")

    n_pre = len(merged)
    # ComCat rows are revised over time: keep the LAST occurrence of each id
    # (the freshly downloaded one wins over the cached one).
    merged = merged.drop_duplicates(subset="id", keep="last")
    n_dupes = n_pre - len(merged)
    merged = merged.sort_values("time", kind="mergesort").reset_index(drop=True)
    merged = merged[(merged.time >= pd.Timestamp(CATALOG_START, tz="UTC")) &
                    (merged.time <= now)]
    merged = merged.reset_index(drop=True)

    # ---- assertions demanded by the protocol
    assert merged["id"].is_unique, "duplicate event ids survived dedup"
    assert merged["time"].is_monotonic_increasing, "times not monotonic after sort"
    assert len(merged) > 10000, f"implausibly small catalog: {len(merged)} rows"
    assert merged.mag.min() >= MINMAG - 1e-9
    assert merged.latitude.between(LAT_MIN, LAT_MAX).all()
    assert merged.longitude.between(LON_MIN, LON_MAX).all()

    merged.to_csv(CACHE, index=False)
    log_lines.append(f"# merged rows before dedup={n_pre}, duplicate ids removed="
                     f"{n_dupes}, final rows={len(merged)}")
    log_lines.append(f"# range {merged.time.min().isoformat()} -> "
                     f"{merged.time.max().isoformat()}")
    log_lines.append(f"# cache written: {CACHE}")

    logf = DATA / f"comcat_log_{now.strftime('%Y%m%d')}.txt"
    logf.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print(f"[log] {logf}")
    return merged, {"windows": len(wins), "window_counts": counts,
                    "rows_before_dedup": int(n_pre), "duplicate_ids_removed": int(n_dupes),
                    "log_file": str(logf), "cache_file": str(CACHE)}


# ================================================================== ETAS pieces
def omori_integral(dt0, dt1, c, p):
    """int_{dt0}^{dt1} (s + c)^(-p) ds, elementwise, p != 1."""
    s = 1.0 - p
    return ((dt1 + c) ** s - (dt0 + c) ** s) / s


def lambda_at(t_now, t_hist, m_hist, par, trunc=TRUNC_DAYS):
    """Intensity at t_now (days) from history, with per-source contributions."""
    dt = t_now - t_hist
    use = (dt > 0) & (dt <= trunc)
    w = 10.0 ** (par["alpha"] * (m_hist[use] - par["M0"]))
    contrib = par["K"] * w * (dt[use] + par["c"]) ** (-par["p"])
    return par["mu"] + contrib.sum(), contrib, np.flatnonzero(use)


def expected_from_history(t_now, t_hist, m_hist, par, T=HORIZON_DAYS, trunc=TRUNC_DAYS):
    """Integral over [now, now+T] of mu + (kernel of EXISTING events only).

    This is the 'no future triggering' expected count of M >= M0 events: it counts
    background plus direct aftershocks of what is already in the catalog, and
    deliberately omits every event triggered by events that have not happened yet.
    Hence a lower bound on the true ETAS expectation.
    """
    dt = t_now - t_hist
    use = (dt > 0) & (dt <= trunc)
    w = 10.0 ** (par["alpha"] * (m_hist[use] - par["M0"]))
    trig = par["K"] * (w * omori_integral(dt[use], dt[use] + T, par["c"], par["p"])).sum()
    return par["mu"] * T + trig, par["mu"] * T, trig


def gr_frac(m, par):
    """Fraction of GR-distributed magnitudes (M0..MAG_CAP, b) that are >= m."""
    if m <= par["M0"]:
        return 1.0
    if m >= MAG_CAP:
        return 0.0
    b = par["b"]
    num = 10.0 ** (-b * (m - par["M0"])) - 10.0 ** (-b * (MAG_CAP - par["M0"]))
    den = 1.0 - 10.0 ** (-b * (MAG_CAP - par["M0"]))
    return num / den


def sample_mags(rng, n, par):
    """Truncated Gutenberg-Richter on [M0, MAG_CAP], inverse CDF."""
    u = rng.random(n)
    b = par["b"]
    top = 1.0 - 10.0 ** (-b * (MAG_CAP - par["M0"]))
    return par["M0"] - np.log10(1.0 - u * top) / b


def sample_omori_times(rng, n, t_src, T_end, par):
    """Sample n offspring times in (t_src, T_end] from the normalised Omori density.

    Inverse-CDF on F(t) proportional to int_0^{t-t_src} (s+c)^-p ds.
    t_src, T_end are simulation-clock days; t_src may be negative (history event).
    """
    c, p = par["c"], par["p"]
    s = 1.0 - p
    lo = np.maximum(0.0, -t_src)          # elapsed time already used up
    hi = T_end - t_src
    A = (lo + c) ** s
    B = (hi + c) ** s
    u = rng.random(n)
    val = A + u * (B - A)
    return t_src + val ** (1.0 / s) - c


def simulate(t_now, t_hist, m_hist, par, n_sims=N_SIMS, T=HORIZON_DAYS,
             seed=SEED, trunc=TRUNC_DAYS):
    """Forward ETAS by CLUSTER (branching) generation, not thinning.

    Each source (a history event, or a simulated one) spawns an independent Poisson
    number of direct offspring inside the horizon with mean
        K * 10^(alpha (M_i - M0)) * int over the remaining window of (t - t_i + c)^-p,
    offspring times drawn from the normalised Omori density restricted to that window,
    magnitudes i.i.d. truncated GR(b, M0, MAG_CAP). Recursed generation by generation
    until a generation is empty. This is exact for the ETAS intensity above (an ETAS
    process is a Poisson cluster process), and avoids thinning's rejection overhead.
    """
    rng = np.random.default_rng(seed)
    dt = t_now - t_hist
    use = (dt > 0) & (dt <= trunc)
    hist_t = -dt[use]                      # simulation clock: now = 0
    hist_m = m_hist[use]
    hist_w = 10.0 ** (par["alpha"] * (hist_m - par["M0"]))
    # expected direct offspring of each history event inside [0, T]
    hist_mean = par["K"] * hist_w * omori_integral(-hist_t, -hist_t + T, par["c"], par["p"])

    mu_mean = par["mu"] * T
    counts, max_mags, capped = np.zeros(n_sims, dtype=np.int64), np.full(n_sims, -np.inf), 0

    for i in range(n_sims):
        n_tot = 0
        mmax = -np.inf
        # ---- generation 0: background + direct offspring of the real history
        n_bg = rng.poisson(mu_mean)
        t_new = list(rng.random(n_bg) * T)
        n_off = rng.poisson(hist_mean)
        nz = np.flatnonzero(n_off)
        for j in nz:
            t_new.extend(sample_omori_times(rng, int(n_off[j]), hist_t[j], T, par))
        t_new = np.asarray(t_new, dtype=float)
        while t_new.size:
            n_tot += t_new.size
            if n_tot > MAX_EVENTS_PER_SIM:
                capped += 1
                break
            m_new = sample_mags(rng, t_new.size, par)
            mmax = max(mmax, float(m_new.max()))
            w = 10.0 ** (par["alpha"] * (m_new - par["M0"]))
            mean = par["K"] * w * omori_integral(0.0, T - t_new, par["c"], par["p"])
            k = rng.poisson(mean)
            nz = np.flatnonzero(k)
            if nz.size == 0:
                break
            child = []
            for j in nz:
                child.extend(sample_omori_times(rng, int(k[j]), t_new[j], T, par))
            t_new = np.asarray(child, dtype=float)
        counts[i] = n_tot
        max_mags[i] = mmax
    return counts, max_mags, capped


# ========================================================================= main
def main():
    ap = argparse.ArgumentParser(description="EXP-L live 7-day ETAS forecast")
    ap.add_argument("--now", default=None, help="ISO UTC override, default = current UTC")
    ap.add_argument("--refresh", action="store_true", help="ignore cache, full re-download")
    ap.add_argument("--nsims", type=int, default=N_SIMS)
    args = ap.parse_args()

    now = (pd.Timestamp(args.now).tz_localize("UTC")
           if args.now and pd.Timestamp(args.now).tz is None
           else pd.Timestamp(args.now) if args.now
           else pd.Timestamp.now(tz="UTC"))
    now = now.tz_convert("UTC")
    run_utc = pd.Timestamp.now(tz="UTC")

    par, exp_h = load_frozen()
    print("=" * 78)
    print("EXP-L  live 7-day forecast from FROZEN ETAS (no refitting)")
    print(f"  now (forecast origin) = {now.isoformat()}")
    print(f"  frozen params: mu={par['mu']:.6f}/d K={par['K']:.6f} alpha={par['alpha']:.4f} "
          f"c={par['c']:.6f} d p={par['p']:.4f} M0={par['M0']} b={par['b']:.4f}")
    print("=" * 78)

    cat, dl = download_catalog(now, refresh=args.refresh)
    n_all = len(cat)
    non_eq = cat[cat["type"] != "earthquake"]
    cat = cat[cat["type"] == "earthquake"].reset_index(drop=True)
    print(f"[catalog] {n_all} rows, {len(non_eq)} non-earthquake "
          f"({sorted(non_eq['type'].unique().tolist())}) removed -> {len(cat)} events")
    print(f"[catalog] {cat.time.min().isoformat()} -> {cat.time.max().isoformat()}")

    epoch = pd.Timestamp("2010-01-01", tz="UTC")
    t_hist = (cat.time - epoch).dt.total_seconds().to_numpy() / 86400.0
    m_hist = cat.mag.to_numpy(dtype=float)
    t_now = (now - epoch).total_seconds() / 86400.0

    # ---- current state
    lam_now, contrib, idx = lambda_at(t_now, t_hist, m_hist, par)
    ratio = lam_now / par["mu"]
    order = np.argsort(contrib)[::-1][:10]
    top = []
    for k in order:
        r = cat.iloc[int(idx[k])]
        top.append({
            "id": r["id"], "time": r["time"].isoformat(), "mag": float(r["mag"]),
            "place": str(r.get("place", "")),
            "days_before_now": float(t_now - t_hist[idx[k]]),
            "contribution_per_day": float(contrib[k]),
            "pct_of_lambda": float(100.0 * contrib[k] / lam_now),
        })
    # empirical mean rate of the driving catalog: the honest yardstick, since mu is
    # only the BACKGROUND part of the model and is far below the observed mean rate.
    span_days = (cat.time.max() - cat.time.min()).total_seconds() / 86400.0
    emp_rate = len(cat) / span_days
    ratio_emp = lam_now / emp_rate
    print(f"\n[state] lambda(now) = {lam_now:.4f} /day   mu = {par['mu']:.4f} /day   "
          f"ratio = {ratio:.3f}x")
    print(f"[state] for scale: the observed 2010-now ComCat mean rate is "
          f"{emp_rate:.3f} /day, so lambda(now) is {ratio_emp:.2f}x the OBSERVED mean. "
          f"mu is only {100 * par['mu'] / emp_rate:.0f}% of the mean rate (the rest is "
          f"triggering), so 'x mu' overstates how elevated the region is.")
    print(f"[state] history truncation = {TRUNC_DAYS:.0f} d "
          f"({int((t_hist > t_now - TRUNC_DAYS).sum())} of {len(cat)} events contribute)")
    print("[state] top contributors to lambda(now):")
    for e in top:
        print(f"   M{e['mag']:.2f}  {e['time'][:19]}  {e['days_before_now']:8.2f} d ago  "
              f"{e['contribution_per_day']:.5f}/d ({e['pct_of_lambda']:.1f}%)  {e['place'][:40]}")

    # ---- EXP-I(ii) alarm state
    alarm_mask = (t_hist > t_now - ALARM_WINDOW_DAYS) & (m_hist >= ALARM_MAG)
    alarm_on = bool(alarm_mask.any())
    alarm_events = [{"id": cat.iloc[int(i)]["id"], "time": cat.iloc[int(i)]["time"].isoformat(),
                     "mag": float(m_hist[i]), "place": str(cat.iloc[int(i)].get("place", ""))}
                    for i in np.flatnonzero(alarm_mask)]
    print(f"\n[alarm] EXP-I(ii) M>={ALARM_MAG} in last {ALARM_WINDOW_DAYS:.0f} d: "
          f"{'ON' if alarm_on else 'OFF (quiet state)'}")
    for e in alarm_events:
        print(f"   M{e['mag']:.2f} {e['time'][:19]} {e['place'][:50]}")

    # ---- analytic no-future-triggering lower bound
    Lam0, Lam_bg, Lam_trig = expected_from_history(t_now, t_hist, m_hist, par)
    bound = {}
    for mth in (4.0, 5.0, 6.0):
        lam_m = Lam0 * gr_frac(mth, par)
        bound[f"P_any_M{int(mth)}"] = float(1.0 - np.exp(-lam_m))
        bound[f"expected_N_M{int(mth)}"] = float(lam_m)
    bound["expected_N_M2.5"] = float(Lam0)
    bound["background_part"] = float(Lam_bg)
    bound["existing_aftershock_part"] = float(Lam_trig)
    print(f"\n[analytic lower bound, no future triggering] expected M>=2.5 count over 7 d "
          f"= {Lam0:.2f} (background {Lam_bg:.2f} + existing-aftershocks {Lam_trig:.2f})")
    for mth in (4, 5, 6):
        print(f"   P(any M>={mth}) >= {bound[f'P_any_M{mth}']:.5f}")

    # ---- simulations
    print(f"\n[sim] {args.nsims} forward ETAS cluster simulations over [now, now+7 d], "
          f"seed={SEED} ...")
    t0 = time.time()
    counts, mmax, capped = simulate(t_now, t_hist, m_hist, par, n_sims=args.nsims, seed=SEED)
    print(f"[sim] done in {time.time() - t0:.1f} s; sims hitting the "
          f"{MAX_EVENTS_PER_SIM} event cap: {capped}")

    sim = {
        "n_sims": int(args.nsims), "seed": SEED, "capped_sims": int(capped),
        "expected_N_M2.5_mean": float(counts.mean()),
        "N_M2.5_p05": float(np.percentile(counts, 5)),
        "N_M2.5_p50": float(np.percentile(counts, 50)),
        "N_M2.5_p95": float(np.percentile(counts, 95)),
        "N_M2.5_max": int(counts.max()),
    }
    for mth in (4.0, 5.0, 6.0):
        p = float((mmax >= mth).mean())
        sim[f"P_any_M{int(mth)}"] = p
        # binomial 95% Wald-ish interval on the MC estimate itself
        se = float(np.sqrt(max(p * (1 - p), 1e-12) / args.nsims))
        sim[f"P_any_M{int(mth)}_mc_se"] = se
    print(f"[sim] expected M>=2.5 in 7 d: mean {sim['expected_N_M2.5_mean']:.1f}  "
          f"(5th-95th pct {sim['N_M2.5_p05']:.0f}-{sim['N_M2.5_p95']:.0f}, "
          f"median {sim['N_M2.5_p50']:.0f})")
    # ---- consistency check: for a subcritical horizon the branching process inflates
    # the seed expectation by 1/(1-n_eff); recover n_eff and compare to theory.
    # E[10^(alpha(M-M0))] under truncated GR on [M0, MAG_CAP]
    b, al = par["b"], par["alpha"]
    top_t = 1.0 - 10.0 ** (-b * (MAG_CAP - par["M0"]))
    e_w = (b / (b - al)) * (1.0 - 10.0 ** (-(b - al) * (MAG_CAP - par["M0"]))) / top_t
    n_eff_theory = float(par["K"] * e_w
                         * omori_integral(0.0, HORIZON_DAYS, par["c"], par["p"]))
    inflation_theory = 1.0 / (1.0 - n_eff_theory) if n_eff_theory < 1 else float("inf")
    inflation_sim = float(counts.mean() / Lam0)
    sim["effective_7d_branching_ratio"] = n_eff_theory
    sim["cascade_inflation_theory_upper"] = float(inflation_theory)
    sim["cascade_inflation_observed"] = inflation_sim
    print(f"[sim] check: 7-day effective branching ratio n_eff = {n_eff_theory:.3f} "
          f"(<1, so simulations terminate); mean/analytic-seed inflation = "
          f"{inflation_sim:.2f}x vs geometric-cascade ceiling {inflation_theory:.2f}x")

    mc_noise_flags = []
    for mth in (4, 5, 6):
        p = sim[f"P_any_M{mth}"]
        lo = bound[f"P_any_M{mth}"]
        note = ""
        if p < lo:
            note = "  <-- BELOW its own analytic floor: pure MC noise, see flags"
            mc_noise_flags.append(
                f"P(any M>={mth}): the {args.nsims}-simulation estimate ({p:.4f}) came out "
                f"BELOW the analytic no-future-triggering lower bound ({lo:.4f}), which is "
                f"impossible in expectation. This is Monte Carlo noise: with {args.nsims} "
                f"sims the resolution is {1.0 / args.nsims:.4f} and only "
                f"{int(round(p * args.nsims))} simulation(s) produced such an event. Use the "
                f"analytic bound as the floor and treat the simulated value as uninformative "
                f"for this magnitude class.")
        print(f"[sim] P(any M>={mth}) = {p:.4f} "
              f"(+/- {1.96 * sim[f'P_any_M{mth}_mc_se']:.4f} MC)"
              f"  [analytic floor {lo:.4f}]{note}")

    # ---- output json
    out = {
        "experiment": "EXP-L",
        "kind": "application (forecast), NOT a hypothesis test",
        "protocol": "PATTERN_PROTOCOL.md :: EXP-L",
        "state_class": "first-run" if not sorted(HERE.glob("results_forecast_*.json"))
                       else "repeat-run (daily reusable script)",
        "run_utc": run_utc.isoformat(),
        "forecast_origin_utc": now.isoformat(),
        "horizon_days": HORIZON_DAYS,
        "refit": False,
        "frozen_params": par,
        "frozen_params_source": "results_exp_h.json :: train_fit.frozen_params "
                                "(+ b_value_train_aki); fit on SCSN 1982-2010, "
                                "scored once on 2010-2018",
        "intensity_convention": "lambda(t) = mu + sum_{t_i<t} K*10^(alpha*(M_i-M0))"
                                "*(t-t_i+c)^(-p), t in days",
        "catalog": {
            "source": "USGS ComCat FDSN event/1/query, format=csv",
            "box": {"lat": [LAT_MIN, LAT_MAX], "lon": [LON_MIN, LON_MAX]},
            "minmagnitude": MINMAG,
            "start": CATALOG_START,
            "end": now.isoformat(),
            "n_rows_downloaded": int(n_all),
            "n_non_earthquake_removed": int(len(non_eq)),
            "non_earthquake_types": sorted(non_eq["type"].unique().tolist()),
            "n_events_used": int(len(cat)),
            "t_first": cat.time.min().isoformat(),
            "t_last": cat.time.max().isoformat(),
            "max_mag": float(cat.mag.max()),
            **dl,
        },
        "current_state": {
            "lambda_now_per_day": float(lam_now),
            "mu_per_day": par["mu"],
            "lambda_over_mu": float(ratio),
            "observed_mean_rate_per_day_2010_now": float(emp_rate),
            "lambda_over_observed_mean_rate": float(ratio_emp),
            "mu_as_fraction_of_observed_mean_rate": float(par["mu"] / emp_rate),
            "ratio_caveat": "lambda/mu is NOT 'x above normal': mu is the background term "
                            "only and is far below the observed mean rate, most of which is "
                            "triggering. Compare lambda(now) to the observed mean rate.",
            "history_truncation_days": TRUNC_DAYS,
            "n_events_within_truncation": int((t_hist > t_now - TRUNC_DAYS).sum()),
            "truncation_note": "contributions from events older than 1000 d are negligible "
                               "at these params (same truncation EXP-H used for the fit; "
                               "EXP-H measured delta_LL/event = -0.0093 vs untruncated)",
            "top_contributors": top,
        },
        "alarm_state_exp_i_ii": {
            "definition": f"any M>={ALARM_MAG} in the trailing {ALARM_WINDOW_DAYS:.0f} d",
            "alarm_on": alarm_on,
            "events": alarm_events,
            "exp_i_ii_test_result_for_reference": {
                "n_eligible": 30, "n_hits": 18, "p_observed": 0.6,
                "poisson_base_rate": 0.0619,
                "note": "from results_exp_i.json; the conditional rate AFTER an M>=5, "
                        "not a claim about the quiet state",
            },
        },
        "analytic_lower_bound_no_future_triggering": bound,
        "simulation": sim,
        "simulation_method": "exact ETAS cluster/branching generation (not thinning): "
                             "background Poisson(mu*T) + Poisson offspring of every history "
                             "event and of every simulated event, Omori-density times, "
                             "i.i.d. truncated-GR magnitudes (b, M0=2.5, cap M8), recursed "
                             "until a generation is empty",
        "flags": [
            "APPLICATION, NOT A TEST: no hypothesis is being evaluated here; nothing about "
            "these numbers confirms or refutes EXP-H.",
            "Catalog source differs from the fitting catalog: params were fit on SCSN "
            "(1982-2010); this forecast is driven by ComCat (ANSS composite, includes NCSN/"
            "NN authoritative regions inside the box). Completeness and magnitude scales are "
            "not identical to SCSN, which is an unquantified error source.",
            "Recent-event magnitudes and locations are preliminary/automatic and get revised.",
            "Temporal ETAS only - no spatial component, so this is a box-wide rate, not a map.",
            "Stationary GR with b frozen at the 1982-2010 train value; no time-varying b.",
            "EXP-H's fitted branching ratio n = 1.161 is formally supercritical over infinite "
            "time; over a 7-day horizon the effective branching ratio is well below 1, so the "
            "simulations terminate, but the parameter set would not be usable for long horizons.",
            "MC uncertainty on each probability is +/- ~1.96*sqrt(p(1-p)/n_sims); for rare "
            "classes (M>=6) the simulated estimate is dominated by MC noise and the analytic "
            "bound is the more informative number.",
            f"CALIBRATION: the simulated 7-day mean of {counts.mean():.1f} M2.5+ events "
            f"({counts.mean() / HORIZON_DAYS:.2f}/day) sits against an observed ComCat mean "
            f"of {emp_rate:.2f}/day for 2010-now. From a quiet state the frozen model runs "
            f"{'below' if counts.mean() / HORIZON_DAYS < emp_rate else 'above'} the "
            f"long-run observed rate; the model's own unconditional mean is undefined "
            f"(supercritical n), so its stationary rate cannot be checked directly. Treat "
            f"the absolute counts as less trustworthy than the relative/conditional signal.",
        ] + mc_noise_flags,
    }
    outf = HERE / f"results_forecast_{now.strftime('%Y%m%d')}.json"
    with open(outf, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\n[out] {outf}")

    # ---- plain-language statement
    st = "ELEVATED (an M>=5 occurred in the last 7 days)" if alarm_on else \
         "QUIET (no M>=5 in the last 7 days)"
    print("\n" + "=" * 78)
    print("PLAIN-LANGUAGE STATEMENT")
    print("=" * 78)
    print(f"""
This is a PROBABILITY FORECAST, not an earthquake prediction. Nobody here is
saying an earthquake will or will not happen. What follows is what a specific,
pre-registered statistical model - an ETAS aftershock model whose five parameters
were fixed on 1982-2010 data and scored once, without adjustment, on 2010-2018 -
says about the next 7 days in southern California, given every M2.5+ event USGS
has recorded up to {now.strftime('%Y-%m-%d %H:%M')} UTC.

Current state: the model's event rate right now is {lam_now:.2f} M2.5+ per day.
That is {ratio:.1f}x the model's BACKGROUND term (mu = {par['mu']:.2f}/day), but that
comparison flatters it - background is only {100 * par['mu'] / emp_rate:.0f}% of what this region
actually does. Against the observed 2010-now average of {emp_rate:.2f}/day, today
is {ratio_emp:.2f}x, i.e. {'somewhat above' if ratio_emp > 1.15 else 'somewhat below' if ratio_emp < 0.87 else 'about'} normal. The region is in the
{st.split(' (')[0]} state by the EXP-I(ii) definition: {st.split('(')[1][:-1]}.

Over the next 7 days the model expects about {sim['expected_N_M2.5_mean']:.0f} M2.5+
events (90% of simulations fall between {sim['N_M2.5_p05']:.0f} and {sim['N_M2.5_p95']:.0f}),
and gives:

    P(at least one M>=4) = {sim['P_any_M4'] * 100:.1f}%   (analytic floor {bound['P_any_M4'] * 100:.1f}%)
    P(at least one M>=5) = {sim['P_any_M5'] * 100:.2f}%   (analytic floor {bound['P_any_M5'] * 100:.2f}%)
    P(at least one M>=6) = {sim['P_any_M6'] * 100:.2f}%   (analytic floor {bound['P_any_M6'] * 100:.2f}%)

The "analytic floor" is what you get counting only background seismicity plus
aftershocks of earthquakes that have ALREADY happened - it ignores everything
that future earthquakes would themselves trigger, so the model's true value is
always at or above it. The gap between the two IS the cascade: the chance that
something happens in the next week and then sets off its own aftershocks. Where
the simulated number falls below its floor (it can, for the rarest class), that
is Monte Carlo noise from only {args.nsims} simulations, not a real result - take
the floor.

What this forecast is good for: the model demonstrated real, measured skill on
held-out data (EXP-H) - but that skill is CONDITIONAL and it is concentrated in
the hours-to-days after a large event, when aftershock physics dominates. In a
quiet state like a typical week, these numbers are close to the long-run base
rate, and the model is telling you very little you could not get from a
climatological table. Do not read a low number as safety, and do not read an
elevated number as a warning of a specific event: a 1-in-1000 weekly chance of
an M6 is both small and non-zero, every week, forever.

Standing caveats: ETAS is temporal only (no location), magnitudes are drawn from
a fixed Gutenberg-Richter law with b={par['b']:.3f}, the driving catalog is ComCat
while the parameters were fit on SCSN, and the most recent events carry
preliminary magnitudes that will be revised.
""".rstrip())
    print("=" * 78)


if __name__ == "__main__":
    main()
