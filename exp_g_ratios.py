"""EXP-G: ratio-sequence miner (the Fibonacci question, disciplined).
Implements PATTERN_PROTOCOL.md EXP-G section (frozen 2026-08-09).

Pipeline:
  1. Load SCSN_original_catalog.txt (raw column order, auto-detected), SoCal box.
  2. Sequences: mainshocks M>=5.0; aftershocks = events within 100 d & 50 km (haversine) of
     the mainshock, mag < mainshock mag, claimed by the LARGEST qualifying mainshock (ties:
     earliest) rather than the nearest -- i.e. never "already claimed by a larger earlier
     mainshock whose window covers them". Sequences split train/test by MAINSHOCK time
     (train < 2010-01-01 UTC); a train sequence's aftershock window is capped at
     2009-12-31 23:59:59 UTC to avoid test leakage. Keep sequences with >= 20 aftershocks
     (post-cap for train).
  3. Ratios: sorted aftershock times -> interevent Δt (drop Δt=0) -> r_i = Δt_{i+1}/Δt_i;
     pool log10(r) across train sequences.
  4. Null: per-sequence Omori (K,c,p) MLE fit (bounds c in [1e-4,10] d, p in [0.5,2.5]);
     500 synthetic sequences per train sequence via inverse-CDF sampling of n events from the
     fitted, window-normalized Omori density; pool synthetic log10(r), keeping per-sim pools.
  5. Candidate bands: log10(r) histogram, width-0.1 bands over [-3,3]; candidate if observed
     density > 99th percentile of the 500 per-band simulation densities. Golden-ratio band
     [0.2,0.3) (log10(1.618)=0.209) evaluated explicitly regardless of candidacy.
  6. TEST: candidate bands (and the golden band, always) rescored on test-period sequences
     (>=20 aftershocks, own 500-sim fitted-Omori nulls); confirmed if p_test < 0.05.
Outputs: results_exp_g.json + maps/exp_g_ratio_hist.png. Raw inputs never modified.
"""
import json
import time
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import minimize

HERE = Path(__file__).parent
DATA = HERE / "data" / "xue_lu_zenodo"
MAPS = HERE / "maps"
MAPS.mkdir(exist_ok=True)

EARTH_R = 6371.0
BOX = dict(lat=(31.5, 38.0), lon=(-122.0, -113.5))
MS_MIN_MAG = 5.0
AFTERSHOCK_WINDOW_D = 100.0
AFTERSHOCK_DIST_KM = 50.0
MIN_AFTERSHOCKS = 20
TRAIN_CUTOFF = pd.Timestamp("2010-01-01", tz="UTC")
TRAIN_CAP = pd.Timestamp("2009-12-31 23:59:59", tz="UTC")
N_SIMS = 500
BAND_WIDTH = 0.1
BAND_LO, BAND_HI = -3.0, 3.0
GOLDEN_LOG10 = float(np.log10(1.618033988749895))  # 0.20898...
SEED = 20260809

RNG = np.random.default_rng(SEED)


# ---------------------------------------------------------------- catalog
def load_catalog(fname):
    """Auto-detect raw vs declustered column order (xue_lu_crosstest.py::load_catalog
    pattern). SCSN_original_catalog.txt is raw: yr mo dy hr mi sec EID lat lon depth mag."""
    raw_cols = ["yr", "mo", "dy", "hr", "mi", "sec", "eid", "lat", "lon", "depth", "mag"]
    dec_cols = ["yr", "mo", "dy", "hr", "mi", "sec", "lat", "lon", "depth", "mag", "eid"]
    probe = pd.read_csv(DATA / fname, sep=r"\s+", header=None, nrows=1000)
    cols = raw_cols if probe[6].abs().max() > 90 else dec_cols
    df = pd.read_csv(DATA / fname, sep=r"\s+", header=None, names=cols)
    assert df.lat.abs().max() <= 90 and df.lon.abs().max() <= 180, f"column detection failed for {fname}"
    sec = df["sec"].astype(float)
    ts = pd.to_datetime(
        dict(year=df.yr, month=df.mo, day=df.dy, hour=df.hr, minute=df.mi, second=0), utc=True
    ) + pd.to_timedelta(sec, unit="s")
    df["time"] = ts
    df["t_unix"] = (ts - pd.Timestamp(0, tz="UTC")).dt.total_seconds().to_numpy()
    return df[["time", "t_unix", "lat", "lon", "depth", "mag", "eid"]].reset_index(drop=True)


def haversine_km(lat1, lon1, lat2, lon2):
    la1, lo1, la2, lo2 = map(np.radians, [lat1, lon1, lat2, lon2])
    a = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * EARTH_R * np.arcsin(np.sqrt(a))


# ---------------------------------------------------------------- sequence construction
def build_sequences(df):
    """Claim aftershocks by the LARGEST qualifying mainshock (ties: earliest), using the
    full uncapped 100 d / 50 km / mag<ms.mag window on the full (train+test) catalog -- this
    determines true physical ownership. Train-vs-test split and the train leakage cap are
    applied afterward, per-sequence, when building each sequence's used aftershock list."""
    ms = df[df.mag >= MS_MIN_MAG].sort_values(["mag", "time"], ascending=[False, True]).copy()
    ms_order = ms.index.to_numpy()  # claim order: largest mag first, ties earliest first
    t_all = df.t_unix.to_numpy()
    lat_all = df.lat.to_numpy()
    lon_all = df.lon.to_numpy()
    mag_all = df.mag.to_numpy()

    claimed_by = np.full(len(df), -1, dtype=np.int64)
    max_dt_s = AFTERSHOCK_WINDOW_D * 86400.0
    for ms_idx in ms_order:
        mt, mlat, mlon, mmag = t_all[ms_idx], lat_all[ms_idx], lon_all[ms_idx], mag_all[ms_idx]
        dt = t_all - mt
        cand = (claimed_by < 0) & (dt > 0) & (dt <= max_dt_s) & (mag_all < mmag)
        if not cand.any():
            continue
        idx = np.where(cand)[0]
        d = haversine_km(mlat, mlon, lat_all[idx], lon_all[idx])
        ok = idx[d <= AFTERSHOCK_DIST_KM]
        claimed_by[ok] = ms_idx

    seqs = []
    for ms_idx in ms.index:
        aft_idx = np.where(claimed_by == ms_idx)[0]
        if len(aft_idx) == 0:
            continue
        mt = t_all[ms_idx]
        ms_time = df.time.iloc[ms_idx]
        is_train = ms_time < TRAIN_CUTOFF
        full_window_end = mt + max_dt_s
        if is_train:
            cap_ts = (TRAIN_CAP - pd.Timestamp(0, tz="UTC")).total_seconds()
            window_end = min(full_window_end, cap_ts)
        else:
            window_end = full_window_end
        window_T_d = (window_end - mt) / 86400.0
        if window_T_d <= 0:
            continue
        keep = aft_idx[t_all[aft_idx] <= window_end]
        if len(keep) < MIN_AFTERSHOCKS:
            continue
        t_days = np.sort((t_all[keep] - mt) / 86400.0)
        seqs.append(dict(
            ms_idx=int(ms_idx), ms_time=ms_time, ms_mag=float(mag_all[ms_idx]),
            ms_lat=float(lat_all[ms_idx]), ms_lon=float(lon_all[ms_idx]),
            is_train=bool(is_train), window_T=float(window_T_d),
            n_aftershocks=int(len(keep)), t_days=t_days,
        ))
    return seqs


# ---------------------------------------------------------------- Omori MLE + simulation
def _omori_I(c, p, T):
    if abs(p - 1.0) < 1e-6:
        return np.log((T + c) / c)
    return ((T + c) ** (1 - p) - c ** (1 - p)) / (1 - p)


def _omori_negloglik(params, t, T):
    c, p = params
    I = _omori_I(c, p, T)
    if not np.isfinite(I) or I <= 0:
        return 1e12
    n = len(t)
    K = n / I
    if K <= 0:
        return 1e12
    ll = n * np.log(K) - p * np.sum(np.log(t + c)) - K * I
    if not np.isfinite(ll):
        return 1e12
    return -ll


def fit_omori(t_days, T):
    best = None
    for c0, p0 in [(0.1, 1.0), (0.01, 1.1), (1.0, 0.9), (0.5, 1.5)]:
        res = minimize(_omori_negloglik, x0=[c0, p0], args=(t_days, T),
                        method="L-BFGS-B", bounds=[(1e-4, 10.0), (0.5, 2.5)])
        if best is None or res.fun < best.fun:
            best = res
    c, p = best.x
    I = _omori_I(c, p, T)
    K = len(t_days) / I
    return float(K), float(c), float(p)


def omori_inverse_cdf(u, c, p, T):
    I = _omori_I(c, p, T)
    if abs(p - 1.0) < 1e-6:
        t = c * (np.exp(u * I) - 1.0)
    else:
        val = c ** (1 - p) + u * I * (1 - p)
        val = np.clip(val, 1e-300, None)
        t = val ** (1.0 / (1.0 - p)) - c
    return np.clip(t, 0.0, T)


def simulate_sequences(K, c, p, T, n, n_sims, rng):
    """n_sims synthetic sequences of n events each, drawn from the fitted normalized Omori
    density over [0,T] via inverse-CDF (K is not needed for sampling shape; kept for record)."""
    u = rng.uniform(0.0, 1.0, size=(n_sims, n))
    t_sim = omori_inverse_cdf(u, c, p, T)
    t_sim.sort(axis=1)
    return t_sim


# ---------------------------------------------------------------- ratios / bands
def sequence_ratios(t_days):
    dt = np.diff(t_days)
    dt = dt[dt > 0]
    if len(dt) < 2:
        return np.array([])
    r = dt[1:] / dt[:-1]
    r = r[(r > 0) & np.isfinite(r)]
    return np.log10(r)


def pooled_log10r(seqs_t_days_list):
    parts = [sequence_ratios(t) for t in seqs_t_days_list]
    parts = [p for p in parts if len(p)]
    if not parts:
        return np.array([])
    return np.concatenate(parts)


def band_edges():
    n_bands = int(round((BAND_HI - BAND_LO) / BAND_WIDTH))
    edges = BAND_LO + np.arange(n_bands + 1) * BAND_WIDTH
    return edges


def band_densities(log10r):
    edges = band_edges()
    n_bands = len(edges) - 1
    if len(log10r) == 0:
        return np.zeros(n_bands)
    counts, _ = np.histogram(log10r, bins=edges)
    return counts / len(log10r)


def band_label(i, edges):
    return f"[{edges[i]:.1f},{edges[i+1]:.1f})"


def band_index_for(value, edges):
    idx = np.searchsorted(edges, value, side="right") - 1
    if idx < 0 or idx >= len(edges) - 1:
        return None
    return int(idx)


# ---------------------------------------------------------------- per-sequence null pooling
def build_null_densities(seqs, n_sims, rng, label):
    """Fit Omori per sequence, simulate n_sims synthetic sequences each, and return
    (n_sims, n_bands) array of per-simulation POOLED densities across all sequences
    (each sim index k pools synthetic sequence #k from every real sequence)."""
    edges = band_edges()
    n_bands = len(edges) - 1
    sim_log10r_by_k = [[] for _ in range(n_sims)]
    fits = []
    t0 = time.time()
    for si, seq in enumerate(seqs):
        K, c, p = fit_omori(seq["t_days"], seq["window_T"])
        fits.append(dict(ms_idx=seq["ms_idx"], K=K, c=c, p=p))
        t_sim = simulate_sequences(K, c, p, seq["window_T"], seq["n_aftershocks"], n_sims, rng)
        for k in range(n_sims):
            lr = sequence_ratios(t_sim[k])
            if len(lr):
                sim_log10r_by_k[k].append(lr)
    dens = np.zeros((n_sims, n_bands))
    for k in range(n_sims):
        pooled = np.concatenate(sim_log10r_by_k[k]) if sim_log10r_by_k[k] else np.array([])
        dens[k] = band_densities(pooled)
    elapsed = time.time() - t0
    print(f"  [{label}] Omori fit+sim for {len(seqs)} sequences, {n_sims} sims: {elapsed:.1f}s")
    return dens, fits, elapsed


# ---------------------------------------------------------------- main
def main():
    t_start = time.time()
    print("Loading SCSN_original_catalog.txt ...")
    df = load_catalog("SCSN_original_catalog.txt")
    df = df[df.lat.between(*BOX["lat"]) & df.lon.between(*BOX["lon"])].reset_index(drop=True)
    print(f"SoCal box: {len(df)} events, {df.time.min()} .. {df.time.max()}, "
          f"M>=5.0 count: {(df.mag >= MS_MIN_MAG).sum()}")

    seqs = build_sequences(df)
    train_seqs = [s for s in seqs if s["is_train"]]
    test_seqs = [s for s in seqs if not s["is_train"]]
    print(f"Sequences with >={MIN_AFTERSHOCKS} aftershocks: train={len(train_seqs)}, test={len(test_seqs)}")
    for s in train_seqs:
        print(f"  TRAIN ms_idx={s['ms_idx']} M{s['ms_mag']:.2f} {s['ms_time']} n_aft={s['n_aftershocks']} "
              f"T={s['window_T']:.1f}d")
    for s in test_seqs:
        print(f"  TEST  ms_idx={s['ms_idx']} M{s['ms_mag']:.2f} {s['ms_time']} n_aft={s['n_aftershocks']} "
              f"T={s['window_T']:.1f}d")

    n_sims = N_SIMS
    # runtime guard: project total time from a timed subset, halve sims if projected > 90 min
    guard_note = None
    if train_seqs:
        rng_probe = np.random.default_rng(SEED - 1)
        probe_n = min(50, n_sims)
        _, _, probe_elapsed = build_null_densities(train_seqs, probe_n, rng_probe, "runtime-probe(train)")
        per_sim_s = probe_elapsed / max(probe_n, 1)
        projected_s = per_sim_s * n_sims * (1 + (len(test_seqs) / max(len(train_seqs), 1)))
        if projected_s > 90 * 60:
            n_sims = n_sims // 2
            guard_note = (f"projected runtime {projected_s/60:.1f} min > 90 min guard; "
                           f"halved N_SIMS {N_SIMS} -> {n_sims}")
            print(f"RUNTIME GUARD: {guard_note}")

    edges = band_edges()
    n_bands = len(edges) - 1
    golden_idx = band_index_for(GOLDEN_LOG10, edges)
    assert edges[golden_idx] <= 0.2 + 1e-9 and abs(edges[golden_idx] - 0.2) < 1e-9, \
        f"golden band edge mismatch: {edges[golden_idx]}"

    results = {
        "protocol": "PATTERN_PROTOCOL.md EXP-G",
        "n_train_sequences": len(train_seqs),
        "n_test_sequences": len(test_seqs),
        "band_width": BAND_WIDTH, "band_range": [BAND_LO, BAND_HI],
        "n_sims": n_sims, "runtime_guard": guard_note,
    }

    if not train_seqs:
        results["verdict"] = "NO TRAIN SEQUENCES with >=20 aftershocks -- protocol cannot run."
        (HERE / "results_exp_g.json").write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
        print(results["verdict"])
        return

    # ---- TRAIN: observed pooled densities + null ----
    train_log10r = pooled_log10r([s["t_days"] for s in train_seqs])
    obs_dens_train = band_densities(train_log10r)
    rng_train = np.random.default_rng(SEED)
    null_dens_train, train_fits, _ = build_null_densities(train_seqs, n_sims, rng_train, "train")
    null_p99_train = np.percentile(null_dens_train, 99, axis=0)

    candidate_mask = obs_dens_train > null_p99_train
    candidate_bands = []
    for i in range(n_bands):
        if candidate_mask[i]:
            null_k = null_dens_train[:, i]
            p_train = (np.sum(null_k >= obs_dens_train[i]) + 1) / (n_sims + 1)
            candidate_bands.append(dict(
                band=band_label(i, edges), lo=float(edges[i]), hi=float(edges[i + 1]),
                observed_density_train=float(obs_dens_train[i]),
                null_p99_train=float(null_p99_train[i]),
                null_mean_train=float(np.mean(null_k)),
                p_train=float(p_train),
            ))
    print(f"Train candidate bands (obs density > null p99): {len(candidate_bands)}")
    for cb in candidate_bands:
        print(f"  {cb['band']}: obs={cb['observed_density_train']:.4f} "
              f"null_p99={cb['null_p99_train']:.4f} p={cb['p_train']:.4f}")

    # ---- golden band (train), regardless of candidacy ----
    gi = golden_idx
    golden_null_k = null_dens_train[:, gi]
    golden_p_train = (np.sum(golden_null_k >= obs_dens_train[gi]) + 1) / (n_sims + 1)
    golden = {
        "band": band_label(gi, edges), "lo": float(edges[gi]), "hi": float(edges[gi + 1]),
        "log10_phi": GOLDEN_LOG10,
        "train": {
            "observed_density": float(obs_dens_train[gi]),
            "null_p99": float(null_p99_train[gi]),
            "null_mean": float(np.mean(golden_null_k)),
            "p_train": float(golden_p_train),
            "is_candidate": bool(candidate_mask[gi]),
        },
    }
    print(f"Golden band {golden['band']} (log10 phi={GOLDEN_LOG10:.4f}): "
          f"obs={obs_dens_train[gi]:.4f} null_p99={null_p99_train[gi]:.4f} "
          f"p_train={golden_p_train:.4f} candidate={bool(candidate_mask[gi])}")

    # ---- TEST ----
    confirmed_bands = []
    if test_seqs:
        test_log10r = pooled_log10r([s["t_days"] for s in test_seqs])
        obs_dens_test = band_densities(test_log10r)
        rng_test = np.random.default_rng(SEED + 1)
        null_dens_test, test_fits, _ = build_null_densities(test_seqs, n_sims, rng_test, "test")

        for cb in candidate_bands:
            i = band_index_for((cb["lo"] + cb["hi"]) / 2, edges)
            null_k = null_dens_test[:, i]
            p_test = (np.sum(null_k >= obs_dens_test[i]) + 1) / (n_sims + 1)
            cb["test"] = dict(observed_density_test=float(obs_dens_test[i]),
                               null_mean_test=float(np.mean(null_k)),
                               p_test=float(p_test), confirmed=bool(p_test < 0.05))
            if cb["test"]["confirmed"]:
                confirmed_bands.append(cb["band"])
            print(f"  TEST {cb['band']}: obs={obs_dens_test[i]:.4f} p_test={p_test:.4f} "
                  f"confirmed={cb['test']['confirmed']}")

        golden_null_test_k = null_dens_test[:, gi]
        golden_p_test = (np.sum(golden_null_test_k >= obs_dens_test[gi]) + 1) / (n_sims + 1)
        golden["test"] = {
            "n_test_sequences": len(test_seqs),
            "observed_density": float(obs_dens_test[gi]),
            "null_mean": float(np.mean(golden_null_test_k)),
            "p_test": float(golden_p_test),
            "confirmed": bool(golden_p_test < 0.05),
        }
        print(f"Golden band TEST: obs={obs_dens_test[gi]:.4f} p_test={golden_p_test:.4f} "
              f"confirmed={golden['test']['confirmed']}")
    else:
        for cb in candidate_bands:
            cb["test"] = {"note": "no test sequences available (>=20 aftershocks, ms.time>=2010)"}
        golden["test"] = {"note": "no test sequences available (>=20 aftershocks, ms.time>=2010)"}
        print("No test sequences available -- TEST stage skipped.")

    if not train_seqs:
        verdict = "no train sequences"
    elif not candidate_bands:
        verdict = ("No candidate bands survived TRAIN screening (observed density never exceeded "
                   "the 99th-percentile Omori-simulation null in any 0.1-wide log10(r) band). "
                   "Golden ratio band shows nothing special "
                   f"(p_train={golden_p_train:.4f}). Consistent with the honest prior.")
    elif not confirmed_bands:
        verdict = (f"{len(candidate_bands)} candidate band(s) survived TRAIN screening but NONE "
                   f"confirmed on TEST (p_test<0.05). Golden ratio band "
                   f"{'was' if candidate_mask[gi] else 'was not'} a train candidate; "
                   f"test result: {golden.get('test', {}).get('p_test', 'n/a')}.")
    else:
        verdict = (f"{len(confirmed_bands)} band(s) confirmed on TEST: {confirmed_bands}. "
                   f"Flagged for replication on an independent region per protocol.")

    results.update({
        "mainshock_min_mag": MS_MIN_MAG,
        "aftershock_window_days": AFTERSHOCK_WINDOW_D,
        "aftershock_dist_km": AFTERSHOCK_DIST_KM,
        "min_aftershocks": MIN_AFTERSHOCKS,
        "train_sequences": [{k: v for k, v in s.items() if k != "t_days"} for s in train_seqs],
        "test_sequences": [{k: v for k, v in s.items() if k != "t_days"} for s in test_seqs],
        "train_omori_fits": train_fits,
        "candidate_bands": candidate_bands,
        "golden_band": golden,
        "confirmed_bands": confirmed_bands,
        "verdict": verdict,
        "runtime_seconds": round(time.time() - t_start, 1),
    })
    (HERE / "results_exp_g.json").write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
    print(f"\nVERDICT: {verdict}")
    print("-> results_exp_g.json")

    # ---- plot ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        centers = (edges[:-1] + edges[1:]) / 2
        lo_env = np.percentile(null_dens_train, 0.5, axis=0)
        hi_env = np.percentile(null_dens_train, 99.5, axis=0)
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.fill_between(centers, lo_env, hi_env, color="tab:blue", alpha=0.25,
                         label="Omori-sim envelope (0.5-99.5 pct, 500 sims)")
        ax.plot(centers, np.percentile(null_dens_train, 99, axis=0), color="tab:blue",
                lw=1, ls="--", label="sim 99th pct")
        ax.plot(centers, obs_dens_train, color="black", lw=1.5, label="observed (train)")
        ax.axvline(GOLDEN_LOG10, color="tab:red", ls=":", lw=1.5,
                   label=f"golden ratio band (log10 phi={GOLDEN_LOG10:.3f})")
        ax.set_xlabel("log10(r), r = Δt_{i+1}/Δt_i")
        ax.set_ylabel("density (fraction of pooled log10(r) in 0.1-wide band)")
        ax.set_title("EXP-G: interevent ratio distribution vs Omori-simulation null (train)")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(MAPS / "exp_g_ratio_hist.png", dpi=150)
        plt.close(fig)
        print("-> maps/exp_g_ratio_hist.png")
    except Exception as e:
        print(f"plot generation failed (non-fatal): {e}")


if __name__ == "__main__":
    main()
