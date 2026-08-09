"""EXP-F: periodicity comb - implements PATTERN_PROTOCOL.md EXP-F (frozen 2026-08-09).

Declustered SoCal catalog, Rayleigh-statistic periodicity search over a log-spaced comb of
periods (0.25 d - 3650 d) plus a named list of candidate periods (lunar/solar/tidal-adjacent).
Null distributions are 500 surrogate catalogs per period, drawn as an inhomogeneous Poisson
process whose rate is a Gaussian-kernel smoothing (sigma = 5*P, clipped) of the observed event
time density - this preserves slow rate structure (secular trends, catalog completeness
changes) while destroying any phase-locking at period P. TRAIN periods with p < 0.01 are
"detected"; each detection is re-scored ONCE on TEST (p < 0.05 + phase agreement within 60 deg
= "confirmed"). 0.5 d / 1 d / 7 d are pre-labeled artifacts (detector cadence / blasting) and
are never interpreted as physical detections even if they fire - they instead serve as the
positive-control method-check, run on the ORIGINAL (non-declustered) catalog, train period only.

Outputs: results_exp_f.json + maps/exp_f_comb.png + console summary. Raw inputs are never
modified.
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

HERE = Path(__file__).parent
DATA = HERE / "data" / "xue_lu_zenodo"
MAPS = HERE / "maps"

BOX = dict(lat=(31.5, 38.0), lon=(-122.0, -113.5))
SPLIT = pd.Timestamp("2010-01-01", tz="UTC")
MAG_FLOOR = 1.5

SEED = 20260809
N_SURR_DEFAULT = 500
RUNTIME_BUDGET_S = 90 * 60
MAX_ELEMS_PER_CHUNK = 8_000_000  # bound peak memory of surrogate arrays

NAMED_PERIODS = [0.5, 1, 7, 13.66, 14.77, 27.32, 27.55, 29.53, 182.62, 365.25]
ARTIFACT_PERIODS = {0.5, 1, 7}
N_LOG_SPACED = 60
PERIOD_MIN_D, PERIOD_MAX_D = 0.25, 3650.0


# ------------------------------------------------------------------ catalog loading
def load_catalog(fname, mag_floor=MAG_FLOOR, box=BOX):
    """Column order differs between the declustered and original SCSN files; auto-detect
    by validating the latitude column, same approach as xue_lu_crosstest.py::load_catalog."""
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
    df = df[
        (df.lat >= box["lat"][0]) & (df.lat <= box["lat"][1])
        & (df.lon >= box["lon"][0]) & (df.lon <= box["lon"][1])
        & (df.mag >= mag_floor)
    ].reset_index(drop=True)
    return df[["t_unix", "time", "lat", "lon", "depth", "mag", "eid"]]


# ------------------------------------------------------------------ period comb
def build_period_comb():
    log_grid = np.geomspace(PERIOD_MIN_D, PERIOD_MAX_D, N_LOG_SPACED)
    named = np.array(sorted(NAMED_PERIODS))
    keep_mask = np.ones(len(log_grid), dtype=bool)
    for p in named:
        close = np.abs(log_grid - p) / p <= 0.01
        keep_mask &= ~close
    periods = np.concatenate([log_grid[keep_mask], named])
    periods = np.unique(np.round(periods, 8))
    periods.sort()
    artifact = np.array([any(abs(p - a) / a <= 0.01 for a in ARTIFACT_PERIODS) for p in periods])
    return periods, artifact


# ------------------------------------------------------------------ Rayleigh statistic
def rayleigh_obs(t_unix, period_d):
    theta = (2 * np.pi / period_d) * np.mod(t_unix / 86400.0, period_d)
    z = np.exp(1j * theta).mean()
    return float(np.abs(z)), float(np.degrees(np.angle(z)))


def circular_diff_deg(a, b):
    d = (a - b + 180.0) % 360.0 - 180.0
    return abs(d)


def build_rate_probs(t_events, t0, t1, period_d):
    """Gaussian-kernel-smoothed event-time density on a regular grid; sigma = 5*P clipped
    to <= span/8. Daily grid for P >= 1 d, hourly grid for sub-day P (protocol spec)."""
    span_d = (t1 - t0) / 86400.0
    sigma_d = min(5.0 * period_d, span_d / 8.0)
    dt_grid = 86400.0 if period_d >= 1.0 else 3600.0
    n_bins = max(int(np.ceil((t1 - t0) / dt_grid)), 10)
    edges = t0 + np.arange(n_bins + 1) * dt_grid
    counts, _ = np.histogram(t_events, bins=edges)
    sigma_bins = max(sigma_d * 86400.0 / dt_grid, 1e-6)
    smoothed = gaussian_filter1d(counts.astype(float), sigma=sigma_bins, mode="nearest")
    smoothed = np.clip(smoothed, 1e-12, None)
    probs = smoothed / smoothed.sum()
    cump = np.cumsum(probs)
    cump[-1] = 1.0  # guard against fp roundoff leaving searchsorted out of range
    return edges, cump, dt_grid, n_bins


def surrogate_R(t_events, t0, t1, period_d, n_surr, rng):
    """500 (or n_surr) surrogate Rayleigh R's: inhomogeneous-Poisson-via-smoothed-rate
    surrogate catalogs, same expected n as observed, chunked to bound peak memory."""
    edges, cump, dt_grid, n_bins = build_rate_probs(t_events, t0, t1, period_d)
    n = len(t_events)
    Rs = np.empty(n_surr)
    chunk = max(1, MAX_ELEMS_PER_CHUNK // max(n, 1))
    done = 0
    while done < n_surr:
        m = min(chunk, n_surr - done)
        r = rng.random((m, n))
        idx = np.clip(np.searchsorted(cump, r), 0, n_bins - 1)
        jitter = rng.random((m, n)) * dt_grid
        t_sur = edges[idx] + jitter
        theta = (2 * np.pi / period_d) * np.mod(t_sur / 86400.0, period_d)
        Rs[done:done + m] = np.abs(np.exp(1j * theta).mean(axis=1))
        done += m
    return Rs


def rayleigh_test(t_events, t0, t1, period_d, n_surr, rng):
    R_obs, phase = rayleigh_obs(t_events, period_d)
    Rs = surrogate_R(t_events, t0, t1, period_d, n_surr, rng)
    p = (int(np.sum(Rs >= R_obs)) + 1) / (n_surr + 1)
    return R_obs, phase, p


# ------------------------------------------------------------------ main
def main():
    t_start = time.time()
    rng = np.random.default_rng(SEED)

    print("Loading declustered SoCal catalog (M>=%.1f)..." % MAG_FLOOR)
    cat = load_catalog("SCSN_decluster_m1.5.txt")
    train = cat[cat.time < SPLIT]
    test = cat[cat.time >= SPLIT]
    t_train = train.t_unix.to_numpy()
    t_test = test.t_unix.to_numpy()
    train_t0, train_t1 = t_train.min(), t_train.max()
    test_t0, test_t1 = t_test.min(), t_test.max()
    print(f"  train n={len(train)} span={train_t1-train_t0:.0f}s ({(train_t1-train_t0)/86400:.0f} d)")
    print(f"  test  n={len(test)} span={test_t1-test_t0:.0f}s ({(test_t1-test_t0)/86400:.0f} d)")

    periods, artifact_flags = build_period_comb()
    print(f"Period comb: {len(periods)} periods "
          f"({N_LOG_SPACED} log-spaced UNION {len(NAMED_PERIODS)} named, deduped)")

    n_surr = N_SURR_DEFAULT
    guard_applied = False
    guard_note = None

    comb_results = []
    per_period_times = []
    for i, (P, is_artifact) in enumerate(zip(periods, artifact_flags)):
        t0p = time.time()
        R_train, phase_train, p_train = rayleigh_test(t_train, train_t0, train_t1, P, n_surr, rng)
        detected = bool(p_train < 0.01)

        R_test = p_test = phase_test = None
        confirmed = False
        if detected:
            R_test, phase_test, p_test = rayleigh_test(t_test, test_t0, test_t1, P, n_surr, rng)
            confirmed = bool(p_test < 0.05 and circular_diff_deg(phase_train, phase_test) <= 60.0)

        comb_results.append(dict(
            period_days=float(P), artifact_label=bool(is_artifact),
            R_train=R_train, p_train=p_train, phase_train_deg=phase_train, detected=detected,
            R_test=R_test, p_test=p_test, phase_test_deg=phase_test, confirmed=confirmed,
        ))
        per_period_times.append(time.time() - t0p)

        # runtime guard: after a handful of periods, extrapolate total wall time incl.
        # test-confirm passes and halve remaining surrogate count if projected > 90 min
        if not guard_applied and i == 4:
            avg = float(np.mean(per_period_times))
            n_detected_so_far = sum(r["detected"] for r in comb_results)
            detect_rate = max(n_detected_so_far / (i + 1), 0.05)
            remaining = len(periods) - (i + 1)
            projected_remaining = avg * remaining * (1.0 + detect_rate)
            projected_total = (time.time() - t_start) + projected_remaining + 120  # + method-check buffer
            if projected_total > RUNTIME_BUDGET_S:
                n_surr = N_SURR_DEFAULT // 2
                guard_applied = True
                guard_note = (f"halved surrogates to {n_surr} after period {i+1}/{len(periods)}: "
                               f"projected total {projected_total/60:.1f} min > 90 min budget")
                print(f"  [runtime guard] {guard_note}")

        if (i + 1) % 10 == 0 or i == len(periods) - 1:
            print(f"  [{i+1}/{len(periods)}] P={P:.4f}d R_train={R_train:.4f} p_train={p_train:.4f} "
                  f"detected={detected} confirmed={confirmed}")

    # -------------------------------------------------------------- method check
    print("\nMethod check: original (non-declustered) catalog, train period, M>=%.1f..." % MAG_FLOOR)
    orig = load_catalog("SCSN_original_catalog.txt")
    orig_train = orig[orig.time < SPLIT]
    t_orig = orig_train.t_unix.to_numpy()
    o_t0, o_t1 = t_orig.min(), t_orig.max()
    print(f"  original-catalog train n={len(orig_train)} span={(o_t1-o_t0)/86400:.0f} d")

    method_check = {}
    mc_n_surr = n_surr  # honor same runtime guard state
    for P in (0.5, 1.0, 7.0):
        R, phase, p = rayleigh_test(t_orig, o_t0, o_t1, P, mc_n_surr, rng)
        method_check[str(P)] = dict(period_days=P, R=R, p=p, phase_deg=phase, fires=bool(p < 0.01))
        print(f"  P={P}d: R={R:.4f} p={p:.4f} fires={p < 0.01}")

    method_check_pass = method_check["1.0"]["fires"] and method_check["7.0"]["fires"]
    method_check["n_surr_used"] = mc_n_surr
    method_check["pass"] = bool(method_check_pass)
    if not method_check_pass:
        print("\n" + "=" * 70)
        print("METHOD CHECK FAILED: 1 d and/or 7 d artifact periods did NOT fire "
              "(p < 0.01) on the original catalog train period. Per protocol, the comb "
              "results below must NOT be interpreted as reliable positive evidence.")
        print("=" * 70 + "\n")
    else:
        print("Method check PASSED: 1 d and 7 d artifacts fire as expected (positive control OK).\n")

    confirmed_nonartifact = [r for r in comb_results if r["confirmed"] and not r["artifact_label"]]
    detections = [r for r in comb_results if r["detected"]]

    elapsed_s = time.time() - t_start
    out = dict(
        protocol="PATTERN_PROTOCOL.md EXP-F",
        generated_at_unix=time.time(),
        seed=SEED,
        n_surr_default=N_SURR_DEFAULT,
        runtime_guard=dict(applied=guard_applied, note=guard_note, n_surr_used_comb=n_surr,
                            elapsed_seconds=elapsed_s),
        train_n=len(train), test_n=len(test),
        train_span_days=float((train_t1 - train_t0) / 86400.0),
        test_span_days=float((test_t1 - test_t0) / 86400.0),
        method_check=method_check,
        comb=comb_results,
        confirmed_nonartifact=confirmed_nonartifact,
    )
    out_path = HERE / "results_exp_f.json"
    out_path.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"Wrote {out_path}")

    # -------------------------------------------------------------- figure
    MAPS.mkdir(exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))
    P_arr = np.array([r["period_days"] for r in comb_results])
    negp = -np.log10(np.array([r["p_train"] for r in comb_results]))
    is_art = np.array([r["artifact_label"] for r in comb_results])
    is_conf = np.array([r["confirmed"] for r in comb_results])
    is_det = np.array([r["detected"] for r in comb_results])

    ax.scatter(P_arr[~is_art & ~is_det], negp[~is_art & ~is_det], c="steelblue", s=18,
               label="not detected", zorder=2)
    ax.scatter(P_arr[is_det & ~is_art & ~is_conf], negp[is_det & ~is_art & ~is_conf], c="orange",
               s=45, label="train-detected, not confirmed", zorder=3)
    ax.scatter(P_arr[is_conf & ~is_art], negp[is_conf & ~is_art], c="crimson", s=70, marker="*",
               label="confirmed (non-artifact)", zorder=4)
    ax.scatter(P_arr[is_art], negp[is_art], c="grey", s=45, marker="x",
               label="artifact-labeled (0.5/1/7 d)", zorder=3)
    ax.axhline(-np.log10(0.01), color="black", linestyle="--", linewidth=1,
               label="detection threshold p=0.01")
    ax.set_xscale("log")
    ax.set_xlabel("Period (days)")
    ax.set_ylabel("-log10(p_train)")
    ax.set_title("EXP-F periodicity comb: TRAIN Rayleigh significance vs period")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig_path = MAPS / "exp_f_comb.png"
    fig.savefig(fig_path, dpi=150)
    print(f"Wrote {fig_path}")

    # -------------------------------------------------------------- console summary
    print("\n" + "=" * 70)
    print("EXP-F SUMMARY")
    print("=" * 70)
    print(f"Method check: {'PASSED' if method_check_pass else 'FAILED'}")
    print(f"TRAIN detections (p<0.01): {len(detections)}/{len(comb_results)}")
    for r in detections:
        tag = " [ARTIFACT]" if r["artifact_label"] else ""
        conf = "CONFIRMED" if r["confirmed"] else ("not confirmed" if r["p_test"] is not None else "")
        print(f"  P={r['period_days']:.4f} d  R_train={r['R_train']:.4f}  p_train={r['p_train']:.4f}"
              f"{tag}  test: p_test={r['p_test']}  {conf}")
    print(f"Confirmed non-artifact periods: {len(confirmed_nonartifact)}")
    for r in confirmed_nonartifact:
        print(f"  P={r['period_days']:.4f} d  phase_train={r['phase_train_deg']:.1f}deg "
              f"phase_test={r['phase_test_deg']:.1f}deg  p_train={r['p_train']:.4f} p_test={r['p_test']:.4f}")
    if guard_applied:
        print(f"Runtime guard applied: {guard_note}")
    print(f"Elapsed: {elapsed_s/60:.1f} min")


if __name__ == "__main__":
    main()
