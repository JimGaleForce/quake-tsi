"""THE AMPLITUDE ARM: does a BIGGER tide trigger more earthquakes? Priced 72.

Fable audit finding 2, and it is the largest hole in the program.

WHAT WAS NEVER TESTED. The 208-test world scan and the 221-test fault-relative scan
between them declare 16 statistics. Read them: twelve are ANGLES (principal bearing,
lunar azimuth, solar azimuth, at harmonic orders 1 to 4) and four are SIGN OR THRESHOLD
FLAGS (`areal < 0`, `rate < 0`, their conjunction, `|sin elev| > 0.5`). Every one is
invariant to the MAGNITUDE of the tidal stress. Multiply the entire tidal field by ten
and not one of the 16 numbers moves.

So the single most physically obvious hypothesis in the whole subject -- A BIGGER
TIDAL STRESS TRIGGERS MORE EARTHQUAKES -- has never been tested by this program at all.
The K-092/D-1 line tested tidal PHASE. This tests tidal AMPLITUDE, which is a different
question and the one a mechanism would actually predict: if the tide adds stress to a
near-critical fault, what matters is how much stress it adds, and that varies by a
factor of roughly three between neap and spring.

THE NULL, AND WHY ITS WINDOW CHANGED. The dwell-time null is waveform-matched exactly as
everywhere else in this program: each event's null draws come from its own neighbourhood
in time, so the null inherits the tide's own amplitude distribution and a statistic
cannot win merely by rediscovering the waveform. But the harmonics scan used a +/-10 day
window, which spans 20 days against a fortnightly period of 14.7653 days -- 1.35 cycles,
so the envelope is sampled UNEVENLY and the null is not amplitude-neutral. That does not
matter when every statistic is amplitude-blind. It matters entirely here.

So this arm sets the null half-window to exactly one fortnightly period, 14.7653 days,
sampling exactly two full cycles per event. The uniformity is not assumed: the null
amplitude marginal is compared against the long-run dwell marginal per region and the
worst-case discrepancy is reported in the artifact. If that discrepancy is large the
result is not to be read.

BOTH MULTIPLICITY CORRECTIONS ARE PAID, because the program now knows it needs both:
  - CELL max |z| over regions x statistics, for an effect concentrated in one place;
  - POOLED max |z| over statistics, for a weak effect present everywhere at once,
    which the cell max is provably blind to (see exp_world_pooled_calibration.py).

WHAT A POSITIVE WOULD MEAN, stated before the run. A survivor in E1/E2 says the rate of
M5+ earthquakes depends on tidal stress amplitude, which is a spring-neap forecast with a
14.77-day period and a known phase -- usable, checkable, and falsifiable on the holdout.
A survivor in E4 says loading and unloading are not symmetric, which is Jim's rope: pulled
from one direction then eased is not the same as never pulled. E5/E6 are the threshold
form: not a smooth dependence but a noise floor the tide has to clear.

THE HOLDOUT IS NOT TOUCHED. Exploration split only, via the same loader as every other
arm, so a survivor here has somewhere to be confirmed.
"""

from __future__ import annotations

import datetime as _dt
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import exp_world_harmonics as W
from engine import tidal_tensor as TT

OUT_JSON = HERE / "results_world_amplitude.json"

FORTNIGHT_DAYS = 14.765294       # Msf, the spring-neap period
NULL_HALF_WINDOW_DAYS = FORTNIGHT_DAYS      # exactly two full cycles per event
S_NULL_PER_EVENT = 200
N_REPLICATES = 5000
RNG_SEED = 20260822
RATE_HALF_MINUTES = 10.0

STATS = ("E1_areal_amp", "E2_shear_amp", "E3_rate_amp",
         "E4_signed_rate", "E5_amp_top_decile", "E6_rate_top_decile")

DECILE = 0.90


def amplitude_features(t_days, lat, lon):
    """Magnitude-bearing tidal quantities at arbitrary (time, site) triples."""
    jd = t_days + W.UNIX_EPOCH_JD
    st = TT.stress_tensor(jd, lat, lon, 0.0)
    areal = st["areal_strain"]
    shear = np.hypot(0.5 * (st["s_NN"] - st["s_EE"]), st["s_NE"])

    dt_d = RATE_HALF_MINUTES / 1440.0
    hi = TT.strain_tensor(jd + dt_d, lat, lon, 0.0)["areal_strain"]
    lo = TT.strain_tensor(jd - dt_d, lat, lon, 0.0)["areal_strain"]
    rate = (hi - lo) / (2.0 * RATE_HALF_MINUTES / 60.0)     # per hour

    return {"areal_abs": np.abs(areal), "shear": shear,
            "rate_abs": np.abs(rate), "rate_signed": rate}


def battery(f, thr_amp, thr_rate, axis=None):
    """The six declared statistics. Thresholds come from the region's own dwell pool,
    never from the events, so E5/E6 cannot be tuned by what the events happen to do."""
    m = (lambda v: float(v.mean())) if axis is None else (lambda v: v.mean(axis=axis))
    return {
        "E1_areal_amp": m(f["areal_abs"]),
        "E2_shear_amp": m(f["shear"]),
        "E3_rate_amp": m(f["rate_abs"]),
        "E4_signed_rate": m(f["rate_signed"]),
        "E5_amp_top_decile": m((f["areal_abs"] > thr_amp).astype(np.float64)),
        "E6_rate_top_decile": m((f["rate_abs"] > thr_rate).astype(np.float64)),
    }


def dwell_thresholds(t, lat, lon, rng):
    """Top-decile cutoffs from a dense uniform sweep of the region's own span."""
    n = min(4000, max(400, t.size * 4))
    tt = rng.uniform(t.min(), t.max(), size=n)
    # Sites are sampled JOINTLY. Drawing lat and lon independently would pair a
    # Turkish latitude with an unrelated longitude and score the tide at a site
    # that does not exist; the first version of this function did exactly that.
    si = rng.integers(0, t.size, size=n)
    f = amplitude_features(tt, lat[si], lon[si])
    return (float(np.quantile(f["areal_abs"], DECILE)),
            float(np.quantile(f["rate_abs"], DECILE)),
            f["areal_abs"])


def run_region(name, t, lat, lon, rng, verbose=True):
    n = t.size
    if n < 30:
        return None
    thr_amp, thr_rate, dwell_amp = dwell_thresholds(t, lat, lon, rng)

    obs = battery(amplitude_features(t, lat, lon), thr_amp, thr_rate)

    off = rng.uniform(-NULL_HALF_WINDOW_DAYS, NULL_HALF_WINDOW_DAYS,
                      size=(n, S_NULL_PER_EVENT))
    nf = amplitude_features((t[:, None] + off).ravel(),
                            np.repeat(lat, S_NULL_PER_EVENT),
                            np.repeat(lon, S_NULL_PER_EVENT))
    nf = {k: v.reshape(n, S_NULL_PER_EVENT) for k, v in nf.items()}

    # Is the null amplitude-neutral? Compare its marginal to the long-run dwell marginal.
    q = np.linspace(0.05, 0.95, 19)
    disc = float(np.max(np.abs(np.quantile(nf["areal_abs"].ravel(), q)
                               - np.quantile(dwell_amp, q))
                        / np.quantile(dwell_amp, q)))

    idx = rng.integers(0, S_NULL_PER_EVENT, size=(n, N_REPLICATES))
    sel = {k: np.take_along_axis(v, idx, axis=1) for k, v in nf.items()}
    nb = battery(sel, thr_amp, thr_rate, axis=0)

    z, cols = {}, {}
    for k in STATS:
        arr = np.asarray(nb[k], dtype=np.float64)
        mu, sd = float(arr.mean()), float(arr.std(ddof=1))
        z[k] = {"observed": float(obs[k]), "null_mean": mu, "null_sd": sd,
                "z": (float(obs[k]) - mu) / sd if sd > 0 else float("nan")}
        cols[k] = (arr - mu) / (sd if sd > 0 else 1e-300)
    if verbose:
        top = max(STATS, key=lambda k: abs(z[k]["z"]))
        print("  %-18s n=%5d  largest |z| = %5.2f (%-18s)  null-neutrality %.3f"
              % (name, n, abs(z[top]["z"]), top, disc), flush=True)
    return {"n": int(n), "per_statistic": z, "_null_cols": cols,
            "null_amplitude_neutrality_worst_rel_dev": disc,
            "threshold_areal_abs": thr_amp, "threshold_rate_abs": thr_rate}


def main():
    declustered = os.environ.get("WORLD_DECLUSTER", "1") == "1"
    rng = np.random.default_rng(RNG_SEED)
    regions = {}
    for path in sorted(glob.glob(str(HERE / "data" / "comcat_world" / "*.csv"))):
        name = os.path.basename(path)[:-4]
        t, la, lo, _h, _f = W.load_region(path, declustered=declustered)
        r = run_region(name, t, la, lo, rng)
        if r is not None:
            regions[name] = r

    names = sorted(regions)
    cells = [(rn, k) for rn in names for k in STATS]

    # ---- correction 1: cell max over regions x statistics
    cellcols = np.stack([np.abs(regions[rn]["_null_cols"][k]) for rn, k in cells])
    cell_null = cellcols.max(axis=0)
    cell_obs = max(abs(regions[rn]["per_statistic"][k]["z"]) for rn, k in cells)
    cell_p = float((np.sum(cell_null >= cell_obs) + 1) / (N_REPLICATES + 1))
    cell_where = max(cells, key=lambda c: abs(regions[c[0]]["per_statistic"][c[1]]["z"]))

    # ---- correction 2: pooled max over statistics
    w = np.array([np.sqrt(regions[rn]["n"]) for rn in names])
    den = np.sqrt((w * w).sum())
    pooled_obs, pooled_null = {}, {}
    for k in STATS:
        zo = np.array([regions[rn]["per_statistic"][k]["z"] for rn in names])
        pooled_obs[k] = float((w * zo).sum() / den)
        zn = np.stack([regions[rn]["_null_cols"][k] for rn in names])
        pooled_null[k] = (w[:, None] * zn).sum(axis=0) / den
    pool_null_max = np.stack([np.abs(pooled_null[k]) for k in STATS]).max(axis=0)
    pool_obs = max(abs(pooled_obs[k]) for k in STATS)
    pool_p = float((np.sum(pool_null_max >= pool_obs) + 1) / (N_REPLICATES + 1))
    pool_where = max(STATS, key=lambda k: abs(pooled_obs[k]))

    worst_neutral = max(regions[rn]["null_amplitude_neutrality_worst_rel_dev"]
                        for rn in names)
    for rn in names:
        regions[rn].pop("_null_cols", None)

    out = {
        "arm": "amplitude / envelope / rate (Fable audit finding 2)",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "declustered": declustered,
        "priced_tests": len(cells) + len(STATS),
        "the_gap": ("all 16 statistics in the committed world scans are angles or sign "
                    "flags, hence invariant to tidal stress MAGNITUDE. Scale the tidal "
                    "field by 10 and not one of them moves. Amplitude was never tested."),
        "null": {"kind": "waveform-matched dwell-time",
                 "half_window_days": NULL_HALF_WINDOW_DAYS,
                 "why_this_window": ("exactly one fortnightly period each side, so each "
                                     "event's null samples two full spring-neap cycles "
                                     "and the null is amplitude-neutral; the +/-10 day "
                                     "window of the harmonics scan spans only 1.35 "
                                     "cycles, which is irrelevant to angle statistics "
                                     "and fatal to amplitude ones"),
                 "worst_amplitude_neutrality_rel_dev": worst_neutral},
        "n_regions": len(names), "regions": names,
        "n_events": int(sum(regions[rn]["n"] for rn in names)),
        "per_region": regions,
        "cell_max_over_regions_x_statistics": {
            "observed_max_abs_z": cell_obs,
            "null_max_p95": float(np.quantile(cell_null, 0.95)),
            "p": cell_p, "where": "%s / %s" % cell_where, "n_cells": len(cells)},
        "pooled_max_over_statistics": {
            "observed_max_abs_pooled_z": pool_obs,
            "null_max_p95": float(np.quantile(pool_null_max, 0.95)),
            "p": pool_p, "where": pool_where, "n_statistics": len(STATS),
            "per_statistic": {k: pooled_obs[k] for k in STATS}},
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    print("\n" + "=" * 76)
    print("AMPLITUDE ARM  (%s, %d regions, %d events, %d priced tests)"
          % ("declustered" if declustered else "FULL", len(names),
             out["n_events"], out["priced_tests"]))
    print("  null amplitude-neutrality, worst relative deviation: %.4f" % worst_neutral)
    print("  CELL   max |z| = %.3f at %s ; null 95th %.3f ; p = %.4f"
          % (cell_obs, "%s / %s" % cell_where,
             out["cell_max_over_regions_x_statistics"]["null_max_p95"], cell_p))
    print("  POOLED max |z| = %.3f at %s ; null 95th %.3f ; p = %.4f"
          % (pool_obs, pool_where,
             out["pooled_max_over_statistics"]["null_max_p95"], pool_p))
    print("\n  pooled z by statistic:")
    for k in STATS:
        print("    %-20s %+7.3f" % (k, pooled_obs[k]))
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
