"""THE MASS SCREEN: thousands of untested candidates, ONE declared family. Priced 1.

Jim's directive: "can we start testing a TON of untested possibilities already?" Yes,
and this is the honest way to do it.

WHY THIS IS AFFORDABLE, WHICH IS THE WHOLE POINT. A declared test is expensive because
it spends multiplicity budget. But a max-statistic correction over a SHARED null ensemble
costs the same whether the family holds ten cells or ten thousand: the null distribution
of `max |z|` is computed over exactly the same cells the observed max ranges over, so the
correction is exact by construction and the family p is honest at any width. Ten thousand
candidates therefore cost ONE declared family. What they cost instead is POWER -- a wider
family raises the null maximum, so a real effect has to be bigger to clear it. That is the
real trade and it is stated here rather than hidden:

    null max |z| over ~60 cells    95th percentile ~ 3.3
    null max |z| over ~3000 cells  95th percentile ~ 4.2

So a screen this wide costs roughly 0.9 sigma of sensitivity against a single-cell effect,
and buys coverage of a hypothesis space nobody has looked at. Against a concentrated
effect that is a bad trade; against "we do not know where to look", it is the only trade.

THE RULE, FROZEN BEFORE THE RUN AND NOT NEGOTIABLE AFTERWARDS.
  1. The family p is `P(max |z| over EVERY cell)`. Every cell. No cell is dropped after
     the fact for being implausible, and no sub-family is quoted as if it had been
     declared separately.
  2. A cell whose |z| clears the null 95th percentile is a CANDIDATE, not a finding. Its
     only privilege is the right to be re-run on the holdout.
  3. The holdout is not touched here. Exploration split only, same loader, same cutoff.
  4. A pooled-across-regions max over statistics is reported alongside, because the cell
     max is provably blind to a weak effect present everywhere at once.

WHAT IS IN THE SPACE, and it is deliberately full of things this program has never scored.
Twelve scalar series and six angular ones, crossed with three scalar statistics and four
harmonic orders, crossed with five subpopulations, crossed with every region:

  SCALARS   areal strain (signed and absolute); horizontal shear magnitude; d(areal)/dt
            (signed and absolute); lunar and solar sin(elevation), signed and absolute;
            lunar distance; the SPRING-NEAP ENVELOPE cos(elongation) and |cos(elongation)|
  ANGLES    principal stress bearing; lunar and solar azimuth; LUNAR AND SOLAR HOUR ANGLE
            -- local lunar time, the natural clock of the semidiurnal tide, and never once
            scored by this program -- and lunar elongation, the synodic phase
  STATS     mean; top-decile occupancy; bottom-decile occupancy; harmonic resultant R_m
  SUBPOPS   all; shallow (< 70 km); deep (>= 70 km); M >= 5.5; M < 5.5

The subpopulations are the point as much as the features are. Every scan this program has
run tested MAIN EFFECTS on whole catalogues, and a conditional effect -- present only in a
near-critical or a particular-depth subpopulation -- is invisible to a main-effect test by
construction while being exactly what a threshold mechanism predicts.

The null is the program's standing waveform-matched dwell-time construction, with the
half-window set to one fortnightly period so the spring-neap envelope is sampled evenly
(the reason is in exp_world_amplitude.py and it applies to every amplitude-like feature
here). Subpopulation cells reuse the parent region's draws restricted to their own rows,
so the null is matched event-by-event within every subpopulation too.

PRICED 1: one declared family, one max-statistic, one p.
"""

from __future__ import annotations

import csv
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

OUT_JSON = HERE / "results_mass_screen.json"

FORTNIGHT_DAYS = 14.765294
NULL_HALF_WINDOW_DAYS = FORTNIGHT_DAYS
S_NULL_PER_EVENT = 200
N_REPLICATES = 4000
RNG_SEED = 20260822
RATE_HALF_MINUTES = 10.0
ORDERS = (1, 2, 3, 4)
MIN_SUBPOP = 30
DEEP_KM = 70.0
MAG_SPLIT = 5.5

SCALARS = ("areal", "areal_abs", "shear", "rate", "rate_abs",
           "moon_sinel", "moon_sinel_abs", "sun_sinel", "sun_sinel_abs",
           "moon_dist", "elong_cos", "elong_cos_abs")
ANGLES = ("bearing", "moon_az", "sun_az", "moon_hourangle", "sun_hourangle", "elongation")
SCALAR_STATS = ("mean", "top_decile", "bot_decile")


def load_region_full(path):
    """Same freeze as exp_world_harmonics.load_region, but keeps depth and magnitude.

    The cutoff, magnitude floor and declustering rule are imported rather than restated
    so this screen cannot silently drift from the arms it will be compared against.
    """
    cut = W.explore_cutoff()
    t, la, lo, mg, dp = [], [], [], [], []
    n_holdout = 0
    for r in csv.DictReader(open(path, newline="", encoding="utf-8")):
        try:
            m = float(r["mag"])
            if m < W.MAG_MIN:
                continue
            ts = _dt.datetime.fromisoformat(r["time"].replace("Z", "+00:00"))
            if ts >= cut:
                n_holdout += 1
                continue
            # EPOCH INVARIANT: days since 1970-01-01Z (see raw_series below,
            # which adds W.UNIX_EPOCH_JD). Identical to W.load_region.
            t.append(ts.timestamp() / 86400.0)
            la.append(float(r["latitude"]))
            lo.append(float(r["longitude"]))
            mg.append(m)
            dp.append(float(r["depth"]) if r.get("depth") not in (None, "") else 0.0)
        except (ValueError, KeyError, TypeError):
            continue
    t = np.asarray(t); la = np.asarray(la); lo = np.asarray(lo)
    mg = np.asarray(mg); dp = np.asarray(dp)
    keep = W.decluster(t, la, lo, mg)
    return t[keep], la[keep], lo[keep], mg[keep], dp[keep], n_holdout


EPOCH_1970 = _dt.datetime(1970, 1, 1, tzinfo=_dt.timezone.utc)


def iso_from_days(t_days):
    """ISO-8601 UTC string for a days-since-1970 value (the EPOCH INVARIANT base)."""
    return (EPOCH_1970 + _dt.timedelta(days=float(t_days))).isoformat()


def assert_epoch(t_days, expect_year, label):
    """FAILURE-FIRST epoch guard, run in every arm that feeds `raw_series`.

    Reconstructs the first event date from `t` and refuses to continue if it is not the
    year the catalogue actually starts in. The 9,131-day defect corrected on 2026-09-02
    would have shown up here as, e.g., QTM starting in 1982 instead of 2008.
    """
    iso = iso_from_days(np.min(t_days))
    print("EPOCH CHECK %-22s first event %s (expect year %d)"
          % (label, iso, expect_year), flush=True)
    got = int(iso[:4])
    if got != expect_year:
        raise SystemExit("EPOCH INVARIANT VIOLATED: %s first event reconstructs to %s, "
                         "expected year %d" % (label, iso, expect_year))


def raw_series(t_days, lat, lon):
    """Every base series at arbitrary (time, site) triples, computed once."""
    jd = t_days + W.UNIX_EPOCH_JD
    st = TT.stress_tensor(jd, lat, lon, 0.0)
    dt_d = RATE_HALF_MINUTES / 1440.0
    hi = TT.strain_tensor(jd + dt_d, lat, lon, 0.0)["areal_strain"]
    lo_ = TT.strain_tensor(jd - dt_d, lat, lon, 0.0)["areal_strain"]
    moon = TT.body_direction(jd, lat, lon, 0.0, "moon")
    sun = TT.body_direction(jd, lat, lon, 0.0, "sun")

    sm, ss = moon["sin_elevation"], sun["sin_elevation"]
    cm = np.sqrt(np.clip(1.0 - sm * sm, 0.0, 1.0))
    cs = np.sqrt(np.clip(1.0 - ss * ss, 0.0, 1.0))
    d2r = np.pi / 180.0
    daz = (moon["azimuth_deg"] - sun["azimuth_deg"]) * d2r
    cos_elong = np.clip(sm * ss + cm * cs * np.cos(daz), -1.0, 1.0)

    areal = st["areal_strain"]
    rate = (hi - lo_) / (2.0 * RATE_HALF_MINUTES / 60.0)
    return {
        "areal": areal, "areal_abs": np.abs(areal),
        "shear": np.hypot(0.5 * (st["s_NN"] - st["s_EE"]), st["s_NE"]),
        "rate": rate, "rate_abs": np.abs(rate),
        "moon_sinel": sm, "moon_sinel_abs": np.abs(sm),
        "sun_sinel": ss, "sun_sinel_abs": np.abs(ss),
        "moon_dist": moon["distance_m"],
        "elong_cos": cos_elong, "elong_cos_abs": np.abs(cos_elong),
        # angles, radians
        "bearing": np.mod(TT.principal_bearing(st["e_NN"], st["e_EE"], st["e_NE"]) * d2r,
                          2.0 * np.pi),
        "moon_az": moon["azimuth_deg"] * d2r,
        "sun_az": sun["azimuth_deg"] * d2r,
        "moon_hourangle": moon["hour_angle_rad"],
        "sun_hourangle": sun["hour_angle_rad"],
        # elongation runs 0..pi (0 = new moon, pi = full); DOUBLED so the two spring
        # states are identified, which is the physically correct axial folding
        "elongation": np.mod(2.0 * np.arccos(cos_elong), 2.0 * np.pi),
    }


def subpops(mag, dep):
    yield "all", np.ones(mag.shape, dtype=bool)
    yield "shallow_lt70km", dep < DEEP_KM
    yield "deep_ge70km", dep >= DEEP_KM
    yield "M_ge5.5", mag >= MAG_SPLIT
    yield "M_lt5.5", mag < MAG_SPLIT


def main():
    rng = np.random.default_rng(RNG_SEED)
    running_max = np.zeros(N_REPLICATES)
    records = []
    pooled_obs, pooled_null, pooled_w = {}, {}, {}
    n_events_total = 0
    region_names = []

    for path in sorted(glob.glob(str(HERE / "data" / "comcat_world" / "*.csv"))):
        name = os.path.basename(path)[:-4]
        t, la, lo, mg, dp, _hold = load_region_full(path)
        if t.size:
            print("EPOCH CHECK %-18s first event %s"
                  % (name, iso_from_days(np.min(t))), flush=True)
        n = t.size
        if n < MIN_SUBPOP:
            print("  %-18s n=%4d  SKIPPED (below %d)" % (name, n, MIN_SUBPOP), flush=True)
            continue
        region_names.append(name)
        n_events_total += n

        obs_raw = raw_series(t, la, lo)
        off = rng.uniform(-NULL_HALF_WINDOW_DAYS, NULL_HALF_WINDOW_DAYS,
                          size=(n, S_NULL_PER_EVENT))
        pool = raw_series((t[:, None] + off).ravel(),
                          np.repeat(la, S_NULL_PER_EVENT),
                          np.repeat(lo, S_NULL_PER_EVENT))
        pool = {k: v.reshape(n, S_NULL_PER_EVENT) for k, v in pool.items()}
        idx = rng.integers(0, S_NULL_PER_EVENT, size=(n, N_REPLICATES))

        masks = [(sn, m) for sn, m in subpops(mg, dp) if int(m.sum()) >= MIN_SUBPOP]

        def emit(cell, sub_name, mask, obs_val, null_col):
            mu, sd = float(null_col.mean()), float(null_col.std(ddof=1))
            if not np.isfinite(sd) or sd <= 0:
                return
            z = (float(obs_val) - mu) / sd
            zn = (null_col - mu) / sd
            np.maximum(running_max, np.abs(zn), out=running_max)
            key = "%s|%s" % (sub_name, cell)
            records.append({"region": name, "subpop": sub_name, "feature": cell,
                            "n": int(mask.sum()), "z": z})
            pooled_obs.setdefault(key, []).append(z)
            pooled_null.setdefault(key, []).append(zn)
            pooled_w.setdefault(key, []).append(np.sqrt(float(mask.sum())))

        for s in SCALARS:
            g = np.take_along_axis(pool[s], idx, axis=1)          # n x R
            thr_hi = np.quantile(pool[s], 0.90)
            thr_lo = np.quantile(pool[s], 0.10)
            hi_p = (pool[s] > thr_hi).astype(np.float64)
            lo_p = (pool[s] < thr_lo).astype(np.float64)
            g_hi = np.take_along_axis(hi_p, idx, axis=1)
            g_lo = np.take_along_axis(lo_p, idx, axis=1)
            o = obs_raw[s]
            for sn, m in masks:
                emit("%s.mean" % s, sn, m, o[m].mean(), g[m].mean(axis=0))
                emit("%s.top_decile" % s, sn, m,
                     (o[m] > thr_hi).mean(), g_hi[m].mean(axis=0))
                emit("%s.bot_decile" % s, sn, m,
                     (o[m] < thr_lo).mean(), g_lo[m].mean(axis=0))
            del g, g_hi, g_lo, hi_p, lo_p

        for a in ANGLES:
            th_p, th_o = pool[a], obs_raw[a]
            for mm in ORDERS:
                gc = np.take_along_axis(np.cos(mm * th_p), idx, axis=1)
                gs = np.take_along_axis(np.sin(mm * th_p), idx, axis=1)
                for sn, m in masks:
                    emit("%s.R%d" % (a, mm), sn, m,
                         np.hypot(np.cos(mm * th_o[m]).mean(),
                                  np.sin(mm * th_o[m]).mean()),
                         np.hypot(gc[m].mean(axis=0), gs[m].mean(axis=0)))
                del gc, gs
        print("  %-18s n=%5d  subpops=%d  cells so far %d"
              % (name, n, len(masks), len(records)), flush=True)

    # ---------------- family verdict
    obs_max = max(abs(r["z"]) for r in records)
    fam_p = float((np.sum(running_max >= obs_max) + 1) / (N_REPLICATES + 1))
    records.sort(key=lambda r: -abs(r["z"]))
    crit95 = float(np.quantile(running_max, 0.95))
    candidates = [r for r in records if abs(r["z"]) > crit95]

    # ---------------- pooled across regions, per (subpop, feature)
    pool_running = np.zeros(N_REPLICATES)
    pool_obs_best, pool_where = 0.0, None
    pooled_rows = []
    for key, zs in pooled_obs.items():
        if len(zs) < 3:
            continue
        w = np.asarray(pooled_w[key])
        den = np.sqrt((w * w).sum())
        po = float((w * np.asarray(zs)).sum() / den)
        pn = (w[:, None] * np.stack(pooled_null[key])).sum(axis=0) / den
        np.maximum(pool_running, np.abs(pn), out=pool_running)
        pooled_rows.append({"key": key, "n_regions": len(zs), "pooled_z": po})
        if abs(po) > abs(pool_obs_best):
            pool_obs_best, pool_where = po, key
    pool_p = float((np.sum(pool_running >= abs(pool_obs_best)) + 1) / (N_REPLICATES + 1))
    pooled_rows.sort(key=lambda r: -abs(r["pooled_z"]))

    out = {
        "arm": "MASS SCREEN: wide candidate space, one declared family",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 1,
        "why_priced_one": ("a max-statistic over a shared null ensemble is exact at any "
                           "family width; the cost of width is POWER, not budget, and "
                           "the power cost is reported as the null 95th percentile"),
        "declustered": True, "holdout_touched": False,
        "n_regions": len(region_names), "regions": region_names,
        "n_events": n_events_total,
        "n_cells": len(records),
        "null": {"kind": "waveform-matched dwell-time",
                 "half_window_days": NULL_HALF_WINDOW_DAYS,
                 "draws_per_event": S_NULL_PER_EVENT,
                 "replicates": N_REPLICATES},
        "space": {"scalars": list(SCALARS), "angles": list(ANGLES),
                  "scalar_stats": list(SCALAR_STATS), "harmonic_orders": list(ORDERS),
                  "subpopulations": ["all", "shallow_lt70km", "deep_ge70km",
                                     "M_ge5.5", "M_lt5.5"]},
        "cell_max_family": {
            "observed_max_abs_z": obs_max,
            "null_max_p95": crit95,
            "null_max_p50": float(np.quantile(running_max, 0.50)),
            "p": fam_p,
            "where": "%s / %s / %s" % (records[0]["region"], records[0]["subpop"],
                                       records[0]["feature"]),
            "n_candidates_over_p95": len(candidates)},
        "pooled_max_family": {
            "observed_max_abs_pooled_z": abs(pool_obs_best),
            "null_max_p95": float(np.quantile(pool_running, 0.95)),
            "p": pool_p, "where": pool_where, "n_keys": len(pooled_rows)},
        "top_30_cells": records[:30],
        "top_15_pooled": pooled_rows[:15],
        "candidates_for_holdout": candidates[:50],
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    print("\n" + "=" * 78)
    print("MASS SCREEN: %d cells, %d regions, %d declustered events, 1 priced family"
          % (len(records), len(region_names), n_events_total))
    print("  CELL   max |z| = %.3f at %s" % (obs_max, out["cell_max_family"]["where"]))
    print("         null max: median %.3f, 95th %.3f   ->  FAMILY p = %.4f"
          % (out["cell_max_family"]["null_max_p50"], crit95, fam_p))
    print("  POOLED max |z| = %.3f at %s ; null 95th %.3f  ->  p = %.4f"
          % (abs(pool_obs_best), pool_where,
             out["pooled_max_family"]["null_max_p95"], pool_p))
    print("\n  top 10 cells by |z|:")
    for r in records[:10]:
        print("    %+6.2f  %-18s %-14s %-24s n=%d"
              % (r["z"], r["region"], r["subpop"], r["feature"], r["n"]))
    print("\n  top 5 pooled:")
    for r in pooled_rows[:5]:
        print("    %+6.2f  %-40s (%d regions)"
              % (r["pooled_z"], r["key"], r["n_regions"]))
    print("\n  candidates over the null 95th (%.3f): %d" % (crit95, len(candidates)))
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
