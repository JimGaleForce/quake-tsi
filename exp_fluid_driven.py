"""P-1.3 THE FLUID-DRIVEN ARM: test where the physics says the effect must be LARGEST.

Priced 1. 137,480 events across 14 regions, on disk since the k034 round and never once
scored for tidal response.

WHY THIS IS THE HIGHEST-PRIOR ARM IN THE PROGRAM, and the argument is borrowed from
observations nobody disputes rather than from hope. Tidal modulation of TREMOR and
LOW-FREQUENCY EARTHQUAKES is established and replicated -- Cascadia, Japan, Parkfield,
multiple independent groups -- and it is not a one-percent effect there, it is order
unity. So the question was never "does the tide move faults". It demonstrably does. The
question is WHERE ALONG THE SPECTRUM the sensitivity turns off, and why there.

The physical parameter that separates tremor from ordinary seismicity is EFFECTIVE NORMAL
STRESS. Tremor sources sit at near-lithostatic pore pressure, so effective stress is small,
the rate-and-state product A*sigma is small, and a kilopascal of tidal stress is a large
perturbation. Ordinary tectonic seismicity sits at high effective stress where the same
kilopascal is negligible. **Geothermal and volcanic fields are the regular-seismicity
populations closest to the tremor end of that axis**: high pore pressure, active fluid
circulation, and at The Geysers, direct industrial injection.

So this arm does not ask "is there a tidal effect". It asks a sharper question with a
pre-declared answer shape: **does tidal sensitivity track effective stress?**

THE DIFFERENTIAL, DECLARED BEFORE THE RUN. The 14 k034 regions split into two groups of
seven on physical grounds fixed in advance, not on anything measured here:

  FLUID-DRIVEN   geysers, coso, long_valley, salton_brawley, lassen, yellowstone,
                 mono_west_nv_mina        (geothermal fields, volcanic centres,
                                           active hydrothermal systems)
  TECTONIC       san_jacinto, parkfield, mendocino, wasatch_slc, smith_valley_nv,
                 cedar_city_ut, little_skull_mtn   (ordinary tectonic faulting)

PREDICTION: tidal SENSITIVITY is larger in the fluid-driven group. The declared statistic
is the difference in MEAN z^2 -- a sign-free magnitude of sensitivity -- between the two
groups, and the null PERMUTES WHICH REGIONS CARRY THE FLUID LABEL, holding the 7/7 split
size fixed. That null is immune to every common-mode artifact by construction, because
relabelling cannot change the data; and the split itself cannot be tuned, because it is
published geology fixed before the run.

A SPECIFICATION ERROR IN THE FIRST VERSION, corrected and left on the record. The first
run declared this prediction but TESTED the difference of pooled SIGNED z, and returned
p = 0.0007. Those are different claims: a signed difference is a statement about preferred
PHASE, not about susceptibility, and the p was driven by the two groups leaning in
opposite directions (fluid +2.45, tectonic -3.57 on elong_cos.mean) rather than by one
being more sensitive. The signed-difference numbers are retained below as
`superseded_signed_difference` and must not be quoted.

THE ARTIFACT THAT WILL CERTAINLY BE PRESENT, NAMED IN ADVANCE. At M >= 1.5 the day/night
detection asymmetry is real and large: `exp_diurnal_discriminator.py` measured roughly 2
percent of events displaced into local night, centred on local midnight, dying above M2.5.
It WILL appear here. So the solar features are carried through the whole battery as an
ARTIFACT MONITOR -- they are expected to fire, and their firing validates that the arm is
sensitive -- while **the declared test is over the LUNAR AND TIDAL features only**. That
separation is fixed here, before the run, and the solar cells are excluded from the
declared family rather than excluded after seeing them.

NULL DISCIPLINE. These catalogues are NOT declustered, and `exp_nearcritical.py` showed
what a per-event null does to clustered data: max |z| = 26.27 of pure dependence, falling
to 2.51 when the null respected sequences. So the primary arm is DECLUSTERED with a
per-event dwell null, and a secondary arm keeps every event with a SEQUENCE-LEVEL null,
one rigid offset per sequence. Both are reported.

A SECOND ERROR, also corrected. The first run declustered at the program's standing
30 d / 150 km, which is calibrated for GLOBAL M>=5 tectonic catalogues. These regions span
20-40 km, so a 150 km radius deletes nearly every event that is not the first of its
swarm: 96,602 events collapsed to 6,321. The radius is now 25 km / 7 d, matched to the
size of the catalogues rather than to a global one, and the primary arm retains 14,405
events. The radius is a property of the REGION, not of the hypothesis, so setting it
correctly is not a researcher degree of freedom.

RESULT: the prediction is REFUTED, not merely unsupported. Mean z^2 is 1.084 (fluid) vs
1.250 (tectonic) in the primary and 1.122 vs 1.407 in the secondary -- fluid-driven regions
are if anything slightly LESS tidally sensitive, p = 0.814 and 0.831 in the predicted
direction. Both groups sit near the null expectation of z^2 = 1. Tidal sensitivity in
regular seismicity at M >= 1.5 does not track the effective-stress axis that separates
tremor from ordinary earthquakes.

Null half-window is one full fortnightly period so the spring-neap envelope is evenly
sampled. Exploration split only.

PRICED 1: one declared family over the lunar/tidal cells, one max-statistic, one p, plus
the declared fluid-minus-tectonic difference.
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

import exp_mass_screen as MS
import exp_world_harmonics as W

OUT_JSON = HERE / "results_fluid_driven.json"

FLUID = ("geysers", "coso", "long_valley", "salton_brawley", "lassen",
         "yellowstone", "mono_west_nv_mina")
TECTONIC = ("san_jacinto", "parkfield", "mendocino", "wasatch_slc",
            "smith_valley_nv", "cedar_city_ut", "little_skull_mtn")

MAG_MIN = 1.5
EXPLORE_FRAC = 0.70
# The program's standing declustering radius is 150 km, calibrated for GLOBAL M>=5
# tectonic catalogues. These regions span roughly 20-40 km, so 150 km deletes nearly
# every event that is not the first of its swarm: it collapsed 96,602 events to 6,321
# on the first run of this arm. A local radius is used instead, and it is a parameter
# of the REGION SIZE, not of the hypothesis.
DECLUSTER_KM_LOCAL = 25.0
DECLUSTER_DAYS_LOCAL = 7.0
NULL_HALF_WINDOW_DAYS = MS.FORTNIGHT_DAYS
S_NULL = 120
N_REPLICATES = 3000
REP_CHUNK = 500
RNG_SEED = 20260822
MIN_N = 200

# solar features are the ARTIFACT MONITOR and are excluded from the declared family
SOLAR = ("sun_sinel", "sun_sinel_abs", "sun_az", "sun_hourangle")
SCALARS = MS.SCALARS
ANGLES = MS.ANGLES


def is_solar(feature):
    return feature.split(".")[0] in SOLAR


def load(path):
    t, la, lo, mg = [], [], [], []
    for r in csv.DictReader(open(path, newline="", encoding="utf-8", errors="replace")):
        try:
            m = float(r["mag"])
            if m < MAG_MIN:
                continue
            ts = _dt.datetime.fromisoformat(r["time"].replace("Z", "+00:00"))
            # EPOCH INVARIANT: days since 1970-01-01Z, as MS.raw_series requires.
            t.append(ts.timestamp() / 86400.0)
            la.append(float(r["latitude"])); lo.append(float(r["longitude"]))
            mg.append(m)
        except (ValueError, KeyError, TypeError):
            continue
    t = np.asarray(t); o = np.argsort(t)
    t, la, lo, mg = t[o], np.asarray(la)[o], np.asarray(lo)[o], np.asarray(mg)[o]
    cut = t.min() + (t.max() - t.min()) * EXPLORE_FRAC
    k = t < cut
    return t[k], la[k], lo[k], mg[k], int((~k).sum())


def score_region(name, t, la, lo, rng, sequence_null, seq_idx=None):
    n = t.size
    n_units = (np.unique(seq_idx).size if sequence_null else n)
    obs_raw = MS.raw_series(t, la, lo)
    if sequence_null:
        off_u = rng.uniform(-NULL_HALF_WINDOW_DAYS, NULL_HALF_WINDOW_DAYS,
                            size=(n_units, S_NULL))
        off = off_u[seq_idx]
        idx_u = rng.integers(0, S_NULL, size=(n_units, N_REPLICATES))
        idx = idx_u[seq_idx]
    else:
        off = rng.uniform(-NULL_HALF_WINDOW_DAYS, NULL_HALF_WINDOW_DAYS, size=(n, S_NULL))
        idx = rng.integers(0, S_NULL, size=(n, N_REPLICATES))
    pool = MS.raw_series((t[:, None] + off).ravel(),
                         np.repeat(la, S_NULL), np.repeat(lo, S_NULL))
    pool = {k: v.reshape(n, S_NULL) for k, v in pool.items()}
    del off

    def nullmean(arr):
        out = np.empty(N_REPLICATES)
        for a in range(0, N_REPLICATES, REP_CHUNK):
            b = min(a + REP_CHUNK, N_REPLICATES)
            out[a:b] = np.take_along_axis(arr, idx[:, a:b], axis=1).mean(axis=0)
        return out

    res = {}
    for s in SCALARS:
        thr_hi = float(np.quantile(pool[s], 0.90))
        thr_lo = float(np.quantile(pool[s], 0.10))
        res["%s.mean" % s] = (obs_raw[s].mean(), nullmean(pool[s]))
        res["%s.top_decile" % s] = ((obs_raw[s] > thr_hi).mean(),
                                    nullmean((pool[s] > thr_hi).astype(np.float64)))
        res["%s.bot_decile" % s] = ((obs_raw[s] < thr_lo).mean(),
                                    nullmean((pool[s] < thr_lo).astype(np.float64)))
    for a in ANGLES:
        for mm in MS.ORDERS:
            gc = nullmean(np.cos(mm * pool[a]))
            gs = nullmean(np.sin(mm * pool[a]))
            res["%s.R%d" % (a, mm)] = (
                float(np.hypot(np.cos(mm * obs_raw[a]).mean(),
                               np.sin(mm * obs_raw[a]).mean())),
                np.hypot(gc, gs))
    out = {}
    for k, (o, nc) in res.items():
        mu, sd = float(nc.mean()), float(nc.std(ddof=1))
        if np.isfinite(sd) and sd > 0:
            out[k] = {"z": (float(o) - mu) / sd, "zn": (nc - mu) / sd}
    return out, n_units


def run_arm(declustered):
    rng = np.random.default_rng(RNG_SEED)
    regions, meta = {}, {}
    for path in sorted(glob.glob(str(HERE / "data" / "k034" / "*.csv"))):
        name = os.path.basename(path)[:-4]
        if name not in FLUID and name not in TECTONIC:
            continue
        t, la, lo, mg, n_hold = load(path)
        if t.size:
            print("EPOCH CHECK %-18s first event %s"
                  % (name, MS.iso_from_days(np.min(t))), flush=True)
        seq_idx = None
        if declustered:
            keep = _decluster_local(t, la, lo, mg)
            t, la, lo, mg = t[keep], la[keep], lo[keep], mg[keep]
        else:
            _, seq_idx = np.unique(_sequence_ids(t, la, lo, mg), return_inverse=True)
        if t.size < MIN_N:
            print("    %-22s n=%6d  SKIPPED" % (name, t.size), flush=True)
            continue
        r, nu = score_region(name, t, la, lo, rng, not declustered, seq_idx)
        regions[name] = r
        meta[name] = {"n": int(t.size), "units": int(nu),
                      "group": "fluid" if name in FLUID else "tectonic",
                      "n_holdout_reserved": n_hold}
        print("    %-22s n=%6d  units=%6d  %s" % (name, t.size, nu, meta[name]["group"]),
              flush=True)
    return regions, meta


def _decluster_local(t, lat, lon, mag):
    """Same rule as the program's standing declustering, at a radius matched to the
    size of these catalogues rather than to a global tectonic one."""
    keep, recent = [], []
    for ii in np.argsort(t):
        while recent and t[ii] - t[recent[0]] > DECLUSTER_DAYS_LOCAL:
            recent.pop(0)
        dep = False
        for jj in reversed(recent):
            if mag[jj] >= mag[ii] and W._great_circle_km(
                    lat[jj], lon[jj], lat[ii], lon[ii]) <= DECLUSTER_KM_LOCAL:
                dep = True
                break
        if not dep:
            keep.append(ii)
        recent.append(ii)
    return np.asarray(sorted(keep), dtype=int)


def group_label_permutation(regions, meta, n_perm=5000, seed=RNG_SEED):
    """THE CORRECTED DECLARED TEST, and the one that actually matches the prediction.

    The prediction was that tidal SENSITIVITY is larger in fluid-driven regions. The
    first version of this arm tested the difference of pooled SIGNED z, which is a
    statement about preferred PHASE, not about susceptibility -- its p = 0.0007 was
    driven by the two groups leaning in opposite directions rather than by one being
    more sensitive. That was a specification error and it is corrected here.

    Statistic: mean over lunar/tidal features of z^2 (a sign-free magnitude of
    sensitivity), fluid group minus tectonic group.

    Null: permute WHICH REGIONS are labelled fluid, holding the 7/7 split size fixed.
    This is the sharpest available null because it tests exactly the claim -- that the
    fluid/tectonic distinction is what organises the sensitivities -- and it is immune
    to any common-mode artifact, since relabelling cannot change the data.
    """
    names = sorted(regions)
    feats = sorted({f for r in regions.values() for f in r if not is_solar(f)})
    Z = np.array([[regions[n][f]["z"] if f in regions[n] else np.nan for f in feats]
                  for n in names])
    fluid_mask = np.array([meta[n]["group"] == "fluid" for n in names])
    k = int(fluid_mask.sum())
    if k < 2 or (len(names) - k) < 2:
        return None

    def stat(mask):
        a = np.nanmean(Z[mask] ** 2)
        b = np.nanmean(Z[~mask] ** 2)
        return float(a - b), float(a), float(b)

    obs, a_obs, b_obs = stat(fluid_mask)
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    idx = np.arange(len(names))
    for i in range(n_perm):
        m = np.zeros(len(names), dtype=bool)
        m[rng.choice(idx, size=k, replace=False)] = True
        null[i] = stat(m)[0]
    p_greater = float((np.sum(null >= obs) + 1) / (n_perm + 1))
    return {
        "statistic": "mean z^2 (fluid) - mean z^2 (tectonic), lunar/tidal features only",
        "mean_z2_fluid": a_obs, "mean_z2_tectonic": b_obs, "observed_difference": obs,
        "null_mean": float(null.mean()), "null_sd": float(null.std(ddof=1)),
        "null_p95": float(np.quantile(null, 0.95)),
        "p_one_sided_fluid_greater": p_greater,
        "n_permutations": n_perm, "n_regions": len(names), "n_features": len(feats),
        "reading": ("this tests the DECLARED prediction directly: is sensitivity larger "
                    "in the fluid-driven group? A signed-difference test cannot answer "
                    "that and the first version of this arm wrongly used one."),
    }


def _sequence_ids(t, lat, lon, mag):
    n = t.size
    seq = np.arange(n)
    recent = []
    for ii in np.argsort(t):
        while recent and t[ii] - t[recent[0]] > W.DECLUSTER_DAYS:
            recent.pop(0)
        for jj in reversed(recent):
            if mag[jj] >= mag[ii] and W._great_circle_km(
                    lat[jj], lon[jj], lat[ii], lon[ii]) <= W.DECLUSTER_KM:
                seq[ii] = seq[jj]
                break
        recent.append(ii)
    return seq


def analyse(regions, meta):
    feats = sorted({k for r in regions.values() for k in r})
    lunar = [f for f in feats if not is_solar(f)]
    solar = [f for f in feats if is_solar(f)]

    # cell max over the DECLARED (lunar/tidal) family only
    run = np.zeros(N_REPLICATES)
    cells = []
    for rn, r in regions.items():
        for f in lunar:
            if f in r:
                np.maximum(run, np.abs(r[f]["zn"]), out=run)
                cells.append({"region": rn, "group": meta[rn]["group"],
                              "feature": f, "n": meta[rn]["n"], "z": r[f]["z"]})
    cells.sort(key=lambda c: -abs(c["z"]))
    obs = abs(cells[0]["z"])
    p = float((np.sum(run >= obs) + 1) / (N_REPLICATES + 1))

    def pooled(group, flist):
        rows, prun = [], np.zeros(N_REPLICATES)
        best, where = 0.0, None
        store = {}
        for f in flist:
            zs, zns, ws = [], [], []
            for rn, r in regions.items():
                if meta[rn]["group"] != group or f not in r:
                    continue
                zs.append(r[f]["z"]); zns.append(r[f]["zn"])
                ws.append(np.sqrt(meta[rn]["units"]))
            if len(zs) < 3:
                continue
            w = np.asarray(ws); den = np.sqrt((w * w).sum())
            po = float((w * np.asarray(zs)).sum() / den)
            pn = (w[:, None] * np.stack(zns)).sum(axis=0) / den
            store[f] = (po, pn)
            np.maximum(prun, np.abs(pn), out=prun)
            rows.append({"feature": f, "pooled_z": po, "n_regions": len(zs)})
            if abs(po) > abs(best):
                best, where = po, f
        rows.sort(key=lambda r: -abs(r["pooled_z"]))
        return rows, prun, best, where, store

    frows, frun, fbest, fwhere, fstore = pooled("fluid", lunar)
    trows, trun, tbest, twhere, tstore = pooled("tectonic", lunar)
    srows, _, sbest, swhere, _ = pooled("fluid", solar)

    # THE DECLARED DIFFERENCE: fluid pooled z minus tectonic pooled z, per feature
    drun = np.zeros(N_REPLICATES)
    drows = []
    for f in lunar:
        if f not in fstore or f not in tstore:
            continue
        dn = fstore[f][1] - tstore[f][1]
        sd = float(dn.std(ddof=1))
        if not np.isfinite(sd) or sd <= 0:
            continue
        dz = ((fstore[f][0] - tstore[f][0]) - float(dn.mean())) / sd
        np.maximum(drun, np.abs((dn - dn.mean()) / sd), out=drun)
        drows.append({"feature": f, "fluid_z": fstore[f][0],
                      "tectonic_z": tstore[f][0], "diff_z": dz})
    drows.sort(key=lambda r: -abs(r["diff_z"]))
    dobs = abs(drows[0]["diff_z"]) if drows else float("nan")
    dp = (float((np.sum(drun >= dobs) + 1) / (N_REPLICATES + 1)) if drows
          else float("nan"))

    return {
        "DECLARED_TEST_group_label_permutation": group_label_permutation(regions, meta),
        "n_regions": len(regions), "n_events": sum(m["n"] for m in meta.values()),
        "n_units": sum(m["units"] for m in meta.values()),
        "per_region": meta,
        "declared_family_cells": len(cells),
        "cell_max_LUNAR_ONLY": {"observed_max_abs_z": obs,
                                "null_max_p95": float(np.quantile(run, 0.95)),
                                "p": p,
                                "where": "%s / %s" % (cells[0]["region"],
                                                      cells[0]["feature"])},
        "pooled_fluid": {"max_abs_pooled_z": abs(fbest), "where": fwhere,
                         "null_max_p95": float(np.quantile(frun, 0.95)),
                         "p": float((np.sum(frun >= abs(fbest)) + 1)
                                    / (N_REPLICATES + 1)), "top": frows[:8]},
        "pooled_tectonic": {"max_abs_pooled_z": abs(tbest), "where": twhere,
                            "null_max_p95": float(np.quantile(trun, 0.95)),
                            "p": float((np.sum(trun >= abs(tbest)) + 1)
                                       / (N_REPLICATES + 1)), "top": trows[:8]},
        "superseded_signed_difference": {
            "observed_max_abs_diff_z": dobs,
            "null_max_p95": float(np.quantile(drun, 0.95)), "p": dp,
            "top": drows[:10]},
        "artifact_monitor_SOLAR_pooled_fluid": {"max_abs_pooled_z": abs(sbest),
                                                "where": swhere, "top": srows[:5]},
        "top_15_cells": cells[:15],
    }


def main():
    out = {"arm": "P-1.3 fluid-driven vs tectonic tidal sensitivity",
           "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
           "priced_tests": 1, "holdout_touched": False,
           "declared_prediction": ("pooled tidal concentration LARGER in fluid-driven "
                                   "than in tectonic regions; the statistic is the "
                                   "difference of pooled z, calibrated against its own "
                                   "null; the 7/7 split is published geology and was "
                                   "fixed before the run"),
           "fluid_regions": list(FLUID), "tectonic_regions": list(TECTONIC),
           "solar_excluded_from_declared_family": list(SOLAR),
           "why_solar_excluded": ("the day/night detection asymmetry is real and was "
                                  "measured at ~2% of events in "
                                  "exp_diurnal_discriminator.py; solar features are "
                                  "carried as an ARTIFACT MONITOR whose firing shows the "
                                  "arm is sensitive, and are excluded from the declared "
                                  "family here BEFORE the run rather than after")}
    for label, dec in (("declustered_PRIMARY", True), ("full_sequence_null_SECONDARY", False)):
        print("\n=== %s ===" % label, flush=True)
        regions, meta = run_arm(dec)
        out[label] = analyse(regions, meta)
        a = out[label]
        print("  events %d  units %d  regions %d"
              % (a["n_events"], a["n_units"], a["n_regions"]))
        print("  CELL max |z| (lunar only) = %.3f at %s ; null 95th %.3f -> p = %.4f"
              % (a["cell_max_LUNAR_ONLY"]["observed_max_abs_z"],
                 a["cell_max_LUNAR_ONLY"]["where"],
                 a["cell_max_LUNAR_ONLY"]["null_max_p95"],
                 a["cell_max_LUNAR_ONLY"]["p"]))
        print("  pooled FLUID    %.3f at %-22s p = %.4f"
              % (a["pooled_fluid"]["max_abs_pooled_z"], a["pooled_fluid"]["where"],
                 a["pooled_fluid"]["p"]))
        print("  pooled TECTONIC %.3f at %-22s p = %.4f"
              % (a["pooled_tectonic"]["max_abs_pooled_z"],
                 a["pooled_tectonic"]["where"], a["pooled_tectonic"]["p"]))
        g = a.get("DECLARED_TEST_group_label_permutation")
        if g:
            print("  *** DECLARED TEST (group-label permutation, sign-free) ***")
            print("      mean z^2  fluid %.3f   tectonic %.3f   difference %+.3f"
                  % (g["mean_z2_fluid"], g["mean_z2_tectonic"],
                     g["observed_difference"]))
            print("      null mean %+.3f sd %.3f 95th %+.3f  ->  p(fluid greater) = %.4f"
                  % (g["null_mean"], g["null_sd"], g["null_p95"],
                     g["p_one_sided_fluid_greater"]))
        print("  [superseded, signed-difference] max |z| = %.3f ; null 95th %.3f -> p = %.4f"
              % (a["superseded_signed_difference"]["observed_max_abs_diff_z"],
                 a["superseded_signed_difference"]["null_max_p95"],
                 a["superseded_signed_difference"]["p"]))
        print("  [artifact monitor] solar pooled fluid %.3f at %s"
              % (a["artifact_monitor_SOLAR_pooled_fluid"]["max_abs_pooled_z"],
                 a["artifact_monitor_SOLAR_pooled_fluid"]["where"]))

    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
