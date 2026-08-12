"""K-069-min -- THE PROSPECTIVE COMMITMENT LOG (frozen emitter; commits, scores nothing).

K-069 pre-registers, with a hash, forward-looking (place, time) windows and scores the
realised seismicity in them against a frozen ETAS. Its full form needs K-059's
*phase-resolved* surface-wave ephemeris T(x,t) -- which does not exist on disk. What
DOES exist, and what this file freezes, is the entry's own reduced form: K-059-min,
"an amplitude-only ephemeris (no phase)" (Popper, P4-5 rank 3 / P5-10 rank 4).

WHAT IS COMMITTED TODAY (and it is a rule + code, not a list of windows, because the
windows cannot exist before their source events do -- that is the entry's own design):

  TRIGGER      any global M >= 7.0.
  DEADLINE     the window set is emitted and hash-appended within 1 h of origin time,
               and always before the first committed window opens (the nearest
               admissible window opens at 1000 km / 4.5 km/s = 62 min after origin).
  EMISSION     the top TOP_N = 100 (0.5 deg cell x 3 h window) cells of the ENVELOPE
               exposure field over the following 10 days, ranked by peak nominal
               dynamic stress, computed by `envelope_exposure()` in this file.
  EXCLUSION    cells within MIN_DELTA_KM = 1000 km of ANY contributing source are
               struck out before ranking -- the source's own aftershock zone and the
               near field are not what this family is about.
  OVERLAY      contributing sources = the trigger plus every catalogued M >= 6.0 in
               the preceding 10 days; exposures ADD. This is the ENVELOPE leg of
               K-060's coherent-minus-envelope comparison, and it is the only leg an
               amplitude-only ephemeris can produce.
  AMPLITUDE    van der Elst & Brodsky (2010) eq.(6), exactly as already used by this
               repository in exp_k034_landers_control.py:
                   log10 A20[um] = M - 1.66 log10(D[deg]) - 2
                   V = 2*pi*A20/T_SW ;  sigma = G_SHEAR * V / C_PHASE
               with G_SHEAR = 30 GPa, C_PHASE = 3500 m/s, T_SW = 20 s. Nominal: it is
               a ranking functional, and the ranking is what is committed.
  TIMING       group-velocity band U in [3.0, 4.5] km/s, peak at 3.7 km/s. A cell's
               exposure is assigned to the 3 h bins overlapping [D/4.5, D/3.0],
               weighted by overlap fraction.
  DOMAIN       the 0.5 deg active-cell domain of the engine design over
               data/comcat_world -- identical to the K-080 census domain, one file,
               one discipline, 13 boxes and NOT the globe.
  SCORING      24 months from the commitment date. Observed M >= 4.5 with epicentre in
               a committed cell and origin time in the committed 3 h bin, versus the
               expectation of the frozen ETAS (the K-080 census parameters) integrated
               over the identical cell-bins with history up to each bin's start.
               Statistic: rate ratio, and bits/event. Unit = committed window;
               N_eff = number of distinct trigger events, NOT the number of windows.
  PASS/FAIL    declared here and nowhere later: rate-ratio CI (on N_eff) excluding 1
               from above = the envelope arm carries prospective information;
               CI containing 1 = a clean prospective null with the committed windows
               behind it. Failure is publishable and is the likelier outcome.

WHAT IS **NOT** COMMITTED, AND WHY (NEEDS-INFRASTRUCTURE, stated so no later reading
can pretend it was committed):
  * The CONSTRUCTIVE-maxima arm. "Constructive" is a statement about phase; an
    amplitude-only ephemeris has none. K-059's phase-resolved T(x,t) is unbuilt.
  * The paired DESTRUCTIVE control set (K-061). It is defined by phase and is
    therefore unbuilt for the same reason. K-069's stated success rule requires BOTH
    ("rate ratio CI excluding 1 in the constructive set AND <= 1 in the destructive
    set") and that rule CANNOT be scored on this log. What this log can be scored on
    is the envelope arm alone, as a one-armed prospective test.
  * Interpretation remains gated on K-034, whose PASS is certified only at >= ~34 kPa
    (k034_report.py). Most committed windows will sit far below that. The commitment
    is still worth making today -- committing costs nothing and un-committing is
    impossible -- but a positive here is not licensed as dynamic triggering until the
    licence reaches the amplitude, and a null is recorded provisional-pending-K-034.

Usage:
    python -u k069_emit.py --init                     # write/refresh the commitment header
    python -u k069_emit.py --demo                     # retrospective DEMO, NOT SCORED
    python -u k069_emit.py --emit --time <iso> --lat <f> --lon <f> --mag <f> [--id <s>]

Deterministic: no randomness. SEED declared and unused. Ties broken by index.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import platform
import sys
import time

import numpy as np
import pandas as pd

from engine import __version__ as ENGINE_VERSION
from engine import design

SEED = 20260811                      # declared; no stochastic step exists in this file

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = "data/comcat_world"
LOG = os.path.join(HERE, "results_k069_prospective_log.json")
DEMO = os.path.join(HERE, "results_k069_demo_emission.json")

# ---- frozen protocol constants (S-9: one value each, declared before any trigger) --
TRIGGER_MAG = 7.0
CONTRIB_MAG = 6.0
HORIZON_DAYS = 10.0
BIN_HOURS = 3.0
TOP_N = 100
MIN_DELTA_KM = 1000.0
U_MIN, U_REF, U_MAX = 3.0, 3.7, 4.5          # km/s group-velocity band
G_SHEAR = 30e9                                # Pa   (vdE&B 2010)
C_PHASE = 3500.0                              # m/s
T_SW = 20.0                                   # s
DLAT = DLON = 0.5
EXPLORE_FRAC = 0.70
SCORING_HORIZON_MONTHS = 24
EMISSION_DEADLINE_MINUTES = 60


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def gcdist_km(la1, lo1, la2, lo2):
    p = np.pi / 180.0
    a = (np.sin((la2 - la1) * p / 2) ** 2
         + np.cos(la1 * p) * np.cos(la2 * p) * np.sin((lo2 - lo1) * p / 2) ** 2)
    return 2 * 6371.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def sigma_primary_Pa(M, r_km):
    """vdE&B eq.(6) surface-wave amplitude -> nominal peak dynamic stress, Pa."""
    Ddeg = np.maximum(r_km, 1.0) / 111.195
    A20_um = 10.0 ** (M - 1.66 * np.log10(Ddeg) - 2.0)
    V = 2 * np.pi * (A20_um * 1e-6) / T_SW
    return G_SHEAR * V / C_PHASE


def cell_domain():
    ctx = design.build_design(data_dir=DATA_DIR, dlat=DLAT, dlon=DLON,
                             explore_frac=EXPLORE_FRAC, verbose=False)
    return (np.asarray(ctx.grid.cell_lat, dtype=float),
            np.asarray(ctx.grid.cell_lon, dtype=float),
            np.asarray(ctx.grid.cell_key, dtype=np.int64))


def envelope_exposure(sources, t_trigger, cell_lat, cell_lon):
    """FROZEN. Returns (E, n_bins) with E[cell, bin] = nominal peak dynamic stress, Pa.

    `sources` = list of dicts with keys time (pd.Timestamp, UTC), lat, lon, mag.
    Bins are BIN_HOURS wide, bin 0 starting at `t_trigger`, covering HORIZON_DAYS.
    Cells within MIN_DELTA_KM of ANY source are set to zero (struck out).
    """
    n_bins = int(round(HORIZON_DAYS * 24.0 / BIN_HOURS))
    E = np.zeros((cell_lat.size, n_bins), dtype=np.float64)
    struck = np.zeros(cell_lat.size, dtype=bool)
    bin_lo = np.arange(n_bins) * BIN_HOURS * 3600.0            # seconds from trigger
    bin_hi = bin_lo + BIN_HOURS * 3600.0

    for s in sources:
        d = gcdist_km(float(s["lat"]), float(s["lon"]), cell_lat, cell_lon)
        struck |= (d < MIN_DELTA_KM)
        amp = sigma_primary_Pa(float(s["mag"]), d)             # Pa, per cell
        off = (pd.Timestamp(s["time"]) - t_trigger).total_seconds()
        t_lo = off + d / U_MAX                     # s (d in km, U in km/s)
        t_hi = off + d / U_MIN
        width = np.maximum(t_hi - t_lo, 1e-9)
        ov = (np.minimum(t_hi[:, None], bin_hi[None, :])
              - np.maximum(t_lo[:, None], bin_lo[None, :]))
        ov = np.clip(ov, 0.0, None) / width[:, None]
        E += amp[:, None] * ov

    E[struck, :] = 0.0
    return E, n_bins, struck


def rank_top(E, cell_lat, cell_lon, cell_key, t_trigger, n_bins):
    flat = E.ravel()
    n_pos = int((flat > 0).sum())
    k = min(TOP_N, n_pos)
    if k == 0:
        return [], 0
    # deterministic: sort by (-exposure, cell_index, bin_index)
    idx = np.argsort(flat, kind="stable")[::-1][: max(k * 4, k)]
    recs = []
    for f in idx:
        c, b = divmod(int(f), n_bins)
        recs.append((-float(flat[f]), c, b))
    recs.sort()
    recs = recs[:k]
    out = []
    for rank, (negE, c, b) in enumerate(recs, start=1):
        w0 = t_trigger + pd.Timedelta(hours=BIN_HOURS * b)
        out.append({
            "rank": rank,
            "cell_index": int(c),
            "cell_key": int(cell_key[c]),
            "lat_center": round(float(cell_lat[c]), 4),
            "lon_center": round(float(cell_lon[c]), 4),
            "bin_index": int(b),
            "window_start_utc": w0.isoformat(),
            "window_end_utc": (w0 + pd.Timedelta(hours=BIN_HOURS)).isoformat(),
            "exposure_kPa": round(-negE / 1e3, 6),
        })
    return out, n_pos


def protocol_block():
    src_files = sorted(glob.glob(os.path.join(HERE, DATA_DIR, "*.csv")))
    return {
        "entry": "K-069 (committed as K-069-min: envelope arm only)",
        "artifact": "prospective commitment log (commits; scores nothing)",
        "committed_by": "worker session, 2026-08-11",
        "commitment_opened_utc": "2026-08-11T00:00:00+00:00",
        "append_only": ("Window sets are APPENDED to `emissions`. No emission is ever "
                        "edited or removed. Each emission carries the sha256 of this "
                        "emitter as it stood when the emission was made."),
        "trigger_rule": {
            "trigger_magnitude": TRIGGER_MAG,
            "catalog": "USGS ComCat, global, first available solution",
            "emission_deadline_minutes": EMISSION_DEADLINE_MINUTES,
            "earliest_admissible_window_opens_minutes":
                round(MIN_DELTA_KM / U_MAX / 60.0, 1),
        },
        "emission_rule": {
            "top_n": TOP_N,
            "cell_deg": DLAT,
            "bin_hours": BIN_HOURS,
            "horizon_days": HORIZON_DAYS,
            "min_delta_km": MIN_DELTA_KM,
            "contributing_source_magnitude": CONTRIB_MAG,
            "contributing_source_lookback_days": HORIZON_DAYS,
            "group_velocity_km_s": [U_MIN, U_REF, U_MAX],
            "amplitude_model": ("van der Elst & Brodsky (2010) eq.(6): "
                                "log10 A20[um] = M - 1.66 log10(D[deg]) - 2; "
                                "V = 2*pi*A20/T_SW; sigma = G*V/c. "
                                "G=30 GPa, c=3500 m/s, T_SW=20 s. Nominal amplitude, "
                                "frozen ranking functional."),
            "overlay": "exposures of all contributing sources ADD (envelope leg)",
            "tie_break": "(-exposure, cell_index, bin_index) ascending; fully deterministic",
            "time_concentration_note": ("Declared in advance, from the 2026-07-17 demo "
                                        "emission: with a single contributing source the "
                                        "top-100 set collapses into the earliest 3 h bins "
                                        "-- an amplitude-only ephemeris has one wave "
                                        "passage and no phase, so the committed set is the "
                                        "wave-passage set. Time diversity appears only "
                                        "under overlay (two or more contributing sources). "
                                        "This makes K-069-min a PROSPECTIVE version of the "
                                        "classical instantaneous remote-triggering test, "
                                        "which is precisely the niche Merton records as "
                                        "empty, and it is NOT a test of interference."),
        },
        "domain": {
            "data_dir": DATA_DIR,
            "cell_domain": ("engine 0.5 deg active-cell domain over data/comcat_world "
                            "(cells with >=1 event in the exploration window); "
                            "identical to the K-080 census domain"),
            "boxes": [os.path.basename(p) for p in src_files],
            "coverage_note": "13 ComCat regional boxes, NOT the globe.",
        },
        "scoring": {
            "horizon_months": SCORING_HORIZON_MONTHS,
            "score_at_utc": "2028-08-11T00:00:00+00:00",
            "target": "M >= 4.5, epicentre in a committed cell, origin time in the committed 3 h bin",
            "expectation": ("frozen ETAS (K-080 census parameters, "
                            "engine/out/cache/etas_params_k080_0p5deg.json), integrated "
                            "over the identical cell-bins with history to each bin's start"),
            "statistic": "rate ratio observed/expected, and bits/event",
            "unit": "committed window",
            "N_eff": "number of distinct trigger events, NOT the number of windows",
            "pass": "rate-ratio CI on N_eff excluding 1 from above",
            "fail": "CI contains 1 -> clean prospective null, published as such",
            "power_note": ("K-069's POWER-STATE: if the expected total ETAS count in the "
                           "committed set is below ~300 at the horizon, the HORIZON "
                           "EXTENDS rather than the claim weakening. Recorded before any "
                           "emission exists."),
        },
        "not_committed_needs_infrastructure": [
            "CONSTRUCTIVE-maxima arm: requires K-059's phase-resolved surface-wave "
            "ephemeris T(x,t). NOT ON DISK. engine/ephemeris.py is an astronomical "
            "(solar/lunar) ephemeris for the tidal miner and is unrelated.",
            "Paired DESTRUCTIVE control set (K-061): phase-defined, unbuilt for the "
            "same reason.",
            "K-069's stated success rule (constructive CI > 1 AND destructive CI <= 1) "
            "therefore CANNOT be scored on this log. The envelope arm is a one-armed "
            "prospective test and is labelled as one.",
            "Interpretation gate K-034 is certified only at >= ~34 kPa; most committed "
            "windows will sit far below it. A null here is provisional-pending-K-034.",
            "When K-059 exists, the constructive/destructive arms open a SECOND, "
            "separately hashed log. They do not retro-fit into this one.",
        ],
        "provenance": {
            "emitter": "k069_emit.py",
            "emitter_sha256": sha256_file(os.path.join(HERE, "k069_emit.py")),
            "engine_version": ENGINE_VERSION,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "seed_declared": SEED,
            "input_sha256": {os.path.basename(p): sha256_file(p) for p in src_files},
        },
    }


def do_emit(args, out_path, demo=False):
    lat, lon, key = cell_domain()
    t_trig = pd.Timestamp(args.time)
    if t_trig.tz is None:
        t_trig = t_trig.tz_localize("UTC")

    sources = [{"time": t_trig, "lat": args.lat, "lon": args.lon, "mag": args.mag,
                "id": args.id or "TRIGGER"}]
    # contributing sources: catalogued M >= CONTRIB_MAG in the preceding HORIZON_DAYS
    frames = []
    for p in sorted(glob.glob(os.path.join(HERE, DATA_DIR, "*.csv"))):
        d = pd.read_csv(p, usecols=["time", "latitude", "longitude", "mag", "id"])
        frames.append(d)
    cat = pd.concat(frames, ignore_index=True).drop_duplicates(subset="id")
    cat["time"] = pd.to_datetime(cat["time"], utc=True, format="mixed")
    sel = ((cat["mag"] >= CONTRIB_MAG)
           & (cat["time"] < t_trig)
           & (cat["time"] >= t_trig - pd.Timedelta(days=HORIZON_DAYS))
           & (cat["id"] != (args.id or "")))
    for _, r in cat[sel].iterrows():
        sources.append({"time": r["time"], "lat": float(r["latitude"]),
                        "lon": float(r["longitude"]), "mag": float(r["mag"]),
                        "id": str(r["id"])})

    E, n_bins, struck = envelope_exposure(sources, t_trig, lat, lon)
    windows, n_pos = rank_top(E, lat, lon, key, t_trig, n_bins)

    emission = {
        "emission_id": f"{'DEMO-' if demo else ''}{args.id or t_trig.strftime('%Y%m%dT%H%M%S')}",
        "scored": (not demo),
        "demo_notice": (None if not demo else
                        "RETROSPECTIVE DEMONSTRATION OF THE EMITTER. NOT SCORED, NOT "
                        "COMMITTED, NOT ADMISSIBLE AS EVIDENCE IN EITHER DIRECTION. It "
                        "exists to prove the frozen emitter runs and produces windows "
                        "within the deadline."),
        "emitted_at_utc": pd.Timestamp.now(tz="UTC").isoformat(),
        "trigger": {"time_utc": t_trig.isoformat(), "lat": args.lat, "lon": args.lon,
                    "mag": args.mag, "id": args.id},
        "contributing_sources": [
            {"id": s["id"], "time_utc": pd.Timestamp(s["time"]).isoformat(),
             "lat": s["lat"], "lon": s["lon"], "mag": s["mag"]} for s in sources],
        "n_cells_domain": int(lat.size),
        "n_cells_struck_out": int(struck.sum()),
        "n_cell_bins_positive": int(n_pos),
        "n_windows_committed": len(windows),
        "exposure_kPa_max": (windows[0]["exposure_kPa"] if windows else None),
        "exposure_kPa_min_committed": (windows[-1]["exposure_kPa"] if windows else None),
        "emitter_sha256": sha256_file(os.path.join(HERE, "k069_emit.py")),
        "windows": windows,
    }

    if demo:
        with open(out_path, "w", encoding="utf-8") as fh:
            json.dump({"protocol": protocol_block(), "emission": emission}, fh, indent=1)
            fh.write("\n")
    else:
        doc = json.load(open(out_path, encoding="utf-8")) if os.path.exists(out_path) \
            else {"protocol": protocol_block(), "emissions": []}
        doc["emissions"].append(emission)
        with open(out_path, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, indent=1)
            fh.write("\n")
    return emission


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--init", action="store_true")
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--emit", action="store_true")
    ap.add_argument("--time")
    ap.add_argument("--lat", type=float)
    ap.add_argument("--lon", type=float)
    ap.add_argument("--mag", type=float)
    ap.add_argument("--id")
    args = ap.parse_args()

    if args.init:
        doc = {"protocol": protocol_block(), "emissions": []}
        with open(LOG, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, indent=1)
            fh.write("\n")
        print(f"wrote {LOG}")
        print(f"SHA256 {sha256_file(LOG)}")
        print(f"emitter SHA256 {sha256_file(os.path.join(HERE, 'k069_emit.py'))}")
        return 0

    if args.demo:
        # most recent M >= TRIGGER_MAG in the catalogue on disk
        frames = []
        for p in sorted(glob.glob(os.path.join(HERE, DATA_DIR, "*.csv"))):
            frames.append(pd.read_csv(p, usecols=["time", "latitude", "longitude",
                                                  "mag", "id"]))
        cat = pd.concat(frames, ignore_index=True).drop_duplicates(subset="id")
        cat["time"] = pd.to_datetime(cat["time"], utc=True, format="mixed")
        big = cat[cat["mag"] >= TRIGGER_MAG].sort_values("time")
        r = big.iloc[-1]
        args.time, args.lat = r["time"].isoformat(), float(r["latitude"])
        args.lon, args.mag, args.id = float(r["longitude"]), float(r["mag"]), str(r["id"])
        print(f"DEMO trigger: {args.id} M{args.mag} {args.time} "
              f"({args.lat:.3f},{args.lon:.3f})")
        t0 = time.time()
        em = do_emit(args, DEMO, demo=True)
        print(f"  windows={em['n_windows_committed']} "
              f"max={em['exposure_kPa_max']} kPa min={em['exposure_kPa_min_committed']} kPa")
        print(f"  emit wall = {time.time() - t0:.1f}s "
              f"(deadline {EMISSION_DEADLINE_MINUTES} min)")
        print(f"wrote {DEMO}")
        print(f"SHA256 {sha256_file(DEMO)}")
        return 0

    if args.emit:
        for k in ("time", "lat", "lon", "mag"):
            if getattr(args, k) is None:
                ap.error(f"--emit requires --{k}")
        if args.mag < TRIGGER_MAG:
            ap.error(f"trigger magnitude {args.mag} < frozen threshold {TRIGGER_MAG}")
        em = do_emit(args, LOG, demo=False)
        print(f"appended emission {em['emission_id']} "
              f"({em['n_windows_committed']} windows) to {LOG}")
        print(f"SHA256 {sha256_file(LOG)}")
        return 0

    ap.error("one of --init / --demo / --emit is required")


if __name__ == "__main__":
    sys.exit(main())
