"""K-080 -- THE OASIS CENSUS (frozen commitment, no claim until scored).

Enumerates every currently-live "oasis" (quiescence region) in the catalogue on disk,
as of the latest catalogue date, under K-076's frozen detector definition, and writes
the census to results_k080_census.json for hash-commitment in HYPOTHESIS_LEDGER.md.

This file computes a DENOMINATOR. It scores nothing. Scoring is at +5 years on a
Molchan trajectory (K-080's committed horizon); +1y and +2y readouts are descriptive
only, per the entry.

FROZEN DEFINITION (K-076 "Construction", as quoted in the ledger; deviations declared
in the JSON under `deviations_from_k076`):
  grid              0.5 deg x 0.5 deg
  neighbourhood     1.5 deg (166.79 km, great-circle, cell centre to cell centre)
  trailing window   90 days (the K-080/Colombia window; K-076's {30,90,365} family)
  Mc                4.5 (the comcat_world floor, fixed, never adaptive)
  baseline          Lambda = int lambda_hat_ETAS dt from the engine's frozen
                    space-time ETAS (etas-v1: nu*mu(c) + K*sum 10^(alpha*(m-4.5)) *
                    (t-t_i+c)^(-p) * Gaussian(sigma)), parameters fitted by Poisson
                    MLE on the EXPLORATION window only and never refit.
                    NOT own-cell climatology, NOT declustered background, NOT regional
                    mean rate -- K-076 rules out all three by name.
  oasis             Lambda >= 4  AND  z <= -2
  UNMEASURABLE      Lambda < 4 (S-15 class; scored neither way, reported in headline)

SIGN CONVENTION, stated because the ledger contains both. K-076 writes
z = (Lambda - N)/sqrt(Lambda) (positive = quiet); K-080 writes the oasis condition as
"z <= -2". Those are inconsistent as literals. The intended object in both is
"unusually quiet", so this file freezes the standardized-residual convention
    z = (N_obs - Lambda) / sqrt(Lambda)      (negative = quiet, oasis at z <= -2)
which makes K-080's stated threshold operative as written. Both signs are emitted per
cell (`z` and `z_k076_sign`) so no later reading can turn on the convention.

Run:  python -u k080_census.py
Deterministic: no randomness anywhere in this file (SEED declared and unused).
"""

from __future__ import annotations

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
from engine import baseline as bl
from engine import design, splits

SEED = 20260811                     # declared; no stochastic step exists in this file

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = "data/comcat_world"
OUT = os.path.join(HERE, "results_k080_census.json")

# ---- frozen detector constants (S-9: one value each, declared before the run) ----
DLAT = DLON = 0.5
NBR_DEG = 1.5
NBR_KM = 166.79                     # 1.5 deg of great circle
WINDOW_DAYS = 90
MC = 4.5
LAMBDA_FLOOR = 4.0                  # K-076's power floor: f=0 at z=2 needs Lambda>=4
Z_THRESHOLD = -2.0
EXPLORE_FRAC = 0.70
ETAS_CACHE = os.path.join("engine", "out", "cache", "etas_params_k080_0p5deg.json")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    t_start = time.time()
    print("=" * 78)
    print("K-080 OASIS CENSUS -- frozen commitment, no claim until scored")
    print("=" * 78)

    ctx = design.build_design(data_dir=DATA_DIR, dlat=DLAT, dlon=DLON,
                             explore_frac=EXPLORE_FRAC, verbose=True)
    explore, _holdout = splits.temporal_split(ctx.n_days, EXPLORE_FRAC)

    counts = ctx.day_counts(MC)                      # (n_cells, n_days)

    print("-" * 78)
    print("fitting/loading frozen ETAS (exploration window only, never refit)")
    base = bl.EtasV1(alpha=0.5, cache_path=ETAS_CACHE, verbose=True,
                     mag_target=MC)
    base.fit(ctx, counts, explore)
    lam = base.lam                                   # (n_cells, n_days) float32
    print(f"  etas params          = {base.params}")
    print(f"  params source        = {base.fit_info.get('source')}")

    # ---- census instant: the last day of the record -------------------------
    t0 = pd.Timestamp(ctx.meta["t0"])
    if t0.tz is None:
        t0 = t0.tz_localize("UTC")
    last_day = ctx.n_days - 1                        # 0-based index of final day
    lo = last_day - WINDOW_DAYS + 1
    assert lo >= 0
    census_time = t0 + pd.Timedelta(days=float(ctx.n_days))   # end of final day, UTC
    win_start = t0 + pd.Timedelta(days=float(lo))
    print("-" * 78)
    print(f"census instant (UTC)   = {census_time.isoformat()}")
    print(f"trailing window        = {win_start.isoformat()} .. {census_time.isoformat()}"
          f"  ({WINDOW_DAYS} d, day idx [{lo},{last_day}])")

    # ---- neighbourhood aggregation -----------------------------------------
    w = ctx.neighbour_matrix(NBR_KM)                 # (n_cells, n_cells) 0/1
    n_obs = w @ counts[:, lo:last_day + 1].sum(axis=1).astype(np.float64)
    Lam = w @ lam[:, lo:last_day + 1].sum(axis=1).astype(np.float64)
    n_nbr = w.sum(axis=1).astype(np.int64)

    z = (n_obs - Lam) / np.sqrt(Lam)                 # negative = quiet (frozen)
    measurable = Lam >= LAMBDA_FLOOR
    oasis = measurable & (z <= Z_THRESHOLD)

    area = ctx.grid.cell_area_km2()
    cell_lat = ctx.grid.cell_lat
    cell_lon = ctx.grid.cell_lon
    cell_key = ctx.grid.cell_key

    n_cells = int(ctx.n_cells)
    n_meas = int(measurable.sum())
    n_unmeas = n_cells - n_meas
    n_oasis = int(oasis.sum())
    a_tot = float(area.sum())
    a_meas = float(area[measurable].sum())
    a_oasis = float(area[oasis].sum())

    # days since the most recent M>=Mc anywhere in each cell's neighbourhood
    ev_last_day = np.full(n_cells, -1, dtype=np.int64)
    np.maximum.at(ev_last_day, ctx.ev_cell.astype(np.int64), ctx.ev_day.astype(np.int64))
    nbr_last = np.full(n_cells, -1, dtype=np.int64)
    for i in range(n_cells):
        idx = np.nonzero(w[i])[0]
        if idx.size:
            nbr_last[i] = ev_last_day[idx].max()
    quiet_days = np.where(nbr_last >= 0, last_day - nbr_last, -1)

    order = np.argsort(z)                            # quietest first
    oasis_idx = [int(i) for i in order if bool(oasis[i])]
    cells = [{
        "rank": r + 1,
        "cell_index": int(i),
        "cell_key": int(cell_key[i]),
        "lat_center": round(float(cell_lat[i]), 4),
        "lon_center": round(float(cell_lon[i]), 4),
        "n_obs_90d_nbhd": float(n_obs[i]),
        "lambda_etas_90d_nbhd": round(float(Lam[i]), 4),
        "z": round(float(z[i]), 4),
        "z_k076_sign": round(float(-z[i]), 4),
        "n_neighbour_cells": int(n_nbr[i]),
        "days_since_last_M4.5_in_nbhd": int(quiet_days[i]),
        "cell_area_km2": round(float(area[i]), 2),
    } for r, i in enumerate(oasis_idx)]

    src_files = sorted(glob.glob(os.path.join(HERE, DATA_DIR, "*.csv")))
    out = {
        "entry": "K-080",
        "artifact": "oasis census (frozen commitment; scores nothing)",
        "committed_by": "worker session, 2026-08-11",
        "no_claim_notice": ("This census makes no claim. It is a denominator committed "
                            "before the numerator exists. Scoring is a Molchan "
                            "trajectory at +5 years against frozen ETAS and against a "
                            "random-alarm reference of identical space-time volume; "
                            "+1y and +2y readouts are DESCRIPTIVE ONLY per the entry."),
        "census_instant_utc": census_time.isoformat(),
        "catalog_latest_event_utc": ctx.meta["load_report"]["t_max"],
        "catalog_earliest_event_utc": ctx.meta["load_report"]["t_min"],
        "trailing_window_utc": [win_start.isoformat(), census_time.isoformat()],
        "definition": {
            "grid_deg": [DLAT, DLON],
            "neighbourhood_deg": NBR_DEG,
            "neighbourhood_km": NBR_KM,
            "trailing_window_days": WINDOW_DAYS,
            "Mc": MC,
            "z": "(N_obs - Lambda) / sqrt(Lambda), negative = quiet",
            "z_threshold_oasis": Z_THRESHOLD,
            "lambda_floor_measurable": LAMBDA_FLOOR,
            "baseline": "frozen space-time ETAS (engine etas-v1), exploration-window fit",
            "unmeasurable_class": "Lambda < 4 -> UNMEASURABLE (S-15), scored neither way",
        },
        "etas": {
            "model": "etas-v1",
            "params": {k: float(v) for k, v in base.params.items()},
            "params_source": base.fit_info.get("source"),
            "fit_window_day_index": base.fit_info.get("fit_window_days"),
            "omori_trunc_days": base.fit_info.get("omori_trunc_days"),
            "bits_per_event_vs_climatology_train":
                base.fit_info.get("bits_per_event_vs_climatology"),
            "cache_path": ETAS_CACHE.replace("\\", "/"),
        },
        "domain": {
            "data_dir": DATA_DIR,
            "coverage_note": ("13 ComCat regional boxes, NOT the whole Earth. The "
                              "northern Andes (Colombia) is a hole in this catalogue; "
                              "K-076's global fill-in download does not exist yet. The "
                              "census is therefore complete for the domain it declares "
                              "and silent outside it."),
            "boxes": [os.path.basename(p) for p in src_files],
            "active_cell_definition": ("engine v1 domain restriction: cells containing "
                                       ">=1 catalogue event during the EXPLORATION "
                                       "window (first 70% of the span). Cells that "
                                       "first activated later are outside the domain."),
            "n_cells_active": n_cells,
            "n_cells_measurable": n_meas,
            "n_cells_unmeasurable": n_unmeas,
            "area_active_km2": a_tot,
            "area_measurable_km2": a_meas,
        },
        "census": {
            "n_oases": n_oasis,
            "area_oases_km2": a_oasis,
            "area_fraction_of_measurable": (a_oasis / a_meas) if a_meas else None,
            "area_fraction_of_active_domain": (a_oasis / a_tot) if a_tot else None,
            "area_fraction_of_earth_surface": a_oasis / 5.10072e8,
            "cell_fraction_of_measurable": (n_oasis / n_meas) if n_meas else None,
            "z_min": float(z.min()),
            "z_median_measurable": float(np.median(z[measurable])) if n_meas else None,
        },
        "oases": cells,
        "provenance": {
            "script": "k080_census.py",
            "script_sha256": sha256_file(os.path.join(HERE, "k080_census.py")),
            "engine_version": ENGINE_VERSION,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "seed_declared": SEED,
            "input_sha256": {os.path.basename(p): sha256_file(p) for p in src_files},
            "catalog_load_report": ctx.meta["load_report"],
            "wall_seconds": round(time.time() - t_start, 1),
        },
        "scoring_protocol": {
            "horizon_years": 5,
            "score_at_utc": (census_time + pd.Timedelta(days=365.25 * 5)).isoformat(),
            "descriptive_readouts_utc": [
                (census_time + pd.Timedelta(days=365.25)).isoformat(),
                (census_time + pd.Timedelta(days=365.25 * 2)).isoformat(),
            ],
            "target": "M >= 6.5, epicentre inside a committed oasis cell",
            "statistic": "Molchan trajectory (alarm-space fraction vs miss rate) + area skill score; ROC secondary; calibration curve tertiary",
            "references": ["frozen ETAS (this file's params)",
                           "random alarms of matched space-time volume",
                           "ETAS-sim alarms through the identical detector"],
            "pass": "Molchan point significantly below the diagonal",
            "fail": "it is not, and we publish the diagonal",
            "holdout_note": ("S-10 / P4-4(b): the engine's 70/30 split has already been "
                             "spent on a quiescence covariate. This census's confirmatory "
                             "scoring runs on events occurring AFTER the census instant, "
                             "which is a fresh window by construction and does not re-spend "
                             "that hash."),
        },
        "deviations_from_k076": [
            "K-076 does not exist yet as a built instrument; this file implements the "
            "minimal detector its 'Construction' paragraph specifies, and nothing more.",
            "Q0 (the injection-recovery licence) has NOT been run. The census makes no "
            "claim and needs no licence to be committed; the SCORING at +5y remains "
            "gated on Q0 exactly as the family's standing mandates require.",
            "Domain is the 13 comcat_world boxes, not the globe (fill-in download absent).",
            "Neighbourhood implemented as a great-circle radius between 0.5 deg cell "
            "centres (166.79 km), not as a lat/lon box.",
            "Trailing window frozen at 90 d, the single member of K-076's "
            "{30,90,365 d} family that K-080's text uses.",
            "Q1 (length-bias): the confirmatory statistic here is the ETAS-expected "
            "count, as mandated; no time-sampled or event-sampled gap frame is used.",
            "Q2 (Mc/station-count rerun) is NOT run here and is owed before scoring.",
        ],
    }

    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, sort_keys=False)
        fh.write("\n")

    digest = sha256_file(OUT)
    print("-" * 78)
    print(f"n_cells active         = {n_cells}")
    print(f"n_cells measurable     = {n_meas}  (Lambda >= {LAMBDA_FLOOR})")
    print(f"n_cells UNMEASURABLE   = {n_unmeas}")
    print(f"LIVE OASES             = {n_oasis}")
    print(f"  area                 = {a_oasis:,.0f} km^2")
    print(f"  frac of measurable   = {a_oasis / a_meas if a_meas else float('nan'):.4f}")
    print(f"  frac of active domain= {a_oasis / a_tot:.4f}")
    print(f"  frac of Earth surface= {a_oasis / 5.10072e8:.5f}")
    print(f"wrote {OUT}")
    print(f"SHA256 {digest}")
    print(f"wall {time.time() - t_start:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
