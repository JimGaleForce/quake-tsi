"""Region partition + the 2R-df phase-incoherent regional sum statistic.

RULE TRACEABILITY (HYPOTHESIS_LEDGER.md)
----------------------------------------
* §P6-4 Rule 4.1  -- the region set is derived from EXPLORATION-WINDOW data only.
  The K-080 census cell list may be reported as a cross-check overlap fraction and
  may NEVER be the selector. `build_regions` takes `ctx` + the mining `window` and
  touches nothing else; there is no code path here that can read the census.
* §P6-4 Rule 4.2  -- the primary blind-spot kill is a PHASE-INCOHERENT statistic,
  not a bigger battery: per-region 2-df score statistics SUMMED to 2R df, one test
  per (feature, lag). The sum is a quadratic form that adds regardless of relative
  phase, so the domain-sum cancellation named in §K87-0(d)(i) cannot hide a signal
  from it.
* §P6-4 Rule 4.3 as AMENDED by §P7-1(d) -- per-region amplitudes are UNRESOLVED by
  default and are not quoted; only the summed statistic is quotable. A region may
  quote an amplitude only if it individually clears its own declared floor under
  the §P7-1(b) formula at its own N and the tranche's own alpha.
* §P6-4 Rule 4.5  -- the partition RULE is frozen in the config hash; the realised
  partition is a deterministic function of that rule and the frozen catalogue, and
  its digest is recorded in the checkpoint and asserted on resume.
* §P7-1(b)        -- the S-15 floor is the FORMULA
      A_min = sqrt(VIF) * (z_alpha + z_0.80) * sqrt(2/N)
  with alpha declared per tranche and VIF measured (F4-58). `a_min` implements it.
* S-9             -- the partition rule is written down ONE way. There is exactly
  one rule id in this module and no alternative is run.

THE DECLARED PARTITION RULE, ONE WAY (S-9)
------------------------------------------
`R2b-lon6-active`:

  1. A cell is ELIGIBLE if it holds at least one catalogue event at the session's
     magnitude target inside the mining window (the exploration window after ETAS
     burn-in). Exploration-window data only.
  2. Region index of an eligible cell = floor(((lon_center + 180) mod 360) / 60),
     i.e. six 60-degree longitude sectors. Longitude is the declared
     phase-coherence axis because solar and lunar tidal phase at fixed UT is a
     function of local time, i.e. of longitude -- which is exactly the
     region-dependent phase that the domain sum cancels.
  3. A sector is RETAINED if it holds at least `min_event_fraction` of the
     in-window events. Cells in a dropped sector are excluded from the statistic
     and both the dropped sector count and its excluded event count are reported.
  4. R is the number of retained sectors. The statistic has 2R df for a 2-df
     (phase) feature and R df for a 1-df (linear) feature.

WHY R = 6 AND NOT §P6-4's ANTICIPATED 12-24: see `results_f4_58_vif.json`. F4-58
measured VIF ~= 24.1 for 2-df phase features, not the ~3.94 §P7-1(a) inferred.
`max_R_for_sum` recomputes the choice from whatever VIF is declared, so the number
is a consequence of a measurement rather than a convention.
"""

from __future__ import annotations

import hashlib
import json
import math

import numpy as np
from scipy import stats

from . import floors

# ------------------------------------------------------------ declared rule --
REGION_RULE_ID = "R2b-lon6-active"
N_SECTORS = 6
MIN_EVENT_FRACTION = 0.01
SECTOR_WIDTH_DEG = 360.0 / N_SECTORS

# The floor's power convention is fixed at 80% by §P7-1(b).
Z_POWER_80 = float(stats.norm.ppf(0.80))

UNRESOLVED = "UNRESOLVED"
UNRESOLVED_RULE = (
    "HYPOTHESIS_LEDGER.md §P6-4 Rule 4.3 as AMENDED by §P7-1(d): per-region "
    "amplitudes are UNRESOLVED by default and are not quoted; the quotable object "
    "is the summed 2R-df statistic. A per-region amplitude may be quoted only for a "
    "region that individually clears its own declared floor under the §P7-1(b) "
    "formula at that region's own N and the tranche's own alpha."
)


# ------------------------------------------------------------- the S-15 floor -
def a_min(vif, alpha, n):
    """§P7-1(b): A_min = sqrt(VIF) * (z_alpha + z_0.80) * sqrt(2/N).

    `alpha` is the DECLARED operating threshold of the tranche (two-sided, so
    z_alpha = Phi^-1(1 - alpha/2) -- the convention §P7-1(a) uses in its own
    arithmetic, where alpha = 0.1/259 gives z = 3.549).

    Delegates to `engine.floors.a_min` so the formula has exactly ONE
    implementation in the build (§P7-8(d) added `engine/floors.py` as its home;
    this stays as the region module's own name for it).
    """
    return floors.a_min(vif, alpha, n)


def max_R_for_sum(vif, alpha, n_total, target_amplitude, r_cap=64):
    """Largest R for which the summed 2R-df statistic still has 80% power.

    A common-amplitude, arbitrary-phase regional signal of amplitude A gives the
    sum a non-centrality of lambda = N*A^2/(2*VIF) REGARDLESS of the phases -- that
    is the whole point of the phase-incoherent sum. Power is 80% when the 20th
    percentile of ncx2(2R, lambda) reaches the alpha critical value of chi2(2R).
    Adding regions costs df and buys phase resolution, so this returns the largest
    R that still clears the floor.
    """
    lam = float(n_total) * float(target_amplitude) ** 2 / (2.0 * float(vif))
    best = 0
    for R in range(1, int(r_cap) + 1):
        crit = float(stats.chi2.isf(float(alpha), 2 * R))
        if float(stats.ncx2.ppf(0.20, 2 * R, lam)) >= crit:
            best = R
        else:
            break
    return best, lam


def max_R_for_per_region_battery(vif, alpha, n_total, target_amplitude):
    """Largest R whose PER-REGION floor still resolves `target_amplitude`.

    A_min(region) = sqrt(VIF)*(z_alpha+z_0.80)*sqrt(2R/N)  <=  A_target
      =>  R <= (N/2) * (A_target / (sqrt(VIF)*(z_alpha+z_0.80)))^2
    Under the F4-58 measurement this is < 2 at any credible target amplitude, which
    is the §P7-1(d) outcome ("essentially no region will clear it").
    """
    z = float(stats.norm.isf(float(alpha) / 2.0))
    k = math.sqrt(float(vif)) * (z + Z_POWER_80)
    r_real = (float(n_total) / 2.0) * (float(target_amplitude) / k) ** 2
    return int(math.floor(r_real)), float(r_real)


# ------------------------------------------------------------ the partition --
def sector_of_lon(lon, n_sectors=N_SECTORS):
    """Region index from cell-centre longitude. Pure, vectorised, deterministic."""
    lon = np.asarray(lon, dtype=np.float64)
    w = 360.0 / int(n_sectors)
    idx = np.floor(((lon + 180.0) % 360.0) / w).astype(np.int64)
    return np.clip(idx, 0, int(n_sectors) - 1)


def build_regions(ctx, y, window, n_sectors=N_SECTORS,
                  min_event_fraction=MIN_EVENT_FRACTION):
    """The declared partition, from EXPLORATION-WINDOW data only (§P6-4 Rule 4.1).

    `y` is the (n_cells, n_days) day-count array at the session's magnitude target
    and `window` the mining window slice. Nothing outside `y[:, window]` is read,
    so the result cannot depend on holdout-window seismicity -- which is Finding A's
    entire complaint against the K-080 census list.

    Returns a dict; `region_of_cell` is -1 for every cell excluded by rule 1 or 3.
    """
    yw = np.asarray(y)[:, window]
    per_cell = np.asarray(yw.sum(axis=1), dtype=np.float64)
    lon = np.asarray(ctx.grid.cell_lon, dtype=np.float64)
    if per_cell.size != lon.size:
        raise ValueError(
            f"day-count array has {per_cell.size} cells but the grid has {lon.size}")

    eligible = per_cell > 0.0                       # rule 1
    sector = sector_of_lon(lon, n_sectors)          # rule 2
    n_total = float(per_cell[eligible].sum())
    if n_total <= 0:
        raise ValueError("no events in the mining window; cannot build a partition")

    sectors = []
    for s in range(int(n_sectors)):
        m = eligible & (sector == s)
        n_ev = float(per_cell[m].sum())
        sectors.append({
            "sector": int(s),
            "lon_lo": -180.0 + s * (360.0 / n_sectors),
            "lon_hi": -180.0 + (s + 1) * (360.0 / n_sectors),
            "n_cells": int(m.sum()),
            "n_events": n_ev,
            "event_fraction": n_ev / n_total,
        })

    retained = [s for s in sectors                  # rule 3
                if s["n_cells"] > 0 and s["event_fraction"] >= float(min_event_fraction)]
    dropped = [s for s in sectors if s not in retained]
    if not retained:
        raise ValueError("every longitude sector was dropped by the activity threshold")

    region_of_cell = np.full(lon.size, -1, dtype=np.int64)
    regions = []
    for r, s in enumerate(retained):                # rule 4
        m = eligible & (sector == s["sector"])
        region_of_cell[m] = r
        regions.append(dict(s, region=int(r)))

    out = {
        "rule_id": REGION_RULE_ID,
        "n_sectors": int(n_sectors),
        "min_event_fraction": float(min_event_fraction),
        "window": [int(window.start), int(window.stop)],
        "selector": "exploration-window active cells only (§P6-4 Rule 4.1)",
        "R": len(regions),
        "regions": regions,
        "dropped_sectors": dropped,
        "n_cells_assigned": int((region_of_cell >= 0).sum()),
        "n_cells_total": int(lon.size),
        "n_events_assigned": float(sum(r["n_events"] for r in regions)),
        "n_events_excluded": float(sum(s["n_events"] for s in dropped)),
        "region_of_cell": region_of_cell,
    }
    out["digest"] = partition_digest(out)
    return out


def partition_digest(part):
    """Stable digest of the REALISED partition, for the checkpoint + report.

    The RULE is what lives in the config hash (§P6-4 Rule 4.5); this pins the
    realisation so that a resumed session cannot silently re-partition.
    """
    blob = json.dumps(
        {
            "rule_id": part["rule_id"],
            "n_sectors": part["n_sectors"],
            "min_event_fraction": part["min_event_fraction"],
            "window": part["window"],
            "regions": [[r["region"], r["lon_lo"], r["lon_hi"], r["n_cells"],
                         round(float(r["n_events"]), 6)] for r in part["regions"]],
            "region_of_cell": np.asarray(part["region_of_cell"]).tolist(),
        },
        sort_keys=True,
    )
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


def regional_series(y, rate, region_of_cell, window, R):
    """(R, n_days) observed counts and ETAS expectations, one row per region.

    The two are aggregated with the SAME cell->region map, so each row is a
    self-consistent (counts, offset) pair and the per-region score test is the
    identical statistic the global test runs, restricted to that region's cells.
    """
    reg = np.asarray(region_of_cell)
    yw = np.asarray(y)[:, window]
    ow = np.asarray(rate)
    if ow.shape != yw.shape:
        raise ValueError(f"offset shape {ow.shape} != counts shape {yw.shape}")
    n_days = yw.shape[1]
    C = np.zeros((int(R), n_days), dtype=np.float64)
    O = np.zeros((int(R), n_days), dtype=np.float64)
    for r in range(int(R)):
        m = reg == r
        C[r] = yw[m].sum(axis=0)
        O[r] = ow[m].sum(axis=0)
    if not (np.isfinite(O).all() and (O > 0).all()):
        raise ValueError("a region has a non-positive ETAS expectation on some day")
    return C, O


def region_floor_table(part, vif, alpha, target_amplitude):
    """S-15 per region: A_min from the §P7-1(b) FORMULA, and measurable/unmeasurable.

    This is the table the banner prints for the optional per-region battery
    (§P6-4 Rule 4.7 item 2, per stratum).
    """
    rows = []
    for r in part["regions"]:
        n = float(r["n_events"])
        am = a_min(vif, alpha, n)
        rows.append({
            "region": int(r["region"]),
            "lon_lo": r["lon_lo"], "lon_hi": r["lon_hi"],
            "n_events": n,
            "a_min": am,
            "a_min_pct": 100.0 * am,
            "target_amplitude": float(target_amplitude),
            "s15": "MEASURABLE" if am <= float(target_amplitude) else "UNMEASURABLE",
        })
    n_meas = sum(1 for x in rows if x["s15"] == "MEASURABLE")
    return {
        "formula": "A_min = sqrt(VIF)*(z_alpha + z_0.80)*sqrt(2/N)",
        "vif": float(vif), "alpha": float(alpha),
        "z_alpha": float(stats.norm.isf(float(alpha) / 2.0)),
        "z_0.80": Z_POWER_80,
        "target_amplitude": float(target_amplitude),
        "rows": rows,
        "n_measurable": n_meas,
        "n_unmeasurable": len(rows) - n_meas,
        "fraction_unmeasurable": (len(rows) - n_meas) / float(len(rows)) if rows else 0.0,
    }


# --------------------------------------------- planted regional-phase signal --
def plant_regional_phase(C, O, values, amplitude, phases, rng=None):
    """Multiply each region's expectation by 1 + A*cos(2*pi*values + phase_r).

    The AMPLITUDE is the same in every region and the PHASE differs, which is the
    §K87-0(d)(i) blind spot exactly: the domain sum of R equally-spaced phases
    cancels to (numerically) zero while every region individually carries a
    full-amplitude modulation. §P6-4 Rule 4.7 item 5 demands that the 2R-df sum
    recover this and the global sum miss it.

    Returns new (C_planted, O) -- the offset is deliberately NOT modified, so the
    planted structure appears as a rate modulation against the declared baseline,
    which is what the score test is testing for.
    """
    rng = rng or np.random.default_rng(0)
    C = np.array(C, dtype=np.float64, copy=True)
    O = np.asarray(O, dtype=np.float64)
    v = np.asarray(values, dtype=np.float64)
    R = C.shape[0]
    if len(phases) != R:
        raise ValueError(f"{len(phases)} phases for {R} regions")
    for r in range(R):
        mod = 1.0 + float(amplitude) * np.cos(2.0 * np.pi * v + float(phases[r]))
        C[r] = rng.poisson(np.maximum(O[r] * mod, 1e-12)).astype(np.float64)
    return C


def equally_spaced_phases(R):
    """R phases whose unit vectors sum to (numerically) zero -- maximal cancellation."""
    return [2.0 * math.pi * r / float(R) for r in range(int(R))]
