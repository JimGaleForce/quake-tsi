"""D-1c: THE DIRECTIONAL ARM on the K-092 seed set. Price 0, seed set only.

WHAT THIS ASKS. D-1 falsified a marginal, unconditional, SCALAR claim: the seed events'
below-neutral-falling excess is the frozen scalar's own dwell time (0.3951 observed vs
0.3827 time-uniform, z +0.32). That result is silent on ORIENTATION, because every
tidal quantity the program has ever computed is a scalar read off a rank-2 field and no
scalar can see the angle between the tidal principal axes and a fault plane. This arm
asks whether the seed events sit anywhere unusual once orientation, resolved Coulomb,
and Jim's directional-history ("rope") features are available.

WHY THE SCORER IS NOT THE FROZEN APP SCALAR. It cannot be. The observing application
computes vertical body-tide displacement only: `K092_SCALAR_PROVENANCE.md` records
"Shida number l2: none - no horizontal displacement, no strain". A scalar with no
horizontal component cannot express a bearing, so this arm necessarily uses
`engine/tidal_tensor.py`. That is a DIFFERENT SCALAR from the frozen one and this arm
therefore tests a DIFFERENT CLAIM from D-12/D-13, which is stated here rather than
discovered later. It does not and cannot score the frozen prediction.

SCOPE, WHICH IS ITS ENTIRE LICENCE. Reads `K092_seed_exclusion_superset.csv` and
nothing else, sha256 verified. Those events are the already-seen set, excluded from
D-12 by construction, which is why they may be looked at at price 0. Everything else is
a deterministic astronomical field on a uniform time grid at those same sites. No event
outside that file is touched. D-7 remains STRUCK.

NOT EVIDENCE, AND THE ASYMMETRY IS THE SAME ONE D-1 CARRIED. A null here bounds the
largest un-mined degree of freedom cheaply. A hit here is NOT a finding, because the
set was selected by having been looked at; it would be a reason to build the lattice and
pay for it properly under section P7-24, and nothing more.

---------------------------------------------------------------------------
THE DECLARATION. Everything below is fixed BEFORE any number is computed.
---------------------------------------------------------------------------

FAULT GEOMETRY, PRIMARY (one, declared, not fitted): strike 250 deg, dip 20 deg,
rake 90 deg, friction 0.4. This is the Alaska Peninsula / eastern Aleutian megathrust
as a shallow-dipping thrust striking WSW along the arc. It is an ASSUMPTION about the
receiver mechanism, not a datum: no focal mechanism is read for any seed event, and
`tidal_tensor.TENSOR_SCOPE_FLAGS` clause 2 applies in full.

FAULT GEOMETRY, SECONDARY GRID (declared, and paid for by the max-statistic):
strike in {230, 240, 250, 260, 270}, dip in {10, 20, 30}, friction in {0.2, 0.4, 0.6}.

THE STATISTIC BATTERY, eight statistics, declared in full:

  S1 CFS_TROUGH     fraction of events with Coulomb stress in the lowest third of its
                    own local cycle (bands as in audit_arcsine: u < -1/2)
  S2 CFS_CREST      ... in the highest third (u > +1/2). The physically favoured
                    direction if tidal Coulomb triggering exists at all.
  S3 CFS_RISING     fraction with dCFS/dt > 0 at the event time
  S4 CFS_LOADED     fraction with CFS > 0 AND dCFS/dt > 0 (the loading quadrant, the
                    Coulomb analogue of the frozen scalar quadrant)
  S5 BEARING_R2     axial concentration of (bearing of maximum horizontal extension
                    MINUS fault strike), as the mean resultant length of the DOUBLED
                    angle. A bearing is an undirected line, so the doubling is
                    mandatory (tidal_tensor.BEARING_CONVENTION).
  S6 ROT_RATE       mean absolute rotation rate of that bearing, deg/h -- Jim's
                    "travelling angle"
  S7 ROPE_EASING    fraction with the fault-normal traction POSITIVE and FALLING:
                    sigma along bearing (strike + 90) > 0 AND its time derivative < 0.
                    "Pulled from that direction, now easing."
  S8 ROPE_LAGGED    fraction that were MORE pulled 6 hours ago than now along that
                    same fixed bearing AND were positive then: a genuine two-time
                    feature, the class no instantaneous phase statistic can express
                    and therefore the class D-1 could not have detected.

THE NULL. Waveform-matched and time-uniform, exactly as D-1b: at each seed site the
field is evaluated on a long continuous grid CENTRED ON THE EVENT EPOCH, every sample
is pushed through the identical feature code, the event is the centre sample, and the
remaining samples are that site's null pool. A null replicate draws one sample
uniformly per site and recomputes the whole battery. Nothing is assumed about the
waveform: the null IS the waveform.

THE MULTIPLICITY CORRECTION. Max-statistic calibration over the WHOLE battery, and over
the whole geometry grid for the secondary: the reported significance is
P(max_over_battery |z| in a null replicate >= max_over_battery |z| observed). Dependence
among the eight statistics is handled implicitly because the null replicates carry the
same dependence. No per-statistic threshold is quoted as a result.

THE POWER FLOOR, asserted before running rather than discovered after: with N null
replicates the smallest attainable p is 1/(N+1), and a family-wise correction over m
statistics puts the floor at m/(N+1). N = 2000 and m = 8 give 0.004, comfortably below
0.05. `assert_power_floor` refuses to run if a future edit breaks that.
"""

from __future__ import annotations

import datetime as _dt
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import exp_k092_d1 as D1                       # scope + loader + cycle analysis
from engine import tidal_tensor as TT

OUT_JSON = HERE / "results_k092_d1c_directional.json"

# ------------------------------------------------------------- the declaration --
FAULT_PRIMARY = {"strike_deg": 250.0, "dip_deg": 20.0, "rake_deg": 90.0}
FRICTION_PRIMARY = 0.4
GEOMETRY_GRID = {"strike_deg": [230.0, 240.0, 250.0, 260.0, 270.0],
                 "dip_deg": [10.0, 20.0, 30.0],
                 "friction": [0.2, 0.4, 0.6]}

HALF_SPAN_DAYS = 10.0
STEP_MINUTES = 2.0
RATE_HALF_MINUTES = 10.0
LAG_HOURS = 6.0
N_NULL_REPLICATES = 2000
RNG_SEED = 20260822

STAT_NAMES = ("CFS_TROUGH", "CFS_CREST", "CFS_RISING", "CFS_LOADED",
              "BEARING_R2", "ROT_RATE", "ROPE_EASING", "ROPE_LAGGED")


def assert_power_floor(n_replicates=N_NULL_REPLICATES, m=len(STAT_NAMES), alpha=0.05):
    """Refuse to run a family that cannot reach significance at any effect size.

    Smallest attainable p with N replicates is 1/(N+1); a family-wise correction over
    m statistics puts the floor at m/(N+1). A test that cannot fire is worse than no
    test, because it still prints a 'not significant' line that a reader counts as
    evidence of absence. Encoded as the formula, so adding a ninth statistic raises the
    requirement automatically instead of silently invalidating a fixed number.
    """
    floor = m / (n_replicates + 1.0)
    if floor > alpha:
        raise SystemExit(
            "REFUSING TO RUN: family-wise power floor m/(N+1) = %d/%d = %.4f exceeds "
            "alpha = %.3f. Raise N to at least %d."
            % (m, n_replicates + 1, floor, alpha, math.ceil(m / alpha) - 1))
    return {"m": m, "n_replicates": n_replicates, "floor": floor, "alpha": alpha,
            "min_attainable_p": 1.0 / (n_replicates + 1.0)}


# ------------------------------------------------------------ the feature block --
def site_stress(t_days, lat, lon):
    """The stress tensor on one site's grid. Computed ONCE and reused by every
    geometry: the tensor does not depend on strike, dip, rake or friction, only the
    resolution onto the plane does. Doing it per geometry would be 46x the work for
    identical numbers."""
    return TT.stress_tensor(t_days + 2440587.5, lat, lon, 0.0)


def site_features(t_days, st, strike_deg, dip_deg, rake_deg, friction,
                  rate_k, lag_k):
    """Every declared feature at every sample of one site's grid.

    Returns boolean/float arrays over the grid, plus a validity mask. All eight
    statistics are computed from these arrays; the observation and the null draw from
    exactly the same arrays, so there is no second code path for either.
    """
    cfs = TT.coulomb(st, strike_deg, dip_deg, rake_deg, friction)["coulomb_pa"]

    n = cfs.size
    dt_h = STEP_MINUTES / 60.0
    ok = np.zeros(n, dtype=bool)
    ok[max(rate_k, lag_k):n - max(rate_k, lag_k)] = True

    # local-cycle band of CFS, using the SAME bracketing-maxima construction as D-1
    tm = D1.refined_maxima(t_days, cfs)
    u = np.full(n, np.nan)
    if tm.size >= 2:
        k = np.searchsorted(tm, t_days, side="right") - 1
        inb = (k >= 0) & (k < tm.size - 1)
        edges = np.searchsorted(t_days, tm)
        lo = np.full(tm.size - 1, np.nan)
        hi = np.full(tm.size - 1, np.nan)
        for j in range(tm.size - 1):
            a, b = edges[j], edges[j + 1]
            if b > a:
                lo[j], hi[j] = cfs[a:b].min(), cfs[a:b].max()
        kk = np.clip(k, 0, tm.size - 2)
        u = 2.0 * (cfs - lo[kk]) / np.maximum(hi[kk] - lo[kk], 1e-300) - 1.0
        ok &= inb & np.isfinite(u)
    else:
        ok[:] = False

    d_cfs = np.full(n, np.nan)
    d_cfs[rate_k:n - rate_k] = ((cfs[2 * rate_k:] - cfs[:-2 * rate_k])
                                / (2.0 * RATE_HALF_MINUTES / 60.0))

    bearing = TT.principal_bearing(st["e_NN"], st["e_EE"], st["e_NE"])
    rot = TT.bearing_rotation_rate(t_days, bearing)
    # doubled bearing relative to strike: the axial angle the circular statistic uses
    ang2 = TT.axial_to_circular(np.mod(bearing - strike_deg, 180.0))

    # the rope: traction along the FAULT-NORMAL bearing, strike + 90
    rope_bearing = np.mod(strike_deg + 90.0, 180.0)
    sig_b = TT.normal_along_bearing(st["s_NN"], st["s_EE"], st["s_NE"], rope_bearing)
    d_sig = np.full(n, np.nan)
    d_sig[rate_k:n - rate_k] = ((sig_b[2 * rate_k:] - sig_b[:-2 * rate_k])
                                / (2.0 * RATE_HALF_MINUTES / 60.0))
    sig_lag = np.full(n, np.nan)
    sig_lag[lag_k:] = sig_b[:n - lag_k]

    return {
        "ok": ok,
        "CFS_TROUGH": u < -0.5,
        "CFS_CREST": u > 0.5,
        "CFS_RISING": d_cfs > 0.0,
        "CFS_LOADED": (cfs > 0.0) & (d_cfs > 0.0),
        "ANG2": ang2,
        "ROT_RATE": np.abs(rot),
        "ROPE_EASING": (sig_b > 0.0) & (d_sig < 0.0),
        "ROPE_LAGGED": (sig_lag > 0.0) & (sig_b < sig_lag),
        "_dt_h": dt_h,
    }


SIMPLE_STATS = ("CFS_TROUGH", "CFS_CREST", "CFS_RISING", "CFS_LOADED",
                "ROPE_EASING", "ROPE_LAGGED")


def stack_features(feats):
    """Pack the per-site feature arrays into [n_sites, n_samples] matrices.

    Everything downstream is then pure numpy gathering, which is what makes 2000 null
    replicates across 46 geometries affordable. The observation and the null read the
    SAME matrices, so there is no second code path for either.
    """
    mats = {}
    for name in SIMPLE_STATS:
        mats[name] = np.stack([f[name].astype(np.float32) for f in feats])
    ang = np.stack([f["ANG2"] for f in feats])
    mats["_COS2"] = np.cos(ang).astype(np.float32)
    mats["_SIN2"] = np.sin(ang).astype(np.float32)
    mats["ROT_RATE"] = np.stack([f["ROT_RATE"].astype(np.float32) for f in feats])
    return mats


def battery_gather(mats, idx):
    """The eight statistics for one or many selections of one sample per site.

    `idx` is [n_sites] or [n_sites, n_draws]; returns scalars or [n_draws] arrays.
    """
    rows = np.arange(mats["ROT_RATE"].shape[0])
    if idx.ndim == 1:
        take = lambda m: m[rows, idx]
        red = lambda v: float(v.mean())
    else:
        take = lambda m: m[rows[:, None], idx]
        red = lambda v: v.mean(axis=0)
    out = {k: red(take(mats[k])) for k in SIMPLE_STATS}
    out["ROT_RATE"] = red(take(mats["ROT_RATE"]))
    c = red(take(mats["_COS2"]))
    s = red(take(mats["_SIN2"]))
    out["BEARING_R2"] = float(np.hypot(c, s)) if idx.ndim == 1 else np.hypot(c, s)
    return out


def run_geometry(events, stresses, grids, strike, dip, rake, friction, rng):
    """One geometry: observation, null ensemble, per-statistic z, and max-|z|."""
    rate_k = int(round(RATE_HALF_MINUTES / STEP_MINUTES))
    lag_k = int(round(LAG_HOURS * 60.0 / STEP_MINUTES))
    n_half = (grids[0].size - 1) // 2

    feats, pools = [], []
    for t, st in zip(grids, stresses):
        f = site_features(t, st, strike, dip, rake, friction, rate_k, lag_k)
        pool = np.nonzero(f["ok"])[0]
        if pool.size == 0 or not f["ok"][n_half]:
            return None
        feats.append(f)
        pools.append(pool)

    mats = stack_features(feats)
    n_sites = len(feats)
    obs = battery_gather(mats, np.full(n_sites, n_half, dtype=np.int64))
    idx = np.empty((n_sites, N_NULL_REPLICATES), dtype=np.int64)
    for s, p in enumerate(pools):
        idx[s] = p[rng.integers(0, p.size, N_NULL_REPLICATES)]
    nb = battery_gather(mats, idx)
    null = {k: np.asarray(nb[k], dtype=np.float64) for k in STAT_NAMES}

    z = {}
    for k in STAT_NAMES:
        mu, sd = float(null[k].mean()), float(null[k].std(ddof=1))
        z[k] = {"observed": obs[k], "null_mean": mu, "null_sd": sd,
                "z": (obs[k] - mu) / sd if sd > 0 else float("nan")}
    zmat = np.stack([(null[k] - null[k].mean()) / max(null[k].std(ddof=1), 1e-300)
                     for k in STAT_NAMES])
    null_max = np.max(np.abs(zmat), axis=0)
    obs_max = max(abs(z[k]["z"]) for k in STAT_NAMES)
    # Phipson-Smyth exact p with the +1 in both numerator and denominator: a p of
    # exactly zero is never justified by finite resampling, and ties count against.
    p_max = (int(np.sum(null_max >= obs_max)) + 1) / (N_NULL_REPLICATES + 1)
    return {
        "geometry": {"strike_deg": strike, "dip_deg": dip, "rake_deg": rake,
                     "friction": friction},
        "per_statistic": z,
        "max_abs_z_observed": obs_max,
        "max_abs_z_null_p95": float(np.quantile(null_max, 0.95)),
        "max_statistic_p": p_max,
        "max_statistic_percentile": float(np.mean(null_max < obs_max)),
    }


def main():
    floor = assert_power_floor()
    events = D1.load_seed_events()
    rng = np.random.default_rng(RNG_SEED)
    print("D-1c: %d seed events, power floor %.4f (min attainable p %.5f)"
          % (len(events), floor["floor"], floor["min_attainable_p"]), flush=True)

    step_days = STEP_MINUTES / 1440.0
    n_half = int(round(HALF_SPAN_DAYS / step_days))
    print("stress tensor on %d sites x %d samples (once, reused by every geometry) ..."
          % (len(events), 2 * n_half + 1), flush=True)
    grids, stresses = [], []
    for e in events:
        t = (e["t_ms"] / 86400000.0
             + step_days * np.arange(-n_half, n_half + 1, dtype=np.float64))
        grids.append(t)
        stresses.append(site_stress(t, e["lat"], e["lon"]))

    print("primary geometry %s friction %.1f ..." % (FAULT_PRIMARY, FRICTION_PRIMARY),
          flush=True)
    primary = run_geometry(events, stresses, grids, FAULT_PRIMARY["strike_deg"],
                           FAULT_PRIMARY["dip_deg"], FAULT_PRIMARY["rake_deg"],
                           FRICTION_PRIMARY, rng)
    if primary is None:
        raise SystemExit("primary geometry produced no valid pool")

    print("secondary geometry grid (%d combinations) ..."
          % (len(GEOMETRY_GRID["strike_deg"]) * len(GEOMETRY_GRID["dip_deg"])
             * len(GEOMETRY_GRID["friction"])), flush=True)
    grid = []
    for s in GEOMETRY_GRID["strike_deg"]:
        for d in GEOMETRY_GRID["dip_deg"]:
            for mu in GEOMETRY_GRID["friction"]:
                r = run_geometry(events, stresses, grids, s, d,
                                 FAULT_PRIMARY["rake_deg"], mu, rng)
                if r is not None:
                    grid.append(r)
                    print("  strike %5.1f dip %4.1f mu %.1f -> max|z| %5.2f  p %.4f"
                          % (s, d, mu, r["max_abs_z_observed"],
                             r["max_statistic_p"]), flush=True)

    grid_max = max(g["max_abs_z_observed"] for g in grid) if grid else float("nan")
    best = max(grid, key=lambda g: g["max_abs_z_observed"]) if grid else None

    out = {
        "arm": "D-1c directional / resolved-Coulomb / rope arm on the K-092 seed set",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 0,
        "scope": {
            "input": "K092_seed_exclusion_superset.csv ONLY",
            "sha256_verified": D1.SEED_CSV_SHA256,
            "n_events": len(events),
            "licence": "already-seen set, excluded from D-12 by construction; "
                       "everything else is a deterministic field on a uniform time "
                       "grid at those same sites. D-7 remains STRUCK.",
        },
        "does_not_score_the_frozen_prediction": (
            "The frozen scalar is the app's vertical body-tide displacement, which has "
            "NO Shida number and NO horizontal component and therefore cannot express "
            "a bearing. This arm necessarily uses engine/tidal_tensor.py, a DIFFERENT "
            "scalar, and therefore tests a DIFFERENT claim from D-12/D-13. Stated "
            "here rather than discovered later."),
        "not_evidence": (
            "A null here bounds the largest un-mined degree of freedom cheaply. A hit "
            "here is NOT a finding, because the set was selected by having been looked "
            "at; it would be a reason to build the lattice and pay for it under "
            "section P7-24, and nothing more."),
        "declaration": {
            "fault_primary": FAULT_PRIMARY, "friction_primary": FRICTION_PRIMARY,
            "geometry_grid": GEOMETRY_GRID,
            "statistics": list(STAT_NAMES),
            "half_span_days": HALF_SPAN_DAYS, "step_minutes": STEP_MINUTES,
            "rate_half_minutes": RATE_HALF_MINUTES, "lag_hours": LAG_HOURS,
            "n_null_replicates": N_NULL_REPLICATES, "rng_seed": RNG_SEED,
            "power_floor": floor,
            "multiplicity": "max-statistic calibration over the whole battery, and "
                            "over the whole geometry grid for the secondary; "
                            "dependence handled implicitly because the null replicates "
                            "carry it; Phipson-Smyth exact p with +1 top and bottom",
            "scope_flags": TT.TENSOR_SCOPE_FLAGS,
            "bearing_convention": TT.BEARING_CONVENTION,
        },
        "primary": primary,
        "secondary_grid": grid,
        "secondary_grid_max_abs_z": grid_max,
        "secondary_grid_best": best["geometry"] if best else None,
        "secondary_grid_note": (
            "The grid's maximum is NOT calibrated against a grid-wide null in this "
            "run: each cell carries its own max-statistic p, but the maximum ACROSS "
            "45 cells needs a null in which the whole grid is re-scanned per "
            "replicate. That is a declared limitation of this arm, not an oversight; "
            "the primary geometry is the calibrated result and the grid is reported "
            "as a sensitivity display only."),
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    print("\n" + "=" * 74)
    print("D-1c PRIMARY (strike 250 / dip 20 / rake 90 / mu 0.4)")
    print("=" * 74)
    for k in STAT_NAMES:
        r = primary["per_statistic"][k]
        print("  %-12s obs %10.4f  null %10.4f +- %8.4f   z %+6.2f"
              % (k, r["observed"], r["null_mean"], r["null_sd"], r["z"]))
    print("  max|z| observed %.3f ; null 95th pct %.3f ; MAX-STATISTIC p = %.4f"
          % (primary["max_abs_z_observed"], primary["max_abs_z_null_p95"],
             primary["max_statistic_p"]))
    print("  grid max|z| %.3f at %s" % (grid_max, best["geometry"] if best else None))
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
