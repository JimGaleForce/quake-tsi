"""D-1e: THE FULL ANGULAR HARMONIC SPECTRUM on the K-092 seed set. Price 0.

WHY THIS ARM EXISTS. Jim's correction: if the mechanism is symmetric -- pull from the
left or from the right both separate the fault -- then a SIGNED directional test
reports null exactly when the effect is total. Ten events from the left and ten from
the right is twenty out of twenty, and it gives a first-moment resultant of zero.

Every angular statistic this program has run so far has been either a first moment or
ONE hand-chosen doubled one. Nobody asked which harmonic carries the concentration.
This arm asks, for every angle it has, at every order from 1 to 4, and reports all of
them.

WHAT THE EXISTING RUNS DID AND DID NOT COVER, stated so this arm is not oversold.
D-1c's BEARING_R2 and D-1d's BEARING_RESID_R2 were BOTH already order-2 statistics on
the bearing relative to fault strike -- i.e. they were correctly axial, and their nulls
do stand for the symmetric case. What was never tested is order 1 (a genuinely
directional effect), order 3 and order 4 (conjugate shear geometry sits at 4), the body
azimuths at any order, and the elevation PARITY split. This arm covers those.

---------------------------------------------------------------------------
THE DECLARATION. Fixed BEFORE any number is computed.
---------------------------------------------------------------------------

ANGLES (each tested at orders m = 1, 2, 3, 4):

  A1  tidal principal-extension bearing MINUS fault strike. Axial by nature: order 2
      is where fault-normal extension lives.
  A2  lunar azimuth MINUS fault strike. Directional by nature: order 1 is meaningful
      here in a way it is not for A1, because the Moon's direction is a true bearing.
  A3  solar azimuth MINUS fault strike.

PARITY (the same argument applied to a linear coordinate):

  P1  EVEN coordinate of lunar sin(elevation), |sin(elev)|, in the top third of its own
      local cycle. This is the coordinate the DEGREE-2 potential actually lives in,
      because P2(cos z) is even and cannot distinguish overhead from underfoot.
  P2  ODD coordinate: lunar elevation above the horizon rather than below. Only a
      degree-3-like asymmetry can carry this one.

Total declared statistics: 3 angles x 4 orders + 2 parity = 14. Power floor 14/2001 =
0.0070, asserted before the run.

THE NULL: waveform-matched and time-uniform, identical in construction to D-1b, D-1c
and D-1d. A bearing's own dwell distribution is strongly non-uniform (axial R2 from
0.14 at Alaska latitudes to 0.59 at the equator), so a uniform-circle null would
manufacture an effect. Every R_m is scored against the R_m the SAME construction
produces on uniformly-random times at the same sites and epochs.

MULTIPLICITY: max-statistic over all 14, calibrated against the same null.

NOT EVIDENCE. The seed set was selected by having been looked at. A null extends the
bound into the symmetric and higher-order directions; a hit is a reason to build the
lattice, not a finding.
"""

from __future__ import annotations

import datetime as _dt
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import exp_k092_d1 as D1
from engine import circular_symmetry as CS
from engine import dwell_null as DN
from engine import tidal_tensor as TT

OUT_JSON = HERE / "results_k092_d1e_harmonics.json"

STRIKE_DEG = 250.0
HALF_SPAN_DAYS = 10.0
STEP_MINUTES = 2.0
RATE_HALF_MINUTES = 10.0
ORDERS = (1, 2, 3, 4)
ANGLES = ("A1_bearing_minus_strike", "A2_moon_az_minus_strike",
          "A3_sun_az_minus_strike")
PARITY = ("P1_elev_even_high", "P2_elev_odd_above")
STAT_NAMES = tuple("%s_m%d" % (a, m) for a in ANGLES for m in ORDERS) + PARITY
N_NULL_REPLICATES = 2000
RNG_SEED = 20260822


def local_cycle_u(t, x):
    tm = D1.refined_maxima(t, x)
    n = x.size
    if tm.size < 2:
        return np.full(n, np.nan)
    k = np.searchsorted(tm, t, side="right") - 1
    ok = (k >= 0) & (k < tm.size - 1)
    edges = np.searchsorted(t, tm)
    lo = np.full(tm.size - 1, np.nan)
    hi = np.full(tm.size - 1, np.nan)
    for j in range(tm.size - 1):
        a, b = edges[j], edges[j + 1]
        if b > a:
            lo[j], hi[j] = x[a:b].min(), x[a:b].max()
    kk = np.clip(k, 0, tm.size - 2)
    u = 2.0 * (x - lo[kk]) / np.maximum(hi[kk] - lo[kk], 1e-300) - 1.0
    return np.where(ok, u, np.nan)


def site_features(t, lat, lon, rate_k):
    jd = t + 2440587.5
    st = TT.stress_tensor(jd, lat, lon, 0.0)
    bearing = TT.principal_bearing(st["e_NN"], st["e_EE"], st["e_NE"])
    moon = TT.body_direction(jd, lat, lon, 0.0, "moon")
    sun = TT.body_direction(jd, lat, lon, 0.0, "sun")

    d2r = np.pi / 180.0
    a1 = (bearing - STRIKE_DEG) * d2r
    a2 = (moon["azimuth_deg"] - STRIKE_DEG) * d2r
    a3 = (sun["azimuth_deg"] - STRIKE_DEG) * d2r

    par = CS.fold_parity(moon["sin_elevation"])
    u_even = local_cycle_u(t, par["even"])

    n = bearing.size
    ok = np.zeros(n, dtype=bool)
    margin = rate_k + int(round(1.5 * 1440.0 / STEP_MINUTES))
    ok[margin:n - margin] = True
    ok &= np.isfinite(u_even) & np.isfinite(a1) & np.isfinite(a2) & np.isfinite(a3)
    return {
        "ok": ok,
        "A1_bearing_minus_strike": np.mod(a1, 2.0 * np.pi),
        "A2_moon_az_minus_strike": np.mod(a2, 2.0 * np.pi),
        "A3_sun_az_minus_strike": np.mod(a3, 2.0 * np.pi),
        "P1_elev_even_high": u_even > (1.0 / 3.0),
        "P2_elev_odd_above": par["sign"] > 0.0,
    }


def battery(mats, idx):
    """R_m for every angle and order, plus the two parity fractions."""
    rows = np.arange(next(iter(mats.values())).shape[0])
    if idx.ndim == 1:
        take = lambda m: m[rows, idx]
        red_mean = lambda v: float(v.mean())
        red_abs = lambda c, s: float(np.hypot(c, s))
    else:
        take = lambda m: m[rows[:, None], idx]
        red_mean = lambda v: v.mean(axis=0)
        red_abs = lambda c, s: np.hypot(c, s)
    out = {}
    for a in ANGLES:
        th = take(mats[a]).astype(np.float64)
        for m in ORDERS:
            out["%s_m%d" % (a, m)] = red_abs(red_mean(np.cos(m * th)),
                                             red_mean(np.sin(m * th)))
    for p in PARITY:
        out[p] = red_mean(take(mats[p]).astype(np.float64))
    return out


def main():
    floor = DN.assert_power_floor(len(STAT_NAMES), N_NULL_REPLICATES, alpha=0.05)
    events = D1.load_seed_events()
    rng = np.random.default_rng(RNG_SEED)
    step_days = STEP_MINUTES / 1440.0
    n_half = int(round(HALF_SPAN_DAYS / step_days))
    rate_k = int(round(RATE_HALF_MINUTES / STEP_MINUTES))
    print("D-1e: %d events, %d declared statistics, power floor %.4f"
          % (len(events), len(STAT_NAMES), floor["floor"]), flush=True)

    feats = []
    for e in events:
        t = (e["t_ms"] / 86400000.0
             + step_days * np.arange(-n_half, n_half + 1, dtype=np.float64))
        f = site_features(t, e["lat"], e["lon"], rate_k)
        if not f["ok"][n_half]:
            raise SystemExit("event %s not scoreable at grid centre" % e["id"])
        feats.append(f)

    mats = {k: np.stack([f[k].astype(np.float32) for f in feats])
            for k in list(ANGLES) + list(PARITY)}
    n_sites = len(feats)
    obs = battery(mats, np.full(n_sites, n_half, dtype=np.int64))
    pool = DN.time_uniform_pool(feats[0]["ok"].size,
                                rate_k + int(round(1.5 * 1440.0 / STEP_MINUTES)))
    idx = DN.draw_replicates(pool, n_sites, N_NULL_REPLICATES, rng)
    bad = sum(int((~f["ok"][idx[s]]).sum()) for s, f in enumerate(feats))
    if bad:
        raise SystemExit("null pool hit %d invalid samples" % bad)
    nb = battery(mats, idx)
    cal = DN.calibrate(obs, {k: np.asarray(nb[k], float) for k in STAT_NAMES},
                       STAT_NAMES)

    # the observed spectrum per angle, as a spectrum rather than as loose numbers
    spectra = {}
    for a in ANGLES:
        th = np.array([float(mats[a][s, n_half]) for s in range(n_sites)])
        sp = CS.harmonic_spectrum(th, max_order=max(ORDERS))
        spectra[a] = {"R": sp["R"],
                      "preferred_direction_deg": {m: np.degrees(v) for m, v
                                                  in sp["preferred_direction_rad"].items()},
                      "dominant": CS.dominant_order(sp)}

    out = {
        "arm": "D-1e full angular harmonic spectrum on the K-092 seed set",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 0,
        "scope": {"input": "K092_seed_exclusion_superset.csv ONLY",
                  "sha256_verified": D1.SEED_CSV_SHA256, "n_events": len(events),
                  "licence": "K-096 first-move licence; already-seen set; D-7 STRUCK"},
        "why": ("a signed directional test reports null exactly when a SYMMETRIC "
                "mechanism is total: ten from the left and ten from the right gives "
                "R1 = 0 and R2 = 1. This arm measures which harmonic carries the "
                "concentration instead of assuming one."),
        "what_was_already_covered": (
            "D-1c BEARING_R2 and D-1d BEARING_RESID_R2 were BOTH already order-2 "
            "statistics on the bearing relative to strike, so those nulls stand for "
            "the symmetric case. Never tested before now: order 1, order 3, order 4 "
            "(conjugate shear), the body azimuths at any order, and the elevation "
            "parity split."),
        "declaration": {
            "strike_deg": STRIKE_DEG, "angles": list(ANGLES), "orders": list(ORDERS),
            "parity": list(PARITY), "statistics": list(STAT_NAMES),
            "n_declared_tests": len(STAT_NAMES),
            "half_span_days": HALF_SPAN_DAYS, "step_minutes": STEP_MINUTES,
            "n_null_replicates": N_NULL_REPLICATES, "rng_seed": RNG_SEED,
            "power_floor": floor,
            "null": "waveform-matched time-uniform (engine/dwell_null.py); a "
                    "uniform-circle null would manufacture an effect",
        },
        "observed_spectra": spectra,
        "calibration": cal,
        "not_evidence": ("selected set. A null extends the bound into the symmetric "
                         "and higher-order directions; a hit is a reason to build the "
                         "lattice, not a finding."),
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    print("\n" + "=" * 78)
    print("D-1e: WHICH HARMONIC CARRIES IT?")
    print("=" * 78)
    for a in ANGLES:
        print("  %s" % a)
        for m in ORDERS:
            r = cal["per_statistic"]["%s_m%d" % (a, m)]
            print("     m=%d  R=%.4f  null %.4f +- %.4f   z %+6.2f   [%s]"
                  % (m, r["observed"], r["null_mean"], r["null_sd"], r["z"],
                     CS.SYMMETRY_MEANING[m].split(":")[0]))
        print("     dominant order (observed): %d" % spectra[a]["dominant"]["order"])
    for p in PARITY:
        r = cal["per_statistic"][p]
        print("  %-22s obs %.4f  null %.4f +- %.4f   z %+6.2f"
              % (p, r["observed"], r["null_mean"], r["null_sd"], r["z"]))
    print("\n  max|z| = %.3f ; null 95th pct %.3f ; MAX-STATISTIC p = %.4f"
          % (cal["max_abs_z_observed"], cal["max_abs_z_null_p95"],
             cal["max_statistic_p"]))
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
