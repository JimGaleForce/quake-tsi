"""WORLD SCAN, FAULT-RELATIVE: every angle measured against the event's OWN interface
geometry from Slab2, not against north. EXPLORATION SPLIT ONLY. PRICED.

WHY. The geographic world scan (exp_world_harmonics.py) declared its own largest
weakness: "an effect organised relative to the FAULT rather than to geography is
DILUTED here." A fault plane does not care which way is north, so if tidal orientation
matters at all it should matter relative to the plane. Slab2 now supplies a strike and
dip for every subduction event with zero free parameters and zero catalogue input into
the geometry, so the diluted coordinate can be replaced with the physical one.

---------------------------------------------------------------------------
THE DECLARATION. Fixed BEFORE the first number.
---------------------------------------------------------------------------

EVENTS: the same world catalogue, M >= 5.0, exploration split (pre-2017-02-14), the
holdout untouched -- PLUS a Slab2 interface assignment with DEPTH MISFIT <= 20 km. The
misfit cut is declared at 20 km following K-101, and its purpose is to keep the analysis
to plausible INTERFACE events: an event on a splay, outer-rise or intraslab plane would
otherwise inherit a geometry that is not its own.

DECLUSTERED IS PRIMARY (30 days / 150 km, the rule section P7-25 froze), full catalogue
SECONDARY. Both scored, neither substituting. The geographic scan's only survivor died
to exactly this correction, so it leads here rather than trailing.

GEOMETRY PER EVENT: strike and dip from Slab2; rake 90 (thrust) DECLARED, because the
interface is a thrust and Slab2 does not supply a rake; friction 0.4 DECLARED. Per
K-111, friction and dip are algebraically inert at the free surface, so no grid over
them appears -- gridding them would repeat a known-vacuous design.

STATISTICS, 13 per region:

  B1  AXIAL angle between the tidal principal-extension axis and the fault trace,
      (bearing - strike) mod 180. BOTH quantities are undirected lines in map view, so
      the relation is pi-periodic and its FUNDAMENTAL order is 2 in the original
      coordinate. It is tested by doubling into the full circle and taking harmonics
      f = 1..4, which correspond to ORIGINAL orders 2, 4, 6, 8. Order 1 in the original
      is not tested here because it is not defined for two axial quantities.
      **f = 1 (original order 2) is Jim's case: pull from either side of the fault
      separates it, and a signed test is blind to exactly that.**

  B2  DIRECTIONAL angle between lunar azimuth and the fault strike, (moon_az - strike)
      mod 360, harmonics m = 1..4. Here order 1 IS meaningful and physical, because the
      interface has a DIP DIRECTION: updip and downdip are genuinely different, and a
      strike is a true direction rather than an axis.

  C1  resolved Coulomb on the event's own plane < 0
  C2  d(Coulomb)/dt > 0
  C3  Coulomb < 0 AND falling -- the fault-relative analogue of the frozen quadrant
  P1  EVEN parity of lunar elevation, |sin(elev)| > 0.5
  P2  ODD parity, sin(elev) > 0

NULL: waveform-matched per event, 200 uniform times within +-10 days at the SAME SITE
with the SAME interface geometry, 5000 replicates. The geometry is held fixed while the
time varies, which is what makes this a test of the tide and not of the tectonics.

MULTIPLICITY: max-statistic across every region and statistic on one shared null
ensemble; per-statistic Stouffer pooling with COCHRAN'S Q, because "some regions use
this and others do not" predicts disagreement.

COVERAGE IS A SUBSET AND THAT IS REPORTED, NOT ABSORBED. Slab2 models subduction zones
only, so Iceland, Turkey, Iran, California and much of the Himalaya drop out. The
regions actually scored are listed in the artifact beside the ones that were lost.

NOT A PROMOTION. A survivor is a candidate for Popper.
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

from engine import circular_symmetry as CS
from engine import dwell_null as DN
from engine import slab2 as SL
from engine import tidal_tensor as TT

DECLUSTERED = os.environ.get("WORLD_DECLUSTER", "1") == "1"
OUT_JSON = HERE / ("results_world_faultrel_declustered.json" if DECLUSTERED
                   else "results_world_faultrel_full.json")

MAG_MIN = 5.0
MISFIT_MAX_KM = 20.0
RAKE_DEG = 90.0
FRICTION = 0.4
EXPLORE_FRAC = 0.70
SPAN_START = _dt.datetime(1995, 1, 1, tzinfo=_dt.timezone.utc)
SPAN_END = _dt.datetime(2026, 8, 9, tzinfo=_dt.timezone.utc)
S_NULL_PER_EVENT = 200
NULL_HALF_WINDOW_DAYS = 10.0
RATE_HALF_MINUTES = 10.0
N_REPLICATES = 5000
RNG_SEED = 20260822
DECLUSTER_DAYS, DECLUSTER_KM = 30.0, 150.0
ORDERS = (1, 2, 3, 4)
FLAGS = ("C1_cfs_below", "C2_cfs_rising", "C3_cfs_quadrant",
         "P1_elev_even", "P2_elev_odd")
STATS = (tuple("B1_axial_f%d" % f for f in ORDERS)
         + tuple("B2_moonaz_m%d" % m for m in ORDERS) + FLAGS)
UNIX_EPOCH_JD = 2440587.5
MIN_EVENTS = 30


def explore_cutoff():
    return SPAN_START + (SPAN_END - SPAN_START) * EXPLORE_FRAC


def _gc_km(a, b, c, d):
    r, p = 6371.0, np.pi / 180.0
    h = (np.sin((c - a) * p / 2) ** 2
         + np.cos(a * p) * np.cos(c * p) * np.sin((d - b) * p / 2) ** 2)
    return 2.0 * r * np.arcsin(np.minimum(1.0, np.sqrt(h)))


def decluster(t, lat, lon, mag):
    order = np.argsort(t)
    keep = []
    for ii in order:
        dep = False
        for jj in reversed(keep):
            if t[ii] - t[jj] > DECLUSTER_DAYS:
                break
            if mag[jj] >= mag[ii] and _gc_km(lat[jj], lon[jj],
                                             lat[ii], lon[ii]) <= DECLUSTER_KM:
                dep = True
                break
        if not dep:
            keep.append(ii)
    return np.asarray(sorted(keep), dtype=int)


def load_region(path):
    cut = explore_cutoff()
    t, la, lo, dp, mg = [], [], [], [], []
    for r in csv.DictReader(open(path, newline="", encoding="utf-8")):
        try:
            m = float(r["mag"])
            if m < MAG_MIN:
                continue
            ts = _dt.datetime.fromisoformat(r["time"].replace("Z", "+00:00"))
            if ts >= cut:
                continue
            d = float(r["depth"])
        except (ValueError, TypeError, KeyError):
            continue
        t.append(ts.timestamp() / 86400.0)
        la.append(float(r["latitude"]))
        lo.append(float(r["longitude"]))
        dp.append(d)
        mg.append(m)
    return tuple(np.asarray(x) for x in (t, la, lo, dp, mg))


def features(t_days, lat, lon, strike, dip):
    """Fault-relative features. Geometry is per event and held fixed as time varies."""
    jd = t_days + UNIX_EPOCH_JD
    st = TT.stress_tensor(jd, lat, lon, 0.0)
    bearing = TT.principal_bearing(st["e_NN"], st["e_EE"], st["e_NE"])
    moon = TT.body_direction(jd, lat, lon, 0.0, "moon")
    cfs = TT.coulomb(st, strike, dip, RAKE_DEG, FRICTION)["coulomb_pa"]
    dt_d = RATE_HALF_MINUTES / 1440.0
    hi = TT.coulomb(TT.stress_tensor(jd + dt_d, lat, lon, 0.0),
                    strike, dip, RAKE_DEG, FRICTION)["coulomb_pa"]
    lo_ = TT.coulomb(TT.stress_tensor(jd - dt_d, lat, lon, 0.0),
                     strike, dip, RAKE_DEG, FRICTION)["coulomb_pa"]
    dcfs = (hi - lo_) / (2.0 * RATE_HALF_MINUTES / 60.0)
    d2r = np.pi / 180.0
    # B1: both are undirected lines -> pi-periodic -> double into the full circle
    axial = np.mod(bearing - strike, 180.0) * d2r * 2.0
    b2 = np.mod(moon["azimuth_deg"] - strike, 360.0) * d2r
    return {
        "B1_axial": axial,
        "B2_moonaz": b2,
        "C1_cfs_below": (cfs < 0.0).astype(np.float32),
        "C2_cfs_rising": (dcfs > 0.0).astype(np.float32),
        "C3_cfs_quadrant": ((cfs < 0.0) & (dcfs < 0.0)).astype(np.float32),
        "P1_elev_even": (np.abs(moon["sin_elevation"]) > 0.5).astype(np.float32),
        "P2_elev_odd": (moon["sin_elevation"] > 0.0).astype(np.float32),
    }


def battery(f, idx=None):
    out = {}
    for key, pref, orders in (("B1_axial", "B1_axial_f", ORDERS),
                              ("B2_moonaz", "B2_moonaz_m", ORDERS)):
        th = f[key] if idx is None else np.take_along_axis(f[key], idx, axis=1)
        for m in orders:
            c, s = np.cos(m * th), np.sin(m * th)
            if idx is None:
                out["%s%d" % (pref, m)] = float(np.hypot(c.mean(), s.mean()))
            else:
                out["%s%d" % (pref, m)] = np.hypot(c.mean(axis=0), s.mean(axis=0))
    for k in FLAGS:
        v = f[k] if idx is None else np.take_along_axis(f[k], idx, axis=1)
        out[k] = float(v.mean()) if idx is None else v.mean(axis=0)
    return out


def main():
    rng = np.random.default_rng(RNG_SEED)
    paths = sorted(glob.glob(str(HERE / "data" / "comcat_world" / "*.csv")))
    prepared, lost = {}, {}
    for path in paths:
        name = os.path.basename(path)[:-4]
        t, la, lo, dp, mg = load_region(path)
        if t.size == 0:
            lost[name] = "no events"
            continue
        a = SL.assign(la, lo, dp)
        ok = a["assigned"] & (a["depth_misfit_km"] <= MISFIT_MAX_KM)
        n_before = int(t.size)
        t, la, lo, dp, mg = t[ok], la[ok], lo[ok], dp[ok], mg[ok]
        strike, dip = a["strike_deg"][ok], a["dip_deg"][ok]
        if DECLUSTERED and t.size:
            k = decluster(t, la, lo, mg)
            t, la, lo, strike, dip = t[k], la[k], lo[k], strike[k], dip[k]
        if t.size < MIN_EVENTS:
            lost[name] = ("%d of %d events survive the Slab2 + misfit + decluster cuts"
                          % (int(t.size), n_before))
            continue
        prepared[name] = (t, la, lo, strike, dip, n_before)

    n_stats = len(prepared) * len(STATS)
    floor = DN.assert_power_floor(n_stats, N_REPLICATES, alpha=0.05)
    print("FAULT-RELATIVE WORLD SCAN%s. %d regions x %d statistics = %d declared "
          "tests ; power floor %.4f"
          % (" (DECLUSTERED)" if DECLUSTERED else "", len(prepared), len(STATS),
             n_stats, floor["floor"]), flush=True)

    regions = {}
    for name, (t, la, lo, strike, dip, n_before) in prepared.items():
        n = t.size
        obs_f = features(t, la, lo, strike, dip)
        obs = battery(obs_f)
        off = rng.uniform(-NULL_HALF_WINDOW_DAYS, NULL_HALF_WINDOW_DAYS,
                          size=(n, S_NULL_PER_EVENT))
        nf = features((t[:, None] + off).ravel(),
                      np.repeat(la, S_NULL_PER_EVENT), np.repeat(lo, S_NULL_PER_EVENT),
                      np.repeat(strike, S_NULL_PER_EVENT),
                      np.repeat(dip, S_NULL_PER_EVENT))
        nf = {k: v.reshape(n, S_NULL_PER_EVENT) for k, v in nf.items()}
        idx = rng.integers(0, S_NULL_PER_EVENT, size=(n, N_REPLICATES))
        nb = battery(nf, idx=idx)
        z, cols = {}, {}
        for k in STATS:
            arr = np.asarray(nb[k], dtype=np.float64)
            mu, sd = float(arr.mean()), float(arr.std(ddof=1))
            z[k] = {"observed": float(obs[k]), "null_mean": mu, "null_sd": sd,
                    "z": (float(obs[k]) - mu) / sd if sd > 0 else float("nan")}
            cols[k] = (arr - mu) / (sd if sd > 0 else 1e-300)
        top = max(STATS, key=lambda k: abs(z[k]["z"]))
        print("  %-18s n=%5d (of %5d)  largest |z| = %5.2f  (%s)"
              % (name, n, n_before, abs(z[top]["z"]), top), flush=True)
        regions[name] = {"n": int(n), "n_before_cuts": int(n_before),
                         "per_statistic": z, "_cols": cols}

    allcols = [regions[r]["_cols"][k] for r in regions for k in STATS]
    gw = np.max(np.abs(np.stack(allcols)), axis=0)
    obs_max = max(abs(regions[r]["per_statistic"][k]["z"])
                  for r in regions for k in STATS)
    p_gw = (int(np.sum(gw >= obs_max)) + 1) / (N_REPLICATES + 1)
    worst = max(((r, k, regions[r]["per_statistic"][k]["z"])
                 for r in regions for k in STATS), key=lambda x: abs(x[2]))
    combined = {k: CS.combine_regions(
        {r: regions[r]["per_statistic"][k]["z"] for r in regions},
        {r: regions[r]["n"] for r in regions}) for k in STATS}
    for r in regions:
        regions[r].pop("_cols", None)

    out = {
        "arm": "WORLD SCAN, FAULT-RELATIVE (Slab2 geometry per event)",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "PRICED": True, "declustered": DECLUSTERED,
        "n_declared_tests": n_stats,
        "declaration": {
            "regions_scored": sorted(regions), "regions_lost": lost,
            "mag_min": MAG_MIN, "misfit_max_km": MISFIT_MAX_KM,
            "rake_deg": RAKE_DEG, "friction": FRICTION,
            "geometry_source": SL.CITATION, "geometry_scope": SL.SCOPE,
            "statistics": list(STATS), "orders": list(ORDERS),
            "B1_note": ("axial: both the tidal extension axis and the fault trace are "
                        "undirected lines, so the relation is pi-periodic and is "
                        "doubled into the full circle. f = 1, 2, 3, 4 here correspond "
                        "to ORIGINAL symmetry orders 2, 4, 6, 8. f = 1 is the case "
                        "where pull from either side separates the fault, which a "
                        "signed test cannot see."),
            "B2_note": ("directional: the interface has a DIP DIRECTION, so updip and "
                        "downdip differ and order 1 is physically meaningful here in a "
                        "way it is not for B1."),
            "explore_cutoff_utc": explore_cutoff().isoformat(),
            "holdout": "NOT READ",
            "n_replicates": N_REPLICATES, "rng_seed": RNG_SEED,
            "power_floor": floor,
        },
        "per_region": regions,
        "combined_across_regions": combined,
        "max_statistic": {
            "observed_max_abs_z": float(obs_max),
            "null_max_p95": float(np.quantile(gw, 0.95)),
            "p": float(p_gw),
            "where": {"region": worst[0], "statistic": worst[1],
                      "z": float(worst[2])},
        },
        "not_a_promotion": "a survivor is a candidate for Popper, not a finding",
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    print("\n" + "=" * 78)
    print("FAULT-RELATIVE RESULT%s" % (" (DECLUSTERED, PRIMARY)" if DECLUSTERED
                                       else " (FULL, SECONDARY)"))
    print("=" * 78)
    print("  events %d across %d regions ; regions lost to Slab2 coverage: %s"
          % (sum(r["n"] for r in regions.values()), len(regions),
             ", ".join(sorted(lost)) or "none"))
    print("  MAX-STATISTIC |z| = %.3f at %s / %s ; null 95th %.3f ; p = %.4f"
          % (obs_max, worst[0], worst[1], np.quantile(gw, 0.95), p_gw))
    print("\n  strongest pooled effects:")
    for k in sorted(STATS, key=lambda k: -abs(combined[k]["z_pooled"]))[:6]:
        c = combined[k]
        print("    %-18s pooled z=%+6.2f   Q p=%.4g %s"
              % (k, c["z_pooled"], c["Q_p"],
                 "HETEROGENEOUS" if c["heterogeneous"] else ""))
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
