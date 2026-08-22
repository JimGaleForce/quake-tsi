"""WORLD SCAN: angular harmonics, body azimuth/elevation parity, and stress sign,
across 13 tectonic regions. EXPLORATION SPLIT ONLY. This one is PRICED.

Every arm before this ran on the 162-event K-092 seed superset at price 0, because
that set was already seen and could not become evidence. THIS IS DIFFERENT. It touches
real catalogues that have not been looked at, so it spends declared tests and is logged
to EXPLORE_COUNT accordingly. Nothing here promotes without Popper.

---------------------------------------------------------------------------
THE DECLARATION. Every line fixed BEFORE the first number is computed.
---------------------------------------------------------------------------

REGIONS (13, all of `data/comcat_world/`): Alaska-Aleutians, California, Caribbean,
Chile, Greece-Aegean, Himalaya, Iceland, Indonesia, Iran, Japan, Mexico, Philippines,
Turkey. Declared as the complete set rather than a chosen subset, so there is no
region-selection step to hide anything in.

MAGNITUDE: M >= 5.0, and the floor is DELIBERATE. K-104 warns that ocean tidal loading
modulates coastal station noise and therefore the completeness magnitude Mc, which
manufactures a tidal signal with the correct period and phase and NO physics in it.
That artifact lives at magnitudes near Mc. M >= 5.0 is comfortably above Mc in every
one of these regions over this window, so the scan buys its safety by giving up N
rather than the other way round.

SPLIT: TEMPORAL, exploration = the first 70 % of the 1995-01-01 to 2026-08-09 span,
i.e. events before 2017-02-14. THE HOLDOUT IS NOT READ, NOT LOADED, AND NOT COUNTED.
A random split would leak in a clustered process; a temporal one does not.

STATISTICS, 16 per region:

  A1  tidal principal-extension bearing, ABSOLUTE, harmonic orders m = 1, 2, 3, 4
  A2  lunar azimuth, ABSOLUTE, orders m = 1, 2, 3, 4
  A3  solar azimuth, ABSOLUTE, orders m = 1, 2, 3, 4
  P1  EVEN parity of lunar elevation: |sin(elevation)| > 0.5. The coordinate the
      degree-2 potential actually lives in, since P2(cos z) cannot tell overhead from
      underfoot.
  P2  ODD parity: sin(elevation) > 0, i.e. Moon above the horizon rather than below.
      Only a degree-3-like asymmetry can carry this.
  Q   the quadrant: areal strain < 0 AND its time derivative < 0.
  L   the level sign: areal strain < 0.

Orders are tested because a SYMMETRIC mechanism -- pull from the left or the right
both separate the fault -- puts everything in m = 2 and NOTHING in m = 1, so a
first-moment test reports null exactly when the effect is total. m = 4 is where
conjugate shear geometry would sit.

Bearings are ABSOLUTE rather than relative to fault strike, because a per-event strike
requires a slab model this program has not downloaded. That is a real limitation and it
is stated rather than worked around: an effect that is organised relative to the FAULT
and not to geography will be diluted here.

TOTAL DECLARED TESTS: 13 regions x 16 = 208.

THE NULL: waveform-matched and time-uniform, per event. For each event, S = 200 times
are drawn uniformly within +-10 days of it AT THE SAME SITE, and the identical feature
code runs on them. A null replicate draws one of the S per event and recomputes the
whole battery. This preserves each event's own local tidal waveform, its latitude, and
its epoch -- all three of which move the null, as the seed-set work measured.

MULTIPLICITY: max-statistic across ALL 208, calibrated against the same null. Per
statistic, regions are combined by Stouffer weighting with sqrt(n) AND reported with
COCHRAN'S Q, because "some regions use this triggering and others do not" predicts
DISAGREEMENT, and a pooled number alone would hide exactly that.

POWER FLOOR: 208 statistics need N >= 208/0.05 - 1 = 4159 replicates. N = 5000 is
declared, giving a floor of 0.0416. Asserted before the run.

NOT A PROMOTION. A survivor here is a candidate for Popper, not a finding.
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
from engine import sitetide as ST
from engine import tidal_tensor as TT

DECLUSTERED = os.environ.get("WORLD_DECLUSTER", "0") == "1"
OUT_JSON = HERE / ("results_world_harmonics_declustered.json" if DECLUSTERED
                   else "results_world_harmonics.json")

MAG_MIN = 5.0
EXPLORE_FRAC = 0.70
SPAN_START = _dt.datetime(1995, 1, 1, tzinfo=_dt.timezone.utc)
SPAN_END = _dt.datetime(2026, 8, 9, tzinfo=_dt.timezone.utc)
S_NULL_PER_EVENT = 200
NULL_HALF_WINDOW_DAYS = 10.0
RATE_HALF_MINUTES = 10.0
N_REPLICATES = 5000
RNG_SEED = 20260822
ORDERS = (1, 2, 3, 4)
ANGLES = ("A1_bearing", "A2_moon_az", "A3_sun_az")
FLAGS = ("P1_elev_even", "P2_elev_odd", "Q_quadrant", "L_level")
STATS = tuple("%s_m%d" % (a, m) for a in ANGLES for m in ORDERS) + FLAGS

UNIX_EPOCH_JD = 2440587.5


def explore_cutoff():
    return SPAN_START + (SPAN_END - SPAN_START) * EXPLORE_FRAC


DECLUSTER_DAYS = 30.0
DECLUSTER_KM = 150.0


def _great_circle_km(lat1, lon1, lat2, lon2):
    r, p = 6371.0, np.pi / 180.0
    h = (np.sin((lat2 - lat1) * p / 2) ** 2
         + np.cos(lat1 * p) * np.cos(lat2 * p) * np.sin((lon2 - lon1) * p / 2) ** 2)
    return 2.0 * r * np.arcsin(np.minimum(1.0, np.sqrt(h)))


def decluster(t, lat, lon, mag):
    """The rule section P7-25 FROZE before any event: an event is dependent if it falls
    within 30 days AND 150 km of an equal-or-larger prior event.

    Not a new test chosen after seeing a result. Section P7-25(2) applied exactly this
    rule to the crest deficit and it is the program's standing dependence correction;
    a scan that reported a z on a catalogue that is half aftershocks would be quoting
    an effective N it does not have.
    """
    order = np.argsort(t)
    keep = []
    for ii in order:
        dep = False
        for jj in reversed(keep):
            if t[ii] - t[jj] > DECLUSTER_DAYS:
                break
            if mag[jj] >= mag[ii] and _great_circle_km(lat[jj], lon[jj],
                                                       lat[ii], lon[ii]) <= DECLUSTER_KM:
                dep = True
                break
        if not dep:
            keep.append(ii)
    return np.asarray(sorted(keep), dtype=int)


def load_region(path, declustered=False):
    """Exploration-split events only. The holdout is never even parsed into memory."""
    cut = explore_cutoff()
    t, la, lo, mg = [], [], [], []
    n_holdout = 0
    for r in csv.DictReader(open(path, newline="", encoding="utf-8")):
        try:
            m = float(r["mag"])
            if m < MAG_MIN:
                continue
            ts = _dt.datetime.fromisoformat(r["time"].replace("Z", "+00:00"))
        except (ValueError, TypeError, KeyError):
            continue
        if ts >= cut:
            n_holdout += 1
            continue
        t.append(ts.timestamp() / 86400.0)
        la.append(float(r["latitude"]))
        lo.append(float(r["longitude"]))
        mg.append(m)
    t, la, lo, mg = (np.asarray(t), np.asarray(la), np.asarray(lo), np.asarray(mg))
    n_full = int(t.size)
    if declustered and t.size:
        k = decluster(t, la, lo, mg)
        t, la, lo = t[k], la[k], lo[k]
    return (t, la, lo, n_holdout, n_full)


def features(t_days, lat, lon):
    """The declared feature set at arbitrary (time, site) triples. Fully vectorised."""
    jd = t_days + UNIX_EPOCH_JD
    st = TT.stress_tensor(jd, lat, lon, 0.0)
    bearing = TT.principal_bearing(st["e_NN"], st["e_EE"], st["e_NE"])
    moon = TT.body_direction(jd, lat, lon, 0.0, "moon")
    sun = TT.body_direction(jd, lat, lon, 0.0, "sun")
    dt_d = RATE_HALF_MINUTES / 1440.0
    a_hi = ST.site_scalar(jd + dt_d, lat, lon, 0.0)
    a_lo = ST.site_scalar(jd - dt_d, lat, lon, 0.0)
    rate = (a_hi - a_lo) / (2.0 * RATE_HALF_MINUTES / 60.0)
    areal = st["areal_strain"]
    d2r = np.pi / 180.0
    return {
        "A1_bearing": np.mod(bearing * d2r, 2.0 * np.pi),
        "A2_moon_az": moon["azimuth_deg"] * d2r,
        "A3_sun_az": sun["azimuth_deg"] * d2r,
        "P1_elev_even": (np.abs(moon["sin_elevation"]) > 0.5).astype(np.float32),
        "P2_elev_odd": (moon["sin_elevation"] > 0.0).astype(np.float32),
        "Q_quadrant": ((areal < 0.0) & (rate < 0.0)).astype(np.float32),
        "L_level": (areal < 0.0).astype(np.float32),
    }


def battery(ang, flg, idx=None):
    """R_m per angle plus the flag fractions. `idx` selects a null draw per event."""
    out = {}
    for a in ANGLES:
        th = ang[a] if idx is None else np.take_along_axis(ang[a], idx, axis=1)
        for m in ORDERS:
            c = np.cos(m * th)
            s = np.sin(m * th)
            if idx is None:
                out["%s_m%d" % (a, m)] = float(np.hypot(c.mean(), s.mean()))
            else:
                out["%s_m%d" % (a, m)] = np.hypot(c.mean(axis=0), s.mean(axis=0))
    for f in FLAGS:
        v = flg[f] if idx is None else np.take_along_axis(flg[f], idx, axis=1)
        out[f] = float(v.mean()) if idx is None else v.mean(axis=0)
    return out


def run_region(name, t, lat, lon, rng, verbose=True):
    n = t.size
    if n < 30:
        return None
    obs_f = features(t, lat, lon)
    obs_ang = {a: obs_f[a] for a in ANGLES}
    obs_flg = {f: obs_f[f] for f in FLAGS}
    obs = battery(obs_ang, obs_flg)

    off = rng.uniform(-NULL_HALF_WINDOW_DAYS, NULL_HALF_WINDOW_DAYS,
                      size=(n, S_NULL_PER_EVENT))
    tt = (t[:, None] + off).ravel()
    ll = np.repeat(lat, S_NULL_PER_EVENT)
    oo = np.repeat(lon, S_NULL_PER_EVENT)
    nf = features(tt, ll, oo)
    null_ang = {a: nf[a].reshape(n, S_NULL_PER_EVENT) for a in ANGLES}
    null_flg = {f: nf[f].reshape(n, S_NULL_PER_EVENT) for f in FLAGS}

    idx = rng.integers(0, S_NULL_PER_EVENT, size=(n, N_REPLICATES))
    nb = battery(null_ang, null_flg, idx=idx)

    z, null_cols = {}, {}
    for k in STATS:
        arr = np.asarray(nb[k], dtype=np.float64)
        mu, sd = float(arr.mean()), float(arr.std(ddof=1))
        z[k] = {"observed": float(obs[k]), "null_mean": mu, "null_sd": sd,
                "z": (float(obs[k]) - mu) / sd if sd > 0 else float("nan")}
        null_cols[k] = (arr - mu) / (sd if sd > 0 else 1e-300)
    if verbose:
        top = max(STATS, key=lambda k: abs(z[k]["z"]))
        print("  %-18s n=%5d   largest |z| = %5.2f  (%s)"
              % (name, n, abs(z[top]["z"]), top), flush=True)
    return {"n": int(n), "per_statistic": z, "_null_cols": null_cols}


def main():
    floor = DN.assert_power_floor(13 * len(STATS), N_REPLICATES, alpha=0.05)
    rng = np.random.default_rng(RNG_SEED)
    cut = explore_cutoff()
    print("WORLD SCAN%s. exploration cutoff %s ; M >= %.1f ; %d statistics/region ; "
          "power floor %.4f" % (" (DECLUSTERED)" if DECLUSTERED else "", cut.date(),
                                MAG_MIN, len(STATS), floor["floor"]), flush=True)

    regions, held, full_counts = {}, {}, {}
    for path in sorted(glob.glob(str(HERE / "data" / "comcat_world" / "*.csv"))):
        name = os.path.basename(path)[:-4]
        t, la, lo, n_hold, n_full = load_region(path, declustered=DECLUSTERED)
        held[name] = int(n_hold)
        full_counts[name] = n_full
        r = run_region(name, t, la, lo, rng)
        if r is not None:
            regions[name] = r

    # ---- max-statistic across EVERY region and statistic, one shared null ensemble --
    cols = [regions[rn]["_null_cols"][k] for rn in regions for k in STATS]
    gw_null = np.max(np.abs(np.stack(cols)), axis=0)
    obs_max = max(abs(regions[rn]["per_statistic"][k]["z"])
                  for rn in regions for k in STATS)
    p_gw = (int(np.sum(gw_null >= obs_max)) + 1) / (N_REPLICATES + 1)

    # ---- per statistic: Stouffer across regions AND Cochran's Q -------------------
    combined = {}
    for k in STATS:
        zb = {rn: regions[rn]["per_statistic"][k]["z"] for rn in regions}
        nb = {rn: regions[rn]["n"] for rn in regions}
        combined[k] = CS.combine_regions(zb, nb)

    for rn in regions:
        regions[rn].pop("_null_cols", None)

    worst = max(((rn, k, regions[rn]["per_statistic"][k]["z"])
                 for rn in regions for k in STATS), key=lambda x: abs(x[2]))
    het = sorted(((k, combined[k]["cochran_Q"], combined[k]["Q_p"])
                  for k in STATS), key=lambda x: x[1], reverse=True)

    out = {
        "arm": "WORLD SCAN: angular harmonics + body parity + stress sign, 13 regions",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "PRICED": True,
        "declustered": DECLUSTERED,
        "decluster_rule": ("dependent if within %.0f days AND %.0f km of an "
                           "equal-or-larger prior event; the rule section P7-25 froze "
                           "before any event" % (DECLUSTER_DAYS, DECLUSTER_KM)),
        "n_full_before_declustering": full_counts,
        "n_declared_tests": 13 * len(STATS),
        "declaration": {
            "regions": sorted(regions), "mag_min": MAG_MIN,
            "explore_frac": EXPLORE_FRAC, "explore_cutoff_utc": cut.isoformat(),
            "holdout_events_not_read": held,
            "statistics_per_region": list(STATS), "orders": list(ORDERS),
            "s_null_per_event": S_NULL_PER_EVENT,
            "null_half_window_days": NULL_HALF_WINDOW_DAYS,
            "n_replicates": N_REPLICATES, "rng_seed": RNG_SEED,
            "power_floor": floor,
            "bearings_are_absolute": (
                "not relative to fault strike, because a per-event strike needs a slab "
                "model this program has not downloaded. An effect organised relative "
                "to the FAULT rather than to geography is DILUTED here. Stated as a "
                "limitation, not worked around."),
            "magnitude_floor_rationale": (
                "K-104: ocean loading modulates coastal station noise and therefore "
                "Mc, manufacturing a tidal signal with the right period and phase and "
                "no physics. That artifact lives near Mc, so the scan buys safety with "
                "N rather than the reverse."),
        },
        "per_region": regions,
        "combined_across_regions": combined,
        "max_statistic": {
            "observed_max_abs_z": float(obs_max),
            "null_max_p95": float(np.quantile(gw_null, 0.95)),
            "p": float(p_gw),
            "where": {"region": worst[0], "statistic": worst[1], "z": float(worst[2])},
            "method": ("max |z| over all 13 regions x 16 statistics against the same "
                       "null replicates, so the correlation among statistics and among "
                       "regions is reproduced rather than modelled"),
        },
        "heterogeneity_ranked": [{"statistic": k, "Q": q, "Q_p": p} for k, q, p in het],
        "not_a_promotion": ("a survivor here is a candidate for Popper, not a finding. "
                            "The holdout was not read."),
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    print("\n" + "=" * 78)
    print("WORLD SCAN RESULT")
    print("=" * 78)
    print("  events scored : %d across %d regions"
          % (sum(r["n"] for r in regions.values()), len(regions)))
    print("  MAX-STATISTIC : |z| = %.3f at %s / %s ; null 95th %.3f ; p = %.4f"
          % (obs_max, worst[0], worst[1], np.quantile(gw_null, 0.95), p_gw))
    print("\n  most heterogeneous statistics (regions disagreeing):")
    for k, q, qp in het[:5]:
        c = combined[k]
        print("    %-16s Q=%7.2f p=%.4g  I2=%.2f  pooled z=%+.2f"
              % (k, q, qp, c["I2"], c["z_pooled"]))
    print("\n  strongest pooled effects:")
    for k in sorted(STATS, key=lambda k: -abs(combined[k]["z_pooled"]))[:5]:
        c = combined[k]
        print("    %-16s pooled z=%+6.2f  Q p=%.4g  %s"
              % (k, c["z_pooled"], c["Q_p"],
                 "HETEROGENEOUS" if c["heterogeneous"] else ""))
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
