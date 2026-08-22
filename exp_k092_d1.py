"""D-1: the declared arcsine / trough-vs-mid-slope readout of the K-092 seed set.

WHAT THIS IS (HYPOTHESIS_LEDGER.md P7-23(C), K092_FREEZE.md pre-scoring gate 2):

  > D-1's declared readout must, for each seed event, report BOTH the level percentile
  > within the local tidal range AND the phase angle, and classify the seed set as
  > TROUGH-CONCENTRATED (artifact-consistent) or MID-SLOPE-CONCENTRATED
  > (artifact-inconsistent), against a classification rule declared in advance --
  > BEFORE any Alaska event outside the seed set is scored.

The rule was declared in `engine/audit_arcsine.py` (bands u < -1/2 / |u| <= 1/2 /
u > +1/2, exactly equiprobable at 1/3 under uniform phase; Z_CLASSIFY = 2.0; quadrant
null 1/4) and committed BEFORE this script existed. This script does not restate the
rule; it imports it, so it cannot drift from it.

SCOPE, WHICH IS ITS ENTIRE LICENCE. D-1 reads `K092_seed_exclusion_superset.csv` and
NOTHING ELSE. Those 162 events are the ALREADY-SEEN set -- the superset of what Jim
scrolled past before the freeze (P7-23(A) "scrolled past is seen"). They are EXCLUDED
from D-12 scoring by construction, which is exactly why they may be looked at here at
price 0. No event outside this file is touched, downloaded or evaluated.

WHAT IT IS NOT. Not evidence. P7-22 Ratification 1 is explicit: if the seed's
observation is quantitatively consistent with the exact null fraction, THE MOTIVATING
OBSERVATION IS EXPLAINED and the entry's priority collapses -- and that is the result,
not a preliminary. The reverse reading is NOT symmetric: a mid-slope concentration on
a set selected by having been looked at is not a positive finding, it only fails to
explain the seed away. D-12 (held-out) and D-13 (prospective, clock running) are the
arms that can produce evidence.

SCORER. The FROZEN scalar is the observing application's, so the primary scorer is the
app's own code, unmodified, through `exp_k092_d1_bridge.mjs` -- the exactness path
recommended by K092_SCALAR_PROVENANCE.md section 5(3). `engine/sitetide.py` runs the
identical downstream analysis alongside as the declared robustness check. The bridge
refuses to run if it cannot reproduce the frozen 1991 readout to the digit.
"""

from __future__ import annotations

import csv
import datetime as _dt
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from engine import audit_arcsine as AA          # noqa: E402  the DECLARED rule
from engine import sitetide as ST               # noqa: E402  the robustness scorer

# ---------------------------------------------------------------- frozen inputs --
SEED_CSV = HERE / "K092_seed_exclusion_superset.csv"
SEED_CSV_SHA256 = "35b5cf0fe9145ec830fba782550f7d8b76dc2708db105eb6003a0f92b68bfcb6"
ASTRO_URL = "file:///D:/CODE/git/earth-tides-globe/src/utils/astro.ts"
BRIDGE = HERE / "exp_k092_d1_bridge.mjs"
OUT_JSON = HERE / "results_k092_d1.json"

# Local-cycle window. MEASURED, not assumed: at +-1 day, 3 of the 162 seed events
# could not be bracketed by two local maxima -- BOTH scorers dropped the SAME three,
# which identifies it as a waveform property and not an implementation quirk. These
# are high-latitude sites where the diurnal inequality stretches the cycle to as much
# as 25.4 h (measured across the seed set: 9.71 h to 25.41 h, mean 15.51 h), so a
# +-1 day half-window can fail to contain a maximum on each side. +-2 days brackets
# every event. 1 minute is 1.3e-3 rad of a 12.4 h cycle before the parabolic
# refinement, which removes it to O(dt^2/T).
HALF_WINDOW_DAYS = 2.0
STEP_MINUTES = 1.0
RATE_HALF_MINUTES = 10.0        # the app's own +-10 min central difference


# ------------------------------------------------------------ the convention note --
# FOUND WHILE BUILDING D-1, RECORDED RATHER THAN SILENTLY RESOLVED.
#
# `K092_FREEZE.md` names the quadrant twice and the two namings disagree in LABEL:
#
#   (a) the two SIGN CONDITIONS: "(tide level below zero) AND (tide level decreasing)"
#   (b) the angle restatement:   "i.e. theta in (pi, 3*pi/2) with theta = 0 at maximum"
#
# With theta = 0 at the MAXIMUM the scalar reads x = A cos(theta), so
#   below neutral  -> cos(theta) < 0    -> theta in (pi/2, 3pi/2)
#   falling        -> -A sin(theta) < 0 -> theta in (0, pi)
#   intersection   -> theta in (pi/2, pi)          NOT (pi, 3pi/2).
# The interval (pi, 3pi/2) is correct for the OTHER convention, x = A sin(theta) with
# theta = 0 at the ascending zero -- which is the convention `audit_arcsine.py`'s
# QUADRANT_NOTE states and the one its `seed_readout` tests.
#
# Both namings pick out THE SAME PHYSICAL QUARTER-CYCLE (below neutral and falling);
# only the angular label differs, and the uniform-phase null probability is exactly
# 1/4 under either. Nothing about the claim, the freeze or the null moves.
#
# HOW THIS SCRIPT HANDLES IT, declared here: the PRIMARY quadrant readout is the two
# SIGN CONDITIONS evaluated directly on the frozen scalar and its +-10 min central
# difference -- the form that is convention-free, the form the freeze leads with, and
# the form K092_SCALAR_PROVENANCE.md section 4 already measured across implementations.
# The angle-interval readout is reported ALONGSIDE, with phases converted into
# audit_arcsine's sin convention (theta_sin = theta_max0 + pi/2) so that its committed
# (pi, 3pi/2) test is applied as written. The two agree only up to the real waveform's
# asymmetry, and their disagreement count is reported rather than assumed to be zero.
CONVENTION_NOTE = (
    "K092_FREEZE.md names the quadrant both as the two sign conditions (level < 0 AND "
    "level falling) and as 'theta in (pi, 3pi/2) with theta = 0 at maximum'. Those two "
    "labels are inconsistent: with theta = 0 at the maximum the scalar is A cos(theta) "
    "and the below-neutral-falling quarter is theta in (pi/2, pi); (pi, 3pi/2) is that "
    "same physical quarter-cycle under the x = A sin(theta) convention that "
    "engine/audit_arcsine.py uses. Same quarter cycle, same exact 1/4 null, different "
    "angular label. D-1 scores the quadrant PRIMARILY by the sign conditions "
    "(convention-free) and reports the angle-interval form alongside with phases "
    "mapped into the sin convention so audit_arcsine's committed test applies as "
    "written.")


# --------------------------------------------------------------------- utilities --
def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_iso_utc(s: str) -> _dt.datetime:
    return _dt.datetime.fromisoformat(s.replace("Z", "+00:00"))


def load_seed_events():
    """The 162 seed-superset events. Verifies the frozen sha256 before reading."""
    got = sha256_of(SEED_CSV)
    if got != SEED_CSV_SHA256:
        raise SystemExit(
            "REFUSING TO RUN: K092_seed_exclusion_superset.csv sha256 %s does not "
            "match the frozen %s. The seed set is part of the freeze; a changed file "
            "is a different arm." % (got, SEED_CSV_SHA256))
    evs = []
    with SEED_CSV.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            t = parse_iso_utc(row["time"])
            evs.append({
                "id": row["id"],
                "time_utc": row["time"],
                "t_ms": int(round(t.timestamp() * 1000.0)),
                "lat": float(row["latitude"]),
                "lon": float(row["longitude"]),
                "depth_km": float(row["depth"]) if row["depth"] not in ("", None) else 0.0,
                "mag": float(row["mag"]),
                "place": row["place"],
            })
    evs.sort(key=lambda e: e["t_ms"])
    return evs


def refined_maxima(t, x):
    """Times of interior local maxima of `x(t)`, parabolic sub-sample refined.

    Same construction as `sitetide.tidal_maxima`, deliberately mirrored so the app arm
    and the sitetide arm differ ONLY in the scalar and never in the analysis.
    """
    i = np.nonzero((x[1:-1] > x[:-2]) & (x[1:-1] >= x[2:]))[0] + 1
    if i.size < 2:
        return np.zeros(0, dtype=np.float64)
    y0, y1, y2 = x[i - 1], x[i], x[i + 1]
    den = y0 - 2.0 * y1 + y2
    shift = np.where(den != 0.0, 0.5 * (y0 - y2) / np.where(den == 0.0, 1.0, den), 0.0)
    dt = float(t[1] - t[0])
    return t[i] + np.clip(shift, -0.5, 0.5) * dt


def local_cycle_readout(t, x, t_event):
    """Per-event readout on one local cycle of the scalar.

    The cycle is the interval between the two refined local maxima that BRACKET the
    event. Within it:
      hi, lo = the cycle's own max and min, so the percentile is "within the local
               tidal range" -- P7-23(C)'s wording, and also what a person reading a
               tidal display can actually see;
      u      = 2 (x - lo)/(hi - lo) - 1, audit_arcsine's normalisation;
      theta  = 2 pi (t - t_before)/(t_after - t_before), zero AT THE MAXIMUM -- the
               D-0 / Tanaka construction (sitetide.phase_from_maxima), reported in
               that convention and converted to the sin convention for the committed
               quadrant test.
    """
    tm = refined_maxima(t, x)
    if tm.size < 2 or not (tm[0] <= t_event <= tm[-1]):
        return None
    k = int(np.searchsorted(tm, t_event, side="right") - 1)
    if k < 0 or k >= tm.size - 1:
        return None
    tb, ta = float(tm[k]), float(tm[k + 1])
    m = (t >= tb) & (t <= ta)
    if not np.any(m):
        return None
    lo, hi = float(x[m].min()), float(x[m].max())
    xe = float(np.interp(t_event, t, x))
    u = 2.0 * (xe - lo) / max(hi - lo, 1e-300) - 1.0
    theta_max0 = 2.0 * math.pi * (t_event - tb) / max(ta - tb, 1e-300)
    return {
        "level": xe,
        "cycle_min": lo, "cycle_max": hi,
        "cycle_hours": (ta - tb) * 24.0,
        "normalised_level_u": u,
        "level_percentile": 0.5 * (u + 1.0),
        "phase_rad_max0": float(np.mod(theta_max0, 2.0 * math.pi)),
        "phase_rad_sin": float(np.mod(theta_max0 + 0.5 * math.pi, 2.0 * math.pi)),
    }


def binom_p_greater(k, n, p0):
    """One-sided exact binomial P(X >= k | n, p0). The D-13 readout's own test form."""
    from math import comb
    return float(sum(comb(n, j) * p0 ** j * (1.0 - p0) ** (n - j)
                     for j in range(k, n + 1)))


def quadrant_block(levels, rates, phases_sin, label):
    """Both quadrant readouts plus the exact binomial, for one scorer."""
    lv = np.asarray(levels, float)
    rt = np.asarray(rates, float)
    ph = np.asarray(phases_sin, float)
    sign_q = (lv < 0.0) & (rt < 0.0)
    angle_q = (ph > math.pi) & (ph < 1.5 * math.pi)
    n = int(lv.size)
    k = int(sign_q.sum())
    return {
        "scorer": label,
        "n": n,
        "sign_condition_quadrant": {
            "definition": "level < 0 AND d(level)/dt < 0 (+-10 min central difference)",
            "count": k, "fraction": k / max(n, 1),
            "null_fraction": AA.QUADRANT_NULL_FRACTION,
            "binomial_p_one_sided_ge": binom_p_greater(k, n, AA.QUADRANT_NULL_FRACTION),
        },
        "angle_interval_quadrant": {
            "definition": "phase in (pi, 3pi/2) under x = A sin(theta) "
                          "(audit_arcsine's committed test, phases converted)",
            "count": int(angle_q.sum()),
            "fraction": float(angle_q.mean()) if n else float("nan"),
            "null_fraction": AA.QUADRANT_NULL_FRACTION,
        },
        "disagreement_between_the_two": {
            "count": int(np.sum(sign_q != angle_q)),
            "note": "waveform asymmetry: the interpolated-phase quarter and the "
                    "sign-condition quarter coincide only for a pure sinusoid.",
        },
        "convention_note": CONVENTION_NOTE,
    }


# ------------------------------------------------------------------ the two arms --
def app_arm(events):
    """PRIMARY: the frozen scalar, scored by the observing application's own code."""
    payload = {
        "astro_url": ASTRO_URL,
        "half_window_days": HALF_WINDOW_DAYS,
        "step_minutes": STEP_MINUTES,
        "events": [{"id": e["id"], "t_ms": e["t_ms"], "lat": e["lat"], "lon": e["lon"]}
                   for e in events],
    }
    tmp_in = HERE / "engine" / "out" / "_k092_d1_bridge_in.json"
    tmp_out = HERE / "engine" / "out" / "_k092_d1_bridge_out.json"
    tmp_in.parent.mkdir(parents=True, exist_ok=True)
    tmp_in.write_text(json.dumps(payload), encoding="utf-8")
    r = subprocess.run(["node", str(BRIDGE), str(tmp_in), str(tmp_out)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("node bridge failed (rc=%d):\n%s" % (r.returncode, r.stderr))
    raw = json.loads(tmp_out.read_text(encoding="utf-8"))
    if not raw["selftest"]["exact_match"]:
        raise SystemExit("bridge provenance self-test did not match exactly")

    per, drops = [], []
    for e, g in zip(events, raw["events"]):
        assert e["id"] == g["id"]
        gr = g["grid"]
        t = (gr["t_ms0"] + gr["step_ms"]
             * np.arange(gr["n"], dtype=np.float64)) / 86400000.0
        x = np.asarray(gr["level"], dtype=np.float64)
        rec = local_cycle_readout(t, x, e["t_ms"] / 86400000.0)
        if rec is None:
            drops.append(e["id"])
            continue
        rec.update({"id": e["id"], "time_utc": e["time_utc"], "mag": e["mag"],
                    "lat": e["lat"], "lon": e["lon"], "place": e["place"],
                    "level_cm_at_event": g["level_cm"],
                    "rate_cm_per_h_at_event": g["rate_cm_per_h"]})
        per.append(rec)
    return raw["selftest"], per, drops


def sitetide_arm(events):
    """ROBUSTNESS: the identical downstream analysis on engine/sitetide.py's scalar."""
    step_days = STEP_MINUTES / 1440.0
    n_half = int(round(HALF_WINDOW_DAYS / step_days))
    rate_k = int(round(RATE_HALF_MINUTES / STEP_MINUTES))
    t0 = _dt.datetime(1900, 1, 1)
    epoch = _dt.datetime(1970, 1, 1, tzinfo=_dt.timezone.utc)
    per, drops = [], []
    for e in events:
        ev_dt = epoch + _dt.timedelta(milliseconds=e["t_ms"])
        d_ev = (ev_dt.replace(tzinfo=None) - t0).total_seconds() / 86400.0
        t = d_ev + step_days * np.arange(-n_half, n_half + 1, dtype=np.float64)
        x = ST.site_scalar_at(t0, t, e["lat"], e["lon"], 0.0)   # areal_strain
        rec = local_cycle_readout(t, x, d_ev)
        if rec is None:
            drops.append(e["id"])
            continue
        c = n_half
        rate = (x[c + rate_k] - x[c - rate_k]) / (2.0 * RATE_HALF_MINUTES / 60.0)
        rec.update({"id": e["id"], "time_utc": e["time_utc"], "mag": e["mag"],
                    "rate_per_h_at_event": float(rate)})
        per.append(rec)
    return per, drops


# -------------------------------------------------------------------------- main --
def main():
    events = load_seed_events()
    print("D-1: %d seed-superset events, sha256 verified against the freeze."
          % len(events), flush=True)

    print("app arm (node bridge, astro.ts unmodified) ...", flush=True)
    selftest, app_per, app_drops = app_arm(events)
    print("  %d scored, %d dropped" % (len(app_per), len(app_drops)), flush=True)

    print("sitetide arm (areal_strain, identical downstream analysis) ...", flush=True)
    st_per, st_drops = sitetide_arm(events)
    print("  %d scored, %d dropped" % (len(st_per), len(st_drops)), flush=True)

    # THE DECLARED READOUT -- imported, not restated. P7-23(C) wants the percentile
    # within the LOCAL cycle range, so the per-event u (already on [-1, +1]) is handed
    # in with lo/hi pinned to -1/+1 rather than letting seed_readout renormalise
    # against the set's own extremes.
    app_read = AA.seed_readout([r["normalised_level_u"] for r in app_per],
                               [r["phase_rad_sin"] for r in app_per], lo=-1.0, hi=1.0)
    st_read = AA.seed_readout([r["normalised_level_u"] for r in st_per],
                              [r["phase_rad_sin"] for r in st_per], lo=-1.0, hi=1.0)

    app_quad = quadrant_block([r["level_cm_at_event"] for r in app_per],
                              [r["rate_cm_per_h_at_event"] for r in app_per],
                              [r["phase_rad_sin"] for r in app_per],
                              "app astro.ts solidTideDisplacementCm (FROZEN, PRIMARY)")
    st_quad = quadrant_block([r["level"] for r in st_per],
                             [r["rate_per_h_at_event"] for r in st_per],
                             [r["phase_rad_sin"] for r in st_per],
                             "engine/sitetide.py areal_strain (ROBUSTNESS)")

    ai = {r["id"]: r for r in app_per}
    si = {r["id"]: r for r in st_per}
    both = sorted(set(ai) & set(si))

    def band(u):
        return "TROUGH" if u < -0.5 else ("CREST" if u > 0.5 else "MID_SLOPE")

    agree_band = sum(1 for i in both
                     if band(ai[i]["normalised_level_u"])
                     == band(si[i]["normalised_level_u"]))
    agree_quad = sum(
        1 for i in both
        if ((ai[i]["level_cm_at_event"] < 0.0 and ai[i]["rate_cm_per_h_at_event"] < 0.0)
            == (si[i]["level"] < 0.0 and si[i]["rate_per_h_at_event"] < 0.0)))

    out = {
        "arm": "D-1 (K092_FREEZE.md pre-scoring gate 2; P7-23(C) declared readout)",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 0,
        "scope": {
            "input": "K092_seed_exclusion_superset.csv ONLY",
            "sha256_verified": SEED_CSV_SHA256,
            "n_events": len(events),
            "licence": ("these are the ALREADY-SEEN events, excluded from D-12 by "
                        "construction (P7-23(A)); no event outside this file was "
                        "touched, downloaded or evaluated"),
            "explore_count": "NOT LOGGED: the seed set is not a scoring window.",
        },
        "not_evidence": (
            "P7-22 Ratification 1: a seed-consistent-with-null result EXPLAINS the "
            "motivating observation and collapses the entry's priority, and that is "
            "the result. The reverse is NOT symmetric -- a mid-slope concentration on "
            "a set selected by having been looked at is not a positive finding; only "
            "D-12 (held-out) and D-13 (prospective) can produce evidence."),
        "declared_rule_source": "engine/audit_arcsine.py (committed before this script)",
        "declared_rule": AA.classify_level_set(np.zeros(0))["rule"],
        "density_note": AA.ARCSINE_DENSITY_NOTE,
        "convention_note": CONVENTION_NOTE,
        "scorers": {
            "primary": "app astro.ts (frozen scalar) via exp_k092_d1_bridge.mjs",
            "robustness": "engine/sitetide.py areal_strain",
            "provenance_selftest": selftest,
            "recommendation_followed": "K092_SCALAR_PROVENANCE.md section 5(3)",
        },
        "analysis": {
            "half_window_days": HALF_WINDOW_DAYS, "step_minutes": STEP_MINUTES,
            "rate_half_minutes": RATE_HALF_MINUTES,
            "local_cycle": "interval between the two refined local maxima bracketing "
                           "the event; percentile taken within that cycle's own range",
        },
        "app_arm": {
            "n_scored": len(app_per), "dropped_ids": app_drops,
            "classification": app_read["classification"],
            "quadrant": app_quad,
            "per_event": app_per,
        },
        "sitetide_arm": {
            "n_scored": len(st_per), "dropped_ids": st_drops,
            "classification": st_read["classification"],
            "quadrant": st_quad,
            "per_event": st_per,
        },
        "cross_scorer_agreement": {
            "n_common": len(both),
            "band_agreement": agree_band / max(len(both), 1),
            "quadrant_agreement": agree_quad / max(len(both), 1),
            "expected": "K092_SCALAR_PROVENANCE.md section 4 measured 99.4% quadrant "
                        "agreement on a time-uniform grid, all disagreement "
                        "boundary-confined.",
        },
        "pure_simulation_controls": AA.run(seed=0),
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    c = app_read["classification"]
    print("\n" + "=" * 72)
    print("D-1 VERDICT (primary scorer = the app's own frozen scalar)")
    print("=" * 72)
    print("  n scored           : %d of %d" % (len(app_per), len(events)))
    for b in AA.BAND_NAMES:
        print("  %-10s       : %3d  (%.4f)  z = %+.2f  [null 1/3]"
              % (b, c["counts"][b], c["fractions"][b], c["z_vs_null"][b]))
    print("  VERDICT            : %s" % c["verdict"])
    q = app_quad["sign_condition_quadrant"]
    print("  quadrant (signs)   : %d/%d = %.4f  [null 0.25]  binom p(>=k) = %.4g"
          % (q["count"], app_quad["n"], q["fraction"], q["binomial_p_one_sided_ge"]))
    print("  quadrant (angles)  : %.4f"
          % app_quad["angle_interval_quadrant"]["fraction"])
    print("  sitetide verdict   : %s" % st_read["classification"]["verdict"])
    print("  cross-scorer bands : %.4f agree ; quadrant %.4f agree"
          % (out["cross_scorer_agreement"]["band_agreement"],
             out["cross_scorer_agreement"]["quadrant_agreement"]))
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
