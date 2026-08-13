"""D-1 -- THE ARCSINE / LEVEL-VS-PHASE CONTROL. The cheapest thing that can end K-092.

HYPOTHESIS_LEDGER.md §P7-22, Ratification 1:

  > **D-1 IS THE MANDATORY FIRST ARM, and it may end the entry.** Verified exact:
  > uniform phase puts **exactly one third** of events in the lowest quarter of the
  > tidal range. **A level plot of a perfectly uniform-phase null looks like
  > below-neutral clustering to any eye, including a trained one.** If the seed's
  > observation is quantitatively consistent with 1/3, **the motivating observation is
  > explained and the entry's priority collapses accordingly** -- and that must be
  > reported as the result, not as a preliminary.

THE ARITHMETIC, DERIVED RATHER THAN SIMULATED
----------------------------------------------
Render a sinusoidal tide as a LEVEL, `L = cos(theta)`, and ask where events fall
within the RANGE `[-1, +1]`. For theta uniform on [0, 2 pi),

    P(L < x) = 1 - arccos(x) / pi                       (the arcsine law)

The lowest quarter of the range is `L < -1 + 2*(1/4) = -1/2`, and

    P(L < -1/2) = 1 - arccos(-1/2)/pi = 1 - (2 pi / 3)/pi = **1/3, exactly.**

So a uniform-phase null -- no tidal triggering whatsoever -- puts a third of its
events in the bottom quarter of the tidal range, and two thirds below the midpoint
would be a 50/50 split only in PHASE. The excess is entirely the `arccos` of the
rendering, i.e. the dwell time of a sinusoid near its turning points. `LOWEST_FRACTION`
gives the general answer for any quantile of the range.

WHY THIS IS ALSO THE D-4 CONTROL HARNESS
-----------------------------------------
The build brief attaches a second job: this module is the **level-vs-phase comparison**
that lets D-1 run later as a DECLARED ARM rather than as a one-off. So
`level_vs_phase_report` reports the two readings of the SAME synthetic catalogue side
by side --

  * the LEVEL reading: fraction of events in the lowest quarter of range, against the
    exact 1/3 the null predicts;
  * the PHASE reading: the same events' phases, against uniform, by the event-path
    statistic D-4 builds (`circstat_event.event_kuiper_watson`) --

and the point of putting them next to each other is that the first LOOKS like a strong
effect and the second is flat, on identical data. `real_waveform_control` does it again
with the actual `sitetide` waveform instead of a pure cosine, because the seed was read
off a real tidal display and a pure cosine could be answered with "the real waveform is
not a sinusoid".

Pure simulation. No catalogue, no region, no window, no holdout, no EXPLORE_COUNT line.
Nothing here is evidence.
"""

from __future__ import annotations

import datetime as _dt
import math

import numpy as np

from . import circstat_event as CE
from . import sitetide as ST

# The exact answer §P7-22 verified (simulated 0.33332 at 1e7 draws, ledger).
LOWEST_QUARTER_FRACTION_EXACT = 1.0 / 3.0

ARCSINE_RULE = (
    "§P7-22 Ratification 1 / §K92-0(7): for phase uniform on [0, 2pi) and a level "
    "rendering L = cos(phase), P(L < -1/2) = 1 - arccos(-1/2)/pi = 1/3 EXACTLY. A "
    "third of a perfectly null catalogue lands in the lowest quarter of the tidal "
    "range, and it looks like below-neutral clustering to any eye.")


def lowest_fraction(q: float) -> float:
    """Exact null probability that a uniform-phase event lands in the lowest `q` of range.

    `1 - arccos(2q - 1)/pi`. At q = 1/4 this is 1/3; at q = 1/2 it is 1/2 (the
    midpoint split IS unbiased, which is why the seed's quartile reading and a median
    reading disagree).
    """
    q = float(q)
    if not (0.0 <= q <= 1.0):
        raise ValueError("q must lie in [0, 1]")
    return 1.0 - math.acos(2.0 * q - 1.0) / math.pi


LOWEST_FRACTION = lowest_fraction        # public alias used in the docstring above


def level_from_phase(theta):
    """The level rendering the control is about: L = cos(theta), range [-1, +1]."""
    return np.cos(np.asarray(theta, dtype=np.float64))


def fraction_in_lowest_quarter(level, q=0.25):
    """Fraction of `level` below the lowest `q` of its OWN observed range.

    Deliberately referenced to the observed min/max rather than to the theoretical
    +-1: that is what a person reading a tidal display does, and it is the only
    version of the statistic that could be computed from a plot.
    """
    x = np.asarray(level, dtype=np.float64)
    lo, hi = float(x.min()), float(x.max())
    return float(np.mean(x < lo + float(q) * (hi - lo)))


# ----------------------------------------------------------------- the control --
def arcsine_control(n=200000, q=0.25, seed=0):
    """Uniform phases -> level rendering -> the 1/3 that has nothing to do with tides."""
    rng = np.random.default_rng(int(seed))
    theta = rng.random(int(n)) * 2.0 * np.pi
    lev = level_from_phase(theta)
    frac = fraction_in_lowest_quarter(lev, q)
    exact = lowest_fraction(q)
    se = math.sqrt(exact * (1.0 - exact) / int(n))
    return {
        "n": int(n), "q": float(q),
        "measured_fraction": frac,
        "exact_null_fraction": exact,
        "binomial_se": se,
        "z": (frac - exact) / se if se > 0 else float("nan"),
        "rule": ARCSINE_RULE,
    }


def level_vs_phase_report(n=20000, q=0.25, seed=0, n_reference=200000):
    """THE DECLARED-ARM OBJECT: the same events read as LEVEL and read as PHASE.

    The level reading reproduces the 1/3 artifact; the phase reading is flat under the
    event-path Kuiper/Watson statistic D-4 builds. Both on identical data, which is
    the whole demonstration: the artifact is in the RENDERING, not in the catalogue.
    """
    rng = np.random.default_rng(int(seed))
    theta = rng.random(int(n)) * 2.0 * np.pi
    ref = rng.random(int(n_reference)) * 2.0 * np.pi
    lev = level_from_phase(theta)
    kw = CE.event_kuiper_watson(theta, ref)
    mom = CE.event_moments(theta)
    return {
        "arm": "D-1 arcsine / level-vs-phase control",
        "n_events": int(n),
        "level_reading": {
            "statistic": "fraction of events in the lowest quarter of tidal range",
            "measured": fraction_in_lowest_quarter(lev, q),
            "exact_uniform_phase_null": lowest_fraction(q),
            "reads_as": ("STRONG below-neutral clustering to an eye, and it is "
                         "EXACTLY what a null produces"),
        },
        "phase_reading": {
            "statistic": "event-path Kuiper V* / Watson U2* vs a uniform pool",
            "V_star": kw["V_star"], "U2_star": kw["U2_star"],
            "R1": mom["R1"], "R2": mom["R2"],
            "reads_as": "FLAT -- the same events carry no phase concentration",
        },
        "verdict_rule": (
            "§P7-22 Ratification 1: if the seeding observation is quantitatively "
            "consistent with the exact null fraction here, THE MOTIVATING "
            "OBSERVATION IS EXPLAINED and that is the result, not a preliminary."),
        "rule": ARCSINE_RULE,
    }


def real_waveform_control(lat_deg=55.34, lon_deg=-160.50, depth_km=10.0,
                          t0=None, n_days=365.0, sample_minutes=10.0, q=0.25):
    """The same control with the ACTUAL `sitetide` waveform instead of a pure cosine.

    A pure-cosine demonstration invites "but the real tide is multi-constituent". It
    is: the beat between M2/S2/K1/O1 broadens the level histogram and moves the number
    off the exact 1/3 -- and the direction and size of that move is a MEASUREMENT this
    module makes rather than a hand-wave.

    The default site is the Sand Point coordinates the seed came from. **This is not a
    look at Alaska in the §P7-22 D-7 sense and cannot become one**: no catalogue, no
    event times, no phase statistic on data -- it evaluates a deterministic
    astronomical waveform on a uniform time grid. Change `lat_deg`/`lon_deg` and the
    conclusion does not move, which `test_audit_arcsine.py` checks at a second site.
    """
    t0 = _dt.datetime(2020, 1, 1) if t0 is None else t0
    dt_days = float(sample_minutes) / (60.0 * 24.0)
    n = int(round(float(n_days) / dt_days))
    t = dt_days * np.arange(n, dtype=np.float64)
    lev = ST.site_scalar_at(t0, t, lat_deg, lon_deg, depth_km)
    return {
        "arm": "D-1 control on the real multi-constituent sitetide waveform",
        "site": {"lat": float(lat_deg), "lon": float(lon_deg),
                 "depth_km": float(depth_km)},
        "note": ("uniform-in-TIME sampling of a deterministic astronomical waveform. "
                 "No catalogue, no event times, no phase statistic on data -- this is "
                 "not a D-7 look."),
        "n_samples": int(n),
        "q": float(q),
        "measured_fraction_lowest_quarter": fraction_in_lowest_quarter(lev, q),
        "pure_cosine_exact": lowest_fraction(q),
        "scalar": ST.SCALAR_FOR_PHASE,
        "scope_flags": ST.SCOPE_FLAGS,
    }


# ------------------------- §P7-23(C): THE DECLARED TROUGH-vs-MID-SLOPE READOUT --
# §P7-23(C) landed after this module's first version and amends D-1:
#
#   > **RULED. D-1's declared readout must, for each seed event, report BOTH the level
#   > percentile within the local tidal range AND the phase angle, and classify the
#   > seed set as TROUGH-CONCENTRATED (artifact-consistent) or MID-SLOPE-CONCENTRATED
#   > (artifact-inconsistent), against a classification rule declared in advance --
#   > BEFORE any Alaska event outside the seed set is scored.**
#
# The physics behind the two-sided reading, verbatim from the ruling: for x = A sin θ
# under uniform phase the level density is arcsine, `f(x) = 1/(pi sqrt(1 - x^2))` --
# **0.318 at mid-level, 1.019 at |x| = 0.95**. Uniform phase piles events at the
# EXTREMES and puts the FEWEST at mid-slope. So trough concentration is what the
# artifact manufactures, and mid-slope concentration is evidence AGAINST the artifact
# explanation -- and mid-fall is also where rate-state sensitivity peaks, so artifact
# and physics separate cleanly and in opposite directions.
#
# THE CLASSIFICATION RULE, DECLARED HERE AND IN ADVANCE (S-9, one way). On the
# normalised level u = 2 (x - min) / (max - min) - 1 in [-1, +1] within the event's
# own local cycle, three bands:
#
#     TROUGH     u <  -1/2
#     MID-SLOPE  -1/2 <= u <= +1/2
#     CREST      u >  +1/2
#
# Under uniform phase these are **exactly equiprobable at 1/3 each**: P(u < -1/2) =
# 1 - arccos(-1/2)/pi = 1/3, P(u > +1/2) = 1/3 by symmetry, so the middle band is 1/3
# too. That equiprobability is why these bands and not others -- the null is a flat
# multinomial and needs no simulation to state.
#
# VERDICT RULE: the set is TROUGH-CONCENTRATED if the trough band holds the largest
# share AND that share exceeds 1/3 by more than `Z_CLASSIFY` binomial standard errors;
# MID-SLOPE-CONCENTRATED under the mirror condition on the middle band; UNCLASSIFIED
# otherwise. `Z_CLASSIFY = 2.0` is declared here, before any data.
BAND_NAMES = ("TROUGH", "MID_SLOPE", "CREST")
BAND_NULL_PROBABILITY = 1.0 / 3.0
Z_CLASSIFY = 2.0

ARCSINE_DENSITY_NOTE = (
    "§P7-23(C): for x = A sin(theta) under uniform phase the level density is "
    "arcsine, f(x) = 1/(pi sqrt(1 - x^2)) -- 0.318 at mid-level and 1.019 at "
    "|x| = 0.95. Uniform phase piles events at the EXTREMES and puts the FEWEST at "
    "mid-slope, so TROUGH concentration is artifact-CONSISTENT and MID-SLOPE "
    "concentration is artifact-INCONSISTENT.")

# 'Below neutral AND falling' is the intersection of theta in (pi, 2pi) with
# theta in (pi/2, 3pi/2), i.e. theta in (pi, 3pi/2): exactly one quadrant.
QUADRANT_NULL_FRACTION = 0.25
QUADRANT_NOTE = (
    "§P7-23(D): 'below neutral and falling' names theta in (pi, 3pi/2) for "
    "x = A sin(theta) -- exactly a quarter cycle, so its uniform-phase null "
    "probability is EXACTLY 1/4 and the specification is rotation-free ONCE THE "
    "SCALAR IS FIXED. Which scalar is a provenance question, not a convention "
    "question (engine/conventions_d.py::assert_scalar_provenance).")


def arcsine_density(x):
    """The level density f(x) = 1/(pi sqrt(1 - x^2)) for x = A sin(theta), |x| < 1."""
    u = np.asarray(x, dtype=np.float64)
    return 1.0 / (np.pi * np.sqrt(np.maximum(1.0 - u * u, 1e-300)))


def normalised_level(level, lo=None, hi=None):
    """u in [-1, +1] within the local tidal range: u = 2(x - lo)/(hi - lo) - 1."""
    x = np.asarray(level, dtype=np.float64)
    a = float(x.min()) if lo is None else float(lo)
    b = float(x.max()) if hi is None else float(hi)
    return 2.0 * (x - a) / max(b - a, 1e-300) - 1.0


def level_percentile(level, lo=None, hi=None):
    """The level's percentile within its local range, in [0, 1]. The (C) readout."""
    return 0.5 * (normalised_level(level, lo, hi) + 1.0)


def classify_level_set(u):
    """§P7-23(C)'s declared classification. `u` is the normalised level in [-1, +1]."""
    x = np.asarray(u, dtype=np.float64)
    n = int(x.size)
    counts = {"TROUGH": int(np.sum(x < -0.5)),
              "MID_SLOPE": int(np.sum((x >= -0.5) & (x <= 0.5))),
              "CREST": int(np.sum(x > 0.5))}
    p = BAND_NULL_PROBABILITY
    se = math.sqrt(p * (1.0 - p) / max(n, 1))
    frac = {k: v / max(n, 1) for k, v in counts.items()}
    z = {k: (frac[k] - p) / se if se > 0 else float("nan") for k in counts}
    top = max(frac, key=frac.get)
    if top == "TROUGH" and z["TROUGH"] > Z_CLASSIFY:
        verdict = "TROUGH-CONCENTRATED (artifact-CONSISTENT)"
    elif top == "MID_SLOPE" and z["MID_SLOPE"] > Z_CLASSIFY:
        verdict = "MID-SLOPE-CONCENTRATED (artifact-INCONSISTENT)"
    else:
        verdict = "UNCLASSIFIED"
    return {
        "n": n, "counts": counts, "fractions": frac, "z_vs_null": z,
        "null_probability_each_band": p,
        "z_threshold_declared": Z_CLASSIFY,
        "verdict": verdict,
        "rule": ("§P7-23(C), declared in advance: on u = 2(x-min)/(max-min) - 1, "
                 "the bands u < -1/2 / |u| <= 1/2 / u > +1/2 are EXACTLY "
                 "equiprobable at 1/3 under uniform phase. TROUGH-CONCENTRATED if "
                 "the trough band leads and exceeds 1/3 by > %.1f binomial SE; "
                 "MID-SLOPE-CONCENTRATED under the mirror condition; else "
                 "UNCLASSIFIED." % Z_CLASSIFY),
        "density_note": ARCSINE_DENSITY_NOTE,
    }


def seed_readout(levels, phases, lo=None, hi=None):
    """§P7-23(C)'s per-event readout: BOTH the level percentile AND the phase angle.

    Returns the per-event pairs and the declared classification of the set. It takes
    levels and phases as ARRAYS and reads nothing from disk -- the caller supplies the
    seed set, and this module never learns which region it came from.
    """
    u = normalised_level(levels, lo, hi)
    pct = 0.5 * (u + 1.0)
    ph = np.mod(np.asarray(phases, dtype=np.float64), 2.0 * np.pi)
    return {
        "per_event": [{"level_percentile": float(a), "normalised_level": float(b),
                       "phase_rad": float(c), "phase_deg": float(np.degrees(c))}
                      for a, b, c in zip(pct, u, ph)],
        "classification": classify_level_set(u),
        "quadrant_below_neutral_and_falling": {
            "fraction": float(np.mean((ph > np.pi) & (ph < 1.5 * np.pi))),
            "null_fraction": QUADRANT_NULL_FRACTION,
            "note": QUADRANT_NOTE,
        },
        "timing_rule": ("§P7-23(C): this classification is a FACT ABOUT THE SEED, "
                        "priced 0, and must be produced BEFORE any event outside the "
                        "seed set is scored."),
    }


def run(seed=0):
    """Everything D-1 owes, as one JSON-able record."""
    return {
        "control": arcsine_control(seed=seed),
        "level_vs_phase": level_vs_phase_report(seed=seed),
        "real_waveform": real_waveform_control(),
        "declared_classification_rule": classify_level_set(np.zeros(0))["rule"],
        "density_note": ARCSINE_DENSITY_NOTE,
        "quadrant_note": QUADRANT_NOTE,
        "priced_tests": 0,
        "explore_count": ("NOT LOGGED: pure simulation, no catalogue, no window, "
                          "no holdout -- nothing to spend."),
    }


if __name__ == "__main__":          # pragma: no cover - operator convenience
    import json
    print(json.dumps(run(), indent=2, default=float))
