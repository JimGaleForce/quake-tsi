"""B1 / §S2 -- THE PER-EVENT PROPERTY JOIN. The object the program did not have.

The miner works on DAILY BINNED COUNTS with covariate columns. The SEARCHER needs
**per-event property values**, evaluated at the event's own timestamp and coordinates,
never at a bin centre -- SEARCHER.md §S2.2(1), and the K-092 lesson behind it: a daily
bin has exactly one hour phase, which is why `observer.obs_utc_hour_phase` is zero on
the count path BY CONSTRUCTION. The SEARCHER lives entirely on the sub-daily path and
inherits `observer.assert_subdaily_gate` as a hard entry gate.

THE FOUR RULES THIS MODULE ENFORCES IN CODE RATHER THAN IN PROSE
-----------------------------------------------------------------
1. **`property_class` is a REQUIRED field and a column without one cannot be built.**
   §P7-24 SP-2 organises the null-validity layer BY PROPERTY CLASS, built once and
   reused forever, and §S11(4) says B1 must carry the class as a required column and
   refuse to emit a property without one. `PropertyClassMissing` is that refusal.
   The classes and their mandatory null layers are SP-2's table, transcribed into
   `PROPERTY_CLASSES` and checked by `may_promote`.

2. **A property whose class has no null layer ATTACHED may appear in the ranked list
   and may NOT promote.** SP-2's last sentence, in code: `may_promote()` returns
   `(False, reason)` and the searcher carries the reason onto the row. This is the
   difference between a checkable precondition and a review-time judgement.

3. **Level properties carry their DWELL-TIME density with them, and a level
   concentration statistic that does not divide by it is not emitted at all**
   (§S2.2(2)). §P7-23(C) verified the arithmetic exactly: for `x = A sin(theta)` under
   uniform phase the level density is arcsine, `f(x) = 1 / (pi sqrt(1 - x^2))` -- 0.318
   at mid-level, 1.019 at |x| = 0.95, and exactly **1/3 of all events in the lowest
   quarter of the range**. `dwell_pit` is the transform; `arcsine_cdf` is the analytic
   case the transform is checked against. A level column built without a dwell measure
   raises `DwellDensityMissing` the moment a concentration statistic asks for its PIT
   values. **This is the single most likely way the SEARCHER would manufacture Jim's
   own Sand Point observation out of nothing, and it is disarmed at the data layer.**

4. **Provenance is a field, not a memory** (§S2.2(3)). Every column carries
   `{source, version, url_or_module, convention, scalar_definition}`. §P7-23(D): the
   quadrant is convention-free only ONCE THE SCALAR IS FIXED, and a property without
   its convention attached is an S-18 clause 1 defect waiting to happen.

THE SITE RULE, declared here because it is a construction choice and S-9 permits one
--------------------------------------------------------------------------------------
Site-dependent properties (the solid tide, local solar hour) are evaluated at the
**REGION'S DECLARED CENTROID**, one site per region, not at each event's own
coordinates. Two reasons, both operative:

  * it makes the property a function of TIME ALONE within a cell, which is what lets
    the ETAS null traverse the *identical* code path with nothing but times -- §P7-22(a)'s
    common-mode requirement taken literally (`sitetide.tidal_maxima`'s own docstring
    makes the same argument for the same reason); and
  * it is the construction K-092 already froze: a region box and one scalar.

The alternative (per-event coordinates) would require the null to carry a location
model as well as a time model, and a location model fitted to the observed events is
exactly the leak the ETAS null exists to prevent. **One value, declared, no
alternatives run (S-9).**

Nothing in this module is evidence. It builds coordinates, not conclusions.
"""

from __future__ import annotations

import datetime as _dt

import numpy as np

from . import (circstat as C, ephemeris as E, marks_ext as MX, sitetide as ST,
               splits)

PROPERTY_RULE_ID = "S2-properties-v1"

# --------------------------------------------------------- SP-2's class table --
# Transcribed from HYPOTHESIS_LEDGER.md §P7-24 SP-2. `layer` is the MANDATORY null
# layer that must EXIST before any cell of that class may promote; `attach` names the
# artifact the searcher must hand the column for the layer to count as built.
PROPERTY_CLASSES = {
    "human-schedule": {
        "examples": "day-of-week, hour-of-day, month, calendar date",
        "layer": ("F7 observer controls AT THE MAGNITUDE IN QUESTION -- weekday/hour "
                  "completeness and reporting-schedule bias must be MEASURED, not "
                  "assumed. Plus anthropogenic seismicity screens (blasting, "
                  "reservoir, injection), which genuinely have weekly cycles."),
        "attach": "observer_features",
    },
    "level-waveform-phase": {
        "examples": "tidal level, tidal phase, lunar distance, any sinusoidal scalar",
        "layer": ("dwell-time-corrected null -- the arcsine density piles "
                  "uniform-phase events at the extremes (§P7-23(C): 0.318 at "
                  "mid-level vs 1.019 at |x| = 0.95). A level-threshold statistic "
                  "against a uniform-phase null is wrong by construction."),
        "attach": "dwell",
    },
    "clustering-derived": {
        "examples": "inter-event time, distance-to-prior, cluster membership",
        "layer": ("FULL ETAS WITH TRIGGERING, never background-only, never a phase "
                  "permutation (§P7-23(A.3)) -- aftershocks inherit their mainshock's "
                  "property neighbourhood and a background-only null scores that "
                  "inheritance as signal."),
        "attach": "mark_null",
    },
    "count-path-periodic": {
        "examples": "any daily-binned cyclic feature",
        "layer": ("MEASURED VIF for that statistic and that null (§P7-12(a)); the "
                  "count-path 24.08 does not transfer (§P7-22 Q1)."),
        "attach": "measured_vif",
    },
    "catalogue-endogenous": {
        "examples": "b-value, depth statistics, magnitude marks",
        "layer": ("composition-preserving null, and the §P7-20(1) same-quantity "
                  "exclusion -- never pair a feature with a mark derived from the "
                  "same quantity."),
        "attach": "mark_null",
    },
}

SP2_RULE = (
    "§P7-24 SP-2: every property belongs to a declared PROPERTY CLASS, and each class "
    "carries a mandatory null-validity layer that must exist before any cell of that "
    "class may promote. A property whose class has no null layer built CAN appear in "
    "the ranked list and CANNOT be promoted. The layer is built once per class and "
    "reused by every property in it forever -- which is what makes 'promote quickly' "
    "affordable.")

DWELL_RULE = (
    "SEARCHER.md §S2.2(2) / §P7-23(C): a level property carries its DWELL-TIME "
    "occupancy measure with it, and a level-concentration statistic that does not "
    "divide by it is NOT EMITTED. For x = A sin(theta) under uniform phase the level "
    "density is arcsine, f(x) = 1/(pi sqrt(1-x^2)): 0.318 at mid-level, 1.019 at "
    "|x| = 0.95, and exactly 1/3 of all events in the lowest quarter of range. A "
    "concentration read off an uncorrected level is that arithmetic, not the Earth.")

SITE_RULE = (
    "S2-properties-v1 site rule: site-dependent properties are evaluated at the "
    "REGION'S DECLARED CENTROID, one site per region -- so the property is a function "
    "of TIME ALONE within a cell and the ETAS null can traverse the IDENTICAL code "
    "path with nothing but times (§P7-22(a) common-mode). One declared value, no "
    "alternatives run (S-9).")

SUBDAILY_RULE = (
    "SEARCHER.md §S2.2(1): every property is evaluated at the event's own timestamp, "
    "never at a daily bin centre. A daily bin has exactly one hour phase, which is "
    "why obs_utc_hour_phase is zero on the count path BY CONSTRUCTION. Columns whose "
    "`subdaily` flag is False have NOT escaped the day-binning sinc and say so per "
    "column (§P7-3(3)).")

PROPERTY_TYPES = ("phase", "categorical", "level")


# ------------------------------------------------------------------ refusals --
class PropertyClassMissing(AssertionError):
    """A property column was built without a declared SP-2 property class."""


class UnknownPropertyClass(AssertionError):
    """A property declared a class that is not in SP-2's table."""


class DwellDensityMissing(AssertionError):
    """A level concentration statistic was asked for without a dwell measure."""


class ProvenanceMissing(AssertionError):
    """A property column was built without its §S2.2(3) provenance fields."""


PROVENANCE_FIELDS = ("source", "version", "url_or_module", "convention",
                     "scalar_definition")


# ------------------------------------------------------ the dwell-time layer --
def arcsine_cdf(x, amplitude=1.0):
    """CDF of `A sin(theta)` for theta ~ U(0, 2pi): the §P7-23(C) analytic case.

    F(x) = 1/2 + arcsin(x/A)/pi. Present so the EMPIRICAL dwell transform below can be
    checked against a case whose answer is known in closed form rather than trusted.
    """
    a = float(amplitude)
    return 0.5 + np.arcsin(np.clip(np.asarray(x, dtype=np.float64) / a, -1.0, 1.0)) / np.pi


def empirical_dwell(values):
    """The occupancy measure of a level property, as a sorted sample of its own values.

    `values` are the property evaluated on a DENSE UNIFORM TIME GRID over the cell's
    window -- i.e. how much TIME the property spends at each level, which is precisely
    the density the arcsine argument is about. Returned as a sorted array; `dwell_pit`
    turns it into the probability-integral transform.
    """
    v = np.asarray(values, dtype=np.float64).ravel()
    v = v[np.isfinite(v)]
    if v.size < 2:
        raise ValueError("a dwell measure needs at least two finite samples")
    return np.sort(v)


def dwell_pit(values, dwell):
    """PIT-transform level values through their own dwell measure -> [0, 1).

    Mid-rank on ties, so a property with plateaus (a clipped index, a quantised Kp)
    does not acquire a spurious pile-up at the plateau's edge. Under the null the
    transformed values are uniform by construction, which is what makes a Kuiper
    statistic on them a DWELL-CORRECTED concentration statistic rather than a
    measurement of the occupancy measure itself.
    """
    if dwell is None:
        raise DwellDensityMissing(DWELL_RULE)
    d = np.asarray(dwell, dtype=np.float64)
    v = np.asarray(values, dtype=np.float64)
    lo = np.searchsorted(d, v, side="left")
    hi = np.searchsorted(d, v, side="right")
    return (0.5 * (lo + hi)) / float(d.size)


def dwell_lowest_quarter_fraction(dwell):
    """Fraction of the occupancy measure in the lowest quarter of the RANGE.

    §P7-23(C)'s headline number: for a pure sinusoid this is 1/3, not 1/4. Reported on
    every level column so the reader sees the size of the correction being applied.
    """
    d = np.asarray(dwell, dtype=np.float64)
    lo, hi = float(d.min()), float(d.max())
    cut = lo + 0.25 * (hi - lo)
    return float(np.count_nonzero(d <= cut)) / float(d.size)


# -------------------------------------------------------------- the column ----
class PropertyColumn:
    """One property, bound to one site, with its class, provenance and null layer.

    `evaluate(day_float)` is the IDENTICAL code path the ETAS null traverses. It takes
    times and nothing else, which is the whole reason for the site rule above: a
    signature that cannot accept a permuted property cannot be handed one.
    """

    __slots__ = ("name", "family", "ptype", "pclass", "evaluate", "provenance",
                 "dwell", "categories", "subdaily", "attached", "notes")

    def __init__(self, name, family, ptype, pclass, evaluate, provenance,
                 dwell=None, categories=None, subdaily=True, attached=None,
                 notes=None):
        if not pclass:
            raise PropertyClassMissing(
                "property %r was built without a declared SP-2 property class. %s"
                % (name, SP2_RULE))
        if pclass not in PROPERTY_CLASSES:
            raise UnknownPropertyClass(
                "property %r declares class %r, which is not one of SP-2's %r. %s"
                % (name, pclass, tuple(PROPERTY_CLASSES), SP2_RULE))
        if ptype not in PROPERTY_TYPES:
            raise ValueError("property type must be one of %r, got %r"
                             % (PROPERTY_TYPES, ptype))
        missing = [f for f in PROVENANCE_FIELDS if f not in (provenance or {})]
        if missing:
            raise ProvenanceMissing(
                "property %r is missing provenance fields %s. §S2.2(3): provenance is "
                "a FIELD, not a memory -- a property without its convention attached "
                "is an S-18 clause 1 defect waiting to happen."
                % (name, ", ".join(missing)))
        if ptype == "level" and dwell is None:
            raise DwellDensityMissing(
                "level property %r was built without a dwell measure. %s"
                % (name, DWELL_RULE))
        self.name = str(name)
        self.family = str(family)
        self.ptype = ptype
        self.pclass = pclass
        self.evaluate = evaluate
        self.provenance = dict(provenance)
        self.dwell = (None if dwell is None else np.asarray(dwell, dtype=np.float64))
        self.categories = (None if categories is None else tuple(categories))
        self.subdaily = bool(subdaily)
        self.attached = set(attached or ())
        self.notes = list(notes or [])
        # AUTO-ATTACHED, and only this one. For class `level-waveform-phase` the
        # mandatory SP-2 layer is the dwell-time correction, and a column of that
        # class carries it BY CONSTRUCTION: a level column holds its occupancy
        # measure (checked above -- it cannot be built without one), and a phase
        # column of that class is compared against the EMPIRICAL reference CDF of the
        # ETAS null pool, which is the same correction estimated rather than
        # tabulated (`circstat_event`'s two-sample construction). The layer is
        # therefore built the moment the column exists, and recording that here is a
        # statement of fact rather than a waiver.
        #
        # `observer_features`, `mark_null`, `measured_vif` are NEVER auto-attached:
        # each is a MEASUREMENT somebody has to make, and SP-2's whole point is that
        # a property whose class layer has not been built cannot promote.
        if pclass == "level-waveform-phase":
            self.attached.add("dwell")

    # -- the SP-2 precondition, as a query rather than a judgement --------------
    def may_promote(self):
        """(bool, reason). SP-2: no null layer built -> ranked list only, never promote."""
        need = PROPERTY_CLASSES[self.pclass]["attach"]
        if need not in self.attached:
            return False, (
                "SP-2: property class %r requires its null layer (%s) to be BUILT and "
                "ATTACHED before any cell of that class may promote; %r is not "
                "attached. The cell may appear in the ranked list and may NOT be "
                "promoted. Mandatory layer: %s"
                % (self.pclass, need, need, PROPERTY_CLASSES[self.pclass]["layer"]))
        return True, "SP-2 null layer %r attached for class %r" % (need, self.pclass)

    def attach(self, what):
        """Record that this class's mandatory null layer has been built for this cell."""
        self.attached.add(str(what))
        return self

    # -- the statistic-facing values -------------------------------------------
    def statistic_values(self, day_float):
        """The values a concentration statistic consumes, per property TYPE.

        phase       -> the angle in [0, 2pi)
        level       -> the DWELL-CORRECTED PIT value in [0, 1), mapped to an angle so
                       the same rotation-invariant Kuiper machinery applies. §S3.1's
                       'Kuiper on the PIT-transformed level', literally.
        categorical -> the integer category index
        """
        v = np.asarray(self.evaluate(day_float), dtype=np.float64)
        if self.ptype == "phase":
            return C.wrap_phase(v)
        if self.ptype == "level":
            if self.dwell is None:                  # unreachable via __init__, kept
                raise DwellDensityMissing(DWELL_RULE)   # pragma: no cover
            return 2.0 * np.pi * dwell_pit(v, self.dwell)
        return v

    def record(self):
        """The column's declaration row -- hash-affecting, provenance attached."""
        ok, why = self.may_promote()
        rec = {
            "property": self.name,
            "family": self.family,
            "type": self.ptype,
            "property_class": self.pclass,
            "class_null_layer": PROPERTY_CLASSES[self.pclass]["layer"],
            "null_layer_attached": sorted(self.attached),
            "may_promote": bool(ok),
            "may_promote_reason": why,
            "subdaily": self.subdaily,
            "categories": (list(self.categories) if self.categories else None),
            "provenance": dict(self.provenance),
            "dwell_time_corrected": bool(self.ptype != "level" or self.dwell is not None),
            "notes": list(self.notes),
        }
        if self.ptype == "level" and self.dwell is not None:
            rec["dwell_n_samples"] = int(self.dwell.size)
            rec["dwell_lowest_quarter_fraction"] = dwell_lowest_quarter_fraction(
                self.dwell)
            rec["dwell_note"] = DWELL_RULE
        if self.ptype == "level":
            rec["arcsine_reference_lowest_quarter"] = 1.0 / 3.0
        return rec


# ------------------------------------------------------------- the families ---
FAMILY_ORDER = ("solid_tide", "ephemeris", "human_clock", "geomagnetic",
                "earth_rotation", "season", "marks", "clocks")

DECLARED_MAG_STRATA = (4.5, 5.0, 5.5, 6.0)
MAG_STRATA_RULE = (
    "SEARCHER.md §S3.2 / §P7-24 SP-6.7: the magnitude strata are ENUMERATED in the "
    "scan declaration {M>=4.5, 5.0, 5.5, 6.0}. A claim at an unenumerated threshold "
    "is a NEW SEED, not a result -- otherwise 'which magnitude counts as a major' is "
    "an undeclared search over thresholds. Per §P7-5(4) the headline effect is the "
    "MINIMUM over the declared Mc set: a worst-case reading, demote-only, priced 0.")

HUMAN_SCHEDULE_CARVE_OUT = (
    "SEARCHER.md §S3.4 (F7-d): the human-schedule arm is scientifically LIVE ONLY at "
    "M >= 6.0 and above, where the catalogue is complete and no M6 goes undetected on "
    "a Sunday; below it the arm is a pure OBSERVER measurement. Three residual "
    "channels survive the carve-out and are named rather than waved past: (i) "
    "magnitude assignment can drift with analyst schedule near the M6 boundary, so an "
    "M >= 6.0 result must be re-run at M >= 6.5 as a demote-only audit; (ii) "
    "origin-time is waveform-derived and schedule-free, so the TIME is clean where "
    "the MAGNITUDE is not; (iii) historical completeness is a declared, hash-affecting "
    "start year set from published completeness studies, not from where our CSV "
    "happens to begin.")

HUMAN_SCHEDULE_MC_FLOOR = 6.0


def _prov(source, version, url_or_module, convention, scalar_definition):
    return {"source": source, "version": version, "url_or_module": url_or_module,
            "convention": convention, "scalar_definition": scalar_definition}


def _dense_grid(day_lo, day_hi, sample_minutes=30.0, cap=400000):
    """A dense uniform TIME grid -- the thing a dwell (occupancy) measure is built on."""
    step = float(sample_minutes) / (60.0 * 24.0)
    n = int(np.ceil((float(day_hi) - float(day_lo)) / step)) + 1
    if n > int(cap):
        step = (float(day_hi) - float(day_lo)) / float(cap)
        n = int(cap) + 1
    return float(day_lo) + step * np.arange(n, dtype=np.float64)


# ---- family A: solid tide (K-092 family) -------------------------------------
def solid_tide_columns(t0: _dt.datetime, day_lo, day_hi, site_lat, site_lon,
                       site_depth_km=0.0, sample_minutes=30.0):
    """Family A. Phase (D-0 Tanaka convention), level, and quadrant, at one site.

    The phase map is `sitetide.tanaka_phase`'s maxima-referenced construction, and the
    MAXIMA ARE PRECOMPUTED ONCE for the site and span so that signal and null share a
    bit-identical map -- `sitetide.tidal_maxima`'s own reason, adopted here.
    """
    pad = 1.0
    t_max = ST.tidal_maxima(t0, float(day_lo) - pad, float(day_hi) + pad,
                            site_lat, site_lon, site_depth_km)

    def phase_fn(day_float):
        return ST.phase_from_maxima(t_max, np.asarray(day_float, dtype=np.float64))

    def level_fn(day_float):
        return ST.site_scalar_at(t0, np.asarray(day_float, dtype=np.float64),
                                 site_lat, site_lon, site_depth_km)

    def quadrant_fn(day_float):
        return np.floor(phase_fn(day_float) / (np.pi / 2.0))

    grid = _dense_grid(day_lo, day_hi, sample_minutes)
    dwell = empirical_dwell(level_fn(grid))

    prov = _prov(
        source="engine/sitetide.py degree-2 Love-number body tide (Wahr/Tanaka)",
        version="D3-sitetide-v1",
        url_or_module="engine.sitetide.site_scalar_at / tanaka_phase",
        convention="D0-tanaka-stressmax-v1 (phase 0 = local scalar MAXIMUM)",
        scalar_definition=ST.SCALAR_FOR_PHASE)
    return [
        PropertyColumn("tide_phase", "solid_tide", "phase", "level-waveform-phase",
                       phase_fn, prov, subdaily=True,
                       notes=[ST.SCOPE_FLAGS, SITE_RULE]),
        PropertyColumn("tide_level", "solid_tide", "level", "level-waveform-phase",
                       level_fn, prov, dwell=dwell, subdaily=True,
                       notes=[DWELL_RULE, SITE_RULE]),
        PropertyColumn("tide_quadrant", "solid_tide", "categorical",
                       "level-waveform-phase", quadrant_fn, prov,
                       categories=(0, 1, 2, 3), subdaily=True,
                       notes=["§P7-23(D): the quadrant is convention-free ONLY once "
                              "the scalar is fixed; the scalar is fixed in "
                              "provenance.scalar_definition.", SITE_RULE]),
    ]


# ---- family B: lunar/solar ephemeris -----------------------------------------
_EPH_PHASES = (
    ("lunar_synodic_phase", "synodic_rad", "moon longitude - sun longitude, 0 = new moon"),
    ("lunar_anomalistic_phase", "anomalistic_rad", "lunar mean anomaly M', 0 = perigee"),
    ("lunar_draconic_phase", "draconic_rad", "argument of latitude F, 0 = ascending node"),
    ("solar_annual_phase", "annual_rad", "solar longitude, 0 = March equinox"),
)


def ephemeris_columns(t0: _dt.datetime, day_lo, day_hi, sample_minutes=30.0):
    """Family B, on the PER-EVENT path `julian_day_at` already supports (§S2.1 B)."""
    cols = []
    prov = _prov(
        source="engine/ephemeris.py, Meeus low-precision series",
        version="ephemeris-v1",
        url_or_module="engine.ephemeris.ephemeris_table_at",
        convention="cycle phases in radians in [0, 2pi), definitions in "
                   "ephemeris.ephemeris_table's docstring",
        scalar_definition="geocentric solar/lunar positions; ~0.2 deg lunar longitude")
    for name, key, meaning in _EPH_PHASES:
        def fn(day_float, _k=key):
            return E.ephemeris_table_at(t0, np.asarray(day_float, dtype=np.float64))[_k]
        p = dict(prov)
        p["scalar_definition"] = meaning
        cols.append(PropertyColumn(name, "ephemeris", "phase", "level-waveform-phase",
                                   fn, p, subdaily=True))

    grid = _dense_grid(day_lo, day_hi, sample_minutes)
    for name, key, meaning in (
            ("lunar_distance", "moon_dist_km", "geocentric lunar distance, km"),
            ("lunar_declination", "moon_dec_deg", "lunar declination, degrees"),
            ("sun_moon_elongation", "elongation_deg",
             "Sun-Moon elongation, degrees; 0 = syzygy (new), 180 = syzygy (full)")):
        def fn(day_float, _k=key):
            return E.ephemeris_table_at(t0, np.asarray(day_float, dtype=np.float64))[_k]
        p = dict(prov)
        p["scalar_definition"] = meaning
        cols.append(PropertyColumn(name, "ephemeris", "level", "level-waveform-phase",
                                   fn, p, dwell=empirical_dwell(fn(grid)),
                                   subdaily=True, notes=[DWELL_RULE]))
    return cols


# ---- family C: human clock ----------------------------------------------------
DOW_NAMES = ("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday",
             "Sunday")


def human_clock_columns(t0: _dt.datetime, site_lon):
    """Family C, PROMOTED from control feature to scannable property under F7 controls.

    Every column here is class `human-schedule`, so under SP-2 none of them can
    promote until `observer.observer_features` has been computed for the region and
    stratum and attached; and under F7-d none of them can promote below M >= 6.0 at
    all. Both refusals are enforced by the searcher, from the class recorded here.
    """
    jd0 = E.julian_day(t0.replace(hour=0, minute=0, second=0, microsecond=0))
    # The Julian Day NUMBER is JDN = floor(jd + 0.5), and `JDN mod 7` is the day of
    # the week with MONDAY = 0 exactly (JDN 0 was a Monday). Checked against a known
    # date rather than trusted: JD 2451545.0 is 2000-01-01 12:00 UT, JDN = 2451545,
    # 2451545 mod 7 = 5 = Saturday, and 2000-01-01 was a Saturday.
    #
    # Written as arithmetic on the Julian Day rather than as a datetime round-trip so
    # that a NULL time array is treated by exactly the same code as an observed one --
    # §P7-22(a)'s common-mode requirement applies to the calendar as much as to the
    # tide, and a `datetime` path for the observed events and an array path for the
    # nulls is precisely how a one-day offset gets in and never gets out.
    def dow_fn(day_float):
        jd = jd0 + np.asarray(day_float, dtype=np.float64)
        return np.mod(np.floor(jd + 0.5).astype(np.int64), 7).astype(np.float64)

    def utc_hour_fn(day_float):
        return MX.utc_hour_phase(np.asarray(day_float, dtype=np.float64))

    def local_hour_fn(day_float):
        return MX.local_solar_hour_phase(np.asarray(day_float, dtype=np.float64),
                                         float(site_lon))

    def month_fn(day_float):
        # month index 0..11 from the Julian Day, inverse of ephemeris.julian_day
        jd = jd0 + np.asarray(day_float, dtype=np.float64) + 0.5
        z = np.floor(jd).astype(np.int64)
        alpha = np.floor((z - 1867216.25) / 36524.25).astype(np.int64)
        a = z + 1 + alpha - (alpha // 4)
        b = a + 1524
        c = np.floor((b - 122.1) / 365.25).astype(np.int64)
        d = np.floor(365.25 * c).astype(np.int64)
        e = np.floor((b - d) / 30.6001).astype(np.int64)
        month = np.where(e < 14, e - 1, e - 13)
        return (month - 1).astype(np.float64)

    prov = _prov(
        source="UTC origin time from the ComCat catalogue; calendar arithmetic on the "
               "Julian Day",
        version="human-clock-v1",
        url_or_module="engine.marks_ext.utc_hour_phase / local_solar_hour_phase; "
                      "engine.ephemeris.julian_day",
        convention="day-of-week Monday = 0; UTC hour as a phase with 0 = 00:00 UTC; "
                   "local solar hour leads UTC by lon/15 h; month index January = 0",
        scalar_definition="calendar position of the event's UTC origin time")
    cols = [
        PropertyColumn("day_of_week", "human_clock", "categorical", "human-schedule",
                       dow_fn, prov, categories=tuple(range(7)), subdaily=False,
                       notes=[HUMAN_SCHEDULE_CARVE_OUT]),
        PropertyColumn("utc_hour_phase", "human_clock", "phase", "human-schedule",
                       utc_hour_fn, prov, subdaily=True,
                       notes=[HUMAN_SCHEDULE_CARVE_OUT]),
        PropertyColumn("local_solar_hour_phase", "human_clock", "phase",
                       "human-schedule", local_hour_fn, prov, subdaily=True,
                       notes=[HUMAN_SCHEDULE_CARVE_OUT, SITE_RULE]),
        PropertyColumn("month_of_year", "human_clock", "categorical",
                       "human-schedule", month_fn, prov,
                       categories=tuple(range(12)), subdaily=False,
                       notes=[HUMAN_SCHEDULE_CARVE_OUT]),
    ]
    return cols


# ---- family F: season ---------------------------------------------------------
def season_columns(t0: _dt.datetime):
    """Family F. Day-of-year as a phase; the hydrological index is NOT BUILT.

    §S2.1 F calls the family PARTIAL: day-of-year is free, local hydrological loading
    is a build. It is not built here, and it is ABSENT rather than approximated -- a
    property family silently running one column under a declaration that says two is
    exactly the defect `marks_ext.build_marks` returns an audit to prevent.
    """
    prov = _prov(
        source="engine/ephemeris.py solar longitude",
        version="season-v1",
        url_or_module="engine.ephemeris.ephemeris_table_at['annual_rad']",
        convention="0 = March equinox (solar longitude), radians in [0, 2pi)",
        scalar_definition="tropical-year phase of the event's UTC origin time")

    def fn(day_float):
        return E.ephemeris_table_at(t0, np.asarray(day_float, dtype=np.float64))[
            "annual_rad"]

    return [PropertyColumn("day_of_year_phase", "season", "phase",
                           "level-waveform-phase", fn, prov, subdaily=True,
                           notes=["hydrological-season index NOT BUILT; absent rather "
                                  "than approximated (§S2.1 F)"])]


# ---- family G: the event's own marks ------------------------------------------
def mark_columns(events):
    """Family G. Per-event marks -- the ONE family that is not a function of time.

    These columns are built from the event table itself, so `evaluate` takes an INDEX
    array rather than times, and the searcher must supply a `mark_null` (a
    composition-preserving or full-ETAS-with-triggering null, per SP-2) for any of
    them to promote. Without it they rank and do not promote, which is SP-2's own
    disposition rather than an omission.

    `dt_prior_days` and `cluster_member` carry the TRANCHE-B VIF_mark warning: both
    were measured ABOVE the 4.575 pooled fallback (§S2.1 G).
    """
    marks, audit = MX.build_marks(events)
    n = int(np.asarray(events["day_float"]).size)
    prov = _prov(
        source="USGS ComCat catalogue marks",
        version="F9-10-marks-v1",
        url_or_module="engine.marks_ext.build_marks",
        convention="marks as declared in marks_ext.MARK_DEFINITIONS",
        scalar_definition="per-event catalogue quantity")

    spec = {
        "depth_km": ("depth", "level", "catalogue-endogenous"),
        "dt_prior_days": ("dt_prior_days", "level", "clustering-derived"),
        "dist_nearest_prior_km": ("dist_nearest_prior_km", "level",
                                  "clustering-derived"),
        "cluster_member": ("cluster_member", "categorical", "clustering-derived"),
    }
    cols = []
    for name, (key, ptype, pclass) in spec.items():
        if key not in marks:
            continue
        vals = np.asarray(marks[key], dtype=np.float64)

        def fn(idx, _v=vals):
            i = np.asarray(idx, dtype=np.int64)
            return _v[i]

        p = dict(prov)
        p["scalar_definition"] = MX.MARK_DEFINITIONS.get(key, key)
        notes = ["MARK COLUMN: `evaluate` takes EVENT INDICES, not times -- this "
                 "family is not a function of time and its null is a mark null."]
        if key in ("dt_prior_days", "cluster_member"):
            notes.append("§S2.1 G: carries the tranche-B VIF_mark warning -- measured "
                         "ABOVE the 4.575 pooled fallback.")
        kw = {}
        if ptype == "level":
            kw["dwell"] = empirical_dwell(vals)
            notes.append("dwell measure for a MARK is its own empirical distribution, "
                         "not a time-occupancy measure: there is no time grid on which "
                         "a catalogue mark dwells. Stated so it is not mistaken for "
                         "the §P7-23(C) arcsine construction.")
        if ptype == "categorical":
            kw["categories"] = (0, 1)
        cols.append(PropertyColumn(name, "marks", ptype, pclass, fn, p,
                                   subdaily=False, notes=notes, **kw))
    return cols, audit


# ---- families D/E: geomagnetic and Earth rotation -----------------------------
def spaceweather_columns(t0: _dt.datetime, day_lo, day_hi, sample_minutes=180.0,
                         data_dir=None):
    """Families D and E, built ONLY if B4's downloads are on disk. Never faked.

    Returns (columns, audit). If a source is UNFETCHED the family is ABSENT and the
    audit says UNFETCHED with the reason -- the build brief's rule, and the same
    discipline `marks_ext.build_marks` applies to a missing lat/lon.
    """
    from . import spaceweather as SW
    series, audit = SW.load_series(data_dir=data_dir)
    cols = []
    grid = _dense_grid(day_lo, day_hi, sample_minutes)
    for key, meta in SW.SERIES.items():
        s = series.get(key)
        if s is None:
            continue

        def fn(day_float, _s=s):
            return SW.interpolate_at(_s, t0, np.asarray(day_float, dtype=np.float64))

        vals = fn(grid)
        if not np.isfinite(vals).any():
            continue
        prov = _prov(source=meta["source"], version=meta["version"],
                     url_or_module=meta["url"], convention=meta["convention"],
                     scalar_definition=meta["definition"])
        cols.append(PropertyColumn(
            meta["property"], meta["family"], "level", "level-waveform-phase",
            fn, prov, dwell=empirical_dwell(vals[np.isfinite(vals)]),
            subdaily=bool(meta.get("subdaily", True)),
            notes=[DWELL_RULE, meta.get("caveat", "")]))
    audit["n_columns"] = len(cols)
    return cols, audit


# --------------------------------------------------------- the matrix build ---
def build_property_matrix(t0: _dt.datetime, day_lo, day_hi, site_lat, site_lon,
                          families=("solid_tide", "ephemeris", "human_clock",
                                    "season"),
                          site_depth_km=0.0, events=None, data_dir=None,
                          sample_minutes=30.0):
    """The declared property lattice for ONE region site. Returns (columns, audit).

    `families` is the DECLARED family list and is hash-affecting: adding a property
    means a NEW scan declaration with a new `m` (SP-6.5), so a family list that
    changes between runs is two scans, not one.
    """
    cols, notes = [], {}
    fams = tuple(families)
    for f in fams:
        if f not in FAMILY_ORDER:
            raise ValueError("unknown property family %r; declared: %r"
                             % (f, FAMILY_ORDER))
    if "solid_tide" in fams:
        cols += solid_tide_columns(t0, day_lo, day_hi, site_lat, site_lon,
                                   site_depth_km, sample_minutes)
    if "ephemeris" in fams:
        cols += ephemeris_columns(t0, day_lo, day_hi, sample_minutes)
    if "human_clock" in fams:
        cols += human_clock_columns(t0, site_lon)
    if "season" in fams:
        cols += season_columns(t0)
    if "marks" in fams:
        if events is None:
            raise ValueError("family 'marks' needs the event table")
        mc, ma = mark_columns(events)
        cols += mc
        notes["marks"] = ma
    if "geomagnetic" in fams or "earth_rotation" in fams:
        sc, sa = spaceweather_columns(t0, day_lo, day_hi, data_dir=data_dir)
        keep = [c for c in sc if c.family in fams]
        cols += keep
        notes["spaceweather"] = sa
    if "clocks" in fams:
        notes["clocks"] = (
            "family 'clocks' (§S2.1 H) is NOT BUILT: engine/clocks.py has the "
            "machinery and the mandatory F8-15 random-clock control, but a per-event "
            "clock evaluation is a build, and each clock is PRICED AS A NEW TEST per "
            "§P7-5(2). Absent rather than approximated.")

    audit = {
        "property_rule_id": PROPERTY_RULE_ID,
        "declared_families": list(fams),
        "n_columns": len(cols),
        "columns": [c.record() for c in cols],
        "site": {"lat": float(site_lat), "lon": float(site_lon),
                 "depth_km": float(site_depth_km), "rule": SITE_RULE},
        "window_days": [float(day_lo), float(day_hi)],
        "sp2_rule": SP2_RULE,
        "dwell_rule": DWELL_RULE,
        "subdaily_rule": SUBDAILY_RULE,
        "mag_strata_rule": MAG_STRATA_RULE,
        "family_notes": notes,
    }
    return cols, audit


def redundancy_collapse(columns, day_float, threshold=0.90):
    """§S3.5: collapse properties with |r| > threshold into ONE declared representative.

    The representative is chosen by a PRE-DECLARED rule -- earliest in the family
    table, then earliest in build order -- never by which one scored better, which is
    the only choice rule that cannot be tuned after seeing a scan. Circular columns are
    compared on the circular-circular correlation of the two angles; level columns on
    their dwell-corrected PIT values; a phase and a level are compared on the level's
    PIT mapped to an angle, which is exactly how the searcher will consume it anyway.

    Returns (kept_columns, audit). This REDUCES the declared multiplicity honestly
    rather than inflating it for appearance.
    """
    vals, keys = {}, []
    for c in columns:
        if c.ptype == "categorical":
            continue                     # no correlation is defined; never collapsed
        try:
            vals[c.name] = np.asarray(c.statistic_values(day_float),
                                      dtype=np.float64)
        except Exception:                # a column that cannot evaluate cannot collapse
            continue
        keys.append(c.name)
    order = {c.name: (FAMILY_ORDER.index(c.family), i)
             for i, c in enumerate(columns)}
    dropped, pairs = {}, []
    for i, a in enumerate(keys):
        if a in dropped:
            continue
        for b in keys[i + 1:]:
            if b in dropped:
                continue
            r = _circ_circ_corr(vals[a], vals[b])
            if abs(r) > float(threshold):
                lose = b if order[a] <= order[b] else a
                keep = a if lose == b else b
                dropped[lose] = keep
                pairs.append({"a": a, "b": b, "r": float(r), "kept": keep,
                              "dropped": lose})
    kept = [c for c in columns if c.name not in dropped]
    return kept, {
        "threshold": float(threshold),
        "n_before": len(columns), "n_after": len(kept),
        "collapsed": pairs,
        "rule": ("SEARCHER.md §S3.5: properties with pairwise |r| > %.2f are collapsed "
                 "into ONE declared representative, chosen by a PRE-DECLARED rule "
                 "(earliest in the family table), never by which scored better."
                 % float(threshold)),
    }


def _circ_circ_corr(a, b):
    """Jammalamadaka-Sarma circular-circular correlation of two angle series."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if a.size < 3:
        return 0.0
    abar = np.arctan2(np.sin(a).mean(), np.cos(a).mean())
    bbar = np.arctan2(np.sin(b).mean(), np.cos(b).mean())
    sa, sb = np.sin(a - abar), np.sin(b - bbar)
    den = np.sqrt((sa ** 2).sum() * (sb ** 2).sum())
    return float((sa * sb).sum() / den) if den > 0 else 0.0


def property_declaration_digest(audit) -> str:
    """A hash over the columns' DECLARATION rows -- values never enter the hash."""
    block = {"property_rule_id": PROPERTY_RULE_ID,
             "columns": audit["columns"],
             "declared_families": audit["declared_families"]}
    return splits.config_hash(block)[:12]


if __name__ == "__main__":          # pragma: no cover - operator convenience
    import json
    t0 = _dt.datetime(2000, 1, 1)
    cols, audit = build_property_matrix(t0, 0.0, 365.0, 38.0, 142.0)
    audit["digest"] = property_declaration_digest(audit)
    print(json.dumps(audit, indent=2, default=str))
