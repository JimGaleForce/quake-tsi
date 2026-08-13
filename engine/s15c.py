"""S-15(c) UNMEASURABLE-BY-WINDOW -- the catalog-wide sweep §P7-14(c) asks for.

§P7-14(c), after two control rows in a battery designed to produce nothing came back
at p < 0.05:

  > *"`nc_jupiter_saturn_synodic` (1.06 cycles in window) and `nc_metonic` (1.11
  > cycles) produced p < 0.05 in a battery designed to produce nothing. That is not a
  > control failure; it is the controls doing their job and surfacing a gap in my own
  > floor formula. §P7-1(b) prices MULTIPLICITY and DISPERSION and is SILENT ON
  > IDENTIFIABILITY."*
  >
  > **Catalog-wide handling: sweep all 281 entries and flag before B is declared.**

The clause itself lives in `engine/floors.py` -- deliberately, so that FUTURE
FEATURES INHERIT IT on the same call path they already use for their power floor.
This module is only the sweep: it enumerates every periodicity this programme has
declared anywhere (the engine's own features, the F9-20 named controls, the
period-scan grid, and the catalog's named long-period entries) and applies the
clause to all of them at once.

WHAT THE SWEEP IS EXPECTED TO FIND, WRITTEN DOWN BEFORE IT RUNS
---------------------------------------------------------------
§P7-14(c) names three in advance: Metonic (1.11 cycles), Jupiter-Saturn synodic
(1.06), and -- *"importantly"* -- the 11-year solar cycle at 1.92. That last one is
not hypothetical: `mine.py`'s own docstring records an earlier build handing F10.7 a
p at the resolution floor with **z = 32**, an artifact of exactly this class. **The
clause retroactively explains a known corpse, which is the best evidence a new rule
can have.** `EXPECTED_CATCHES` records that prediction in code; `sweep()` checks it
and says so, because a prediction that is not checked is a prediction that was never
made.

THE CATALOG PERIOD TABLE IS DECLARED, NOT SCRAPED
-------------------------------------------------
Every row below carries its catalog entry ID and its period as the catalog states it.
Scraping the markdown would be fragile and unauditable; this table can be diffed
against the catalog by eye, and a row that is wrong is wrong visibly.
"""

from __future__ import annotations

import numpy as np

from . import ephemeris as eph, floors, mine as M

DAYS_PER_YEAR = 365.25

# ------------------------------------------------------------ declared periods --
# (catalog_id, name, period_days, note). Only entries with a DECLARED period appear:
# S-15(c) is a statement about periodic features and says nothing about the rest.
CATALOG_PERIODS = (
    ("F1-01", "moon_synodic_phase", eph.SYNODIC_MONTH, "lunar synodic month"),
    ("F1-02", "moon_anomalistic_phase", eph.ANOMALISTIC_MONTH, "anomalistic month"),
    ("F1-03", "moon_draconic_phase", eph.DRACONIC_MONTH, "draconic month"),
    ("F1-04", "annual_phase", eph.TROPICAL_YEAR, "tropical year"),
    ("F1-05", "lunar_nodal_phase_18_6y", 6798.4, "lunar nodal regression, 18.61 y"),
    ("F1-07", "apsidal_precession_phase_8_85y", 3232.6, "lunar perigee precession"),
    ("F1-19", "lunar_semidiurnal_phase_M2", 0.5175, "M2; structurally notched"),
    ("F1-20", "solar_semidiurnal_phase_S2", 0.5, "S2; exact day-binning zero"),
    ("F1-36", "planetary_synodic_jupiter_saturn", M.JUPITER_SATURN_SYNODIC_DAYS,
     "Jupiter-Saturn synodic; a named F9-20 mechanism-free control"),
    ("F1-37", "planetary_synodic_venus", 583.9, "Earth-Venus synodic"),
    ("F1-39", "solar_rotation_phase_27d", 27.2753, "solar rotation (Carrington)"),
    ("F1-40", "solar_cycle_phase_11y", 11.0 * DAYS_PER_YEAR,
     "Schwabe solar cycle -- the F10.7 z=32 corpse's own band"),
    ("F2-01", "spring_neap_phase", eph.SYNODIC_MONTH / 2, "spring-neap beat"),
    ("F2-02", "half_draconic_phase", eph.DRACONIC_MONTH / 2, "declination-tide beat"),
    ("F2-03", "perigean_spring_beat", 411.8, "perigean-spring cycle"),
    ("F2-04", "eclipse_year_beat", 173.3, "eclipse half-year"),
    ("F2-05", "beat_annual_draconic", 205.9, "annual-draconic beat"),
    ("F2-06", "beat_annual_anomalistic", 182.6, "annual-anomalistic beat"),
    ("F2-08", "saros_phase", 6585.32, "saros cycle"),
    ("F2-09", "metonic_phase", M.METONIC_DAYS,
     "Metonic 19 y; a named F9-20 mechanism-free control"),
    ("F3-19", "LOD_annual_amplitude", 400.0, "trailing-window LOD annual amplitude"),
    ("F3-22", "chandler_wobble_phase", 433.0, "Chandler wobble"),
    ("F3-23", "chandler_annual_beat_6_4y", 2350.0, "Chandler-annual beat, 6.4 y"),
    ("F6-01", "seasonal_phase_by_hemisphere", eph.TROPICAL_YEAR, "seasonal"),
    ("F6-05", "ENSO_ONI_index", 3.5 * DAYS_PER_YEAR, "ENSO quasi-period, ~3.5 y"),
    ("F7-01", "diurnal_detection_amplitude", 1.0,
     "the observer's diurnal cycle; exact day-binning zero on the count path"),
    ("F7-03", "weekly_detection_amplitude", 7.0, "human-operations week"),
)

# §P7-14(c)'s own advance list, checked by `sweep()`.
EXPECTED_CATCHES = ("metonic_phase", "planetary_synodic_jupiter_saturn",
                    "solar_cycle_phase_11y")

CORPSE_NOTE = (
    "mine.py's own docstring records an earlier build handing F10.7 (the 11 y solar "
    "cycle, 1.92 cycles in this window) a p at the resolution floor with z = 32. "
    "S-15(c) retroactively explains that corpse, which is the best evidence a new "
    "rule can have (§P7-14(c)).")

ORTHOGONALITY_NOTE = (
    "S-15(c) is ORTHOGONAL to the S-15 power floor. A feature can pass a_min() and "
    "fail this, which is exactly what happened to the two Tranche A control rows. "
    "Amplitude floors are about N; this is about identifiability, and no N repairs "
    "it.")


def engine_feature_periods(t0=None, n_days=2000):
    """The engine's OWN declared periods, read off the features rather than retyped."""
    import datetime as _dt
    t0 = t0 or _dt.datetime(2000, 1, 1)
    return [(f.name, float(f.period_hint))
            for f in M.ephemeris_features(t0, int(n_days))
            if f.period_hint]


def period_scan_grid_check(record_days=floors.EXPLORATION_RECORD_DAYS,
                           period_max=None):
    """Does the period scan's own upper cap respect the clause? Reported, not fixed.

    `mine_session.PERIOD_MAX` is 4,000 d against a cut at record/3 = 2,572 d, so the
    top ~36% of the scan's declared grid sits in territory S-15(c) says is not
    identifiable. That is a live finding about a declared grid, and it belongs in the
    sweep's output rather than in a silent clamp: changing a declared grid is a new
    declaration (§P6-3 rule 5), not a build fix.
    """
    from . import mine_session as ms
    pmax = float(period_max if period_max is not None else ms.PERIOD_MAX)
    cut = floors.max_identifiable_period(record_days)
    return {
        "period_scan_max_days": pmax,
        "s15c_cut_days": cut,
        "grid_exceeds_cut": bool(pmax > cut),
        "fraction_of_log_grid_above_cut": (
            float(np.log(pmax / cut) / np.log(pmax / ms.PERIOD_MIN))
            if pmax > cut else 0.0),
        "disposition": ("REPORTED, NOT CLAMPED. The period grid is a DECLARED "
                        "config value; narrowing it is a new declaration and a new "
                        "config hash (§P6-3 rule 5), not a build-time fix. Any peak "
                        "the scan reports above %.0f d must carry the "
                        "UNMEASURABLE-BY-WINDOW verdict." % cut),
    }


def sweep(record_days=floors.EXPLORATION_RECORD_DAYS,
          min_cycles=floors.MIN_CYCLES_IN_WINDOW, include_engine=True):
    """The catalog-wide S-15(c) sweep. Returns the table, the list and the check."""
    items = [(name, p) for _cid, name, p, _n in CATALOG_PERIODS]
    seen = {n for n, _ in items}
    if include_engine:
        items += [(n, p) for n, p in engine_feature_periods() if n not in seen]
    tab = floors.window_sweep(items, record_days, min_cycles)

    by_name = {r["feature"]: r for r in tab["rows"]}
    meta = {name: (cid, note) for cid, name, _p, note in CATALOG_PERIODS}
    for r in tab["rows"]:
        cid, note = meta.get(r["feature"], (None, "engine feature (period_hint)"))
        r["catalog_id"] = cid
        r["catalog_note"] = note

    caught = [n for n in EXPECTED_CATCHES
              if by_name.get(n, {}).get("verdict") == floors.UNMEASURABLE_BY_WINDOW]
    tab.update({
        "expected_catches": list(EXPECTED_CATCHES),
        "expected_catches_confirmed": caught,
        "expected_catches_missed": [n for n in EXPECTED_CATCHES if n not in caught],
        "prediction_check": ("§P7-14(c) named %d entries in advance; the sweep "
                             "caught %d of them"
                             % (len(EXPECTED_CATCHES), len(caught))),
        "corpse_note": CORPSE_NOTE,
        "orthogonality_note": ORTHOGONALITY_NOTE,
        "period_scan_grid": period_scan_grid_check(record_days),
        "what_this_does_not_do": (
            "an UNMEASURABLE-BY-WINDOW feature is scored NEITHER WAY. It is not a "
            "null, it may not be counted as one, and it may not be counted as a "
            "detection either. It is reported in the headline fraction and removed "
            "from both numerators."),
    })
    return tab


def table_lines(rep):
    """The sweep as report-ready markdown rows, longest period first."""
    out = ["| catalog | feature | period (d) | cycles in window | verdict |",
           "| --- | --- | ---: | ---: | --- |"]
    for r in rep["rows"]:
        out.append("| %s | `%s` | %.4g | %.3g | %s |"
                   % (r.get("catalog_id") or "-", r["feature"], r["period_days"],
                      r["cycles_in_window"], r["verdict"]))
    return out


def main(argv=None):
    import argparse
    import json
    ap = argparse.ArgumentParser(description="S-15(c) catalog-wide sweep (§P7-14(c))")
    ap.add_argument("--record-days", type=float,
                    default=floors.EXPLORATION_RECORD_DAYS)
    ap.add_argument("--json", default=None)
    a = ap.parse_args(argv)
    rep = sweep(record_days=a.record_days)
    print("S-15(c) UNMEASURABLE-BY-WINDOW sweep -- %s" % floors.WINDOW_CLAUSE_SOURCE)
    print("record %.0f d; minimum %g cycles; cut at period > %.0f d"
          % (rep["record_days"], rep["min_cycles_required"], rep["cut_period_days"]))
    print("")
    for line in table_lines(rep):
        print(line)
    print("")
    print("UNMEASURABLE-BY-WINDOW: %d of %d examined"
          % (rep["n_unmeasurable_by_window"], rep["n_examined"]))
    for r in rep["unmeasurable"]:
        print("  %-34s %10.4g d  %.3g cycles" % (r["feature"], r["period_days"],
                                                 r["cycles_in_window"]))
    print("\n" + rep["prediction_check"])
    print(rep["corpse_note"])
    g = rep["period_scan_grid"]
    print("\nperiod-scan grid: max %.0f d vs cut %.0f d -> exceeds = %s"
          % (g["period_scan_max_days"], g["s15c_cut_days"], g["grid_exceeds_cut"]))
    if a.json:
        with open(a.json, "w", encoding="utf-8") as fh:
            json.dump(rep, fh, indent=1, default=float)
        print("wrote %s" % a.json)
    return rep


if __name__ == "__main__":
    main()
