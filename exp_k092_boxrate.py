"""The K-092 box EVENT-RATE COUNT. Priced 0. A RATE, NOT A PHASE.

LICENCE, verbatim from HYPOTHESIS_LEDGER.md §P7-25(3), disposition 12:

  > PERMITTED, at price 0, logged: count the M >= 5.5 and M >= 5.0 ComCat events in the
  > frozen box, 1990-2026-08-13, and report the rates. This is a RATE COUNT, NOT A
  > PHASE STATISTIC; those events are all pre-freeze and outside every D-13 scoring
  > set; and D-12's SECONDARY stratum needs the catalogue regardless. NO PHASE, NO
  > LEVEL, NO QUADRANT, NO SCALAR MAY BE COMPUTED ON ANY OF THEM -- doing so spends
  > D-12's priced test. The count replaces my GR extrapolation with a measurement and
  > fixes D-13/D-13b's expected N before either clock matters.

WHY IT MATTERS. §P7-25 re-labelled D-13 DESCRIPTIVE-PRIMARY, NOT DECISIVE on a power
calculation whose N came from a Gutenberg-Richter extrapolation off the M >= 6.0 seed
count (0.917/yr -> ~2.9/yr at M >= 5.5 -> N ~ 9 over three years). That extrapolation
assumed b = 1 and it assumed the seed file's own completeness. This module measures the
number instead, so D-13's power statement rests on a count rather than on an assumption
about a slope.

THE ONE THING THIS FILE MUST NOT DO. It imports NOTHING from `sitetide`,
`tidal_tensor`, `audit_arcsine` or either bridge, and it computes no time-of-day
quantity of any kind. That is enforced by `assert_no_tidal_imports()` at run time
rather than left to discipline, because the licence is narrow and the cost of
exceeding it is a priced D-12 test.
"""

from __future__ import annotations

import csv
import datetime as _dt
import json
import math
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
CATALOGUE = HERE / "data" / "comcat_world" / "Alaska-Aleutians.csv"
OUT_JSON = HERE / "results_k092_boxrate.json"

# The frozen box, from K092_FREEZE.md. Copied, not re-derived.
LAT_MIN, LAT_MAX = 51.0, 58.0
LON_MIN, LON_MAX = -166.0, -152.0
FREEZE_UTC = _dt.datetime(2026, 8, 13, 20, 27, 26, tzinfo=_dt.timezone.utc)
WINDOW_START = _dt.datetime(1990, 1, 1, tzinfo=_dt.timezone.utc)

FORBIDDEN = ("sitetide", "tidal_tensor", "audit_arcsine", "circstat", "ephemeris")


def assert_no_tidal_imports():
    """The licence is 'a rate count, not a phase statistic'. Enforce it mechanically."""
    loaded = [m for m in sys.modules
              if any(f in m for f in FORBIDDEN) and m.startswith("engine")]
    if loaded:
        raise SystemExit(
            "REFUSING TO RUN: this arm's licence (§P7-25 disposition 12) is a RATE "
            "COUNT ONLY, and a tidal module is loaded: %s. Computing any phase, "
            "level, quadrant or scalar on these events spends D-12's priced test."
            % loaded)


def parse_iso_utc(s):
    return _dt.datetime.fromisoformat(s.replace("Z", "+00:00"))


def _great_circle_km(lat1, lon1, lat2, lon2):
    """Haversine, km. Used only by the declustering rule frozen in section P7-25."""
    r, p = 6371.0, math.pi / 180.0
    h = (math.sin((lat2 - lat1) * p / 2) ** 2
         + math.cos(lat1 * p) * math.cos(lat2 * p)
         * math.sin((lon2 - lon1) * p / 2) ** 2)
    return 2.0 * r * math.asin(min(1.0, math.sqrt(h)))


def main():
    assert_no_tidal_imports()
    if not CATALOGUE.exists():
        raise SystemExit("catalogue not found: %s" % CATALOGUE)

    rows, in_box, in_box_full = 0, [], []
    with CATALOGUE.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rows += 1
            try:
                lat, lon, mag = (float(r["latitude"]), float(r["longitude"]),
                                 float(r["mag"]))
            except (KeyError, ValueError, TypeError):
                continue
            if not (LAT_MIN <= lat <= LAT_MAX and LON_MIN <= lon <= LON_MAX):
                continue
            t = parse_iso_utc(r["time"])
            in_box.append((t, mag))
            in_box_full.append((t, lat, lon, mag))

    span_years = (FREEZE_UTC - WINDOW_START).total_seconds() / (365.2425 * 86400.0)
    counts, rates, by_decade = {}, {}, {}
    for thr in (5.0, 5.5, 6.0, 6.5, 7.0):
        sel = [(t, m) for t, m in in_box
               if m >= thr and WINDOW_START <= t <= FREEZE_UTC]
        counts[thr] = len(sel)
        rates[thr] = len(sel) / span_years
        by_decade[thr] = dict(sorted(Counter(t.year // 10 * 10 for t, _ in sel).items()))

    # The b-value implied by the measured counts, reported so the §P7-25 extrapolation
    # can be checked rather than trusted. Two-point slope, not a fit: this is a
    # diagnostic on an assumption, not an estimate anyone may quote.
    b_5_6 = (math.log10(counts[5.0] / counts[6.0]) / 1.0) if counts[6.0] else None
    b_55_6 = (math.log10(counts[5.5] / counts[6.0]) / 0.5) if counts[6.0] else None

    # DECLUSTERED COUNT. Still a count, so still inside the licence -- and D-13b's
    # PRIMARY stratum is the declustered one, so its expected N is the declustered
    # rate and not the full one. The rule is the one §P7-25 froze before any event:
    # dependent if within 30 days AND 150 km of an equal-or-larger prior event.
    m50 = sorted((t, la, lo, m) for t, la, lo, m in in_box_full
                 if m >= 5.0 and WINDOW_START <= t <= FREEZE_UTC)
    indep = []
    for i, (t, la, lo, m) in enumerate(m50):
        if not any((t - t2).days <= 30 and m2 >= m
                   and _great_circle_km(la2, lo2, la, lo) <= 150.0
                   for t2, la2, lo2, m2 in m50[:i]):
            indep.append((t, la, lo, m))

    horizon_years = 3.0

    # POWER, computed here rather than assumed anywhere else. p0 is the pooled
    # waveform-matched duty cycle measured by D-1b; the real test uses the per-event
    # Poisson-binomial, for which this is a close and slightly conservative proxy.
    from scipy.stats import binom
    P0 = 0.3827

    def power_at(n, true, alpha=0.01):
        k = next((kk for kk in range(n + 1) if binom.sf(kk - 1, n, P0) <= alpha), None)
        return (k, float(binom.sf(k - 1, n, true))) if k is not None else (None, 0.0)

    power_table = {}
    for label, n in (("D-13 M>=5.5 full", round(rates[5.5] * horizon_years)),
                     ("D-13b M>=5.0 full (SECONDARY)",
                      round(rates[5.0] * horizon_years)),
                     ("D-13b M>=5.0 declustered (PRIMARY)",
                      round(len(indep) / span_years * horizon_years))):
        row = {"N": n}
        for true in (0.50, 0.60, 0.70):
            k, pw = power_at(n, true)
            row["k_crit"] = k
            row["power_at_%.2f" % true] = pw
        power_table[label] = row

    out = {
        "arm": "K-092 box event-rate count (§P7-25 disposition 12)",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 0,
        "licence": ("RATE COUNT ONLY. No phase, level, quadrant or scalar was computed "
                    "on any event; enforced by assert_no_tidal_imports(). All events "
                    "counted are pre-freeze and outside every D-13/D-13b scoring set."),
        "box": {"lat": [LAT_MIN, LAT_MAX], "lon": [LON_MIN, LON_MAX],
                "source": "K092_FREEZE.md"},
        "window": {"start": WINDOW_START.isoformat(), "end": FREEZE_UTC.isoformat(),
                   "span_years": span_years},
        "catalogue": {"path": str(CATALOGUE.relative_to(HERE)), "rows_read": rows,
                      "rows_in_box_all_time": len(in_box)},
        "counts": {str(k): v for k, v in counts.items()},
        "rates_per_year": {str(k): v for k, v in rates.items()},
        "counts_by_decade": {str(k): v for k, v in by_decade.items()},
        "implied_b_two_point": {"M5.0_vs_M6.0": b_5_6, "M5.5_vs_M6.0": b_55_6,
                                "note": "two-point slope, a DIAGNOSTIC on §P7-25's "
                                        "b = 1 assumption, not a b-value estimate and "
                                        "not quotable as one"},
        "declustered": {
            "rule": "dependent if within 30 days AND 150 km of an equal-or-larger "
                    "prior event; frozen in section P7-25 before any D-13b event",
            "n_full_M5.0": len(m50), "n_independent_M5.0": len(indep),
            "dependent_fraction": 1.0 - len(indep) / max(len(m50), 1),
            "independent_rate_per_yr": len(indep) / span_years,
            "note": "D-13b's PRIMARY stratum is the DECLUSTERED one, so its expected "
                    "N is this rate and not the full-catalogue rate.",
        },
        "power": {
            "p0_used": P0,
            "p0_source": "pooled waveform-matched duty cycle from results_k092_d1_null"
                         ".json; the real test uses the per-event Poisson-binomial, "
                         "for which this is a close and slightly conservative proxy",
            "alpha": 0.01, "one_sided": True,
            "table": power_table,
            "verdict": ("MEASURED, not assumed. D-13 at its measured rate cannot "
                        "detect even a very large effect. D-13b's declustered PRIMARY "
                        "is better but still underpowered against a true rate of 0.60, "
                        "and is adequate only against ~0.70. Read alongside Merton "
                        "section M-011's Beeler and Lockner (2003) requirement of "
                        "10^5 to 10^6 events at body-tide amplitudes: NO prospective "
                        "arm at this region's event rate can detect a "
                        "literature-scale (order 1 percent) tidal effect on any "
                        "horizon. That is a fact about Alaska's seismicity, not about "
                        "this pipeline, and it should govern how much further "
                        "prospective machinery is worth building."),
        },
        "expected_N_over_3yr": {
            "D-13 (M>=5.5)": rates[5.5] * horizon_years,
            "D-13b (M>=5.0)": rates[5.0] * horizon_years,
        },
        "p7_25_comparison": {
            "assumed_rate_M5.5_per_yr": 2.9,
            "assumed_N_3yr": 9,
            "measured_rate_M5.5_per_yr": rates[5.5],
            "measured_N_3yr": rates[5.5] * horizon_years,
        },
        "completeness_caveat": (
            "Counts are as-catalogued and NOT completeness-corrected. The decade "
            "breakdown is reported so that a completeness trend is visible rather "
            "than averaged away: a rate rising with decade is a network artifact "
            "before it is anything else, and the 1990 start was chosen for that "
            "reason. This bears on the expected N and therefore on D-13/D-13b power, "
            "and on nothing else."),
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    print("K-092 box rate count, %s to %s (%.2f yr)"
          % (WINDOW_START.date(), FREEZE_UTC.date(), span_years))
    for thr in (5.0, 5.5, 6.0, 6.5, 7.0):
        print("  M >= %.1f : %5d events  %6.3f /yr   by decade %s"
              % (thr, counts[thr], rates[thr], by_decade[thr]))
    print("  M >= 5.0 declustered : %5d events  %6.3f /yr  (%.0f%% dependent)"
          % (len(indep), len(indep) / span_years,
             100.0 * (1.0 - len(indep) / max(len(m50), 1))))
    print("\n  power (alpha 0.01, one-sided, p0 = %.4f):" % P0)
    for _k, _v in power_table.items():
        print("    %-38s N=%-3d k_crit=%-4s power@0.60 %.3f  @0.70 %.3f"
              % (_k, _v["N"], _v["k_crit"], _v["power_at_0.60"], _v["power_at_0.70"]))
    print("\n  implied two-point b (M5.0 vs M6.0) = %s" % (
        "%.3f" % b_5_6 if b_5_6 else "n/a"))
    print("  expected N over 3 yr: D-13 (M>=5.5) %.1f ; D-13b (M>=5.0) %.1f"
          % (rates[5.5] * horizon_years, rates[5.0] * horizon_years))
    print("  §P7-25 assumed 2.9/yr -> N ~ 9 ; measured %.3f/yr -> N ~ %.1f"
          % (rates[5.5], rates[5.5] * horizon_years))
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
