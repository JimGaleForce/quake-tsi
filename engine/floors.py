"""The S-15 declaration floor, in one place, with its measured VIF attached.

HYPOTHESIS_LEDGER.md §P7-1(b), scale RETRACTED and re-measured by §P7-8(a):

    A_min = sqrt(VIF) * (z_alpha + z_0.80) * sqrt(2 / N)

`A_min` is the smallest SINUSOIDAL RATE-MODULATION AMPLITUDE (as a fraction of
the baseline rate) that the count path can declare at 80% power, at a declared
two-sided operating threshold `alpha`, over `N` events.

WHY THIS MODULE EXISTS -- §P7-8(d), stated so nobody has to reconstruct it
---------------------------------------------------------------------------
§P7-8(d) ruled: *"F9-19 / G-M1 arm (ii) plants signals and demands recovery at
Ahat/A in [0.8, 1.2]: a signal planted below ~15% at global aggregation cannot be
recovered at that tolerance, and the gate would fail for POWER reasons while
reading as an INSTRUMENT failure. Plant at >= 2x the operative floor ... This
ruling exists to prevent a false G-M1 failure that would otherwise have been
recorded against the pipeline."*

Read that twice, because it is the whole point. A plant below the floor produces
a test that fails. The failure is real. But its CAUSE is that the experiment had
no power to see what it planted -- and the failure would be written down as
"the recovery instrument does not work", which is a different claim about a
different object, and a false one. The floor is not a formality here: it is the
line between "the pipeline is broken" and "this test could never have worked".
Every harness that plants a signal must therefore compute its floor from THIS
module at the N and alpha it actually runs at, and plant at least `PLANT_FACTOR`
times it -- never at a hardcoded amplitude chosen before the VIF was measured.

CONVENTIONS, CARRIED (S-18, §P7-8(e): a number may not enter an arithmetic
unless it carries its convention)
-------------------------------------------------------------------------
  * `alpha`  -- TWO-SIDED. z_alpha = Phi^-1(1 - alpha/2). This is the convention
                §P7-1(a) used in its own arithmetic (alpha = 0.1/259 -> z = 3.549).
  * power    -- fixed at 80%; z_0.80 = Phi^-1(0.80) = 0.8416. Not adjustable,
                because §P7-1(b) fixed it.
  * VIF      -- MEASURED, not assumed: 24.0818 over the 2-df phase features of
                `session_20260811T022953` (results_f4_58_vif.json,
                `overall_vif_df2_phase.median`). Dispersion is INSIDE the floor
                and OUTSIDE every p-value -- §P7-8(a) is explicit that no p-value,
                BH threshold or max-statistic may be "corrected" by the VIF, since
                the surrogate nulls already carry the dispersion by construction.
  * N        -- the event count the statistic is actually computed over. For a
                per-region statistic that is the REGION's N, not the domain's.

PROVISIONAL, AND WHY THE DIRECTION IS SAFE
------------------------------------------
§P7-8(b) adopts VIF = 24.08 PROVISIONALLY and requires item (ii) -- the
true-VIF=1 simulation control (`f4_58_vif_control.py`,
results_f4_58_vif_control.json) -- before the floor may retire any catalog entry
on floor grounds or be quoted externally. A larger VIF declares MORE things
unmeasurable, so over-stating it costs scope and cannot manufacture a claim.
Planting is the one place where the direction is also safe: a floor that is too
big makes the plant bigger, which cannot turn a broken instrument into a passing
test.
"""

from __future__ import annotations

import json
import math
import os

import numpy as np
from scipy import stats

# ------------------------------------------------------------------ constants --
Z_POWER_80 = float(stats.norm.ppf(0.80))
POWER = 0.80

# §P7-8(a): the F4-58 MEASUREMENT over 2-df phase features, stable to three
# significant figures across all four sessions on disk (24.08 / 24.09 / 24.66 /
# 24.09). NOT §P7-1(a)'s inferred 3.94 (falsified) and not Kepler's 9.7.
MEASURED_VIF_DF2_PHASE = 24.081827389301417
MEASURED_VIF_SOURCE = (
    "results_f4_58_vif.json :: sessions[session_20260811T022953]"
    ".overall_vif_df2_phase.median (HYPOTHESIS_LEDGER.md §P7-8(a))")
MEASURED_VIF_STATUS = (
    "PROVISIONAL (§P7-8(b)). Item (ii), the true-VIF=1 simulation control, has "
    "now RUN (f4_58_vif_control.py -> results_f4_58_vif_control.json): on 24 "
    "catalogs Poisson-simulated from the fitted ETAS lambda the same estimator "
    "returns VIF = 1.00 flat from 30 d to 800 d blocks (log-log slope "
    "-0.002 +/- 0.002, p = 0.47), reproducing 4.1% of the 24.08 measured on the "
    "Earth. That is Reading A and it supports the flat 24.08; PROMOTION FROM "
    "PROVISIONAL IS THE POPPER SEAT'S CALL, not this module's, so the label "
    "stands until it is made.")

# Declared operating thresholds, BH q = 0.10 at the most conservative rung.
FDR_Q = 0.10
M_DECLARED_TRANCHE_A = 713
M_DECLARED_20260811 = 259
ALPHA_TRANCHE_A = FDR_Q / M_DECLARED_TRANCHE_A
ALPHA_20260811 = FDR_Q / M_DECLARED_20260811

# §P7-8(d): "Plant at >= 2x the operative floor".
PLANT_FACTOR = 2.0

RESULTS_JSON = "results_f4_58_vif.json"


# ---------------------------------------------------------- §P7-14(c): S-15(c) --
# UNMEASURABLE-BY-WINDOW. A periodic feature with fewer than 3 full cycles in the
# analysis window is UNMEASURABLE-BY-WINDOW *regardless of N*, declared before the
# run, reported in the headline fraction, and scored neither way.
#
# THIS IS AN IDENTIFIABILITY LIMIT AND IS ORTHOGONAL TO THE POWER FLOOR ABOVE.
# `a_min` prices multiplicity and dispersion and is silent on identifiability: a
# feature with ~1 cycle in the window is not distinguishable from a trend or an
# offset, so the fitted "period" is reporting record length, not periodicity. No
# amount of N fixes it, which is precisely why it cannot live inside a formula whose
# only lever is N. A feature can PASS the amplitude floor and FAIL this, and that is
# exactly what happened to `nc_jupiter_saturn_synodic` (1.06 cycles) and `nc_metonic`
# (1.11 cycles) in the Tranche A control battery -- two p < 0.05 rows in a battery
# designed to produce nothing.
#
# THE THRESHOLD IS NOT NEW AND IS NOT MINE. §P7-14(c) adopts the rule this programme
# already implemented and justified in `mine.py:harmonic_ladder`, which refuses any
# rung longer than record/3: *"with fewer than three observed cycles an epoch fold is
# not measuring a period, and a ladder that 'wins' at 5827 d on 7716 days is
# reporting the record length."* It is stricter than the catalog's `< 2 cycles` and
# it retroactively explains a known corpse -- `mine.py`'s own docstring records an
# earlier build handing F10.7 (11 y solar cycle, 1.92 cycles in window) a p at the
# resolution floor with z = 32.
MIN_CYCLES_IN_WINDOW = 3.0
WINDOW_CLAUSE = "S-15(c)"
WINDOW_CLAUSE_SOURCE = (
    "HYPOTHESIS_LEDGER.md §P7-14(c), threshold inherited from "
    "mine.py:harmonic_ladder's hi_cap = record/3")
UNMEASURABLE_BY_WINDOW = "UNMEASURABLE-BY-WINDOW"
MEASURABLE_BY_WINDOW = "MEASURABLE-BY-WINDOW"

# The v2 exploration window, stated so a caller that omits `record_days` gets the
# declared record rather than a guess. 7,716 days -> the cut falls at 2,572 d.
EXPLORATION_RECORD_DAYS = 7716.0


class UnmeasurableByWindow(AssertionError):
    """A period the analysis window cannot identify, at any N (S-15(c))."""


class PlantBelowFloor(AssertionError):
    """A planted amplitude that the test could not have recovered anyway."""


# ------------------------------------------------------------------- formula --
def z_alpha(alpha: float) -> float:
    """Two-sided critical z at the declared threshold: Phi^-1(1 - alpha/2)."""
    return float(stats.norm.isf(float(alpha) / 2.0))


def a_min(vif: float, alpha: float, n: float) -> float:
    """§P7-1(b): A_min = sqrt(VIF) * (z_alpha + z_0.80) * sqrt(2/N)."""
    n = float(n)
    if n <= 0:
        return float("inf")
    return float(math.sqrt(float(vif)) * (z_alpha(alpha) + Z_POWER_80)
                 * math.sqrt(2.0 / n))


def floor_constant(vif: float, alpha: float) -> float:
    """The `c` in A_min = c / sqrt(N). 32.3 at Tranche A under the measured VIF."""
    return float(math.sqrt(float(vif)) * (z_alpha(alpha) + Z_POWER_80)
                 * math.sqrt(2.0))


def measured_vif(path: str | None = None, default: float = MEASURED_VIF_DF2_PHASE):
    """Re-read the measured VIF from results_f4_58_vif.json; fall back to the constant.

    The constant is the authority for reproducibility (it is what every recorded
    number in the ledger was computed with); this reader exists so a re-run of
    F4-58 that moves the measurement is picked up by anything that asks for it.
    """
    p = path or os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), RESULTS_JSON)
    try:
        with open(p, "r", encoding="utf-8") as fh:
            d = json.load(fh)
        for s in d.get("sessions", []):
            if s.get("session") == "session_20260811T022953":
                v = s.get("overall_vif_df2_phase", {}).get("median")
                if v:
                    return float(v)
    except (OSError, json.JSONDecodeError, KeyError, TypeError):
        pass
    return float(default)


def cycles_in_window(period_days, record_days=EXPLORATION_RECORD_DAYS):
    """How many full cycles of `period_days` the analysis window contains."""
    p = float(period_days)
    if not (p > 0) or not np.isfinite(p):
        return float("nan")
    return float(record_days) / p


def max_identifiable_period(record_days=EXPLORATION_RECORD_DAYS,
                            min_cycles=MIN_CYCLES_IN_WINDOW):
    """The longest period the window can identify. 2,572 d over the v2 window."""
    return float(record_days) / float(min_cycles)


def window_report(period_days, record_days=EXPLORATION_RECORD_DAYS,
                  min_cycles=MIN_CYCLES_IN_WINDOW, name=None):
    """The S-15(c) line for one periodic feature. Verdict, not a p-value.

    Every field a report needs to print this without re-deriving anything, including
    the sentence that must accompany the verdict: a feature scored neither way is
    not a null and may not be counted as one.
    """
    cyc = cycles_in_window(period_days, record_days)
    bad = bool(np.isfinite(cyc) and cyc < float(min_cycles))
    return {
        "feature": name,
        "period_days": float(period_days),
        "record_days": float(record_days),
        "cycles_in_window": cyc,
        "min_cycles_required": float(min_cycles),
        "max_identifiable_period_days": max_identifiable_period(record_days,
                                                               min_cycles),
        "verdict": UNMEASURABLE_BY_WINDOW if bad else MEASURABLE_BY_WINDOW,
        "clause": WINDOW_CLAUSE,
        "clause_source": WINDOW_CLAUSE_SOURCE,
        "scored": not bad,
        "note": (("%.3g cycles in a %.0f d window is an IDENTIFIABILITY failure, "
                  "not a power failure: the fitted period is reporting record "
                  "length. Scored neither way -- this is NOT a null and may not be "
                  "counted as one, and no amount of N repairs it."
                  % (cyc, float(record_days))) if bad else
                 ("%.3g cycles in a %.0f d window: identifiable. The power floor "
                  "a_min() still applies and is a separate question."
                  % (cyc, float(record_days)))),
    }


def unmeasurable_by_window(period_days, record_days=EXPLORATION_RECORD_DAYS,
                           min_cycles=MIN_CYCLES_IN_WINDOW):
    """True when S-15(c) bites. The one-line form for a filter expression."""
    return window_report(period_days, record_days,
                         min_cycles)["verdict"] == UNMEASURABLE_BY_WINDOW


def assert_measurable_by_window(period_days, record_days=EXPLORATION_RECORD_DAYS,
                                min_cycles=MIN_CYCLES_IN_WINDOW, name=None):
    """Raise unless the window can identify this period. Returns the report.

    Wired here rather than in each caller so that FUTURE FEATURES INHERIT THE CLAUSE
    (§P7-14(c) asks for exactly that): anything that plants, scores or declares a
    periodic feature goes through this module for its floor already, and now gets
    the identifiability check on the same call path.
    """
    rep = window_report(period_days, record_days, min_cycles, name=name)
    if rep["verdict"] == UNMEASURABLE_BY_WINDOW:
        raise UnmeasurableByWindow(
            f"{name or 'feature'}: period {rep['period_days']:.4g} d gives "
            f"{rep['cycles_in_window']:.3g} cycles in a {rep['record_days']:.0f} d "
            f"window, below the S-15(c) minimum of {min_cycles:g} "
            f"(max identifiable period {rep['max_identifiable_period_days']:.0f} d). "
            f"This is an IDENTIFIABILITY limit and no N repairs it "
            f"(HYPOTHESIS_LEDGER.md §P7-14(c)). Score it neither way; do not plant "
            f"here and do not report a bound here.")
    return rep


def window_sweep(items, record_days=EXPLORATION_RECORD_DAYS,
                 min_cycles=MIN_CYCLES_IN_WINDOW):
    """S-15(c) over an iterable of (name, period_days). Returns the full table."""
    rows = [window_report(p, record_days, min_cycles, name=n) for n, p in items
            if p is not None and float(p) > 0]
    unm = [r for r in rows if r["verdict"] == UNMEASURABLE_BY_WINDOW]
    return {
        "clause": WINDOW_CLAUSE, "clause_source": WINDOW_CLAUSE_SOURCE,
        "record_days": float(record_days),
        "min_cycles_required": float(min_cycles),
        "cut_period_days": max_identifiable_period(record_days, min_cycles),
        "n_examined": len(rows), "n_unmeasurable_by_window": len(unm),
        "unmeasurable": sorted(unm, key=lambda r: -r["period_days"]),
        "rows": sorted(rows, key=lambda r: -r["period_days"]),
    }


# ---------------------------------- §P7-10(c)/§P7-13: the MARK-AXIS floor -----
# The count-path floor above is in RATE-MODULATION units and does not apply to a
# rank or circular-linear correlation. §P7-8(d) names that category error explicitly
# ("the mark axis needs its own floor, in its own units") and §P7-10(c) supplies it:
#
#     rho_min = sqrt(VIF_mark) * (z_alpha + z_0.80) / sqrt(n - 1)
#             = sqrt(VIF_mark) * 4.732 / sqrt(n - 1)      at alpha = 1.0e-4
#
# At n = 46,585 that is 0.0219 at VIF_mark = 1, 0.0439 at 4, 0.0658 at 9, 0.1074 at
# 24. §P7-10(c) is binding on F9-10: *"F9-10 declares its floor from the measured
# VIF_mark before it runs, or it does not run."*
#
# THE FALLBACK IS 4.575 AND THE REASON IS CENSORING, NOT CONSERVATISM (§P7-13(b)).
# Three of the 43 scored mark rows are censored at the Monte Carlo floor, and they
# are censored PRECISELY BECAUSE THEY ARE EXTREME -- their bounds 13.87 / 22.70 /
# 31.12 are 3x to 7x the median 4.345. Dropping them is not dropping data at random,
# it is dropping the upper tail because it is the upper tail, which biases a pooled
# median downward -- the anti-conservative direction for a fallback applied to a
# feature of unknown VIF. So the pooled fallback re-admits them at their bounds:
# 4.575, giving rho_min = 0.0469.
#
# AND THE DISTINCTION THAT MAKES THAT LEGITIMATE, carried into code because it is
# the kind of thing that gets lost: re-admitting a bound is legitimate for a
# POPULATION summary and illegitimate for a PER-FEATURE floor. The three censored
# features remain UNMEASURABLE BY DECLARATION (§P7-11(c)) -- their own VIF is still a
# bound, and a bound cannot produce a floor for the feature it belongs to.
ALPHA_TRANCHE_B = 1.0e-4
VIF_MARK_FALLBACK = 4.575
VIF_MARK_FALLBACK_SOURCE = (
    "HYPOTHESIS_LEDGER.md §P7-13(b) -- pooled mark-axis fallback, bounds-readmitted "
    "(the measurements-only median is 4.345 and is retained beside it)")
VIF_MARK_MEASUREMENTS_ONLY_MEDIAN = 4.345
MARK_Z_SUM = 4.732          # z_{alpha=1e-4, two-sided} + z_0.80 = 3.891 + 0.842


class MarkVifIsABound(AssertionError):
    """A per-feature mark floor was asked for from a CENSORED (bound) VIF."""


def rho_min(n_events, vif_mark=None, alpha=ALPHA_TRANCHE_B):
    """§P7-10(c): the smallest mark-axis correlation declarable at 80% power."""
    n = float(n_events)
    if n <= 1:
        return float("inf")
    v = VIF_MARK_FALLBACK if vif_mark is None else float(vif_mark)
    z = MARK_Z_SUM if abs(float(alpha) - ALPHA_TRANCHE_B) < 1e-15 else (
        z_alpha(alpha) + Z_POWER_80)
    return float(math.sqrt(v) * z / math.sqrt(n - 1.0))


def mark_floor_report(n_events, vif_mark=None, vif_is_bound=False, feature=None,
                      alpha=ALPHA_TRANCHE_B):
    """Everything §P7-10(c) requires printed next to a mark-axis result.

    `vif_is_bound=True` marks a VIF that saturated at the Monte Carlo floor. Per
    §P7-13(a) such a value is an UPPER BOUND, over-stated, and per §P7-11(c) it may
    not produce a per-feature floor -- so the report says UNMEASURABLE BY DECLARATION
    and `assert_mark_floor_declarable` refuses.
    """
    used_fallback = vif_mark is None
    v = VIF_MARK_FALLBACK if used_fallback else float(vif_mark)
    fl = rho_min(n_events, v, alpha)
    return {
        "feature": feature,
        "n_events": float(n_events),
        "alpha": float(alpha),
        "vif_mark": v,
        "vif_mark_is_bound": bool(vif_is_bound),
        "vif_mark_source": (VIF_MARK_FALLBACK_SOURCE if used_fallback
                            else "measured per-feature (F4-58M)"),
        "vif_mark_pooled_fallback": VIF_MARK_FALLBACK,
        "vif_mark_measurements_only_median": VIF_MARK_MEASUREMENTS_ONLY_MEDIAN,
        "rho_min": fl,
        "declarable": not bool(vif_is_bound),
        "verdict": ("UNMEASURABLE BY DECLARATION (§P7-11(c): this feature's own "
                    "VIF_mark is a censored BOUND, and a bound cannot produce a "
                    "floor for the feature it belongs to)" if vif_is_bound
                    else "floor declared"),
        "formula": "rho_min = sqrt(VIF_mark) * 4.732 / sqrt(n - 1)  (§P7-10(c))",
        "population_vs_per_feature": (
            "§P7-13(b): re-admitting censored rows at their bounds is legitimate for "
            "a POPULATION summary (the pooled fallback) and illegitimate for a "
            "PER-FEATURE floor. Two different uses of the same number."),
    }


def assert_mark_floor_declarable(n_events, vif_mark=None, vif_is_bound=False,
                                 feature=None, alpha=ALPHA_TRANCHE_B):
    """§P7-10(c): *declare the floor from the measured VIF_mark, or do not run.*"""
    rep = mark_floor_report(n_events, vif_mark, vif_is_bound, feature, alpha)
    if not rep["declarable"]:
        raise MarkVifIsABound(
            f"{feature or 'mark test'}: VIF_mark is a censored BOUND "
            f"({rep['vif_mark']:.4g}), so no per-feature mark floor may be derived "
            f"from it (§P7-11(c), §P7-13(a)). The feature stays UNMEASURABLE BY "
            f"DECLARATION. Use the pooled fallback only for a POPULATION summary, "
            f"never as this feature's own floor.")
    return rep


def assert_mark_plant_above_floor(rho, n_events, vif_mark=None,
                                  factor=PLANT_FACTOR, feature=None,
                                  alpha=ALPHA_TRANCHE_B):
    """The mark-axis analogue of `assert_plant_above_floor`, in correlation units."""
    rep = mark_floor_report(n_events, vif_mark, False, feature, alpha)
    need = float(factor) * rep["rho_min"]
    rep.update({"planted_rho": float(rho), "required_plant_factor": float(factor),
                "required_min_plant": need,
                "rho_over_floor": float(rho) / max(rep["rho_min"], 1e-300),
                "compliant": bool(abs(float(rho)) >= need)})
    if not rep["compliant"]:
        raise PlantBelowFloor(
            f"{feature or 'mark plant'}: planted rho {float(rho):.4f} is below "
            f"{factor:g}x the §P7-10(c) mark floor at n = {float(n_events):.0f} "
            f"(rho_min = {rep['rho_min']:.4f}, required >= {need:.4f}; "
            f"VIF_mark = {rep['vif_mark']:.4g}, alpha = {float(alpha):.3e}). "
            f"Raise the plant; do not lower the floor.")
    return rep


# ------------------------------------------------------- the planting contract --
def min_plant_amplitude(n, alpha: float = ALPHA_TRANCHE_A,
                        vif: float | None = None,
                        factor: float = PLANT_FACTOR) -> float:
    """§P7-8(d): the smallest amplitude a harness is allowed to plant at this N.

    `factor * a_min(...)`. Plant at or above this, or the test is measuring its own
    power and calling the answer an instrument verdict.
    """
    v = MEASURED_VIF_DF2_PHASE if vif is None else float(vif)
    return float(factor) * a_min(v, alpha, n)


def plant_report(amplitude, n, alpha: float = ALPHA_TRANCHE_A,
                 vif: float | None = None, factor: float = PLANT_FACTOR,
                 what: str = "plant") -> dict:
    """Everything a harness must state next to a planted amplitude (§P7-8(d))."""
    v = MEASURED_VIF_DF2_PHASE if vif is None else float(vif)
    fl = a_min(v, alpha, n)
    need = float(factor) * fl
    return {
        "what": what,
        "planted_amplitude": float(amplitude),
        "N": float(n),
        "alpha": float(alpha),
        "vif": v,
        "vif_source": MEASURED_VIF_SOURCE,
        "vif_status": MEASURED_VIF_STATUS,
        "operative_floor_A_min": fl,
        "required_plant_factor": float(factor),
        "required_min_plant": need,
        "margin_over_required": float(amplitude) - need,
        "amplitude_over_floor": (float(amplitude) / fl) if fl > 0 else float("inf"),
        "compliant": bool(float(amplitude) >= need),
        "rule": ("HYPOTHESIS_LEDGER.md §P7-8(d): plants must sit at >= 2x the "
                 "operative floor at the relevant N, so that a failure to recover "
                 "is an INSTRUMENT verdict and not a POWER verdict wearing an "
                 "instrument's clothes."),
    }


def assert_plant_above_floor(amplitude, n, alpha: float = ALPHA_TRANCHE_A,
                             vif: float | None = None,
                             factor: float = PLANT_FACTOR,
                             what: str = "plant") -> dict:
    """Raise unless `amplitude` clears `factor` x the floor at `n`. Returns the report."""
    rep = plant_report(amplitude, n, alpha, vif, factor, what)
    if not rep["compliant"]:
        raise PlantBelowFloor(
            f"{what}: planted amplitude {rep['planted_amplitude']:.4f} is below "
            f"{factor:g}x the operative S-15 floor at N = {rep['N']:.0f} "
            f"(A_min = {rep['operative_floor_A_min']:.4f}, required "
            f">= {rep['required_min_plant']:.4f}; VIF = {rep['vif']:.4f}, "
            f"alpha = {rep['alpha']:.3e}). "
            f"A plant below the floor fails for POWER reasons while reading as an "
            f"INSTRUMENT failure -- HYPOTHESIS_LEDGER.md §P7-8(d) exists to prevent "
            f"that false negative being recorded against the pipeline. Raise the "
            f"plant; do not lower the floor."
        )
    return rep
