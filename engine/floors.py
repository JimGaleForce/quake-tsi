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
