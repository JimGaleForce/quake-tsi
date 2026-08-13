"""engine/floors.py -- the S-15 declaration floor and the §P7-8(d) planting contract.

The floor exists so that a harness which plants a signal cannot accidentally test
its own power and record the answer as an instrument verdict.
"""

from __future__ import annotations

import math

import pytest
from scipy import stats

from engine import floors, regions as R


def test_formula_matches_the_ledgers_own_arithmetic():
    """§P7-1(a): alpha = 0.1/259 -> z_alpha = 3.549, two-sided."""
    assert floors.z_alpha(0.10 / 259) == pytest.approx(3.549, abs=5e-4)
    assert floors.Z_POWER_80 == pytest.approx(float(stats.norm.ppf(0.80)))


def test_floor_constant_reproduces_the_P7_8d_table():
    """§P7-8(d): 32.3/sqrt(N) at Tranche A, i.e. 14.9% at N = 46,585."""
    c = floors.floor_constant(floors.MEASURED_VIF_DF2_PHASE, floors.ALPHA_TRANCHE_A)
    assert c == pytest.approx(32.3, abs=0.05)
    am = floors.a_min(floors.MEASURED_VIF_DF2_PHASE, floors.ALPHA_TRANCHE_A, 46585)
    assert 100.0 * am == pytest.approx(14.9, abs=0.1)
    # the 2026-08-11 session's own row: m = 259 -> 30.5/sqrt(N) -> 14.1%
    c259 = floors.floor_constant(floors.MEASURED_VIF_DF2_PHASE, floors.ALPHA_20260811)
    assert c259 == pytest.approx(30.5, abs=0.05)


def test_measured_vif_is_the_measurement_not_the_falsified_guess():
    assert floors.MEASURED_VIF_DF2_PHASE == pytest.approx(24.08, abs=0.01)
    assert floors.measured_vif() == pytest.approx(floors.MEASURED_VIF_DF2_PHASE,
                                                  rel=1e-9)
    # explicitly NOT Popper's falsified 3.94 nor Kepler's 9.7
    assert abs(floors.MEASURED_VIF_DF2_PHASE - 3.94) > 1.0
    assert abs(floors.MEASURED_VIF_DF2_PHASE - 9.7) > 1.0


def test_measured_vif_falls_back_when_the_results_file_is_missing(tmp_path):
    assert floors.measured_vif(path=str(tmp_path / "nope.json")) == pytest.approx(
        floors.MEASURED_VIF_DF2_PHASE)


def test_floor_scales_as_one_over_sqrt_N():
    a1 = floors.a_min(24.0, 1e-4, 10_000)
    a4 = floors.a_min(24.0, 1e-4, 40_000)
    assert a1 / a4 == pytest.approx(2.0, rel=1e-12)


def test_floor_scales_as_sqrt_VIF():
    assert (floors.a_min(4.0, 1e-4, 1000) / floors.a_min(1.0, 1e-4, 1000)
            == pytest.approx(2.0, rel=1e-12))


def test_regions_a_min_is_the_same_implementation():
    """One formula, one implementation -- regions.a_min delegates."""
    for vif, alpha, n in ((24.08, 0.10 / 713, 32000), (1.0, 0.05, 500),
                          (9.7, 1e-5, 1_000_000)):
        assert R.a_min(vif, alpha, n) == floors.a_min(vif, alpha, n)


def test_zero_or_negative_N_is_infinite_not_a_crash():
    assert math.isinf(floors.a_min(24.0, 1e-4, 0))
    assert math.isinf(floors.a_min(24.0, 1e-4, -5))


# ---------------------------------------------- §P7-8(d): the planting contract --
def test_plant_below_two_times_the_floor_raises():
    n = 32_000.0
    fl = floors.a_min(floors.MEASURED_VIF_DF2_PHASE, floors.ALPHA_TRANCHE_A, n)
    with pytest.raises(floors.PlantBelowFloor):
        floors.assert_plant_above_floor(1.999 * fl, n)
    ok = floors.assert_plant_above_floor(2.001 * fl, n)
    assert ok["compliant"]


def test_the_pre_P7_8d_plant_of_0_35_would_have_been_refused():
    """The concrete case §P7-8(d) was written about: 0.35 at one region's N."""
    with pytest.raises(floors.PlantBelowFloor):
        floors.assert_plant_above_floor(0.35, 32_000.0)
    # and the raised plant passes
    assert floors.assert_plant_above_floor(0.40, 32_000.0)["compliant"]


def test_the_refusal_explains_power_versus_instrument():
    """The message is the point: a reader must not record this as a pipeline bug."""
    with pytest.raises(floors.PlantBelowFloor) as e:
        floors.assert_plant_above_floor(0.01, 32_000.0, what="G-M1 arm (ii)")
    msg = str(e.value)
    assert "POWER" in msg and "INSTRUMENT" in msg
    assert "P7-8(d)" in msg
    assert "G-M1 arm (ii)" in msg
    assert "Raise the plant; do not lower the floor" in msg


def test_plant_report_carries_its_conventions():
    """S-18: a number may not enter an arithmetic unless it carries its convention."""
    rep = floors.plant_report(0.40, 32_000.0)
    for k in ("alpha", "vif", "vif_source", "vif_status", "N",
              "operative_floor_A_min", "required_min_plant", "rule"):
        assert k in rep and rep[k] not in (None, "")
    assert "PROVISIONAL" in rep["vif_status"]
    assert rep["amplitude_over_floor"] == pytest.approx(
        0.40 / rep["operative_floor_A_min"])


def test_smaller_N_means_a_higher_floor_so_per_region_binds():
    """The binding aggregation is the smallest N, which is why per-region binds."""
    per_region = floors.min_plant_amplitude(32_000.0)
    globally = floors.min_plant_amplitude(6 * 32_000.0)
    assert per_region > globally
    assert per_region / globally == pytest.approx(math.sqrt(6.0), rel=1e-12)
