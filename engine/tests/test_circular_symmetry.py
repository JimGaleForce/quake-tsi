"""Acceptance tests for engine/circular_symmetry.py.

The load-bearing one is `test_ten_left_ten_right_is_not_a_null`, which encodes the
exact scenario the module exists for: a totally effective symmetric mechanism that
every first-moment circular test calls uniform. If that test ever fails, the program
has gone back to reporting a null on R1 as a null on the hypothesis.
"""

from __future__ import annotations

import numpy as np
import pytest

from engine import circular_symmetry as CS


# ------------------------------------------------------- THE ONE THAT MATTERS --
def test_ten_left_ten_right_is_not_a_null():
    """20 of 20 fired on horizontal pull. R1 says uniform; R2 says total."""
    angles = np.concatenate([np.zeros(10), np.full(10, np.pi)])
    sp = CS.harmonic_spectrum(angles, max_order=4)
    assert sp["R"][1] == pytest.approx(0.0, abs=1e-12), "first moment cancels"
    assert sp["R"][2] == pytest.approx(1.0, abs=1e-12), "second moment is total"
    dom = CS.dominant_order(sp)
    assert dom["order"] == 2
    assert "AXIAL" in dom["meaning"]
    # and the analytic uniform test would have called R1 non-significant
    assert CS.rayleigh_p_uniform(sp["R"][1], 20) > 0.9


def test_the_same_thing_with_noise_still_reads_axial():
    """Threshold taken from theory, not guessed.

    For a wrapped normal of width sigma the m-th resultant is exp(-m^2 sigma^2 / 2).
    At sigma = 0.35 that is 0.783 for m = 2, so the first version of this test, which
    asserted R2 > 0.85, was demanding more than the construction can produce. The
    expectation was wrong, not the code.
    """
    rng = np.random.default_rng(0)
    n, sigma = 400, 0.35
    lobe = rng.normal(0.0, sigma, n // 2)
    angles = np.concatenate([lobe, lobe + np.pi])
    sp = CS.harmonic_spectrum(angles, max_order=4)
    expected_r2 = np.exp(-(2.0 ** 2) * sigma ** 2 / 2.0)
    assert sp["R"][2] == pytest.approx(expected_r2, abs=0.06), (sp["R"][2], expected_r2)
    assert sp["R"][1] < 0.15
    assert CS.dominant_order(sp)["order"] == 2


def test_a_genuinely_directional_effect_reads_order_one():
    """The module must not simply always say 'axial'."""
    rng = np.random.default_rng(1)
    angles = rng.normal(1.0, 0.4, 500)
    sp = CS.harmonic_spectrum(angles, max_order=4)
    assert CS.dominant_order(sp)["order"] == 1
    assert sp["R"][1] > sp["R"][2]


def test_four_lobed_shear_geometry_reads_order_four():
    """Conjugate shear planes at 90 degree spacing."""
    base = np.array([0.0, np.pi / 2, np.pi, 3 * np.pi / 2])
    angles = np.repeat(base, 50)
    sp = CS.harmonic_spectrum(angles, max_order=6)
    assert sp["R"][4] == pytest.approx(1.0, abs=1e-12)
    assert sp["R"][1] == pytest.approx(0.0, abs=1e-12)
    assert sp["R"][2] == pytest.approx(0.0, abs=1e-12)
    assert CS.dominant_order(sp)["order"] == 4


# ------------------------------------------------------------------ mechanics --
def test_uniform_angles_are_flat_at_every_order():
    a = np.linspace(0.0, 2.0 * np.pi, 200001)[:-1]
    sp = CS.harmonic_spectrum(a, max_order=6)
    for m in range(1, 7):
        assert sp["R"][m] < 1e-3, (m, sp["R"][m])


def test_rayleigh_matches_the_known_uniform_distribution():
    """Cross-check of the analytic form on the one case where it is valid."""
    rng = np.random.default_rng(7)
    ps = []
    for _ in range(400):
        a = rng.uniform(0.0, 2.0 * np.pi, 60)
        ps.append(CS.rayleigh_p_uniform(CS.harmonic_resultant(a, 1), 60))
    ps = np.array(ps)
    assert 0.02 < np.mean(ps < 0.05) < 0.10, np.mean(ps < 0.05)


def test_fold_maps_opposite_directions_together_at_order_two():
    a = np.array([0.1, 0.1 + np.pi, 2.0, 2.0 + np.pi])
    f = CS.fold(a, 2)
    assert f[0] == pytest.approx(f[1])
    assert f[2] == pytest.approx(f[3])


def test_fold_at_order_four_also_merges_forward_and_back():
    a = np.array([0.3, 0.3 + np.pi / 2, 0.3 + np.pi, 0.3 + 3 * np.pi / 2])
    f = CS.fold(a, 4)
    assert np.allclose(f, f[0])


def test_harmonic_phase_recovers_the_preferred_axis():
    axis = 0.7
    angles = np.concatenate([np.full(100, axis), np.full(100, axis + np.pi)])
    got = CS.harmonic_phase(angles, 2)
    assert got == pytest.approx(axis, abs=1e-9)
    assert 0.0 <= got < np.pi


def test_relative_angle_folds_against_a_reference():
    theta = np.array([0.0, np.pi])
    strike = np.full(2, 0.3)
    r = CS.relative_angle(theta, strike, order=2)
    assert r[0] == pytest.approx(r[1])


def test_parity_split_is_actually_a_split():
    """The recorded bug: `odd` was written sign(v)*v, which is |v|."""
    v = np.array([-3.0, -1.0, 0.0, 2.0, 5.0])
    p = CS.fold_parity(v)
    assert np.allclose(p["even"], np.abs(v))
    assert np.allclose(p["odd"], v)
    assert not np.allclose(p["even"], p["odd"])
    assert np.allclose(p["sign"], np.sign(v))


def test_parity_even_coordinate_is_symmetric_and_odd_is_not():
    v = np.array([-2.0, -1.0, 1.0, 2.0])
    p = CS.fold_parity(v)
    pm = CS.fold_parity(-v)
    assert np.allclose(p["even"], pm["even"])          # even survives reflection
    assert np.allclose(p["odd"], -pm["odd"])           # odd flips


def test_empty_and_nan_input_do_not_crash():
    assert np.isnan(CS.harmonic_resultant(np.array([]), 2))
    a = np.array([0.0, np.nan, np.pi])
    assert CS.harmonic_resultant(a, 2) == pytest.approx(1.0, abs=1e-12)


def test_dominant_order_reports_the_LOWEST_order_not_the_largest_R():
    """An order-2 symmetry is automatically order-4, and the fundamental one wins.

    The exact ten-left/ten-right set gives R2 = R4 = 1.0. Taking the argmax would name
    a four-lobed shear mechanism where the data shows a two-lobed axial one; this is
    the bug the first version of dominant_order had.
    """
    angles = np.concatenate([np.zeros(10), np.full(10, np.pi)])
    sp = CS.harmonic_spectrum(angles, max_order=6)
    assert sp["R"][2] == pytest.approx(1.0, abs=1e-12)
    assert sp["R"][4] == pytest.approx(1.0, abs=1e-12)
    assert sp["R"][6] == pytest.approx(1.0, abs=1e-12)
    dom = CS.dominant_order(sp)
    assert dom["order"] == 2, "the fundamental order, not a multiple of it"
    assert dom["also_high"] == [4, 6]
    assert dom["expected_also_high"] == [4, 6], "multiples of 2 within the range"
    assert dom["margin_over_next_independent"] > 0.5
    assert "waveform-matched null" in dom["caution"]


def test_a_true_order_four_pattern_is_not_demoted_to_order_two():
    """The lowest-order rule must not simply always answer 2."""
    base = np.array([0.0, np.pi / 2, np.pi, 3 * np.pi / 2])
    sp = CS.harmonic_spectrum(np.repeat(base, 50), max_order=6)
    dom = CS.dominant_order(sp)
    assert dom["order"] == 4
    assert sp["R"][2] == pytest.approx(0.0, abs=1e-12)


# ------------------------------- world vs regions: the pooling trap, one level up --
def test_two_regions_with_perpendicular_axes_pool_to_exactly_zero():
    """THE SECOND TRAP, and it is the ten-left/ten-right problem one level up.

    Region A concentrates on the 0 axis, region B on the 90 axis. Each has R2 = 1.0,
    i.e. each is a TOTAL effect. Pool the raw angles and R2 = 0.0 exactly, because the
    two axes are antipodal on the doubled circle. A world-level test reports nothing
    while every region individually is maximal.
    """
    a = np.concatenate([np.zeros(50), np.full(50, np.pi)])
    b = np.concatenate([np.full(50, np.pi / 2), np.full(50, 3 * np.pi / 2)])
    assert CS.harmonic_resultant(a, 2) == pytest.approx(1.0, abs=1e-12)
    assert CS.harmonic_resultant(b, 2) == pytest.approx(1.0, abs=1e-12)
    pooled = CS.harmonic_resultant(np.concatenate([a, b]), 2)
    assert pooled == pytest.approx(0.0, abs=1e-12), "the pooled test sees nothing"

    per = CS.per_region_spectra({"A": a, "B": b}, max_order=4)
    assert per["A"]["R"][2] == pytest.approx(1.0, abs=1e-12)
    assert per["B"]["R"][2] == pytest.approx(1.0, abs=1e-12)
    ag = CS.axis_agreement(per, order=2)
    assert ag["agreement_R"] == pytest.approx(0.0, abs=1e-9), \
        "the regions are real and they disagree"


def test_regions_that_agree_show_high_axis_agreement():
    a = np.concatenate([np.zeros(50), np.full(50, np.pi)])
    b = np.concatenate([np.full(50, 0.05), np.full(50, 0.05 + np.pi)])
    per = CS.per_region_spectra({"A": a, "B": b}, max_order=4)
    ag = CS.axis_agreement(per, order=2)
    assert ag["agreement_R"] > 0.99
    assert abs(ag["common_axis_rad"] - 0.025) < 0.05


def test_cochran_q_flags_disagreeing_regions():
    z = {"A": 4.0, "B": -4.0, "C": 3.5, "D": -3.5}
    n = {k: 100 for k in z}
    rep = CS.combine_regions(z, n)
    assert abs(rep["z_pooled"]) < 1.0, "pooled looks like nothing"
    assert rep["heterogeneous"] is True, "but Q says the regions disagree"
    assert rep["I2"] > 0.8


def test_cochran_q_is_quiet_when_regions_agree():
    z = {"A": 2.0, "B": 2.2, "C": 1.8, "D": 2.1}
    n = {k: 100 for k in z}
    rep = CS.combine_regions(z, n)
    assert rep["z_pooled"] > 3.0
    assert rep["heterogeneous"] is False
    assert rep["I2"] < 0.5


def test_combine_regions_weights_by_sqrt_n():
    """Combining agreeing regions must RAISE the pooled z, not preserve it.

    The first version of this test asserted z_pooled == 2.0, which misunderstands
    Stouffer combination: two independent z = 2 results are stronger evidence than
    one. With w = sqrt(n) the answer is 220/sqrt(10100) = 2.189, and the big region
    dominates as it should.
    """
    rep = CS.combine_regions({"big": 2.0, "small": 2.0},
                             {"big": 10000, "small": 100})
    assert rep["z_pooled"] == pytest.approx(2.189, abs=1e-3)
    assert rep["z_pooled"] > 2.0
    assert rep["n_regions"] == 2
    # a single region reproduces its own z exactly
    solo = CS.combine_regions({"big": 2.0}, {"big": 10000})
    assert solo["z_pooled"] == pytest.approx(2.0, rel=1e-12)


def test_combine_regions_survives_empty_and_bad_input():
    assert CS.combine_regions({}, {})["n_regions"] == 0
    rep = CS.combine_regions({"A": float("nan"), "B": 2.0}, {"A": 10, "B": 10})
    assert rep["n_regions"] == 1
