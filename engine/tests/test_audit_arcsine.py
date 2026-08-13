"""D-1: the arcsine / level-vs-phase control. §P7-22 Ratification 1's mandatory arm.

The exact arithmetic Popper verified before ruling -- *"uniform phase puts
(2pi/3)/(2pi) = 1/3 of events in the lowest quarter of the tidal range, simulated
0.33332 at 1e7 draws"* -- reproduced on demand, plus the level-vs-phase comparison that
lets D-1 run later as a declared arm.
"""

from __future__ import annotations

import numpy as np
import pytest

from engine import audit_arcsine as A

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def test_lowest_quarter_fraction_is_exactly_one_third():
    """The closed form, exactly: 1 - arccos(-1/2)/pi = 1 - (2pi/3)/pi = 1/3."""
    assert A.lowest_fraction(0.25) == pytest.approx(1.0 / 3.0, abs=1e-15)
    assert A.LOWEST_QUARTER_FRACTION_EXACT == pytest.approx(1.0 / 3.0, abs=1e-15)


def test_the_median_split_is_unbiased_and_the_quartile_split_is_not():
    """Why the seed's quartile reading and a median reading disagree."""
    assert A.lowest_fraction(0.5) == pytest.approx(0.5, abs=1e-15)
    assert A.lowest_fraction(0.25) > 0.25
    assert A.lowest_fraction(0.75) == pytest.approx(2.0 / 3.0, abs=1e-15)


def test_simulation_reproduces_one_third_to_monte_carlo_precision():
    """§P7-22's own check, run small: uniform phase -> 1/3 in the lowest quarter."""
    r = A.arcsine_control(n=400000, seed=1)
    assert abs(r["z"]) < 4.0, r
    assert r["measured_fraction"] == pytest.approx(1.0 / 3.0, abs=0.004)


def test_level_reading_looks_like_clustering_while_phase_reading_is_flat():
    """THE demonstration: identical data, two readings, only one of them is a signal."""
    r = A.level_vs_phase_report(n=40000, seed=2)
    lev = r["level_reading"]
    ph = r["phase_reading"]
    assert lev["measured"] == pytest.approx(1.0 / 3.0, abs=0.01)
    # a third in the bottom quarter is a 33% excess over the 25% an eye expects
    assert lev["measured"] > 1.25 * 0.25
    # and the same events carry no phase concentration at all
    assert ph["R1"] < 0.02
    assert ph["V_star"] < 2.5           # ~1.6-1.75 is the uniform Kuiper scale


def test_the_artifact_survives_the_real_multi_constituent_waveform():
    """A pure cosine invites 'the real tide is not a sinusoid'. It is not, and it is worse.

    The degree-2 potential goes as P2(cos z), whose range is [-1/2, +1]: the scalar is
    SKEWED, so the level rendering puts even more than a third of a null catalogue in
    the lowest quarter of range. The measurement is what is asserted here; the number
    is not rounded into the docstring of the module under test.
    """
    r = A.real_waveform_control(n_days=200.0, sample_minutes=10.0)
    assert r["measured_fraction_lowest_quarter"] > 1.0 / 3.0
    assert r["measured_fraction_lowest_quarter"] < 0.9


def test_the_conclusion_does_not_depend_on_the_site():
    """Two sites, same conclusion -- which is what makes this not a look at a region."""
    a = A.real_waveform_control(lat_deg=55.34, lon_deg=-160.50, n_days=200.0,
                                sample_minutes=10.0)
    b = A.real_waveform_control(lat_deg=-33.0, lon_deg=-71.6, n_days=200.0,
                                sample_minutes=10.0)
    assert a["measured_fraction_lowest_quarter"] > 1.0 / 3.0
    assert b["measured_fraction_lowest_quarter"] > 1.0 / 3.0


def test_the_arm_record_states_that_it_may_end_the_entry():
    r = A.level_vs_phase_report(n=2000, seed=0)
    assert "MOTIVATING OBSERVATION IS EXPLAINED" in r["verdict_rule"]


def test_fraction_uses_the_observed_range_not_the_theoretical_one():
    """The statistic must be computable from a plot, which is what the seed was."""
    lev = np.array([-2.0, -1.9, 0.0, 1.0, 2.0])
    # lowest quarter of [-2, 2] is < -1: two of five
    assert A.fraction_in_lowest_quarter(lev, 0.25) == pytest.approx(0.4)


# ------------------------- §P7-23(C): the declared trough-vs-mid-slope readout --
def test_the_three_declared_bands_are_exactly_equiprobable_under_the_null():
    """Why THESE bands: the null is a flat multinomial and needs no simulation."""
    assert A.lowest_fraction(0.25) == pytest.approx(1.0 / 3.0, abs=1e-15)
    assert 1.0 - A.lowest_fraction(0.75) == pytest.approx(1.0 / 3.0, abs=1e-15)
    rng = np.random.default_rng(7)
    c = A.classify_level_set(np.sin(rng.random(400000) * 2 * np.pi))
    for band in A.BAND_NAMES:
        assert c["fractions"][band] == pytest.approx(1.0 / 3.0, abs=0.005), band
    # and at a seed-set-sized n a null draw is not classified either way
    n_unclassified = sum(
        A.classify_level_set(np.sin(rng.random(300) * 2 * np.pi))["verdict"]
        == "UNCLASSIFIED" for _ in range(40))
    assert n_unclassified >= 34, n_unclassified      # ~2-sided z=2 -> ~95% each way


def test_the_arcsine_density_numbers_the_ruling_verified():
    """§P7-23(C): 0.318 at mid-level, 1.019 at |x| = 0.95."""
    assert float(A.arcsine_density(0.0)) == pytest.approx(0.3183, abs=1e-3)
    assert float(A.arcsine_density(0.95)) == pytest.approx(1.019, abs=1e-3)


def test_trough_and_mid_slope_sets_classify_in_opposite_directions():
    """The two-sided diagnostic: the artifact and the physics separate cleanly."""
    rng = np.random.default_rng(8)
    trough = -0.5 - 0.5 * rng.random(400)          # all in the trough band
    mid = rng.uniform(-0.4, 0.4, 400)              # all mid-slope
    assert A.classify_level_set(trough)["verdict"].startswith("TROUGH-CONCENTRATED")
    assert "artifact-CONSISTENT" in A.classify_level_set(trough)["verdict"]
    assert A.classify_level_set(mid)["verdict"].startswith("MID-SLOPE-CONCENTRATED")
    assert "artifact-INCONSISTENT" in A.classify_level_set(mid)["verdict"]


def test_the_seed_readout_reports_both_percentile_and_phase_per_event():
    """§P7-23(C) requires BOTH, per event, and the set-level classification."""
    ph = np.array([0.2, 3.5, 4.0, 4.4, 5.0])
    lev = np.sin(ph)
    r = A.seed_readout(lev, ph)
    assert len(r["per_event"]) == 5
    for e in r["per_event"]:
        assert 0.0 <= e["level_percentile"] <= 1.0
        assert 0.0 <= e["phase_rad"] < 2 * np.pi
        assert e["phase_deg"] == pytest.approx(np.degrees(e["phase_rad"]))
    assert "verdict" in r["classification"]
    assert "BEFORE any event outside the seed set is scored" in r["timing_rule"]


def test_below_neutral_and_falling_is_exactly_a_quarter_cycle():
    """§P7-23(D): the quadrant is rotation-free once the scalar is fixed."""
    assert A.QUADRANT_NULL_FRACTION == 0.25
    rng = np.random.default_rng(9)
    ph = rng.random(400000) * 2 * np.pi
    r = A.seed_readout(np.sin(ph), ph)["quadrant_below_neutral_and_falling"]
    assert r["fraction"] == pytest.approx(0.25, abs=0.005)
    assert "provenance question" in r["note"]
