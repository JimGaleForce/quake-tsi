"""F9-01 and F9-04: the statistics, their nulls, and the two controls §P7-3 demands.

The counted invariants here are the S-17 recovery demands themselves, run small:

  * F9-01 recovers a planted ANTIPODAL two-lobed structure whose first moment is ~0;
  * F9-01 does NOT fire on a planted PURE SINUSOID -- §P7-3(1), the negative control
    Kepler did not name, and the one without which a second-moment "detection" is
    unfalsifiably confounded with the first moment;
  * F9-04 recovers a CONCENTRATED-PHASE plant that the Rayleigh-form 2-df statistic
    provably misses -- the money demonstration;
  * F9-04 does NOT fire on the day-binning lattice -- §P7-3(2)'s negative control.
"""

from __future__ import annotations

import numpy as np
import pytest

from engine import circstat as C, floors, recovery_b as R

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


# --------------------------------------------------------------- moments -----
def test_second_moment_sees_antipodal_that_first_moment_cannot():
    """The shape §K87-0(c) says no v1 statistic could have produced."""
    th = np.linspace(0, 2 * np.pi, 20000, endpoint=False)
    w = 1.0 + 0.6 * np.cos(2 * th)              # two lobes, first moment exactly 0
    r1, _ = C.circular_moment(th, w, k=1)
    r2, _ = C.circular_moment(th, w, k=2)
    assert r1 < 1e-3, r1
    assert r2 > 0.25, r2
    rep = C.moments_report(th, w)
    assert rep["reading"].startswith("AXIAL")


def test_first_moment_plant_reads_as_harmonic_ambiguous():
    th = np.linspace(0, 2 * np.pi, 20000, endpoint=False)
    rep = C.moments_report(th, 1.0 + 0.6 * np.cos(th))
    assert rep["R1"] > 0.25 and rep["R2"] < 1e-3
    assert rep["reading"].startswith("HARMONIC-OR-AMBIGUOUS")


def test_doubled_angle_design_is_orthogonal_to_the_fundamental():
    """Why the sinusoid negative control holds, checked rather than asserted."""
    th = np.linspace(0, 2 * np.pi, 4096, endpoint=False)
    X1 = np.column_stack([np.sin(th), np.cos(th)])
    X2 = C.doubled_angle_design(th)
    assert np.abs(X1.T @ X2).max() < 1e-8


# ------------------------------------------------- Kuiper / Watson mechanics --
def test_kuiper_watson_zero_when_observed_matches_expected():
    th = np.linspace(0, 2 * np.pi, 500, endpoint=False)
    lam = 3.0 + np.cos(th)
    out = C.kuiper_watson(th, lam, lam)
    assert out["V"] == pytest.approx(0.0, abs=1e-12)
    assert out["U2"] == pytest.approx(0.0, abs=1e-12)


def test_wrap_phase_closes_the_seam_and_merges_float_noise_ties():
    """`np.mod(2*pi*t, 2*pi)` lands on BOTH sides of the seam; that is the hazard."""
    noisy = np.mod(2 * np.pi * np.arange(3000, dtype=float), 2 * np.pi)
    assert np.unique(noisy).size > 1                  # they are NOT exactly equal
    assert noisy.max() > 6.28                         # and some are near 2pi, not 0
    # grouping WITHOUT closing the seam sees two phases where there is one
    assert C.phase_group_ends(np.sort(noisy)).size == 2
    assert C.phase_group_ends(np.sort(C.wrap_phase(noisy))).size == 1


def test_day_lattice_phase_has_exactly_one_distinct_value():
    for P in (1.0, 0.5, 0.25):
        th = C.day_lattice_phase(2000, period_days=P)
        assert C.phase_group_ends(np.sort(th)).size == 1, P


def test_day_lattice_negative_control_cannot_fire():
    """§P7-3(2): the statistic must not fire on a construction that is only the grid."""
    rng = np.random.default_rng(3)
    n = 2000
    counts = rng.poisson(8.0, n).astype(float)
    # a strong TREND, which is exactly what a tie-blind implementation would report
    counts = counts * np.linspace(0.5, 1.5, n)
    offset = np.full(n, 8.0)
    th = C.day_lattice_phase(n, 1.0)
    rows = C.omnibus_test(th, counts, offset, 200, rng, block_days=30.0,
                          periodic=True)
    assert rows[0]["statistic"] == pytest.approx(0.0, abs=1e-12)
    assert rows[1]["statistic"] == pytest.approx(0.0, abs=1e-12)
    assert min(r["p_raw"] for r in rows) > 0.10


# ------------------------------------------------------- the plant algebra ----
def test_narrow_arc_is_mean_preserving():
    th = np.linspace(0, 2 * np.pi, 100000, endpoint=False)
    for duty in (0.02, 0.10, 0.25):
        m = C.narrow_arc_intensity(th, 0.8, duty=duty)
        assert m.mean() == pytest.approx(1.0, abs=2e-3), duty


def test_fundamental_coefficient_matches_a_measured_fft():
    th = np.linspace(0, 2 * np.pi, 200000, endpoint=False)
    A, duty = 0.7, 0.10
    m = C.narrow_arc_intensity(th, A, duty=duty) - 1.0
    a1 = 2.0 * np.abs(np.fft.rfft(m)[1]) / m.size
    assert a1 == pytest.approx(C.fundamental_coefficient(A, duty), rel=2e-3)


def test_k_fold_arcs_kill_harmonics_1_and_2_exactly():
    """The algebra the decisive demonstration rests on, measured on the FFT."""
    th = np.linspace(0, 2 * np.pi, 300000, endpoint=False)
    m = C.k_fold_arc_intensity(th, 1.5, n_arcs=3, duty_each=0.05) - 1.0
    sp = 2.0 * np.abs(np.fft.rfft(m)) / m.size
    assert sp[1] < 1e-6, sp[1]          # fundamental: the Rayleigh-form is blind
    assert sp[2] < 1e-6, sp[2]          # second harmonic: F9-01 is blind
    assert sp[3] > 0.05, sp[3]          # third: Kuiper/Watson see the CDF excursion


def test_single_arc_kuiper_edge_is_bounded_by_pi_over_2():
    for d in (0.001, 0.01, 0.05, 0.10, 0.25):
        assert C.single_arc_kuiper_edge(d) <= C.SINGLE_ARC_KUIPER_EDGE_CEILING + 1e-9
    assert C.single_arc_kuiper_edge(1e-6) == pytest.approx(np.pi / 2, rel=1e-4)


def test_kuiper_equivalent_amplitude_round_trips():
    a = C.kuiper_equivalent_amplitude(2.7, 0.05)
    assert C.arc_amplitude_for_kuiper_equivalent(a, 0.05) == pytest.approx(2.7)


# ============================== THE S-17 RECOVERY DEMANDS, RUN SMALL ==========
BAND = 14.765
N_DAYS = 3000
RATE = 8.0
N_SURR = 300
ALPHA = 0.01


def _floor(n_days=N_DAYS, rate=RATE):
    return floors.a_min(floors.MEASURED_VIF_DF2_PHASE, floors.ALPHA_TRANCHE_B,
                        rate * n_days)


def test_f9_01_positive_control_two_lobed_plant_recovered():
    """§P7-3(1) positive control: a two-lobed plant with first moment ~ 0."""
    fl = _floor()
    amp = floors.PLANT_FACTOR * fl
    floors.assert_plant_above_floor(amp, RATE * N_DAYS, floors.ALPHA_TRANCHE_B,
                                    what="F9-01 positive control")
    th = np.mod(2 * np.pi * np.arange(N_DAYS) / BAND, 2 * np.pi)
    hits = 0
    for rep in range(3):
        _t, c, o, _m = R.simulate_phase_catalog(
            N_DAYS, BAND, R.antipodal_modulation(th, amp), seed=1200 + rep,
            rate=RATE)
        m2, _r1 = R._second_moment_and_rayleigh(th, c, o, BAND, N_SURR, 1200 + rep)
        hits += int(m2["p_raw"] <= ALPHA)
    assert hits == 3, hits


def test_f9_01_negative_control_pure_sinusoid_must_not_fire():
    """§P7-3(1)'s NEGATIVE control -- the one the whole statistic's meaning rests on.

    A planted pure sinusoid at the SAME amplitude as the positive control. The
    first-moment statistic must see it (otherwise the plant is too small and the
    control proves nothing); the SECOND-moment statistic must not.
    """
    fl = _floor()
    amp = floors.PLANT_FACTOR * fl
    th = np.mod(2 * np.pi * np.arange(N_DAYS) / BAND, 2 * np.pi)
    fired_second, fired_first = 0, 0
    for rep in range(4):
        s = 7700 + rep
        _t, c, o, _m = R.simulate_phase_catalog(
            N_DAYS, BAND, R.sinusoid_modulation(th, amp), seed=s, rate=RATE)
        m2, r1 = R._second_moment_and_rayleigh(th, c, o, BAND, N_SURR, s)
        fired_second += int(m2["p_raw"] <= ALPHA)
        fired_first += int(r1["p_raw"] <= ALPHA)
    assert fired_first >= 3, ("the plant is too small to make the negative control "
                             "meaningful: the FIRST moment did not see it either "
                             "(%d/4)" % fired_first)
    assert fired_second == 0, ("§P7-3(1) VIOLATED: the second-moment statistic fired "
                               "on a pure sinusoid %d/4 times" % fired_second)


def test_kuiper_sees_what_rayleigh_misses():
    """THE MONEY DEMONSTRATION (§P7-3(2), MINING_CATALOG F9-04).

    A concentrated-phase plant -- three sharp arcs of 5% duty each -- sized so its
    KUIPER-equivalent amplitude is 2x the operative floor. Three equally spaced arcs
    put exactly zero power in harmonics 1 and 2, so the Rayleigh-form 2-df statistic
    and the F9-01 second moment are blind BY ALGEBRA STATED IN ADVANCE (checked
    independently in `test_k_fold_arcs_kill_harmonics_1_and_2_exactly`), while Kuiper
    and Watson read the whole CDF.

    This is what F9-04 was admitted for, and the F9-01 column in the same assertion
    is what shows the two new statistics are COMPLEMENTARY rather than two names for
    one repair.
    """
    fl = _floor()
    a_eq = floors.PLANT_FACTOR * fl
    floors.assert_plant_above_floor(a_eq, RATE * N_DAYS, floors.ALPHA_TRANCHE_B,
                                    what="F9-04 decisive plant")
    amp = C.arc_amplitude_for_kuiper_equivalent(a_eq, 0.05)
    th = np.mod(2 * np.pi * np.arange(N_DAYS) / BAND, 2 * np.pi)
    k_hits = r_hits = m2_hits = 0
    reps = 4
    for rep in range(reps):
        s = 5500 + rep
        _t, c, o, _m = R.simulate_phase_catalog(
            N_DAYS, BAND, C.k_fold_arc_intensity(th, amp, 3, 0.05), seed=s,
            rate=RATE)
        omni, p_ray = R._omnibus_and_rayleigh(th, c, o, BAND, N_SURR, s)
        m2, _ = R._second_moment_and_rayleigh(th, c, o, BAND, N_SURR, s)
        k_hits += int(min(omni[0]["p_raw"], omni[1]["p_raw"]) <= ALPHA)
        r_hits += int(p_ray <= ALPHA)
        m2_hits += int(m2["p_raw"] <= ALPHA)
    assert k_hits >= reps - 1, ("Kuiper/Watson failed to recover a concentrated "
                                "plant at 2x its own floor: %d/%d" % (k_hits, reps))
    assert r_hits == 0, ("the Rayleigh-form statistic fired %d/%d times on a plant "
                         "with zero fundamental -- the demonstration is not what it "
                         "claims to be" % (r_hits, reps))
    assert m2_hits == 0, ("the F9-01 second moment fired %d/%d times on a plant with "
                          "zero second harmonic" % (m2_hits, reps))


def test_omnibus_p_raw_follows_the_engines_periodic_null_rule():
    """A periodic feature's p_raw is the bootstrap alone; otherwise the max."""
    rng = np.random.default_rng(11)
    n = 800
    th = np.mod(2 * np.pi * np.arange(n) / 27.3, 2 * np.pi)
    c = rng.poisson(6.0, n).astype(float)
    o = np.full(n, 6.0)
    per = C.omnibus_test(th, c, o, 150, np.random.default_rng(2), periodic=True)
    non = C.omnibus_test(th, c, o, 150, np.random.default_rng(2), periodic=False)
    for a, b in zip(per, non):
        assert a["p_raw"] == a["p_block_bootstrap"]
        assert b["p_raw"] == max(b["p_circular_shift"], b["p_block_bootstrap"])
