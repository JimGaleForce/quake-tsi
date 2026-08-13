"""D-4: the EVENT-PATH Kuiper/Watson/R1/2nd-moment port and its ETAS event-time null.

The demands §P7-22 attaches to this build, one test each:

  * a planted arc concentration at **2x the S-15 floor** is RECOVERED (§P7-8(d): a
    plant below the floor fails for power reasons while reading as an instrument
    failure);
  * the uniform-phase null is FLAT -- p-values uniform, gated by the one-sided KS in
    the anti-conservative direction only (§P6-9(a));
  * the §P7-21(c) fix is INHERITED: no circular-shift null for a periodic phase, and
    the p is mid-rank tie-tolerant;
  * the null is ETAS EVENT TIMES through the identical phase map, and a phase
    permutation cannot be passed (§P7-22(a)'s common-mode requirement).
"""

from __future__ import annotations

import datetime as _dt

import numpy as np
import pytest

from engine import circstat as C, circstat_event as CE, floors, sitetide as ST

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

T0 = _dt.datetime(2010, 1, 1)
N_DAYS = 1200
SITE = dict(lat_deg=38.0, lon_deg=142.0, depth_km=25.0)


def _intensity(n_days=N_DAYS):
    d = np.arange(n_days, dtype=np.float64)
    lam = 1.0 + 0.2 * np.cos(2 * np.pi * d / 365.2425)
    return np.arange(n_days + 1, dtype=np.float64), lam


def _phase_fn():
    """The D-0 waveform phase, ONE precomputed maxima array shared by all callers."""
    tmax = ST.tidal_maxima(T0, -1.0, N_DAYS + 1.0, **SITE)
    return lambda t: ST.phase_from_maxima(tmax, t)


# ------------------------------------------------------ the two-sample kernel --
def test_two_sample_kuiper_is_zero_for_identical_samples():
    th = np.linspace(0, 2 * np.pi, 500, endpoint=False)
    kw = CE.event_kuiper_watson(th, th)
    assert kw["V"] == pytest.approx(0.0, abs=1e-12)
    assert kw["U2"] == pytest.approx(0.0, abs=1e-12)


def test_two_sample_kuiper_is_rotation_invariant():
    """§P7-22(a) turns on this exact property: a rigid rotation cannot move V."""
    rng = np.random.default_rng(3)
    a = np.mod(rng.normal(1.0, 0.6, 4000), 2 * np.pi)
    b = rng.random(40000) * 2 * np.pi
    v0 = CE.event_kuiper_watson(a, b)["V"]
    for shift in (0.3, 1.7, 4.4):
        v = CE.event_kuiper_watson(a + shift, b + shift)["V"]
        assert v == pytest.approx(v0, rel=1e-9)


# ------------------------------------------------------------ planted arm ----
def test_planted_arc_at_twice_the_floor_is_recovered():
    """§P7-8(d): plant at 2x the operative floor, then require recovery.

    The floor is the S-15 formula at this N with VIF = 1 (declared for the harness;
    D-2 is what measures the operative event-path VIF). The sinusoidal-equivalent
    amplitude is converted to an ARC height with
    `circstat.arc_amplitude_for_kuiper_equivalent`, which is the units step §P7-8(d)
    exists to catch -- sizing by the arc height directly would compare a CDF excursion
    against a Fourier coefficient.
    """
    rng = np.random.default_rng(11)
    edges, lam = _intensity()
    fn = _phase_fn()
    n_events = 3000
    duty = 0.10

    a_eq = floors.min_plant_amplitude(n_events, alpha=floors.ALPHA_TRANCHE_B,
                                      vif=1.0, factor=floors.PLANT_FACTOR)
    rep = floors.assert_plant_above_floor(a_eq, n_events,
                                          alpha=floors.ALPHA_TRANCHE_B, vif=1.0,
                                          what="D-4 event-path arc plant")
    assert rep["compliant"]
    arc = C.arc_amplitude_for_kuiper_equivalent(a_eq, duty)

    cand = CE.resample_event_times(edges, lam, 60000, rng)
    obs_t = CE.thin_by_phase_intensity(
        cand, fn(cand), lambda th: C.narrow_arc_intensity(th, arc, duty=duty),
        n_events, rng)

    nulls = [CE.resample_event_times(edges, lam, n_events, rng) for _ in range(199)]
    ref = CE.resample_event_times(edges, lam, 30000, rng)
    rows = CE.event_omnibus(fn(obs_t), nulls, fn, reference_times=ref)

    by = {r["test"]: r for r in rows}
    assert by["kuiper_V"]["primary"] is True
    assert by["kuiper_V"]["p_raw"] <= 0.01, by["kuiper_V"]["p_raw"]
    assert by["watson_U2"]["p_raw"] <= 0.01


def test_uniform_phase_null_is_flat_ks_one_sided():
    """§P6-9(a): p-values uniform, gated on D+ only (D- is permitted conservatism)."""
    rng = np.random.default_rng(5)
    edges, lam = _intensity()
    fn = _phase_fn()
    n_events = 800
    ref = CE.resample_event_times(edges, lam, 20000, rng)
    ref_ph = fn(ref)

    ps = []
    for _ in range(60):
        obs = CE.resample_event_times(edges, lam, n_events, rng)
        nulls = [CE.resample_event_times(edges, lam, n_events, rng)
                 for _ in range(79)]
        rows = CE.event_omnibus(fn(obs), nulls, fn, reference_phases=ref_ph)
        ps.append({r["test"]: r["p_raw"] for r in rows}["kuiper_V"])

    ks = CE.ks_uniform_one_sided(ps, a=0.01)
    assert ks["pass"], ks
    assert 0.3 < ks["mean_p"] < 0.7, ks["mean_p"]


# ------------------------------------------------- the inherited §P7-21 fix ---
def test_no_circular_shift_null_for_a_periodic_phase():
    """§P7-21(c), inherited: the shift null is ABSENT with its reason on the row."""
    rng = np.random.default_rng(7)
    edges, lam = _intensity()
    fn = _phase_fn()
    obs = CE.resample_event_times(edges, lam, 800, rng)
    nulls = [CE.resample_event_times(edges, lam, 800, rng) for _ in range(49)]
    ref = CE.resample_event_times(edges, lam, 20000, rng)
    for r in CE.event_omnibus(fn(obs), nulls, fn, reference_times=ref):
        assert r["p_circular_shift"] is None
        assert r["p_circular_shift_note"] is C.PERIODIC_SHIFT_OMITTED
        assert r["p_raw"] == r["p_etas_event_null"]


def test_p_is_mid_rank_tie_tolerant_and_continuous():
    """The §P7-21 repair itself: jitter in the observed statistic moves p by O(1/B)."""
    draws = np.array([1.0] * 50 + [2.0] * 50)
    p_lo = C.tie_tolerant_p(draws, 1.0 - 1e-12, draws.size)
    p_hi = C.tie_tolerant_p(draws, 1.0 + 1e-12, draws.size)
    assert abs(p_lo - p_hi) < 1e-9
    # ties count as half: 50 above, 50 tied -> (1 + 50 + 25) / 101
    assert C.tie_tolerant_p(draws, 1.0, draws.size) == pytest.approx(76.0 / 101.0)


# ------------------------------------------------------ the common-mode null --
def test_the_null_takes_event_times_and_a_phase_permutation_cannot_be_passed():
    """§P7-22(a): the design is saved by refusing a phase permutation. By signature.

    `event_omnibus` calls `phase_fn` on every null element itself, so what a caller
    supplies is TIMES. Handing it permuted phases and an identity map is possible only
    by lying about the map, and the test records that the ARGUMENT is times -- the
    deterministic warp of the phase computation then cancels common-mode.
    """
    import inspect
    sig = inspect.signature(CE.event_omnibus)
    assert "null_times" in sig.parameters
    assert "phase_fn" in sig.parameters
    assert "null_phases" not in sig.parameters
    assert "S-8" in CE.COMMON_MODE_NOTE or "identical" in CE.COMMON_MODE_NOTE


def test_a_deterministic_phase_warp_cancels_between_signal_and_null():
    """The cancellation §P7-22(a) rests on, demonstrated rather than asserted.

    Apply a strong deterministic warp to the phase map. Under a true null the p-value
    distribution must be unmoved, because the warp is common-mode to observed, null
    and reference.
    """
    rng = np.random.default_rng(13)
    edges, lam = _intensity()
    base = _phase_fn()

    def warped(t):
        th = base(t)
        return np.mod(th + 0.9 * np.sin(2.0 * th) + 0.4 * np.cos(3.0 * th), 2 * np.pi)

    ref = CE.resample_event_times(edges, lam, 20000, rng)
    ps = []
    for _ in range(30):
        obs = CE.resample_event_times(edges, lam, 800, rng)
        nulls = [CE.resample_event_times(edges, lam, 800, rng) for _ in range(79)]
        rows = CE.event_omnibus(warped(obs), nulls, warped, reference_times=ref)
        ps.append({r["test"]: r["p_raw"] for r in rows}["kuiper_V"])
    ks = CE.ks_uniform_one_sided(ps, a=0.01)
    assert ks["pass"], ks


# --------------------------------------------------------------- the record ---
def test_rows_carry_the_scope_clauses_they_are_required_to_carry():
    rng = np.random.default_rng(2)
    edges, lam = _intensity()
    fn = _phase_fn()
    obs = CE.resample_event_times(edges, lam, 500, rng)
    nulls = [CE.resample_event_times(edges, lam, 500, rng) for _ in range(29)]
    ref = CE.resample_event_times(edges, lam, 10000, rng)
    rows = CE.event_omnibus(fn(obs), nulls, fn, reference_times=ref)
    assert {r["test"] for r in rows} == {"kuiper_V", "watson_U2", "R1_contrast",
                                         "second_moment_contrast"}
    for r in rows:
        assert "may not be quoted" in r["fail_is_not_a_bound"]
        assert "psi" in r["psi_unreportable"]
        assert r["path"] == "event"
    assert sum(1 for r in rows if r["primary"]) == 1
