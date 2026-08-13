"""D-3: `engine/sitetide.py` -- the degree-2 Love-number body tide at a site.

The acceptance tests §K92-1 D-3 declares, plus the three the build brief adds:

  * the M2 / S2 / K1 / O1 constituent periods are PRESENT in a spectral check of a
    synthetic year at a mid-latitude site;
  * the amplitude is order 1e-8 areal strain (and ~1 kPa in the F-016 stress
    convention, which is the quantity this module states it emits);
  * phase is CONTINUOUS across day boundaries -- the property that separates a genuine
    sub-daily feature from a day-binned one evaluated at an event (§P7-3(3));
  * the Doodson constant comes out at its published 2.628 m^2/s^2, which checks the
    whole potential path in one number.
"""

from __future__ import annotations

import datetime as _dt

import numpy as np
import pytest

from engine import sitetide as ST

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

T0 = _dt.datetime(2020, 1, 1)
MID_LAT_SITE = dict(lat_deg=38.0, lon_deg=142.0, depth_km=20.0)   # mid-latitude


# ------------------------------------------------------------- the constants --
def test_doodson_constant_matches_published_value():
    """(3/4) GM_moon a^2 / d^3 = 2.628 m^2/s^2 -- the whole potential path in one number."""
    assert ST.doodson_constant() == pytest.approx(2.628, abs=0.002)


def test_love_number_derived_constants():
    """The two combinations an implementation is most likely to get wrong."""
    # areal-strain factor 2 h2 - n(n+1) l2 at n = 2
    assert ST.AREAL_FACTOR == pytest.approx(2 * 0.6078 - 6 * 0.0847, abs=1e-12)
    assert ST.AREAL_FACTOR == pytest.approx(0.7074, abs=1e-4)
    # gravimetric factor delta = 1 + h2 - 1.5 k2, nominal ~1.155
    assert ST.GRAVIMETRIC_DELTA == pytest.approx(1.155, abs=0.002)


# ------------------------------------------------------------- the spectrum --
def test_m2_s2_k1_o1_present_in_a_synthetic_year():
    """§K92-1 D-3's acceptance test: recover the known constituent periods."""
    f, a = ST.constituent_spectrum(T0, n_days=365.0, sample_minutes=30.0,
                                   **MID_LAT_SITE)
    noise = float(np.median(a))
    for name in ("M2", "S2", "K1", "O1"):
        period = ST.CONSTITUENT_PERIODS_DAYS[name]
        f_peak, a_peak = ST.peak_near(f, a, period, half_width_cpd=0.02)
        # the peak sits at the constituent frequency to within one FFT bin
        assert abs(f_peak - 1.0 / period) < 2.0 / 365.0, name
        # and it is a peak, not a ripple: >= 100x the median spectral level
        assert a_peak > 100.0 * noise, (name, a_peak, noise)


def test_semidiurnal_and_diurnal_split_between_moon_and_sun():
    """S2 is SOLAR: dropping the Sun must kill it and leave M2 essentially intact."""
    f, a_both = ST.constituent_spectrum(T0, n_days=365.0, sample_minutes=30.0,
                                        **MID_LAT_SITE)
    _f, a_moon = ST.constituent_spectrum(T0, n_days=365.0, sample_minutes=30.0,
                                         bodies=("moon",), **MID_LAT_SITE)
    _fm, m_both = ST.peak_near(f, a_both, ST.CONSTITUENT_PERIODS_DAYS["M2"], 0.005)
    _fm2, m_moon = ST.peak_near(f, a_moon, ST.CONSTITUENT_PERIODS_DAYS["M2"], 0.005)
    _fs, s_both = ST.peak_near(f, a_both, ST.CONSTITUENT_PERIODS_DAYS["S2"], 0.005)
    _fs2, s_moon = ST.peak_near(f, a_moon, ST.CONSTITUENT_PERIODS_DAYS["S2"], 0.005)
    assert m_moon == pytest.approx(m_both, rel=0.05)
    assert s_moon < 0.05 * s_both


# ------------------------------------------------------------ the amplitudes --
def test_areal_strain_is_order_1e_minus_8_and_the_stress_equivalent_is_kPa():
    """The declared emitted quantity, at its declared order of magnitude."""
    jd = ST.E.julian_day_at(T0, np.arange(0.0, 60.0, 1.0 / 96.0))
    rec = ST.site_tide(jd, **MID_LAT_SITE)
    peak = float(np.max(np.abs(rec["areal_strain"])))
    assert 1e-9 < peak < 1e-7, peak                       # order 1e-8
    # F-016 convention: sigma_m = E * e_areal / (3 (1 - nu)), E = 75 GPa, nu = 0.25
    s = float(np.max(np.abs(rec["mean_stress_pa"])))
    assert s == pytest.approx(peak * 75.0e9 / (3.0 * 0.75), rel=1e-12)
    assert 200.0 < s < 5000.0, s                          # ~1 kPa body-tide scale
    # and the equilibrium tide is the familiar few tens of cm
    assert 0.1 < float(np.max(np.abs(rec["equilibrium_m"]))) < 0.8


def test_radial_displacement_scales_with_h2():
    jd = ST.E.julian_day_at(T0, np.arange(0.0, 30.0, 0.01))
    rec = ST.site_tide(jd, **MID_LAT_SITE)
    np.testing.assert_allclose(rec["radial_disp_m"],
                               ST.H2 * rec["equilibrium_m"], rtol=1e-12)


def test_depth_reduces_the_potential_by_the_r_squared_scaling():
    """Depth enters as (r/a)^2 and nothing else -- stated in the docstring, checked here."""
    jd = ST.E.julian_day_at(T0, np.arange(0.0, 5.0, 0.01))
    shallow = ST.site_scalar(jd, 38.0, 142.0, 0.0)
    deep = ST.site_scalar(jd, 38.0, 142.0, 50.0)
    phi, _lon, r0 = ST.geocentric_site(38.0, 142.0, 0.0)
    _phi, _l, r1 = ST.geocentric_site(38.0, 142.0, 50.0)
    np.testing.assert_allclose(deep, shallow * (r1 / r0) ** 2, rtol=1e-10)
    assert 0.98 < float(np.mean(deep / shallow)) < 1.0


def test_latitude_dependence_semidiurnal_vanishes_at_the_pole():
    """cos^2(lat) semidiurnal geometry: the M2 line must die at high latitude."""
    f, a_mid = ST.constituent_spectrum(T0, 38.0, 142.0, 10.0, n_days=200.0,
                                       sample_minutes=30.0)
    _f, a_pol = ST.constituent_spectrum(T0, 88.0, 142.0, 10.0, n_days=200.0,
                                        sample_minutes=30.0)
    _fm, m_mid = ST.peak_near(f, a_mid, ST.CONSTITUENT_PERIODS_DAYS["M2"], 0.01)
    _fp, m_pol = ST.peak_near(f, a_pol, ST.CONSTITUENT_PERIODS_DAYS["M2"], 0.01)
    assert m_pol < 0.02 * m_mid


# ------------------------------------------------------------ the continuity --
def test_scalar_is_continuous_across_day_boundaries():
    """No day-binning seam: the scalar is a smooth function of continuous time.

    The failure this catches is a feature that is really day-binned wearing an
    event-time signature -- §P7-3(3)/F9-10's whole distinction. A day-binned scalar
    would step at every integer `day_float`; a genuine one does not.
    """
    eps = 1e-6
    at_boundary = np.arange(1.0, 40.0)                 # integer day_float
    off_boundary = at_boundary + 0.37                  # nothing special about these
    def step(x):
        return np.abs(ST.site_scalar_at(T0, x + eps, **MID_LAT_SITE)
                      - ST.site_scalar_at(T0, x - eps, **MID_LAT_SITE))
    # The increment across an integer day_float must be no larger than the increment
    # across an ordinary interior time. A day-binned scalar would step by a full
    # cycle's worth here and by exactly zero there.
    assert float(step(at_boundary).max()) < 3.0 * float(step(off_boundary).max())
    assert float(step(off_boundary).max()) > 0.0


def test_tanaka_phase_is_continuous_and_uniform_in_time():
    """D-0's phase: 0 at a maximum, wraps once per cycle, no day-boundary seam."""
    t = np.arange(0.0, 60.0, 1.0 / 288.0)          # 5-minute event times, 60 days
    ph, info = ST.tanaka_phase(T0, t, **MID_LAT_SITE)
    assert info["n_nan"] == 0
    assert np.all((ph >= 0.0) & (ph < 2 * np.pi))
    # the mean cycle is semidiurnal-ish at a mid-latitude site
    assert 10.0 < info["mean_cycle_hours"] < 14.0
    # phase advances monotonically within a cycle: exactly one wrap per maximum
    # (the maxima array spans one pad day beyond the events at each end, so a few
    # cycles carry no events and produce no wrap)
    n_wraps = int(np.sum(np.diff(ph) < -np.pi))
    assert n_wraps <= info["n_maxima"] - 1
    assert n_wraps >= info["n_maxima"] - 6
    # uniform-in-time sampling gives a near-uniform phase histogram, by construction
    h, _e = np.histogram(ph, bins=12, range=(0.0, 2 * np.pi))
    assert h.min() > 0.5 * h.mean()
    # the grid error is stated and small
    assert info["phase_grid_error_rad"] < 1e-3


def test_phase_from_precomputed_maxima_matches_the_one_shot_call():
    """The shared-maxima fast path is BIT-comparable to the direct call.

    This is §P7-22(a)'s common-mode requirement as an invariant: the D-2/D-4 harnesses
    phase signal and null against ONE precomputed maxima array, and that must be the
    same map the one-shot call uses or the cancellation argument does not hold.
    """
    t = np.arange(0.5, 30.0, 0.017)
    ph_direct, _i = ST.tanaka_phase(T0, t, **MID_LAT_SITE)
    tmax = ST.tidal_maxima(T0, -1.0, 31.0, **MID_LAT_SITE)
    ph_shared = ST.phase_from_maxima(tmax, t)
    np.testing.assert_allclose(ph_direct, ph_shared, atol=1e-12)


def test_constituent_phase_is_a_pure_clock_and_not_the_d0_convention():
    """The two phase functions are different objects and must not be confused."""
    t = np.arange(0.0, 10.0, 0.01)
    p = ST.constituent_phase(T0, t, ST.CONSTITUENT_PERIODS_DAYS["M2"])
    assert np.all((p >= 0) & (p < 2 * np.pi))
    # exactly periodic in the constituent period
    np.testing.assert_allclose(
        p, ST.constituent_phase(T0, t + ST.CONSTITUENT_PERIODS_DAYS["M2"],
                                ST.CONSTITUENT_PERIODS_DAYS["M2"]), atol=1e-9)
    w, _i = ST.tanaka_phase(T0, t, **MID_LAT_SITE)
    assert float(np.mean(np.abs(np.mod(w - p + np.pi, 2 * np.pi) - np.pi))) > 0.2


# ------------------------------------------------------------------ scope ----
def test_scope_flags_are_attached_to_every_record():
    """A record that could be quoted must carry what it is not entitled to claim."""
    jd = ST.E.julian_day_at(T0, np.array([0.0, 1.0]))
    rec = ST.site_tide(jd, **MID_LAT_SITE)
    assert "ocean loading" in rec["scope_flags"].lower()
    assert "K-090(c)" in rec["scope_flags"]
    assert rec["scalar_for_phase"] == ST.SCALAR_FOR_PHASE
    assert rec["love_numbers"]["source"].startswith("IERS Conventions (2010)")
