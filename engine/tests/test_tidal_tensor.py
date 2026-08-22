"""Acceptance tests for engine/tidal_tensor.py.

The load-bearing one is `test_trace_equals_sitetide_areal_strain`: the tensor is built
from the Shida number and the analytic angular derivatives, by a completely different
route from `sitetide.areal_strain`, and its trace must reproduce that scalar to
floating-point. That is an independent re-derivation of AREAL_FACTOR = 2h2 - 6l2, so if
the tensor were wrong the identity would break. Everything else here is guarding the
pieces that identity cannot see: the shear component, the bearing algebra, the axial
doubling, and the Coulomb resolution.
"""

from __future__ import annotations

import datetime as _dt

import numpy as np
import pytest

from engine import ephemeris as E
from engine import sitetide as ST
from engine import tidal_tensor as TT

# A site with no special symmetry, and one that is not the Alaska seed site: these
# tests must not be a look at anything.
LAT, LON, DEPTH = 35.7, 139.7, 0.0
T0 = _dt.datetime(2020, 6, 1)
DAYS = np.linspace(0.0, 3.0, 433)
JD = E.julian_day_at(T0, DAYS)


def test_trace_equals_sitetide_areal_strain():
    """(l2/ga) Laplacian(W) + 2 h2 W/ga must equal (2h2 - 6l2) W/ga, exactly."""
    s = TT.strain_tensor(JD, LAT, LON, DEPTH)
    ref = ST.site_scalar(JD, LAT, LON, DEPTH, which="areal_strain")
    assert np.allclose(s["areal_strain"], ref, rtol=1e-11, atol=1e-18)
    # and the constant it re-derives is the one sitetide hard-codes
    assert abs((2.0 * TT.H2 - 6.0 * TT.L2) - ST.AREAL_FACTOR) < 1e-15


def test_radial_strain_matches_sitetide():
    s = TT.strain_tensor(JD, LAT, LON, DEPTH)
    ref = ST.site_tide(JD, LAT, LON, DEPTH)["radial_strain"]
    assert np.allclose(s["radial_strain"], ref, rtol=1e-11, atol=1e-18)


@pytest.mark.parametrize("lat,lon", [(35.7, 139.7), (-22.0, -70.0), (58.0, -152.0)])
def test_analytic_derivatives_against_finite_differences(lat, lon):
    """W_p, W_pp, W_l, W_ll, W_pl vs central differences of sitetide.tidal_potential.

    The finite differences are taken in GEOCENTRIC latitude, so the geodetic input is
    inverted first -- doing it in geodetic latitude would silently compare two
    different derivatives and pass to within a few percent, which is exactly the kind
    of near-miss the program's failure catalogue warns about.
    """
    jd = float(E.julian_day_at(T0, 1.234))
    phi_gc, lam, _ = ST.geocentric_site(lat, lon, 0.0)
    d = TT.potential_and_derivatives(jd, lat, lon, 0.0)

    def w_of(phi, lm):
        """W at a geocentric latitude/longitude, holding the radius fixed."""
        _, _, r_s = ST.geocentric_site(lat, lon, 0.0)
        tot = 0.0
        for name, gm, pos, dist in (
            ("moon", ST.GM_MOON, E.moon_position(jd),
             np.asarray(E.moon_position(jd)["dist_km"]) * 1000.0),
            ("sun", ST.GM_SUN, E.sun_position(jd),
             np.asarray(E.sun_position(jd)["dist_au"]) * 1.495978707e11),
        ):
            ha = ST.gmst_rad(jd) + lm - pos["ra_deg"] * TT.DEG
            dec = pos["dec_deg"] * TT.DEG
            c = (np.sin(phi) * np.sin(dec)
                 + np.cos(phi) * np.cos(dec) * np.cos(ha))
            tot = tot + (gm / dist) * (r_s / dist) ** 2 * 0.5 * (3 * c * c - 1)
        return float(tot)

    h = 1e-5
    fd_p = (w_of(phi_gc + h, lam) - w_of(phi_gc - h, lam)) / (2 * h)
    fd_pp = (w_of(phi_gc + h, lam) - 2 * w_of(phi_gc, lam)
             + w_of(phi_gc - h, lam)) / (h * h)
    fd_l = (w_of(phi_gc, lam + h) - w_of(phi_gc, lam - h)) / (2 * h)
    fd_ll = (w_of(phi_gc, lam + h) - 2 * w_of(phi_gc, lam)
             + w_of(phi_gc, lam - h)) / (h * h)
    fd_pl = ((w_of(phi_gc + h, lam + h) - w_of(phi_gc + h, lam - h)
              - w_of(phi_gc - h, lam + h) + w_of(phi_gc - h, lam - h))
             / (4 * h * h))

    scale = max(abs(float(d["W"])), 1e-12)
    assert abs(float(d["W_p"]) - fd_p) < 1e-5 * scale
    assert abs(float(d["W_pp"]) - fd_pp) < 1e-4 * scale
    assert abs(float(d["W_l"]) - fd_l) < 1e-5 * scale
    assert abs(float(d["W_ll"]) - fd_ll) < 1e-4 * scale
    assert abs(float(d["W_pl"]) - fd_pl) < 1e-4 * scale


def test_bearing_is_the_argmax_of_normal_strain():
    """principal_bearing must agree with a brute-force scan of normal_along_bearing."""
    s = TT.strain_tensor(JD[:37], LAT, LON, DEPTH)
    b = TT.principal_bearing(s["e_NN"], s["e_EE"], s["e_NE"])
    grid = np.linspace(0.0, 180.0, 3601)[:-1]
    for i in range(len(b)):
        vals = TT.normal_along_bearing(s["e_NN"][i], s["e_EE"][i], s["e_NE"][i], grid)
        brute = grid[int(np.argmax(vals))]
        d = abs(brute - b[i])
        assert min(d, 180.0 - d) < 0.2, (i, brute, b[i])


def test_bearing_value_equals_larger_eigenvalue():
    s = TT.strain_tensor(JD[:23], LAT, LON, DEPTH)
    b = TT.principal_bearing(s["e_NN"], s["e_EE"], s["e_NE"])
    got = TT.normal_along_bearing(s["e_NN"], s["e_EE"], s["e_NE"], b)
    for i in range(len(b)):
        m = np.array([[s["e_NN"][i], s["e_NE"][i]], [s["e_NE"][i], s["e_EE"][i]]])
        assert abs(got[i] - np.linalg.eigvalsh(m).max()) < 1e-18 + 1e-9 * abs(got[i])


def test_normal_along_bearing_hand_cases():
    """Pure North extension: bearing 0 sees it all, bearing 90 sees none."""
    assert TT.normal_along_bearing(1.0, 0.0, 0.0, 0.0) == pytest.approx(1.0)
    assert TT.normal_along_bearing(1.0, 0.0, 0.0, 90.0) == pytest.approx(0.0, abs=1e-15)
    assert TT.normal_along_bearing(1.0, 0.0, 0.0, 45.0) == pytest.approx(0.5)
    # pure shear: principal axes at +-45
    assert TT.principal_bearing(0.0, 0.0, 1.0) == pytest.approx(45.0)


def test_axial_doubling_identifies_antipodal_bearings():
    assert TT.axial_to_circular(10.0) == pytest.approx(TT.axial_to_circular(190.0))
    assert TT.axial_to_circular(0.0) == pytest.approx(0.0)


def test_bearing_rotation_rate_on_a_known_sweep():
    """A bearing sweeping at a constant 5 deg/h must read back 5 deg/h through 180."""
    t_days = np.linspace(0.0, 2.0, 2001)
    b = np.mod(5.0 * t_days * 24.0, 180.0)
    r = TT.bearing_rotation_rate(t_days, b)
    assert np.allclose(r[5:-5], 5.0, atol=1e-6)


def test_fault_vectors_are_orthonormal():
    for strike in (0.0, 37.0, 250.0, 359.0):
        for dip in (5.0, 20.0, 60.0, 89.0):
            for rake in (-120.0, 0.0, 90.0, 170.0):
                n, u = TT.fault_vectors(strike, dip, rake)
                assert abs(np.linalg.norm(n) - 1.0) < 1e-12
                assert abs(np.linalg.norm(u) - 1.0) < 1e-12
                assert abs(float(np.dot(n, u))) < 1e-12


def test_fault_vectors_known_geometry():
    """Strike 0, dip 90, rake 0: a vertical N-S left-lateral plane."""
    n, u = TT.fault_vectors(0.0, 90.0, 0.0)
    assert np.allclose(n, [0.0, 1.0, 0.0], atol=1e-12)     # normal points East
    assert np.allclose(u, [1.0, 0.0, 0.0], atol=1e-12)     # slip points North
    # horizontal thrust plane, dip 0: normal is straight up (Down component -1)
    n2, _ = TT.fault_vectors(0.0, 0.0, 90.0)
    assert np.allclose(n2, [0.0, 0.0, -1.0], atol=1e-12)


def test_coulomb_hand_computed_vertical_strike_slip():
    """Uniaxial N-S tension on a vertical N-S plane: no shear, no normal traction."""
    stress = {"s_NN": np.array([1.0e6]), "s_EE": np.array([0.0]),
              "s_NE": np.array([0.0])}
    c = TT.coulomb(stress, 0.0, 90.0, 0.0, friction=0.4)
    assert float(c["normal_pa"][0]) == pytest.approx(0.0, abs=1e-6)
    assert float(c["shear_pa"][0]) == pytest.approx(0.0, abs=1e-6)
    # rotate the plane to E-W: now it feels the full normal traction and no shear
    c2 = TT.coulomb(stress, 90.0, 90.0, 0.0, friction=0.4)
    assert float(c2["normal_pa"][0]) == pytest.approx(1.0e6, rel=1e-12)
    assert float(c2["shear_pa"][0]) == pytest.approx(0.0, abs=1e-6)
    # CFS = tau + mu * sigma_n = 0 + 0.4 * 1e6. Not sigma_n + mu * sigma_n: the
    # normal traction enters ONLY through the friction term.
    assert float(c2["coulomb_pa"][0]) == pytest.approx(0.4 * 1.0e6, rel=1e-12)


def test_coulomb_pure_shear_on_45_degree_plane():
    """Pure NE shear: a vertical plane striking 45 feels pure normal traction."""
    stress = {"s_NN": np.array([0.0]), "s_EE": np.array([0.0]),
              "s_NE": np.array([1.0e6])}
    c = TT.coulomb(stress, 45.0, 90.0, 0.0, friction=0.0)
    assert abs(float(c["normal_pa"][0])) == pytest.approx(1.0e6, rel=1e-9)
    assert float(c["shear_pa"][0]) == pytest.approx(0.0, abs=1.0)


def test_stress_tensor_is_plane_stress():
    s = TT.stress_tensor(JD[:11], LAT, LON, DEPTH)
    assert np.all(s["s_DD"] == 0.0)
    # plane-stress inversion must return the strains it came from
    ey, nu = TT.YOUNGS_MODULUS_PA, TT.POISSON_RATIO
    e_nn = (s["s_NN"] - nu * s["s_EE"]) / ey
    assert np.allclose(e_nn, s["e_NN"], rtol=1e-10)


def test_shear_along_bearing_vanishes_on_principal_axes():
    s = TT.strain_tensor(JD[:29], LAT, LON, DEPTH)
    b = TT.principal_bearing(s["e_NN"], s["e_EE"], s["e_NE"])
    sh = TT.shear_along_bearing(s["e_NN"], s["e_EE"], s["e_NE"], b)
    assert np.max(np.abs(sh)) < 1e-9 * np.max(np.abs(s["areal_strain"]))


def test_poles_raise_rather_than_return_nonsense():
    with pytest.raises(ValueError):
        TT.strain_tensor(JD[:3], 90.0, 0.0, 0.0)


def test_scope_flags_travel():
    s = TT.strain_tensor(JD[:3], LAT, LON, DEPTH)
    assert "OCEAN" in s["scope_flags"].upper()
    assert "PLANE STRESS" in s["scope_flags"].upper()
    c = TT.coulomb(TT.stress_tensor(JD[:3], LAT, LON, DEPTH), 250.0, 20.0, 90.0)
    assert "scope_flags" in c and "DECLARED ASSUMPTION" in c["scope_flags"].upper()


# ------------------------------------------- degree 3 (§P7-25 disposition 22) --
def test_degree3_is_off_by_default_and_changes_nothing():
    """The repair must be strictly additive: default-off reproduces the old numbers."""
    a = ST.site_tide(JD[:64], LAT, LON, DEPTH)
    b = ST.site_tide(JD[:64], LAT, LON, DEPTH, degree3=False)
    assert np.array_equal(a["areal_strain"], b["areal_strain"])
    assert np.array_equal(a["radial_disp_m"], b["radial_disp_m"])
    assert a["degree3"] is False and a["potential_moon_degree3"] is None


def test_degree3_amplitude_scaling():
    """W3/W2 = (r/d) * P3(cos z)/P2(cos z).

    The COEFFICIENT ratio is r/d ~ 6371/384400 ~ 1/60; the realised amplitude ratio
    also carries |P3|/|P2|, which averages well below one, so the measured figure is
    ~1/127 and not ~1/60. This test's first version asserted the naive 1/60 band and
    failed at 1/127 -- the bound was wrong, not the code, and it is recorded here so
    nobody re-tightens it. The coefficient scaling itself is asserted exactly below.
    """
    rec = ST.site_tide(JD, LAT, LON, DEPTH, degree3=True)
    w2 = np.abs(np.asarray(rec["potential_moon"], dtype=float))
    w3 = np.abs(np.asarray(rec["potential_moon_degree3"], dtype=float))
    hi = w2 > np.percentile(w2, 50)
    ratio = np.mean(w3[hi]) / np.mean(w2[hi])
    assert 1.0 / 300.0 < ratio < 1.0 / 25.0, ratio


def test_degree3_coefficient_scaling_is_exactly_r_over_d():
    """W3/W2 divided by P3/P2 must equal r_site/d_moon, to floating point."""
    jd = float(E.julian_day_at(T0, 0.77))
    rec = ST.site_tide(jd, LAT, LON, DEPTH, degree3=True)
    phi, lam, r_s = ST.geocentric_site(LAT, LON, DEPTH)
    m = E.moon_position(jd)
    d = float(np.asarray(m["dist_km"])) * 1000.0
    cz = ST._cos_zenith(phi, lam, jd, m["ra_deg"], m["dec_deg"])
    got = (float(rec["potential_moon_degree3"]) / float(rec["potential_moon"]))
    want = (r_s / d) * float(ST._p3(cz)) / float(ST._p2(cz))
    assert abs(got - want) < 1e-12 * abs(want)


def test_degree3_uses_its_own_love_numbers():
    """2 h3 - 12 l3, not the degree-2 factor. Mixing them is the likely error here."""
    assert abs(ST.AREAL_FACTOR_3 - (2.0 * ST.H3 - 12.0 * ST.L3)) < 1e-15
    assert abs(ST.AREAL_FACTOR_3 - 0.4016) < 1e-4
    assert ST.AREAL_FACTOR_3 != ST.AREAL_FACTOR


def test_degree3_perturbs_but_does_not_dominate():
    """Enabling it must move the scalar by a few percent, not reorder it."""
    off = ST.site_scalar(JD, LAT, LON, DEPTH)
    on = ST.site_scalar(JD, LAT, LON, DEPTH, degree3=True)
    rel = np.max(np.abs(on - off)) / np.max(np.abs(off))
    assert 1e-4 < rel < 0.10, rel
    assert np.corrcoef(off, on)[0, 1] > 0.999


# --------------- the free-surface degeneracy (Kepler K-111 / K-108, verified) --
def test_free_surface_coulomb_is_a_scalar_multiple_of_the_fault_normal_bearing():
    """AT THE FREE SURFACE, FRICTION AND DIP ARE ALGEBRAICALLY INERT. Assert it.

    With sigma_DD = sigma_rN = sigma_rE = 0 the fault normal's only surviving part is
    its horizontal component sin(dip) * m, where m is the unit vector at azimuth
    strike + 90. Hence for rake 90:

        sigma_n = sin^2(d) * N(s+90)
        tau     = -sin(d) cos(d) * N(s+90)
        CFS     = sin(d) (mu sin(d) - cos(d)) * N(s+90)

    Both tractions are the SAME time function scaled by a geometry constant, so
    gridding friction or dip at depth 0 produces cells that are identical up to a
    positive scale and cannot be distinct tests.

    This was found AFTER exp_k092_d1c_directional.py had already gridded 5 strikes x
    3 dips x 3 frictions and read the result as 45 cells. It is 5. The assertion
    exists so that a future result which DOES vary with mu is immediately recognised
    as evidence that the vertical stress components have become non-zero (depth, or
    an ocean load) rather than as a finding about friction.
    """
    jd = E.julian_day_at(T0, np.linspace(0.0, 20.0, 14401))
    st = TT.stress_tensor(jd, 55.35, -160.5, 0.0)
    for strike in (230.0, 250.0, 270.0):
        n_perp = TT.normal_along_bearing(st["s_NN"], st["s_EE"], st["s_NE"],
                                         strike + 90.0)
        ref = None
        for dip in (10.0, 20.0, 30.0):
            for mu in (0.2, 0.4, 0.6):
                c = TT.coulomb(st, strike, dip, 90.0, mu)["coulomb_pa"]
                d = np.radians(dip)
                want = np.sin(d) * (mu * np.sin(d) - np.cos(d)) * n_perp
                assert np.allclose(c, want, rtol=1e-12, atol=1e-9 * np.abs(c).max())
                # every cell at this strike is the same series up to positive scale
                if ref is None:
                    ref = c
                else:
                    assert abs(abs(np.corrcoef(ref, c)[0, 1]) - 1.0) < 1e-9


def test_stressing_rate_is_orthogonal_to_the_level_by_identity():
    """d/dt of a stationary field is uncorrelated with the field at lag zero.

    This is why the rate axis is the exact orthogonal complement of the level axis
    that D-1 bounded, and it is a theorem rather than a measurement.
    """
    t = np.linspace(0.0, 30.0, 43201)
    jd = E.julian_day_at(T0, t)
    s = ST.site_scalar(jd, 55.35, -160.5, 0.0)
    ds = np.gradient(s, t * 24.0)
    assert abs(float(np.corrcoef(s, ds)[0, 1])) < 1e-3


# ------------------------------- body azimuth / elevation (Jim's original list) --
def test_sin_elevation_is_exactly_the_potentials_own_cos_zenith():
    """The identity that keeps these angles consistent with the tide beside them.

    If these were computed from a second, better ephemeris they would not correspond
    to the potential they are supposed to explain. Exact equality, not approximate.
    """
    jd = E.julian_day_at(T0, np.arange(0.0, 3.0, 1.0 / 1440.0))
    phi, lon, _ = ST.geocentric_site(55.35, -160.5, 0.0)
    for body, posfn in (("moon", E.moon_position), ("sun", E.sun_position)):
        d = TT.body_direction(jd, 55.35, -160.5, 0.0, body)
        pos = posfn(jd)
        cz = ST._cos_zenith(phi, lon, jd, pos["ra_deg"], pos["dec_deg"])
        assert np.max(np.abs(d["sin_elevation"] - cz)) == 0.0


def test_the_sun_transits_south_in_the_north_and_north_in_the_south():
    jd = E.julian_day_at(_dt.datetime(2020, 6, 21), np.arange(0.0, 2.0, 1.0 / 1440.0))
    north = TT.body_direction(jd, 55.35, -160.5, 0.0, "sun")
    i = int(np.argmax(north["elevation_deg"]))
    assert abs(north["azimuth_deg"][i] - 180.0) < 1.0
    south = TT.body_direction(jd, -40.0, -160.5, 0.0, "sun")
    j = int(np.argmax(south["elevation_deg"]))
    az = south["azimuth_deg"][j]
    assert min(az, 360.0 - az) < 1.0


def test_elevation_and_azimuth_stay_in_range():
    jd = E.julian_day_at(T0, np.arange(0.0, 5.0, 1.0 / 720.0))
    for body in ("moon", "sun"):
        d = TT.body_direction(jd, 12.0, 77.0, 0.0, body)
        assert np.all(d["elevation_deg"] >= -90.0) and np.all(d["elevation_deg"] <= 90.0)
        assert np.all(d["azimuth_deg"] >= 0.0) and np.all(d["azimuth_deg"] < 360.0)


def test_body_azimuth_is_a_true_direction_not_an_axis():
    """It sweeps the full circle, so doubling it -- as a principal-axis bearing must
    be doubled -- would be wrong. The convention string says so and this proves it."""
    jd = E.julian_day_at(T0, np.arange(0.0, 3.0, 1.0 / 1440.0))
    d = TT.body_direction(jd, 20.0, 0.0, 0.0, "moon")
    assert d["azimuth_deg"].max() - d["azimuth_deg"].min() > 300.0
    assert "do NOT double" in d["convention"]


def test_p2_is_blind_to_the_sign_of_cos_zenith_and_p3_is_not():
    """Why lunar elevation is orthogonal to everything the program has bounded.

    The degree-2 potential goes as P2(cos z), which is EVEN, so it is identical for
    the Moon overhead and the Moon underfoot. Elevation distinguishes them. The
    degree-3 term is ODD and does too, which is what gives the elevation axis a
    physical carrier rather than merely a free coordinate.
    """
    cz = np.linspace(-1.0, 1.0, 4001)
    p2 = 0.5 * (3.0 * cz ** 2 - 1.0)
    p3 = 0.5 * (5.0 * cz ** 3 - 3.0 * cz)
    assert np.allclose(p2, p2[::-1])            # even
    assert np.allclose(p3, -p3[::-1])           # odd


def test_lunar_elevation_is_nearly_orthogonal_to_the_bounded_scalar():
    """The measurement that answers whether this axis is worth having at all."""
    from engine import residual as RS
    t = np.arange(0.0, 30.0, 1.0 / 1440.0)
    jd = E.julian_day_at(T0, t)
    scal = ST.site_scalar(jd, 55.35, -160.5, 0.0)
    d = TT.body_direction(jd, 55.35, -160.5, 0.0, "moon")
    assert RS.new_information_fraction(d["sin_elevation"], [scal]) > 0.99
    # and the solar one is NOT, so the test is not passing on a triviality
    ds = TT.body_direction(jd, 55.35, -160.5, 0.0, "sun")
    assert RS.new_information_fraction(ds["sin_elevation"], [scal]) < 0.75


def test_body_direction_rates_are_physically_sane():
    t = np.arange(0.0, 5.0, 1.0 / 1440.0)
    jd = E.julian_day_at(T0, t)
    r = TT.body_direction_rates(t, jd, 20.0, 0.0, 0.0, "moon")
    med = float(np.median(r["azimuth_rate_deg_per_h"]))
    assert 5.0 < med < 25.0, med          # ~360 deg per ~24.8 h lunar day
    assert np.all(np.isfinite(r["elevation_rate_deg_per_h"]))


def test_free_surface_coulomb_annihilates_at_mu_equals_cot_dip():
    """K-313: CFS = sin(d)(mu sin d - cos d) N(s+90) is IDENTICALLY ZERO at mu = cot(d).

    At the default friction 0.4 the annihilating dip is 68.2 degrees. A receiver there
    is scored with a Coulomb series that is exactly zero, so any statistic built on it
    is meaningless -- and it fails SILENTLY, because a zero series produces a
    well-formed null rather than an error.

    MEASURED CONSEQUENCE FOR WORK ALREADY COMMITTED: the fault-relative world scan used
    Slab2 INTERFACE dips, median 16.9 degrees, and only 96 of 7,004 events (1.37 %) fell
    within |geometry factor| < 0.05. So the defect is real and immaterial there. It is
    NOT immaterial for the vertical strike-slip populations K-127 needs, which is
    exactly where a future arm would put its receivers, and this test exists so that
    trap is caught before that arm is written rather than after.
    """
    import math
    st = TT.stress_tensor(JD[:64], LAT, LON, DEPTH)
    for mu in (0.2, 0.4, 0.6, 0.8):
        d_ann = math.degrees(math.atan(1.0 / mu))
        c = TT.coulomb(st, 250.0, d_ann, 90.0, mu)["coulomb_pa"]
        scale = np.abs(TT.normal_along_bearing(st["s_NN"], st["s_EE"], st["s_NE"],
                                               340.0)).max()
        assert np.abs(c).max() < 1e-12 * scale, (mu, d_ann, np.abs(c).max())
    # and away from the annihilating dip the series is emphatically non-zero
    c = TT.coulomb(st, 250.0, 20.0, 90.0, 0.4)["coulomb_pa"]
    assert np.abs(c).max() > 0.0


def test_coulomb_geometry_factor_matches_the_closed_form():
    """The factor itself, so a future reader can see where the zero comes from."""
    import math
    for dip in (5.0, 20.0, 45.0, 68.2, 80.0):
        for mu in (0.2, 0.4, 0.6):
            r = math.radians(dip)
            want = math.sin(r) * (mu * math.sin(r) - math.cos(r))
            st = TT.stress_tensor(JD[:16], LAT, LON, DEPTH)
            n_perp = TT.normal_along_bearing(st["s_NN"], st["s_EE"], st["s_NE"], 340.0)
            got = TT.coulomb(st, 250.0, dip, 90.0, mu)["coulomb_pa"]
            assert np.allclose(got, want * n_perp, rtol=1e-10,
                               atol=1e-9 * np.abs(n_perp).max())
