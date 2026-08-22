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
