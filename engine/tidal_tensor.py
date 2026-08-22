"""THE TIDAL STRAIN AND STRESS TENSOR, AND THE DIRECTIONAL AXES OFF IT.

WHY THIS MODULE EXISTS. Everything the program has computed from the tide until now
is a SCALAR read off a rank-2 field: `sitetide.areal_strain` (the declared
SCALAR_FOR_PHASE), `radial_disp_m`, and the observing application's vertical body-tide
displacement. Every one of them discards ORIENTATION. The tide is a tensor whose
principal axes sweep in bearing through the day; a fault has a fixed strike, dip and
rake. No scalar statistic can see the angle between them, so no scalar statistic the
program has ever run could have detected an orientation-dependent effect. D-1 falsified
a marginal, unconditional, scalar claim; it is silent on this class by construction.

`K092_SCALAR_PROVENANCE.md` section 5(1) already named the gap and marked it open:
"Coulomb-on-a-thrust still would [differ by a phase offset]; that translation remains
un-attempted and separately declared." This module attempts it.

WHAT IT COMPUTES

  1. The horizontal tidal STRAIN tensor in local (North, East) at the site, from the
     degree-2 potential and the IERS-2010 nominal Love/Shida pair (h2, l2). The Shida
     number l2 has never been used by this program for anything except the areal
     factor; it is what carries the horizontal information.
  2. The STRESS tensor under the free-surface plane-stress condition, in the F-016
     elastic convention (E = 75 GPa, nu = 0.25) reused verbatim so numbers stay
     comparable to everything else the program has published.
  3. RESOLVED COULOMB stress on a declared fault plane (strike, dip, rake, friction),
     Aki-Richards geometry, tension-positive.
  4. THE DIRECTIONAL AXES ("the rope"): the bearing of maximum horizontal extension,
     its rotation rate, and the traction resolved along an arbitrary FIXED bearing so
     that pull-then-ease along one direction is expressible as a signed history rather
     than as an instantaneous scalar.

THE VALIDATION THAT MAKES IT TRUSTWORTHY. The trace of the horizontal strain tensor
computed here must equal `sitetide.areal_strain` EXACTLY, because
areal = (l2/ga) * Laplacian_Omega(W) + 2 h2 W/ga and for a degree-2 surface harmonic
Laplacian_Omega(W) = -6 W, giving (2 h2 - 6 l2) W/ga = `sitetide.AREAL_FACTOR` * W/ga.
That identity is an independent derivation of a constant the program already trusts, so
if the tensor is wrong the identity breaks. `engine/tests/test_tidal_tensor.py` asserts
it, asserts the analytic angular derivatives against finite differences of
`sitetide.tidal_potential`, and asserts the Coulomb resolution against a hand-computed
case.

SCOPE FLAGS, WHICH ARE STRICTLY WIDER THAN sitetide's AND MUST TRAVEL WITH ANY NUMBER
FROM HERE. Everything in `sitetide.SCOPE_FLAGS` applies (body tide only, NO OCEAN
LOADING, spherical elastic laterally-homogeneous Earth, nominal degree-2 Love numbers
with no frequency dependence), and this module adds three of its own:

  * FREE-SURFACE PLANE STRESS. The stress tensor is formed under sigma_rr = sigma_rN =
    sigma_rE = 0, which is exact AT the free surface and approximate at hypocentral
    depth. Depth enters only through the (r/a)^2 potential scaling, as in sitetide. Any
    absolute Coulomb amplitude from this module is therefore NOT a measurement; the
    licensed use is a DETECTION statistic against an identically-constructed null,
    which is the same licence P7-22(a) grants the scalar path.
  * THE FAULT MECHANISM IS AN ASSUMPTION, NOT A DATUM. sitetide's flag says "no
    receiver fault mechanism, so no Coulomb resolution -- do not quote signs you cannot
    defend". This module does not lift that flag; it makes the mechanism an EXPLICIT,
    DECLARED input so the assumption is visible and can be gridded and priced. A
    Coulomb number here is conditional on the strike/dip/rake/friction handed in.
  * OCEAN LOADING IS LARGEST EXACTLY WHERE ORIENTATION MATTERS MOST. At a coastal
    subduction site the missing ocean load is both large and DIRECTIONAL (it is
    dominated by the coast-normal direction), so it is a worse omission for a bearing
    statistic than it is for a scalar one. Stated here so it cannot be forgotten
    downstream.

Nothing in this module is evidence. It is machinery.
"""

from __future__ import annotations

import numpy as np

from . import ephemeris as E
from . import sitetide as ST

DEG = np.pi / 180.0

# Reused verbatim from sitetide so the two paths cannot drift apart.
H2 = ST.H2
L2 = ST.L2
G_SURFACE = ST.G_SURFACE
R_MEAN = ST.R_MEAN
GM_MOON = ST.GM_MOON
GM_SUN = ST.GM_SUN
YOUNGS_MODULUS_PA = ST.YOUNGS_MODULUS_PA
POISSON_RATIO = ST.POISSON_RATIO

GA = G_SURFACE * R_MEAN          # the m^2/s^2 normalisation W is divided by

TENSOR_SCOPE_FLAGS = (
    "ALL OF sitetide.SCOPE_FLAGS, PLUS: (1) FREE-SURFACE PLANE STRESS -- the stress "
    "tensor is formed under sigma_rr = sigma_rN = sigma_rE = 0, exact at the free "
    "surface and approximate at hypocentral depth; no absolute Coulomb amplitude from "
    "this module is a measurement, only a detection statistic against an "
    "identically-constructed null. (2) THE FAULT MECHANISM IS A DECLARED ASSUMPTION, "
    "not a datum -- sitetide's 'no receiver mechanism, so no Coulomb resolution' flag "
    "is NOT lifted, it is made explicit so the assumption can be gridded and priced. "
    "(3) The absent ocean load is DIRECTIONAL and largest at coastal subduction "
    "sites, so it is a worse omission for a bearing statistic than for a scalar one.")

BEARING_CONVENTION = (
    "Bearings are AZIMUTHS IN DEGREES CLOCKWISE FROM NORTH, in [0, 180) because a "
    "principal axis is an undirected line and 'pull along 030' and 'pull along 210' "
    "are the same state. Any circular statistic on a bearing is therefore an AXIAL "
    "statistic on a pi-periodic variable and must be doubled before a standard "
    "circular test is applied -- see axial_to_circular().")


# --------------------------------------------------- the potential and its slopes --
def _body_terms(jd, phi_gc, lon, r_s, gm, ra_deg, dec_deg, dist_m):
    """W and its first and second angular derivatives at the site, for one body.

    W = K * P2(c) with c = cos(zenith angle) and K = (GM/d) (r/d)^2, exactly as
    `sitetide.tidal_potential` builds it. Derivatives are w.r.t. GEOCENTRIC LATITUDE
    phi and longitude lambda at fixed radius, taken analytically:

        c        = sin(phi) sin(dec) + cos(phi) cos(dec) cos(H),  H = GMST + lon - RA
        c_phi    = cos(phi) sin(dec) - sin(phi) cos(dec) cos(H)
        c_phiphi = -c
        c_lam    = -cos(phi) cos(dec) sin(H)
        c_lamlam = -cos(phi) cos(dec) cos(H)
        c_philam =  sin(phi) cos(dec) sin(H)

    and W = K(3c^2 - 1)/2 gives W_phi = 3 K c c_phi, W_phiphi = 3K(c_phi^2 - c^2),
    W_lam = 3 K c c_lam, W_lamlam = 3K(c_lam^2 + c c_lamlam),
    W_philam = 3K(c_phi c_lam + c c_philam).
    """
    ha = ST.gmst_rad(jd) + lon - np.asarray(ra_deg, dtype=np.float64) * DEG
    dec = np.asarray(dec_deg, dtype=np.float64) * DEG
    sphi, cphi = np.sin(phi_gc), np.cos(phi_gc)
    sdec, cdec = np.sin(dec), np.cos(dec)
    sha, cha = np.sin(ha), np.cos(ha)

    c = sphi * sdec + cphi * cdec * cha
    c_p = cphi * sdec - sphi * cdec * cha
    c_pp = -c
    c_l = -cphi * cdec * sha
    c_ll = -cphi * cdec * cha
    c_pl = sphi * cdec * sha

    k = (gm / dist_m) * (r_s / dist_m) ** 2
    return {
        "W": k * 0.5 * (3.0 * c * c - 1.0),
        "W_p": 3.0 * k * c * c_p,
        "W_pp": 3.0 * k * (c_p * c_p + c * c_pp),
        "W_l": 3.0 * k * c * c_l,
        "W_ll": 3.0 * k * (c_l * c_l + c * c_ll),
        "W_pl": 3.0 * k * (c_p * c_l + c * c_pl),
    }


def potential_and_derivatives(jd, lat_deg, lon_deg, depth_km=0.0,
                              bodies=("moon", "sun")):
    """W, W_phi, W_phiphi, W_lambda, W_lambdalambda, W_philambda summed over bodies."""
    jd = np.asarray(jd, dtype=np.float64)
    phi_gc, lon, r_s = ST.geocentric_site(lat_deg, lon_deg, depth_km)
    keys = ("W", "W_p", "W_pp", "W_l", "W_ll", "W_pl")
    shape = np.broadcast(jd, np.asarray(r_s)).shape
    out = {k: np.zeros(shape, dtype=np.float64) for k in keys}
    if "moon" in bodies:
        m = E.moon_position(jd)
        t = _body_terms(jd, phi_gc, lon, r_s, GM_MOON, m["ra_deg"], m["dec_deg"],
                        np.asarray(m["dist_km"], dtype=np.float64) * 1000.0)
        for k in keys:
            out[k] = out[k] + t[k]
    if "sun" in bodies:
        s = E.sun_position(jd)
        t = _body_terms(jd, phi_gc, lon, r_s, GM_SUN, s["ra_deg"], s["dec_deg"],
                        np.asarray(s["dist_au"], dtype=np.float64) * 1.495978707e11)
        for k in keys:
            out[k] = out[k] + t[k]
    out["phi_gc"] = phi_gc
    return out


# ------------------------------------------------------------- the strain tensor --
def strain_tensor(jd, lat_deg, lon_deg, depth_km=0.0, bodies=("moon", "sun")):
    """Horizontal tidal strain in local (North, East), EXTENSION POSITIVE.

    From the classical body-tide displacement field u_r = h2 W/g,
    u_theta = (l2/g) dW/dtheta, u_lambda = (l2/(g sin theta)) dW/dlambda, the surface
    strains are

        e_thth = (l2/ga) W_thth                 + (h2/ga) W
        e_lala = (l2/ga) [ W_lala / sin^2(th)
                           + cot(th) W_th ]      + (h2/ga) W
        e_thla = (l2/ga) [ W_thla / sin(th)
                           - cot(th) W_la / sin(th) ]

    Rewritten in geocentric LATITUDE phi (theta = pi/2 - phi, so d/dtheta = -d/dphi,
    d2/dtheta2 = d2/dphi2, sin theta = cos phi, cot theta = tan phi) and mapped to
    local North = -theta:

        e_NN = (l2/ga) W_pp + (h2/ga) W
        e_EE = (l2/ga) [ W_ll / cos^2(phi) - tan(phi) W_p ] + (h2/ga) W
        e_NE = (l2/ga) [ W_pl + tan(phi) W_l ] / cos(phi)

    `e_rr` is carried through from `sitetide` unchanged so the two modules cannot
    disagree about the vertical: it is the free-surface value -nu/(1-nu) * areal.
    """
    d = potential_and_derivatives(jd, lat_deg, lon_deg, depth_km, bodies)
    phi = d["phi_gc"]
    cphi = np.cos(phi)
    if np.any(np.abs(cphi) < 1e-6):
        raise ValueError("strain_tensor is singular at the geographic poles")
    tphi = np.tan(phi)

    e_nn = (L2 / GA) * d["W_pp"] + (H2 / GA) * d["W"]
    e_ee = (L2 / GA) * (d["W_ll"] / (cphi * cphi) - tphi * d["W_p"]) + (H2 / GA) * d["W"]
    e_ne = (L2 / GA) * (d["W_pl"] + tphi * d["W_l"]) / cphi

    areal = e_nn + e_ee
    e_rr = -(POISSON_RATIO / (1.0 - POISSON_RATIO)) * areal
    return {
        "e_NN": e_nn, "e_EE": e_ee, "e_NE": e_ne,
        "areal_strain": areal,
        "radial_strain": e_rr,
        "potential": d["W"],
        "sign_convention": "EXTENSION POSITIVE",
        "scope_flags": TENSOR_SCOPE_FLAGS,
    }


def stress_tensor(jd, lat_deg, lon_deg, depth_km=0.0, bodies=("moon", "sun")):
    """Horizontal stress in local (North, East), Pa, TENSION POSITIVE.

    Free-surface PLANE STRESS: sigma_rr = sigma_rN = sigma_rE = 0, so

        sigma_NN = E/(1-nu^2) (e_NN + nu e_EE)
        sigma_EE = E/(1-nu^2) (e_EE + nu e_NN)
        sigma_NE = E/(1+nu) e_NE

    with the F-016 elastic convention (E = 75 GPa, nu = 0.25) reused verbatim.
    """
    s = strain_tensor(jd, lat_deg, lon_deg, depth_km, bodies)
    ey, nu = YOUNGS_MODULUS_PA, POISSON_RATIO
    a = ey / (1.0 - nu * nu)
    out = dict(s)
    out.update({
        "s_NN": a * (s["e_NN"] + nu * s["e_EE"]),
        "s_EE": a * (s["e_EE"] + nu * s["e_NN"]),
        "s_NE": (ey / (1.0 + nu)) * s["e_NE"],
        "s_DD": np.zeros_like(s["e_NN"]),
        "stress_units": "Pa, TENSION POSITIVE, free-surface plane stress",
    })
    return out


# ------------------------------------------------------- the directional axes --
def principal_bearing(e_nn, e_ee, e_ne):
    """Bearing (deg CW from North, in [0, 180)) of MAXIMUM horizontal EXTENSION.

    Normal strain along azimuth alpha (unit vector (E, N) = (sin a, cos a)):

        e(a) = (e_EE + e_NN)/2 + cos(2a) (e_NN - e_EE)/2 + sin(2a) e_NE

    which is maximised at 2a = atan2(e_NE, (e_NN - e_EE)/2). Returned mod 180 because a
    principal axis is an undirected line: see BEARING_CONVENTION.
    """
    two_a = np.arctan2(np.asarray(e_ne, dtype=np.float64),
                       0.5 * (np.asarray(e_nn, dtype=np.float64)
                              - np.asarray(e_ee, dtype=np.float64)))
    return np.mod(0.5 * two_a / DEG, 180.0)


def normal_along_bearing(e_nn, e_ee, e_ne, bearing_deg):
    """Normal strain (or stress, same algebra) resolved along a FIXED bearing.

    This is the quantity Jim's rope image names: how hard the site is being pulled
    ALONG ONE DIRECTION, as a signed number with a sign history, rather than how hard
    it is being pulled in some direction-free average sense.
    """
    a = np.asarray(bearing_deg, dtype=np.float64) * DEG
    return (np.asarray(e_ee, dtype=np.float64) * np.sin(a) ** 2
            + 2.0 * np.asarray(e_ne, dtype=np.float64) * np.sin(a) * np.cos(a)
            + np.asarray(e_nn, dtype=np.float64) * np.cos(a) ** 2)


def shear_along_bearing(e_nn, e_ee, e_ne, bearing_deg):
    """Shear component in the frame whose first axis lies along `bearing_deg`.

    With p = (cos a, sin a) in (North, East) and q the perpendicular,

        tau(a) = p . T . q = -(e_NN - e_EE)/2 * sin(2a) + e_NE * cos(2a)

    The leading MINUS is not cosmetic: with the opposite sign the expression fails to
    vanish on the principal axes, which is the identity
    `test_shear_along_bearing_vanishes_on_principal_axes` exists to enforce. It caught
    exactly that error on this function's first version.
    """
    a = np.asarray(bearing_deg, dtype=np.float64) * DEG
    return (-0.5 * (np.asarray(e_nn, dtype=np.float64)
                    - np.asarray(e_ee, dtype=np.float64)) * np.sin(2.0 * a)
            + np.asarray(e_ne, dtype=np.float64) * np.cos(2.0 * a))


def axial_to_circular(bearing_deg):
    """Map an axial bearing in [0, 180) to a circular angle in [0, 2pi).

    Doubling is mandatory before any standard circular statistic: a pi-periodic
    variable tested with a 2pi-periodic statistic is a silently wrong test, and the
    program's Kuiper/Watson machinery is 2pi by construction.
    """
    return np.mod(2.0 * np.asarray(bearing_deg, dtype=np.float64) * DEG, 2.0 * np.pi)


def bearing_rotation_rate(t_days, bearing_deg):
    """d(bearing)/dt in degrees per hour, unwrapped on the DOUBLED angle.

    Jim's "travelling angle". Unwrapping has to happen on the doubled angle or every
    passage through 180 degrees registers as a spurious jump; the result is halved back
    into bearing units afterwards.
    """
    t = np.asarray(t_days, dtype=np.float64)
    two = np.unwrap(2.0 * np.asarray(bearing_deg, dtype=np.float64) * DEG)
    return np.gradient(0.5 * two / DEG, t * 24.0)


# --------------------------------------------------------------- fault geometry --
def fault_vectors(strike_deg, dip_deg, rake_deg):
    """Unit normal and slip vectors in (North, East, Down), Aki & Richards convention.

    strike is measured clockwise from North with the fault dipping to the RIGHT of the
    strike direction; dip from horizontal; rake in the fault plane from the strike
    direction, positive anticlockwise looking at the hanging wall (90 deg = pure
    thrust, -90 = pure normal, 0 = left-lateral strike slip).
    """
    s = np.asarray(strike_deg, dtype=np.float64) * DEG
    d = np.asarray(dip_deg, dtype=np.float64) * DEG
    r = np.asarray(rake_deg, dtype=np.float64) * DEG
    n = np.stack([-np.sin(d) * np.sin(s),
                  np.sin(d) * np.cos(s),
                  -np.cos(d)], axis=-1)
    u = np.stack([np.cos(r) * np.cos(s) + np.cos(d) * np.sin(r) * np.sin(s),
                  np.cos(r) * np.sin(s) - np.cos(d) * np.sin(r) * np.cos(s),
                  -np.sin(r) * np.sin(d)], axis=-1)
    return n, u


def coulomb(stress, strike_deg, dip_deg, rake_deg, friction=0.4):
    """Resolved Coulomb failure stress on a declared plane, Pa, from `stress_tensor`.

        CFS = tau + mu * sigma_n

    with tau the shear traction along the SLIP direction and sigma_n the normal
    traction, TENSION POSITIVE, so that unclamping raises CFS. The stress tensor is
    the free-surface plane-stress one, i.e. the only non-zero components are NN, EE
    and NE; the Down row and column are exactly zero.

    `friction` is a FREE PARAMETER and this function does not choose it. Anything that
    scans it owes the multiplicity, per the module docstring's second scope flag.
    """
    n, u = fault_vectors(strike_deg, dip_deg, rake_deg)
    s_nn, s_ee, s_ne = stress["s_NN"], stress["s_EE"], stress["s_NE"]
    # traction t = sigma . n, with sigma having zero Down row/column
    t_n = s_nn * n[..., 0] + s_ne * n[..., 1]
    t_e = s_ne * n[..., 0] + s_ee * n[..., 1]
    t_d = np.zeros_like(np.asarray(t_n))
    sigma_n = t_n * n[..., 0] + t_e * n[..., 1] + t_d * n[..., 2]
    tau = t_n * u[..., 0] + t_e * u[..., 1] + t_d * u[..., 2]
    return {
        "shear_pa": tau,
        "normal_pa": sigma_n,
        "coulomb_pa": tau + float(friction) * sigma_n,
        "geometry": {"strike_deg": float(np.asarray(strike_deg).ravel()[0]),
                     "dip_deg": float(np.asarray(dip_deg).ravel()[0]),
                     "rake_deg": float(np.asarray(rake_deg).ravel()[0]),
                     "friction": float(friction)},
        "convention": ("Aki-Richards (N, E, Down); TENSION POSITIVE so unclamping "
                       "raises CFS; free-surface plane stress"),
        "scope_flags": TENSOR_SCOPE_FLAGS,
    }


def site_directional(jd, lat_deg, lon_deg, depth_km=0.0, bodies=("moon", "sun"),
                     fault=None, friction=0.4, fixed_bearing_deg=None):
    """Everything directional at one site and arbitrary times, as one record."""
    st = stress_tensor(jd, lat_deg, lon_deg, depth_km, bodies)
    out = dict(st)
    out["bearing_max_extension_deg"] = principal_bearing(st["e_NN"], st["e_EE"],
                                                         st["e_NE"])
    out["bearing_convention"] = BEARING_CONVENTION
    if fixed_bearing_deg is not None:
        out["strain_along_fixed_bearing"] = normal_along_bearing(
            st["e_NN"], st["e_EE"], st["e_NE"], fixed_bearing_deg)
        out["stress_along_fixed_bearing_pa"] = normal_along_bearing(
            st["s_NN"], st["s_EE"], st["s_NE"], fixed_bearing_deg)
        out["fixed_bearing_deg"] = float(fixed_bearing_deg)
    if fault is not None:
        out["coulomb"] = coulomb(st, fault["strike_deg"], fault["dip_deg"],
                                 fault["rake_deg"], friction)
    return out


# ------------------------------------------------- WHERE THE BODY ACTUALLY IS --
# Jim's original question named "azimuth, elevation ... pull/push/angular/speed/
# directional". The BEARING functions above answer the stress half of that: which way
# the tide is pulling the ground. They do NOT answer the source half: where the Moon
# and Sun actually are in the sky from the epicentre.
#
# §P7-25(5) REFRAMED K-096 item 4 (body azimuth/elevation) into item 2's tidal-traction
# azimuth, on the grounds that it "is closer to the forcing than to the stress and
# invites an unconstrained space". That reframing may well be right, but it was an
# argument rather than a measurement, and `engine/residual.py` now makes it a
# measurement: compute the axis and read its new-information fraction against the axis
# D-1 already bounded. Kepler's proposed rule 6 says to audit the projection BEFORE
# proposing the axis, and that is what these functions exist to make possible.
#
# NOTHING HERE IS A STATISTIC. It is the coordinate, so the question can be settled by
# arithmetic instead of by preference.

def body_direction(jd, lat_deg, lon_deg, depth_km=0.0, body="moon"):
    """Azimuth and elevation of the Moon or Sun as seen from the site.

    Geometric and GEOCENTRIC: no refraction, no parallax correction, no light-time --
    the same approximations `sitetide` already declares, so these angles are consistent
    with the potential computed beside them rather than being a second, better
    ephemeris that would not correspond to it.

        sin(elevation) = sin(phi) sin(dec) + cos(phi) cos(dec) cos(H)
        azimuth        = atan2(-sin(H) cos(dec),
                               sin(dec) cos(phi) - cos(dec) sin(phi) cos(H))

    with H the local hour angle. Azimuth is degrees CLOCKWISE FROM NORTH in [0, 360)
    and is a TRUE direction, not an axial one -- unlike a principal-axis bearing it is
    2-pi periodic, so it must NOT be doubled before a circular statistic.
    """
    jd = np.asarray(jd, dtype=np.float64)
    phi_gc, lon, _ = ST.geocentric_site(lat_deg, lon_deg, depth_km)
    pos = E.moon_position(jd) if body == "moon" else E.sun_position(jd)
    ha = ST.gmst_rad(jd) + lon - np.asarray(pos["ra_deg"], dtype=np.float64) * DEG
    dec = np.asarray(pos["dec_deg"], dtype=np.float64) * DEG
    sin_el = (np.sin(phi_gc) * np.sin(dec)
              + np.cos(phi_gc) * np.cos(dec) * np.cos(ha))
    sin_el = np.clip(sin_el, -1.0, 1.0)
    az = np.arctan2(-np.sin(ha) * np.cos(dec),
                    np.sin(dec) * np.cos(phi_gc)
                    - np.cos(dec) * np.sin(phi_gc) * np.cos(ha))
    dist = (np.asarray(pos["dist_km"], dtype=np.float64) * 1000.0 if body == "moon"
            else np.asarray(pos["dist_au"], dtype=np.float64) * 1.495978707e11)
    return {
        "body": body,
        "elevation_deg": np.degrees(np.arcsin(sin_el)),
        "azimuth_deg": np.mod(np.degrees(az), 360.0),
        "sin_elevation": sin_el,                 # == cos(zenith), the potential's own
        "distance_m": dist,
        "hour_angle_rad": np.mod(ha, 2.0 * np.pi),
        "convention": ("azimuth degrees CLOCKWISE FROM NORTH in [0, 360), a TRUE "
                       "direction and 2-pi periodic -- do NOT double it the way an "
                       "axial principal-axis bearing must be doubled; elevation "
                       "degrees above the horizon, geometric and geocentric, no "
                       "refraction and no parallax"),
    }


def body_direction_rates(t_days, jd, lat_deg, lon_deg, depth_km=0.0, body="moon"):
    """d(elevation)/dt and d(azimuth)/dt in degrees per hour.

    Jim's "speed" and "travelling angle" applied to the SOURCE rather than to the
    stress. The azimuth derivative is unwrapped on the full 2-pi circle (not the
    doubled angle, because a body azimuth is a true direction and not an axis).
    """
    d = body_direction(jd, lat_deg, lon_deg, depth_km, body)
    t_h = np.asarray(t_days, dtype=np.float64) * 24.0
    az_unwrapped = np.degrees(np.unwrap(np.radians(d["azimuth_deg"])))
    return {
        "elevation_rate_deg_per_h": np.gradient(d["elevation_deg"], t_h),
        "azimuth_rate_deg_per_h": np.gradient(az_unwrapped, t_h),
    }
