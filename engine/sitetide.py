"""D-3 -- SOLID-EARTH (BODY) TIDE AT A SITE, degree 2, from the existing ephemeris.

MINING_CATALOG F1-30 Tier 1, and §K92-0(5) item 2: *"the Love-number closed form
(`sitetide.py`) -- no download, body tide only"*. Zero downloads, zero new
dependencies: the lunar and solar positions are already computed by
`engine/ephemeris.py` and this module does the site projection and the Love-number
response. Nothing here is evidence.

WHAT IT EMITS, STATED ONCE AND UNAMBIGUOUSLY
---------------------------------------------
The PRIMARY quantity is **AREAL STRAIN**, `e_areal = e_theta_theta + e_phi_phi`,
dimensionless, **EXTENSION POSITIVE**. Everything else in this module is derived from
it or reported beside it:

  * `potential`        degree-2 tidal potential W2, m^2/s^2 (the generating field)
  * `equilibrium_m`    W2 / g, m -- the equilibrium tide height (the ~0.36 m lunar
                       scale that D-3's acceptance test names)
  * `radial_disp_m`    h2 * W2 / g, m -- the vertical body-tide displacement
  * `areal_strain`     the PRIMARY scalar (see below for the formula)
  * `volumetric_strain` areal strain closed with the free-surface condition
  * `mean_stress_pa`   the F-016 stress convention, E = 75 GPa, nu = 0.25 (below)

`SCALAR_FOR_PHASE = "areal_strain"` and `engine/conventions_d.py` (D-0) attaches the
Tanaka phase origin to exactly that name. One scalar, one sign, one unit (S-9).

THE FORMULAE, AND WHY EACH CONSTANT IS THE ONE IT IS
-----------------------------------------------------
Degree-2 tidal potential of a body of mass-parameter GM at geocentric distance d, at a
site at geocentric radius r_s whose geocentric zenith angle to the body is z:

    W2 = (GM / d) * (r_s / d)^2 * P2(cos z),      P2(x) = (3 x^2 - 1) / 2

which is the leading term of the standard expansion of the perturbing potential
(Melchior, *The Tides of the Planet Earth*, 2nd ed., ch. 2; Agnew, "Earth Tides",
*Treatise on Geophysics* vol. 3, eq. 3-5). The Doodson constant follows from it:
D = (3/4) GM_moon a^2 / d_moon^3 = 2.628 m^2/s^2 at mean lunar distance, which
`doodson_constant()` computes rather than hard-codes so the test can check it.

LOVE / SHIDA NUMBERS -- nominal degree-2 elastic values, **IERS Conventions (2010),
Petit & Luzum (eds.), IERS Technical Note 36, chapter 6, Table 6.3**:

    h2 = 0.6078       (vertical displacement)
    l2 = 0.0847       (horizontal displacement -- Shida number)
    k2 = 0.30190      (potential / gravimetric; carried for completeness and used
                       only in the gravimetric factor, not in the strains)

These are the nominal values for the elastic Earth at degree 2; the IERS frequency
dependence (the FCN resonance in the diurnal band, the anelastic corrections) is NOT
applied. That is a stated approximation: it moves diurnal amplitudes by order 1 % and
phases by well under a degree, which is far below the ocean-loading error this build
is already declared unable to carry (§P7-22(a), K-090(c)).

SURFACE STRAINS from a degree-n potential, with u_r = (h_n/g) W_n and horizontal
displacement u_h = (l_n/g) * a * grad_h W_n (grad_h in angle / a):

    e_theta_theta = (1/a) (du_theta/dtheta + u_r)
    e_phi_phi     = (1/a) ((1/sin th) du_phi/dphi + u_theta cot th + u_r)
    e_areal = e_tt + e_pp = [ 2 h_n + l_n * Lambda ] W_n / (g a)

where Lambda is the angular Laplacian eigenvalue, Lambda = -n(n+1) = -6 at n = 2. So

    e_areal = (2 h2 - 6 l2) * W2 / (g a) = 0.7074 * W2 / (g a)

(Agnew, *Treatise on Geophysics* 3.06, eq. 44-46; Melchior ch. 6.) The combination
(2h2 - 6l2) is the standard degree-2 areal-strain factor, and writing it as a named
constant `AREAL_FACTOR` is deliberate -- it is the single number an implementation is
most likely to get wrong, and the unit test checks the amplitude it implies.

FREE SURFACE and the VOLUMETRIC / STRESS PROXY. At a traction-free surface
sigma_rr = 0, so e_rr = -(nu / (1 - nu)) * e_areal and

    e_vol = e_areal + e_rr = e_areal * (1 - 2 nu) / (1 - nu)

Under the same plane-stress closure the horizontal stress trace is
sigma_tt + sigma_pp = E / (1 - nu) * e_areal, and with sigma_rr = 0 the MEAN NORMAL
STRESS is

    sigma_m = E * e_areal / (3 (1 - nu))

**E = 75 GPa and nu = 0.25 are the F-016 protocol's declared convention** (see
`exp_f016_phase_clock_null.py` / COSO_FIG4C_PROTOCOL.md), reused verbatim so that a
number from this module can be compared with an F-016 number without a units
argument. At the lunar peak this gives |sigma_m| ~ 1.3 kPa, the familiar body-tide
stress scale.

WHAT IS DELIBERATELY *NOT* MODELLED -- SCOPE FLAGS, not caveats
----------------------------------------------------------------
1. **Ocean loading. Absent.** §P7-22(a) scopes K-090(c) precisely around this: a
   rotation-invariant DETECTION statistic against an identically-warped null is
   permitted; **an amplitude or an absolute psi from this module is not a
   measurement** in a coastal or subduction setting. `SCOPE_FLAGS` says so and every
   `site_tide` record carries it.
2. **No receiver mechanism**, so no Coulomb resolution. F1-30's Pit (M-006): *"do not
   quote signs you cannot defend."*
3. **Spherical, non-rotating, elastic, laterally homogeneous Earth**; geodetic
   latitude is converted to geocentric with the WGS-84 flattening, but the site radius
   uses the WGS-84 ellipsoid radius rather than a geoid.
4. **Ephemeris accuracy is `engine/ephemeris.py`'s** -- ~0.2 deg lunar longitude,
   ~1000 km lunar distance. In the semidiurnal band 0.2 deg of lunar longitude is
   ~24 s of M2 phase (0.2 deg of phase at 12.42 h), i.e. ~0.006 rad of tidal phase.
   Amplitude error from the distance term is ~3 * 1000/384400 ~ 0.8 %.
   `EPHEMERIS_ERROR_NOTE` states it; it is small compared with everything in 1-3.
5. **Depth** enters only through the (r_s/a)^2 scaling of the potential; the radial
   dependence of h2/l2 inside the Earth is NOT modelled. At 50 km the potential
   scaling is a 1.6 % effect and the neglected Love-number depth dependence is of the
   same order -- both far below (1).
"""

from __future__ import annotations

import datetime as _dt

import numpy as np

from . import ephemeris as E

DEG = np.pi / 180.0

# ------------------------------------------------------------------ constants --
# Gravitational parameters, m^3/s^2. IAU 2015 / DE430 nominal values.
GM_MOON = 4.9028001e12
GM_SUN = 1.32712440018e20

# WGS-84 ellipsoid + standard gravity.
A_EARTH = 6378137.0                 # equatorial radius, m
F_EARTH = 1.0 / 298.257223563       # flattening
G_SURFACE = 9.80665                 # standard gravity, m/s^2
R_MEAN = 6371000.0                  # mean radius, m (used for the g*a normalisation)

# Degree-2 Love/Shida numbers, nominal elastic, IERS Conventions (2010) TN36 Table 6.3.
H2 = 0.6078
L2 = 0.0847
K2 = 0.30190

# The degree-2 areal-strain factor 2 h2 - n(n+1) l2 with n = 2. Named because it is
# the single constant an implementation is most likely to get wrong.
AREAL_FACTOR = 2.0 * H2 - 6.0 * L2      # = 0.7074

# The gravimetric factor delta = 1 + h2 - 3/2 k2, reported for cross-checking against
# any published tidal-gravity number (nominal 1.1562).
GRAVIMETRIC_DELTA = 1.0 + H2 - 1.5 * K2

# F-016's declared elastic convention, reused verbatim so numbers are comparable.
YOUNGS_MODULUS_PA = 75.0e9
POISSON_RATIO = 0.25

# The scalar D-0's convention is attached to (see engine/conventions_d.py).
SCALAR_FOR_PHASE = "areal_strain"

# The mean periods of the four constituents the acceptance test looks for, in DAYS.
# Standard values (Doodson); used for the spectral check and for the D-2 band scan.
CONSTITUENT_PERIODS_DAYS = {
    "M2": 12.420601 / 24.0,      # principal lunar semidiurnal
    "S2": 12.000000 / 24.0,      # principal solar semidiurnal
    "N2": 12.658348 / 24.0,      # larger lunar elliptic semidiurnal
    "K1": 23.934470 / 24.0,      # lunisolar diurnal
    "O1": 25.819342 / 24.0,      # principal lunar diurnal
    "P1": 24.065890 / 24.0,      # principal solar diurnal
    "Mf": 13.660791,             # lunar fortnightly
    "Mm": 27.554550,             # lunar monthly (anomalistic)
    "Ssa": 182.621095,           # solar semiannual
}

SCOPE_FLAGS = (
    "BODY TIDE ONLY. No ocean loading (§P7-22(a) scopes K-090(c): a "
    "rotation-invariant DETECTION statistic against an identically-warped null is "
    "permitted; an amplitude or an absolute psi from this module is NOT a "
    "measurement in a coastal or subduction setting). No receiver fault mechanism, "
    "so no Coulomb resolution -- F1-30 Pit / M-006, 'do not quote signs you cannot "
    "defend'. Spherical elastic laterally-homogeneous Earth, nominal IERS-2010 "
    "degree-2 Love numbers with NO frequency dependence (no FCN resonance, no "
    "anelasticity). Depth enters only through the (r/a)^2 potential scaling.")

EPHEMERIS_ERROR_NOTE = (
    "Ephemeris is engine/ephemeris.py (Meeus truncated): ~0.2 deg lunar longitude, "
    "~1000 km lunar distance. That is ~0.006 rad of M2 tidal phase and ~0.8 % of "
    "amplitude -- both negligible beside the absent ocean loading.")


# ------------------------------------------------------------- site geometry --
def geocentric_site(lat_deg, lon_deg, depth_km=0.0):
    """(geocentric latitude rad, longitude rad, geocentric radius m) for a site.

    Geodetic -> geocentric latitude with the WGS-84 flattening:
    tan(phi_gc) = (1 - f)^2 tan(phi_gd). The radius is the ellipsoid radius at that
    latitude, minus the depth. Depth is in KM because that is the unit the catalogue
    carries it in (`mine.load_event_marks` -> `marks['depth']`).
    """
    lat = np.asarray(lat_deg, dtype=np.float64) * DEG
    lon = np.asarray(lon_deg, dtype=np.float64) * DEG
    dep = np.asarray(depth_km, dtype=np.float64) * 1000.0
    phi_gc = np.arctan((1.0 - F_EARTH) ** 2 * np.tan(lat))
    # ellipsoid radius at geodetic latitude (exact for an ellipse of revolution)
    b = A_EARTH * (1.0 - F_EARTH)
    r_ell = np.sqrt(((A_EARTH ** 2 * np.cos(lat)) ** 2 + (b ** 2 * np.sin(lat)) ** 2)
                    / ((A_EARTH * np.cos(lat)) ** 2 + (b * np.sin(lat)) ** 2))
    return phi_gc, lon, r_ell - dep


def gmst_rad(jd):
    """Greenwich Mean Sidereal Time, radians. IAU 1982 series (Meeus 12.4).

    theta = 280.46061837 + 360.98564736629 (JD - J2000)
            + 0.000387933 T^2 - T^3 / 38710000,  degrees.
    """
    jd = np.asarray(jd, dtype=np.float64)
    d = jd - E.J2000
    t = d / 36525.0
    th = (280.46061837 + 360.98564736629 * d
          + 0.000387933 * t * t - t ** 3 / 38710000.0)
    return np.mod(th, 360.0) * DEG


def _cos_zenith(phi_gc, lon, jd, ra_deg, dec_deg):
    """cos of the geocentric zenith angle from the site to a body."""
    ha = gmst_rad(jd) + lon - np.asarray(ra_deg, dtype=np.float64) * DEG
    dec = np.asarray(dec_deg, dtype=np.float64) * DEG
    return (np.sin(phi_gc) * np.sin(dec)
            + np.cos(phi_gc) * np.cos(dec) * np.cos(ha))


def _p2(x):
    return 0.5 * (3.0 * np.asarray(x, dtype=np.float64) ** 2 - 1.0)


def doodson_constant(d_moon_m=3.844e8):
    """(3/4) GM_moon a^2 / d^3 -- 2.628 m^2/s^2 at mean lunar distance.

    Computed, never hard-coded: this is the amplitude scale the acceptance test
    checks the whole potential path against.
    """
    return 0.75 * GM_MOON * R_MEAN ** 2 / float(d_moon_m) ** 3


# ---------------------------------------------------------------- the potential --
def tidal_potential(jd, lat_deg, lon_deg, depth_km=0.0, bodies=("moon", "sun")):
    """Degree-2 tidal potential W2 at the site, m^2/s^2, summed over bodies.

    Returns a dict with the total and the per-body contributions, because the
    lunar/solar split is what makes the M2-vs-S2 spectral check interpretable.
    """
    jd = np.asarray(jd, dtype=np.float64)
    phi_gc, lon, r_s = geocentric_site(lat_deg, lon_deg, depth_km)
    out = {}
    total = np.zeros(np.broadcast(jd, np.asarray(r_s)).shape, dtype=np.float64)
    if "moon" in bodies:
        m = E.moon_position(jd)
        d = np.asarray(m["dist_km"], dtype=np.float64) * 1000.0
        cz = _cos_zenith(phi_gc, lon, jd, m["ra_deg"], m["dec_deg"])
        w = (GM_MOON / d) * (r_s / d) ** 2 * _p2(cz)
        out["moon"] = w
        total = total + w
    if "sun" in bodies:
        s = E.sun_position(jd)
        d = np.asarray(s["dist_au"], dtype=np.float64) * 1.495978707e11
        cz = _cos_zenith(phi_gc, lon, jd, s["ra_deg"], s["dec_deg"])
        w = (GM_SUN / d) * (r_s / d) ** 2 * _p2(cz)
        out["sun"] = w
        total = total + w
    out["total"] = total
    return out


# ------------------------------------------------------------ the site record --
def site_tide(jd, lat_deg, lon_deg, depth_km=0.0, bodies=("moon", "sun")):
    """Every declared body-tide quantity at (lat, lon, depth) and arbitrary times.

    `jd` may be any shape; all outputs share it. Times are Julian Days -- use
    `engine.ephemeris.julian_day_at(t0, day_float)` to get them from the catalogue's
    own sub-daily `day_float`, which is what makes this an EVENT-TIME quantity and not
    a day-binned one (F9-10 / §P7-3(3)).
    """
    w = tidal_potential(jd, lat_deg, lon_deg, depth_km, bodies)
    w2 = w["total"]
    ga = G_SURFACE * R_MEAN
    areal = AREAL_FACTOR * w2 / ga
    e_rr = -(POISSON_RATIO / (1.0 - POISSON_RATIO)) * areal
    vol = areal + e_rr
    sigma_m = YOUNGS_MODULUS_PA * areal / (3.0 * (1.0 - POISSON_RATIO))
    return {
        "jd": np.asarray(jd, dtype=np.float64),
        "potential": w2,
        "potential_moon": w.get("moon"),
        "potential_sun": w.get("sun"),
        "equilibrium_m": w2 / G_SURFACE,
        "radial_disp_m": H2 * w2 / G_SURFACE,
        "areal_strain": areal,
        "radial_strain": e_rr,
        "volumetric_strain": vol,
        "mean_stress_pa": sigma_m,
        "scalar_for_phase": SCALAR_FOR_PHASE,
        "sign_convention": "EXTENSION POSITIVE",
        "elastic_convention": {"E_Pa": YOUNGS_MODULUS_PA, "nu": POISSON_RATIO,
                               "source": "F-016 protocol convention, reused verbatim"},
        "love_numbers": {"h2": H2, "l2": L2, "k2": K2,
                         "source": "IERS Conventions (2010) TN36 Table 6.3, nominal "
                                   "elastic degree 2",
                         "areal_factor_2h2_minus_6l2": AREAL_FACTOR,
                         "gravimetric_delta": GRAVIMETRIC_DELTA},
        "scope_flags": SCOPE_FLAGS,
        "ephemeris_note": EPHEMERIS_ERROR_NOTE,
    }


def site_scalar(jd, lat_deg, lon_deg, depth_km=0.0, which=SCALAR_FOR_PHASE,
                bodies=("moon", "sun")):
    """Just the declared scalar as a bare array -- the hot path for phase work."""
    rec = site_tide(jd, lat_deg, lon_deg, depth_km, bodies)
    if which not in rec:
        raise KeyError("unknown site-tide scalar %r" % (which,))
    return np.asarray(rec[which], dtype=np.float64)


def site_scalar_at(t0: _dt.datetime, day_float, lat_deg, lon_deg, depth_km=0.0,
                   which=SCALAR_FOR_PHASE, bodies=("moon", "sun")):
    """The declared scalar at the catalogue's own continuous `day_float` times.

    Routed through `ephemeris.julian_day_at` -- the one function §P7-3(3) names as
    the difference between a genuinely sub-daily feature and a day-binned one
    evaluated at an event.
    """
    jd = E.julian_day_at(t0, day_float)
    return site_scalar(jd, lat_deg, lon_deg, depth_km, which=which, bodies=bodies)


# --------------------------------------------------- D-0's phase, in the code --
def tidal_maxima(t0: _dt.datetime, day_lo, day_hi, lat_deg, lon_deg, depth_km=0.0,
                 which=SCALAR_FOR_PHASE, bodies=("moon", "sun"),
                 grid_minutes=5.0):
    """Sub-sample-refined TIMES of the site scalar's local maxima over [lo, hi], days.

    Split out from `tanaka_phase` because the maxima are a property of the SITE AND
    THE SPAN ALONE: every catalogue, every ETAS null replicate and every planted arm
    on the same site share them. Computing them once and phasing many event sets
    against them is what keeps the D-2/D-4 null ensembles affordable, and -- more
    importantly -- it makes the phase map BIT-IDENTICAL across signal and null, which
    is §P7-22(a)'s common-mode requirement taken literally rather than trusted.
    """
    dt_days = float(grid_minutes) / (60.0 * 24.0)
    lo, hi = float(day_lo), float(day_hi)
    n = int(np.ceil((hi - lo) / dt_days)) + 1
    tg = lo + dt_days * np.arange(n, dtype=np.float64)
    s = site_scalar_at(t0, tg, lat_deg, lon_deg, depth_km, which=which, bodies=bodies)
    i = np.nonzero((s[1:-1] > s[:-2]) & (s[1:-1] >= s[2:]))[0] + 1
    if i.size < 2:
        raise ValueError("fewer than two tidal maxima on the grid; widen the span")
    # parabolic sub-sample refinement of each maximum's TIME: error O(dt^2/T), not
    # O(dt/T) -- ~1e-4 rad at 5 minutes, four orders below the ephemeris error.
    y0, y1, y2 = s[i - 1], s[i], s[i + 1]
    den = (y0 - 2.0 * y1 + y2)
    shift = np.where(den != 0.0, 0.5 * (y0 - y2) / np.where(den == 0, 1.0, den), 0.0)
    return tg[i] + np.clip(shift, -0.5, 0.5) * dt_days


def phase_from_maxima(t_max, day_float):
    """D-0's phase for event times against precomputed maxima. NaN outside the span."""
    tm = np.asarray(t_max, dtype=np.float64)
    d = np.asarray(day_float, dtype=np.float64)
    k = np.searchsorted(tm, d, side="right") - 1
    ok = (k >= 0) & (k < tm.size - 1)
    ph = np.full(d.shape, np.nan, dtype=np.float64)
    kk = k[ok]
    t0m, t1m = tm[kk], tm[kk + 1]
    ph[ok] = 2.0 * np.pi * (d[ok] - t0m) / np.maximum(t1m - t0m, 1e-12)
    return np.mod(ph, 2.0 * np.pi)


def tanaka_phase(t0: _dt.datetime, day_float, lat_deg, lon_deg, depth_km=0.0,
                 which=SCALAR_FOR_PHASE, bodies=("moon", "sun"),
                 grid_minutes=5.0, pad_days=1.0):
    """Phase in [0, 2pi) with 0 at the local MAXIMUM of the scalar (D-0 / Tanaka).

    Construction, exactly as `engine/conventions_d.py` declares it and as Tanaka
    (2002) defines it: evaluate the scalar on a fine uniform time grid spanning the
    events (plus a pad), locate the interior local maxima, and for each event LINEARLY
    INTERPOLATE IN TIME between the bracketing maxima:

        phase = 2 pi * (t - t_max_before) / (t_max_after - t_max_before)

    So phase 0 and 2pi are successive maxima and phase pi is (to within the waveform's
    asymmetry) the intervening minimum.

    `grid_minutes = 5` resolves the ~12.4 h semidiurnal cycle with ~149 samples; the
    extremum locations are refined by a parabolic fit through the three grid samples
    around each discrete maximum, so the phase error from the grid is O(dt^2 / T) and
    not O(dt / T) -- about 1e-4 rad at 5 minutes, four orders below the ~0.006 rad
    ephemeris error. `phase_grid_error_rad` is returned so the number is stated rather
    than assumed.

    Events outside the interior maxima (before the first or after the last) get NaN;
    `pad_days` exists to make that set empty in normal use.
    """
    d = np.asarray(day_float, dtype=np.float64)
    if d.size == 0:
        return np.zeros(0, dtype=np.float64), {"n_maxima": 0, "n_nan": 0}
    dt_days = float(grid_minutes) / (60.0 * 24.0)
    t_max = tidal_maxima(t0, float(d.min()) - float(pad_days),
                         float(d.max()) + float(pad_days), lat_deg, lon_deg,
                         depth_km, which=which, bodies=bodies,
                         grid_minutes=grid_minutes)
    ph = phase_from_maxima(t_max, d)
    ok = np.isfinite(ph)
    info = {
        "convention": "D0-tanaka-stressmax-v1 (phase 0 = local scalar MAXIMUM)",
        "scalar": which,
        "n_maxima": int(t_max.size),
        "n_nan": int((~ok).sum()),
        "grid_minutes": float(grid_minutes),
        "mean_cycle_hours": float(np.mean(np.diff(t_max)) * 24.0),
        "phase_grid_error_rad": float(2.0 * np.pi * (dt_days ** 2)
                                      / max(float(np.mean(np.diff(t_max))), 1e-12)),
        "scope_flags": SCOPE_FLAGS,
    }
    return ph, info


def constituent_phase(t0: _dt.datetime, day_float, period_days, phase0=0.0):
    """Phase at a FIXED constituent frequency -- the per-band axis D-2 scans.

    This is NOT the D-0 convention (which is waveform-referenced); it is the
    single-frequency clock a band-by-band variance measurement needs, and it is kept
    in a separate function so the two can never be confused at a call site.
    """
    d = np.asarray(day_float, dtype=np.float64)
    return np.mod(2.0 * np.pi * d / float(period_days) + float(phase0), 2.0 * np.pi)


def constituent_spectrum(t0: _dt.datetime, lat_deg, lon_deg, depth_km=0.0,
                         n_days=365.0, sample_minutes=30.0,
                         which=SCALAR_FOR_PHASE, bodies=("moon", "sun")):
    """Amplitude spectrum of the site scalar over `n_days`, cycles/day.

    D-3's acceptance test: M2, S2, K1 and O1 must be present as resolved peaks.
    Returns (freq_cpd, amplitude), amplitude being 2|X_k|/N so a pure cosine of
    amplitude A reads A.
    """
    dt_days = float(sample_minutes) / (60.0 * 24.0)
    n = int(round(float(n_days) / dt_days))
    t = dt_days * np.arange(n, dtype=np.float64)
    s = site_scalar_at(t0, t, lat_deg, lon_deg, depth_km, which=which, bodies=bodies)
    s = s - s.mean()
    X = np.fft.rfft(s)
    amp = 2.0 * np.abs(X) / n
    freq = np.fft.rfftfreq(n, d=dt_days)
    return freq, amp


def peak_near(freq, amp, period_days, half_width_cpd=0.02):
    """(peak frequency, peak amplitude) in a window around 1/period. Test helper."""
    f0 = 1.0 / float(period_days)
    m = np.abs(freq - f0) <= float(half_width_cpd)
    if not m.any():
        return float("nan"), 0.0
    j = np.nonzero(m)[0]
    k = j[int(np.argmax(amp[j]))]
    return float(freq[k]), float(amp[k])


if __name__ == "__main__":          # pragma: no cover - operator convenience
    import json
    t0 = _dt.datetime(2020, 1, 1)
    rec = site_tide(E.julian_day_at(t0, np.arange(0.0, 2.0, 0.01)), 45.0, -125.0, 20.0)
    print(json.dumps({
        "doodson_constant_m2_s2": doodson_constant(),
        "areal_factor": AREAL_FACTOR,
        "gravimetric_delta": GRAVIMETRIC_DELTA,
        "peak_areal_strain": float(np.max(np.abs(rec["areal_strain"]))),
        "peak_mean_stress_Pa": float(np.max(np.abs(rec["mean_stress_pa"]))),
        "peak_radial_disp_m": float(np.max(np.abs(rec["radial_disp_m"]))),
        "peak_equilibrium_m": float(np.max(np.abs(rec["equilibrium_m"]))),
    }, indent=2))
