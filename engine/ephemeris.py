"""Zero-download astronomical ephemeris (Meeus-style low precision, plain numpy).

Everything the miner's family-1 features need is a deterministic function of UTC time:
solar and lunar longitude/declination, lunar distance, and the three lunar cycle
phases (synodic, anomalistic, draconic). No packages beyond numpy, no network.

Accuracy, stated honestly and tested in engine/tests/test_mine.py:

  * solar longitude   ~ 0.01 deg      (Meeus ch.25 "low accuracy" series)
  * lunar longitude   ~ 0.2  deg      (Meeus ch.47 truncated to the leading terms)
  * lunar distance    ~ 1000 km       (parallax series, 4 terms)
  * lunar declination ~ 0.3  deg

That is far finer than the one-day sampling of the engine's design, which is the only
resolution any of these features are ever used at: a 0.2 deg error in lunar longitude
is 0.016 days of synodic phase. These are NOT ephemerides for eclipse prediction.

All angles are returned in degrees unless the name says `_rad`. Time is the Julian Day
(TT is not distinguished from UTC: the ~70 s difference is 1e-5 of a day).
"""

from __future__ import annotations

import datetime as _dt

import numpy as np

J2000 = 2451545.0
DEG = np.pi / 180.0

# Mean cycle lengths, days (for documentation, tests and the period-scan ladder).
SYNODIC_MONTH = 29.530588
ANOMALISTIC_MONTH = 27.554550
DRACONIC_MONTH = 27.212221
TROPICAL_YEAR = 365.242190


def julian_day(dt: _dt.datetime) -> float:
    """Julian Day of a (timezone-aware or naive-UTC) datetime."""
    y, m = dt.year, dt.month
    d = (dt.day + (dt.hour + (dt.minute + dt.second / 60.0) / 60.0) / 24.0)
    if m <= 2:
        y, m = y - 1, m + 12
    a = y // 100
    b = 2 - a + a // 4                      # Gregorian calendar
    return (np.floor(365.25 * (y + 4716)) + np.floor(30.6001 * (m + 1))
            + d + b - 1524.5)


def day_grid(t0: _dt.datetime, n_days: int, hour: float = 12.0) -> np.ndarray:
    """Julian Days for day indices 0..n_days-1, sampled at `hour` UTC of each day.

    The engine's day `t` is the 24 h window [t0 + t, t0 + t + 1); a feature for that
    day is evaluated at its MIDPOINT (12:00 UTC) by default, which is the smallest
    possible representation error for a daily feature and, critically, uses nothing
    from any other day.
    """
    base = julian_day(t0.replace(hour=0, minute=0, second=0, microsecond=0))
    return base + hour / 24.0 + np.arange(int(n_days), dtype=np.float64)


def _wrap360(x):
    return np.mod(x, 360.0)


def obliquity_deg(jd):
    """Mean obliquity of the ecliptic (Meeus 22.2, truncated)."""
    t = (np.asarray(jd, dtype=np.float64) - J2000) / 36525.0
    return 23.439291 - 0.0130042 * t - 1.64e-7 * t * t + 5.04e-7 * t ** 3


def sun_position(jd):
    """Geometric solar position. Meeus ch.25, 'low accuracy' (~0.01 deg)."""
    jd = np.asarray(jd, dtype=np.float64)
    t = (jd - J2000) / 36525.0
    l0 = _wrap360(280.46646 + 36000.76983 * t + 0.0003032 * t * t)     # mean longitude
    m = _wrap360(357.52911 + 35999.05029 * t - 0.0001537 * t * t)      # mean anomaly
    mr = m * DEG
    c = ((1.914602 - 0.004817 * t - 1.4e-5 * t * t) * np.sin(mr)
         + (0.019993 - 0.000101 * t) * np.sin(2 * mr)
         + 0.000289 * np.sin(3 * mr))
    true_lon = _wrap360(l0 + c)
    e = 0.016708634 - 0.000042037 * t - 1.267e-7 * t * t
    v = mr + c * DEG
    r_au = 1.000001018 * (1 - e * e) / (1 + e * np.cos(v))
    eps = obliquity_deg(jd) * DEG
    lam = true_lon * DEG
    dec = np.arcsin(np.sin(eps) * np.sin(lam)) / DEG
    ra = _wrap360(np.arctan2(np.cos(eps) * np.sin(lam), np.cos(lam)) / DEG)
    return {"lon_deg": true_lon, "mean_anom_deg": m, "dec_deg": dec, "ra_deg": ra,
            "dist_au": r_au}


def moon_position(jd):
    """Geocentric lunar position. Meeus ch.47 truncated to the leading terms.

    Returns ecliptic longitude/latitude, distance (km), declination, and the four
    fundamental arguments (D, M, M', F) that the miner's cycle phases are built from.
    """
    jd = np.asarray(jd, dtype=np.float64)
    t = (jd - J2000) / 36525.0
    lp = _wrap360(218.3164477 + 481267.88123421 * t - 0.0015786 * t * t)   # L'
    d = _wrap360(297.8501921 + 445267.1114034 * t - 0.0018819 * t * t)     # elongation
    m = _wrap360(357.5291092 + 35999.0502909 * t - 0.0001536 * t * t)      # sun anomaly
    mp = _wrap360(134.9633964 + 477198.8675055 * t + 0.0087414 * t * t)    # moon anomaly
    f = _wrap360(93.2720950 + 483202.0175233 * t - 0.0036539 * t * t)      # arg latitude

    dr, mr, mpr, fr = d * DEG, m * DEG, mp * DEG, f * DEG

    lon = (lp
           + 6.288774 * np.sin(mpr)
           + 1.274027 * np.sin(2 * dr - mpr)
           + 0.658314 * np.sin(2 * dr)
           + 0.213618 * np.sin(2 * mpr)
           - 0.185116 * np.sin(mr)
           - 0.114332 * np.sin(2 * fr)
           + 0.058793 * np.sin(2 * dr - 2 * mpr)
           + 0.057066 * np.sin(2 * dr - mr - mpr)
           + 0.053322 * np.sin(2 * dr + mpr)
           + 0.045758 * np.sin(2 * dr - mr)
           - 0.040923 * np.sin(mr - mpr)
           - 0.034720 * np.sin(dr)
           - 0.030383 * np.sin(mr + mpr)
           + 0.015327 * np.sin(2 * dr - 2 * fr)
           - 0.012528 * np.sin(mpr + 2 * fr)
           + 0.010980 * np.sin(mpr - 2 * fr))
    lat = (5.128122 * np.sin(fr)
           + 0.280602 * np.sin(mpr + fr)
           + 0.277693 * np.sin(mpr - fr)
           + 0.173237 * np.sin(2 * dr - fr)
           + 0.055413 * np.sin(2 * dr - mpr + fr)
           + 0.046271 * np.sin(2 * dr - mpr - fr)
           + 0.032573 * np.sin(2 * dr + fr)
           + 0.017198 * np.sin(2 * mpr + fr)
           + 0.009266 * np.sin(2 * dr + mpr - fr))
    # equatorial horizontal parallax (deg) -> distance
    par = (0.950724
           + 0.051818 * np.cos(mpr)
           + 0.009531 * np.cos(2 * dr - mpr)
           + 0.007843 * np.cos(2 * dr)
           + 0.002824 * np.cos(2 * mpr)
           + 0.000857 * np.cos(2 * dr + mpr))
    dist_km = 6378.14 / np.sin(par * DEG)

    eps = obliquity_deg(jd) * DEG
    lam, bet = _wrap360(lon) * DEG, lat * DEG
    dec = np.arcsin(np.sin(bet) * np.cos(eps)
                    + np.cos(bet) * np.sin(eps) * np.sin(lam)) / DEG
    ra = _wrap360(np.arctan2(np.sin(lam) * np.cos(eps) - np.tan(bet) * np.sin(eps),
                             np.cos(lam)) / DEG)
    return {"lon_deg": _wrap360(lon), "lat_deg": lat, "dist_km": dist_km,
            "dec_deg": dec, "ra_deg": ra,
            "D_deg": d, "M_deg": m, "Mprime_deg": mp, "F_deg": f}


def ephemeris_table(t0: _dt.datetime, n_days: int, hour: float = 12.0):
    """All family-1 quantities for day indices 0..n_days-1, as a dict of arrays.

    Cycle phases are returned in RADIANS in [0, 2pi) and are the canonical
    definitions used by the miner:

      synodic      = moon longitude - sun longitude   (0 = new moon)
      anomalistic  = lunar mean anomaly M'            (0 = perigee)
      draconic     = argument of latitude F           (0 = ascending node)
      annual       = solar longitude                  (0 = March equinox)
    """
    jd = day_grid(t0, n_days, hour=hour)
    s = sun_position(jd)
    m = moon_position(jd)
    elong = np.arccos(np.clip(
        np.cos(m["lat_deg"] * DEG) * np.cos((m["lon_deg"] - s["lon_deg"]) * DEG),
        -1.0, 1.0)) / DEG
    return {
        "jd": jd,
        "sun_lon_deg": s["lon_deg"],
        "sun_dec_deg": s["dec_deg"],
        "sun_dist_au": s["dist_au"],
        "moon_lon_deg": m["lon_deg"],
        "moon_lat_deg": m["lat_deg"],
        "moon_dec_deg": m["dec_deg"],
        "moon_dist_km": m["dist_km"],
        "elongation_deg": elong,
        "synodic_rad": np.mod((m["lon_deg"] - s["lon_deg"]) * DEG, 2 * np.pi),
        "anomalistic_rad": np.mod(m["Mprime_deg"] * DEG, 2 * np.pi),
        "draconic_rad": np.mod(m["F_deg"] * DEG, 2 * np.pi),
        "annual_rad": np.mod(s["lon_deg"] * DEG, 2 * np.pi),
    }
