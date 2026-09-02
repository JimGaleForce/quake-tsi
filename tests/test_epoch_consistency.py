"""EPOCH INVARIANT regression test (supervisor correction, 2026-09-02).

THE DEFECT THIS PINS DOWN
-------------------------
The shared feature builders -- `exp_world_harmonics.features` and
`exp_mass_screen.raw_series` -- convert an event time column to a Julian date with

    jd = t_days + W.UNIX_EPOCH_JD          # UNIX_EPOCH_JD = JD of 1970-01-01T00:00Z

That conversion is only correct if `t_days` is DAYS SINCE THE UNIX EPOCH.
`exp_world_harmonics.load_region` always satisfied that. Four other loaders --
`exp_highn.load_zenodo`, `exp_mass_screen.load_region_full`, `exp_fluid_driven.load`
and `exp_nearcritical.load` -- returned days since `W.SPAN_START` (1995-01-01)
instead, so every tidal/lunar feature in the arms they feed was evaluated
9,131 days (25.0 yr) BEFORE the event it belonged to.

The invariant is stated ONCE here and asserted for EVERY loader, so a new loader that
gets it wrong fails this file rather than silently returning a null against a
scrambled lunar covariate.
"""

from __future__ import annotations

import datetime as _dt
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import exp_fluid_driven as FD          # noqa: E402
import exp_highn as HN                 # noqa: E402
import exp_mass_screen as MS           # noqa: E402
import exp_nearcritical as NC          # noqa: E402
import exp_world_harmonics as W        # noqa: E402

UNIX_EPOCH = _dt.datetime(1970, 1, 1, tzinfo=_dt.timezone.utc)
TOL = 1e-6                              # days

# A common fixture instant, comfortably inside every arm exploration split
# (the earliest cutoff is exp_nearcritical / exp_fluid_driven 70 % of their own span,
# and exp_world_harmonics.explore_cutoff() is 2017-02-14).
COMMON = _dt.datetime(2012, 3, 4, 5, 6, 7, tzinfo=_dt.timezone.utc)


def _known_jd(ts):
    """Julian date of a UTC datetime, computed independently of the code under test.

    Fliegel-Van Flandern Gregorian-calendar JDN, so this is NOT a restatement of
    `t + UNIX_EPOCH_JD`; it is an external check on it.
    """
    y, m, d = ts.year, ts.month, ts.day
    a = (14 - m) // 12
    yy = y + 4800 - a
    mm = m + 12 * a - 3
    jdn = d + (153 * mm + 2) // 5 + 365 * yy + yy // 4 - yy // 100 + yy // 400 - 32045
    frac = (ts.hour - 12) / 24.0 + ts.minute / 1440.0 + ts.second / 86400.0
    return jdn + frac


def _zenodo_line(ts, lat, lon, dep, mag, eid):
    return "%d %d %d %d %d %.3f %.4f %.4f %.2f %.2f %s\n" % (
        ts.year, ts.month, ts.day, ts.hour, ts.minute,
        ts.second + ts.microsecond / 1e6, lat, lon, dep, mag, eid)


def _csv_text(rows, header="time,latitude,longitude,depth,mag"):
    out = [header]
    for ts, lat, lon, dep, mag in rows:
        out.append("%s,%.4f,%.4f,%.2f,%.2f"
                   % (ts.strftime("%Y-%m-%dT%H:%M:%S.000Z"), lat, lon, dep, mag))
    return "\n".join(out) + "\n"


# Three widely separated, widely spaced sites so no arm declustering rule
# (30 days AND 150 km of an equal-or-larger prior event) removes any of them.
SITES = ((33.0, -116.0), (36.5, -121.0), (40.0, -124.5))


def _three(ts_list, mags):
    return [(ts, SITES[i][0], SITES[i][1], 5.0, mags[i])
            for i, ts in enumerate(ts_list)]


# --------------------------------------------------------------------------- fixtures
# Each loader gets dates matching the catalogue it actually reads, so the "first event
# ISO date" sanity check in each arm console output is exercised by the same values.
TIMES = {
    # QTM starts 2008 (exp_highn / exp_learned / exp_diurnal_discriminator)
    "zenodo": [_dt.datetime(2008, 1, 2, 3, 4, 5, tzinfo=_dt.timezone.utc),
               COMMON,
               _dt.datetime(2016, 6, 30, 23, 59, 59, tzinfo=_dt.timezone.utc)],
    # ComCat world regions, M >= 5 (exp_mass_screen)
    "region_full": [_dt.datetime(1996, 5, 6, 7, 8, 9, tzinfo=_dt.timezone.utc),
                    COMMON,
                    _dt.datetime(2016, 11, 12, 13, 14, 15, tzinfo=_dt.timezone.utc)],
    # k034 / fluid-driven regional catalogues start 1985
    "fluid": [_dt.datetime(1985, 7, 8, 9, 10, 11, tzinfo=_dt.timezone.utc),
              COMMON,
              _dt.datetime(2026, 2, 3, 4, 5, 6, tzinfo=_dt.timezone.utc)],
    # comcat_socal_m25 (exp_nearcritical)
    "nearcritical": [_dt.datetime(2010, 9, 8, 7, 6, 5, tzinfo=_dt.timezone.utc),
                     COMMON,
                     _dt.datetime(2020, 4, 5, 6, 7, 8, tzinfo=_dt.timezone.utc)],
}


def _load_zenodo(tmp_path, monkeypatch=None):
    p = tmp_path / "zen.txt"
    p.write_text("".join(_zenodo_line(ts, SITES[i][0], SITES[i][1], 5.0, 2.0 + i,
                                      "e%d" % i)
                         for i, ts in enumerate(TIMES["zenodo"])), encoding="utf-8")
    return HN.load_zenodo(p)[0]


def _load_region_full(tmp_path, monkeypatch=None):
    p = tmp_path / "region.csv"
    p.write_text(_csv_text(_three(TIMES["region_full"], (5.5, 5.6, 5.7))),
                 encoding="utf-8")
    return MS.load_region_full(p)[0]


def _load_fluid(tmp_path, monkeypatch=None):
    p = tmp_path / "fluid.csv"
    p.write_text(_csv_text(_three(TIMES["fluid"], (3.0, 3.1, 3.2))), encoding="utf-8")
    return FD.load(p)[0]


def _load_nearcritical(tmp_path, monkeypatch=None):
    p = tmp_path / "socal.csv"
    p.write_text(_csv_text(_three(TIMES["nearcritical"], (3.0, 3.1, 3.2))),
                 encoding="utf-8")
    monkeypatch.setattr(NC, "CATALOG", p)
    return NC.load()[0]


def _load_world(tmp_path):
    p = tmp_path / "world.csv"
    p.write_text(_csv_text(_three(TIMES["region_full"], (5.5, 5.6, 5.7))),
                 encoding="utf-8")
    return W.load_region(p)[0]


LOADERS = {
    "exp_fluid_driven.load": ("fluid", _load_fluid),
    "exp_highn.load_zenodo": ("zenodo", _load_zenodo),
    "exp_mass_screen.load_region_full": ("region_full", _load_region_full),
    "exp_nearcritical.load": ("nearcritical", _load_nearcritical),
}


@pytest.mark.parametrize("name", sorted(LOADERS))
def test_loader_returns_days_since_unix_epoch(name, tmp_path, monkeypatch):
    """`t + W.UNIX_EPOCH_JD` must be the TRUE Julian date of the fixture instant."""
    key, fn = LOADERS[name]
    t = fn(tmp_path, monkeypatch)
    assert t.size >= 1, "%s dropped the whole fixture" % name

    # The first row survives every arm magnitude floor, cutoff and declustering rule.
    want = _known_jd(TIMES[key][0])
    got = float(np.min(t)) + W.UNIX_EPOCH_JD
    assert got == pytest.approx(want, abs=TOL), (
        "%s: JD off by %.3f d -- the epoch invariant is broken" % (name, got - want))


@pytest.mark.parametrize("name", sorted(LOADERS))
def test_loader_agrees_with_load_region_on_a_common_instant(name, tmp_path,
                                                            monkeypatch):
    """All four loaders and `W.load_region` must place the SAME instant identically."""
    key, fn = LOADERS[name]
    t = fn(tmp_path, monkeypatch)
    idx = TIMES[key].index(COMMON)
    mine = float(np.sort(t)[idx])

    tw = _load_world(tmp_path)
    theirs = float(np.sort(tw)[TIMES["region_full"].index(COMMON)])

    assert mine == pytest.approx(theirs, abs=TOL), (
        "%s disagrees with exp_world_harmonics.load_region by %.3f d"
        % (name, mine - theirs))
    assert mine == pytest.approx(
        (COMMON - UNIX_EPOCH).total_seconds() / 86400.0, abs=TOL)


def test_the_9131_day_defect_would_be_caught():
    """The exact failure mode: SPAN_START-based days are 9131 d short of the invariant."""
    off = (W.SPAN_START - UNIX_EPOCH).total_seconds() / 86400.0
    assert off == pytest.approx(9131.0, abs=TOL)
    bad = (COMMON - W.SPAN_START).total_seconds() / 86400.0
    assert bad + W.UNIX_EPOCH_JD == pytest.approx(_known_jd(COMMON) - 9131.0, abs=TOL)


def test_unix_epoch_jd_constant_is_right():
    assert W.UNIX_EPOCH_JD == pytest.approx(_known_jd(UNIX_EPOCH), abs=TOL)


def test_day_of_year_anchor_is_a_january_first():
    """`exp_learned.build_features` takes `np.mod(t_days, 365.2422)` as day-of-year.

    That is only meaningful if t = 0 is a 1 January. It is, once the invariant holds.
    """
    assert UNIX_EPOCH.month == 1 and UNIX_EPOCH.day == 1
