"""Acceptance tests for engine/slab2.py.

Validated against KNOWN TRENCH GEOMETRY rather than against a reimplementation: the
Aleutian, Japan, Sunda and Peru-Chile trenches have strikes that are textbook facts,
and a reader that got the longitude convention or the row/column order wrong would not
reproduce all four.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from engine import slab2

pytestmark = pytest.mark.skipif(
    not os.path.isdir(slab2.SLAB2_DIR),
    reason="Slab2 distribution not present (data/ is gitignored)")


def test_all_27_slabs_are_present():
    codes = slab2.available_slabs()
    assert len(codes) == 27, codes
    for expect in ("alu", "sam", "kur", "sum", "cas", "hel"):
        assert expect in codes


@pytest.mark.parametrize("name,lat,lon,depth,strike_lo,strike_hi,code", [
    # Aleutian megathrust at Sand Point: arc strikes WSW
    ("sandpoint", 54.6, -160.5, 40.0, 240.0, 262.0, "alu"),
    # Japan Trench off Tohoku: strikes roughly N-S
    ("tohoku", 38.3, 142.4, 25.0, 180.0, 205.0, "kur"),
    # Sunda Trench off Sumatra: strikes NW-SE
    ("sumatra", 3.3, 95.9, 30.0, 295.0, 325.0, "sum"),
    # Peru-Chile Trench off central Chile: strikes roughly N-S
    ("chile", -33.0, -72.5, 30.0, 350.0, 20.0, "sam"),
])
def test_known_trench_geometry(name, lat, lon, depth, strike_lo, strike_hi, code):
    a = slab2.assign([lat], [lon], [depth])
    assert a["assigned"][0], name
    assert a["slab_code"][0] == code
    s = a["strike_deg"][0]
    if strike_lo <= strike_hi:
        assert strike_lo <= s <= strike_hi, (name, s)
    else:                                    # wraps through 0
        assert s >= strike_lo or s <= strike_hi, (name, s)
    assert 0.0 < a["dip_deg"][0] < 60.0, (name, a["dip_deg"][0])


def test_the_declared_alaska_geometry_was_right():
    """D-1c declared strike 250 / dip 20 from published knowledge, before Slab2.

    Slab2 puts the Sand Point interface at strike ~251 and dip ~15. Recorded as a
    standing check that the hand-declared geometry those arms used was not far off.
    """
    a = slab2.assign([54.6], [-160.5], [40.0])
    assert abs(a["strike_deg"][0] - 250.0) < 6.0
    assert abs(a["dip_deg"][0] - 20.0) < 8.0


def test_longitude_convention_is_handled():
    """Slab2 grids are 0..360, ComCat is -180..180. Easiest thing to get wrong here."""
    a = slab2.assign([54.6], [-160.5], [40.0])
    b = slab2.assign([54.6], [199.5], [40.0])
    assert a["strike_deg"][0] == b["strike_deg"][0]
    assert a["dip_deg"][0] == b["dip_deg"][0]


def test_non_subduction_points_are_unassigned_not_extrapolated():
    a = slab2.assign([0.0, 64.0, 40.0], [-140.0, -19.0, -100.0], [10.0, 10.0, 10.0])
    assert not np.any(a["assigned"])
    assert np.all(np.isnan(a["strike_deg"]))


def test_depth_misfit_is_returned_and_is_a_real_number():
    a = slab2.assign([54.6], [-160.5], [40.0])
    assert np.isfinite(a["depth_misfit_km"][0])
    assert a["depth_misfit_km"][0] == pytest.approx(
        abs(40.0 - a["interface_depth_km"][0]), abs=1e-9)


def test_depth_disambiguates_overlapping_slabs():
    """The tie-break is declared: closest interface IN DEPTH when a depth is given."""
    shallow = slab2.assign([54.6], [-160.5], [10.0])
    deep = slab2.assign([54.6], [-160.5], [200.0])
    assert shallow["assigned"][0] and deep["assigned"][0]
    # both resolve; the misfit differs because the event depth differs
    assert shallow["depth_misfit_km"][0] != deep["depth_misfit_km"][0]


def test_vectorised_and_scalar_agree():
    lat = [54.6, 38.3, 3.3]
    lon = [-160.5, 142.4, 95.9]
    dep = [40.0, 25.0, 30.0]
    vec = slab2.assign(lat, lon, dep)
    for i in range(3):
        one = slab2.assign([lat[i]], [lon[i]], [dep[i]])
        assert one["strike_deg"][0] == vec["strike_deg"][i]
        assert one["slab_code"][0] == vec["slab_code"][i]


def test_citation_and_scope_travel_with_the_answer():
    a = slab2.assign([54.6], [-160.5], [40.0])
    assert "Hayes" in a["citation"] and "10.5066" in a["citation"]
    assert "INTERFACE" in a["scope"] and "NOT modelled" in a["scope"]
