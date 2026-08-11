import os

import pytest

from engine import datasets

DATA = datasets.DEFAULT_DATA_DIR
have_data = os.path.isdir(DATA) and bool(os.listdir(DATA))
pytestmark = pytest.mark.skipif(not have_data, reason="comcat_world data not present")


def test_load_catalog_invariants():
    cat, rep = datasets.load_catalog(DATA)
    assert 50_000 < rep["n_events"] < 150_000
    assert rep["dup_removed"] > 0, "overlapping boxes must share events"
    assert cat["mag"].min() >= datasets.CATALOG_MAG_FLOOR
    assert cat["t"].is_monotonic_increasing
    assert cat["id"].is_unique
    assert cat["t"].dt.tz is not None, "timestamps must be tz-aware UTC"
    assert len(rep["files"]) >= 10


def test_no_subfloor_rows_survive():
    cat, rep = datasets.load_catalog(DATA, mag_floor=5.0)
    assert cat["mag"].min() >= 5.0
    assert rep["dropped_floor"] > 0
