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


# ------------------------------- §P7-8(c)(4): REFUSE, never silently substitute --
# The bug: engine/datasets.py:CATALOG_MAG_FLOOR = 4.5 meant a run declaring
# --mag-target 4.0 selected exactly the same events as one declaring 4.5, and
# session_20260812T021707 came out BITWISE IDENTICAL to session_20260812T004857
# while declaring 550 further tests. The ruling is that an unsupported --mag must
# RAISE with a clear message, never silently substitute.
@pytest.mark.parametrize("mag", [4.4999, 4.0, 3.0, 0.0, -1.0])
def test_sub_floor_magnitude_is_refused(mag):
    with pytest.raises(datasets.UnsupportedMagnitude):
        datasets.assert_mag_supported(mag)


@pytest.mark.parametrize("mag", [4.5, 4.6, 5.0, 6.5])
def test_supported_magnitude_is_returned_unchanged(mag):
    assert datasets.assert_mag_supported(mag) == pytest.approx(mag)


def test_refusal_message_names_the_constant_the_value_and_the_ruling():
    """A refusal a reader cannot act on is only half a fix."""
    with pytest.raises(datasets.UnsupportedMagnitude) as e:
        datasets.assert_mag_supported(4.0, what="--mag-target")
    msg = str(e.value)
    assert "--mag-target=4" in msg
    assert "CATALOG_MAG_FLOOR" in msg
    assert "4.5" in msg
    assert "SILENTLY CLAMPED" in msg
    assert "P7-8(c)" in msg


def test_loader_itself_refuses_a_sub_floor_request():
    """Not just the CLI: going straight at the loader must not reopen the clamp."""
    with pytest.raises(datasets.UnsupportedMagnitude):
        datasets.load_catalog(DATA, mag_floor=4.0)


def test_unsupported_magnitude_is_a_value_error():
    """Callers that catch ValueError keep working; nothing silently passes."""
    assert issubclass(datasets.UnsupportedMagnitude, ValueError)


def test_cli_refuses_sub_floor_mag_target_before_doing_any_work():
    """The refusal fires at argument-parse time -- no design, no session dir."""
    from engine import cli
    for argv in (["mine", "--mag-target", "4.0"],
                 ["fit-etas", "--mag-target", "4.0"],
                 ["run", "--covariate", "recent_rate", "--mag-target", "4.0"]):
        with pytest.raises(datasets.UnsupportedMagnitude):
            cli.main(argv)
