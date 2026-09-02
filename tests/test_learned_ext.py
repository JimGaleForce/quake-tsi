"""Failure-first tests for `engine_ext_forcing.py`.

These are COUNTED INVARIANTS against the real files on disk, not mocks. The numbers below
were MEASURED on 2026-09-02 and cross-checked against INVENTORY.md section 2. A parser
that silently drops rows, shifts a column index, or fills a gap with zero fails here
rather than in a result.

Run only this file:  python -u -m pytest tests/test_learned_ext.py -q
"""

from __future__ import annotations

import datetime as _dt
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

import engine_ext_forcing as F      # noqa: E402

DATA = os.path.join(HERE, "data", "spaceweather")

# ---- MEASURED 2026-09-02 against the files on disk -------------------------------
EOP_ROWS = 23381
EOP_JD_MIN, EOP_JD_MAX = 2437665.5, 2461045.5          # 1962-01-01 .. 2026-01-05
KP_ROWS = 276536
KP_JD_MIN, KP_JD_MAX = 2426707.5, 2461274.375          # 1932-01-01 .. 2026-08-21T21Z
OMNI_FILES = 32
OMNI_ROWS = 280512
OMNI_JD_MIN = 2449718.5                                # 1995-01-01T00Z


@pytest.fixture(scope="module")
def series():
    s, _a = F.load_forcing(DATA)
    return s


# ------------------------------------------------------------------ row counts --
def test_eop_row_count_and_span():
    jd, cols = F.parse_iers_eop(os.path.join(DATA, F.EOP_FILE))
    assert jd.size == EOP_ROWS, "IERS C04 row count changed: %d" % jd.size
    assert jd[0] == pytest.approx(EOP_JD_MIN)
    assert jd[-1] == pytest.approx(EOP_JD_MAX)
    assert set(cols) == {"x_pole", "y_pole", "lod"}
    # daily, strictly increasing, no duplicate or skipped day
    d = np.diff(jd)
    assert np.all(d == pytest.approx(1.0)), "IERS C04 is not a gapless daily grid"


def test_kp_row_count_span_and_cadence():
    jd, kp, ap = F.parse_kp_ap(os.path.join(DATA, F.KP_FILE))
    assert jd.size == KP_ROWS, "Kp row count changed: %d" % jd.size
    assert jd[0] == pytest.approx(KP_JD_MIN)
    assert jd[-1] == pytest.approx(KP_JD_MAX)
    assert np.allclose(np.diff(jd), 0.125), "Kp is not a gapless 3-hourly grid"
    # header lines must NOT have leaked in as data
    assert kp.size == jd.size and ap.size == jd.size


def test_omni_file_count_row_count_and_grid():
    paths = [os.path.join(DATA, f) for f in sorted(os.listdir(DATA))
             if f.startswith("omni2_") and f.endswith(".dat")]
    assert len(paths) == OMNI_FILES, "expected %d OMNI2 year files" % OMNI_FILES
    jd, cols, stat = F.parse_omni2(paths)
    assert stat["n_rows"] == OMNI_ROWS, "OMNI2 row count changed: %d" % stat["n_rows"]
    assert stat["n_short_rows"] == 0, "a short OMNI2 row would shift every column index"
    assert jd[0] == pytest.approx(OMNI_JD_MIN)
    assert np.allclose(np.diff(jd), 1.0 / 24.0), "OMNI2 is not a gapless hourly grid"
    assert set(cols) == set(F.OMNI_COLS)


# ------------------------------------------------------------- column identity --
def test_omni_column_indices_are_the_declared_quantities(series):
    """A shifted index would give a plausible-looking series of the wrong quantity.

    These are physical range checks on the finite values, plus the OMNI2 first-row
    literal for 2015-01-01T00Z which was read off the file by hand.
    """
    o = series["omni"]
    fin = lambda k: o[k][np.isfinite(o[k])]      # noqa: E731
    assert 0.0 < fin("imf_b").mean() < 20.0
    assert abs(fin("bz_gsm").mean()) < 2.0        # near zero on average, by symmetry
    assert 200.0 < fin("sw_speed").mean() < 700.0
    assert 0.5 < fin("proton_density").mean() < 30.0
    assert -60.0 < fin("dst").mean() < 10.0       # ring current is negative on average

    i = int(np.searchsorted(o["jd"], 2457023.5))  # 2015-01-01T00:00Z
    assert o["jd"][i] == pytest.approx(2457023.5)
    assert o["imf_b"][i] == pytest.approx(6.1)
    assert o["bz_gsm"][i] == pytest.approx(1.1)
    assert o["proton_density"][i] == pytest.approx(2.8)
    assert o["sw_speed"][i] == pytest.approx(568.0)
    assert o["dst"][i] == pytest.approx(-18.0)


def test_eop_first_row_literal(series):
    """First data line of eopc04_14_IAU2000.62-now.csv, read off the file by hand."""
    e = series["eop"]
    assert e["jd"][0] == pytest.approx(2437665.5)          # MJD 37665
    assert e["x_pole"][0] == pytest.approx(-0.012700)
    assert e["y_pole"][0] == pytest.approx(0.213000)
    assert e["lod"][0] == pytest.approx(0.0017230)


def test_kp_first_row_literal(series):
    """First data line of Kp_ap_since_1932.txt: 1932 01 01 00.0 ... 3.333 18 1."""
    k = series["kp"]
    assert k["jd"][0] == pytest.approx(2426707.5)
    assert k["kp"][0] == pytest.approx(3.333)
    assert k["ap"][0] == pytest.approx(18.0)


# ------------------------------------------------------------- interpolation --
def test_hold_returns_the_interval_value_not_the_neighbour(series):
    """Kp is a step function over its own 3-hour interval; a hold must not average."""
    k = series["kp"]
    i = 5000
    jd0 = k["jd"][i]
    for frac in (0.0, 0.25, 0.5, 0.999):
        q = jd0 + frac * 0.125
        got = F.sample_hold(k["jd"], k["kp"], np.array([q]))[0]
        assert got == pytest.approx(k["kp"][i]), "hold leaked across the interval edge"
    nxt = F.sample_hold(k["jd"], k["kp"], np.array([jd0 + 0.125]))[0]
    assert nxt == pytest.approx(k["kp"][i + 1])


def test_linear_hits_samples_exactly_and_midpoints_correctly(series):
    e = series["eop"]
    i = 10000
    exact = F.sample_linear(e["jd"], e["lod"], np.array([e["jd"][i]]))[0]
    assert exact == pytest.approx(e["lod"][i], rel=0, abs=1e-15)
    mid = F.sample_linear(e["jd"], e["lod"], np.array([e["jd"][i] + 0.5]))[0]
    assert mid == pytest.approx(0.5 * (e["lod"][i] + e["lod"][i + 1]))


def test_off_coverage_is_nan_not_zero(series):
    """The whole fill policy in one assertion: off-coverage is NaN, never 0.0."""
    for src, val in (("eop", "lod"), ("kp", "kp"), ("omni", "sw_speed")):
        s = series[src]
        before = np.array([s["jd"][0] - 10.0])
        after = np.array([s["jd"][-1] + 10.0])
        f = F.sample_linear if src == "eop" else F.sample_hold
        assert np.isnan(f(s["jd"], s[val], before)[0])
        assert np.isnan(f(s["jd"], s[val], after)[0])


def test_trailing_mean_is_nan_aware_and_trailing():
    v = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
    m = F.trailing_mean(v, 3)
    assert m[0] == pytest.approx(1.0)
    assert m[1] == pytest.approx(1.5)
    assert m[2] == pytest.approx(1.5)             # NaN ignored, not treated as 0
    assert m[3] == pytest.approx(3.0)             # mean(2, nan, 4) = 3
    assert m[4] == pytest.approx(4.5)             # mean(nan, 4, 5) = 4.5
    assert np.all(np.isnan(F.trailing_mean(np.array([np.nan, np.nan]), 2)))


# ----------------------------------------------------------------- the block --
def test_jd_from_days_since_is_explicit_about_epoch():
    """The 9131-day trap: days-since-1995 read as days-since-1970 (see the EPOCH AUDIT)."""
    e95 = _dt.datetime(1995, 1, 1, tzinfo=_dt.timezone.utc)
    e70 = _dt.datetime(1970, 1, 1, tzinfo=_dt.timezone.utc)
    assert F.jd_from_days_since(0.0, e70) == pytest.approx(F.UNIX_EPOCH_JD)
    assert F.jd_from_days_since(0.0, e95) == pytest.approx(2449718.5)
    assert (F.jd_from_days_since(0.0, e95)
            - F.jd_from_days_since(0.0, e70)) == pytest.approx(9131.0)


def test_event_block_shapes_flags_and_no_silent_fill(series):
    jd = F.jd_from_days_since(np.array([5000.0, 6000.0, 7000.0, 8000.0]),
                              _dt.datetime(2000, 1, 1, tzinfo=_dt.timezone.utc))
    cols, fill = F.event_block(jd, series)
    assert set(cols) == set(F.BLOCK3)
    for c in F.BLOCK3:
        assert cols[c].shape == jd.shape
        assert not np.isinf(cols[c]).any()
    assert fill["n_rows"] == 4
    for c in F.BLOCK3_VALUES:
        assert fill["per_feature_nan"][c] == int((~np.isfinite(cols[c])).sum())
    # flags are exactly the NaN mask of their source's anchor column
    assert np.array_equal(cols["miss_eop"], (~np.isfinite(cols["lod"])).astype(float))
    assert np.array_equal(cols["miss_kp"], (~np.isfinite(cols["kp"])).astype(float))
    assert np.array_equal(cols["miss_omni"],
                          (~np.isfinite(cols["sw_speed"])).astype(float))


def test_event_block_off_coverage_flags_rather_than_zero_fills(series):
    """1980 is inside EOP and Kp coverage but a decade before OMNI2 begins."""
    jd = F.jd_from_days_since(np.array([0.0]),
                              _dt.datetime(1980, 6, 1, tzinfo=_dt.timezone.utc))
    cols, fill = F.event_block(jd, series)
    assert cols["miss_omni"][0] == 1.0
    assert cols["miss_eop"][0] == 0.0
    assert cols["miss_kp"][0] == 0.0
    for c in F.OMNI_INSTANT + F.OMNI_MEAN24:
        assert np.isnan(cols[c][0]), "%s was filled instead of flagged" % c
        assert cols[c][0] != 0.0
    assert np.isfinite(cols["lod"][0]) and np.isfinite(cols["kp"][0])


def test_lod_differences_are_actual_differences(series):
    jd = F.jd_from_days_since(np.array([100.0, 200.0]),
                              _dt.datetime(2010, 1, 1, tzinfo=_dt.timezone.utc))
    cols, _f = F.event_block(jd, series)
    e = series["eop"]
    direct1 = (F.sample_linear(e["jd"], e["lod"], jd)
               - F.sample_linear(e["jd"], e["lod"], jd - 1.0))
    direct7 = (F.sample_linear(e["jd"], e["lod"], jd)
               - F.sample_linear(e["jd"], e["lod"], jd - 7.0))
    assert np.allclose(cols["lod_d1"], direct1)
    assert np.allclose(cols["lod_d7"], direct7)
    assert np.allclose(cols["pm_radius"], np.hypot(cols["x_pole"], cols["y_pole"]))


def test_schedule_audit_measures_the_named_artifacts(series):
    """ARTIFACT 1 must be MEASURED, not asserted away -- the audit has to return it."""
    lo = F.jd_from_days_since(0.0, _dt.datetime(2008, 1, 1, tzinfo=_dt.timezone.utc))
    hi = F.jd_from_days_since(0.0, _dt.datetime(2015, 1, 1, tzinfo=_dt.timezone.utc))
    a = F.schedule_audit(series, float(lo), float(hi))
    for nm in ("kp", "ap", "sw_speed", "dst"):
        assert nm in a
        assert a[nm]["n"] > 1000
        assert a[nm]["ut_diurnal_amplitude_over_sd"] is not None
        assert a[nm]["weekday_peak_to_trough_over_sd"] is not None
    # ARTIFACT 2: Kp really is quantised, so the audit's distinct-value count is small
    assert a["kp"]["n_distinct_values"] < 40
    assert a["sw_speed"]["n_distinct_values"] > 500
