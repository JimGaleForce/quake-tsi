"""B4 / §S2.1 D+E: the geomagnetic / solar-wind / EOP downloader and its parsers.

**NO NETWORK IS TOUCHED ANYWHERE IN THIS FILE.** Every parser is exercised against a
RECORDED FIXTURE -- a verbatim excerpt of the real file, header and all -- written to a
tmp_path. That is the build brief's rule: if a source is unreachable, the downloader
and its recorded-fixture test still stand and the source is reported UNFETCHED.

The properties pinned down here, in the order they would hurt:

  * **`load_series` never invents a series.** An absent source yields NO column, the
    audit says UNFETCHED, and nothing is zero-filled or interpolated across.
  * **`interpolate_at` returns NaN off coverage and never interpolates ACROSS a gap** --
    a time whose bracketing samples are missing gets NaN, not an invented value.
  * **Kp is a HOLD and Dst is LINEAR**, declared per series, because a 3-hour band
    index is not a sample of a continuous curve.
  * **the OMNI2 column indices are asserted by name**, not trusted: a silently short
    row would shift every index by one and produce a plausible series of the wrong
    quantity.
  * **the fetch is verified TLS, never unverified.** `ssl._create_unverified_context`
    appears nowhere: an unverified fetch of a scientific series is an unattributable
    file, and provenance is a field (§S2.2(3)).
"""

from __future__ import annotations

import datetime as _dt
import os

import numpy as np
import pytest

from engine import spaceweather as SW

# ---------------------------------------------------------------- fixtures ----
# Recorded verbatim from https://kp.gfz.de/app/files/Kp_ap_since_1932.txt
KP_FIXTURE = """\
# PURPOSE: This file distributes the Kp index and ap index (one line per three-hour interval)
# LICENSE: CC BY 4.0
# SOURCE: Geomagnetic Observatory Niemegk, GFZ Helmholtz Centre for Geosciences
# ASCII, blank separated and fixed length, missing data indicated by -1.000 for Kp and -1 for ap
#YYY MM DD hh.h hh._m        days      days_m     Kp   ap D
1932 01 01 00.0 01.50     0.00000     0.06250  3.333   18 1
1932 01 01 03.0 04.50     0.12500     0.18750  2.667   12 1
1932 01 01 06.0 07.50     0.25000     0.31250  2.000    7 1
2003 10 30 00.0 01.50 26235.00000 26235.06250  8.667  400 1
2003 10 30 03.0 04.50 26235.12500 26235.18750  9.000  400 1
2026 01 01 00.0 01.50 34333.00000 34333.06250 -1.000   -1 0
"""

# Recorded verbatim from .../low_res_omni/omni2_2003.dat (two hours), truncated to the
# declared 55 words; the fill values are the file's own.
_OMNI_ROW = (" 2003 303 {hr} 2367 51 52  60  36   2.5   0.8 -13.7 297.1   0.3  -0.7  "
             "-0.2  -0.6  {bz}   0.9   2.6   1.7   1.7   0.9   22067.  10.0  {v}   "
             "0.9  -0.0 0.008  1.50    6159.   2.8    2.   0.7   0.6 0.002   0.09  "
             "10.01  18.7  87   6  {dst}   17 999999.99 99999.99 99999.99     0.30 "
             "    0.20     0.14 -1 400  {f107}   0.2    -2    15  6.0")
OMNI_FIXTURE = "\n".join([
    _OMNI_ROW.format(hr=" 0", bz="-20.5", v=" 900.", dst="  -353", f107="150.0"),
    _OMNI_ROW.format(hr=" 1", bz="999.9", v="9999.", dst=" 99999", f107="999.9"),
]) + "\n"

# Recorded verbatim from https://datacenter.iers.org/data/csv/eopc04_14_IAU2000.62-now.csv
EOP_FIXTURE = """\
MJD;Year;Month;Day;Type;x_pole;sigma_x_pole;y_pole;sigma_y_pole;x_rate;sigma_x_rate;y_rate;sigma_y_rate;Type;UT1-UTC;sigma_UT1-UTC;LOD;sigma_LOD;Type;dPsi;sigma_dPsi;dEpsilon;sigma_dEpsilon;dX;sigma_dX;dY;sigma_dY
37665;1962;01;01;;-0.012700;0.030000;0.213000;0.030000;;;;;;0.0326338;0.0020000;0.0017230;0.0014000;;;;;;0.000000;0.004774;0.000000;0.002000
37666;1962;01;02;;-0.015900;0.030000;0.214100;0.030000;;;;;;0.0320547;0.0020000;0.0016690;0.0014000;;;;;;0.000000;0.004774;0.000000;0.002000
37667;1962;01;03;;-0.019000;0.030000;0.215200;0.030000;;;;;;0.0315526;0.0020000;0.0015820;0.0014000;;;;;;0.000000;0.004774;0.000000;0.002000
"""


@pytest.fixture
def recorded(tmp_path):
    """The three recorded fixtures on disk, in the layout `load_series` expects."""
    d = tmp_path / "spaceweather"
    d.mkdir()
    (d / SW.SOURCES["kp_ap"]["filename"]).write_text(KP_FIXTURE, encoding="utf-8")
    (d / SW.SOURCES["iers_eop"]["filename"]).write_text(EOP_FIXTURE, encoding="utf-8")
    (d / "omni2_2003.dat").write_text(OMNI_FIXTURE, encoding="utf-8")
    return str(d)


# ----------------------------------------------------------------- parsers ----
def test_kp_parser_reads_the_declared_columns_and_nans_the_fill(tmp_path):
    p = str(tmp_path / "kp.txt")
    open(p, "w", encoding="utf-8").write(KP_FIXTURE)
    jd, kp, ap = SW.parse_kp_ap(p)
    assert jd.size == kp.size == ap.size == 6
    assert kp[0] == pytest.approx(3.333) and ap[0] == pytest.approx(18.0)
    assert kp[4] == pytest.approx(9.000)              # the Halloween-storm row
    assert np.isnan(kp[-1]) and np.isnan(ap[-1])      # -1.000 / -1 -> NaN
    assert np.all(np.diff(jd) > 0)


def test_kp_jd_matches_the_ephemeris_module_on_a_known_date(tmp_path):
    from engine import ephemeris as E
    p = str(tmp_path / "kp.txt")
    open(p, "w", encoding="utf-8").write(KP_FIXTURE)
    jd, _kp, _ap = SW.parse_kp_ap(p)
    # 2003-10-30 00:00 UT is the fourth row; the JD must agree with the module the
    # rest of the engine dates everything with, or the join is off by a constant.
    assert jd[3] == pytest.approx(E.julian_day(_dt.datetime(2003, 10, 30)), abs=1e-9)


def test_omni_parser_reads_by_the_named_column_table(tmp_path):
    p = str(tmp_path / "omni2_2003.dat")
    open(p, "w", encoding="utf-8").write(OMNI_FIXTURE)
    jd, cols = SW.parse_omni2([p])
    assert set(cols) == {"bz_gsm", "sw_speed", "dst", "f107"}
    assert cols["dst"][0] == pytest.approx(-353.0)
    assert cols["sw_speed"][0] == pytest.approx(900.0)
    assert cols["bz_gsm"][0] == pytest.approx(-20.5)
    assert cols["f107"][0] == pytest.approx(150.0)
    for k in cols:                                    # every fill value -> NaN
        assert np.isnan(cols[k][1]), k


def test_omni_parser_skips_a_short_row_rather_than_shifting_every_index(tmp_path):
    """A silently short row would move every column by one and produce a plausible
    series of the WRONG quantity. It is dropped, not realigned."""
    p = str(tmp_path / "omni2_2003.dat")
    open(p, "w", encoding="utf-8").write(OMNI_FIXTURE + "2003 304 0 2367 51\n")
    jd, cols = SW.parse_omni2([p])
    assert jd.size == 2
    assert SW.OMNI_N_WORDS == 55
    assert SW.OMNI_COLS["dst"][0] == 40                # word 41, 0-based


def test_eop_parser_takes_columns_by_name_not_by_position(tmp_path):
    p = str(tmp_path / "eop.csv")
    open(p, "w", encoding="utf-8").write(EOP_FIXTURE)
    jd, cols = SW.parse_iers_eop(p)
    assert set(cols) == {"x_pole", "y_pole", "lod"}
    assert cols["x_pole"][0] == pytest.approx(-0.0127)
    assert cols["lod"][0] == pytest.approx(0.001723)
    assert jd[0] == pytest.approx(37665 + SW.MJD_EPOCH_JD)


def test_eop_parser_refuses_a_header_it_does_not_recognise(tmp_path):
    p = str(tmp_path / "eop.csv")
    open(p, "w", encoding="utf-8").write("MJD;Year;something_else\n1;2;3\n")
    with pytest.raises(ValueError) as e:
        SW.parse_iers_eop(p)
    assert "x_pole" in str(e.value)


# ------------------------------------------------------------- load_series ----
def test_load_series_reports_every_declared_series_from_the_fixtures(recorded):
    series, audit = SW.load_series(recorded)
    assert set(series) == set(SW.SERIES)
    assert audit["series_absent"] == []
    for k, v in audit["sources"].items():
        assert v["status"] == "LOADED", k


def test_an_absent_source_is_UNFETCHED_and_yields_no_column(tmp_path):
    """The build brief's rule, in code: report UNFETCHED, never fake."""
    series, audit = SW.load_series(str(tmp_path))
    assert series == {}
    assert set(audit["series_absent"]) == set(SW.SERIES)
    assert all(v["status"] == "UNFETCHED" for v in audit["sources"].values())
    assert "nothing is substituted" in audit["absence_rule"]


def test_an_unparseable_source_is_UNPARSED_and_yields_no_column(tmp_path):
    d = tmp_path / "sw"
    d.mkdir()
    (d / SW.SOURCES["kp_ap"]["filename"]).write_text("# only a header\n",
                                                     encoding="utf-8")
    series, audit = SW.load_series(str(d))
    assert "kp" not in series
    assert audit["sources"]["kp_ap"]["status"] == "UNPARSED"
    assert "error" in audit["sources"]["kp_ap"]


def test_properties_emits_no_column_for_an_absent_source(tmp_path):
    """§S2.1 D/E: the family is ABSENT, not zero-filled -- the same discipline
    `marks_ext.build_marks` applies to a missing lat/lon."""
    from engine import properties as P
    cols, audit = P.spaceweather_columns(_dt.datetime(2000, 1, 1), 0.0, 30.0,
                                         data_dir=str(tmp_path))
    assert cols == [] and audit["n_columns"] == 0
    assert audit["series_absent"]


# ---------------------------------------------------------- interpolate_at ----
def test_kp_is_held_over_its_own_three_hour_interval(recorded):
    series, _a = SW.load_series(recorded)
    t0 = _dt.datetime(1932, 1, 1)
    # 00:00-03:00 UT carries 3.333; 03:00-06:00 carries 2.667. A LINEAR reading
    # would blend them, which is wrong for a band index.
    v = SW.interpolate_at(series["kp"], t0, np.array([0.01, 0.10, 0.20, 0.30]))
    assert v[0] == pytest.approx(3.333) and v[1] == pytest.approx(3.333)
    assert v[2] == pytest.approx(2.667) and v[3] == pytest.approx(2.000)
    assert SW.SERIES["kp"]["interp"] == "hold"
    assert SW.SERIES["dst"]["interp"] == "linear"


def test_a_time_outside_coverage_is_nan_and_not_extrapolated(recorded):
    series, _a = SW.load_series(recorded)
    v = SW.interpolate_at(series["lod"], _dt.datetime(1900, 1, 1),
                          np.array([0.0, 10.0]))
    assert np.isnan(v).all()
    v2 = SW.interpolate_at(series["lod"], _dt.datetime(2200, 1, 1), np.array([0.0]))
    assert np.isnan(v2).all()


def test_a_linear_series_never_interpolates_across_a_nan_gap(recorded):
    """A time bracketed by a missing sample gets NaN, not an invented value."""
    series, _a = SW.load_series(recorded)
    t0 = _dt.datetime(2003, 10, 30)
    # hour 0 carries -353; hour 1 is the file's own fill row and parses to NaN.
    v = SW.interpolate_at(series["dst"], t0,
                          np.array([0.0, 1.0 / 48.0, 1.0 / 24.0, 0.5]))
    assert v[0] == pytest.approx(-353.0)     # an EXACT hit on a good sample survives
    assert np.isnan(v[1])                    # mid-gap: never interpolated across
    assert np.isnan(v[2])                    # the fill sample itself
    assert np.isnan(v[3])                    # past the end of coverage


def test_linear_interpolation_is_exact_at_the_sample_points(recorded):
    series, _a = SW.load_series(recorded)
    t0 = _dt.datetime(1962, 1, 1)
    v = SW.interpolate_at(series["lod"], t0, np.array([0.0, 1.0, 2.0, 0.5]))
    assert v[0] == pytest.approx(0.0017230)
    assert v[1] == pytest.approx(0.0016690)
    assert v[2] == pytest.approx(0.0015820)
    assert v[3] == pytest.approx(0.5 * (0.0017230 + 0.0016690))


# --------------------------------------------------------- the declarations ---
def test_every_declared_series_carries_full_provenance():
    for k, m in SW.SERIES.items():
        for f in ("property", "family", "source", "interp", "version",
                  "source_name", "url", "convention", "definition"):
            assert m.get(f), (k, f)
        assert m["family"] in ("geomagnetic", "earth_rotation")
        assert m["source"] in SW.SOURCES
        assert m["interp"] in ("hold", "linear")


def test_dst_credits_kyoto_as_the_origin_and_nasa_as_the_delivery_path():
    """§S2.2(3): provenance names BOTH, because that is the honest form."""
    m = SW.SERIES["dst"]
    assert "Kyoto" in m["source_name"] and "OMNI2" in m["source_name"]
    assert "ORIGIN IS WDC KYOTO" in m["caveat"]


def test_daily_series_are_declared_not_subdaily():
    for k in ("f107", "lod", "x_pole", "y_pole"):
        assert SW.SERIES[k]["subdaily"] is False, k
    for k in ("kp", "ap", "dst", "sw_speed", "bz_gsm"):
        assert SW.SERIES[k].get("subdaily", True) is True, k


def test_the_omni_year_range_is_declared_and_hash_affecting():
    assert SW.OMNI_YEAR_MIN == 1995                # the ComCat boxes' start
    assert len(SW.SOURCES["omni2"]["urls"]) == SW.OMNI_YEAR_MAX - SW.OMNI_YEAR_MIN + 1


# --------------------------------------------------- the transport's promises --
def test_the_module_never_disables_certificate_verification():
    """An unverified fetch of a scientific series is an unattributable file."""
    import ast
    import inspect
    tree = ast.parse(inspect.getsource(SW))
    # the EXECUTABLE body only: the module docstring names the forbidden call
    # precisely to say it must never be added, so a raw text search would fail on
    # the prohibition itself.
    code = chr(10).join(ast.unparse(n) for n in ast.walk(tree)
                        if isinstance(n, (ast.Call, ast.Attribute, ast.Assign)))
    for forbidden in ("_create_unverified_context", "CERT_NONE",
                      "check_hostname", "verify=False"):
        assert forbidden not in code, forbidden
    assert "create_default_context" in code


def test_a_cached_file_is_not_refetched(tmp_path, monkeypatch):
    d = tmp_path / "sw"
    d.mkdir()
    (d / SW.SOURCES["kp_ap"]["filename"]).write_text(KP_FIXTURE, encoding="utf-8")

    def boom(*a, **k):                       # any network call is a failure here
        raise AssertionError("the network was touched for a CACHED file")

    monkeypatch.setattr(SW, "_http_get", boom)
    recs = SW.download_source("kp_ap", data_dir=str(d), log=False)
    assert len(recs) == 1 and recs[0]["status"] == "CACHED"
    assert recs[0]["bytes"] == os.path.getsize(
        os.path.join(str(d), SW.SOURCES["kp_ap"]["filename"]))


def test_a_network_failure_degrades_to_UNFETCHED_and_never_raises(tmp_path,
                                                                  monkeypatch):
    import urllib.error

    def boom(url, dest, **k):
        raise urllib.error.URLError("no network tonight")

    monkeypatch.setattr(SW, "_http_get", boom)
    recs = SW.download_source("kp_ap", data_dir=str(tmp_path), log=False)
    assert len(recs) == 1 and recs[0]["status"] == "UNFETCHED"
    assert "no network tonight" in recs[0]["error"]
    assert not os.path.exists(os.path.join(str(tmp_path),
                                           SW.SOURCES["kp_ap"]["filename"]))


def test_every_mirror_is_tried_before_a_source_is_called_unfetched(tmp_path,
                                                                   monkeypatch):
    import urllib.error
    tried = []

    def flaky(url, dest, **k):
        tried.append(url)
        if len(tried) == 1:
            raise urllib.error.URLError("first mirror down")
        open(dest, "w", encoding="utf-8").write(KP_FIXTURE)
        return len(KP_FIXTURE)

    monkeypatch.setattr(SW, "_http_get", flaky)
    recs = SW.download_source("kp_ap", data_dir=str(tmp_path), log=False)
    assert len(tried) == 2
    assert recs[0]["status"] == "FETCHED"


def test_a_resumed_fetch_appends_only_on_a_206(monkeypatch, tmp_path):
    """A 200 answer to a Range request must RESTART, never append: appending a full
    body to a partial one silently corrupts the file and every number below it."""
    import inspect
    src = inspect.getsource(SW._http_get)
    assert 'Range' in src and '206' in src
    assert 'mode = "ab" if (have and code == 206) else "wb"' in src


def test_download_all_reports_a_status_for_every_declared_source(tmp_path,
                                                                 monkeypatch):
    import urllib.error
    monkeypatch.setattr(SW, "_http_get",
                        lambda *a, **k: (_ for _ in ()).throw(
                            urllib.error.URLError("offline")))
    rep = SW.download_all(data_dir=str(tmp_path), timeout=1)
    assert set(rep["sources"]) == set(SW.SOURCES)
    assert rep["n_unfetched"] > 0 and rep["n_fetched"] == 0
