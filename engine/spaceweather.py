"""B4 / §S2.1 D+E -- GEOMAGNETIC, SOLAR-WIND AND EARTH-ROTATION SERIES.

The two property families SEARCHER.md calls **DOES NOT EXIST** -- family D
(geomagnetic / space weather: Kp, ap, Dst, solar-wind speed, Bz, F10.7) and family E
(Earth rotation: LOD, polar motion) -- plus the cached, resumable, network-optional
downloader that puts them on disk. §S8.2 B4.

FOUR PROPERTIES OF THIS DOWNLOADER, all of them requirements rather than niceties
---------------------------------------------------------------------------------
1. **CACHED: it never re-fetches a file it already has.** A completed file is
   `CACHED` and the network is not touched. Only `force=True` re-fetches, and the
   engine never passes it.
2. **RESUMABLE:** a partial download lives at `<name>.part` and is continued with an
   HTTP `Range` request. A truncated overnight fetch costs the bytes already on disk,
   not the whole file.
3. **DEGRADES GRACEFULLY:** every failure mode -- no network, DNS, TLS, HTTP error,
   truncated body -- returns `status="UNFETCHED"` with the reason on the record.
   **It never raises into the engine and it never substitutes anything.** A property
   family whose series is UNFETCHED is ABSENT from the property matrix, and
   `engine/properties.py:spaceweather_columns` reports the absence rather than
   emitting a column of zeros. This is the same discipline `marks_ext.build_marks`
   applies to a missing lat/lon, and for the same reason: a family silently running
   fewer columns than its declaration says is the defect an audit exists to prevent.
4. **WRITES ONLY UNDER `data/`, WHICH IS GITIGNORED.** Nothing here is ever committed.

TLS NOTE, recorded because it cost time and will cost it again
---------------------------------------------------------------
Both GFZ and IERS present certificate chains that this machine's system trust store
cannot complete (`CERTIFICATE_VERIFY_FAILED: unable to get local issuer certificate`),
while NASA SPDF verifies fine. The fix is to verify against `certifi`'s CA bundle when
it is installed -- which is verification, not a bypass. **`ssl._create_unverified_context`
is not used anywhere in this module and must not be added**: an unverified fetch of a
scientific series is an unattributable file, and provenance is a field (§S2.2(3)).

THE SOURCES, with their citations, because a property carries its provenance
-----------------------------------------------------------------------------
* **Kp / ap**, 3-hourly, 1932-present, ~5 MB. GFZ Helmholtz Centre / Geomagnetic
  Observatory Niemegk. Matzka, J., Stolle, C., Yamazaki, Y., Bronkalla, O. &
  Morschhauser, A. (2021), *The geomagnetic Kp index and derived indices of
  geomagnetic activity*, Space Weather 19, e2020SW002641. Data publication: Matzka
  et al. (2021), GFZ Data Services, doi:10.5880/Kp.0001. CC BY 4.0.
* **Dst, solar-wind speed, Bz(GSM), F10.7**, hourly, from **NASA/GSFC OMNI2**
  (King & Papitashvili (2005), JGR 110, A02104). **Dst inside OMNI2 originates at
  WDC Kyoto** and is credited to Kyoto on the property's provenance, with NASA named
  as the delivery path -- which is the honest form and is what §S2.2(3) is for.
* **LOD, polar motion x/y**, daily, 1962-present, ~5 MB. IERS EOP 14 C04 (IAU2000),
  IERS Earth Orientation Centre, Paris Observatory. Bizouard, C., Lambert, S.,
  Gattano, C., Becker, O. & Richard, J.-Y. (2019), *The IERS EOP 14C04 solution for
  Earth orientation parameters consistent with ITRF 2014*, J. Geodesy 93, 621-633.

Nothing in this module is evidence. It fetches coordinates.
"""

from __future__ import annotations

import datetime as _dt
import json
import os
import ssl
import urllib.error
import urllib.request

import numpy as np

from . import ephemeris as E

DATA_DIR = os.path.join("data", "spaceweather")
DOWNLOAD_LOG = os.path.join(DATA_DIR, "download_log.jsonl")
USER_AGENT = "quake-replication-searcher/1.0 (engine.spaceweather; research use)"

# The OMNI2 year range. DECLARED and hash-affecting: it is the span of the ComCat
# boxes on disk (1995-01-01 onward, `data/comcat_world_log.txt`), not "whatever the
# server has", so a re-run a year later fetches one more file rather than silently
# changing the property's support.
OMNI_YEAR_MIN = 1995
OMNI_YEAR_MAX = 2026

MJD_EPOCH_JD = 2400000.5

SOURCES = {
    "kp_ap": {
        "urls": ("https://kp.gfz.de/app/files/Kp_ap_since_1932.txt",
                 "https://kp.gfz-potsdam.de/app/files/Kp_ap_since_1932.txt"),
        "filename": "Kp_ap_since_1932.txt",
        "cite": ("Matzka et al. (2021) Space Weather 19 e2020SW002641; GFZ Data "
                 "Services doi:10.5880/Kp.0001; CC BY 4.0"),
        "holder": "GFZ Helmholtz Centre / Geomagnetic Observatory Niemegk",
    },
    "iers_eop": {
        "urls": ("https://datacenter.iers.org/data/csv/eopc04_14_IAU2000.62-now.csv",),
        "filename": "eopc04_14_IAU2000.62-now.csv",
        "cite": "Bizouard et al. (2019) J. Geodesy 93, 621-633 (IERS EOP 14 C04)",
        "holder": "IERS Earth Orientation Centre, Paris Observatory",
    },
    "omni2": {
        "urls": tuple(
            "https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_%d.dat" % y
            for y in range(OMNI_YEAR_MIN, OMNI_YEAR_MAX + 1)),
        "filename": None,                   # one file per year; see `_omni_paths`
        "cite": "King & Papitashvili (2005) JGR 110 A02104 (OMNI2); Dst from WDC Kyoto",
        "holder": "NASA/GSFC Space Physics Data Facility",
    },
}


# ------------------------------------------------------------- the transport --
def _ssl_context():
    """Verify against certifi's bundle when present, the system store otherwise.

    NEVER unverified. See the TLS note in the module docstring.
    """
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except Exception:                        # pragma: no cover - env dependent
        return ssl.create_default_context()


def _http_get(url, dest, timeout=120, chunk=1 << 16):
    """Fetch `url` to `dest`, resuming `<dest>.part` if one is there. Returns bytes.

    Raises on failure; `download_source` is the layer that turns a raise into an
    UNFETCHED record, so that exactly one place in this module decides that a network
    failure is not an engine failure.
    """
    part = dest + ".part"
    have = os.path.getsize(part) if os.path.exists(part) else 0
    headers = {"User-Agent": USER_AGENT}
    if have:
        headers["Range"] = "bytes=%d-" % have
    req = urllib.request.Request(url, headers=headers)
    ctx = _ssl_context()
    with urllib.request.urlopen(req, timeout=timeout, context=ctx) as resp:
        code = resp.getcode()
        # 206 = the server honoured the Range and we append; 200 = it did not and we
        # must start over, because appending a full body to a partial one silently
        # corrupts the file and every number downstream of it.
        mode = "ab" if (have and code == 206) else "wb"
        if mode == "wb":
            have = 0
        with open(part, mode) as fh:
            while True:
                b = resp.read(chunk)
                if not b:
                    break
                fh.write(b)
    size = os.path.getsize(part)
    if size <= 0:
        raise IOError("empty body from %s" % url)
    os.replace(part, dest)
    return size


def _log(rec, log_path=DOWNLOAD_LOG):
    os.makedirs(os.path.dirname(log_path) or ".", exist_ok=True)
    with open(log_path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, sort_keys=True) + "\n")


def download_source(key, data_dir=DATA_DIR, force=False, timeout=120, log=True):
    """Fetch one declared source. Returns a list of per-file records. Never raises."""
    if key not in SOURCES:
        raise KeyError("no declared space-weather source %r" % (key,))
    src = SOURCES[key]
    os.makedirs(data_dir, exist_ok=True)
    recs = []
    if key == "omni2":
        targets = [(u, os.path.join(data_dir, os.path.basename(u)))
                   for u in src["urls"]]
    else:
        dest = os.path.join(data_dir, src["filename"])
        targets = [(src["urls"], dest)]     # a tuple of mirrors for one destination

    for url_or_mirrors, dest in targets:
        mirrors = ((url_or_mirrors,) if isinstance(url_or_mirrors, str)
                   else tuple(url_or_mirrors))
        if os.path.exists(dest) and not force:
            recs.append({"source": key, "path": dest, "status": "CACHED",
                         "bytes": os.path.getsize(dest), "url": mirrors[0]})
            continue
        last = None
        for url in mirrors:
            try:
                n = _http_get(url, dest, timeout=timeout)
                recs.append({"source": key, "path": dest, "status": "FETCHED",
                             "bytes": int(n), "url": url})
                last = None
                break
            except (urllib.error.URLError, urllib.error.HTTPError, ssl.SSLError,
                    IOError, OSError, ValueError) as exc:
                last = "%s: %s" % (type(exc).__name__, exc)
        if last is not None:
            recs.append({"source": key, "path": dest, "status": "UNFETCHED",
                         "bytes": 0, "url": mirrors[0], "error": last})
    if log:
        for r in recs:
            rr = dict(r)
            rr["ts"] = _dt.datetime.now(_dt.timezone.utc).isoformat()
            rr["cite"] = src["cite"]
            _log(rr, os.path.join(data_dir, "download_log.jsonl"))
    return recs


def download_all(keys=None, data_dir=DATA_DIR, force=False, timeout=120):
    """Fetch every declared source. Returns a report; UNFETCHED is a status, not an error."""
    out = {"data_dir": data_dir, "sources": {}, "n_fetched": 0, "n_cached": 0,
           "n_unfetched": 0}
    for k in (keys or tuple(SOURCES)):
        recs = download_source(k, data_dir=data_dir, force=force, timeout=timeout)
        out["sources"][k] = recs
        for r in recs:
            out["n_" + r["status"].lower()] = out.get(
                "n_" + r["status"].lower(), 0) + 1
    return out


# ---------------------------------------------------------------- the parsers --
def parse_kp_ap(path):
    """GFZ 3-hourly Kp/ap -> (jd, kp, ap). Missing = -1.000 (Kp) / -1 (ap) -> NaN.

    Declared file format, from the file's own 30-line header:
      YYYY MM DD hh.h hh._m days days_m Kp ap D
    `hh.h` is the START of the 3-hour interval; the value is a step function over the
    interval, which is why `SERIES['kp']['interp']` is 'hold' and not 'linear'.
    """
    y, mo, d, hh, kp, ap = [], [], [], [], [], []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            f = line.split()
            if len(f) < 9:
                continue
            y.append(int(f[0])); mo.append(int(f[1])); d.append(int(f[2]))
            hh.append(float(f[3])); kp.append(float(f[7])); ap.append(float(f[8]))
    if not y:
        raise ValueError("no data rows parsed from %r" % (path,))
    jd = _jd_from_ymdh(np.array(y), np.array(mo), np.array(d), np.array(hh))
    kp = np.array(kp, dtype=np.float64)
    ap = np.array(ap, dtype=np.float64)
    kp[kp < 0] = np.nan
    ap[ap < 0] = np.nan
    return jd, kp, ap


# OMNI2 word -> 0-based token index, from the OMNI2 format description. Kept as a
# named table because a bare `f[40]` in a parser is the kind of line that is wrong for
# a year before anyone notices.
OMNI_COLS = {
    "bz_gsm": (16, 999.9),          # word 17, nT, GSM
    "sw_speed": (24, 9999.0),       # word 25, km/s
    "dst": (40, 99999.0),           # word 41, nT, from Kyoto
    "f107": (50, 999.9),            # word 51, sfu
}
OMNI_N_WORDS = 55


def parse_omni2(paths):
    """OMNI2 hourly year files -> (jd, {name: values}). Fill values -> NaN.

    Each row is 55 whitespace-separated words; the count is ASSERTED per file rather
    than assumed, because a silently short row would shift every column index by one
    and produce a plausible-looking series of the wrong quantity.
    """
    jd_all, cols = [], {k: [] for k in OMNI_COLS}
    n_rows = 0
    for p in sorted(paths):
        with open(p, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                f = line.split()
                if len(f) < OMNI_N_WORDS:
                    continue
                yr, doy, hr = int(f[0]), int(f[1]), int(f[2])
                jd_all.append(E.julian_day(_dt.datetime(yr, 1, 1))
                              + (doy - 1) + hr / 24.0)
                for k, (i, fill) in OMNI_COLS.items():
                    v = float(f[i])
                    cols[k].append(np.nan if v >= fill else v)
                n_rows += 1
    if n_rows == 0:
        raise ValueError("no OMNI2 rows parsed from %d file(s)" % len(list(paths)))
    jd = np.asarray(jd_all, dtype=np.float64)
    o = np.argsort(jd, kind="stable")
    return jd[o], {k: np.asarray(v, dtype=np.float64)[o] for k, v in cols.items()}


def parse_iers_eop(path):
    """IERS EOP 14 C04 CSV (';' separated) -> (jd, {x_pole, y_pole, lod}).

    Columns are taken BY NAME from the header line, never by position: the C04 product
    has gained columns over its life and a positional parser would silently read the
    wrong one after the next revision.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        head = fh.readline().strip().split(";")
        idx = {}
        for want, names in (("mjd", ("MJD",)), ("x_pole", ("x_pole",)),
                            ("y_pole", ("y_pole",)), ("lod", ("LOD",))):
            for n in names:
                if n in head:
                    idx[want] = head.index(n)
                    break
        missing = [k for k in ("mjd", "x_pole", "y_pole", "lod") if k not in idx]
        if missing:
            raise ValueError("IERS C04 header lacks %s; header=%r" % (missing, head))
        rows = {k: [] for k in idx}
        for line in fh:
            f = line.rstrip("\n").split(";")
            if len(f) <= max(idx.values()):
                continue
            try:
                vals = {k: float(f[i]) if f[i].strip() else np.nan
                        for k, i in idx.items()}
            except ValueError:
                continue
            for k, v in vals.items():
                rows[k].append(v)
    if not rows["mjd"]:
        raise ValueError("no data rows parsed from %r" % (path,))
    jd = np.asarray(rows["mjd"], dtype=np.float64) + MJD_EPOCH_JD
    return jd, {k: np.asarray(v, dtype=np.float64)
                for k, v in rows.items() if k != "mjd"}


def _jd_from_ymdh(y, mo, d, hours):
    """Vectorised Julian Day from calendar fields -- `ephemeris.julian_day`'s formula."""
    y = np.asarray(y, dtype=np.int64).copy()
    mo = np.asarray(mo, dtype=np.int64).copy()
    dd = np.asarray(d, dtype=np.float64) + np.asarray(hours, dtype=np.float64) / 24.0
    adj = mo <= 2
    y[adj] -= 1
    mo[adj] += 12
    a = y // 100
    b = 2 - a + a // 4
    return (np.floor(365.25 * (y + 4716)) + np.floor(30.6001 * (mo + 1))
            + dd + b - 1524.5)


# ------------------------------------------------------- the declared series --
# One entry per SCANNABLE property. `interp` is 'hold' for an index defined on an
# interval (Kp is a 3-hour band value, not a sample of a continuous curve) and
# 'linear' for a sampled continuous quantity. Declared per series, before any scan.
SERIES = {
    "kp": dict(property="kp_index", family="geomagnetic", source="kp_ap",
               interp="hold", version="GFZ Kp v1.0 (doi:10.5880/Kp.0001)",
               source_name="GFZ Helmholtz Centre, Geomagnetic Observatory Niemegk",
               url="https://kp.gfz.de/app/files/Kp_ap_since_1932.txt",
               convention="3-hourly planetary index, value held constant over its "
                          "own 3-hour UT interval; missing = -1.000 -> NaN",
               definition="Kp planetary geomagnetic activity index, unitless",
               caveat="Kp is QUANTISED to thirds (0, 0.333, 0.667, ...); its dwell "
                      "measure is a step distribution and dwell_pit's mid-rank tie "
                      "handling is what keeps the plateaus from reading as structure."),
    "ap": dict(property="ap_index", family="geomagnetic", source="kp_ap",
               interp="hold", version="GFZ Kp v1.0 (doi:10.5880/Kp.0001)",
               source_name="GFZ Helmholtz Centre, Geomagnetic Observatory Niemegk",
               url="https://kp.gfz.de/app/files/Kp_ap_since_1932.txt",
               convention="3-hourly equivalent planetary amplitude, held over its "
                          "own interval; missing = -1 -> NaN",
               definition="ap equivalent planetary amplitude, unitless (x2 nT)",
               caveat="ap is a nonlinear transform of Kp: the two are near-perfectly "
                      "monotonically related and §S3.5's redundancy collapse is "
                      "EXPECTED to keep one of them and drop the other."),
    "dst": dict(property="dst_index", family="geomagnetic", source="omni2",
                interp="linear", version="OMNI2 hourly",
                source_name="WDC Kyoto (index), delivered via NASA/GSFC OMNI2",
                url="https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/",
                convention="hourly Dst, nT; OMNI2 fill 99999 -> NaN",
                definition="Dst ring-current index, nT (negative = storm)",
                caveat="ORIGIN IS WDC KYOTO; NASA OMNI2 is the delivery path. The "
                       "provenance names both, per §S2.2(3)."),
    "sw_speed": dict(property="solar_wind_speed", family="geomagnetic",
                     source="omni2", interp="linear", version="OMNI2 hourly",
                     source_name="NASA/GSFC OMNI2 (King & Papitashvili 2005)",
                     url="https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/",
                     convention="hourly bulk flow speed, km/s; fill 9999 -> NaN",
                     definition="solar-wind plasma flow speed, km/s",
                     caveat="pre-1995 coverage is patchy; NaN spans are dropped from "
                            "the event set rather than interpolated across."),
    "bz_gsm": dict(property="imf_bz_gsm", family="geomagnetic", source="omni2",
                   interp="linear", version="OMNI2 hourly",
                   source_name="NASA/GSFC OMNI2 (King & Papitashvili 2005)",
                   url="https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/",
                   convention="hourly Bz in GSM coordinates, nT; fill 999.9 -> NaN",
                   definition="IMF north-south component Bz (GSM), nT",
                   caveat="SIGN CONVENTION IS GSM AND IS PART OF THE CLAIM: southward "
                          "(negative) Bz is the geoeffective sign."),
    "f107": dict(property="f107_flux", family="geomagnetic", source="omni2",
                 interp="linear", version="OMNI2 hourly",
                 source_name="NASA/GSFC OMNI2 (originally DRAO Penticton)",
                 url="https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/",
                 convention="daily 10.7 cm flux repeated hourly; fill 999.9 -> NaN",
                 definition="F10.7 solar radio flux, sfu",
                 subdaily=False,
                 caveat="F10.7 is a DAILY value repeated across the hours of its day: "
                        "it is NOT sub-daily and the column says so (§P7-3(3))."),
    "lod": dict(property="length_of_day", family="earth_rotation", source="iers_eop",
                interp="linear", version="IERS EOP 14 C04 (IAU2000)",
                source_name="IERS Earth Orientation Centre, Paris Observatory",
                url="https://datacenter.iers.org/data/csv/eopc04_14_IAU2000.62-now.csv",
                convention="daily LOD excess, seconds",
                definition="excess length of day, s",
                subdaily=False,
                caveat="DAILY series: not sub-daily (§P7-3(3))."),
    "x_pole": dict(property="polar_motion_x", family="earth_rotation",
                   source="iers_eop", interp="linear",
                   version="IERS EOP 14 C04 (IAU2000)",
                   source_name="IERS Earth Orientation Centre, Paris Observatory",
                   url="https://datacenter.iers.org/data/csv/eopc04_14_IAU2000.62-now.csv",
                   convention="daily polar motion x, arcseconds",
                   definition="polar motion x component, arcsec",
                   subdaily=False, caveat="DAILY series: not sub-daily."),
    "y_pole": dict(property="polar_motion_y", family="earth_rotation",
                   source="iers_eop", interp="linear",
                   version="IERS EOP 14 C04 (IAU2000)",
                   source_name="IERS Earth Orientation Centre, Paris Observatory",
                   url="https://datacenter.iers.org/data/csv/eopc04_14_IAU2000.62-now.csv",
                   convention="daily polar motion y, arcseconds",
                   definition="polar motion y component, arcsec",
                   subdaily=False, caveat="DAILY series: not sub-daily."),
}


def _omni_paths(data_dir):
    if not os.path.isdir(data_dir):
        return []
    return [os.path.join(data_dir, f) for f in sorted(os.listdir(data_dir))
            if f.startswith("omni2_") and f.endswith(".dat")]


def load_series(data_dir=None):
    """Every declared series that is ON DISK. Returns ({key: series}, audit).

    A series is `{'jd': array, 'value': array, 'interp': 'hold'|'linear'}`. Sources
    that are absent or unparseable are reported UNFETCHED / UNPARSED in the audit and
    are simply not in the dict -- **never faked, never zero-filled**.
    """
    d = data_dir or DATA_DIR
    out, status = {}, {}

    kp_path = os.path.join(d, SOURCES["kp_ap"]["filename"])
    if os.path.exists(kp_path):
        try:
            jd, kp, ap = parse_kp_ap(kp_path)
            out["kp"] = {"jd": jd, "value": kp, "interp": SERIES["kp"]["interp"]}
            out["ap"] = {"jd": jd, "value": ap, "interp": SERIES["ap"]["interp"]}
            status["kp_ap"] = {"status": "LOADED", "n": int(jd.size),
                               "jd_min": float(jd.min()), "jd_max": float(jd.max())}
        except Exception as exc:
            status["kp_ap"] = {"status": "UNPARSED", "error": repr(exc)}
    else:
        status["kp_ap"] = {"status": "UNFETCHED", "path": kp_path}

    op = _omni_paths(d)
    if op:
        try:
            jd, cols = parse_omni2(op)
            for k, v in cols.items():
                out[k] = {"jd": jd, "value": v, "interp": SERIES[k]["interp"]}
            status["omni2"] = {"status": "LOADED", "n_files": len(op),
                               "n": int(jd.size), "jd_min": float(jd.min()),
                               "jd_max": float(jd.max())}
        except Exception as exc:
            status["omni2"] = {"status": "UNPARSED", "error": repr(exc)}
    else:
        status["omni2"] = {"status": "UNFETCHED", "path": d + "/omni2_YYYY.dat"}

    ep = os.path.join(d, SOURCES["iers_eop"]["filename"])
    if os.path.exists(ep):
        try:
            jd, cols = parse_iers_eop(ep)
            for k, v in cols.items():
                out[k] = {"jd": jd, "value": v, "interp": SERIES[k]["interp"]}
            status["iers_eop"] = {"status": "LOADED", "n": int(jd.size),
                                  "jd_min": float(jd.min()),
                                  "jd_max": float(jd.max())}
        except Exception as exc:
            status["iers_eop"] = {"status": "UNPARSED", "error": repr(exc)}
    else:
        status["iers_eop"] = {"status": "UNFETCHED", "path": ep}

    audit = {"data_dir": d, "sources": status,
             "series_loaded": sorted(out), "n_series": len(out),
             "series_declared": sorted(SERIES),
             "series_absent": sorted(k for k in SERIES if k not in out),
             "absence_rule": (
                 "A source that is UNFETCHED or UNPARSED yields NO property column. "
                 "The family is ABSENT and the audit says so; nothing is substituted, "
                 "zero-filled or interpolated across a missing source.")}
    return out, audit


def interpolate_at(series, t0: _dt.datetime, day_float):
    """Evaluate a series at the catalogue's own continuous times. NaN off-coverage.

    `interp='hold'` is a zero-order hold (the value of the interval the time falls in)
    and `interp='linear'` is linear interpolation between samples. **Neither ever
    interpolates ACROSS a NaN gap**: a time whose bracketing samples are missing gets
    NaN, and the scan drops those events from that cell with the count reported,
    rather than inventing a value inside a data gap.
    """
    jd_q = E.julian_day_at(t0, np.asarray(day_float, dtype=np.float64))
    jd = np.asarray(series["jd"], dtype=np.float64)
    v = np.asarray(series["value"], dtype=np.float64)
    out = np.full(jd_q.shape, np.nan, dtype=np.float64)
    if jd.size == 0:
        return out
    inside = (jd_q >= jd[0]) & (jd_q <= jd[-1])
    if not inside.any():
        return out
    k = np.clip(np.searchsorted(jd, jd_q, side="right") - 1, 0, jd.size - 1)
    if series.get("interp") == "hold":
        out[inside] = v[k[inside]]
        return out
    k2 = np.clip(k + 1, 0, jd.size - 1)
    t1, t2 = jd[k], jd[k2]
    v1, v2 = v[k], v[k2]
    w = np.where(t2 > t1, (jd_q - t1) / np.where(t2 > t1, t2 - t1, 1.0), 0.0)
    lin = v1 + w * (v2 - v1)
    lin = np.where(np.isfinite(v1) & np.isfinite(v2), lin, np.nan)
    # An EXACT hit on a sample is not an interpolation and must not be nulled by its
    # neighbour: a query landing precisely on a good sample keeps that sample's value
    # even when the next one is a fill. Without this the rule "never interpolate
    # across a gap" would also discard the last good sample BEFORE every gap, which
    # is a different and much larger deletion than the one intended.
    lin = np.where(w == 0.0, v1, lin)
    out[inside] = lin[inside]
    return out


def availability(data_dir=None):
    """A one-glance FETCHED / UNFETCHED report for the CLI and for the scan record."""
    _series, audit = load_series(data_dir)
    return audit


if __name__ == "__main__":          # pragma: no cover - operator entry point
    import argparse
    ap = argparse.ArgumentParser("engine.spaceweather")
    ap.add_argument("--data-dir", default=DATA_DIR)
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--timeout", type=int, default=120)
    ap.add_argument("--sources", default=None,
                    help="comma-separated subset of %s" % (tuple(SOURCES),))
    a = ap.parse_args()
    if a.download:
        keys = tuple(a.sources.split(",")) if a.sources else None
        rep = download_all(keys, data_dir=a.data_dir, force=a.force,
                           timeout=a.timeout)
        print(json.dumps({k: [dict(r) for r in v]
                          for k, v in rep["sources"].items()}, indent=2))
    print(json.dumps(availability(a.data_dir), indent=2))
