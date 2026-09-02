"""NON-TIDAL FORCING AT EVENT TIME -- the block-3 covariate loader. Priced 0.

WHAT THIS IS. A per-event covariate block built from the three dated forcing series that
have been on disk since 2026-08-22 and used by nothing: IERS Earth-orientation (LOD and
polar motion, daily, 1962->now), GFZ Kp/ap (3-hourly, 1932->now) and NASA/GSFC OMNI2
(hourly solar wind and indices, 1995->2026). It exists so that `exp_learned_ext.py` can
fit a nested model on a covariate class OTHER than the tides, on exactly the same matched
case-control design.

WHY THESE TWO VERY DIFFERENT THINGS SHARE A BLOCK, and why the block is honest about it:

  * LOD and polar motion are a GENUINE, SMALL mechanical channel. Earth-rotation variations
    load the lithosphere; the Bendick & Bilham decadal claim is a specific published
    hypothesis this program holds the data to test.
  * Kp / ap / OMNI are a MATCHED PLACEBO. The physical coupling from geomagnetic activity
    to fault stress is orders of magnitude below tidal. Their defensible use is to measure
    how often this pipeline manufactures survivors on a property class that has REAL
    temporal structure and a near-zero physical prior (playbook P-2.4). A survivor here is
    a false positive by construction, and its rate is the calibration number.

FILL POLICY, WHICH IS THE WHOLE REASON THIS IS A MODULE AND NOT SIX LINES INLINE.
Nothing is ever zero-filled. A time outside a series' coverage, or inside one of its data
gaps, yields NaN in the value column AND a 1 in that source's explicit `miss_*` indicator
column. The indicator is a FEATURE handed to the model, so the model cannot exploit
missingness without the audit seeing it, and `event_block` returns a per-series count of
exactly how many event rows were filled. `HistGradientBoostingClassifier` consumes NaN
natively, so no imputation happens anywhere in the chain.

TIME CONVENTION, WHICH IS NOT COSMETIC HERE. Every function takes JULIAN DAY directly.
The forcing series are DATED: an epoch error does not degrade gracefully into a weaker
signal, it silently walks the query off the end of the data (OMNI2 begins 1995). Callers
converting from a days-since-epoch catalogue column must pass the epoch explicitly through
`jd_from_days_since`. See the EPOCH AUDIT in `exp_learned_ext.py`.

ARTIFACTS THAT COULD FAKE A RESULT FROM THIS BLOCK, named before the run:
  1. HUMAN-SCHEDULE LEAKAGE. Kp is derived from magnetometer observatories on a human
     operating schedule and is reported on a UT 3-hour grid; if it carries weekly or
     diurnal structure of its own it can proxy the catalogue's own cultural-noise
     day/night artifact. `schedule_audit()` MEASURES this rather than assuming it away.
  2. COUNT-PATH PERIODICITY. Kp and ap are quantised (Kp to thirds), so their empirical
     distributions are step functions and any statistic sensitive to ties can read the
     plateaus as structure.
  3. LONG-PERIOD TREND CONFOUNDING. LOD has decadal and interannual trends; catalogue
     completeness also drifts on those scales. In the matched design each stratum's
     controls come from the event's OWN +/- fortnight, so epoch is matched to a fortnight
     and the decadal channel is largely differenced away -- which is also why this design
     CANNOT test the decadal LOD claim, only a sub-fortnight one. Residual exposure: LOD
     changes measurably across a fortnight (seasonal and tidal-in-LOD terms), so a trend
     that is steep at the scale of two weeks is not removed.
  4. THE DAY/NIGHT DETECTION ARTIFACT. Roughly 2 percent of catalogued events are displaced
     into local night by cultural noise (`exp_diurnal_discriminator.py`). Model A is handed
     local solar hour, day of week and day of year for free, so any block-3 gain must be
     beyond it. OMNI/Kp are UT-referenced and SoCal is at a fixed UT offset, so a UT-grid
     periodicity can alias onto local solar hour: that is exactly what artifact 1 measures.
"""

from __future__ import annotations

import datetime as _dt
import os

import numpy as np

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "data", "spaceweather")

MJD_EPOCH_JD = 2400000.5
UNIX_EPOCH_JD = 2440587.5          # JD of 1970-01-01T00:00:00Z

EOP_FILE = "eopc04_14_IAU2000.62-now.csv"
KP_FILE = "Kp_ap_since_1932.txt"

# OMNI2 hourly: 0-BASED token index -> (name, fill sentinel). Word numbers are from the
# OMNI2 format description; a bare integer index in a parser is the kind of line that is
# wrong for a year before anyone notices, so both are written down.
#   word  9 -> idx  8  scalar field magnitude average |B|, nT
#   word 17 -> idx 16  Bz GSM, nT
#   word 24 -> idx 23  proton density, n/cc
#   word 25 -> idx 24  plasma flow speed, km/s
#   word 41 -> idx 40  Dst, nT
OMNI_COLS = {
    "imf_b":          (8,  999.9),
    "bz_gsm":         (16, 999.9),
    "proton_density": (23, 999.9),
    "sw_speed":       (24, 9999.0),
    "dst":            (40, 99999.0),
}
OMNI_N_WORDS = 55

OMNI_INSTANT = tuple(OMNI_COLS)                       # value at the prior hour
OMNI_MEAN24 = tuple(k + "_24h" for k in OMNI_COLS)    # trailing 24 h running mean

EOP_FEATURES = ("lod", "lod_d1", "lod_d7",
                "x_pole", "y_pole", "pm_radius", "pm_rate")
KP_FEATURES = ("kp", "ap")
BLOCK3_VALUES = EOP_FEATURES + KP_FEATURES + OMNI_INSTANT + OMNI_MEAN24
BLOCK3_FLAGS = ("miss_eop", "miss_kp", "miss_omni")
BLOCK3 = BLOCK3_VALUES + BLOCK3_FLAGS


# ------------------------------------------------------------------ time helpers --
def jd_from_days_since(t_days, epoch):
    """JD from a days-since-`epoch` column. The epoch is EXPLICIT and mandatory."""
    e = epoch if epoch.tzinfo else epoch.replace(tzinfo=_dt.timezone.utc)
    d70 = (e - _dt.datetime(1970, 1, 1, tzinfo=_dt.timezone.utc)).total_seconds() / 86400.0
    return np.asarray(t_days, dtype=np.float64) + d70 + UNIX_EPOCH_JD


def _jd_from_ymdh(y, mo, d, hours):
    """Vectorised Julian Day from calendar fields (same formula as engine.ephemeris)."""
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


# ---------------------------------------------------------------------- parsers --
def parse_iers_eop(path):
    """IERS EOP 14 C04 CSV (';' separated) -> (jd, {x_pole, y_pole, lod}), by NAME.

    Columns are taken by header name, never by position: the C04 product has gained
    columns over its life and a positional parser would silently read the wrong one.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        head = fh.readline().strip().split(";")
        idx = {}
        for want, key in (("MJD", "mjd"), ("x_pole", "x_pole"),
                          ("y_pole", "y_pole"), ("LOD", "lod")):
            if want not in head:
                raise ValueError("IERS C04 header lacks %r; header=%r" % (want, head))
            idx[key] = head.index(want)
        rows = {k: [] for k in idx}
        for line in fh:
            f = line.rstrip("\n").split(";")
            if len(f) <= max(idx.values()):
                continue
            try:
                vals = {k: (float(f[i]) if f[i].strip() else np.nan)
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


def parse_kp_ap(path):
    """GFZ 3-hourly Kp/ap -> (jd, kp, ap). Missing = -1.000 / -1 -> NaN.

    File's own declared format: `YYYY MM DD hh.h hh._m days days_m Kp ap D`, where `hh.h`
    is the START of the 3-hour interval. The value is a step function over that interval,
    which is why the sampler is a zero-order hold and not a linear interpolation.
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
    kp = np.asarray(kp, dtype=np.float64)
    ap = np.asarray(ap, dtype=np.float64)
    kp[kp < 0] = np.nan
    ap[ap < 0] = np.nan
    return jd, kp, ap


def parse_omni2(paths):
    """OMNI2 hourly year files -> (jd, {name: values}, stats). Fill sentinels -> NaN.

    Each row carries 55 whitespace-separated words. The count is CHECKED per row rather
    than assumed, because a silently short row would shift every column index by one and
    produce a plausible-looking series of the wrong quantity.
    """
    paths = sorted(paths)
    jd_all, cols = [], {k: [] for k in OMNI_COLS}
    n_rows, n_short = 0, 0
    jan1 = {}
    for p in paths:
        with open(p, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                f = line.split()
                if len(f) < OMNI_N_WORDS:
                    if line.strip():
                        n_short += 1
                    continue
                yr, doy, hr = int(f[0]), int(f[1]), int(f[2])
                if yr not in jan1:
                    jan1[yr] = float(_jd_from_ymdh(np.array([yr]), np.array([1]),
                                                   np.array([1]), np.array([0.0]))[0])
                jd_all.append(jan1[yr] + (doy - 1) + hr / 24.0)
                for k, (i, fill) in OMNI_COLS.items():
                    v = float(f[i])
                    cols[k].append(np.nan if v >= fill else v)
                n_rows += 1
    if n_rows == 0:
        raise ValueError("no OMNI2 rows parsed from %d file(s)" % len(paths))
    jd = np.asarray(jd_all, dtype=np.float64)
    o = np.argsort(jd, kind="stable")
    return (jd[o], {k: np.asarray(v, dtype=np.float64)[o] for k, v in cols.items()},
            {"n_files": len(paths), "n_rows": n_rows, "n_short_rows": n_short})


# ------------------------------------------------------------------- the sampler --
def sample_hold(jd_grid, values, jd_q):
    """Zero-order hold: the value of the sample interval each query falls in.

    NaN outside coverage. NEVER interpolates across a gap and never substitutes.
    """
    jd_q = np.asarray(jd_q, dtype=np.float64)
    out = np.full(jd_q.shape, np.nan, dtype=np.float64)
    if jd_grid.size == 0:
        return out
    inside = (jd_q >= jd_grid[0]) & (jd_q <= jd_grid[-1])
    if not inside.any():
        return out
    k = np.clip(np.searchsorted(jd_grid, jd_q, side="right") - 1, 0, jd_grid.size - 1)
    out[inside] = values[k[inside]]
    return out


def sample_linear(jd_grid, values, jd_q):
    """Linear interpolation between bracketing samples; NaN outside coverage.

    A query landing EXACTLY on a good sample keeps that sample even when its neighbour is
    a fill -- otherwise the rule "never interpolate across a gap" would also discard the
    last good sample before every gap, a much larger deletion than the one intended.
    """
    jd_q = np.asarray(jd_q, dtype=np.float64)
    out = np.full(jd_q.shape, np.nan, dtype=np.float64)
    if jd_grid.size == 0:
        return out
    inside = (jd_q >= jd_grid[0]) & (jd_q <= jd_grid[-1])
    if not inside.any():
        return out
    k = np.clip(np.searchsorted(jd_grid, jd_q, side="right") - 1, 0, jd_grid.size - 1)
    k2 = np.clip(k + 1, 0, jd_grid.size - 1)
    t1, t2 = jd_grid[k], jd_grid[k2]
    v1, v2 = values[k], values[k2]
    w = np.where(t2 > t1, (jd_q - t1) / np.where(t2 > t1, t2 - t1, 1.0), 0.0)
    lin = v1 + w * (v2 - v1)
    lin = np.where(np.isfinite(v1) & np.isfinite(v2), lin, np.nan)
    lin = np.where(w == 0.0, v1, lin)
    out[inside] = lin[inside]
    return out


def trailing_mean(values, window):
    """NaN-aware trailing mean over `window` samples INCLUDING the current one.

    Returns NaN where the window holds no finite sample. Requires a regular grid; the
    caller asserts regularity (`load_forcing` measures the step and raises otherwise).
    """
    v = np.asarray(values, dtype=np.float64)
    ok = np.isfinite(v)
    s = np.concatenate([[0.0], np.cumsum(np.where(ok, v, 0.0))])
    c = np.concatenate([[0.0], np.cumsum(ok.astype(np.float64))])
    i = np.arange(v.size) + 1
    lo = np.maximum(0, i - window)
    num = s[i] - s[lo]
    den = c[i] - c[lo]
    return np.where(den > 0, num / np.maximum(den, 1.0), np.nan)


def jd_to_iso(jd):
    return (_dt.datetime(1970, 1, 1, tzinfo=_dt.timezone.utc)
            + _dt.timedelta(days=float(jd) - UNIX_EPOCH_JD)).isoformat()


# --------------------------------------------------------------------- the load --
def load_forcing(data_dir=None):
    """Parse all three sources. Returns (series, audit) with MEASURED counts.

    A source that is absent or unparseable raises; it is never substituted, zero-filled,
    or silently skipped.
    """
    d = data_dir or DATA_DIR
    series, audit = {}, {"data_dir": d, "sources": {}}

    p = os.path.join(d, EOP_FILE)
    jd, cols = parse_iers_eop(p)
    step = np.unique(np.round(np.diff(jd), 9))
    series["eop"] = dict(jd=jd, **cols)
    audit["sources"]["iers_eop"] = {
        "file": EOP_FILE, "n_rows": int(jd.size),
        "jd_min": float(jd.min()), "jd_max": float(jd.max()),
        "utc_min": jd_to_iso(jd.min()), "utc_max": jd_to_iso(jd.max()),
        "step_days_unique": [float(x) for x in step[:5]],
        "n_nan": {k: int(np.isnan(v).sum()) for k, v in cols.items()}}

    p = os.path.join(d, KP_FILE)
    jdk, kp, ap = parse_kp_ap(p)
    stepk = np.unique(np.round(np.diff(jdk), 6))
    series["kp"] = {"jd": jdk, "kp": kp, "ap": ap}
    audit["sources"]["kp_ap"] = {
        "file": KP_FILE, "n_rows": int(jdk.size),
        "jd_min": float(jdk.min()), "jd_max": float(jdk.max()),
        "utc_min": jd_to_iso(jdk.min()), "utc_max": jd_to_iso(jdk.max()),
        "step_days_unique": [float(x) for x in stepk[:5]],
        "n_nan": {"kp": int(np.isnan(kp).sum()), "ap": int(np.isnan(ap).sum())}}

    paths = [os.path.join(d, f) for f in sorted(os.listdir(d))
             if f.startswith("omni2_") and f.endswith(".dat")]
    jdo, ocols, ostat = parse_omni2(paths)
    stepo = np.unique(np.round(np.diff(jdo), 6))
    regular = bool(stepo.size == 1 and abs(stepo[0] - 1.0 / 24.0) < 1e-6)
    means = {k + "_24h": trailing_mean(v, 24) for k, v in ocols.items()}
    series["omni"] = dict(jd=jdo, **ocols, **means)
    audit["sources"]["omni2"] = {
        "n_files": ostat["n_files"], "n_rows": ostat["n_rows"],
        "n_short_rows": ostat["n_short_rows"],
        "jd_min": float(jdo.min()), "jd_max": float(jdo.max()),
        "utc_min": jd_to_iso(jdo.min()), "utc_max": jd_to_iso(jdo.max()),
        "step_days_unique": [float(x) for x in stepo[:5]],
        "grid_is_regular_hourly": regular,
        "n_nan": {k: int(np.isnan(v).sum()) for k, v in ocols.items()}}
    if not regular:
        raise ValueError("OMNI2 grid is not regular hourly; trailing means would be "
                         "misaligned. measured steps=%r" % (stepo[:10],))
    return series, audit


# ------------------------------------------------------------------ the block --
def event_block(jd_q, series):
    """Block 3 at arbitrary Julian days. Returns (cols, fill_report).

    `cols` maps every name in BLOCK3 to a float array. Values are NaN where the source
    does not cover the time; the three `miss_*` flags are 1.0 exactly there. Nothing is
    imputed. `fill_report` counts, per feature, how many query rows are NaN.
    """
    jd_q = np.asarray(jd_q, dtype=np.float64)
    e, k, o = series["eop"], series["kp"], series["omni"]
    cols = {}

    lod = sample_linear(e["jd"], e["lod"], jd_q)
    cols["lod"] = lod
    cols["lod_d1"] = lod - sample_linear(e["jd"], e["lod"], jd_q - 1.0)
    cols["lod_d7"] = lod - sample_linear(e["jd"], e["lod"], jd_q - 7.0)
    xp = sample_linear(e["jd"], e["x_pole"], jd_q)
    yp = sample_linear(e["jd"], e["y_pole"], jd_q)
    cols["x_pole"], cols["y_pole"] = xp, yp
    r_now = np.hypot(xp, yp)
    r_prev = np.hypot(sample_linear(e["jd"], e["x_pole"], jd_q - 1.0),
                      sample_linear(e["jd"], e["y_pole"], jd_q - 1.0))
    cols["pm_radius"] = r_now
    cols["pm_rate"] = r_now - r_prev

    cols["kp"] = sample_hold(k["jd"], k["kp"], jd_q)
    cols["ap"] = sample_hold(k["jd"], k["ap"], jd_q)

    for name in OMNI_INSTANT + OMNI_MEAN24:
        cols[name] = sample_hold(o["jd"], o[name], jd_q)

    cols["miss_eop"] = (~np.isfinite(cols["lod"])).astype(np.float64)
    cols["miss_kp"] = (~np.isfinite(cols["kp"])).astype(np.float64)
    cols["miss_omni"] = (~np.isfinite(cols["sw_speed"])).astype(np.float64)

    for c in BLOCK3:
        if np.isinf(cols[c]).any():
            raise ValueError("infinite value produced in column %r" % c)

    fill = {"n_rows": int(jd_q.size),
            "per_feature_nan": {c: int((~np.isfinite(cols[c])).sum())
                                for c in BLOCK3_VALUES},
            "n_rows_any_source_missing":
                int((cols["miss_eop"] + cols["miss_kp"] + cols["miss_omni"] > 0).sum()),
            "policy": ("NaN is carried into the model (HistGradientBoosting consumes it "
                       "natively) and is flagged by an explicit miss_* indicator column. "
                       "No zero-fill, no imputation, no interpolation across a gap.")}
    return cols, fill


# ------------------------------------------------------------- artifact measures --
def schedule_audit(series, jd_lo, jd_hi):
    """ARTIFACTS 1 and 2, MEASURED rather than assumed, over the arm's own time span.

    Reports (a) the amplitude of the UT-diurnal and weekly cycles in Kp, ap and the OMNI
    instants as a fraction of the series sd, and (b) the number of distinct values, which
    exposes the Kp/ap quantisation. A UT-diurnal cycle in a covariate is the mechanism by
    which it could proxy the catalogue's own local-solar-hour detection artifact.
    """
    out = {}
    for src, names in (("kp", KP_FEATURES), ("omni", OMNI_INSTANT)):
        s = series[src]
        m = (s["jd"] >= jd_lo) & (s["jd"] <= jd_hi)
        jd = s["jd"][m]
        frac = np.mod(jd + 0.5, 1.0)                       # UT fraction of day
        dow = np.mod(np.floor(jd + 0.5).astype(np.int64), 7)
        for nm in names:
            v = s[nm][m]
            ok = np.isfinite(v)
            if ok.sum() < 100:
                continue
            vv, ff, dd = v[ok], frac[ok], dow[ok]
            mu, sd = float(np.mean(vv)), float(np.std(vv))
            c, sn = np.cos(2 * np.pi * ff), np.sin(2 * np.pi * ff)
            amp_ut = float(2.0 * np.hypot(np.mean((vv - mu) * c),
                                          np.mean((vv - mu) * sn)))
            wk = np.array([vv[dd == j].mean() if (dd == j).sum() else np.nan
                           for j in range(7)])
            amp_wk = float(np.nanmax(wk) - np.nanmin(wk))
            out[nm] = {
                "n": int(ok.sum()), "mean": mu, "sd": sd,
                "ut_diurnal_amplitude": amp_ut,
                "ut_diurnal_amplitude_over_sd": (amp_ut / sd) if sd > 0 else None,
                "weekday_peak_to_trough": amp_wk,
                "weekday_peak_to_trough_over_sd": (amp_wk / sd) if sd > 0 else None,
                "n_distinct_values": int(np.unique(vv).size)}
    return out
