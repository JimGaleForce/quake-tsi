"""`mine` mode -- the EQ-24 pattern miner. A GENERATOR, NEVER EVIDENCE.

What it does: sweeps many candidate features x lags against the ETAS-residualized
global daily target, applies multiplicity discipline (surrogate nulls, BH-FDR, a
harmonic ladder and an aliasing audit), and emits ranked hypothesis stubs.

What it does NOT do: produce a result anybody may quote. Every output file says so,
in the first ten lines. A stub that survives here has earned exactly one thing: the
right to be written down as a hypothesis, which then walks the ordinary
Kepler -> Merton floor -> Popper ruling -> EQ-23 frozen-holdout path like anything
else. Nothing here spends a holdout hash.

Design (ruled 2026-08-11, EQ-24):

  target      daily GLOBAL (domain-wide) event occurrence vs the ETAS baseline:
              r_t = sum_c y(c,t) - sum_c lambda_etas(c,t), over the EXPLORATION
              window minus the 365-day ETAS burn-in. Plus per-event marks
              (magnitude, depth) for feature-vs-mark tests. v1 is temporal/global
              only -- no per-cell spatial mining.

  features    (1) ephemeris, zero-download, deterministic in t;
              (2) beats/interactions of those cycles (the moire rule: the beat
                  nobody computed is a first-class feature, not a nuisance);
              (3) optional downloads (F10.7 / Kp / Ap / LOD), graceful skip;
              (4) causal catalogue-derived global marks (trailing b-value, deep
                  fraction, mean depth).

  scoring     (a) Poisson GLM of daily domain counts on the feature with
                  log(sum lambda_etas) as offset;
              (b) Lomb-Scargle period scan of r_t, 2..4000 d, with a harmonic
                  ladder {P/3, P/2, P, 2P, 3P};
              (c) mark tests (Spearman; circular-linear for phase features).

  discipline  circular-shift surrogate nulls, BH-FDR at q=0.1 across the whole
              session, the harmonic ladder, and an aliasing audit at 2 d / 7 d
              binning plus unbinned event times.

One honest deviation from the ruling, flagged rather than hidden, and MEASURED
(engine/tests/test_mine.py::test_block_bootstrap_null_is_calibrated_and_shift_null_is_not):

  A circular time shift of an evenly sampled series does not change its power
  spectrum at all -- a shift is a pure phase rotation. Two of this miner's three
  scoring paths are, at bottom, measurements of that spectrum at a fixed frequency:
  the Lomb-Scargle period scan explicitly, and the 2-df (sin, cos) GLM test of a
  DETERMINISTIC CYCLE implicitly. For those, circular-shift surrogates inherit the
  very signal they are meant to destroy. Measured on a synthetic 8% modulation at
  29.53 d over 7716 days: observed chi2 = 152, median circular-shift surrogate = 97,
  and the shift null returns p = 0.001 only because the observed value happens to
  top a distribution it already contaminated; at a 2% modulation it returns a
  meaningless p = 0.22 for a real effect.

  So the null is chosen per test rather than uniformly:
    * periodic features (families 1-2) -> circular MOVING-BLOCK BOOTSTRAP of the
      (counts, offset) pair, block length >= 2 cycles of the feature under test;
    * aperiodic features (families 3-4) -> the MORE CONSERVATIVE of circular shift
      (as ruled) and the same block bootstrap;
    * mark tests -> the more conservative of circular shift and a block bootstrap
      of the mark series along the event sequence;
    * period scan -> an AR(1) red-noise null matched to the residual's own lag-1
      autocorrelation and a permutation (white) null, more conservative of the two.
  Both p-values are recorded for every test, so the choice is auditable.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import math
import os
import time

import numpy as np

from . import baseline as bl, datasets, design, ephemeris as eph, splits
from . import gpd_tail

MINE_DIR = os.path.join("engine", "out", "mine")
DATA_LOG = os.path.join(MINE_DIR, "data_log.jsonl")

FDR_Q = 0.10
LN2 = math.log(2.0)

GENERATOR_NOT_EVIDENCE = (
    "GENERATOR, NOT EVIDENCE. Nothing in this file is a result. It is a ranked list "
    "of hypotheses produced by an automated sweep over the EXPLORATION window only, "
    "with in-sample effect sizes. No line here may be quoted, cited, or reported as "
    "a finding. A candidate becomes evidence only by walking the ordinary path: "
    "written up as a K-entry, given a Popper ruling, and then tested ONCE against "
    "the frozen holdout (which costs a hash). This run spent no holdout hash."
)

AMPLITUDE_HONESTY = (
    "EXPECTED-AMPLITUDE FRAMING. The tidal corpse is the calibration: kPa-scale "
    "forcing on a fault that is already critically stressed yields BOUNDS, not "
    "detections. Almost every exotic forcing in this feature set is weaker than the "
    "solid-earth tide. The expected product of this hypothesis is therefore an upper "
    "bound, or at most a state-dependent interaction (the K-077 shape) -- not a "
    "standalone predictive signal. Read the quoted rate modulation as the size of "
    "the thing that must be excluded, not as a forecast."
)

# ------------------------------------------------ K-089-R (offset/lag scan) ---
# Popper's §P5-5 ruling, executed. The free-scan exception is per FEATURE and per
# AXIS, never per feature-kind:
#   * the OFFSET axis is exactly free for all 9 kind='phase' features (a phase offset
#     is an orthogonal rotation of [sin, cos] and the 2-df quadratic form is exactly
#     invariant), so the scan is complete at offset 0 at zero multiplicity cost;
#   * the LAG axis is free only where theta is LINEAR in t. Four phase features are
#     built from mean arguments (or close enough); the other five are built on true
#     solar/lunar longitudes (Meeus, 6.29 deg equation of centre) and a lag is NOT a
#     rotation for them -- measured min R^2 down to 0.905 at lag 15 d;
#   * the LADDER-RUNG axis is never free (a harmonic is a new feature).
# Tranche 1 is therefore the lag scan over the 13 lag-unscanned cyclic features:
# 8 kind='linear' (240 new tests) + 5 non-lag-free kind='phase' (150 new tests).
TRANCHE1_LABEL = "K-089-R tranche 1"
TRANCHE1_LAGS = tuple(range(0, 31))          # 0..30 d inclusive -> 30 NEW lags each

# The four phase features whose lag axis is provably free (min R^2 >= 1 - 1e-4
# against the lag-0 column space); see LAG_FREE_TOL and lag_invariance_audit().
LAG_FREE_PHASE_FEATURES = (
    "moon_anomalistic_phase", "moon_draconic_phase", "half_draconic_phase",
    "annual_phase",
)
LAG_FREE_TOL = 1e-4

AXIS_AUDIT_LINE = (
    "AXIS AUDIT (K-089, the clause that survives unamended). Axes scanned in this "
    "session, axes provably invariant, and axes neither -- stated so no reader has to "
    "infer it from the config."
)

MIGNAN_BROCCARDO = (
    "the known ML-mining failure mode from M-008 (DeVries Nature net beaten by "
    "2-parameter logistic -- Mignan & Broccardo 2019) quoted in the README as the "
    "standing warning: prefer the boring model; a boosted-tree find that a 2-param "
    "GLM can't reproduce after conditioning is treated as leakage until proven "
    "otherwise."
)


# ============================================================== features ===
class Feature:
    """One candidate feature.

    kind='phase'  -> `values` is an angle in radians; the design is [sin, cos] (2 df)
    kind='linear' -> `values` is a real series; the design is one standardised column
    """

    def __init__(self, name, family, kind, values, describe="", lags=(0,),
                 causality_exempt=False, period_hint=None):
        self.name = name
        self.family = int(family)
        self.kind = kind
        self.values = np.asarray(values, dtype=np.float64)
        self.describe = describe
        self.lags = tuple(int(l) for l in lags)
        self.causality_exempt = bool(causality_exempt)
        self.period_hint = period_hint

    @property
    def periodic(self):
        """True for the deterministic ephemeris cycles (families 1 and 2).

        For these the circular-shift null is provably powerless (a shift preserves
        the target's power spectrum, which is exactly what a fixed-frequency test
        measures), so the block bootstrap is the primary null. See
        `score_stat_block_bootstrap`.
        """
        return self.family in (1, 2)

    @property
    def block_days(self):
        """Bootstrap block length, set by the FEATURE'S OWN timescale.

        The block bootstrap is only a valid null if its blocks are long compared to
        the structure being tested: too short and it destroys the feature's own
        timescale in the surrogates, which makes the null too narrow and the test
        ANTI-CONSERVATIVE. That is not hypothetical -- an earlier build used a flat
        90-day block and handed F10.7 (an 11-year solar cycle) a p at the resolution
        floor with z = 32, purely because 90-day blocks cannot represent an 11-year
        cycle. So the length is measured, not guessed:

          * cyclic features       -> 2 x the known period;
          * everything else       -> 4 x the e-folding time of the feature's own
                                     autocorrelation function.

        Clipped to [30, 800] days: 800 d over a 7716-day window is only ~10 blocks,
        which is exactly as brutal as a decadal claim on 21 years of data deserves.
        """
        if self.period_hint:
            base = 2.0 * float(self.period_hint)
        else:
            v = self.values - self.values.mean()
            denom = float(np.dot(v, v))
            if denom <= 0:
                return 90.0
            tau = 1.0
            for lag in range(1, min(1200, v.size // 4)):
                if float(np.dot(v[lag:], v[:-lag])) / denom < math.exp(-1.0):
                    tau = float(lag)
                    break
                tau = float(lag)
            base = 4.0 * tau
        return float(np.clip(base, 30.0, 800.0))

    def design(self, sl: slice, lag: int = 0) -> np.ndarray:
        """(n_window, k) design columns for a day window, shifted back by `lag` days."""
        start, stop = sl.start - lag, sl.stop - lag
        if start < 0:
            raise ValueError(f"lag {lag} runs off the head of the record")
        v = self.values[start:stop]
        if self.kind == "phase":
            return np.column_stack([np.sin(v), np.cos(v)])
        sd = v.std()
        if sd <= 0:
            raise ValueError(f"{self.name}: constant over the window (0 variance)")
        return ((v - v.mean()) / sd)[:, None]

    @property
    def df(self):
        return 2 if self.kind == "phase" else 1


def ephemeris_features(t0, n_days, lags=(0,), lag_features=None):
    """Families 1 and 2: deterministic functions of t. No downloads, ever.

    `lags` defaults to (0,) -- the historical behaviour, byte-for-byte -- so every
    existing caller keeps a single lag-0 test per cyclic feature.

    `lag_features` selects WHICH features receive the grid (K-089-R: the exception is
    per feature and per axis, so a blanket grid would pay multiplicity for tests that
    are provably identical):
      * None        -> every family-1/2 feature gets `lags`;
      * "tranche1"  -> only the 13 lag-unscanned cyclic features (the 8 kind='linear'
                       plus the 5 kind='phase' whose lag is not a rotation). The four
                       features in LAG_FREE_PHASE_FEATURES stay at lag 0: their
                       invariance is a theorem, and scanning them is waste, not rigour;
      * an iterable of feature names -> exactly those.
    """
    e = eph.ephemeris_table(t0, n_days)
    syn, ano, dra, ann = (e["synodic_rad"], e["anomalistic_rad"],
                          e["draconic_rad"], e["annual_rad"])
    f = []
    add = f.append

    # ---- family 1: the cycles themselves -------------------------------
    add(Feature("moon_synodic_phase", 1, "phase", syn,
                "lunar phase angle (0 = new moon)", period_hint=eph.SYNODIC_MONTH,
                causality_exempt=True))
    add(Feature("moon_anomalistic_phase", 1, "phase", ano,
                "lunar anomalistic phase (0 = perigee)",
                period_hint=eph.ANOMALISTIC_MONTH, causality_exempt=True))
    add(Feature("moon_draconic_phase", 1, "phase", dra,
                "lunar draconic phase (0 = ascending node)",
                period_hint=eph.DRACONIC_MONTH, causality_exempt=True))
    add(Feature("annual_phase", 1, "phase", ann,
                "solar ecliptic longitude (0 = March equinox)",
                period_hint=eph.TROPICAL_YEAR, causality_exempt=True))
    add(Feature("moon_distance", 1, "linear", e["moon_dist_km"],
                "lunar geocentric distance (km)", causality_exempt=True,
                period_hint=eph.ANOMALISTIC_MONTH))
    add(Feature("moon_declination", 1, "linear", e["moon_dec_deg"],
                "lunar declination (deg)", causality_exempt=True,
                period_hint=eph.DRACONIC_MONTH))
    add(Feature("moon_abs_declination", 1, "linear", np.abs(e["moon_dec_deg"]),
                "|lunar declination| (deg; the 13.66 d declination-tide envelope)",
                causality_exempt=True, period_hint=eph.DRACONIC_MONTH / 2))
    add(Feature("sun_declination", 1, "linear", e["sun_dec_deg"],
                "solar declination (deg)", causality_exempt=True,
                period_hint=eph.TROPICAL_YEAR))
    add(Feature("sun_moon_elongation", 1, "linear", e["elongation_deg"],
                "sun-moon elongation (deg)", causality_exempt=True,
                period_hint=eph.SYNODIC_MONTH))

    # ---- family 2: beats and interactions (the moire rule) --------------
    add(Feature("spring_neap_phase", 2, "phase", np.mod(2 * syn, 2 * np.pi),
                "spring-neap beat = 2 x synodic (14.765 d)",
                period_hint=eph.SYNODIC_MONTH / 2, causality_exempt=True))
    add(Feature("half_draconic_phase", 2, "phase", np.mod(2 * dra, 2 * np.pi),
                "2 x draconic (13.606 d; the declination-tide beat)",
                period_hint=eph.DRACONIC_MONTH / 2, causality_exempt=True))
    add(Feature("perigean_spring_beat", 2, "phase", np.mod(syn - ano, 2 * np.pi),
                "synodic - anomalistic beat (411.8 d perigean-spring cycle)",
                period_hint=1.0 / abs(1 / eph.SYNODIC_MONTH - 1 / eph.ANOMALISTIC_MONTH),
                causality_exempt=True))
    add(Feature("eclipse_year_beat", 2, "phase", np.mod(syn - dra, 2 * np.pi),
                "synodic - draconic beat (the 173.3 d eclipse half-year)",
                period_hint=1.0 / abs(1 / eph.SYNODIC_MONTH - 1 / eph.DRACONIC_MONTH),
                causality_exempt=True))
    add(Feature("annual_synodic_beat", 2, "phase",
                np.mod(ann - syn, 2 * np.pi),
                "annual - synodic beat", causality_exempt=True,
                period_hint=1.0 / abs(1 / eph.TROPICAL_YEAR - 1 / eph.SYNODIC_MONTH)))
    # interaction terms: products, entered as linear covariates
    add(Feature("perigee_syzygy", 2, "linear",
                np.cos(2 * syn) * (400000.0 / e["moon_dist_km"]) ** 3,
                "perigean-spring tidal proxy: cos(2*synodic) x (r0/r)^3",
                causality_exempt=True,
                period_hint=1.0 / abs(1 / eph.SYNODIC_MONTH
                                      - 1 / eph.ANOMALISTIC_MONTH)))
    add(Feature("declination_product", 2, "linear",
                e["moon_dec_deg"] * e["sun_dec_deg"],
                "lunar x solar declination product (the declination beat)",
                causality_exempt=True, period_hint=eph.TROPICAL_YEAR))
    add(Feature("tidal_potential_proxy", 2, "linear",
                (400000.0 / e["moon_dist_km"]) ** 3
                * (0.5 - 1.5 * np.sin(np.radians(e["moon_dec_deg"])) ** 2),
                "zonal lunar tidal-potential proxy (r^-3 x P2(sin dec))",
                causality_exempt=True,
                period_hint=1.0 / abs(1 / eph.DRACONIC_MONTH
                                      - 1 / eph.ANOMALISTIC_MONTH)))

    lags = tuple(int(l) for l in lags)
    if lags != (0,):
        if lag_features is None:
            sel = {ft.name for ft in f}
        elif lag_features == "tranche1":
            sel = set(tranche1_lag_features(f))
        else:
            sel = set(lag_features)
        unknown = sel - {ft.name for ft in f}
        if unknown:
            raise ValueError(f"lag_features names no such feature: {sorted(unknown)}")
        for ft in f:
            if ft.name in sel:
                ft.lags = lags
    return f


def tranche1_lag_features(feats):
    """The 13 K-089-R tranche-1 names, derived from the feature list, not hardcoded.

    8 kind='linear' cyclic features (the axis §K87-0(b') proved was never tested) plus
    the 5 kind='phase' features whose lag is not an exact rotation. Returns names in
    feature order.
    """
    return tuple(ft.name for ft in feats
                 if ft.family in (1, 2)
                 and (ft.kind == "linear"
                      or ft.name not in LAG_FREE_PHASE_FEATURES))


def lag_invariance_audit(feats, window, lags=(1, 3, 7, 15, 30)):
    """Clause 3 of K-089, restated by §P5-5: prove FREE numerically, before the scan.

    For every kind='phase' feature, regress the lag-L design [sin, cos] onto the lag-0
    design and take the MINIMUM R^2 over the two columns. A lag that is an exact
    rotation gives 1.0. Declaration comes from LAG_FREE_PHASE_FEATURES; the audit
    checks the declaration and reports `ok=False` on either failure mode Popper names:
    a feature declared FREE that measures below tolerance is a BUG; a feature declared
    PRICED that measures at 1.0 is a MISSED SAVING.
    """
    rows = []
    for ft in feats:
        if ft.kind != "phase" or ft.family not in (1, 2):
            continue
        A = ft.design(window, 0)
        A1 = np.column_stack([A, np.ones(A.shape[0])])
        pinv = np.linalg.pinv(A1)
        r2 = {}
        for L in lags:
            B = ft.design(window, int(L))
            fitb = A1 @ (pinv @ B)
            resid = B - fitb
            sst = ((B - B.mean(axis=0)) ** 2).sum(axis=0)
            r2[int(L)] = float(np.min(1.0 - (resid ** 2).sum(axis=0)
                                      / np.maximum(sst, 1e-300)))
        worst = min(r2.values())
        declared_free = ft.name in LAG_FREE_PHASE_FEATURES
        measured_free = worst >= 1.0 - LAG_FREE_TOL
        rows.append({
            "feature": ft.name, "declared": "FREE" if declared_free else "PRICED",
            "min_r2_by_lag": r2, "worst_min_r2": worst,
            "measured_free": measured_free,
            "ok": bool(declared_free == measured_free),
            "verdict": ("OK" if declared_free == measured_free else
                        ("BUG: declared FREE, measures PRICED" if declared_free else
                         "MISSED SAVING: declared PRICED, measures FREE")),
        })
    return rows


def scan_axis_audit(feats, tranche1):
    """The K-089 audit line's content: which axes were scanned, invariant, neither."""
    cyc = [f for f in feats if f.family in (1, 2)]
    scanned = [f.name for f in cyc if len(f.lags) > 1]
    free_lag = [f.name for f in cyc
                if len(f.lags) == 1 and f.name in LAG_FREE_PHASE_FEATURES]
    unscanned = [f.name for f in cyc
                 if len(f.lags) == 1 and f.name not in LAG_FREE_PHASE_FEATURES]
    n_new = sum(len(f.lags) - 1 for f in cyc)
    return {
        "tranche": TRANCHE1_LABEL if tranche1 else "none (lag 0 only for cyclic features)",
        "offset_axis": ("PROVABLY INVARIANT for all kind='phase' features (orthogonal "
                        "rotation of [sin, cos]; complete at offset 0, zero "
                        "multiplicity cost). NOT applicable to kind='linear'."),
        "lag_axis_scanned": scanned,
        "lag_axis_provably_free_not_scanned": free_lag,
        "lag_axis_neither_scanned_nor_free": unscanned,
        "ladder_rung_axis": ("NOT SCANNED for cyclic features (never free: a harmonic "
                             "is a new feature, not a rotation). Tranche 3, declared "
                             "at 68 tests, not run here."),
        "n_new_lag_tests_cyclic": int(n_new),
    }


def catalog_features(ctx, marks, lags):
    """Family 4: causal (strictly trailing) catalogue-derived global marks."""
    n_days = ctx.n_days
    day = marks["day"]
    mag = marks["mag"]
    dep = marks["depth"]

    def per_day(w):
        a = np.zeros((1, n_days), dtype=np.float64)
        np.add.at(a[0], day, w)
        return a

    n_d = per_day(np.ones(day.size))
    m_d = per_day(mag)
    z_d = per_day(dep)
    deep_d = per_day((dep > 70.0).astype(float))

    def trail(a, w):
        return design.EngineContext.trailing_sum(a.astype(np.float32), w)[0].astype(
            np.float64)

    def safe(num, den, fill):
        out = np.full(num.size, float(fill))
        ok = den > 0
        out[ok] = num[ok] / den[ok]
        return out

    mc, dm = float(ctx.meta.get("mag_floor", 4.5)), 0.1
    n90 = trail(n_d, 90)
    mean_m90 = safe(trail(m_d, 90), n90, mc + 0.3)
    # Aki (1965) MLE b-value; guard the degenerate mean == Mc - dM/2
    b90 = 1.0 / (math.log(10.0) * np.maximum(mean_m90 - (mc - dm / 2.0), 1e-3))
    b90 = np.clip(b90, 0.2, 3.0)

    n30 = trail(n_d, 30)
    deep30 = safe(trail(deep_d, 30), n30, 0.0)
    depth30 = safe(trail(z_d, 30), n30, 30.0)

    return [
        Feature("b_value_90d", 4, "linear", b90,
                "trailing 90 d global Aki b-value (STRICTLY causal)", lags=lags),
        Feature("deep_fraction_30d", 4, "linear", deep30,
                "trailing 30 d global share of events with z > 70 km (causal)",
                lags=lags),
        Feature("mean_depth_30d", 4, "linear", depth30,
                "trailing 30 d global mean hypocentral depth, km (causal)", lags=lags),
    ]


# =============================================== optional data downloads ===
DOWNLOADS = [
    {
        "key": "gfz_kp_ap_f107",
        "url": "https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_Ap_SN_F107_since_1932.txt",
        "gives": ["Ap_geomagnetic", "F107_solar_flux"],
    },
    {
        "key": "iers_lod",
        "url": "https://datacenter.iers.org/data/latestVersion/finals.all.iau2000.txt",
        "gives": ["length_of_day"],
    },
]


def _fetch(url, timeout=25):
    """Plain HTTPS GET with a certifi SSL context. Returns bytes or raises."""
    import ssl
    import urllib.request

    try:
        import certifi
        ctxt = ssl.create_default_context(cafile=certifi.where())
    except Exception:                                   # no certifi -> system store
        ctxt = ssl.create_default_context()
    req = urllib.request.Request(url, headers={"User-Agent": "EQ-23-engine/mine"})
    with urllib.request.urlopen(req, timeout=timeout, context=ctxt) as fh:
        return fh.read()


def _log_data(rec):
    os.makedirs(os.path.dirname(DATA_LOG) or ".", exist_ok=True)
    with open(DATA_LOG, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, sort_keys=True) + "\n")


def _parse_gfz(text, t0, n_days):
    """GFZ Kp/ap/Ap/SN/F10.7 daily file -> (Ap, F10.7) day-indexed arrays."""
    ap = np.full(n_days, np.nan)
    f107 = np.full(n_days, np.nan)
    for line in text.splitlines():
        if not line or line.lstrip().startswith("#"):
            continue
        f = line.split()
        if len(f) < 27:
            continue
        try:
            d = _dt.date(int(f[0]), int(f[1]), int(f[2]))
            i = (d - t0.date()).days
            if not (0 <= i < n_days):
                continue
            a, fx = float(f[23]), float(f[25])
        except (ValueError, IndexError):
            continue
        if a >= 0:
            ap[i] = a
        if fx > 0:
            f107[i] = fx
    return ap, f107


def _parse_iers_lod(text, t0, n_days):
    """IERS finals.all (IAU2000) -> length-of-day excess (ms), day-indexed."""
    lod = np.full(n_days, np.nan)
    base_mjd = eph.julian_day(t0.replace(hour=0, minute=0, second=0,
                                         microsecond=0)) - 2400000.5
    for line in text.splitlines():
        if len(line) < 86:
            continue
        try:
            mjd = float(line[7:15])
            v = line[79:86].strip()
            if not v:
                continue
            i = int(round(mjd - base_mjd))
            if 0 <= i < n_days:
                lod[i] = float(v)
        except ValueError:
            continue
    return lod


def _fill_gaps(a, name):
    """Forward/back fill a day-indexed series; returns (series, coverage_fraction)."""
    a = np.asarray(a, dtype=np.float64)
    ok = np.isfinite(a)
    cov = float(ok.mean())
    if cov == 0.0:
        raise ValueError(f"{name}: no usable values")
    idx = np.where(ok, np.arange(a.size), 0)
    np.maximum.accumulate(idx, out=idx)
    out = a[idx]
    first = int(np.argmax(ok))
    out[:first] = a[first]
    return out, cov


def download_features(t0, n_days, lags, enabled=True, timeout=25, verbose=True):
    """Family 3. Any failure logs a reason and returns fewer features -- never raises."""
    feats, log = [], []
    if not enabled:
        log.append({"key": "*", "status": "skipped", "reason": "--no-download"})
        if verbose:
            print("  family 3 (downloads): SKIPPED (--no-download)")
        return feats, log

    for src in DOWNLOADS:
        rec = {"key": src["key"], "url": src["url"],
               "ts": _dt.datetime.now(_dt.timezone.utc).isoformat()}
        t_a = time.time()
        try:
            raw = _fetch(src["url"], timeout=timeout)
            rec.update(status="ok", n_bytes=len(raw), seconds=round(time.time() - t_a, 2),
                       sha256=hashlib.sha256(raw).hexdigest())
            text = raw.decode("utf-8", errors="replace")
            frozen = os.path.join(MINE_DIR, "frozen", src["key"] + ".txt")
            os.makedirs(os.path.dirname(frozen), exist_ok=True)
            with open(frozen, "w", encoding="utf-8") as fh:
                fh.write(text)
            rec["frozen_copy"] = frozen.replace("\\", "/")

            if src["key"] == "gfz_kp_ap_f107":
                ap, f107 = _parse_gfz(text, t0, n_days)
                for nm, arr, desc in (
                    ("Ap_geomagnetic", ap, "GFZ daily Ap geomagnetic index"),
                    ("F107_solar_flux", f107, "Penticton F10.7 cm solar radio flux"),
                ):
                    s, cov = _fill_gaps(arr, nm)
                    s = np.roll(s, 1)            # strictly-before-t: use day t-1
                    s[0] = s[1]
                    feats.append(Feature(nm, 3, "linear", s,
                                         desc + " (lagged 1 d: strictly causal)",
                                         lags=lags))
                    rec.setdefault("coverage", {})[nm] = round(cov, 4)
            elif src["key"] == "iers_lod":
                lod = _parse_iers_lod(text, t0, n_days)
                s, cov = _fill_gaps(lod, "length_of_day")
                s = np.roll(s, 1)
                s[0] = s[1]
                feats.append(Feature("length_of_day", 3, "linear", s,
                                     "IERS excess length-of-day, ms "
                                     "(lagged 1 d: strictly causal)", lags=lags))
                rec.setdefault("coverage", {})["length_of_day"] = round(cov, 4)
            if verbose:
                print(f"  family 3 {src['key']}: OK ({len(raw)/1e6:.2f} MB, "
                      f"sha256 {rec['sha256'][:12]}..., coverage "
                      f"{rec.get('coverage')})")
        except Exception as exc:                        # noqa: BLE001 -- graceful skip
            rec.update(status="failed", seconds=round(time.time() - t_a, 2),
                       reason=f"{type(exc).__name__}: {exc}")
            if verbose:
                print(f"  family 3 {src['key']}: SKIPPED ({rec['reason']})")
        _log_data(rec)
        log.append(rec)
    return feats, log


# ================================================================ target ===
def build_target(ctx, base, y, window: slice):
    """(counts_t, offset_t) over the mining window: domain-summed observed/expected."""
    counts = y[:, window].sum(axis=0).astype(np.float64)
    offset = base.rate(window).sum(axis=0).astype(np.float64)
    assert counts.size == offset.size
    assert np.isfinite(offset).all() and (offset > 0).all()
    if counts.sum() <= 0:
        raise ValueError("bulk transform produced 0 events in the mining window")
    return counts, offset


def load_event_marks(ctx, data_dir, mag_floor):
    """Per-event (day, day_float, mag, depth) for the domain, aligned to ctx.ev_*.

    `mine` needs depth, which the design cache does not carry (adding it would
    invalidate the cache and force an ETAS refit). It goes back through the SAME
    loader the design was built from -- never a raw CSV -- and then asserts the
    reconstruction is identical to ctx.ev_day / ctx.ev_mag before using it.
    """
    cat, _rep = datasets.load_catalog(data_dir, mag_floor)
    days_f, lat, lon, mag, _t0 = datasets.catalog_arrays(cat)
    from . import grid as gridmod
    day = gridmod.day_index(days_f, ctx.n_days)
    cell = ctx.grid.cell_index(lat, lon)
    keep = cell >= 0
    ev_day = day[keep]
    ev_mag = mag[keep]
    assert ev_day.size == ctx.ev_day.size, (
        f"mark reconstruction has {ev_day.size} events, design has {ctx.ev_day.size}")
    assert np.array_equal(ev_day.astype(np.int64), ctx.ev_day), \
        "mark reconstruction does not align with the cached design (day index)"
    assert np.allclose(ev_mag, ctx.ev_mag), \
        "mark reconstruction does not align with the cached design (magnitude)"
    return {"day": ev_day.astype(np.int64), "day_float": days_f[keep],
            "mag": ev_mag.astype(np.float64),
            "depth": cat["depth"].to_numpy(dtype="float64")[keep]}


# ============================================================== GLM tests ===
def glm_fit(X, y, offset, max_iter=60, tol=1e-10):
    """Newton fit of log lambda = a + X b + log(offset). Returns coefficients + LL."""
    X = np.asarray(X, dtype=np.float64)
    n, k = X.shape
    D = np.column_stack([np.ones(n), X])
    y = np.asarray(y, dtype=np.float64)
    off = np.asarray(offset, dtype=np.float64)
    beta = np.zeros(k + 1)
    beta[0] = math.log(y.sum() / off.sum())
    ll = -np.inf
    converged = False
    for it in range(1, max_iter + 1):
        lam = off * np.exp(D @ beta)
        ll_new = float((y * np.log(lam) - lam).sum())
        g = D.T @ (y - lam)
        H = (D * lam[:, None]).T @ D
        try:
            step = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            break
        beta = beta + step
        if abs(ll_new - ll) < tol * max(1.0, abs(ll_new)):
            ll = ll_new
            converged = True
            break
        ll = ll_new
    lam = off * np.exp(D @ beta)
    ll = float((y * np.log(lam) - lam).sum())
    H = (D * lam[:, None]).T @ D
    cov = np.linalg.inv(H)
    a0 = math.log(y.sum() / off.sum())
    ll0 = float((y * np.log(off * math.exp(a0)) - off * math.exp(a0)).sum())
    n_ev = float(y.sum())
    return {
        "a": float(beta[0]), "beta": beta[1:].tolist(),
        "se": np.sqrt(np.maximum(np.diag(cov)[1:], 0.0)).tolist(),
        "ll_model": ll, "ll_null": ll0, "n_events": n_ev,
        "bits_per_event": (ll - ll0) / (n_ev * LN2),
        "lr_chi2": 2.0 * (ll - ll0), "df": k, "converged": bool(converged),
        "n_iter": it,
    }


def _corr_all_shifts(x, s):
    """c[k] = sum_t x[t] * s[t-k] for every circular shift k, via one FFT pair."""
    n = x.size
    return np.fft.irfft(np.fft.rfft(x) * np.conj(np.fft.rfft(s)), n)


def score_stat_all_shifts(X, y, offset):
    """Rao score statistic for beta=0 under EVERY circular shift of the target.

    Shifting the (counts, offset) PAIR keeps the null intercept exactly invariant
    (both totals are preserved), so S[0] is the observed statistic and S[k] for
    k != 0 is an exact surrogate draw -- no sampling error, no iteration. This is
    why the miner can afford a full-enumeration null instead of 1000 draws.
    """
    X = np.asarray(X, dtype=np.float64)
    n, k = X.shape
    y = np.asarray(y, dtype=np.float64)
    off = np.asarray(offset, dtype=np.float64)
    lam0 = off * (y.sum() / off.sum())
    d = y - lam0
    W = float(lam0.sum())

    U = np.column_stack([_corr_all_shifts(X[:, j], d) for j in range(k)])
    A = np.column_stack([_corr_all_shifts(X[:, j], lam0) for j in range(k)])
    B = np.empty((n, k, k))
    for i in range(k):
        for j in range(i, k):
            b = _corr_all_shifts(X[:, i] * X[:, j], lam0)
            B[:, i, j] = b
            B[:, j, i] = b
    V = B - A[:, :, None] * A[:, None, :] / W
    # ridge only against exact singularity (never reached for these designs)
    V = V + np.eye(k)[None, :, :] * (1e-12 * np.abs(B[:, np.arange(k), np.arange(k)]
                                                    ).max())
    return _quadform(V, U)


def block_bootstrap_idx(n, n_boot, block, rng):
    """Circular moving-block bootstrap indices, (n_boot, n) int32.

    Fixed block length `block`, blocks drawn with replacement from all n circular
    starting positions, concatenated and truncated to length n.
    """
    L = int(max(1, min(int(block), n)))
    nb = int(math.ceil(n / L))
    starts = rng.integers(0, n, size=(int(n_boot), nb))
    idx = (starts[:, :, None] + np.arange(L)[None, None, :]) % n
    return idx.reshape(int(n_boot), nb * L)[:, :n].astype(np.int32)


def score_stat_block_bootstrap(X, y, offset, n_boot, rng, mean_block=90.0,
                               chunk=200):
    """Score statistics under a circular moving-block bootstrap of the target.

    Why this exists (and why it is the PRIMARY null for periodic features): a
    circular time shift of the target does not change the target's power spectrum,
    and the 2-df score test on (sin, cos) of a fixed-frequency feature IS a
    measurement of that spectrum at that frequency. So circular-shift surrogates
    inherit the very signal they are supposed to destroy. Measured on a synthetic
    8%-amplitude 29.53 d modulation over 7716 days: observed chi2 = 152, and the
    MEDIAN circular-shift surrogate = 97. The shift null is not wrong, it is
    powerless -- for periodic features only.

    The block bootstrap resamples the (counts, offset) pair in whole blocks of
    `mean_block` days, so it preserves the target's short-range
    clustering while destroying its alignment with a deterministic cycle. The block
    length must be at least a couple of cycles of the feature under test, or the
    surrogates cannot represent the structure being tested; the caller sets it from
    the feature's own timescale.
    """
    X = np.asarray(X, dtype=np.float64)
    n, _k = X.shape
    y = np.asarray(y, dtype=np.float64)
    off = np.asarray(offset, dtype=np.float64)
    pre = _score_rows_precompute(X)
    out = np.empty(int(n_boot))
    done = 0
    while done < n_boot:
        b = int(min(chunk, n_boot - done))
        idx = block_bootstrap_idx(n, b, mean_block, rng)
        out[done:done + b] = _score_rows(y[idx], off[idx], pre)
        done += b
    return out


def _score_rows_precompute(X):
    """The design-only part of the batched score statistic, hoisted out of the loop."""
    X = np.asarray(X, dtype=np.float64)
    k = X.shape[1]
    pairs = [(i, j) for i in range(k) for j in range(i, k)]
    return {"X": X, "k": k, "pairs": pairs,
            "XX": {(i, j): X[:, i] * X[:, j] for i, j in pairs}}


def _score_rows(ys, os_, pre):
    """Rao score statistic for beta=0, one value per ROW of a (b, n) batch.

    Factored out of `score_stat_block_bootstrap` unchanged so that the per-region
    statistic (§P6-4 Rule 4.2) is computed by literally the same arithmetic as the
    global one -- the regional sum must not be a second, subtly different, test.
    """
    X, k, pairs, XX = pre["X"], pre["k"], pre["pairs"], pre["XX"]
    ys = np.atleast_2d(np.asarray(ys, dtype=np.float64))
    os_ = np.atleast_2d(np.asarray(os_, dtype=np.float64))
    b = ys.shape[0]
    lam0 = os_ * (ys.sum(axis=1) / os_.sum(axis=1))[:, None]
    d = ys - lam0
    U = d @ X
    A = lam0 @ X
    W = lam0.sum(axis=1)
    V = np.empty((b, k, k))
    for i, j in pairs:
        v = lam0 @ XX[(i, j)] - A[:, i] * A[:, j] / W
        V[:, i, j] = v
        V[:, j, i] = v
    return _quadform(V, U)


# ------------------------------- §P6-4 Rule 4.2: the 2R-df phase-incoherent sum -
# The blind spot §K87-0(d)(i) named is that a coherent regional signal with
# REGION-DEPENDENT PHASE is cancelled by the domain sum before any test runs. The
# cure is not more tests, it is a statistic that cannot cancel: a sum of per-region
# 2-df quadratic forms. Each term is |U_r|^2 in the V_r metric, which is
# non-negative and phase-blind, so R regions in antiphase ADD rather than cancel.
# One test per (feature, lag); the sum is the tested statistic and every surrogate
# goes through the same sum.
def score_stat_regions(X, C, O):
    """Per-region score statistics, (R,). C and O are (R, n_days)."""
    return _score_rows(C, O, _score_rows_precompute(X))


def score_stat_regsum(X, C, O):
    """The tested statistic: sum of the per-region score statistics, ~ chi2(R*df)."""
    return float(np.sum(score_stat_regions(X, C, O)))


def score_stat_regsum_all_shifts(X, C, O):
    """The regional sum under EVERY common circular shift of the (C, O) pair.

    A single shift applied to all regions at once is one surrogate world: it
    destroys alignment with the deterministic feature while preserving every
    region's own serial structure AND the cross-region dependence. Index 0 is the
    observed sum. Exhaustive, exact, and it goes through the same sum as the
    observed statistic -- which is the requirement, not a convenience.
    """
    C = np.atleast_2d(np.asarray(C, dtype=np.float64))
    O = np.atleast_2d(np.asarray(O, dtype=np.float64))
    total = None
    for r in range(C.shape[0]):
        s = score_stat_all_shifts(X, C[r], O[r])
        total = s if total is None else total + s
    return total


def score_stat_regsum_block_bootstrap(X, C, O, n_boot, rng, mean_block=90.0,
                                      chunk=200):
    """Block-bootstrap null for the regional SUM, (n_boot,).

    ONE index draw per surrogate, shared by every region. That is deliberate: the
    regions are spatially dependent, and resampling them on independent indices
    would manufacture a narrower null than the data supports. Sharing the index
    keeps the cross-region covariance of the real series in every surrogate world.
    """
    C = np.atleast_2d(np.asarray(C, dtype=np.float64))
    O = np.atleast_2d(np.asarray(O, dtype=np.float64))
    R, n = C.shape
    pre = _score_rows_precompute(X)
    out = np.zeros(int(n_boot))
    done = 0
    while done < n_boot:
        b = int(min(chunk, n_boot - done))
        idx = block_bootstrap_idx(n, b, mean_block, rng)
        acc = np.zeros(b)
        for r in range(R):
            acc += _score_rows(C[r][idx], O[r][idx], pre)
        out[done:done + b] = acc
        done += b
    return out


def _quadform(V, U):
    """U^T V^-1 U per row, falling back to the pseudo-inverse if V is singular.

    V goes singular exactly when a design column has been flattened -- which is what
    the aliasing audit does on purpose when it re-bins a 2-day alternation at 2-day
    resolution. That is a real answer ("the effect does not survive this lattice"),
    not an error, so it must return 0 rather than raise.
    """
    try:
        return np.einsum("bi,bi->b", U, np.linalg.solve(V, U[:, :, None])[:, :, 0])
    except np.linalg.LinAlgError:
        return np.einsum("bi,bi->b", U, np.einsum("bij,bj->bi", np.linalg.pinv(V), U))


def bootstrap_p(S_obs, S_boot):
    ge = int((np.asarray(S_boot) >= float(S_obs)).sum())
    return (1.0 + ge) / (1.0 + np.asarray(S_boot).size)


# =============================================== canonical test-key schema ===
# Every random stream in the miner is addressed by a COMPLETE test key hashed as
# canonical JSON -- not by a raw tuple, so the schema can gain fields without
# silently re-keying the fields that already exist. `region` and `bin_width` are
# present and defaulted NOW, unused today, precisely so that adding per-region or
# per-binning strata later does not change the key schema of anything else.
TEST_KEY_FIELDS = ("master_seed", "feature", "kind", "lag", "rung", "null_type",
                   "region", "bin_width")


def test_key(master_seed, feature, kind, lag=None, rung=None, null_type=None,
             region=None, bin_width=None):
    """The complete, canonical identity of one test (or one rung of one test)."""
    return {
        "master_seed": int(master_seed),
        "feature": str(feature),
        "kind": str(kind),
        "lag": None if lag is None else int(lag),
        "rung": None if rung is None else int(rung),
        "null_type": None if null_type is None else str(null_type),
        "region": None if region is None else str(region),
        "bin_width": None if bin_width is None else float(bin_width),
    }


def test_key_digest(key):
    """sha256 of the key's canonical JSON. Stable across processes and versions."""
    missing = set(TEST_KEY_FIELDS) - set(key)
    if missing:
        raise KeyError(f"test key missing fields {sorted(missing)}")
    blob = json.dumps({k: key[k] for k in TEST_KEY_FIELDS}, sort_keys=True,
                      separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def seed_sequence_for(key):
    """SeedSequence for a test key. Depends on the key and nothing else."""
    d = test_key_digest(key)
    return np.random.SeedSequence(entropy=int(key["master_seed"]),
                                  spawn_key=(int(d[:32], 16),))


def rung_seed_sequence(parent, rung):
    """Child SeedSequence for index `rung`, addressable BY INDEX.

    `SeedSequence.spawn(k)` is stateful -- it hands out children in call order, so
    the stream a rung gets would depend on how many rungs were spawned before it,
    which is exactly the process-boundary dependence we must not have. Constructing
    the child by index reproduces what `spawn` would have produced for that index
    while being addressable: rung 2 is rung 2 whether rung 1 ran in this process,
    another process, or not at all.
    """
    return np.random.SeedSequence(entropy=parent.entropy,
                                  spawn_key=tuple(parent.spawn_key) + (int(rung),))


def surrogate_seed_sequence(parent, index):
    """Child SeedSequence for SURROGATE `index` of one Monte Carlo null.

    Deliberately the SAME construction as `rung_seed_sequence` -- an index-addressable
    child, never `SeedSequence.spawn` -- and deliberately keyed on the SURROGATE
    index rather than on a chunk index. That is the whole design of the period
    scan's decomposition:

      * surrogate i is surrogate i whatever chunk it lands in, whatever process
        drew it, and whether or not surrogates 0..i-1 were drawn at all;
      * therefore the CHUNK SIZE is a pure execution knob -- it changes how the work
        is packed for the pool and nothing else. `--jobs 1` and `--jobs 32` produce
        identical numbers, and so do chunk=10 and chunk=200
        (engine/tests/test_mine_parallel.py::test_period_chunk_size_invariance).

    Note what is NOT done here: the chunk index is NOT added to TEST_KEY_FIELDS.
    Adding a field to the canonical key would re-key every existing digest in the
    schema, which the phase-1 contract forbids. The spawn-key child path is exactly
    the extensibility axis phase 1 built for sub-indices, and it leaves every
    already-declared test digest byte-identical.
    """
    return rung_seed_sequence(parent, index)


# ================================================ adaptive surrogate ladder ===
# OFF by default. When enabled (`--ladder`) its parameters enter the config hash,
# because a stopping rule is a STATISTICAL choice, not an execution detail.
#
# `chunk` is an IMPLEMENTATION DETAIL -- the number of surrogates generated per
# vectorised call -- and is NOT a statistical stage. The rule below is defined on
# single draws and recovers the exact stopping index inside whatever chunk trips
# it, so changing `chunk` changes speed and nothing else. It is nevertheless
# hashed, because a claim that it changes nothing should be checkable.
LADDER_DEFAULTS = {
    "rule": "besag_clifford_1991",
    "h": 25,              # stop once this many surrogates reach the observed stat
    "chunk": 200,         # vectorisation chunk size. NOT a statistical stage.
}


def ladder_n_max(ladder_cfg, key):
    """N_max for one test's declared stratum. FUNCTION OF THE TEST KEY ONLY.

    Fixed before the run. There is deliberately no path by which an observed
    statistic, a p-value, a rank or a survivor list can raise or lower a test's
    budget: rank-based escalation would invalidate the sequential p-value, and no
    amount of care elsewhere would repair it. Strata may key on any field of the
    test key; the default is a single stratum at the preset's surrogate count.
    """
    strata = ladder_cfg.get("n_max_by_kind") or {}
    return int(strata.get(key.get("kind"), ladder_cfg["n_max"]))


def besag_clifford_p(draw_chunk, s_obs, n_max, rung_rng, h=25, chunk=200,
                     keep_draws=False):
    """Sequential Monte Carlo test with an EXACTLY valid p under optional stopping.

    Besag & Clifford (1991), "Sequential Monte Carlo p-values", Biometrika 78(2)
    301-304. Surrogates are drawn until `h` of them reach the observed statistic.
    If the h-th exceedance lands on draw l (counted inclusive), the p-value is

        p = h / l

    and if the cap `n_max` is reached with only c < h exceedances it falls back to
    the ordinary fixed-sample estimate p = (1 + c) / (1 + n_max). Both branches are
    valid under the null -- which is the whole point, since an ordinary Monte Carlo
    p-value computed after peeking would not be.

    CHUNKING IS NOT A STAGE. Draws are generated `chunk` at a time for
    vectorisation. The chunk that trips the rule is generated in full and `l` is
    then recovered from the ORDER of exceedances INSIDE it. Setting l to the chunk
    boundary instead would report a smaller p than the rule permits
    (anti-conservative); recovering the exact l is the rule exactly as published,
    and makes the answer independent of `chunk`.

    `draw_chunk(b, rng)` returns b surrogate statistics. `rung_rng(i)` returns the
    generator for chunk i; it must be index-addressable (see `rung_seed_sequence`)
    so that a given chunk's draws do not depend on which process drew the ones
    before it.

    This is deliberately the ONLY place the stopping rule lives, so an adjudication
    that changes h, the chunk schedule, or the rule itself changes one function.

    `keep_draws` retains the surrogate statistics in the returned dict under the key
    `draws`, for §P6-2 tail extrapolation, and does NOT change a single number this
    function computes. It is off by default because those arrays are per-test
    vectors of up to N_max floats and this dict is checkpointed as JSON: the caller
    that asks for them must POP them before the row is stored.
    """
    s_obs = float(s_obs)
    n_max, h = int(n_max), int(h)
    chunk = int(max(1, chunk))
    if h < 1:
        raise ValueError("besag_clifford_p: h must be >= 1")
    total, c, l_stop, rung = 0, 0, None, 0
    kept = [] if keep_draws else None
    while total < n_max:
        b = int(min(chunk, n_max - total))
        s = np.asarray(draw_chunk(b, rung_rng(rung)), dtype=np.float64).ravel()
        if s.size != b:
            raise ValueError(f"draw_chunk({b}) returned {s.size} statistics")
        if kept is not None:
            kept.append(s.copy())
        hit = s >= s_obs
        n_hit = int(hit.sum())
        if c + n_hit >= h:
            k = int(np.nonzero(np.cumsum(hit) == (h - c))[0][0])
            l_stop = total + k + 1
            c += n_hit
            total += b
            rung += 1
            break
        c += n_hit
        total += b
        rung += 1
    base = {"rule": "besag_clifford_1991", "h": h, "chunk": chunk, "n_max": n_max,
            "n_drawn": int(total), "n_rungs": int(rung), "n_exceed": int(c)}
    if kept is not None:
        base["draws"] = np.concatenate(kept) if kept else np.empty(0)
    if l_stop is not None:
        base.update({"p": h / float(l_stop), "l_stop": int(l_stop),
                     "stopped_early": True,
                     "p_basis": "h/l (Besag-Clifford sequential)"})
    else:
        base.update({"p": (1.0 + c) / (1.0 + n_max), "l_stop": None,
                     "stopped_early": False,
                     "p_basis": "(1+c)/(1+n_max) (cap reached)"})
    return base


def empirical_p(S, n_surr, rng, min_shift=30):
    """p = rank of S[0] among circular-shift surrogates. Returns (p, n_used)."""
    n = S.size
    ok = np.arange(n)
    ok = ok[(ok >= min_shift) & (ok <= n - min_shift)]
    if ok.size == 0:
        raise ValueError("no admissible circular shifts")
    if n_surr < ok.size:
        ok = rng.choice(ok, size=int(n_surr), replace=False)
    ge = int((S[ok] >= S[0]).sum())
    return (1.0 + ge) / (1.0 + ok.size), int(ok.size)


def chi2_sf(x, df):
    """Upper tail of chi2 without importing scipy at module scope."""
    from scipy import stats
    return float(stats.chi2.sf(max(float(x), 0.0), int(df)))


# ============================================================ period scan ===
# ------------------------------------------------------------------------------
# THE CACHED-BASIS LOMB-SCARGLE KERNEL (v2 phase 1c). ONE KERNEL, NO FLAG.
#
# Phase 1b measured the ceiling and named the cause: `scipy.signal.lombscargle` is
# pure numpy and allocates ~8 float64 temporaries of shape (n_time, n_freq) per call
# -- about 49 MB EACH at 7716 x 800 -- so a periodogram is memory-bandwidth bound
# and aggregate throughput saturates a few x above serial however many cores are
# thrown at it (proven in phase 1b with zero engine code in the loop).
#
# Everything scipy computes from (x, freqs) alone is IDENTICAL for the observed
# series and for all ~2 x n_mc surrogates of one scan. Reading scipy's own source
# (scipy/signal/_spectral_py.py, 1.18.0) with uniform weights = 1/N and
# floating_mean=False, the computation is:
#
#     CC0 = sum_j cos^2(w t_j)/N        SS0 = 1 - CC0      CS0 = sum_j cos sin /N
#     tau = 0.5 * atan2(2 CS0, CC0 - SS0)
#     C[k,j] = cos(w_k t_j - tau_k)     S[k,j] = sin(w_k t_j - tau_k)
#     CC  = sum_j C^2/N (clipped >= epsneg)      SS = 1 - CC (clipped >= epsneg)
#     YC  = sum_j (x_j/N) C[k,j]        YS = sum_j (x_j/N) S[k,j]
#     YY  = sum_j (x_j/N) x_j
#     power = 2 (YC^2/CC + YS^2/SS) * 0.5/YY
#
# The ONLY x-dependent quantities are YC, YS and YY. So tau, C, S, CC and SS are
# built ONCE per (t, periods) pair and every subsequent series costs exactly one
# BLAS matrix-vector product against the stacked basis plus one dot product. Per
# periodogram that is one streaming read of the basis (~99 MB at the quick preset)
# instead of ~8 allocate-write-read passes over the same footprint.
#
# WHY THERE IS NO --kernel FLAG, AND WHY THAT IS THE SAFER CHOICE. A per-call
# kernel choice is exactly the structure that could let an observed statistic and
# its surrogates be scored by different code, which would corrupt a rank while every
# individual number stayed defensible. There is therefore ONE implementation of
# `lomb_scargle_power` and no selector of any kind: `period_observed` and
# `period_surrogate_maxima` cannot disagree because there is nothing to disagree
# about. The licence for the change is §P6-5 (a deliberate kernel change) plus the
# validation in engine/tests/test_mine_kernel.py: max relative deviation from
# scipy.signal.lombscargle below 1e-10 over the full period grid on real and
# synthetic series, and the quick-preset p-values unchanged.
#
# WHY THE BASIS IS NOT BATCHED OVER SURROGATES. Evaluating a whole chunk of
# surrogates as one (n_time, n_chunk) matrix product would cut the basis traffic
# further, but BLAS is free to sum a dgemm differently from a dgemv, so the answer
# would acquire a dependence on the CHUNK SIZE at the last ulp -- and chunk size is
# declared a pure execution knob, pinned by
# test_mine_parallel.py::test_period_chunk_size_invariance. One series at a time is
# bit-stable under every chunking, so that is what this does.
# ------------------------------------------------------------------------------

# Per-worker ceiling on ONE cached basis. The quick preset (7716 days x 800 trial
# periods) needs 2 x 7716 x 800 x 8 B = 98.8 MB, so 30 workers hold ~2.9 GB, which
# this machine has. The assertion exists for the grid nobody has run yet: at
# n_periods=3000 (the `full` preset) the basis is 370 MB per worker and 30 workers
# would want 11 GB. That must fail LOUDLY at the first build, naming the fix, rather
# than swap the machine to death halfway through a sweep.
LS_BASIS_MAX_BYTES = 512 * 1024 * 1024
LS_BASIS_CACHE_MAX = 4               # entries kept per process
LS_BASIS_BUILD_BLOCK = 64            # frequency rows built at once (transient only)
LS_BASIS_ANNOUNCE_BYTES = 1 << 20    # print the footprint for bases above 1 MB

_LS_BASIS_CACHE = {}                 # digest -> basis dict. PER PROCESS.
_LS_BASIS_LAST = None                # (t_obj, periods_obj, basis) identity fast path


def ls_basis_footprint_bytes(n_time, n_freq):
    """Bytes ONE cached basis occupies: the stacked (2 n_freq, n_time) float64 array.

    The (2 n_freq,) CC/SS vectors are four orders of magnitude smaller and are not
    counted; this number is the one that decides whether jobs x workers fits in RAM.
    """
    return int(2 * int(n_time) * int(n_freq) * 8)


def ls_basis_digest(t, periods):
    """Content hash of the (time grid, period grid) pair. Stable across processes."""
    t = np.ascontiguousarray(np.asarray(t, dtype=np.float64).ravel())
    p = np.ascontiguousarray(np.atleast_1d(np.asarray(periods, dtype=np.float64)
                                           ).ravel())
    h = hashlib.sha256()
    h.update(np.asarray([t.size, p.size], dtype=np.int64).tobytes())
    h.update(t.tobytes())
    h.update(p.tobytes())
    return h.hexdigest()


def build_ls_basis(t, periods, block=LS_BASIS_BUILD_BLOCK, announce=True):
    """Build the x-independent half of the Lomb-Scargle periodogram, once.

    Returns {"B": (n_time, 2 n_freq) float64 with columns [cos-cols | sin-cols]
    already rotated by tau, "CC", "SS", "n_time", "n_freq", "bytes"}. Built in
    frequency blocks so the TRANSIENT working set is `block` columns, not another
    full copy.

    THE LAYOUT IS SCIPY'S, DELIBERATELY, AND IT COST 1.2 ms PER PERIODOGRAM TO KEEP.
    The obvious layout is (2 n_freq, n_time) with a row-wise matvec, which measured
    3.08 ms against this one's 4.31 ms. It was rejected: it reduces every accumulation
    to a different summation order from `scipy.signal.lombscargle`, and at the
    degenerate Nyquist endpoint P = 2 d on an integer day grid -- where sin(pi j) is
    identically zero, SS is clipped to machine `epsneg` and CS is pure round-off --
    that turned a 1e-16 difference in CS into a 1.2e-8 RELATIVE difference in the
    reported power. Harmless in fact (the point is a grid endpoint, never a local
    maximum, five orders below the peak) but it would have shipped a validation
    number of 1e-8 against a stated bar of 1e-10, i.e. an excuse. Keeping scipy's own
    (n_time, n_freq) orientation and its exact dot-product sequence measures
    9.2e-15 instead, with no excuse attached.
    """
    t = np.ascontiguousarray(np.asarray(t, dtype=np.float64).ravel())
    periods = np.ascontiguousarray(
        np.atleast_1d(np.asarray(periods, dtype=np.float64)).ravel())
    n, nf = int(t.size), int(periods.size)
    if n == 0 or nf == 0:
        raise ValueError("Lomb-Scargle basis needs a non-empty time and period grid")
    nbytes = ls_basis_footprint_bytes(n, nf)
    if nbytes > LS_BASIS_MAX_BYTES:
        raise MemoryError(
            f"Lomb-Scargle cached basis would need {nbytes / 2**20:.0f} MiB per "
            f"WORKER PROCESS ({n} time samples x {nf} trial periods x 2 x 8 B), "
            f"over the {LS_BASIS_MAX_BYTES / 2**20:.0f} MiB per-worker ceiling "
            f"(engine.mine.LS_BASIS_MAX_BYTES). Every worker holds its own copy, so "
            f"the real bill is that number times --jobs. Fix by running fewer "
            f"workers, cutting n_periods, or raising the ceiling DELIBERATELY after "
            f"checking total RAM against --jobs.")

    freqs = (2.0 * np.pi / periods).reshape(1, -1)
    tcol = t.reshape(-1, 1)
    weights = np.full((n, 1), 1.0 / n)         # scipy's normalised uniform weights
    B = np.empty((n, 2 * nf), dtype=np.float64)
    CC = np.empty(nf, dtype=np.float64)
    SS = np.empty(nf, dtype=np.float64)
    epsneg = float(np.finfo(np.float64).epsneg)
    for a in range(0, nf, int(max(1, block))):
        b = min(a + int(max(1, block)), nf)
        # every line below is scipy/signal/_spectral_py.py's lombscargle body for
        # uniform weights and floating_mean=False, in its order, on a column block
        freqst = freqs[:, a:b] * tcol
        coswt, sinwt = np.cos(freqst), np.sin(freqst)
        cc0 = np.dot(weights.T, coswt * coswt)
        ss0 = 1.0 - cc0
        cs0 = np.dot(weights.T, coswt * sinwt)
        tau = 0.5 * np.arctan2(2.0 * cs0, cc0 - ss0)
        freqst_tau = freqst - tau
        ct = np.cos(freqst_tau, out=coswt)
        st = np.sin(freqst_tau, out=sinwt)
        cc = np.dot(weights.T, ct * ct)
        ss = 1.0 - cc
        # scipy clips BOTH after forming SS = 1 - CC; replicate exactly, in order
        cc[cc < epsneg] = epsneg
        ss[ss < epsneg] = epsneg
        B[:, a:b] = ct
        B[:, nf + a:nf + b] = st
        CC[a:b] = np.squeeze(cc, axis=0)
        SS[a:b] = np.squeeze(ss, axis=0)
    basis = {"B": B, "CC": CC, "SS": SS, "n_time": n, "n_freq": nf,
             "bytes": nbytes, "digest": None}
    if announce and nbytes >= LS_BASIS_ANNOUNCE_BYTES:
        print(f"  [pid {os.getpid()}] Lomb-Scargle basis cached: {n} x {nf} "
              f"= {nbytes / 2**20:.1f} MiB per worker process "
              f"(ceiling {LS_BASIS_MAX_BYTES / 2**20:.0f} MiB)", flush=True)
    return basis


def ls_basis(t, periods):
    """The cached basis for this (t, periods) pair, built at most once per process.

    Two lookups, cheapest first: an IDENTITY check against the last basis used (which
    is the hit in the inner surrogate loop, where the same two arrays are passed
    hundreds of times), then a content digest. The identity path keeps references to
    the arrays it matched, so an id() cannot be recycled underneath it.
    """
    global _LS_BASIS_LAST
    last = _LS_BASIS_LAST
    if last is not None and last[0] is t and last[1] is periods:
        return last[2]
    d = ls_basis_digest(t, periods)
    basis = _LS_BASIS_CACHE.get(d)
    if basis is None:
        basis = build_ls_basis(t, periods)
        basis["digest"] = d
        if len(_LS_BASIS_CACHE) >= LS_BASIS_CACHE_MAX:
            _LS_BASIS_CACHE.pop(next(iter(_LS_BASIS_CACHE)))
        _LS_BASIS_CACHE[d] = basis
    _LS_BASIS_LAST = (t, periods, basis)
    return basis


def lomb_scargle_power(t, x, periods):
    """Normalised Lomb-Scargle power at the given trial periods.

    THE one kernel. Identical code for the observed statistic and for every
    surrogate of every null, by construction: there is no selector to get wrong.
    Numerically equivalent to `scipy.signal.lombscargle(t, x - mean(x), 2 pi/periods,
    normalize=True)` to better than 1e-10 relative (validated, not asserted).
    """
    x = np.asarray(x, dtype=np.float64)
    x = x - x.mean()
    bas = ls_basis(t, periods)
    if x.size != bas["n_time"]:
        raise ValueError(f"series has {x.size} samples, time grid has "
                         f"{bas['n_time']}")
    nf = bas["n_freq"]
    y = x.reshape(-1, 1)
    wy = y * (1.0 / bas["n_time"])             # scipy's weights_y, uniform weights
    z = np.dot(wy.T, bas["B"])                 # [YC | YS], one streaming pass
    yc, ys = z[0, :nf], z[0, nf:]
    yy = float(np.dot(wy.T, y)[0, 0])
    pgram = 2.0 * (yc * yc / bas["CC"] + ys * ys / bas["SS"])
    return pgram * (0.5 / yy)


def lomb_scargle_power_scipy(t, x, periods):
    """The pre-1c kernel, kept ONLY as the validation reference. Never called in a run.

    This is deliberately not reachable from any scan: it exists so
    engine/tests/test_mine_kernel.py can measure the deviation of the shipped kernel
    against the library it replaced, and so a future reader can re-measure it.
    """
    from scipy import signal
    x = np.asarray(x, dtype=np.float64)
    x = x - x.mean()
    w = np.atleast_1d(2.0 * np.pi / np.asarray(periods, dtype=np.float64))
    return np.atleast_1d(
        signal.lombscargle(np.asarray(t, dtype=np.float64), x, w, normalize=True))


# --------------------------------------------------------------------------
# THE PERIOD SCAN, DECOMPOSED (v2 phase 1b)
#
# The scan is 2 x n_mc INDEPENDENT full Lomb-Scargle periodograms plus one for the
# observed series. Before phase 1b it was a single indivisible task and, at the
# quick preset, 90.7% of the whole sweep -- an Amdahl ceiling of 1.10x however many
# workers were thrown at it (measured, engine/bench_v2.py, phase-1 report).
#
# It is now four separable pieces:
#   period_observed          the observed periodogram + peak picking (deterministic)
#   period_surrogate_maxima  the max power of surrogates [start, stop) of one null
#   period_assemble          the p-values, given ALL the maxima (counted invariant)
#   period_scan              the three of them wired up in one process
#
# THE SEEDING RULE, which is the load-bearing part. Surrogate i of null `nt` draws
# from `surrogate_seed_sequence(parent[nt], i)` -- keyed on the SURROGATE INDEX, not
# on the chunk. So the chunk boundaries are not a statistical stage: any chunking of
# [0, n_mc) reproduces the same 2 x n_mc surrogates in the same order, and the pool's
# scheduling cannot reach the answer. The old code drew all n_mc AR(1) series from
# one (n_mc, n) normal draw off a shared threaded generator, which was reproducible
# only for a single process running the tasks in one fixed order. That stream is
# DELIBERATELY superseded (§P6-5): determinism independent of scheduling is worth
# more than bit-identity with the pre-decomposition build.
# --------------------------------------------------------------------------
PERIOD_NULLS = ("ar1", "permutation")

# The max-power MC cap. The 95th percentile of a max-power distribution converges
# fast, and each draw costs a full periodogram. The cap sets the period scan's own
# p floor at 1/(cap+1), which the report states.
PERIOD_N_MC_CAP = 500

# Surrogates per subtask. PURE EXECUTION DETAIL, like --jobs: provably result-
# invariant by the per-surrogate seeding above, pinned by a test, and therefore
# deliberately NOT in the config hash. Small enough that the largest period subtask
# is a fraction of the sweep (the whole point of phase 1b), large enough that
# per-task pickling/dispatch stays in the noise.
PERIOD_SURROGATE_CHUNK = 10


def period_n_mc(n_surr):
    """Surrogates per null for the period scan. Capped -- see PERIOD_N_MC_CAP."""
    return int(min(int(n_surr), PERIOD_N_MC_CAP))


def period_chunks(n_mc, chunk=None):
    """[(start, stop), ...] partitioning [0, n_mc). Execution only, never statistics."""
    c = int(chunk if chunk else PERIOD_SURROGATE_CHUNK)
    if c < 1:
        raise ValueError("period chunk size must be >= 1")
    n_mc = int(n_mc)
    return [(s, min(s + c, n_mc)) for s in range(0, n_mc, c)]


def period_seed_seqs(master_seed):
    """The parent SeedSequence of each period-scan null, from the canonical key.

    Two keys where phase 1 declared one (`null_type="ar1_and_permutation"`): the two
    nulls are now separate task families and must not share a stream. Every OTHER
    key digest in the schema is untouched -- no field was added to TEST_KEY_FIELDS.
    """
    return {nt: seed_sequence_for(test_key(int(master_seed), "period_scan", "period",
                                           null_type=nt))
            for nt in PERIOD_NULLS}


def ar1_params(x):
    """(phi, sd) of the AR(1) matched to x's own lag-1 autocorrelation."""
    x = np.asarray(x, dtype=np.float64)
    x = x - x.mean()
    phi = float(np.dot(x[1:], x[:-1]) / max(np.dot(x[:-1], x[:-1]), 1e-30))
    phi = float(np.clip(phi, -0.98, 0.98))
    sd = float(x.std() * math.sqrt(max(1.0 - phi * phi, 1e-6)))
    return phi, sd


def ar1_draw(n, phi, sd, rng):
    """ONE AR(1) realisation of length n, stationary at the first sample."""
    n = int(n)
    e = rng.normal(0.0, sd, size=n)
    out = np.empty(n)
    out[0] = e[0] / math.sqrt(max(1.0 - phi * phi, 1e-6))
    for i in range(1, n):
        out[i] = phi * out[i - 1] + e[i]
    return out


def period_surrogate_maxima(t, resid, periods, null_type, start, stop, parent_seq):
    """Max LS power of surrogates [start, stop) of one null. The parallel unit.

    Pure: no state, no shared generator, no dependence on which other surrogates
    exist. Surrogate i's stream is `surrogate_seed_sequence(parent_seq, i)`.
    """
    if null_type not in PERIOD_NULLS:
        raise ValueError(f"unknown period null {null_type!r}; expected {PERIOD_NULLS}")
    resid = np.asarray(resid, dtype=np.float64)
    start, stop = int(start), int(stop)
    if stop < start:
        raise ValueError(f"period chunk [{start}, {stop}) is empty-backwards")
    phi, sd = ar1_params(resid) if null_type == "ar1" else (0.0, 0.0)
    out = np.empty(stop - start)
    for j, i in enumerate(range(start, stop)):
        g = np.random.default_rng(surrogate_seed_sequence(parent_seq, i))
        s = (ar1_draw(resid.size, phi, sd, g) if null_type == "ar1"
             else g.permutation(resid))
        out[j] = float(lomb_scargle_power(t, s, periods).max())
    return out


def period_observed(t, resid, periods, n_peaks=8):
    """The observed periodogram and its peak indices. Deterministic, no randomness."""
    power = lomb_scargle_power(t, resid, periods)
    # local maxima, strongest first, with a 5% fractional-period exclusion zone
    is_pk = np.r_[False, (power[1:-1] > power[:-2]) & (power[1:-1] >= power[2:]), False]
    cand = np.argsort(power)[::-1]
    peaks = []
    for i in cand:
        if not is_pk[i]:
            continue
        if any(abs(math.log(periods[i] / periods[j])) < 0.05 for j in peaks):
            continue
        peaks.append(int(i))
        if len(peaks) >= n_peaks:
            break
    return power, peaks


def period_assemble(power, periods, peaks, maxima, phi, n_mc, gpd=None,
                    gpd_rng=None):
    """p-values from the collected maxima. COUNTS THE SURROGATES BEFORE USING THEM.

    `maxima` is {null_type: array}. A silently short null -- one dropped chunk, one
    cancelled future -- would shift EVERY p in the scan by a fraction of a rank, so
    the count is an assertion, not an assumption.
    """
    n_mc = int(n_mc)
    got = {}
    for nt in PERIOD_NULLS:
        if nt not in maxima:
            raise RuntimeError(
                f"period scan: null {nt!r} produced NO surrogate maxima at all "
                f"(expected {n_mc}). Refusing to compute a p-value.")
        a = np.asarray(maxima[nt], dtype=np.float64).ravel()
        if a.size != n_mc:
            raise RuntimeError(
                f"period scan: null {nt!r} collected {a.size} surrogate maxima, "
                f"declared {n_mc}. A dropped or duplicated chunk shifts every "
                f"p-value in this scan. Refusing to compute a p-value.")
        if not np.isfinite(a).all():
            raise RuntimeError(
                f"period scan: null {nt!r} returned non-finite surrogate maxima")
        got[nt] = a
    max_red, perm = got["ar1"], got["permutation"]
    power = np.asarray(power, dtype=np.float64)
    periods = np.asarray(periods, dtype=np.float64)
    enabled = bool(gpd and gpd.get("enabled"))
    out = []
    for i in peaks:
        i = int(i)
        p_red = (1.0 + int((max_red >= power[i]).sum())) / (1.0 + n_mc)
        p_perm = (1.0 + int((perm >= power[i]).sum())) / (1.0 + n_mc)
        row = {
            "period_days": float(periods[i]), "power": float(power[i]),
            "p_red_noise": p_red, "p_permutation": p_perm,
            "p_raw": max(p_red, p_perm),          # the CONSERVATIVE one
            "n_surrogates": n_mc, "ar1_phi": float(phi),
        }
        # §P6-2: BOTH period nulls are independently generated (AR(1) draws and
        # permutations), so both are permitted substrates. Each is priced on its
        # own and the CONSERVATIVE (larger) of the two binds -- the same rule the
        # un-extrapolated p_raw already follows, so extrapolation cannot make the
        # answer less conservative than the pair of nulls jointly allows.
        pr_red = gpd_tail.price_row(p_red, n_mc, power[i],
                                    max_red if enabled else None, "ar1", gpd_rng,
                                    enabled=enabled)
        pr_perm = gpd_tail.price_row(p_perm, n_mc, power[i],
                                     perm if enabled else None, "permutation",
                                     gpd_rng, enabled=enabled)
        binding = pr_red if pr_red["p_bh"] >= pr_perm["p_bh"] else pr_perm
        row.update(binding)
        row["gpd_red_noise"] = pr_red["gpd"]
        row["gpd_permutation"] = pr_perm["gpd"]
        row["p_method_binding_null"] = ("ar1" if binding is pr_red
                                        else "permutation")
        out.append(row)
    return out, {"ar1_phi": float(phi), "n_mc": n_mc,
                 "max_power_red_p95": float(np.quantile(max_red, 0.95)),
                 "max_power_perm_p95": float(np.quantile(perm, 0.95))}


def period_scan(t, resid, periods, n_surr, seeds, n_peaks=8, verbose=True,
                chunk=None):
    """Lomb-Scargle scan + max-power nulls (AR(1) red noise AND permutation white).

    `seeds` is either a master seed (int) or the dict `period_seed_seqs` returns.
    It is NOT a Generator: a threaded generator would make the answer depend on how
    the surrogates were split across chunks and processes, which is exactly what
    phase 1b removes. This single-process entry point exists for tests, small
    callers and as the UNCHUNKED REFERENCE the chunked session path is checked
    against; the session driver runs the same pieces as separate tasks.
    """
    if not isinstance(seeds, dict):
        seeds = period_seed_seqs(seeds)
    resid = np.asarray(resid, dtype=np.float64)
    n_mc = period_n_mc(n_surr)
    power, peaks = period_observed(t, resid, periods, n_peaks=n_peaks)
    phi, _sd = ar1_params(resid)
    t0 = time.time()
    maxima = {}
    for nt in PERIOD_NULLS:
        parts = [period_surrogate_maxima(t, resid, periods, nt, a, b, seeds[nt])
                 for a, b in period_chunks(n_mc, chunk)]
        maxima[nt] = np.concatenate(parts) if parts else np.empty(0)
    if verbose:
        print(f"  period scan: {len(periods)} trial periods, {n_mc} AR(1) + {n_mc} "
              f"permutation surrogates (AR(1) phi = {phi:.3f}) in "
              f"{time.time()-t0:.1f}s")
    return period_assemble(power, periods, peaks, maxima, phi, n_mc)


def fold_lr(counts, lam0, days, period, n_bins=8):
    """Epoch-folding concentration: LR of a per-phase-bin rate vs a flat rate."""
    ph = np.mod(np.asarray(days, dtype=np.float64), float(period)) / float(period)
    b = np.minimum((ph * n_bins).astype(int), n_bins - 1)
    Y = np.bincount(b, weights=counts, minlength=n_bins)
    L = np.bincount(b, weights=lam0, minlength=n_bins)
    L = L * (Y.sum() / max(L.sum(), 1e-30))
    m = (Y > 0) & (L > 0)
    return float(2.0 * ((Y[m] * np.log(Y[m] / L[m])).sum() - (Y.sum() - L.sum())))


def harmonic_ladder(counts, lam0, days, period, n_bins=8, max_extend=3,
                    max_period=None):
    """Score {P/3, P/2, P, 2P, 3P}, extending at the winning edge while it improves.

    Rungs longer than a third of the record are not scored at all: with fewer than
    three observed cycles an epoch fold is not measuring a period, and a ladder that
    "wins" at 5827 d on 7716 days of data is reporting the record length.
    """
    hi_cap = float(max_period if max_period else max(len(counts) / 3.0, 4.0))
    rungs = {}
    mults = [1 / 3.0, 1 / 2.0, 1.0, 2.0, 3.0]
    for m in mults:
        p = period * m
        if 2.0 <= p <= hi_cap:
            rungs[round(p, 4)] = {"multiplier": m, "fold_lr": fold_lr(counts, lam0,
                                                                     days, p, n_bins)}
    for _ in range(max_extend):
        ks = sorted(rungs, key=lambda k: rungs[k]["fold_lr"], reverse=True)
        best = ks[0]
        lo, hi = min(rungs), max(rungs)
        if best == hi:
            p = hi * 2.0
        elif best == lo:
            p = lo / 2.0
        else:
            break
        if not (2.0 <= p <= hi_cap) or round(p, 4) in rungs:
            break
        rungs[round(p, 4)] = {"multiplier": p / period,
                              "fold_lr": fold_lr(counts, lam0, days, p, n_bins)}
    if not rungs:
        return {"rungs": [], "winning_period_days": float(period),
                "winning_multiplier": 1.0, "n_bins": n_bins,
                "verdict": f"NOT LADDERED: P exceeds record/3 = {hi_cap:.0f} d",
                "max_period_days": hi_cap}
    best = max(rungs, key=lambda k: rungs[k]["fold_lr"])
    return {
        "max_period_days": hi_cap,
        "rungs": [{"period_days": k, **rungs[k]} for k in sorted(rungs)],
        "winning_period_days": float(best),
        "winning_multiplier": float(rungs[best]["multiplier"]),
        "verdict": ("rung P (as scanned)" if abs(rungs[best]["multiplier"] - 1) < 1e-9
                    else f"rung {rungs[best]['multiplier']:.4g}xP WINS over P"),
        "n_bins": n_bins,
    }


# =============================================================== mark tests ===
def _ranks(a):
    a = np.asarray(a, dtype=np.float64)
    order = np.argsort(a, kind="stable")
    r = np.empty(a.size, dtype=np.float64)
    r[order] = np.arange(1, a.size + 1, dtype=np.float64)
    # average ties
    _, inv, cnt = np.unique(a, return_inverse=True, return_counts=True)
    if (cnt > 1).any():
        sums = np.bincount(inv, weights=r)
        r = (sums / cnt)[inv]
    return r


def _pearson_all_shifts(x, s):
    """Pearson r of x against every circular shift of s (both are 1-D)."""
    n = x.size
    xc = x - x.mean()
    sc = s - s.mean()
    num = _corr_all_shifts(xc, sc)
    den = math.sqrt(float(np.dot(xc, xc)) * float(np.dot(sc, sc)))
    return num / max(den, 1e-30)


def _mark_stat_matrix(feature_at_event, R, kind):
    """Association statistic between a fixed feature and each row of rank matrix R."""
    R = np.atleast_2d(np.asarray(R, dtype=np.float64))
    Rc = R - R.mean(axis=1, keepdims=True)
    Rn = np.sqrt((Rc * Rc).sum(axis=1))

    def corr(v):
        vc = v - v.mean()
        return (Rc @ vc) / np.maximum(Rn * math.sqrt(float(np.dot(vc, vc))), 1e-30)

    if kind == "phase":
        th = np.asarray(feature_at_event, dtype=np.float64)
        c, s = np.cos(th), np.sin(th)
        rcs = float(np.corrcoef(c, s)[0, 1])
        r_xc, r_xs = corr(c), corr(s)
        v = (r_xc ** 2 + r_xs ** 2 - 2 * r_xc * r_xs * rcs) / max(1 - rcs ** 2, 1e-12)
        return np.sqrt(np.maximum(v, 0.0)), None       # circular-linear correlation
    rho = corr(_ranks(feature_at_event))
    return np.abs(rho), rho


def mark_test(feature_at_event, mark, kind, n_surr, rng, min_shift=100,
              block_events=500, ladder=None, rung_rng=None, gpd=None,
              gpd_rng=None):
    """Spearman (linear) or circular-linear (phase) feature-vs-mark association.

    Two nulls, both resampling the MARK series along the time-ordered event
    sequence (which destroys the feature-mark pairing while preserving both
    marginal distributions and the autocorrelation of the marks):

      * circular shift -- all admissible shifts evaluated exactly, because ranks
        are shift-equivariant and the rank correlation is a cross-correlation;
      * circular moving-block bootstrap in blocks of `block_events` events, which
        is the null that still has power when the feature is a deterministic cycle
        (see `score_stat_block_bootstrap` for why a shift alone does not).

    `p_raw` is the LARGER (more conservative) of the two.
    """
    rm = _ranks(mark)
    if kind == "phase":
        th = np.asarray(feature_at_event, dtype=np.float64)
        c, s = np.cos(th), np.sin(th)
        rcs = float(np.corrcoef(c, s)[0, 1])
        r_xc = _pearson_all_shifts(c, rm)
        r_xs = _pearson_all_shifts(s, rm)
        stat = (r_xc ** 2 + r_xs ** 2 - 2 * r_xc * r_xs * rcs) / max(1 - rcs ** 2, 1e-12)
        stat = np.sqrt(np.maximum(stat, 0.0))
        effect, df = float(stat[0]), 2
    else:
        rf = _ranks(feature_at_event)
        rho = _pearson_all_shifts(rf, rm)
        stat = np.abs(rho)
        effect, df = float(rho[0]), 1
    p_shift, n_used = empirical_p(stat, n_surr, rng, min_shift=min_shift)

    n_ev = rm.size

    def _draw(b, g):
        idx = block_bootstrap_idx(n_ev, int(b), block_events, g)
        v, _ = _mark_stat_matrix(feature_at_event, rm[idx], kind)
        return v

    lad, boot = None, None
    want_draws = bool(gpd and gpd.get("enabled"))
    if ladder:
        # `ladder` is a dict of stopping-rule parameters; see besag_clifford_p.
        # NOTE: the ladder touches ONLY this Monte Carlo block-bootstrap null.
        # `p_circular_shift` above comes from empirical_p, an EXHAUSTIVE enumeration
        # of admissible shifts whose floor is a property of the record length, not
        # of a surrogate budget; sequential stopping is meaningless there and is
        # deliberately not applied.
        if rung_rng is None:
            raise ValueError("mark_test: ladder requires rung_rng")
        lad = besag_clifford_p(_draw, stat[0],
                               n_max=int(ladder.get("n_max", n_surr)),
                               rung_rng=rung_rng,
                               h=int(ladder.get("h", 25)),
                               chunk=int(ladder.get("chunk", 200)),
                               keep_draws=want_draws)
        # POPPED, never stored: `ladder_stop` is checkpointed as JSON.
        boot = lad.pop("draws", None)
        p_boot = lad["p"]
        mc_n_max = int(lad["n_max"])
    else:
        boot = np.empty(int(n_surr))
        done, chunk = 0, 100
        while done < n_surr:
            b = int(min(chunk, n_surr - done))
            boot[done:done + b] = _draw(b, rng)
            done += b
        p_boot = bootstrap_p(stat[0], boot)
        mc_n_max = int(n_surr)
    out = {"effect": effect, "statistic": float(stat[0]),
           "p_circular_shift": p_shift, "p_block_bootstrap": p_boot,
           "p_raw": max(p_shift, p_boot), "n_surrogates": min(n_used, int(n_surr)),
           "df": df, "null": "max(circular-shift, block-bootstrap)",
           "test": "circular-linear" if kind == "phase" else "spearman"}
    if lad is not None:
        out["ladder_stop"] = lad
    # §P6-2(5): the p_method label goes on EVERY row, whether or not extrapolation
    # is enabled. The block bootstrap is the only permitted substrate here -- the
    # circular-shift enumeration is forbidden by §P6-2(4) and is passed as
    # `other_p`, which the conservative max() keeps binding where it is larger.
    out.update(gpd_tail.price_row(
        p_boot, mc_n_max, stat[0], boot if want_draws else None,
        "block_bootstrap", gpd_rng, other_p=p_shift,
        other_floor=1.0 / (n_used + 1.0), enabled=want_draws))
    return out


# ============================================================ multiplicity ===
def benjamini_hochberg(pvals, q=FDR_Q):
    """BH step-up. Returns (adjusted q-values, boolean pass mask) in input order."""
    p = np.asarray(pvals, dtype=np.float64)
    n = p.size
    if n == 0:
        return np.array([]), np.array([], dtype=bool)
    order = np.argsort(p)
    ps = p[order]
    adj = ps * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    out = np.empty(n)
    out[order] = adj
    return out, out <= q


# ========================================================== aliasing audit ===
def rebin(counts, offset, X, days, width):
    """Aggregate the target over `width`-day bins; average the design within a bin."""
    n = counts.size
    nb = n // width
    cut = nb * width
    c = counts[:cut].reshape(nb, width).sum(axis=1)
    o = offset[:cut].reshape(nb, width).sum(axis=1)
    x = X[:cut].reshape(nb, width, X.shape[1]).mean(axis=1)
    d = days[:cut].reshape(nb, width).mean(axis=1)
    return c, o, x, d


def aliasing_audit_glm(feature, lag, window, counts, offset, n_surr, rng):
    """Re-test a GLM survivor at 1 d (reference), 2 d and 7 d binning."""
    days = np.arange(counts.size, dtype=np.float64)
    X1 = feature.design(window, lag)
    out = {}
    for w in (1, 2, 7):
        if w == 1:
            c, o, x = counts, offset, X1
        else:
            c, o, x, _d = rebin(counts, offset, X1, days, w)
        if (x.std(axis=0) < 1e-9 * max(np.abs(x).max(), 1e-12)).any():
            # the re-binning flattened a design column: the pattern IS the lattice
            out[f"{w}d"] = {"chi2": 0.0, "z_equivalent": 0.0, "p_raw": 1.0,
                            "n_bins": int(c.size),
                            "note": "design column is constant at this binning -- "
                                    "the pattern lives entirely on the 1-day lattice"}
            continue
        S = score_stat_all_shifts(x, c, o)
        p_shift, n_used = empirical_p(S, n_surr, rng, min_shift=max(3, 30 // w))
        Sb = score_stat_block_bootstrap(x, c, o, min(n_surr, 500), rng,
                                        mean_block=max(2, feature.block_days / w))
        p_boot = bootstrap_p(S[0], Sb)
        p = p_boot if feature.periodic else max(p_shift, p_boot)
        out[f"{w}d"] = {"chi2": float(S[0]), "z_equivalent": float(math.sqrt(max(S[0], 0))),
                        "p_raw": p, "p_circular_shift": p_shift,
                        "p_block_bootstrap": p_boot,
                        "n_surrogates": min(n_used, n_surr), "n_bins": int(c.size)}
    z1 = out["1d"]["z_equivalent"]
    ratios = {k: (v["z_equivalent"] / z1 if z1 > 0 else 0.0) for k, v in out.items()}
    suspect = any(r < 0.5 for k, r in ratios.items() if k != "1d")
    out["z_ratio_vs_1d"] = {k: round(v, 3) for k, v in ratios.items()}
    out["verdict"] = "LATTICE-SUSPECT" if suspect else "STABLE"
    out["rule"] = ("LATTICE-SUSPECT if any alternative binning drops the "
                   "surrogate-calibrated effect below 50% of its 1-day value")
    return out


def aliasing_audit_period(period, counts, offset, ev_times, n_surr, rng, day0):
    """Re-test a period claim at 2 d / 7 d binning AND on unbinned event times."""
    days = np.arange(counts.size, dtype=np.float64)
    lam0 = offset * (counts.sum() / offset.sum())
    out = {"fold_lr": {}}
    for w in (1, 2, 7):
        if w == 1:
            c, l0, d = counts, lam0, days
        else:
            X = np.zeros((counts.size, 1))
            c, o, _x, d = rebin(counts, offset, X, days, w)
            l0 = o * (c.sum() / o.sum())
        out["fold_lr"][f"{w}d"] = fold_lr(c, l0, d, period)

    # ---- unbinned: Rayleigh statistic on real event times, ETAS-simulated null
    t = np.asarray(ev_times, dtype=np.float64) - day0
    ph = 2 * np.pi * np.mod(t, period) / period
    R = abs(complex(np.cos(ph).sum(), np.sin(ph).sum())) / math.sqrt(max(t.size, 1))
    n_mc = int(min(n_surr, 1000))
    sim = np.empty(n_mc)
    for i in range(n_mc):
        k = rng.poisson(lam0)
        idx = np.repeat(np.arange(k.size), k)
        tt = idx + rng.random(idx.size)
        q = 2 * np.pi * np.mod(tt, period) / period
        sim[i] = abs(complex(np.cos(q).sum(), np.sin(q).sum())) / math.sqrt(
            max(tt.size, 1))
    p_unbinned = (1.0 + int((sim >= R).sum())) / (1.0 + n_mc)
    z1 = out["fold_lr"]["1d"]
    ratios = {k: (v / z1 if z1 > 0 else 0.0) for k, v in out["fold_lr"].items()}
    suspect = any(r < 0.5 for k, r in ratios.items() if k != "1d") or p_unbinned > 0.05
    out.update({
        "fold_lr_ratio_vs_1d": {k: round(v, 3) for k, v in ratios.items()},
        "unbinned_rayleigh": {"R": float(R), "p": p_unbinned, "n_surrogates": n_mc,
                              "n_events": int(t.size)},
        "verdict": "LATTICE-SUSPECT" if suspect else "STABLE",
        "rule": ("LATTICE-SUSPECT if a 2 d or 7 d binning halves the folded "
                 "concentration, or the unbinned event-time Rayleigh test fails "
                 "at p > 0.05 against an ETAS-simulated null"),
    })
    return out
