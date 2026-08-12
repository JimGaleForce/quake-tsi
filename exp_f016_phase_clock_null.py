"""F-016 - null calibration of anchor-based tidal PHASE CLOCKS. FROZEN PROTOCOL.

WHAT THIS MEASURES, IN ONE SENTENCE. Anchor-based phase conventions force each
segment of a tidal cycle to occupy a fixed span of phase; a mixed diurnal/semidiurnal
tide does NOT spend equal time in those segments; therefore uniform-in-time events do
not emerge uniform in phase. The induced first-harmonic amplitude is a property of the
INSTRUMENT, not of the Earth, and this script measures it for two conventions:

  quarter_anchored  "quarter-anchored (Tanaka-style) phase assignment": anchors the
                    trough / ascending zero / peak / descending zero of the stress
                    series at 0 / 90 / 180 / 270 deg and interpolates linearly in time
                    between anchors, then shifts by 180 deg so peak stress = 0 and
                    wraps to [-180, +180]. Four anchors per cycle.
  two_anchor        THIS repository's own convention, imported unmodified from
                    `coso_positive_control.phase_series`: trough -> peak -> trough,
                    -180 at the preceding trough, 0 at the peak, +180 at the following
                    trough, piecewise linear in time. Two anchors per cycle.

THE TWO QUANTITIES, AND WHY THEY ARE REPORTED SEPARATELY (F-016 gap #2, binding).
Conflating them is the exact error this freeze exists to prevent:

  instrument_response  Deterministic. Every sample of a UNIFORM TIME GRID is pushed
                       through the convention and the first-harmonic amplitude of the
                       resulting phase distribution is computed. No random numbers are
                       involved at all, so there is NO sampling noise in it. This is
                       the finding.
  sampling_floor       Monte Carlo. N uniform-random event times are drawn and the
                       same amplitude is computed. A finite N carries a Rayleigh floor
                       of E|R| = sqrt(pi/N) even against a perfectly flat clock. This
                       is BOOKKEEPING - it tells a re-runner how much of a small
                       measured amplitude is their own N, and nothing else. It is not
                       an instrument property and may never be quoted as one.

Every number in results_phase_clock_null.json carries a "quantity" field with one of
those two values. Nothing is reported without it.

DISTRIBUTIONS, NOT HEADLINES (F-016 requirement). The response depends on the phasing
of the constituents relative to each other, which is a free parameter of any synthetic
mix and an accident of epoch in any real one. Every synthetic mix is therefore run
over >= N_REALIZATIONS independent constituent-phase draws and reported as
min / median / max (plus quartiles), never as a single number. Every real series is
run over >= N_REALIZATIONS disjoint record segments, which is the same axis measured
the only way a fixed record allows.

DETERMINISM. All seeds are fixed constants below. There is no time-based, hostname-
based or path-based nondeterminism anywhere in this file; re-running it reproduces
every digit.

PROVENANCE / SCOPE, stated so it cannot be over-read.
  * The quarter-anchored convention is a shared convention of the tidal-triggering
    literature. This script characterises the CONVENTION. It is not a criticism of any
    paper: a surrogate or occupancy correction cancels this response identically, and
    the careful work in this literature carries one.
  * `two_anchor` is OUR convention and OUR defect. It is included precisely so that
    this measurement is not something we point only at other people's code.

USAGE
  python -u exp_f016_phase_clock_null.py            # full battery (~20-40 min)
  python -u exp_f016_phase_clock_null.py --smoke    # reduced, end-to-end (~1-2 min)
Output: results_phase_clock_null.json (next to this file).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results_phase_clock_null.json")

# ============================================================ frozen protocol ===
# Changing anything in this block changes PROTOCOL_HASH and therefore the identity of
# the protocol. That is the point of freezing it.
PROTOCOL = {
    "id": "F-016-phase-clock-null",
    "version": "1.0",
    "frozen": "2026-08-11",
    "conventions": ["quarter_anchored", "two_anchor"],
    "quantities": ["instrument_response", "sampling_floor"],
    "seed_constituent_phases": 20260811,
    "seed_monte_carlo": 20260812,
    "seed_segments": 20260813,
    "n_realizations": 24,            # >= 20 constituent-phase draws per mix (F-016)
    "n_segments_real": 24,           # >= 20 disjoint segments per real series
    "years": 10.0,                   # synthetic record length
    "dt_hours": 1.0,                 # synthetic sampling
    "edge_guard_hours": 48.0,        # events are drawn away from the record ends
    "mc_event_counts": [10_000, 100_000, 1_000_000],
    "mc_repeats": 8,                 # MC draws per (mix, N) for the floor distribution
    "robustness_dt_hours": [1.0, 0.5, 0.25, 0.1],
    "robustness_years": [3.0, 10.0, 30.0],
    "robustness_n_realizations": 6,
    "mixes": {
        # (relative amplitude, period in hours). Periods are the standard constituents.
        "M2_only": [[1.00, 12.4206]],
        "M2_S2_semidiurnal": [[1.00, 12.4206], [0.46, 12.0000]],
        "SoCal_like_M2_S2_K1_O1": [[1.00, 12.4206], [0.46, 12.0000],
                                   [0.58, 23.9345], [0.41, 25.8193]],
        "diurnal_dominant": [[0.40, 12.4206], [0.18, 12.0000],
                             [1.00, 23.9345], [0.70, 25.8193]],
    },
    "real_series": {
        # Both are REAL-SITE series (a specific latitude/longitude and epoch), not
        # synthetic constituent sums. Neither is an instrumental strainmeter record;
        # that gap is named in the output.
        "coso_mean_stress_native_6000s": {
            "path": "data/xue_lu_zenodo/Tidal_Vol.txt",
            "companions": ["data/xue_lu_zenodo/Tidal_N_0.txt",
                           "data/xue_lu_zenodo/Tidal_N_90.txt"],
            "dt_seconds": 6000.0,
            "note": "Coso / SoCal tidal stress, K-035-era Zenodo series, NATIVE "
                    "sampling (no spline upsampling)",
        },
        "long_valley_volumetric_600s": {
            "path": "data/lv_tidal_vol.npz",
            "dt_seconds": 600.0,
            "note": "Long Valley caldera volumetric body-tide series computed from "
                    "JPL ephemerides for the real site, K-035/LV era",
        },
    },
    "elastic": {"E_young_Pa": 75e9, "nu": 0.25},
}
PROTOCOL_HASH = hashlib.sha256(
    json.dumps(PROTOCOL, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

SMOKE_OVERRIDES = {
    "n_realizations": 3,
    "n_segments_real": 3,
    "years": 2.0,
    "mc_event_counts": [10_000, 100_000],
    "mc_repeats": 2,
    "robustness_dt_hours": [1.0, 0.25],
    "robustness_years": [3.0],
    "robustness_n_realizations": 2,
}

DEG = np.pi / 180.0


def seed_of(*parts):
    """Deterministic seed from a label. NOT python's hash(): that is randomised per
    process by PYTHONHASHSEED and would make this file irreproducible."""
    h = hashlib.sha256("|".join(str(p) for p in parts).encode()).digest()
    return int.from_bytes(h[:8], "big") % (2 ** 32)


# ================================================================ conventions ===
def _interp1_at_zero(x, y):
    """Zero-crossing time by linear interpolation on a sorted sample set.

    Sample points are sorted before interpolation (the reference implementation's
    gridded-interpolant semantics); returns NaN if 0 lies outside the sampled range.
    """
    o = np.argsort(x, kind="stable")
    xs, ys = x[o], y[o]
    if 0.0 < xs[0] or 0.0 > xs[-1]:
        return np.nan
    return float(np.interp(0.0, xs, ys))


def quarter_anchored_phase(t_grid, stress, t_ev):
    """Quarter-anchored (Tanaka-style) phase assignment.

    Anchors trough / ascending zero / peak / descending zero at 0/90/180/270 deg,
    interpolates phase piecewise-linearly in time between anchors, shifts by 180 deg
    so that peak stress sits at 0, and wraps to [-180, +180].

    Faithful to the reference implementation in three details that matter and that a
    from-scratch rewrite would get wrong:
      (i)  single-sample negative intervals are removed before anchoring;
      (ii) if the cycle's minimum is found AFTER its maximum, the minimum is re-taken
           over the leading part of the segment (trough-before-peak guard);
      (iii) a segment whose ascending or descending zero cannot be interpolated is
           skipped rather than patched.
    Returns NaN outside the anchored span.
    """
    ids = (stress < 0).astype(int)
    d = np.diff(ids)
    op = np.where(d == -1)[0]           # end of a negative interval
    op2 = np.where(d == 1)[0] + 1       # start of a negative interval
    common = np.intersect1d(op, op2)    # single-sample negative intervals
    op = op[~np.isin(op, common)]
    op2 = op2[~np.isin(op2, common)]
    if len(op) and len(op2) and op[0] < op2[0]:
        op = op[1:]
    anchors_t, anchors_p = [], []
    k = 0
    n = len(t_grid)
    for i in range(len(op2) - 1):
        hi = min(op2[i + 1] + 1, n - 1)
        icut = np.arange(op2[i], hi + 1)
        sc, tc = stress[icut], t_grid[icut]
        imax = int(np.argmax(sc))
        imin = int(np.argmin(sc))
        if imin > imax:
            imin = int(np.argmin(sc[:imax + 1]))
        seg1_s, seg1_t = sc[imin:imax + 1], tc[imin:imax + 1]
        seg2_s, seg2_t = sc[imax:], tc[imax:]
        if len(seg1_s) < 2 or len(seg2_s) < 2:
            continue
        t_asc = _interp1_at_zero(seg1_s, seg1_t)
        t_dsc = _interp1_at_zero(seg2_s, seg2_t)
        if not (np.isfinite(t_asc) and np.isfinite(t_dsc)):
            continue
        base = k * 360.0
        anchors_t += [tc[imin], t_asc, tc[imax], t_dsc]
        anchors_p += [base + 0.0, base + 90.0, base + 180.0, base + 270.0]
        k += 1
    at = np.asarray(anchors_t, dtype=np.float64)
    ap = np.asarray(anchors_p, dtype=np.float64)
    if at.size < 4:
        return np.full(np.shape(t_ev), np.nan)
    o = np.argsort(at, kind="stable")
    pha = np.interp(t_ev, at[o], ap[o], left=np.nan, right=np.nan)
    return (pha - 180.0 + 180.0) % 360.0 - 180.0


def two_anchor_phase_series(stress):
    """This repository's own convention, imported unmodified.

    `coso_positive_control.phase_series` returns a per-SAMPLE phase array (NaN where
    the sample falls outside a complete trough-peak-trough cycle). Importing it rather
    than re-implementing it is deliberate: the point of the exercise is to calibrate
    the function this program actually runs, not a clean-room copy of it.
    """
    global _PHASE_SERIES
    if _PHASE_SERIES is None:
        if HERE not in sys.path:
            sys.path.insert(0, HERE)
        from coso_positive_control import phase_series  # noqa: E402 (deliberate)
        _PHASE_SERIES = phase_series
    return _PHASE_SERIES(np.asarray(stress, dtype=np.float64))


_PHASE_SERIES = None


# ================================================================= statistics ===
def harmonic_amplitude(phase_deg):
    """First-harmonic amplitude (%) and peak-aligned cosine component (%).

    a = 2*<cos(theta)>, b = 2*<sin(theta)>; amplitude = hypot(a, b). Under a fair
    clock both vanish. Bin-count ratios are deliberately NOT computed: max/min bin
    ratio is dominated by sampling discretisation of the anchors and diverges at
    coarse dt, while the fitted amplitude is stable. Amplitudes only.
    """
    ph = np.asarray(phase_deg, dtype=np.float64)
    ph = ph[np.isfinite(ph)]
    if ph.size == 0:
        return dict(amplitude_pct=float("nan"), cos_component_pct=float("nan"),
                    sin_component_pct=float("nan"), n=0)
    a = 2.0 * float(np.mean(np.cos(ph * DEG)))
    b = 2.0 * float(np.mean(np.sin(ph * DEG)))
    return dict(amplitude_pct=float(np.hypot(a, b) * 100.0),
                cos_component_pct=float(a * 100.0),
                sin_component_pct=float(b * 100.0), n=int(ph.size))


def _dist(vals, quantity, label, **extra):
    """min / median / max (+ quartiles) over realizations. Never a headline alone."""
    v = np.asarray([x for x in vals if np.isfinite(x)], dtype=np.float64)
    out = {"quantity": quantity, "label": label, "n_realizations": int(v.size)}
    if v.size:
        out.update({
            "min_pct": float(v.min()), "q1_pct": float(np.percentile(v, 25)),
            "median_pct": float(np.median(v)), "q3_pct": float(np.percentile(v, 75)),
            "max_pct": float(v.max()), "values_pct": [float(x) for x in v],
        })
    out.update(extra)
    return out


def rayleigh_floor_pct(n_events):
    """Analytic Monte Carlo floor: E|R| = sqrt(pi/N) for a perfectly fair clock."""
    return float(100.0 * np.sqrt(np.pi / float(n_events)))


# ================================================================== synthetic ===
def synth_series(mix, years, dt_hours, rng):
    t = np.arange(0.0, 24.0 * 365.0 * years, dt_hours)
    s = np.zeros_like(t)
    for amp, period in mix:
        s += amp * np.cos(2.0 * np.pi * t / period + rng.uniform(0.0, 2.0 * np.pi))
    return t, s


def instrument_response_on_grid(t, s, convention, edge_guard_h):
    """DETERMINISTIC: the uniform time grid itself is the event set. No RNG."""
    if convention == "quarter_anchored":
        keep = (t >= t[0] + edge_guard_h) & (t <= t[-1] - edge_guard_h)
        return harmonic_amplitude(quarter_anchored_phase(t, s, t[keep]))
    ph = two_anchor_phase_series(s)
    keep = (t >= t[0] + edge_guard_h) & (t <= t[-1] - edge_guard_h)
    return harmonic_amplitude(ph[keep])


def sampling_floor_on_grid(t, s, convention, n_events, rng, edge_guard_h):
    """MONTE CARLO: N uniform-random event times. Bookkeeping, not a finding."""
    if convention == "quarter_anchored":
        t_ev = rng.uniform(t[0] + edge_guard_h, t[-1] - edge_guard_h, int(n_events))
        return harmonic_amplitude(quarter_anchored_phase(t, s, t_ev))
    ph = two_anchor_phase_series(s)
    ok = np.flatnonzero(np.isfinite(ph)
                        & (t >= t[0] + edge_guard_h) & (t <= t[-1] - edge_guard_h))
    if ok.size == 0:
        return harmonic_amplitude(np.array([]))
    idx = rng.choice(ok, size=int(n_events), replace=True)
    return harmonic_amplitude(ph[idx])


# ======================================================== real tidal series ====
def _load_real(key, spec, cfg):
    """Load a real-site tidal series -> (t_hours, stress, provenance dict) or None."""
    path = os.path.join(HERE, spec["path"])
    if not os.path.exists(path):
        return None
    dt_h = spec["dt_seconds"] / 3600.0
    if path.endswith(".npz"):
        z = np.load(path)
        s = np.asarray(z["eps"], dtype=np.float64)
        t = (np.asarray(z["t_unix"], dtype=np.float64) - float(z["t_unix"][0])) / 3600.0
    else:
        # Mean (volumetric) tidal stress from the strain components, exactly as
        # coso_positive_control.load_stress builds it, but at NATIVE sampling.
        E, nu = cfg["elastic"]["E_young_Pa"], cfg["elastic"]["nu"]
        G = E / (2.0 * (1.0 + nu))
        lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
        exx = np.loadtxt(os.path.join(HERE, spec["companions"][0])) * 1e-9
        eyy = np.loadtxt(os.path.join(HERE, spec["companions"][1])) * 1e-9
        evol = np.loadtxt(path) * 1e-9
        ezz = evol - exx - eyy
        sxx = (lam + 2 * G) * exx + lam * eyy + lam * ezz
        syy = lam * exx + (lam + 2 * G) * eyy + lam * ezz
        szz = lam * exx + lam * eyy + (lam + 2 * G) * ezz
        s = (sxx + syy + szz) / 3.0
        t = np.arange(s.size, dtype=np.float64) * dt_h
    s = s - float(np.mean(s))
    prov = {
        "path": spec["path"], "dt_seconds": spec["dt_seconds"],
        "n_samples": int(s.size), "record_days": float((t[-1] - t[0]) / 24.0),
        "sha256_first_1e6_bytes": _sha_head(path),
        "note": spec["note"],
    }
    return t, s, prov


def _sha_head(path, n_bytes=1_000_000):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read(n_bytes))
    return h.hexdigest()


# ======================================================================= run ====
def run(cfg, verbose=True):
    t_start = time.time()
    res = {
        "protocol": PROTOCOL, "protocol_hash": PROTOCOL_HASH,
        "effective_config": cfg,
        "quantity_definitions": {
            "instrument_response": (
                "DETERMINISTIC. First-harmonic amplitude of the phase distribution of "
                "a UNIFORM TIME GRID pushed through the convention. No random numbers "
                "are involved; there is no sampling noise in it. THIS IS THE FINDING."),
            "sampling_floor": (
                "MONTE CARLO BOOKKEEPING. Same amplitude computed from N uniform-"
                "random event times. Carries a Rayleigh floor of sqrt(pi/N) even "
                "against a perfectly fair clock. NOT an instrument property and may "
                "never be quoted as one."),
        },
        "conventions": {
            "quarter_anchored": "quarter-anchored (Tanaka-style) phase assignment: "
                                "trough/ascending-zero/peak/descending-zero at "
                                "0/90/180/270 deg, linear in time between anchors, "
                                "shifted so peak stress = 0 deg.",
            "two_anchor": "this repository's own convention, imported unmodified from "
                          "coso_positive_control.phase_series: trough -180, peak 0, "
                          "trough +180, linear in time.",
        },
        "synthetic": {}, "real": {}, "robustness": {}, "gaps": [],
    }
    convs = PROTOCOL["conventions"]
    eg = PROTOCOL["edge_guard_hours"]

    # ---------------------------------------------------------- synthetic ----
    for mix_name, mix in PROTOCOL["mixes"].items():
        res["synthetic"][mix_name] = {"mix_amp_periodhours": mix}
        for conv in convs:
            t_mix = time.time()
            rng_c = np.random.default_rng(
                seed_of(PROTOCOL["seed_constituent_phases"], mix_name, conv))
            amps, coss = [], []
            grids = []
            for r in range(cfg["n_realizations"]):
                t, s = synth_series(mix, cfg["years"], PROTOCOL["dt_hours"], rng_c)
                m = instrument_response_on_grid(t, s, conv, eg)
                amps.append(m["amplitude_pct"])
                coss.append(m["cos_component_pct"])
                if r < 2:
                    grids.append((t, s))
            block = {
                "instrument_response": _dist(
                    amps, "instrument_response",
                    f"{conv} first-harmonic amplitude, deterministic uniform grid, "
                    f"{cfg['years']:.0f} yr @ {PROTOCOL['dt_hours']} h",
                    deterministic=True, monte_carlo=False),
                "instrument_response_cos_component": _dist(
                    coss, "instrument_response",
                    f"{conv} peak-aligned cosine component (SIGN is waveform-"
                    f"dependent), deterministic uniform grid",
                    deterministic=True, monte_carlo=False),
                "sampling_floor": {},
            }
            # ---- Monte Carlo floor at several N (bookkeeping, several repeats) --
            rng_mc = np.random.default_rng(
                seed_of(PROTOCOL["seed_monte_carlo"], mix_name, conv))
            t, s = grids[0]
            for n_ev in cfg["mc_event_counts"]:
                vals = [sampling_floor_on_grid(t, s, conv, n_ev, rng_mc, eg)
                        ["amplitude_pct"] for _ in range(cfg["mc_repeats"])]
                block["sampling_floor"][f"N={n_ev}"] = _dist(
                    vals, "sampling_floor",
                    f"{conv} amplitude from {n_ev} uniform-random event times "
                    f"(realization 0 of the constituent-phase draw). This is "
                    f"instrument response PLUS event-sampling noise; the noise part "
                    f"is the quantity being reported, and it is what a re-runner's "
                    f"own N buys them for free.",
                    deterministic=False, monte_carlo=True, n_events=int(n_ev),
                    analytic_rayleigh_floor_pct=rayleigh_floor_pct(n_ev),
                    instrument_response_same_grid_pct=float(amps[0]))
            res["synthetic"][mix_name][conv] = block
            if verbose:
                ir = block["instrument_response"]
                print(f"  [{time.time()-t_start:7.1f}s] synthetic {mix_name:<24s} "
                      f"{conv:<16s} instrument_response "
                      f"{ir['min_pct']:.2f}..{ir['max_pct']:.2f}% "
                      f"(median {ir['median_pct']:.2f}%) over "
                      f"{ir['n_realizations']} phasings  [{time.time()-t_mix:.1f}s]")

    # --------------------------------------------------------------- real ----
    found_any = False
    for key, spec in PROTOCOL["real_series"].items():
        loaded = _load_real(key, spec, PROTOCOL)
        if loaded is None:
            res["gaps"].append(f"real series {key} ({spec['path']}) NOT FOUND on disk; "
                               f"not run")
            if verbose:
                print(f"  real series {key}: MISSING ({spec['path']})")
            continue
        found_any = True
        t, s, prov = loaded
        entry = {"provenance": prov}
        n_seg = cfg["n_segments_real"]
        # disjoint, equal-length segments: the only realization axis a fixed record has
        seg_len = s.size // n_seg
        for conv in convs:
            t_c = time.time()
            amps, coss = [], []
            for i in range(n_seg):
                sl = slice(i * seg_len, (i + 1) * seg_len)
                ts, ss = t[sl] - t[sl][0], s[sl]
                m = instrument_response_on_grid(ts, ss, conv, eg)
                amps.append(m["amplitude_pct"])
                coss.append(m["cos_component_pct"])
            entry[conv] = {
                "instrument_response": _dist(
                    amps, "instrument_response",
                    f"{conv} first-harmonic amplitude on the REAL series, "
                    f"deterministic uniform grid, {n_seg} disjoint segments of "
                    f"{seg_len * spec['dt_seconds'] / 86400.0:.0f} d",
                    deterministic=True, monte_carlo=False,
                    segment_days=float(seg_len * spec["dt_seconds"] / 86400.0)),
                "instrument_response_cos_component": _dist(
                    coss, "instrument_response",
                    f"{conv} peak-aligned cosine component on the REAL series",
                    deterministic=True, monte_carlo=False),
                "sampling_floor": {},
            }
            rng_mc = np.random.default_rng(
                seed_of(PROTOCOL["seed_monte_carlo"], key, conv))
            sl = slice(0, seg_len)
            ts, ss = t[sl] - t[sl][0], s[sl]
            for n_ev in cfg["mc_event_counts"]:
                vals = [sampling_floor_on_grid(ts, ss, conv, n_ev, rng_mc, eg)
                        ["amplitude_pct"] for _ in range(cfg["mc_repeats"])]
                entry[conv]["sampling_floor"][f"N={n_ev}"] = _dist(
                    vals, "sampling_floor",
                    f"{conv} amplitude from {n_ev} uniform-random event times on "
                    f"segment 0 of the REAL series. Instrument response PLUS "
                    f"event-sampling noise; only the noise part is this quantity.",
                    deterministic=False, monte_carlo=True, n_events=int(n_ev),
                    analytic_rayleigh_floor_pct=rayleigh_floor_pct(n_ev),
                    instrument_response_same_grid_pct=float(amps[0]))
            if verbose:
                ir = entry[conv]["instrument_response"]
                print(f"  [{time.time()-t_start:7.1f}s] REAL      {key:<32s} "
                      f"{conv:<16s} instrument_response "
                      f"{ir['min_pct']:.2f}..{ir['max_pct']:.2f}% "
                      f"(median {ir['median_pct']:.2f}%) over {ir['n_realizations']} "
                      f"segments  [{time.time()-t_c:.1f}s]")
        res["real"][key] = entry
    if not found_any:
        res["gaps"].append(
            "NO real tidal series was available: every number in this file is "
            "SYNTHETIC-ONLY and says nothing about any specific study region.")
    else:
        res["gaps"].append(
            "The real series here are modelled body-tide / tidal-stress series for "
            "real sites and epochs, NOT instrumental strainmeter records, and neither "
            "carries ocean loading. A borrowed instrumental record from a published "
            "study's own site remains an open gap.")

    # --------------------------------------------------------- robustness ----
    mix_name = "SoCal_like_M2_S2_K1_O1"
    mix = PROTOCOL["mixes"][mix_name]
    for conv in convs:
        rb = {"sampling_refinement": {}, "record_length": {}}
        for dt_h in cfg["robustness_dt_hours"]:
            rng_r = np.random.default_rng(
                seed_of(PROTOCOL["seed_constituent_phases"], "dt", conv, dt_h))
            vals = []
            for _ in range(cfg["robustness_n_realizations"]):
                t, s = synth_series(mix, min(cfg["years"], 3.0), dt_h, rng_r)
                vals.append(instrument_response_on_grid(t, s, conv, eg)
                            ["amplitude_pct"])
            rb["sampling_refinement"][f"dt={dt_h}h"] = _dist(
                vals, "instrument_response",
                f"{conv} amplitude vs sampling refinement ({mix_name})",
                deterministic=True, monte_carlo=False)
            if verbose:
                d = rb["sampling_refinement"][f"dt={dt_h}h"]
                print(f"  [{time.time()-t_start:7.1f}s] robustness dt={dt_h:<6}h "
                      f"{conv:<16s} median {d['median_pct']:.2f}%")
        for yrs in cfg["robustness_years"]:
            rng_r = np.random.default_rng(
                seed_of(PROTOCOL["seed_constituent_phases"], "yr", conv, yrs))
            vals = []
            for _ in range(cfg["robustness_n_realizations"]):
                t, s = synth_series(mix, yrs, PROTOCOL["dt_hours"], rng_r)
                vals.append(instrument_response_on_grid(t, s, conv, eg)
                            ["amplitude_pct"])
            rb["record_length"][f"years={yrs}"] = _dist(
                vals, "instrument_response",
                f"{conv} amplitude vs record length ({mix_name})",
                deterministic=True, monte_carlo=False)
        res["robustness"][conv] = rb

    res["elapsed_seconds"] = round(time.time() - t_start, 1)
    return res


def summarize(res, verbose=True):
    """Print the two quantities side by side, which is the whole point of F-016."""
    print()
    print("=" * 78)
    print("F-016 SUMMARY -- two separately labelled quantities, never merged")
    print("=" * 78)
    for conv in PROTOCOL["conventions"]:
        print(f"\nconvention: {conv}")
        print(f"  {'mix / series':<34s} {'INSTRUMENT RESPONSE (det.)':<30s} "
              f"{'SAMPLING FLOOR (MC)':<24s}")
        for mix_name in PROTOCOL["mixes"]:
            b = res["synthetic"][mix_name][conv]
            ir = b["instrument_response"]
            n_max = max(int(k.split("=")[1]) for k in b["sampling_floor"])
            sf = b["sampling_floor"][f"N={n_max}"]
            print(f"  {mix_name:<34s} "
                  f"{ir['min_pct']:5.2f}-{ir['max_pct']:5.2f}% "
                  f"(med {ir['median_pct']:5.2f}%)      "
                  f"MC {sf['median_pct']:6.3f}% @ N={n_max} vs det "
                  f"{sf['instrument_response_same_grid_pct']:6.3f}% on the same grid; "
                  f"analytic floor {sf['analytic_rayleigh_floor_pct']:.3f}%")
        for key, entry in res.get("real", {}).items():
            if conv not in entry:
                continue
            ir = entry[conv]["instrument_response"]
            n_max = max(int(k.split("=")[1]) for k in entry[conv]["sampling_floor"])
            sf = entry[conv]["sampling_floor"][f"N={n_max}"]
            print(f"  REAL {key:<29s} "
                  f"{ir['min_pct']:5.2f}-{ir['max_pct']:5.2f}% "
                  f"(med {ir['median_pct']:5.2f}%)      "
                  f"MC {sf['median_pct']:6.3f}% @ N={n_max} vs det "
                  f"{sf['instrument_response_same_grid_pct']:6.3f}% on the same grid; "
                  f"analytic floor {sf['analytic_rayleigh_floor_pct']:.3f}%")
    print()
    for g in res["gaps"]:
        print("GAP: " + g)
    print()
    print("READ THIS BEFORE QUOTING ANY NUMBER ABOVE. The left column is the "
          "INSTRUMENT RESPONSE:\ndeterministic, no random numbers, and it is the "
          "finding. The right column is a MONTE\nCARLO estimate on the same grid -- "
          "instrument response PLUS event-sampling noise --\nprinted next to the "
          "deterministic value and the analytic floor sqrt(pi/N) so the two\n"
          "contributions can be told apart by inspection. Where the instrument is "
          "flat (the\nsemidiurnal controls), the MC column IS the sampling floor and "
          "nothing else: a 0.1-\n0.9% reading there is the re-runner's own N, not a "
          "bias. They are not the same\nquantity and neither substitutes for the "
          "other.")


def main(argv=None):
    ap = argparse.ArgumentParser(description="F-016 phase-clock null calibration "
                                             "(frozen protocol)")
    ap.add_argument("--smoke", action="store_true",
                    help="reduced realizations/record length; end-to-end check only")
    ap.add_argument("--out", default=None, help="output JSON path")
    args = ap.parse_args(argv)

    cfg = {k: PROTOCOL[k] for k in
           ("n_realizations", "n_segments_real", "years", "mc_event_counts",
            "mc_repeats", "robustness_dt_hours", "robustness_years",
            "robustness_n_realizations")}
    cfg["mode"] = "full"
    if args.smoke:
        cfg.update(SMOKE_OVERRIDES)
        cfg["mode"] = "smoke"

    print("=" * 78)
    print("F-016 -- null calibration of anchor-based tidal phase clocks "
          "(FROZEN PROTOCOL)")
    print(f"protocol {PROTOCOL['id']} v{PROTOCOL['version']}  "
          f"hash {PROTOCOL_HASH[:16]}...")
    print(f"mode = {cfg['mode']}; conventions = "
          f"{', '.join(PROTOCOL['conventions'])}")
    print(f"realizations per synthetic mix = {cfg['n_realizations']} "
          f"(constituent-phase draws); segments per real series = "
          f"{cfg['n_segments_real']}")
    print("quantities reported SEPARATELY: instrument_response (deterministic) vs "
          "sampling_floor (Monte Carlo)")
    print("=" * 78)

    res = run(cfg)
    out = args.out or (OUT if not args.smoke
                       else os.path.join(HERE, "results_phase_clock_null_smoke.json"))
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    summarize(res)
    print(f"\nelapsed {res['elapsed_seconds']} s  ->  {out}")
    if args.smoke:
        print("SMOKE MODE: reduced realizations and record length. These numbers are "
              "an end-to-end\nplumbing check, NOT the F-016 calibration. Run without "
              "--smoke for the frozen result.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
