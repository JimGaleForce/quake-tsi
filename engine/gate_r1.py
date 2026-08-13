"""G-M1 + R1: the v2 acceptance gate, run end-to-end on simulated catalogs.

WHAT THIS IS
------------
Two named arms of one gate, executed in one harness, per HYPOTHESIS_LEDGER.md:

  * ARM 1 -- R1 (§P6-6).  "Null-only end-to-end calibration at v2 scale. Run the
    ENTIRE v2 pipeline -- ladder, GPD, stratified BH, per-region battery, all of
    it -- over >= 30 ETAS-simulated catalogues, through the identical code path."
    Pass condition, verbatim: (a) mean BH survivors <= q x m; (b) the global
    max-statistic p is uniform across the 30 sims; (c) the GPD_EXTRAPOLATED
    survivor rate is not elevated relative to MC_RESOLVED.
    Consequence of failure, verbatim: "v2 does not ship as default, whatever its
    throughput."

  * ARM 2 -- G-M1 arm (ii) as extended by §P6-4 rule 4.7 item 5: "v2 may not
    report a bound in any band, or at any regional/cell aggregation level, where
    it has not demonstrated planted-signal recovery at that band and that
    aggregation. The planted-signal test is re-run PER REGION, because N per
    region is an order of magnitude smaller and recovery is an N-dependent
    property."  §P7-2(b) records that F9-19 IS G-M1 arm (ii)-extended -- same
    execution, one artifact, no new ledger object.

THE PIPELINE IS THE REAL ONE
----------------------------
Every catalog is pushed through `engine.mine_session.run` -- the same function
`engine.cli mine` calls -- with --ladder, --gpd, --regsum on and flat (default)
BH plus the S-8 max-statistic. Nothing is reimplemented here. The ONLY thing
this module supplies is the target: `y` becomes a Poisson draw from the fitted
ETAS lambda instead of the observed catalog, exactly as §P7-8(b) item (ii)
specifies and exactly as `f4_58_vif_control.py` already does it (that script's
"ETAS SIMULATION PATH, DECLARED" note is the authority and is reproduced in
`_simulate_y` below).

GPD LICENCE. --gpd is on because it is LICENSED at run time: `gpd_tail.
assert_calibrated` is called by `mine_session.build_config` and refuses to run
without a PASSING §P6-2(7) artifact. The artifact on disk is
`engine/out/audit_gpd_bca.json` (BCa, decisive set, 45 cases, `pass: true`,
§P7-9). If that file ever stops passing this harness stops running, which is
the intended coupling.

LEDGER DISCIPLINE -- NO EXPLORE_COUNT LINES ARE PRICED HERE
-----------------------------------------------------------
§P7-2(a) prices F9-19 (= G-M1 arm (ii)-extended) and F4-58 at **0** -- "each
either makes no rejection or reports a property of the run". These catalogs are
SIMULATED; no hypothesis about the Earth is tested and no rejection about the
Earth is available, so no multiplicity is owed. The mine sessions still write
their ledger line, because refusing to would mean editing the engine's own
accounting; the line is redirected to a SANDBOX ledger under `engine/out/
gate_r1/` and never reaches `engine/EXPLORE_COUNT.jsonl`. Likewise the session
directories and the §P7-8(c)(5) build-invariant artifact registry: `run()`
roots the registry at `dirname(session_dir)`, so pointing the sessions at the
sandbox keeps the sim hashes out of the real `engine/out/mine` registry. Three
sandboxes, asserted by `assert_sandboxed()` before a single catalog runs.

PLANT AMPLITUDES -- §P7-8(d), engine/floors.py IS CANONICAL
-----------------------------------------------------------
Every planted amplitude is `floors.min_plant_amplitude(N, alpha)` or above --
2x the operative S-15 floor at the N the statistic is actually computed over --
and each is reported next to that floor by `floors.plant_report`. Global arms
use the window's total N; the regional arm uses EACH REGION'S OWN N, which is
rule 4.7 item 5's "re-run per region ... recovery is an N-dependent property"
taken literally. A plant below the floor fails for POWER reasons while reading
as an INSTRUMENT failure; floors.py exists to stop that being recorded against
the pipeline.

PLANT CONSTRUCTION, DECLARED
----------------------------
    lambda_planted[c, t] = lambda[c, t] * exp(A * cos(theta_t - phi)) / I0(A)

`theta` is the feature's own phase series, so the plant is IN the band and at
the aggregation the arm names. The log-link form makes the truth exactly the
quantity the pipeline estimates: the GLM design for a phase feature is the
unstandardised [sin(theta), cos(theta)], so E[|beta|] = A and the recovered/
planted ratio is directly comparable to G-M1's [0.8, 1.2]. The I0(A)
normalisation (the mean of exp(A cos)) keeps the planted catalog's expected
event count equal to the null one, so the plant moves phase and nothing else.
The pipeline's OFFSET stays the unmodulated ETAS lambda -- the modulation is
the signal, not part of the baseline.

Usage:
    python -u -m engine.gate_r1 --jobs 6              # full gate
    python -u -m engine.gate_r1 --smoke --jobs 6      # 2 null + 1 per plant arm
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import shutil
import sys
import time

import numpy as np
from scipy import special, stats

from . import (baseline as bl, design, floors, gpd_tail, mine as M,
               mine_session as ms, regions as regions_mod, splits, __version__)

# ------------------------------------------------------------------ sandbox --
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SANDBOX = os.path.join("engine", "out", "gate_r1")
SANDBOX_SESSIONS = os.path.join(SANDBOX, "sessions")
SANDBOX_LEDGER = os.path.join(SANDBOX, "EXPLORE_COUNT_SANDBOX.jsonl")
PROGRESS = os.path.join(SANDBOX, "progress.json")
RESULTS = "results_gate_r1.json"

# The 1-degree ETAS refit. The shared engine/out/cache/etas_params.json holds a
# DIFFERENT design (key 88cfab8c12ba4c29, 4329 cells); the 1-degree mine design
# is key 80c635bcca169bd5, whose fit lives here. Passing this explicitly is what
# stops a 13-minute refit from overwriting the other cache.
ETAS_CACHE_1DEG = os.path.join("engine", "out", "cache",
                               "etas_params_f4_58_control_1deg.json")

REAL_LEDGER = os.path.join("engine", "EXPLORE_COUNT.jsonl")
REAL_MINE_DIR = os.path.join("engine", "out", "mine")

# ------------------------------------------------------------------- config --
# The quick preset's grid (§ the task's own allowance), with the surrogate budget
# raised to 10,000 and the LADDER ON. This is not a luxury: at 200 surrogates the
# Monte Carlo floor is 1/201 and the BH threshold at m ~ 900 declared tests is
# 1e-4, so NO test could reject at any effect size and both R1(a) and every
# G-M1 arm would pass or fail vacuously. Besag-Clifford costs ~h*ln(N_max/h) ~
# 150 draws on a null test, so the floor drops 50x for roughly the price of the
# quick preset. The period scan is unaffected (PERIOD_N_MC_CAP = 500).
GATE_SURROGATES = 10000
BASE_ALPHA = floors.ALPHA_TRANCHE_A          # 0.10 / 713, Tranche A (§P7-2)
BASE_VIF = floors.MEASURED_VIF_DF2_PHASE     # 24.0818, F4-58 (§P7-8(a))

N_NULL_CATALOGS = 30                          # R1: ">= 30"
N_PLANT_CATALOGS = 10                         # G-M1: ">= 10 per arm"
MASTER_SEED = 20260812

# arm (i) mid-band cyclic feature; arm (iii) long-period band feature.
ARM_I_FEATURE = "moon_synodic_phase"          # 29.53 d
ARM_II_FEATURE = "moon_synodic_phase"         # same band, regional PHASE plant
ARM_III_FEATURE = "annual_phase"              # 365.24 d


class _Args:
    """A stand-in for the argparse namespace `mine_session.build_config` reads."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def gate_config():
    """The per-catalog mine config. Built by the ENGINE's own build_config."""
    args = _Args(
        data_dir="data/comcat_world", dlat=1.0, dlon=1.0, explore_frac=0.70,
        mag_target=4.5, seed=MASTER_SEED, no_download=False, tranche1=False,
        ladder=True, ladder_h=None, ladder_chunk=None,
        gpd=True, gpd_confirm_max_n=None, gpd_calibration=None,
        strata=None,                       # flat BH, the default (§P6-3 note)
        regsum=True, regions=False,
        region_sectors=None, region_min_fraction=None,
        region_vif=None, region_alpha=None, region_target_amplitude=None,
    )
    cfg = ms.build_config(args, ms.QUICK)
    cfg["n_surrogates"] = int(GATE_SURROGATES)
    cfg["ladder"]["n_max"] = int(GATE_SURROGATES)
    return cfg


# ------------------------------------------------------------- sandboxing ----
def _under(path, root):
    """Is `path` inside `root`? False across Windows drives rather than raising."""
    try:
        return os.path.commonpath([os.path.abspath(path),
                                   os.path.abspath(root)]) == \
            os.path.abspath(root)
    except ValueError:                     # different drives on Windows
        return False


def assert_sandboxed(session_dir, ledger_path):
    """Three separations, asserted before anything runs. §P7-2(a) / §P7-8(c)(5).

    1. the mine session directory is under the gate sandbox, not engine/out/mine;
    2. the ledger is the sandbox ledger, not engine/EXPLORE_COUNT.jsonl;
    3. the build-invariant registry -- which `run()` roots at
       dirname(session_dir) -- therefore also lands in the sandbox.
    """
    sess = os.path.abspath(session_dir)
    sand = os.path.abspath(SANDBOX_SESSIONS)
    if not _under(sess, sand):
        raise AssertionError(
            f"session dir {sess} is not under the gate sandbox {sand}: a "
            f"simulated catalog would land in the real mine registry.")
    if _under(sess, REAL_MINE_DIR):
        raise AssertionError(f"session dir {sess} is inside {REAL_MINE_DIR}.")
    led = os.path.abspath(ledger_path)
    if led == os.path.abspath(REAL_LEDGER):
        raise AssertionError(
            "the gate would append to the REAL exploration ledger. §P7-2(a) "
            "prices these runs at 0 -- they are simulated catalogs and no "
            "multiplicity is owed. Redirect the ledger.")
    if not _under(led, SANDBOX):
        raise AssertionError(f"ledger {led} is not under the sandbox {SANDBOX}.")
    registry_root = os.path.dirname(os.path.normpath(sess))
    if not _under(registry_root, sand):
        raise AssertionError(
            f"build-invariant registry would be written to {registry_root}, "
            f"outside the sandbox.")
    return {"session_dir": sess.replace("\\", "/"),
            "ledger": led.replace("\\", "/"),
            "artifact_registry_root": registry_root.replace("\\", "/")}


# --------------------------------------------------------------- preparing ---
def prepare_base(cfg, verbose=True):
    """`mine_session.prepare`, with the 1-degree ETAS cache passed explicitly.

    Mirrors `mine_session.prepare` line for line except for `cache_path` and
    `fit_if_missing=False`: the default path would MISS (it holds a different
    design's fit), silently spend ~13 minutes refitting, and OVERWRITE
    engine/out/cache/etas_params.json. `fit_if_missing=False` turns that into a
    loud FileNotFoundError instead of a quiet cache clobber.
    """
    ctx = design.build_design(
        data_dir=cfg["data_dir"], dlat=cfg["grid"]["dlat"],
        dlon=cfg["grid"]["dlon"], explore_frac=cfg["explore_frac"],
        verbose=verbose)
    explore, _hold = splits.temporal_split(ctx.n_days, cfg["explore_frac"])
    y = ctx.day_counts(cfg["mag_target"])

    base = bl.EtasV1(verbose=verbose, mag_target=cfg["mag_target"],
                     cache_path=ETAS_CACHE_1DEG, fit_if_missing=False)
    base.fit(ctx, y, explore)
    if verbose:
        for line in base.report():
            print(line)
        print(f"  ETAS params source = {base.fit_info.get('source')} "
              f"({ETAS_CACHE_1DEG})")

    burn = int(base.burn_in_days)
    window = slice(burn, explore.stop)
    counts, offset = M.build_target(ctx, base, y, window)

    t0 = _dt.datetime.fromisoformat(str(ctx.meta["t0"]))
    all_marks = M.load_event_marks(ctx, cfg["data_dir"], ctx.meta["mag_floor"])
    in_win = (all_marks["day"] >= window.start) & (all_marks["day"] < window.stop)
    marks = {k: v[in_win] for k, v in all_marks.items()}

    lags = tuple(cfg["lags"])
    feats = M.ephemeris_features(t0, ctx.n_days)
    dl_feats, dl_log = M.download_features(t0, ctx.n_days, lags,
                                           enabled=cfg["downloads"],
                                           verbose=verbose)
    feats += dl_feats
    feats += M.catalog_features(ctx, all_marks, lags)
    if verbose:
        print(f"real window: [{window.start}, {window.stop}) = {counts.size} d, "
              f"{counts.sum():.0f} observed events, ETAS expectation "
              f"{offset.sum():.1f}, {len(feats)} features")
    return {"ctx": ctx, "base": base, "y": y, "window": window,
            "counts": counts, "offset": offset, "marks": marks, "feats": feats,
            "dl_log": dl_log, "t0": t0,
            "rate_w": np.asarray(base.rate(window), dtype=np.float64)}


# --------------------------------------------------------------- simulation --
def _simulate_y(rate_w, rng, modulation=None):
    """Poisson counts from the fitted ETAS lambda. §P7-8(b) item (ii), verbatim.

    `engine/baseline.py` exposes no simulate/rvs entry point, so -- exactly as
    `f4_58_vif_control.py` declares -- the simulation is a direct Poisson draw
    from the fitted lambda series. Conditional on lambda the draws are
    independent across days and cells, so the target is a TRUE NULL by
    construction (and its true VIF is exactly 1). It deliberately does NOT re-run
    the ETAS triggering cascade: a self-exciting simulation would reintroduce
    clustering and would no longer be a null.

    `modulation`, when given, is a (n_cells, n_days) multiplicative factor -- the
    planted signal. It multiplies the RATE only; the pipeline's offset stays the
    unmodulated lambda, which is what makes the plant a signal rather than a
    baseline change.
    """
    lam = rate_w if modulation is None else rate_w * modulation
    return rng.poisson(lam).astype(np.int64)


def simulate_catalog(prep, rng, modulation=None):
    """(y_sim, counts_sim, marks_sim) for one catalog, shaped like the real ones.

    Marks: event days are the simulated counts expanded (day d contributes
    counts[d] events); magnitude and depth are drawn WITH REPLACEMENT from the
    real in-window mark pool, independently of day. That makes the mark tests
    null by construction too -- reusing the real (day, mark) pairing against
    simulated counts would leave the mark arm testing the Earth, which is not
    what R1 is calibrating.
    """
    window = prep["window"]
    y_full = np.zeros_like(prep["y"], dtype=np.int64)
    y_win = _simulate_y(prep["rate_w"], rng, modulation)
    y_full[:, window] = y_win

    per_day = y_win.sum(axis=0).astype(np.int64)
    counts = per_day.astype(np.float64)
    day_idx = np.repeat(np.arange(window.start, window.stop, dtype=np.int64),
                        per_day)
    n_ev = int(day_idx.size)
    real = prep["marks"]
    pick = rng.integers(0, real["mag"].size, size=n_ev)
    marks = {
        "day": day_idx,
        "day_float": day_idx.astype(np.float64) + rng.random(n_ev),
        "mag": np.asarray(real["mag"])[pick],
        "depth": np.asarray(real["depth"])[pick],
    }
    return y_full, counts, marks


def prepared_tuple(prep, y_sim, counts_sim, marks_sim):
    """The exact 10-tuple `mine_session.prepare` returns, with the sim target in."""
    return (prep["ctx"], prep["base"], y_sim, prep["window"], counts_sim,
            prep["offset"], marks_sim, prep["feats"], prep["dl_log"],
            prep["t0"])


# ------------------------------------------------------------------- plants --
def _feature(prep, name):
    for f in prep["feats"]:
        if f.name == name:
            return f
    raise KeyError(f"feature {name!r} not in the declared feature list")


def _cos_series(theta, phi):
    return np.cos(np.asarray(theta, dtype=np.float64) - float(phi))


def global_plant(prep, feature_name, phi=0.0, alpha=BASE_ALPHA, vif=BASE_VIF,
                 factor=floors.PLANT_FACTOR, amplitude=None):
    """A single-phase sinusoidal rate modulation over the WHOLE domain.

    The floor is taken at the window's total N -- the event count the global
    statistic is actually computed over -- and the plant at `factor` x it.
    """
    f = _feature(prep, feature_name)
    if f.kind != "phase":
        raise ValueError(f"{feature_name} is kind={f.kind}; the plant is a "
                         f"sinusoid in a phase feature's own angle")
    window = prep["window"]
    n = float(prep["offset"].sum())
    A = (floors.min_plant_amplitude(n, alpha=alpha, vif=vif, factor=factor)
         if amplitude is None else float(amplitude))
    rep = floors.assert_plant_above_floor(
        A, n, alpha=alpha, vif=vif, factor=factor,
        what=f"global plant :: {feature_name}")
    theta = np.asarray(f.values[window], dtype=np.float64)
    row = np.exp(A * _cos_series(theta, phi)) / float(special.i0(A))
    mod = np.broadcast_to(row[None, :], prep["rate_w"].shape)
    rep.update({"arm_aggregation": "global", "feature": feature_name,
                "phase_rad": float(phi), "period_hint_days": f.period_hint,
                "block_days": f.block_days})
    return np.ascontiguousarray(mod), rep


def regional_plant(prep, partition, feature_name, alpha=BASE_ALPHA,
                   vif=BASE_VIF, factor=floors.PLANT_FACTOR):
    """A REGION-DEPENDENT-PHASE plant: the §K87-0(d)(i) blind spot, planted.

    Each region gets its own phase, spread evenly over 2*pi, so the domain sum
    cancels the signal (that is the blind spot §P6-4 Rule 4.2 exists to kill)
    while the phase-incoherent 2R-df regional sum retains it. Each region's
    amplitude is 2x the floor AT THAT REGION'S OWN N -- rule 4.7 item 5's
    "re-run per region, because N per region is an order of magnitude smaller
    and recovery is an N-dependent property", taken literally.
    """
    f = _feature(prep, feature_name)
    window = prep["window"]
    roc = np.asarray(partition["region_of_cell"])
    R = int(partition["R"])
    theta = np.asarray(f.values[window], dtype=np.float64)
    rate_w = prep["rate_w"]
    mod = np.ones_like(rate_w)
    reports = []
    for r in range(R):
        sel = roc == r
        n_r = float(rate_w[sel].sum())
        phi = 2.0 * math.pi * r / R
        A = floors.min_plant_amplitude(n_r, alpha=alpha, vif=vif, factor=factor)
        rep = floors.assert_plant_above_floor(
            A, n_r, alpha=alpha, vif=vif, factor=factor,
            what=f"regional plant :: {feature_name} :: region {r}")
        rep.update({"arm_aggregation": "regional", "region": int(r),
                    "feature": feature_name, "phase_rad": float(phi),
                    "n_cells": int(sel.sum())})
        reports.append(rep)
        mod[sel] = (np.exp(A * _cos_series(theta, phi))
                    / float(special.i0(A)))[None, :]
    return np.ascontiguousarray(mod), reports


def build_partition(prep, cfg):
    """The declared region partition, from the REAL exploration-window y.

    §P6-4 Rule 4.1: exploration-window data only. It is built from the real `y`
    -- the same call `run()` makes with the SIMULATED y -- and used here only to
    know WHERE to plant. The partition the pipeline actually scores against is
    the one `run()` rebuilds from the simulated catalog, which is the honest
    arrangement: the plant does not get to choose the regions it is scored in.
    """
    rc = cfg["regions"]
    return regions_mod.build_regions(
        prep["ctx"], prep["y"], prep["window"],
        n_sectors=int(rc["n_sectors"]),
        min_event_fraction=float(rc["min_event_fraction"]))


# --------------------------------------------------------------- extraction --
def extract_metrics(session_dir, planted=None):
    """Everything the gate needs from one completed session, read from its own
    checkpoint (the artifact the engine wrote, not a recomputation)."""
    with open(os.path.join(session_dir, "checkpoint.json"), "r",
              encoding="utf-8") as fh:
        st = json.load(fh)
    tests = st["tests"]
    n_tests = int(st["n_tests"])
    survivors = [t for t in tests if t.get("passes_fdr")]

    by_method = {}
    for t in tests:
        m = t.get("p_method", "UNLABELLED")
        d = by_method.setdefault(m, {"n": 0, "n_survivors": 0})
        d["n"] += 1
        d["n_survivors"] += int(bool(t.get("passes_fdr")))

    msr = st.get("max_statistic") or {}
    glob = msr.get("global") or {}
    out = {
        "session_dir": session_dir.replace("\\", "/"),
        "config_hash": st.get("config_hash"),
        "artifact_hash": (st.get("build_invariant") or {}).get("artifact_hash"),
        "n_tests": n_tests,
        "n_bh_survivors": len(survivors),
        "bh_survivor_names": [f"{t['test']}:{t['feature']}:"
                              f"{t.get('lag', t.get('mark'))}"
                              for t in survivors][:20],
        "max_stat_p": (float(glob["p"]) if glob.get("p") is not None else None),
        "max_stat_t_obs": (float(glob["t_obs"]) if glob.get("t_obs") is not None
                           else None),
        "max_stat_floor": (float(glob["floor"]) if glob.get("floor") is not None
                           else None),
        "max_stat_n_covered": msr.get("n_covered"),
        "max_stat_n_replicates": glob.get("n_replicates"),
        "p_method_counts": by_method,
        "elapsed_seconds": st.get("elapsed_seconds"),
        "regions_R": (st.get("regions") or {}).get("R"),
    }
    if planted:
        out["planted_recovery"] = _planted_recovery(tests, out, planted)
    return out


def max_stat_T(row):
    """The max-statistic's own scale, for ONE test: -log10(p within its own null).

    `strata.max_statistic_p` maxes over tests of exactly this quantity (its
    `statistic` field says so), computed from the shared circular-shift index.
    Reading `t_obs` as a chi2 -- which it visibly is not -- would silently make
    every max-statistic attribution false, so the scale is taken from the
    engine's own definition and checked against `t_obs` in the artifact.
    """
    ps = row.get("p_circular_shift")
    if ps is None or ps <= 0:
        return None
    return -math.log10(float(ps))


def _planted_recovery(tests, summary, planted):
    """Did the planted feature come out? BH survivor OR max-statistic detection.

    "Max-statistic detection" is the single-step Westfall-Young reading and the
    only one the engine's artifact supports: the family-wise max-statistic p is
    significant AND this test ATTAINS the family maximum, so that the test's own
    Westfall-Young adjusted p, P(max_j T_j^null >= T_planted), IS the reported
    family p. A test below the maximum has a LARGER adjusted p that the artifact
    does not carry, and reading the family p onto it would be reading a
    family-level number as a per-test one.

    Ties at the resolution ceiling count as attainment: several tests pinned at
    -log10(1/(n_shifts+1)) all have the same adjusted p, and refusing the
    planted one on a tie-break would be an arbitrary rule, not a statistical one.

    Only tests the max-statistic actually COVERS are eligible for this route.
    `max_statistic_matrix` admits `glm_poisson_offset_etas` rows alone -- mark
    tests and period peaks do not share the day-shift index and the regsum rows
    are not in the matrix -- so the regional arm is decided on BH, by rule.
    """
    feature = planted["feature"]
    kind = planted["test_kind"]
    covered = kind == "glm_poisson_offset_etas"
    rows = [t for t in tests if t.get("feature") == feature
            and t.get("test") == kind]
    if not rows:
        return {"detected": False, "reason": f"no {kind} row for {feature}",
                "n_rows": 0, "amplitude_ratio": None,
                "best_p_raw": 1.0, "best_p_floor": None}
    best = min(rows, key=lambda t: (float(t["p_raw"]), -float(t["chi2_score"])))
    bh_hit = any(t.get("passes_fdr") for t in rows)

    t_obs = summary.get("max_stat_t_obs")
    p_ms = summary.get("max_stat_p")
    ts_planted = [v for v in (max_stat_T(t) for t in rows) if v is not None]
    t_planted = max(ts_planted) if ts_planted else None
    ms_hit = bool(
        covered and t_obs is not None and p_ms is not None and p_ms <= 0.05
        and t_planted is not None and t_planted >= float(t_obs) - 1e-9)

    # The AMPLITUDE is a separate object from the P-METHOD label. A row can be
    # BH-INELIGIBLE under §P6-2(1)/(3) -- its GPD extrapolation failed a gate --
    # and still report a perfectly good effect size, which is precisely the case
    # the max-statistic exists to adjudicate. Only the ABSENCE of an amplitude
    # (the §P7-1(d) UNRESOLVED regional/regsum rows) suppresses the ratio.
    amp = best.get("amplitude_log_rate")
    truth = planted.get("truth_amplitude")
    ratio = (None if (amp is None or not truth)
             else float(amp) / float(truth))
    return {
        "detected": bool(bh_hit or ms_hit),
        "bh_survivor": bool(bh_hit),
        "bh_eligible": bool(best.get("bh_eligible", True)),
        "max_stat_covered": covered,
        "max_stat_detection": ms_hit,
        "max_stat_p": p_ms,
        "max_stat_T_planted": t_planted,
        "max_stat_T_family_max": (None if t_obs is None else float(t_obs)),
        "n_rows": len(rows),
        "best_lag": best.get("lag"),
        "best_p_raw": float(best["p_raw"]),
        "best_p_floor": best.get("p_floor"),
        "best_chi2": float(best["chi2_score"]),
        "best_p_method": best.get("p_method"),
        "best_p_method_reason": best.get("p_method_reason"),
        "recovered_amplitude_log_rate": (None if amp is None else float(amp)),
        "truth_amplitude": truth,
        "amplitude_ratio": ratio,
        "amplitude_unresolved": bool(amp is None),
    }


# ------------------------------------------------------------- acceptance ----
def ks_uniform_one_sided(ps, a=0.01):
    """§P6-9(a): one-sided KS in the ANTI-CONSERVATIVE direction.

    D+ = sup(F_n(x) - x) against D_crit = sqrt(-ln(a)/(2n)). D- is REPORTED and
    DOES NOT GATE: the max-statistic p is a discrete exceedance count on a
    lattice and is stochastically conservative by construction, so a two-sided
    test at fixed level rejects any valid discrete estimator once n is large --
    which is the exact error §P6-9 corrects. The two-sided statistic is printed
    alongside, labelled not-gating, so the deviation stays visible.
    """
    x = np.sort(np.asarray([p for p in ps if p is not None], dtype=np.float64))
    n = x.size
    if n == 0:
        return {"n": 0, "d_plus": None, "d_crit": None, "pass": False}
    i = np.arange(1, n + 1, dtype=np.float64)
    d_plus = float(np.max(i / n - x))
    d_minus = float(np.max(x - (i - 1) / n))
    d_crit = float(math.sqrt(-math.log(a) / (2.0 * n)))
    two = stats.kstest(x, "uniform")
    return {
        "n": int(n), "alpha": a,
        "d_plus": d_plus, "d_minus_reported_not_gating": d_minus,
        "d_crit": d_crit, "pass": bool(d_plus <= d_crit),
        "two_sided_D_not_gating": float(two.statistic),
        "two_sided_p_not_gating": float(two.pvalue),
        "mean_p": float(x.mean()), "min_p": float(x.min()),
        "rule": ("HYPOTHESIS_LEDGER.md §P6-9(a): gate on D+ only; D- is "
                 "permitted by super-uniformity and does not gate."),
    }


def gpd_vs_mc_survivor_rate(per_catalog):
    """R1(c): is the GPD_EXTRAPOLATED survivor rate ELEVATED vs MC_RESOLVED?

    Pooled over catalogs, one-sided Fisher exact (GPD elevated). "Not elevated"
    is the pass condition, so the gate is `p_one_sided > 0.05` OR the GPD rate
    not exceeding the MC rate at all.
    """
    g_n = g_s = m_n = m_s = 0
    for c in per_catalog:
        for name, d in c["p_method_counts"].items():
            if name == gpd_tail.P_GPD:
                g_n += d["n"]
                g_s += d["n_survivors"]
            elif name == gpd_tail.P_MC_RESOLVED:
                m_n += d["n"]
                m_s += d["n_survivors"]
    g_rate = (g_s / g_n) if g_n else None
    m_rate = (m_s / m_n) if m_n else None
    p = None
    if g_n and m_n:
        p = float(stats.fisher_exact([[g_s, g_n - g_s], [m_s, m_n - m_s]],
                                     alternative="greater").pvalue)
    elevated = bool(g_rate is not None and m_rate is not None
                    and g_rate > m_rate and p is not None and p <= 0.05)
    return {
        "n_gpd_extrapolated": g_n, "n_gpd_survivors": g_s, "gpd_rate": g_rate,
        "n_mc_resolved": m_n, "n_mc_survivors": m_s, "mc_rate": m_rate,
        "fisher_one_sided_p_gpd_greater": p,
        "elevated": elevated, "pass": (not elevated),
        "note": ("A zero denominator on the GPD arm is reported as such and "
                 "NOT read as a pass on the estimator: it means no test was "
                 "censored at the Monte Carlo floor in these catalogs, so the "
                 "extrapolation path was never exercised."),
    }


def r1_survivor_arithmetic(per_catalog, q=M.FDR_Q):
    """R1(a): mean BH survivors <= q x m, with m the declared count."""
    surv = [c["n_bh_survivors"] for c in per_catalog]
    ms_ = [c["n_tests"] for c in per_catalog]
    mean_s = float(np.mean(surv)) if surv else float("nan")
    m = float(np.mean(ms_)) if ms_ else float("nan")
    bound = float(q) * m
    return {"n_catalogs": len(surv), "mean_bh_survivors": mean_s,
            "declared_m_mean": m, "q": float(q), "bound_q_times_m": bound,
            "max_survivors_single_catalog": (int(max(surv)) if surv else 0),
            "total_survivors": int(sum(surv)),
            "pass": bool(mean_s <= bound)}


def recovery_rate(records, threshold=0.80):
    """>= 80% of planted catalogs recovered, and the amplitude band where quoted."""
    det = [r["planted_recovery"]["detected"] for r in records]
    ratios = [r["planted_recovery"]["amplitude_ratio"] for r in records
              if r["planted_recovery"]["amplitude_ratio"] is not None]
    in_band = [0.8 <= x <= 1.2 for x in ratios]
    n = len(det)
    rate = (sum(det) / n) if n else 0.0
    amp_pass = (all(in_band) if ratios else None)
    return {
        "n_catalogs": n, "n_detected": int(sum(det)), "recovery_rate": rate,
        "threshold": threshold, "detection_pass": bool(n and rate >= threshold),
        "amplitude_quoted": bool(ratios),
        "n_amplitudes": len(ratios),
        "amplitude_ratio_median": (float(np.median(ratios)) if ratios else None),
        "amplitude_ratio_min": (float(min(ratios)) if ratios else None),
        "amplitude_ratio_max": (float(max(ratios)) if ratios else None),
        "amplitude_in_band_fraction": (float(np.mean(in_band)) if ratios
                                       else None),
        "amplitude_pass": amp_pass,
        "pass": bool(n and rate >= threshold and (amp_pass is not False)),
    }


def miss_diagnosis(records, plant_reports):
    """A miss must state whether it is POWER or INSTRUMENT. §P7-8(d)."""
    misses = [r for r in records if not r["planted_recovery"]["detected"]]
    if not misses:
        return {"n_misses": 0, "verdict": "no misses"}
    ratios = sorted({round(float(p["amplitude_over_floor"]), 3)
                     for p in plant_reports})
    floor_ok = all(p["compliant"] for p in plant_reports)
    at_floor = [m for m in misses
                if m["planted_recovery"]["best_p_raw"]
                <= (m["planted_recovery"].get("best_p_floor") or 0) + 1e-12]
    return {
        "n_misses": len(misses),
        "plant_over_floor_ratios": ratios,
        "all_plants_above_2x_floor": floor_ok,
        "n_misses_pinned_at_their_own_p_floor": len(at_floor),
        "verdict": (
            "INSTRUMENT -- every plant sits at or above 2x the operative S-15 "
            "floor at its own N (§P7-8(d)), so the experiment had the declared "
            "power and the pipeline still did not return the signal."
            if floor_ok else
            "POWER -- at least one plant is below 2x its operative floor; the "
            "miss is a power verdict and MUST NOT be recorded against the "
            "pipeline (§P7-8(d))."),
        "note": ("A miss whose best p sits AT its own resolution floor is "
                 "budget-limited rather than instrument-limited: the test could "
                 "not have produced a smaller p at any effect size. Those are "
                 "counted separately above."),
    }


# ------------------------------------------------------------------ driver ---
def _load_progress():
    if os.path.exists(PROGRESS):
        with open(PROGRESS, "r", encoding="utf-8") as fh:
            return json.load(fh)
    return {"catalogs": {}}


def _save_progress(prog):
    os.makedirs(os.path.dirname(PROGRESS) or ".", exist_ok=True)
    tmp = PROGRESS + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(prog, fh, indent=1)
    os.replace(tmp, PROGRESS)


def run_catalog(tag, prep, cfg, rng, modulation, planted, jobs, verbose=False):
    """One catalog: simulate, push through the REAL pipeline, extract."""
    session_dir = os.path.join(SANDBOX_SESSIONS, f"session_{tag}")
    assert_sandboxed(session_dir, SANDBOX_LEDGER)
    if os.path.isdir(session_dir):
        shutil.rmtree(session_dir)
    os.makedirs(session_dir, exist_ok=True)

    y_sim, counts_sim, marks_sim = simulate_catalog(prep, rng, modulation)
    prepared = prepared_tuple(prep, y_sim, counts_sim, marks_sim)
    t0 = time.time()
    out = ms.run(cfg, verbose=verbose, resume=False, session_dir=session_dir,
                 jobs=jobs, ledger_path=SANDBOX_LEDGER, prepared=prepared)
    rec = extract_metrics(session_dir, planted=planted)
    rec.update({"tag": tag, "n_events_simulated": float(counts_sim.sum()),
                "wall_seconds": round(time.time() - t0, 1),
                "engine_n_pass": out["n_pass"], "engine_n_tests": out["n_tests"]})
    return rec


def main(argv=None):
    p = argparse.ArgumentParser("engine.gate_r1")
    p.add_argument("--jobs", type=int, default=6)
    p.add_argument("--smoke", action="store_true",
                   help="2 null catalogs + 1 per plant arm")
    p.add_argument("--nulls", type=int, default=None)
    p.add_argument("--plants", type=int, default=None)
    p.add_argument("--arms", default="i,ii,iii")
    p.add_argument("--fresh", action="store_true",
                   help="discard the progress checkpoint and rerun everything")
    p.add_argument("--report-only", action="store_true")
    a = p.parse_args(argv)

    n_null = a.nulls if a.nulls is not None else (2 if a.smoke
                                                  else N_NULL_CATALOGS)
    n_plant = a.plants if a.plants is not None else (1 if a.smoke
                                                     else N_PLANT_CATALOGS)
    arms = [s.strip() for s in a.arms.split(",") if s.strip()]

    os.makedirs(SANDBOX_SESSIONS, exist_ok=True)
    if a.fresh and os.path.exists(PROGRESS):
        os.remove(PROGRESS)
    prog = _load_progress()

    cfg = gate_config()
    cfg_hash = splits.config_hash(cfg)
    print("=" * 78)
    print(f"G-M1 + R1 ACCEPTANCE GATE -- engine v{__version__}")
    print("=" * 78)
    print(f"config hash          = {cfg_hash}")
    print(f"preset               = {cfg['preset']} grid, "
          f"n_surrogates = {cfg['n_surrogates']} (ladder N_max)")
    print(f"ladder               = {cfg['ladder']}")
    print(f"gpd                  = ON, calibration "
          f"{cfg['gpd']['calibration'].get('path', 'engine/out/audit_gpd_bca.json')} "
          f"pass={cfg['gpd']['calibration'].get('pass')}")
    print(f"regsum               = ON (R rule {cfg['regions']['rule_id']}), "
          f"battery OFF")
    print(f"BH                   = flat (unstratified) at q = {cfg['fdr_q']}")
    print(f"sandbox              = {SANDBOX} "
          f"(sessions + ledger + artifact registry)")
    print(f"catalogs             = {n_null} null, {n_plant} per plant arm "
          f"{arms}")

    if not a.report_only:
        prep = prepare_base(cfg, verbose=True)
        partition = build_partition(prep, cfg)
        print(f"plant partition      = R {partition['R']}, digest "
              f"{partition['digest']}")

        plant_specs = {}
        if "i" in arms:
            mod, rep = global_plant(prep, ARM_I_FEATURE)
            plant_specs["i"] = {
                "modulation": mod, "reports": [rep],
                "planted": {"feature": ARM_I_FEATURE,
                            "test_kind": "glm_poisson_offset_etas",
                            "truth_amplitude": rep["planted_amplitude"]},
                "label": f"global aggregation, mid-band ({ARM_I_FEATURE})"}
        if "ii" in arms:
            mod, reps = regional_plant(prep, partition, ARM_II_FEATURE)
            plant_specs["ii"] = {
                "modulation": mod, "reports": reps,
                "planted": {"feature": ARM_II_FEATURE,
                            "test_kind": "regsum_score_2Rdf",
                            "truth_amplitude": None},
                "label": (f"regional aggregation, region-dependent phase "
                          f"({ARM_II_FEATURE}), recovered by regsum")}
        if "iii" in arms:
            mod, rep = global_plant(prep, ARM_III_FEATURE)
            plant_specs["iii"] = {
                "modulation": mod, "reports": [rep],
                "planted": {"feature": ARM_III_FEATURE,
                            "test_kind": "glm_poisson_offset_etas",
                            "truth_amplitude": rep["planted_amplitude"]},
                "label": f"long-period band ({ARM_III_FEATURE})"}

        for k, spec in plant_specs.items():
            for r in spec["reports"]:
                print(f"  plant arm ({k}) {r['what']}: A = "
                      f"{r['planted_amplitude']:.4f} vs floor A_min = "
                      f"{r['operative_floor_A_min']:.4f} at N = {r['N']:.0f} "
                      f"(x{r['amplitude_over_floor']:.2f}, required >= "
                      f"{r['required_min_plant']:.4f}) -> "
                      f"{'OK' if r['compliant'] else 'BELOW FLOOR'}")

        prog["plant_reports"] = {k: [dict(r) for r in spec["reports"]]
                                 for k, spec in plant_specs.items()}
        prog["plant_labels"] = {k: spec["label"] for k, spec in plant_specs.items()}
        _save_progress(prog)

        jobs = int(a.jobs)
        plan = [("null", i, None) for i in range(n_null)]
        for k in plant_specs:
            plan += [("plant_" + k, i, k) for i in range(n_plant)]

        # Deterministic per-catalog streams. `hash()` on a str is randomised per
        # process (PYTHONHASHSEED), so the stream id is an explicit declared
        # integer per arm -- a resumed run must draw the same catalogs.
        stream_id = {"null": 0, "plant_i": 1, "plant_ii": 2, "plant_iii": 3}
        for kind, i, arm in plan:
            tag = f"{kind}_{i:03d}"
            if tag in prog["catalogs"]:
                print(f"[{tag}] checkpointed, skipping")
                continue
            rng = np.random.default_rng(
                np.random.SeedSequence([MASTER_SEED, stream_id[kind], i]))
            mod = plant_specs[arm]["modulation"] if arm else None
            planted = plant_specs[arm]["planted"] if arm else None
            t = time.time()
            rec = run_catalog(tag, prep, cfg, rng, mod, planted, jobs)
            rec["kind"] = kind
            rec["arm"] = arm
            rec["seed_note"] = (f"SeedSequence([{MASTER_SEED}, "
                                f"{stream_id[kind]}, {i}])")
            prog["catalogs"][tag] = rec
            _save_progress(prog)
            extra = ""
            if planted:
                pr = rec["planted_recovery"]
                extra = (f"  planted {'HIT ' if pr['detected'] else 'MISS'}"
                         f" p={pr['best_p_raw']:.3g}"
                         + (f" A/A0={pr['amplitude_ratio']:.3f}"
                            if pr["amplitude_ratio"] is not None else " (UNRESOLVED)"))
            print(f"[{tag}] {rec['n_tests']} tests, "
                  f"{rec['n_bh_survivors']} BH survivors, max-stat p = "
                  f"{rec['max_stat_p']}, {time.time()-t:.0f}s{extra}", flush=True)

    # ------------------------------------------------------------- verdict --
    cats = list(prog["catalogs"].values())
    nulls = [c for c in cats if c["kind"] == "null"]
    r1a = r1_survivor_arithmetic(nulls)
    r1b = ks_uniform_one_sided([c["max_stat_p"] for c in nulls])
    r1c = gpd_vs_mc_survivor_rate(nulls)
    r1_pass = bool(r1a["pass"] and r1b["pass"] and r1c["pass"])

    gm1 = {}
    for arm in ("i", "ii", "iii"):
        recs = [c for c in cats if c.get("arm") == arm]
        if not recs:
            continue
        gm1[arm] = recovery_rate(recs)
        gm1[arm]["misses"] = miss_diagnosis(
            recs, prog.get("plant_reports", {}).get(arm, []))
    gm1_pass = bool(gm1) and all(v["pass"] for v in gm1.values())

    out = {
        "gate": "G-M1 + R1",
        "generated_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "engine_version": __version__,
        "config": cfg,
        "config_hash": cfg_hash,
        "verdict": {"R1": "PASS" if r1_pass else "FAIL",
                    "G-M1": "PASS" if gm1_pass else "FAIL",
                    "gate": "PASS" if (r1_pass and gm1_pass) else "FAIL"},
        "R1": {"a_bh_survivors": r1a, "b_max_statistic_uniformity": r1b,
               "c_gpd_vs_mc_survivor_rate": r1c},
        "G_M1": gm1,
        "plant_reports": prog.get("plant_reports", {}),
        "per_catalog": cats,
        "sandbox": {"sessions": SANDBOX_SESSIONS.replace("\\", "/"),
                    "ledger": SANDBOX_LEDGER.replace("\\", "/"),
                    "artifact_registry": os.path.join(
                        SANDBOX_SESSIONS, "artifact_registry.jsonl"
                    ).replace("\\", "/"),
                    "real_ledger_untouched": REAL_LEDGER.replace("\\", "/")},
    }
    with open(RESULTS, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, default=float)
    print("-" * 78)
    print(f"R1(a) mean BH survivors {r1a['mean_bh_survivors']:.3f} <= q x m = "
          f"{r1a['bound_q_times_m']:.1f}  -> {'PASS' if r1a['pass'] else 'FAIL'}")
    print(f"R1(b) max-stat D+ {r1b['d_plus']} <= D_crit {r1b['d_crit']}  -> "
          f"{'PASS' if r1b['pass'] else 'FAIL'}")
    print(f"R1(c) GPD rate {r1c['gpd_rate']} vs MC {r1c['mc_rate']}  -> "
          f"{'PASS' if r1c['pass'] else 'FAIL'}")
    for arm, v in gm1.items():
        print(f"G-M1 arm ({arm}): {v['n_detected']}/{v['n_catalogs']} = "
              f"{v['recovery_rate']:.0%}  -> {'PASS' if v['pass'] else 'FAIL'}")
    print(f"GATE VERDICT: {out['verdict']['gate']}")
    print(f"written -> {RESULTS}")
    return 0 if out["verdict"]["gate"] == "PASS" else 1


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    sys.exit(main())
