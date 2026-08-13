"""Session driver for `mine` mode: checkpointed sweep -> report.md + stubs.json.

Split out of engine/mine.py (which holds the statistics) so the orchestration --
resume logic, task ordering, ledger accounting and report writing -- reads in one
sitting. Everything this module writes carries the generator-not-evidence banner.


RNG POLICY (v2) -- READ THIS BEFORE COMPARING TWO RUNS
------------------------------------------------------
There are two randomness regimes and they do NOT produce the same numbers. This is
deliberate and the difference is the price of parallelism.

* `--jobs 1` (the default) is the AUDIT BASELINE for the GLM and MARK tests. One
  `default_rng(cfg["seed"])` is created in `run()` and threaded through those tasks
  in a fixed order -- GLM sweep in feature order, then the mark tests in feature
  order. Every draw depends on every draw before it. Output is bit-identical to the
  pre-parallel engine for the same config and seed.

* THE PERIOD SCAN IS THE EXCEPTION, AND IT IS DELIBERATE (v2 phase 1b, §P6-5).
  It no longer touches the shared stream in EITHER mode. Its 2 x n_mc surrogates
  are subtasks, and surrogate i draws from an index-addressable child of that
  null's key sequence (`engine.mine.surrogate_seed_sequence`), so:
    - `--jobs 1` and `--jobs 32` give bit-identical period results;
    - the chunk size is an execution knob and provably changes nothing;
    - and the pre-1b shared-stream period numbers are SUPERSEDED. That was the
      price, paid knowingly: the old stream was reproducible only for one process
      running one fixed task order, and the scan was 90.7% of the sweep as a single
      indivisible task (Amdahl ceiling 1.10x). Determinism independent of scheduling
      is worth more than bit-identity with a build that could not be parallelised.
  Second-order consequence, stated rather than discovered later: in `--jobs 1` the
  shared stream is no longer advanced by the period scan, so the post-FDR aliasing
  audit (which runs after it and threads the same stream) also draws different
  surrogates than it did before phase 1b.

* `--jobs != 1` derives an INDEPENDENT `numpy.random.SeedSequence` per task from
  the sha256 of the CANONICAL JSON of the complete test key -- master seed, feature
  name, test kind, lag, rung, null type, region, bin width (`engine.mine.test_key`,
  `engine.mine.seed_sequence_for`). Not a raw tuple: a hashed schema can gain a
  field without re-keying the fields that already exist, and `region`/`bin_width`
  are present and defaulted NOW so that per-region or per-binning strata later do
  not change anybody else's stream. The parent asserts every declared stream has a
  distinct digest and refuses to run on a collision (`assert_task_keys_unique`).

  Nothing is threaded, so nothing depends on the order in which the pool happens
  to finish: `--jobs 2` and `--jobs 4` and `--jobs 32` all return the same numbers.
  They are NOT the same numbers as `--jobs 1`, because a per-task stream cannot
  reproduce a single shared stream. Both are valid nulls; only one is the baseline.

* LADDER RUNGS cross process boundaries. When `--ladder` is on, each chunk (rung)
  of a test's sequential draw takes an index-addressable child SeedSequence of that
  test's sequence (`engine.mine.rung_seed_sequence`), never a continuing generator.
  Rung 2 is therefore the same draws whether rungs 1 and 2 ran in one process or
  two, or whether rung 1 ran at all in this invocation. `SeedSequence.spawn` is
  deliberately NOT used: it hands out children in call order, which is exactly the
  process-boundary dependence we must not have.

* ORDER-DETERMINISTIC REDUCTION. Results are sorted by `order_key` -- a function of
  the test key alone -- before any aggregation, and `order_key` is the final
  tiebreak of every downstream sort, so BH ties and table rank ties break by test
  key and never by arrival order.

  The post-FDR aliasing audit is also given per-task derived streams in parallel
  mode (keyed on the audit's own checkpoint key), so it does not depend on how many
  tests happened to survive FDR ahead of it.

`--jobs` is an EXECUTION detail and is deliberately absent from the config hash:
the same sweep run on 1 core or 30 cores is the same sweep, and a resumable session
must not be orphaned by a core count. `--ladder`, by contrast, changes the sampling
rule, so its parameters DO enter the config hash (see `build_config`).

Parallelism is process-based (`spawn`, so Windows-safe): worker entry points are
module-level, arguments are picklable, and workers NEVER write. The parent owns
every checkpoint write, the report, the stubs and the append-only ledger.
"""

from __future__ import annotations

import concurrent.futures as _cf
from concurrent.futures.process import BrokenProcessPool
import datetime as _dt
import hashlib
import json
import math
import multiprocessing as _mp
import os
import time
import traceback

import numpy as np

from . import (baseline as bl, datasets, design, mine as M, splits, __version__)
from scipy import stats

from . import gpd_tail, floors, regions as regions_mod, strata as strata_mod
from . import circstat, dispositions as disp, marks_ext, observer as observer_mod
# `write_report` binds a LOCAL named `floors` (a list of p-value floors),
# which shadows the module for the whole function body. This alias is how
# the S-15(c) section reaches the module without renaming a local that
# predates it.
from . import floors as floors_mod

QUICK = {
    "n_surrogates": 200, "n_periods": 800, "n_peaks": 5,
    "lags": (0, 1, 3, 7, 14, 30), "label": "quick",
}
OVERNIGHT = {
    "n_surrogates": 10000, "n_periods": 3000, "n_peaks": 10,
    "lags": tuple(range(0, 31)), "label": "overnight",
}
DEFAULT = {
    "n_surrogates": 1000, "n_periods": 1500, "n_peaks": 8,
    "lags": (0, 1, 2, 3, 5, 7, 10, 14, 21, 30), "label": "default",
}

PERIOD_MIN, PERIOD_MAX = 2.0, 4000.0


# ---------------------------------------------------------------- RNG policy ---
def task_seed_sequence(master_seed, feature, kind, lag=None, null_type=None,
                       region=None, bin_width=None):
    """Deterministic SeedSequence for a task. See module docstring + M.test_key.

    Derived from the sha256 of the COMPLETE test key's canonical JSON -- master
    seed, feature, kind, lag, rung, null type, region, bin width -- and nothing
    else. In particular NOT from task index, submission order or completion order,
    so the stream a task sees is independent of how the pool schedules it.
    """
    return M.seed_sequence_for(M.test_key(master_seed, feature, kind, lag=lag,
                                          null_type=null_type, region=region,
                                          bin_width=bin_width))


def task_rng(master_seed, feature, kind, lag=None, **kw):
    return np.random.default_rng(
        task_seed_sequence(master_seed, feature, kind, lag, **kw))


def rung_rng_factory(master_seed, feature, kind, lag=None, null_type=None,
                     region=None, bin_width=None):
    """Index-addressable per-rung generators for one test's ladder.

    Rung i is reproducible on its own: it does not matter whether rungs 0..i-1 were
    drawn in this process, another process, or (in a resumed run) never at all.
    """
    parent = task_seed_sequence(master_seed, feature, kind, lag,
                               null_type=null_type, region=region,
                               bin_width=bin_width)
    return lambda i: np.random.default_rng(M.rung_seed_sequence(parent, i))


def ledger_record(cfg, n_tests, session_dir):
    """The one line a completed session appends to the exploration ledger.

    `n_declared_tests` is stated explicitly and by name alongside `n_tests` so a
    reader never has to infer whether the count meant DECLARED or SURVIVING. They
    are the same number here -- the declared family size -- and saying so costs one
    field and removes an ambiguity that would otherwise be resolved by guessing.

    Parent-only: workers never see this function, let alone the ledger path.
    """
    return {
        "ts": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "kind": "mine",
        "n_tests": int(n_tests),
        "n_declared_tests": int(n_tests),
        "hash": splits.config_hash(cfg),
        "config": cfg,
        "session_dir": session_dir.replace("\\", "/"),
    }


def test_floor(t, n_surr):
    """Smallest p this test could ever have reported, and where that floor comes from.

    Two DIFFERENT floors, reported separately because they have different causes and
    different cures:

      * the ENUMERATION floor, 1/(n_admissible_shifts + 1). The circular-shift null
        is an exhaustive enumeration, so this is a property of the RECORD LENGTH.
        More surrogates cannot lower it, and the ladder never touches it.
      * the MONTE CARLO floor, 1/(N_max + 1), a property of the surrogate BUDGET.
        This is the one the ladder economises on -- and note it is set by N_max, the
        pre-declared cap, not by however many draws a laddered test actually made:
        a test that stopped at rung 1 could still have gone to the cap.

    Where p_raw is max(shift, bootstrap), the effective floor is the LARGER of the
    two. Where the block bootstrap is the sole null (periodic features) it is the
    Monte Carlo floor alone. Side effect: annotates `t`.
    """
    lad = t.get("ladder_stop")
    mc_n = int(lad["n_max"]) if lad else int(n_surr)
    floor_mc = 1.0 / (mc_n + 1.0)
    floor_enum = None
    if t["test"] == "lomb_scargle_peak":
        # the period scan has no enumeration null at all -- see the DEVIATION note
        floor_mc = 1.0 / (int(t.get("n_surrogates", n_surr)) + 1.0)
    elif t.get("p_circular_shift") is not None:
        floor_enum = 1.0 / (int(t.get("n_surrogates", n_surr)) + 1.0)
    t["floor_enumeration"] = floor_enum
    t["floor_monte_carlo"] = float(floor_mc)
    periodic_only = str(t.get("null", "")).startswith("block-bootstrap")
    if floor_enum is None or periodic_only:
        return float(floor_mc)
    return float(max(floor_enum, floor_mc))


def assert_task_keys_unique(keys):
    """Fail loudly on a test-key collision rather than silently sharing a stream."""
    digests = {}
    for k in keys:
        d = M.test_key_digest(k)
        if d in digests:
            raise RuntimeError(
                f"test-key digest collision: {digests[d]} and {k} hash to {d}. "
                f"Two tests would share a random stream. Refusing to run.")
        digests[d] = k
    return digests


def resolve_jobs(jobs):
    """--jobs N; 0 means auto = max(1, cpu_count() - 2). Never returns < 1."""
    j = int(jobs)
    if j == 0:
        return max(1, (os.cpu_count() or 1) - 2)
    return max(1, j)


# ----------------------------------------------------------------- task bodies -
# One body per test kind, called by BOTH drivers. The sequential driver hands them
# the single shared rng (which is what makes --jobs 1 bit-identical to v1); the
# parallel driver hands them per-task derived streams. Nothing else differs.
def glm_task(f, window, counts, offset, n_surr, rng_for_lag, ladder=None, seed=0,
             gpd=None):
    rows = []
    want_draws = bool(gpd and gpd.get("enabled"))
    for lag in f.lags:
        rng = rng_for_lag(lag)
        X = f.design(window, lag)
        fit = M.glm_fit(X, counts, offset)
        S = M.score_stat_all_shifts(X, counts, offset)
        # EXHAUSTIVE enumeration of admissible circular shifts. Never laddered:
        # its resolution floor is a property of the record length, not of a
        # surrogate budget, so sequential stopping has nothing to stop.
        p_shift, n_used = M.empirical_p(S, n_surr, rng)
        lad = None
        if ladder:
            key = M.test_key(seed, f.name, "glm", lag=lag,
                             null_type="block_bootstrap")
            rr = rung_rng_factory(seed, f.name, "glm", lag,
                                  null_type="block_bootstrap")
            lad = M.besag_clifford_p(
                lambda b, g: M.score_stat_block_bootstrap(
                    X, counts, offset, int(b), g, mean_block=f.block_days),
                S[0], n_max=M.ladder_n_max(ladder, key), rung_rng=rr,
                h=int(ladder.get("h", 25)), chunk=int(ladder.get("chunk", 200)),
                keep_draws=want_draws)
            # POPPED, never checkpointed: see besag_clifford_p(keep_draws=...).
            Sb = lad.pop("draws", None)
            p_boot = lad["p"]
            mc_n_max = int(lad["n_max"])
        else:
            Sb = M.score_stat_block_bootstrap(X, counts, offset, n_surr, rng,
                                              mean_block=f.block_days)
            p_boot = M.bootstrap_p(S[0], Sb)
            mc_n_max = int(n_surr)
        # circular shifts are provably powerless for a deterministic cycle, so
        # for periodic features the block bootstrap IS the null; elsewhere both
        # are valid and the more conservative one is reported.
        p = p_boot if f.periodic else max(p_shift, p_boot)
        b = np.asarray(fit["beta"])
        amp = float(np.hypot(*b)) if f.kind == "phase" else float(abs(b[0]))
        row = {
            "feature": f.name, "family": f.family, "kind": f.kind, "lag": int(lag),
            "test": "glm_poisson_offset_etas", "df": f.df,
            "beta": fit["beta"], "se": fit["se"], "amplitude_log_rate": amp,
            "pct_rate_modulation": 100.0 * (math.exp(amp) - 1.0),
            "bits_per_event": fit["bits_per_event"],
            "chi2_score": float(S[0]), "p_parametric": M.chi2_sf(S[0], f.df),
            "p_circular_shift": p_shift, "p_block_bootstrap": p_boot,
            "block_days": f.block_days,
            "null": ("block-bootstrap (periodic feature)" if f.periodic
                     else "max(circular-shift, block-bootstrap)"),
            "p_raw": p, "n_surrogates": min(n_used, n_surr),
            "converged": fit["converged"],
        }
        if lad is not None:
            row["ladder_stop"] = lad
        # §P6-2(4)+(5). The block bootstrap is the ONLY permitted substrate; the
        # exhaustive circular-shift enumeration is forbidden and enters only as
        # `other_p`, which stays binding wherever it is the larger (= more
        # conservative) of the two. For a PERIODIC feature the shift null is
        # provably powerless and is not part of p_raw at all, so it is not passed.
        row.update(M.gpd_tail.price_row(
            p_boot, mc_n_max, S[0], Sb if want_draws else None,
            "block_bootstrap",
            task_rng(seed, f.name, "gpd", lag) if want_draws else None,
            other_p=(0.0 if f.periodic else p_shift),
            other_floor=(0.0 if f.periodic else 1.0 / (n_used + 1.0)),
            enabled=want_draws))
        rows.append(row)
    return rows


# --------------------- §P6-4 Rule 4.2: the 2R-df phase-incoherent regional sum -
def _region_unresolved_rows(f, X, C, O, part, reg_cfg):
    """Per-region amplitudes, COMPUTED and reported UNRESOLVED.

    §P6-4 Rule 4.3 as AMENDED by §P7-1(d): these are not quoted and never enter the
    ranked list. They are carried so that a future lower-Mc regional catalogue can
    pass the conditional clause on the arithmetic instead of on a re-argument, and
    each row states, in itself, the rule that makes it unresolved and the floor it
    failed.
    """
    per_region = M.score_stat_regions(X, C, O)
    rows = []
    for r, meta in enumerate(part["regions"]):
        fit = M.glm_fit(X, C[r], O[r])
        b = np.asarray(fit["beta"])
        amp = float(np.hypot(*b)) if f.kind == "phase" else float(abs(b[0]))
        n_r = float(C[r].sum())
        floor = regions_mod.a_min(reg_cfg["vif"], reg_cfg["alpha"], n_r)
        clears = amp <= 0.0 or False  # a per-region amplitude clears only via `floor`
        clears = bool(amp >= 0.0 and floor <= float(reg_cfg["target_amplitude"]))
        rows.append({
            "region": int(r),
            "lon_lo": meta["lon_lo"], "lon_hi": meta["lon_hi"],
            "n_events_in_window": n_r,
            "chi2_score_region": float(per_region[r]),
            "amplitude_log_rate": amp,
            "pct_rate_modulation": 100.0 * (math.exp(amp) - 1.0),
            "beta": fit["beta"], "se": fit["se"],
            "a_min_formula": floor,
            "a_min_pct": 100.0 * floor,
            "s15": "MEASURABLE" if clears else "UNMEASURABLE",
            # The label is on the row itself, not inferred by the reader.
            "p_method": regions_mod.UNRESOLVED,
            "p_method_reason": regions_mod.UNRESOLVED_RULE,
            "quotable": bool(clears),
        })
    return rows


def regsum_task(f, window, C, O, n_surr, rng_for_lag, part, reg_cfg, seed=0):
    """One 2R-df regional-sum test per (feature, lag). Rule 4.2, priced 1 test each.

    The tested statistic is the SUM of the per-region score statistics, and every
    surrogate -- circular shift and block bootstrap alike -- is pushed through the
    identical sum (`M.score_stat_regsum_*`). A null built region-by-region and then
    summed after the fact would not be the null of this statistic.
    """
    R = int(part["R"])
    rows = []
    for lag in f.lags:
        rng = rng_for_lag(lag)
        X = f.design(window, lag)
        S = M.score_stat_regsum_all_shifts(X, C, O)
        p_shift, n_used = M.empirical_p(S, n_surr, rng)
        Sb = M.score_stat_regsum_block_bootstrap(X, C, O, n_surr, rng,
                                                 mean_block=f.block_days)
        p_boot = M.bootstrap_p(S[0], Sb)
        p = p_boot if f.periodic else max(p_shift, p_boot)
        df = R * int(f.df)
        row = {
            "feature": f.name, "family": f.family, "kind": f.kind, "lag": int(lag),
            "test": "regsum_score_2Rdf", "df": df, "R": R,
            "region": None,
            "region_rule_id": part["rule_id"], "region_digest": part["digest"],
            "chi2_score": float(S[0]),
            "p_parametric": M.chi2_sf(S[0], df),
            "p_circular_shift": p_shift, "p_block_bootstrap": p_boot,
            "block_days": f.block_days,
            "null": ("block-bootstrap (periodic feature), regional sum"
                     if f.periodic else
                     "max(circular-shift, block-bootstrap), regional sum"),
            "p_raw": p, "n_surrogates": min(n_used, n_surr),
            # The sum is a DETECTION statistic over the full event count, not a
            # per-region estimator -- §P7-1(d). It has no single amplitude and does
            # not get one invented for the ranked list.
            "amplitude_log_rate": None, "pct_rate_modulation": None,
            "amplitude_note": (
                "NO SINGLE AMPLITUDE. The 2R-df sum is a detection statistic over "
                "the full event count, not an estimator (§P7-1(d)). Per-region "
                "amplitudes are carried in `regions_unresolved` and are UNRESOLVED."),
            "regions_unresolved": _region_unresolved_rows(f, X, C, O, part, reg_cfg),
        }
        rows.append(row)
    return rows


def region_battery_task(f, window, C, O, r, n_surr, rng_for_lag, part, reg_cfg,
                        seed=0):
    """Optional per-region battery (§P6-4 Rule 4.3/4.4), one test per (feature, lag).

    Priced at R x the declared count and stratified on `region` (§P6-4 Rule 4.5;
    engine.strata already carries the region axis). Every row states its own
    §P7-1(b) floor and its S-15 verdict, and a row below its floor is scored
    neither way (Rule 4.4).
    """
    rows = []
    meta = part["regions"][int(r)]
    n_r = float(C[int(r)].sum())
    floor = regions_mod.a_min(reg_cfg["vif"], reg_cfg["alpha"], n_r)
    measurable = floor <= float(reg_cfg["target_amplitude"])
    for lag in f.lags:
        rng = rng_for_lag(lag)
        X = f.design(window, lag)
        c, o = C[int(r)], O[int(r)]
        fit = M.glm_fit(X, c, o)
        S = M.score_stat_all_shifts(X, c, o)
        p_shift, n_used = M.empirical_p(S, n_surr, rng)
        Sb = M.score_stat_block_bootstrap(X, c, o, n_surr, rng,
                                          mean_block=f.block_days)
        p_boot = M.bootstrap_p(S[0], Sb)
        p = p_boot if f.periodic else max(p_shift, p_boot)
        b = np.asarray(fit["beta"])
        amp = float(np.hypot(*b)) if f.kind == "phase" else float(abs(b[0]))
        rows.append({
            "feature": f.name, "family": f.family, "kind": f.kind, "lag": int(lag),
            "test": "glm_poisson_offset_etas_region", "df": int(f.df),
            "region": int(r),
            "lon_lo": meta["lon_lo"], "lon_hi": meta["lon_hi"],
            "region_rule_id": part["rule_id"], "region_digest": part["digest"],
            "n_events_in_window": n_r,
            "beta": fit["beta"], "se": fit["se"],
            "amplitude_log_rate": amp,
            "pct_rate_modulation": 100.0 * (math.exp(amp) - 1.0),
            "bits_per_event": fit["bits_per_event"],
            "chi2_score": float(S[0]), "p_parametric": M.chi2_sf(S[0], f.df),
            "p_circular_shift": p_shift, "p_block_bootstrap": p_boot,
            "block_days": f.block_days,
            "null": ("block-bootstrap (periodic feature)" if f.periodic
                     else "max(circular-shift, block-bootstrap)"),
            "p_raw": p, "n_surrogates": min(n_used, n_surr),
            "converged": fit["converged"],
            "a_min_formula": floor, "a_min_pct": 100.0 * floor,
            "s15": "MEASURABLE" if measurable else "UNMEASURABLE",
            # Rule 4.4: an unmeasurable cell is scored NEITHER way and may not emit
            # a stub. `bh_eligible` False keeps it out of the BH vector entirely.
            "bh_eligible": bool(measurable),
            "p_method": (gpd_tail.P_MC_RESOLVED if measurable
                         else regions_mod.UNRESOLVED),
            "p_method_reason": (
                "region clears its own §P7-1(b) floor at the declared tranche alpha"
                if measurable else
                "§P6-4 Rule 4.4 + §P7-1(d): region does not clear its own §P7-1(b) "
                "floor (A_min = %.4f > target %.4f); scored neither way, no stub."
                % (floor, float(reg_cfg["target_amplitude"]))),
            "amplitude_label": ("SELECTION-BIASED UPPER BOUND"
                                if measurable else regions_mod.UNRESOLVED),
        })
    return rows


def marks_task(f, window, fe_day, marks, n_surr, rng, ladder=None, seed=0,
                gpd=None):
    vals = f.values[window][fe_day]
    rows = []
    for mk in ("mag", "depth"):
        kk = f"mark_{mk}"
        rr, lad_cfg = None, None
        if ladder:
            key = M.test_key(seed, f.name, kk, null_type="block_bootstrap")
            # N_max is resolved from the DECLARED stratum of this test key, before
            # any statistic is looked at (M.ladder_n_max enforces that).
            lad_cfg = dict(ladder, n_max=M.ladder_n_max(ladder, key))
            rr = rung_rng_factory(seed, f.name, kk, null_type="block_bootstrap")
        r = M.mark_test(vals, marks[mk], f.kind, n_surr, rng,
                        ladder=lad_cfg, rung_rng=rr, gpd=gpd,
                        gpd_rng=(task_rng(seed, f.name, "gpd:" + kk)
                                 if (gpd and gpd.get("enabled")) else None))
        r.update({"feature": f.name, "family": f.family, "kind": f.kind,
                  "mark": mk, "n_events": int(marks[mk].size)})
        rows.append(r)
    return rows


# =========================================================== TRANCHE B kinds ===
# Three NEW TEST KINDS, each with the same task/key/checkpoint/strata discipline as
# the three that already exist: one body per kind, called by both drivers, addressed
# by a canonical `M.test_key`, checkpointed by the parent, and routed to a declared
# §P6-3 stratum by `strata._test_kind`. They are inserted ONLY when
# `cfg["tranche_b"]["enabled"]` is set, so a default session's config, hash, task
# list, declared streams and resumable sessions are byte-identical to phase 2a.
#
# NOTHING HERE RUNS ON REAL DATA IN THIS BUILD. §P7-14(d)/§P7-15(b): B's BUILD is
# authorized; B's RUN is gated on its per-statistic G-M1 clearance, which is what
# `engine/recovery_b.py` produces on SIMULATED catalogues only.
TRANCHE_B_RULE_ID = "TRANCHE-B-v1"


def feature_phase(f, window, lag=0):
    """The phase vector a circular statistic is computed on, for ONE feature.

    Declared construction, one branch each, no alternatives (S-9):
      * `kind='phase'`  -> the feature's own angle. This is the definition.
      * `kind='linear'` with a `period_hint` -> the DAY INDEX folded on the declared
        period, `2 pi t / P`. Not the feature's values: F9-01 and F9-04 are
        statistics of a PHASE DISTRIBUTION, and a linear cyclic feature's phase is
        the position in its own declared cycle, which is what "17 cyclic features"
        in F9-01's price means.
      * anything else -> refused. A circular statistic on a feature with no declared
        cycle would be measuring the record, not a cycle.
    """
    if f.kind == "phase":
        v = f.values[window.start - int(lag):window.stop - int(lag)]
        return np.mod(np.asarray(v, dtype=np.float64), 2 * np.pi)
    if f.period_hint:
        t = np.arange(window.start, window.stop, dtype=np.float64) - float(lag)
        return np.mod(2 * np.pi * t / float(f.period_hint), 2 * np.pi)
    raise ValueError(
        "%s: a circular statistic needs a declared cycle; this feature is "
        "kind=%r with no period_hint. Refusing to invent one." % (f.name, f.kind))


def moment2_task(f, window, counts, offset, n_surr, rng, seed=0):
    """F9-01: the second circular moment, as the 2-df score on the doubled angle."""
    th = feature_phase(f, window)
    r = circstat.second_moment_test(th, counts, offset, n_surr, rng,
                                    block_days=f.block_days, periodic=f.periodic)
    r.update({"feature": f.name, "family": f.family, "kind": f.kind, "lag": 0,
              "mark": None,
              "period_days": (float(f.period_hint) if f.period_hint else None),
              "tranche": "B", "catalog_entry": "F9-01"})
    r.update(_window_clause(f))
    return [r]


def omnibus_task(f, window, counts, offset, n_surr, rng, seed=0):
    """F9-04: Kuiper V and Watson U^2, two rows, under the engine's two nulls."""
    th = feature_phase(f, window)
    rows = circstat.omnibus_test(th, counts, offset, n_surr, rng,
                                 block_days=f.block_days, periodic=f.periodic)
    wc = _window_clause(f)
    for r in rows:
        r.update({"feature": f.name, "family": f.family, "kind": f.kind, "lag": 0,
                  "mark": None,
                  "period_days": (float(f.period_hint) if f.period_hint else None),
                  "tranche": "B", "catalog_entry": "F9-04"})
        r.update(wc)
    return rows


def marksx_task(f, window, fe_day, marks, ext_marks, n_surr, rng, seed=0,
                subdaily_values=None, subdaily=False):
    """F9-10: the mark axis over the 7 declared marks, optionally at EVENT TIMES.

    `subdaily_values` is the feature RE-DERIVED at each event's own `day_float`
    (`marks_ext.event_time_feature_values`). When it is present the row is genuinely
    sub-daily and says so; when it is absent the feature is read at its day value and
    the row says THAT, because a day-binned feature evaluated at an event has not
    escaped the sinc and a row that implied otherwise would be the single most
    misleading thing this tranche could emit.
    """
    if subdaily and subdaily_values is not None:
        vals = np.asarray(subdaily_values, dtype=np.float64)
        if f.kind == "phase":
            vals = np.mod(vals, 2 * np.pi)
    else:
        vals = f.values[window][fe_day]
        subdaily = False
    rows = []
    for mk in marks_ext.SCORED_MARK_NAMES:
        if mk not in ext_marks:
            continue
        r = M.mark_test(vals, ext_marks[mk], f.kind, n_surr, rng)
        n_ev = int(np.asarray(ext_marks[mk]).size)
        # §P7-10(c): F9-10 declares its floor from VIF_mark BEFORE it runs.
        fl = floors.mark_floor_report(n_ev, None, feature="%s x %s" % (f.name, mk))
        r.update({
            "feature": f.name, "family": f.family, "kind": f.kind, "mark": mk,
            "mark_axis": "F9-10", "n_events": n_ev, "lag": None,
            "subdaily": bool(subdaily),
            "subdaily_note": marks_ext.SUBDAILY_NOTE,
            "mark_definition": marks_ext.MARK_DEFINITIONS.get(mk),
            "rho_min": fl["rho_min"], "mark_floor": fl,
            "statistic_over_floor": (float(abs(r["statistic"])) / fl["rho_min"]
                                     if fl["rho_min"] > 0 else None),
            "tranche": "B", "catalog_entry": "F9-10",
        })
        r.update(_window_clause(f))
        rows.append(r)
    return rows


def window_clause_census(tests, record_days):
    """S-15(c) over every row that carries a period -- LABEL, never a grid change.

    §P7-16 settled the period-scan grid question in the direction that costs nothing
    and hides nothing: **`PERIOD_MAX` stays at 4,000 d and peaks between record/3 and
    4,000 d are REPORTED as UNMEASURABLE-BY-WINDOW.** Clamping the grid would have
    been a change to a declared config value (a new declaration under §P6-3 rule 5)
    AND it would have made the affected band invisible rather than labelled -- a scan
    that silently cannot look somewhere is worse than one that looks and says the
    answer is not identifiable there.

    So this annotates and counts; it does not filter, and it does not touch any
    p-value. An UNMEASURABLE-BY-WINDOW row is scored NEITHER WAY in the S-15 headline
    fraction -- it is not a null and it is not a detection -- and the row keeps its
    place in the BH vector because the declared denominator is the declaration.
    """
    cut = floors.max_identifiable_period(float(record_days))
    flagged = []
    for t in tests:
        p = t.get("period_days")
        if not p or not np.isfinite(float(p)) or float(p) <= 0:
            continue
        rep = floors.window_report(float(p), float(record_days),
                                   name=str(t.get("feature")))
        t["s15c"] = rep
        t["s15c_verdict"] = rep["verdict"]
        t["s15c_scored"] = rep["scored"]
        if rep["verdict"] == floors.UNMEASURABLE_BY_WINDOW:
            flagged.append({"test": t["test"], "feature": t.get("feature"),
                            "period_days": float(p),
                            "cycles_in_window": rep["cycles_in_window"],
                            "p_raw": float(t.get("p_raw", float("nan")))})
    n_per = sum(1 for t in tests if t.get("period_days"))
    return {
        "clause": floors.WINDOW_CLAUSE,
        "clause_source": floors.WINDOW_CLAUSE_SOURCE,
        "record_days": float(record_days),
        "cut_period_days": cut,
        "period_scan_max_days": float(PERIOD_MAX),
        "grid_unchanged": True,
        "n_rows_with_a_period": int(n_per),
        "n_unmeasurable_by_window": len(flagged),
        "rows": sorted(flagged, key=lambda r: -r["period_days"]),
        "disposition": (
            "§P7-16: the period grid is UNCHANGED (PERIOD_MAX = %.0f d) and every "
            "peak above %.0f d is REPORTED UNMEASURABLE-BY-WINDOW. Labelled, not "
            "clamped: a scan that silently cannot look somewhere is worse than one "
            "that looks and says the answer is not identifiable there."
            % (float(PERIOD_MAX), cut)),
    }


def _window_clause(f):
    """S-15(c) on every Tranche B row that has a declared period (§P7-14(c))."""
    if not f.period_hint:
        return {"s15c": None}
    rep = floors.window_report(float(f.period_hint), name=f.name)
    return {"s15c": rep, "s15c_verdict": rep["verdict"],
            "s15c_scored": rep["scored"]}


def max_statistic_matrix(feats, window, counts, offset, tests, min_shift=30):
    """The joint null for S-8's max-statistic: one column per test, ONE SHARED index.

    §P6-3(4) makes the max-statistic -- not BH -- the anti-repartition guarantee, and
    it only IS one if the columns share a replicate index: the statistic has to be
    "the largest of these tests ON THE SAME SURROGATE WORLD", not "the largest of
    these tests each on its own surrogates". The circular-shift enumeration provides
    exactly that for the GLM family: shift k is one surrogate world, every GLM test
    is evaluated in it, and the enumeration is exhaustive rather than sampled, so
    the joint null is EXACT under arbitrary dependence between the tests -- no PRDS
    assumption and no correction factor.

    COVERAGE IS DECLARED, NOT ASSUMED. Mark tests resample along the EVENT sequence
    and period peaks resample the residual SERIES; neither shares the day-shift
    index, so neither can be entered into this matrix without inventing a
    correspondence that does not exist. They are reported as NOT COVERED rather than
    quietly folded in, and `n_covered` is printed next to the declared count.
    """
    fmap = {f.name: f for f in feats}
    cols, obs = [], []
    idx = []
    admissible = None
    for i, t in enumerate(tests):
        if t.get("test") != "glm_poisson_offset_etas":
            continue
        f = fmap.get(t["feature"])
        if f is None:
            continue
        S = M.score_stat_all_shifts(f.design(window, int(t["lag"])), counts, offset)
        if admissible is None:
            n = S.size
            a = np.arange(n)
            admissible = a[(a >= min_shift) & (a <= n - min_shift)]
        cols.append(S[admissible])
        obs.append(float(S[0]))
        idx.append(i)
    if not cols:
        return None, None, []
    return np.column_stack(cols), np.array(obs), idx


def confirm_gpd_candidates(tests, feats, window, counts, offset, marks, fe_day,
                           cfg, verbose=True):
    """§P6-2(6), the load-bearing rule. A GPD survivor is a CANDIDATE, not a stub.

    BH has already run with the GPD rows entered at their CI-upper. Every survivor
    still labelled GPD_EXTRAPOLATED is marked CANDIDATE-REQUIRES-BRUTE-FORCE and gets
    a TARGETED single-test Monte Carlo at N >= 10/p_gpd, resolving its p directly.
    ONLY that brute-force p may be written to stubs.json -- which is affordable
    precisely because selection has reduced the extrapolated survivors to a handful,
    and is the entire economic argument for spending the budget here.

    A candidate whose required N exceeds the declared ceiling stays a CANDIDATE and
    emits NO stub. That is the honest outcome: the confirmation was not affordable,
    so the claim was not made.
    """
    cap = int((cfg.get("gpd") or {}).get("confirm_max_n", 0) or 0)
    fmap = {f.name: f for f in feats}
    out = []
    for t in tests:
        if not t.get("passes_fdr") or t.get("p_method") != gpd_tail.P_GPD:
            continue
        t["candidate_label"] = gpd_tail.CANDIDATE_LABEL
        p_gpd = float(t.get("p_bh"))
        n_need = gpd_tail.brute_force_n(p_gpd)
        rec = {"feature": t["feature"], "test": t["test"], "p_gpd_ci_upper": p_gpd,
               "n_required": n_need, "cap": cap}
        if t["test"] == "lomb_scargle_peak":
            rec.update({"status": "NOT CONFIRMED (not attempted)",
                        "reason": ("each period-scan surrogate is a FULL "
                                   "periodogram; a targeted brute force at "
                                   f"N = {n_need} is not affordable and is not "
                                   "attempted. No stub is emitted.")})
        elif n_need > cap:
            rec.update({"status": "NOT CONFIRMED (over the declared ceiling)",
                        "reason": (f"N >= 10/p_gpd = {n_need} exceeds the declared "
                                   f"confirmation ceiling {cap}. No stub is "
                                   f"emitted; the row stays a candidate.")})
        else:
            t0 = time.perf_counter()
            rng = task_rng(int(cfg["seed"]), t["feature"], "gpd_confirm",
                           t.get("lag"))
            if t["test"] == "glm_poisson_offset_etas":
                f = fmap[t["feature"]]
                X = f.design(window, int(t["lag"]))
                Sb = M.score_stat_block_bootstrap(X, counts, offset, n_need, rng,
                                                  mean_block=f.block_days)
                p_bf = M.bootstrap_p(t["chi2_score"], Sb)
            else:
                f = fmap[t["feature"]]
                vals = f.values[window][fe_day]
                r = M.mark_test(vals, marks[t["mark"]], f.kind, n_need, rng)
                p_bf = float(r["p_block_bootstrap"])
            rec.update({"status": "CONFIRMED (brute force)", "p_brute_force": p_bf,
                        "n_drawn": n_need,
                        "seconds": round(time.perf_counter() - t0, 1),
                        "floor": 1.0 / (n_need + 1.0)})
            t["p_brute_force"] = p_bf
            t["p_method"] = gpd_tail.P_MC_RESOLVED
            t["p_method_reason"] = (f"GPD candidate CONFIRMED by a targeted "
                                    f"single-test Monte Carlo at N = {n_need} "
                                    f">= 10/p_gpd (§P6-2(6)); the brute-force p is "
                                    f"the only number written to stubs.json")
        t["gpd_confirmation"] = rec
        out.append(rec)
        if verbose:
            print(f"  §P6-2(6) {gpd_tail.CANDIDATE_LABEL}: {t['feature']} "
                  f"({t['test']}) p_gpd={p_gpd:.3g}, N>=10/p={n_need} -> "
                  f"{rec['status']}"
                  + (f", p_brute={rec['p_brute_force']:.4g}"
                     if "p_brute_force" in rec else ""))
    return out


def period_grid(cfg):
    """The trial-period grid. A pure function of the config, so every subtask that
    needs it can rebuild it instead of shipping it."""
    return np.exp(np.linspace(math.log(PERIOD_MIN), math.log(PERIOD_MAX),
                              int(cfg["n_periods"])))


# ---- the period scan as N subtasks (phase 1b) ------------------------------
# Task keys:
#   "period_obs"            observed periodogram + peak picking (no randomness)
#   "psurr:<null>:<start>"  max power of surrogates [start, start+chunk) of <null>
# The parent then assembles. `period_scan` (the checkpoint key downstream code and
# the report already know) is written by the parent from those pieces, so nothing
# past the assembly changed shape.
def period_task_keys(cfg, n_surr, chunk=None):
    n_mc = M.period_n_mc(n_surr)
    keys = ["period_obs"]
    for nt in M.PERIOD_NULLS:
        for a, _b in M.period_chunks(n_mc, chunk):
            keys.append(f"psurr:{nt}:{a}")
    return keys


def period_obs_task(days, resid, cfg):
    power, peaks = M.period_observed(days, resid, period_grid(cfg),
                                     n_peaks=int(cfg["n_peaks"]))
    return {"power": [float(v) for v in power], "peaks": [int(i) for i in peaks],
            "ar1_phi": float(M.ar1_params(resid)[0]),
            "n_trial_periods": int(cfg["n_periods"])}


def period_surrogate_task(days, resid, cfg, null_type, start, stop, seed):
    seqs = M.period_seed_seqs(int(seed))
    vals = M.period_surrogate_maxima(days, resid, period_grid(cfg), null_type,
                                     start, stop, seqs[null_type])
    return {"null_type": null_type, "start": int(start), "stop": int(stop),
            "maxima": [float(v) for v in vals]}


def collect_period_maxima(results, cfg, n_surr, chunk=None):
    """Reassemble the per-null surrogate maxima IN SURROGATE ORDER, counting them.

    FAILURE-FIRST. A missing chunk, a short chunk or a chunk that returned the wrong
    slice is named here -- null type and chunk start -- and stops the run. It cannot
    be allowed to become a quietly shifted p-value: every peak's p is a rank against
    exactly n_mc draws, so losing one chunk moves every p in the scan.
    """
    n_mc = M.period_n_mc(n_surr)
    out = {}
    for nt in M.PERIOD_NULLS:
        parts = []
        for a, b in M.period_chunks(n_mc, chunk):
            key = f"psurr:{nt}:{a}"
            r = results.get(key)
            if r is None:
                raise RuntimeError(
                    f"period scan: chunk {key} (null {nt!r}, surrogates "
                    f"[{a}, {b})) has no result. Refusing to price the scan.")
            got = len(r["maxima"])
            if int(r["start"]) != a or int(r["stop"]) != b or got != b - a:
                raise RuntimeError(
                    f"period scan: chunk {key} (null {nt!r}) claims "
                    f"[{r['start']}, {r['stop']}) with {got} maxima but was "
                    f"declared [{a}, {b}). Refusing to price the scan.")
            parts.append(np.asarray(r["maxima"], dtype=np.float64))
        out[nt] = np.concatenate(parts) if parts else np.empty(0)
    return out, n_mc


def period_finish(counts, offset, days, cfg, obs, maxima, n_mc):
    """Parent-side assembly: p-values, harmonic ladder, folded amplitudes."""
    periods = period_grid(cfg)
    # The GPD stream for the period scan is derived from the test key like every
    # other stream (feature "period_scan", kind "period"), so the extrapolation is
    # reproducible and independent of how the surrogate chunks were scheduled.
    peaks, meta = M.period_assemble(obs["power"], periods, obs["peaks"], maxima,
                                    obs["ar1_phi"], n_mc,
                                    gpd=cfg.get("gpd"),
                                    gpd_rng=task_rng(cfg["seed"], "period_scan",
                                                     "gpd"))
    return _period_decorate(counts, offset, days, cfg, peaks, meta)


def _period_decorate(counts, offset, days, cfg, peaks, meta):
    lam0 = offset * (counts.sum() / offset.sum())
    for pk in peaks:
        pk["ladder"] = M.harmonic_ladder(counts, lam0, days, pk["period_days"],
                                         max_period=counts.size / 3.0)
        pk["feature"] = f"period_{pk['period_days']:.4g}d"
        pk["family"] = 0
        pk["test"] = "lomb_scargle_peak"
        wp = pk["ladder"]["winning_period_days"]
        ph = np.mod(days, wp) / wp
        b = np.minimum((ph * 8).astype(int), 7)
        Y = np.bincount(b, weights=counts, minlength=8)
        L = np.bincount(b, weights=lam0, minlength=8)
        ratio = Y / np.maximum(L * (Y.sum() / L.sum()), 1e-9)
        pk["pct_rate_modulation"] = float(100.0 * (ratio.max() - ratio.min()) / 2.0)
        pk["fold_ratio_by_phase_bin"] = [round(float(v), 4) for v in ratio]
    return {"peaks": peaks, "scan": meta,
            "n_trial_periods": int(cfg["n_periods"]),
            "period_range_days": [PERIOD_MIN, PERIOD_MAX]}


def period_task(counts, offset, days, resid, cfg, n_surr, seeds, verbose=False,
                chunk=None):
    """The whole period scan in ONE process. The unchunked reference.

    Not used by `run()` (which runs the pieces as separate tasks) but kept as the
    thing the chunked path is checked against: same seeds -> byte-identical peaks.
    """
    peaks, meta = M.period_scan(days, resid, period_grid(cfg), n_surr, seeds,
                                n_peaks=int(cfg["n_peaks"]), verbose=verbose,
                                chunk=chunk)
    return _period_decorate(counts, offset, days, cfg, peaks, meta)


# ------------------------------------------------------------------ workers ---
# Windows uses the `spawn` start method, so everything below must be importable at
# module level and every argument must be picklable. The bulk payload (counts,
# offset, features, marks) is shipped ONCE per worker through the pool initializer;
# per-task messages carry only a string key.
_WORKER_PAYLOAD = {}


def _worker_init(payload):
    global _WORKER_PAYLOAD
    _WORKER_PAYLOAD = payload


def dispatch_task(key, W):
    """Run one named task against a payload dict. Pure: no file writes, ever."""
    kind, _, name = key.partition(":")
    seed, n_surr, ladder = W["seed"], W["n_surr"], W["ladder"]
    gpd = W.get("gpd")
    if kind == "glm":
        f = W["feats"][name]
        return glm_task(f, W["window"], W["counts"], W["offset"], n_surr,
                        lambda lag: task_rng(seed, name, "glm", lag,
                                             null_type="circular_shift"),
                        ladder, seed=seed, gpd=gpd)
    if kind == "marks":
        f = W["feats"][name]
        return marks_task(f, W["window"], W["fe_day"], W["marks"], n_surr,
                          task_rng(seed, name, "marks", None,
                                   null_type="circular_shift"),
                          ladder, seed=seed, gpd=gpd)
    if kind == "moment2":
        f = W["feats"][name]
        return moment2_task(f, W["window"], W["counts"], W["offset"], n_surr,
                            task_rng(seed, name, "moment2", None,
                                     null_type="circular_shift"), seed=seed)
    if kind == "omnibus":
        f = W["feats"][name]
        return omnibus_task(f, W["window"], W["counts"], W["offset"], n_surr,
                            task_rng(seed, name, "omnibus", None,
                                     null_type="circular_shift"), seed=seed)
    if kind == "markx":
        f = W["feats"][name]
        tb = W.get("tranche_b") or {}
        return marksx_task(f, W["window"], W["fe_day"], W["marks"],
                           W["ext_marks"], n_surr,
                           task_rng(seed, name, "markx", None,
                                    null_type="circular_shift"), seed=seed,
                           subdaily_values=(W.get("subdaily_values") or {}).get(name),
                           subdaily=bool(tb.get("subdaily")))
    if kind == "regsum":
        f = W["feats"][name]
        return regsum_task(f, W["window"], W["reg_counts"], W["reg_offset"], n_surr,
                           lambda lag: task_rng(seed, name, "regsum", lag,
                                                null_type="circular_shift",
                                                region="ALL"),
                           W["partition"], W["regions_cfg"], seed=seed)
    if kind == "region":
        r_txt, _, fname = name.partition(":")
        r = int(r_txt)
        f = W["feats"][fname]
        return region_battery_task(
            f, W["window"], W["reg_counts"], W["reg_offset"], r, n_surr,
            lambda lag: task_rng(seed, fname, "region", lag,
                                 null_type="circular_shift", region=r),
            W["partition"], W["regions_cfg"], seed=seed)
    if key == "period_obs":
        return period_obs_task(W["days"], W["resid"], W["cfg"])
    if kind == "psurr":
        null_type, _, start = name.partition(":")
        a = int(start)
        b = min(a + int(W.get("period_chunk") or M.PERIOD_SURROGATE_CHUNK),
                M.period_n_mc(n_surr))
        return period_surrogate_task(W["days"], W["resid"], W["cfg"], null_type,
                                     a, b, seed)
    raise KeyError(f"unknown mine task key {key!r}")


def worker_run_task(key):
    """Pool entry point. Returns (key, result, error_text) -- never raises upward.

    A worker that dies outright (OOM, segfault, taskkill) cannot return anything;
    that case is detected by the PARENT via BrokenProcessPool, which is why the
    parent keeps the future->key map. This function covers the other case: an
    ordinary Python exception, whose traceback must reach the parent still carrying
    the feature/test name, because "one of 90 tasks failed" is not a bug report.
    """
    t = time.perf_counter()
    try:
        return key, dispatch_task(key, _WORKER_PAYLOAD), None, time.perf_counter() - t
    except BaseException:                                    # noqa: BLE001
        return key, None, traceback.format_exc(), time.perf_counter() - t


def _run_tasks_parallel(ckpt, keys, payload, jobs, verbose=True):
    """Run the not-yet-checkpointed `keys` across a spawn pool. Parent writes."""
    todo = [k for k in keys if not ckpt.done(k)]
    if verbose:
        print(f"parallel: {len(todo)} tasks over {jobs} worker processes "
              f"({len(keys) - len(todo)} already checkpointed)")
    if not todo:
        return 0
    # One BLAS thread per worker: 30 processes each spawning 32 OpenMP threads is
    # slower than sequential. Set in the PARENT because spawned children inherit
    # os.environ at creation, and by then it is too late to set it from inside.
    saved = {}
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS"):
        saved[var] = os.environ.get(var)
        os.environ[var] = "1"
    n_done = 0
    try:
        ctx = _mp.get_context("spawn")
        with _cf.ProcessPoolExecutor(max_workers=int(jobs), mp_context=ctx,
                                     initializer=_worker_init,
                                     initargs=(payload,)) as ex:
            futs = {ex.submit(worker_run_task, k): k for k in todo}
            try:
                for fut in _cf.as_completed(futs):
                    key = futs[fut]
                    try:
                        rkey, result, err, secs = fut.result()
                    except BrokenProcessPool as exc:
                        pending = sorted(futs[f] for f in futs if not f.done())
                        raise RuntimeError(
                            f"mine worker pool died while running {key!r}; "
                            f"{len(pending)} task(s) unfinished: "
                            f"{', '.join(pending[:10])}"
                            f"{' ...' if len(pending) > 10 else ''}. "
                            f"Completed tasks are checkpointed -- rerun to resume."
                        ) from exc
                    if err is not None:
                        raise RuntimeError(
                            f"mine task {key!r} raised in a worker process:\n{err}")
                    # PARENT-ONLY WRITE. A kill between two puts loses at most the
                    # one task in flight; the ledger is appended once, later, by the
                    # parent, so a partial run cannot double-count multiplicity.
                    ckpt.record_seconds(rkey, secs)
                    ckpt.put(rkey, result)
                    n_done += 1
                    if verbose:
                        print(f"  [{n_done}/{len(todo)}] done {rkey} "
                              f"({secs:.1f}s in-worker)")
            except BaseException:
                for f in futs:
                    f.cancel()
                raise
    finally:
        for var, val in saved.items():
            if val is None:
                os.environ.pop(var, None)
            else:
                os.environ[var] = val
    return n_done


# --------------------------------------------- §P7-8(c)(5) BUILD INVARIANT ---
# "Bitwise identity between two runs that differ in a declared parameter is a
#  build invariant worth asserting -- add it to the run harness: two sessions
#  whose configs differ must not produce identical artifact hashes."
#
# The incident: `session_20260812T021707` (--mag-target 4.0) and
# `session_20260812T004857` (--mag-target 4.5) are bitwise identical, because
# engine/datasets.py:CATALOG_MAG_FLOOR silently clamped the 4.0 request. An
# overnight run was believed to be an independent replicate and was not; had the
# two ever been cited as agreeing, the agreement would have been vacuous.
#
# This check WARNS. It does not delete, quarantine, rename or refuse anything: a
# collision is evidence about the build, and destroying either artifact would
# destroy the evidence. The clamp that caused this particular collision is now
# refused outright (engine/datasets.assert_mag_supported); this invariant is the
# backstop for the NEXT parameter that turns out not to bind.
ARTIFACT_REGISTRY = "artifact_hashes.jsonl"


def artifact_content_hash(tests):
    """sha256 over a session's RESULTS, canonically serialised.

    Hashes the test rows only -- not timings, not the session id, not the
    creation timestamp -- so the hash answers exactly one question: "did these
    two runs compute the same numbers?".
    """
    blob = json.dumps(tests, sort_keys=True, default=repr, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def check_build_invariant(session_dir, config_hash, tests, root=None,
                          register=True, n_tests=None):
    """Assert that configs which differ produce artifacts which differ.

    Returns a dict with `ok`, `artifact_hash`, `collisions` and (when violated) a
    ready-to-print `message`. A prior session with the SAME config hash is not a
    collision -- that is a reproduction, which is the desired behaviour.
    """
    root = root or M.MINE_DIR
    path = os.path.join(root, ARTIFACT_REGISTRY)
    art = artifact_content_hash(tests)
    me = os.path.basename(os.path.normpath(session_dir))

    prior = []
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    prior.append(json.loads(line))
                except json.JSONDecodeError:
                    continue

    collisions = [r for r in prior
                  if r.get("artifact_hash") == art
                  and r.get("config_hash") != config_hash
                  and r.get("session") != me]

    rec = {
        "session": me,
        "session_dir": str(session_dir).replace("\\", "/"),
        "config_hash": config_hash,
        "artifact_hash": art,
        "n_tests": (len(tests) if n_tests is None else int(n_tests)),
        "recorded": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "collides_with": [r["session"] for r in collisions],
    }
    already = any(r.get("session") == me and r.get("artifact_hash") == art
                  for r in prior)
    if register and not already:
        os.makedirs(root, exist_ok=True)
        splits._append_jsonl(path, rec)

    msg = None
    if collisions:
        lines = [
            "!" * 78,
            "BUILD INVARIANT VIOLATED (HYPOTHESIS_LEDGER.md §P7-8(c)(5)):",
            "  two sessions whose CONFIGS DIFFER produced IDENTICAL artifacts.",
            f"  this session : {me}  config_hash={config_hash}",
        ]
        for r in collisions:
            lines.append(f"  collides with: {r.get('session')}  "
                         f"config_hash={r.get('config_hash')}  "
                         f"({r.get('session_dir')})")
        lines += [
            f"  shared artifact content hash: {art[:16]}...",
            "  MEANING: some declared parameter did not bind. These runs are NOT",
            "  independent replicates and must never be cited as agreeing -- the",
            "  agreement would be vacuous. Their EXPLORE_COUNT.jsonl lines declare",
            "  the same tests more than once; per §P7-8(c)(1) those lines are NOT",
            "  edited or reduced, and the over-count stands (the conservative",
            "  direction). Find the parameter that did not bind before quoting",
            "  either session.",
            "  Nothing has been deleted. This is a warning, by rule.",
            "!" * 78,
        ]
        msg = "\n".join(lines)

    return {"ok": not collisions, "artifact_hash": art, "collisions": collisions,
            "registry": path.replace("\\", "/"), "record": rec, "message": msg}


# ------------------------------------------------------------ checkpointing ---
def _cfg_hash(cfg):
    return splits.config_hash(cfg)[:12]


def find_resumable(cfg_hash, root=M.MINE_DIR):
    """The newest incomplete session directory for this exact configuration."""
    if not os.path.isdir(root):
        return None
    best = None
    for name in sorted(os.listdir(root)):
        path = os.path.join(root, name, "checkpoint.json")
        if not (name.startswith("session_") and os.path.exists(path)):
            continue
        try:
            with open(path, "r", encoding="utf-8") as fh:
                ck = json.load(fh)
        except (OSError, json.JSONDecodeError):
            continue
        if ck.get("config_hash") == cfg_hash and not ck.get("complete"):
            best = (os.path.join(root, name), ck)
    return best


class Checkpoint:
    def __init__(self, path, cfg, cfg_hash):
        self.path = path
        self.state = {
            "engine_version": __version__, "kind": "mine",
            "config": cfg, "config_hash": cfg_hash,
            "created": _dt.datetime.now(_dt.timezone.utc).isoformat(),
            "complete": False, "completed_tasks": [], "results": {},
            "banner": M.GENERATOR_NOT_EVIDENCE,
        }

    @classmethod
    def load(cls, path):
        with open(path, "r", encoding="utf-8") as fh:
            state = json.load(fh)
        obj = cls.__new__(cls)
        obj.path = path
        obj.state = state
        return obj

    def done(self, key):
        return key in self.state["results"]

    def record_seconds(self, key, seconds):
        """Per-task cost, so a slow sweep can be diagnosed without re-running it.

        Recorded for BOTH drivers and always measured around the task body itself,
        so parallel numbers are per-task compute and not wall-clock-with-queueing.
        """
        self.state.setdefault("task_seconds", {})[key] = round(float(seconds), 3)

    def put(self, key, value):
        self.state["results"][key] = value
        if key not in self.state["completed_tasks"]:
            self.state["completed_tasks"].append(key)
        self.save()

    def save(self):
        os.makedirs(os.path.dirname(self.path) or ".", exist_ok=True)
        tmp = self.path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(self.state, fh, indent=1)
        os.replace(tmp, self.path)


# ------------------------------------------------------------------- driver ---
# ---- declared inputs to the §P7-1(b) floor, all three MEASURED or DECLARED -----
# VIF: F4-58 as re-specified in §P7-1(c), measured over existing session output --
# median over 2-df phase-feature tests with a clearly-null surrogate p. See
# `results_f4_58_vif.json` and `f4_58_vif.py`. NOT the 3.94 §P7-1(a) inferred.
REGION_VIF_DEFAULT = 24.081827389301417
# alpha: the DECLARED operating threshold of the tranche this runs in -- Tranche A
# (§P7-2, ~713 declared tests) under BH at q = 0.10, most conservative rung.
REGION_ALPHA_DEFAULT = 0.10 / 713.0
# target amplitude: §P6-4 Finding B's own reference rate modulation (20%), not a
# number chosen after seeing what R it produces.
REGION_TARGET_AMPLITUDE_DEFAULT = 0.20


def build_config(args, preset):
    cfg = {
        "engine_version": __version__,
        "mode": "mine",
        "preset": preset["label"],
        "n_surrogates": int(preset["n_surrogates"]),
        "n_periods": int(preset["n_periods"]),
        "n_peaks": int(preset["n_peaks"]),
        "lags": list(preset["lags"]),
        # K-089-R tranche 1: the lag grid for the 13 lag-unscanned CYCLIC features.
        # Off by default -- the historical single lag-0 test per cyclic feature.
        "tranche1": bool(getattr(args, "tranche1", False)),
        "tranche1_lags": (list(M.TRANCHE1_LAGS)
                          if getattr(args, "tranche1", False) else []),
        "fdr_q": M.FDR_Q,
        "baseline": "etas",
        # §P7-8(c)(4): refused at config-build time, so a programmatic caller that
        # never touches the CLI still cannot declare a magnitude the catalogue
        # cannot supply. A clamped --mag-target is what made two sessions with
        # different declared configs bitwise identical.
        "mag_target": float(datasets.assert_mag_supported(args.mag_target)),
        "grid": {"dlat": float(args.dlat), "dlon": float(args.dlon)},
        "explore_frac": float(args.explore_frac),
        "data_dir": args.data_dir.replace("\\", "/"),
        "downloads": bool(not args.no_download),
        "seed": int(args.seed),
    }
    # --ladder changes the SAMPLING RULE, so its parameters must be part of the
    # config hash: two runs that stop drawing under different rules are not the
    # same experiment. The keys are inserted ONLY when the ladder is on, so a
    # default run's config -- and therefore its hash, and therefore its resumable
    # sessions -- is byte-identical to the pre-v2 engine.
    #
    # --jobs is NOT here, on purpose. See the module docstring.
    if getattr(args, "ladder", False):
        cfg["ladder"] = dict(M.LADDER_DEFAULTS)
        for k in ("h", "chunk"):
            v = getattr(args, "ladder_" + k, None)
            if v is not None:
                cfg["ladder"][k] = int(v)
        # N_max is fixed HERE, before the run, from the preset alone. Nothing
        # downstream may raise it for an interesting-looking test.
        cfg["ladder"]["n_max"] = int(cfg["n_surrogates"])
    # --gpd changes the ESTIMATOR (a censored p becomes an extrapolated one) and
    # --strata changes the MULTIPLICITY PROCEDURE, so both are hash-affecting and
    # both are OFF by default. Inserted only when on, so a default run's config --
    # and therefore its hash and its resumable sessions -- is byte-identical to the
    # pre-§P6-2/§P6-3 engine.
    if getattr(args, "gpd", False):
        cfg["gpd"] = {
            "enabled": True,
            "threshold_q": gpd_tail.GPD_THRESHOLD_Q,
            "stability_q": list(gpd_tail.GPD_STABILITY_QS),
            "ad_alpha": gpd_tail.GPD_AD_ALPHA,
            "ci_level": gpd_tail.GPD_CI_LEVEL,
            "n_boot": gpd_tail.GPD_N_BOOT,
            "n_ad_boot": gpd_tail.GPD_N_AD_BOOT,
            "decades": gpd_tail.GPD_DECADES,
            "confirm_factor": gpd_tail.GPD_CONFIRM_FACTOR,
            "confirm_max_n": int(getattr(args, "gpd_confirm_max_n", None)
                                 or 500000),
            # §P6-2(7) is BINDING, not advisory: no calibration (or a REJECT one)
            # and the engine refuses to run with extrapolation on. The verdict is
            # part of the config -- and therefore of the hash -- because a
            # re-calibration is a new licence and a new experiment.
            "calibration": gpd_tail.assert_calibrated(
                getattr(args, "gpd_calibration", None)),
        }
    # §P6-4 Rule 4.2/4.5. The partition RULE, the declared VIF, the declared tranche
    # alpha and the declared target amplitude are all hash-affecting: a run under a
    # different partition or a different floor is a different experiment and may not
    # resume an existing session. Inserted ONLY when the regional statistic is on,
    # so a default run's config -- and therefore its hash, and therefore its
    # resumable sessions -- stays byte-identical to phase 2a.
    if getattr(args, "regsum", False) or getattr(args, "regions", False):
        cfg["regions"] = {
            "enabled": True,
            "rule_id": regions_mod.REGION_RULE_ID,
            "n_sectors": int(getattr(args, "region_sectors", None)
                             or regions_mod.N_SECTORS),
            "min_event_fraction": float(
                getattr(args, "region_min_fraction", None)
                if getattr(args, "region_min_fraction", None) is not None
                else regions_mod.MIN_EVENT_FRACTION),
            "battery": bool(getattr(args, "regions", False)),
            "vif": float(getattr(args, "region_vif", None) or REGION_VIF_DEFAULT),
            "alpha": float(getattr(args, "region_alpha", None)
                           or REGION_ALPHA_DEFAULT),
            "target_amplitude": float(getattr(args, "region_target_amplitude", None)
                                      or REGION_TARGET_AMPLITUDE_DEFAULT),
        }
    # F9-20 arm 1 (§P7-2). The control battery CHANGES THE DECLARED FAMILY -- 713
    # priced tests enter the BH vector -- so it is hash-affecting, and inserted only
    # when on, so a default run's config, hash and resumable sessions are unchanged.
    if getattr(args, "controls", False):
        cfg["controls"] = {
            "enabled": True,
            "rule_id": M.F9_20_RULE_ID,
            "n_named": M.F9_20_N_NAMED,
            "n_matched": M.F9_20_N_MATCHED,
            "lags": list(M.F9_20_LAGS),
            "n_declared_tests": int(M.F9_20_N_DECLARED_TESTS),
            "family": M.F9_20_FAMILY,
        }
    # TRANCHE B (§P7-3, §P7-14(d), §P7-15(b)). Three NEW STATISTICS enter the
    # declared family, so the block is hash-affecting and is inserted ONLY when on:
    # a default run's config, hash and resumable sessions stay byte-identical to
    # phase 2a. `subdaily` is separate from `mark_axis` because §P7-3(3) gates ONLY
    # the sub-daily arm -- the fortnightly-and-longer mark arm may proceed inside B
    # without waiting, and collapsing the two flags would either over-gate the arm
    # that is free or under-gate the arm that is not.
    if getattr(args, "tranche_b", False):
        cfg["tranche_b"] = {
            "enabled": True,
            "rule_id": TRANCHE_B_RULE_ID,
            "second_moment": bool(getattr(args, "tb_second_moment", True)),
            "omnibus": bool(getattr(args, "tb_omnibus", True)),
            "mark_axis": bool(getattr(args, "tb_mark_axis", True)),
            "subdaily": bool(getattr(args, "tb_subdaily", False)),
            "marks_declared": list(marks_ext.MARK_NAMES),
            "marks_scored": list(marks_ext.SCORED_MARK_NAMES),
            "marks_deduplicated": list(marks_ext.DEDUPLICATED_MARKS),
            "alpha": floors.ALPHA_TRANCHE_B,
            "vif_mark_fallback": floors.VIF_MARK_FALLBACK,
            "n_declared_tests": (int(getattr(args, "tb_declared", 0) or 0) or None),
        }
    path = getattr(args, "strata", None)
    if path:
        # §P6-3(3)+(5): the partition FILE CONTENT is hash-affecting (its sha256 is
        # in the config), so re-partitioning cannot resume an existing session --
        # it is a new config hash and a new ledger line, which is rule 5 exactly.
        cfg["strata"] = strata_mod.load_partition(path, q=cfg["fdr_q"])
        strata_mod.assert_budget_identity(cfg["strata"]["strata"], cfg["fdr_q"])
    return cfg


def prepare(cfg, verbose=True):
    """Design + ETAS baseline + target series + the full feature list."""
    ctx = design.build_design(
        data_dir=cfg["data_dir"], dlat=cfg["grid"]["dlat"], dlon=cfg["grid"]["dlon"],
        explore_frac=cfg["explore_frac"], verbose=verbose)
    explore, _hold = splits.temporal_split(ctx.n_days, cfg["explore_frac"])
    y = ctx.day_counts(cfg["mag_target"])

    base = bl.EtasV1(verbose=verbose, mag_target=cfg["mag_target"])
    base.fit(ctx, y, explore)
    if verbose:
        for line in base.report():
            print(line)

    burn = int(base.burn_in_days)
    window = slice(burn, explore.stop)
    counts, offset = M.build_target(ctx, base, y, window)
    if verbose:
        print(f"mining window (days) = [{window.start}, {window.stop}) "
              f"= {counts.size} days")
        print(f"  observed events    = {counts.sum():.0f}")
        print(f"  ETAS expectation   = {offset.sum():.1f}")
        print(f"  residual mean/sd   = {(counts-offset).mean():+.4f} / "
              f"{(counts-offset).std():.4f} events/day")

    t0 = _dt.datetime.fromisoformat(str(ctx.meta["t0"]))
    all_marks = M.load_event_marks(ctx, cfg["data_dir"], ctx.meta["mag_floor"])
    in_win = (all_marks["day"] >= window.start) & (all_marks["day"] < window.stop)
    marks = {k: v[in_win] for k, v in all_marks.items()}
    if verbose:
        print(f"  marks in window    = {marks['day'].size} events "
              f"(magnitude + depth)")

    lags = tuple(cfg["lags"])
    if cfg.get("tranche1"):
        feats = M.ephemeris_features(t0, ctx.n_days,
                                     lags=tuple(cfg["tranche1_lags"]),
                                     lag_features="tranche1")
        if verbose:
            n_scan = sum(1 for f in feats if len(f.lags) > 1)
            print(f"{M.TRANCHE1_LABEL}: lag grid "
                  f"{cfg['tranche1_lags'][0]}..{cfg['tranche1_lags'][-1]} d "
                  f"({len(cfg['tranche1_lags'])} lags) on {n_scan} lag-unscanned "
                  f"cyclic features; {len(M.LAG_FREE_PHASE_FEATURES)} lag-free phase "
                  f"features stay at lag 0 (invariance is a theorem, not a test)")
    else:
        feats = M.ephemeris_features(t0, ctx.n_days)
    dl_feats, dl_log = M.download_features(t0, ctx.n_days, lags,
                                           enabled=cfg["downloads"], verbose=verbose)
    feats += dl_feats
    feats += M.catalog_features(ctx, all_marks, lags)
    # F9-20 arm 1, LAST, so the matched controls take the REAL feature list in its
    # declared order as donors and no control can ever become a donor for another.
    ctl_cfg = cfg.get("controls")
    if ctl_cfg and ctl_cfg.get("enabled"):
        ctls = M.negative_control_features(t0, ctx.n_days, feats,
                                           lags=tuple(ctl_cfg["lags"]),
                                           seed=int(cfg["seed"]))
        feats += ctls
        if verbose:
            print(f"{M.F9_20_LABEL}: {len(ctls)} controls "
                  f"({ctl_cfg['n_named']} named mechanism-free + "
                  f"{ctl_cfg['n_matched']} matched surrogates) x "
                  f"{len(ctl_cfg['lags'])} lags = "
                  f"{ctl_cfg['n_declared_tests']} PRICED tests in their own "
                  f"declared stratum (§P7-2(a)); GLM only -- no mark axis, no "
                  f"regional axis")
    if verbose:
        by_fam = {}
        for f in feats:
            by_fam.setdefault(f.family, []).append(f.name)
        for fam in sorted(by_fam):
            print(f"  family {fam}: {len(by_fam[fam])} features -> "
                  f"{', '.join(by_fam[fam])}")
    return ctx, base, y, window, counts, offset, marks, feats, dl_log, t0


def run(cfg, verbose=True, resume=True, session_dir=None, jobs=1,
        ledger_path=None, prepared=None):
    """`prepared` reuses an already-built context; `ledger_path` redirects the ledger.

    `prepared` is the tuple `prepare(cfg)` returns. It exists for the throughput
    benchmark, which must time the TEST LOOP -- the only thing the parallel layer
    changes -- and not the design build and ETAS fit, which are identical work in
    every cell and would otherwise dominate the wall clock and flatten the speedup
    into meaninglessness. Passing a context prepared under a DIFFERENT config would
    silently mismatch the report, so ordinary callers leave it None.

    It exists for ONE reason: the throughput benchmark runs the same sweep four
    times to time it, and four benchmark lines in `engine/EXPLORE_COUNT.jsonl`
    would inflate the reported multiplicity of real science with runs that were
    never about a hypothesis. Benchmarks point this at scratch. Real runs leave it
    None and get the real ledger. Either way, only the PARENT ever appends.
    """
    ledger_path = ledger_path or splits.EXPLORE_COUNT
    os.makedirs(M.MINE_DIR, exist_ok=True)
    ch = _cfg_hash(cfg)
    ckpt = None
    if session_dir is None and resume:
        found = find_resumable(ch)
        if found:
            session_dir, _state = found
            ckpt = Checkpoint.load(os.path.join(session_dir, "checkpoint.json"))
            if verbose:
                print(f"RESUMING session {session_dir} "
                      f"({len(ckpt.state['results'])} tasks already complete)")
    if session_dir is None:
        stamp = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
        session_dir = os.path.join(M.MINE_DIR, f"session_{stamp}")
        os.makedirs(session_dir, exist_ok=True)
    if ckpt is None:
        ckpt = Checkpoint(os.path.join(session_dir, "checkpoint.json"), cfg, ch)
        ckpt.save()

    t_start = time.time()
    (ctx, base, y, window, counts, offset, marks, feats,
     dl_log, t0) = prepared if prepared is not None else prepare(cfg,
                                                                 verbose=verbose)
    ckpt.state["data_log"] = dl_log
    ckpt.state["window"] = [window.start, window.stop]

    # ---- §P6-4 Rule 4.1/4.2/4.5: the region partition -----------------------
    # Built from EXPLORATION-WINDOW data only (`y[:, window]`); the K-080 census
    # list is not reachable from here. The RULE is in the config hash; the realised
    # partition's digest is pinned in the checkpoint so a resumed session cannot
    # silently re-partition under it.
    regpart = reg_counts = reg_offset = None
    reg_cfg = cfg.get("regions")
    if reg_cfg and reg_cfg.get("enabled"):
        regpart = regions_mod.build_regions(
            ctx, y, window,
            n_sectors=int(reg_cfg["n_sectors"]),
            min_event_fraction=float(reg_cfg["min_event_fraction"]))
        prev = ckpt.state.get("regions", {}).get("digest")
        if prev and prev != regpart["digest"]:
            raise RuntimeError(
                f"region partition digest changed on resume ({prev} -> "
                f"{regpart['digest']}). The declared partition is frozen (§P6-4 Rule "
                f"4.5); refusing to continue a session under a different one.")
        reg_counts, reg_offset = regions_mod.regional_series(
            y, base.rate(window), regpart["region_of_cell"], window, regpart["R"])
        floors = regions_mod.region_floor_table(
            regpart, reg_cfg["vif"], reg_cfg["alpha"], reg_cfg["target_amplitude"])
        r_sum, lam = regions_mod.max_R_for_sum(
            reg_cfg["vif"], reg_cfg["alpha"], float(counts.sum()),
            reg_cfg["target_amplitude"])
        r_bat, r_bat_real = regions_mod.max_R_for_per_region_battery(
            reg_cfg["vif"], reg_cfg["alpha"], float(counts.sum()),
            reg_cfg["target_amplitude"])
        state_part = {k: v for k, v in regpart.items() if k != "region_of_cell"}
        state_part["floors"] = floors
        state_part["R_arithmetic"] = {
            "vif": reg_cfg["vif"], "alpha": reg_cfg["alpha"],
            "target_amplitude": reg_cfg["target_amplitude"],
            "N_window": float(counts.sum()),
            "noncentrality_at_target": lam,
            "max_R_for_2Rdf_sum_at_80pct_power": int(r_sum),
            "max_R_for_per_region_battery": int(r_bat),
            "max_R_for_per_region_battery_real": r_bat_real,
            "R_declared": int(regpart["R"]),
            "sum_measurable": bool(regpart["R"] <= r_sum),
        }
        ckpt.state["regions"] = state_part
        if verbose:
            print(f"regions: rule {regpart['rule_id']} -> R = {regpart['R']} "
                  f"({regpart['n_cells_assigned']}/{regpart['n_cells_total']} cells, "
                  f"{regpart['n_events_assigned']:.0f} events assigned, "
                  f"{regpart['n_events_excluded']:.0f} excluded), digest "
                  f"{regpart['digest']}")
            print(f"  §P7-1(b) floor inputs: VIF = {reg_cfg['vif']:.4g} (F4-58), "
                  f"alpha = {reg_cfg['alpha']:.4g}, target amplitude = "
                  f"{reg_cfg['target_amplitude']:.3g}")
            print(f"  2R-df SUM at R = {regpart['R']}: "
                  f"{'MEASURABLE' if regpart['R'] <= r_sum else 'UNMEASURABLE'} "
                  f"(max R with 80% power = {r_sum})")
            print(f"  PER-REGION battery: {floors['n_unmeasurable']}/"
                  f"{len(floors['rows'])} regions UNMEASURABLE per S-15 "
                  f"(max R whose per-region floor resolves the target = {r_bat})")

    # ---- K-089 clause 3, restated by §P5-5: prove FREE numerically BEFORE the scan --
    axis_audit = M.scan_axis_audit(feats, bool(cfg.get("tranche1")))
    inv = M.lag_invariance_audit(feats, window)
    axis_audit["lag_invariance"] = inv
    ckpt.state["axis_audit"] = axis_audit
    if verbose:
        print(f"axis audit: {axis_audit['tranche']}; "
              f"{len(axis_audit['lag_axis_scanned'])} cyclic features lag-scanned, "
              f"{len(axis_audit['lag_axis_provably_free_not_scanned'])} provably "
              f"lag-free, {len(axis_audit['lag_axis_neither_scanned_nor_free'])} "
              f"neither; {axis_audit['n_new_lag_tests_cyclic']} NEW cyclic lag tests")
        bad = [r for r in inv if not r["ok"]]
        print(f"  lag-invariance audit: {len(inv)} phase features checked, "
              f"{len(bad)} declaration mismatches"
              + ("" if not bad else
                 " -> " + "; ".join(f"{r['feature']}: {r['verdict']}" for r in bad)))
    ckpt.state["session_dir"] = session_dir.replace("\\", "/")
    ckpt.save()

    rng = np.random.default_rng(cfg["seed"])
    n_surr = int(cfg["n_surrogates"])
    ladder = cfg.get("ladder") or None
    resid = counts - offset
    days = np.arange(counts.size, dtype=np.float64)
    fe_day = marks["day"] - window.start

    # The declared task list. Order is the v1 order (GLM sweep, mark tests, period
    # scan) because the sequential driver's shared rng stream depends on it -- for
    # the GLM and mark tasks. The period-scan subtasks do NOT touch that stream:
    # each surrogate draws from its own index-addressable child sequence, so their
    # position in this list is presentational only (see M.surrogate_seed_sequence).
    period_keys = period_task_keys(cfg, n_surr)
    # §P6-4 Rule 4.2 costs ONE task per feature (all its lags inside); the optional
    # battery (Rule 4.3/4.4) costs R x that and is priced as such in the banner.
    # F9-20 prices 23 controls x 31 lags = 713 GLM tests AND NOTHING ELSE. A control
    # is therefore excluded from the mark axis and the regional axis by name: a
    # battery that quietly grew a mark test per control would be running 759 tests
    # under a declaration that says 713.
    sci_feats = [f for f in feats if not getattr(f, "control", False)]
    # A SUB-DAILY-ONLY control (observer.obs_utc_hour_phase) is constant across a
    # daily bin BY CONSTRUCTION, so `Feature.design` would raise on its zero-variance
    # column. It must still be DECLARED -- §P7-3(3)'s gate reads the declared feature
    # set, and it is the one control that is live exactly where the sub-daily mark arm
    # is -- so it is carried in `feats` and excluded from the count path here.
    glm_feats = [f for f in feats if not getattr(f, "subdaily_only", False)]
    regsum_keys = [f"regsum:{f.name}" for f in sci_feats] if regpart else []
    region_keys = ([f"region:{r}:{f.name}"
                    for r in range(regpart["R"]) for f in sci_feats]
                   if (regpart and reg_cfg.get("battery")) else [])
    # ---- TRANCHE B's three new kinds (§P7-3, §P7-14(d)) ---------------------
    # Cyclic science features only: F9-01 and F9-04 are statistics of a PHASE
    # DISTRIBUTION, so a feature with no declared cycle has nothing for them to be
    # computed on, and a CONTROL is excluded from the science axes by name exactly
    # as F9-20 excludes it from the mark and regional axes.
    tb_cfg = cfg.get("tranche_b") or {}
    tb_on = bool(tb_cfg.get("enabled"))
    cyc_feats = ([f for f in sci_feats if f.kind == "phase" or f.period_hint]
                 if tb_on else [])
    moment2_keys = ([f"moment2:{f.name}" for f in cyc_feats]
                    if tb_cfg.get("second_moment") else [])
    omnibus_keys = ([f"omnibus:{f.name}" for f in cyc_feats]
                    if tb_cfg.get("omnibus") else [])
    markx_keys = ([f"markx:{f.name}" for f in sci_feats]
                  if tb_cfg.get("mark_axis") else [])
    ext_marks, subdaily_values = None, {}
    if markx_keys:
        ext_marks, ext_marks_audit = marks_ext.build_marks(marks)
        ckpt.state["f9_10_marks"] = ext_marks_audit
        ckpt.state["f9_10_redundancy_audit"] = marks_ext.redundancy_audit(ext_marks)
        if tb_cfg.get("subdaily"):
            # §P7-3(3): the sub-daily arm is GATED on the observer controls, and the
            # gate is checked HERE -- before a single sub-daily surrogate is drawn --
            # because a sub-daily result computed without them is indistinguishable
            # from one computed with them except in what it is allowed to mean.
            ckpt.state["subdaily_gate"] = observer_mod.assert_subdaily_gate(
                [f.name for f in feats])
            subdaily_values, sub_flags = marks_ext.event_time_feature_values(
                t0, marks["day_float"])
            ckpt.state["subdaily_features"] = {
                k: bool(v) for k, v in sub_flags.items()}
        elif verbose:
            print("tranche B: mark axis FORTNIGHTLY-AND-LONGER only "
                  "(sub-daily arm not requested; §P7-3(3) gate not exercised)")

    task_keys = ([f"glm:{f.name}" for f in glm_feats]
                 + [f"marks:{f.name}" for f in sci_feats]
                 + period_keys + regsum_keys + region_keys
                 + moment2_keys + omnibus_keys + markx_keys)

    # Every random stream in the session is addressed by a canonical test key.
    # Two tests sharing a digest would share a stream, which is a silent
    # correlation between "independent" nulls -- so it is checked, not assumed.
    all_keys = []
    for f in glm_feats:
        for lag in f.lags:
            all_keys.append(M.test_key(cfg["seed"], f.name, "glm", lag=lag,
                                       null_type="block_bootstrap"))
        if getattr(f, "control", False):
            continue
        for mk in ("mag", "depth"):
            all_keys.append(M.test_key(cfg["seed"], f.name, f"mark_{mk}",
                                       null_type="block_bootstrap"))
    # The period scan's two nulls are now separate task families with separate
    # streams (phase 1b). No field was added to TEST_KEY_FIELDS, so every OTHER
    # declared digest in this session is byte-identical to phase 1.
    for _nt in M.PERIOD_NULLS:
        all_keys.append(M.test_key(cfg["seed"], "period_scan", "period",
                                   null_type=_nt))
    # The regional streams use the `region` field that phase 1 reserved in
    # TEST_KEY_FIELDS -- no field was added, so every OTHER declared digest in this
    # session is byte-identical to phase 2a. The regsum test is keyed region="ALL"
    # (it consumes every region at once) and each battery test by its integer id,
    # so a battery test and its sum can never share a stream.
    if regpart:
        for f in sci_feats:
            for lag in f.lags:
                all_keys.append(M.test_key(cfg["seed"], f.name, "regsum", lag=lag,
                                           null_type="circular_shift",
                                           region="ALL"))
        if reg_cfg.get("battery"):
            for r in range(regpart["R"]):
                for f in sci_feats:
                    for lag in f.lags:
                        all_keys.append(M.test_key(
                            cfg["seed"], f.name, "region", lag=lag,
                            null_type="circular_shift", region=r))
    for f in cyc_feats:
        if moment2_keys:
            all_keys.append(M.test_key(cfg["seed"], f.name, "moment2",
                                       null_type="block_bootstrap"))
        if omnibus_keys:
            all_keys.append(M.test_key(cfg["seed"], f.name, "omnibus",
                                       null_type="block_bootstrap"))
    if markx_keys:
        for f in sci_feats:
            all_keys.append(M.test_key(cfg["seed"], f.name, "markx",
                                       null_type="block_bootstrap"))
    digests = assert_task_keys_unique(all_keys)
    n_declared_streams = len(all_keys)
    if len(digests) != n_declared_streams:            # belt and braces
        raise RuntimeError(f"test-key digest count {len(digests)} != "
                           f"{n_declared_streams} declared streams")
    if verbose:
        print(f"test-key schema: {len(M.TEST_KEY_FIELDS)} fields "
              f"{M.TEST_KEY_FIELDS}; {n_declared_streams} declared streams, "
              f"{len(digests)} distinct digests (no collisions)")

    n_jobs = resolve_jobs(jobs)
    if verbose:
        print(f"execution: --jobs {jobs} -> {n_jobs} process(es); "
              + ("SEQUENTIAL, shared rng stream for glm/marks (v1-bit-identical "
                 "audit baseline); period scan on per-surrogate streams"
                 if n_jobs == 1 else
                 "PARALLEL, per-task derived SeedSequence (order-independent; "
                 "glm/marks NOT bit-identical to --jobs 1, period scan IS)")
              + (f"; ladder ON {ladder}" if ladder else "; ladder off"))
    # PRE-FLIGHT MEMORY STATEMENT (v2 phase 1c). Each worker builds and holds its
    # OWN cached Lomb-Scargle basis, so the bill is per-worker x --jobs. The
    # per-worker ceiling is enforced inside `M.build_ls_basis`, which can only see
    # one process; the AGGREGATE is only knowable here, so it is printed here,
    # BEFORE the pool is created rather than after the machine starts swapping.
    _basis_bytes = M.ls_basis_footprint_bytes(counts.size, int(cfg["n_periods"]))
    if verbose:
        print(f"cached Lomb-Scargle basis: {counts.size} days x "
              f"{int(cfg['n_periods'])} trial periods = "
              f"{_basis_bytes / 2**20:.1f} MiB per worker process "
              f"(ceiling {M.LS_BASIS_MAX_BYTES / 2**20:.0f} MiB), "
              f"{n_jobs * _basis_bytes / 2**30:.2f} GiB total at {n_jobs} "
              f"process(es)")
    if _basis_bytes > M.LS_BASIS_MAX_BYTES:
        raise MemoryError(
            f"the period grid for this config needs {_basis_bytes / 2**20:.0f} MiB "
            f"of cached Lomb-Scargle basis PER WORKER, over the "
            f"{M.LS_BASIS_MAX_BYTES / 2**20:.0f} MiB ceiling; at --jobs {n_jobs} "
            f"that is {n_jobs * _basis_bytes / 2**30:.1f} GiB. Refusing before the "
            f"pool is created. Cut n_periods, cut --jobs, or raise "
            f"engine.mine.LS_BASIS_MAX_BYTES deliberately.")

    # The worker payload. Built for BOTH drivers: the sequential path runs the
    # period subtasks through the very same `dispatch_task`, which is what makes
    # --jobs 1 and --jobs N bit-identical there.
    payload = {
        "feats": {f.name: f for f in feats}, "window": window,
        "counts": counts, "offset": offset, "days": days, "resid": resid,
        "marks": marks, "fe_day": fe_day, "seed": int(cfg["seed"]),
        "n_surr": n_surr, "ladder": ladder, "cfg": cfg,
        "gpd": cfg.get("gpd"),
        "period_chunk": M.PERIOD_SURROGATE_CHUNK,
        "partition": regpart, "regions_cfg": reg_cfg,
        "reg_counts": reg_counts, "reg_offset": reg_offset,
        "tranche_b": tb_cfg, "ext_marks": ext_marks,
        "subdaily_values": subdaily_values,
    }

    if n_jobs == 1:
        # -------------------------------------------- (a) GLM sweep --------
        for f in glm_feats:
            key = f"glm:{f.name}"
            if ckpt.done(key):
                if verbose:
                    print(f"  [skip, checkpointed] {key}")
                continue
            t_f = time.perf_counter()
            rows = glm_task(f, window, counts, offset, n_surr,
                            lambda lag: rng, ladder, seed=int(cfg["seed"]),
                            gpd=cfg.get("gpd"))
            ckpt.record_seconds(key, time.perf_counter() - t_f)
            ckpt.put(key, rows)
            if verbose:
                best = min(rows, key=lambda r: r["p_raw"])
                print(f"  glm {f.name:<24s} {len(rows):>2d} lags in "
                      f"{time.perf_counter()-t_f:5.1f}s   best p={best['p_raw']:.4g} "
                      f"(lag {best['lag']}, {best['pct_rate_modulation']:+.2f}%/sd)")

        # ------------------------------------------ (c) mark tests ---------
        for f in sci_feats:
            key = f"marks:{f.name}"
            if ckpt.done(key):
                if verbose:
                    print(f"  [skip, checkpointed] {key}")
                continue
            t_f = time.perf_counter()
            rows = marks_task(f, window, fe_day, marks, n_surr, rng, ladder,
                              seed=int(cfg["seed"]), gpd=cfg.get("gpd"))
            ckpt.record_seconds(key, time.perf_counter() - t_f)
            ckpt.put(key, rows)
            if verbose:
                print(f"  marks {f.name:<22s} mag p={rows[0]['p_raw']:.4g} "
                      f"depth p={rows[1]['p_raw']:.4g}")

        # ---------------------------------------- (b) period scan ----------
        # SAME SUBTASKS AS THE PARALLEL PATH, same per-surrogate streams. This is
        # the design-B choice: --jobs 1 and --jobs N give bit-identical period
        # results, at the price of no longer reproducing the pre-1b shared-stream
        # numbers. See the module docstring.
        for key in period_keys:
            if ckpt.done(key):
                if verbose:
                    print(f"  [skip, checkpointed] {key}")
                continue
            t_f = time.perf_counter()
            res = dispatch_task(key, payload)
            ckpt.record_seconds(key, time.perf_counter() - t_f)
            ckpt.put(key, res)
        if verbose:
            print(f"  period scan: {len(period_keys)} subtasks "
                  f"({M.period_n_mc(n_surr)} AR(1) + {M.period_n_mc(n_surr)} "
                  f"permutation surrogates in chunks of "
                  f"{M.PERIOD_SURROGATE_CHUNK})")

        # ------------------------------ (d) regional sum + battery ---------
        # These run through `dispatch_task` on per-task derived streams in BOTH
        # drivers, so --jobs 1 and --jobs N are bit-identical here (the design-B
        # choice already made for the period scan).
        # ------------------ (e) TRANCHE B: moment2 / omnibus / markx -------
        # Through `dispatch_task` on per-task derived streams in BOTH drivers, the
        # design-B choice already made for the period scan and the regional axes:
        # these kinds have no v1 shared-stream baseline to reproduce, so buying
        # scheduling-independent determinism costs nothing here.
        for key in regsum_keys + region_keys + moment2_keys + omnibus_keys + markx_keys:
            if ckpt.done(key):
                if verbose:
                    print(f"  [skip, checkpointed] {key}")
                continue
            t_f = time.perf_counter()
            res = dispatch_task(key, payload)
            ckpt.record_seconds(key, time.perf_counter() - t_f)
            ckpt.put(key, res)
            if verbose:
                best = min(res, key=lambda r: r["p_raw"])
                print(f"  {key:<40s} {len(res):>2d} lag(s) in "
                      f"{time.perf_counter()-t_f:5.1f}s   best p={best['p_raw']:.4g}")
    else:
        _run_tasks_parallel(ckpt, task_keys, payload, n_jobs, verbose=verbose)

    # BULK INVARIANT. Every declared task must have a result before anything is
    # reported. A silently short pool (cancelled future, swallowed error) would
    # otherwise produce a report whose multiplicity denominator is a lie.
    missing = [k for k in task_keys if not ckpt.done(k)]
    if missing:
        raise RuntimeError(
            f"mine session incomplete: {len(missing)} of {len(task_keys)} declared "
            f"tasks have no result ({', '.join(missing[:10])}"
            f"{' ...' if len(missing) > 10 else ''}). Refusing to write a report.")

    # ------------------------------------------ period-scan assembly -------
    # PARENT ONLY, and after the bulk invariant. `collect_period_maxima` re-counts
    # the surrogates chunk by chunk and `M.period_assemble` re-counts the total
    # before any rank is taken: a dropped chunk must be a crash, never a p-value.
    maxima, n_mc = collect_period_maxima(ckpt.state["results"], cfg, n_surr)
    per = period_finish(counts, offset, days, cfg,
                        ckpt.state["results"]["period_obs"], maxima, n_mc)
    per["n_subtasks"] = len(period_keys)
    per["surrogate_chunk"] = int(M.PERIOD_SURROGATE_CHUNK)
    ckpt.state["results"]["period_scan"] = per
    if verbose:
        print(f"period scan assembled from {len(period_keys)} subtasks: "
              f"{n_mc} AR(1) + {n_mc} permutation surrogates counted exactly, "
              f"{len(per['peaks'])} peaks priced")

    # --------------------------------------------------- multiplicity ------
    # ORDER-DETERMINISTIC REDUCTION. Every row carries `order_key`, a function of
    # its test key alone (declared feature ordinal, test kind, lag, mark). The list
    # is sorted by it before ANY aggregation, and every downstream sort uses it as
    # the final tiebreak, so nothing -- not BH, not the ranked table, not the stub
    # list -- can depend on the order in which a pool happened to return results.
    tests = []
    for i, f in enumerate(feats):
        for t in ckpt.state["results"].get(f"glm:{f.name}", []):
            t["order_key"] = [0, i, int(t["lag"]), 0, f.name]
            tests.append(t)
        for j, t in enumerate(ckpt.state["results"].get(f"marks:{f.name}", [])):
            t["order_key"] = [0, i, -1, 1 + j, f.name]
            tests.append(t)
    for j, t in enumerate(ckpt.state["results"]["period_scan"]["peaks"]):
        t["order_key"] = [1, j, -1, 0, str(t["feature"])]
        tests.append(t)
    # Regional rows sort AFTER the phase-2a families and before nothing, on a
    # leading ordinal of their own, so adding them cannot perturb the order -- and
    # therefore the BH tie-breaks -- of any pre-existing row.
    for i, f in enumerate(feats):
        for t in ckpt.state["results"].get(f"regsum:{f.name}", []):
            t["order_key"] = [2, i, int(t["lag"]), 0, f.name]
            tests.append(t)
    if regpart and reg_cfg.get("battery"):
        for r in range(regpart["R"]):
            for i, f in enumerate(feats):
                for t in ckpt.state["results"].get(f"region:{r}:{f.name}", []):
                    t["order_key"] = [3, i, int(t["lag"]), 1 + r, f.name]
                    tests.append(t)
    # TRANCHE B rows sort on leading ordinals 4/5/6, AFTER every pre-existing
    # family, so switching the tranche on cannot perturb the order -- and therefore
    # the BH tie-breaks -- of any row that existed before it.
    for i, f in enumerate(feats):
        for t in ckpt.state["results"].get(f"moment2:{f.name}", []):
            t["order_key"] = [4, i, -1, 0, f.name]
            tests.append(t)
        for j, t in enumerate(ckpt.state["results"].get(f"omnibus:{f.name}", [])):
            t["order_key"] = [5, i, -1, j, f.name]
            tests.append(t)
        for j, t in enumerate(ckpt.state["results"].get(f"markx:{f.name}", [])):
            t["order_key"] = [6, i, -1, j, f.name]
            tests.append(t)
    tests.sort(key=lambda t: t["order_key"])

    # ---- §P7-17: exactly one disposition per executed row --------------------
    # Tagged BEFORE the multiplicity block, because the tag is what decides whether
    # a row is entitled to a rejection at all. A COMPONENT-OF row is half of another
    # row's claim; it is attached to its parent and it never appears standalone --
    # `write_stubs` raises if one reaches `stubs.json`.
    if tb_on:
        cyc_names = [f.name for f in cyc_feats]
        prior_keys = set()
        disp.tag_rows(tests, prior_keys=prior_keys, cyclic_features=cyc_names)
        disp.assert_one_disposition(tests)
        tests, moved = disp.attach_components(tests)
        ckpt.state["dispositions"] = {
            "counts": disp.counts_by_disposition(tests + moved),
            "n_component_rows_attached_to_parents": len(moved),
            "rule": dict(disp.DISPOSITION_RULE),
            "note": ("§P7-17: COMPONENT-OF rows are removed from the standalone "
                     "test list and attached inside their parents' records. They "
                     "are NOT deleted -- they are half of a claim and the claim "
                     "needs them -- they are simply never readable alone."),
        }
        if verbose:
            print("§P7-17 dispositions: "
                  + ", ".join(f"{k} = {v}" for k, v in
                              sorted(ckpt.state["dispositions"]["counts"].items()))
                  + f"; {len(moved)} component row(s) attached to parents")

    n_tests = len(tests)
    # §P6-2(5): every row carries a p_method. Rows produced by a path that predates
    # the labelling (or by a resumed checkpoint written before it) are labelled here
    # rather than left blank, because "no label" is exactly what rule 5 forbids.
    for t in tests:
        # §P6-4 Rule 4.4 + §P7-1(d): an UNRESOLVED row already carries its label and
        # its reason and is INELIGIBLE for rejection. Re-labelling it MC_RESOLVED
        # here would silently put an unmeasurable region back into the BH vector,
        # which is precisely what the rule forbids -- so it is excluded by name.
        if t.get("p_method") == regions_mod.UNRESOLVED:
            t["p_bh"] = float(t["p_raw"])
            t["bh_eligible"] = False
            continue
        if t.get("p_method") not in gpd_tail.P_METHODS:
            t["p_method"] = gpd_tail.P_MC_RESOLVED
            t["p_bh"] = float(t["p_raw"])
            t["bh_eligible"] = True
            t["p_method_reason"] = "labelled at assembly (no extrapolation path)"

    # §P6-3. Strata are a DECLARATION read from the partition file, whose content is
    # hash-affecting; the default is one stratum holding the whole declared family,
    # which reproduces flat BH exactly (engine/tests/test_strata.py pins that).
    part = cfg.get("strata")
    if part:
        strata_decl = [dict(s) for s in part["strata"]]
        stratum_names = [strata_mod.stratum_of(t, part) for t in tests]
        for t, nm in zip(tests, stratum_names):
            t["stratum"] = nm
    else:
        strata_decl = strata_mod.flat_partition(n_tests, M.FDR_Q)
        stratum_names = [strata_mod.UNSTRATIFIED] * n_tests
        for t in tests:
            t["stratum"] = strata_mod.UNSTRATIFIED

    # §P6-2(2)+(6): GPD rows enter BH at their CI-UPPER, never at the point
    # estimate; §P6-2(1): a test whose extrapolation failed its gates is INELIGIBLE
    # for rejection but still counts against its stratum's declared m_s.
    p = np.array([float(t.get("p_bh", t["p_raw"])) for t in tests])
    eligible = np.array([bool(t.get("bh_eligible", True)) for t in tests])
    q, passed, bh_meta = strata_mod.stratified_bh(
        p, stratum_names, strata_decl, M.FDR_Q, eligible=eligible)
    for t, qq, pa in zip(tests, q, passed):
        t["bh_q"] = float(qq)
        t["passes_fdr"] = bool(pa)
    ckpt.state["bh"] = bh_meta
    ckpt.state["p_method_census"] = gpd_tail.census(tests)

    # ---- §P6-6 R3 / §P6-4(4.7) item 3: the winner's-curse LABEL, on EVERY row --
    # Not only on survivors, and not only on the ranked table. At this scale every
    # quoted amplitude is the maximum of a search and is upward-biased by selection;
    # a label that appears only where somebody remembered to put it is not a label.
    for t in tests:
        if t.get("amplitude_label"):
            continue
        if ("amplitude_log_rate" in t or "pct_rate_modulation" in t
                or t.get("test") == "lomb_scargle_peak"):
            t["amplitude_label"] = SELECTION_BIASED

    # ---- Tranche A read-outs. All priced at 0: each reports a property of the
    # run and makes no rejection (§P7-2(a)).
    tranche_a = None
    # ---- S-15(c), on EVERY session (§P7-14(c) + §P7-16) ---------------------
    # Unconditional: the clause is a property of the analysis WINDOW, not of a
    # tranche, so a session that never heard of Tranche B still gets its
    # long-period peaks labelled rather than quietly scored.
    s15c_census = window_clause_census(tests, float(counts.size))
    ckpt.state["s15c_window_clause"] = s15c_census
    if verbose and s15c_census["n_rows_with_a_period"]:
        print(f"S-15(c): {s15c_census['n_unmeasurable_by_window']}/"
              f"{s15c_census['n_rows_with_a_period']} rows with a period are "
              f"UNMEASURABLE-BY-WINDOW (cut at "
              f"{s15c_census['cut_period_days']:.0f} d over a "
              f"{s15c_census['record_days']:.0f} d window; grid UNCHANGED at "
              f"{s15c_census['period_scan_max_days']:.0f} d)")

    if cfg.get("controls", {}).get("enabled"):
        tranche_a = {
            "label": M.F9_20_LABEL,
            "declared": dict(cfg["controls"]),
            "F10_25_control_calibration": control_calibration(tests, strata_decl),
            "R3_F9_17_winners_curse": winners_curse_report(tests),
            "S15_by_stratum": s15_by_stratum(tests, float(counts.sum()),
                                             strata_decl),
        }
        ckpt.state["tranche_a"] = tranche_a
        if verbose:
            cc = tranche_a["F10_25_control_calibration"]
            wc = tranche_a["R3_F9_17_winners_curse"]
            print(f"  F9-20 control arm: {cc['control_arm']['n_survivors']} "
                  f"survivor(s) of {cc['control_arm']['n_executed']} executed "
                  f"({cc['control_arm']['n_declared']} declared), vs "
                  f"q*m = {cc['control_arm']['expected_false_by_chance_q_times_m']:.1f} "
                  f"expected by chance")
            print(f"  real arm:          {cc['real_arm']['n_survivors']} "
                  f"survivor(s) of {cc['real_arm']['n_executed']} executed "
                  f"({cc['real_arm']['n_declared']} declared), vs "
                  f"q*m = {cc['real_arm']['expected_false_by_chance_q_times_m']:.1f}")
            print(f"  R3 winner's curse: top-{wc['selection_size_k']} median "
                  f"amplitude real {wc['real_arm_median_amplitude_top_k']:.4f} vs "
                  f"null {wc['null_arm_median_amplitude_top_k']:.4f} "
                  f"(excess {wc['excess_over_null_top_k']:+.4f}); every amplitude "
                  f"labelled {SELECTION_BIASED}")
            for r in tranche_a["S15_by_stratum"]["rows"]:
                print(f"  S-15 {r['stratum']}: {r['n_measurable']} MEASURABLE / "
                      f"{r['n_unmeasurable']} UNMEASURABLE of "
                      f"{r['n_count_path_rows']} count-path rows at target "
                      f"{tranche_a['S15_by_stratum']['target_amplitude']:.0%} "
                      f"(A_min {r['A_min_min']:.4f}..{r['A_min_max']:.4f}, "
                      f"alpha_s = {r['alpha_s']:.4g}); max obs/floor "
                      f"{r['max_obs_over_floor']:.3f}")

    # Per-test resolution floors, which are HETEROGENEOUS once the ladder is on
    # (and already heterogeneous without it, because the enumeration floor is set
    # by the record length while the Monte Carlo floor is set by the budget).
    floors = [test_floor(t, n_surr) for t in tests]
    for t, fl in zip(tests, floors):
        t["p_floor"] = float(fl)
    bh_thresh = M.FDR_Q / max(n_tests, 1)
    n_above = sum(1 for fl in floors if fl > bh_thresh)
    enum_floors = sorted({t["floor_enumeration"] for t in tests
                          if t.get("floor_enumeration") is not None})
    mc_floors = sorted({t["floor_monte_carlo"] for t in tests
                        if t.get("floor_monte_carlo") is not None})
    if verbose:
        print(f"multiplicity: {n_tests} tests, BH-FDR at q={M.FDR_Q} -> "
              f"{int(passed.sum())} survive")
        print(f"  BH threshold at the DECLARED count ({n_tests} tests): "
              f"p <= {bh_thresh:.3g} for the smallest.")
        print(f"  RESOLUTION: **{n_above} of {n_tests} declared tests have a "
              f"resolution floor ABOVE the BH threshold** and therefore cannot "
              f"attain it at any effect size.")
        print(f"    enumeration floors (exhaustive circular shifts, set by RECORD "
              f"LENGTH, never laddered): "
              + (", ".join(f"{v:.3g}" for v in enum_floors[:6]) or "none")
              + (" ..." if len(enum_floors) > 6 else ""))
        print(f"    Monte Carlo floors (surrogate budget, laddered when --ladder): "
              + (", ".join(f"{v:.3g}" for v in mc_floors[:6]) or "none")
              + (" ..." if len(mc_floors) > 6 else ""))
        if n_above:
            print(f"  For those {n_above}, S-8's sim-calibrated max-statistic over "
                  f"the declared family -- not BH -- is the confirmatory "
                  f"instrument; the BH line for them is descriptive only.")
        n_at_floor = sum(1 for t, fl in zip(tests, floors)
                         if t["p_raw"] <= fl + 1e-12)
        k_min = int(math.ceil(min(floors) * max(n_tests, 1) / M.FDR_Q))
        if k_min > 1:
            print(f"  WARNING: BH can reject only if >= {k_min} tests tie at the "
                  f"smallest floor {min(floors):.5f} ({n_at_floor} tests sit at "
                  f"their own floor). Survivors are provisional; rerun with more "
                  f"surrogates before trusting the ordering.")

    # ------------------------------------- S-8 max-statistic (§P6-3(4)) ----
    t_ms = time.perf_counter()
    ms_mat, ms_obs, ms_idx = max_statistic_matrix(feats, window, counts, offset,
                                                  tests)
    if ms_mat is not None:
        ms = strata_mod.max_statistic_report(
            ms_mat, ms_obs, [tests[i]["stratum"] for i in ms_idx], strata_decl)
        ms["n_covered"] = len(ms_idx)
        ms["n_declared"] = n_tests
        ms["coverage_note"] = (
            f"{len(ms_idx)} of {n_tests} declared tests carry a SHARED surrogate "
            f"index (the exhaustive circular-shift enumeration) and are therefore "
            f"in the joint max-statistic null. Mark tests resample along the event "
            f"sequence and period peaks resample the residual series; neither "
            f"shares that index, so they are NOT COVERED rather than folded in.")
        ms["seconds"] = round(time.perf_counter() - t_ms, 2)
        ckpt.state["max_statistic"] = ms
        if verbose:
            print(f"  S-8 max-statistic (§P6-3(4)), GLOBAL and "
                  f"partition-invariant: p = {ms['global']['p']:.4g} over "
                  f"{ms['global']['n_tests']} tests x "
                  f"{ms['global']['n_replicates']} shared surrogate worlds "
                  f"(floor {ms['global']['floor']:.3g}); {ms['seconds']}s")
            for r in ms["per_stratum"]:
                print(f"    stratum {r['stratum']}: p = {r['p']:.4g} "
                      f"({r['n_tests']} covered tests) | GLOBAL p = "
                      f"{r['global_p']:.4g} (printed adjacent, §P6-3(4))")
    else:
        ckpt.state["max_statistic"] = None

    # --------------------------- §P6-2(6) GPD candidate confirmation ------
    ckpt.state["gpd_confirmations"] = confirm_gpd_candidates(
        tests, feats, window, counts, offset, marks, fe_day, cfg, verbose=verbose)

    if verbose:
        cen = gpd_tail.census(tests)
        print("  " + gpd_tail.census_line(cen, n_tests))

    # ------------------------------------------------- aliasing audit ------
    fmap = {f.name: f for f in feats}
    for t in tests:
        if not t["passes_fdr"]:
            continue
        key = ("audit:" + t["test"] + ":" + t["feature"] + ":"
               + str(t.get("lag", t.get("mark", ""))))
        if ckpt.done(key):
            t["aliasing"] = ckpt.state["results"][key]
            continue
        # Sequential mode keeps threading the shared stream (v1 behaviour). Parallel
        # mode derives a stream from the audit's own key, so the audit result does
        # not depend on how many survivors preceded it.
        arng = rng if n_jobs == 1 else task_rng(cfg["seed"], key, "aliasing", None)
        if t["test"] == "glm_poisson_offset_etas":
            aud = M.aliasing_audit_glm(fmap[t["feature"]], t["lag"], window,
                                       counts, offset, n_surr, arng)
        elif t["test"] == "lomb_scargle_peak":
            aud = M.aliasing_audit_period(t["ladder"]["winning_period_days"], counts,
                                          offset, marks["day_float"], n_surr, arng,
                                          float(window.start))
        else:
            aud = {"verdict": "N/A (mark test: no time lattice claim)"}
        ckpt.put(key, aud)
        t["aliasing"] = aud
        if verbose:
            print(f"  aliasing audit {t['feature']} ({t['test']}): {aud['verdict']}")

    # ------------------------------------------------------- outputs -------
    ckpt.state["tests"] = tests
    ckpt.state["n_tests"] = n_tests
    ckpt.state["elapsed_seconds"] = round(time.time() - t_start, 1)
    ckpt.state["complete"] = True

    # §P7-8(c)(5) build invariant, checked at session end and BEFORE the report is
    # written so a violation is inside the artifact, not only in the console.
    # The registry lives beside the sessions it indexes (engine/out/mine in
    # production; the tmp dir under test), so a test session never writes into the
    # real registry and a redirected run keeps its own.
    invariant = check_build_invariant(
        session_dir, ch, tests, n_tests=n_tests,
        root=(os.path.dirname(os.path.normpath(session_dir)) or M.MINE_DIR))
    ckpt.state["build_invariant"] = {
        k: v for k, v in invariant.items() if k != "collisions"}
    ckpt.state["build_invariant"]["collides_with"] = [
        {"session": c.get("session"), "config_hash": c.get("config_hash")}
        for c in invariant["collisions"]]
    ckpt.save()
    if invariant["message"]:
        print(invariant["message"], flush=True)
    elif verbose:
        print(f"build invariant      -> OK (artifact "
              f"{invariant['artifact_hash'][:12]}..., no differing-config "
              f"collision in {invariant['registry']})")

    rep = write_report(session_dir, cfg, ckpt.state, tests, feats, counts, offset,
                       marks, base, window, build_invariant=invariant)
    stubs = write_stubs(session_dir, cfg, tests)

    # One ledger line per SESSION, not per invocation: a resumed session is the same
    # sweep continued, so counting it twice would inflate the reported multiplicity.
    if not ckpt.state.get("ledger_logged"):
        splits._append_jsonl(ledger_path,
                             ledger_record(cfg, n_tests, session_dir))
        ckpt.state["ledger_logged"] = True
        ckpt.save()
        if verbose:
            print(f"session ledger       -> {ledger_path} "
                  f"(kind=mine, n_tests={n_tests}, "
                  f"n_declared_tests={n_tests})")
    elif verbose:
        print(f"session ledger       -> already logged for this session "
              f"(resumed run; multiplicity is not double-counted)")
    return {"session_dir": session_dir, "report": rep, "stubs": stubs,
            "n_tests": n_tests, "n_pass": int(sum(t["passes_fdr"] for t in tests)),
            "elapsed": ckpt.state["elapsed_seconds"]}


# =================== TRANCHE A: F10-25, R3/F9-17, S-15 per stratum ==========
# All three are PRICED AT 0 (§P7-2(a)): each reports a property of the run rather
# than making a rejection. They are computed here, in the parent, from rows that
# already exist -- no new surrogate, no new test, no new declaration.

SELECTION_BIASED = "SELECTION-BIASED UPPER BOUND"
# §P6-4 Finding B's own reference rate modulation. Declared BEFORE the run, not
# chosen after seeing which floors it makes measurable.
S15_TARGET_AMPLITUDE = 0.20
# The programme's own standing bound on sinusoidal rate modulation, carried as a
# SECOND and stricter S-15 reference so the headline cannot be read as "measurable"
# at a target this programme would actually care about.
S15_TARGET_STANDING_BOUND = 0.056


def is_control_row(t):
    """A row belonging to the F9-20 negative-control arm."""
    return int(t.get("family", -1)) == M.F9_20_FAMILY


def measured_vif_table(path=None, session="session_20260812T004857"):
    """Per-feature MEASURED VIF (F4-58), §P7-12(a)(1)'s primary.

    Returns `{feature: {...}}` for the named session, or `{}` if the measurement
    file is unreadable. §P7-12(a)(3) reinstates the pooled flat 24.08 as the
    acceptable fallback wherever a per-feature measurement does not exist, on its
    measured <=8% accuracy across the whole 30-800 d block ladder -- so a missing
    file DEGRADES the floor's precision and never blocks the run.
    """
    p = path or os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), floors.RESULTS_JSON)
    try:
        with open(p, "r", encoding="utf-8") as fh:
            doc = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return {}
    for s in doc.get("sessions", []):
        if s.get("session") != session:
            continue
        out = {}
        for r in s.get("per_feature", []):
            v = r.get("vif_median")
            if v:
                out[str(r["feature"])] = {
                    "vif": float(v), "block_days": r.get("block_days"),
                    "df": r.get("df"), "n_usable": r.get("n_usable")}
        return out
    return {}


def s15_by_stratum(tests, n_events, strata_decl, vif_table=None,
                   target=S15_TARGET_AMPLITUDE):
    """§P6-4(4.7) item 2: the measurable/unmeasurable count and fraction, PER STRATUM.

    The floor is §P7-1(b) as re-measured by §P7-8(a) and finalised by §P7-12:

        A_min = sqrt(VIF) * (z_alpha + z_0.80) * sqrt(2/N)

    with the FEATURE'S OWN measured VIF where F4-58 has one (§P7-12(a)(1)) and the
    pooled flat 24.08 where it does not (§P7-12(a)(3)), and with `alpha` the DECLARED
    operating threshold OF THE STRATUM THE ROW RUNS IN -- q_s / m_s, the most
    conservative BH rung inside that stratum. That last point is §P7-1(b)'s whole
    reason for being a formula instead of a scalar: the multiplicity term is a
    property of the declared count, so two strata of different size do not share a
    floor and a scalar carried between them is wrong in both.
    """
    vif_table = measured_vif_table() if vif_table is None else vif_table
    rows = []
    for s in strata_decl:
        nm = s["name"]
        alpha = float(s["q_s"]) / max(int(s["m_s"]), 1)
        sel = [t for t in tests if t.get("stratum") == nm
               and t.get("test") == "glm_poisson_offset_etas"]
        meas = unmeas = 0
        best_ratio, best_row, fl_all = 0.0, None, []
        for t in sel:
            rec = vif_table.get(t["feature"])
            vif = float(rec["vif"]) if rec else floors.MEASURED_VIF_DF2_PHASE
            fl = floors.a_min(vif, alpha, n_events)
            obs = math.expm1(float(t.get("amplitude_log_rate", 0.0)))
            t["a_min_formula"] = fl
            t["a_min_vif"] = vif
            t["a_min_vif_source"] = ("F4-58 per-feature measurement (§P7-12(a)(1))"
                                     if rec else
                                     "pooled flat 24.08 fallback (§P7-12(a)(3))")
            t["a_min_alpha"] = alpha
            t["obs_amplitude_fraction"] = obs
            t["obs_over_floor"] = (obs / fl) if fl > 0 else float("inf")
            ok = fl <= float(target)
            t["s15"] = "MEASURABLE" if ok else "UNMEASURABLE"
            meas += int(ok)
            unmeas += int(not ok)
            fl_all.append(fl)
            if t["obs_over_floor"] > best_ratio:
                best_ratio, best_row = t["obs_over_floor"], t
        n = meas + unmeas
        rows.append({
            "stratum": nm, "m_s": int(s["m_s"]), "q_s": float(s["q_s"]),
            "alpha_s": alpha, "n_count_path_rows": n,
            "n_measurable": meas, "n_unmeasurable": unmeas,
            "fraction_unmeasurable": (unmeas / n) if n else float("nan"),
            "A_min_min": (min(fl_all) if fl_all else None),
            "A_min_max": (max(fl_all) if fl_all else None),
            "n_measurable_at_standing_bound": sum(
                1 for f in fl_all if f <= S15_TARGET_STANDING_BOUND),
            "max_obs_over_floor": best_ratio,
            "max_obs_over_floor_row": (None if best_row is None else
                                       "%s lag %s" % (best_row["feature"],
                                                      best_row.get("lag"))),
        })
    return {
        "target_amplitude": float(target),
        "second_reference_amplitude": S15_TARGET_STANDING_BOUND,
        "conventions": (
            "S-18: every number here carries its convention. alpha is TWO-SIDED and "
            "is q_s/m_s, the stratum's most conservative BH rung; power is fixed at "
            "80%; A_min is a SINUSOIDAL rate-modulation amplitude as a fraction of "
            "the baseline rate; N is the domain event count in the mining window. "
            "The second reference (5.6%) is the programme's own standing bound, "
            "carried so a MEASURABLE verdict at the declared 20% cannot be read as "
            "measurable at an amplitude this programme would care about."),
        "N_events": float(n_events),
        "vif_rule": ("§P7-12(a): the feature's OWN measured VIF where F4-58 has one "
                     "(a direct read, not a fit); the pooled flat 24.08 as the "
                     "acceptable fallback where it does not, accurate to <=8% "
                     "across the 30-800 d block ladder (+0.5% at 30 d, -7.5% at "
                     "800 d). The §P7-10(b) curve is SUPERSEDED and is not used."),
        "n_features_with_own_vif": len(vif_table),
        "rows": rows,
    }


def _decl(by_name, names):
    return int(sum(by_name[n]["m_s"] for n in names if n in by_name))


def control_calibration(tests, strata_decl):
    """F10-25 + F9-20's own read-out: survivors per arm against expectation.

    F10-25 is "the ratio of survivors among real features to survivors among the
    F9-20 negative controls", and the catalog is explicit that it means nothing
    unless the two arms ran in the same tranche, at the same lags, under the same
    nulls -- which is why the controls are in THIS session and THIS vector rather
    than in a run of their own.
    """
    ctl = [t for t in tests if is_control_row(t)]
    real = [t for t in tests if not is_control_row(t)]
    by_name = {s["name"]: s for s in strata_decl}
    ctl_strata = sorted({t.get("stratum") for t in ctl})
    real_strata = sorted({t.get("stratum") for t in real})

    def _exp(names):
        return float(sum(by_name[n]["m_s"] * by_name[n]["q_s"]
                         for n in names if n in by_name))

    n_ctl_s = sum(1 for t in ctl if t["passes_fdr"])
    n_real_s = sum(1 for t in real if t["passes_fdr"])
    return {
        "control_arm": {
            "strata": ctl_strata, "n_declared": _decl(by_name, ctl_strata),
            "n_executed": len(ctl), "n_survivors": n_ctl_s,
            "expected_false_by_chance_q_times_m": _exp(ctl_strata),
            "survivors_are": ("MEASURED FALSE POSITIVES (§P7-2(a)). A survivor in "
                              "the control stratum is never reported as a finding."),
        },
        "real_arm": {
            "strata": real_strata, "n_declared": _decl(by_name, real_strata),
            "n_executed": len(real), "n_survivors": n_real_s,
            "expected_false_by_chance_q_times_m": _exp(real_strata),
        },
        "F10_25_survivor_ratio": ((n_real_s / n_ctl_s) if n_ctl_s else None),
        "F10_25_note": (
            "UNDEFINED when the control arm has zero survivors, and that is the "
            "expected outcome of an honest machine at this scale rather than a "
            "failure of the statistic. The interpretable object at zero survivors is "
            "the PAIR of counts against the pair of q*m expectations, printed above "
            "rather than collapsed into a ratio with a zero denominator."),
        "what_a_control_survivor_would_mean": (
            "a measured false positive under the ACTUAL dependence structure of the "
            "ACTUAL data -- the number BH's theorem assumes rather than measures, "
            "and the reason §P7-2(a) priced this battery instead of exempting it."),
    }


def winners_curse_report(tests, top_k=10):
    """R3 / F9-17 (§P6-4(4.7) item 3, §P6-6 R3). The label, and the null comparison.

    R3 accepts EITHER a selection-debiased estimate OR the `SELECTION-BIASED UPPER
    BOUND` label plus the median effect among survivors of an equal-sized null run.
    This build takes the second branch, and the F9-20 control arm IS the null run:
    null features, real data, same lags, same nulls, same session, same vector.

    THE SIZE MISMATCH IS STATED WITH ITS DIRECTION DERIVED, not asserted (S-18
    clause 4). The control arm declares 713 count-path tests against the real arm's
    500, so the null arm is the LARGER search; a fixed upper order statistic of a
    larger sample from the same null is stochastically larger, so the null arm's
    selected median is if anything an OVER-estimate of the selection inflation a
    size-matched null would show. The comparison therefore UNDER-states the real
    arm's excess and cannot manufacture one.
    """
    def _amp(t):
        return abs(math.expm1(float(t.get("amplitude_log_rate", 0.0))))

    def _p(t):
        return float(t.get("p_bh", t["p_raw"]))

    rows = [t for t in tests if t.get("test") == "glm_poisson_offset_etas"]
    ctl = sorted((t for t in rows if is_control_row(t)),
                 key=lambda t: (_p(t), t.get("order_key", [])))
    real = sorted((t for t in rows if not is_control_row(t)),
                  key=lambda t: (_p(t), t.get("order_key", [])))
    n_real_surv = sum(1 for t in real if t["passes_fdr"])

    # DEDUPLICATE THE LAG AXIS BEFORE SELECTING, and this is not cosmetic. At long
    # period the 31 lags of one feature are near-identical columns: a raw top-10 by
    # p can be ten lags of a SINGLE feature, and a "median amplitude among the top
    # 10" computed on that is one feature's amplitude quoted ten times. The primary
    # comparison therefore takes each feature's own best row and then the top k
    # FEATURES. The row-level figure is kept beside it, labelled, because it is what
    # a naive reader would compute and the difference is worth seeing.
    def _best_per_feature(seq):
        seen, out = set(), []
        for t in seq:                       # already sorted by p
            if t["feature"] in seen:
                continue
            seen.add(t["feature"])
            out.append(t)
        return out

    ctl_f, real_f = _best_per_feature(ctl), _best_per_feature(real)
    k = min(max(int(n_real_surv), int(top_k)), len(ctl_f), len(real_f))
    k_rows = min(max(int(n_real_surv), int(top_k)), len(ctl), len(real))

    def _med(seq):
        v = [_amp(t) for t in seq]
        return float(np.median(v)) if v else float("nan")

    return {
        "rule": "§P6-6 R3 / §P6-4(4.7) item 3, implemented as F9-17.",
        "branch_taken": ("LABEL + NULL-RUN MEDIAN. Every amplitude this session "
                         "reports carries `" + SELECTION_BIASED + "`; the median "
                         "effect among an equal-sized selection from the F9-20 null "
                         "arm is printed beside it."),
        "label_applied_to_every_amplitude": SELECTION_BIASED,
        "n_real_survivors": int(n_real_surv),
        "selection_size_k": int(k),
        "selection_rule": ("each feature's own smallest p_bh, then the k features "
                           "with the smallest of those, k = max(number of real "
                           "survivors, %d). Declared here rather than chosen after "
                           "seeing the numbers, and IDENTICAL in both arms -- an "
                           "unequal selection rule compares a maximum against a "
                           "mean and calls the difference physics." % int(top_k)),
        "real_arm_median_amplitude_top_k": _med(real_f[:k]),
        "null_arm_median_amplitude_top_k": _med(ctl_f[:k]),
        "real_arm_max_amplitude_top_k": max((_amp(t) for t in real_f[:k]),
                                            default=float("nan")),
        "null_arm_max_amplitude_top_k": max((_amp(t) for t in ctl_f[:k]),
                                            default=float("nan")),
        "real_arm_top_k_features": [t["feature"] for t in real_f[:k]],
        "null_arm_top_k_features": [t["feature"] for t in ctl_f[:k]],
        "real_arm_median_amplitude_all": _med(real),
        "null_arm_median_amplitude_all": _med(ctl),
        "excess_over_null_top_k": _med(real_f[:k]) - _med(ctl_f[:k]),
        "row_level_secondary": {
            "note": ("the same statistic WITHOUT lag-axis deduplication -- what a "
                     "naive top-k by p returns. Reported so the deduplication is "
                     "visible as a choice rather than hidden as a default."),
            "k": int(k_rows),
            "real_arm_median_amplitude": _med(real[:k_rows]),
            "null_arm_median_amplitude": _med(ctl[:k_rows]),
            "n_distinct_real_features_in_top_k": len({t["feature"]
                                                      for t in real[:k_rows]}),
            "n_distinct_null_features_in_top_k": len({t["feature"]
                                                      for t in ctl[:k_rows]}),
        },
        "kind_balance": (
            "the 20 matched controls copy their donor's `kind`, so the null arm "
            "carries the same mix of 2-df phase and 1-df linear columns as the real "
            "arm; the two medians are not comparing a phase amplitude against a "
            "per-sd slope."),
        "size_mismatch": {
            "n_real_count_path": len(real), "n_null_count_path": len(ctl),
            "direction_of_bias": (
                "the null arm is the LARGER search (%d vs %d count-path rows), so "
                "its selected median is stochastically the larger; the comparison "
                "therefore UNDER-states the real arm's excess and cannot "
                "manufacture one." % (len(ctl), len(real))),
            "derivation": (
                "a fixed upper order statistic of a larger sample drawn from the "
                "same distribution is stochastically larger; both arms are drawn "
                "under the same null by construction, so n is the only asymmetry "
                "left. Derived, not asserted -- S-18 candidate clause 4."),
        },
        "what_this_does_not_buy": (
            "a debiased estimate. The label is a LABEL: it says the number is an "
            "upper bound, not how much of it is selection. Nothing here licenses "
            "quoting any amplitude as an effect size."),
    }


def rank_stability(ref_tests, other_tests, top_k=100, label=""):
    """R4 (reseeding) and F10-24 (data resampling): top-k rank correlation.

    §P7-2(b) is explicit that these are ADJACENT, BOTH REQUIRED, and must never be
    conflated: R4 asks whether the RNG moves the ranking, F10-24 whether the DATA
    moves it. Same statistic, two different failure modes, so the caller supplies
    the label and it is carried into the output.
    """
    def _key(t):
        return (t["test"], t["feature"], t.get("lag"), t.get("mark"))

    def _rank(ts):
        s = sorted(ts, key=lambda t: (float(t.get("p_bh", t["p_raw"])),
                                      t.get("order_key", [])))
        return {_key(t): i for i, t in enumerate(s)}

    ra, rb = _rank(ref_tests), _rank(other_tests)
    top = [k for k, v in sorted(ra.items(), key=lambda kv: kv[1])[:top_k]]
    common = [k for k in top if k in rb]
    if len(common) < 3:
        return {"label": label, "n_top": top_k, "n_common": len(common),
                "spearman_rho": None,
                "verdict": "NOT COMPUTABLE (fewer than 3 shared rows)"}
    a = np.array([ra[k] for k in common], dtype=float)
    b = np.array([rb[k] for k in common], dtype=float)
    rho = float(stats.spearmanr(a, b).statistic)
    top_b = {k for k, v in sorted(rb.items(), key=lambda kv: kv[1])[:top_k]}
    overlap = len(set(top) & top_b)
    return {
        "label": label, "n_top": top_k, "n_common": len(common),
        "spearman_rho": rho,
        "top_k_set_overlap": overlap,
        "top_k_set_overlap_fraction": overlap / max(len(top), 1),
        "verdict": ("STABLE" if rho >= 0.9 else
                    "PARTIALLY STABLE" if rho >= 0.5 else
                    "UNSTABLE -- the ranking is noise and the banner must say so"),
    }


# ------------------------------------------------------------------ report ---
def _row_effect(t):
    if t["test"] in ("glm_poisson_offset_etas", "glm_poisson_offset_etas_region"):
        if t["kind"] == "phase":
            return f"amp {t['amplitude_log_rate']:.4f} log-rate"
        return f"beta {t['beta'][0]:+.4f} +/- {t['se'][0]:.4f} /sd"
    if t["test"] == "regsum_score_2Rdf":
        # §P7-1(d): the sum has no amplitude and is not given one.
        return (f"chi2 {t['chi2_score']:.2f} on {t['df']} df "
                f"(R = {t['R']}); NO AMPLITUDE -- detection statistic, not estimator")
    if t["test"] == "lomb_scargle_peak":
        return f"LS power {t['power']:.4f}"
    if t["test"] == "second_circular_moment_score":
        # F9-01's Pit: R2 is never quoted without R1 beside it, because the reading
        # is AMBIGUOUS between an axial and a harmonic response and the two numbers
        # together are what disambiguates it.
        return (f"R2 {t['R2']:.4f} (R1 {t['R1']:.4f}); chi2 {t['chi2_score']:.2f} "
                f"on 2 df at the doubled angle")
    if t["test"] in ("kuiper_V", "watson_U2"):
        return f"V* {t['V_star']:.4f} / U2* {t['U2_star']:.4f}"
    return f"{t['test']} {t['effect']:+.4f}"


def _row_where(t):
    if t["test"] == "regsum_score_2Rdf":
        return f"lag {t['lag']} d, all {t['R']} regions summed ({t['region_rule_id']})"
    if t["test"] == "glm_poisson_offset_etas_region":
        return (f"lag {t['lag']} d, region {t['region']} "
                f"[{t['lon_lo']:+.0f}, {t['lon_hi']:+.0f}) deg lon")
    if t["test"] == "glm_poisson_offset_etas":
        return f"lag {t['lag']} d"
    if t["test"] == "lomb_scargle_peak":
        return f"P = {t['period_days']:.4g} d"
    if t["test"] in ("second_circular_moment_score", "kuiper_V", "watson_U2"):
        p = t.get("period_days")
        return ("phase distribution of `%s`%s" % (t["feature"],
                "" if not p else " (P = %.4g d)" % p))
    if t.get("mark_axis") == "F9-10":
        return ("mark %s, %s" % (t["mark"],
                "SUB-DAILY (event times)" if t.get("subdaily") else "day-binned"))
    return f"mark {t.get('mark')}"


def _amp_note(t):
    if t["test"] == "regsum_score_2Rdf":
        return t["amplitude_note"]
    if t["test"] == "lomb_scargle_peak":
        return f"+/-{t['pct_rate_modulation']:.1f}% folded rate"
    if t["test"] in ("glm_poisson_offset_etas", "glm_poisson_offset_etas_region"):
        return f"+/-{t['pct_rate_modulation']:.2f}% rate"
    if t["test"] == "second_circular_moment_score":
        return f"+/-{t['pct_rate_modulation']:.2f}% rate at the SECOND harmonic"
    if t["test"] in ("kuiper_V", "watson_U2"):
        # F9-04 is an omnibus: it has no amplitude parameter at all, and inventing
        # one for a report column is exactly the §P7-1(d) error the regsum row
        # already refuses to make.
        return "OMNIBUS -- no amplitude parameter (detection statistic)"
    if t.get("mark_axis") == "F9-10":
        return (f"rank correlation; floor rho_min = {t['rho_min']:.4f} "
                f"(§P7-10(c), VIF_mark = {t['mark_floor']['vif_mark']:.3f})")
    return "rank correlation (no rate amplitude)"


def write_report(session_dir, cfg, state, tests, feats, counts, offset, marks,
                 base, window, build_invariant=None):
    path = os.path.join(session_dir, "report.md")
    order = sorted(tests, key=lambda t: (t["bh_q"], t["p_raw"], t.get("order_key", [])))
    n_pass = sum(t["passes_fdr"] for t in tests)
    L = []
    A = L.append

    tranche = bool(cfg.get("tranche1"))
    A(f"# mine session report -- {os.path.basename(session_dir)}"
      + (f" -- **{M.TRANCHE1_LABEL}**" if tranche else ""))
    A("")
    if tranche:
        tl = cfg["tranche1_lags"]
        A(f"> **RUN LABEL: {M.TRANCHE1_LABEL}.** Lag grid {tl[0]}..{tl[-1]} d "
          f"({len(tl)} lags) applied to the 13 lag-unscanned cyclic features "
          f"(8 `kind='linear'` + 5 `kind='phase'` whose lag is not a rotation); the "
          f"4 provably lag-free phase features "
          f"({', '.join('`%s`' % n for n in M.LAG_FREE_PHASE_FEATURES)}) stay at "
          f"lag 0 by theorem. See the axis audit below.")
        A("")
    A("> **" + M.GENERATOR_NOT_EVIDENCE.split(".")[0] + ".**")
    A(">")
    for line in _wrap(M.GENERATOR_NOT_EVIDENCE):
        A("> " + line)
    A("")
    A("> **Standing warning (EQ-24, verbatim):** " + M.MIGNAN_BROCCARDO)
    A("")
    A("## Configuration")
    A("")
    A("```json")
    A(json.dumps(cfg, indent=1))
    A("```")
    A("")
    A(f"- engine v{state['engine_version']}, preset `{cfg['preset']}`, "
      f"elapsed {state.get('elapsed_seconds')} s")
    A(f"- baseline: **{base.name}** -- {base.caveat}")
    A(f"- mining window: exploration days [{window.start}, {window.stop}) = "
      f"{counts.size} days (365 d ETAS burn-in dropped)")
    A(f"- target: daily domain-wide counts vs sum of lambda_etas; "
      f"{counts.sum():.0f} observed events, {offset.sum():.1f} expected")
    A(f"- marks: {marks['mag'].size} events (magnitude, depth)")
    if build_invariant is not None:
        # §P7-8(c)(5). Printed on every run, pass or fail: a build invariant that
        # is only visible when it fires is not auditable.
        if build_invariant["ok"]:
            A(f"- build invariant (§P7-8(c)(5)): **OK** -- artifact content hash "
              f"`{build_invariant['artifact_hash'][:16]}`, no session with a "
              f"DIFFERENT config hash shares it "
              f"(registry `{build_invariant['registry']}`)")
        else:
            A("")
            A(f"> ### BUILD INVARIANT VIOLATED (§P7-8(c)(5))")
            A(">")
            A(f"> This session (`{os.path.basename(os.path.normpath(session_dir))}`, "
              f"config hash `{state.get('config_hash')}`) produced an artifact "
              f"IDENTICAL to a session with a different config hash:")
            for c in build_invariant["collisions"]:
                A(f"> - `{c.get('session')}` (config hash `{c.get('config_hash')}`)")
            A(">")
            A(f"> Shared artifact content hash: "
              f"`{build_invariant['artifact_hash'][:16]}`.")
            A(">")
            for line in _wrap(
                "Some declared parameter did not bind. These runs are NOT "
                "independent replicates and must never be cited as agreeing -- the "
                "agreement would be vacuous. Their EXPLORE_COUNT.jsonl lines "
                "declare the same tests more than once; per §P7-8(c)(1) those "
                "lines are not edited or reduced and the over-count stands, which "
                "is the conservative direction. Nothing has been deleted: a "
                "collision is evidence about the build and destroying either "
                "artifact would destroy the evidence."
            ):
                A("> " + line)
            A("")
    A(f"- features: {len(feats)} "
      + ", ".join(f"family {k}: {v}" for k, v in sorted(
          _count_by(feats, lambda f: f.family).items())))
    A(f"- **{state['n_tests']} tests** in this session; BH-FDR at q = {cfg['fdr_q']} "
      f"-> **{n_pass} survive**")
    # ---- §P6-4 rule 4.7 item 1: the declared count AND what chance alone buys ----
    _m_decl = int(state["n_tests"])
    _exp_false = float(cfg["fdr_q"]) * _m_decl
    A(f"- **DECLARED TEST COUNT {_m_decl}, AND WHAT CHANCE ALONE BUYS AT THIS "
      f"THRESHOLD (§P6-4(4.7) item 1).** Operating at q = {cfg['fdr_q']}, the "
      f"expected number of FALSE discoveries among the survivors is "
      f"**q x m = {_exp_false:.1f}**. Read that number before the ranked list: with "
      f"{n_pass} survivor(s) reported, chance alone is expected to supply "
      f"{_exp_false:.1f} of them. BH controls the EXPECTED PROPORTION of false "
      f"discoveries, not their number, and it makes no statement at all about any "
      f"individual row.")
    # ---- §P6-4 rule 4.7 item 4: the p-method census + the resolvability count ----
    _cen = state.get("p_method_census") or gpd_tail.census(tests)
    A(f"- **P-METHOD CENSUS (§P6-2(5), §P6-4(4.7) item 4).** "
      f"{_cen.get(gpd_tail.P_MC_RESOLVED, 0)} `MC_RESOLVED` / "
      f"{_cen.get(gpd_tail.P_GPD, 0)} `GPD_EXTRAPOLATED` / "
      f"{_cen.get(gpd_tail.P_UNRESOLVED, 0)} `UNRESOLVED`, of {_m_decl} declared "
      f"tests. `MC_RESOLVED` means the Monte Carlo actually resolved the p; "
      f"`GPD_EXTRAPOLATED` means it did not and a gated GPD tail fit was quoted at "
      f"the UPPER end of its 95% bootstrap CI; `UNRESOLVED` means the p is the "
      f"floor and the true p is somewhere at or below it. A reader must never have "
      f"to infer which kind of number they are looking at."
      + (f" **{_cen.get(gpd_tail.P_GPD, 0)} extrapolated survivor(s) are labelled "
         f"{gpd_tail.CANDIDATE_LABEL} and emit no stub** (§P6-2(6))."
         if _cen.get(gpd_tail.P_GPD, 0) else ""))
    # ---- S-15(c) UNMEASURABLE-BY-WINDOW (§P7-14(c), §P7-16) ------------------
    _wc = state.get("s15c_window_clause")
    if _wc and _wc["n_rows_with_a_period"]:
        A("")
        A("## S-15(c) UNMEASURABLE-BY-WINDOW -- identifiability, not power")
        A("")
        A(f"*A periodic feature with fewer than "
          f"{floors_mod.MIN_CYCLES_IN_WINDOW:g} full cycles in the analysis window is "
          f"UNMEASURABLE-BY-WINDOW **regardless of N**. Over this "
          f"{_wc['record_days']:.0f} d window the cut falls at "
          f"**period > {_wc['cut_period_days']:.0f} d**. This is an "
          f"IDENTIFIABILITY limit and is orthogonal to the S-15 power floor: a "
          f"feature can pass `a_min` and fail this, and no amount of N repairs it. "
          f"Threshold inherited from `mine.py:harmonic_ladder`'s "
          f"`hi_cap = record/3`. ({floors_mod.WINDOW_CLAUSE_SOURCE})*")
        A("")
        A(f"**The period grid is UNCHANGED** at "
          f"{_wc['period_scan_max_days']:.0f} d (§P7-16). Peaks between "
          f"{_wc['cut_period_days']:.0f} d and "
          f"{_wc['period_scan_max_days']:.0f} d are REPORTED here rather than "
          f"removed from the scan: a scan that silently cannot look somewhere is "
          f"worse than one that looks and says the answer is not identifiable "
          f"there.")
        A("")
        A(f"**{_wc['n_unmeasurable_by_window']} of "
          f"{_wc['n_rows_with_a_period']}** rows carrying a period are "
          f"UNMEASURABLE-BY-WINDOW. They are **scored NEITHER WAY**: not a null, "
          f"not a detection, and removed from both numerators of the S-15 headline "
          f"fraction.")
        if _wc["rows"]:
            A("")
            A("| test | feature | period (d) | cycles in window | p_raw | verdict |")
            A("| --- | --- | ---: | ---: | ---: | --- |")
            for r in _wc["rows"][:40]:
                A(f"| `{r['test']}` | `{r['feature']}` | {r['period_days']:.4g} | "
                  f"{r['cycles_in_window']:.3g} | {r['p_raw']:.4g} | "
                  f"{floors_mod.UNMEASURABLE_BY_WINDOW} |")
            if len(_wc["rows"]) > 40:
                A(f"| ... | *{len(_wc['rows']) - 40} more* | | | | |")

    # ---- TRANCHE A (§P7-2): the control arm, R3, and S-15 per stratum --------
    _ta = state.get("tranche_a")
    if _ta:
        A("")
        A(f"## TRANCHE A -- {_ta['label']}")
        A("")
        for line in _wrap(
            "This session runs TWO ARMS IN ONE DECLARED VECTOR, at the same lags, "
            "under the same nulls, against the same real data, with one config "
            "hash. The REAL arm is the already-declared K-089-R tranche-1 sweep, "
            "re-occupying its slots in this session's BH denominator and presented "
            "as NO NEW SCOPE. The CONTROL arm is F9-20's negative-control battery: "
            "declared certain-zero-effect features, PRICED at 713 tests in their "
            "own declared stratum per HYPOTHESIS_LEDGER.md §P7-2(a), because an "
            "unpriced control is an unaudited channel and because F10-25's ratio "
            "means nothing unless both arms face the same threshold."
        ):
            A(line)
        A("")
        cc = _ta["F10_25_control_calibration"]
        A("### F10-25 -- the survivor ratio, and what chance alone buys per arm")
        A("")
        A("| arm | declared m_s | executed | survivors | expected by chance (q*m) |")
        A("| --- | ---: | ---: | ---: | ---: |")
        for nm, d in (("REAL (already declared)", cc["real_arm"]),
                      ("F9-20 CONTROL (new price)", cc["control_arm"])):
            A(f"| {nm} | {d['n_declared']} | {d['n_executed']} | "
              f"**{d['n_survivors']}** | "
              f"{d['expected_false_by_chance_q_times_m']:.1f} |")
        A("")
        rat = cc["F10_25_survivor_ratio"]
        A(f"- **F10-25 survivor ratio (real : control) = "
          f"{'UNDEFINED' if rat is None else format(rat, '.3g')}**")
        for line in _wrap(cc["F10_25_note"]):
            A("  " + line)
        A("")
        for line in _wrap("A survivor in the control stratum is "
                          + cc["what_a_control_survivor_would_mean"]
                          + " It is a MEASURED FALSE POSITIVE (\u00a7P7-2(a)) and is "
                          + "never reported as a finding."):
            A(line)
        A("")
        wc = _ta["R3_F9_17_winners_curse"]
        A("### R3 / F9-17 -- winner's curse (§P6-4(4.7) item 3)")
        A("")
        for line in _wrap(wc["branch_taken"]):
            A(line)
        A("")
        A(f"| quantity | REAL arm | F9-20 NULL arm |")
        A("| --- | ---: | ---: |")
        A(f"| count-path rows | {wc['size_mismatch']['n_real_count_path']} | "
          f"{wc['size_mismatch']['n_null_count_path']} |")
        A(f"| median amplitude, top-{wc['selection_size_k']} by p | "
          f"{wc['real_arm_median_amplitude_top_k']:.4f} | "
          f"{wc['null_arm_median_amplitude_top_k']:.4f} |")
        A(f"| max amplitude, top-{wc['selection_size_k']} by p | "
          f"{wc['real_arm_max_amplitude_top_k']:.4f} | "
          f"{wc['null_arm_max_amplitude_top_k']:.4f} |")
        A(f"| median amplitude, all rows | "
          f"{wc['real_arm_median_amplitude_all']:.4f} | "
          f"{wc['null_arm_median_amplitude_all']:.4f} |")
        A("")
        A(f"- selection rule: {wc['selection_rule']}")
        A(f"- excess of the real arm over the null arm at top-"
          f"{wc['selection_size_k']}: **{wc['excess_over_null_top_k']:+.4f}** "
          f"(fractional rate modulation)")
        A(f"- size mismatch, direction DERIVED not asserted: "
          f"{wc['size_mismatch']['direction_of_bias']} "
          f"{wc['size_mismatch']['derivation']}")
        A(f"- **WHAT THIS DOES NOT BUY:** {wc['what_this_does_not_buy']}")
        A(f"- {wc['kind_balance']}")
        _rl = wc["row_level_secondary"]
        A(f"- without lag-axis deduplication the same statistic reads real "
          f"{_rl['real_arm_median_amplitude']:.4f} vs null "
          f"{_rl['null_arm_median_amplitude']:.4f} over the top-{_rl['k']} ROWS, "
          f"which span only {_rl['n_distinct_real_features_in_top_k']} real and "
          f"{_rl['n_distinct_null_features_in_top_k']} null distinct features. "
          f"{_rl['note']}")
        A("")
        s15 = _ta["S15_by_stratum"]
        A("### S-15 per stratum (§P6-4(4.7) item 2) under the §P7-12 floors")
        A("")
        A(f"| stratum | m_s | alpha_s | count-path rows | MEASURABLE | "
          f"UNMEASURABLE | A_min range | measurable at {s15['second_reference_amplitude']:.1%} | max obs/floor |")
        A("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
        for r in s15["rows"]:
            A(f"| `{r['stratum']}` | {r['m_s']} | {r['alpha_s']:.4g} | "
              f"{r['n_count_path_rows']} | {r['n_measurable']} | "
              f"{r['n_unmeasurable']} | "
              + (f"{r['A_min_min']:.4f}..{r['A_min_max']:.4f}"
                 if r["A_min_min"] is not None else "n/a")
              + f" | {r['n_measurable_at_standing_bound']} | "
              f"{r['max_obs_over_floor']:.3f} |")
        A("")
        A(f"- declared S-15 target amplitude: "
          f"**{s15['target_amplitude']:.0%}** rate modulation; second reference "
          f"**{s15['second_reference_amplitude']:.1%}** (the programme's own "
          f"standing bound), N = {s15['N_events']:.0f} events.")
        for line in _wrap(s15["conventions"]):
            A("  " + line)
        for line in _wrap(s15["vif_rule"]):
            A("  " + line)
        A(f"- per-feature measured VIF available for "
          f"{s15['n_features_with_own_vif']} features; every other row uses the "
          f"pooled flat 24.08 fallback.")
        for line in _wrap(
            "ASYMMETRY, STATED RATHER THAN AVERAGED AWAY. F4-58 has a per-feature "
            "VIF for every REAL feature and none for any CONTROL, so the real arm "
            "runs on measured floors spanning the full 30-800 d block ladder while "
            "the control arm runs on the single pooled fallback. The control "
            "stratum's floor is therefore uniform where the real stratum's is not, "
            "and any cross-arm comparison of obs/floor inherits that. Direction, "
            "DERIVED rather than asserted: the real arm's measured VIFs run 15.8 to "
            "120.9 with a median above the pooled 24.08, so the real arm's floors "
            "are on the whole HIGHER and its obs/floor ratios on the whole LOWER "
            "than they would be under a common VIF. The asymmetry therefore "
            "flatters the CONTROL arm, which is the safe direction for a "
            "calibration whose whole purpose is to catch the real arm claiming too "
            "much."
        ):
            A("  " + line)
        A("")
    m = max(state["n_tests"], 1)
    bh_thresh = cfg["fdr_q"] / m
    floors = [t.get("p_floor", 1.0 / (cfg["n_surrogates"] + 1.0)) for t in tests]
    floor = min(floors) if floors else 1.0 / (cfg["n_surrogates"] + 1.0)
    n_above = sum(1 for fl in floors if fl > bh_thresh)
    enum_f = sorted({t["floor_enumeration"] for t in tests
                     if t.get("floor_enumeration") is not None})
    mc_f = sorted({t["floor_monte_carlo"] for t in tests
                   if t.get("floor_monte_carlo") is not None})
    A(f"- **MULTIPLICITY, PRICED AT THE DECLARED COUNT, AND THE FLOORS THAT LIMIT IT "
      f"(stated together, S-8 caveat).** Declared family = **{m} tests**; BH at "
      f"q = {cfg['fdr_q']} demands **p <= {bh_thresh:.3g}** for the smallest. "
      f"**{n_above} of {m} declared tests have a resolution floor ABOVE that "
      f"threshold** and therefore cannot attain it at any effect size. "
      + ("Every test in this session is resolvable against the BH line."
         if n_above == 0 else
         f"For those {n_above}, the BH line is DESCRIPTIVE ONLY. Per S-8 -- this "
         f"program's multiplicity standard since round 1, and not BH -- the "
         f"confirmatory statistic for a scanned family is the MAXIMUM absolute "
         f"effect over the whole declared family compared to that same maximum "
         f"computed over sim catalogues through the identical code path, whose "
         f"family-wise p is resolvable at 1/(N+1) regardless of family size. Scan "
         f"size costs POWER, not RESOLUTION."))
    A(f"- **the two floors are different things and are reported separately.** "
      f"The ENUMERATION floor comes from the exhaustive circular-shift null and is "
      f"a property of the RECORD LENGTH -- no surrogate budget can lower it, and "
      f"the adaptive ladder never touches it: "
      + (", ".join(f"{v:.3g}" for v in enum_f[:8]) or "n/a")
      + (" ..." if len(enum_f) > 8 else "")
      + f". The MONTE CARLO floor comes from the surrogate BUDGET (1/(N_max+1), "
        f"where N_max is the PRE-DECLARED cap, not the number of draws a laddered "
        f"test actually made): "
      + (", ".join(f"{v:.3g}" for v in mc_f[:8]) or "n/a")
      + (" ..." if len(mc_f) > 8 else "") + ".")
    k_min = int(math.ceil(floor * m / cfg["fdr_q"]))
    n_at_floor = sum(1 for t, fl in zip(tests, floors)
                     if t["p_raw"] <= fl + 1e-12)
    A(f"- **resolution (read this before the table):** the smallest attainable "
      f"empirical p anywhere in this session is {floor:.5f}. Every p at its own "
      f"floor is CENSORED -- the true p is somewhere at or below it -- so the BH q "
      f"attached to it is an upper bound computed from a tie, not a measurement. "
      + (f"**No test in this session can pass FDR at all** ({k_min} > "
         f"{m} tests would have to tie at the floor)."
         if k_min > m else
         f"At this resolution BH can only reject if at least {k_min} "
         f"{'test ties' if k_min == 1 else 'tests tie'} at the floor "
         f"simultaneously; {n_at_floor} sit at their own floor."
         + ("" if k_min == 1 else
            " A survivor list produced from a tie at the floor is provisional --"
            " rerun with more surrogates before believing the ordering.")))
    A(f"- the period scan carries its own, coarser floor: its max-power Monte Carlo "
      f"is capped at {state['results']['period_scan']['scan']['n_mc']} draws "
      f"(each draw is a full periodogram), so no period peak can report p below "
      f"{1.0 / (state['results']['period_scan']['scan']['n_mc'] + 1):.5f}.")
    A(f"- ledger: one line appended to `engine/EXPLORE_COUNT.jsonl` "
      f"(kind=mine, n_tests={state['n_tests']}), so holdout multiplicity reporting "
      f"includes this sweep. No holdout hash was spent.")
    A("")

    # ------------------------- §P6-4 Rule 4.7: the region banner, five lines ----
    rp = state.get("regions")
    if rp:
        rc = cfg.get("regions", {})
        ar = rp["R_arithmetic"]
        fl = rp["floors"]
        n_regsum = sum(1 for t in tests if t["test"] == "regsum_score_2Rdf")
        n_batt = sum(1 for t in tests
                     if t["test"] == "glm_poisson_offset_etas_region")
        A("## Region battery (§P6-4 Rules 4.1-4.7, §P7-1(b)/(d))")
        A("")
        A(f"- **Partition rule `{rp['rule_id']}`, R = {rp['R']}, digest "
          f"`{rp['digest']}`.** Selector: {rp['selector']}. The K-080 census cell "
          f"list is NOT the selector (§P6-4 Finding A); "
          f"{rp['n_cells_assigned']}/{rp['n_cells_total']} cells assigned, "
          f"{rp['n_events_assigned']:.0f} events in, {rp['n_events_excluded']:.0f} "
          f"excluded by the activity threshold.")
        A(f"- **Declared test count from this axis: {n_regsum} regional-sum tests"
          + (f" + {n_batt} per-region battery tests (R x the declared count)"
             if n_batt else " (battery OFF)")
          + f"**; at q = {M.FDR_Q} over {len(tests)} declared tests, "
            f"{M.FDR_Q * len(tests):.1f} survivors are expected BY CHANCE.")
        A(f"- **S-15, from the FORMULA (§P7-1(b)) not a scalar:** "
          f"`A_min = sqrt(VIF)*(z_alpha + z_0.80)*sqrt(2/N)` with **VIF = "
          f"{fl['vif']:.4g}** (MEASURED, F4-58 per §P7-1(c) -- not the 3.94 "
          f"§P7-1(a) inferred), **alpha = {fl['alpha']:.4g}** (declared tranche "
          f"threshold, z = {fl['z_alpha']:.3f}), target amplitude "
          f"{fl['target_amplitude']:.3g}. Per-region: "
          f"**{fl['n_unmeasurable']}/{len(fl['rows'])} regions UNMEASURABLE "
          f"({100.0*fl['fraction_unmeasurable']:.0f}%)**.")
        A(f"- **The 2R-df SUM at R = {rp['R']} is "
          f"{'MEASURABLE' if ar['sum_measurable'] else 'UNMEASURABLE'}**: at the "
          f"declared target amplitude the sum has 80% power up to R = "
          f"{ar['max_R_for_2Rdf_sum_at_80pct_power']} "
          f"(non-centrality {ar['noncentrality_at_target']:.1f} on "
          f"N = {ar['N_window']:.0f} events). The per-region battery's own floor "
          f"resolves the target only up to R = "
          f"{ar['max_R_for_per_region_battery']} "
          f"({ar['max_R_for_per_region_battery_real']:.2f} before flooring) -- "
          f"which is §P7-1(d)'s anticipated outcome, stated rather than hidden.")
        A("- **Per-region amplitudes are UNRESOLVED and are not quoted** (§P6-4 "
          "Rule 4.3 as amended by §P7-1(d)); they appear in `stubs.json` under "
          "`region_amplitudes_unresolved`, labelled, and never in the ranked list. "
          "The quotable object is the summed statistic alone.")
        A("- **G-M1 restated (§P6-4 Rule 4.7 item 5):** no output here may be "
          "entered for or against any ledger entry until G-M1 clears, and no bound "
          "may be reported at any aggregation level without demonstrated "
          "planted-signal recovery AT THAT AGGREGATION. Recovery at the regional "
          "aggregation is demonstrated at engine-test level in "
          "`engine/tests/test_regsum.py`; that is NOT the full G-M1 run.")
        A("")
        A("| region | lon | N events | A_min (formula) | S-15 |")
        A("|---|---|---:|---:|---|")
        for r in fl["rows"]:
            A(f"| {r['region']} | [{r['lon_lo']:+.0f}, {r['lon_hi']:+.0f}) | "
              f"{r['n_events']:.0f} | {r['a_min_pct']:.1f}% | {r['s15']} |")
        A("")

    # ---------------------------------- §P6-3: the stratum table + S-8 ---------
    bh_meta = state.get("bh")
    ms = state.get("max_statistic")
    if bh_meta:
        strat_on = bool(cfg.get("strata"))
        A("## Multiplicity partition and the max-statistic (§P6-3)")
        A("")
        A(f"- **mode: {'STRATIFIED (weighted BH)' if strat_on else 'UNSTRATIFIED (flat BH)'}**"
          + (f", partition `{cfg['strata']['path']}` sha256 "
             f"`{cfg['strata']['sha256'][:16]}` -- the FILE CONTENT is "
             f"hash-affecting, so re-partitioning is a new session with a new "
             f"config hash and a new ledger line (§P6-3(5))."
             if strat_on else
             ". Flat BH over one stratum holding the whole declared family, which "
             "is the default; `--strata <partition.json>` switches this on."))
        A("")
        for line in strata_mod.budget_table_lines(bh_meta["table"],
                                                  bh_meta["identity"]):
            A(line)
        A("")
        for line in _wrap(strata_mod.HONEST_TRADE):
            A(line)
        A("")
        A("**Per-stratum BH at a flat q would NOT control global FDR.** Applying BH "
          "at level q independently within each of S families controls the AVERAGE "
          "OVER FAMILIES of FDR (Benjamini & Bogomolov 2014), not the overall FDR "
          "across all tests. The identity above is what converts it into weighted "
          "BH with pre-specified weights, which controls global FDR at q exactly, "
          "and the engine refuses to run without it.")
        A("")
        A("Errored, ineligible and never-executed tests count as NON-REJECTIONS "
          "against their declared `m_s` (§P6-3(3)): the denominator is the "
          "DECLARATION, not the execution. The `ineligible/errored` column is that "
          "count.")
        A("")
    if ms:
        A("### S-8 sim-calibrated max-statistic -- the anti-repartition guarantee")
        A("")
        for line in _wrap(strata_mod.MAXSTAT_NOTE):
            A(line)
        A("")
        A(f"- **GLOBAL max-statistic p = {ms['global']['p']:.4g}** over "
          f"{ms['global']['n_tests']} covered tests x "
          f"{ms['global']['n_replicates']} shared surrogate worlds "
          f"(floor {ms['global']['floor']:.3g}). Statistic: "
          f"{ms['global'].get('statistic')}.")
        for line in _wrap(ms.get("coverage_note", "")):
            A("  " + line)
        A("")
        A("| stratum | covered tests | max-statistic p | **GLOBAL p (adjacent, "
          "§P6-3(4))** |")
        A("| --- | ---: | ---: | ---: |")
        for r in ms["per_stratum"]:
            pv = "n/a (no covered test)" if not np.isfinite(r["p"]) else f"{r['p']:.4g}"
            A(f"| `{r['stratum']}` | {r['n_tests']} | {pv} | "
              f"**{r['global_p']:.4g}** |")
        A("")
        A("No stratum is reported without the global figure on the same row, per "
          "§P6-3(4). The global number is invariant to how the family is "
          "partitioned -- that is exactly what makes stratification a discipline "
          "rather than a knob.")
        A("")

    # ------------------------------------- §P6-2: the p-method table -----------
    gconf = state.get("gpd_confirmations") or []
    if cfg.get("gpd") or any(t.get("p_method") == gpd_tail.P_GPD for t in tests):
        A("## GPD tail extrapolation (§P6-2)")
        A("")
        A(gpd_tail.census_line(_cen, _m_decl))
        A("")
        A("| test | p_method | p_raw (MC) | p entering BH | p_AD | xi stability | "
          "reason |")
        A("| --- | --- | ---: | ---: | ---: | --- | --- |")
        for t in sorted(tests, key=lambda t: (t.get("p_bh", t["p_raw"]),
                                              t.get("order_key", []))):
            if t.get("p_method") == gpd_tail.P_MC_RESOLVED and not t.get("gpd"):
                continue
            g = t.get("gpd") or {}
            pad = g.get("p_ad")
            A(f"| `{t['feature']}` {_row_where(t)} | `{t.get('p_method')}` | "
              f"{t['p_raw']:.4g} | {t.get('p_bh', t['p_raw']):.4g} | "
              f"{('%.3f' % pad) if pad is not None else 'n/a'} | "
              f"{g.get('xi_stability_pass')} | "
              f"{(t.get('p_method_reason') or '')[:160]} |")
        A("")
        if gconf:
            A("### §P6-2(6) confirmation of GPD candidates")
            A("")
            A("A survivor labelled `GPD_EXTRAPOLATED` **is not a survivor**. BH ran "
              "with it at its CI-upper; surviving makes it a "
              f"`{gpd_tail.CANDIDATE_LABEL}`, and it becomes a stub only after a "
              "targeted single-test Monte Carlo at N >= 10/p_gpd resolves its p "
              "directly. **Only the brute-force p is written to `stubs.json`.**")
            A("")
            A("| test | p_gpd (CI-upper) | N required | status | p_brute_force |")
            A("| --- | ---: | ---: | --- | ---: |")
            for r in gconf:
                pbf = r.get("p_brute_force")
                A(f"| `{r['feature']}` ({r['test']}) | {r['p_gpd_ci_upper']:.4g} | "
                  f"{r['n_required']} | {r['status']} | "
                  f"{('%.4g' % pbf) if pbf is not None else '--'} |")
            A("")
    dl = state.get("data_log", [])
    if dl:
        A("### Optional downloads (family 3)")
        A("")
        A("| source | status | detail |")
        A("| --- | --- | --- |")
        for r in dl:
            det = (f"{r.get('n_bytes','?')} B, sha256 `{r.get('sha256','')[:16]}`, "
                   f"coverage {r.get('coverage')}" if r.get("status") == "ok"
                   else str(r.get("reason", "")))
            A(f"| `{r['key']}` | {r['status']} | {det} |")
        A("")
        A("Frozen copies and hashes are in `engine/out/mine/data_log.jsonl` "
          "(sniff-grade hygiene). This is deliberately NOT `download_log.md`, which "
          "belongs to frozen-protocol runs only.")
        A("")

    aa = state.get("axis_audit")
    if aa:
        A("## Axis audit (K-089, the standing rule)")
        A("")
        for line in _wrap(M.AXIS_AUDIT_LINE):
            A(line)
        A("")
        A(f"- **tranche:** {aa['tranche']}")
        A(f"- **offset axis:** {aa['offset_axis']}")
        A(f"- **lag axis, SCANNED ({len(aa['lag_axis_scanned'])} features, "
          f"{aa['n_new_lag_tests_cyclic']} new cyclic tests):** "
          + (", ".join(f"`{n}`" for n in aa["lag_axis_scanned"]) or "none"))
        A(f"- **lag axis, PROVABLY FREE and therefore not scanned "
          f"({len(aa['lag_axis_provably_free_not_scanned'])}):** "
          + (", ".join(f"`{n}`" for n in aa["lag_axis_provably_free_not_scanned"])
             or "none"))
        A(f"- **lag axis, NEITHER scanned nor provably free "
          f"({len(aa['lag_axis_neither_scanned_nor_free'])}) -- the honest gap:** "
          + (", ".join(f"`{n}`" for n in aa["lag_axis_neither_scanned_nor_free"])
             or "none"))
        A(f"- **ladder-rung axis:** {aa['ladder_rung_axis']}")
        A("")
        A("### Lag-invariance proof, executed before the scan (min R^2 of the "
          "lagged [sin, cos] design on the lag-0 column space)")
        A("")
        A("| phase feature | declared | lag 1 | 3 | 7 | 15 | 30 | worst | verdict |")
        A("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |")
        for r in aa.get("lag_invariance", []):
            byl = r["min_r2_by_lag"]
            cells = " | ".join(f"{byl[str(L)] if str(L) in byl else byl[L]:.5f}"
                               for L in (1, 3, 7, 15, 30) if (L in byl or str(L) in byl))
            A(f"| `{r['feature']}` | {r['declared']} | {cells} | "
              f"{r['worst_min_r2']:.5f} | {r['verdict']} |")
        A("")
        A(f"Tolerance: FREE means min R^2 >= 1 - {M.LAG_FREE_TOL:g}. A feature "
          "declared FREE that measures below tolerance is a BUG; one declared PRICED "
          "that measures at 1.0 is a MISSED SAVING. Both are reported, neither is "
          "silently corrected.")
        A("")

    A("## Ranked candidates (top 25 by BH q, then raw p)")
    A("")
    A("| # | feature | lag/period | effect | raw p | p_method | BH q | stratum | "
      "surrogates | ladder | aliasing | amplitude honesty |")
    A("| ---: | --- | --- | --- | ---: | --- | ---: | --- | ---: | --- | --- | --- |")
    for i, t in enumerate(order[:25], 1):
        lad = (t["ladder"]["verdict"] if t["test"] == "lomb_scargle_peak" else "n/a")
        ali = t.get("aliasing", {}).get("verdict", "-- (did not pass FDR)")
        # §P6-2(5): p_method is on EVERY row of the table, not only the survivors.
        A(f"| {i} | `{t['feature']}` | {_row_where(t)} | {_row_effect(t)} | "
          f"{t['p_raw']:.4g} | `{t.get('p_method')}` | {t['bh_q']:.4g} | "
          f"{t.get('stratum', '')} | {t.get('n_surrogates', '')} | "
          f"{lad} | {ali} | {_amp_note(t)} |")
    A("")

    A("## Survivors (post-FDR) in full")
    A("")
    if n_pass == 0:
        A("**Nothing survived BH-FDR at q = %.2f.**" % cfg["fdr_q"])
        A("")
        A("That is a clean result for the instrument and is reported as such: with "
          "the ETAS baseline absorbing Omori-Utsu clustering, the global daily "
          "residual carries no ephemeris-, space-weather- or catalogue-mark "
          "association that survives this session's multiplicity. It is NOT evidence "
          "of absence -- see the power note below.")
    for t in order:
        if not t["passes_fdr"]:
            continue
        A(f"### `{t['feature']}` -- {_row_where(t)}")
        A("")
        A(f"- test: `{t['test']}`, {t.get('df', '?')} df")
        A(f"- effect: {_row_effect(t)}; {_amp_note(t)}")
        A(f"- raw p = {t['p_raw']:.4g} (empirical, {t.get('n_surrogates')} "
          f"circular-shift surrogates); BH q = {t['bh_q']:.4g}")
        if t["test"] == "glm_poisson_offset_etas":
            A(f"- in-sample bits/event vs etas-v1: {t['bits_per_event']:+.6f}")
        if t["test"] == "lomb_scargle_peak":
            lad = t["ladder"]
            A(f"- harmonic ladder ({lad['n_bins']} phase bins, epoch-folding "
              f"likelihood ratio): **{lad['verdict']}**")
            A("")
            A("  | rung | period (d) | fold LR |")
            A("  | --- | ---: | ---: |")
            for r in lad["rungs"]:
                A(f"  | {r['multiplier']:.4g} x P | {r['period_days']:.4g} | "
                  f"{r['fold_lr']:.2f} |")
            A("")
        ali = t.get("aliasing")
        if ali:
            A(f"- aliasing audit: **{ali.get('verdict')}** -- {ali.get('rule','')}")
            if "z_ratio_vs_1d" in ali:
                A(f"  - surrogate-calibrated effect ratio vs 1 d binning: "
                  f"{ali['z_ratio_vs_1d']}")
            if "fold_lr_ratio_vs_1d" in ali:
                A(f"  - folded-concentration ratio vs 1 d binning: "
                  f"{ali['fold_lr_ratio_vs_1d']}")
                u = ali["unbinned_rayleigh"]
                A(f"  - unbinned event-time Rayleigh (n = {u['n_events']} events, "
                  f"ETAS-simulated null): R = {u['R']:.3f}, p = {u['p']:.4g}")
        A("")
        A("- " + M.AMPLITUDE_HONESTY)
        A("")

    A("## Instrument validation -- what the miner rediscovered")
    A("")
    A("A miner run against a clean residual should find NOTHING at strong "
      "significance. Anything it does find that is already known physics or a known "
      "catalogue artifact is **validation of the instrument, not discovery**. The "
      "candidates below are pre-labelled as such wherever they match a known "
      "structure:")
    A("")
    A("| candidate | raw p | known-structure label |")
    A("| --- | ---: | --- |")
    for t in order[:25]:
        A(f"| `{t['feature']}` {_row_where(t)} | {t['p_raw']:.4g} | "
          f"{_known_label(t)} |")
    A("")

    A("## Method notes, including one honest deviation")
    A("")
    A("1. **Surrogate nulls.** Every GLM and mark test is calibrated against "
      "CIRCULAR TIME SHIFTS of the target (the counts/offset pair for GLM tests, "
      "the mark series along the time-ordered event sequence for mark tests). "
      "Because the score statistic and the rank correlations are exact "
      "cross-correlations, ALL admissible shifts are evaluated in closed form "
      "rather than sampled; when the requested surrogate count exceeds the number "
      "of admissible shifts the run uses every one of them and reports the actual "
      "count in the `surrogates` column. Shifts within 30 days of zero are excluded.")
    A(f"2. **DEVIATION, flagged not hidden.** A circular shift of an evenly sampled "
      f"series does not change its Lomb-Scargle power at all -- a shift is a pure "
      f"phase rotation -- so circular-shift surrogates are mathematically vacuous "
      f"as a null for the PERIOD SCAN. The period scan therefore uses two Monte "
      f"Carlo nulls instead: an AR(1) red-noise null matched to the residual's own "
      f"lag-1 autocorrelation (phi = "
      f"{state['results']['period_scan']['scan']['ar1_phi']:.3f}) and a permutation "
      f"(white) null, and reports the MORE CONSERVATIVE of the two p-values. "
      f"Circular shifts remain the null everywhere else.")
    A("2b. **Block length is measured, not guessed.** A block bootstrap is a valid "
      "null only if its blocks are long compared to the structure being tested; too "
      "short and the surrogates cannot represent the feature's own timescale, which "
      "makes the null too narrow and the test ANTI-CONSERVATIVE. Cyclic features use "
      "2x their known period; everything else uses 4x the e-folding time of its own "
      "autocorrelation, clipped to [30, 800] days. The per-test value is recorded as "
      "`block_days`. Residual risk, stated: a feature whose autocorrelation is fast "
      "but whose ENVELOPE is slow (Ap is the example -- storm-scale ACF, solar-cycle "
      "envelope) gets a short block and its p is therefore optimistic.")
    A("3. **Harmonic ladder.** Every candidate period P is scored by epoch-folding "
      "likelihood ratio at {P/3, P/2, P, 2P, 3P}, extended at the winning edge "
      "while the score keeps improving. The reported period is the winning rung and "
      "all rung scores are printed. Rungs longer than a third of the record "
      f"({counts.size / 3.0:.0f} d here) are not scored at all: with fewer than "
      "three observed cycles an epoch fold measures the record length, not a "
      "period.")
    A("3b. **The two period-scan nulls disagree on purpose, and the report takes "
      "the loser.** The global daily residual is red (AR(1) phi = "
      f"{state['results']['period_scan']['scan']['ar1_phi']:.3f}), so a permutation "
      "(white) null calls essentially every peak significant. The AR(1) null is the "
      "one that knows the residual is autocorrelated. Reporting the more "
      "conservative of the two is what keeps the period scan from manufacturing a "
      "candidate list out of red noise -- and the gap between the two columns is "
      "the size of the mistake that would have been made.")
    A("4. **Aliasing audit.** Every post-FDR survivor is re-tested at 2-day and "
      "7-day binning, and period claims are additionally re-tested on UNBINNED "
      "event times (real catalogue timestamps, sub-day precision) against an "
      "ETAS-simulated inhomogeneous-Poisson null. A pattern whose effect halves "
      "under re-binning, or fails unbinned, is flagged LATTICE-SUSPECT: a pattern "
      "that moves when the lattice moves is the lattice.")
    A("5. **Causality.** Family-1 and family-2 features are deterministic functions "
      "of t (ephemeris) and are therefore EXEMPT from the causality-shuffle test by "
      "construction: no rearrangement of the catalogue can change them. Family-3 "
      "(downloaded indices) are lagged one day so day t uses only values published "
      "strictly before t. Family-4 (catalogue-derived) are strictly trailing "
      "windows, exclusive of today, and ARE included in "
      "`engine/tests/test_causality.py`.")
    A("6. **In-sample.** Every effect size here is fitted and evaluated on the "
      "exploration window. They are upper bounds, exactly as `--mode explore` is.")
    lad_cfg = cfg.get("ladder")
    if lad_cfg:
        A(f"7. **ADAPTIVE SURROGATE LADDER WAS ON for this session** "
          f"(`--ladder`, rule `{lad_cfg.get('rule')}`, h = {lad_cfg.get('h')}, "
          f"N_max = {lad_cfg.get('n_max')}). Each test's Monte Carlo null was drawn "
          f"until h surrogates reached its observed statistic; its p is then h/l, "
          f"where l is the draw index of the h-th exceedance, or (1+c)/(1+N_max) if "
          f"the cap was reached with c < h. That is the Besag & Clifford (1991) "
          f"sequential Monte Carlo p-value, EXACTLY valid under optional stopping -- "
          f"an ordinary MC p computed after peeking would not be. Uninteresting "
          f"tests stop after {lad_cfg.get('h')} exceedances instead of consuming the "
          f"full {lad_cfg.get('n_max')}; interesting ones still pay the full cap, so "
          f"the Monte Carlo floor at the sharp end is unchanged.")
        A(f"   - **N_max is fixed per declared stratum BEFORE the run** and depends "
          f"only on the test key (feature, kind, lag). There is no rank-based "
          f"escalation: no observed statistic, p-value or survivor list can raise a "
          f"test's budget, because that would invalidate the sequential p-value.")
        A(f"   - **The ladder applies to MONTE CARLO nulls only** (block bootstrap; "
          f"AR(1) and permutation where used). It does NOT touch the exhaustive "
          f"circular-shift enumeration, whose floor is a property of the record "
          f"length. The two floors are reported separately in the configuration "
          f"section above.")
        A(f"   - Chunk size {lad_cfg.get('chunk')} is a VECTORISATION DETAIL, not a "
          f"statistical stage: the exact stopping index is recovered inside "
          f"whichever chunk trips the rule, so the answer does not depend on it. "
          f"Per-test stopping detail is in the `ladder_stop` field of "
          f"`checkpoint.json`. All stopping-rule parameters are part of the config "
          f"hash: a run under a different rule is a different experiment.")
    else:
        A("7. **Adaptive surrogate ladder: OFF** (the default). Every test drew the "
          "full surrogate budget; no optional stopping was applied, so the p-values "
          "here are ordinary fixed-sample Monte Carlo p-values.")
    A("")
    A("## Power / bounds note")
    A("")
    sd = float((counts - offset).std())
    mean = float(counts.mean())
    A(f"The global daily residual has mean count {mean:.2f} events/day and residual "
      f"sd {sd:.2f} over {counts.size} days. A sinusoidal rate modulation of "
      f"amplitude A is detectable at this session's thresholds only above roughly "
      f"{100 * 2.5 * sd / (mean * math.sqrt(counts.size / 2.0)):.2f}% "
      f"(a 2.5-sigma-equivalent rule of thumb, not a formal power calculation). "
      f"Absence of a survivor is a BOUND at about that amplitude, not an absence of "
      f"effect.")
    A("")
    A("That rule of thumb assumes independent days, so it is OPTIMISTIC -- and most "
      "optimistic exactly where the block bootstrap is longest. A feature tested "
      "with 700-day blocks has an effective sample of ~11 independent pieces, not "
      "7716, and its true detection threshold is several times the number above. "
      "Read the per-test `block_days` before quoting any bound: the longer the "
      "block, the weaker the bound.")
    A("")
    A("---")
    A("")
    A(M.AMPLITUDE_HONESTY)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return path


KNOWN = [
    (("annual_phase", "sun_declination"), "KNOWN: annual/seasonal cycle -- seasonal "
     "catalogue completeness, hydrological loading and reporting cadence all live "
     "here. Validation of the instrument, not discovery."),
    (("b_value_90d", "deep_fraction_30d", "mean_depth_30d"),
     "KNOWN: catalogue-composition drift (network eras, Mc drift, deep-slab "
     "sequences). Expected to be found; validation, not discovery."),
    (("moon_synodic_phase", "moon_anomalistic_phase", "moon_draconic_phase",
      "spring_neap_phase", "half_draconic_phase", "perigean_spring_beat",
      "eclipse_year_beat", "annual_synodic_beat", "perigee_syzygy",
      "tidal_potential_proxy", "moon_distance", "sun_moon_elongation",
      "moon_declination", "moon_abs_declination", "declination_product"),
     "KNOWN-CLASS: solid-earth tidal forcing. The tidal corpse says expect a bound "
     "at the sub-percent level, not a detection."),
    (("F107_solar_flux", "Ap_geomagnetic", "length_of_day"),
     "OPEN-CLASS: space weather / rotation. Contested literature; no accepted "
     "mechanism at this amplitude."),
]


def _known_label(t):
    name = t["feature"]
    for names, label in KNOWN:
        if name in names:
            return label
    if t["test"] == "lomb_scargle_peak":
        p = t["period_days"]
        for ref, lab in ((365.25, "annual"), (182.6, "semi-annual"),
                         (29.53, "synodic month"), (14.77, "spring-neap"),
                         (27.55, "anomalistic month"), (13.66, "declination tide"),
                         (7.0, "WEEKLY -- catalogue/reporting cadence, an artifact"),
                         (1.0, "diurnal")):
            if abs(math.log(p / ref)) < 0.03:
                return f"KNOWN period: {lab}. Validation of the instrument."
        return "UNLABELLED period -- no known structure matched; treat as candidate."
    return "unlabelled"


def _count_by(items, keyfn):
    out = {}
    for i in items:
        out[keyfn(i)] = out.get(keyfn(i), 0) + 1
    return out


def _wrap(text, width=88):
    words, line, out = text.split(), "", []
    for w in words:
        if len(line) + len(w) + 1 > width:
            out.append(line)
            line = w
        else:
            line = (line + " " + w).strip()
    if line:
        out.append(line)
    return out


# ------------------------------------------------------------------- stubs ---
def write_stubs(session_dir, cfg, tests):
    path = os.path.join(session_dir, "stubs.json")
    # §P7-17, enforced at the one place it matters most. `stubs.json` is the file a
    # human reads as "things worth looking at", so a COMPONENT-OF row arriving here
    # is precisely the failure the ruling names: a half-claim being read as a claim,
    # against a threshold priced for a different test.
    disp.assert_no_component_standalone(tests, "stubs.json")
    entries, candidates = [], []
    for t in sorted(tests, key=lambda t: (t["bh_q"], t["p_raw"], t.get("order_key", []))):
        if not t["passes_fdr"]:
            continue
        # §P6-2(6). A SURVIVOR LABELLED GPD_EXTRAPOLATED IS NOT A SURVIVOR: it is a
        # candidate, and it becomes a stub only after the targeted brute-force MC
        # has resolved its p directly. `confirm_gpd_candidates` relabels a confirmed
        # row to MC_RESOLVED and attaches `p_brute_force`; anything still carrying
        # the GPD label at this point was not confirmed and is listed separately,
        # never as a stub.
        if t.get("p_method") == gpd_tail.P_GPD:
            candidates.append({
                "status": gpd_tail.CANDIDATE_LABEL,
                "observable": _observable(t),
                "feature": t["feature"], "test": t["test"],
                "p_method": gpd_tail.P_GPD,
                "p_gpd_ci_upper": t.get("p_bh"),
                "bh_q": t["bh_q"],
                "confirmation": t.get("gpd_confirmation"),
                "why_not_a_stub": (
                    "§P6-2(6): BH ran with this row at its 95%-CI upper end, and it "
                    "survived -- which makes it a CANDIDATE, not a stub. Only a "
                    "targeted single-test Monte Carlo at N >= 10/p_gpd may resolve "
                    "it, and only that brute-force p may be written here."),
            })
            continue
        entries.append({
            "status": "DRAFT STUB -- generator output, not a K-entry",
            "observable": _observable(t),
            "feature": t["feature"],
            "family": t["family"],
            "test": t["test"],
            "where": _row_where(t),
            "effect_size": _row_effect(t),
            "rate_modulation": _amp_note(t),
            # §P6-6 R3: the label is on the STUB, not only in the report prose. A
            # stub is the object that travels; an unlabelled amplitude inside one
            # is exactly how a selection maximum gets read as an effect size.
            "amplitude_label": t.get("amplitude_label"),
            "a_min_formula": t.get("a_min_formula"),
            "a_min_vif": t.get("a_min_vif"),
            "a_min_vif_source": t.get("a_min_vif_source"),
            "s15": t.get("s15"),
            "obs_over_floor": t.get("obs_over_floor"),
            "is_negative_control": is_control_row(t),
            "p_raw": t["p_raw"],
            # §P6-2(5): the label is on EVERY row, so a reader never has to infer
            # which kind of number they are looking at.
            "p_method": t.get("p_method"),
            "p_method_reason": t.get("p_method_reason"),
            "p_entered_into_bh": t.get("p_bh", t["p_raw"]),
            "p_brute_force": t.get("p_brute_force"),
            "stratum": t.get("stratum"),
            "bh_q": t["bh_q"],
            "n_surrogates": t.get("n_surrogates"),
            "ladder": t.get("ladder", {}).get("verdict"),
            "aliasing_verdict": t.get("aliasing", {}).get("verdict"),
            "known_structure_label": _known_label(t),
            "caveats": [M.GENERATOR_NOT_EVIDENCE, M.AMPLITUDE_HONESTY],
            "next_step": ("Write up as a K-entry with a Popper ruling, then -- and "
                          "only then -- freeze a config and spend ONE holdout hash "
                          "via `python -u -m engine.cli run --mode holdout`."),
        })
    # §P6-4 Rule 4.3 as AMENDED by §P7-1(d). Per-region amplitudes are COMPUTED and
    # carried here, labelled UNRESOLVED with the rule that makes them unresolved.
    # They are deliberately NOT in `stubs` and NOT in the report's ranked list: this
    # block is the whole of their existence in the output.
    unresolved = []
    for t in sorted(tests, key=lambda t: t.get("order_key", [])):
        for r in t.get("regions_unresolved", []) or []:
            if r.get("quotable"):
                continue
            unresolved.append(dict(
                r, feature=t["feature"], lag=t["lag"], test=t["test"],
                region_rule_id=t.get("region_rule_id"),
                region_digest=t.get("region_digest"),
                why_not_a_stub=regions_mod.UNRESOLVED_RULE,
            ))
    for t in tests:
        if (t.get("test") == "glm_poisson_offset_etas_region"
                and t.get("p_method") == regions_mod.UNRESOLVED):
            unresolved.append({
                "feature": t["feature"], "lag": t["lag"], "test": t["test"],
                "region": t["region"], "n_events_in_window": t["n_events_in_window"],
                "amplitude_log_rate": t["amplitude_log_rate"],
                "pct_rate_modulation": t["pct_rate_modulation"],
                "a_min_formula": t["a_min_formula"], "s15": t["s15"],
                "p_method": regions_mod.UNRESOLVED,
                "p_method_reason": t["p_method_reason"],
                "why_not_a_stub": regions_mod.UNRESOLVED_RULE,
            })

    payload = {
        "banner": M.GENERATOR_NOT_EVIDENCE,
        "standing_warning_eq24": M.MIGNAN_BROCCARDO,
        "session": os.path.basename(session_dir),
        "config": cfg,
        "n_tests": len(tests),
        "n_stubs": len(entries),
        "n_region_amplitudes_unresolved": len(unresolved),
        "region_amplitudes_unresolved_rule": regions_mod.UNRESOLVED_RULE,
        "region_amplitudes_unresolved": unresolved,
        "p_method_census": gpd_tail.census(tests),
        "p_method_census_line": gpd_tail.census_line(gpd_tail.census(tests),
                                                     len(tests)),
        "n_gpd_candidates": len(candidates),
        "gpd_candidates_requiring_brute_force": candidates,
        "stubs": entries,
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    return path


def _observable(t):
    if t["test"] == "regsum_score_2Rdf":
        return (f"Daily M>=4.5 occurrence rate depends on `{t['feature']}` at lag "
                f"{t['lag']} d in AT LEAST ONE of {t['R']} declared regions, with "
                f"ARBITRARY per-region phase, after ETAS residualization "
                f"(phase-incoherent 2R-df sum, §P6-4 Rule 4.2).")
    if t["test"] == "glm_poisson_offset_etas_region":
        return (f"Daily M>=4.5 occurrence rate in region {t['region']} "
                f"[{t['lon_lo']:+.0f}, {t['lon_hi']:+.0f}) deg longitude depends on "
                f"`{t['feature']}` at lag {t['lag']} d, after ETAS residualization.")
    if t["test"] == "lomb_scargle_peak":
        return (f"Global daily M>={4.5} occurrence residual against etas-v1 carries a "
                f"periodicity at {t['ladder']['winning_period_days']:.4g} d.")
    if t["test"] == "glm_poisson_offset_etas":
        return (f"Global daily M>={4.5} occurrence rate depends on `{t['feature']}` "
                f"at lag {t['lag']} d, after ETAS residualization.")
    return (f"Event {t['mark']} is associated with `{t['feature']}` at the time of "
            f"the event, after ETAS residualization of occurrence.")
