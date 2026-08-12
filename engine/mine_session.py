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
import json
import math
import multiprocessing as _mp
import os
import time
import traceback

import numpy as np

from . import (baseline as bl, design, mine as M, splits, __version__)

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
def glm_task(f, window, counts, offset, n_surr, rng_for_lag, ladder=None, seed=0):
    rows = []
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
                h=int(ladder.get("h", 25)), chunk=int(ladder.get("chunk", 200)))
            p_boot = lad["p"]
        else:
            Sb = M.score_stat_block_bootstrap(X, counts, offset, n_surr, rng,
                                              mean_block=f.block_days)
            p_boot = M.bootstrap_p(S[0], Sb)
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
        rows.append(row)
    return rows


def marks_task(f, window, fe_day, marks, n_surr, rng, ladder=None, seed=0):
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
                        ladder=lad_cfg, rung_rng=rr)
        r.update({"feature": f.name, "family": f.family, "kind": f.kind,
                  "mark": mk, "n_events": int(marks[mk].size)})
        rows.append(r)
    return rows


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
    peaks, meta = M.period_assemble(obs["power"], periods, obs["peaks"], maxima,
                                    obs["ar1_phi"], n_mc)
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
    if kind == "glm":
        f = W["feats"][name]
        return glm_task(f, W["window"], W["counts"], W["offset"], n_surr,
                        lambda lag: task_rng(seed, name, "glm", lag,
                                             null_type="circular_shift"),
                        ladder, seed=seed)
    if kind == "marks":
        f = W["feats"][name]
        return marks_task(f, W["window"], W["fe_day"], W["marks"], n_surr,
                          task_rng(seed, name, "marks", None,
                                   null_type="circular_shift"),
                          ladder, seed=seed)
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
        "mag_target": float(args.mag_target),
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
    task_keys = ([f"glm:{f.name}" for f in feats]
                 + [f"marks:{f.name}" for f in feats]
                 + period_keys)

    # Every random stream in the session is addressed by a canonical test key.
    # Two tests sharing a digest would share a stream, which is a silent
    # correlation between "independent" nulls -- so it is checked, not assumed.
    all_keys = []
    for f in feats:
        for lag in f.lags:
            all_keys.append(M.test_key(cfg["seed"], f.name, "glm", lag=lag,
                                       null_type="block_bootstrap"))
        for mk in ("mag", "depth"):
            all_keys.append(M.test_key(cfg["seed"], f.name, f"mark_{mk}",
                                       null_type="block_bootstrap"))
    # The period scan's two nulls are now separate task families with separate
    # streams (phase 1b). No field was added to TEST_KEY_FIELDS, so every OTHER
    # declared digest in this session is byte-identical to phase 1.
    for _nt in M.PERIOD_NULLS:
        all_keys.append(M.test_key(cfg["seed"], "period_scan", "period",
                                   null_type=_nt))
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

    # The worker payload. Built for BOTH drivers: the sequential path runs the
    # period subtasks through the very same `dispatch_task`, which is what makes
    # --jobs 1 and --jobs N bit-identical there.
    payload = {
        "feats": {f.name: f for f in feats}, "window": window,
        "counts": counts, "offset": offset, "days": days, "resid": resid,
        "marks": marks, "fe_day": fe_day, "seed": int(cfg["seed"]),
        "n_surr": n_surr, "ladder": ladder, "cfg": cfg,
        "period_chunk": M.PERIOD_SURROGATE_CHUNK,
    }

    if n_jobs == 1:
        # -------------------------------------------- (a) GLM sweep --------
        for f in feats:
            key = f"glm:{f.name}"
            if ckpt.done(key):
                if verbose:
                    print(f"  [skip, checkpointed] {key}")
                continue
            t_f = time.perf_counter()
            rows = glm_task(f, window, counts, offset, n_surr,
                            lambda lag: rng, ladder, seed=int(cfg["seed"]))
            ckpt.record_seconds(key, time.perf_counter() - t_f)
            ckpt.put(key, rows)
            if verbose:
                best = min(rows, key=lambda r: r["p_raw"])
                print(f"  glm {f.name:<24s} {len(rows):>2d} lags in "
                      f"{time.perf_counter()-t_f:5.1f}s   best p={best['p_raw']:.4g} "
                      f"(lag {best['lag']}, {best['pct_rate_modulation']:+.2f}%/sd)")

        # ------------------------------------------ (c) mark tests ---------
        for f in feats:
            key = f"marks:{f.name}"
            if ckpt.done(key):
                if verbose:
                    print(f"  [skip, checkpointed] {key}")
                continue
            t_f = time.perf_counter()
            rows = marks_task(f, window, fe_day, marks, n_surr, rng, ladder,
                              seed=int(cfg["seed"]))
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
        for t in ckpt.state["results"][f"glm:{f.name}"]:
            t["order_key"] = [0, i, int(t["lag"]), 0, f.name]
            tests.append(t)
        for j, t in enumerate(ckpt.state["results"][f"marks:{f.name}"]):
            t["order_key"] = [0, i, -1, 1 + j, f.name]
            tests.append(t)
    for j, t in enumerate(ckpt.state["results"]["period_scan"]["peaks"]):
        t["order_key"] = [1, j, -1, 0, str(t["feature"])]
        tests.append(t)
    tests.sort(key=lambda t: t["order_key"])

    p = np.array([t["p_raw"] for t in tests])
    q, passed = M.benjamini_hochberg(p, M.FDR_Q)
    for t, qq, pa in zip(tests, q, passed):
        t["bh_q"] = float(qq)
        t["passes_fdr"] = bool(pa)
    n_tests = len(tests)

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
    ckpt.save()

    rep = write_report(session_dir, cfg, ckpt.state, tests, feats, counts, offset,
                       marks, base, window)
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


# ------------------------------------------------------------------ report ---
def _row_effect(t):
    if t["test"] == "glm_poisson_offset_etas":
        if t["kind"] == "phase":
            return f"amp {t['amplitude_log_rate']:.4f} log-rate"
        return f"beta {t['beta'][0]:+.4f} +/- {t['se'][0]:.4f} /sd"
    if t["test"] == "lomb_scargle_peak":
        return f"LS power {t['power']:.4f}"
    return f"{t['test']} {t['effect']:+.4f}"


def _row_where(t):
    if t["test"] == "glm_poisson_offset_etas":
        return f"lag {t['lag']} d"
    if t["test"] == "lomb_scargle_peak":
        return f"P = {t['period_days']:.4g} d"
    return f"mark {t['mark']}"


def _amp_note(t):
    if t["test"] == "lomb_scargle_peak":
        return f"+/-{t['pct_rate_modulation']:.1f}% folded rate"
    if t["test"] == "glm_poisson_offset_etas":
        return f"+/-{t['pct_rate_modulation']:.2f}% rate"
    return "rank correlation (no rate amplitude)"


def write_report(session_dir, cfg, state, tests, feats, counts, offset, marks,
                 base, window):
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
    A(f"- features: {len(feats)} "
      + ", ".join(f"family {k}: {v}" for k, v in sorted(
          _count_by(feats, lambda f: f.family).items())))
    A(f"- **{state['n_tests']} tests** in this session; BH-FDR at q = {cfg['fdr_q']} "
      f"-> **{n_pass} survive**")
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
    A("| # | feature | lag/period | effect | raw p | BH q | surrogates | ladder | "
      "aliasing | amplitude honesty |")
    A("| ---: | --- | --- | --- | ---: | ---: | ---: | --- | --- | --- |")
    for i, t in enumerate(order[:25], 1):
        lad = (t["ladder"]["verdict"] if t["test"] == "lomb_scargle_peak" else "n/a")
        ali = t.get("aliasing", {}).get("verdict", "-- (did not pass FDR)")
        A(f"| {i} | `{t['feature']}` | {_row_where(t)} | {_row_effect(t)} | "
          f"{t['p_raw']:.4g} | {t['bh_q']:.4g} | {t.get('n_surrogates', '')} | "
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
    entries = []
    for t in sorted(tests, key=lambda t: (t["bh_q"], t["p_raw"], t.get("order_key", []))):
        if not t["passes_fdr"]:
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
            "p_raw": t["p_raw"],
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
    payload = {
        "banner": M.GENERATOR_NOT_EVIDENCE,
        "standing_warning_eq24": M.MIGNAN_BROCCARDO,
        "session": os.path.basename(session_dir),
        "config": cfg,
        "n_tests": len(tests),
        "n_stubs": len(entries),
        "stubs": entries,
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    return path


def _observable(t):
    if t["test"] == "lomb_scargle_peak":
        return (f"Global daily M>={4.5} occurrence residual against etas-v1 carries a "
                f"periodicity at {t['ladder']['winning_period_days']:.4g} d.")
    if t["test"] == "glm_poisson_offset_etas":
        return (f"Global daily M>={4.5} occurrence rate depends on `{t['feature']}` "
                f"at lag {t['lag']} d, after ETAS residualization.")
    return (f"Event {t['mark']} is associated with `{t['feature']}` at the time of "
            f"the event, after ETAS residualization of occurrence.")
