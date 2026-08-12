"""Wall-clock benchmark for the EQ-24 miner v2 throughput layer.

Runs the SAME quick-preset-sized sweep four ways on real cached data and prints a
speedup table:

    --jobs 1 --no-ladder      the v1 audit baseline
    --jobs 0 --no-ladder      auto parallel = max(1, cpu_count() - 2)
    --jobs 1 --ladder         sequential, Besag-Clifford sequential stopping
    --jobs 0 --ladder         both

Every number printed is measured in this process by `time.perf_counter` around a
real `mine_session.run` call. Nothing here is extrapolated, modelled or reused from
a previous run.

WHAT IS TIMED. The design build and ETAS baseline fit are done ONCE, before any
cell, and are EXCLUDED from every row -- they are identical work in every cell, the
parallel layer does not touch them, and at the quick preset they are several times
larger than the sweep itself, so including them would flatten every speedup toward
1.0x and say nothing about the thing being measured. Their cost is printed
separately; end-to-end wall clock for a cell is its row plus that number.

Run:
    python -u -m engine.bench_v2 --no-download                    # all four cells
    python -u -m engine.bench_v2 --no-download --cells 1seq,1par  # a subset
    python -u -m engine.bench_v2 --no-download --repeat 3         # medians

Run from replication/. Requires the design cache to already exist -- build it with
`python -u -m engine.cli cache` first. Each cell writes a throwaway session under
engine/out/mine/ (gitignored scratch) and is forced to start fresh, because timing
a resumed session would time the checkpoint reader.

HISTORY, AND WHY THE TABLE STILL PRINTS THE CEILING. In phase 1 the sweep was
dominated by ONE task: `period_scan` cost ~101 s of a ~112 s sweep (90.7%), because
it computed 2 x n_mc full Lomb-Scargle periodograms inside a single task. The pool
boundary was the per-feature/per-test loop and `period_scan` was ONE member of it,
so the Amdahl ceiling was 1.10x however many workers were used.

Phase 1b decomposed it: the surrogate periodograms are now their own subtasks
("period_obs" plus "psurr:<null>:<start>"), seeded PER SURROGATE INDEX off the
canonical-key machinery, so chunking is provably not a statistical stage. The
ceiling line below is kept -- printed from the measured per-task seconds on every
run -- because a headline speedup should never be readable without the granularity
that bounds it. It bought 3.41x, which was 88-93% of a HARDWARE ceiling: the kernel
itself was the wall.

Phase 1c took the kernel. The cos/sin basis over (time grid, trial-period grid) is
identical for the observed series and for every surrogate in a scan, so it is built
once per worker process and each series then costs one BLAS pass
(engine.mine.lomb_scargle_power). Measured on this machine: 261 ms -> 2.9 ms per
7716 x 800 periodogram, a whole quick-preset period scan 102.5 s -> 1.9 s, with the
quick preset's period p-values EXACTLY unchanged.

Which moved the bottleneck, and the table below is now the only honest guide to
where. Read `--scaling` before choosing a job count: at the quick preset the sweep
is small enough that per-worker warm-up (including one ~0.66 s basis build per
process that touches a period task) can cost more than the parallelism saves.

WHY THE CELLS ARE NOT DIRECTLY COMPARABLE AS SCIENCE. The ladder cells do LESS
WORK, on purpose -- that is the speedup being measured -- and their p-values come
from a different (also valid) rule. The parallel cells use per-task RNG streams and
so do not reproduce the sequential numbers. This script measures TIME. It makes no
claim that the four cells produce the same statistics, and prints the per-cell test
counts and survivor counts so a reader can see what each cell actually did.
"""

from __future__ import annotations

import argparse
import json
import os
import statistics
import sys
import time

from . import datasets, mine as M, mine_session as ms, splits, __version__


class _Args:
    """The subset of the CLI namespace `build_config` reads."""

    def __init__(self, **kw):
        self.mag_target = 4.5
        self.dlat = 0.5
        self.dlon = 0.5
        self.explore_frac = 0.7
        self.data_dir = datasets.DEFAULT_DATA_DIR
        self.no_download = True
        self.seed = 20260811
        self.tranche1 = False
        self.ladder = False
        self.ladder_h = None
        self.ladder_chunk = None
        for k, v in kw.items():
            setattr(self, k, v)


# Scratch ledger. See the note in `_one`.
BENCH_LEDGER = os.path.join(M.MINE_DIR, "bench_explore_count.jsonl")

CELLS = {
    "1seq": {"jobs": 1, "ladder": False, "label": "--jobs 1  (v1 baseline)"},
    "1par": {"jobs": 0, "ladder": False, "label": "--jobs 0  (auto parallel)"},
    "1seqL": {"jobs": 1, "ladder": True, "label": "--jobs 1 --ladder"},
    "1parL": {"jobs": 0, "ladder": True, "label": "--jobs 0 --ladder"},
}


def _cfg_for(cell, preset, args_kw):
    return ms.build_config(_Args(ladder=CELLS[cell]["ladder"], **args_kw), preset)


def _one(cell, preset, args_kw, prepared, verbose=False):
    spec = CELLS[cell]
    cfg = _cfg_for(cell, preset, args_kw)
    t = time.perf_counter()
    # resume=False: a resumed session would measure the checkpoint reader, not the
    # sweep. Every cell therefore starts from zero completed tasks.
    #
    # `prepared` is passed in so the design build and ETAS fit -- byte-identical
    # work in every cell, and far larger than the sweep at this preset -- happen
    # ONCE and are excluded from the timed region. Including them would have every
    # cell measure the same large constant and flatten every speedup toward 1.0x.
    #
    # The benchmark must NEVER touch engine/EXPLORE_COUNT.jsonl: timing the same
    # sweep four times is not four hypotheses, and four ledger lines would inflate
    # the declared multiplicity of real science with runs nobody proposed.
    out = ms.run(cfg, verbose=verbose, resume=False, jobs=spec["jobs"],
                 ledger_path=BENCH_LEDGER, prepared=prepared)
    elapsed = time.perf_counter() - t
    return {
        "cell": cell, "label": spec["label"], "jobs_flag": spec["jobs"],
        "jobs_actual": ms.resolve_jobs(spec["jobs"]), "ladder": spec["ladder"],
        "seconds": elapsed, "n_tests": out["n_tests"], "n_pass": out["n_pass"],
        "config_hash": splits.config_hash(cfg)[:12],
        "session_dir": out["session_dir"],
        "task_seconds": _task_seconds(out["session_dir"]),
    }


def _task_seconds(session_dir):
    with open(os.path.join(session_dir, "checkpoint.json"), encoding="utf-8") as fh:
        return json.load(fh).get("task_seconds", {})


def _stages(task_seconds):
    """Per-stage compute, and the Amdahl ceiling the task granularity allows.

    This is printed because a speedup number without it is actively misleading. If
    one task is most of the sweep, no number of workers helps, and the honest
    reading of a 1.1x is "the pool worked and the granularity is wrong" rather than
    "parallelism does not help here".
    """
    glm = sum(v for k, v in task_seconds.items() if k.startswith("glm:"))
    marks = sum(v for k, v in task_seconds.items() if k.startswith("marks:"))
    # phase 1b: the period scan is now "period_obs" + N "psurr:<null>:<start>"
    # subtasks. The old single "period_scan" key is still summed so a checkpoint
    # written by the phase-1 build reads correctly in this table.
    per = (task_seconds.get("period_scan", 0.0)
           + task_seconds.get("period_obs", 0.0)
           + sum(v for k, v in task_seconds.items() if k.startswith("psurr:")))
    total = glm + marks + per
    biggest = max(task_seconds.values()) if task_seconds else 0.0
    return {"glm": glm, "marks": marks, "period_scan": per, "total": total,
            "largest_single_task": biggest,
            "amdahl_ceiling": (total / biggest) if biggest > 0 else float("inf"),
            "n_tasks": len(task_seconds)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--no-download", action="store_true", default=True,
                    help="use cached data only (the default here, and required "
                         "for a benchmark: a network fetch is not compute)")
    ap.add_argument("--download", dest="no_download", action="store_false")
    ap.add_argument("--cells", default="1seq,1par,1seqL,1parL",
                    help="comma-separated subset of " + ",".join(CELLS))
    ap.add_argument("--repeat", type=int, default=1,
                    help="repetitions per cell; the MEDIAN is reported")
    ap.add_argument("--surrogates", type=int, default=None,
                    help="override the quick preset's surrogate count")
    ap.add_argument("--data-dir", default=datasets.DEFAULT_DATA_DIR)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--scaling", default=None,
                    help="comma-separated worker counts (e.g. 1,2,4,8,16,30). "
                         "Runs the no-ladder sweep at each and prints the scaling "
                         "curve. A headline speedup at one job count says nothing "
                         "about where the curve bends, and the bend is the finding")
    a = ap.parse_args(argv)

    cells = [c.strip() for c in a.cells.split(",") if c.strip()]
    bad = [c for c in cells if c not in CELLS]
    if bad:
        ap.error(f"unknown cells {bad}; choose from {sorted(CELLS)}")

    preset = dict(ms.QUICK)
    if a.surrogates:
        preset["n_surrogates"] = int(a.surrogates)
    args_kw = {"data_dir": a.data_dir, "no_download": a.no_download}

    print("=" * 78)
    print(f"EQ-24 miner v2 throughput benchmark -- engine v{__version__}")
    print("=" * 78)
    print(f"preset          = {preset['label']} "
          f"(n_surrogates={preset['n_surrogates']}, "
          f"n_periods={preset['n_periods']}, lags={len(preset['lags'])})")
    print(f"cpu_count       = {os.cpu_count()}, auto --jobs 0 -> "
          f"{ms.resolve_jobs(0)} processes")
    print(f"data            = {a.data_dir}, downloads "
          f"{'DISABLED (cached only)' if a.no_download else 'enabled'}")
    print(f"repeats         = {a.repeat} (median reported)")
    print(f"ledger          = {BENCH_LEDGER} (SCRATCH -- the real "
          f"engine/EXPLORE_COUNT.jsonl is deliberately not touched)")
    print()
    print("NOTE: the four cells do not do the same amount of statistical work. "
          "The ladder cells stop early on purpose; the parallel cells use per-task "
          "RNG streams. This table measures WALL CLOCK only.")
    print()

    # ---- prepare ONCE, outside every timed region ----
    print("preparing design + ETAS baseline once (excluded from all timings) ...",
          flush=True)
    t_prep = time.perf_counter()
    prepared = ms.prepare(_cfg_for(cells[0], preset, args_kw), verbose=a.verbose)
    prep_s = time.perf_counter() - t_prep
    n_feats = len(prepared[7])
    print(f"  prepared in {prep_s:.1f} s: {n_feats} features, "
          f"{prepared[4].size} days in the mining window")
    print()

    rows = []
    for cell in cells:
        runs = []
        for r in range(a.repeat):
            print(f"  running {cell} ({CELLS[cell]['label']}) "
                  f"rep {r + 1}/{a.repeat} ...", flush=True)
            runs.append(_one(cell, preset, args_kw, prepared, verbose=a.verbose))
        med = statistics.median(x["seconds"] for x in runs)
        row = dict(runs[0])
        row["seconds"] = med
        row["all_seconds"] = [round(x["seconds"], 2) for x in runs]
        rows.append(row)
        print(f"    -> {med:.1f} s, {row['n_tests']} tests, "
              f"{row['n_pass']} survive FDR")

    base = next((r for r in rows if r["cell"] == "1seq"), None)
    print()
    print("=" * 78)
    print("WALL CLOCK")
    print("=" * 78)
    print(f"{'cell':<7} {'configuration':<26} {'procs':>5} {'ladder':>7} "
          f"{'seconds':>9} {'speedup':>8} {'tests':>6} {'pass':>5}")
    print("-" * 78)
    for r in rows:
        sp = (f"{base['seconds'] / r['seconds']:.2f}x"
              if base and r["seconds"] > 0 else "--")
        print(f"{r['cell']:<7} {r['label']:<26} {r['jobs_actual']:>5} "
              f"{('on' if r['ladder'] else 'off'):>7} {r['seconds']:>9.1f} "
              f"{sp:>8} {r['n_tests']:>6} {r['n_pass']:>5}")
    print("-" * 78)
    if base:
        print(f"speedup is relative to cell 1seq ({base['seconds']:.1f} s), "
              f"which is the v1 sequential audit baseline.")
    print(f"design + ETAS preparation took {prep_s:.1f} s and is EXCLUDED from every "
          f"row: it is identical work in every cell and the parallel layer does not "
          f"touch it. End-to-end wall clock for one cell is its row PLUS that "
          f"{prep_s:.1f} s.")

    st = _stages(rows[0]["task_seconds"]) if rows else None
    if st and st["total"] > 0:
        print()
        print("=" * 78)
        print("WHERE THE TIME GOES, AND WHAT THAT DOES TO THE CEILING")
        print("=" * 78)
        print(f"  measured in cell {rows[0]['cell']}, {st['n_tasks']} tasks, "
              f"per-task compute:")
        for name in ("glm", "marks", "period_scan"):
            print(f"    {name:<14s} {st[name]:8.1f} s  "
                  f"({100 * st[name] / st['total']:5.1f}% of the sweep)")
        print(f"    {'TOTAL':<14s} {st['total']:8.1f} s")
        print(f"  largest SINGLE task = {st['largest_single_task']:.1f} s "
              f"({100 * st['largest_single_task'] / st['total']:.1f}% of the sweep).")
        print(f"  AMDAHL CEILING at this task granularity = "
              f"{st['amdahl_ceiling']:.2f}x, however many workers are used, "
              f"because the pool cannot split one task.")
        if base:
            par = next((r for r in rows if r["jobs_actual"] > 1), None)
            if par and par["seconds"] > 0:
                got = base["seconds"] / par["seconds"]
                print(f"  ACHIEVED {got:.2f}x against a granularity ceiling of "
                      f"{st['amdahl_ceiling']:.2f}x on {par['jobs_actual']} "
                      f"workers.")
                if got < 0.25 * st["amdahl_ceiling"]:
                    print("  >>> The gap is NOT task granularity; look at the "
                          "stage table above for which stage now dominates. "
                          "HISTORY: through phase 1b the binding limit was the "
                          "Lomb-Scargle kernel -- scipy 1.18's `lombscargle` is "
                          "pure numpy and allocates several (n_time x n_freq) "
                          "float64 temporaries plus 2 x n_time x n_freq "
                          "transcendentals PER CALL, so the scan was "
                          "memory-bandwidth bound and aggregate throughput "
                          "saturated a few x above serial however many cores were "
                          "used. Phase 1c removed that: the cos/sin basis over "
                          "(t, period grid) is identical for the observed series "
                          "and every surrogate in a scan, so it is now built once "
                          "per worker process and each series costs one BLAS pass "
                          "(engine.mine.lomb_scargle_power). If the period stage "
                          "is no longer the largest line above, the next lever is "
                          "wherever the table says it is.")
        if st["amdahl_ceiling"] < 2.0:
            print(f"  >>> The measured parallel speedup above should be read "
                  f"against {st['amdahl_ceiling']:.2f}x, NOT against the worker "
                  f"count. The pool is not the bottleneck; the task granularity "
                  f"is. Splitting the largest task is a SCOPE DECISION, not a "
                  f"tuning knob -- see the note in this module's docstring.")
    if a.repeat > 1:
        print()
        for r in rows:
            print(f"  {r['cell']:<7} individual seconds: {r['all_seconds']}")
    print()
    print("config hashes (note --jobs is absent from all of them; --ladder is not):")
    for r in rows:
        print(f"  {r['cell']:<7} {r['config_hash']}  ladder="
              f"{'on' if r['ladder'] else 'off'}")

    # ---- the jobs scaling curve -------------------------------------------
    # Deliberately the SAME sweep at every point, ladder off, so the only thing
    # varying is the worker count. Reported with parallel efficiency, because a
    # speedup number without efficiency hides exactly where the curve stops paying.
    scaling = []
    if a.scaling:
        js = [int(x) for x in a.scaling.split(",") if x.strip()]
        print()
        print("=" * 78)
        print(f"JOBS SCALING CURVE (no ladder, {a.repeat} rep(s) per point, "
              f"median reported)")
        print("=" * 78)
        for j in js:
            runs = []
            for r in range(a.repeat):
                print(f"  jobs={j} rep {r + 1}/{a.repeat} ...", flush=True)
                CELLS["_scale"] = {"jobs": j, "ladder": False,
                                   "label": f"--jobs {j}"}
                try:
                    runs.append(_one("_scale", preset, args_kw, prepared,
                                     verbose=a.verbose))
                finally:
                    CELLS.pop("_scale", None)
            med = statistics.median(x["seconds"] for x in runs)
            scaling.append({"jobs_flag": j,
                            "jobs_actual": ms.resolve_jobs(j),
                            "seconds": med,
                            "all_seconds": [round(x["seconds"], 2) for x in runs],
                            "n_tests": runs[0]["n_tests"],
                            "n_pass": runs[0]["n_pass"],
                            "task_seconds_total": sum(
                                runs[0]["task_seconds"].values())})
            print(f"    -> {med:.2f} s")
        s1 = next((r["seconds"] for r in scaling if r["jobs_actual"] == 1), None)
        print()
        print(f"{'jobs':>6} {'seconds':>10} {'speedup':>9} {'efficiency':>11} "
              f"{'compute-s':>10} {'tests':>6}")
        print("-" * 60)
        for r in scaling:
            sp = (s1 / r["seconds"]) if s1 and r["seconds"] > 0 else float("nan")
            eff = sp / r["jobs_actual"] if r["jobs_actual"] else float("nan")
            print(f"{r['jobs_actual']:>6} {r['seconds']:>10.2f} {sp:>8.2f}x "
                  f"{100 * eff:>10.1f}% {r['task_seconds_total']:>10.1f} "
                  f"{r['n_tests']:>6}")
        print("-" * 60)
        print("compute-s is the SUM of per-task seconds measured inside the "
              "workers; wall clock below it is the parallel gain, wall clock "
              "above it is dispatch, pickling and per-process basis construction.")

    payload = {"engine_version": __version__, "preset": preset,
               "scaling": scaling,
               "cpu_count": os.cpu_count(), "repeat": a.repeat,
               "prepare_seconds_excluded": round(prep_s, 2),
               "n_features": n_feats,
               "banner": M.GENERATOR_NOT_EVIDENCE, "rows": rows}
    out_path = os.path.join(M.MINE_DIR, "bench_v2.json")
    os.makedirs(M.MINE_DIR, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    print(f"\nwritten -> {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
