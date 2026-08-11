"""engine CLI:  python -u -m engine.cli {cache,list,run}   (run from replication/)"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
import time

import numpy as np

from . import BASELINE_CAVEAT, __version__, baseline as bl, covariates, datasets, design, score, splits

OUT_DIR = os.path.join("engine", "out")


def _banner(title):
    print("=" * 78)
    print(title)
    print("=" * 78)


def _caveat(baseline_name):
    print(f"baseline={baseline_name} (clustering NOT absorbed; sniff-grade only)"
          if baseline_name == "climatology-v1" else f"baseline={baseline_name}")


# ------------------------------------------------------------------ cache ---
def cmd_cache(args):
    _banner(f"EQ-23 engine v{__version__} -- design cache")
    t0 = time.time()
    ctx = design.build_design(
        data_dir=args.data_dir, dlat=args.dlat, dlon=args.dlon,
        explore_frac=args.explore_frac, rebuild=args.rebuild, verbose=True,
    )
    print(f"cache ready in {time.time()-t0:.1f}s")
    print(BASELINE_CAVEAT)
    return ctx


# ------------------------------------------------------------------- list ---
def cmd_list(args):
    _banner("registered covariates")
    for name, fn in covariates.available().items():
        print(f"  {name:<14s} {fn.describe}")
        print(f"                 defaults={fn.defaults}")
    print()
    print("baselines: " + ", ".join(sorted(bl.BASELINES)) + "   (etas = NotImplementedError)")
    print(BASELINE_CAVEAT)


# -------------------------------------------------------------------- run ---
def _build_config(args):
    return {
        "engine_version": __version__,
        "covariate": args.covariate,
        "params": json.loads(args.params) if args.params else {},
        "baseline": args.baseline,
        "grid": {"dlat": args.dlat, "dlon": args.dlon},
        "explore_frac": args.explore_frac,
        "mag_target": args.mag_target,
        "alpha": args.alpha,
        "data_dir": args.data_dir.replace("\\", "/"),
    }


def _run_core(cfg, mode, verbose=True):
    ctx = design.build_design(
        data_dir=cfg["data_dir"], dlat=cfg["grid"]["dlat"], dlon=cfg["grid"]["dlon"],
        explore_frac=cfg["explore_frac"], verbose=verbose,
    )
    explore, holdout = splits.temporal_split(ctx.n_days, cfg["explore_frac"])

    y = ctx.day_counts(cfg["mag_target"])
    cls = bl.get_baseline(cfg["baseline"])
    base = cls(alpha=cfg.get("alpha", 0.5)) if cls is bl.ClimatologyV1 else cls()
    base.fit(ctx, y, explore)
    if verbose:
        for line in base.report():
            print(line)

    t0 = time.time()
    z, zstats = covariates.compute(cfg["covariate"], ctx, cfg["params"])
    burn = zstats["burn_in"]
    if verbose:
        print(f"covariate            = {cfg['covariate']} {zstats['params']}")
        print(f"  computed in        = {time.time()-t0:.1f}s")
        print(f"  burn-in days       = {burn} (dropped from both fit and eval)")
        print(f"  standardised range = {zstats['min']:.3f} .. {zstats['max']:.3f} "
              f"(raw mean {zstats['raw_mean']:.4f}, sd {zstats['raw_sd']:.4f})")

    fit_sl = slice(burn, explore.stop)
    eval_sl = slice(burn, explore.stop) if mode == "explore" else holdout

    fit = score.fit_poisson(base.mu, z[:, fit_sl], y[:, fit_sl])
    ig = score.information_gain(base.mu, z[:, eval_sl], y[:, eval_sl], fit["a"], fit["beta"])
    mol = score.molchan(base.mu, z[:, eval_sl], y[:, eval_sl], fit["a"], fit["beta"])
    gain, nev = score.green_red_cells(
        ctx, base.mu, z[:, eval_sl], y[:, eval_sl], fit["a"], fit["beta"], ig["a0_baseline"]
    )
    results = {
        "mode": mode,
        "fit": fit,
        "covariate_stats": {k: v for k, v in zstats.items() if k != "params"},
        "fit_window_days": [fit_sl.start, fit_sl.stop],
        "eval_window_days": [eval_sl.start, eval_sl.stop],
        "n_events_fit": float(y[:, fit_sl].sum()),
        "information_gain": {k: v for k, v in ig.items()},
        "molchan": {k: v for k, v in mol.items() if k not in ("tau", "nu")},
        "molchan_trajectory": {"tau": mol["tau"], "nu": mol["nu"]},
        "in_sample": mode == "explore",
    }
    return ctx, results, gain, nev


def _print_results(cfg, res, n_green, n_red, extra=()):
    print("-" * 78)
    print(f"RESULT  covariate={cfg['covariate']}  mode={res['mode']}")
    f = res["fit"]
    ig = res["information_gain"]
    print(f"  beta                 = {f['beta']:+.5f} +/- {f['se_beta']:.5f} "
          f"(z = {f['z_stat']:+.2f})")
    print(f"  converged            = {f['converged']} in {f['n_iter']} Newton steps")
    print(f"  intercept a          = {f['a']:+.5f}   (nuisance re-normalisation)")
    print(f"  fit window (days)    = {res['fit_window_days']}  n_events={res['n_events_fit']:.0f}")
    print(f"  eval window (days)   = {res['eval_window_days']}  n_events={ig['n_events_eval']:.0f}")
    print(f"  BITS/EVENT           = {ig['bits_per_event']:+.6f}   "
          f"(total {ig['bits_total']:+.1f} bits)")
    m = res["molchan"]
    print(f"  Molchan skill        = {m['molchan_skill']:+.4f} "
          f"(0 = no skill; miss rate at tau=0.10 is {m['nu_at_tau_0.10']:.3f})")
    print(f"  green/red cells      = {n_green} green / {n_red} red")
    if res["in_sample"]:
        print("  NOTE: mode=explore -> beta is fitted AND evaluated on the exploration")
        print("        period. This bits/event is IN-SAMPLE and is an upper bound.")
    for line in extra:
        print(line)
    print(BASELINE_CAVEAT)
    print("-" * 78)


def cmd_run(args):
    cfg = _build_config(args)
    runid = _dt.datetime.now().strftime("%Y%m%dT%H%M%S") + "_" + splits.config_hash(cfg)[:8]
    _banner(f"EQ-23 engine v{__version__} -- run {runid}  mode={args.mode}")
    if args.mode != "holdout":
        print(f"config hash          = {splits.config_hash(cfg)}")

    if args.mode == "holdout":
        if not args.config:
            print("ERROR: --mode holdout requires --config <file.json>. The holdout is spent")
            print("       on first contact; the configuration must be frozen in a file first.")
            return 2
        with open(args.config, "r", encoding="utf-8") as fh:
            cfg = json.load(fh)
        h = splits.config_hash(cfg)
        print(f"frozen config        = {args.config}")
        print(f"config hash          = {h}")
        try:
            splits.assert_holdout_unused(cfg)
        except splits.HoldoutRefused as e:
            print()
            print(str(e))
            return 3
        n_prior = splits.n_explore_runs_since_last_holdout()
        print(f"n_explore_runs since last holdout = {n_prior}  (reported multiplicity)")

    ctx, res, gain, nev = _run_core(cfg, args.mode)

    os.makedirs(OUT_DIR, exist_ok=True)
    cells_path = os.path.join(OUT_DIR, f"cells_{runid}.json")
    n_green, n_red = score.export_cells(
        cells_path, ctx, gain, nev, runid,
        {"config": cfg, "baseline_caveat": BASELINE_CAVEAT, "mode": args.mode},
    )

    extra = []
    if cfg["covariate"] == "anisotropy":
        from .covariates.anisotropy import CAVEAT
        extra.append("  " + CAVEAT)

    if args.mode == "holdout":
        summary = {
            "beta": res["fit"]["beta"], "se_beta": res["fit"]["se_beta"],
            "bits_per_event": res["information_gain"]["bits_per_event"],
            "n_events_eval": res["information_gain"]["n_events_eval"],
            "molchan_skill": res["molchan"]["molchan_skill"],
            "runid": runid,
        }
        rec = splits.record_holdout(cfg, summary)
        extra.append(f"  HOLDOUT SPENT: appended to {splits.HOLDOUT_LOG}")
        extra.append(f"  n_explore_runs since last holdout = {rec['n_explore_runs']} "
                     f"(reported multiplicity)")
        res["n_explore_runs"] = rec["n_explore_runs"]
    else:
        splits.log_explore_run(cfg)
        extra.append(f"  explore run logged to {splits.EXPLORE_COUNT}")

    _print_results(cfg, res, n_green, n_red, extra)
    print(f"  green/red export     = {cells_path}")

    res_path = os.path.join(OUT_DIR, f"run_{runid}.json")
    with open(res_path, "w", encoding="utf-8") as fh:
        json.dump({"runid": runid, "config": cfg, "results": res,
                   "baseline_caveat": BASELINE_CAVEAT}, fh, indent=1)
    print(f"  full result JSON     = {res_path}")
    return 0


def main(argv=None):
    p = argparse.ArgumentParser("engine.cli", description=f"EQ-23 engine v{__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--data-dir", default=datasets.DEFAULT_DATA_DIR)
        sp.add_argument("--dlat", type=float, default=1.0)
        sp.add_argument("--dlon", type=float, default=1.0)
        sp.add_argument("--explore-frac", type=float, default=0.70)

    c = sub.add_parser("cache", help="build/refresh the design cache")
    common(c)
    c.add_argument("--rebuild", action="store_true")
    c.set_defaults(fn=cmd_cache)

    l = sub.add_parser("list", help="list registered covariates")
    l.set_defaults(fn=cmd_list)

    r = sub.add_parser("run", help="fit + score one covariate")
    common(r)
    r.add_argument("--covariate", required=True)
    r.add_argument("--params", default=None, help="JSON dict of covariate params")
    r.add_argument("--baseline", default="climatology-v1")
    r.add_argument("--mode", choices=["explore", "holdout"], default="explore")
    r.add_argument("--config", default=None, help="frozen config JSON (required for holdout)")
    r.add_argument("--mag-target", type=float, default=4.5)
    r.add_argument("--alpha", type=float, default=0.5)
    r.set_defaults(fn=cmd_run)

    args = p.parse_args(argv)
    rc = args.fn(args)
    return rc if isinstance(rc, int) else 0


if __name__ == "__main__":
    sys.exit(main())
