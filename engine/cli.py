"""engine CLI:  python -u -m engine.cli {cache,list,run}   (run from replication/)"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
import time

import numpy as np

from . import (BASELINE_CAVEAT, __version__, baseline as bl, baseline_caveat,
               canonical_baseline, covariates, datasets, design, score, splits)

OUT_DIR = os.path.join("engine", "out")

BASELINE_CHOICES = ["climatology", "climatology-v1", "etas", "etas-v1"]


def _banner(title):
    print("=" * 78)
    print(title)
    print("=" * 78)


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
    print("baselines:")
    print(f"  climatology-v1  {BASELINE_CAVEAT}")
    print(f"  etas            {baseline_caveat('etas')}")
    print("                  params cached in engine/out/cache/etas_params.json; "
          "(re)fit with `cli fit-etas`")


# -------------------------------------------------------------------- run ---
def _build_config(args):
    return {
        "engine_version": __version__,
        "covariate": args.covariate,
        "params": json.loads(args.params) if args.params else {},
        "baseline": canonical_baseline(args.baseline),
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
    base = bl.make_baseline(cfg["baseline"], alpha=cfg.get("alpha", 0.5),
                            mag_target=cfg["mag_target"], verbose=verbose)
    t_bl = time.time()
    base.fit(ctx, y, explore)
    if verbose:
        for line in base.report():
            print(line)
        print(f"  baseline ready in  = {time.time()-t_bl:.1f}s")

    t0 = time.time()
    z, zstats = covariates.compute(cfg["covariate"], ctx, cfg["params"])
    burn = max(zstats["burn_in"], getattr(base, "burn_in_days", 0))
    if verbose:
        print(f"covariate            = {cfg['covariate']} {zstats['params']}")
        print(f"  computed in        = {time.time()-t0:.1f}s")
        print(f"  burn-in days       = {burn} (dropped from both fit and eval; "
              f"covariate {zstats['burn_in']}, baseline "
              f"{getattr(base, 'burn_in_days', 0)})")
        print(f"  standardised range = {zstats['min']:.3f} .. {zstats['max']:.3f} "
              f"(raw mean {zstats['raw_mean']:.4f}, sd {zstats['raw_sd']:.4f})")

    fit_sl = slice(burn, explore.stop)
    eval_sl = slice(burn, explore.stop) if mode == "explore" else holdout

    mu_fit = base.rate(fit_sl)          # 1-D for climatology, 2-D for etas
    mu_eval = base.rate(eval_sl)
    fit = score.fit_poisson(mu_fit, z[:, fit_sl], y[:, fit_sl])
    ig = score.information_gain(mu_eval, z[:, eval_sl], y[:, eval_sl], fit["a"], fit["beta"])
    mol = score.molchan(mu_eval, z[:, eval_sl], y[:, eval_sl], fit["a"], fit["beta"])
    gain, nev = score.green_red_cells(
        ctx, mu_eval, z[:, eval_sl], y[:, eval_sl], fit["a"], fit["beta"], ig["a0_baseline"]
    )
    results = {
        "mode": mode,
        "baseline": {"name": base.name, "caveat": base.caveat,
                     **({"params": base.params, "fit_info":
                         {k: v for k, v in base.fit_info.items() if k != "bounds"}}
                        if isinstance(base, bl.EtasV1) else {})},
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
    print(f"RESULT  covariate={cfg['covariate']}  mode={res['mode']}  "
          f"baseline={res.get('baseline', {}).get('name', cfg['baseline'])}")
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
    print(baseline_caveat(cfg["baseline"]))
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
        {"config": cfg, "baseline_caveat": baseline_caveat(cfg["baseline"]),
         "mode": args.mode},
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
                   "baseline_caveat": baseline_caveat(cfg["baseline"])}, fh, indent=1)
    print(f"  full result JSON     = {res_path}")
    return 0


# --------------------------------------------------------------- fit-etas ---
def cmd_fit_etas(args):
    _banner(f"EQ-23 engine v{__version__} -- fit ETAS baseline")
    ctx = design.build_design(
        data_dir=args.data_dir, dlat=args.dlat, dlon=args.dlon,
        explore_frac=args.explore_frac, verbose=True,
    )
    explore, _hold = splits.temporal_split(ctx.n_days, args.explore_frac)
    y = ctx.day_counts(args.mag_target)

    base = bl.EtasV1(alpha=args.alpha, verbose=True, polish=not args.no_polish,
                     trunc_days=args.omori_trunc, mag_target=args.mag_target)
    if args.refit and os.path.exists(base.cache_path):
        print(f"--refit: discarding cached params in {base.cache_path}")
        os.remove(base.cache_path)

    print("fitting ETAS by Poisson maximum likelihood on the EXPLORATION window only.")
    print(f"  Omori truncation     = {args.omori_trunc} d (declared approximation)")
    print(f"  bounds               = "
          + ", ".join(f"{k}{list(bl.PARAM_BOUNDS[k])}" for k in bl.PARAM_NAMES))
    t0 = time.time()
    base.fit(ctx, y, explore)
    print(f"  wall clock           = {time.time()-t0:.1f}s "
          f"(params {base.fit_info.get('source')})")
    print("-" * 78)
    for line in base.report():
        print(line)
    info = base.fit_info
    if "n_objective_evals" in info:
        print(f"  objective evals      = {info['n_objective_evals']} "
              f"({', '.join(info['methods'])})")
    print(f"  params cache         = {base.cache_path}")
    print(baseline_caveat("etas"))
    print("-" * 78)
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
    r.add_argument("--baseline", default="climatology-v1", choices=BASELINE_CHOICES)
    r.add_argument("--mode", choices=["explore", "holdout"], default="explore")
    r.add_argument("--config", default=None, help="frozen config JSON (required for holdout)")
    r.add_argument("--mag-target", type=float, default=4.5)
    r.add_argument("--alpha", type=float, default=0.5)
    r.set_defaults(fn=cmd_run)

    f = sub.add_parser("fit-etas", help="fit (or reload) the ETAS baseline parameters")
    common(f)
    f.add_argument("--mag-target", type=float, default=4.5)
    f.add_argument("--alpha", type=float, default=0.5,
                   help="climatology smoothing used for the ETAS background shape")
    f.add_argument("--omori-trunc", type=int, default=bl.OMORI_TRUNC_DAYS)
    f.add_argument("--refit", action="store_true", help="discard the cached params first")
    f.add_argument("--no-polish", action="store_true",
                   help="skip the Nelder-Mead polish after L-BFGS-B")
    f.set_defaults(fn=cmd_fit_etas)

    args = p.parse_args(argv)
    rc = args.fn(args)
    return rc if isinstance(rc, int) else 0


if __name__ == "__main__":
    sys.exit(main())
