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


# ------------------------------------------------------------------- mine ---
def cmd_mine(args):
    from . import mine as M, mine_session as ms

    preset = ms.QUICK if args.quick else (ms.OVERNIGHT if args.overnight else ms.DEFAULT)
    cfg = ms.build_config(args, preset)
    if args.surrogates:
        cfg["n_surrogates"] = int(args.surrogates)
    _banner(f"EQ-23 engine v{__version__} -- MINE (EQ-24 pattern miner) "
            f"preset={cfg['preset']}"
            + (f" -- {M.TRANCHE1_LABEL}" if cfg["tranche1"] else ""))
    print("*** " + M.GENERATOR_NOT_EVIDENCE)
    print()
    print("standing warning (EQ-24, verbatim):")
    print("  " + M.MIGNAN_BROCCARDO)
    print()
    print(f"config hash          = {splits.config_hash(cfg)}")
    print(f"surrogates           = {cfg['n_surrogates']} (circular shifts; capped at "
          f"the number of admissible shifts and reported per test)")
    print(f"lag grid             = {cfg['lags'][0]}..{cfg['lags'][-1]} "
          f"({len(cfg['lags'])} lags, aperiodic features only)")
    if cfg["tranche1"]:
        tl = cfg["tranche1_lags"]
        print(f"tranche              = {M.TRANCHE1_LABEL}: cyclic lag grid "
              f"{tl[0]}..{tl[-1]} d ({len(tl)} lags) on the 13 lag-unscanned cyclic "
              f"features")
        print(f"                       = {13 * (len(tl) - 1)} NEW tests "
              f"(8 linear x {len(tl)-1} = {8*(len(tl)-1)}; 5 non-lag-free phase x "
              f"{len(tl)-1} = {5*(len(tl)-1)}); the 4 provably lag-free phase "
              f"features stay at lag 0")
    print(f"FDR                  = Benjamini-Hochberg at q={cfg['fdr_q']}")
    n_jobs = ms.resolve_jobs(args.jobs)
    print(f"jobs                 = {n_jobs} process(es) "
          + ("(sequential; shared rng stream = the v1 audit baseline)" if n_jobs == 1
             else "(parallel; per-task derived SeedSequence, order-independent, "
                  "NOT bit-identical to --jobs 1). Not in the config hash."))
    if cfg.get("ladder"):
        L = cfg["ladder"]
        print(f"ladder               = ON, {L['rule']} h={L['h']} "
              f"N_max={L['n_max']} (chunk={L['chunk']}, an implementation detail) "
              f"-- IN the config hash (statistical choice). Applies to Monte Carlo "
              f"nulls only, never to the exhaustive shift enumeration.")
    else:
        print(f"ladder               = off (full fixed surrogate budget per test)")

    out = ms.run(cfg, verbose=True, resume=not args.new_session, jobs=args.jobs)
    print("-" * 78)
    print(f"MINE SESSION COMPLETE  {out['session_dir']}")
    print(f"  tests                = {out['n_tests']}")
    print(f"  survive BH-FDR       = {out['n_pass']}")
    print(f"  report               = {out['report']}")
    print(f"  stubs                = {out['stubs']}")
    print(f"  elapsed              = {out['elapsed']} s")
    print(baseline_caveat("etas"))
    print("*** " + M.GENERATOR_NOT_EVIDENCE)
    print("-" * 78)
    return 0


# ----------------------------------------------------------------- search ---
# B5 / SEARCHER.md §S8.2: `search` mode alongside `explore` / `holdout` / `mine`,
# with the same GENERATOR_NOT_EVIDENCE banner `mine.py:84` already prints AND a
# second banner naming §S0's sentence. ADDITIVE ONLY: nothing above this block was
# modified, and no existing subcommand's behaviour changes by its existence.
#
# THE REFUSAL THIS SUBCOMMAND EXISTS TO CARRY. §P7-24 SP-7 is BINDING: *"No real scan
# runs before this passes."* `--run` therefore reads the SP-7 gate artifact and
# REFUSES when it is absent or FAILED. There is no flag that bypasses it, on the same
# principle as `splits.assert_holdout_unused`: deleting the artifact is the human's
# act, not the engine's.
def cmd_search(args):
    from . import (freeze_gen as FG, lattice_s1 as LAT, properties as PR,
                   searcher as SE, searcher_gate as SG, spaceweather as SW)

    _banner(f"EQ-23 engine v{__version__} -- SEARCH (the SEARCHER, §S0..§S11)")
    SE.print_banners()

    l3 = None
    print(f"lattice rule id      = {LAT.LATTICE_RULE_ID}")
    print(f"lattice digest       = {LAT.declaration_digest(l3)}")
    print(f"  L1 named regions   = {LAT.N_L1}")
    print(f"  L2 tectonic classes= {LAT.N_L2} {LAT.BIRD_CLASSES}")
    print(f"  L3 grid cells      = {len(LAT.l3_grid_cells())} declared at "
          f"{LAT.L3_CELL_DEG:g} deg (measurability filter N >= {LAT.L3_MIN_EVENTS} "
          f"applied on the EXPLORATION WINDOW ONLY)")
    print(f"  seeded exclusions  = {LAT.SEEDED_REGIONS}")
    print()

    av = SW.availability(args.data_dir_sw)
    print("space-weather / EOP series (families D and E):")
    for k, v in sorted(av["sources"].items()):
        print(f"  {k:<10s} {v['status']}"
              + (f"  n={v['n']}" if "n" in v else "")
              + (f"  err={v.get('error')}" if v["status"] != "LOADED" else ""))
    if av["series_absent"]:
        print(f"  ABSENT series      = {av['series_absent']} "
              f"(no column is emitted for an absent source; nothing is substituted)")
    print()

    decl = SE.declare_scan(
        [LAT.region(n) for n in SG.GATE_REGIONS], SG.GATE_PROPERTIES,
        SG.GATE_MAG_STRATA, SG.GATE_FAMILIES, l3_cells=l3,
        scan_id=args.scan_id) if args.declare_only or not args.run else None
    if decl is not None:
        print(f"DECLARATION (dry run -- nothing scanned)")
        print(f"  scan id            = {decl['scan_id']}")
        print(f"  m (FULL cell count)= {decl['m']}  -- ASSERTED at scan time (SP-6.1)")
        print(f"  alpha = q/m        = {decl['alpha']:.6e}  (q = {decl['q']}) "
              f"-- ONE rule, the OR-limb is REFUSED (SP-3)")
        print(f"  K cap              = {decl['k_cap']}   ranking = {decl['ranking_function']}")
        print(f"  n_events floor     = {decl['n_events_floor']}")
        print(f"  magnitude strata   = {decl['mag_strata']} (SP-6.7)")
        print(f"  config hash        = {decl['config_hash']}")
        print()

    gate_path = args.gate_artifact or SG.GATE_ARTIFACT
    gate = None
    if os.path.exists(gate_path):
        with open(gate_path, "r", encoding="utf-8") as fh:
            gate = json.load(fh)
    print("SP-7 GATE (BINDING -- no real scan runs before it passes):")
    if gate is None:
        print(f"  NOT RUN. Artifact absent at {gate_path}.")
        print(f"  Run it with:  python -u -m engine.searcher_gate")
    else:
        print(f"  verdict            = {gate['verdict']}")
        print(f"  promotions         = {gate['n_promotions_total']} over "
              f"{gate['n_catalogs']} catalogues (allowed <= {gate['max_total_allowed']})")
        print(f"  binomial tail p    = {gate['binomial_upper_tail_p']:.4f}")
        print(f"  vacuity control    = {gate['vacuity_control']['verdict']}")
        print(f"  artifact           = {gate_path}")
    print()

    if args.run:
        if gate is None or not gate.get("passed"):
            print("REFUSED: a real scan requires a PASSING SP-7 gate artifact.")
            print("  §P7-24 SP-7: \"No real scan runs before this passes. Same "
                  "standing as R1 for v2; same reason.\"")
            print("  A failure there means the nominal p's are wrong -- an SP-2 null "
                  "layer is invalid -- which is precisely the failure mode the "
                  "multiplicity arithmetic cannot catch.")
            print("  No flag bypasses this refusal.")
            return 2
        print("REFUSED: `--run` is not wired to a catalogue in this build.")
        print("  The scan runner (engine/searcher.py:run_scan) takes prepared cell "
              "specs; building them from a real catalogue is the operator step that "
              "must be declared and hashed first, and this CLI will not improvise it.")
        return 2

    print("*** " + SE.THE_SENTENCE)
    print("*** " + FG.PROMOTION_BUYS)
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

    mi = sub.add_parser("mine", help="EQ-24 pattern miner (GENERATOR, NOT EVIDENCE)")
    common(mi)
    mi.add_argument("--mag-target", type=float, default=4.5)
    mi.add_argument("--quick", action="store_true",
                    help="small grid, 200 surrogates (<10 min)")
    mi.add_argument("--overnight", action="store_true",
                    help="full lag grid + 10k surrogates, checkpointed per feature")
    mi.add_argument("--tranche1", action="store_true",
                    help="K-089-R tranche 1: scan lags 0..30 d for the 13 "
                         "lag-unscanned cyclic features (8 linear + 5 non-lag-free "
                         "phase); +390 tests, priced at the declared count")
    mi.add_argument("--controls", action="store_true",
                    help="F9-20 arm 1 (Tranche A, §P7-2): add the negative-control "
                         "feature battery -- 3 named mechanism-free lines + 20 "
                         "matched null surrogates x 31 lags = 713 PRICED tests in "
                         "their own declared stratum. GLM only; no mark or regional "
                         "axis. Hash-affecting.")
    mi.add_argument("--surrogates", type=int, default=None,
                    help="override the preset surrogate count")
    mi.add_argument("--no-download", action="store_true",
                    help="skip family-3 optional downloads entirely")
    mi.add_argument("--new-session", action="store_true",
                    help="do not resume an incomplete session with the same config")
    mi.add_argument("--seed", type=int, default=20260811)
    mi.add_argument("--jobs", type=int, default=1,
                    help="worker processes for the per-feature test loop. 1 (the "
                         "default) is sequential and bit-identical to v1; 0 is auto "
                         "= max(1, cpu_count()-2). EXECUTION ONLY -- not part of the "
                         "config hash, and parallel runs use per-task derived RNG "
                         "streams (see engine/mine_session.py docstring)")
    mi.add_argument("--ladder", action="store_true",
                    help="adaptive surrogate ladder: Besag-Clifford sequential "
                         "stopping on the block-bootstrap null. Off by default. "
                         "CHANGES THE CONFIG HASH when on, because it is a "
                         "statistical choice")
    mi.add_argument("--ladder-h", type=int, default=None,
                    help="exceedances at which a laddered test stops (default 25)")
    mi.add_argument("--ladder-chunk", type=int, default=None,
                    help="vectorisation chunk size for ladder draws (default 200). "
                         "An IMPLEMENTATION DETAIL, not a statistical stage: the "
                         "exact stopping index is recovered inside the chunk, so "
                         "changing this changes speed and nothing else")
    mi.add_argument("--gpd", action="store_true",
                    help="§P6-2 GPD tail extrapolation for p-values CENSORED at the "
                         "Monte Carlo floor: MLE fit to the exceedances over the "
                         "declared top 10%%, gated on Anderson-Darling p >= 0.05 AND "
                         "xi-stability across the top {5,10,20}%%, quoted at the "
                         "UPPER end of the 95%% bootstrap CI, capped one decade "
                         "below 1/(N_max+1), and FORBIDDEN on the all-shifts "
                         "enumeration null. Off by default. CHANGES THE CONFIG HASH")
    mi.add_argument("--gpd-confirm-max-n", type=int, default=None,
                    help="ceiling on the targeted brute-force confirmation of a GPD "
                         "candidate (§P6-2(6) demands N >= 10/p_gpd). Default "
                         "500000. A candidate needing more stays a candidate and "
                         "emits NO stub")
    mi.add_argument("--gpd-calibration", default=None,
                    help="path to the §P6-2(7) calibration artifact (default "
                         "engine/out/audit_gpd_bca.json). --gpd REFUSES TO RUN without "
                         "a PASSING one: an instrument that has not demonstrated "
                         "recovery in the range it reports in may not report there")
    mi.add_argument("--strata", default=None,
                    help="§P6-3 stratified (weighted) BH: path to a declared "
                         "partition JSON, (feature_family x test_kind [x region]) "
                         "with a declared m_s and q_s per stratum. The engine "
                         "ASSERTS sum_s m_s q_s == m q and REFUSES TO RUN on "
                         "violation. Default UNSTRATIFIED (flat BH). The FILE "
                         "CONTENT is hash-affecting")
    mi.add_argument("--regsum", action="store_true",
                    help="§P6-4 Rule 4.2: the PHASE-INCOHERENT regional statistic. "
                         "Partitions the domain into R declared regions from "
                         "EXPLORATION-WINDOW data only (Rule 4.1), computes the "
                         "per-region 2-df score statistic and SUMS it to 2R df -- "
                         "one test per (feature, lag), which kills the "
                         "§K87-0(d)(i) domain-sum cancellation blind spot at zero "
                         "multiplicity cost. Per-region amplitudes are reported "
                         "UNRESOLVED (§P7-1(d)). Hash-affecting")
    mi.add_argument("--regions", action="store_true",
                    help="§P6-4 Rule 4.3/4.4: ALSO run the per-region battery, one "
                         "test per (feature, lag, region), priced at R x the "
                         "declared count, with region as a stratum axis and an "
                         "S-15 measurable/unmeasurable verdict per region from the "
                         "§P7-1(b) FORMULA floor. Implies --regsum. Hash-affecting")
    mi.add_argument("--region-sectors", type=int, default=None,
                    help="R before the activity threshold (default 6, chosen from "
                         "the F4-58 VIF measurement -- see results_f4_58_vif.json)")
    mi.add_argument("--region-min-fraction", type=float, default=None,
                    help="a longitude sector is retained only if it holds at least "
                         "this fraction of in-window events (default 0.01)")
    mi.add_argument("--region-vif", type=float, default=None,
                    help="declared VIF for the §P7-1(b) floor formula. Default is "
                         "the F4-58 MEASUREMENT (24.08, 2-df phase features), not "
                         "the 3.94 §P7-1(a) inferred. Hash-affecting")
    mi.add_argument("--region-alpha", type=float, default=None,
                    help="declared operating threshold of the tranche this runs in "
                         "(default 0.10/713 = Tranche A at BH q=0.10). "
                         "Hash-affecting")
    mi.add_argument("--region-target-amplitude", type=float, default=None,
                    help="declared effect of interest, as a rate-modulation "
                         "fraction (default 0.20 -- §P6-4 Finding B's own reference "
                         "amplitude). Hash-affecting")
    mi.set_defaults(fn=cmd_mine)

    # B5: the SEARCHER's operator entry point. Purely additive.
    se = sub.add_parser("search", help="the SEARCHER concentration scan "
                                       "(GENERATOR, NOT EVIDENCE; SP-7-gated)")
    se.add_argument("--data-dir-sw", default=None,
                    help="space-weather / EOP directory (default data/spaceweather)")
    se.add_argument("--scan-id", default=None)
    se.add_argument("--gate-artifact", default=None,
                    help="path to the SP-7 gate artifact (default "
                         "engine/out/searcher_sp7_gate.json)")
    se.add_argument("--declare-only", action="store_true",
                    help="print the frozen declaration and stop (the default "
                         "behaviour without --run)")
    se.add_argument("--run", action="store_true",
                    help="attempt a real scan. REFUSED without a PASSING SP-7 gate "
                         "artifact (§P7-24 SP-7 is binding); no flag bypasses it")
    se.set_defaults(fn=cmd_search)

    args = p.parse_args(argv)
    # §P7-8(c)(4): an unsupported --mag-target RAISES here, before any design is
    # built or any session directory is created. It must never be silently clamped
    # to the catalogue floor -- see engine/datasets.assert_mag_supported for the
    # session-duplication incident that earned this check.
    if getattr(args, "mag_target", None) is not None:
        datasets.assert_mag_supported(args.mag_target, what="--mag-target")
    rc = args.fn(args)
    return rc if isinstance(rc, int) else 0


if __name__ == "__main__":
    sys.exit(main())
