"""F4-58 item (ii): the VIF ESTIMATOR CONTROL, required by HYPOTHESIS_LEDGER.md §P7-8(b).

WHAT THIS IS
------------
§P7-8(a) measured VIF = 24.08 over the 2-df phase features of
`session_20260811T022953` using the §P7-1(c) identity

    VIF(test) = chi2_obs / chi2.ppf(1 - p_boot, df)

§P7-8(b) then ruled that the measurement, on its own, **cannot distinguish two
readings with opposite consequences**:

  * Reading A -- the EARTH is overdispersed. Daily global M>=4.5 counts carry
    residual clustering the ETAS offset does not absorb. VIF = 24 is real, the
    S-15 declaration floor is real, and the programme is simply less powerful
    than it believed.
  * Reading B -- OUR OWN NULL is too wide. `Feature.block_days` clips to
    [30, 800] days; at 800 d over a 7,716-day window the circular moving-block
    bootstrap has ~10 blocks, and a null built from ten blocks is enormously
    variable BY CONSTRUCTION. Then the inflation is our null and not the Earth,
    the fix is a better null rather than a bigger floor, and the miner has been
    UNDER-POWERED rather than over-confident.

This script is the positive control that separates them. It simulates catalogs
from the engine's OWN fitted ETAS lambda where **the true VIF is 1 by
construction** (independent Poisson draws given lambda: Var = mean, exactly, and
no residual clustering of any kind survives), and runs **the same estimator,
through the same engine code path, at the same block lengths**. Any VIF above 1
that comes back is manufactured by the estimator + null, not by the Earth.

    VIF_control(block) ~ 1        => Reading A: the inflation is the Earth.
    VIF_control(block) ~ VIF_real => Reading B: the inflation is the null.

PRICING -- ZERO. THIS IS A MEASUREMENT, NOT A TRANCHE.
------------------------------------------------------
Priced tests: 0. No EXPLORE_COUNT.jsonl line is written by this script, on the
same §P7-1(c) logic that priced F4-58 itself at zero: this makes **no rejection
about the Earth and enters no BH vector**. Every statistic here is computed on
SIMULATED data whose generating process is known and stated, and the object under
test is the ESTIMATOR (S-17 applied to the newest estimator in the building), not
a hypothesis about the world. There is no multiplicity to declare because there
is no claim that could be inflated by it -- a false positive here would be a
statement about our own arithmetic, which is exactly what we are trying to
provoke.

CODE-PATH IDENTITY (what is reused, verbatim)
---------------------------------------------
  * lambda      -- `engine.baseline.EtasV1` fitted by the engine's own ML fit on
                   the exploration window, then `engine.mine.build_target` for
                   the domain-summed (counts, offset) pair. The offset IS the
                   fitted lambda series; it is the simulation's ground truth.
  * features    -- `engine.mine.ephemeris_features` / `catalog_features` and the
                   family-3 parsers `_parse_gfz` / `_parse_iers_lod` run over the
                   FROZEN copies the sessions downloaded (engine/out/mine/frozen),
                   so the designs and every `Feature.block_days` are the real ones.
  * statistic   -- `engine.mine.score_stat_all_shifts` (observed chi2).
  * null        -- `engine.mine.score_stat_block_bootstrap` at the feature's own
                   `block_days`, and `engine.mine.empirical_p` for the
                   circular-shift arm, combined by the miner's own rule
                   (p = p_boot for periodic features, max(p_shift, p_boot) else).
  * VIF         -- `f4_58_vif.vif_of`, IMPORTED, not reimplemented.

The ONLY thing that changes is the target: `counts` becomes a Poisson draw from
`offset` instead of the observed catalog.

ETAS SIMULATION PATH, DECLARED
------------------------------
`engine/baseline.py` exposes no simulate/rvs entry point (checked: EtasV1 has
fit / rate / mu / report and no sampler). So the simulation is a direct Poisson
draw from the fitted lambda series -- `rng.poisson(offset)` on the same
domain-summed daily lambda the miner scores against. This is the construction
§P7-8(b) item (ii) specifies verbatim ("Simulate Poisson counts from the fitted
ETAS lambda"), and it is the construction that makes the true VIF exactly 1:
conditional on lambda the draws are independent across days, so the only
dispersion in the simulated target is Poisson dispersion, which is what the
chi2 reference distribution assumes. It deliberately does NOT re-run the ETAS
triggering cascade -- a self-exciting simulation would reintroduce clustering
and would no longer be a VIF = 1 control.

Usage:
    python -u f4_58_vif_control.py                # full run
    python -u f4_58_vif_control.py --smoke        # 3 catalogs, 300 boots
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy import stats

REPO = os.path.dirname(os.path.abspath(__file__))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from engine import baseline as bl              # noqa: E402
from engine import design as design_mod        # noqa: E402
from engine import mine as M                   # noqa: E402
from engine import splits                      # noqa: E402
from f4_58_vif import vif_of                   # noqa: E402  -- the SAME estimator

OUT_PATH = os.path.join(REPO, "results_f4_58_vif_control.json")
PREP_CACHE = os.path.join(REPO, "engine", "out", "cache", "f4_58_control_prep.npz")
# A DEDICATED ETAS parameter cache. The shared engine/out/cache/etas_params.json
# currently holds the 0.5-degree fit; refitting the 1-degree configuration through
# the default path would overwrite it. Nothing under engine/out/cache is deleted
# or modified by this script -- it only adds these two files.
ETAS_CACHE_1DEG = os.path.join(REPO, "engine", "out", "cache",
                               "etas_params_f4_58_control_1deg.json")

FROZEN = os.path.join(REPO, "engine", "out", "mine", "frozen")
REAL_RESULTS = os.path.join(REPO, "results_f4_58_vif.json")
PRIMARY_SESSION = "session_20260811T022953"

# The primary session's configuration, from
# engine/out/mine/session_20260811T022953/checkpoint.json:config.
CFG = {
    "data_dir": "data/comcat_world",
    "dlat": 1.0,
    "dlon": 1.0,
    "explore_frac": 0.70,
    "mag_target": 4.5,
    "seed": 20260811,
}

# The forced block-length ladder. Spans the clip range of `Feature.block_days`
# ([30, 800] days) and includes every block length the real features actually
# used, so the control curve is evaluated AT the real curve's abscissae.
BLOCK_LADDER = [30.0, 60.0, 64.256, 104.0, 108.0, 150.0, 200.0, 300.0, 404.0,
                500.0, 600.0, 693.24, 730.484, 800.0]

# ---------------------------------------------------------------- worker state
_W = {}


# ======================================================================= prep ==
def _family3_from_frozen(t0, n_days, lags=(0,)):
    """Family 3 covariates rebuilt from the FROZEN download copies.

    Same parsers, same gap-filling, same 1-day causal roll as
    `engine.mine.download_features` -- the network fetch is replaced by the
    frozen text the sessions themselves wrote. Returns [] if a file is absent.
    """
    feats = []
    gfz = os.path.join(FROZEN, "gfz_kp_ap_f107.txt")
    iers = os.path.join(FROZEN, "iers_lod.txt")
    if os.path.exists(gfz):
        with open(gfz, "r", encoding="utf-8", errors="replace") as fh:
            ap, f107 = M._parse_gfz(fh.read(), t0, n_days)
        for nm, arr in (("Ap_geomagnetic", ap), ("F107_solar_flux", f107)):
            s, _cov = M._fill_gaps(arr, nm)
            s = np.roll(s, 1)
            s[0] = s[1]
            feats.append(M.Feature(nm, 3, "linear", s, "frozen copy", lags=lags))
    if os.path.exists(iers):
        with open(iers, "r", encoding="utf-8", errors="replace") as fh:
            lod = M._parse_iers_lod(fh.read(), t0, n_days)
        s, _cov = M._fill_gaps(lod, "length_of_day")
        s = np.roll(s, 1)
        s[0] = s[1]
        feats.append(M.Feature("length_of_day", 3, "linear", s, "frozen copy",
                               lags=lags))
    return feats


def prepare(verbose=True):
    """(offset, feature_records, meta). `offset` IS the fitted ETAS lambda series."""
    if os.path.exists(PREP_CACHE):
        z = np.load(PREP_CACHE, allow_pickle=True)
        meta = json.loads(str(z["meta"]))
        recs = json.loads(str(z["feature_meta"]))
        for i, r in enumerate(recs):
            r["X"] = z["X_%d" % i]
        if verbose:
            print("prep cache HIT  %s" % PREP_CACHE)
        return z["offset"], z["counts_real"], recs, meta

    if verbose:
        print("prep cache MISS -- building design + fitting ETAS (this is the "
              "expensive step; it is cached afterwards)")
    ctx = design_mod.build_design(
        data_dir=CFG["data_dir"], dlat=CFG["dlat"], dlon=CFG["dlon"],
        explore_frac=CFG["explore_frac"], verbose=verbose)
    explore, _hold = splits.temporal_split(ctx.n_days, CFG["explore_frac"])
    y = ctx.day_counts(CFG["mag_target"])
    base = bl.EtasV1(verbose=verbose, mag_target=CFG["mag_target"],
                     cache_path=ETAS_CACHE_1DEG)
    base.fit(ctx, y, explore)
    burn = int(base.burn_in_days)
    window = slice(burn, explore.stop)
    counts, offset = M.build_target(ctx, base, y, window)

    t0 = _dt.datetime.fromisoformat(str(ctx.meta["t0"]))
    feats = list(M.ephemeris_features(t0, ctx.n_days))
    feats += _family3_from_frozen(t0, ctx.n_days)
    marks = M.load_event_marks(ctx, CFG["data_dir"], ctx.meta["mag_floor"])
    feats += M.catalog_features(ctx, marks, (0,))

    recs = []
    for f in feats:
        recs.append({
            "name": f.name, "family": int(f.family), "kind": f.kind,
            "df": int(f.df), "periodic": bool(f.periodic),
            "block_days": float(f.block_days),
            "X": np.asarray(f.design(window, 0), dtype=np.float64),
        })

    meta = {
        "n_days_window": int(counts.size),
        "window": [int(window.start), int(window.stop)],
        "n_events_observed": float(counts.sum()),
        "etas_expected": float(offset.sum()),
        "etas_params": base.params,
        "etas_cache": ETAS_CACHE_1DEG.replace("\\", "/"),
        "design_n_cells": int(ctx.n_cells),
        "design_n_days": int(ctx.n_days),
        "t0": str(ctx.meta["t0"]),
        "config": CFG,
    }

    payload = {"offset": offset, "counts_real": counts,
               "meta": json.dumps(meta),
               "feature_meta": json.dumps([{k: v for k, v in r.items() if k != "X"}
                                           for r in recs])}
    for i, r in enumerate(recs):
        payload["X_%d" % i] = r["X"]
    os.makedirs(os.path.dirname(PREP_CACHE), exist_ok=True)
    np.savez_compressed(PREP_CACHE, **payload)
    if verbose:
        print("prep cache WROTE %s" % PREP_CACHE)
    return offset, counts, recs, meta


# ===================================================================== worker ==
def _init(offset, designs, dfs, periodics, n_boot, seed_base):
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    _W["offset"] = offset
    _W["designs"] = designs
    _W["dfs"] = dfs
    _W["periodic"] = periodics
    _W["n_boot"] = int(n_boot)
    _W["seed_base"] = int(seed_base)


def _simulate(cat):
    """The VIF = 1 target: independent Poisson draws from the fitted lambda."""
    rng = np.random.default_rng([_W["seed_base"], 1, int(cat)])
    return rng.poisson(_W["offset"]).astype(np.float64)


def _one(task):
    """(arm, catalog, feature index, block length) -> one VIF measurement."""
    arm, cat, fi, block = task
    off = _W["offset"]
    X = _W["designs"][fi]
    df = _W["dfs"][fi]
    periodic = _W["periodic"][fi]
    y = _simulate(cat)
    rng = np.random.default_rng([_W["seed_base"], 2, int(cat), int(fi),
                                 int(round(block * 1000))])
    S = M.score_stat_all_shifts(X, y, off)
    chi2_obs = float(S[0])
    Sb = M.score_stat_block_bootstrap(X, y, off, _W["n_boot"], rng,
                                      mean_block=float(block))
    p_boot = M.bootstrap_p(chi2_obs, Sb)
    if periodic:
        p_shift = None
        p = p_boot
    else:
        p_shift, _n_used = M.empirical_p(S, 10 ** 9, rng)
        p = max(p_shift, p_boot)
    vif, status = vif_of(chi2_obs, p, df)
    return {"arm": arm, "catalog": int(cat), "feature_idx": int(fi),
            "block_days": float(block), "df": int(df), "periodic": bool(periodic),
            "chi2_obs": chi2_obs, "p_boot": float(p_boot),
            "p_shift": (None if p_shift is None else float(p_shift)),
            "p_used": float(p), "vif": vif, "status": status,
            "p_at_boot_floor": bool(p_boot <= 1.0 / (_W["n_boot"] + 1.0) + 1e-12)}


# =================================================================== analysis ==
def _summ(vals):
    a = np.asarray([v for v in vals if v is not None and np.isfinite(v)], float)
    if a.size == 0:
        return None
    return {"n": int(a.size), "median": float(np.median(a)),
            "mean": float(np.mean(a)), "q25": float(np.percentile(a, 25)),
            "q75": float(np.percentile(a, 75)), "min": float(a.min()),
            "max": float(a.max())}


def _loglog(xs, ys):
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    ok = np.isfinite(xs) & np.isfinite(ys) & (xs > 0) & (ys > 0)
    if ok.sum() < 4:
        return {"n": int(ok.sum()), "note": "too few points"}
    lr = stats.linregress(np.log10(xs[ok]), np.log10(ys[ok]))
    rho, prho = stats.spearmanr(xs[ok], ys[ok])
    return {"n": int(ok.sum()), "loglog_slope": float(lr.slope),
            "loglog_slope_stderr": float(lr.stderr), "loglog_p": float(lr.pvalue),
            "spearman_rho": float(rho), "spearman_p": float(prho)}


def load_real():
    """The real-data curve: per-feature median VIF and block length, primary session."""
    with open(REAL_RESULTS, "r", encoding="utf-8") as fh:
        d = json.load(fh)
    prim = next(s for s in d["sessions"] if s["session"] == PRIMARY_SESSION)
    per_feat = {f["feature"]: f for f in prim["per_feature"]}
    return d, prim, per_feat


# ======================================================================= main ==
def main(argv=None):
    ap = argparse.ArgumentParser("f4_58_vif_control")
    ap.add_argument("--catalogs", type=int, default=24,
                    help=">= 20 required by §P7-8(b) item (ii)")
    ap.add_argument("--boots", type=int, default=2000)
    ap.add_argument("--jobs", type=int, default=0)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args(argv)
    if args.smoke:
        args.catalogs, args.boots = 3, 300

    t_start = time.time()
    offset, counts_real, recs, meta = prepare(verbose=True)
    offset = np.asarray(offset, dtype=np.float64)
    designs = [np.asarray(r["X"], dtype=np.float64) for r in recs]
    dfs = [int(r["df"]) for r in recs]
    periodics = [bool(r["periodic"]) for r in recs]
    names = [r["name"] for r in recs]
    blocks = [float(r["block_days"]) for r in recs]

    print("window days = %d, observed events = %.0f, ETAS expected = %.1f"
          % (offset.size, counts_real.sum(), offset.sum()))
    print("features = %d  (block_days span %.1f .. %.1f d)"
          % (len(recs), min(blocks), max(blocks)))

    # ---- task list ----------------------------------------------------------
    # arm "natural": every real feature at ITS OWN block_days -- matches the real
    #                measurement feature-for-feature.
    # arm "ladder" : the df=2 phase features at every block length on the ladder --
    #                isolates block length from feature identity, which is the
    #                variable §P7-8(b) says the two readings disagree about.
    tasks = []
    for cat in range(args.catalogs):
        for fi in range(len(recs)):
            tasks.append(("natural", cat, fi, blocks[fi]))
    phase_idx = [i for i, r in enumerate(recs) if r["df"] == 2 and r["periodic"]]
    for cat in range(args.catalogs):
        for fi in phase_idx:
            for b in BLOCK_LADDER:
                tasks.append(("ladder", cat, fi, b))
    print("tasks = %d  (natural %d, ladder %d over %d phase features)"
          % (len(tasks),
             args.catalogs * len(recs),
             args.catalogs * len(phase_idx) * len(BLOCK_LADDER), len(phase_idx)))

    jobs = args.jobs or max(1, (os.cpu_count() or 1) - 2)
    rows = []
    with ProcessPoolExecutor(
            max_workers=jobs,
            initializer=_init,
            initargs=(offset, designs, dfs, periodics, args.boots, CFG["seed"]),
    ) as ex:
        for i, row in enumerate(ex.map(_one, tasks, chunksize=4)):
            rows.append(row)
            if (i + 1) % 200 == 0:
                print("  %d/%d  (%.0fs)" % (i + 1, len(tasks), time.time() - t_start),
                      flush=True)
    print("sweep done in %.0fs" % (time.time() - t_start))

    for r in rows:
        r["feature"] = names[r["feature_idx"]]

    ok = [r for r in rows if r["status"] == "OK"]
    nat = [r for r in ok if r["arm"] == "natural"]
    lad = [r for r in ok if r["arm"] == "ladder"]

    # ---- control curve: VIF vs block length ---------------------------------
    ladder_curve = {}
    for b in BLOCK_LADDER:
        s = _summ([r["vif"] for r in lad if abs(r["block_days"] - b) < 1e-6])
        if s:
            ladder_curve["%.3f" % b] = s

    natural_by_feature = {}
    for nm in sorted(set(names)):
        s = _summ([r["vif"] for r in nat if r["feature"] == nm])
        if s:
            natural_by_feature[nm] = {
                "block_days": blocks[names.index(nm)],
                "df": dfs[names.index(nm)],
                "periodic": periodics[names.index(nm)],
                "vif_control": s,
            }

    # ---- real curve and the ratio -------------------------------------------
    real_doc, real_prim, real_feat = load_real()
    comparison = []
    for nm, rec in natural_by_feature.items():
        rf = real_feat.get(nm)
        if not rf or rf.get("vif_median") is None:
            continue
        c = rec["vif_control"]["median"]
        rv = float(rf["vif_median"])
        comparison.append({
            "feature": nm,
            "df": rec["df"],
            "block_days_control": rec["block_days"],
            "block_days_real": rf["block_days"],
            "vif_control_median": c,
            "vif_real_median": rv,
            "ratio_real_over_control": (rv / c) if c > 0 else None,
            "share_of_real_explained_by_null": (c / rv) if rv > 0 else None,
        })
    comparison.sort(key=lambda r: r["block_days_control"])

    ctrl_slope = _loglog([r["block_days_control"] for r in comparison],
                         [r["vif_control_median"] for r in comparison])
    ladder_slope = _loglog([float(k) for k in ladder_curve],
                           [v["median"] for v in ladder_curve.values()])
    real_slope = real_prim["block_length_dependence"]

    # ---- headline: the df=2 phase objects, which are what the floor is built on
    df2_ctrl_nat = _summ([r["vif"] for r in nat if r["df"] == 2])
    real_df2 = real_prim["overall_vif_df2_phase"]["median"]
    long_blocks = [b for b in BLOCK_LADDER if b >= 700.0]
    ctrl_long = _summ([r["vif"] for r in lad
                       if any(abs(r["block_days"] - b) < 1e-6 for b in long_blocks)])
    ctrl_short = _summ([r["vif"] for r in lad if r["block_days"] <= 60.0 + 1e-6])

    share_long = (ctrl_long["median"] / real_df2) if ctrl_long else None

    if share_long is None:
        verdict = "UNRESOLVED (no long-block control measurements)"
        reading = "UNRESOLVED"
    elif share_long >= 0.50:
        reading = "B"
        verdict = ("READING B OPERATING. On data whose TRUE VIF is 1 by "
                   "construction the same estimator returns a median VIF of "
                   "%.2f at block lengths >= 700 d, i.e. %.0f%% of the 24.08 "
                   "measured on the Earth. The inflation is substantially "
                   "manufactured by the ~10-block bootstrap null, the miner is "
                   "UNDER-POWERED rather than over-confident, and the flat "
                   "VIF = 24.08 floor is to that extent an artifact of the null."
                   % (ctrl_long["median"], 100.0 * share_long))
    elif share_long >= 0.20:
        reading = "A+B (mixed)"
        verdict = ("MIXED. The VIF = 1 control returns median %.2f at block "
                   "lengths >= 700 d, %.0f%% of the 24.08 measured on the Earth. "
                   "A material part of the inflation is the null's own width and "
                   "a material part is not. The VIF term must be re-specified as "
                   "BLOCK-LENGTH-CONDITIONAL against this simulated baseline "
                   "rather than adopted as a flat 24."
                   % (ctrl_long["median"], 100.0 * share_long))
    else:
        reading = "A"
        verdict = ("READING A STANDS. On true-VIF=1 data the estimator returns "
                   "median %.2f at block lengths >= 700 d, only %.0f%% of the "
                   "24.08 measured on the Earth. The estimator is not "
                   "manufacturing the inflation; the daily global M>=4.5 count "
                   "series is overdispersed relative to the fitted ETAS lambda, "
                   "the floor is real, and the programme is less powerful than "
                   "it believed."
                   % (ctrl_long["median"], 100.0 * share_long))

    out = {
        "id": "F4-58 item (ii) -- VIF ESTIMATOR CONTROL",
        "title": ("Measured VIF on catalogs simulated from the fitted ETAS lambda, "
                  "where the true VIF is 1 by construction"),
        "ruling": ("HYPOTHESIS_LEDGER.md §P7-8(b) ITEM (ii), SPECIFIED. Required "
                   "before the S-15 floor may retire any catalog entry on floor "
                   "grounds, before it may be quoted in any external artifact, and "
                   "before Tranche B is priced."),
        "priced_tests": 0,
        "priced_tests_note": (
            "Zero. No EXPLORE_COUNT.jsonl line is written. §P7-1(c) pricing-at-zero "
            "logic applies verbatim: this makes no rejection about the Earth and "
            "enters no BH vector. Every number here is computed on SIMULATED data "
            "with a known generating process; the object under test is the VIF "
            "ESTIMATOR, not a hypothesis about the world, so there is no "
            "multiplicity that could inflate any claim."),
        "generated_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "construction": {
            "true_vif": 1.0,
            "why": ("counts[t] ~ Poisson(lambda[t]) independently across days, "
                    "lambda = the engine's own fitted ETAS intensity summed over "
                    "the domain. Conditional on lambda the variance equals the "
                    "mean exactly, which is what the chi2 reference distribution "
                    "in the VIF identity assumes."),
            "etas_simulation_path": (
                "DIRECT POISSON DRAW from the fitted lambda series. engine/baseline.py "
                "exposes no simulate/rvs entry point (EtasV1: fit / rate / mu / report "
                "only), so no engine sampler was bypassed. A self-exciting re-simulation "
                "of the triggering cascade was deliberately NOT used: it would "
                "reintroduce clustering and would no longer be a true-VIF=1 control."),
            "estimator": "f4_58_vif.vif_of, imported (not reimplemented)",
            "null": ("engine.mine.score_stat_block_bootstrap at the same block "
                     "lengths, plus engine.mine.empirical_p for the circular-shift "
                     "arm; combined by the miner's own rule (p_boot for periodic "
                     "features, max(p_shift, p_boot) otherwise)"),
            "statistic": "engine.mine.score_stat_all_shifts",
            "n_catalogs": args.catalogs,
            "n_bootstrap_draws_per_test": args.boots,
            "seed_base": CFG["seed"],
        },
        "design": meta,
        "block_ladder": BLOCK_LADDER,
        "control_curve_ladder_df2_phase": ladder_curve,
        "control_by_feature_at_natural_block": natural_by_feature,
        "real_vs_control": comparison,
        "block_length_dependence": {
            "control_ladder": ladder_slope,
            "control_natural_across_features": ctrl_slope,
            "real_primary_session": real_slope,
        },
        "headline": {
            "real_vif_df2_phase_median": real_df2,
            "control_vif_df2_phase_median_at_natural_blocks": df2_ctrl_nat,
            "control_vif_at_long_blocks_ge_700d": ctrl_long,
            "control_vif_at_short_blocks_le_60d": ctrl_short,
            "share_of_real_vif_reproduced_at_long_blocks": share_long,
            "reading": reading,
        },
        "verdict": verdict,
        "caveats": [
            "The control uses lag 0 only. The real per-feature medians for "
            "families 3 and 4 are taken over 31 lags; lag changes the design's "
            "phase, not its block length or its dispersion, so the comparison is "
            "made at matched block_days and the lag reduction is declared here "
            "rather than corrected for.",
            "n_bootstrap draws per test is %d here against 50,000 in the real "
            "session. That sets the RESOLUTION of p_boot (and so of a single VIF "
            "value), not its expectation; the reported quantity is a median over "
            "%d catalogs x features, and tests whose p_boot sat at the Monte Carlo "
            "floor are flagged per row (`p_at_boot_floor`)."
            % (args.boots, args.catalogs),
            "The simulated target has the same MEAN function as the Earth's but "
            "not its clustering. That is the point of the control, and it is also "
            "its limit: it bounds how much VIF the ESTIMATOR manufactures, it does "
            "not estimate how much the Earth carries.",
            "engine/out/cache is only ADDED to by this script (a prep npz and a "
            "dedicated 1-degree ETAS parameter file). Nothing there is deleted or "
            "overwritten.",
        ],
        "rows": rows,
    }

    with open(OUT_PATH, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print("\nwrote %s" % OUT_PATH)

    # ---------------------------------------------------------------- printout
    print("\n=== CONTROL CURVE (true VIF = 1), df=2 phase features, forced ladder ===")
    print("%10s %8s %10s %10s %10s" % ("block_d", "n", "VIF_med", "q25", "q75"))
    for k, v in ladder_curve.items():
        print("%10.1f %8d %10.3f %10.3f %10.3f"
              % (float(k), v["n"], v["median"], v["q25"], v["q75"]))
    print("\ncontrol ladder block-length dependence:", ladder_slope)
    print("real  primary block-length dependence:", real_slope)

    print("\n=== REAL vs CONTROL at each feature's own block length ===")
    print("%-26s %8s %10s %10s %10s" %
          ("feature", "block_d", "VIF_ctrl", "VIF_real", "real/ctrl"))
    for c in comparison:
        print("%-26s %8.1f %10.3f %10.3f %10.2f"
              % (c["feature"], c["block_days_control"], c["vif_control_median"],
                 c["vif_real_median"],
                 -1 if c["ratio_real_over_control"] is None
                 else c["ratio_real_over_control"]))

    print("\n=== VERDICT ===")
    print("reading:", reading)
    print(verdict)
    return 0


if __name__ == "__main__":
    sys.exit(main())
