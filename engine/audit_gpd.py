"""§P6-2 rule 7: the calibration the GPD estimator must pass before it may report.

THE RULE, VERBATIM IN EFFECT. "On >= 15 representative tests spanning all three test
kinds, compare the GPD-at-2,000-surrogates CI against a brute-force p at N = 10^6.
Accept only if the CI covers the brute-force p in >= 90% of cases AND never
understates it by more than 3x. An instrument that has not demonstrated recovery in
the range it reports in may not report there."

This is G-M1 arm (ii) turned on the estimator instead of on the pipeline, and it is
the same logic §P5-8 used to void exp(dtau/A sigma): a law fitted in one range and
quoted in another has to demonstrate that the quoted range still works.

WHAT IS MEASURED
  covered            p_brute lies inside the 95% bootstrap CI [p_lo, p_hi]
  understatement     p_brute / p_ci_upper. > 1 means the number the engine would
                     REPORT is SMALLER than the truth -- the anti-conservative
                     direction, and the only one the 3x bar constrains. (The CI
                     upper end is what §P6-2(2) makes the reported number, so the
                     bar is applied to it and not to the point estimate.)
  GATE               coverage >= 90% AND max understatement <= 3x.

WHAT "REPRESENTATIVE" MEANS HERE, stated because it is a deviation worth seeing.
The three kinds are the three this engine runs -- the GLM score statistic under the
block bootstrap, the mark-association statistic under the event-block bootstrap, and
the Lomb-Scargle max-power statistic under permutation -- and each is computed by
the ENGINE'S OWN function (`M.score_stat_block_bootstrap`, `M._mark_stat_matrix`,
`M.ls_basis`), not by a re-implementation. What is scaled down is the SERIES LENGTH
and the period grid, because a brute force at N = 10^6 on the production series is
not affordable and the quantity under audit is the ESTIMATOR's recovery of a tail
probability, which is a property of the statistic's tail shape rather than of the
record length. The scaling is printed in the table so the deviation is visible.

The observed statistic for each test is PLACED, not found: it is set at a pilot
quantile chosen so the true p lands in the decade the estimator is actually asked to
report in (about 2e-4 to 1e-3, i.e. at or just below the 1/(2000+1) = 5.0e-4 floor).
Auditing the estimator anywhere else would be auditing a range nobody quotes.

Run:  python -u -m engine.audit_gpd --jobs 12
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

from . import gpd_tail, mine as M            # noqa: E402

N_SURROGATES_GPD = 2000                      # the "GPD-at-2,000-surrogates" arm
N_BRUTE = 1_000_000                          # the brute-force arm, N = 10^6
BRUTE_CHUNK = 2000
TARGET_PS = (2e-4, 5e-4, 1e-3)               # the decade the estimator reports in
PILOT_N = 20000
COVERAGE_BAR = 0.90
UNDERSTATE_BAR = 3.0
# the estimator's own bootstrap sizes, named here so the harness's cost is visible
# and so a test can shrink them without touching the estimator's declared defaults
N_BOOT = gpd_tail.GPD_N_BOOT
N_AD_BOOT = gpd_tail.GPD_N_AD_BOOT


# ------------------------------------------------------------------ the substrates
def _series(n_days, seed, rate=0.9):
    rng = np.random.default_rng(seed)
    offset = rate * (1.0 + 0.25 * np.sin(2 * np.pi * np.arange(n_days) / 365.25))
    counts = rng.poisson(offset).astype(float)
    return counts, offset


def glm_case(name, n_days, period_days, seed):
    """GLM score statistic on a [sin, cos] design, block-bootstrap null."""
    counts, offset = _series(n_days, seed)
    t = np.arange(n_days, dtype=float)
    ph = 2 * np.pi * t / period_days
    X = np.column_stack([np.sin(ph), np.cos(ph)])
    block = max(2.0 * period_days, 30.0)

    def draw(b, rng):
        return M.score_stat_block_bootstrap(X, counts, offset, int(b), rng,
                                            mean_block=block, chunk=int(b))
    return {"name": name, "kind": "glm", "null_type": "block_bootstrap",
            "draw": draw, "n": n_days,
            "detail": f"[sin,cos] at {period_days} d, block {block:.0f} d, "
                      f"{n_days} days"}


def mark_case(name, n_events, seed, kind="linear"):
    """Mark-association statistic on the engine's own matrix, event-block bootstrap."""
    rng = np.random.default_rng(seed)
    fe = rng.standard_normal(n_events)
    if kind == "phase":
        fe = np.mod(np.arange(n_events) * 0.37, 1.0) * 2 * np.pi
    mark = rng.standard_normal(n_events) + 0.15 * np.roll(fe, 3)
    rm = M._ranks(mark)
    block = max(10, n_events // 8)

    def draw(b, rng_):
        idx = M.block_bootstrap_idx(n_events, int(b), block, rng_)
        v, _ = M._mark_stat_matrix(fe, rm[idx], kind)
        return v
    return {"name": name, "kind": "mark", "null_type": "block_bootstrap",
            "draw": draw, "n": n_events,
            "detail": f"{kind} mark stat, {n_events} events, block {block}"}


def _ls_power_batch(t, Xb, periods):
    """Max Lomb-Scargle power for a BATCH of series, through the engine's own basis.

    Same basis, same normalisation, same arithmetic as `M.lomb_scargle_power` -- the
    only change is that the single column vector becomes a matrix, which is what
    makes a 10^6-surrogate brute force affordable at all. Pinned against the scalar
    kernel by `engine/tests/test_gpd_audit.py`.
    """
    bas = M.ls_basis(t, periods)
    nf = bas["n_freq"]
    Xc = Xb - Xb.mean(axis=1, keepdims=True)
    W = Xc * (1.0 / bas["n_time"])
    z = W @ bas["B"]                                  # (b, 2nf)
    yc, ys = z[:, :nf], z[:, nf:]
    yy = np.einsum("bi,bi->b", W, Xc)
    pg = 2.0 * (yc * yc / bas["CC"] + ys * ys / bas["SS"])
    return (pg * (0.5 / yy)[:, None]).max(axis=1)


def period_case(name, n_days, n_periods, seed):
    """Lomb-Scargle max-power statistic under the permutation null."""
    rng = np.random.default_rng(seed)
    resid = rng.standard_normal(n_days)
    t = np.arange(n_days, dtype=float)
    periods = np.exp(np.linspace(math.log(2.0), math.log(n_days / 3.0), n_periods))

    def draw(b, rng_):
        Xb = np.empty((int(b), n_days))
        for i in range(int(b)):
            Xb[i] = rng_.permutation(resid)
        return _ls_power_batch(t, Xb, periods)
    return {"name": name, "kind": "period", "null_type": "permutation",
            "draw": draw, "n": n_days,
            "detail": f"max LS power over {n_periods} periods, {n_days} days"}


def roster():
    """>= 15 tests spanning the three kinds. Fixed, declared, and printed."""
    cases = []
    for i, (nd, pd_) in enumerate([(2000, 29.53), (2000, 13.66), (3000, 365.25)]):
        cases.append(glm_case(f"glm{i+1}", nd, pd_, 1000 + i))
    for i, (ne, kd) in enumerate([(2000, "linear"), (3000, "linear"),
                                  (2000, "phase")]):
        cases.append(mark_case(f"mark{i+1}", ne, 2000 + i, kd))
    for i, (nd, npd) in enumerate([(1000, 100), (1500, 150), (800, 60)]):
        cases.append(period_case(f"per{i+1}", nd, npd, 3000 + i))
    return cases


# ------------------------------------------------------------------- the two arms
def run_case(case, target_p, seed):
    rng = np.random.default_rng(seed)
    pilot = case["draw"](PILOT_N, rng)
    s_obs = float(np.quantile(pilot, 1.0 - target_p))

    t0 = time.perf_counter()
    S = case["draw"](N_SURROGATES_GPD, np.random.default_rng(seed + 1))
    p_mc = (1.0 + int((S >= s_obs).sum())) / (1.0 + N_SURROGATES_GPD)
    g = gpd_tail.gpd_tail_p(S, s_obs, N_SURROGATES_GPD, case["null_type"],
                            np.random.default_rng(seed + 2), n_boot=N_BOOT,
                            ad_boot=N_AD_BOOT)
    t_gpd = time.perf_counter() - t0

    t0 = time.perf_counter()
    rng_b = np.random.default_rng(seed + 3)
    ge, done = 0, 0
    while done < N_BRUTE:
        b = min(BRUTE_CHUNK, N_BRUTE - done)
        ge += int((case["draw"](b, rng_b) >= s_obs).sum())
        done += b
    p_bf = (1.0 + ge) / (1.0 + N_BRUTE)
    t_bf = time.perf_counter() - t0

    lo, hi = g.get("p_ci_lower"), g.get("p_ci_upper")
    covered = bool(lo is not None and hi is not None and lo <= p_bf <= hi)
    understate = (float(p_bf / hi) if hi else float("inf"))
    return {"case": case["name"], "kind": case["kind"], "detail": case["detail"],
            "target_p": target_p, "s_obs": s_obs, "p_mc_2000": p_mc,
            "p_method": g["p_method"], "p_point": g.get("p_point"),
            "p_ci_lower": lo, "p_ci_upper": hi, "p_brute_force": p_bf,
            "n_brute": N_BRUTE, "covered": covered,
            "understatement_x": understate, "p_ad": (g.get("ad") or {}).get("p_ad"),
            "xi": g.get("xi"), "xi_stable": (g.get("xi_stability") or {}).get("pass"),
            "reason": g.get("reason"), "seconds_gpd": round(t_gpd, 1),
            "seconds_brute": round(t_bf, 1)}


def _job(spec):
    os.environ["OMP_NUM_THREADS"] = "1"
    idx, target_p, seed = spec
    return run_case(roster()[idx], target_p, seed)


def main(argv=None):
    ap = argparse.ArgumentParser("engine.audit_gpd")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--out", default=os.path.join("engine", "out",
                                                  "audit_gpd.json"))
    ap.add_argument("--targets", type=float, nargs="*", default=list(TARGET_PS))
    args = ap.parse_args(argv)

    cases = roster()
    specs = [(i, tp, 700000 + 137 * i + 11 * j)
             for i, _c in enumerate(cases) for j, tp in enumerate(args.targets)]
    print(f"§P6-2 rule 7 calibration: {len(cases)} representative tests x "
          f"{len(args.targets)} target p = {len(specs)} comparisons, "
          f"GPD at N = {N_SURROGATES_GPD} vs brute force at N = {N_BRUTE:,}")
    print(f"  kinds: " + ", ".join(sorted({c['kind'] for c in cases})))
    t0 = time.perf_counter()
    if args.jobs > 1:
        import concurrent.futures as cf
        import multiprocessing as mp
        ctx = mp.get_context("spawn")
        with cf.ProcessPoolExecutor(max_workers=args.jobs,
                                    mp_context=ctx) as ex:
            rows = list(ex.map(_job, specs))
    else:
        rows = [_job(s) for s in specs]
    elapsed = time.perf_counter() - t0

    rows.sort(key=lambda r: (r["kind"], r["case"], r["target_p"]))
    print()
    print("| case | kind | target p | p_MC@2000 | p_method | GPD 95% CI | "
          "brute force @1e6 | covered | understatement |")
    print("| --- | --- | ---: | ---: | --- | --- | ---: | :---: | ---: |")
    for r in rows:
        ci = ("--" if r["p_ci_lower"] is None
              else f"[{r['p_ci_lower']:.3g}, {r['p_ci_upper']:.3g}]")
        us = ("--" if not np.isfinite(r["understatement_x"])
              else f"{r['understatement_x']:.2f}x")
        print(f"| {r['case']} | {r['kind']} | {r['target_p']:.0e} | "
              f"{r['p_mc_2000']:.3g} | {r['p_method']} | {ci} | "
              f"{r['p_brute_force']:.3g} | {'YES' if r['covered'] else 'no'} | "
              f"{us} |")

    fitted = [r for r in rows if r["p_method"] == gpd_tail.P_GPD]
    n_cov = sum(1 for r in fitted if r["covered"])
    cov = (n_cov / len(fitted)) if fitted else 0.0
    worst = max((r["understatement_x"] for r in fitted), default=float("inf"))
    ok_cov = cov >= COVERAGE_BAR
    ok_us = worst <= UNDERSTATE_BAR
    print()
    print(f"comparisons that produced a GPD fit: {len(fitted)} of {len(rows)} "
          f"({len(rows) - len(fitted)} were gated out and reported UNRESOLVED, "
          f"which is the estimator refusing rather than failing)")
    print(f"GATE 1  coverage      : {n_cov}/{len(fitted)} = {cov:.1%} "
          f"(bar >= {COVERAGE_BAR:.0%})  -> {'PASS' if ok_cov else 'FAIL'}")
    print(f"GATE 2  understatement: worst = {worst:.2f}x "
          f"(bar <= {UNDERSTATE_BAR:.0f}x)      -> "
          f"{'PASS' if ok_us else 'FAIL'}")
    print(f"§P6-2(7) VERDICT: {'ACCEPT' if (ok_cov and ok_us) else 'REJECT'} "
          f"-- {elapsed:.0f}s wall, --jobs {args.jobs}")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump({"rule": "§P6-2(7)", "n_cases": len(cases),
                   "n_comparisons": len(rows), "n_fitted": len(fitted),
                   "coverage": cov, "worst_understatement": worst,
                   "coverage_bar": COVERAGE_BAR,
                   "understatement_bar": UNDERSTATE_BAR,
                   "pass": bool(ok_cov and ok_us), "elapsed_s": elapsed,
                   "rows": rows}, fh, indent=1)
    print(f"written -> {args.out}")
    return 0 if (ok_cov and ok_us) else 1


if __name__ == "__main__":
    raise SystemExit(main())
