"""Popper's binding acceptance test for the adaptive surrogate ladder.

The ladder (`--ladder`, Besag & Clifford 1991 sequential Monte Carlo) may not
become a default, and may not be trusted in a report, unless it passes BOTH of the
checks below. They are separate claims and they fail for different reasons:

  1. UNIFORMITY. Under a true null the ladder's p-values must be U(0,1). This is
     what "exactly valid under optional stopping" means, and it is the whole
     justification for stopping early at all.

  2. AGREEMENT. For at least 99% of tests, the ladder's p must agree with the
     p that the SAME draws would have produced had they run to the full N_max,
     within 3 * sqrt(p(1-p)/n) -- three standard errors of a binomial proportion
     at n = N_max. This catches a ladder that is valid but useless. The two
     p-values are computed from the IDENTICAL surrogate stream, so any
     disagreement is the stopping rule's doing and nothing else's.


TWO FLAGGED DEVIATIONS FROM THE CRITERIA AS LITERALLY WORDED
------------------------------------------------------------
Both criteria are measured and printed exactly as worded. Neither is quietly
relaxed. But as worded, both are unsatisfiable by ANY correct Besag-Clifford rule,
for reasons that are properties of sequential Monte Carlo and not of this
implementation, so the exit code is gated on stated repairs. The adjudicator, not
this script, should decide whether the repairs are accepted.

  (1) "KS distance from U(0,1) not significant at 0.01" cannot be met at n >= 5000
      because a sequential Monte Carlo p-value is DISCRETE: it takes the values
      h/l for integer l, plus (1+c)/(1+N_max). At N_max=1000, h=25 there are only
      ~800 attainable values in [0,1], so the ECDF is a staircase and a two-sided
      KS test against a CONTINUOUS uniform rejects on the step structure alone, at
      any sample size large enough to see it. The rejection is entirely in the
      CONSERVATIVE direction (measured D- = 0.041, D+ = 0.001 at n = 20000):
      discreteness makes the p-value stochastically LARGER than uniform, which is
      validity, not a defect.

      REPAIR, gated: the one-sided statistic D+ = sup_a [ECDF(a) - a], which tests
      the only direction that can hurt anyone -- P(p <= a) > a, i.e.
      anti-conservatism. The rule PASSES when D+ is below the one-sided KS critical
      value at 0.01. This is a strictly stronger requirement in the direction that
      matters and a strictly weaker one in the direction that does not.

  (2) "agreement within 3*sqrt(p(1-p)/n) for >= 99%" uses only the FULL estimator's
      standard error. The ladder's own estimator h/l has standard deviation about
      q*sqrt((1-q)/h) -- roughly 1/sqrt(h) = 20% RELATIVE at h = 25 -- which is not
      error, it is the deliberate price of stopping early and is the entire source
      of the saving. At q = 0.5, N_max = 1000 the tolerance is +/-0.047 while the
      ladder's own SD is 0.071, so the criterion must fail in the bulk however
      correct the rule is. Measured: 57% at h = 25.

      REPAIR, gated: the same 3-sigma comparison against the COMBINED standard
      error sqrt(p(1-p)/N_max + p^2(1-p)/h) of the two estimators being compared.
      Measured: 99.6%, which clears the 99% bar as stated.

      Reported alongside, ungated, because it is the number a reader actually
      wants: agreement under the LITERAL tolerance restricted to the
      decision-relevant tail (p_full <= 0.1), where the ladder mostly runs to the
      cap and the two estimators mostly coincide exactly. Measured: 97%.

§P6-9, THE AMENDED CRITERIA (2026-08-12). Both repairs above were GRANTED IN FULL,
with three additions, and this module now implements the amended rule rather than
proposing it. What gates, and what does not:

  GATE 1  one-sided KS D+ = sup_a (F_n(a) - a) against D_crit = sqrt(-ln(a)/(2n))
          at a = 0.01. D- is REPORTED AND DOES NOT GATE: super-uniformity, not
          uniformity, is the property that licenses the ladder, and discreteness
          guarantees D- > 0 for any correct estimator. The two-sided statistic and
          its p-value are still printed, labelled as the superseded literal rule
          and as NOT GATING, so the deviation stays visible instead of becoming
          folklore.
  GATE 2  agreement within 3 x the COMBINED standard error of the two estimators
          being compared, sqrt(p^2/h + p(1-p)/N), at >= 99.0% of tests. Note the
          ladder term is p^2/h exactly as §P6-9 words it (SE_ladder = p/sqrt(h)),
          not the p^2(1-p)/h this file used when it was still proposing the repair.
  GATE 3  (§P6-9 addition 2) signed-bias: |median[(p_ladder - p_full)/SE_combined]|
          <= 0.25. A two-sided count bar is blind to a small systematic downward
          shift of the ladder p's, which is exactly the anti-conservative failure
          this audit exists to catch; a centred cloud and a shifted cloud can have
          identical exceedance counts and only the signed statistic separates them.
  GATE 4  the cap-branch invariant: where the cap was reached with c < h the two
          estimators are the same formula on the same draws and must agree to the
          bit.

  REPORTED, NOT GATING (§P6-9 addition 1): the realised conservatism
  P(p <= a)/a at a in {0.1, 0.01, 0.001} -- the operating points where the ladder
  is actually spent -- alongside the ladder's own relative standard error
  1/sqrt(h). A sup-distance somewhere in the middle of the unit interval is not
  what costs the miner power; conservatism at the rejection thresholds is.

  §P6-9 addition 3: THE AUDIT IS RE-RUN AT EACH DECLARED PRODUCTION N_max, and a
  stratum may not run at an N_max the audit has not covered. The estimator's
  discreteness and its whole cap branch (1+c)/(1+N_max) are functions of N_max, so
  an audit at 1000 does not license production at 10^4. That is S-17 applied to the
  ladder itself. `--n-max` exists for exactly this, and `--jobs` exists because at
  production N_max the honest audit is 20,000 x N_max surrogate draws for the
  full-run comparison arm -- 200 million at N_max = 10^4 -- which is worth
  parallelising and is NOT worth shortcutting.

Run:
    python -u -m engine.audit_ladder --fast                    # ~5k tests, seconds
    python -u -m engine.audit_ladder                           # 20k tests, N_max=1e3
    python -u -m engine.audit_ladder --n-max 10000 --jobs 16   # production N_max

Run from replication/. Writes nothing; prints a verdict and exits non-zero on
failure so it can gate a build.
"""

from __future__ import annotations

import argparse
import concurrent.futures as _cf
import math
import multiprocessing as _mp
import os
import sys
import time

import numpy as np

from . import mine as M

# §P6-9 addition 1: the rejection thresholds where the ladder is actually spent.
OPERATING_POINTS = (0.1, 0.01, 0.001)
# §P6-9 addition 2: the signed-bias gate.
SIGNED_BIAS_MAX = 0.25

# The null the audit draws from. A test statistic whose null distribution is known
# in closed form, so "the p-value should be uniform" is a statement about the
# STOPPING RULE and not about some approximation elsewhere. Exponential rather than
# uniform on purpose: a heavy-ish right tail is where a sequential rule is most
# likely to misbehave, and the exceedance probability is still exactly computable.
def _null_draw(b, rng):
    return rng.exponential(1.0, size=int(b))


def _observed_from_null(rng):
    """An observed statistic drawn from the SAME null -> p must be uniform."""
    return float(rng.exponential(1.0))


def audit_range(start, stop, n_max, h, chunk, seed):
    """Run audit tests [start, stop). PURE, and addressable by test index.

    Every test's streams come from its own canonical key, keyed on the test INDEX,
    so this slices identically however the range is cut up: `--jobs 1` and
    `--jobs 30` audit the same 20,000 tests with the same draws in the same order.
    Same design as the period scan's per-surrogate seeding (§P6-5), for the same
    reason -- a parallel audit whose answer depends on the scheduler is not an audit.
    """
    n = int(stop) - int(start)
    p_ladder = np.empty(n)
    p_full = np.empty(n)
    n_drawn = np.empty(n, dtype=np.int64)
    l_stop_none = np.zeros(n, dtype=bool)

    for j, i in enumerate(range(int(start), int(stop))):
        # Each audit test gets its own key-derived stream, exactly as a real test
        # does -- so this audits the production seeding path, not a toy one.
        key = M.test_key(seed, f"audit_{i}", "ladder_audit", lag=0,
                         null_type="uniformity_audit")
        ss = M.seed_sequence_for(key)
        # The observed statistic gets its own key (a different null_type), never a
        # rung of the surrogate stream -- if it shared one, "observed" and
        # "surrogate" would be correlated and the audit would be measuring itself.
        obs = _observed_from_null(np.random.default_rng(M.seed_sequence_for(
            M.test_key(seed, f"audit_{i}", "ladder_audit", lag=0,
                       null_type="observed_statistic"))))

        def rung_rng(r, _ss=ss):
            return np.random.default_rng(M.rung_seed_sequence(_ss, r))

        res = M.besag_clifford_p(_null_draw, obs, n_max=n_max, rung_rng=rung_rng,
                                 h=h, chunk=chunk)
        p_ladder[j] = res["p"]
        n_drawn[j] = res["n_drawn"]
        l_stop_none[j] = res["l_stop"] is None

        # The same stream, run to the cap, scored the ordinary fixed-sample way.
        # Identical draws => any difference is the stopping rule's alone. This arm
        # is what costs N_max draws per test and is why --jobs exists; it is NOT
        # shortcut-able, because the whole agreement criterion is defined against
        # the counterfactual full run of THESE draws.
        c, total, r = 0, 0, 0
        while total < n_max:
            b = min(chunk, n_max - total)
            c += int((_null_draw(b, rung_rng(r)) >= obs).sum())
            total += b
            r += 1
        p_full[j] = (1.0 + c) / (1.0 + n_max)
    return p_ladder, p_full, n_drawn, l_stop_none


def _audit_range_star(a):
    return audit_range(*a)


def run_audit(n_tests=20000, n_max=1000, h=25, chunk=200, seed=20260812,
              ks_alpha=0.01, agree_frac=0.99, verbose=True, jobs=1):
    from scipy import stats

    t0 = time.time()
    jobs = max(1, int(jobs))
    if jobs == 1:
        p_ladder, p_full, n_drawn, l_stop_none = audit_range(
            0, n_tests, n_max, h, chunk, seed)
    else:
        # ~4 slices per worker so a straggler cannot hold the whole audit.
        n_slices = min(n_tests, jobs * 4)
        edges = np.linspace(0, n_tests, n_slices + 1).astype(int)
        args = [(int(edges[k]), int(edges[k + 1]), n_max, h, chunk, seed)
                for k in range(n_slices) if edges[k + 1] > edges[k]]
        saved = {}
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS"):
            saved[var] = os.environ.get(var)
            os.environ[var] = "1"
        try:
            ctx = _mp.get_context("spawn")
            with _cf.ProcessPoolExecutor(max_workers=jobs, mp_context=ctx) as ex:
                parts = list(ex.map(_audit_range_star, args))
        finally:
            for var, val in saved.items():
                if val is None:
                    os.environ.pop(var, None)
                else:
                    os.environ[var] = val
        # concatenated IN SLICE ORDER, so the arrays are indexed by test index
        p_ladder = np.concatenate([p[0] for p in parts])
        p_full = np.concatenate([p[1] for p in parts])
        n_drawn = np.concatenate([p[2] for p in parts])
        l_stop_none = np.concatenate([p[3] for p in parts])
    if p_ladder.size != n_tests:
        raise RuntimeError(f"audit assembled {p_ladder.size} tests, declared "
                           f"{n_tests}: a dropped slice would move every number "
                           f"below. Refusing to report.")

    # ---- criterion 1, as literally worded: two-sided KS vs continuous U(0,1) ----
    ks = stats.kstest(p_ladder, "uniform")
    ks_pass_literal = bool(ks.pvalue > ks_alpha)

    # ---- criterion 1, gated repair: one-sided D+, the anti-conservative side ----
    s = np.sort(p_ladder)
    idx = np.arange(1, n_tests + 1)
    d_plus = float(np.max(idx / n_tests - s))          # ECDF above the diagonal
    d_minus = float(np.max(s - (idx - 1) / n_tests))   # ECDF below (= discreteness)
    # one-sided KS critical value at alpha, from P(D+ > d) = exp(-2 n d^2)
    d_crit = float(math.sqrt(-0.5 * math.log(ks_alpha) / n_tests))
    valid_pass = bool(d_plus <= d_crit)

    # ---- §P6-9 addition 1: realised conservatism AT THE OPERATING POINTS -------
    # Reported, never gating. D+ and D- are sups over the whole unit interval; what
    # costs the miner power is the ratio P(p <= a)/a at the thresholds where a
    # rejection is actually declared. A ratio below 1 is conservatism (power given
    # up); above 1 is anti-conservatism, which is what GATE 1 exists to forbid.
    op_points = []
    for a in OPERATING_POINTS:
        frac = float(np.mean(p_ladder <= a))
        frac_full = float(np.mean(p_full <= a))
        op_points.append({
            "alpha": float(a), "p_le_alpha": frac, "ratio": frac / a,
            "p_le_alpha_full_run": frac_full, "ratio_full_run": frac_full / a,
            "n_at_or_below": int((p_ladder <= a).sum()),
        })
    ladder_rel_se = 1.0 / math.sqrt(h)          # the ladder's own relative SE

    # ---- criterion 2, as literally worded (SUPERSEDED, not gating) ------------
    q = np.clip(p_full * (1.0 - p_full), 0.0, None)
    resolution = 1.0 / (n_max + 1.0)
    tol_literal = np.maximum(3.0 * np.sqrt(q / n_max), resolution)
    agree_literal = float(np.mean(np.abs(p_ladder - p_full) <= tol_literal))

    # ---- GATE 2, §P6-1(6b) AMENDED: 3 x the COMBINED SE of both estimators ----
    # sqrt(p^2/h + p(1-p)/N), the root-sum-square of the two estimators' standard
    # errors, each under its own sampling scheme. The ladder term is p^2/h exactly
    # as §P6-9 words it (SE_ladder = p/sqrt(h)); this file used p^2(1-p)/h while it
    # was still PROPOSING the repair, which was very slightly the more permissive
    # of the two, and the granted wording is what ships.
    se_comb = np.sqrt(p_full ** 2 / h + q / n_max)
    tol_comb = np.maximum(3.0 * se_comb, resolution)
    agree_combined = float(np.mean(np.abs(p_ladder - p_full) <= tol_comb))

    # ---- GATE 3, §P6-9 addition 2: the SIGNED bias ----------------------------
    # 99.5% two-sided agreement is fully compatible with a small systematic
    # DOWNWARD shift of the ladder p's -- the anti-conservative failure this audit
    # exists to catch. A centred cloud and a shifted cloud can have identical
    # exceedance counts; only the signed statistic separates them.
    signed = (p_ladder - p_full) / np.maximum(se_comb, 1e-300)
    signed_median = float(np.median(signed))
    signed_mean = float(np.mean(signed))
    bias_pass = bool(abs(signed_median) <= SIGNED_BIAS_MAX)

    # ---- reported alongside: the literal tolerance in the tail that matters ----
    tail = p_full <= 0.10
    agree_tail = (float(np.mean((np.abs(p_ladder - p_full) <= tol_literal)[tail]))
                  if tail.any() else float("nan"))
    # Where the cap was reached with c < h the two estimators are the SAME formula
    # on the SAME draws, so they must agree to the bit. Anything else is a bug.
    capped = np.asarray(l_stop_none, dtype=bool)
    exact_on_cap = (float(np.mean(p_ladder[capped] == p_full[capped]))
                    if capped.any() else float("nan"))

    saving = 1.0 - float(n_drawn.mean()) / n_max
    agree_pass = bool(agree_combined >= agree_frac)
    exact_pass = bool(not capped.any() or exact_on_cap == 1.0)
    out = {
        "n_tests": int(n_tests), "n_max": int(n_max), "h": int(h),
        "chunk": int(chunk), "seed": int(seed),
        "ks_stat": float(ks.statistic), "ks_p": float(ks.pvalue),
        "ks_alpha": ks_alpha, "ks_pass_literal": ks_pass_literal,
        "d_plus": d_plus, "d_minus": d_minus, "d_crit_one_sided": d_crit,
        "validity_pass": valid_pass,
        "n_distinct_p": int(np.unique(p_ladder).size),
        "agree_frac_literal": agree_literal,
        "agree_frac_combined": agree_combined,
        "agree_frac_tail_literal": agree_tail,
        "agree_required": agree_frac, "agree_pass": agree_pass,
        "frac_capped": float(capped.mean()), "exact_on_cap": exact_on_cap,
        "exact_on_cap_pass": exact_pass,
        "mean_draws": float(n_drawn.mean()), "draw_saving": saving,
        # §P6-9 addition 1 (reported, not gating)
        "operating_points": op_points, "ladder_relative_se": ladder_rel_se,
        # §P6-9 addition 2 (GATE 3)
        "signed_bias_median": signed_median, "signed_bias_mean": signed_mean,
        "signed_bias_max": SIGNED_BIAS_MAX, "signed_bias_pass": bias_pass,
        "jobs": int(jobs),
        "pass": bool(valid_pass and agree_pass and bias_pass and exact_pass),
        "pass_literal_as_worded": bool(ks_pass_literal
                                       and agree_literal >= agree_frac),
        "elapsed_s": round(time.time() - t0, 1),
    }
    if verbose:
        def _v(ok):
            return "PASS" if ok else "FAIL"
        print("=" * 78)
        print(f"LADDER UNIFORMITY AUDIT -- criteria per HYPOTHESIS_LEDGER.md "
              f"§P6-9 (amended)")
        print("=" * 78)
        print(f"  {n_tests} null tests, N_max={n_max}, h={h}, chunk={chunk}, "
              f"seed={seed}, jobs={jobs} ({out['elapsed_s']} s)")
        print(f"  full-run comparison arm drew {n_tests * n_max / 1e6:.2f} M "
              f"surrogates; the ladder arm drew "
              f"{n_drawn.sum() / 1e6:.2f} M.")
        print()
        print("  GATE 1 -- validity (§P6-1(6a) AMENDED: one-sided, D+ only)")
        print(f"    D+ (anti-conservative side) = {d_plus:.6f} vs D_crit = "
              f"sqrt(-ln(a)/(2n)) = {d_crit:.6f} at a={ks_alpha} -> "
              f"{_v(valid_pass)}")
        print(f"    D- (conservative side)      = {d_minus:.6f}   REPORTED, "
              f"DOES NOT GATE -- super-uniformity is the licensed property and "
              f"discreteness guarantees D- > 0")
        print(f"    [SUPERSEDED, NOT GATING] two-sided KS vs continuous U(0,1): "
              f"D = {out['ks_stat']:.6f}, p = {out['ks_p']:.4g} "
              f"({'would REJECT' if not ks_pass_literal else 'would accept'} "
              f"at {ks_alpha})")
        print(f"        p takes only {out['n_distinct_p']} distinct values over "
              f"{n_tests} tests: a sequential MC p is discrete, so the two-sided "
              f"test rejects on the staircase alone. The rejection is "
              f"{'CONSERVATIVE' if d_minus > d_plus else 'ANTI-CONSERVATIVE'}.")
        print()
        print("  GATE 2 -- agreement with the SAME draws run to N_max")
        print(f"    within 3 x combined SE sqrt(p^2/h + p(1-p)/N): "
              f"{agree_combined * 100:.2f}% -> {_v(agree_pass)} "
              f"(needs >= {agree_frac * 100:.1f}%)")
        print(f"    [SUPERSEDED, NOT GATING] within 3*sqrt(p(1-p)/N_max) only: "
              f"{agree_literal * 100:.2f}% -- counts only the full estimator's "
              f"error, so it cannot be met in the bulk however correct the rule is")
        print(f"    [reported] superseded tolerance in the decision-relevant tail "
              f"p_full <= 0.1: {agree_tail * 100:.2f}%")
        print()
        print("  GATE 3 -- signed bias (§P6-9 addition 2)")
        print(f"    median[(p_ladder - p_full)/SE_combined] = {signed_median:+.4f} "
              f"vs |median| <= {SIGNED_BIAS_MAX} -> {_v(bias_pass)}")
        print(f"    (mean {signed_mean:+.4f}; a NEGATIVE median is the "
              f"anti-conservative direction -- ladder p's systematically below the "
              f"full-run p's)")
        print()
        print("  GATE 4 -- cap-branch invariant")
        print(f"    cap reached with c < h in {out['frac_capped'] * 100:.1f}% of "
              f"tests; on those the two estimators are the same formula on the "
              f"same draws and agree bit-for-bit in {exact_on_cap * 100:.2f}% -> "
              f"{_v(exact_pass)}")
        print()
        print("  REPORTED, NOT GATING -- realised conservatism at the operating "
              "points (§P6-9 addition 1)")
        print(f"    {'alpha':>8} {'P(p<=a)':>10} {'ratio':>8} {'n':>7}   "
              f"{'full-run P':>11} {'ratio':>8}")
        for op in op_points:
            print(f"    {op['alpha']:>8.3g} {op['p_le_alpha']:>10.5f} "
                  f"{op['ratio']:>8.3f} {op['n_at_or_below']:>7d}   "
                  f"{op['p_le_alpha_full_run']:>11.5f} "
                  f"{op['ratio_full_run']:>8.3f}")
        print(f"    ladder's own relative standard error 1/sqrt(h) = "
              f"{ladder_rel_se:.4f} ({ladder_rel_se * 100:.1f}%) at h={h}")
        print(f"    ratio < 1 is conservatism (power given up); ratio > 1 is the "
              f"anti-conservatism GATE 1 forbids. These are what the ladder COST, "
              f"which D- alone does not tell you.")
        print()
        print(f"  saving: mean {out['mean_draws']:.1f} draws of {n_max} "
              f"({saving * 100:.1f}% of surrogate work saved)")
        print(f"  VERDICT (§P6-9 amended gates 1-4): {_v(out['pass'])}")
        print(f"  VERDICT (original §P6-1(6) exactly as first worded, SUPERSEDED): "
              f"{_v(out['pass_literal_as_worded'])}")
        print(f"  This run licenses the ladder AT N_max = {n_max} ONLY (§P6-9 "
              f"addition 3 / S-17): a stratum may not run at an N_max the audit "
              f"has not covered.")
        print("=" * 78)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fast", action="store_true",
                    help="5000 tests at N_max=500 (the suite-sized run)")
    ap.add_argument("--n-tests", type=int, default=None)
    ap.add_argument("--n-max", type=int, default=None)
    ap.add_argument("--h", type=int, default=M.LADDER_DEFAULTS["h"])
    ap.add_argument("--chunk", type=int, default=M.LADDER_DEFAULTS["chunk"])
    ap.add_argument("--seed", type=int, default=20260812)
    ap.add_argument("--jobs", type=int, default=1,
                    help="worker processes; 0 = max(1, cpu_count() - 2). The "
                         "audit is sliced by TEST INDEX and every stream is "
                         "key-derived, so the result is identical at any --jobs")
    ap.add_argument("--json", default=None,
                    help="also write the full result dict here")
    a = ap.parse_args(argv)
    n_tests = a.n_tests if a.n_tests else (5000 if a.fast else 20000)
    n_max = a.n_max if a.n_max else (500 if a.fast else 1000)
    jobs = a.jobs if a.jobs else max(1, (os.cpu_count() or 2) - 2)
    out = run_audit(n_tests=n_tests, n_max=n_max, h=a.h, chunk=a.chunk,
                    seed=a.seed, jobs=jobs)
    if a.json:
        import json
        with open(a.json, "w", encoding="utf-8") as fh:
            json.dump(out, fh, indent=1)
        print(f"written -> {a.json}")
    return 0 if out["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
