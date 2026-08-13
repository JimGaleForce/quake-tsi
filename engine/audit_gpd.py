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

ATTEMPT 2, AND WHY THE DESIGN CHANGED WHILE THE BARS DID NOT (§P7-7). Attempt 1 ran
the percentile interval on 27 comparisons and scored 22/26 = 84.6% coverage. That is
NOT a rejection: the exact 95% interval on 22/26 is [0.675, 0.946] and it CONTAINS
the 90% bar, so the design could not resolve its own threshold and the correct label
is UNRESOLVED-AT-THE-DESIGN'S-RESOLUTION. §P7-7 therefore pre-declared exactly one
estimator change (percentile -> BCa; see `gpd_tail`) and exactly one design change:

  * the BARS ARE UNTOUCHED -- >= 90% coverage, <= 3x worst understatement, as
    §P6-2(7) wrote them. Moving a bar after seeing a result would void the harness.
  * the DECISIVE SET is >= 125 FRESH comparisons, from n >= 0.9*0.1/((0.90-0.846)/2)^2
    = 125, the size that separates the observed 0.846 from the bar at 2 SE. Here it
    is 45 fresh cases x 3 targets = 135.
  * the LEGACY 27 are retained and re-scored under BCa, and reported SEPARATELY as a
    frozen subset that LICENSES NOTHING -- they are the data that motivated the fix,
    so a pass on them would be a fit to its own training set. `gpd_tail.
    assert_calibrated` refuses any artifact under 125 comparisons, which makes that
    separation structural rather than editorial.
  * ONE estimator is scored. The percentile endpoints ride along in every row as a
    labelled reference column so the swap's effect is visible, and they are not
    eligible to be the verdict. Reporting whichever of the two passed is the forking
    path §P7-7(a) names as inadmissible.

Run:  python -u -m engine.audit_gpd --jobs 16
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
    """THE LEGACY SET. The 9 cases x 3 targets = 27 comparisons of attempt 1.

    FROZEN, and it LICENSES NOTHING (§P7-7(b)). It is retained and re-scored under
    BCa on exactly the seeds it ran on before, so the change in the interval can be
    read comparison by comparison -- but it is the data that MOTIVATED the BCa swap,
    so scoring the fix on it would be scoring a fix on its own training set. It is
    reported separately, and `assert_calibrated` refuses any artifact with fewer
    than 125 comparisons precisely so this set can never quietly become the licence.
    """
    cases = []
    for i, (nd, pd_) in enumerate([(2000, 29.53), (2000, 13.66), (3000, 365.25)]):
        cases.append(glm_case(f"glm{i+1}", nd, pd_, 1000 + i))
    for i, (ne, kd) in enumerate([(2000, "linear"), (3000, "linear"),
                                  (2000, "phase")]):
        cases.append(mark_case(f"mark{i+1}", ne, 2000 + i, kd))
    for i, (nd, npd) in enumerate([(1000, 100), (1500, 150), (800, 60)]):
        cases.append(period_case(f"per{i+1}", nd, npd, 3000 + i))
    return cases


# --------------------------------------------------------------- the decisive set
# §P7-7(b), stated in numbers. Separating an observed 0.846 from the 0.90 bar at
# 2 SE needs n >= 0.9*0.1/((0.90-0.846)/2)^2 = 125 comparisons, so the decisive set
# is 45 FRESH cases x 3 target p = 135 -- the same three test kinds and the same
# three target levels as before, on new seeds and new substrate sizes that no
# previous run has touched.
#
# THE SUBSTRATE SIZES ARE SCALED, AND HERE IS THE ARITHMETIC. The brute-force arm
# costs, from the phase-2a artifact, about 0.116 s per GLM day, 0.078 s per mark
# event and 0.06-0.09 s per period case at N = 10^6. 135 comparisons at the legacy
# sizes would be ~6.5 CPU-hours; the sizes below hold the decisive set to roughly
# 4.5 CPU-hours, which at --jobs 16 is well inside one run. What is scaled is the
# SERIES LENGTH and the period grid -- never the estimator, never N_SURROGATES_GPD,
# never N_BRUTE, and never a bar.
FRESH_CASE_SEED = 51000                  # case construction; disjoint from 1000-3999
FRESH_RUN_SEED = 900017                  # run seeds; disjoint from the 700000 block

_FRESH_GLM = [(1000, 7.0), (1200, 13.66), (1400, 27.55), (1600, 29.53),
              (1800, 91.31), (1000, 14.77), (1200, 45.0), (1400, 182.6),
              (1600, 9.13), (1800, 27.55), (1100, 29.53), (1300, 60.0),
              (1500, 13.66), (1700, 121.75), (1900, 7.0)]
_FRESH_MARK = [(1200, "linear"), (1400, "linear"), (1600, "linear"),
               (1800, "linear"), (2000, "linear"), (2200, "linear"),
               (2400, "linear"), (1200, "phase"), (1400, "phase"),
               (1600, "phase"), (1800, "phase"), (2000, "phase"),
               (2200, "phase"), (2400, "phase"), (2600, "phase")]
_FRESH_PERIOD = [(800, 60), (900, 70), (1000, 80), (1100, 90), (1200, 100),
                 (1300, 110), (1400, 120), (1500, 130), (1600, 140),
                 (800, 100), (1000, 60), (1200, 150), (1400, 80),
                 (900, 120), (1100, 70)]


def fresh_roster():
    """THE DECISIVE SET. 45 fresh cases, 15 per kind, spanning the three kinds."""
    cases = []
    for i, (nd, pd_) in enumerate(_FRESH_GLM):
        cases.append(glm_case(f"fglm{i+1:02d}", nd, pd_, FRESH_CASE_SEED + i))
    for i, (ne, kd) in enumerate(_FRESH_MARK):
        cases.append(mark_case(f"fmark{i+1:02d}", ne,
                               FRESH_CASE_SEED + 100 + i, kd))
    for i, (nd, npd) in enumerate(_FRESH_PERIOD):
        cases.append(period_case(f"fper{i+1:02d}", nd, npd,
                                 FRESH_CASE_SEED + 200 + i))
    return cases


SET_DECISIVE = "decisive"
SET_LEGACY = "legacy"


# --------------------------------------------------------------- the top-up (§P7-9)
# The first BCa run scored 122 of 135 -- three short of the 125 the design requires,
# once §P7-9(1) resolved that count to mean SCORED comparisons. The remedy is to add
# fresh batches until the decisive set is the size it was always supposed to be.
#
# THE STOPPING RULE, DECLARED BEFORE THE BATCHES RUN AND IMPLEMENTED BELOW:
#   run whole batches in the fixed order declared here; after EACH COMPLETE batch,
#   stop iff (scored comparisons) >= 125.
# It reads ONE quantity: how many rows produced a GPD fit. It never reads coverage
# and never reads understatement. That is the whole point -- a rule that looked at
# the bars and stopped when they were clear would be optional stopping on the audit
# itself, which manufactures a pass out of sampling noise no matter how honest each
# individual comparison was. For the same reason the batches are run WHOLE: no
# peeking after each comparison and quitting on the good ones.
#
# The batch is 3 cases (one per kind) x 3 targets = 9 run comparisons. It is
# deliberately the smallest unit that still spans all three kinds, because
# worst-understatement is a MAXIMUM over the pooled set and drifts upward with n:
# every extra comparison is another draw at the 3x ceiling, so the top-up should
# overshoot 125 by as little as the design allows. Five batches are pre-declared so
# the sequence is fixed and public even though one is expected to suffice.
TOPUP_CASE_SEED = 61000
TOPUP_RUN_SEED = 1_300_000
_TOPUP_BATCHES = [
    [("glm", (1300, 21.0)), ("mark", (1900, "linear")), ("period", (1050, 85))],
    [("glm", (1550, 36.5)), ("mark", (2100, "phase")), ("period", (1250, 95))],
    [("glm", (1150, 73.0)), ("mark", (1700, "linear")), ("period", (950, 105))],
    [("glm", (1750, 11.0)), ("mark", (2300, "phase")), ("period", (1350, 65))],
    [("glm", (1450, 182.6)), ("mark", (1500, "linear")), ("period", (850, 115))],
]
N_TOPUP_BATCHES = len(_TOPUP_BATCHES)


def topup_set_name(batch):
    return f"topup{batch}"


def topup_roster(batch):
    """Batch `batch` (1-based) of the top-up. Same generators, new names and seeds."""
    spec = _TOPUP_BATCHES[int(batch) - 1]
    base = TOPUP_CASE_SEED + 1000 * int(batch)
    cases = []
    for i, (kind, params) in enumerate(spec):
        name = f"t{batch}{kind}{i+1}"
        if kind == "glm":
            cases.append(glm_case(name, params[0], params[1], base + i))
        elif kind == "mark":
            cases.append(mark_case(name, params[0], base + i, params[1]))
        else:
            cases.append(period_case(name, params[0], params[1], base + i))
    return cases


def roster_for(set_name):
    if set_name == SET_DECISIVE:
        return fresh_roster()
    if str(set_name).startswith("topup"):
        return topup_roster(int(str(set_name)[5:]))
    return roster()


def specs_for(set_name, targets):
    """The (set, case index, target p, seed) grid. Seeds are declared, not drawn."""
    cases = roster_for(set_name)
    if set_name == SET_DECISIVE:
        return [(set_name, i, tp, FRESH_RUN_SEED + 271 * i + 13 * j)
                for i, _c in enumerate(cases) for j, tp in enumerate(targets)]
    if str(set_name).startswith("topup"):
        b = int(str(set_name)[5:])
        return [(set_name, i, tp, TOPUP_RUN_SEED + 10000 * b + 271 * i + 13 * j)
                for i, _c in enumerate(cases) for j, tp in enumerate(targets)]
    # the legacy seeds, byte for byte the ones attempt 1 ran on, so the re-score is
    # the SAME 27 comparisons through a different interval and nothing else
    return [(set_name, i, tp, 700000 + 137 * i + 11 * j)
            for i, _c in enumerate(cases) for j, tp in enumerate(targets)]


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

    # THE SCORED INTERVAL IS BCa AND ONLY BCa. The percentile endpoints travel in
    # the row under their own names so the table can print what the superseded
    # estimator would have said, and they are scored into `covered_percentile` /
    # `understatement_x_percentile` -- REFERENCE COLUMNS. The gate, the verdict and
    # the artifact's `pass` read the BCa columns and nothing else (§P7-7(a)).
    lo, hi = g.get("p_ci_lower"), g.get("p_ci_upper")
    covered = bool(lo is not None and hi is not None and lo <= p_bf <= hi)
    understate = (float(p_bf / hi) if hi else float("inf"))
    plo, phi = g.get("p_ci_lower_percentile"), g.get("p_ci_upper_percentile")
    p_cov = bool(plo is not None and phi is not None and plo <= p_bf <= phi)
    p_us = (float(p_bf / phi) if phi else float("inf"))
    bca = g.get("bca") or {}
    return {"case": case["name"], "kind": case["kind"], "detail": case["detail"],
            "target_p": target_p, "s_obs": s_obs, "p_mc_2000": p_mc,
            "p_method": g["p_method"], "p_point": g.get("p_point"),
            "p_ci_lower": lo, "p_ci_upper": hi, "p_brute_force": p_bf,
            "n_brute": N_BRUTE, "covered": covered,
            "understatement_x": understate,
            "ci_method": ("BCa" if (bca and not bca.get("fallback"))
                          else ("percentile-fallback" if bca else None)),
            "bca_z0": bca.get("z0"), "bca_a": bca.get("a"),
            "bca_alpha_lo": bca.get("alpha_lo"), "bca_alpha_hi": bca.get("alpha_hi"),
            "bca_jackknife": bca.get("jackknife"), "bca_n_jack": bca.get("n_jack"),
            "bca_fallback": bca.get("fallback"), "bca_notes": bca.get("notes"),
            "p_ci_lower_percentile": plo, "p_ci_upper_percentile": phi,
            "covered_percentile": p_cov, "understatement_x_percentile": p_us,
            "p_ad": (g.get("ad") or {}).get("p_ad"),
            "xi": g.get("xi"), "xi_stable": (g.get("xi_stability") or {}).get("pass"),
            "reason": g.get("reason"), "seconds_gpd": round(t_gpd, 1),
            "seconds_brute": round(t_bf, 1)}


def _job(spec):
    os.environ["OMP_NUM_THREADS"] = "1"
    set_name, idx, target_p, seed = spec
    r = run_case(roster_for(set_name)[idx], target_p, seed)
    r["set"] = set_name
    return r


def score(rows, column="understatement_x", cov_key="covered"):
    """The gate, computed once and applied to whichever column is named.

    Only rows that produced a GPD FIT are scored. A gated-out row is the estimator
    REFUSING to extrapolate, which is the behaviour §P6-2(1) asks for, and scoring a
    refusal as a coverage failure would penalise exactly the right answer.
    """
    fitted = [r for r in rows if r["p_method"] == gpd_tail.P_GPD]
    n_cov = sum(1 for r in fitted if r[cov_key])
    cov = (n_cov / len(fitted)) if fitted else 0.0
    worst = max((r[column] for r in fitted), default=float("inf"))
    ok_cov = bool(cov >= COVERAGE_BAR)
    ok_us = bool(worst <= UNDERSTATE_BAR)
    # n_scored IS n_fitted -- the rows that entered the two bars. It is named
    # separately because §P7-9(1) makes it the quantity the design is stated in, and
    # a number that licenses an instrument should not have to be inferred from a
    # synonym.
    return {"n_comparisons": len(rows), "n_fitted": len(fitted),
            "n_scored": len(fitted),
            "n_covered": n_cov, "coverage": cov, "worst_understatement": worst,
            "coverage_bar": COVERAGE_BAR, "understatement_bar": UNDERSTATE_BAR,
            "coverage_pass": ok_cov, "understatement_pass": ok_us,
            "pass": bool(ok_cov and ok_us)}


def by_kind(rows, cov_key="covered", column="understatement_x"):
    out = {}
    for k in sorted({r["kind"] for r in rows}):
        out[k] = score([r for r in rows if r["kind"] == k], column, cov_key)
    return out


def _print_table(rows, title):
    print()
    print(f"### {title}")
    print("| case | kind | target p | p_MC@2000 | p_method | GPD 95% BCa CI | "
          "brute force @1e6 | covered | understatement | (ref) percentile CI | "
          "(ref) cov |")
    print("| --- | --- | ---: | ---: | --- | --- | ---: | :---: | ---: | --- "
          "| :---: |")
    for r in rows:
        ci = ("--" if r["p_ci_lower"] is None
              else f"[{r['p_ci_lower']:.3g}, {r['p_ci_upper']:.3g}]")
        pci = ("--" if r.get("p_ci_lower_percentile") is None
               else f"[{r['p_ci_lower_percentile']:.3g}, "
                    f"{r['p_ci_upper_percentile']:.3g}]")
        us = ("--" if not np.isfinite(r["understatement_x"])
              else f"{r['understatement_x']:.2f}x")
        flag = "*" if r.get("bca_fallback") else ""
        print(f"| {r['case']}{flag} | {r['kind']} | {r['target_p']:.0e} | "
              f"{r['p_mc_2000']:.3g} | {r['p_method']} | {ci} | "
              f"{r['p_brute_force']:.3g} | {'YES' if r['covered'] else 'no'} | "
              f"{us} | {pci} | "
              f"{'YES' if r.get('covered_percentile') else 'no'} |")


def _print_verdict(label, s, extra=""):
    # §P7-9(3): understatement leads, WITH its n, because it is a MAXIMUM over the
    # scored set and therefore drifts upward as the set grows -- a number quoted
    # without the n it was taken over is not interpretable. Coverage is a mean and
    # comes second.
    print()
    print(f"[{label}] SCORED comparisons: {s['n_scored']} of {s['n_comparisons']} "
          f"run ({s['n_comparisons'] - s['n_scored']} gated out and reported "
          f"UNRESOLVED, which is the estimator refusing rather than failing; those "
          f"rows entered NEITHER bar)")
    print(f"[{label}] GATE 1  worst understatement: "
          f"{s['worst_understatement']:.2f}x over n = {s['n_scored']} scored "
          f"(bar <= {UNDERSTATE_BAR:.0f}x)  -> "
          f"{'PASS' if s['understatement_pass'] else 'FAIL'}")
    print(f"[{label}] GATE 2  coverage            : {s['n_covered']}/"
          f"{s['n_scored']} = {s['coverage']:.1%} (bar >= {COVERAGE_BAR:.0%})"
          f"       -> {'PASS' if s['coverage_pass'] else 'FAIL'}")
    print(f"[{label}] VERDICT: {'ACCEPT' if s['pass'] else 'REJECT'} {extra}")


def main(argv=None):
    ap = argparse.ArgumentParser("engine.audit_gpd")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--out", default=os.path.join("engine", "out",
                                                  "audit_gpd_bca.json"))
    ap.add_argument("--targets", type=float, nargs="*", default=list(TARGET_PS))
    ap.add_argument("--sets", nargs="*", default=[SET_DECISIVE, SET_LEGACY],
                    choices=[SET_DECISIVE, SET_LEGACY])
    ap.add_argument("--no-topup", action="store_true",
                    help="skip the §P7-9(2) top-up batches (diagnostics only; the "
                         "result cannot license)")
    args = ap.parse_args(argv)

    min_scored = gpd_tail.GPD_MIN_CALIBRATION_COMPARISONS

    def _run(specs, jobs):
        if jobs > 1:
            import concurrent.futures as cf
            import multiprocessing as mp
            ctx = mp.get_context("spawn")
            with cf.ProcessPoolExecutor(max_workers=jobs, mp_context=ctx) as ex:
                return list(ex.map(_job, specs))
        return [_job(s) for s in specs]

    specs = []
    for sn in args.sets:
        specs.extend(specs_for(sn, args.targets))
    n_dec = sum(1 for s in specs if s[0] == SET_DECISIVE)
    print(f"§P6-2(7) calibration, attempt 2 (BCa, §P7-7(a)) -- "
          f"GPD at N = {N_SURROGATES_GPD} vs brute force at N = {N_BRUTE:,}")
    print(f"  DECISIVE set  : {len(fresh_roster())} fresh cases x "
          f"{len(args.targets)} target p = {n_dec} comparisons RUN, plus "
          f"pre-declared top-up batches until >= {min_scored} are SCORED "
          f"(§P7-7(b) as resolved by §P7-9(1))")
    print(f"  LEGACY set    : {len(specs) - n_dec} comparisons, re-scored under "
          f"BCa on the original seeds, REPORTED SEPARATELY AND LICENSING NOTHING")
    print(f"  bars UNCHANGED: coverage >= {COVERAGE_BAR:.0%}, worst "
          f"understatement <= {UNDERSTATE_BAR:.0f}x")
    print(f"  kinds: " + ", ".join(sorted({c['kind'] for c in fresh_roster()})))
    t0 = time.perf_counter()
    rows = _run(specs, args.jobs)

    # §P7-9(2): the top-up, under the stopping rule declared at _TOPUP_BATCHES.
    # The condition reads n_scored and NOTHING else -- not coverage, not
    # understatement -- and each batch is run to completion before it is consulted.
    topup_log = []
    if SET_DECISIVE in args.sets and not args.no_topup:
        for b in range(1, N_TOPUP_BATCHES + 1):
            dec_now = [r for r in rows if r["set"] != SET_LEGACY]
            n_scored = score(dec_now)["n_scored"]
            topup_log.append({"before_batch": b, "n_scored": n_scored,
                              "min_required": min_scored,
                              "ran": bool(n_scored < min_scored)})
            if n_scored >= min_scored:
                print(f"\nSTOPPING RULE: {n_scored} scored >= {min_scored} "
                      f"required -- no further batches. (The rule read the scored "
                      f"count only; coverage and understatement were not "
                      f"consulted.)")
                break
            sn = topup_set_name(b)
            bs = specs_for(sn, args.targets)
            print(f"\nSTOPPING RULE: {n_scored} scored < {min_scored} required "
                  f"-> running top-up batch {b} WHOLE ({len(bs)} comparisons, "
                  f"{len(roster_for(sn))} fresh cases, one per kind)")
            rows.extend(_run(bs, args.jobs))
    elapsed = time.perf_counter() - t0

    rows.sort(key=lambda r: (r["kind"], r["case"], r["target_p"]))
    # the decisive set is the POOLED set: the fresh 45-case run plus every top-up
    # batch the stopping rule called for. Both arms are scored over the pool.
    dec = [r for r in rows if r["set"] != SET_LEGACY]
    leg = [r for r in rows if r["set"] == SET_LEGACY]

    if dec:
        _print_table(dec, "DECISIVE SET -- fresh comparisons, BCa. This is the "
                          "licensing evidence. (* = BCa degenerate, percentile "
                          "fallback, flagged)")
    if leg:
        _print_table(leg, "LEGACY SET -- the 27 attempt-1 comparisons re-scored "
                          "under BCa. FROZEN SUBSET, LICENSES NOTHING (§P7-7(b)): "
                          "it is the data that motivated the fix.")

    s_dec = score(dec) if dec else None
    s_leg = score(leg) if leg else None
    # REFERENCE ONLY: what the superseded percentile interval scores on the same
    # comparisons. Printed so the swap's effect is visible; it is NOT a candidate
    # verdict and it cannot license anything (§P7-7(a)).
    s_dec_ref = (score(dec, "understatement_x_percentile", "covered_percentile")
                 if dec else None)

    if s_dec:
        _print_verdict("DECISIVE", s_dec,
                       f"-- {elapsed:.0f}s wall, --jobs {args.jobs}")
        print("[DECISIVE] per kind: " + "; ".join(
            f"{k}: {v['n_covered']}/{v['n_fitted']} = {v['coverage']:.1%}, "
            f"worst {v['worst_understatement']:.2f}x"
            for k, v in by_kind(dec).items()))
        print(f"[DECISIVE] (reference, NOT a verdict) the superseded percentile "
              f"interval on the same comparisons: "
              f"{s_dec_ref['n_covered']}/{s_dec_ref['n_fitted']} = "
              f"{s_dec_ref['coverage']:.1%} coverage, worst "
              f"{s_dec_ref['worst_understatement']:.2f}x")
        nfb = sum(1 for r in dec if r.get("bca_fallback"))
        print(f"[DECISIVE] BCa degenerate -> percentile fallback on {nfb} of "
              f"{len(dec)} rows")
    if s_leg:
        _print_verdict("LEGACY (licenses nothing)", s_leg)
        print("[LEGACY] per kind: " + "; ".join(
            f"{k}: {v['n_covered']}/{v['n_fitted']} = {v['coverage']:.1%}, "
            f"worst {v['worst_understatement']:.2f}x"
            for k, v in by_kind(leg).items()))

    # The artifact's TOP-LEVEL verdict is the DECISIVE set's and nothing else --
    # `load_calibration` reads these keys, so this is where "the legacy subset
    # licenses nothing" has to be true in the bytes, not only in the prose.
    # §P7-9(1): SCORED, not run.
    ok = bool(s_dec and s_dec["pass"] and s_dec["n_scored"] >= min_scored)
    doc = {"rule": "§P6-2(7)", "ruling": "§P7-7", "attempt": 2,
           "ci_method": "BCa", "set": SET_DECISIVE,
           "estimator_note": ("BCa replaces the percentile interval everywhere; the "
                              "percentile columns in these rows are a labelled "
                              "reference and were not eligible to be the verdict"),
           "n_cases": len(fresh_roster()), "elapsed_s": elapsed,
           "min_scored_required": min_scored,
           "scored_count_basis": ("SCORED comparisons (rows that produced a GPD "
                                  "fit and entered both bars), per §P7-9(1); "
                                  "gated-out rows are excluded"),
           "topup": {"rule": ("run whole pre-declared batches until n_scored >= "
                              "125; the condition reads the scored count only and "
                              "never coverage or understatement (§P7-9(2))"),
                     "batches_declared": N_TOPUP_BATCHES,
                     "batches_run": sum(1 for t in topup_log if t["ran"]),
                     "log": topup_log},
           "pass": ok, "rows": dec,
           "decisive_by_kind": by_kind(dec) if dec else {},
           "percentile_reference": s_dec_ref,
           "legacy": {"note": ("the 27 attempt-1 comparisons re-scored under BCa; "
                               "frozen subset, licenses nothing (§P7-7(b))"),
                      "summary": s_leg, "by_kind": by_kind(leg) if leg else {},
                      "rows": leg}}
    doc.update(s_dec or {"n_comparisons": 0, "n_fitted": 0, "n_scored": 0,
                         "coverage": 0.0,
                         "worst_understatement": float("inf"),
                         "coverage_bar": COVERAGE_BAR,
                         "understatement_bar": UNDERSTATE_BAR})
    doc["pass"] = ok
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1)
    print()
    if ok:
        print(f"§P7-7 OUTCOME: LICENSED -- --gpd may run. Worst understatement "
              f"{s_dec['worst_understatement']:.2f}x over n = {s_dec['n_scored']} "
              f"scored; coverage {s_dec['coverage']:.1%}.")
    elif s_dec and s_dec["n_scored"] < min_scored:
        print(f"§P7-7 OUTCOME: INCOMPLETE -- only {s_dec['n_scored']} scored of "
              f"{min_scored} required, and every declared top-up batch has been "
              f"spent. --gpd stays blocked.")
    else:
        failed = ([] + (["worst understatement "
                         f"{s_dec['worst_understatement']:.2f}x > "
                         f"{UNDERSTATE_BAR:.0f}x over n = {s_dec['n_scored']}"]
                        if s_dec and not s_dec["understatement_pass"] else [])
                  + ([f"coverage {s_dec['coverage']:.1%} < {COVERAGE_BAR:.0%}"]
                     if s_dec and not s_dec["coverage_pass"] else []))
        print(f"§P7-7(c) OUTCOME: REJECTED-STANDING -- --gpd stays blocked. "
              f"Failed on: {'; '.join(failed) or 'no decisive set'}. This was the "
              f"last estimator change; a third requires a §P6-2 amendment argued "
              f"on its own merits.")
    print(f"written -> {args.out}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
