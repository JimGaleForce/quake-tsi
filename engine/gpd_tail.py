"""GPD tail extrapolation for censored Monte Carlo p-values -- §P6-2, all 7 rules.

WHAT THIS IS FOR. A surrogate-based p-value is censored at 1/(N+1): once zero of N
surrogates reach the observed statistic, the Monte Carlo has stopped measuring and
is only reporting its own budget. Knijnenburg et al. (2009) is the field's standard
answer -- fit a generalised Pareto distribution to the exceedances over a high
threshold and read the tail probability off the fit.

WHY IT IS FENCED. §P6-2 grants this PASS-WITH-CONDITIONS and the conditions are
severe, because this program has already voided one estimator for the identical sin
(§P5-8 killed exp(dtau/A sigma) for being fitted at 7-347 kPa and quoted at 0.3 kPa).
A GPD p quoted a decade below the data that fitted it is the same manufacture in a
different coordinate. Every one of the seven rules is implemented here or, where it
is a pipeline rule rather than an estimator rule, is implemented at the single call
site and named in this docstring so the trace is not lost:

  RULE 1  goodness-of-fit gate, declared before the run: threshold u = the declared
          top 10% of surrogate statistics, GPD by MLE, and BOTH of
            (a) Anderson-Darling of the exceedances against the FITTED GPD,
                p_AD >= 0.05                                  -> `ad_gate`
            (b) shape stability: refit at top {5%, 10%, 20%} and require xi to
                agree across all three within its own bootstrap CI
                                                              -> `xi_stability`
          Fail either -> NO extrapolation; p is reported UNRESOLVED and the test is
          INELIGIBLE FOR BH REJECTION (enforced at the call site, see
          `mine_session.apply_p_methods` / the BH assembly).
  RULE 2  never a point estimate: the reported number is `p_ci_upper`, the upper end
          of the 95% bootstrap CI. `p_point` is carried for diagnosis ONLY and no
          downstream consumer may read it.                    -> `gpd_tail_p`
          The interval is BCa (bias-corrected and accelerated) -- see THE INTERVAL,
          below. There is exactly ONE interval estimator.
  RULE 3  one-decade cap: nothing below 0.1/(N_max+1). Below it -> UNRESOLVED,
          never "very significant".                           -> `gpd_tail_p`
  RULE 4  FORBIDDEN on the all-shifts enumeration null -- ~8,000 circular-shift
          statistics are strongly dependent across neighbouring shifts and the
          iid-exceedances assumption fails hardest exactly there. Structural: the
          null type is a REQUIRED argument and a forbidden one RAISES.
                                                              -> `assert_null_permitted`
  RULE 5  mandatory label p_method in {MC_RESOLVED, GPD_EXTRAPOLATED, UNRESOLVED}
          on every row of stubs.json and report.md plus a three-number census in the
          banner.                                             -> `P_METHODS`, `census`
  RULE 6  the confirmation rule: BH runs with GPD rows at CI-upper; a GPD survivor
          is a CANDIDATE-REQUIRES-BRUTE-FORCE, not a stub, until a targeted
          single-test MC at N >= 10/p_gpd resolves it directly. Only the brute-force
          p may be written to stubs.json.                     -> `brute_force_n`,
          and `mine_session.confirm_gpd_candidates`
  RULE 7  calibration before default: >= 15 representative tests spanning the three
          test kinds, GPD-at-2,000-surrogates CI vs brute force at N = 1e6; accept
          only if the CI covers in >= 90% of cases and never understates by more
          than 3x.                                            -> `engine/audit_gpd.py`
          §P7-7(b) tightened the DESIGN (not the bars): the licensing set must be
          >= 125 FRESH comparisons, because at n = 26 the harness could not resolve
          its own 90% threshold. §P7-9(1) resolved that count to mean SCORED
          comparisons -- rows the gates threw out entered neither bar and carry no
          power.                                              -> `assert_calibrated`

THE INTERVAL: BCa, AND WHY IT IS THE ONLY ONE (§P7-7(a)). Attempt 1 used the
percentile bootstrap and scored UNRESOLVED-AT-DESIGN-RESOLUTION: 22/26 = 84.6%
coverage, whose exact 95% interval [0.675, 0.946] contains the 90% bar, so 26
comparisons could not tell a broken estimator from a working one. Every one of the
four misses was upper-edge under-coverage of a percentile interval on a strongly
right-skewed tail quantity, which is exactly the failure mode bias-correction and
acceleration exist to repair, so §P7-7 pre-declared ONE swap:

    p_lo = G^-1(alpha1), p_hi = G^-1(alpha2),  G = the bootstrap ECDF of p*,
    alpha_j = Phi( z0 + (z0 + z_j) / (1 - a (z0 + z_j)) ),   z_j = the normal
                                                             quantiles of the level
    z0 = Phi^-1( #{p*_b < p_hat} / B )        bias correction
    a  = sum d_i^3 / (6 (sum d_i^2)^{3/2}),   d_i = mean(p_(.)) - p_(i)
                                              over the JACKKNIFE of the surrogate set

BCa REPLACES the percentile interval everywhere -- `p_ci_lower`/`p_ci_upper`, and
therefore `p_report` and everything BH sees, are BCa. The percentile endpoints
survive ONLY as `p_ci_lower_percentile`/`p_ci_upper_percentile`, a labelled
reference column for the audit table. They are NOT selectable and no consumer may
read them: running BCa and percentile and studentised and reporting whichever
passes is the forking path §P7-7(a) names as what would make the swap inadmissible.
The one exception is DEGENERATE ARITHMETIC -- z0 or a non-finite, the BCa
denominator 1 - a(z0+z) non-positive, or alpha1 >= alpha2 -- where there is no BCa
interval to report at all; there the percentile endpoints are used as the defined
fallback and the row is FLAGGED (`bca.fallback`), so a fallback can never be
mistaken for a BCa result or chosen because it looked better.

BOTH resamplings re-derive u, zeta, xi and beta FROM SCRATCH. The bootstrap
resamples the whole surrogate set and the jackknife deletes from it, and both route
through the single function `_tail_p_from_sample`, so neither can condition on a
quantity estimated off the full sample. That conditioning was the phase-2a defect
(the old interval priced only GPD parameter uncertainty, holding u and zeta fixed);
its repair was a CORRECTNESS fix and it stands, and the shared helper is what makes
the jackknife's estimand identical to the bootstrap's by construction rather than by
inspection.

THE ANDERSON-DARLING METHOD, STATED (the gotcha `scipy.stats.anderson` does not
cover). scipy's `anderson` supports norm/expon/logistic/gumbel/extreme1 and has no
GPD case, and there is no closed-form null distribution for A^2 when the GPD's two
parameters are ESTIMATED from the same sample -- the standard tables assume known
parameters and are badly anti-conservative otherwise. So A^2 is computed here
directly from the fitted CDF,

    A^2 = -n - (1/n) * sum_i (2i-1) * [ ln z_(i) + ln(1 - z_(n+1-i)) ],
    z_(i) = F_gpd(x_(i); xi_hat, beta_hat) on the sorted exceedances,

and its p-value comes from a PARAMETRIC BOOTSTRAP that reproduces the estimation:
each replicate draws n points from the FITTED GPD, REFITS by MLE, and recomputes
A^2, so the null distribution carries the same parameter-estimation shrinkage as the
observed statistic. p_AD = (1 + #{A^2_b >= A^2_obs}) / (1 + B). That is a Stephens
(1976)-style parametric-bootstrap GOF test, not a table lookup.
"""

import hashlib
import json
import math
import os

import numpy as np
from scipy import stats

# ------------------------------------------------------------ declared constants
# DECLARED BEFORE THE RUN (§P6-2(1)). These are not tunables: changing any of them
# changes the estimator, and the estimator is calibrated as a whole by audit_gpd.
GPD_THRESHOLD_Q = 0.90          # u = the declared top 10% of surrogate statistics
GPD_STABILITY_QS = (0.95, 0.90, 0.80)   # top 5%, 10%, 20% -- rule 1(b)
GPD_AD_ALPHA = 0.05             # rule 1(a): p_AD >= 0.05
GPD_CI_LEVEL = 0.95             # rule 2: the 95% bootstrap CI
GPD_N_BOOT = 400                # bootstrap replicates for the CI and for xi's CI
GPD_N_AD_BOOT = 200             # parametric-bootstrap replicates for p_AD
GPD_MIN_EXCEEDANCES = 25        # below this an MLE tail fit is not an estimate
GPD_DECADES = 1.0               # rule 3: one decade below the MC floor and no more
GPD_CONFIRM_FACTOR = 10.0       # rule 6: brute force at N >= 10 / p_gpd

# BCa (§P7-7(a)). The jackknife that supplies the acceleration is PER-POINT --
# delete one surrogate, re-derive u, zeta, xi, beta, recompute p -- up to
# GPD_JACKKNIFE_MAX_POINTS. It is not an approximation and at the engine's declared
# N = 2,000 it costs ~8 s per test against a brute-force reference arm of 130-350 s,
# so the exact version is affordable and is what runs. Above the cap the per-point
# cost is linear in n and would start to dominate, so the jackknife becomes a
# GROUPED (delete-d) jackknife over GPD_JACKKNIFE_GROUPS contiguous blocks of the
# surrogate index -- the surrogates are iid draws, so contiguous blocks are an
# arbitrary partition and carry no ordering. Which one ran is recorded per row in
# `bca.jackknife`, so the table never has to be trusted on this point.
GPD_JACKKNIFE_MAX_POINTS = 2500
GPD_JACKKNIFE_GROUPS = 250
# §P7-7(b), as resolved by §P7-9: the licensing set must be >= 125 SCORED
# comparisons -- not 125 run. A comparison that the gates threw out contributed
# nothing to either bar, so counting it toward the design's power would be counting
# an absence as evidence. The first BCa run made exactly this mistake in my favour:
# 135 run, 13 gated out, 122 scored, and it was three short. Bars unchanged, design
# enlarged -- 26 comparisons cannot resolve a 90% threshold.
GPD_MIN_CALIBRATION_COMPARISONS = 125

P_MC_RESOLVED = "MC_RESOLVED"
P_GPD = "GPD_EXTRAPOLATED"
P_UNRESOLVED = "UNRESOLVED"
P_METHODS = (P_MC_RESOLVED, P_GPD, P_UNRESOLVED)

CANDIDATE_LABEL = "CANDIDATE-REQUIRES-BRUTE-FORCE"

# RULE 4, structurally. The all-shifts enumeration null is FORBIDDEN as a GPD
# substrate; so is anything else whose "surrogates" are a deterministic sweep of a
# single record. Only INDEPENDENTLY GENERATED surrogates are permitted.
FORBIDDEN_NULL_TYPES = frozenset({
    "all_shifts", "circular_shift", "circular-shift", "enumeration",
    "all-shifts", "shift_enumeration",
})
PERMITTED_NULL_TYPES = frozenset({
    "block_bootstrap", "block-bootstrap", "ar1", "permutation", "etas_sim",
    "etas-sim", "parametric",
})


class ForbiddenNullError(ValueError):
    """Raised when a GPD fit is attempted on a null §P6-2(4) forbids."""


def assert_null_permitted(null_type):
    """RULE 4. Refuse to extrapolate off a dependent-surrogate null. Structural.

    Not a warning and not a flag: neighbouring circular shifts of one record share
    almost all of their data, the exceedances are the opposite of iid, and the GPD's
    entire licence is iid exceedances. An unknown null type is ALSO refused -- the
    permitted list is a declaration, so a null nobody has thought about defaults to
    forbidden rather than to permitted.
    """
    if null_type is None:
        raise ForbiddenNullError(
            "GPD tail extrapolation requires an explicit null_type (§P6-2(4)); "
            "None is not a declaration.")
    key = str(null_type).strip().lower()
    if key in FORBIDDEN_NULL_TYPES:
        raise ForbiddenNullError(
            f"§P6-2(4): GPD tail extrapolation is FORBIDDEN on the '{null_type}' "
            f"null. Circular-shift statistics are strongly dependent across "
            f"neighbouring shifts and the GPD's iid-exceedances assumption fails "
            f"hardest exactly there. Permitted substrates: "
            f"{', '.join(sorted(PERMITTED_NULL_TYPES))}.")
    if key not in PERMITTED_NULL_TYPES:
        raise ForbiddenNullError(
            f"§P6-2(4): null type '{null_type}' is not on the declared permitted "
            f"list {sorted(PERMITTED_NULL_TYPES)}. An undeclared null is refused, "
            f"not assumed independent.")
    return key


# ------------------------------------------------------------------- GPD fitting
def fit_gpd(exceedances):
    """MLE fit of GPD(xi, beta) to exceedances over a threshold, location fixed at 0.

    Returns (xi, beta). scipy parameterises genpareto as (c, loc, scale) with c = xi;
    `floc=0` is not an option here but a requirement -- the exceedances are already
    measured from u, and letting loc float would re-estimate the threshold that rule
    1 declared.
    """
    x = np.asarray(exceedances, dtype=np.float64).ravel()
    x = x[np.isfinite(x) & (x > 0)]
    if x.size < 2:
        raise ValueError("fit_gpd: need at least 2 positive exceedances")
    xi, _loc, beta = stats.genpareto.fit(x, floc=0.0)
    return float(xi), float(beta)


def ad_statistic(exceedances, xi, beta):
    """Anderson-Darling A^2 of exceedances against the FITTED GPD. See module docstring."""
    x = np.sort(np.asarray(exceedances, dtype=np.float64).ravel())
    n = x.size
    if n < 2:
        return float("nan")
    z = stats.genpareto.cdf(x, xi, loc=0.0, scale=beta)
    # clip only against exact 0/1, which would send a finite statistic to infinity
    # purely through floating point; the clip is 1e-12, far below any real signal.
    z = np.clip(z, 1e-12, 1.0 - 1e-12)
    i = np.arange(1, n + 1)
    return float(-n - np.sum((2 * i - 1) * (np.log(z) + np.log(1.0 - z[::-1]))) / n)


def ad_gate(exceedances, xi, beta, rng, n_boot=GPD_N_AD_BOOT,
            alpha=GPD_AD_ALPHA):
    """RULE 1(a). p_AD by parametric bootstrap that REFITS on every replicate.

    Refitting matters: A^2 computed against a fitted GPD is systematically smaller
    than A^2 against the true one, so a null distribution simulated WITHOUT refitting
    would be shifted up and the gate would pass almost anything.
    """
    x = np.asarray(exceedances, dtype=np.float64).ravel()
    n = x.size
    a2_obs = ad_statistic(x, xi, beta)
    if not np.isfinite(a2_obs):
        return {"a2": float("nan"), "p_ad": 0.0, "n_boot": 0, "pass": False,
                "alpha": float(alpha),
                "method": "parametric bootstrap, refit per replicate"}
    ge = 0
    n_ok = 0
    for _ in range(int(n_boot)):
        xb = stats.genpareto.rvs(xi, loc=0.0, scale=beta, size=n,
                                 random_state=rng)
        xb = xb[np.isfinite(xb) & (xb > 0)]
        if xb.size < 2:
            continue
        try:
            xib, _l, betab = stats.genpareto.fit(xb, floc=0.0)
        except Exception:
            continue
        a2b = ad_statistic(xb, xib, betab)
        if not np.isfinite(a2b):
            continue
        n_ok += 1
        ge += int(a2b >= a2_obs)
    p_ad = (1.0 + ge) / (1.0 + n_ok) if n_ok else 0.0
    return {"a2": float(a2_obs), "p_ad": float(p_ad), "n_boot": int(n_ok),
            "pass": bool(p_ad >= alpha), "alpha": float(alpha),
            "method": ("Anderson-Darling against the fitted GPD; p by parametric "
                       "bootstrap with MLE refit on every replicate (scipy.stats."
                       "anderson has no GPD case and its tables assume known "
                       "parameters)")}


def _xi_ci(exceedances, rng, n_boot=GPD_N_BOOT, level=GPD_CI_LEVEL):
    """Nonparametric bootstrap CI for xi at one threshold."""
    x = np.asarray(exceedances, dtype=np.float64).ravel()
    n = x.size
    xis = []
    for _ in range(int(n_boot)):
        xb = x[rng.integers(0, n, size=n)]
        try:
            xib, _l, _b = stats.genpareto.fit(xb, floc=0.0)
        except Exception:
            continue
        if np.isfinite(xib):
            xis.append(float(xib))
    if len(xis) < 10:
        return (float("nan"), float("nan"))
    lo = (1.0 - level) / 2.0 * 100.0
    return (float(np.percentile(xis, lo)), float(np.percentile(xis, 100.0 - lo)))


def xi_stability(S, rng, quantiles=GPD_STABILITY_QS, n_boot=GPD_N_BOOT,
                 level=GPD_CI_LEVEL, min_exceedances=GPD_MIN_EXCEEDANCES):
    """RULE 1(b). Refit at top {5%, 10%, 20%} and require xi to agree within its CI.

    "Agree within its own bootstrap CI" is read as: the three CIs have a point in
    common (equivalently max(lo) <= min(hi)). That is the weakest reading that still
    means agreement, and it is the one a reader would check by eye off the printed
    table. A shape parameter that walks with the threshold means the tail is not GPD
    over the range being used, which is precisely the failure that manufactures a
    p-value out of the fit's curvature.
    """
    S = np.asarray(S, dtype=np.float64).ravel()
    rows = []
    for qq in quantiles:
        u = float(np.quantile(S, qq))
        exc = S[S > u] - u
        if exc.size < min_exceedances:
            rows.append({"q": float(qq), "u": u, "n_exceed": int(exc.size),
                         "xi": float("nan"), "ci": [float("nan"), float("nan")],
                         "ok": False})
            continue
        try:
            xi, beta = fit_gpd(exc)
        except Exception:
            rows.append({"q": float(qq), "u": u, "n_exceed": int(exc.size),
                         "xi": float("nan"), "ci": [float("nan"), float("nan")],
                         "ok": False})
            continue
        lo, hi = _xi_ci(exc, rng, n_boot=n_boot, level=level)
        rows.append({"q": float(qq), "u": u, "n_exceed": int(exc.size),
                     "xi": float(xi), "beta": float(beta), "ci": [lo, hi],
                     "ok": bool(np.isfinite(lo) and np.isfinite(hi))})
    ok = all(r["ok"] for r in rows) and len(rows) == len(quantiles)
    if ok:
        lo_max = max(r["ci"][0] for r in rows)
        hi_min = min(r["ci"][1] for r in rows)
        agree = bool(lo_max <= hi_min)
    else:
        lo_max = hi_min = float("nan")
        agree = False
    return {"rows": rows, "pass": bool(agree), "ci_overlap_lo": lo_max,
            "ci_overlap_hi": hi_min,
            "rule": "xi CIs at top {5,10,20}% must share a common point"}


# ------------------------------------------------- the estimand, one recipe only
def _tail_p_from_sample(S, s_obs, threshold_q, min_exceedances=2):
    """p = zeta * SF_GPD(s_obs - u), with u, zeta, xi, beta ALL re-derived from S.

    THE single recipe. The point estimate, every bootstrap replicate and every
    jackknife replicate call this and nothing else, so no replicate can condition on
    a quantity estimated off the full sample -- the estimand is the extrapolated p
    and each replicate must re-estimate the whole chain that produces it. Returns
    None (never raises, never a partial answer) whenever the chain cannot be run:
    too few exceedances, an observed statistic that no longer clears the replicate's
    own threshold, a failed MLE, or a non-finite p.
    """
    S = np.asarray(S, dtype=np.float64).ravel()
    n = S.size
    if n < 1:
        return None
    u = float(np.quantile(S, threshold_q))
    exc = S[S > u] - u
    if exc.size < int(min_exceedances) or float(s_obs) <= u:
        return None
    try:
        xi, _loc, beta = stats.genpareto.fit(exc, floc=0.0)
    except Exception:
        return None
    zeta = exc.size / float(n)
    p = float(zeta * stats.genpareto.sf(float(s_obs) - u, xi, loc=0.0, scale=beta))
    if not np.isfinite(p):
        return None
    return {"p": p, "u": u, "zeta": float(zeta), "xi": float(xi),
            "beta": float(beta), "n_exceed": int(exc.size)}


def _bootstrap_ps(S, s_obs, rng, threshold_q=GPD_THRESHOLD_Q, n_boot=GPD_N_BOOT):
    """Nonparametric bootstrap of the extrapolated p, resampling the WHOLE set."""
    S = np.asarray(S, dtype=np.float64).ravel()
    n = S.size
    ps = []
    for _ in range(int(n_boot)):
        r = _tail_p_from_sample(S[rng.integers(0, n, size=n)], s_obs, threshold_q)
        if r is not None:
            ps.append(r["p"])
    return np.asarray(ps, dtype=np.float64)


def _jackknife_ps(S, s_obs, threshold_q=GPD_THRESHOLD_Q,
                  max_points=GPD_JACKKNIFE_MAX_POINTS,
                  groups=GPD_JACKKNIFE_GROUPS):
    """Leave-one-out (or grouped delete-d) jackknife of the extrapolated p.

    Returns (values, mode). `mode` is reported so the acceleration's provenance is
    on the row: "leave-one-out" is the exact jackknife, "grouped-<g>" is the
    delete-d jackknife used only when n exceeds the declared cap.
    """
    S = np.asarray(S, dtype=np.float64).ravel()
    n = S.size
    if n < 3:
        return np.asarray([], dtype=np.float64), "none"
    if n <= int(max_points):
        blocks = [(i, i + 1) for i in range(n)]
        mode = "leave-one-out"
    else:
        g = int(min(groups, n))
        edges = np.linspace(0, n, g + 1).astype(int)
        blocks = [(int(edges[k]), int(edges[k + 1])) for k in range(g)
                  if edges[k + 1] > edges[k]]
        mode = f"grouped-{len(blocks)}"
    vals = []
    for lo, hi in blocks:
        sub = np.concatenate([S[:lo], S[hi:]])
        r = _tail_p_from_sample(sub, s_obs, threshold_q)
        if r is not None:
            vals.append(r["p"])
    return np.asarray(vals, dtype=np.float64), mode


def bca_interval(theta_hat, boot, jack, level=GPD_CI_LEVEL):
    """§P7-7(a). The BCa interval, plus the percentile endpoints as a REFERENCE only.

    Returns a dict:
      lo, hi          the interval the estimator reports. BCa, unless `fallback`.
      percentile      [lo, hi] of the plain percentile interval -- diagnosis and the
                      audit's reference column. NEVER an alternative the caller may
                      select; see the module docstring.
      z0, a           the two corrections, so a reader can recompute the endpoints
      alpha_lo/_hi    the adjusted probability levels actually read off the ECDF
      fallback        True iff the BCa arithmetic was degenerate and the percentile
                      endpoints were used instead
      notes           every guard that fired, in words

    The degenerate cases are enumerated rather than smoothed over:
      * fewer than 20 valid bootstrap replicates -> no ECDF worth inverting
      * every replicate on one side of the point estimate -> z0 = +-inf, so the
        fraction is clamped to +-1/(2B) (the standard continuity guard) and flagged
      * a jackknife with no spread -> sum d^2 = 0 -> a := 0, which degrades BCa to
        BC, the correct limit rather than a division by zero
      * 1 - a(z0 + z) <= 0 -> the BCa map is undefined at this level; fall back.
        This is a GUARD, not a routine outcome: |sum d^3| <= (sum d^2)^{3/2} for any
        real vector, so |a| <= 1/6 always, and the denominator can only go
        non-positive if |z0 + z| > 6 -- which the z0 clamp keeps out of reach at any
        B this estimator will ever run at.
      * alpha_lo >= alpha_hi -> the map inverted the interval; fall back
    """
    boot = np.asarray(boot, dtype=np.float64).ravel()
    boot = boot[np.isfinite(boot)]
    B = int(boot.size)
    lo_pct = (1.0 - float(level)) / 2.0 * 100.0
    out = {"lo": None, "hi": None, "percentile": None, "z0": None, "a": None,
           "alpha_lo": None, "alpha_hi": None, "n_boot": B, "n_jack": 0,
           "jackknife": None, "fallback": False, "level": float(level),
           "notes": [], "method": "BCa (bias-corrected and accelerated) bootstrap"}
    if B < 20:
        out["notes"].append(f"only {B} valid bootstrap replicates (< 20); no "
                            f"interval can be formed")
        return out
    perc = [float(np.percentile(boot, lo_pct)),
            float(np.percentile(boot, 100.0 - lo_pct))]
    out["percentile"] = perc

    # ---- bias correction z0. Ties with the point estimate are split, so that a
    # discrete pile-up at theta_hat does not push z0 down by its whole mass.
    below = float(np.count_nonzero(boot < float(theta_hat)))
    ties = float(np.count_nonzero(boot == float(theta_hat)))
    frac = (below + 0.5 * ties) / B
    eps = 0.5 / B
    if frac <= 0.0 or frac >= 1.0:
        out["notes"].append(
            f"every bootstrap replicate lies on one side of the point estimate "
            f"(fraction below = {frac:.3g}); z0 would be infinite, so the fraction "
            f"is clamped to {eps:.3g} from the boundary")
        frac = min(max(frac, eps), 1.0 - eps)
    z0 = float(stats.norm.ppf(frac))
    if not np.isfinite(z0):                                  # pragma: no cover
        out["notes"].append("z0 is not finite; falling back to percentile")
        out["fallback"] = True
        out["lo"], out["hi"] = perc
        return out
    out["z0"] = z0

    # ---- acceleration a, from the jackknife of the estimand
    jack = np.asarray(jack, dtype=np.float64).ravel()
    jack = jack[np.isfinite(jack)]
    out["n_jack"] = int(jack.size)
    a = 0.0
    if jack.size < 3:
        out["notes"].append(f"jackknife produced only {jack.size} valid replicates "
                            f"(< 3); acceleration set to 0, which degrades BCa to "
                            f"the bias-corrected interval")
    else:
        d = float(jack.mean()) - jack
        s2 = float(np.sum(d * d))
        s3 = float(np.sum(d ** 3))
        if s2 <= 0.0:
            out["notes"].append("jackknife replicates are identical (sum d^2 = 0); "
                                "acceleration set to 0")
        else:
            a = s3 / (6.0 * s2 ** 1.5)
    if not np.isfinite(a):                                   # pragma: no cover
        out["notes"].append("acceleration is not finite; set to 0")
        a = 0.0
    out["a"] = float(a)

    # ---- the adjusted levels
    z_lo = float(stats.norm.ppf((1.0 - float(level)) / 2.0))
    z_hi = float(stats.norm.ppf(1.0 - (1.0 - float(level)) / 2.0))

    def _adj(z):
        den = 1.0 - a * (z0 + z)
        if not np.isfinite(den) or den <= 0.0:
            return None
        v = float(stats.norm.cdf(z0 + (z0 + z) / den))
        return v if np.isfinite(v) else None

    al, ah = _adj(z_lo), _adj(z_hi)
    if al is None or ah is None or not (al < ah):
        out["notes"].append(
            f"the BCa map is degenerate at z0 = {z0:.4g}, a = {a:.4g} "
            f"(alpha_lo = {al}, alpha_hi = {ah}); the percentile endpoints are used "
            f"and this row is FLAGGED as a fallback, not reported as BCa")
        out["fallback"] = True
        out["lo"], out["hi"] = perc
        return out
    if al < eps or ah > 1.0 - eps:
        out["notes"].append(
            f"adjusted levels ({al:.4g}, {ah:.4g}) run past the outermost bootstrap "
            f"order statistic; clamped to [{eps:.3g}, {1 - eps:.3g}] -- the ECDF "
            f"holds no information beyond its own extremes")
        al = max(al, eps)
        ah = min(ah, 1.0 - eps)
    out["alpha_lo"], out["alpha_hi"] = float(al), float(ah)
    out["lo"] = float(np.percentile(boot, 100.0 * al))
    out["hi"] = float(np.percentile(boot, 100.0 * ah))
    return out


# --------------------------------------------------------------- the estimator --
def gpd_floor(n_max, decades=GPD_DECADES):
    """RULE 3. The hard floor: one decade below the Monte Carlo floor 1/(N_max+1)."""
    return float(10.0 ** (-decades) / (float(n_max) + 1.0))


def brute_force_n(p_gpd, factor=GPD_CONFIRM_FACTOR):
    """RULE 6. Surrogate count for the targeted single-test confirmation: N >= 10/p."""
    p = float(p_gpd)
    if not (p > 0.0):
        raise ValueError("brute_force_n: p must be positive")
    return int(math.ceil(factor / p))


def gpd_tail_p(S, s_obs, n_max, null_type, rng, threshold_q=GPD_THRESHOLD_Q,
               n_boot=GPD_N_BOOT, level=GPD_CI_LEVEL, ad_boot=GPD_N_AD_BOOT,
               min_exceedances=GPD_MIN_EXCEEDANCES, decades=GPD_DECADES):
    """The whole of §P6-2(1)-(4) in one call. Returns a dict; never raises on a bad fit.

    The ONE thing it does raise on is a forbidden null (rule 4), because that is a
    programming error at the call site rather than a data outcome.

    Returned keys a downstream consumer may read:
      p_method   one of P_METHODS
      p_report   the number that enters BH -- CI-upper for GPD rows, and for
                 UNRESOLVED rows the censored MC floor itself (never the point
                 estimate, never the extrapolation)
      p_ci_upper / p_ci_lower / p_point   the fit, for the report table
      reason     why a non-GPD label was assigned, in words
    """
    null_key = assert_null_permitted(null_type)
    S = np.asarray(S, dtype=np.float64).ravel()
    S = S[np.isfinite(S)]
    s_obs = float(s_obs)
    n = S.size
    floor = gpd_floor(n_max, decades)
    mc_floor = 1.0 / (float(n_max) + 1.0)
    out = {"p_method": P_UNRESOLVED, "p_report": mc_floor, "p_point": None,
           "p_ci_lower": None, "p_ci_upper": None, "null_type": null_key,
           "n_surrogates": int(n), "threshold_q": float(threshold_q),
           "floor": floor, "mc_floor": mc_floor, "reason": None,
           "ad": None, "xi_stability": None, "xi": None, "beta": None, "u": None}

    if n < 1:
        out["reason"] = "no surrogates"
        return out
    u = float(np.quantile(S, threshold_q))
    out["u"] = u
    exc = S[S > u] - u
    out["n_exceed"] = int(exc.size)
    if exc.size < min_exceedances:
        out["reason"] = (f"only {exc.size} exceedances over the declared top "
                         f"{100 * (1 - threshold_q):.0f}% threshold "
                         f"(< {min_exceedances}); an MLE tail fit here is not an "
                         f"estimate")
        return out
    if s_obs <= u:
        out["reason"] = ("observed statistic is below the tail threshold; the "
                         "Monte Carlo resolves it directly and extrapolation is "
                         "neither needed nor licensed")
        return out

    try:
        xi, beta = fit_gpd(exc)
    except Exception as e:                                   # pragma: no cover
        out["reason"] = f"GPD MLE failed: {e}"
        return out
    out["xi"], out["beta"] = float(xi), float(beta)

    # RULE 1(a)
    ad = ad_gate(exc, xi, beta, rng, n_boot=ad_boot)
    out["ad"] = ad
    # RULE 1(b)
    xs = xi_stability(S, rng, n_boot=n_boot, level=level,
                      min_exceedances=min_exceedances)
    out["xi_stability"] = xs
    if not ad["pass"]:
        out["reason"] = (f"§P6-2(1a) Anderson-Darling gate FAILED: p_AD = "
                         f"{ad['p_ad']:.4f} < {GPD_AD_ALPHA}. The exceedances are "
                         f"not GPD, so no extrapolation.")
        return out
    if not xs["pass"]:
        out["reason"] = ("§P6-2(1b) xi-stability gate FAILED: the shape parameter "
                         "does not agree across the top {5,10,20}% thresholds "
                         "within its own bootstrap CI, so the tail is not GPD over "
                         "the range being used.")
        return out

    point = _tail_p_from_sample(S, s_obs, threshold_q)
    if point is None:                                        # pragma: no cover
        out["reason"] = "the point extrapolation could not be formed"
        return out
    p_point = point["p"]
    out["p_point"] = p_point
    out["zeta"] = point["zeta"]                   # P(X > u), the exceedance rate

    # RULE 2: the 95% bootstrap CI, and the UPPER end of it is the answer.
    #
    # THE BOOTSTRAP RESAMPLES THE WHOLE SURROGATE SET, NOT JUST THE EXCEEDANCES.
    # Resampling only the exceedances holds BOTH the threshold u and the exceedance
    # rate zeta = n_exc/n fixed at their sample values, so the interval prices only
    # the GPD parameter uncertainty and is too narrow by the two other sources. The
    # §P6-2(7) calibration MEASURED that: with the exceedances-only bootstrap the CI
    # covered the brute-force p in 22 of 26 cases (84.6%). Resampling S itself and
    # re-deriving u, zeta and (xi, beta) on every replicate propagates all three,
    # which is what a bootstrap of this estimator has to do. That was a CORRECTNESS
    # repair and it stands; the JACKKNIFE below is held to the identical standard by
    # sharing `_tail_p_from_sample` with the bootstrap.
    ps = _bootstrap_ps(S, s_obs, rng, threshold_q=threshold_q, n_boot=n_boot)
    if ps.size < 20:
        out["reason"] = "bootstrap CI could not be formed (too few valid refits)"
        return out
    out["n_boot_ok"] = int(ps.size)

    # §P7-7(a): BCa. Bias correction from the bootstrap distribution's position
    # relative to the point estimate; acceleration from the jackknife of the SAME
    # estimand, the extrapolated p.
    jack, jmode = _jackknife_ps(S, s_obs, threshold_q=threshold_q)
    bca = bca_interval(p_point, ps, jack, level=level)
    bca["jackknife"] = jmode
    out["bca"] = bca
    if bca["lo"] is None or bca["hi"] is None:               # pragma: no cover
        out["reason"] = ("bootstrap CI could not be formed ("
                         + "; ".join(bca["notes"]) + ")")
        return out
    p_lo, p_hi = bca["lo"], bca["hi"]
    out["p_ci_lower"], out["p_ci_upper"] = p_lo, p_hi
    # REFERENCE ONLY (§P7-7(a)): the percentile endpoints are kept so the audit can
    # print what the superseded estimator would have said. Nothing reads them.
    out["p_ci_lower_percentile"] = bca["percentile"][0]
    out["p_ci_upper_percentile"] = bca["percentile"][1]

    # RULE 3: one decade and no further. Note the test is on the number we would
    # REPORT (the CI upper end), not on the point estimate -- reporting a CI-upper
    # above the floor while the point estimate is below it is exactly the
    # conservative behaviour rule 2 asks for.
    if p_hi < floor:
        out["reason"] = (f"§P6-2(3) one-decade cap: the 95% CI upper end "
                         f"{p_hi:.3g} is below the floor {floor:.3g} = "
                         f"0.1/(N_max+1). Beyond one decade the shape parameter "
                         f"dominates the answer and the CI stops meaning anything.")
        return out
    out["p_method"] = P_GPD
    out["p_report"] = p_hi                       # RULE 2, enforced here
    out["reason"] = None
    return out


def compact(g):
    """The subset of a `gpd_tail_p` result that is small enough to store per row.

    The bootstrap replicate vectors and the AD replicates are NOT stored: they are
    per-test arrays of hundreds of floats and they would go into the checkpoint JSON
    and stubs.json, where nobody reads them and everybody pays for them.
    """
    if g is None:
        return None
    keep = ("p_method", "p_point", "p_ci_lower", "p_ci_upper", "u", "xi", "beta",
            "n_exceed", "n_surrogates", "threshold_q", "floor", "mc_floor",
            "null_type", "reason",
            # reference only, never read by a consumer (§P7-7(a))
            "p_ci_lower_percentile", "p_ci_upper_percentile")
    out = {k: g.get(k) for k in keep}
    b = g.get("bca") or {}
    out["ci_method"] = "BCa" if (b and not b.get("fallback")) else (
        "percentile-fallback" if b else None)
    out["bca_z0"] = b.get("z0")
    out["bca_a"] = b.get("a")
    out["bca_alpha_lo"] = b.get("alpha_lo")
    out["bca_alpha_hi"] = b.get("alpha_hi")
    out["bca_jackknife"] = b.get("jackknife")
    out["bca_n_jack"] = b.get("n_jack")
    out["bca_fallback"] = b.get("fallback")
    out["bca_notes"] = b.get("notes")
    ad = g.get("ad") or {}
    out["p_ad"] = ad.get("p_ad")
    out["ad_a2"] = ad.get("a2")
    out["ad_pass"] = ad.get("pass")
    xs = g.get("xi_stability") or {}
    out["xi_stability_pass"] = xs.get("pass")
    out["xi_by_threshold"] = [
        {"q": r.get("q"), "xi": r.get("xi"), "ci": r.get("ci"),
         "n_exceed": r.get("n_exceed")} for r in xs.get("rows", [])]
    return out


def price_row(p_mc, mc_n_max, s_obs, S_surr, null_type, rng, other_p=0.0,
              other_floor=0.0, enabled=False, **kw):
    """Label one test's p and decide what number enters BH. Rules 1, 2, 3, 5, 6.

    Arguments are deliberately the raw ingredients rather than a result row, so that
    this is testable without building a pipeline:

      p_mc      the Monte Carlo p from the INDEPENDENT-surrogate null
      mc_n_max  the DECLARED surrogate cap (not the number of draws actually made)
      s_obs     the observed statistic
      S_surr    that null's surrogate statistics, or None if they were not retained
      other_p   any other component of the reported p_raw that GPD may not touch --
                in this engine the exhaustive circular-shift enumeration, on which
                extrapolation is FORBIDDEN (rule 4). The reported p stays the
                conservative max of the components.
      enabled   GPD extrapolation is OPT-IN (`--gpd`). Off, the label is still
                emitted on every row (rule 5) and describes the censoring state; on,
                the full §P6-2 path runs and a failed gate makes the test
                INELIGIBLE FOR BH REJECTION (rule 1).

    Returns a dict to merge into the row: p_method, p_bh, bh_eligible, gpd.
    """
    mc_floor = 1.0 / (float(mc_n_max) + 1.0)
    p_mc = float(p_mc)
    other_p = float(other_p or 0.0)
    p_report = max(p_mc, other_p)
    censored = (p_mc <= mc_floor + 1e-15) and (p_mc >= other_p)
    if not censored:
        # The OTHER component can be binding and be sitting on ITS OWN floor -- the
        # exhaustive circular-shift enumeration, whose floor is set by the record
        # length. That p is censored too, and rule 5 says say so; but §P6-2(4)
        # FORBIDS extrapolating off it, so the honest label is UNRESOLVED with the
        # test still ELIGIBLE (nothing about it failed a gate -- there is simply no
        # licensed instrument that could resolve it further).
        if other_p > p_mc and other_p <= float(other_floor) + 1e-15:
            return {"p_method": P_UNRESOLVED, "p_bh": p_report,
                    "bh_eligible": True, "gpd": None,
                    "p_method_reason": (
                        f"the binding component is the circular-shift enumeration "
                        f"at its own floor {other_floor:.3g}, a property of the "
                        f"RECORD LENGTH; §P6-2(4) forbids GPD extrapolation on that "
                        f"null, so the p stays censored")}
        return {"p_method": P_MC_RESOLVED, "p_bh": p_report, "bh_eligible": True,
                "gpd": None,
                "p_method_reason": "the Monte Carlo resolved this p directly"}
    if not enabled:
        return {"p_method": P_UNRESOLVED, "p_bh": p_report, "bh_eligible": True,
                "gpd": None,
                "p_method_reason": (f"censored at the Monte Carlo floor "
                                    f"{mc_floor:.3g}; GPD extrapolation is OFF "
                                    f"(--gpd not passed), so the p is the floor and "
                                    f"the true p is somewhere at or below it")}
    if S_surr is None:
        return {"p_method": P_UNRESOLVED, "p_bh": p_report, "bh_eligible": False,
                "gpd": None,
                "p_method_reason": ("censored, and this null's surrogate statistics "
                                    "were not retained, so no tail fit is possible")}
    g = gpd_tail_p(S_surr, s_obs, mc_n_max, null_type, rng, **kw)
    if g["p_method"] == P_GPD:
        return {"p_method": P_GPD, "p_bh": max(g["p_ci_upper"], other_p),
                "bh_eligible": True, "gpd": compact(g),
                "p_method_reason": ("extrapolated; the number entering BH is the "
                                    "UPPER end of the 95% bootstrap CI (§P6-2(2)), "
                                    "never the point estimate")}
    # RULE 1: a failed gate means NO extrapolation and NO BH rejection.
    return {"p_method": P_UNRESOLVED, "p_bh": p_report, "bh_eligible": False,
            "gpd": compact(g),
            "p_method_reason": g.get("reason") or "extrapolation not licensed"}


# ------------------------------------------- rule 7: calibration before default -
# The licensing artifact is the BCa one. `audit_gpd.json` -- the phase-2a percentile
# run -- is NOT it and can never become it: §P7-7(c) scored that attempt UNRESOLVED
# at its design's resolution, and §P7-7(b) requires the licensing set to be >= 125
# FRESH comparisons, which the 26-comparison legacy set is not. Pointing the default
# here is what makes "the legacy subset licenses nothing" structural.
CALIBRATION_PATH = os.path.join("engine", "out", "audit_gpd_bca.json")
LEGACY_CALIBRATION_PATH = os.path.join("engine", "out", "audit_gpd.json")


class NotCalibratedError(RuntimeError):
    """Raised when extrapolation is requested and §P6-2(7) has not been cleared."""


def load_calibration(path=None):
    """The §P6-2(7) verdict artifact, or None if the audit has never been run."""
    p = path or CALIBRATION_PATH
    if not os.path.exists(p):
        return None
    with open(p, "rb") as fh:
        raw = fh.read()
    doc = json.loads(raw.decode("utf-8"))
    return {"path": str(p).replace("\\", "/"),
            "sha256": hashlib.sha256(raw).hexdigest(),
            "pass": bool(doc.get("pass")),
            "coverage": doc.get("coverage"),
            "coverage_bar": doc.get("coverage_bar"),
            "worst_understatement": doc.get("worst_understatement"),
            "understatement_bar": doc.get("understatement_bar"),
            "n_comparisons": doc.get("n_comparisons"),
            "n_fitted": doc.get("n_fitted"),
            "n_scored": doc.get("n_scored"),
            "ci_method": doc.get("ci_method"),
            "set": doc.get("set")}


def assert_calibrated(path=None):
    """RULE 7, made binding rather than decorative.

    "An instrument that has not demonstrated recovery in the range it reports in may
    not report there." So `--gpd` does not merely PREFER a calibration -- it requires
    one, and a REJECT verdict is as disqualifying as no verdict at all. The
    calibration summary is returned so it can be written into the config hash: a
    re-calibration is a new licence and therefore a new experiment.
    """
    cal = load_calibration(path)
    p = path or CALIBRATION_PATH
    if cal is None:
        raise NotCalibratedError(
            f"§P6-2(7): GPD tail extrapolation was requested but no calibration "
            f"artifact exists at {p}. An instrument that has not demonstrated "
            f"recovery in the range it reports in may not report there. Run "
            f"`python -u -m engine.audit_gpd --jobs N` first. Refusing to run.")
    if not cal["pass"]:
        raise NotCalibratedError(
            f"§P6-2(7): the calibration at {p} is a REJECT -- coverage "
            f"{cal['coverage']} against a bar of {cal['coverage_bar']}, worst "
            f"understatement {cal['worst_understatement']}x against a bar of "
            f"{cal['understatement_bar']}x. The estimator may not report until it "
            f"clears. Refusing to run.")
    # §P7-9(1): the count is SCORED comparisons. `n_comparisons` is how many were
    # RUN, and the two differ by every row the gates threw out -- rows that entered
    # neither bar and so cannot have contributed power to either. An artifact that
    # does not state its scored count cannot be checked against the design and is
    # refused rather than given the benefit of the doubt.
    n_scored = cal.get("n_scored")
    if n_scored is None:
        raise NotCalibratedError(
            f"§P7-9(1): the calibration at {p} does not report `n_scored`, so the "
            f"size of the set that actually entered the two bars is unknown. The "
            f"design is stated in SCORED comparisons and cannot be checked against "
            f"a run count. Re-run `python -u -m engine.audit_gpd`. Refusing to run.")
    if int(n_scored) < GPD_MIN_CALIBRATION_COMPARISONS:
        raise NotCalibratedError(
            f"§P7-9(1): the calibration at {p} SCORED only {n_scored} comparisons "
            f"(of {cal.get('n_comparisons')} run), below the required decisive set "
            f"of {GPD_MIN_CALIBRATION_COMPARISONS}. Gated-out rows entered neither "
            f"bar, so counting them toward the design's power would count an "
            f"absence as evidence. A set this size cannot resolve its own 90% "
            f"coverage bar -- attempt 1's 22/26 ran an exact 95% interval of "
            f"[0.675, 0.946], which CONTAINS the bar. A pass on an underpowered "
            f"design is not a pass. Refusing to run.")
    if cal.get("ci_method") not in (None, "BCa"):
        raise NotCalibratedError(
            f"§P7-7(a): the calibration at {p} was produced with interval method "
            f"'{cal['ci_method']}', not the declared BCa estimator. A licence is "
            f"issued to the estimator that was calibrated and to no other. "
            f"Refusing to run.")
    return cal


# ------------------------------------------------------------- labels and census
def census(rows, key="p_method"):
    """RULE 5. The three-number census that goes in the banner. Always three keys."""
    out = {m: 0 for m in P_METHODS}
    for r in rows:
        m = r.get(key) if isinstance(r, dict) else r
        if m in out:
            out[m] += 1
        else:
            out.setdefault("UNLABELLED", 0)
            out["UNLABELLED"] += 1
    return out


def census_line(counts, total=None):
    """One printable line. Named numbers, so a reader never has to infer the kind."""
    n = sum(counts.get(m, 0) for m in P_METHODS) if total is None else int(total)
    return (f"p-method census (§P6-2(5)): {counts.get(P_MC_RESOLVED, 0)} "
            f"{P_MC_RESOLVED} / {counts.get(P_GPD, 0)} {P_GPD} / "
            f"{counts.get(P_UNRESOLVED, 0)} {P_UNRESOLVED}  (of {n} declared tests)")
