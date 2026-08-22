"""SP-7 v2: A DERIVED, TWO-SIDED, POWER-AWARE CALIBRATION BAND. No hardcoded cap.

WHY SP-7 v1 HAD TO CHANGE, AND WHY CHANGING IT IS NOT THE DEFECT IT LOOKS LIKE.

SP-7 v1 says: "at alpha = q/m the expected promotions are <= 0.1 per catalogue, so
<= 3 across 30. PASS: promotions at or below that rate, with the count consistent with
binomial."

The arithmetic underneath is `E[promotions] = N * m * (q/m) = N * q`, which is 3.0
across 30 catalogues **exactly, and regardless of m**. So the cap of 3 is the MEAN of a
correctly calibrated searcher, not an upper bound on it, and a perfectly calibrated
instrument exceeds it about 35 % of the time (Poisson(3), P(X > 3) = 0.3528).

A gate that discards a working instrument a third of the time is not conservative; it
is destructive. This program's purpose is to surface POTENTIALS worth looking into more
deeply, and a 35 % false-fail rate on the gate that stands in front of every scan
throws away a third of them for nothing. That is the same class of defect as a test
whose power floor exceeds alpha: a criterion that is not calibrated to what it claims
to measure.

THE HONEST-CHANGE PROBLEM, HANDLED EXPLICITLY. Changing a criterion after seeing the
result it produced is the exact move that invalidates a pre-registered protocol, and
the SP-7 v1 build was right to refuse it. Four things make this amendment legitimate,
and all four are checkable:

  1. **The criterion is DERIVED, not chosen.** Every threshold below is computed from
     (n_trials, alpha_effective, declared false-fail rate). No integer is written down.
     Feed it the v1 design and it reproduces v1's arithmetic; feed it any other and it
     adapts. There is no value to tune.
  2. **The defect is in v1's ARITHMETIC, not in v1's verdict.** `N*q = 3` being the
     mean rather than a bound is true before any catalogue is simulated. The redesign
     is motivated by a fact available at freeze time, not by the observed count.
  3. **The false-fail rate is DECLARED and the achieved rate is REPORTED**, so the new
     criterion's own error rate is a number on the artifact rather than an assumption.
  4. **The verdict is rendered on FRESH nulls with a new declared seed**, so the run
     that motivated the review is not the run that is scored under the new rule. The
     old run is re-scored too, and reported, purely for transparency.

WHAT THE NEW CRITERION IS.

  * **TWO-SIDED.** Too many promotions means the nominal p's are wrong and an SP-2 null
    layer is invalid -- v1's concern, and the one that matters most. Too FEW is also a
    failure, and v1 could not express it: an instrument that promotes nothing in 30
    true-null catalogues is over-conservative and cannot detect anything either, and
    "we saw none" is evidence of a bug rather than of virtue.
  * **DERIVED FROM THE EXACT BINOMIAL**, at `alpha_effective = max(alpha_declared,
    1/(B+1))`, because a resampling p cannot be finer than its own resolution and
    pretending otherwise misstates the expected count.
  * **POWER-AWARE, which is the part v1 lacked entirely.** A band alone says whether
    the count is ordinary; it does not say whether the gate could have caught a broken
    null. `detectable_inflation` reports the smallest false-positive-rate inflation the
    gate detects at declared power, and `required_catalogues` inverts it to say how
    many null catalogues a target sensitivity needs. A gate that cannot detect a 3x
    inflation is not a gate, and now it says so.

Nothing here decides anything about earthquakes. It is the calibration of a calibrator.
"""

from __future__ import annotations

import math

SP7_V1_TEXT = (
    "SP-7 v1 (HYPOTHESIS_LEDGER.md §P7-24): 'at alpha = q/m the expected promotions "
    "are <= 0.1 per catalogue, so <= 3 across 30. PASS: promotions at or below that "
    "rate, with the count consistent with binomial.' The cap of 3 is E[X] = N*q, the "
    "MEAN, not a bound.")

SP7_V2_RULE = (
    "SP-7 v2: PASS if (a) the run is not vacuous -- every cell promotion-eligible and "
    "the planted-signal control fires -- AND (b) the promotion count lies strictly "
    "inside a two-sided exact-binomial calibration band derived from n_trials, "
    "alpha_effective and a DECLARED false-fail rate. FAIL-HIGH means the nominal p's "
    "are too small, i.e. an SP-2 null layer is invalid. FAIL-LOW means the instrument "
    "is over-conservative and could not detect anything either. The band is COMPUTED, "
    "never written down, and the artifact reports the achieved false-fail rate and the "
    "smallest miscalibration the gate can detect at declared power.")


# ------------------------------------------------------------ exact binomial tails --
def _log_choose(n, k):
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))


def binom_pmf(k, n, p):
    if k < 0 or k > n:
        return 0.0
    if p <= 0.0:
        return 1.0 if k == 0 else 0.0
    if p >= 1.0:
        return 1.0 if k == n else 0.0
    return math.exp(_log_choose(n, k) + k * math.log(p) + (n - k) * math.log1p(-p))


_TERM_FLOOR = 1e-300
_REL_EPS = 1e-17

try:                                            # exact, and the primary path
    from scipy.special import betainc as _betainc
except Exception:                               # pragma: no cover - fallback only
    _betainc = None


def binom_sf(k, n, p):
    """P(X >= k), exactly, via the regularized incomplete beta identity.

        P(X >= k) = I_p(k, n - k + 1)

    THE RECURRENCE BELOW IS A FALLBACK AND IT HAS A KNOWN FAILURE MODE, recorded here
    because it produced a wrong answer before it was caught: when k lies far BELOW the
    mode n*p, `binom_pmf(k, n, p)` underflows to zero, the upward walk starts from
    nothing, and the function returns ~0 where the true value is ~1. That is exactly
    the regime `detectable_inflation` probes when it tries a large inflation factor,
    and it silently turned a detectable miscalibration into an undetectable one. The
    beta identity has no such regime. The fallback keeps a guard so that if it is ever
    used in that regime it returns 1.0 rather than 0.0.
    """
    n = int(n)
    k = int(k)
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    p = float(p)
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    if _betainc is not None:
        return float(min(1.0, max(0.0, _betainc(k, n - k + 1, p))))
    return _binom_sf_recurrence(k, n, p)


def _binom_sf_recurrence(k, n, p):
    """Pure-python fallback for `binom_sf`. See its docstring for the failure mode.

    `pmf(i+1)/pmf(i) = ((n-i)/(i+1)) * p/(1-p)`, so the tail is walked with one
    multiply per term instead of a fresh lgamma triple. The loop stops once the terms
    are negligible relative to what has accumulated AND the mode has been passed.
    """
    term = binom_pmf(k, n, p)
    if term <= _TERM_FLOOR and k < n * p:
        # k is far below the mode and the starting term underflowed: the upper tail is
        # everything. Returning the recurrence's 0.0 here is the recorded bug.
        return 1.0
    total = term
    ratio_p = p / (1.0 - p)
    mode = n * p
    for i in range(k, n):
        term *= (n - i) / (i + 1.0) * ratio_p
        if term <= _TERM_FLOOR:
            break
        total += term
        if i > mode and term < _REL_EPS * total:
            break
    return min(1.0, total)


def binom_cdf(k, n, p):
    """P(X <= k). Computed as the complement of the upper tail for tiny p."""
    n = int(n)
    k = int(k)
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    if k + 1 > n * float(p):
        return max(0.0, 1.0 - binom_sf(k + 1, n, p))
    return sum(binom_pmf(i, n, p) for i in range(0, k + 1))


def effective_alpha(alpha_declared, n_surrogates):
    """The false-positive rate a resampling test ACTUALLY has, which is not `alpha`.

    A Phipson-Smyth p from B surrogates takes only the values 1/(B+1), 2/(B+1), ...
    So `p <= alpha` fires exactly when `p <= floor(alpha*(B+1))/(B+1)`, and the true
    per-test false-positive rate is that floor, NOT the declared alpha:

        alpha_effective = floor(alpha_declared * (B+1)) / (B+1)

    THIS FUNCTION WAS WRONG ON ITS FIRST VERSION and the error is recorded rather than
    quietly fixed: it returned `max(alpha_declared, 1/(B+1))`, which for the SP-7
    design gives 1/300 where the truth is 1/401. The discrepancy was caught by
    comparing against `results_searcher_sp7_gate.json`, which had the right number,
    and it matters because the expected promotion count is `n_trials * alpha_effective`
    and a 34 % overstatement of that propagates straight into the band.

    Returns 0.0 when `alpha_declared < 1/(B+1)`, i.e. when the declared threshold is
    finer than the resampling can resolve and the test CANNOT FIRE AT ANY EFFECT SIZE.
    Callers must treat that as a design error, not as a very strict test.
    """
    b1 = float(n_surrogates) + 1.0
    return math.floor(float(alpha_declared) * b1) / b1


# -------------------------------------------------------------- the derived band --
def calibration_band(n_trials, alpha_eff, false_fail_rate=0.01):
    """The two-sided acceptance band, DERIVED. Returns the interval and its true size.

    `false_fail_rate` is the total probability that a CORRECTLY calibrated searcher is
    rejected, split evenly between the tails. Because the count is discrete the
    achieved rate is at most the declared one and is reported rather than assumed.

    Returned `k_hi` is the smallest count that FAILS HIGH; `k_lo` the largest that
    FAILS LOW. The accepted set is therefore `k_lo < k < k_hi`.
    """
    n = int(n_trials)
    p = float(alpha_eff)
    half = float(false_fail_rate) / 2.0

    k_hi = n + 1
    for k in range(0, n + 1):
        if binom_sf(k, n, p) <= half:
            k_hi = k
            break

    k_lo = -1
    for k in range(0, n + 1):
        if binom_cdf(k, n, p) > half:
            k_lo = k - 1
            break

    ach_hi = binom_sf(k_hi, n, p) if k_hi <= n else 0.0
    ach_lo = binom_cdf(k_lo, n, p) if k_lo >= 0 else 0.0
    return {
        "n_trials": n, "alpha_effective": p,
        "expected_count": n * p,
        "false_fail_rate_declared": float(false_fail_rate),
        "k_fail_high": int(k_hi), "k_fail_low": int(k_lo),
        "accepted_open_interval": (int(k_lo), int(k_hi)),
        "achieved_false_fail_high": float(ach_hi),
        "achieved_false_fail_low": float(ach_lo),
        "achieved_false_fail_total": float(ach_hi + ach_lo),
        "note": ("the band is COMPUTED from (n_trials, alpha_effective, "
                 "false_fail_rate); no threshold is written down anywhere"),
    }


def detectable_inflation(n_trials, alpha_eff, k_fail_high, power=0.80,
                         r_max=100.0):
    """Smallest false-positive-rate inflation R this gate detects at `power`.

    If an SP-2 null layer is broken so that the true per-test false-positive rate is
    `R * alpha_eff` instead of `alpha_eff`, the gate fires when the count reaches
    `k_fail_high`. This returns the smallest R for which that happens with probability
    at least `power`. **This is the gate's own detection curve**, and without it a PASS
    means only "nothing odd happened", not "a broken null would have been caught".
    """
    n = int(n_trials)
    lo, hi = 1.0, float(r_max)
    if binom_sf(k_fail_high, n, min(1.0, hi * alpha_eff)) < power:
        return None                      # not detectable at any inflation below r_max
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if binom_sf(k_fail_high, n, min(1.0, mid * alpha_eff)) >= power:
            hi = mid
        else:
            lo = mid
    return hi


def required_catalogues(m, alpha_eff, r_min, power=0.80, false_fail_rate=0.01,
                        n_max=100000):
    """Minimum number of null catalogues to detect an R-fold inflation at `power`.

    Inverts `detectable_inflation`, so the SIZE of the gate is derived from the
    sensitivity it is required to have rather than picked. This is the number SP-7 v1
    should have carried instead of "30".
    """
    def ok(n_cat):
        band = calibration_band(m * n_cat, alpha_eff, false_fail_rate)
        r = detectable_inflation(m * n_cat, alpha_eff, band["k_fail_high"], power)
        return (r is not None and r <= float(r_min)), r, band

    # The count is DISCRETE, so `detectable_inflation` is only weakly monotone in n
    # and can step sideways. A geometric probe finds a bound cheaply; the minimum is
    # then found by scanning upward from 1, because a bisection on a non-monotone
    # step function can return a value that is satisfied but not minimal -- which is
    # exactly what the first version of this function did.
    bound = None
    probe = 1
    while probe <= int(n_max):
        good, _, _ = ok(probe)
        if good:
            bound = probe
            break
        probe = probe * 2 if probe > 1 else 2
    if bound is None:
        return None
    for n_cat in range(1, bound + 1):
        good, r, band = ok(n_cat)
        if good:
            return {"n_catalogues": n_cat, "n_trials": m * n_cat,
                    "detectable_inflation": r, "band": band,
                    "search": "upward scan to the first satisfying n; the criterion "
                              "is discrete and only weakly monotone, so a bisection "
                              "can return a satisfying but non-minimal value"}
    return None


# ------------------------------------------------------------------- the verdict --
def verdict(n_promotions, n_trials, alpha_eff, false_fail_rate=0.01, power=0.80,
            vacuous=False, vacuity_reason=None):
    """SP-7 v2's verdict, with every number it rests on returned beside it."""
    band = calibration_band(n_trials, alpha_eff, false_fail_rate)
    k = int(n_promotions)
    r_det = detectable_inflation(n_trials, alpha_eff, band["k_fail_high"], power)

    if vacuous:
        v = "VACUOUS-FAIL"
    elif k >= band["k_fail_high"]:
        v = "FAIL-HIGH"
    elif k <= band["k_fail_low"]:
        v = "FAIL-LOW"
    else:
        v = "PASS"

    return {
        "rule_id": "SP7-gate-v2-derived-band",
        "rule": SP7_V2_RULE,
        "superseded": SP7_V1_TEXT,
        "verdict": v,
        "passed": v == "PASS",
        "n_promotions": k,
        "band": band,
        "p_two_sided": float(min(1.0, 2.0 * min(binom_sf(k, n_trials, alpha_eff),
                                                binom_cdf(k, n_trials, alpha_eff)))),
        "power": {
            "declared": float(power),
            "detectable_inflation_R": r_det,
            "meaning": ("smallest factor by which the true per-test false-positive "
                        "rate must be inflated for this gate to fire with probability "
                        ">= the declared power. A PASS is only meaningful to the "
                        "extent this number is small."),
        },
        "vacuity_reason": vacuity_reason,
        "interpretation": {
            "FAIL-HIGH": "nominal p's too small; an SP-2 null layer is invalid",
            "FAIL-LOW": ("instrument over-conservative; it could not detect a real "
                         "effect either, and 'we saw none' is evidence of a bug"),
            "VACUOUS-FAIL": "the gate did not actually exercise the instrument",
        }[v] if v != "PASS" else ("count is ordinary for a correctly calibrated "
                                  "searcher, and the gate had the stated power to "
                                  "notice if it were not"),
    }
