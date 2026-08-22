"""THE DWELL-TIME NULL LAYER, generalised. §P7-24 SP-2 as a module instead of a habit.

WHY THIS EXISTS, AND WHY IT IS NOT A REFACTOR.

§P7-24 SP-2 makes one null layer mandatory for the property class *level / waveform
phase*: "dwell-time-corrected null -- a level-threshold statistic against a
uniform-phase null is wrong by construction." D-1b built that layer for one statistic
(the below-neutral-falling quadrant) and D-1c built it again for eight more. Building it
a third time by hand is how the third one gets it subtly wrong.

More importantly, §P7-25(5) K-096 condition (ii) makes it a PROMOTION GATE for a class
that does not have one yet:

  > Every axis gets its waveform-matched time-uniform null, built by the same
  > machinery, before its number is looked at. A directional axis needs its own
  > dwell-time null just as badly as a level did. A bearing distribution from a
  > rotating principal axis is NOT uniform on the circle, and reading an azimuth
  > concentration against a uniform-circle null would be D-1's exact error in a new
  > coordinate. ... A property whose class has no null layer cannot promote.

This module is that layer. It is deliberately generic about the field: it takes SAMPLED
FEATURE ARRAYS and knows nothing about tides, so the same machinery serves a level, a
sign quadrant, a bearing, a rotation rate, or anything else sampled from a deterministic
waveform.

THE ONE IDEA. If a quantity is a deterministic function of time at a site, then the null
for "an event happened at a uniformly random time" is the quantity's own OCCUPANCY
DISTRIBUTION over time at that site. Not a uniform distribution over the quantity's
range, not a uniform distribution over phase, and for a bearing not a uniform
distribution over the circle. You do not model the occupancy; you sample it, from the
same array the observation was read from, so the observation and the null cannot diverge
through a second code path.

WHAT THE CALLER STILL OWES, because this module cannot check it:
  * The sampling pool must be GEOMETRY- AND HYPOTHESIS-INDEPENDENT. If the pool depends
    on the thing being tested, the null is contaminated. `time_uniform_pool` builds one
    from margins alone for exactly this reason.
  * The declared statistic list must be fixed BEFORE the numbers are seen, and its
    length is the `m` that `assert_power_floor` prices.
  * A null is only as good as the span it is drawn over. `span_sensitivity` exists so
    that dependence is measured rather than assumed.
"""

from __future__ import annotations

import math

import numpy as np

SP2_RULE = (
    "§P7-24 SP-2: property class 'level / waveform phase' carries the mandatory null "
    "layer 'dwell-time-corrected null -- a level-threshold statistic against a "
    "uniform-phase null is wrong by construction'. §P7-25(5) K-096(ii) extends the "
    "same requirement to DIRECTIONAL axes: a bearing distribution from a rotating "
    "principal axis is not uniform on the circle.")


# ----------------------------------------------------------------- power floor --
def assert_power_floor(m, n_replicates, alpha=0.05):
    """Refuse a family that cannot reach significance at any effect size.

    Smallest attainable p with N replicates is 1/(N+1); a family-wise correction over
    m statistics puts the floor at m/(N+1). A test that cannot fire is worse than no
    test, because it still prints a 'not significant' line that a reader counts as
    evidence of absence. Encoded as the FORMULA, so adding a statistic raises the
    requirement automatically instead of silently invalidating a fixed number.
    """
    floor = m / (n_replicates + 1.0)
    if floor > alpha:
        raise ValueError(
            "family-wise power floor m/(N+1) = %d/%d = %.4f exceeds alpha = %.3f; "
            "raise N to at least %d" % (m, n_replicates + 1, floor, alpha,
                                        math.ceil(m / alpha) - 1))
    return {"m": int(m), "n_replicates": int(n_replicates), "floor": floor,
            "alpha": alpha, "min_attainable_p": 1.0 / (n_replicates + 1.0)}


# ------------------------------------------------------------ the sampling pool --
def time_uniform_pool(n_samples, margin):
    """Indices eligible as null draws: the interior, by MARGIN ALONE.

    Deliberately depends on nothing but the grid length and the margin the derived
    features need. A pool that depended on the hypothesis under test would contaminate
    the null with the hypothesis, and the failure is silent.
    """
    margin = int(margin)
    if 2 * margin >= n_samples:
        raise ValueError("margin %d leaves no interior in %d samples"
                         % (margin, n_samples))
    return np.arange(margin, n_samples - margin, dtype=np.int64)


def draw_replicates(pool, n_sites, n_replicates, rng):
    """One uniformly random sample index per site per replicate: [n_sites, n_rep]."""
    idx = np.empty((int(n_sites), int(n_replicates)), dtype=np.int64)
    for s in range(int(n_sites)):
        idx[s] = pool[rng.integers(0, pool.size, int(n_replicates))]
    return idx


# ------------------------------------------------------- occupancy diagnostics --
def occupancy(values, bins):
    """The fraction of TIME a sampled quantity spends in each bin. The null itself.

    `bins` is a sequence of (name, boolean-mask-callable). Returned fractions are the
    dwell-time probabilities that replace whatever closed form a caller was about to
    assume.
    """
    v = np.asarray(values)
    ok = np.isfinite(v) if v.dtype.kind == "f" else np.ones(v.shape, dtype=bool)
    n = int(ok.sum())
    return {name: float(np.sum(ok & fn(v)) / max(n, 1)) for name, fn in bins}


def axial_occupancy(bearings_deg, n_bins=36):
    """How far a bearing's dwell distribution is from uniform on the circle.

    THE POINT OF THIS FUNCTION IS TO STOP SOMEONE ASSUMING IT IS UNIFORM. It returns
    the histogram, the axial mean resultant length R2 of the DOUBLED angle, and a
    chi-square against uniformity. R2 materially above zero means a uniform-circle null
    would have been wrong, and by how much.
    """
    b = np.mod(np.asarray(bearings_deg, dtype=np.float64), 180.0)
    b = b[np.isfinite(b)]
    if b.size == 0:
        raise ValueError("no finite bearings")
    hist, _ = np.histogram(b, bins=n_bins, range=(0.0, 180.0))
    exp = b.size / n_bins
    chi2 = float(np.sum((hist - exp) ** 2 / exp))
    two = 2.0 * b * np.pi / 180.0
    r2 = float(np.hypot(np.mean(np.cos(two)), np.mean(np.sin(two))))
    return {
        "n": int(b.size), "n_bins": int(n_bins),
        "histogram_fraction": (hist / b.size).tolist(),
        "axial_R2": r2,
        "chi2_vs_uniform": chi2, "chi2_df": n_bins - 1,
        "uniform_circle_null_would_be_wrong_by": r2,
        "note": ("R2 is the axial concentration of the NULL itself. If it is not ~0, "
                 "a uniform-circle null is wrong and §P7-25(5) K-096(ii) applies."),
        "rule": SP2_RULE,
    }


# ------------------------------------------------- max-statistic calibration --
def calibrate(observed, null_draws_stats, stat_names):
    """Max-statistic calibration of a declared battery against a null ensemble.

    `observed` maps stat name -> scalar. `null_draws_stats` maps stat name -> array of
    per-replicate values. Returns per-statistic z, the observed max |z|, and the exact
    p of that maximum against the null's own max |z| distribution.

    Dependence among the statistics is handled IMPLICITLY: correlated statistics
    inflate the null maximum exactly as they inflate the observed one, so the
    correlation is reproduced rather than modelled. That is the whole reason this is
    the right correction for a collinear lattice, and it is why §P7-25(5) disposition
    19 makes the simulated max-statistic PRIMARY and q/m only a conservative bracket.
    """
    z, cols = {}, []
    for k in stat_names:
        arr = np.asarray(null_draws_stats[k], dtype=np.float64)
        mu, sd = float(arr.mean()), float(arr.std(ddof=1))
        z[k] = {"observed": float(observed[k]), "null_mean": mu, "null_sd": sd,
                "z": (float(observed[k]) - mu) / sd if sd > 0 else float("nan")}
        cols.append((arr - mu) / (sd if sd > 0 else 1e-300))
    null_max = np.max(np.abs(np.stack(cols)), axis=0)
    obs_max = max(abs(z[k]["z"]) for k in stat_names)
    n = null_max.size
    # Phipson-Smyth: ties count AGAINST significance, and the +1 top and bottom
    # prevents a p of exactly zero, which finite resampling never justifies.
    p = (int(np.sum(null_max >= obs_max)) + 1) / (n + 1)
    return {
        "per_statistic": z,
        "max_abs_z_observed": obs_max,
        "max_abs_z_null_p95": float(np.quantile(null_max, 0.95)),
        "max_statistic_p": p,
        "percentile_of_observed": float(np.mean(null_max < obs_max)),
        "n_replicates": int(n),
        "p_convention": ("Phipson-Smyth exact: (#{null >= observed} + 1)/(N + 1). "
                         "Ties count against significance; p is never exactly zero."),
        "dependence": ("handled implicitly -- the null replicates carry the same "
                       "correlation among statistics that the observation does, so it "
                       "is reproduced rather than modelled"),
    }


def combine_max_statistic(per_cell_null_max, observed_max_over_cells):
    """Calibrate a maximum taken ACROSS cells of a grid, not just within one.

    Requires every cell to have been scored on the SAME replicates, so the
    per-replicate max-|z| vectors align and their element-wise maximum is the null for
    "the largest z anywhere in the grid". Scanning a grid and quoting the best cell's
    own p is the error this function exists to prevent.
    """
    stack = np.stack([np.asarray(a, dtype=np.float64) for a in per_cell_null_max])
    if stack.ndim != 2:
        raise ValueError("per_cell_null_max must be a sequence of equal-length arrays")
    gw = np.max(stack, axis=0)
    n = gw.size
    p = (int(np.sum(gw >= float(observed_max_over_cells))) + 1) / (n + 1)
    return {
        "n_cells": int(stack.shape[0]), "n_replicates": int(n),
        "observed_max_over_cells": float(observed_max_over_cells),
        "null_max_over_cells_p95": float(np.quantile(gw, 0.95)),
        "grid_wide_max_statistic_p": p,
        "requires": ("every cell scored on the SAME replicates; otherwise the "
                     "element-wise maximum is not a null for anything"),
    }


def span_sensitivity(values_by_span, bins):
    """Occupancy as a function of the span it was measured over.

    §P7-25(1)(i) flagged that D-1b's +-10-day span was not pre-declared and that the
    resulting constant is span-dependent even though the verdict was not. This makes
    that dependence a MEASUREMENT rather than a caveat: pass the same feature sampled
    over several spans and read off how much the null constant moves.
    """
    out = {}
    for span, vals in values_by_span.items():
        out[str(span)] = occupancy(vals, bins)
    names = list(next(iter(out.values())).keys())
    out["_range"] = {nm: (min(o[nm] for k, o in out.items() if k != "_range"),
                          max(o[nm] for k, o in out.items() if k != "_range"))
                     for nm in names}
    return out


# ------------------------------------------ SENSITIVITY: does the dial move at all? --
# Kepler's proposed rule 7, adopted here as machinery. `assert_power_floor` above asks
# whether a p CAN be small. Nothing asked whether the statistic RESPONDS to the effect
# it is supposed to measure. A statistic that cannot move is not a null result, it is a
# broken instrument, and quoting its null as evidence of absence is the same class of
# error as reading a concentration against a uniform-circle null.
#
# MEASURED, NOT ASSUMED, and the min-over-nuisance clause is the whole point: a bound
# must hold at the adversary's choice of unknown phase, not at the analyst's best one.

def circular_sensitivity(pool_angles, n, order, eps=0.30, n_phase=8, n_rep=60,
                         n_null=800, rng=None):
    """Power-relevant sensitivity of |R_m| and of the 2-dof complex statistic.

    Plants a density modulation `1 + eps*cos(m*(theta - theta0))` on the site's own
    dwell pool, sweeps the unknown phase theta0, and reports BOTH estimators.

    THE MEASURED FACT this exists to expose: `d|R_m|/d(eps)` is SIGNED and passes
    through zero, because a perturbation off-axis from an already-long dwell resultant
    ROTATES it instead of lengthening it, and |R| discards the rotation. The 2-dof
    Mahalanobis statistic on the complex resultant keeps the rotation and has a
    non-vanishing minimum.
    """
    rng = np.random.default_rng() if rng is None else rng
    pool = np.asarray(pool_angles, dtype=np.float64)
    zs = np.empty(int(n_null), dtype=complex)
    for i in range(int(n_null)):
        zs[i] = np.exp(1j * order * rng.choice(pool, size=int(n))).mean()
    r_null = np.abs(zs)
    mu = np.array([zs.real.mean(), zs.imag.mean()])
    cov = np.cov(np.vstack([zs.real, zs.imag]))
    inv = np.linalg.inv(cov)

    def maha(z):
        d = np.array([z.real, z.imag]) - mu
        return float(np.sqrt(d @ inv @ d))

    m_null = np.mean([maha(z) for z in zs[:200]])
    dR, dM = [], []
    for th0 in np.linspace(0.0, 2.0 * np.pi, int(n_phase), endpoint=False):
        w = np.clip(1.0 + eps * np.cos(order * (pool - th0)), 0.0, None)
        w = w / w.sum()
        rr, mm = [], []
        for _ in range(int(n_rep)):
            s = rng.choice(pool, size=int(n), p=w)
            z = np.exp(1j * order * s).mean()
            rr.append(abs(z))
            mm.append(maha(z))
        dR.append((np.mean(rr) - r_null.mean()) / eps)
        dM.append((np.mean(mm) - m_null) / eps)
    return {
        "order": int(order), "n": int(n), "eps": float(eps),
        "dR_by_phase": [float(x) for x in dR],
        "dM_by_phase": [float(x) for x in dM],
        "dR_min_abs": float(np.min(np.abs(dR))),
        "dR_max_abs": float(np.max(np.abs(dR))),
        "dR_changes_sign": bool(np.min(dR) < 0 < np.max(dR)),
        "dM_min": float(np.min(dM)),
        "two_sigma_on_R": float(2.0 * r_null.std(ddof=1)),
        "note": ("|R| sensitivity is signed and crosses zero, so it is blind at some "
                 "unknown phases; the 2-dof statistic has a non-vanishing minimum. "
                 "A bound must be quoted at the MINIMUM over phase, never the "
                 "maximum, because the phase is a nuisance parameter the adversary "
                 "chooses."),
    }


def assert_sensitivity_floor(sensitivity, min_dM=0.5):
    """Refuse a statistic whose worst-case response is below a declared floor.

    The companion to `assert_power_floor`: that one refuses a family that cannot reach
    significance, this one refuses a statistic that cannot register the effect.
    """
    got = float(sensitivity["dM_min"])
    if got < float(min_dM):
        raise ValueError(
            "sensitivity floor: worst-case 2-dof response %.4f per unit effect is "
            "below the declared floor %.4f -- this statistic cannot see what it "
            "claims to measure" % (got, min_dM))
    return {"dM_min": got, "floor": float(min_dM), "passed": True}
