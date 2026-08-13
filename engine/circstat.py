"""TRANCHE B's new circular statistics: F9-01 (second moment) and F9-04 (Kuiper/Watson).

Both are NEW ESTIMATORS, and §P7-3 attaches a demonstrated-recovery demand to each
one before any result from it means anything (S-17 candidate, §P6-7). The harnesses
that discharge those demands live in `engine/recovery_b.py`; this module is only the
statistics and their nulls. Nothing here is evidence.

WHAT F9-01 IS, AND THE ONE THING IT MUST NOT DO
-----------------------------------------------
MINING_CATALOG F9-01: *"for a candidate cycle, compute the second circular moment of
the event-phase distribution, |mean(exp(2 i theta))|, in addition to the first."*
§K87-0(c) is the reason: every statistic in the v1 session was a first-moment,
single-phase statistic, and a two-stage unlock/release process -- unlock at phase X
in cycle N, release at a dispersed phase in cycle N+k -- decays geometrically in the
first moment while surviving in the second.

On the count path the second moment is EXACTLY the engine's existing 2-df score
statistic evaluated on the DOUBLED ANGLE: the design becomes [sin 2t, cos 2t] instead
of [sin t, cos t]. That is not a convenience, it is the whole reason the negative
control below is decisive:

  > **§P7-3(1), the control Kepler did not name.** *A planted PURE SINUSOID -- the
  > first-harmonic alternative the programme has already bounded -- must NOT fire the
  > second-moment statistic beyond its nominal size.* Without it, "a second-moment
  > detection" is unfalsifiably confounded with the first moment.

The mathematical reason it holds is one line and it is worth writing down, because a
control that passes for a reason nobody can state is a control nobody can trust:
over a uniform phase grid `cos(2t)` and `sin(2t)` are ORTHOGONAL to `cos(t)` and
`sin(t)`, so a first-harmonic modulation projects to zero on the doubled-angle
design. It is not "small", it is zero up to the discreteness of the grid -- and the
test in `engine/tests/test_circstat.py` checks the empirical size, not the algebra,
because the algebra is not what could be wrong in an implementation.

The DESCRIPTIVE moments (the catalog's own `|mean(exp(k i theta))|` at k = 1 and 2)
are reported alongside, because F9-01's own Pit clause requires it: a second-moment
hit is ambiguous between an AXIAL response and a HARMONIC one, and the report must
say so with both numbers next to each other rather than quoting the second alone.

WHAT F9-04 IS, AND WHAT IT DOES *NOT* BUY
------------------------------------------
MINING_CATALOG F9-04: Kuiper's V and Watson's U^2 in place of / alongside the 2-df
quadratic form, because *"the 2-df quadratic form is optimal against a sinusoid and
weak against everything else. A sharp response confined to 10% of the cycle is
nearly invisible to it and obvious to Kuiper."* That is the money demonstration and
it is built as a test (`test_kuiper_sees_what_rayleigh_misses`).

  > **§P7-3(2), stated in the module because it is the thing most likely to be
  > overclaimed downstream.** Kuiper on DAILY-BINNED data inherits the same sinc
  > notches as every other count-path statistic. It buys SHAPE sensitivity. It does
  > NOT buy BAND coverage, and it may not be described as buying both. The
  > day-lattice negative control (`day_lattice_phase`) exists to keep that honest.

THE NULL, AND WHY IT IS THE SAME NULL
--------------------------------------
Kuiper and Watson have no closed-form null under a non-uniform baseline intensity
(F9-04's own Pit), so both are calibrated by the engine's existing two-null
discipline, unchanged and for the same reasons `mine.glm_task` gives:

  * CIRCULAR SHIFT of the (counts, offset) PAIR against a fixed phase vector. Both
    totals are preserved, so the null intercept is exactly invariant, and shift 0 is
    the observed value. Powerless for a deterministic cycle -- see below.
  * CIRCULAR MOVING-BLOCK BOOTSTRAP of the same pair in blocks of the feature's own
    `block_days`, which is the null that still has power when the feature IS a
    deterministic cycle.

`p_raw` is the more conservative of the two, except for a periodic feature where the
shift null is provably powerless and the bootstrap IS the null -- exactly the rule
`glm_task` already applies, carried here rather than re-derived.
"""

from __future__ import annotations

import math

import numpy as np

from . import mine as M

# The chunk width for the shift null: `n_shifts x n_days` float64 materialised at
# once. 256 x 7716 is 15 MB, which is the point of chunking it at all.
SHIFT_CHUNK = 256

# Phases closer than this are ONE phase. See `phase_group_ends` for why a tolerance
# is load-bearing rather than fussy: 1e-9 rad is a million times smaller than the
# ~1e-3 rad by which adjacent days differ in any real ephemeris cycle, and a million
# times LARGER than the 1e-13 float residue of `np.mod(2*pi*t, 2*pi)`.
PHASE_TIE_TOL = 1e-9

# Empirical-vs-expected CDF differences below this are float noise. See
# `_kuiper_watson_core` for why snapping them to exact zero is load-bearing.
CDF_NOISE_FLOOR = 1e-12

AMBIGUITY_NOTE = (
    "MINING_CATALOG F9-01 Pit: the second moment is driven by ANY two-lobed "
    "structure, including a semidiurnal or half-cycle response, so a hit is "
    "ambiguous between an AXIAL response and a HARMONIC one and must be reported as "
    "such. `R1` and `R2` are printed together for that reason -- a large R2 beside a "
    "large R1 is a harmonic reading, a large R2 beside R1 ~ 0 is the axial reading "
    "§K87-0(c) says no v1 statistic could have produced.")

BAND_NOTE = (
    "§P7-3(2): Kuiper/Watson on DAILY-BINNED phases inherit the day-binning sinc "
    "notches exactly as the 2-df form does. They buy SHAPE sensitivity, not BAND "
    "coverage. Any recovery claim from them is per band and must be demonstrated "
    "band by band; the diurnal and semidiurnal bands stay structurally notched.")


# ------------------------------------------------------- descriptive moments --
def circular_moment(theta, weights=None, k=1):
    """|mean(exp(i k theta))| and its mean angle, weight-aware. The catalog's own form.

    `weights` are event counts per phase sample (the count path) or all-ones (one
    row per event, the mark path). Returns (R, mean_angle_rad).
    """
    th = np.asarray(theta, dtype=np.float64)
    w = (np.ones_like(th) if weights is None
         else np.asarray(weights, dtype=np.float64))
    tot = float(w.sum())
    if tot <= 0:
        return 0.0, float("nan")
    z = complex(float((w * np.cos(k * th)).sum()) / tot,
                float((w * np.sin(k * th)).sum()) / tot)
    return float(abs(z)), float(np.angle(z))


def moments_report(theta, weights=None):
    """R1 and R2 side by side, with the ambiguity note attached (F9-01 Pit)."""
    r1, a1 = circular_moment(theta, weights, k=1)
    r2, a2 = circular_moment(theta, weights, k=2)
    return {"R1": r1, "mean_angle_1_rad": a1, "R2": r2, "mean_angle_2_rad": a2,
            "reading": ("AXIAL (R2 >> R1): the structure §K87-0(c) says no "
                        "first-moment statistic could have produced"
                        if r2 > 2.0 * r1 else
                        "HARMONIC-OR-AMBIGUOUS (R2 not dominant over R1)"),
            "ambiguity_note": AMBIGUITY_NOTE}


# ------------------------------------------------- F9-01 on the count path ---
def doubled_angle_design(theta):
    """[sin 2t, cos 2t] -- the second-moment design, column order matching Feature.design."""
    th = np.asarray(theta, dtype=np.float64)
    return np.column_stack([np.sin(2.0 * th), np.cos(2.0 * th)])


def second_moment_test(theta, counts, offset, n_surr, rng, block_days=90.0,
                       periodic=False, min_shift=30):
    """F9-01 on daily counts: the 2-df score statistic on the DOUBLED angle.

    Same two nulls, same conservatism rule and same `p_raw` convention as
    `mine_session.glm_task`, because it is the same statistic on a different design
    matrix -- and running it under a different null discipline would make its result
    incomparable with the first-moment result it exists to be compared against.
    """
    th = np.asarray(theta, dtype=np.float64)
    counts = np.asarray(counts, dtype=np.float64)
    offset = np.asarray(offset, dtype=np.float64)
    X2 = doubled_angle_design(th)
    S = M.score_stat_all_shifts(X2, counts, offset)
    p_shift, n_used = M.empirical_p(S, n_surr, rng, min_shift=min_shift)
    Sb = M.score_stat_block_bootstrap(X2, counts, offset, int(n_surr), rng,
                                      mean_block=float(block_days))
    p_boot = M.bootstrap_p(S[0], Sb)
    p = p_boot if periodic else max(p_shift, p_boot)

    fit = M.glm_fit(X2, counts, offset)
    b = np.asarray(fit["beta"])
    amp = float(np.hypot(*b))
    first = M.glm_fit(np.column_stack([np.sin(th), np.cos(th)]), counts, offset)
    b1 = np.asarray(first["beta"])

    out = {
        "test": "second_circular_moment_score", "df": 2,
        "statistic": float(S[0]), "chi2_score": float(S[0]),
        "p_parametric": M.chi2_sf(S[0], 2),
        "p_circular_shift": p_shift, "p_block_bootstrap": p_boot, "p_raw": p,
        "n_surrogates": min(n_used, int(n_surr)),
        "block_days": float(block_days),
        "null": ("block-bootstrap (periodic feature)" if periodic
                 else "max(circular-shift, block-bootstrap)"),
        "beta": fit["beta"], "se": fit["se"],
        "amplitude_log_rate": amp,
        "pct_rate_modulation": 100.0 * (math.exp(amp) - 1.0),
        "bits_per_event": fit["bits_per_event"],
        "converged": bool(fit["converged"]),
        # the first-moment fit on the SAME window, so the ambiguity reading is
        # available on the row itself and nobody has to join two tables to get it
        "first_moment_amplitude_log_rate": float(np.hypot(*b1)),
        "harmonic_over_fundamental": (
            float(amp / max(float(np.hypot(*b1)), 1e-12))),
        "ambiguity_note": AMBIGUITY_NOTE,
    }
    out.update(moments_report(th, counts))
    return out


# ----------------------------------------- F9-04: Kuiper V and Watson U^2 ----
def wrap_phase(theta, tie_tol=PHASE_TIE_TOL):
    """Reduce to [0, 2pi) AND close the seam at 2pi. Both halves are load-bearing.

    `np.mod(2*pi*t, 2*pi)` for integer `t` does not merely wobble around zero -- it
    lands on BOTH SIDES of the seam, returning ~1e-12 for some `t` and ~2pi - 1e-12
    for others. Sorted, those are the first and last elements of the array with the
    whole circle between them, so a tie-grouping that only compares neighbours sees
    TWO phase groups where the construction has one, and the day-lattice negative
    control fires on the seam. Snapping the near-2pi tail to 0 closes it.
    """
    th = np.mod(np.asarray(theta, dtype=np.float64), 2 * np.pi)
    return np.where(2 * np.pi - th < float(tie_tol), 0.0, th)


def phase_group_ends(theta_sorted, tie_tol=PHASE_TIE_TOL):
    """Indices of the LAST sample in each run of equal phases. Ties are the point.

    THE DAY-LATTICE NEGATIVE CONTROL DEPENDS ENTIRELY ON THIS. A feature whose period
    divides the 1-day sampling lattice has ONE distinct phase value: every day is
    tied. Evaluate the empirical-vs-expected CDF difference at every sample index and
    the "sort order" is just time order, so Kuiper silently becomes a test for DRIFT
    in the record -- and it fires, on a construction whose only structure is the
    grid, which is the exact false positive §P7-3(2) demands it not produce.
    Evaluating only at DISTINCT phase values makes V and U^2 identically zero there,
    which is the correct answer and the one the control checks.

    A TOLERANCE, and why it is not fussiness. `np.mod(2*pi*t, 2*pi)` for integer `t`
    does not return exact zeros -- it returns residuals of order 1e-13 that DIFFER
    per `t`. Exact-equality grouping would then see 7,716 distinct phases where there
    is one, and the day-lattice control would fire on floating-point noise. Real
    ephemeris phases differ by ~1e-3 rad between adjacent days, a million times the
    tolerance, so nothing real is ever merged by it.
    """
    th = np.asarray(theta_sorted, dtype=np.float64)
    n = th.size
    end = np.empty(n, dtype=bool)
    if n:
        end[-1] = True
        if n > 1:
            end[:-1] = np.abs(th[1:] - th[:-1]) > float(tie_tol)
    return np.nonzero(end)[0]


def _kuiper_watson_core(obs, exp, ends):
    """(V, U2) from phase-SORTED weight vectors, evaluated at DISTINCT phases only.

    obs/exp are (..., n): any leading batch dimension is carried, which is what makes
    the surrogate loops one call instead of n_surr calls. `ends` comes from
    `phase_group_ends` and is a property of the phase vector alone, so it is computed
    once and reused across every shift and every bootstrap replicate.
    """
    obs = np.asarray(obs, dtype=np.float64)
    exp = np.asarray(exp, dtype=np.float64)
    so = obs.sum(axis=-1, keepdims=True)
    se = exp.sum(axis=-1, keepdims=True)
    Fn = np.cumsum(obs, axis=-1)[..., ends] / np.maximum(so, 1e-300)
    F0 = np.cumsum(exp, axis=-1)[..., ends] / np.maximum(se, 1e-300)
    d = Fn - F0
    # CDF differences below CDF_NOISE_FLOOR are float noise, not statistics, and
    # snapping them to EXACT zero is what makes the degenerate case (one distinct
    # phase -> V and U^2 identically 0) produce p = 1 instead of p = 1/(B+1). Without
    # it, an observed 1e-19 beats a surrogate 0.0 and the day-lattice control reports
    # the resolution floor on a statistic that is mathematically zero.
    d = np.where(np.abs(d) < CDF_NOISE_FLOOR, 0.0, d)
    # sup and inf include the origin, where both CDFs are 0 -- so D+ and D- are
    # non-negative by construction and V is a genuine Kuiper statistic.
    V = np.maximum(d.max(axis=-1), 0.0) + np.maximum((-d).max(axis=-1), 0.0)
    p = np.diff(F0, axis=-1, prepend=0.0)            # dF0 over distinct phases
    dbar = (d * p).sum(axis=-1, keepdims=True)
    U2 = (((d - dbar) ** 2) * p).sum(axis=-1)        # Watson's U^2 (unscaled)
    return V, U2


def _scaled(V, U2, n_eff):
    """Stephens' n-scalings, so the statistic is comparable across surrogates.

    This matters more than it looks: a moving-block bootstrap resamples DAYS, so the
    total event count moves between replicates. Comparing a raw V computed on 46,585
    events against one computed on 46,402 is comparing two different statistics. The
    scaled forms are monotone in V and U2 at fixed n, so the OBSERVED p-value is
    unchanged by this choice -- it is the surrogates that need it.
    """
    n = np.maximum(np.asarray(n_eff, dtype=np.float64), 1.0)
    rn = np.sqrt(n)
    # Stephens' additive small-sample shift (-0.1 + 0.1/n) is DELIBERATELY OMITTED
    # from U2: it exists to make the asymptotic closed-form tail apply, and this
    # statistic is calibrated by surrogates rather than by that tail. Keeping it
    # would drive a degenerate U2 = 0 (the day-lattice control) to a NEGATIVE
    # statistic, which is meaningless and would have to be explained in every report
    # forever. What is kept is the n-scaling, which is what the surrogates need.
    return (V * (rn + 0.155 + 0.24 / rn), U2 * n * (1.0 + 0.8 / n))


def kuiper_watson(theta, counts, offset, order=None, ends=None):
    """Observed (V*, U2*) for daily phases `theta` with observed/expected day weights.

    The null CDF is the ETAS baseline's own phase distribution -- `offset`, not
    uniform -- which is F9-04's Pit clause ("no closed-form null under a non-uniform
    baseline intensity") handled by construction rather than by simulation alone.
    """
    th = wrap_phase(theta)
    o = np.argsort(th, kind="stable") if order is None else order
    e = phase_group_ends(th[o]) if ends is None else ends
    obs = np.asarray(counts, dtype=np.float64)[o]
    exp = np.asarray(offset, dtype=np.float64)[o]
    V, U2 = _kuiper_watson_core(obs, exp, e)
    Vs, U2s = _scaled(V, U2, obs.sum())
    return {"V": float(V), "U2": float(U2), "V_star": float(Vs),
            "U2_star": float(U2s), "n_events": float(obs.sum()),
            "n_distinct_phases": int(e.size)}


def _kw_batch(counts, offset, idx, order, ends):
    """(V*, U2*) for a batch of index vectors `idx` (b, n) into (counts, offset)."""
    obs = np.asarray(counts, dtype=np.float64)[idx][:, order]
    exp = np.asarray(offset, dtype=np.float64)[idx][:, order]
    V, U2 = _kuiper_watson_core(obs, exp, ends)
    return _scaled(V, U2, obs.sum(axis=-1))


def kuiper_watson_shift_null(theta, counts, offset, n_surr, rng, min_shift=30,
                             chunk=SHIFT_CHUNK):
    """Circular-shift surrogates of the (counts, offset) PAIR against fixed phases.

    Exactly `score_stat_all_shifts`' construction -- shifting the pair preserves both
    totals, so the null intercept is invariant -- but evaluated at a SAMPLE of
    admissible shifts rather than all of them, because Kuiper needs a sort-ordered
    cumulative sum per shift and the full n x n matrix is not affordable at
    production n. The sample is drawn without replacement exactly as
    `mine.empirical_p` draws it, and `n_used` is reported so the resolution floor
    1/(n_used+1) is a stated property of the run.

    Returns (p_V, p_U2, n_used, obs_dict).
    """
    th = wrap_phase(theta)
    n = th.size
    order = np.argsort(th, kind="stable")
    ends = phase_group_ends(th[order])
    obs = kuiper_watson(th, counts, offset, order=order, ends=ends)

    ok = np.arange(n)
    ok = ok[(ok >= int(min_shift)) & (ok <= n - int(min_shift))]
    if ok.size == 0:
        raise ValueError("no admissible circular shifts")
    if int(n_surr) < ok.size:
        ok = rng.choice(ok, size=int(n_surr), replace=False)
    base = np.arange(n)
    ge_v = ge_u = 0
    for a in range(0, ok.size, int(chunk)):
        k = ok[a:a + int(chunk)]
        idx = (base[None, :] - k[:, None]) % n       # roll the pair by k
        Vs, U2s = _kw_batch(counts, offset, idx, order, ends)
        ge_v += int((Vs >= obs["V_star"]).sum())
        ge_u += int((U2s >= obs["U2_star"]).sum())
    return ((1.0 + ge_v) / (1.0 + ok.size), (1.0 + ge_u) / (1.0 + ok.size),
            int(ok.size), obs)


def kuiper_watson_block_bootstrap(theta, counts, offset, n_boot, rng,
                                  mean_block=90.0, chunk=200):
    """Moving-block bootstrap surrogates of the (counts, offset) pair, blocks in days.

    Same resampler as `mine.score_stat_block_bootstrap` (`mine.block_bootstrap_idx`),
    so the two statistics face the same null construction and their p-values are
    comparable. Returns (V*_draws, U2*_draws).
    """
    th = wrap_phase(theta)
    n = th.size
    order = np.argsort(th, kind="stable")
    ends = phase_group_ends(th[order])
    outV = np.empty(int(n_boot))
    outU = np.empty(int(n_boot))
    done = 0
    while done < int(n_boot):
        b = int(min(chunk, int(n_boot) - done))
        idx = M.block_bootstrap_idx(n, b, mean_block, rng)
        Vs, U2s = _kw_batch(counts, offset, idx, order, ends)
        outV[done:done + b] = Vs
        outU[done:done + b] = U2s
        done += b
    return outV, outU


def omnibus_test(theta, counts, offset, n_surr, rng, block_days=90.0,
                 periodic=False, min_shift=30):
    """F9-04: Kuiper V and Watson U^2 as TWO rows, under the engine's two nulls.

    Two rows and not one, because F9-04 prices `2 statistics x 17 cyclic features`
    and because its own Pit clause requires both to be reported: they are less
    powerful than the 2-df form when the truth IS sinusoidal, and a run that quoted
    only whichever won would be a forking path with a statistic for a knob.
    """
    p_v, p_u, n_used, obs = kuiper_watson_shift_null(
        theta, counts, offset, n_surr, rng, min_shift=min_shift)
    bV, bU = kuiper_watson_block_bootstrap(theta, counts, offset, int(n_surr), rng,
                                           mean_block=float(block_days))
    rows = []
    for name, stat, p_shift, draws in (
            ("kuiper_V", obs["V_star"], p_v, bV),
            ("watson_U2", obs["U2_star"], p_u, bU)):
        p_boot = M.bootstrap_p(stat, draws)
        rows.append({
            "test": name, "df": None,
            "statistic": float(stat),
            "V": obs["V"], "U2": obs["U2"],
            "V_star": obs["V_star"], "U2_star": obs["U2_star"],
            "n_events": obs["n_events"],
            "p_circular_shift": float(p_shift), "p_block_bootstrap": float(p_boot),
            "p_raw": float(p_boot if periodic else max(p_shift, p_boot)),
            "n_surrogates": int(min(n_used, int(n_surr))),
            "block_days": float(block_days),
            "null": ("block-bootstrap (periodic feature)" if periodic
                     else "max(circular-shift, block-bootstrap)"),
            "band_note": BAND_NOTE,
        })
    return rows


# ------------------------------------------- the day-lattice negative control -
def day_lattice_phase(n_days, period_days=1.0, phase0=0.0):
    """§P7-3(2)'s negative control: a construction whose ONLY structure is the grid.

    A feature whose period divides the 1-day sampling lattice exactly is CONSTANT
    over daily bins, so any "structure" a statistic reports on it is the lattice
    talking. Kuiper and Watson must not fire on it. Returned as a phase vector so it
    can be pushed through the unmodified statistic rather than special-cased.
    """
    t = np.arange(int(n_days), dtype=np.float64)
    # The fraction is taken BEFORE the 2pi scaling: for a period that divides the
    # lattice exactly, `t / P` is an exact float integer and `mod 1` is exactly 0, so
    # the construction really does have one phase rather than 7,716 phases that
    # differ in the last bit.
    frac = np.mod(t / float(period_days) + float(phase0) / (2 * np.pi), 1.0)
    return 2 * np.pi * frac


def narrow_arc_intensity(theta, amplitude, duty=0.10, centre=0.0):
    """A ~`duty` duty-cycle boxcar response on the unit circle, mean-preserving.

    F9-04's positive control (§P7-3(2)): *"a planted narrow-arc response (~10% duty
    cycle) at a declared amplitude, which is precisely the alternative the 2-df form
    is weak against."* The multiplier is `1 + amplitude` inside the arc and is
    reduced outside so the CYCLE-AVERAGED rate is unchanged -- otherwise the plant
    would also be a rate change and the recovery would be partly trivial.

    Why the 2-df form is weak against it, derived rather than asserted (S-18's
    fourth clause: a stated direction without a derivation is a guess wearing the
    grammar of a result). Write the mean-preserving multiplier as
    `1 + K (1_arc - d)` with `K = A / (1 - d)`. Its first harmonic has magnitude

        a1 = (K / pi) * |1 - exp(-2 pi i d)| = 2 A sin(pi d) / (pi (1 - d))

    which at d = 0.10 is **0.219 A**. So a 60% narrow-arc excursion presents ~13% to
    a statistic that only looks at the first harmonic, and the whole rest of the
    excursion is invisible to it and visible to Kuiper. `fundamental_coefficient`
    returns exactly this number and the demonstration test asserts against it.
    """
    th = np.mod(np.asarray(theta, dtype=np.float64) - float(centre), 2 * np.pi)
    inside = th < (2 * np.pi * float(duty))
    m = np.where(inside, 1.0 + float(amplitude),
                 1.0 - float(amplitude) * float(duty) / (1.0 - float(duty)))
    return m


def k_fold_arc_intensity(theta, amplitude, n_arcs=3, duty_each=0.05):
    """`n_arcs` equally spaced narrow arcs, mean-preserving. THE DECISIVE PLANT.

    For arcs at `theta_j = 2 pi j / n`, the Fourier coefficient at harmonic `m` is
    nonzero ONLY when `m == 0 (mod n)`: the arc positions form a complete residue
    system, so `sum_j exp(-i m theta_j) = 0` unless `n | m`. At `n_arcs = 3` that
    kills harmonics 1 AND 2 exactly.

    That is what makes this the decisive F9-04 demonstration rather than a rhetorical
    one. The plant is a genuinely CONCENTRATED phase response -- three sharp arcs of
    5% duty each -- and:

      * the **2-df Rayleigh-form** statistic reads the FIRST harmonic: blind, exactly;
      * the **F9-01 second-moment** statistic reads the SECOND: blind, exactly;
      * **Kuiper and Watson** read the whole CDF: they see it.

    So the same catalogue separates all three statistics, and it also shows F9-01 and
    F9-04 are COMPLEMENTARY rather than two names for the same repair -- which is
    worth more than either demonstration alone.
    """
    ph = np.mod(np.asarray(theta, dtype=np.float64), 2 * np.pi)
    inside = np.zeros(ph.shape, dtype=bool)
    for j in range(int(n_arcs)):
        inside |= np.mod(ph - 2 * np.pi * j / float(n_arcs), 2 * np.pi) \
            < 2 * np.pi * float(duty_each)
    d = float(n_arcs) * float(duty_each)
    return np.where(inside, 1.0 + float(amplitude),
                    1.0 - float(amplitude) * d / (1.0 - d))


def kuiper_equivalent_amplitude(arc_amplitude, duty_each):
    """The sinusoidal amplitude a narrow arc presents to KUIPER, derived.

    Kuiper's V is a CDF excursion. A mean-preserving arc of height `A` and duty `d`
    accumulates a CDF excess of exactly `A * d` across the arc (density excess
    `A/(2 pi)` over an interval of length `2 pi d`, times the `(1-d)` normalisation
    that cancels), so `V_arc = A d`. A sinusoid of amplitude `a` gives
    `V_sin = 2 * a / (2 pi) = a / pi`. Equating them:

        a_equivalent = pi * A * d

    This is the number a plant must be sized by if "Kuiper recovered it" is to mean
    anything against a floor quoted in sinusoidal-amplitude units -- and sizing by
    `A` instead would compare a CDF excursion against a Fourier coefficient, which is
    a units error of exactly the kind §P7-8(d) exists to catch.
    """
    return float(np.pi * float(arc_amplitude) * float(duty_each))


def arc_amplitude_for_kuiper_equivalent(a_equivalent, duty_each):
    """Invert `kuiper_equivalent_amplitude`: the arc height for a target a_eq."""
    return float(a_equivalent) / (np.pi * float(duty_each))


# The ceiling on what a SINGLE mean-preserving boxcar can buy Kuiper over the 2-df
# form, derived rather than asserted (S-18's fourth clause):
#
#   a_kuiper / a_rayleigh = [pi A d] / [2 A sin(pi d) / (pi (1-d))]
#                         = pi^2 d (1-d) / (2 sin(pi d))   -->  pi/2 as d -> 0.
#
# So a single narrow arc gives Kuiper AT MOST a factor 1.571 in equivalent amplitude,
# which is a real edge and a modest one. MINING_CATALOG F9-04's "nearly invisible to
# it and obvious to Kuiper" overstates the single-arc case; the decisive case is a
# plant whose FUNDAMENTAL IS SUPPRESSED, which is `k_fold_arc_intensity`. This is
# recorded here rather than in a report, because the next person to size a narrow-arc
# plant will read the code and not the report.
SINGLE_ARC_KUIPER_EDGE_CEILING = float(np.pi / 2.0)


def single_arc_kuiper_edge(duty):
    """The measured-in-advance edge ratio at this duty. See the constant above."""
    d = float(duty)
    return float(np.pi ** 2 * d * (1.0 - d) / (2.0 * np.sin(np.pi * d)))


def fundamental_coefficient(amplitude, duty):
    """The first-harmonic amplitude a narrow-arc plant presents to a 2-df test.

    `2 A sin(pi d) / (pi (1 - d))` -- see `narrow_arc_intensity` for the derivation.
    This is the number that decides whether a narrow-arc demonstration is honest: if
    it sits ABOVE the S-15 floor the 2-df form should have found the plant and a
    "Rayleigh missed it" claim would be a power accident, not a shape result.
    """
    d = float(duty)
    return 2.0 * float(amplitude) * math.sin(math.pi * d) / (math.pi * (1.0 - d))
