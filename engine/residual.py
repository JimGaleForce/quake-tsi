"""PROJECTION AND RESIDUAL MACHINERY: turn "that axis is redundant" into a NUMBER.

WHY THIS EXISTS. K-108 established that D-1c's resolved-Coulomb statistics are ~89 %
variance-shared with the scalar axis D-1 had already bounded (`corr = -0.9506`,
measured), so D-1c was not independent evidence. The natural response is to discard the
arm. THAT RESPONSE THROWS AWAY THE 11 %.

A bound on one direction is not a bound on the space. The right operation is not
rejection but PROJECTION: split every candidate statistic into the part the bounded
axis already explains and the part it cannot, and then test only the second. That
converts a binary "redundant, discard" into a graded "this much is new, and here is
what it says", which is what a search for potentials worth looking into more deeply
actually needs.

The same move answers Kepler's proposed SP-9 (declare a battery's EFFECTIVE rank
alongside its nominal statistic count) and Wegener's W-008-P1 (regress the residual on
the artifact and require the slope to be zero). Both are contractions of the machinery
here.

WHAT IS AND IS NOT LEGITIMATE ABOUT THIS.

  * The projection is computed on the DETERMINISTIC FIELD, not on the events. It uses
    no catalogue, so choosing to residualise costs nothing and cannot be data-snooping.
  * The residual is a DIFFERENT STATISTIC from the raw one and needs its OWN
    waveform-matched null. It does not inherit the raw statistic's null, and
    `engine/dwell_null.py` must be re-run on the residual series. A residual scored
    against the raw statistic's null would be a new instance of the exact error D-1
    was killed by.
  * `shared_variance` is a property of the FIELD at a site and epoch, so it varies.
    Report its distribution across sites; a single pooled number hides the variation
    that decides whether a transfer coefficient is meaningful at all (K-098's own
    falsifier).
  * Projecting out an axis REDUCES the degrees of freedom, and the reduction must be
    counted. `effective_rank` is how many independent statistics a battery really
    contains, and it is what the multiplicity correction should be priced against.
"""

from __future__ import annotations

import numpy as np


def _demean(x):
    a = np.asarray(x, dtype=np.float64)
    return a - a.mean()


def project_out(target, bases):
    """Component of `target` orthogonal to the linear span of `bases`.

    `bases` is a sequence of series of the same length. Least squares on the demeaned
    columns, so the residual is the part of `target` that no linear combination of the
    already-bounded axes can produce. Returns the residual and the fraction of
    `target`'s variance the bases explain.
    """
    y = _demean(target)
    cols = [_demean(b) for b in bases]
    if not cols:
        return y, 0.0
    a = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(a, y, rcond=None)
    fit = a @ coef
    resid = y - fit
    vy = float(np.dot(y, y))
    shared = float(np.dot(fit, fit) / vy) if vy > 0 else 0.0
    return resid, min(1.0, max(0.0, shared))


def shared_variance(target, bases):
    """Just the r^2 that `bases` explain of `target`. The K-098 transfer coefficient."""
    return project_out(target, bases)[1]


def new_information_fraction(target, bases):
    """1 - r^2: the share of a candidate axis that a bound does NOT already cover.

    This is the number that should sit beside every proposed new statistic. An axis at
    0.11 is not worthless -- D-1c's Coulomb axis was exactly that -- but it must be
    priced as 0.11 of a new test rather than as a whole one.
    """
    return 1.0 - shared_variance(target, bases)


def effective_rank(series, tol=0.01):
    """How many independent statistics a battery really contains.

    Eigenvalues of the correlation matrix; the rank is the number needed to reach
    `1 - tol` of the total variance. Kepler's SP-9 asks for exactly this to be declared
    alongside `n_declared_tests`, because a nominal count of m badly overstates the
    breadth of a collinear search, and overstating it corrupts any correction computed
    from the count.
    """
    cols = [_demean(s) for s in series]
    a = np.column_stack(cols)
    sd = a.std(axis=0)
    keep = sd > 0
    if not np.any(keep):
        return {"n_nominal": len(cols), "effective_rank": 0,
                "eigenvalues": [], "note": "all columns constant"}
    a = a[:, keep] / sd[keep]
    c = np.corrcoef(a, rowvar=False)
    c = np.atleast_2d(c)
    ev = np.sort(np.linalg.eigvalsh(c))[::-1]
    ev = np.clip(ev, 0.0, None)
    tot = ev.sum()
    cum = np.cumsum(ev) / max(tot, 1e-300)
    rank = int(np.searchsorted(cum, 1.0 - tol) + 1)
    return {
        "n_nominal": len(cols),
        "n_non_constant": int(keep.sum()),
        "effective_rank": rank,
        "eigenvalues": ev.tolist(),
        "variance_explained_by_first": float(ev[0] / max(tot, 1e-300)),
        "tol": tol,
        "note": ("effective rank is the number of principal components needed to "
                 "reach 1 - tol of the total variance; a battery whose effective rank "
                 "is far below its nominal count is not the search it claims to be"),
    }


def axial_shared_variance(bearings_deg, bases):
    """Shared variance for an AXIAL (pi-periodic) quantity, done on the doubled angle.

    A bearing is not a linear quantity, so it cannot be regressed directly. Its
    doubled-angle cosine and sine components can, and the shared variance is the
    variance-weighted mean of the two. Doing this on the raw bearing instead would be
    the same class of error as testing a pi-periodic variable with a 2-pi statistic.
    """
    b = np.mod(np.asarray(bearings_deg, dtype=np.float64), 180.0)
    two = 2.0 * b * np.pi / 180.0
    c, s = np.cos(two), np.sin(two)
    rc, sh_c = project_out(c, bases)
    rs, sh_s = project_out(s, bases)
    vc, vs = float(np.var(c)), float(np.var(s))
    tot = vc + vs
    shared = (sh_c * vc + sh_s * vs) / tot if tot > 0 else 0.0
    return {"shared_variance": float(shared),
            "residual_cos": rc, "residual_sin": rs,
            "component_shared": {"cos2": sh_c, "sin2": sh_s}}


def orthogonality_report(candidates, bounded, names=None):
    """The table that should precede any new arm: what is new, and how new.

    `candidates` maps name -> series; `bounded` is the list of already-bounded axes.
    Returns per-candidate shared variance and new-information fraction, plus the
    battery's effective rank, so an arm can be judged BEFORE it is run rather than
    after it has spent a multiplicity budget on directions that were already closed.
    """
    rows = {}
    for k, v in candidates.items():
        sh = shared_variance(v, bounded)
        rows[k] = {"shared_variance": sh, "new_information_fraction": 1.0 - sh}
    rank = effective_rank(list(candidates.values()))
    order = sorted(rows, key=lambda k: rows[k]["new_information_fraction"],
                   reverse=True)
    return {
        "per_candidate": rows,
        "ranked_by_novelty": order,
        "battery_effective_rank": rank,
        "rule": ("Kepler SP-9 (proposed): declare a battery's EFFECTIVE rank alongside "
                 "its nominal statistic count. An arm whose effective rank is "
                 "materially below its nominal count must say so in its own result "
                 "file."),
    }
