"""SYMMETRY-AWARE ANGULAR TESTING. The harmonic that carries the signal is MEASURED.

THE POINT, IN ONE EXAMPLE, BECAUSE IT IS THE WHOLE MODULE.

Twenty events. Ten fire when the tide pulls from the left, ten when it pulls from the
right. If the mechanism is "horizontal extension across the fault opens it", that is
TWENTY OUT OF TWENTY -- a total effect. Every standard first-order circular test
(Rayleigh, Schuster, the mean resultant) reports it as PERFECTLY UNIFORM, because the
two lobes cancel:

    R1 = |mean(exp(i*theta))|   = 0.0     <- "no effect"
    R2 = |mean(exp(2i*theta))|  = 1.0     <- maximal

**A NULL ON R1 IS NOT A NULL ON THE HYPOTHESIS.** It is a null on the ASYMMETRIC part
of the hypothesis, and if the physics is symmetric that is the part that was never
going to be there. This program has run first-moment statistics and one hand-chosen
doubled one; it has never asked WHICH harmonic carries the concentration.

THE FIX. Do not guess the symmetry. Compute the whole angular spectrum

    R_m = |mean(exp(i * m * theta))|,  m = 1 .. M

and report every order, with the multiplicity of having looked at M of them priced.
Each order corresponds to a physical symmetry, and naming them is what makes the test
interpretable rather than a fishing expedition:

    m = 1   DIRECTIONAL. Pull from the north differs from pull from the south.
            Requires a mechanism that distinguishes a direction from its opposite.
    m = 2   AXIAL (left/right equivalent). The fundamental domain is 180 degrees.
            THIS IS THE ONE A SYMMETRIC MECHANISM LIVES IN, and the one that a
            first-moment test is blind to. Fault-normal extension, principal-axis
            alignment, and "pull from either side opens it" all sit here.
    m = 4   QUADRANTAL (left/right AND forward/back equivalent). The fundamental
            domain is 90 degrees. Shear-driven mechanisms with conjugate planes sit
            here: maximum shear occurs at 45 degrees to the principal axes, so a
            shear-controlled process is naturally four-lobed.
    m >= 3  reported because the data is entitled to say something the list above did
            not anticipate. An unexpected m = 3 or m = 5 is far more likely to be an
            artifact or a sampling lattice than a mechanism, and must be treated that
            way, but it is not this module's job to hide it.

THE SAME ARGUMENT APPLIES TO PARITY, NOT ONLY TO ANGLE. The degree-2 tidal potential
goes as P2(cos z), which is EVEN: it is identical for the Moon overhead and the Moon
underfoot. That was described elsewhere in this program as the degree-2 potential being
"blind" to the sign, i.e. as a limitation. IT MAY BE THE PHYSICS. If a body-tide
mechanism is symmetric under overhead/underfoot then the EVEN coordinate is the correct
one and the odd coordinate is testing something that cannot matter. `fold_parity`
provides both so the question is settled by measurement.

WHAT THIS MODULE DOES NOT DO. It does not supply the null. An angular concentration
must be scored against a WAVEFORM-MATCHED null from `engine/dwell_null.py`, because a
tidal bearing's own dwell distribution is strongly non-uniform (axial R2 measured from
0.14 at Alaska latitudes to 0.59 at the equator) and reading a concentration against a
uniform-circle null would manufacture an effect. The Rayleigh forms here are provided
ONLY as an analytic cross-check on synthetic uniform data and are labelled as such.
"""

from __future__ import annotations

import math

import numpy as np

SYMMETRY_MEANING = {
    1: "DIRECTIONAL: pull from one side differs from pull from the opposite side",
    2: "AXIAL: left/right equivalent, 180 deg fundamental domain -- where a symmetric "
       "mechanism lives, and where a first-moment test is blind",
    3: "three-lobed: no standard tidal mechanism; treat as artifact until shown "
       "otherwise",
    4: "QUADRANTAL: left/right AND forward/back equivalent, 90 deg domain -- conjugate "
       "shear planes sit here",
    5: "five-lobed: no standard tidal mechanism; treat as artifact until shown "
       "otherwise",
    6: "six-lobed: no standard tidal mechanism; treat as artifact until shown "
       "otherwise",
}

THE_TWENTY_EVENT_EXAMPLE = (
    "Ten events with the tide pulling from the left and ten from the right is TWENTY "
    "OUT OF TWENTY for a symmetric mechanism, and it gives R1 = 0 and R2 = 1. A "
    "first-moment circular test calls that perfectly uniform. A null on R1 is a null "
    "on the ASYMMETRIC part of the hypothesis and says nothing about the symmetric "
    "part.")


def harmonic_resultant(angles_rad, order):
    """R_m = |mean(exp(i m theta))|, the m-th trigonometric moment's length.

    R_m = 1 means the angles are perfectly concentrated once folded into 2pi/m;
    R_m = 0 means no concentration AT THAT SYMMETRY ORDER, which is not the same as
    no concentration.
    """
    a = np.asarray(angles_rad, dtype=np.float64)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan")
    z = np.exp(1j * float(order) * a)
    return float(np.abs(z.mean()))


def harmonic_phase(angles_rad, order):
    """The preferred direction at order m, back in the original coordinate.

    Returned in [0, 2pi/m): the fundamental domain of that symmetry. At m = 2 this is
    the preferred AXIS, in [0, 180) degrees.
    """
    a = np.asarray(angles_rad, dtype=np.float64)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan")
    z = np.exp(1j * float(order) * a).mean()
    return float(np.mod(np.angle(z) / float(order), 2.0 * np.pi / float(order)))


def harmonic_spectrum(angles_rad, max_order=6):
    """Every R_m and its preferred direction, m = 1..max_order.

    This is the object that replaces "we computed the mean resultant and it was flat".
    """
    a = np.asarray(angles_rad, dtype=np.float64)
    a = a[np.isfinite(a)]
    return {
        "n": int(a.size),
        "max_order": int(max_order),
        "R": {m: harmonic_resultant(a, m) for m in range(1, int(max_order) + 1)},
        "preferred_direction_rad": {m: harmonic_phase(a, m)
                                    for m in range(1, int(max_order) + 1)},
        "meaning": {m: SYMMETRY_MEANING.get(m, "higher order")
                    for m in range(1, int(max_order) + 1)},
        "multiplicity_note": ("max_order harmonics were examined; the correction owed "
                              "is over max_order, and the max-statistic across them "
                              "must be calibrated against the same waveform-matched "
                              "null the observation used"),
        "the_example": THE_TWENTY_EVENT_EXAMPLE,
    }


def fold(angles_rad, order):
    """Fold into the fundamental domain [0, 2pi/order) of an order-m symmetry.

    At order 2 this maps 'pull from the left' and 'pull from the right' onto the same
    value, which is the operation a symmetric mechanism requires and which a signed
    statistic refuses to perform.
    """
    a = np.asarray(angles_rad, dtype=np.float64)
    return np.mod(a, 2.0 * np.pi / float(order))


def fold_parity(x):
    """Split a signed quantity into its EVEN and ODD parts about zero.

    The angular argument applied to a linear coordinate. For elevation, the even part
    is 'how far from the horizon, either way' -- the overhead/underfoot-symmetric
    coordinate that the degree-2 potential actually lives in -- and the odd part is
    'above or below', which only a degree-3-like asymmetry can carry.
    """
    v = np.asarray(x, dtype=np.float64)
    # `odd` is v itself. The first version wrote sign(v)*v, which is |v| -- it would
    # have made the odd coordinate identical to the even one and quietly turned the
    # parity test into a comparison of a thing with itself.
    return {"even": np.abs(v), "odd": v.copy(),
            "sign": np.sign(v),
            "note": ("EVEN = |x|, the symmetric coordinate; ODD = x itself, which a "
                     "symmetric mechanism cannot use. P2(cos z) is even, so a "
                     "degree-2 body-tide mechanism is symmetric under "
                     "overhead/underfoot by construction and the EVEN coordinate is "
                     "the physically correct one to test.")}


def relative_angle(theta_rad, reference_rad, order=2):
    """Angle of `theta` relative to a `reference` direction, folded at `order`.

    The natural coordinate for "how does the tidal pull sit relative to the fault
    strike". At order 2 a pull along strike and a pull anti-along strike are the same
    state, which is what a fault plane actually sees.
    """
    d = np.asarray(theta_rad, dtype=np.float64) - np.asarray(reference_rad,
                                                             dtype=np.float64)
    return fold(d, order)


# ---------------------------------------------------- analytic cross-check only --
def rayleigh_p_uniform(r, n):
    """Rayleigh p for R_m under a TRULY UNIFORM null. Cross-check only.

    NOT for use on real tidal angles: a tidal bearing's own dwell distribution is far
    from uniform, so this p is meaningless there and would manufacture an effect. It
    exists to validate the harmonic machinery against synthetic uniform data, which is
    the one case where the analytic answer is known.
    """
    n = int(n)
    if n <= 0:
        return float("nan")
    z = n * float(r) ** 2
    return float(math.exp(-z) * (1.0 + (2.0 * z - z * z) / (4.0 * n)))


def dominant_order(spectrum, tol=0.05):
    """The FUNDAMENTAL symmetry order: the LOWEST m carrying the concentration.

    NOT simply the largest R_m, and the distinction is not cosmetic. A pattern with
    order-2 symmetry is automatically order-4, order-6 and every other multiple of 2:
    the exact ten-left/ten-right set gives R2 = R4 = 1.0. Taking the argmax would
    report order 4 or higher at random among the ties and would name a four-lobed shear
    mechanism where the data shows a two-lobed axial one.

    So: find the largest R over the spectrum, then report the SMALLEST order within
    `tol` of it. Divisors of the true order are what a real symmetry produces, and the
    smallest of them is the one with the mechanism attached.

    `also_high` lists the other orders within tolerance, which for a genuine order-k
    symmetry should be exactly its multiples -- a useful consistency check, since an
    order-2 result that does NOT also light up order 4 is incoherent.
    """
    r = spectrum["R"]
    rmax = max(r.values())
    within = sorted(m for m in r if r[m] >= rmax - float(tol))
    best = int(within[0])
    rest = sorted((v for m, v in r.items() if m not in within), reverse=True)
    return {
        "order": best,
        "R": float(r[best]),
        "meaning": SYMMETRY_MEANING.get(best, "higher order"),
        "also_high": [int(m) for m in within if m != best],
        "expected_also_high": [int(m) for m in r if m % best == 0 and m != best],
        "margin_over_next_independent": (float(r[best] - rest[0]) if rest
                                         else float("nan")),
        "why_lowest_not_largest": (
            "an order-k symmetry implies every multiple of k, so the argmax over "
            "orders picks a multiple at random among ties and names the wrong "
            "mechanism; the fundamental order is the smallest one carrying it"),
        "caution": ("a dominant order is a DESCRIPTION until it is calibrated against "
                    "a waveform-matched null and corrected for having examined every "
                    "order; orders 3, 5 and 6 have no standard tidal mechanism and "
                    "should be read as artifact first"),
    }


# ------------------------------- WORLD AND REGIONS: pooling that cannot hide it --
# Jim's second requirement: "test everything, the world AND specific REGIONS as some
# regions may use this triggering and others not."
#
# THE TRAP POOLING SETS. If half the regions carry a real order-2 concentration at one
# preferred axis and the other half carry it at the perpendicular axis, the pooled
# resultant is ZERO and the world-level test reports nothing while every region
# individually is maximal. That is the ten-left/ten-right problem again, one level up:
# a null on the pooled estimate is not a null on the hypothesis.
#
# So a regional analysis owes THREE numbers, never one:
#   1. the per-region effect, each against ITS OWN waveform-matched null,
#   2. the pooled effect, weighted by sqrt(n),
#   3. a HETEROGENEITY statistic, because genuine disagreement between regions is
#      itself the finding when the hypothesis is "some regions use this and others do
#      not". Cochran's Q, and a large Q means the pooling was inappropriate and the
#      pooled number must not be quoted as though it summarised anything.

def combine_regions(z_by_region, n_by_region):
    """Pooled z, Cochran's Q and I^2 across regions.

    Each region supplies a z already scored against ITS OWN null -- pooling raw angles
    across regions would be meaningless, since each site has a different dwell
    distribution. Weights are sqrt(n), so

        z_pooled = sum(w_i z_i) / sqrt(sum(w_i^2))

    which is the standard inverse-variance combination when each z has unit variance.
    """
    keys = sorted(z_by_region)
    z = np.array([float(z_by_region[k]) for k in keys])
    n = np.array([float(n_by_region[k]) for k in keys])
    ok = np.isfinite(z) & (n > 0)
    z, n, keys = z[ok], n[ok], [k for k, m in zip(keys, ok) if m]
    if z.size == 0:
        return {"n_regions": 0, "z_pooled": float("nan")}
    w = np.sqrt(n)
    z_pooled = float(np.sum(w * z) / math.sqrt(np.sum(w * w)))

    # COCHRAN'S Q ON UNIT WEIGHTS, and the reason is not a detail. Each z_i is ALREADY
    # standardised against its own region's null, so it has unit variance whatever n_i
    # is -- the region's size is inside the z, not beside it. Weighting Q by n a
    # second time counts the information twice and manufactures heterogeneity: the
    # first version of this function did exactly that and called four regions at
    # z = 1.8 to 2.2 significantly heterogeneous at p = 0.033.
    #
    # CAVEAT THAT TRAVELS WITH THE NUMBER: on the z scale this tests heterogeneity of
    # DETECTABILITY-SCALED effect. If the regions differ greatly in n, a common true
    # effect still yields different z, which reads as heterogeneity. Where that matters
    # the caller should pass effect sizes and standard errors instead of z.
    z_bar = float(np.mean(z))
    q = float(np.sum((z - z_bar) ** 2))
    df = int(z.size - 1)
    i2 = float(max(0.0, (q - df) / q)) if q > 0 and df > 0 else 0.0
    try:
        from scipy.stats import chi2
        q_p = float(chi2.sf(q, df)) if df > 0 else float("nan")
    except Exception:                                   # pragma: no cover
        q_p = float("nan")
    return {
        "n_regions": int(z.size),
        "regions": keys,
        "z_by_region": {k: float(v) for k, v in zip(keys, z)},
        "z_pooled": z_pooled,
        "z_weighted_mean": z_bar,
        "cochran_Q": q, "Q_df": df, "Q_p": q_p,
        "I2": i2,
        "heterogeneous": bool(q_p == q_p and q_p < 0.05),
        "interpretation": (
            "z_pooled is only a summary if Q is small. A LARGE Q with a small "
            "z_pooled is the signature of regions that disagree -- which is exactly "
            "what 'some regions use this triggering and others do not' predicts, and "
            "is a FINDING rather than a failure. Quoting z_pooled alone in that case "
            "reports nothing while the data says something."),
        "weights": "sqrt(n) per region for the POOLED z; UNIT weights for Q, because each z is already standardised against its own null",
    }


def per_region_spectra(angles_by_region, max_order=6):
    """The harmonic spectrum of every region, plus the axis each one prefers.

    The preferred AXIS is reported per region precisely so that regions concentrating
    on perpendicular axes are visible as such rather than cancelling silently in a
    pooled resultant.
    """
    out = {}
    for k, a in angles_by_region.items():
        sp = harmonic_spectrum(a, max_order=max_order)
        out[k] = {
            "n": sp["n"],
            "R": sp["R"],
            "preferred_direction_rad": sp["preferred_direction_rad"],
            "dominant": dominant_order(sp),
        }
    return out


def axis_agreement(per_region, order=2):
    """Do the regions that show a concentration prefer the SAME axis?

    Returns the resultant of the per-region preferred directions at `order`, weighted
    by each region's own concentration. High agreement means one global axis; low
    agreement with high per-region R means the regions are real but disagree, which is
    the case a pooled test destroys.
    """
    ang, wt = [], []
    for k, v in per_region.items():
        d = v["preferred_direction_rad"].get(order)
        if d is None or not np.isfinite(d):
            continue
        ang.append(float(d) * order)          # back onto the full circle for averaging
        wt.append(float(v["R"].get(order, 0.0)) * math.sqrt(max(v["n"], 1)))
    if not ang:
        return {"n_regions": 0, "agreement_R": float("nan")}
    a = np.asarray(ang)
    w = np.asarray(wt)
    z = np.sum(w * np.exp(1j * a)) / max(np.sum(w), 1e-300)
    return {
        "n_regions": int(a.size),
        "agreement_R": float(np.abs(z)),
        "common_axis_rad": float(np.mod(np.angle(z) / order, 2.0 * np.pi / order)),
        "interpretation": (
            "agreement_R near 1 means the regions prefer a common axis and pooling is "
            "appropriate; near 0 with high per-region R means the regions each have a "
            "real axis and they disagree, and the pooled resultant would be zero for "
            "the wrong reason"),
    }
