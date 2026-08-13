"""D-4 -- KUIPER / WATSON / R1 / 2nd MOMENT ON THE **EVENT PATH**, ETAS event-time null.

The count-path versions live in `engine/circstat.py` and are unchanged. This module is
the port §K92-1 D-4 prices at 0, and it lands under three binding constraints, all of
which are implemented here rather than described:

1. **§P7-22(a), the common-mode requirement -- the design's load-bearing part.**
   The null is **ETAS-simulated EVENT TIMES pushed through the IDENTICAL phase
   computation**, never a permutation of the observed phases. Popper:

     > *"A deterministic warp of a uniform circular distribution is generally not
     > uniform, so this is a live hazard. But Kepler's declared null is ETAS-simulated
     > event times pushed through the same body-tide phase computation, so any
     > deterministic non-uniformity of that mapping is COMMON-MODE to signal and null
     > and cancels exactly. ... His explicit refusal to use a phase permutation is
     > what saves the design."*

   `event_omnibus` therefore takes NULL EVENT TIMES, not null phases, and calls the
   caller-supplied `phase_fn` on them itself. A phase-permutation null is not merely
   discouraged here -- there is no argument through which one can be passed.

2. **§P7-21(c), inherited and BINDING** (§P7-22 ratification 4: *"D-4 may not land
   without inheriting it"*). Two halves, both carried:
   * **No circular-shift null for a periodic phase feature.** A tidal phase IS
     periodic, so the shift null is provably powerless and is NOT COMPUTED -- absent,
     with `circstat.PERIODIC_SHIFT_OMITTED` on the row as its reason. `p_raw` is the
     resampling p alone.
   * **Mid-rank, tie-tolerant p** (`circstat.tie_tolerant_p`) against the null
     ensemble, so float jitter in the observed statistic moves p by O(1/B) instead of
     by a whole plateau.
   Both are USED FROM `circstat`, not re-implemented: a copy could drift, and §P7-21's
   defect was exactly a statistic that looked right and was not.

3. **The statistic is TWO-SAMPLE.** On the count path the expected CDF is the ETAS
   offset's own phase distribution. On the event path there is no offset vector, so
   the reference CDF is the EMPIRICAL phase distribution of a large pool of ETAS-
   simulated events -- the same construction, estimated instead of computed. The core
   is `circstat._kuiper_watson_core` verbatim, which is what carries the phase-tie
   grouping, the CDF noise floor and the degenerate-case guarantee across unchanged.

   Two-sample n-scaling uses `n_eff = n1 n2 / (n1 + n2)`, the standard two-sample
   effective size; with a reference pool much larger than the sample this tends to
   n1, i.e. to the one-sample scaling, which is the sense in which the port is a port.

WHAT THIS MODULE DOES NOT DO
-----------------------------
It does not report psi. §P7-22(a) keeps the concentration ANGLE unreportable under
K-090(c) for a body-tide-only instrument in an ocean-loaded setting; `event_moments`
returns the mean angles because the second moment's axial/harmonic reading needs both
R1 and R2 side by side (F9-01's Pit), and every record carries `PSI_UNREPORTABLE`
saying the angle may not be quoted as a measurement.

It also does not decide anything about the Earth. Nothing here is evidence.
"""

from __future__ import annotations

import numpy as np

from . import circstat as C

# The primary statistic, declared per §P7-22 Q4 ("Declare which is PRIMARY").
PRIMARY_STATISTIC = "kuiper_V"

PSI_UNREPORTABLE = (
    "§P7-22(a): K-090(c) is SCOPED, not overturned. A body-tide-only instrument in "
    "an ocean-loaded setting is an approximately unbiased detector of CONCENTRATION "
    "and a BIASED estimator of the angle. R1/R2 mean angles are returned for the "
    "F9-01 axial-vs-harmonic reading ONLY; psi may not be quoted as a measurement "
    "until D-10/K-093 bounds the loading warp.")

FAIL_IS_NOT_A_BOUND = (
    "§P7-22(a): a PASS is informative (the multi-constituent loading warp is "
    "common-mode with the null and cannot manufacture concentration). A FAIL is NOT "
    "a bound and may not be quoted as one -- the warp can DESTROY real concentration "
    "by an unbounded factor, and D-10/K-093 is the gate on quoting a null.")

COMMON_MODE_NOTE = (
    "§P7-22(a): the null is ETAS-simulated EVENT TIMES pushed through the IDENTICAL "
    "phase computation (S-8's 'through the identical code path'). A phase "
    "PERMUTATION null is refused by construction -- it would destroy the clustering "
    "that carries the entire VIF problem, and it would not cancel the loading warp.")


# ------------------------------------------------------- two-sample statistics --
def _two_sample_weights(theta_obs, theta_ref):
    """Union-sorted (obs, exp) weight vectors and the phase-tie group ends.

    Observed events contribute weight 1 to `obs` and 0 to `exp`; reference events the
    reverse. `circstat._kuiper_watson_core` then forms the two empirical CDFs and
    their difference at DISTINCT phases only -- the tie handling that
    `circstat.phase_group_ends` exists for rides along unchanged.
    """
    a = C.wrap_phase(np.asarray(theta_obs, dtype=np.float64).ravel())
    b = C.wrap_phase(np.asarray(theta_ref, dtype=np.float64).ravel())
    th = np.concatenate([a, b])
    w_obs = np.concatenate([np.ones(a.size), np.zeros(b.size)])
    w_ref = np.concatenate([np.zeros(a.size), np.ones(b.size)])
    o = np.argsort(th, kind="stable")
    ths = th[o]
    return ths, w_obs[o], w_ref[o], C.phase_group_ends(ths), a.size, b.size


def event_kuiper_watson(theta_obs, theta_ref):
    """Two-sample Kuiper V and Watson U^2 for EVENT phases against a reference pool.

    Returns the raw and n-scaled statistics. The scaling is
    `circstat._scaled` at `n_eff = n1 n2 / (n1 + n2)`, so V* and U2* are comparable
    across catalogues of different size -- which the resampling null requires, because
    an ETAS-simulated catalogue does not have exactly the observed event count.
    """
    _th, w_obs, w_ref, ends, n1, n2 = _two_sample_weights(theta_obs, theta_ref)
    if n1 == 0 or n2 == 0:
        raise ValueError("both the event set and the reference pool must be non-empty")
    V, U2 = C._kuiper_watson_core(w_obs, w_ref, ends)
    n_eff = float(n1) * float(n2) / float(n1 + n2)
    Vs, U2s = C._scaled(V, U2, n_eff)
    return {"V": float(V), "U2": float(U2), "V_star": float(Vs),
            "U2_star": float(U2s), "n_events": int(n1), "n_reference": int(n2),
            "n_eff": n_eff, "n_distinct_phases": int(ends.size)}


def event_moments(theta):
    """R1, R2 and their mean angles for unweighted event phases (F9-01's pair)."""
    out = C.moments_report(np.asarray(theta, dtype=np.float64))
    out["psi_unreportable"] = PSI_UNREPORTABLE
    return out


def event_second_moment_stat(theta_obs, theta_ref):
    """A scalar second-moment contrast: |R2_obs - R2_ref| on the doubled angle.

    Priced beside Kuiper/Watson/R1 per §P7-22 Q4 ("keep all four priced"), and
    referenced to the pool rather than to zero so that any deterministic
    second-harmonic structure the phase MAP itself imposes is common-mode and cancels
    -- the same argument §P7-22(a) makes for the whole design.
    """
    r2o, _ = C.circular_moment(np.asarray(theta_obs, dtype=np.float64), k=2)
    r2r, _ = C.circular_moment(np.asarray(theta_ref, dtype=np.float64), k=2)
    return float(abs(r2o - r2r))


def event_r1_stat(theta_obs, theta_ref):
    """|R1_obs - R1_ref|, the first-resultant contrast, same referencing argument."""
    r1o, _ = C.circular_moment(np.asarray(theta_obs, dtype=np.float64), k=1)
    r1r, _ = C.circular_moment(np.asarray(theta_ref, dtype=np.float64), k=1)
    return float(abs(r1o - r1r))


# ------------------------------------------------ the ETAS event-time resampler --
def resample_event_times(day_edges, lam, n_events, rng):
    """Draw `n_events` CONTINUOUS event times from a piecewise-constant intensity.

    `lam[i]` is the intensity on `[day_edges[i], day_edges[i+1])`. Times are drawn by
    inverse-CDF sampling of the normalised intensity, so the draw is an inhomogeneous
    Poisson process CONDITIONED ON N -- conditionally independent given lambda, which
    is what makes a catalogue built this way a TRUE NULL with true VIF exactly 1
    (`engine/gate_r1.py::_simulate_y` gives the same argument for the count path).

    Conditioning on N is deliberate: the statistic is a CDF comparison, so holding the
    sample size fixed removes an N-fluctuation that would otherwise enter the null's
    dispersion and be mistaken for phase structure.
    """
    e = np.asarray(day_edges, dtype=np.float64)
    l = np.asarray(lam, dtype=np.float64)
    if e.size != l.size + 1:
        raise ValueError("day_edges must have one more element than lam")
    w = l * np.diff(e)
    tot = float(w.sum())
    if not (tot > 0):
        raise ValueError("intensity integrates to zero")
    cw = np.concatenate([[0.0], np.cumsum(w)]) / tot
    u = rng.random(int(n_events))
    k = np.clip(np.searchsorted(cw, u, side="right") - 1, 0, l.size - 1)
    frac = (u - cw[k]) / np.maximum(cw[k + 1] - cw[k], 1e-300)
    return np.sort(e[k] + frac * (e[k + 1] - e[k]))


def thin_by_phase_intensity(times, phases, multiplier, n_keep, rng):
    """Rejection-thin an event set so its phases carry a declared modulation.

    The recovery-harness primitive for the event path: draw a large candidate set from
    the null, then ACCEPT each candidate with probability `multiplier(phase) / max`.
    The survivors are a sample from the null intensity TIMES the declared phase
    response -- a genuine rate modulation, not a phase relabelling -- which is the
    only planting construction a rotation-invariant statistic can honestly be tested
    against (§P7-8(d): a plant must be a signal, not a re-drawing of the answer).

    Returns the kept TIMES; the caller re-derives phases through the same map so that
    the planted arm and the null differ in nothing but the acceptance step.
    """
    t = np.asarray(times, dtype=np.float64)
    m = np.asarray(multiplier(np.asarray(phases, dtype=np.float64)), dtype=np.float64)
    keep = rng.random(t.size) < (m / m.max())
    kept = t[keep]
    if kept.size < int(n_keep):
        raise ValueError("candidate pool too small after thinning: %d < %d"
                         % (kept.size, int(n_keep)))
    return np.sort(kept[:int(n_keep)])


# ---------------------------------------------------------------- the omnibus --
def event_omnibus(theta_obs, null_times, phase_fn, reference_times=None,
                  reference_phases=None):
    """The four event-path statistics with p-values from an ETAS EVENT-TIME null.

    Parameters
    ----------
    theta_obs : the observed events' phases (already through `phase_fn`).
    null_times : sequence of B arrays of simulated event TIMES -- not phases. Each is
        pushed through `phase_fn` here, which is what makes the null common-mode with
        the signal (§P7-22(a)). Passing phases is impossible by signature.
    phase_fn : callable(times) -> phases. The IDENTICAL computation used for
        `theta_obs`. In Tranche D this is `sitetide.tanaka_phase` (D-0/D-3) or
        `sitetide.constituent_phase`, bound to one site.
    reference_times / reference_phases : the pool that supplies the reference CDF.
        Supply exactly one. It must be DISJOINT from `null_times`, or every null
        replicate is compared against a pool containing itself and the null's
        dispersion is biased low.

    Returns a list of rows in `circstat.omnibus_test`'s shape, with
    `p_circular_shift = None` and `PERIODIC_SHIFT_OMITTED` as its reason (§P7-21(c)).
    """
    if (reference_times is None) == (reference_phases is None):
        raise ValueError("supply exactly one of reference_times / reference_phases")
    ref = (np.asarray(reference_phases, dtype=np.float64) if reference_times is None
           else np.asarray(phase_fn(reference_times), dtype=np.float64))
    ref = ref[np.isfinite(ref)]

    obs = np.asarray(theta_obs, dtype=np.float64)
    obs = obs[np.isfinite(obs)]

    kw = event_kuiper_watson(obs, ref)
    stats_obs = {
        "kuiper_V": kw["V_star"],
        "watson_U2": kw["U2_star"],
        "R1_contrast": event_r1_stat(obs, ref),
        "second_moment_contrast": event_second_moment_stat(obs, ref),
    }

    draws = {k: [] for k in stats_obs}
    for t in null_times:
        ph = np.asarray(phase_fn(t), dtype=np.float64)
        ph = ph[np.isfinite(ph)]
        if ph.size == 0:
            continue
        k2 = event_kuiper_watson(ph, ref)
        draws["kuiper_V"].append(k2["V_star"])
        draws["watson_U2"].append(k2["U2_star"])
        draws["R1_contrast"].append(event_r1_stat(ph, ref))
        draws["second_moment_contrast"].append(event_second_moment_stat(ph, ref))

    b = len(draws[PRIMARY_STATISTIC])
    if b == 0:
        raise ValueError("the ETAS event-time null produced no usable replicates")

    mom = event_moments(obs)
    rows = []
    for name, stat in stats_obs.items():
        d = np.asarray(draws[name], dtype=np.float64)
        p = C.tie_tolerant_p(d, stat, d.size)
        rows.append({
            "test": name,
            "path": "event",
            "primary": bool(name == PRIMARY_STATISTIC),
            "df": None,
            "statistic": float(stat),
            "V": kw["V"], "U2": kw["U2"],
            "V_star": kw["V_star"], "U2_star": kw["U2_star"],
            "R1": mom["R1"], "R2": mom["R2"],
            "n_events": int(kw["n_events"]),
            "n_reference": int(kw["n_reference"]),
            "p_circular_shift": None,
            "p_circular_shift_note": C.PERIODIC_SHIFT_OMITTED,
            "p_etas_event_null": float(p),
            "p_raw": float(p),
            "p_method_note": "mid-rank, tie-tolerant (§P7-21)",
            "n_surrogates": int(d.size),
            "p_resolution_floor": float(1.0 / (d.size + 1.0)),
            "null": "ETAS event-time resampling through the identical phase map",
            "common_mode_note": COMMON_MODE_NOTE,
            "fail_is_not_a_bound": FAIL_IS_NOT_A_BOUND,
            "psi_unreportable": PSI_UNREPORTABLE,
            "band_note": ("§P7-3(2)'s day-binning sinc note does NOT apply on the "
                          "event path: these phases are computed at each event's own "
                          "sub-daily time. The sub-daily arm's own gate is "
                          "observer.assert_subdaily_gate (§P7-3(3), D-5)."),
        })
    return rows


def ks_uniform_one_sided(ps, a=0.01):
    """§P6-9(a) one-sided KS, imported from `gate_r1` so there is exactly one copy."""
    from . import gate_r1
    return gate_r1.ks_uniform_one_sided(ps, a=a)
