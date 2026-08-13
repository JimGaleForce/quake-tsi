"""D-2 -- THE EVENT-PATH VIF MEASUREMENT. SIM ONLY. Priced 0. Nothing here is evidence.

HYPOTHESIS_LEDGER.md §P7-22 Q1, RULED:

  > **Kepler's reading is correct and I make it binding.** F4-58 measured VIF on the
  > **daily-binned count path against a block-bootstrap null**. The event-phase
  > statistic has a different binning (none), a different null (ETAS event times), and
  > a different dependence structure. **Transferring it would be S-18 clause 1
  > exactly.** §K92-0(6)'s table is a **BRACKET**, quotable only as a bracket, until
  > D-2 reports.
  >
  > **Recorded prediction.** §P7-11(a) measured the residual's power excess as 7.6x at
  > 2-7 d, 19.6x at 7-30 d, 64.0x at 365-1000 d -- monotone decreasing toward higher
  > frequency. A semidiurnal phase statistic samples at ~0.5 d, **below the shortest
  > band ever measured**, so I expect **D-2 to return substantially below 24.1**.
  > **S-17 forbids me from quoting that extrapolation as a number** -- which is
  > exactly why D-2 is a measurement and not an inference.

WHAT IS MEASURED, AND THE IDENTITY IT USES
-------------------------------------------
§P7-1(c)'s identity is `VIF = chi2_obs / chi2_ppf(1 - p_surrogate, df)`: the excess of
the observed statistic over the quantile the surrogate tail implies. Kuiper's V is not
a chi2, so the adaptation is stated rather than assumed. V is a CDF excursion and
scales as `n_eff^(-1/2)`, so the SCALED statistic `V*` is the analogue of a
`sqrt(chi2)`: variance inflation by a factor `VIF` multiplies `V*` by `sqrt(VIF)`.
Hence the PRIMARY estimator here is the second-moment ratio

    VIF_band = E[ V*^2 | catalogue ] / E[ V*^2 | theoretical reference ]

which is exactly the brief's *"the Kuiper statistic's dispersion vs its theoretical
reference"*, is dimensionless, is 1 by construction when the catalogue matches the
reference, and needs no tail extrapolation. The SECONDARY estimator is §P7-1(c)'s own
tail form: map each catalogue's `V*` through the reference distribution to a
probability, convert to a 2-df chi2 equivalent, and take
`median(chi2_eq) / median(chi2_2)`. It is reported beside the primary with §P7-10(a)'s
floor CENSORING applied, because a saturated tail probability is a bound and not a
measurement, and the two agreeing is the check that neither is an artefact of its own
mapping.

THE THEORETICAL REFERENCE is `n` i.i.d. uniform phases at the same `n`, computed by
Monte Carlo rather than from Stephens' asymptotic tail: the statistic here is
TWO-SAMPLE against a finite pool, so the asymptotic one-sample tail is the wrong
reference and using it would put a known bias in the denominator.

THE TWO ARMS, AND WHICH ONE IS THE D-2 NUMBER
----------------------------------------------
**ARM A -- PRIMARY, and the D-2 deliverable. True VIF = 1 by construction.** Events are
drawn from a declared smooth intensity by inverse-CDF sampling
(`circstat_event.resample_event_times`), i.e. conditionally independent given lambda --
the same argument `gate_r1._simulate_y` gives for the count path. Any excess over 1
here is INSTRUMENT-SIDE: the intensity's own phase distribution, the finite record, and
the phase map's deterministic warp. This is the number the §K92-0(6) bracket is waiting
on.

**ARM B -- SECONDARY, SENSITIVITY, and clearly not the D-2 number.** A self-exciting
temporal Hawkes/ETAS with genuine sub-daily aftershock clustering, so the statistic
faces the dependence structure a real catalogue has. **SCOPE FLAG:** its branching
ratio is **DECLARED, NOT FITTED**. The frozen fit
(`engine/out/cache/etas_params_k080_0p5deg.json`) is a SPATIO-temporal model whose `K`
is not separable from the spatial kernel's normalisation, so a temporal-only branching
ratio is not identified from it. The Omori SHAPE (`c`, `p`) is taken from the frozen
fit; the branching ratio is declared and swept. Arm B therefore bounds the effect of
clustering at a declared aftershock fraction; it is not a measurement of the real
catalogue's VIF and may not be quoted as one.

DISCIPLINE
----------
Simulation only: no catalogue is read, no region is named, no window is touched, no
holdout hash is spent, and **no EXPLORE_COUNT line is written** -- there is nothing to
spend. Price 0, per §P7-22's table.
"""

from __future__ import annotations

import datetime as _dt

import numpy as np
from scipy import stats

from . import circstat_event as CE
from . import sitetide as ST

# ------------------------------------------------------------ declared sim ----
# All of these are CONSTRUCTION CHOICES, written one way (S-9) and not tuned to any
# outcome. They describe a plausible regional subduction box, not any real one.
SIM_T0 = _dt.datetime(2005, 1, 1)
SIM_N_DAYS = 3650                   # 10 years
SIM_N_EVENTS = 3000                 # a regional-box-scale N
SIM_SITE = {"lat_deg": 38.0, "lon_deg": 142.0, "depth_km": 30.0}   # mid-latitude
SIM_ANNUAL_AMPLITUDE = 0.15         # seasonal modulation of the background rate
SIM_ERA_RAMP = 0.60                 # network-era growth over the record, fractional

# Omori SHAPE from the frozen ETAS fit; branching ratio DECLARED (see module docstring).
FROZEN_OMORI_C_DAYS = 0.0010000983347727406
FROZEN_OMORI_P = 1.193684024163046
FROZEN_OMORI_TRUNC_DAYS = 365.0
DECLARED_BRANCHING_RATIOS = (0.3, 0.5, 0.7)

# The constituent bands the VIF is reported across.
BANDS = ("M2", "S2", "N2", "K1", "O1", "P1", "Mf", "Mm", "Ssa")

# Popper's recorded prediction, carried so the verdict is scored against the
# prediction as recorded rather than against a remembered version of it.
RECORDED_PREDICTION = {
    "source": "HYPOTHESIS_LEDGER.md §P7-22 Q1",
    "claim": ("the event-phase VIF returns SUBSTANTIALLY BELOW 24.1, because the "
              "residual power excess is monotone decreasing toward high frequency "
              "(7.6x at 2-7 d, 19.6x at 7-30 d, 64.0x at 365-1000 d) and a "
              "semidiurnal statistic samples below the shortest band ever measured"),
    "count_path_value_not_transferable": 24.1,
    "s17_note": ("S-17 forbids quoting the extrapolation as a number; this module "
                 "measures instead of extrapolating"),
}

CENSOR_NOTE = (
    "§P7-10(a): a tail probability sitting on its Monte Carlo resolution floor is a "
    "BOUND, not a measurement, and is excluded from the quantile-form VIF and "
    "reported separately with its count. The moment-form estimator does not use a "
    "tail and is not censored.")

SCOPE_FLAGS = (
    "SIM ONLY, price 0. Arm A's catalogues have TRUE VIF = 1 by construction "
    "(conditionally independent draws from a declared smooth intensity), so its "
    "number is the INSTRUMENT-SIDE inflation of the event path. Arm B's branching "
    "ratio is DECLARED, NOT FITTED -- the frozen spatio-temporal K is not separable "
    "from the spatial kernel's normalisation -- so Arm B bounds the effect of "
    "clustering at a declared aftershock fraction and is NOT a measurement of the "
    "real catalogue's VIF.")


# --------------------------------------------------------- declared intensity --
def declared_intensity(n_days=SIM_N_DAYS):
    """A smooth daily intensity: network-era ramp x annual cycle. No sub-daily shape.

    Having NO sub-daily structure is the point: it means every sub-daily statistic
    computed on Arm A is facing a genuine null, so an excess over 1 cannot be the
    intensity smuggling in the very thing being measured.
    """
    d = np.arange(int(n_days), dtype=np.float64)
    era = 1.0 + SIM_ERA_RAMP * d / max(float(n_days) - 1.0, 1.0)
    ann = 1.0 + SIM_ANNUAL_AMPLITUDE * np.cos(2.0 * np.pi * d / 365.2425)
    lam = era * ann
    return np.arange(int(n_days) + 1, dtype=np.float64), lam / lam.mean()


# ------------------------------------------------------------------- arm A ----
def arm_a_catalog(edges, lam, n_events, rng):
    """One TRUE-VIF-1 catalogue: conditionally independent draws from lambda."""
    return CE.resample_event_times(edges, lam, int(n_events), rng)


# ------------------------------------------------------------------- arm B ----
def _omori_delays(n, rng, c=FROZEN_OMORI_C_DAYS, p=FROZEN_OMORI_P,
                  trunc=FROZEN_OMORI_TRUNC_DAYS):
    """Delays from the normalised truncated Omori density, by inverse CDF.

    g(t) ∝ (1 + t/c)^(-p) on [0, trunc]; G(t) = (1 - (1+t/c)^(1-p)) / Z with
    Z = 1 - (1 + trunc/c)^(1-p).
    """
    z = 1.0 - (1.0 + trunc / c) ** (1.0 - p)
    u = rng.random(int(n)) * z
    return c * ((1.0 - u) ** (1.0 / (1.0 - p)) - 1.0)


def arm_b_catalog(edges, lam, n_target, branching_ratio, rng,
                  max_generations=40):
    """One self-exciting catalogue: background + Omori cascade, CONTINUOUS times.

    The cascade gives real sub-daily clustering -- an aftershock minutes after its
    parent -- which is the dependence a semidiurnal event-phase statistic actually
    faces and which Arm A deliberately does not have. Sub-critical by construction
    (`branching_ratio < 1`); the expected total is `n_bg / (1 - branching_ratio)`, so
    the background count is set to hit `n_target` in expectation.
    """
    n = float(branching_ratio)
    if not (0.0 <= n < 1.0):
        raise ValueError("branching ratio must be sub-critical")
    n_bg = max(int(round(float(n_target) * (1.0 - n))), 1)
    t = CE.resample_event_times(edges, lam, n_bg, rng)
    out = [t]
    t_end = float(edges[-1])
    for _ in range(int(max_generations)):
        k = rng.poisson(n, size=t.size)
        tot = int(k.sum())
        if tot == 0:
            break
        parents = np.repeat(t, k)
        child = parents + _omori_delays(tot, rng)
        child = child[child < t_end]
        if child.size == 0:
            break
        out.append(child)
        t = child
    return np.sort(np.concatenate(out))


# ---------------------------------------------------------------- the phases --
def band_phase_fn(period_days):
    """The fixed-frequency clock for one constituent band."""
    def f(times):
        return ST.constituent_phase(SIM_T0, times, period_days)
    return f


def waveform_phase_fn(site=None, grid_minutes=5.0, n_days=SIM_N_DAYS):
    """The D-0 / Tanaka waveform phase at the declared site -- the composite axis.

    This is the axis D-8 would actually run on; the per-band clocks decompose it. The
    site's tidal maxima are computed ONCE for the whole span and shared by every
    catalogue and every null replicate, so the phase map is bit-identical across
    signal and null -- §P7-22(a)'s common-mode requirement, enforced by construction
    rather than by care.
    """
    s = dict(SIM_SITE if site is None else site)
    t_max = ST.tidal_maxima(SIM_T0, -1.0, float(n_days) + 1.0, s["lat_deg"],
                            s["lon_deg"], s["depth_km"], grid_minutes=grid_minutes)

    def f(times):
        return ST.phase_from_maxima(t_max, times)
    return f


# ------------------------------------------------------------- the estimator --
def _vstar(theta, reference):
    return CE.event_kuiper_watson(theta, reference)["V_star"]


def reference_vstar_draws(n_events, n_reference, n_draws, rng):
    """V* for i.i.d. uniform phases at the same n -- the THEORETICAL reference.

    Monte Carlo rather than Stephens' asymptotic tail, because the statistic is
    two-sample against a finite pool and the one-sample asymptotic tail is the wrong
    denominator.
    """
    ref = rng.random(int(n_reference)) * 2.0 * np.pi
    return np.array([_vstar(rng.random(int(n_events)) * 2.0 * np.pi, ref)
                     for _ in range(int(n_draws))], dtype=np.float64)


def vif_from_draws(sim_vstar, ref_vstar):
    """Both estimators, side by side, with §P7-10(a) censoring on the quantile form."""
    s = np.asarray(sim_vstar, dtype=np.float64)
    r = np.asarray(ref_vstar, dtype=np.float64)

    # PRIMARY -- moment form. V* ~ sqrt(chi2)-like, so the ratio of second moments
    # is the variance-inflation analogue and needs no tail.
    vif_moment = float((s ** 2).mean() / (r ** 2).mean())

    # SECONDARY -- §P7-1(c)'s tail form, with the reference distribution supplying the
    # probability and a 2-df chi2 supplying the equivalent statistic.
    rs = np.sort(r)
    b = rs.size
    # mid-rank survival probability of each sim value in the reference sample
    gt = b - np.searchsorted(rs, s, side="right")
    eq = np.searchsorted(rs, s, side="right") - np.searchsorted(rs, s, side="left")
    p_ref = (1.0 + gt + 0.5 * eq) / (1.0 + b)
    floor = 1.0 / (1.0 + b)
    censored = p_ref <= floor + 1e-15
    keep = ~censored
    if keep.any():
        chi2_eq = stats.chi2.isf(p_ref[keep], 2)
        vif_quantile = float(np.median(chi2_eq) / stats.chi2.ppf(0.5, 2))
    else:
        vif_quantile = None
    # Bootstrap CI over CATALOGUES for the primary estimator, so the table carries an
    # uncertainty and a value near 1 is not read as an exact 1.
    br = np.random.default_rng(12345)
    idx = br.integers(0, s.size, size=(2000, s.size))
    boot = (s[idx] ** 2).mean(axis=1) / (r ** 2).mean()
    lo, hi = np.percentile(boot, [2.5, 97.5])

    return {
        "vif_moment": vif_moment,
        "vif_moment_ci95": [float(lo), float(hi)],
        "vif_quantile_p7_1c": vif_quantile,
        "n_catalogs": int(s.size),
        "n_reference_draws": int(r.size),
        "n_censored_at_floor": int(censored.sum()),
        "censor_note": CENSOR_NOTE,
        "mean_Vstar_sim": float(s.mean()),
        "mean_Vstar_reference": float(r.mean()),
        "sd_Vstar_sim": float(s.std(ddof=1)) if s.size > 1 else None,
        "sd_Vstar_reference": float(r.std(ddof=1)) if r.size > 1 else None,
    }


# ------------------------------------------------------------------- the run --
def run(n_catalogs=200, n_reference_draws=400, n_reference_pool=60000,
        n_events=SIM_N_EVENTS, seed=20260813, bands=BANDS,
        branching_ratios=DECLARED_BRANCHING_RATIOS, include_waveform=True,
        verbose=True):
    """Measure the event-path VIF per constituent band, both arms. Returns a record."""
    rng = np.random.default_rng(int(seed))
    edges, lam = declared_intensity()

    axes = [(b, band_phase_fn(ST.CONSTITUENT_PERIODS_DAYS[b])) for b in bands]
    if include_waveform:
        axes.append(("WAVEFORM_D0", waveform_phase_fn()))

    # Arm A catalogues, drawn once and reused across axes so that the per-band numbers
    # are the SAME catalogues seen through different clocks -- which is what makes the
    # band-to-band comparison a decomposition rather than 9 unrelated runs.
    cats_a = [arm_a_catalog(edges, lam, n_events, rng) for _ in range(int(n_catalogs))]
    pool_a = arm_a_catalog(edges, lam, int(n_reference_pool), rng)

    cats_b = {}
    for n in branching_ratios:
        cats_b[n] = [arm_b_catalog(edges, lam, n_events, n, rng)
                     for _ in range(int(n_catalogs))]

    out = {
        "item": "D-2 event-path VIF (SIM ONLY, price 0)",
        "recorded_prediction": RECORDED_PREDICTION,
        "scope_flags": SCOPE_FLAGS,
        "declared_sim": {
            "t0": SIM_T0.isoformat(), "n_days": SIM_N_DAYS,
            "n_events_per_catalog": int(n_events),
            "n_catalogs": int(n_catalogs),
            "site": SIM_SITE,
            "reference_pool": int(n_reference_pool),
            "reference_draws": int(n_reference_draws),
            "omori_shape_from_frozen_fit": {"c_days": FROZEN_OMORI_C_DAYS,
                                            "p": FROZEN_OMORI_P,
                                            "trunc_days": FROZEN_OMORI_TRUNC_DAYS},
            "branching_ratios_declared": list(branching_ratios),
            "seed": int(seed),
        },
        "arm_a_true_vif_1": {},
        "arm_b_declared_clustering": {n: {} for n in branching_ratios},
        "priced_tests": 0,
        "explore_count": ("NOT LOGGED: simulation only, no catalogue, no window, "
                          "no holdout."),
    }

    for name, fn in axes:
        ref_pool = np.asarray(fn(pool_a), dtype=np.float64)
        ref_pool = ref_pool[np.isfinite(ref_pool)]
        ref_draws = reference_vstar_draws(n_events, min(n_reference_pool, 60000),
                                          int(n_reference_draws), rng)

        v_a = []
        for t in cats_a:
            ph = np.asarray(fn(t), dtype=np.float64)
            ph = ph[np.isfinite(ph)]
            v_a.append(_vstar(ph, ref_pool))
        out["arm_a_true_vif_1"][name] = vif_from_draws(v_a, ref_draws)
        if verbose:
            r = out["arm_a_true_vif_1"][name]
            print("  ARM A %-12s VIF_moment = %6.3f  [%.3f, %.3f]  "
                  "VIF_quantile = %s"
                  % (name, r["vif_moment"], r["vif_moment_ci95"][0],
                     r["vif_moment_ci95"][1],
                     "n/a" if r["vif_quantile_p7_1c"] is None
                     else "%6.3f" % r["vif_quantile_p7_1c"]), flush=True)

        for n, cats in cats_b.items():
            v_b = []
            for t in cats:
                ph = np.asarray(fn(t), dtype=np.float64)
                ph = ph[np.isfinite(ph)]
                v_b.append(_vstar(ph, ref_pool))
            out["arm_b_declared_clustering"][n][name] = vif_from_draws(v_b, ref_draws)
            if verbose:
                r = out["arm_b_declared_clustering"][n][name]
                print("  ARM B n=%.1f %-12s VIF_moment = %6.3f"
                      % (n, name, r["vif_moment"]), flush=True)

    worst_a = max(v["vif_moment"] for v in out["arm_a_true_vif_1"].values())
    worst_b = max(v["vif_moment"]
                  for arm in out["arm_b_declared_clustering"].values()
                  for v in arm.values())
    out["verdict"] = {
        "max_vif_arm_a": float(worst_a),
        "max_vif_arm_b_declared": float(worst_b),
        "count_path_value": 24.1,
        "prediction_upheld_arm_a": bool(worst_a < 24.1),
        "prediction_upheld_arm_b_declared": bool(worst_b < 24.1),
        "reading": ("Popper's §P7-22 Q1 prediction is that the event-path VIF returns "
                    "substantially below the count path's 24.1. The D-2 number is "
                    "ARM A. Arm B is a declared-clustering sensitivity and bounds the "
                    "clustering contribution; it is not a measurement of the real "
                    "catalogue."),
    }
    return out


if __name__ == "__main__":          # pragma: no cover - operator entry point
    import json
    import os
    import sys
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    rec = run(verbose=True)
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "out", "audit_event_vif.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(rec, fh, indent=2, default=float)
    print("wrote", path)
    json.dump(rec["verdict"], sys.stdout, indent=2, default=float)
    print()
