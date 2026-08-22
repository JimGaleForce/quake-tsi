"""B2 / §S3 -- THE REGION-LATTICE RUNNER: the coincidence scan, and its accountant.

> **THE SENTENCE (SEARCHER.md §S0, and it is the whole discipline).** The SEARCHER
> never produces evidence. It produces a RANKING and a CONTROL-CALIBRATED SURVIVOR
> COUNT, and its only licensed output is a mechanically generated FREEZE FILE. Every
> claim this program ever makes from a SEARCHER hit is made by the frozen test that
> the freeze file specifies, run afterwards, on data the scan did not touch.

This module iterates the declared (region x property x magnitude-stratum) lattice,
computes ONE declared statistic per property type, runs the matched control arm at the
IDENTICAL threshold, and emits `(n_survivors_real, n_survivors_control)` plus a ranked
list capped at K. It is priced **EXPLORATORY-UNPRICED** under Rule 4.4 / §P7-4 / SP-1.

WHAT IS FROZEN BEFORE THE SCAN, AND THEREFORE CANNOT BE TUNED AFTER IT
-----------------------------------------------------------------------
`declare_scan` builds the declaration and hashes it. Everything below is inside the
hash; **`--jobs` is not**, on the same rule the miner already follows (execution, not
statistics), and neither is any chunk size.

  * **`m`, the FULL declared cell count** -- SP-6.1: *"Understating `m` is the way to
    cheat this protocol, and it is the only way, so it is the thing that gets an
    assertion."* `assert_cell_count` is that assertion, in the same shape as
    `strata.assert_partition_total`.
  * **`alpha = q / m` with `q = 0.10`** -- SP-3, and **the OR-limb is REFUSED**. There
    is no second, more lenient limb in this module and adding one is a new declaration.
    Jim's own 14-of-14 Monday clears the strict bar by four orders of magnitude, which
    is the demonstration that strong claims do not need a discount.
  * the **declared concentration forms**, the **n >= 8 floor** (§S6.2's cliff: 14-of-14
    survives m = 336,000 at 4.95e-07 while **10-of-14 fails at 0.683**, so the floor is
    load-bearing rather than decorative), the **magnitude strata** (SP-6.7), the
    **property family list**, and the **lattice digest**.

THE NULL IS ALWAYS ETAS-SIMULATED EVENT TIMES THROUGH THE IDENTICAL PROPERTY MAP
---------------------------------------------------------------------------------
Never a permutation of the property. §P7-22(a): a deterministic warp of the property
map (ocean loading, admittance, ephemeris approximation) is common-mode to signal and
null and **cancels exactly** only if the null traverses the same code. `scan_cell`
takes NULL TIMES and calls the column's own `evaluate` on them; a permuted property
cannot be passed through this signature. **Full ETAS with triggering, not
background-only** (§P7-23(A) condition 3) is the caller's obligation and is recorded on
every row as a declared field, because this module cannot verify the provenance of an
array of times and says so rather than implying otherwise.

THE CONTROL ARM (§S3.3(1), F7-b)
---------------------------------
The control arm is a matched scan over **synthetic properties with identical marginal
structure**: each real column's time axis is warped by an **F8-15-STYLE random monotone
clock with a DECLARED coefficient of variation** (`CONTROL_CV`, knots every
`CONTROL_KNOT_DAYS`). A monotone warp preserves the property's marginal/dwell structure
and destroys its alignment with the events, which is exactly what a matched control
needs to be. **It is scanned identically and its survivors are counted in the same
multiplicity at the same alpha** -- Kepler's own F9-20 correction to Popper, ratified
at §P7-6. *A control channel that faces a different threshold than the science channel
is not a control.* **A raw survivor count with no control arm may not be quoted, in
any forum, ever** (§S3.3(1)).

> **READ `CONTROL_DEGENERACY_FOUND` BEFORE CHANGING THE CONTROL.** The first version of
> this arm seeded `clocks.random_monotone_clock` from the scan's IDENTITY time axis.
> F8-15 matches its reference clock's CV, the identity clock's CV is ZERO, and both
> F8-15 modes therefore return the identity map -- so the "control" was BIT-IDENTICAL to
> the science arm and the first SP-7 gate run scored it as such. The defect was found by
> `engine/tests/test_searcher.py`, is recorded rather than quietly repaired, and is the
> reason the warp's lumpiness is a DECLARED constant instead of one matched from data.

WHAT LEAVES THE RUN
--------------------
`(n_survivors_real, n_survivors_control)`; the ranked list capped at K = 30; the S-8
sim-calibrated max-statistic, **reported whether or not anything promotes** (SP-6.4);
and one `EXPLORE_COUNT.jsonl` line carrying the FULL declared cell count as
`n_declared_tests` even though the scan is unpriced (SP-1.4: *unpriced never means
uncounted*). **Zero stubs, zero K-entries, zero findings** (§S3.3(3)).

Nothing in this module is evidence.
"""

from __future__ import annotations

import datetime as _dt
import math
import os

import numpy as np

from . import (circstat as C, circstat_event as CE, clocks as CL,
               lattice_s1 as LAT, observer as OBS, properties as P, splits,
               strata as ST)

SEARCHER_RULE_ID = "S3-searcher-v1"

GENERATOR_NOT_EVIDENCE = (
    "THE SEARCHER IS A GENERATOR, NOT EVIDENCE. Nothing it prints is a finding, a "
    "claim, or a K-entry. Every row below is EXPLORATORY-UNPRICED under Rule 4.4 / "
    "§P7-4 / §P7-24 SP-1.")

THE_SENTENCE = (
    "SEARCHER.md §S0: the SEARCHER never produces evidence. It produces a RANKING and "
    "a CONTROL-CALIBRATED SURVIVOR COUNT, and its only licensed output is a "
    "mechanically generated FREEZE FILE. Every claim this program ever makes from a "
    "SEARCHER hit is made by the frozen test that the freeze file specifies, run "
    "afterwards, on data the scan did not touch.")

CONTROL_MANDATORY = (
    "SEARCHER.md §S3.3(1): only the DIFFERENCE between the real and control survivor "
    "counts is interpretable. A raw survivor count with no control arm may not be "
    "quoted, in any forum, ever.")

OR_LIMB_REFUSED = (
    "§P7-24 SP-3: the promotion rule is ONE rule. alpha = q/m with q = 0.10 and m the "
    "scan's FULL declared cell count. The proposed OR-limb ('or the strong-claim "
    "binomial clears 1e-6 at n >= 10') is REFUSED: a flat 1e-6 limb is 10x more "
    "lenient than q/m at m = 1e6 and 1,000x at 1e8, so it would silently become the "
    "operative rule for every large scan and decouple the bar from the search size. "
    "An OR of two rules is two chances to promote. Honouring small-n does not require "
    "lowering the bar; it requires not confusing 'few events' with 'weak effect.'")

ANALYTIC_IS_NOT_THE_HEADLINE = (
    "SEARCHER.md §S6.3 / §P7-23(A): the reported p is ALWAYS the full-ETAS-null p. "
    "The analytic tail printed beside it is a SANITY CHECK and is never the headline "
    "-- it assumes an independence and a uniformity that the null measures directly, "
    "and fourteen events in one region across a century are not fourteen independent "
    "draws if any of them are aftershocks of each other.")

# ------------------------------------------- the constants, frozen before any run --
Q_DEFAULT = 0.10                     # SP-3
K_CAP = 30                           # §S3.3(2), §S4.1
N_EVENTS_FLOOR = 8                   # §S4.1 / §S6.2's cliff
DECLARED_FORMS = ("quadrant", "single-cell", "arc")
HUMAN_SCHEDULE_MC_FLOOR = P.HUMAN_SCHEDULE_MC_FLOOR      # 6.0, F7-d
RANKING_FUNCTION = "-log10(p_real), descending, within magnitude stratum"


class CellCountMismatch(AssertionError):
    """The number of cells evaluated differs from the declared `m` (SP-6.1)."""


class ControlArmMissing(AssertionError):
    """A scan was run without its matched control arm (§S3.3(1), F7-b)."""


class HoldoutTouched(AssertionError):
    """A scan was handed events outside the exploration window (SP-6.2)."""


def alpha_from_m(m, q=Q_DEFAULT):
    """SP-3's threshold. One rule, no limbs."""
    m = int(m)
    if m <= 0:
        raise ValueError("m must be positive; it is the scan's FULL declared cell count")
    return float(q) / float(m)


def assert_cell_count(n_evaluated, m_declared):
    """SP-6.1, in the shape of `strata.assert_partition_total`. The anti-cheat."""
    if int(n_evaluated) != int(m_declared):
        raise CellCountMismatch(
            "the scan evaluated %d cells and declared m = %d. §P7-24 SP-6.1: the "
            "scan's full cell count is fixed and hashed BEFORE the scan runs and the "
            "searcher ASSERTS that the number of cells evaluated equals m. "
            "Understating m is the way to cheat this protocol, and it is the only "
            "way, so it is the thing that gets an assertion."
            % (int(n_evaluated), int(m_declared)))
    return int(m_declared)


# ------------------------------------------------------------ the declaration --
def declare_scan(lattice_regions, property_names, mag_strata, families,
                 q=Q_DEFAULT, k_cap=K_CAP, n_floor=N_EVENTS_FLOOR,
                 forms=DECLARED_FORMS, l3_cells=None, scan_id=None,
                 window=None, extra=None):
    """The frozen, hash-affecting scan declaration. Built BEFORE any statistic runs.

    `m` is `n_regions x n_properties x n_mag_strata` -- the FULL lattice, including
    every cell that will turn out to be too small to score. A cell that is skipped for
    lack of events is still a declared cell and still counts in `m`; dropping it would
    be understating `m`, which SP-6.1 names as the only way to cheat the protocol.
    """
    regions = [r["name"] for r in lattice_regions]
    props = [str(p) for p in property_names]
    strata_mc = [float(s) for s in mag_strata]
    m = len(regions) * len(props) * len(strata_mc)
    decl = {
        "searcher_rule_id": SEARCHER_RULE_ID,
        "scan_id": str(scan_id or _dt.datetime.now(_dt.timezone.utc)
                       .strftime("scan_%Y%m%dT%H%M%S")),
        "lattice_rule_id": LAT.LATTICE_RULE_ID,
        "lattice_digest": LAT.declaration_digest(l3_cells),
        "property_rule_id": P.PROPERTY_RULE_ID,
        "regions": regions,
        "properties": props,
        "mag_strata": strata_mc,
        "mag_strata_rule": P.MAG_STRATA_RULE,
        "families": list(families),
        "m": int(m),
        "q": float(q),
        "alpha": alpha_from_m(m, q),
        "k_cap": int(k_cap),
        "n_events_floor": int(n_floor),
        "declared_forms": list(forms),
        "ranking_function": RANKING_FUNCTION,
        "human_schedule_mc_floor": float(HUMAN_SCHEDULE_MC_FLOOR),
        "window": (list(window) if window is not None else None),
        "or_limb": OR_LIMB_REFUSED,
        "the_sentence": THE_SENTENCE,
        "control_mandatory": CONTROL_MANDATORY,
    }
    if extra:
        decl.update(extra)
    decl["config_hash"] = splits.config_hash(decl)
    return decl


# ---------------------------------------------------------- the statistics ----
# §S3.1: exactly ONE primary statistic per property TYPE, fixed for all cells, no
# alternatives run (S-9). phase and level share the circular path because a level is
# consumed as its DWELL-CORRECTED PIT mapped to an angle -- §S3.1's "Kuiper on the
# PIT-transformed level", literally, which is why the arcsine correction is a
# transform rather than a caveat.
PRIMARY_STATISTIC = {
    "phase": "kuiper_V_event_path",
    "level": "kuiper_V_on_dwell_corrected_PIT",
    "categorical": "max_cell_binomial_tail_vs_empirical_base_rates",
}


def _circular_stat(theta_obs, theta_ref):
    """Kuiper V*, two-sample, event path. `circstat_event` is used, never re-implemented."""
    return float(CE.event_kuiper_watson(theta_obs, theta_ref)["V_star"])


def _base_rates(cats_ref, categories):
    """The catalogue's OWN empirical base rates, from the null pool. Never 1/k.

    §S3.1: *"scored against the catalogue's own empirical base rates, not against 1/k
    -- the difference between 1/7 and the observed weekday rate costs 4.5 orders of
    magnitude and it is not optional."* §S6.2(1) works the arithmetic: an inflated
    Monday rate of 0.186 moves P(all 14) from 1.47e-12 to 5.93e-11. **The correction
    is 45% of the evidence on a logarithmic scale.**
    """
    ref = np.asarray(cats_ref, dtype=np.float64)
    n = float(ref.size)
    rates = np.array([max(float(np.count_nonzero(ref == c)) / n, 1.0 / (n + 1.0))
                      for c in categories], dtype=np.float64)
    return rates / rates.sum()


def _binom_sf(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p), by the regularised incomplete beta.

    Exact tail, not a normal approximation: the whole small-n instrument depends on
    the tail being right at k = n, where a normal approximation is worthless.
    """
    k, n = int(k), int(n)
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    from scipy import special
    return float(special.betainc(k, n - k + 1, float(p)))


def _categorical_stat(cats_obs, base_rates, categories):
    """The pre-declared max-cell binomial tail, as -log10(min tail). §S3.1."""
    obs = np.asarray(cats_obs, dtype=np.float64)
    n = int(obs.size)
    if n == 0:
        return 0.0, None, 1.0
    best_p, best_c, best_k = 1.0, None, 0
    for c, pr in zip(categories, base_rates):
        k = int(np.count_nonzero(obs == c))
        tail = _binom_sf(k, n, float(pr))
        if tail < best_p:
            best_p, best_c, best_k = tail, c, k
    return -math.log10(max(best_p, 1e-300)), best_c, best_p


CONTROL_CV = 1.0
CONTROL_KNOT_DAYS = 90.0

CONTROL_DEGENERACY_FOUND = (
    "A DEFECT FOUND BY ITS OWN TEST, RECORDED RATHER THAN QUIETLY FIXED. The first "
    "version of this control warped the scan's time axis with "
    "`clocks.random_monotone_clock` seeded from the IDENTITY clock. F8-15 matches its "
    "reference clock's mean increment and CV -- and the identity clock's increments "
    "are all equal, so its CV is ZERO and BOTH F8-15 modes return the identity: "
    "`block` permutes identical blocks, and `iid` draws Gamma with shape 1/CV^2 -> "
    "infinity, i.e. a constant. The 'control' was therefore BIT-IDENTICAL to the "
    "science arm, and the first SP-7 gate run scored it as such. The searcher's "
    "control is not a warp of a clock derived FROM the data; it is a warp with a "
    "DECLARED lumpiness, and CONTROL_CV is that declaration.")


def _control_column(col, rng, mode="declared-cv", cv=CONTROL_CV,
                    knot_days=CONTROL_KNOT_DAYS, span_days=None):
    """A matched synthetic property: the SAME map on a random monotone time warp.

    §S3.3(1)'s *"synthetic properties with identical marginal structure"*, built as an
    **F8-15-style random monotone reparameterisation of the scan's time axis with a
    DECLARED coefficient of variation** (`CONTROL_CV = 1.0`, knots every
    `CONTROL_KNOT_DAYS = 90` days). Read `CONTROL_DEGENERACY_FOUND` before changing
    any of that: seeding the warp from the identity clock makes it the identity map,
    which is how a control can be a duplicate of the thing it controls.

    WHY THIS PRESERVES MARGINAL STRUCTURE AND DESTROYS ALIGNMENT, stated rather than
    assumed. The warp is piecewise-linear with knots 90 days apart, so within a knot
    interval it is a pure RATE CHANGE. Every declared property here is quasi-periodic
    with a period far shorter than 90 days (12.4 h to 29.5 d), so each interval still
    sweeps many complete cycles and the property's occupancy measure over the span is
    essentially unchanged -- while the phase an event lands on is completely
    re-drawn. **`monotone` matters**: a non-monotone scramble would not be a clock,
    and F8-15's whole point is that the control is a CLOCK.

    The construction is validated against `clocks.clock_summary` on the way out, so
    the claim that it is strictly monotone with the declared CV is checked and not
    asserted.
    """
    span = float(span_days if span_days is not None else 3650.0)
    n_knots = max(8, int(round(span / float(knot_days))))
    ref = np.linspace(0.0, span, n_knots + 1)
    # Gamma increments with the DECLARED CV: shape = 1/cv^2, then endpoint-matched so
    # the warp maps the span onto itself exactly (F8-15's endpoint identity).
    shape = 1.0 / max(float(cv) ** 2, 1e-12)
    inc = rng.gamma(shape, 1.0 / shape, size=n_knots)
    inc = inc * (span / float(inc.sum()))
    tau = np.concatenate([[0.0], np.cumsum(inc)])
    CL.clock_increments(tau)                 # raises unless strictly monotone

    def warped(day_float):
        d = np.asarray(day_float, dtype=np.float64)
        return col.evaluate(np.interp(d, ref, tau))

    ctrl = P.PropertyColumn(
        col.name + "__F8-15-control", col.family, col.ptype, col.pclass, warped,
        dict(col.provenance), dwell=col.dwell, categories=col.categories,
        subdaily=col.subdaily, attached=set(col.attached),
        notes=list(col.notes) + [
            "F8-15-STYLE RANDOM-CLOCK CONTROL (%s, CV = %.2f, knots every %.0f d): "
            "the same property map on a random monotone warp of the time axis. "
            "Marginal structure preserved, alignment destroyed. Scanned at the "
            "IDENTICAL alpha and counted in the SAME multiplicity (F7-b). %s | %s"
            % (mode, float(cv), float(knot_days), CL.MANDATORY_NOTE,
               CONTROL_DEGENERACY_FOUND)])
    return ctrl


def concentration_form(col):
    """Which of the DECLARED forms this column's statistic can express (§S4.1)."""
    if col.ptype == "categorical":
        if col.categories is not None and len(col.categories) == 4 and \
                col.name.endswith("quadrant"):
            return "quadrant"
        return "single-cell"
    return "arc"


# ------------------------------------------------------------- the cell scan --
def scan_cell(col, obs_times, null_times, ref_times, region, mc, decl,
              observer_report=None, null_provenance="UNDECLARED"):
    """One lattice cell -> one row. No promotion decision is taken here.

    `null_times` is a sequence of B arrays of SIMULATED EVENT TIMES -- not property
    values. They are pushed through `col.evaluate` inside this function, which is what
    makes the null common-mode with the signal (§P7-22(a)) and is why a permuted
    property cannot be handed to this signature.
    """
    obs = np.asarray(obs_times, dtype=np.float64)
    ref = np.asarray(ref_times, dtype=np.float64)
    n_obs_declared = int(obs.size)

    if col.ptype == "categorical":
        v_obs = np.asarray(col.evaluate(obs), dtype=np.float64)
        v_ref = np.asarray(col.evaluate(ref), dtype=np.float64)
        ok = np.isfinite(v_obs)
        v_obs = v_obs[ok]
        v_ref = v_ref[np.isfinite(v_ref)]
        cats = col.categories or tuple(sorted(set(v_ref.tolist())))
        rates = _base_rates(v_ref, cats)
        stat, cell, tail = _categorical_stat(v_obs, rates, cats)
        draws = []
        for t in null_times:
            vn = np.asarray(col.evaluate(t), dtype=np.float64)
            vn = vn[np.isfinite(vn)]
            if vn.size:
                draws.append(_categorical_stat(vn, rates, cats)[0])
        analytic = {"p_analytic_binomial_max_cell": float(tail),
                    "p_analytic_times_k": float(min(1.0, tail * len(cats))),
                    "max_cell": (None if cell is None else float(cell)),
                    "base_rates": [float(r) for r in rates],
                    "note": ANALYTIC_IS_NOT_THE_HEADLINE}
    else:
        v_obs = np.asarray(col.statistic_values(obs), dtype=np.float64)
        v_ref = np.asarray(col.statistic_values(ref), dtype=np.float64)
        v_obs = v_obs[np.isfinite(v_obs)]
        v_ref = v_ref[np.isfinite(v_ref)]
        stat = _circular_stat(v_obs, v_ref)
        draws = []
        for t in null_times:
            vn = np.asarray(col.statistic_values(t), dtype=np.float64)
            vn = vn[np.isfinite(vn)]
            if vn.size:
                draws.append(_circular_stat(vn, v_ref))
        analytic = {"note": ANALYTIC_IS_NOT_THE_HEADLINE}

    d = np.asarray(draws, dtype=np.float64)
    if d.size == 0:
        raise ValueError("the ETAS event-time null produced no usable replicates for "
                         "cell (%s, %s, Mc=%s)" % (region["name"], col.name, mc))
    p = float(C.tie_tolerant_p(d, float(stat), d.size))

    ok_layer, why_layer = col.may_promote()
    row = {
        "region": region["name"],
        "layer": region.get("layer"),
        "region_class": region.get("class"),
        "property": col.name,
        "family": col.family,
        "type": col.ptype,
        "property_class": col.pclass,
        "mc": float(mc),
        "statistic": PRIMARY_STATISTIC[col.ptype],
        "value": float(stat),
        "p_real": p,
        "n_events": int(v_obs.size),
        "n_events_declared": n_obs_declared,
        "n_reference": int(v_ref.size),
        "n_surrogates": int(d.size),
        "p_resolution_floor": float(1.0 / (d.size + 1.0)),
        "null_draws": d,                      # kept for the S-8 max-statistic matrix
        "null": "ETAS-simulated EVENT TIMES through the IDENTICAL property map",
        "null_provenance": str(null_provenance),
        "common_mode_note": CE.COMMON_MODE_NOTE,
        "concentration_form": concentration_form(col),
        "dwell_time_corrected": bool(col.ptype != "level" or col.dwell is not None),
        "sp2_null_layer_built": bool(ok_layer),
        "sp2_reason": why_layer,
        "subdaily": bool(col.subdaily),
        "human_schedule": bool(col.pclass == "human-schedule"),
        "provenance": dict(col.provenance),
        "analytic_sanity_check": analytic,
        "observer": observer_report,
    }
    return row


def promote(row, decl):
    """§S4.1 as superseded by SP-3. Returns (bool, list_of_failed_conditions).

    Every condition is checked and REPORTED, not short-circuited, so a near-miss row
    says which single condition it failed rather than only that it failed.
    """
    alpha = float(decl["alpha"])
    fails = []
    if not (row["p_real"] <= alpha):
        fails.append("p_real %.3e > alpha %.3e (SP-3, alpha = q/m)"
                     % (row["p_real"], alpha))
    pc = row.get("p_control")
    if pc is None:
        fails.append("no control arm: %s" % CONTROL_MANDATORY)
    elif not (pc > alpha):
        fails.append("control fired at p_control %.3e <= alpha %.3e (F7-b)"
                     % (pc, alpha))
    if row["concentration_form"] not in decl["declared_forms"]:
        fails.append("concentration form %r is not in the declared forms %r"
                     % (row["concentration_form"], decl["declared_forms"]))
    if row["n_events"] < int(decl["n_events_floor"]):
        fails.append("n_events %d < declared floor %d (§S6.2: 10-of-14 fails m = "
                     "336,000 at 0.683 while 14-of-14 clears by four orders -- the "
                     "instrument detects near-total concentration at small n and "
                     "nothing weaker)" % (row["n_events"], decl["n_events_floor"]))
    if not row["dwell_time_corrected"]:
        fails.append("level property with no dwell measure: %s" % P.DWELL_RULE)
    if not row["sp2_null_layer_built"]:
        fails.append(row["sp2_reason"])
    if row["human_schedule"] and row["mc"] < float(decl["human_schedule_mc_floor"]):
        fails.append("F7-d: a human-schedule property at Mc = %.1f < %.1f is a pure "
                     "OBSERVER measurement and is not eligible for promotion. %s"
                     % (row["mc"], decl["human_schedule_mc_floor"],
                        P.HUMAN_SCHEDULE_CARVE_OUT))
    if row["region"] in set(LAT.excluded_names_for_family(row["family"])):
        fails.append("region %r has seeded this property family and is excluded BY "
                     "NAME, permanently (§S1.1(c))" % row["region"])
    return (not fails), fails


# ------------------------------------------------------------------ the scan --
def run_scan(decl, cells, rng=None, control_mode="declared-cv", log_explore=True,
             explore_path=None, verbose=False):
    """The scan. `cells` is an iterable of prepared cell specs; returns the report.

    A cell spec is a dict:
        {'region': lattice region dict, 'mc': float, 'column': PropertyColumn,
         'obs_times': array, 'null_times': [arrays], 'ref_times': array,
         'observer': dict|None, 'null_provenance': str}

    The caller builds the specs, because building them requires a catalogue and
    §P7-22(b)(i)'s discipline is that the module that DECLARES the lattice never
    touches one. `run_scan` asserts that exactly `m` cells arrived (SP-6.1) and that
    the control arm exists (F7-b).

    **The `EXPLORE_COUNT.jsonl` append happens HERE, in the parent, exactly once**, on
    the same rule the miner follows: a ledger line written from a worker is a ledger
    line written an unknown number of times.
    """
    rng = rng or np.random.default_rng(20260813)
    rows, null_cols = [], []
    specs = list(cells)
    assert_cell_count(len(specs), decl["m"])

    for i, spec in enumerate(specs):
        col = spec["column"]
        span = spec.get("span_days")
        row = scan_cell(col, spec["obs_times"], spec["null_times"],
                        spec["ref_times"], spec["region"], spec["mc"], decl,
                        observer_report=spec.get("observer"),
                        null_provenance=spec.get("null_provenance", "UNDECLARED"))
        ctrl_col = spec.get("control_column") or _control_column(
            col, rng, mode=control_mode, span_days=span)
        crow = scan_cell(ctrl_col, spec["obs_times"],
                         spec.get("control_null_times", spec["null_times"]),
                         spec["ref_times"], spec["region"], spec["mc"], decl,
                         observer_report=spec.get("observer"),
                         null_provenance=spec.get("null_provenance", "UNDECLARED"))
        row["p_control"] = float(crow["p_real"])
        row["control_value"] = float(crow["value"])
        row["control_construction"] = (
            "F8-15-style random monotone clock warp (%s, CV = %.2f, knots every "
            "%.0f d), identical alpha, same multiplicity (F7-b). %s"
            % (control_mode, CONTROL_CV, CONTROL_KNOT_DAYS,
               CONTROL_DEGENERACY_FOUND))
        ok, fails = promote(row, decl)
        row["promoted"] = bool(ok)
        row["promotion_failures"] = fails
        cok, cfails = promote({**crow, "p_control": 1.0}, decl)
        row["control_promoted"] = bool(cok)
        row["control_promotion_failures"] = cfails
        rows.append(row)
        null_cols.append((row["null_draws"], crow["null_draws"]))
        if verbose and (i + 1) % 25 == 0:
            print("  cell %d/%d" % (i + 1, len(specs)))

    n_real = int(sum(1 for r in rows if r["promoted"]))
    n_ctrl = int(sum(1 for r in rows if r["control_promoted"]))

    # -- S-8 max-statistic, reported whether or not anything promotes (SP-6.4) ----
    b = min(int(d.size) for d, _c in null_cols)
    null_matrix = np.column_stack([d[:b] for d, _c in null_cols])
    observed = np.array([r["value"] for r in rows], dtype=np.float64)
    stratum_names = [r["layer"] or "L1" for r in rows]
    strata_table = LAT.region_strata(
        l3_cells=[{"name": n} for n in
                  sorted({r["region"] for r in rows if r["layer"] == "L3"})],
        q=decl["q"],
        n_properties=len(decl["properties"]), n_mag_strata=len(decl["mag_strata"]))
    present = {t["name"] for t in strata_table} & set(stratum_names)
    maxstat = ST.max_statistic_report(
        null_matrix, observed, stratum_names,
        strata=[t for t in strata_table if t["name"] in present] or None)

    # -- the ranked list, capped at K, by the PRE-DECLARED ranking function -------
    ranked = sorted(rows, key=lambda r: (-(-math.log10(max(r["p_real"], 1e-300))),
                                         r["region"], r["property"]))
    ranked = ranked[:int(decl["k_cap"])]
    for k, r in enumerate(ranked, start=1):
        r["rank_in_scan"] = k

    report = {
        "searcher_rule_id": SEARCHER_RULE_ID,
        "scan_id": decl["scan_id"],
        "config_hash": decl["config_hash"],
        "m_declared": int(decl["m"]),
        "n_cells_evaluated": len(rows),
        "alpha": float(decl["alpha"]),
        "q": float(decl["q"]),
        "n_survivors_real": n_real,
        "n_survivors_control": n_ctrl,
        "survivor_difference": n_real - n_ctrl,
        "control_calibration_note": CONTROL_MANDATORY,
        "max_statistic": maxstat,
        "ranked": [{k: v for k, v in r.items() if k != "null_draws"}
                   for r in ranked],
        "rows": [{k: v for k, v in r.items() if k != "null_draws"} for r in rows],
        "banner": GENERATOR_NOT_EVIDENCE,
        "the_sentence": THE_SENTENCE,
        "or_limb": OR_LIMB_REFUSED,
        "zero_outputs": ("§S3.3(3): zero stubs, zero K-entries, zero findings. "
                         "Unpriced rows cannot produce findings (§P7-16(4)) and the "
                         "same rows cannot drive a confirmatory statistic (§P7-22 "
                         "Ratification)."),
    }

    if log_explore:
        # SP-1.4: the scan appends its FULL declared cell count even though unpriced.
        # Unpriced never means uncounted. PARENT ONLY, exactly once.
        line = dict(decl)
        line["n_declared_tests"] = int(decl["m"])
        line["priced"] = 0
        line["category"] = "EXPLORATORY-UNPRICED (Rule 4.4 / §P7-4 / SP-1)"
        line["n_survivors_real"] = n_real
        line["n_survivors_control"] = n_ctrl
        splits.log_explore_run(line, explore_path or splits.EXPLORE_COUNT)
        report["explore_count_logged"] = True
    else:
        report["explore_count_logged"] = False
    return report


# ------------------------------------------------------------- entry guards ---
def assert_scan_entry(decl, columns, region_names, observer_feature_names=None,
                      event_days=None, explore_window=None):
    """Every gate that must pass BEFORE the first statistic of a real scan.

    1. **SP-6.2 -- scans never touch the temporal holdout.** Absolute, no exceptions,
       no flag. Checked on the event times against the declared exploration window.
    2. **§S1.1(c) -- the per-property-family by-name exclusions**, delegated to
       `lattice_s1.assert_family_exclusions` which delegates the tidal case to
       `regions_d.assert_alaska_excluded`: one implementation, not three.
    3. **§S2.2(1) -- the sub-daily observer gate**, `observer.assert_subdaily_gate`,
       whenever any declared column is sub-daily. The SEARCHER lives entirely on the
       sub-daily path and therefore inherits it as a HARD ENTRY GATE.
    4. **F8-15 -- the random-clock control is declared**, since the control arm is
       built from one and `clocks.assert_random_clock_control` is the check that a
       clock scan has one.
    """
    out = {}
    if event_days is not None and explore_window is not None:
        d = np.asarray(event_days, dtype=np.float64)
        lo, hi = float(explore_window[0]), float(explore_window[1])
        bad = int(np.count_nonzero((d < lo) | (d > hi)))
        if bad:
            raise HoldoutTouched(
                "%d of %d event times lie outside the declared exploration window "
                "[%g, %g]. §P7-24 SP-6.2: scans never touch the temporal holdout -- "
                "absolute, no exceptions, no flag." % (bad, d.size, lo, hi))
        out["holdout"] = {"window": [lo, hi], "n_events": int(d.size),
                          "n_outside": 0, "rule": "SP-6.2"}
    fams = sorted({c.family for c in columns})
    for f in fams:
        LAT.assert_family_exclusions(f, region_names)
    out["family_exclusions"] = {f: list(LAT.excluded_names_for_family(f))
                                for f in fams}
    if any(c.subdaily for c in columns):
        out["subdaily_gate"] = OBS.assert_subdaily_gate(observer_feature_names or [])
    out["random_clock_control"] = CL.assert_random_clock_control(
        ["scan-time-axis"], ["F8-15-style-declared-cv-%.2f" % CONTROL_CV])
    out["control_degeneracy_note"] = CONTROL_DEGENERACY_FOUND
    out["banner"] = GENERATOR_NOT_EVIDENCE
    out["the_sentence"] = THE_SENTENCE
    return out


def print_banners():
    """The two banners §S8.2 B5 requires the CLI to print."""
    print("*** " + GENERATOR_NOT_EVIDENCE)
    print()
    print("*** " + THE_SENTENCE)
    print()


if __name__ == "__main__":          # pragma: no cover - operator convenience
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    print_banners()
    print("engine.searcher is a library; the operator entry point is "
          "`python -u -m engine.cli search` (B5) and the SP-7 gate is "
          "`python -u -m engine.searcher_gate`.")
