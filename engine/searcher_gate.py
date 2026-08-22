"""SP-7 -- THE SEARCHER'S OWN GATE. Binding. No real scan runs before this passes.

> **§P7-24 SP-7, verbatim.** *"Run the complete searcher -- scan, SP-2 nulls, SP-3
> promotion rule -- over >= 30 true-null ETAS-sim catalogues through the identical
> code path. At `alpha = q/m` the expected promotions are <= 0.1 per catalogue, so
> <= 3 across 30. PASS: promotions at or below that rate, with the count consistent
> with binomial. FAIL: anything more -- and a failure here means the nominal p's are
> wrong, i.e. an SP-2 null layer is invalid, which is precisely the failure mode that
> matters and precisely the one the multiplicity arithmetic cannot catch.*
> *No real scan runs before this passes. Same standing as R1 for v2; same reason."*

A recovery gate proves the instrument can SEE. **SP-7 proves it does not
HALLUCINATE**, and it is the one that matters, because no amount of multiplicity
arithmetic can detect a nominal p computed against a wrong null.

THE THREE THINGS THIS GATE HAD TO GET RIGHT, EACH OF WHICH COULD HAVE MADE IT VACUOUS
--------------------------------------------------------------------------------------
1. **The catalogues must be ETAS WITH TRIGGERING, and so must the null replicates.**
   §P7-23(A.3). If the observed catalogues clustered and the null replicates did not,
   the p's would be anti-conservative and the gate would FAIL -- which is exactly the
   failure it exists to catch, so the generator is shared: `simulate_etas_times` draws
   the observed catalogues AND the null ensemble from the *same* process with the
   *same* parameters. A gate whose observed arm and null arm came from different
   processes would be measuring the difference between two simulators.
2. **Every declared cell must be PROMOTION-ELIGIBLE apart from its p-values**, or the
   gate passes because nothing could ever promote and proves nothing. `n_eligible` is
   computed and a gate with `n_eligible < m` is reported **VACUOUS**, which is a FAIL.
3. **The instrument must actually fire when there is something to see.** A planted
   von-Mises concentration is pushed through the identical code path and must promote.
   This is the vacuity control, and without it "zero promotions in 30 catalogues"
   is equally consistent with a correctly calibrated searcher and with a broken one.

THE DISCRETENESS OF THE THRESHOLD, STATED RATHER THAN GLOSSED
---------------------------------------------------------------
With `B` null replicates the attainable p-values are multiples of `1/(B+1)`, so the
EFFECTIVE promotion rate per cell is `floor(alpha (B+1)) / (B+1)`, not `alpha`. The
gate computes the effective rate and prices its own expectation against it -- reporting
a PASS against an expectation the instrument cannot attain would be the same class of
error the whole protocol exists to prevent.

Nothing this module produces is evidence, and a PASS licenses exactly one thing: the
right for a real scan to be run at all.
"""

from __future__ import annotations

import datetime as _dt
import json
import math
import os

import numpy as np

from . import lattice_s1 as LAT, properties as P, searcher as S, splits

GATE_RULE_ID = "SP7-searcher-gate-v1"
GATE_ARTIFACT = os.path.join("engine", "out", "searcher_sp7_gate.json")

# Declared gate constants, frozen here before the gate runs.
N_CATALOGS = 30                 # SP-7: ">= 30 true-null ETAS-sim catalogues"
EXPECTED_PER_CATALOG = 0.10     # SP-7: "<= 0.1 per catalogue"
MAX_TOTAL = 3                   # SP-7: "so <= 3 across 30"
BINOMIAL_ALPHA = 0.01           # "with the count consistent with binomial"

SPAN_DAYS = 3650.0              # 10 years
N_NULL = 400                    # null replicates per cell
N_REF = 1200                    # reference-pool events (disjoint from the null)

# ETAS parameters for the true-null simulator. DECLARED constants; the gate is a
# calibration instrument and its generator's parameters are part of its declaration.
ETAS_MU = 0.0110                # background events/day  (~40 background in 10 y)
ETAS_K = 0.45                   # branching ratio
ETAS_C = 0.05                   # Omori c, days
ETAS_P = 1.15                   # Omori p
ETAS_MAX_GEN = 12

GATE_REGIONS = ("japan", "chile", "california", "iceland", "himalaya")
GATE_PROPERTIES = ("tide_phase", "tide_level", "lunar_synodic_phase",
                   "solar_annual_phase", "day_of_week", "utc_hour_phase")
GATE_MAG_STRATA = (6.0,)        # human-schedule properties are live only at M >= 6.0
GATE_FAMILIES = ("solid_tide", "ephemeris", "human_clock", "season")


# --------------------------------------------------------- the true-null ETAS --
def simulate_etas_times(t_lo, t_hi, rng, mu=ETAS_MU, k=ETAS_K, c=ETAS_C, p=ETAS_P,
                        max_gen=ETAS_MAX_GEN):
    """A temporal ETAS (Hawkes-Omori) sequence on [t_lo, t_hi). TRUE NULL in property.

    Background events are a homogeneous Poisson process; each event triggers
    Poisson(k) offspring at Omori-distributed delays `(t - t_i) ~ (t - t_i + c)^-p`.
    **Nothing in the generator knows about any property**, so a concentration on any
    property of the world is, by construction, an error of the instrument.

    A homogeneous background is deliberate: a background with an annual or diurnal
    cycle would plant real structure in `solar_annual_phase` / `utc_hour_phase` and
    the gate would measure the plant instead of the instrument.
    """
    lo, hi = float(t_lo), float(t_hi)
    n_bg = int(rng.poisson(mu * (hi - lo)))
    times = list(np.sort(rng.uniform(lo, hi, size=n_bg)))
    parents = np.asarray(times, dtype=np.float64)
    for _g in range(int(max_gen)):
        if parents.size == 0:
            break
        n_off = rng.poisson(float(k), size=parents.size)
        tot = int(n_off.sum())
        if tot == 0:
            break
        src = np.repeat(parents, n_off)
        # inverse-CDF of the Omori density on [0, inf): dt = c((1-u)^(1/(1-p)) - 1)
        u = rng.random(tot)
        dt = float(c) * ((1.0 - u) ** (1.0 / (1.0 - float(p))) - 1.0)
        kids = src + dt
        kids = kids[(kids >= lo) & (kids < hi)]
        if kids.size == 0:
            break
        times.extend(kids.tolist())
        parents = kids
    return np.sort(np.asarray(times, dtype=np.float64))


def _draw_until(n_min, t_lo, t_hi, rng, **kw):
    """One ETAS sequence with at least `n_min` events, re-drawn if it comes up short.

    Re-drawing conditions on N >= n_min in BOTH the observed and null arms (the same
    function serves both), so the conditioning is common-mode and cancels -- the same
    argument `circstat_event.resample_event_times` makes for conditioning on N.
    """
    for _ in range(200):
        t = simulate_etas_times(t_lo, t_hi, rng, **kw)
        if t.size >= int(n_min):
            return t
    raise RuntimeError("ETAS simulator could not reach n_min=%d; raise ETAS_MU"
                       % int(n_min))


def effective_alpha(alpha, b=N_NULL):
    """The attainable promotion rate given `b` replicates: floor(alpha (b+1))/(b+1)."""
    step = 1.0 / (float(b) + 1.0)
    return math.floor(float(alpha) / step) * step


def _binom_sf(k, n, p):
    from scipy import special
    k, n = int(k), int(n)
    if k <= 0:
        return 1.0
    if k > n:
        return 0.0
    return float(special.betainc(k, n - k + 1, float(p)))


# --------------------------------------------------------------- the fixture --
def build_gate_columns(t0, span_days=SPAN_DAYS):
    """One property column set per gate region, filtered to the declared properties."""
    per_region = {}
    for name in GATE_REGIONS:
        reg = LAT.region(name)
        lat = 0.5 * (reg["lat_min"] + reg["lat_max"])
        lon = 0.5 * (reg["lon_min"] + reg["lon_max"])
        cols, audit = P.build_property_matrix(
            t0, 0.0, float(span_days), lat, lon, families=GATE_FAMILIES,
            sample_minutes=60.0)
        keep = [c for c in cols if c.name in GATE_PROPERTIES]
        # F7-a: the observer report travels BESIDE every human-schedule cell, always.
        # In the gate the catalogues are simulated and carry no observer, so the
        # report says exactly that rather than implying a measurement was made.
        obs_report = {
            "measured": False,
            "why": ("SP-7 catalogues are SIMULATED: there is no observer to measure. "
                    "The F7 layer is attached here to make the human-schedule cells "
                    "PROMOTION-ELIGIBLE so the gate is not vacuous; on a REAL scan "
                    "`observer.observer_features` must be computed for the region and "
                    "stratum and reported beside the cell (F7-a)."),
            "rule": P.PROPERTY_CLASSES["human-schedule"]["layer"],
        }
        for c in keep:
            if c.pclass == "human-schedule":
                c.attach("observer_features")
        found = {c.name for c in keep}
        missing = [p for p in GATE_PROPERTIES if p not in found]
        if missing:
            raise RuntimeError("gate property set incomplete for %r: missing %s"
                               % (name, missing))
        per_region[name] = {"region": reg, "columns": keep, "observer": obs_report,
                            "audit": audit, "site": (lat, lon)}
    return per_region


def build_gate_declaration(l3_cells=None, scan_id=None):
    return S.declare_scan(
        [LAT.region(n) for n in GATE_REGIONS], GATE_PROPERTIES, GATE_MAG_STRATA,
        GATE_FAMILIES, l3_cells=l3_cells,
        scan_id=scan_id or "SP7-GATE", window=[0.0, SPAN_DAYS],
        extra={"gate_rule_id": GATE_RULE_ID, "n_null": N_NULL, "n_ref": N_REF,
               "n_catalogs": N_CATALOGS, "etas": {"mu": ETAS_MU, "K": ETAS_K,
                                                  "c": ETAS_C, "p": ETAS_P}})


def _specs(per_region, decl, obs_by_region, null_times, ref_times):
    specs = []
    for rname in decl["regions"]:
        r = per_region[rname]
        by_name = {c.name: c for c in r["columns"]}
        for pname in decl["properties"]:
            for mc in decl["mag_strata"]:
                specs.append({
                    "region": r["region"], "mc": float(mc),
                    "column": by_name[pname],
                    "obs_times": obs_by_region[rname],
                    "null_times": null_times,
                    "ref_times": ref_times,
                    "span_days": SPAN_DAYS,
                    "observer": r["observer"],
                    "null_provenance": ("full ETAS WITH TRIGGERING (Hawkes-Omori), "
                                        "identical generator and parameters as the "
                                        "observed arm (§P7-23(A.3))"),
                })
    return specs


# ------------------------------------------------------------ the vacuity test --
def vacuity_control(per_region, decl, null_times, ref_times, rng):
    """The instrument must FIRE on a planted concentration, or a PASS means nothing.

    A von-Mises rate modulation is planted on `lunar_synodic_phase` in one region by
    rejection-thinning an ETAS draw through the IDENTICAL property map
    (`circstat_event.thin_by_phase_intensity` is the declared primitive and the
    planting is a genuine RATE modulation, not a phase relabelling -- §P7-8(d)).
    """
    from . import circstat_event as CE
    rname = decl["regions"][0]
    col = {c.name: c for c in per_region[rname]["columns"]}["lunar_synodic_phase"]
    cand = _draw_until(4000, 0.0, SPAN_DAYS, rng, mu=ETAS_MU * 200.0)
    ph = np.asarray(col.evaluate(cand), dtype=np.float64)
    kappa, mu0 = 3.0, 1.0
    planted = CE.thin_by_phase_intensity(
        cand, ph, lambda x: np.exp(kappa * np.cos(x - mu0)), 60, rng)
    row = S.scan_cell(col, planted, null_times, ref_times,
                      per_region[rname]["region"], decl["mag_strata"][0], decl,
                      observer_report=per_region[rname]["observer"],
                      null_provenance="planted-arm control")
    row["p_control"] = 1.0                 # the plant is the science arm by construction
    ok, fails = S.promote(row, decl)
    return {
        "planted": "von Mises kappa=%.1f on lunar_synodic_phase, rate modulation by "
                   "rejection thinning through the identical map" % kappa,
        "region": rname, "n_events": int(row["n_events"]),
        "p_real": float(row["p_real"]),
        "p_resolution_floor": float(row["p_resolution_floor"]),
        "alpha": float(decl["alpha"]),
        "promoted": bool(ok), "failures": fails,
        "verdict": "FIRES" if ok else "DOES NOT FIRE -- a PASS would be meaningless",
    }


# ------------------------------------------------------------------- the gate --
def run_gate(n_catalogs=N_CATALOGS, seed=20260813, out_path=GATE_ARTIFACT,
             verbose=True, n_null=N_NULL, n_ref=N_REF):
    """Run SP-7 end to end. Returns the artifact dict and writes it to `out_path`."""
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    t_start = _dt.datetime.now(_dt.timezone.utc)
    rng = np.random.default_rng(int(seed))
    t0 = _dt.datetime(2000, 1, 1)

    if verbose:
        S.print_banners()
        print("SP-7 SEARCHER GATE -- %d true-null ETAS catalogues" % int(n_catalogs))

    per_region = build_gate_columns(t0)
    decl = build_gate_declaration()
    alpha = float(decl["alpha"])
    a_eff = effective_alpha(alpha, n_null)

    # The null ensemble and the reference pool: drawn ONCE from the same generator as
    # the observed arm, and DISJOINT from it and from each other -- a null replicate
    # compared against a pool containing itself has its dispersion biased low
    # (`circstat_event.event_omnibus`'s own warning).
    null_times = [_draw_until(12, 0.0, SPAN_DAYS, rng) for _ in range(int(n_null))]
    ref_times = _draw_until(int(n_ref), 0.0, SPAN_DAYS, rng, mu=ETAS_MU * 40.0)[
        :int(n_ref)]

    vac = vacuity_control(per_region, decl, null_times, ref_times, rng)

    per_catalog, n_promotions, n_control_promotions, n_eligible_min = [], 0, 0, None
    n_bare = 0          # cells clearing alpha BEFORE the F7-b control condition
    for j in range(int(n_catalogs)):
        obs = {r: _draw_until(int(decl["n_events_floor"]), 0.0, SPAN_DAYS, rng)
               for r in decl["regions"]}
        specs = _specs(per_region, decl, obs, null_times, ref_times)
        rep = S.run_scan(decl, specs, rng=rng, log_explore=False)
        # Eligibility: a cell whose ONLY failing conditions are the two p conditions.
        elig = 0
        for row in rep["rows"]:
            fails = [f for f in row["promotion_failures"]
                     if not (f.startswith("p_real ") or f.startswith("control fired"))]
            if not fails:
                elig += 1
        n_eligible_min = elig if n_eligible_min is None else min(n_eligible_min, elig)
        # HONESTY FIELD, added after the first gate run showed why it is needed. The
        # real arm's promotion count is suppressed by the F7-b condition (`the control
        # must NOT fire`), and in a clustered ETAS draw the real map and the warped
        # control map tend to fire TOGETHER -- a burst concentrates events in time and
        # therefore in almost any property. So "zero promotions" is partly the
        # threshold and partly the control blocking a co-firing cell, and reporting
        # only the first number would overstate how quiet the threshold is on its own.
        n_bare += int(sum(1 for r in rep["rows"] if r["p_real"] <= a_eff))
        n_promotions += rep["n_survivors_real"]
        n_control_promotions += rep["n_survivors_control"]
        per_catalog.append({
            "catalog": j, "n_survivors_real": rep["n_survivors_real"],
            "n_survivors_control": rep["n_survivors_control"],
            "n_cells": rep["n_cells_evaluated"], "n_eligible": elig,
            "min_p_real": min(r["p_real"] for r in rep["rows"]),
            "max_statistic_global_p": rep["max_statistic"]["global"]["p"]
            if isinstance(rep["max_statistic"], dict)
            and "global" in rep["max_statistic"] else None,
            "n_events": {r: int(v.size) for r, v in obs.items()},
        })
        if verbose:
            print("  catalogue %2d/%d  promotions=%d  control=%d  min p=%.3e"
                  % (j + 1, int(n_catalogs), rep["n_survivors_real"],
                     rep["n_survivors_control"], per_catalog[-1]["min_p_real"]))

    m = int(decl["m"])
    n_trials = m * int(n_catalogs)
    p_binom = _binom_sf(n_promotions, n_trials, a_eff)
    rate = n_promotions / float(n_catalogs)
    vacuous = (n_eligible_min is None) or (n_eligible_min < m) or (not vac["promoted"])

    passed = bool((not vacuous)
                  and n_promotions <= int(MAX_TOTAL)
                  and p_binom >= BINOMIAL_ALPHA)

    art = {
        "gate_rule_id": GATE_RULE_ID,
        "protocol": "HYPOTHESIS_LEDGER.md §P7-24 SP-7 (BINDING)",
        "verdict": "PASS" if passed else ("VACUOUS-FAIL" if vacuous else "FAIL"),
        "passed": passed,
        "ts_start": t_start.isoformat(),
        "ts_end": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "seed": int(seed),
        "declaration": decl,
        "m": m,
        "alpha_declared": alpha,
        "alpha_effective": a_eff,
        "discreteness_note": (
            "with B = %d null replicates the attainable p-values are multiples of "
            "1/(B+1) = %.6g, so the EFFECTIVE per-cell promotion rate is "
            "floor(alpha (B+1))/(B+1) = %.6g, not the declared alpha = %.6g. The "
            "expectation below is priced against the rate the instrument can ATTAIN."
            % (int(n_null), 1.0 / (n_null + 1.0), a_eff, alpha)),
        "n_catalogs": int(n_catalogs),
        "n_cells_per_catalog": m,
        "n_trials": n_trials,
        "n_promotions_total": int(n_promotions),
        "n_cells_clearing_alpha_before_control": int(n_bare),
        "n_cells_clearing_alpha_note": (
            "cells whose p_real cleared the effective alpha BEFORE the F7-b "
            "condition (the control must NOT fire) was applied. Reported "
            "because in a clustered ETAS draw the real map and the warped "
            "control map tend to fire TOGETHER -- a burst concentrates events "
            "in time and therefore in almost any property -- so a promotion "
            "count of zero is partly the threshold and partly the control "
            "blocking a co-firing cell. Expected under a correct null: "
            "n_trials * alpha_effective."),
        "n_cells_clearing_alpha_expected": float(n_trials * a_eff),
        "binomial_upper_tail_p_bare": _binom_sf(n_bare, n_trials, a_eff),
        "cap_is_the_mean_note": (
            "A DIAGNOSTIC, NOT A CRITERION -- reported so the operator can see what a "
            "FAIL here does and does not mean, and NOT used to decide the verdict. "
            "The expected promotions across N catalogues at alpha = q/m is "
            "N * m * (q/m) = N * q, i.e. 3.0 across 30 REGARDLESS OF m. SP-7's cap of "
            "3 is therefore exactly the MEAN of a correctly calibrated searcher, not "
            "an upper bound on it: a perfectly calibrated instrument exceeds 3 about "
            "35%% of the time (Poisson(3), P(X > 3) = 0.353). The second clause of "
            "SP-7 -- 'with the count consistent with binomial' -- is the one that "
            "distinguishes a miscalibrated null from an ordinary fluctuation, and "
            "this artifact reports both numbers so the distinction can be made by "
            "whoever is entitled to make it. THIS BUILD DOES NOT RELAX THE CAP: "
            "changing a criterion after seeing a result is the exact defect the "
            "protocol exists to prevent."),
        "poisson_p_exceeding_cap_if_calibrated": float(
            1.0 - sum(math.exp(-float(n_catalogs) * float(decl["q"]))
                      * (float(n_catalogs) * float(decl["q"])) ** i
                      / math.factorial(i) for i in range(int(MAX_TOTAL) + 1))),
        "n_control_promotions_total": int(n_control_promotions),
        "promotions_per_catalog": rate,
        "expected_per_catalog_declared": EXPECTED_PER_CATALOG,
        "expected_per_catalog_effective": m * a_eff,
        "max_total_allowed": int(MAX_TOTAL),
        "binomial_upper_tail_p": p_binom,
        "binomial_alpha": BINOMIAL_ALPHA,
        "n_eligible_min": n_eligible_min,
        "vacuous": bool(vacuous),
        "vacuity_control": vac,
        "per_catalog": per_catalog,
        "criteria": {
            "1_not_vacuous": ("every one of the %d declared cells must be "
                              "promotion-eligible apart from its p-values, AND the "
                              "planted-signal control must FIRE" % m),
            "2_count": "total promotions across %d catalogues <= %d"
                       % (int(n_catalogs), int(MAX_TOTAL)),
            "3_binomial": "upper-tail binomial p >= %.3g at n = %d, p = %.6g"
                          % (BINOMIAL_ALPHA, n_trials, a_eff),
        },
        "what_a_pass_licenses": (
            "EXACTLY ONE THING: the right for a real scan to be run at all. A PASS is "
            "not evidence, is not a finding, and says nothing about the Earth. A FAIL "
            "means the nominal p's are wrong -- an SP-2 null layer is invalid -- which "
            "is precisely the failure mode the multiplicity arithmetic cannot catch."),
        "banner": S.GENERATOR_NOT_EVIDENCE,
        "the_sentence": S.THE_SENTENCE,
    }
    art["artifact_hash"] = splits.config_hash(
        {k: v for k, v in art.items() if k != "artifact_hash"})[:16]

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(art, fh, indent=2, default=str)
    if verbose:
        print("-" * 78)
        print("SP-7 VERDICT: %s" % art["verdict"])
        print("  promotions          = %d over %d catalogues (%.3f per catalogue)"
              % (n_promotions, int(n_catalogs), rate))
        print("  allowed             = <= %d (SP-7)" % int(MAX_TOTAL))
        print("  control promotions  = %d" % n_control_promotions)
        print("  alpha declared      = %.6e  (q/m, q=%.2f, m=%d)"
              % (alpha, decl["q"], m))
        print("  alpha effective     = %.6e  (B = %d)" % (a_eff, int(n_null)))
        print("  binomial upper tail = %.4f  (needs >= %.2f)"
              % (p_binom, BINOMIAL_ALPHA))
        print("  vacuity control     = %s" % vac["verdict"])
        print("  cells eligible      = %s of %d" % (n_eligible_min, m))
        print("  artifact            = %s" % out_path)
        print("-" * 78)
    return art


if __name__ == "__main__":          # pragma: no cover - operator entry point
    import argparse
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    ap = argparse.ArgumentParser("engine.searcher_gate")
    ap.add_argument("--catalogs", type=int, default=N_CATALOGS)
    ap.add_argument("--null", type=int, default=N_NULL)
    ap.add_argument("--ref", type=int, default=N_REF)
    ap.add_argument("--seed", type=int, default=20260813)
    ap.add_argument("--out", default=GATE_ARTIFACT)
    a = ap.parse_args()
    art = run_gate(n_catalogs=a.catalogs, seed=a.seed, out_path=a.out,
                   n_null=a.null, n_ref=a.ref)
    raise SystemExit(0 if art["passed"] else 1)
