"""§P7-15: the recovery-versus-amplitude curve. Priced 0, assigned as a diagnostic.

THE QUESTION, AND WHY IT IS NOT COSMETIC
----------------------------------------
§P7-15 flagged an arithmetic that does not close. At the tranche-A operating
threshold the §P7-12(a) floor formula predicts recovery of **0.800 at 1x the
floor, 0.9992 at 1.5x and 1.0000 at 2x**. The G-M1 gate (results_gate_r1.json)
returned **80% at >= 2x** on arm (i) and 90% on arm (iii) -- arm (i)'s number is
what the formula predicts at ONE x, not two. Either

  (a) the plants were mislabelled and effectively sat at ~1x (formula fine), or
  (b) the effective global floor is ~2x larger than the formula gives, i.e. the
      formula is ANTI-CONSERVATIVE at global aggregation, or
  (c) mixed -- one aggregation is fine and another is not.

"It is the difference between a bound quoted at ~15% and one quoted at ~30%
globally" (§P7-15). So it is measured, not argued.

WHAT THIS RUNS
--------------
The SAME machinery as the gate: `engine.gate_r1`'s plant construction,
simulation, sandboxing, pipeline invocation and recovery extraction, at the same
config hash. Only the plant AMPLITUDE moves, across 1x / 1.5x / 2x / 3x the
operative floor, on each of the three aggregations the gate cleared:

    arm (i)   global aggregation, mid-band       (moon_synodic_phase, 29.53 d)
    arm (ii)  regional 2R-df, region-dependent phase (regsum)
    arm (iii) global aggregation, long-period    (annual_phase, 365.24 d)

>= 10 catalogs per (arm, amplitude) cell. The 2x cells are NOT re-run: the gate
already executed exactly that cell -- same config hash, same feature, same
factor, same construction, same seeds -- and re-drawing it would spend four
hours to obtain a second sample of a quantity already measured. They are read
from `engine/out/gate_r1/progress.json` and labelled `source: gate_r1`.

A DELIBERATE DEPARTURE FROM §P7-8(d), DECLARED
----------------------------------------------
§P7-8(d) requires plants at >= 2x the operative floor, precisely so that a miss
is an instrument verdict and not a power verdict. **This sweep deliberately
plants BELOW that line**, at 1x and 1.5x, because measuring power below the
2x line is the entire assignment: the curve cannot be measured without visiting
the amplitudes whose power is in question. Every cell therefore carries
`p7_8d_compliant` and the sub-2x cells are explicitly NOT instrument verdicts --
a miss there is a POWER measurement, which is what it is for. No G-M1 verdict
is taken from this module; G-M1 was decided in `results_gate_r1.json` at >= 2x.

THE PREDICTED CURVE, DERIVED RATHER THAN COPIED
-----------------------------------------------
§P7-12(a)'s floor is `A_min = sqrt(VIF) * (z_alpha + z_0.80) * sqrt(2/N)`, i.e.
the amplitude at which the normal-approximation power is exactly 0.80 at a
two-sided threshold alpha. Inverting the same construction, the power at
`A = k * A_min` is

    power(k) = Phi( k * (z_alpha + z_0.80) - z_alpha )

`predicted_power` below is that function; `test_recovery_curve.py` pins it
against §P7-15's own three published values (0.800 / 0.9992 / 1.0000), so if it
ever stops reproducing them the module is wrong rather than the ledger.

THE ESTIMATOR OF THE ANSWER
---------------------------
The quantity in dispute is a single number per arm: the amplitude at which the
pipeline ACTUALLY reaches 80% recovery, against the formula's `A_min`. So the
measured curve is fitted with the SAME one-parameter family, with the floor
free:

    power(A) = Phi( A * (z_alpha + z_0.80) / A_eff - z_alpha )

`A_eff` is the effective 80%-power floor, estimated by binomial maximum
likelihood over every catalog in the arm. `A_eff / A_min_formula` IS the
inflation factor §P7-15 asks for: 1.0 means the formula is right, ~2 means it is
anti-conservative by 2x at that aggregation. The raw per-cell rates and their
Wilson intervals are reported alongside and are the primary evidence; the fit is
a summary of them, not a substitute.

Usage:
    python -u -m engine.recovery_curve --jobs 8
    python -u -m engine.recovery_curve --jobs 8 --smoke   # 2 catalogs per cell
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import sys
import time

import numpy as np
from scipy import optimize, stats

from . import floors, gate_r1 as G, splits, __version__

PROGRESS = os.path.join(G.SANDBOX, "curve_progress.json")
RESULTS = "results_recovery_curve.json"
GATE_PROGRESS = os.path.join(G.SANDBOX, "progress.json")
GATE_RESULTS = "results_gate_r1.json"

# The declared sweep. 2.0 is present because it is the gate's own cell and the
# curve is meaningless without the point that raised the question.
FACTORS = (1.0, 1.5, 2.0, 3.0)
N_PER_CELL = 10
CURVE_SEED = 20260813

ARMS = {
    "i": {"feature": G.ARM_I_FEATURE, "kind": "global",
          "test_kind": "glm_poisson_offset_etas",
          "label": "global aggregation, mid-band (moon_synodic_phase, 29.53 d)"},
    "ii": {"feature": G.ARM_II_FEATURE, "kind": "regional",
           "test_kind": "regsum_score_2Rdf",
           "label": "regional 2R-df sum, region-dependent phase (regsum)"},
    "iii": {"feature": G.ARM_III_FEATURE, "kind": "global",
            "test_kind": "glm_poisson_offset_etas",
            "label": "global aggregation, long-period (annual_phase, 365.24 d)"},
}
# Stream ids are declared integers, never `hash()` (randomised per process), so
# a resumed sweep redraws the identical catalogs.
STREAM_ID = {("i", 1.0): 11, ("i", 1.5): 12, ("i", 3.0): 14,
             ("ii", 1.0): 21, ("ii", 1.5): 22, ("ii", 3.0): 24,
             ("iii", 1.0): 31, ("iii", 1.5): 32, ("iii", 3.0): 34}


# ----------------------------------------------------------- predicted curve --
def predicted_power(k, alpha=G.BASE_ALPHA):
    """§P7-12(a)'s own power function at `k` x its floor. See the module docstring."""
    za = floors.z_alpha(alpha)
    return float(stats.norm.cdf(float(k) * (za + floors.Z_POWER_80) - za))


def effective_floor_ratio(amplitude, a_eff, alpha=G.BASE_ALPHA):
    """The same family with the floor free: power at `amplitude` given `a_eff`."""
    za = floors.z_alpha(alpha)
    return float(stats.norm.cdf(
        float(amplitude) * (za + floors.Z_POWER_80) / float(a_eff) - za))


def fit_effective_floor(amplitudes, n_trials, n_hits, alpha=G.BASE_ALPHA,
                        a_grid=None):
    """Binomial MLE for `A_eff`, the amplitude at which power is truly 0.80.

    One free parameter, the same functional family as the formula, so the ratio
    `A_eff / A_min_formula` is a like-for-like statement about the formula and
    not a comparison between two different models.

    Returns None for `a_eff` when the data are degenerate (all cells 0% or all
    100%): the likelihood is then monotone and the MLE sits at a boundary, which
    is reported as a BOUND rather than dressed up as a point estimate.
    """
    amplitudes = np.asarray(amplitudes, dtype=np.float64)
    n_trials = np.asarray(n_trials, dtype=np.float64)
    n_hits = np.asarray(n_hits, dtype=np.float64)
    total, hits = float(n_trials.sum()), float(n_hits.sum())
    degenerate = (hits == 0.0) or (hits == total)

    def nll(log_a):
        a = math.exp(log_a)
        p = np.array([effective_floor_ratio(x, a, alpha) for x in amplitudes])
        p = np.clip(p, 1e-12, 1 - 1e-12)
        return -float((n_hits * np.log(p)
                       + (n_trials - n_hits) * np.log1p(-p)).sum())

    lo, hi = math.log(amplitudes.min() / 50.0), math.log(amplitudes.max() * 50.0)
    res = optimize.minimize_scalar(nll, bounds=(lo, hi), method="bounded")
    a_eff = float(math.exp(res.x))
    out = {"a_eff": (None if degenerate else a_eff),
           "a_eff_raw_mle": a_eff,
           "degenerate": bool(degenerate),
           "nll": float(res.fun), "converged": bool(res.success)}
    if degenerate:
        out["bound"] = ("all cells at 100%: A_eff is BELOW the smallest "
                        "amplitude visited" if hits == total else
                        "all cells at 0%: A_eff is ABOVE the largest "
                        "amplitude visited")
    else:
        # profile-likelihood 95% interval (chi2_1 cutoff 3.841 / 2)
        target = res.fun + 1.920729410347062
        try:
            lo_b = optimize.brentq(lambda x: nll(x) - target, lo, res.x)
            hi_b = optimize.brentq(lambda x: nll(x) - target, res.x, hi)
            out["a_eff_ci95"] = [float(math.exp(lo_b)), float(math.exp(hi_b))]
        except (ValueError, RuntimeError):
            out["a_eff_ci95"] = None
    return out


def amplitude_at_power(a_eff, power, alpha=G.BASE_ALPHA):
    """Invert the fitted curve: the amplitude at which recovery reaches `power`."""
    za = floors.z_alpha(alpha)
    return float(a_eff) * (za + float(stats.norm.ppf(power))) / (
        za + floors.Z_POWER_80)


def wilson(k, n, z=1.959963984540054):
    """Wilson score interval. Never a Wald interval at 10/10 or 0/10."""
    if n == 0:
        return [None, None]
    p = k / n
    d = 1.0 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return [float((c - h) / d), float((c + h) / d)]


def interpolate_threshold(points, target):
    """Smallest amplitude reaching `target` by linear interpolation on the RAW rates.

    `points` is [(amplitude, rate), ...]. Returns None (with a reason) when the
    curve never crosses the target inside the swept range, so the answer is a
    stated bound rather than an extrapolation -- S-17's rule applied to this
    module's own estimator.
    """
    pts = sorted(points)
    if not pts:
        return {"amplitude": None, "reason": "no cells"}
    if pts[0][1] >= target:
        return {"amplitude": None, "at_or_below": pts[0][0],
                "reason": (f"already at/above {target:.0%} at the SMALLEST "
                           f"amplitude swept ({pts[0][0]:.4f}); the threshold "
                           f"is a bound, not an interpolation")}
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if y0 < target <= y1:
            if y1 == y0:
                return {"amplitude": float(x1), "reason": "flat segment"}
            return {"amplitude": float(x0 + (target - y0) * (x1 - x0)
                                       / (y1 - y0)),
                    "reason": "linear interpolation between adjacent cells",
                    "bracket": [x0, x1]}
    return {"amplitude": None, "above": pts[-1][0],
            "reason": (f"never reaches {target:.0%} within the swept range "
                       f"(max amplitude {pts[-1][0]:.4f}, max rate "
                       f"{pts[-1][1]:.2f}); the threshold is a lower bound")}


# ------------------------------------------------------------------- driver ---
def _load(path, default):
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as fh:
            return json.load(fh)
    return default


def _save(path, obj):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=1, default=float)
    os.replace(tmp, path)


def build_plant(prep, partition, arm, factor):
    """The plant for one (arm, factor) cell, with its floor report(s)."""
    spec = ARMS[arm]
    if spec["kind"] == "global":
        mod, rep = G.global_plant(prep, spec["feature"], factor=factor)
        reps = [rep]
        truth = rep["planted_amplitude"]
    else:
        mod, reps = G.regional_plant(prep, partition, spec["feature"],
                                     factor=factor)
        truth = None            # §P7-1(d): the 2R-df sum reports UNRESOLVED
    for r in reps:
        r["declared_factor"] = float(factor)
        r["p7_8d_compliant"] = bool(factor >= floors.PLANT_FACTOR)
        r["p7_8d_note"] = (
            "at or above the §P7-8(d) 2x line: a miss here IS an instrument "
            "verdict" if factor >= floors.PLANT_FACTOR else
            "DELIBERATELY BELOW the §P7-8(d) 2x line. A miss here is a POWER "
            "measurement and is the point of this sweep; it is NOT an "
            "instrument verdict and no G-M1 conclusion is taken from it.")
    planted = {"feature": spec["feature"], "test_kind": spec["test_kind"],
               "truth_amplitude": truth}
    return mod, reps, planted


def reuse_gate_cells(prog):
    """The 2x cells, read from the gate's own artifact rather than re-run."""
    gate = _load(GATE_PROGRESS, {"catalogs": {}})
    gate_reports = gate.get("plant_reports", {})
    n = 0
    for tag, rec in gate["catalogs"].items():
        arm = rec.get("arm")
        if not arm or arm not in ARMS:
            continue
        key = f"curve_{arm}_2.0x_{tag.rsplit('_', 1)[1]}"
        if key in prog["cells"]:
            continue
        r = dict(rec)
        r.update({"arm": arm, "factor": 2.0, "source": "gate_r1",
                  "source_tag": tag,
                  "source_note": ("read from engine/out/gate_r1/progress.json: "
                                  "the gate executed exactly this cell at the "
                                  "same config hash and the same construction")})
        prog["cells"][key] = r
        n += 1
    if gate_reports:
        prog.setdefault("plant_reports", {})
        for arm, reps in gate_reports.items():
            prog["plant_reports"][f"{arm}@2.0x"] = reps
    return n


def main(argv=None):
    p = argparse.ArgumentParser("engine.recovery_curve")
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--per-cell", type=int, default=None)
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--arms", default="i,ii,iii")
    p.add_argument("--fresh", action="store_true")
    p.add_argument("--report-only", action="store_true")
    a = p.parse_args(argv)

    per_cell = a.per_cell if a.per_cell is not None else (2 if a.smoke
                                                          else N_PER_CELL)
    arms = [s.strip() for s in a.arms.split(",") if s.strip() in ARMS]

    os.makedirs(G.SANDBOX_SESSIONS, exist_ok=True)
    if a.fresh and os.path.exists(PROGRESS):
        os.remove(PROGRESS)
    prog = _load(PROGRESS, {"cells": {}, "plant_reports": {}})

    cfg = G.gate_config()
    cfg_hash = splits.config_hash(cfg)
    gate = _load(GATE_RESULTS, {})
    print("=" * 78)
    print(f"§P7-15 RECOVERY-VERSUS-AMPLITUDE CURVE -- engine v{__version__}")
    print("=" * 78)
    print(f"config hash          = {cfg_hash}")
    if gate.get("config_hash"):
        same = gate["config_hash"] == cfg_hash
        print(f"gate config hash     = {gate['config_hash']} "
              f"({'IDENTICAL -- cells are comparable' if same else 'DIFFERS'})")
        if not same:
            raise SystemExit(
                "the curve config hash differs from the gate's. The 2x cells "
                "could not then be reused and the curve would not be "
                "comparable to the number §P7-15 is asking about. Refusing.")
    n_reused = reuse_gate_cells(prog)
    print(f"reused from the gate = {n_reused} catalogs at factor 2.0")
    print(f"sweep                = factors {FACTORS}, arms {arms}, "
          f"{per_cell} catalogs per cell")
    print("predicted (§P7-12a):  " + ", ".join(
        f"{k}x -> {predicted_power(k):.4f}" for k in FACTORS))
    print("§P7-8(d) NOTE: the 1.0x and 1.5x cells are DELIBERATELY below the "
          "2x plant line;\n  a miss there is a POWER measurement, not an "
          "instrument verdict, and no G-M1\n  conclusion is taken from this "
          "module.")

    if not a.report_only:
        prep = G.prepare_base(cfg, verbose=True)
        partition = G.build_partition(prep, cfg)
        plants = {}
        for arm in arms:
            for k in FACTORS:
                if k == 2.0:
                    continue
                mod, reps, planted = build_plant(prep, partition, arm, k)
                plants[(arm, k)] = (mod, planted)
                prog["plant_reports"][f"{arm}@{k}x"] = [dict(r) for r in reps]
                for r in reps:
                    print(f"  plant {arm}@{k}x {r['what']}: A = "
                          f"{r['planted_amplitude']:.4f} vs floor "
                          f"{r['operative_floor_A_min']:.4f} at N = "
                          f"{r['N']:.0f} (x{r['amplitude_over_floor']:.2f})"
                          + ("" if r["p7_8d_compliant"] else "  [BELOW 2x -- "
                             "power measurement by design]"))
        _save(PROGRESS, prog)

        for arm in arms:
            for k in FACTORS:
                if k == 2.0:
                    continue
                mod, planted = plants[(arm, k)]
                for i in range(per_cell):
                    tag = f"curve_{arm}_{k}x_{i:03d}"
                    if tag in prog["cells"]:
                        print(f"[{tag}] checkpointed, skipping")
                        continue
                    rng = np.random.default_rng(np.random.SeedSequence(
                        [CURVE_SEED, STREAM_ID[(arm, k)], i]))
                    t = time.time()
                    rec = G.run_catalog(tag, prep, cfg, rng, mod, planted,
                                        int(a.jobs))
                    rec.update({"arm": arm, "factor": float(k),
                                "source": "recovery_curve",
                                "seed_note": (f"SeedSequence([{CURVE_SEED}, "
                                              f"{STREAM_ID[(arm, k)]}, {i}])")})
                    prog["cells"][tag] = rec
                    _save(PROGRESS, prog)
                    pr = rec["planted_recovery"]
                    print(f"[{tag}] {'HIT ' if pr['detected'] else 'MISS'} "
                          f"p={pr['best_p_raw']:.3g} "
                          + (f"A/A0={pr['amplitude_ratio']:.3f} "
                             if pr["amplitude_ratio"] is not None
                             else "(UNRESOLVED) ")
                          + f"{time.time()-t:.0f}s", flush=True)

    # -------------------------------------------------------------- analysis --
    cells = list(prog["cells"].values())
    curve = {}
    for arm in ARMS:
        rows, amps, nn, hh = [], [], [], []
        for k in FACTORS:
            recs = [c for c in cells
                    if c.get("arm") == arm and float(c.get("factor", -1)) == k]
            if not recs:
                continue
            det = [bool(c["planted_recovery"]["detected"]) for c in recs]
            reps = prog["plant_reports"].get(f"{arm}@{k}x", [])
            amp = (float(np.mean([r["planted_amplitude"] for r in reps]))
                   if reps else None)
            fl = (float(np.mean([r["operative_floor_A_min"] for r in reps]))
                  if reps else None)
            ratios = [c["planted_recovery"]["amplitude_ratio"] for c in recs
                      if c["planted_recovery"]["amplitude_ratio"] is not None]
            rows.append({
                "factor": k, "planted_amplitude": amp,
                "operative_floor_A_min": fl,
                "n": len(det), "n_detected": int(sum(det)),
                "recovery_rate": sum(det) / len(det),
                "wilson95": wilson(sum(det), len(det)),
                "predicted_power": predicted_power(k),
                "p7_8d_compliant": bool(k >= floors.PLANT_FACTOR),
                "amplitude_ratio_median": (float(np.median(ratios))
                                           if ratios else None),
                "amplitude_ratio_n": len(ratios),
                "source": sorted({c.get("source", "?") for c in recs}),
            })
            if amp is not None:
                amps.append(amp)
                nn.append(len(det))
                hh.append(int(sum(det)))
        if not rows:
            continue
        fit = fit_effective_floor(amps, nn, hh) if amps else {}
        a_formula = rows[0]["operative_floor_A_min"] if rows[0][
            "operative_floor_A_min"] else None
        entry = {"label": ARMS[arm]["label"], "cells": rows,
                 "a_min_formula": a_formula, "fit": fit}
        if fit.get("a_eff") and a_formula:
            entry["a_eff_over_a_formula"] = fit["a_eff"] / a_formula
            entry["fitted_A80"] = amplitude_at_power(fit["a_eff"], 0.80)
            entry["fitted_A95"] = amplitude_at_power(fit["a_eff"], 0.95)
            entry["fitted_A80_in_floor_units"] = entry["fitted_A80"] / a_formula
            entry["fitted_A95_in_floor_units"] = entry["fitted_A95"] / a_formula
        pts = [(r["planted_amplitude"], r["recovery_rate"]) for r in rows
               if r["planted_amplitude"] is not None]
        entry["raw_A80"] = interpolate_threshold(pts, 0.80)
        entry["raw_A95"] = interpolate_threshold(pts, 0.95)
        curve[arm] = entry

    # ---- the §P7-15 verdict: (a) labelling, (b) formula, or (c) mixed -------
    verdict = classify(curve)

    out = {
        "object": "§P7-15 recovery-versus-amplitude curve",
        "priced": 0,
        "pricing_citation": ("HYPOTHESIS_LEDGER.md §P7-2(a): F9-19 and the "
                             "zero-priced diagnostics make no rejection about "
                             "the Earth; §P7-1(c) / §P7-15 assign this sweep "
                             "at price 0. Simulated catalogs, sandboxed "
                             "ledger, no EXPLORE_COUNT line."),
        "generated_at": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "engine_version": __version__,
        "config_hash": cfg_hash, "config": cfg,
        "predicted_curve_source": ("§P7-12(a) floor formula inverted: "
                                   "power(k) = Phi(k*(z_alpha+z_0.80) - "
                                   "z_alpha); reproduces §P7-15's published "
                                   "0.800 / 0.9992 / 1.0000 at 1/1.5/2x"),
        "alpha": G.BASE_ALPHA, "vif": G.BASE_VIF,
        "p7_8d_departure": ("The 1.0x and 1.5x cells sit BELOW the §P7-8(d) 2x "
                            "plant line, deliberately: the curve cannot be "
                            "measured without visiting the amplitudes whose "
                            "power is in question. Misses there are POWER "
                            "measurements. No G-M1 verdict is taken from this "
                            "module; G-M1 was decided at >= 2x in "
                            "results_gate_r1.json."),
        "verdict": verdict,
        "curve": curve,
        "plant_reports": prog.get("plant_reports", {}),
        "per_catalog": cells,
        "sandbox": {"sessions": G.SANDBOX_SESSIONS.replace("\\", "/"),
                    "ledger": G.SANDBOX_LEDGER.replace("\\", "/"),
                    "real_ledger_untouched": G.REAL_LEDGER.replace("\\", "/")},
    }
    _save(RESULTS, out)

    print("-" * 78)
    for arm, e in curve.items():
        print(f"arm ({arm}) {e['label']}")
        print(f"  A_min formula = {e['a_min_formula']:.4f}"
              if e["a_min_formula"] else "  A_min formula = n/a")
        for r in e["cells"]:
            print(f"    {r['factor']:>4}x  A={r['planted_amplitude']:.4f}  "
                  f"{r['n_detected']:>2}/{r['n']:<2} = {r['recovery_rate']:>5.0%}"
                  f"  [{r['wilson95'][0]:.2f},{r['wilson95'][1]:.2f}]   "
                  f"predicted {r['predicted_power']:.4f}")
        if e.get("a_eff_over_a_formula"):
            print(f"  A_eff = {e['fit']['a_eff']:.4f} = "
                  f"{e['a_eff_over_a_formula']:.2f}x the formula's floor; "
                  f"A80 = {e['fitted_A80']:.4f}, A95 = {e['fitted_A95']:.4f}")
        else:
            print(f"  A_eff: {e['fit'].get('bound', 'not estimable')}")
    print("-" * 78)
    print(f"§P7-15 VERDICT: ({verdict['choice']}) {verdict['headline']}")
    print(f"written -> {RESULTS}")
    return 0


def classify(curve):
    """(a) plant labelling, (b) formula anti-conservative, or (c) mixed.

    The two hypotheses make DIFFERENT predictions about a quantity this sweep
    measures directly, so the classification is a reading of evidence rather
    than a judgement call:

      * (a) says the plants were not where they were labelled. The direct test
        is the RECOVERED/PLANTED amplitude ratio, which the pipeline reports on
        every global arm: if the plants really sat at ~1x, the ratio would be
        ~0.5, not ~1.0.
      * (b) says the plants were where they were labelled and the pipeline
        needed more amplitude than the formula claims. Its signature is
        A_eff / A_min_formula > 1 with the ratio near 1.0.
    """
    per_arm = {}
    for arm, e in curve.items():
        ratios = [r["amplitude_ratio_median"] for r in e["cells"]
                  if r["amplitude_ratio_median"] is not None]
        med = float(np.median(ratios)) if ratios else None
        infl = e.get("a_eff_over_a_formula")
        labelling_off = bool(med is not None and not (0.8 <= med <= 1.2))
        formula_off = bool(infl is not None and infl >= 1.25)
        per_arm[arm] = {
            "amplitude_ratio_median": med,
            "a_eff_over_a_formula": infl,
            "plants_delivered_as_labelled": (None if med is None
                                             else not labelling_off),
            "formula_anti_conservative": (None if infl is None else formula_off),
            "note": ("amplitude UNRESOLVED by §P7-1(d) on this arm; delivery is "
                     "asserted from the plant construction, not measured"
                     if med is None else None),
        }
    off = [a for a, v in per_arm.items() if v["formula_anti_conservative"]]
    mis = [a for a, v in per_arm.items()
           if v["plants_delivered_as_labelled"] is False]
    if mis and not off:
        choice, head = "a", ("plant labelling: the recovered/planted amplitude "
                             "ratio is off, the formula is not implicated")
    elif off and not mis:
        choice = "b" if len(off) == len(
            [a for a, v in per_arm.items()
             if v["formula_anti_conservative"] is not None]) else "c"
        head = (f"the floor formula is ANTI-CONSERVATIVE at arm(s) "
                f"{', '.join(off)}; plants were delivered as labelled")
    elif off and mis:
        choice, head = "c", "mixed: both effects present"
    else:
        choice, head = "none", ("neither: plants delivered as labelled and no "
                                "arm shows an inflated effective floor")
    return {"choice": choice, "headline": head, "per_arm": per_arm,
            "arms_with_inflated_floor": off,
            "arms_with_mislabelled_plants": mis}


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    sys.exit(main())
