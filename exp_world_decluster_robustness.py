"""DECLUSTERING ROBUSTNESS: the world-scan result under FOUR declustering choices,
with residual clustering MEASURED rather than assumed. Priced 0 (re-scoring, not a new
battery).

WHY THIS IS A GATE AND NOT AN IMPROVEMENT. Kepler's C-5, and it is right. The world
scan reported its result under exactly TWO arms -- full and one 30-day/150-km window --
and those two disagreed violently about the only survivor, which went from z = +4.76 to
z = -0.51. That disagreement is precisely the sensitivity the program's own prior-art
record warns about:

  * Luen & Stark (2012), GJI 189:691-700: SoCal M>=3.8 declustered with
    Gardner-Knopoff FAILS the stationary-Poisson hypothesis, and Reasenberg fails it
    even at M>=4.0. Declustered catalogues are NOT Poisson, so an inference that
    assumes they are is unsafe.
  * Merton M-007.4 A3, verbatim: "report the result under >= 3 declustering choices
    including 'none', and pre-register which is primary."

A bound quoted under one window is a bound conditional on that window. This module
reports it under four, and measures how much clustering each one actually leaves.

THE FOUR CHOICES, declared:

  R0  NONE. The full catalogue. Included because M-007.4 A3 requires it and because it
      bounds the effect of the correction from the other end.
  R1  WINDOW 30 d / 150 km. The rule section P7-25 froze before any event and the
      program's standing dependence correction. PRE-REGISTERED PRIMARY.
  R2  WINDOW 60 d / 250 km. A deliberately more aggressive window, to show the
      direction the answer moves when more is removed.
  R3  GARDNER-KNOPOFF (1974) magnitude-dependent windows, the field's most-cited rule:
        L(M) = 10^(0.1238 M + 0.983) km
        T(M) = 10^(0.032 M + 2.7389) d for M >= 6.5, else 10^(0.5409 M - 0.547) d

RESIDUAL CLUSTERING IS MEASURED, NOT ASSUMED. For each rule the INDEX OF DISPERSION
(Fano factor: variance/mean of counts in fixed 30-day bins) is computed per region and
pooled. A homogeneous Poisson process gives 1.0. A value materially above 1 means the
"declustered" catalogue is still clustered and its effective N is below its nominal N,
which is exactly the Luen & Stark result and exactly what makes a z computed on it
optimistic.

WHAT IS RE-SCORED: the same 16-statistic battery, the same waveform-matched null, the
same seed. Nothing is added, so no new tests are declared and this is priced 0 -- it is
the sensitivity analysis of an already-priced scan, not a second scan.
"""

from __future__ import annotations

import datetime as _dt
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import exp_world_harmonics as W
from engine import circular_symmetry as CS

OUT_JSON = HERE / "results_world_decluster_robustness.json"
BIN_DAYS = 30.0


def gk_windows(mag):
    """Gardner-Knopoff (1974) magnitude-dependent space and time windows."""
    m = np.asarray(mag, dtype=np.float64)
    L = 10.0 ** (0.1238 * m + 0.983)
    T = np.where(m >= 6.5, 10.0 ** (0.032 * m + 2.7389),
                 10.0 ** (0.5409 * m - 0.547))
    return L, T


def decluster_window(t, lat, lon, mag, days, km):
    order = np.argsort(t)
    keep = []
    for ii in order:
        dep = False
        for jj in reversed(keep):
            if t[ii] - t[jj] > days:
                break
            if mag[jj] >= mag[ii] and W._great_circle_km(lat[jj], lon[jj],
                                                         lat[ii], lon[ii]) <= km:
                dep = True
                break
        if not dep:
            keep.append(ii)
    return np.asarray(sorted(keep), dtype=int)


def decluster_gk(t, lat, lon, mag):
    """Gardner-Knopoff: the window is set by the PRIOR (larger) event's magnitude."""
    order = np.argsort(t)
    L, T = gk_windows(mag)
    keep = []
    for ii in order:
        dep = False
        for jj in reversed(keep):
            if t[ii] - t[jj] > T[jj]:
                continue          # this prior event's window has closed for us
            if mag[jj] >= mag[ii] and W._great_circle_km(lat[jj], lon[jj],
                                                         lat[ii], lon[ii]) <= L[jj]:
                dep = True
                break
        if not dep:
            keep.append(ii)
    return np.asarray(sorted(keep), dtype=int)


def index_of_dispersion(t, bin_days=BIN_DAYS):
    """Fano factor of counts in fixed bins. 1.0 for a homogeneous Poisson process."""
    if t.size < 10:
        return float("nan")
    lo, hi = t.min(), t.max()
    nb = max(int(np.ceil((hi - lo) / bin_days)), 2)
    counts, _ = np.histogram(t, bins=nb, range=(lo, lo + nb * bin_days))
    m = counts.mean()
    return float(counts.var(ddof=1) / m) if m > 0 else float("nan")


RULES = (
    ("R0_none", None),
    ("R1_window_30d_150km", ("win", 30.0, 150.0)),
    ("R2_window_60d_250km", ("win", 60.0, 250.0)),
    ("R3_gardner_knopoff", ("gk",)),
)


def main():
    cut = W.explore_cutoff()
    paths = sorted(glob.glob(str(HERE / "data" / "comcat_world" / "*.csv")))

    raw = {}
    for path in paths:
        name = os.path.basename(path)[:-4]
        t, la, lo, mg = [], [], [], []
        import csv
        for r in csv.DictReader(open(path, newline="", encoding="utf-8")):
            try:
                m = float(r["mag"])
                if m < W.MAG_MIN:
                    continue
                ts = _dt.datetime.fromisoformat(r["time"].replace("Z", "+00:00"))
                if ts >= cut:
                    continue
            except (ValueError, TypeError, KeyError):
                continue
            t.append(ts.timestamp() / 86400.0)
            la.append(float(r["latitude"]))
            lo.append(float(r["longitude"]))
            mg.append(m)
        if t:
            raw[name] = tuple(np.asarray(x) for x in (t, la, lo, mg))

    out = {
        "arm": "declustering robustness of the world scan (Kepler C-5 gate)",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 0,
        "why_priced_zero": ("re-scoring an already-priced battery under alternative "
                            "dependence corrections; no statistic is added"),
        "mandate": ("Merton M-007.4 A3: report the result under >= 3 declustering "
                    "choices including 'none', and pre-register which is primary. "
                    "Luen & Stark (2012) GJI 189:691-700: declustered catalogues fail "
                    "the stationary-Poisson hypothesis, so a z computed on one is "
                    "optimistic by an unmeasured amount."),
        "pre_registered_primary": "R1_window_30d_150km",
        "rules": {}, "by_rule": {},
    }

    for rule_name, spec in RULES:
        rng = np.random.default_rng(W.RNG_SEED)
        regions, disp, kept = {}, {}, {}
        for name, (t, la, lo, mg) in raw.items():
            if spec is None:
                k = np.arange(t.size)
            elif spec[0] == "win":
                k = decluster_window(t, la, lo, mg, spec[1], spec[2])
            else:
                k = decluster_gk(t, la, lo, mg)
            tt, ll, oo = t[k], la[k], lo[k]
            kept[name] = {"n": int(tt.size), "n_full": int(t.size),
                          "fraction_kept": float(tt.size / max(t.size, 1))}
            disp[name] = index_of_dispersion(tt)
            r = W.run_region(name, tt, ll, oo, rng, verbose=False)
            if r is not None:
                regions[name] = r
        if not regions:
            continue
        cols = [regions[rn]["_null_cols"][k] for rn in regions for k in W.STATS]
        gw = np.max(np.abs(np.stack(cols)), axis=0)
        obs_max = max(abs(regions[rn]["per_statistic"][k]["z"])
                      for rn in regions for k in W.STATS)
        p_gw = (int(np.sum(gw >= obs_max)) + 1) / (W.N_REPLICATES + 1)
        worst = max(((rn, k, regions[rn]["per_statistic"][k]["z"])
                     for rn in regions for k in W.STATS), key=lambda x: abs(x[2]))
        comb = {k: CS.combine_regions(
            {rn: regions[rn]["per_statistic"][k]["z"] for rn in regions},
            {rn: regions[rn]["n"] for rn in regions}) for k in W.STATS}
        n_het = sum(1 for v in comb.values() if v["heterogeneous"])
        for rn in regions:
            regions[rn].pop("_null_cols", None)

        dvals = [v for v in disp.values() if np.isfinite(v)]
        rec = {
            "n_regions": len(regions),
            "n_events": int(sum(r["n"] for r in regions.values())),
            "kept": kept,
            "index_of_dispersion": {"per_region": disp,
                                    "median": float(np.median(dvals)),
                                    "max": float(np.max(dvals))},
            "max_statistic": {"observed_max_abs_z": float(obs_max),
                              "p": float(p_gw),
                              "where": {"region": worst[0], "statistic": worst[1],
                                        "z": float(worst[2])}},
            "n_heterogeneous_of_16": int(n_het),
            "pooled_z": {k: comb[k]["z_pooled"] for k in W.STATS},
            "largest_pooled_abs_z": float(max(abs(comb[k]["z_pooled"])
                                              for k in W.STATS)),
        }
        out["by_rule"][rule_name] = rec
        print("  %-22s n=%5d  disp(med)=%5.2f  max|z|=%5.2f p=%.4f  het %d/16  "
              "largest pooled |z|=%.2f"
              % (rule_name, rec["n_events"], rec["index_of_dispersion"]["median"],
                 obs_max, p_gw, n_het, rec["largest_pooled_abs_z"]), flush=True)

    ms = {k: v["max_statistic"]["p"] for k, v in out["by_rule"].items()}
    out["verdict"] = {
        "max_statistic_p_by_rule": ms,
        "any_rule_clears_0.05": bool(any(p < 0.05 for p in ms.values())),
        "primary_rule_p": ms.get("R1_window_30d_150km"),
        "reading": ("the bound is only as good as the rule it was computed under. "
                    "Reported across four so the sensitivity is visible rather than "
                    "hidden, with residual clustering measured per rule so a rule "
                    "that leaves the catalogue clustered cannot be mistaken for one "
                    "that fixed it."),
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
