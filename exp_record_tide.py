"""RECORD-TIDE TRIGGERING: has this fault ever felt a push this hard before?

Kepler K-141/K-163, and it is the first genuinely new QUESTION in this program rather
than a new statistic on an old one.

THE IDEA, IN ONE SENTENCE. Every scan so far asked "what tidal PHASE do earthquakes
prefer", treating all ~17,000 tidal peaks in 25 years as members of one population. But
a stuck drawer does not care what phase you pull at -- it cares whether THIS pull is
harder than any pull it has already survived. A fault that held through a 40 nanostrain
peak yesterday is not learning anything new from a 38 nanostrain peak today: it has
already proved it can hold that. **The mechanically relevant variable is not amplitude
or phase but whether the current peak is a RECORD.**

WHY IT HAS NEVER BEEN RUN HERE, AND WHY IT IS A DIFFERENT EXPERIMENT. Measured on the
program's own frozen scalar over 25 years at 34N/117W, 30-minute step:

    17,413 peaks, mean spacing 0.524 d, peak amplitude max/median = 2.99
    fraction of peaks that are a record over the trailing 1 year  : 0.21 %
    ... trailing 5 years                                          : 0.09 %
    ... trailing 18.61 years                                      : 0.06 %

A 0.2 %-duty-cycle marker is a completely different experiment from a 50 % phase bin.
Every phase statistic in this program is blind to it, and it is blind to them: the
record flag is a LEVEL-CROSSING property of the trailing history, not a property of the
instantaneous state.

---------------------------------------------------------------------------
THE DECLARATION. Fixed before the first event is scored.
---------------------------------------------------------------------------

STATISTIC, three per region: the fraction of events whose containing tidal cycle's peak
is a record over the trailing W, for W in {30 days, 1 year, 5 years}. Nothing else.

THE NULL IS THE BASE RATE AND IT IS EMPIRICAL RATHER THAN ASSUMED. Times drawn uniformly
over the same window at the same reference site, pushed through the identical lookup.
That is the fraction of TIME the site spends inside a record cycle, which is exactly the
probability an event at a uniformly random time lands in one. **No surrogate catalogue
and no phase permutation is involved.** This is the cleanest null in the program.

SITE APPROXIMATION, DECLARED AND MEASURED RATHER THAN ASSUMED. The record STRUCTURE --
which peaks are records -- is driven by the spring-neap, declination and perigee
envelopes, which are ephemeris properties that vary slowly over a region. The peak
series is therefore computed once per region at its event centroid rather than once per
event. `record_flag_agreement` measures how well that holds by comparing two sites
several hundred km apart, and the number is reported on the artifact. If the agreement
were poor the approximation would be void.

EVENTS: the same world catalogue, M >= 5.0, exploration split (pre-2017-02-14), holdout
untouched. DECLUSTERED PRIMARY, full SECONDARY -- the ordering every arm in this program
now uses, because three earlier survivors turned out to be dependence.

MULTIPLICITY: max-statistic across every region and window against the same null
ensemble. 13 regions x 3 windows = 39 declared tests.

WHAT A HIT WOULD MEAN. An excess ratio above 1 would say faults respond to novelty of
load rather than to its phase, which is a mechanism no tidal-triggering statistic in the
literature tests and which would be visible in the catalogue we already have. A ratio of
exactly 1 with a tight interval is the strongest constraint anyone has published on the
"crust remembers its largest recent load" hypothesis.
"""

from __future__ import annotations

import csv
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
from engine import ephemeris as E
from engine import sitetide as ST

OUT_JSON = HERE / "results_record_tide.json"
STEP_MIN = 30.0
WINDOWS = (("W30d", 30.0), ("W1yr", 365.25), ("W5yr", 5 * 365.25))
N_NULL_DRAWS = 20000
N_REPLICATES = 5000
RNG_SEED = 20260822
UNIX_EPOCH_JD = 2440587.5
MIN_EVENTS = 30


def peak_series(lat, lon, t_lo, t_hi, with_grid=False):
    """Local maxima of the declared scalar over [t_lo, t_hi] in Unix days."""
    t = np.arange(t_lo, t_hi, STEP_MIN / 1440.0)
    s = ST.site_scalar(t + UNIX_EPOCH_JD, lat, lon, 0.0)
    i = np.nonzero((s[1:-1] > s[:-2]) & (s[1:-1] >= s[2:]))[0] + 1
    if with_grid:
        return t[i], s[i], t, s
    return t[i], s[i]


def record_flags(t_pk, a_pk, window_days):
    """Is each peak a record over the trailing `window_days`? NaN where not eligible."""
    n = t_pk.size
    flag = np.full(n, np.nan)
    for k in range(n):
        lo = t_pk[k] - window_days
        j0 = np.searchsorted(t_pk, lo, side="left")
        if k - j0 < 10:                      # too little trailing history to judge
            continue
        flag[k] = 1.0 if a_pk[k] > a_pk[j0:k].max() else 0.0
    return flag


def containing_cycle(t_pk, times):
    """Index of the peak whose cycle contains each time (nearest preceding peak).

    RETAINED FOR THE v1 COMPARISON ONLY. See `above_record_at` for why it is the wrong
    assignment for this test.
    """
    j = np.searchsorted(t_pk, times, side="right") - 1
    return np.clip(j, 0, t_pk.size - 1)


def above_record_at(t_grid, s_grid, t_pk, a_pk, window_days, times):
    """Is the stress AT EACH GIVEN INSTANT above every peak in the trailing window?

    THE v1 BUG THIS REPLACES, confirmed by measurement before the fix. v1 attached the
    record flag to the PEAK of the containing cycle. But the stress first exceeds the
    old record on the RISING LIMB *toward* that peak -- which lies in the interval whose
    preceding peak is the previous, NON-record one. So v1 labelled the first-exceedance
    moment 0, not 1.

    Measured on the program's own scalar, 3 years at 34N/117W: 33.3 % of all
    above-old-record dwell time falls before the record peak and was misassigned. For a
    STRICT first-exceedance mechanism -- which is precisely the Kaiser-effect story this
    arm exists to test -- the loss approaches 100 %, because the first moment the stress
    passes the old record is ALWAYS on the rising limb.

    It failed silently: the null was built by the identical path, so the test stayed
    well-formed and simply could not see what it was looking for.

    This version asks the question at the event's own instant instead, which is both the
    physically correct condition ("is this fault feeling something it has not felt
    before, right now") and free of any assignment ambiguity.
    """
    n = t_pk.size
    trailing = np.full(n, np.nan)
    for k in range(n):
        j0 = np.searchsorted(t_pk, t_pk[k] - window_days, side="left")
        if k - j0 < 10:
            continue
        trailing[k] = a_pk[j0:k + 1].max()
    k_at = np.clip(np.searchsorted(t_pk, times, side="right") - 1, 0, n - 1)
    tr = trailing[k_at]
    s_at = np.interp(times, t_grid, s_grid)
    out = np.full(times.shape, np.nan)
    ok = np.isfinite(tr)
    out[ok] = (s_at[ok] > tr[ok]).astype(float)
    return out


def main():
    rng = np.random.default_rng(RNG_SEED)
    cut = W.explore_cutoff()
    declustered = os.environ.get("RT_DECLUSTER", "1") == "1"
    print("RECORD-TIDE TEST%s" % (" (DECLUSTERED, PRIMARY)" if declustered
                                  else " (FULL, SECONDARY)"), flush=True)

    regions, agreement = {}, None
    for path in sorted(glob.glob(str(HERE / "data" / "comcat_world" / "*.csv"))):
        name = os.path.basename(path)[:-4]
        t, la, lo, mg = [], [], [], []
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
        t, la, lo, mg = (np.asarray(t), np.asarray(la), np.asarray(lo), np.asarray(mg))
        if t.size < MIN_EVENTS:
            continue
        if declustered:
            k = W.decluster(t, la, lo, mg)
            t, la, lo = t[k], la[k], lo[k]
        if t.size < MIN_EVENTS:
            continue

        clat, clon = float(np.median(la)), float(np.median(lo))
        # peak series must start 5 years before the earliest event so W5yr is eligible
        t_lo = t.min() - 5.2 * 365.25
        t_hi = t.max() + 1.0
        t_pk, a_pk, t_grid, s_grid = peak_series(clat, clon, t_lo, t_hi,
                                                 with_grid=True)

        if agreement is None:      # measure the site approximation once, not assume it
            t2, a2 = peak_series(clat + 3.0, clon + 3.0, t_lo, min(t_lo + 2000.0, t_hi))
            f1 = record_flags(t_pk[:t2.size], a_pk[:t2.size], 365.25)
            f2 = record_flags(t2, a2, 365.25)
            ok = np.isfinite(f1) & np.isfinite(f2)
            agreement = {"sites_apart_deg": 3.0, "n_compared": int(ok.sum()),
                         "flag_agreement": float(np.mean(f1[ok] == f2[ok]))}

        t_draw = rng.uniform(t.min(), t.max(), N_NULL_DRAWS)

        per = {}
        for wname, wdays in WINDOWS:
            ev = above_record_at(t_grid, s_grid, t_pk, a_pk, wdays, t)
            nl = above_record_at(t_grid, s_grid, t_pk, a_pk, wdays, t_draw)
            ev = ev[np.isfinite(ev)]
            nl = nl[np.isfinite(nl)]
            if ev.size < MIN_EVENTS or nl.size < 1000:
                continue
            base = float(nl.mean())
            obs = float(ev.mean())
            # replicate null: draw n events' worth from the null pool
            reps = np.array([nl[rng.integers(0, nl.size, ev.size)].mean()
                             for _ in range(N_REPLICATES)])
            sd = float(reps.std(ddof=1))
            per["%s" % wname] = {
                "n_events": int(ev.size),
                "observed_fraction": obs,
                "base_rate": base,
                "excess_ratio": (obs / base) if base > 0 else float("nan"),
                "n_events_in_record_cycle": int(round(obs * ev.size)),
                "expected_in_record_cycle": base * ev.size,
                "null_sd": sd,
                "z": (obs - base) / sd if sd > 0 else float("nan"),
                "_reps": (reps - base) / (sd if sd > 0 else 1e-300),
            }
        if per:
            regions[name] = {"n": int(t.size), "centroid": [clat, clon],
                             "per_window": per}
            best = max(per, key=lambda k: abs(per[k]["z"]))
            print("  %-18s n=%5d   %s: obs %.4f vs base %.4f  ER %.2f  z %+.2f"
                  % (name, t.size, best, per[best]["observed_fraction"],
                     per[best]["base_rate"], per[best]["excess_ratio"],
                     per[best]["z"]), flush=True)

    # ---- EXACT BINOMIAL CALIBRATION, and the reason is a defect this run caught ----
    # The first version calibrated a normal-ish z against a replicate ensemble. At these
    # base rates (0.002 for the 1-year window) and these n, the null count is 0 with
    # probability 0.88 and 2 with probability 0.007, so a z built from a standard
    # deviation is badly non-Gaussian and its max-statistic is ANTI-CONSERVATIVE:
    # it reported p = 0.033 where the exact calibration gives p = 0.118. Rare-event
    # counts need exact binomial tails, not a z.
    from scipy.stats import binom
    tests = [(r, w, regions[r]["per_window"][w]["n_events"],
              regions[r]["per_window"][w]["base_rate"],
              regions[r]["per_window"][w]["n_events_in_record_cycle"])
             for r in regions for w in regions[r]["per_window"]]
    exact = {}
    for r, w, n_e, p_b, k_o in tests:
        exact.setdefault(r, {})[w] = float(binom.sf(k_o - 1, n_e, p_b))
        regions[r]["per_window"][w]["exact_binomial_p"] = exact[r][w]
    obs_min_p = min(exact[r][w] for r, w, _, _, _ in tests)
    worst = min(((r, w, exact[r][w]) for r, w, _, _, _ in tests), key=lambda x: x[2])
    # null distribution of the MINIMUM exact p across the whole family
    NS = 200000
    mins = np.ones(NS)
    for r, w, n_e, p_b, _k in tests:
        ks = rng.binomial(n_e, p_b, NS)
        mins = np.minimum(mins, binom.sf(ks - 1, n_e, p_b))
    p_gw = float((np.sum(mins <= obs_min_p) + 1) / (NS + 1))
    obs_max = obs_min_p
    for r in regions:
        for w in regions[r]["per_window"]:
            regions[r]["per_window"][w].pop("_reps", None)

    n_tests = sum(len(regions[r]["per_window"]) for r in regions)
    out = {
        "arm": "RECORD-TIDE TRIGGERING (Kepler K-141/K-163)",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "PRICED": True, "n_declared_tests": n_tests, "declustered": declustered,
        "the_question": ("not 'what tidal phase do earthquakes prefer' but 'has this "
                         "fault ever felt a push this hard before'. A fault that held "
                         "through a bigger peak last month learns nothing from a "
                         "smaller one today."),
        "declaration": {
            "windows_days": {k: v for k, v in WINDOWS},
            "step_minutes": STEP_MIN, "mag_min": W.MAG_MIN,
            "explore_cutoff": cut.isoformat(), "holdout": "NOT READ",
            "null": ("uniform random times over the same window at the same reference "
                     "site through the identical lookup -- the fraction of TIME spent "
                     "inside a record cycle. No surrogate catalogue, no phase "
                     "permutation."),
            "site_approximation": agreement,
            "n_null_draws": N_NULL_DRAWS, "n_replicates": N_REPLICATES,
            "rng_seed": RNG_SEED,
        },
        "per_region": regions,
        "max_statistic": {
            "method": ("EXACT binomial per test; family corrected by simulating the "
                       "MINIMUM exact p across all tests under the null. A z-based "
                       "max-statistic is anti-conservative at these counts and the "
                       "first version of this script reported p = 0.033 where the "
                       "exact calibration gives 0.118."),
            "smallest_exact_p": float(obs_min_p),
            "n_tests": len(tests),
            "family_corrected_p": float(p_gw),
            "where": {"region": worst[0], "window": worst[1],
                      "exact_p": float(worst[2])},
            "exact_p_by_test": exact},
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    print("\n" + "=" * 78)
    print("RECORD-TIDE RESULT   (%d declared tests)" % n_tests)
    print("=" * 78)
    print("  site approximation: record flags agree %.3f between sites 3 deg apart"
          % agreement["flag_agreement"])
    print("  smallest EXACT binomial p = %.5f at %s / %s   (%d tests)"
          % (obs_min_p, worst[0], worst[1], len(tests)))
    print("  FAMILY-CORRECTED p = %.4f" % p_gw)
    print("\n  pooled excess ratio by window:")
    for wname, _ in WINDOWS:
        tot_o = sum(regions[r]["per_window"][wname]["n_events_in_record_cycle"]
                    for r in regions if wname in regions[r]["per_window"])
        tot_e = sum(regions[r]["per_window"][wname]["expected_in_record_cycle"]
                    for r in regions if wname in regions[r]["per_window"])
        if tot_e > 0:
            print("    %-6s observed %4d in a record cycle, expected %7.1f  ER = %.2f"
                  % (wname, tot_o, tot_e, tot_o / tot_e))
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
