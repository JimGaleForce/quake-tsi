"""WHY DOES A REGION RISE ABOVE THE FOLD? Regress the per-region z on physical
covariates, across every scan already run. PRICED.

THE STRATEGIC INVERSION, which is Jim's and which is the reason this exists.

Every arm so far asked "IS there a tidal effect", pooled across regions, and returned
null. **A field with a ZERO GLOBAL MEAN can still be strongly STRUCTURED.** Our pooled
nulls are entirely compatible with tidal susceptibility varying enormously from place to
place and averaging to nothing. So the question changes:

    not  "does the tide trigger earthquakes"
    but  "what predicts WHERE the tide triggers earthquakes"

and the regional scatter stops being noise and becomes the dependent variable. If one of
fifty factors governs it, finding that factor is worth more than any single regional
result, because it is the thing that would transfer to a region nobody has scored.

WHY THIS IS NOT A TEST WE HAVE ALREADY RUN. Cochran's Q asks whether the between-region
variance exceeds sampling noise. At 11 regions that is a WEAK test, and it uses only the
SPREAD. A regression against a declared covariate uses the ORDERING, which carries far
more information: eleven regions ranked correctly by a physical variable is a strong
statement even when their spread is unremarkable. Q said "no heterogeneity" in the
declustered arms; that is emphatically NOT the same as "the variation is unstructured",
and nobody has asked the second question.

---------------------------------------------------------------------------
THE DECLARATION. Fixed before the first correlation is computed.
---------------------------------------------------------------------------

COVARIATES, eight, all computable from what is on disk, and each with a DIRECTION
declared in advance so a sign flip counts against it:

  C1 mean |latitude|              -- the dwell artifact is a strong monotone function
                                     of latitude, so this is the ARTIFACT covariate and
                                     is included to be ruled OUT, not in
  C2 body-tide amplitude          -- sd of areal strain at the region's own sites. More
                                     forcing should mean more effect. DIRECTION: +
  C3 bearing dwell concentration  -- axial R2 of the tidal bearing. This is a
                                     DETECTABILITY covariate, not a physical one: a
                                     concentrated dwell destroys dynamic range. Included
                                     so it can be separated from physics. DIRECTION: -
  C4 mean event depth             -- shallower should be more tidally sensitive
                                     (lower effective normal stress). DIRECTION: -
  C5 mean interface dip (Slab2)   -- geometry; no strong prior, declared two-sided
  C6 subduction fraction          -- fraction of events with a Slab2 interface within
                                     20 km. Ocean-loaded shallow thrusts are where the
                                     literature's largest effects live. DIRECTION: +
  C7 declustered event rate       -- productivity/loading proxy. DIRECTION: two-sided
  C8 mean magnitude               -- population character. DIRECTION: two-sided

RESPONSES: the per-region z of every statistic in the two DECLUSTERED PRIMARY arms
already committed -- 16 geographic statistics and 13 fault-relative ones.

STATISTIC: Spearman rank correlation between covariate and per-region z. Rank, not
Pearson, because eleven points and one outlier region would otherwise drive everything.

NULL: PERMUTE THE REGION LABELS. This preserves both marginal distributions exactly --
the covariate values and the z values are unchanged, only their pairing is destroyed --
so it tests exactly the ordering claim and nothing else. 20,000 permutations.

MULTIPLICITY: max |rho| across every covariate x statistic pair, calibrated against the
same permutation ensemble. That is the correct correction and it is why this is priced.

WHAT A HIT WOULD AND WOULD NOT MEAN. A hit on C1 or C3 is an ARTIFACT result and says
the scan's regional scatter is instrumental. A hit on C2, C4 or C6 is a MECHANISM HINT
and is the thing worth chasing. A hit on nothing says the regional variation is
unstructured with respect to everything we can currently measure, which is itself worth
knowing and bounds the "some regions use this" hypothesis.

NOT A PROMOTION. A survivor here is a candidate for Popper.
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
from engine import slab2 as SL
from engine import tidal_tensor as TT

OUT_JSON = HERE / "results_region_covariates.json"
N_PERM = 20000
RNG_SEED = 20260822

COVARIATE_DIRECTION = {
    "C1_abs_latitude": "artifact-covariate, two-sided, included to be RULED OUT",
    "C2_bodytide_amplitude": "+",
    "C3_bearing_dwell_R2": "- (detectability, not physics)",
    "C4_mean_depth_km": "-",
    "C5_mean_interface_dip": "two-sided",
    "C6_subduction_fraction": "+",
    "C7_rate_per_year": "two-sided",
    "C8_mean_magnitude": "two-sided",
}


def region_covariates():
    """Eight declared covariates per region, computed from what is on disk."""
    cut = W.explore_cutoff()
    rows = {}
    for path in sorted(glob.glob(str(HERE / "data" / "comcat_world" / "*.csv"))):
        name = os.path.basename(path)[:-4]
        t, la, lo, dp, mg = [], [], [], [], []
        for r in csv.DictReader(open(path, newline="", encoding="utf-8")):
            try:
                m = float(r["mag"])
                if m < W.MAG_MIN:
                    continue
                ts = _dt.datetime.fromisoformat(r["time"].replace("Z", "+00:00"))
                if ts >= cut:
                    continue
                d = float(r["depth"])
            except (ValueError, TypeError, KeyError):
                continue
            t.append(ts.timestamp() / 86400.0)
            la.append(float(r["latitude"]))
            lo.append(float(r["longitude"]))
            dp.append(d)
            mg.append(m)
        if len(t) < 30:
            continue
        t, la, lo, dp, mg = (np.asarray(t), np.asarray(la), np.asarray(lo),
                             np.asarray(dp), np.asarray(mg))
        k = W.decluster(t, la, lo, mg)
        t, la, lo, dp, mg = t[k], la[k], lo[k], dp[k], mg[k]
        if t.size < 30:
            continue

        # tidal quantities at the region's own centroid-ish sample of sites
        sel = np.linspace(0, t.size - 1, min(24, t.size)).astype(int)
        amps, r2s = [], []
        for i in sel:
            tt = np.arange(0.0, 30.0, 2.0 / 1440.0)
            jd = E.julian_day_at(_dt.datetime(2010, 1, 1), tt)
            st = TT.stress_tensor(jd, la[i], lo[i], 0.0)
            amps.append(float(np.std(st["areal_strain"])))
            b = TT.principal_bearing(st["e_NN"], st["e_EE"], st["e_NE"])
            two = 2.0 * np.radians(np.mod(b, 180.0))
            r2s.append(float(np.hypot(np.mean(np.cos(two)), np.mean(np.sin(two)))))

        a = SL.assign(la, lo, dp)
        sub = a["assigned"] & (a["depth_misfit_km"] <= 20.0)
        dip = (float(np.nanmean(a["dip_deg"][sub])) if np.any(sub) else float("nan"))
        span_yr = (t.max() - t.min()) / 365.2425

        rows[name] = {
            "C1_abs_latitude": float(np.mean(np.abs(la))),
            "C2_bodytide_amplitude": float(np.mean(amps)),
            "C3_bearing_dwell_R2": float(np.mean(r2s)),
            "C4_mean_depth_km": float(np.mean(dp)),
            "C5_mean_interface_dip": dip,
            "C6_subduction_fraction": float(np.mean(sub)),
            "C7_rate_per_year": float(t.size / max(span_yr, 1e-9)),
            "C8_mean_magnitude": float(np.mean(mg)),
            "_n": int(t.size),
        }
    return rows


def spearman(x, y):
    def rank(v):
        o = np.argsort(np.argsort(v))
        return o.astype(float)
    rx, ry = rank(np.asarray(x, float)), rank(np.asarray(y, float))
    rx -= rx.mean(); ry -= ry.mean()
    d = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / d) if d > 0 else 0.0


def main():
    rng = np.random.default_rng(RNG_SEED)
    cov = region_covariates()
    print("covariates computed for %d regions" % len(cov), flush=True)

    responses = {}
    for label, fname, statkey in (
        ("geographic", "results_world_harmonics_declustered.json", "per_statistic"),
        ("faultrel", "results_world_faultrel_declustered.json", "per_statistic"),
    ):
        d = json.load(open(HERE / fname, encoding="utf-8"))
        for rn, rec in d["per_region"].items():
            for sk, sv in rec[statkey].items():
                responses.setdefault("%s:%s" % (label, sk), {})[rn] = sv["z"]

    cov_names = [c for c in COVARIATE_DIRECTION]
    pairs, obs = [], []
    for cn in cov_names:
        for sn, zmap in responses.items():
            regions = sorted(set(zmap) & set(cov))
            regions = [r for r in regions if np.isfinite(cov[r][cn])]
            if len(regions) < 6:
                continue
            x = np.array([cov[r][cn] for r in regions])
            y = np.array([zmap[r] for r in regions])
            pairs.append((cn, sn, regions))
            obs.append(spearman(x, y))
    obs = np.asarray(obs)
    print("declared pairs: %d covariates x responses = %d" % (len(cov_names), obs.size),
          flush=True)

    # permutation null: shuffle region labels, recompute every pair on the same shuffle
    null_max = np.empty(N_PERM)
    for p in range(N_PERM):
        perm_cache = {}
        best = 0.0
        for (cn, sn, regions), _o in zip(pairs, obs):
            key = tuple(regions)
            if key not in perm_cache:
                perm_cache[key] = rng.permutation(len(regions))
            pi = perm_cache[key]
            x = np.array([cov[r][cn] for r in regions])
            y = np.array([responses[sn][r] for r in regions])[pi]
            best = max(best, abs(spearman(x, y)))
        null_max[p] = best
        if (p + 1) % 4000 == 0:
            print("  permutations %d/%d" % (p + 1, N_PERM), flush=True)

    obs_max = float(np.max(np.abs(obs)))
    p_max = (int(np.sum(null_max >= obs_max)) + 1) / (N_PERM + 1)
    iw = int(np.argmax(np.abs(obs)))

    ranked = sorted(range(obs.size), key=lambda i: -abs(obs[i]))[:15]
    out = {
        "arm": "regional covariate regression: WHY does a region rise above the fold",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "PRICED": True, "n_declared_tests": int(obs.size),
        "declaration": {
            "covariates": COVARIATE_DIRECTION,
            "statistic": "Spearman rank correlation, region-level",
            "null": "permute REGION LABELS; preserves both marginals exactly and "
                    "destroys only the pairing, so it tests the ordering claim alone",
            "n_permutations": N_PERM, "rng_seed": RNG_SEED,
            "why_not_cochran_Q": ("Q tests whether between-region VARIANCE exceeds "
                                  "noise and is weak at 11 regions; this tests whether "
                                  "the ORDERING is organised by a physical covariate, "
                                  "which uses far more of the information"),
        },
        "region_covariates": cov,
        "max_statistic": {
            "observed_max_abs_rho": obs_max,
            "null_max_p95": float(np.quantile(null_max, 0.95)),
            "p": float(p_max),
            "where": {"covariate": pairs[iw][0], "response": pairs[iw][1],
                      "rho": float(obs[iw]), "n_regions": len(pairs[iw][2])},
        },
        "top_15_pairs": [{"covariate": pairs[i][0], "response": pairs[i][1],
                          "rho": float(obs[i]), "n_regions": len(pairs[i][2])}
                         for i in ranked],
        "interpretation": {
            "artifact_covariates": ["C1_abs_latitude", "C3_bearing_dwell_R2"],
            "mechanism_covariates": ["C2_bodytide_amplitude", "C4_mean_depth_km",
                                     "C6_subduction_fraction"],
            "note": ("a hit on an artifact covariate says the regional scatter is "
                     "instrumental; a hit on a mechanism covariate is the thing worth "
                     "chasing; a hit on nothing bounds the 'some regions use this' "
                     "hypothesis against everything we can currently measure"),
        },
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    print("\n" + "=" * 76)
    print("MAX |rho| = %.3f  (%s vs %s, %d regions)"
          % (obs_max, pairs[iw][0], pairs[iw][1], len(pairs[iw][2])))
    print("null 95th = %.3f   PERMUTATION p = %.4f"
          % (np.quantile(null_max, 0.95), p_max))
    print("\ntop pairs by |rho|:")
    for i in ranked[:10]:
        tag = " [ARTIFACT covariate]" if pairs[i][0] in ("C1_abs_latitude",
                                                         "C3_bearing_dwell_R2") else ""
        print("  %-24s %-34s rho = %+.3f%s"
              % (pairs[i][0], pairs[i][1], obs[i], tag))
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
