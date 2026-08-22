"""IS IT US OR IS IT THE WORLD? A non-circular calibration of the null pipeline. Priced 0.

The declustered per-region z's are under-dispersed: variance 0.641, and once the cells'
mutual correlation is measured rather than assumed (`exp_dispersion_dof.py`) that is 2.15
sigma, p = 0.032. Small, but it points somewhere serious, because there are only two
explanations and they lead opposite ways:

  US    -- the pipeline inflates null standard deviations, so every z in the program is
           too small, every test under-powered, and every published bound too loose. If
           this is it, stop hunting physics and go fix the instrument.
  WORLD -- the pipeline is exact and the real catalogue genuinely sits closer to its own
           null mean than a random null draw does. Odd, but a fact about the data.

THE TEST, AND WHY IT IS NOT CIRCULAR. The obvious check -- take one null replicate as a
pseudo-catalogue and score it against the other replicates -- is worthless, because the
null columns are standardised using those same replicates, so the answer is 1.000 by
construction and proves nothing. Instead this builds a SECOND, INDEPENDENT null ensemble
from a different seed and different offsets, then scores ensemble B's draws as if they
were the observed catalogue against ensemble A's mean and standard deviation.

Under a correct pipeline those pseudo-z's are marginally N(0,1) and their variance across
cells is 1.000. If the pipeline inflates null spread, the pseudo-z variance comes back at
roughly 0.64 -- the same deficit the real data shows -- and the defect is ours. If it comes
back at 1.000, the pipeline is exact and the deficit belongs to the catalogue.

This is the same discriminator logic the program used on the latitude artifact: find a
comparison where the two explanations predict DIFFERENT numbers, rather than arguing about
which is more plausible.
"""

from __future__ import annotations

import datetime as _dt
import glob
import json
import math
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import exp_world_harmonics as W

OUT_JSON = HERE / "results_null_calibration.json"
N_PSEUDO = 200          # independent pseudo-catalogues scored per cell


def ensemble(t, lat, lon, rng, n_rep):
    """One null ensemble: per-event pools, then n_rep independent draws of the battery."""
    n = t.size
    off = rng.uniform(-W.NULL_HALF_WINDOW_DAYS, W.NULL_HALF_WINDOW_DAYS,
                      size=(n, W.S_NULL_PER_EVENT))
    nf = W.features((t[:, None] + off).ravel(),
                    np.repeat(lat, W.S_NULL_PER_EVENT),
                    np.repeat(lon, W.S_NULL_PER_EVENT))
    ang = {a: nf[a].reshape(n, W.S_NULL_PER_EVENT) for a in W.ANGLES}
    flg = {f: nf[f].reshape(n, W.S_NULL_PER_EVENT) for f in W.FLAGS}
    idx = rng.integers(0, W.S_NULL_PER_EVENT, size=(n, n_rep))
    return W.battery(ang, flg, idx=idx)


def main():
    pseudo_z = []
    per_region = {}
    for path in sorted(glob.glob(str(HERE / "data" / "comcat_world" / "*.csv"))):
        name = os.path.basename(path)[:-4]
        t, la, lo, _h, _f = W.load_region(path, declustered=True)
        if t.size < 30:
            continue
        # A: the reference ensemble.  B: an INDEPENDENT ensemble, different seed AND
        # different offsets, whose draws stand in for the observed catalogue.
        A = ensemble(t, la, lo, np.random.default_rng(W.RNG_SEED), W.N_REPLICATES)
        B = ensemble(t, la, lo, np.random.default_rng(W.RNG_SEED + 999983), N_PSEUDO)
        zs = []
        for k in W.STATS:
            a = np.asarray(A[k], dtype=np.float64)
            mu, sd = float(a.mean()), float(a.std(ddof=1))
            if not np.isfinite(sd) or sd <= 0:
                continue
            zs.append((np.asarray(B[k], dtype=np.float64) - mu) / sd)
        zs = np.stack(zs)                       # [cells, N_PSEUDO]
        pseudo_z.append(zs)
        per_region[name] = {"n": int(t.size),
                            "pseudo_z_variance": float(zs.var(ddof=1))}
        print("  %-18s n=%5d   pseudo-z variance %.4f"
              % (name, t.size, zs.var(ddof=1)), flush=True)

    Z = np.concatenate(pseudo_z, axis=0)        # [all cells, N_PSEUDO]
    # variance across CELLS within each pseudo-catalogue, then averaged: this is exactly
    # the quantity the dispersion gate measures on the real catalogue
    per_draw = Z.var(axis=0, ddof=1)
    v_mean = float(per_draw.mean())
    v_se = float(per_draw.std(ddof=1) / math.sqrt(per_draw.size))

    out = {
        "arm": "non-circular calibration of the dwell-null pipeline",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 0,
        "design": ("ensemble B is built with a different seed and different offsets from "
                   "ensemble A, then scored against A's mean and sd; a correct pipeline "
                   "returns pseudo-z variance 1.000"),
        "n_cells": int(Z.shape[0]), "n_pseudo_catalogues": int(Z.shape[1]),
        "pseudo_z_variance_mean": v_mean,
        "pseudo_z_variance_se": v_se,
        "pseudo_z_mean": float(Z.mean()),
        "real_catalogue_variance_for_comparison": 0.641,
        "per_region": per_region,
    }
    lo_, hi_ = v_mean - 2.5 * v_se, v_mean + 2.5 * v_se
    if lo_ <= 1.0 <= hi_:
        out["verdict"] = (
            "PIPELINE IS CALIBRATED. Pseudo-z variance %.4f +/- %.4f brackets 1.000, so "
            "the dwell-null construction does not inflate null spread. The real "
            "catalogue's 0.641 is therefore a property of the DATA, not of the code -- a "
            "2.15 sigma curiosity, not a broken instrument, and not something to build "
            "on at this strength." % (v_mean, v_se))
    elif hi_ < 1.0:
        out["verdict"] = (
            "PIPELINE INFLATES NULL SPREAD. Pseudo-z variance %.4f +/- %.4f is BELOW "
            "1.000 on pure null data, so the deficit is manufactured by the code and "
            "every z in this program is too small by a factor of about %.2f. Every bound "
            "is too loose. Fix the instrument before running anything else."
            % (v_mean, v_se, 1.0 / math.sqrt(v_mean)))
    else:
        out["verdict"] = (
            "PIPELINE UNDERSTATES NULL SPREAD. Pseudo-z variance %.4f +/- %.4f is ABOVE "
            "1.000 on pure null data, which is the dangerous direction: z's are too "
            "large and the program may have been generating false positives."
            % (v_mean, v_se))

    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print("\n" + "=" * 74)
    print("pseudo-z variance on PURE NULL data: %.4f +/- %.4f  (expect 1.000)"
          % (v_mean, v_se))
    print("real catalogue for comparison:       0.641")
    print("\n" + out["verdict"])
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
