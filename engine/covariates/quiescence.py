"""quiescence -- K-076 style seismic-quiescence z-score.

z = -(observed - expected) / sqrt(expected)

`observed` = count of M>=m_min within radius_km over the trailing window W days,
strictly before day t. `expected` = the same neighbourhood's EXPLORATION-period
climatological count over W days. Sign is flipped so that POSITIVE z means "quieter
than climatology" -- i.e. a positive fitted beta would mean quiescence precedes
elevated rate, which is the K-076 claim.

Caveats (v1, declared not solved):
  * `expected` is time-stationary. A cell whose catalogue completeness improved over
    the record will read as chronically "active", not quiescent.
  * The climatology used here is fitted on the exploration period only, so it is a
    legitimate training-time transform, but it is the SAME climatology the baseline
    uses -- this covariate is close to a residual of the baseline by construction.
"""

from __future__ import annotations

import numpy as np

from . import register

DEFAULTS = {"m_min": 4.5, "radius_km": 100.0, "window_days": 365}


@register("quiescence", defaults=DEFAULTS, burn_in=365,
          describe="-(obs-exp)/sqrt(exp) over trailing window vs exploration climatology")
def compute(ctx, p):
    w = int(p["window_days"])
    p.setdefault("burn_in", w)
    obs = ctx.past_counts(p["m_min"], p["radius_km"], w)

    # neighbourhood climatology from the EXPLORATION window only
    counts = ctx.day_counts(p["m_min"])
    nbr = ctx.neighbour_matrix(p["radius_km"])
    train_tot = counts[:, : ctx.n_explore_days].sum(axis=1).astype(np.float64)
    nbr_rate = (nbr @ train_tot) / float(ctx.n_explore_days)      # events/day/nbhd
    exp = np.maximum(nbr_rate * w, 0.1)                           # floor: avoid /0
    z = -(obs - exp[:, None]) / np.sqrt(exp)[:, None]
    return z.astype(np.float32)
