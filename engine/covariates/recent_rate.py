"""recent_rate -- the sanity covariate.

log(1 + count of M>=m_min within radius_km, over the past `days` days).

This one MUST come out positive-skill under the climatology-v1 baseline. That is
simultaneously (a) the engine's own smoke test -- if this is flat, the plumbing is
broken -- and (b) the demonstration of the declared v1 gap: a time-stationary
climatology does not absorb aftershock clustering, so "there were recent nearby
earthquakes" is worth real bits against it. It is NOT a discovery. An ETAS baseline
is expected to eat most or all of it.
"""

from __future__ import annotations

import numpy as np

from . import register

DEFAULTS = {"m_min": 4.5, "radius_km": 100.0, "days": 30}


@register("recent_rate", defaults=DEFAULTS, burn_in=30,
          describe="log1p(count of M>=m_min within radius_km over past `days`)")
def compute(ctx, p):
    p.setdefault("burn_in", int(p["days"]))
    n = ctx.past_counts(p["m_min"], p["radius_km"], p["days"])
    return np.log1p(n)
