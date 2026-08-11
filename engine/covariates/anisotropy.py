"""anisotropy -- K-081 sniff-grade directional concentration of recent large events.

For each cell and day, take the M>=m_min events in the annulus [r0, r1] km over the
past `days` days, compute the azimuths from the cell to those events, and return the
circular resultant length R in [0, 1]: R~0 means the surrounding activity is spread
evenly around the cell, R~1 means it is all on one side.

STRUCTURE-AWARE CAVEAT (v1: printed, not solved):
  * Azimuths are measured to the SOURCE CELL CENTRE, not to the hypocentre. At 1 deg
    cells over a 100-600 km annulus that is a bearing error of order 10 deg; it is
    enough for a coarse "one-sided vs surrounded" contrast and not enough for anything
    fault-resolved.
  * R is undefined when the annulus is empty; those cell-days get R=0, which is the
    same value an isotropically-surrounded cell gets. The count N is NOT carried as a
    second covariate in v1, so "no data" and "no anisotropy" are conflated.
  * No fault geometry, no focal mechanisms, no plate-boundary orientation enters this.
    A real structure-aware version needs all three.
"""

from __future__ import annotations

import numpy as np

from . import register

DEFAULTS = {"m_min": 5.5, "r_inner_km": 100.0, "r_outer_km": 600.0, "days": 365}

CAVEAT = (
    "anisotropy caveat: azimuths taken to source CELL CENTRES (not hypocentres); "
    "empty annulus and isotropic annulus both map to R=0; no fault geometry or focal "
    "mechanisms are used. Structure-aware only in the loosest sense."
)


@register("anisotropy", defaults=DEFAULTS, burn_in=365,
          describe="circular resultant length R of azimuths to prior M>=m_min in an annulus")
def compute(ctx, p):
    p.setdefault("burn_in", int(p["days"]))
    r, _n = ctx.past_azimuth_resultant(
        p["m_min"], p["r_inner_km"], p["r_outer_km"], int(p["days"])
    )
    return np.nan_to_num(r, nan=0.0, posinf=0.0, neginf=0.0)
