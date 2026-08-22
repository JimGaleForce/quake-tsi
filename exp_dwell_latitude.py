"""THE DWELL ARTIFACT AS A FUNCTION OF LATITUDE. A property of the instrument, not of
any catalogue. Priced 0. No events, no region, no window, no holdout.

WHY THIS IS WORTH MEASURING.

Everything killed or wounded in §P7-25 was killed by DWELL TIME: the frozen quadrant
occupies 38.3% of the clock rather than 25%, so a null catalogue looked like a strong
effect. `engine/dwell_null.py` then measured that the same hazard exists in the
directional coordinate, and found the bearing's axial concentration R2 running from
0.14 at Alaska latitudes to 0.59 at the equator. That second number raised the question
this module answers: HOW MUCH of the dwell artifact is a function of LATITUDE alone?

The answer matters for two separate reasons.

  1. AS A HAZARD MAP. Any level, quadrant or bearing statistic compared ACROSS regions
     at different latitudes is comparing quantities with different nulls. A study that
     pools a Chilean segment with an Aleutian one, or that compares them, inherits a
     latitude-dependent artifact that has nothing to do with earthquakes. This module
     says how big that is before anyone pools anything.

  2. AS A DISCRIMINATOR, which is the more interesting use. A pure waveform artifact
     MUST track latitude, because it is a property of the tidal geometry at the site.
     A real triggering effect need not, and generally should not: it should track
     tectonic setting, forcing amplitude and criticality. So the latitude dependence of
     a claimed effect is a cheap, powerful, and to this program's knowledge un-exploited
     test of whether the effect is instrumental. If a "signal" has the same latitude
     profile as the null's dwell structure, it is the dwell structure.

WHAT IS COMPUTED. At a ladder of latitudes, on a deterministic astronomical field over
a fixed span, sampled uniformly in time: the three-band level occupancy under both
normalisations, the below-neutral-and-falling quadrant duty cycle, and the axial
concentration of the principal-extension bearing. Longitude is averaged over so the
result is a function of latitude alone rather than of one meridian's accidents.

Nothing here is evidence about earthquakes. It is a characterisation of the ruler.
"""

from __future__ import annotations

import datetime as _dt
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from engine import dwell_null as DN            # noqa: E402
from engine import ephemeris as E              # noqa: E402
from engine import tidal_tensor as TT          # noqa: E402

OUT_JSON = HERE / "results_dwell_latitude.json"

T0 = _dt.datetime(2020, 1, 1)
SPAN_DAYS = 60.0
STEP_MINUTES = 2.0
RATE_HALF_MINUTES = 10.0
LATITUDES = list(range(0, 76, 5))
LONGITUDES = [0.0, 60.0, 120.0, 180.0, 240.0, 300.0]   # averaged over

BANDS = (("TROUGH", lambda u: u < -0.5),
         ("MID", lambda u: (u >= -0.5) & (u <= 0.5)),
         ("CREST", lambda u: u > 0.5))


def local_cycle_u(t, s):
    """u within each bracketing-maxima cycle. Same construction as D-1b."""
    import exp_k092_d1 as D1
    tm = D1.refined_maxima(t, s)
    if tm.size < 2:
        return np.full(s.shape, np.nan)
    k = np.searchsorted(tm, t, side="right") - 1
    ok = (k >= 0) & (k < tm.size - 1)
    edges = np.searchsorted(t, tm)
    lo = np.full(tm.size - 1, np.nan)
    hi = np.full(tm.size - 1, np.nan)
    for j in range(tm.size - 1):
        a, b = edges[j], edges[j + 1]
        if b > a:
            lo[j], hi[j] = s[a:b].min(), s[a:b].max()
    kk = np.clip(k, 0, tm.size - 2)
    u = 2.0 * (s - lo[kk]) / np.maximum(hi[kk] - lo[kk], 1e-300) - 1.0
    return np.where(ok, u, np.nan)


def one_site(lat, lon, t, jd, rate_k):
    st = TT.strain_tensor(jd, lat, lon, 0.0)
    s = st["areal_strain"]
    rate = np.full(s.shape, np.nan)
    rate[rate_k:-rate_k] = ((s[2 * rate_k:] - s[:-2 * rate_k])
                            / (2.0 * RATE_HALF_MINUTES / 60.0))
    m = np.isfinite(rate)
    quad = float(np.mean((s[m] < 0.0) & (rate[m] < 0.0)))
    g = DN.occupancy(2.0 * (s - s.min()) / (s.max() - s.min()) - 1.0, BANDS)
    lc = DN.occupancy(local_cycle_u(t, s), BANDS)
    bearing = TT.principal_bearing(st["e_NN"], st["e_EE"], st["e_NE"])
    ax = DN.axial_occupancy(bearing)
    return {"quadrant_duty": quad, "global": g, "local_cycle": lc,
            "axial_R2": ax["axial_R2"],
            "axial_chi2_per_df": ax["chi2_vs_uniform"] / ax["chi2_df"]}


def main():
    step_days = STEP_MINUTES / 1440.0
    t = step_days * np.arange(int(SPAN_DAYS / step_days), dtype=np.float64)
    jd = E.julian_day_at(T0, t)
    rate_k = int(round(RATE_HALF_MINUTES / STEP_MINUTES))

    rows = []
    for lat in LATITUDES:
        per_lon = [one_site(float(lat), lon, t, jd, rate_k) for lon in LONGITUDES]
        agg = {
            "lat": lat,
            "quadrant_duty": float(np.mean([p["quadrant_duty"] for p in per_lon])),
            "quadrant_duty_sd": float(np.std([p["quadrant_duty"] for p in per_lon])),
            "axial_R2": float(np.mean([p["axial_R2"] for p in per_lon])),
            "axial_chi2_per_df": float(np.mean([p["axial_chi2_per_df"]
                                                for p in per_lon])),
        }
        for key in ("global", "local_cycle"):
            for b in ("TROUGH", "MID", "CREST"):
                agg["%s_%s" % (key, b)] = float(np.mean([p[key][b] for p in per_lon]))
        rows.append(agg)
        print("  lat %2d  quad %.4f  R2 %.4f  localcyc T/M/C %.3f/%.3f/%.3f"
              % (lat, agg["quadrant_duty"], agg["axial_R2"],
                 agg["local_cycle_TROUGH"], agg["local_cycle_MID"],
                 agg["local_cycle_CREST"]), flush=True)

    lats = np.array([r["lat"] for r in rows], dtype=float)
    r2 = np.array([r["axial_R2"] for r in rows])
    qd = np.array([r["quadrant_duty"] for r in rows])
    out = {
        "arm": "dwell artifact vs latitude (instrument characterisation)",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 0,
        "not_evidence": ("A deterministic astronomical field on a uniform time grid. "
                         "No catalogue, no events, no region, no window, no holdout. "
                         "This characterises the RULER, not any measurement made with "
                         "it."),
        "construction": {"t0": T0.isoformat(), "span_days": SPAN_DAYS,
                         "step_minutes": STEP_MINUTES,
                         "rate_half_minutes": RATE_HALF_MINUTES,
                         "latitudes": LATITUDES, "longitudes_averaged": LONGITUDES,
                         "scalar": "sitetide areal_strain via tidal_tensor",
                         "scope_flags": TT.TENSOR_SCOPE_FLAGS},
        "rows": rows,
        "summary": {
            "axial_R2_range": [float(r2.min()), float(r2.max())],
            "axial_R2_at_equator": float(r2[0]),
            "axial_R2_at_lat55": float(np.interp(55.0, lats, r2)),
            "quadrant_duty_range": [float(qd.min()), float(qd.max())],
            "quadrant_duty_at_lat55": float(np.interp(55.0, lats, qd)),
            "spearman_R2_vs_lat": float(np.corrcoef(
                np.argsort(np.argsort(lats)), np.argsort(np.argsort(r2)))[0, 1]),
        },
        "the_two_uses": {
            "hazard_map": ("any level, quadrant or bearing statistic compared ACROSS "
                           "latitudes is comparing quantities with different nulls; "
                           "the size of that is tabulated here"),
            "discriminator": ("a waveform artifact MUST track latitude because it is a "
                              "property of the tidal geometry at the site; a real "
                              "triggering effect should track tectonic setting, "
                              "forcing amplitude and criticality instead. A claimed "
                              "effect whose latitude profile matches the dwell "
                              "structure IS the dwell structure. Cheap, and to this "
                              "program's knowledge un-exploited."),
        },
        "honest_limit": (
            "This is a BODY-TIDE field with no ocean loading, and ocean loading is "
            "both larger and differently distributed at coastal sites. The latitude "
            "profile of the TOTAL forcing is therefore NOT this profile, and the "
            "discriminator above is only valid where the body tide dominates or where "
            "the ocean term is modelled. Stated here because it is exactly the kind of "
            "caveat that gets dropped when a table is quoted."),
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print("\n  axial R2: %.4f at the equator -> %.4f at lat 55 ; range %.4f to %.4f"
          % (r2[0], out["summary"]["axial_R2_at_lat55"], r2.min(), r2.max()))
    print("  quadrant duty cycle: %.4f to %.4f (never 0.25 anywhere)"
          % (qd.min(), qd.max()))
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
