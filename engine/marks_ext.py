"""F9-10 -- THE MARK AXIS EXTENDED, and the immunity it forfeits by existing.

MINING_CATALOG F9-10: *"extend the existing mark tests (currently 2 marks) to:
magnitude, depth, log-moment, distance to nearest prior event, time since prior
event, cluster membership, region... A forcing that does not change the rate may
still change WHICH events happen, and mark tests see that while count tests cannot.
S15: mark tests use the full event set (46,585): MEASURABLE, and they escape the
day-binning sinc entirely because they are computed at event times. That last
property is underexploited and is worth a tranche on its own."*  Price: 7 marks x
23+ features = 161+ tests.

THE HAZARD KEPLER DID NOT PRICE, AND IT IS THE REASON HALF THIS MODULE EXISTS
-----------------------------------------------------------------------------
§P7-3(3): *"escaping the sinc means escaping the notch that was protecting us."* The
day-binning zero and the global-longitude-summation zero are what let §P5-1 claim
STRUCTURAL immunity to the S1/S2 detection-cycle systematic. The mark path has
neither zero. So:

  * the SUB-DAILY arm is gated on `observer.assert_subdaily_gate` (§P7-3(3));
  * **G-M1 arm (i) must be re-run ON THE MARK PATH SPECIFICALLY**, and there it is a
    different test than it was on the count path. On the count path v1 is *designed*
    blind, so arm (i) demanded `A_hat / A_ref < 0.1` and *"could only destroy, never
    license"*. On the mark path the pipeline is NOT blind, so arm (i) becomes a LIVE
    FALSIFICATION: the mark path is required to RECOVER a planted local-solar-hour
    modulation, and a mark path that cannot recover a known artifact of known
    amplitude has no business reporting a sub-daily mark result. The harness is
    `recovery_b.gm1_arm_i_mark_path`; this module supplies its feature.

WHAT "AT EVENT TIMES" ACTUALLY REQUIRES
---------------------------------------
Evaluating a DAY-BINNED feature at an event does not escape day-binning: every event
on day t still gets the day-t value, and the sinc is exactly where it was. The escape
requires the feature re-derived at the event's own `day_float`, which is what
`event_time_feature_values` does through `ephemeris.ephemeris_table_at` and the
shared `mine.ephemeris_feature_specs` definitions. Features with no closed form in t
(family 4, catalogue-derived trailing statistics) cannot escape and are returned
DAY-BINNED with `subdaily=False` on the row -- stated per feature, so no row can be
read as sub-daily that is not.

ONE PRICING FINDING THIS MODULE MAKES AND DOES NOT ACT ON
---------------------------------------------------------
`log_moment = 1.5 * mag + 9.1` is a strictly increasing function of `mag`, and BOTH
declared mark statistics are rank statistics (`mine.mark_test` ranks the mark before
doing anything else: Spearman on the ranks, circular-linear on the ranks). A strictly
monotone transform does not move ranks. So the `log_moment` mark test is
BIT-IDENTICAL to the `mag` mark test -- 23 of the declared 161 tests are provably the
same test run twice. `redundancy_audit` PROVES that numerically rather than asserting
it (§P5-5(2)'s declare-then-prove discipline, applied to the transform axis exactly
as §P7-5(5) extends it there).

**This module does not resolve it.** Dropping the 23 would lower the declared count,
which is an adjudication, not a build decision; keeping them is the conservative
direction (a larger denominator can only make BH stricter). The audit is reported
with the declaration and the Popper seat decides. Nothing here is evidence.
"""

from __future__ import annotations

import numpy as np

from . import ephemeris as eph, mine as M

# The seven declared marks, in declaration order. `kind` is how the mark enters the
# rank statistic; every one of them is reduced to ranks by `mine.mark_test`.
MARK_NAMES = ("mag", "depth", "log_moment", "dist_nearest_prior_km",
              "dt_prior_days", "cluster_member", "lon_sector")
N_MARKS = len(MARK_NAMES)

# Declared constants. Frozen here, before any run, per S-9: one declared value each,
# no alternatives tried.
NEAREST_PRIOR_WINDOW_DAYS = 30.0
NEAREST_PRIOR_CAP_KM = 20015.0          # antipodal; the value for "no prior in window"
CLUSTER_RADIUS_KM = 100.0
CLUSTER_WINDOW_DAYS = 10.0
LON_SECTORS = 12
EARTH_R_KM = 6371.0
MOMENT_A, MOMENT_B = 1.5, 9.1           # log10 M0 = 1.5 M + 9.1 (Hanks & Kanamori)

MARK_DEFINITIONS = {
    "mag": "event magnitude (the v1 mark, unchanged)",
    "depth": "hypocentral depth, km (the v1 mark, unchanged)",
    "log_moment": ("log10 seismic moment = %.1f M + %.1f. PROVABLY RANK-IDENTICAL "
                   "to `mag` under both declared mark statistics -- see "
                   "`redundancy_audit`" % (MOMENT_A, MOMENT_B)),
    "dist_nearest_prior_km": (
        "great-circle distance to the nearest event in the preceding %.0f d; "
        "events with no prior in that window take the declared cap %.0f km"
        % (NEAREST_PRIOR_WINDOW_DAYS, NEAREST_PRIOR_CAP_KM)),
    "dt_prior_days": "time since the immediately preceding event, days (day_float)",
    "cluster_member": ("1 if some prior event lies within %.0f km and %.0f d, else 0"
                       % (CLUSTER_RADIUS_KM, CLUSTER_WINDOW_DAYS)),
    "lon_sector": ("longitude sector index 0..%d, ORDINAL IN LONGITUDE. A rank "
                   "statistic on it tests a monotone longitude gradient in the "
                   "mark and is NOT an omnibus over regions -- stated because the "
                   "catalog says 'region' and a reader would assume otherwise"
                   % (LON_SECTORS - 1)),
}

SUBDAILY_NOTE = (
    "§P7-3(3): a mark row is SUB-DAILY only if its feature was re-derived at the "
    "event's own day_float. A day-binned feature evaluated at an event is still "
    "day-binned and still notched; `subdaily` on the row says which it is, per row, "
    "so no result can be read as sub-daily that is not.")


# --------------------------------------------------------------- the marks ----
def _haversine_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp = p2 - p1
    dl = np.radians(lon2) - np.radians(lon1)
    a = np.sin(dp / 2.0) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2.0) ** 2
    return 2.0 * EARTH_R_KM * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))


def _nearest_prior(day_float, lat, lon, window=NEAREST_PRIOR_WINDOW_DAYS,
                   cap=NEAREST_PRIOR_CAP_KM, cluster_km=CLUSTER_RADIUS_KM,
                   cluster_days=CLUSTER_WINDOW_DAYS):
    """(distance to nearest prior within `window`, cluster flag) per event.

    STRICTLY CAUSAL: only events with a strictly smaller `day_float` are candidates,
    so a mark can never be built from the future of its own event. The `window`
    truncation is a DECLARED approximation, not an accident -- an untruncated
    nearest-prior over 46,585 events is O(N^2), and a 30 d window at a global rate of
    ~6 events/day already carries ~180 candidates per event.
    """
    d = np.asarray(day_float, dtype=np.float64)
    o = np.argsort(d, kind="stable")
    ds, las, los = d[o], np.asarray(lat)[o], np.asarray(lon)[o]
    n = ds.size
    lo = np.searchsorted(ds, ds - float(window), side="left")
    dist = np.full(n, float(cap))
    clus = np.zeros(n, dtype=np.float64)
    for i in range(n):
        a, b = int(lo[i]), i
        if b <= a:
            continue
        dk = _haversine_km(las[i], los[i], las[a:b], los[a:b])
        j = int(np.argmin(dk))
        dist[i] = float(dk[j])
        near = dk <= float(cluster_km)
        if near.any() and (ds[i] - ds[a:b][near]).min() <= float(cluster_days):
            clus[i] = 1.0
    out_d = np.empty(n)
    out_c = np.empty(n)
    out_d[o] = dist
    out_c[o] = clus
    return out_d, out_c


def build_marks(marks, n_sectors=LON_SECTORS):
    """The 7 declared F9-10 marks from the event table. Returns (dict, audit).

    `marks` is `mine.load_event_marks`'s output. `lat`/`lon` are required for the
    three spatial marks; without them the three are omitted and the audit says so
    rather than substituting something -- a mark axis quietly running 4 marks under a
    declaration that says 7 is the failure this returns an audit to prevent.
    """
    d = np.asarray(marks["day_float"], dtype=np.float64)
    mag = np.asarray(marks["mag"], dtype=np.float64)
    out = {"mag": mag,
           "depth": np.asarray(marks["depth"], dtype=np.float64),
           "log_moment": MOMENT_A * mag + MOMENT_B}

    o = np.argsort(d, kind="stable")
    dt = np.empty(d.size)
    ds = d[o]
    dt_sorted = np.concatenate([[float(NEAREST_PRIOR_WINDOW_DAYS)], np.diff(ds)])
    # The FIRST event has no prior, and its value must be a DECLARED CONSTANT rather
    # than anything computed from the record: a median-of-the-rest would make event
    # 0's mark depend on events that came after it, which is precisely the leak
    # `_nearest_prior`'s strict causality is built to prevent -- and it would silently
    # change every mark when the record is extended.
    dt[o] = dt_sorted
    out["dt_prior_days"] = dt

    have_space = ("lat" in marks and "lon" in marks)
    if have_space:
        lat = np.asarray(marks["lat"], dtype=np.float64)
        lon = np.asarray(marks["lon"], dtype=np.float64)
        dist, clus = _nearest_prior(d, lat, lon)
        out["dist_nearest_prior_km"] = dist
        out["cluster_member"] = clus
        out["lon_sector"] = np.floor(np.mod(lon, 360.0)
                                     / (360.0 / float(n_sectors))).astype(np.float64)
    audit = {
        "declared_marks": list(MARK_NAMES),
        "built_marks": [m for m in MARK_NAMES if m in out],
        "omitted_marks": [m for m in MARK_NAMES if m not in out],
        "omission_reason": (None if have_space else
                            "event lat/lon absent from the mark table: the three "
                            "spatial marks cannot be built and are OMITTED rather "
                            "than substituted"),
        "definitions": dict(MARK_DEFINITIONS),
        "n_events": int(d.size),
    }
    return out, audit


def redundancy_audit(mark_values):
    """§P5-5(2) declare-then-prove, on the TRANSFORM axis (§P7-5(5)).

    Proves numerically which declared marks are rank-identical to which others, i.e.
    which of the declared tests are provably the same test run twice under the two
    rank statistics the mark axis actually uses. Returns the pairs and the count.
    """
    names = [n for n in MARK_NAMES if n in mark_values]
    ranks = {n: M._ranks(np.asarray(mark_values[n], dtype=np.float64))
             for n in names}
    pairs = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if np.array_equal(ranks[a], ranks[b]):
                pairs.append({
                    "a": a, "b": b,
                    "verdict": "RANK-IDENTICAL -> the two mark tests are the same "
                               "test; one of them is provably FREE under §P5-5(2)",
                    "proof": "np.array_equal(rank(%s), rank(%s)) is True over all "
                             "%d events" % (a, b, ranks[a].size)})
    return {
        "n_marks_examined": len(names),
        "n_rank_identical_pairs": len(pairs),
        "pairs": pairs,
        "expected": ("log_moment = %.1f*mag + %.1f is strictly increasing, and both "
                     "declared mark statistics rank the mark first, so this pair is "
                     "predicted before it is measured -- which is what makes the "
                     "measurement a check rather than a discovery"
                     % (MOMENT_A, MOMENT_B)),
        "disposition": ("REPORTED, NOT ACTED ON. Dropping a declared test lowers the "
                        "BH denominator, which is an adjudication and not a build "
                        "decision. Keeping it is the conservative direction."),
    }


# ------------------------------------------- event-time (sub-daily) features ---
def event_time_feature_values(t0, day_float, day_index=None, day_binned=None):
    """Family-1/2 feature values RE-DERIVED AT EVENT TIMES. The escape, made real.

    Returns (values_by_name, subdaily_by_name). Names come from the SAME
    `mine.ephemeris_feature_specs` the daily features are built from, so a sub-daily
    feature and its daily namesake cannot silently diverge.

    `day_binned` (optional: {name: daily array}) supplies the fallback for features
    with no closed form in t -- they are returned at their DAY value with
    `subdaily=False`, because a day-binned feature evaluated at an event has not
    escaped anything.
    """
    e = eph.ephemeris_table_at(t0, np.asarray(day_float, dtype=np.float64))
    vals, sub = {}, {}
    for spec in M.ephemeris_feature_specs(e):
        vals[spec["name"]] = np.asarray(spec["values"], dtype=np.float64)
        sub[spec["name"]] = True
    if day_binned:
        di = (np.asarray(day_index, dtype=np.int64) if day_index is not None
              else np.floor(np.asarray(day_float)).astype(np.int64))
        for name, series in day_binned.items():
            if name in vals:
                continue
            vals[name] = np.asarray(series, dtype=np.float64)[di]
            sub[name] = False
    return vals, sub


def local_solar_hour_phase(day_float, lon):
    """G-M1 arm (i)'s feature: hour of LOCAL solar day, as a phase in [0, 2pi).

    Local solar time leads UTC by lon/15 hours, i.e. by lon/360 of a day. The
    LOCAL-ness is the whole point: §P5-1's second exact zero is the global longitude
    summation, and a UTC-hour feature summed over the globe cancels the S1/S2
    artifact that a LOCAL-hour feature does not. Arm (i) on the mark path is
    therefore run on the local phase, which is the one that can actually fire.
    """
    return np.mod(2 * np.pi * (np.asarray(day_float, dtype=np.float64)
                               + np.asarray(lon, dtype=np.float64) / 360.0),
                  2 * np.pi)


def utc_hour_phase(day_float):
    """The UTC-hour companion to `local_solar_hour_phase`, for the F7 control."""
    return np.mod(2 * np.pi * np.asarray(day_float, dtype=np.float64), 2 * np.pi)


def harmonic_amplitude(theta, values):
    """Least-squares first-harmonic amplitude of `values` against phase `theta`.

    The recovery quantity for G-M1 arm (i) on the mark path: fit
    `values ~ c + a cos(theta) + b sin(theta)` and return (sqrt(a^2 + b^2), phase).
    Reported in the MARK's own units, which for the S1/S2 artifact is magnitude --
    the observable a detection-threshold modulation actually moves.
    """
    th = np.asarray(theta, dtype=np.float64)
    v = np.asarray(values, dtype=np.float64)
    X = np.column_stack([np.ones(th.size), np.cos(th), np.sin(th)])
    beta, *_ = np.linalg.lstsq(X, v, rcond=None)
    return float(np.hypot(beta[1], beta[2])), float(np.arctan2(beta[2], beta[1]))
