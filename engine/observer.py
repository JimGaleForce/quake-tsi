"""F7-01/02/03 -- THE OBSERVER AS PART OF THE SYSTEM, and the gate they hold shut.

MINING_CATALOG family 7: *"the catalog is a measurement made by an instrument that
is itself inside the system, and the instrument's state is a covariate we have never
entered. Family 3's standing warning (a solar hit is an observer hit until proven
otherwise) is unresolvable without this family."*

WHY THIS MODULE IS LOAD-BEARING FOR TRANCHE B AND NOT FOR TRANCHE A
-------------------------------------------------------------------
§P7-2(b) built these in A and named their consumer: *"Their consumer is B's mark
axis, not A itself."* §P7-3(3) says why, and it is the sharpest sentence in the
whole tranche:

  > *"Kepler is right that mark tests are computed at event times and therefore
  > escape the day-binning sinc entirely... BUT THAT PROPERTY IS EXACTLY WHAT
  > FORFEITS THE PROGRAM'S ONLY STRUCTURAL DEFENCE AGAINST THE OBSERVER. §P5-1 let
  > us claim, for the first time, that our fortnightly numbers are structurally
  > immune to the S1/S2 detection-cycle systematic -- and that immunity comes from
  > the two exact zeros of day-binning and global longitude summation. The mark path
  > has neither zero. Escaping the sinc means escaping the notch that was protecting
  > us."*
  >
  > **RULED: F9-10's sub-daily arm is GATED on the F7-01/02/03 observer controls.**

`assert_subdaily_gate` is that ruling in code. It is a hard refusal, not a warning,
because the failure it prevents is silent: a sub-daily mark result computed without
observer controls looks exactly like a sub-daily mark result computed with them, and
the difference only shows up in what it is allowed to mean.

WHAT THE CONTROLS ARE
---------------------
Three MEASURED instrument states, all computed from `day_float` (never from daily
bins -- F7-01's Pit is explicit that the diurnal cycle lives in the notched band and
computing it from daily bins destroys the very thing it measures):

  F7-01 `obs_diurnal_amplitude_365d`   the 24 h detection cycle's amplitude, trailing
                                       365 d, in a fixed band just above global Mc.
                                       K-091's always-on positive control made
                                       continuous: a known, large, physically certain
                                       artifact whose amplitude over time IS a direct
                                       measurement of the network's noise floor.
  F7-02 `obs_diurnal_amplitude_b<k>`   the same in four magnitude bands. The band at
                                       which it reaches zero is the global
                                       completeness magnitude, measured by a method
                                       entirely independent of the
                                       magnitude-frequency distribution.
  F7-03 `obs_weekly_amplitude_365d`    the 7 d human-operations cycle. A second,
                                       independent artifact ruler with a period in
                                       the band the miner IS sensitive to (sinc ~ 1),
                                       unlike the diurnal one -- which is what makes
                                       it diagnostic rather than merely present.

Plus the CATALOG-COMPOSITION NULL FEATURES the sub-daily arm is gated on, which are
not instrument states but composition proxies -- the ways the catalogue's own
make-up drifts under a fixed physical Earth:

  `obs_network_era_365d`      trailing 365 d event count: the network-era proxy.
  `obs_mc_drift_365d`         trailing 365 d completeness proxy (mode of the
                              magnitude histogram, the maximum-curvature estimator).
  `obs_day_of_week_phase`     the operations clock, as a phase.
  `obs_utc_hour_phase`        the UTC-hour clock. DEGENERATE ON THE COUNT PATH BY
                              CONSTRUCTION and that is the point: it is constant
                              across a daily bin, so it can only be tested on the
                              sub-daily mark path -- which makes it the exact control
                              for the arm that escapes day-binning.

ALL OF THEM ARE CONTROL FEATURES (`control=True`, family 7), so they land in a
declared CONTROL STRATUM under §P6-3 and a survivor on one is a MEASURED OBSERVER
ARTIFACT, reported as such, never as a finding. Nothing in this module is evidence.
"""

from __future__ import annotations

import numpy as np

from . import mine as M

OBSERVER_FAMILY = 7
TRAIL_DAYS = 365
MC_BAND_WIDTH = 0.5

# F7-02's four declared bands. The catalog names them; they are not tuned here.
F7_02_BANDS = ((4.5, 5.0), (5.0, 5.5), (5.5, 6.0), (6.0, 10.0))

# The names the sub-daily gate requires. Declared as a frozen tuple so the gate
# cannot be satisfied by whatever happens to be in the feature list on the day.
F7_01_NAME = "obs_diurnal_amplitude_365d"
F7_03_NAME = "obs_weekly_amplitude_365d"
COMPOSITION_NULLS = ("obs_network_era_365d", "obs_mc_drift_365d",
                     "obs_day_of_week_phase", "obs_utc_hour_phase")
REQUIRED_FOR_SUBDAILY = (F7_01_NAME, F7_03_NAME) + COMPOSITION_NULLS

SUBDAILY_GATE_RULE = (
    "HYPOTHESIS_LEDGER.md §P7-3(3): F9-10's SUB-DAILY arm is gated on the "
    "F7-01/02/03 observer controls, because the mark path's escape from the "
    "day-binning sinc is ALSO an escape from the S1/S2 structural immunity §P5-1 "
    "let us claim. The fortnightly-and-longer mark arm may proceed without waiting; "
    "the sub-daily arm may not.")

OBSERVER_CONTROL_BANNER = (
    "F7 OBSERVER CONTROLS -- these are CONTROL features in a declared control "
    "stratum (§P6-3). A survivor here is a MEASURED OBSERVER ARTIFACT and is "
    "reported as one. Family 3's standing warning applies in reverse: a hit on one "
    "of these is the instrument, and a hit on a science feature that co-moves with "
    "one of these is the instrument until shown otherwise.")


class SubDailyGateNotSatisfied(AssertionError):
    """The sub-daily mark arm was asked to run without its observer controls."""


# ------------------------------------------------------------- trailing tools --
def _trailing_resultant(day_float, phase, sel, n_days, trail=TRAIL_DAYS):
    """Trailing-window mean resultant length of `phase` over selected events.

    Returns (R_t, n_t) for each day t, over events in [t + 1 - trail, t + 1). The
    computation is over EVENT TIMES: `phase` is built from `day_float`, so the 24 h
    cycle is measured where it lives rather than where daily bins would put it
    (F7-01's Pit, which is not optional -- a daily bin has exactly one value of the
    diurnal phase and the amplitude it measures is zero by construction).
    """
    d = np.asarray(day_float, dtype=np.float64)[sel]
    ph = np.asarray(phase, dtype=np.float64)[sel]
    o = np.argsort(d, kind="stable")
    d, ph = d[o], ph[o]
    cs_c = np.concatenate([[0.0], np.cumsum(np.cos(ph))])
    cs_s = np.concatenate([[0.0], np.cumsum(np.sin(ph))])
    edge_hi = np.arange(1, int(n_days) + 1, dtype=np.float64)
    hi = np.searchsorted(d, edge_hi, side="left")
    lo = np.searchsorted(d, edge_hi - float(trail), side="left")
    n = (hi - lo).astype(np.float64)
    C = cs_c[hi] - cs_c[lo]
    S = cs_s[hi] - cs_s[lo]
    R = np.hypot(C, S) / np.maximum(n, 1.0)
    R[n < 1] = 0.0
    return R, n


def _fill_head(v, n_valid):
    """Carry the first well-determined value backwards over the burn-in.

    A trailing-365 d feature is not defined on day 3. Rather than leave a ramp that
    every downstream statistic would read as a trend, the head is held flat at the
    first fully-determined value -- and the mining window starts after the baseline's
    own burn-in anyway, so no test ever sees these days.
    """
    v = np.asarray(v, dtype=np.float64).copy()
    k = int(min(max(n_valid, 1), v.size - 1))
    v[:k] = v[k]
    return v


def _amplitude_from_R(R):
    """Fractional sinusoidal amplitude from a mean resultant length: a = 2R.

    For a rate proportional to `1 + a cos(theta)` the mean resultant length of the
    event phases is `a / 2`. Stated rather than assumed, per S-18's fourth clause.
    """
    return 2.0 * np.asarray(R, dtype=np.float64)


def diurnal_amplitude(marks, n_days, mag_lo, mag_hi, trail=TRAIL_DAYS,
                      period_days=1.0):
    """F7-01/F7-02/F7-03's kernel: trailing amplitude of a sub-daily/weekly cycle."""
    d = np.asarray(marks["day_float"], dtype=np.float64)
    mag = np.asarray(marks["mag"], dtype=np.float64)
    sel = (mag >= float(mag_lo)) & (mag < float(mag_hi))
    ph = np.mod(2 * np.pi * d / float(period_days), 2 * np.pi)
    R, n = _trailing_resultant(d, ph, sel, n_days, trail)
    valid = int(np.argmax(n >= 30)) if (n >= 30).any() else int(trail)
    return _fill_head(_amplitude_from_R(R), valid), n


def mc_maximum_curvature(marks, n_days, trail=TRAIL_DAYS, bin_width=0.1,
                         mc_floor=4.5):
    """Trailing-window maximum-curvature Mc: the mode of the magnitude histogram.

    The standard, cheap completeness estimator. It is entered as a CONTROL, not as a
    correction: its job is to say whether a "finding" co-moves with the catalogue's
    own completeness drift, which is the F7 family's whole reason to exist.
    """
    d = np.asarray(marks["day_float"], dtype=np.float64)
    mag = np.asarray(marks["mag"], dtype=np.float64)
    o = np.argsort(d, kind="stable")
    d, mag = d[o], mag[o]
    edges = np.arange(float(mc_floor), float(mag.max()) + 2 * bin_width, bin_width)
    edge_hi = np.arange(1, int(n_days) + 1, dtype=np.float64)
    hi = np.searchsorted(d, edge_hi, side="left")
    lo = np.searchsorted(d, edge_hi - float(trail), side="left")
    out = np.full(int(n_days), float(mc_floor))
    # cumulative histogram over events -> trailing histogram by difference
    binned = np.clip(np.searchsorted(edges, mag, side="right") - 1, 0,
                     edges.size - 2)
    cum = np.zeros((mag.size + 1, edges.size - 1), dtype=np.int32)
    np.add.at(cum, (np.arange(1, mag.size + 1), binned), 1)
    cum = np.cumsum(cum, axis=0)
    h = cum[hi] - cum[lo]
    tot = h.sum(axis=1)
    ok = tot >= 30
    out[ok] = edges[np.argmax(h[ok], axis=1)] + bin_width / 2.0
    return _fill_head(out, int(np.argmax(ok)) if ok.any() else int(trail))


# ----------------------------------------------------------------- the features -
def observer_features(marks, n_days, mc=4.5, lags=(0,), trail=TRAIL_DAYS):
    """F7-01/02/03 plus the catalog-composition nulls, as declared CONTROL features.

    Every one carries `control=True` and `family=7`, so `mine_session.run` keeps them
    off the science axes by name and `strata.stratum_of` routes them to their own
    declared stratum. Lags default to (0,) -- these are controls, and a control that
    quietly grew a 31-lag grid would be running 31x the declaration it was priced at.
    """
    lags = tuple(int(l) for l in lags)
    n_days = int(n_days)
    day = np.arange(n_days, dtype=np.float64)
    feats = []

    a_diurnal, _n = diurnal_amplitude(marks, n_days, mc, mc + MC_BAND_WIDTH,
                                      trail, period_days=1.0)
    feats.append(M.Feature(
        F7_01_NAME, OBSERVER_FAMILY, "linear", a_diurnal,
        "F7-01: trailing %d d amplitude of the 24 h detection cycle, M in "
        "[%.2f, %.2f), computed from day_float and NEVER from daily bins"
        % (trail, mc, mc + MC_BAND_WIDTH),
        lags=lags, control=True, periodic_override=False,
        block_days_override=float(trail)))

    for k, (lo, hi) in enumerate(F7_02_BANDS):
        a_b, n_b = diurnal_amplitude(marks, n_days, lo, hi, trail, period_days=1.0)
        feats.append(M.Feature(
            "obs_diurnal_amplitude_b%d" % k, OBSERVER_FAMILY, "linear", a_b,
            "F7-02: F7-01 in magnitude band [%.1f, %.1f). The band where this "
            "reaches zero IS the global completeness magnitude, measured "
            "independently of the magnitude-frequency distribution (%d events in "
            "the last trailing window)" % (lo, hi, int(n_b[-1])),
            lags=lags, control=True, periodic_override=False,
            block_days_override=float(trail)))

    a_week, _nw = diurnal_amplitude(marks, n_days, mc, 10.0, trail, period_days=7.0)
    feats.append(M.Feature(
        F7_03_NAME, OBSERVER_FAMILY, "linear", a_week,
        "F7-03: trailing %d d amplitude of the 7 d human-operations cycle. Unlike "
        "the diurnal one this sits at sinc ~ 1, in the band the miner IS sensitive "
        "to -- which is what makes it diagnostic" % trail,
        lags=lags, control=True, periodic_override=False,
        block_days_override=float(trail)))

    # ---- catalog-composition nulls ----------------------------------------
    d = np.asarray(marks["day_float"], dtype=np.float64)
    ds = np.sort(d)
    edge_hi = np.arange(1, n_days + 1, dtype=np.float64)
    rate = (np.searchsorted(ds, edge_hi, side="left")
            - np.searchsorted(ds, edge_hi - float(trail), side="left")
            ).astype(np.float64)
    feats.append(M.Feature(
        "obs_network_era_365d", OBSERVER_FAMILY, "linear",
        _fill_head(rate, int(trail)),
        "network-era proxy: trailing %d d event count. Detection capability grows "
        "monotonically with network build-out, so this is the composition drift a "
        "long-period 'finding' is most likely to be" % trail,
        lags=lags, control=True, periodic_override=False,
        block_days_override=float(trail)))

    feats.append(M.Feature(
        "obs_mc_drift_365d", OBSERVER_FAMILY, "linear",
        mc_maximum_curvature(marks, n_days, trail, mc_floor=mc),
        "Mc-drift proxy: trailing %d d maximum-curvature completeness magnitude. "
        "Entered as a CONTROL and never as a correction" % trail,
        lags=lags, control=True, periodic_override=False,
        block_days_override=float(trail)))

    feats.append(M.Feature(
        "obs_day_of_week_phase", OBSERVER_FAMILY, "phase",
        np.mod(2 * np.pi * day / 7.0, 2 * np.pi),
        "observer clock: day-of-week phase. Analysts review on weekdays, so this "
        "can enter the catalogue at the REVIEW stage rather than the detection "
        "stage -- a different magnitude dependence, and that difference is "
        "diagnostic (F7-03 Pit)",
        lags=lags, control=True, period_hint=7.0, periodic_override=True))

    # Constant across a daily bin BY CONSTRUCTION: on the count path its design is
    # degenerate, which is exactly what makes it the right control for the arm that
    # escapes day-binning. `mine.Feature.design` would raise on the zero-variance
    # column, so the count path must never be handed this feature -- `subdaily_only`
    # marks it and `assert_subdaily_gate` is what reads the mark.
    utc = M.Feature(
        "obs_utc_hour_phase", OBSERVER_FAMILY, "phase",
        np.zeros(n_days, dtype=np.float64),
        "observer clock: UTC hour-of-day phase. ZERO ON THE COUNT PATH BY "
        "CONSTRUCTION (a daily bin has exactly one hour phase) and live only on the "
        "sub-daily mark path -- which is precisely the arm §P7-3(3) gates",
        lags=(0,), control=True, period_hint=1.0, periodic_override=True)
    utc.subdaily_only = True
    feats.append(utc)
    return feats


def count_path_features(feats):
    """The subset a COUNT-path session may score: drops the sub-daily-only controls."""
    return [f for f in feats if not getattr(f, "subdaily_only", False)]


# ------------------------------------------------------------------- the gate --
def subdaily_gate_report(feature_names):
    """Which of the required observer controls are present, and the verdict."""
    have = set(str(n) for n in feature_names)
    missing = [n for n in REQUIRED_FOR_SUBDAILY if n not in have]
    return {
        "rule": SUBDAILY_GATE_RULE,
        "required": list(REQUIRED_FOR_SUBDAILY),
        "present": [n for n in REQUIRED_FOR_SUBDAILY if n in have],
        "missing": missing,
        "satisfied": not missing,
        "banner": OBSERVER_CONTROL_BANNER,
        "what_it_does_not_license": (
            "the gate being SATISFIED licenses the sub-daily arm to RUN. It does "
            "not license any sub-daily result to be read as Earth: a sub-daily mark "
            "result that co-moves with any of these controls is the observer, and "
            "the controls exist to make that visible, not to subtract it."),
    }


def assert_subdaily_gate(feature_names):
    """§P7-3(3), enforced. Refuses the sub-daily mark arm without its controls."""
    rep = subdaily_gate_report(feature_names)
    if not rep["satisfied"]:
        raise SubDailyGateNotSatisfied(
            "F9-10's SUB-DAILY mark arm requires the F7-01/02/03 observer controls "
            "and the catalog-composition nulls in the declared feature set; missing "
            + ", ".join(rep["missing"]) + ". " + SUBDAILY_GATE_RULE +
            " The fortnightly-and-longer mark arm is NOT gated and may proceed.")
    return rep
