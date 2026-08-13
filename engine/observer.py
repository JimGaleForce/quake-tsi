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


# --------------------------------------------- D-5: the gate against a SESSION --
# §K92-1 D-5 / §P7-22: *"Clear the F9-10 sub-daily gate: F7-01/02/03 observer
# controls. `observer.assert_subdaily_gate` hard-refuses until this exists. Nothing
# sub-daily runs before it."* The controls themselves were built in Tranche A
# (§P7-2(b)); what was still owed is the EVALUATION -- pointing the gate at a real
# session's F7 rows and reading out what it says. That is what this section is.
#
# TWO STAGES, and the distinction is the whole reason the PASS path is worth writing.
#
#   STAGE 1, THE GATE ITSELF -- PRESENCE. This is the ruled condition and the only
#   thing that can PASS or REFUSE: are the F7-01/02/03 controls and the four
#   catalog-composition nulls in the session's declared, SCORED feature set? §P7-3(3)
#   gates on their existence, not on their values, and widening the gate to "and they
#   must be quiet" would be a new ruling this module is not entitled to make.
#
#   STAGE 2, THE READING -- REPORTED, NEVER GATING. What the controls actually
#   measured in that session. A control sitting at its Monte Carlo resolution floor is
#   a MEASURED OBSERVER ARTIFACT of high significance, and it is exactly what
#   OBSERVER_CONTROL_BANNER says must be visible. It does not change the verdict; it
#   changes what a sub-daily result would be allowed to MEAN, which is the clause
#   `what_it_does_not_license` already carries.

SESSION_GATE_STAGES = (
    "STAGE 1 (GATES): are the required F7 controls present in the session's scored "
    "feature set? §P7-3(3) gates on existence. STAGE 2 (REPORTED, NOT GATING): what "
    "those controls measured. A control at its Monte Carlo resolution floor is a "
    "measured observer artifact; it does not move the verdict, it moves what a "
    "sub-daily result would be permitted to mean.")

COUNT_PATH_DEGENERACY_NOTE = (
    "`obs_utc_hour_phase` is ZERO ON THE COUNT PATH BY CONSTRUCTION -- a daily bin "
    "has exactly one UTC-hour phase -- so `observer_features` marks it "
    "`subdaily_only` and `count_path_features` DROPS it before a count-path session "
    "scores anything. A count-path session therefore CANNOT satisfy this gate, and "
    "that is the design working: the one control that is specific to the arm which "
    "escapes day-binning is the one a day-binned session cannot supply. Clearing the "
    "gate requires a session that scores the SUB-DAILY MARK path, where that feature "
    "is live.")


def _session_tests(checkpoint):
    return list((checkpoint or {}).get("tests") or [])


def session_gate_reading(checkpoint):
    """STAGE 2: what the F7 controls measured in this session. Never gating.

    Reads only fields the checkpoint already records -- no surrogates are drawn and
    nothing is recomputed, exactly as `f4_58_vif.py` reads a session.
    """
    rows = []
    for r in _session_tests(checkpoint):
        name = str(r.get("feature", ""))
        if not name.startswith("obs_"):
            continue
        p_shift = r.get("p_circular_shift")
        p_boot = r.get("p_block_bootstrap")
        p_floor = r.get("p_floor")
        at_floor = []
        if p_shift is not None and p_floor is not None and p_shift <= p_floor * 1.0001:
            at_floor.append("circular_shift")
        if (p_boot is not None and r.get("n_surrogates")
                and p_boot <= 1.0 / (1.0 + float(r["n_surrogates"])) * 1.0001):
            at_floor.append("block_bootstrap")
        rows.append({
            "feature": name,
            "family": r.get("family"),
            "df": r.get("df"),
            "chi2_score": r.get("chi2_score"),
            "pct_rate_modulation": r.get("pct_rate_modulation"),
            "p_circular_shift": p_shift,
            "p_block_bootstrap": p_boot,
            "p_raw": r.get("p_raw"),
            "p_resolution_floor": p_floor,
            "at_resolution_floor": at_floor,
            "disposition": r.get("disposition"),
        })
    rows.sort(key=lambda x: x["feature"])
    n_floor = sum(1 for x in rows if x["at_resolution_floor"])
    return {
        "n_observer_rows": len(rows),
        "n_at_resolution_floor": n_floor,
        "rows": rows,
        "reading": (
            "OBSERVER STRUCTURE IS LIVE AND MEASURED: %d of %d F7 control rows have "
            "at least one null at its Monte Carlo resolution floor. Per the F7 "
            "banner these are MEASURED OBSERVER ARTIFACTS, reported as such and "
            "never as findings -- and per §P7-3(3) they are the reason the sub-daily "
            "arm is gated at all." % (n_floor, len(rows))
            if n_floor else
            "no F7 control row sits at its Monte Carlo resolution floor in this "
            "session"),
        "gating": False,
        "banner": OBSERVER_CONTROL_BANNER,
    }


def session_subdaily_gate(checkpoint, session_id=None):
    """D-5: evaluate the sub-daily gate against a real session. Both stages.

    Returns the full record. `verdict` is `"PASS"` or `"REFUSE"` and is decided by
    STAGE 1 ALONE (§P7-3(3)); STAGE 2's reading is attached and explicitly
    non-gating.
    """
    names = sorted({str(r.get("feature", ""))
                    for r in _session_tests(checkpoint) if r.get("feature")})
    presence = subdaily_gate_report(names)
    reading = session_gate_reading(checkpoint)
    subdaily_only_missing = [n for n in presence["missing"]
                             if n == "obs_utc_hour_phase"]
    return {
        "item": "D-5 sub-daily gate evaluation against a real session",
        "session": session_id,
        "config_hash": (checkpoint or {}).get("config_hash"),
        "session_kind": (checkpoint or {}).get("kind"),
        "n_scored_tests": len(_session_tests(checkpoint)),
        "stages": SESSION_GATE_STAGES,
        "stage_1_presence": presence,
        "stage_2_reading": reading,
        "verdict": "PASS" if presence["satisfied"] else "REFUSE",
        "verdict_basis": ("STAGE 1 only (§P7-3(3) gates on the controls' EXISTENCE). "
                          "STAGE 2 is reported and does not gate."),
        "count_path_degeneracy": (COUNT_PATH_DEGENERACY_NOTE
                                  if subdaily_only_missing else None),
        "rule": SUBDAILY_GATE_RULE,
        "what_a_pass_would_not_license": presence["what_it_does_not_license"],
    }


def assert_session_subdaily_gate(checkpoint, session_id=None):
    """The hard refusal, against a session. Returns the record when it PASSES."""
    rec = session_subdaily_gate(checkpoint, session_id)
    if rec["verdict"] != "PASS":
        raise SubDailyGateNotSatisfied(
            "session %s does not clear the F9-10 sub-daily gate; missing %s. %s%s"
            % (session_id, ", ".join(rec["stage_1_presence"]["missing"]),
               SUBDAILY_GATE_RULE,
               " " + COUNT_PATH_DEGENERACY_NOTE
               if rec["count_path_degeneracy"] else ""))
    return rec
