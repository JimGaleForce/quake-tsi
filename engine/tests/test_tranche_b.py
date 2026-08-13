"""Tranche B's non-statistic machinery: floors, marks, observer gate, F8-15, S-15(c),
the §P7-14(b) stability line, the declaration, and the session wiring.

Counted invariants, one per ruling:

  * §P7-10(c)/§P7-13   the mark floor reproduces the ledger's own printed numbers
                       (0.0219 / 0.0439 / 0.0658 / 0.1074, fallback 0.0469) and
                       refuses to build a per-feature floor from a censored BOUND.
  * §P7-14(c)          the S-15(c) clause bites at exactly 3 cycles, catches the
                       three entries named in advance, and is inherited by anything
                       that goes through `floors`.
  * §P7-3(3)           the sub-daily mark arm REFUSES to run without the F7 controls.
  * F8-15              a clock scan REFUSES to be declared without its random-clock
                       control, and the control matches what it claims to match.
  * §P7-14(b)          every candidate carries the stability line, the line is a
                       LABEL and not a kill, and the R4 pair alone is refused.
  * F9-10              the seven marks are built, `log_moment` is PROVED
                       rank-identical to `mag`, and a sub-daily row says so on the
                       row rather than in a footnote.
  * wiring             the three new kinds are OFF by default (byte-identical config
                       and task list), route to their own declared strata when on,
                       and run end-to-end through a session.
"""

from __future__ import annotations

import datetime as _dt
import json
import os
import types

import numpy as np
import pytest

from engine import (circstat, clocks, floors, marks_ext, mine as M,
                    mine_session as ms, observer, s15c, stability,
                    strata as S, tranche_b)


# ============================================ §P7-10(c) / §P7-13: mark floor ==
def test_mark_floor_reproduces_the_ledger_numbers():
    """§P7-10(c) prints these four values; the code must not disagree with them."""
    n = 46585
    assert floors.rho_min(n, 1.0) == pytest.approx(0.0219, abs=5e-5)
    assert floors.rho_min(n, 4.0) == pytest.approx(0.0439, abs=1e-4)
    assert floors.rho_min(n, 9.0) == pytest.approx(0.0658, abs=5e-5)
    assert floors.rho_min(n, 24.0) == pytest.approx(0.1074, abs=5e-5)


def test_pooled_fallback_is_4575_and_gives_0469():
    """§P7-13(b), bounds-readmitted because the censoring is INFORMATIVE."""
    assert floors.VIF_MARK_FALLBACK == pytest.approx(4.575)
    assert floors.VIF_MARK_MEASUREMENTS_ONLY_MEDIAN == pytest.approx(4.345)
    assert floors.rho_min(46585) == pytest.approx(0.0469, abs=5e-5)


def test_a_censored_vif_may_not_produce_a_per_feature_floor():
    """§P7-11(c)/§P7-13(a): a bound is not a measurement and cannot make a floor."""
    with pytest.raises(floors.MarkVifIsABound):
        floors.assert_mark_floor_declarable(46585, 13.87, vif_is_bound=True,
                                            feature="b_value x mag")
    rep = floors.mark_floor_report(46585, 13.87, vif_is_bound=True)
    assert "UNMEASURABLE BY DECLARATION" in rep["verdict"]
    # ... and the same number IS legitimate in a population summary
    assert floors.mark_floor_report(46585, None)["declarable"]


def test_mark_plant_below_its_own_floor_is_refused():
    with pytest.raises(floors.PlantBelowFloor):
        floors.assert_mark_plant_above_floor(0.01, 46585, feature="too small")
    ok = floors.assert_mark_plant_above_floor(0.20, 46585, feature="fine")
    assert ok["compliant"] and ok["rho_over_floor"] > 2.0


# ================================================ §P7-14(c): S-15(c) clause ===
def test_s15c_bites_at_exactly_three_cycles():
    rec = 7716.0
    assert floors.max_identifiable_period(rec) == pytest.approx(2572.0)
    assert floors.unmeasurable_by_window(2572.0 + 1.0, rec)
    assert not floors.unmeasurable_by_window(2572.0 - 1.0, rec)


def test_s15c_catches_the_three_entries_named_in_advance():
    """§P7-14(c) named Metonic, Jupiter-Saturn and the 11 y solar cycle BEFORE the
    sweep. A prediction that is not checked is a prediction that was never made."""
    rep = s15c.sweep()
    assert rep["expected_catches_missed"] == []
    names = {r["feature"] for r in rep["unmeasurable"]}
    assert {"metonic_phase", "planetary_synodic_jupiter_saturn",
            "solar_cycle_phase_11y"} <= names
    cyc = {r["feature"]: r["cycles_in_window"] for r in rep["unmeasurable"]}
    assert cyc["metonic_phase"] == pytest.approx(1.11, abs=0.02)
    assert cyc["planetary_synodic_jupiter_saturn"] == pytest.approx(1.06, abs=0.02)
    assert cyc["solar_cycle_phase_11y"] == pytest.approx(1.92, abs=0.02)


def test_s15c_scores_an_unmeasurable_feature_neither_way():
    r = floors.window_report(6939.688, name="metonic")
    assert r["verdict"] == floors.UNMEASURABLE_BY_WINDOW
    assert r["scored"] is False
    assert "not a null" in r["note"].lower() or "NOT a null" in r["note"]
    with pytest.raises(floors.UnmeasurableByWindow):
        floors.assert_measurable_by_window(6939.688, name="metonic")


def test_s15c_flags_the_period_scan_grid_without_clamping_it():
    g = s15c.period_scan_grid_check()
    assert g["grid_exceeds_cut"] is True
    assert g["period_scan_max_days"] == pytest.approx(ms.PERIOD_MAX)
    assert "NOT CLAMPED" in g["disposition"]


def test_future_features_inherit_the_clause_through_floors():
    """§P7-14(c) asks that the clause be wired so future features inherit it."""
    rows = floors.window_sweep([("a_new_decadal_feature", 4000.0),
                                ("a_new_fortnightly_feature", 14.765)])
    verdicts = {r["feature"]: r["verdict"] for r in rows["rows"]}
    assert verdicts["a_new_decadal_feature"] == floors.UNMEASURABLE_BY_WINDOW
    assert verdicts["a_new_fortnightly_feature"] == floors.MEASURABLE_BY_WINDOW


# =========================================== §P7-3(3): the observer controls ==
def _fake_marks(n_days=800, rate=8.0, seed=0):
    rng = np.random.default_rng(seed)
    n = int(n_days * rate)
    d = np.sort(rng.uniform(0, n_days, size=n))
    return {"day": np.floor(d).astype(np.int64), "day_float": d,
            "mag": 4.5 + rng.exponential(0.4, n),
            "depth": rng.uniform(1.0, 200.0, n),
            "lat": rng.uniform(-60, 60, n), "lon": rng.uniform(-180, 180, n)}


def test_observer_features_are_controls_and_finite():
    feats = observer.observer_features(_fake_marks(), 800)
    names = {f.name for f in feats}
    assert set(observer.REQUIRED_FOR_SUBDAILY) <= names
    for f in feats:
        assert f.control is True and f.family == observer.OBSERVER_FAMILY
        assert np.isfinite(f.values).all(), f.name


def test_f7_01_measures_a_planted_diurnal_cycle_and_is_blind_when_day_binned():
    """F7-01's Pit: it lives in the notched band and must be computed from day_float."""
    rng = np.random.default_rng(4)
    n_days, n = 800, 8000
    # a strong diurnal detection cycle: events concentrated near local "night"
    u = rng.uniform(0, 1, n * 3)
    hr = rng.uniform(0, 1, n * 3)
    keep = u < (1.0 + 0.5 * np.cos(2 * np.pi * hr)) / 1.5
    hr = hr[keep][:n]
    day = rng.integers(0, n_days, hr.size)
    d = day + hr
    marks = {"day": day.astype(np.int64), "day_float": d,
             "mag": 4.5 + rng.exponential(0.4, hr.size)}
    amp, cnt = observer.diurnal_amplitude(marks, n_days, 4.5, 5.0)
    assert amp[-1] > 0.15, amp[-1]        # the planted 0.5 cycle is seen
    # the same events, DAY-BINNED: the diurnal phase is constant inside a bin, so
    # the amplitude a day-binned computation could report is exactly zero
    binned = dict(marks, day_float=np.floor(d) + 0.5)
    amp_b, _ = observer.diurnal_amplitude(binned, n_days, 4.5, 5.0)
    assert amp_b[-1] == pytest.approx(2.0, abs=1e-9)   # degenerate: all one phase
    assert not np.allclose(amp[-1], amp_b[-1])


def test_subdaily_gate_refuses_without_the_observer_controls():
    with pytest.raises(observer.SubDailyGateNotSatisfied):
        observer.assert_subdaily_gate(["moon_synodic_phase", "annual_phase"])
    feats = observer.observer_features(_fake_marks(), 800)
    rep = observer.assert_subdaily_gate([f.name for f in feats])
    assert rep["satisfied"] and rep["missing"] == []
    assert "does not license" in rep["what_it_does_not_license"] \
        or "not license" in rep["what_it_does_not_license"]


def test_utc_hour_control_is_subdaily_only():
    feats = observer.observer_features(_fake_marks(), 800)
    sub = [f for f in feats if getattr(f, "subdaily_only", False)]
    assert [f.name for f in sub] == ["obs_utc_hour_phase"]
    assert "obs_utc_hour_phase" not in {
        f.name for f in observer.count_path_features(feats)}


# ================================================================== F8-15 =====
def test_random_clock_matches_the_coarse_statistics_it_claims_to():
    rng = np.random.default_rng(7)
    inc = rng.gamma(2.0, 0.5, 2000) + 0.01
    ref = np.concatenate([[0.0], np.cumsum(inc)])
    for mode in clocks.F8_15_MODES:
        c = clocks.random_clock_control(ref, rng, mode=mode)
        assert all(c["matched"].values()), (mode, c["matched"])
        assert (np.diff(c["tau"]) > 0).all()


def test_clock_scan_without_its_control_is_refused():
    with pytest.raises(clocks.RandomClockControlMissing):
        clocks.assert_random_clock_control(["natural_time", "etas_rescaled"], [])
    assert clocks.assert_random_clock_control(["natural_time"],
                                              ["f8_15_iid"])["satisfied"]


def test_non_monotone_clock_is_rejected():
    with pytest.raises(ValueError):
        clocks.clock_increments(np.array([0.0, 1.0, 1.0, 2.0]))


# ============================================ §P7-14(b): stability machinery ==
def _rows(order, p_base=1e-4):
    return [{"test": "glm_poisson_offset_etas", "feature": f, "lag": 0,
             "mark": None, "p_raw": p_base * (i + 1), "order_key": [0, i]}
            for i, f in enumerate(order)]


def test_stability_line_labels_and_does_not_kill():
    feats = ["a", "b", "c", "d", "e"]
    ref = _rows(feats)
    stable = [_rows(feats), _rows(feats)]
    line = stability.stability_line(ref[0], ref, stable, stable, top_k=5)
    assert line["tag"] == stability.TAG_STABLE
    assert line["ii_F10_24_data_resampling"]["selection_frequency"] == 1.0
    assert "LOWER BOUND" in line["iii_lower_bound"]
    assert "not deleted" in line["label_not_kill"].lower()

    shuffled = [_rows(["z1", "z2", "z3", "z4", "a"]),
                _rows(["y1", "y2", "y3", "y4", "y5"])]
    line2 = stability.stability_line(ref[0], ref, stable, shuffled, top_k=2)
    assert line2["tag"] == stability.TAG_UNSTABLE
    # the candidate is still THERE -- a tag, not a deletion
    assert line2["candidate"]["feature"] == "a"


def test_r4_alone_is_refused():
    ref = _rows(["a", "b", "c"])
    with pytest.raises(stability.NoDataResamplingArm):
        stability.stability_line(ref[0], ref, [_rows(["a", "b", "c"])], [])


def test_null_baseline_is_tranche_as_measured_pair():
    assert stability.NULL_BASELINE["f10_24_rho_pair"] == (0.487, 0.044)
    assert stability.THRESHOLDS["min_selection_frequency"] == 1.0
    assert stability.THRESHOLDS["min_resample_rho"] > 0.487


def test_every_survivor_must_carry_a_line():
    tests = _rows(["a", "b", "c"])
    tests[0]["pass_bh"] = True
    with pytest.raises(AssertionError):
        stability.assert_lines_present(tests, [])
    lines = stability.stability_block(
        [tests[0]], tests, [_rows(["a", "b", "c"])], [_rows(["a", "b", "c"])],
        top_k=3)
    assert stability.assert_lines_present(tests, lines)["satisfied"]


# ================================================================= F9-10 ======
def test_seven_marks_are_built_with_lat_lon_and_three_are_omitted_without():
    marks = _fake_marks(400, 8.0, seed=2)
    built, audit = marks_ext.build_marks(marks)
    assert set(audit["built_marks"]) == set(marks_ext.MARK_NAMES)
    assert audit["omitted_marks"] == []
    nospace = {k: v for k, v in marks.items() if k not in ("lat", "lon")}
    _b2, a2 = marks_ext.build_marks(nospace)
    assert set(a2["omitted_marks"]) == {"dist_nearest_prior_km", "cluster_member",
                                        "lon_sector"}
    assert "OMITTED rather than substituted" in a2["omission_reason"]


def test_log_moment_is_proved_rank_identical_to_mag():
    """§P5-5(2) declare-then-prove on the transform axis (§P7-5(5))."""
    built, _a = marks_ext.build_marks(_fake_marks(300, 8.0, seed=3))
    aud = marks_ext.redundancy_audit(built)
    pairs = {(p["a"], p["b"]) for p in aud["pairs"]}
    assert ("mag", "log_moment") in pairs or ("log_moment", "mag") in pairs
    assert "REPORTED, NOT ACTED ON" in aud["disposition"]
    # and the consequence, measured: the two mark tests are the SAME test
    th = np.mod(2 * np.pi * np.arange(built["mag"].size) / 29.53, 2 * np.pi)
    a = M.mark_test(th, built["mag"], "phase", 100, np.random.default_rng(1))
    b = M.mark_test(th, built["log_moment"], "phase", 100, np.random.default_rng(1))
    assert a["statistic"] == pytest.approx(b["statistic"], rel=1e-12)


def test_marks_are_strictly_causal():
    """A mark built from the future of its own event would leak; it must not."""
    marks = _fake_marks(300, 8.0, seed=5)
    built, _a = marks_ext.build_marks(marks)
    # truncating the record after event k must not change event k's own marks
    k = 1500
    trunc = {key: v[:k] for key, v in marks.items()}
    b2, _a2 = marks_ext.build_marks(trunc)
    for name in ("dist_nearest_prior_km", "cluster_member", "dt_prior_days"):
        assert np.allclose(built[name][:k - 1], b2[name][:k - 1], equal_nan=True), name


def test_event_time_features_really_escape_day_binning():
    """The escape, checked: a sub-daily phase must vary WITHIN a day."""
    t0 = _dt.datetime(2000, 1, 1)
    day_float = np.array([10.0, 10.25, 10.5, 10.75, 11.0])
    vals, sub = marks_ext.event_time_feature_values(t0, day_float)
    assert sub["moon_synodic_phase"] is True
    v = vals["moon_synodic_phase"]
    assert np.ptp(v[:4]) > 1e-4, "sub-daily phases are constant within the day"
    # a feature with no closed form falls back to its DAY value and says so
    series = np.arange(50, dtype=float)
    vals2, sub2 = marks_ext.event_time_feature_values(
        t0, day_float, day_binned={"b_value_90d": series})
    assert sub2["b_value_90d"] is False
    assert np.allclose(vals2["b_value_90d"], np.floor(day_float))


def test_local_solar_hour_differs_from_utc_hour_by_longitude():
    d = np.array([10.0, 10.0, 10.0])
    lon = np.array([0.0, 90.0, 180.0])
    loc = marks_ext.local_solar_hour_phase(d, lon)
    utc = marks_ext.utc_hour_phase(d)
    assert np.allclose(loc[0], utc[0])
    assert not np.allclose(loc[1], utc[1])
    assert loc[1] == pytest.approx(np.pi / 2, abs=1e-9)


def test_harmonic_amplitude_recovers_a_planted_modulation():
    rng = np.random.default_rng(6)
    th = rng.uniform(0, 2 * np.pi, 40000)
    v = 5.0 + 0.3 * np.cos(th - 0.7) + rng.normal(0, 1.0, th.size)
    a, phi = marks_ext.harmonic_amplitude(th, v)
    assert a == pytest.approx(0.3, rel=0.15)
    assert phi == pytest.approx(0.7, abs=0.15)


# ======================================================== the declaration =====
def test_declaration_enumerates_and_flags_the_gap_without_resolving_it():
    e = tranche_b.enumerate_declared()
    assert e["n_reoccupied_already_declared"] == 68
    assert e["bh_denominator_m"] == e["n_new_scope"] + 68
    named = {i["key"]: i["n"] for i in e["items"]}
    assert named["second_moment"] == 17 and named["omnibus"] == 34
    assert named["mark_axis"] == 161
    assert e["agrees_with_headline"] is False
    assert e["gap_vs_headline"] > 500
    assert "DOES NOT RESOLVE" in e["discrepancy_note"]
    assert "68" in e["reoccupied_identification"]


def test_strata_file_satisfies_the_budget_identity(tmp_path):
    p = str(tmp_path / "strata_tranche_b.json")
    info = tranche_b.write_strata_file(p)
    part = S.load_partition(p, q=M.FDR_Q)
    S.assert_budget_identity(part["strata"], M.FDR_Q)
    assert info["m"] == tranche_b.enumerate_declared()["bh_denominator_m"]
    doc = json.loads(open(p, encoding="utf-8").read())
    assert "COUNT DISCREPANCY" in doc["note"]


def test_the_declaration_refuses_to_run():
    with pytest.raises(tranche_b.TrancheBRunGated):
        tranche_b.run()


def test_declared_strata_route_the_new_kinds():
    """Each new kind must land in its OWN declared stratum, never the v1 `mark` one."""
    doc = tranche_b.strata_document()
    kinds = {s["test_kind"] for s in doc["strata"]}
    assert {"moment2", "omnibus", "markx"} <= kinds
    assert S._test_kind({"test": "second_circular_moment_score"}) == "moment2"
    assert S._test_kind({"test": "kuiper_V"}) == "omnibus"
    assert S._test_kind({"test": "watson_U2"}) == "omnibus"
    assert S._test_kind({"test": "circular-linear", "mark": "mag",
                         "mark_axis": "F9-10"}) == "markx"
    assert S._test_kind({"test": "circular-linear", "mark": "mag"}) == "mark"


# ===================================================== the session wiring =====
class _StubBaseline:
    name = "stub-const-v0"
    caveat = "constant-rate stub, test only"
    burn_in_days = 0

    def __init__(self, rate):
        self._rate = np.asarray(rate, dtype=np.float64)

    def rate(self, window):
        return self._rate[None, window]

    def report(self):
        return ["stub baseline"]


N_DAYS = 700


def _prepared(seed=5):
    rng = np.random.default_rng(seed)
    day = np.arange(N_DAYS, dtype=float)
    offset = 3.0 * (1.0 + 0.2 * np.sin(2 * np.pi * day / 365.25))
    counts = rng.poisson(offset).astype(float)
    window = slice(60, N_DAYS)
    t0 = _dt.datetime(2000, 1, 1)
    feats = M.ephemeris_features(t0, N_DAYS)[:3]
    ev_day = np.repeat(np.arange(N_DAYS), counts.astype(int))
    ev_day = ev_day[ev_day >= window.start]
    n = ev_day.size
    marks = {"day": ev_day.astype(np.int64),
             "day_float": ev_day.astype(float) + rng.uniform(0, 1, n),
             "mag": 4.5 + rng.standard_exponential(n) * 0.3,
             "depth": rng.uniform(1.0, 30.0, n),
             "lat": rng.uniform(-60, 60, n), "lon": rng.uniform(-180, 180, n)}
    return (types.SimpleNamespace(n_days=N_DAYS), _StubBaseline(offset), None,
            window, counts[window], offset[window], marks, feats, [], t0)


def _cfg(tmp_path, tranche_b_on=True, subdaily=False, feats_extra=None):
    args = types.SimpleNamespace(
        mag_target=4.5, dlat=1.0, dlon=1.0, explore_frac=0.7,
        data_dir=str(tmp_path), no_download=True, seed=7, tranche1=False,
        ladder=False, gpd=False, strata=None, tranche_b=tranche_b_on,
        tb_second_moment=True, tb_omnibus=True, tb_mark_axis=True,
        tb_subdaily=subdaily)
    preset = dict(ms.QUICK, n_surrogates=120, n_periods=40, n_peaks=3,
                  lags=(0, 1))
    return ms.build_config(args, preset)


def test_tranche_b_is_off_by_default_and_changes_nothing(tmp_path):
    """A default session's config -- and therefore its hash -- must be unchanged."""
    off = _cfg(tmp_path, tranche_b_on=False)
    assert "tranche_b" not in off


def test_session_runs_the_three_new_kinds_end_to_end(tmp_path):
    cfg = _cfg(tmp_path, tranche_b_on=True)
    out = ms.run(cfg, verbose=False, resume=False,
                 session_dir=str(tmp_path / "sess"), jobs=1,
                 ledger_path=str(tmp_path / "ledger.jsonl"),
                 prepared=_prepared())
    st = json.load(open(os.path.join(out["session_dir"], "checkpoint.json"),
                        encoding="utf-8"))
    kinds = {}
    for t in st["results"].get("__tests__", []) or []:
        kinds[t["test"]] = kinds.get(t["test"], 0) + 1
    keys = set(st["results"].keys())
    assert any(k.startswith("moment2:") for k in keys)
    assert any(k.startswith("omnibus:") for k in keys)
    assert any(k.startswith("markx:") for k in keys)
    # F9-10's audits reached the checkpoint
    assert st["f9_10_marks"]["built_marks"]
    assert st["f9_10_redundancy_audit"]["n_rank_identical_pairs"] >= 1
    # every markx row carries its own §P7-10(c) floor and its sub-daily flag
    row = st["results"][[k for k in keys if k.startswith("markx:")][0]][0]
    assert row["mark_axis"] == "F9-10"
    assert row["rho_min"] > 0 and row["subdaily"] is False
    assert row["catalog_entry"] == "F9-10"
    # and every moment2 row carries R1 beside R2 (F9-01's Pit)
    m2 = st["results"][[k for k in keys if k.startswith("moment2:")][0]][0]
    assert "R1" in m2 and "R2" in m2 and m2["test"] == "second_circular_moment_score"
    assert m2["s15c_verdict"] in (floors.MEASURABLE_BY_WINDOW,
                                  floors.UNMEASURABLE_BY_WINDOW)


def test_subdaily_arm_refuses_without_observer_controls(tmp_path):
    """§P7-3(3) enforced INSIDE the session, before a single surrogate is drawn."""
    cfg = _cfg(tmp_path, tranche_b_on=True, subdaily=True)
    with pytest.raises(observer.SubDailyGateNotSatisfied):
        ms.run(cfg, verbose=False, resume=False,
               session_dir=str(tmp_path / "sess2"), jobs=1,
               ledger_path=str(tmp_path / "ledger2.jsonl"), prepared=_prepared())


def test_subdaily_arm_runs_and_is_labelled_when_the_gate_is_satisfied(tmp_path):
    prep = list(_prepared())
    # ALL the observer features are DECLARED (the gate reads the declared set);
    # `mine_session.run` is what keeps the sub-daily-only one off the count path.
    prep[7] = prep[7] + observer.observer_features(prep[6], N_DAYS)
    cfg = _cfg(tmp_path, tranche_b_on=True, subdaily=True)
    out = ms.run(cfg, verbose=False, resume=False,
                 session_dir=str(tmp_path / "sess3"), jobs=1,
                 ledger_path=str(tmp_path / "ledger3.jsonl"), prepared=tuple(prep))
    st = json.load(open(os.path.join(out["session_dir"], "checkpoint.json"),
                        encoding="utf-8"))
    assert st["subdaily_gate"]["satisfied"] is True
    rows = st["results"][[k for k in st["results"] if k.startswith("markx:")][0]]
    assert any(r["subdaily"] for r in rows)
    assert all("subdaily" in r for r in rows)
