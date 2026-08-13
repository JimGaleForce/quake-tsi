"""D-5: the sub-daily gate evaluated against a real session's F7 rows.

§K92-1 D-5 / §P7-3(3). The controls were built in Tranche A; the owed debt was the
EVALUATION -- the PASS path wired up and pointed at a session that exists. These tests
cover both outcomes, because a gate that has only ever been observed refusing has not
been shown to be a gate:

  * a feature set carrying all the required controls PASSES;
  * the real Tranche B session REFUSES, for the reason the design says it must;
  * STAGE 2's reading is reported and DOES NOT gate.
"""

from __future__ import annotations

import os

import pytest

from engine import audit_subdaily_gate as G, observer as O


def _rows(names):
    return {"kind": "mine", "config_hash": "test", "tests": [
        {"feature": n, "family": 7, "df": 1, "chi2_score": 1.0,
         "p_circular_shift": 0.5, "p_block_bootstrap": 0.5, "p_raw": 0.5,
         "p_floor": 1e-4, "n_surrogates": 9999, "disposition": "UNPRICED-CONTROL"}
        for n in names]}


# ------------------------------------------------------------- the PASS path --
def test_the_gate_passes_when_every_required_control_is_present():
    """The path that had never been exercised: STAGE 1 satisfied -> PASS."""
    ck = _rows(list(O.REQUIRED_FOR_SUBDAILY) + ["moon_declination"])
    rec = O.session_subdaily_gate(ck, "synthetic_subdaily_session")
    assert rec["verdict"] == "PASS"
    assert rec["stage_1_presence"]["satisfied"] is True
    assert rec["stage_1_presence"]["missing"] == []
    assert rec["count_path_degeneracy"] is None
    # and it returns rather than raises
    assert O.assert_session_subdaily_gate(ck, "synthetic")["verdict"] == "PASS"


def test_a_pass_does_not_license_reading_a_subdaily_result_as_earth():
    ck = _rows(O.REQUIRED_FOR_SUBDAILY)
    rec = O.session_subdaily_gate(ck, "synthetic")
    assert "does not license" in rec["what_a_pass_would_not_license"]
    assert "is the observer" in rec["what_a_pass_would_not_license"]


# ----------------------------------------------------------- the REFUSE path --
def test_the_gate_refuses_a_count_path_feature_set_and_says_why():
    """`obs_utc_hour_phase` is degenerate on the count path and is dropped there.

    The count-path name list is built by the REAL machinery -- `observer_features`
    then `count_path_features` -- so if the drop rule ever changed, this test moves
    with it instead of asserting a frozen list.
    """
    import numpy as np
    n_days = 800
    rng = np.random.default_rng(0)
    marks = {"day_float": np.sort(rng.random(4000) * n_days),
             "mag": 4.5 + rng.exponential(0.4, 4000),
             "depth": rng.random(4000) * 50.0}
    feats = O.count_path_features(O.observer_features(marks, n_days))
    ck = _rows([f.name for f in feats])
    rec = O.session_subdaily_gate(ck, "synthetic_count_path")
    assert rec["verdict"] == "REFUSE"
    assert rec["stage_1_presence"]["missing"] == ["obs_utc_hour_phase"]
    assert "ZERO ON THE COUNT PATH BY CONSTRUCTION" in rec["count_path_degeneracy"]
    with pytest.raises(O.SubDailyGateNotSatisfied):
        O.assert_session_subdaily_gate(ck, "synthetic_count_path")


# --------------------------------------------------- STAGE 2 is not a gate ----
def test_stage_2_reading_is_reported_and_does_not_gate():
    """Controls screaming at the resolution floor must not flip a satisfied STAGE 1."""
    ck = _rows(O.REQUIRED_FOR_SUBDAILY)
    for r in ck["tests"]:
        r["p_circular_shift"] = 1e-4        # exactly the recorded floor
        r["chi2_score"] = 700.0
    rec = O.session_subdaily_gate(ck, "synthetic_loud")
    assert rec["stage_2_reading"]["gating"] is False
    assert rec["stage_2_reading"]["n_at_resolution_floor"] == len(ck["tests"])
    assert "MEASURED OBSERVER ARTIFACT" in rec["stage_2_reading"]["reading"]
    assert rec["verdict"] == "PASS"          # STAGE 1 alone decides


# ------------------------------------------------------- the REAL session -----
_REAL = os.path.join(G.MINE_DIR, G.DEFAULT_SESSION, "checkpoint.json")


@pytest.mark.skipif(not os.path.exists(_REAL),
                    reason="the Tranche B session is not on this disk")
def test_the_real_tranche_b_session_refuses_and_the_reading_is_live():
    """What the gate says TODAY, from the session on disk. Reported honestly."""
    rec = G.run()
    assert rec["verdict"] == "REFUSE"
    assert "obs_utc_hour_phase" in rec["stage_1_presence"]["missing"]
    # the other five required controls ARE there -- the debt is the sub-daily path,
    # not the controls
    assert set(rec["stage_1_presence"]["present"]) == {
        "obs_diurnal_amplitude_365d", "obs_weekly_amplitude_365d",
        "obs_network_era_365d", "obs_mc_drift_365d", "obs_day_of_week_phase"}
    # and the observer structure is measured and live
    rd = rec["stage_2_reading"]
    assert rd["n_observer_rows"] >= 9
    assert rd["n_at_resolution_floor"] >= 1
    floors_hit = {r["feature"] for r in rd["rows"] if r["at_resolution_floor"]}
    assert "obs_weekly_amplitude_365d" in floors_hit
