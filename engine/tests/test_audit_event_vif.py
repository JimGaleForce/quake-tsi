"""D-2: the event-path VIF harness. SIM only, price 0.

§P7-22 Q1 makes D-2 mandatory and §K92-0(6)'s floor table a BRACKET until it lands.
The tests here are invariants of the INSTRUMENT, not of the Earth:

  * the estimator returns ~1 on catalogues whose true VIF is 1 by construction --
    without that, an excess measured anywhere else would be the estimator;
  * the estimator RESPONDS: a self-exciting catalogue with declared clustering returns
    a VIF above 1, so a value of 1 in Arm A is a measurement and not a dead statistic;
  * the record carries its scope flags, its price of 0, and Popper's recorded
    prediction verbatim, so the verdict is scored against the prediction as recorded.

The production numbers are produced by `python -u -m engine.audit_event_vif`; these run
small so the suite stays affordable.
"""

from __future__ import annotations

import numpy as np
import pytest

from engine import audit_event_vif as V

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


# ------------------------------------------------------------- the estimator --
def test_the_estimator_returns_one_on_its_own_reference():
    """Identity check: reference against reference must give exactly the ratio 1."""
    rng = np.random.default_rng(0)
    draws = V.reference_vstar_draws(600, 8000, 120, rng)
    r = V.vif_from_draws(draws, draws)
    assert r["vif_moment"] == pytest.approx(1.0, abs=1e-12)


def test_the_estimator_responds_to_a_known_inflation():
    """Scale V* by sqrt(3) and the moment-form VIF must read 3 -- the stated scaling.

    This is the one place the V*-is-sqrt(chi2)-like adaptation is checkable in
    isolation, and it is what licenses reading the ratio as a variance inflation.
    """
    rng = np.random.default_rng(1)
    ref = V.reference_vstar_draws(400, 4000, 200, rng)
    r = V.vif_from_draws(ref * np.sqrt(3.0), ref)
    assert r["vif_moment"] == pytest.approx(3.0, rel=1e-9)


# ------------------------------------------------------------- the simulators --
def test_arm_a_catalogues_have_true_vif_one_by_construction():
    """Conditionally independent draws from lambda: the gate_r1 argument, event-side."""
    rng = np.random.default_rng(2)
    edges, lam = V.declared_intensity(400)
    t = V.arm_a_catalog(edges, lam, 5000, rng)
    assert t.size == 5000
    assert np.all(np.diff(t) >= 0)
    assert 0.0 <= t.min() and t.max() < 400.0
    # the realised time distribution follows the declared intensity, not a uniform
    first, second = np.mean(t < 200.0), np.mean(t >= 200.0)
    assert second > first                       # the network-era ramp is present


def test_arm_b_catalogues_carry_genuine_sub_daily_clustering():
    """The dependence Arm A deliberately does not have: aftershocks minutes apart."""
    rng = np.random.default_rng(3)
    edges, lam = V.declared_intensity(3650)     # the declared record length
    a = V.arm_a_catalog(edges, lam, 3000, rng)
    b = V.arm_b_catalog(edges, lam, 3000, 0.7, rng)
    # fraction of inter-event gaps under ten minutes: ~0.6 % for a Poisson process at
    # this rate, and several times that once an Omori cascade is present
    f_a = float(np.mean(np.diff(a) < 1.0 / 144.0))
    f_b = float(np.mean(np.diff(b) < 1.0 / 144.0))
    assert f_b > 3.0 * f_a, (f_a, f_b)
    with pytest.raises(ValueError):
        V.arm_b_catalog(edges, lam, 100, 1.0, rng)      # supercritical is refused


# ------------------------------------------------------------------- the run --
def test_small_run_arm_a_is_near_one_and_arm_b_is_above_it():
    """The two arms, small. Arm A ~ 1; Arm B strictly above it at the same band."""
    rec = V.run(n_catalogs=60, n_reference_draws=200, n_reference_pool=20000,
                n_events=1500, bands=("M2", "Mf"), branching_ratios=(0.7,),
                include_waveform=False, verbose=False, seed=99)
    for band in ("M2", "Mf"):
        a = rec["arm_a_true_vif_1"][band]
        b = rec["arm_b_declared_clustering"][0.7][band]
        lo, hi = a["vif_moment_ci95"]
        assert lo < 1.0 < hi, (band, a["vif_moment_ci95"])
        assert b["vif_moment"] > a["vif_moment"], band
    # and both are far below the count path's 24.1
    assert rec["verdict"]["prediction_upheld_arm_a"]
    assert rec["verdict"]["max_vif_arm_a"] < 5.0


def test_the_record_carries_the_prediction_the_price_and_the_scope_flags():
    rec = V.run(n_catalogs=12, n_reference_draws=60, n_reference_pool=4000,
                n_events=600, bands=("M2",), branching_ratios=(0.5,),
                include_waveform=False, verbose=False, seed=4)
    assert rec["priced_tests"] == 0
    assert "NOT LOGGED" in rec["explore_count"]
    assert rec["recorded_prediction"]["count_path_value_not_transferable"] == 24.1
    assert "S-17" in rec["recorded_prediction"]["s17_note"]
    assert "DECLARED, NOT FITTED" in rec["scope_flags"]
    assert "TRUE VIF = 1" in rec["scope_flags"]


def test_the_omori_shape_comes_from_the_frozen_fit():
    """Shape frozen, branching ratio declared -- and the record says which is which."""
    import json
    import os
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(V.__file__))),
                        "engine", "out", "cache", "etas_params_k080_0p5deg.json")
    if not os.path.exists(path):
        pytest.skip("frozen ETAS params not on this disk")
    with open(path, "r", encoding="utf-8") as fh:
        p = json.load(fh)["params"]
    assert V.FROZEN_OMORI_C_DAYS == pytest.approx(p["c"], rel=1e-12)
    assert V.FROZEN_OMORI_P == pytest.approx(p["p"], rel=1e-12)
