"""§P6-2 rule 7's harness: the roster is honest and the batched kernel is the engine's.

The audit itself (`python -u -m engine.audit_gpd --jobs N`) is a minutes-long run and
is not a unit test. What IS unit-tested is everything that would make its table a
lie: that the roster really spans the three test kinds at the declared count, that
the batched Lomb-Scargle used for the 10^6-surrogate arm is numerically the engine's
own scalar kernel, and that the accept/reject bars are the ones §P6-2(7) states.
"""

import numpy as np
import pytest

from engine import audit_gpd as A
from engine import mine as M


def test_roster_spans_the_three_test_kinds_at_the_declared_count():
    cases = A.roster()
    kinds = {c["kind"] for c in cases}
    assert kinds == {"glm", "mark", "period"}
    # >= 15 representative TESTS: the roster x the declared target-p grid
    assert len(cases) * len(A.TARGET_PS) >= 15
    assert all(c["null_type"] in ("block_bootstrap", "ar1", "permutation")
               for c in cases)


def test_no_roster_case_uses_a_forbidden_null():
    """§P6-2(4) again: the calibration may not be run on the all-shifts null either."""
    from engine import gpd_tail as G
    for c in A.roster():
        assert G.assert_null_permitted(c["null_type"])


def test_the_bars_are_the_ones_the_rule_states():
    assert A.COVERAGE_BAR == 0.90
    assert A.UNDERSTATE_BAR == 3.0
    assert A.N_BRUTE == 1_000_000
    assert A.N_SURROGATES_GPD == 2000


def test_batched_lomb_scargle_equals_the_engine_scalar_kernel():
    rng = np.random.default_rng(2)
    n, nf = 250, 30
    t = np.arange(n, dtype=float)
    periods = np.exp(np.linspace(np.log(2.0), np.log(n / 3.0), nf))
    Xb = rng.standard_normal((7, n))
    batch = A._ls_power_batch(t, Xb, periods)
    scalar = np.array([M.lomb_scargle_power(t, Xb[i], periods).max()
                       for i in range(Xb.shape[0])])
    assert np.allclose(batch, scalar, rtol=1e-11, atol=0.0)


def test_each_case_draw_returns_the_requested_number_of_statistics():
    for c in A.roster():
        v = np.asarray(c["draw"](12, np.random.default_rng(1)))
        assert v.shape == (12,) and np.isfinite(v).all()


@pytest.mark.parametrize("kind", ["glm", "mark", "period"])
def test_run_case_is_wired_end_to_end_at_a_toy_brute_force(kind, monkeypatch):
    """The whole comparison path, with the 10^6 arm shrunk so it runs in a test."""
    monkeypatch.setattr(A, "N_BRUTE", 4000)
    monkeypatch.setattr(A, "PILOT_N", 3000)
    monkeypatch.setattr(A, "N_SURROGATES_GPD", 500)
    monkeypatch.setattr(A, "N_BOOT", 40)
    monkeypatch.setattr(A, "N_AD_BOOT", 40)
    case = next(c for c in A.roster() if c["kind"] == kind)
    r = A.run_case(case, 0.01, 12345)
    assert r["kind"] == kind
    assert r["p_method"] in ("MC_RESOLVED", "GPD_EXTRAPOLATED", "UNRESOLVED")
    assert 0.0 < r["p_brute_force"] <= 1.0
    assert isinstance(r["covered"], bool)
