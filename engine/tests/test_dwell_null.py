"""Acceptance tests for engine/dwell_null.py, the §P7-24 SP-2 null layer.

The two that matter most are `test_sinusoid_reproduces_the_exact_arcsine_thirds` --
which anchors the machinery to the one case with a closed form, so a wrong
implementation cannot hide behind "the waveform is complicated" -- and
`test_real_tidal_bearing_null_is_not_uniform_on_the_circle`, which is the measurement
§P7-25(5) K-096(ii) demands before any directional axis may be read.
"""

from __future__ import annotations

import datetime as _dt

import numpy as np
import pytest

from engine import dwell_null as DN
from engine import ephemeris as E
from engine import tidal_tensor as TT

BANDS = (("TROUGH", lambda u: u < -0.5),
         ("MID", lambda u: (u >= -0.5) & (u <= 0.5)),
         ("CREST", lambda u: u > 0.5))

# Non-Alaska site: these tests must not be a look at anything.
LAT, LON = 35.7, 139.7
T0 = _dt.datetime(2020, 6, 1)


# ------------------------------------------------------------- the power floor --
def test_power_floor_refuses_a_family_that_cannot_fire():
    with pytest.raises(ValueError):
        DN.assert_power_floor(m=8, n_replicates=128, alpha=0.05)   # floor 0.062
    ok = DN.assert_power_floor(m=8, n_replicates=2000, alpha=0.05)
    assert ok["floor"] == pytest.approx(8 / 2001.0)
    assert ok["min_attainable_p"] == pytest.approx(1 / 2001.0)


def test_power_floor_is_the_formula_not_a_constant():
    """Adding a statistic must raise the requirement automatically."""
    a = DN.assert_power_floor(m=8, n_replicates=300)
    b = DN.assert_power_floor(m=12, n_replicates=300)
    assert b["floor"] > a["floor"]
    with pytest.raises(ValueError):
        DN.assert_power_floor(m=200, n_replicates=300)


# ---------------------------------------------------------- the closed-form anchor --
def test_sinusoid_reproduces_the_exact_arcsine_thirds():
    """A pure sinusoid sampled uniformly in TIME must give exactly 1/3 per band.

    This is the case §P7-22 Ratification 1 derived in closed form: for uniform phase
    and u = cos(theta) normalised to [-1, 1], P(u < -1/2) = 1 - arccos(-1/2)/pi = 1/3
    exactly, and the other two follow by symmetry. If the machinery cannot reproduce
    the one number that is known analytically, nothing it says about a real waveform
    is worth reading.
    """
    theta = np.linspace(0.0, 2000.0 * np.pi, 4_000_001)
    occ = DN.occupancy(np.cos(theta), BANDS)
    for name in ("TROUGH", "MID", "CREST"):
        assert occ[name] == pytest.approx(1.0 / 3.0, abs=2e-3), (name, occ)


def test_occupancy_bands_partition():
    x = np.random.default_rng(0).normal(size=50_000)
    occ = DN.occupancy(np.clip(x, -1, 1), BANDS)
    assert sum(occ.values()) == pytest.approx(1.0)


def test_occupancy_ignores_nan():
    v = np.array([-0.9, np.nan, 0.0, 0.9, np.nan])
    occ = DN.occupancy(v, BANDS)
    assert sum(occ.values()) == pytest.approx(1.0)
    assert occ["TROUGH"] == pytest.approx(1 / 3.0)


# --------------------------------------------------- the real waveform is not 1/3 --
def test_the_normalisation_choice_changes_the_null_and_neither_choice_is_thirds():
    """MEASURED, and it is the trap this module exists to stop.

    The band null depends on HOW the level is normalised, and the two natural choices
    give materially different answers. Against the GLOBAL min/max of the span the
    distribution is strongly CREST-DEPLETED (~0.09), because the global maximum is set
    by the largest spring tide in the span and the signal spends almost no time near
    it. Against each LOCAL CYCLE's own range both tails are elevated and the middle is
    depleted, which is the classic dwell-at-the-turning-points shape. NEITHER IS 1/3,
    so anyone quoting a band null owes the normalisation along with the number -- the
    "conventions travel with numbers" rule, in the one place it silently matters most.
    This test's first version asserted a crest EXCESS under global normalisation and
    failed at 0.088; the expectation was wrong, not the code.

    The exact figures are site- and span-dependent (§P7-25(1)(i)) and pinning it here would
    make this test a hostage to a span nobody declared.
    """
    t = np.arange(0.0, 60.0, 1.0 / 1440.0)
    jd = E.julian_day_at(T0, t)
    s = TT.strain_tensor(jd, LAT, LON, 0.0)["areal_strain"]
    g = DN.occupancy(2.0 * (s - s.min()) / (s.max() - s.min()) - 1.0, BANDS)
    assert g["CREST"] < 0.15                       # MEASURED ~0.088, see below
    assert g["TROUGH"] > 1.0 / 3.0

    lc = DN.occupancy(_local_cycle_u(t, s), BANDS)
    assert lc["TROUGH"] > 1.0 / 3.0
    assert lc["CREST"] > 1.5 * g["CREST"]          # local-cycle lifts the crest tail
    assert lc["MID"] < 1.0 / 3.0                   # dwell piles at the turning points

    for occ in (g, lc):
        assert abs(occ["TROUGH"] - 1.0 / 3.0) > 0.04


def _local_cycle_u(t, s):
    """u normalised within each bracketing-maxima cycle -- the D-1b construction."""
    import exp_k092_d1 as D1
    tm = D1.refined_maxima(t, s)
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


def test_local_cycle_null_reproduces_the_committed_D1b_constants():
    """Cross-check against a number already committed by a DIFFERENT code path.

    results_k092_d1_null.json measured the local-cycle band probabilities at the
    Alaska seed sites as TROUGH 0.4800 / MID 0.3091 / CREST 0.2109, via
    exp_k092_d1_null.py and the app's own scalar. This module, built independently and
    taking a different route on a different scalar, must land in the same
    neighbourhood at a comparable site. Exact agreement is NOT expected -- different
    sites, different span, different scalar, and §P7-25(1)(i) already recorded that the
    constant is span-dependent -- but a disagreement in SHAPE would mean one of the two
    is wrong, and that is worth a standing test.
    """
    t = np.arange(0.0, 60.0, 1.0 / 1440.0)
    jd = E.julian_day_at(T0, t)
    s = TT.strain_tensor(jd, 56.0, -156.0, 0.0)["areal_strain"]
    occ = DN.occupancy(_local_cycle_u(t, s), BANDS)
    assert occ["TROUGH"] == pytest.approx(0.48, abs=0.06)
    assert occ["MID"] == pytest.approx(0.31, abs=0.06)
    assert occ["CREST"] == pytest.approx(0.21, abs=0.06)


# ------------------------------------------- the directional null, K-096(ii) --
def test_real_tidal_bearing_null_is_not_uniform_on_the_circle():
    """§P7-25(5) K-096(ii), measured rather than asserted.

    The bearing of maximum horizontal extension sweeps, and it does NOT sweep at a
    constant rate, so its dwell distribution is not uniform. Reading an azimuth
    concentration against a uniform-circle null would be D-1's exact error in a new
    coordinate, and this test is the evidence that the error would be material.
    """
    t = np.arange(0.0, 60.0, 1.0 / 1440.0)
    jd = E.julian_day_at(T0, t)
    s = TT.strain_tensor(jd, LAT, LON, 0.0)
    b = TT.principal_bearing(s["e_NN"], s["e_EE"], s["e_NE"])
    rep = DN.axial_occupancy(b)
    assert rep["axial_R2"] > 0.05, rep["axial_R2"]
    assert rep["chi2_vs_uniform"] > 3.0 * rep["chi2_df"]


def test_axial_occupancy_is_flat_for_a_uniform_bearing():
    b = np.linspace(0.0, 180.0, 200_001)[:-1]
    rep = DN.axial_occupancy(b)
    assert rep["axial_R2"] < 1e-3


# ----------------------------------------------------------------- the pool --
def test_pool_is_margin_only_and_rejects_impossible_margins():
    p = DN.time_uniform_pool(1000, 100)
    assert p[0] == 100 and p[-1] == 899 and p.size == 800
    with pytest.raises(ValueError):
        DN.time_uniform_pool(100, 50)


def test_draw_replicates_stays_inside_the_pool():
    rng = np.random.default_rng(7)
    pool = DN.time_uniform_pool(500, 50)
    idx = DN.draw_replicates(pool, n_sites=9, n_replicates=200, rng=rng)
    assert idx.shape == (9, 200)
    assert idx.min() >= 50 and idx.max() <= 449


# ---------------------------------------------------------- max-statistic --
def test_calibrate_recovers_a_planted_effect_and_ignores_a_null_one():
    """Positive and negative control on the calibrator itself."""
    rng = np.random.default_rng(3)
    names = ("a", "b", "c")
    null = {k: rng.normal(size=4000) for k in names}
    quiet = DN.calibrate({k: 0.0 for k in names}, null, names)
    assert quiet["max_statistic_p"] > 0.2
    loud = DN.calibrate({"a": 6.0, "b": 0.0, "c": 0.0}, null, names)
    assert loud["max_statistic_p"] <= 2.0 / 4001.0
    assert loud["max_abs_z_observed"] > 5.0


def test_calibrate_p_is_never_zero_and_counts_ties_against():
    names = ("a",)
    null = {"a": np.zeros(100)}
    r = DN.calibrate({"a": 0.0}, null, names)
    assert r["max_statistic_p"] > 0.0


def test_correlated_statistics_do_not_inflate_significance():
    """Ten copies of the same statistic must not behave like ten independent looks."""
    rng = np.random.default_rng(11)
    base = rng.normal(size=5000)
    names = tuple("s%d" % i for i in range(10))
    null = {k: base + 1e-9 * rng.normal(size=5000) for k in names}
    r = DN.calibrate({k: 2.0 for k in names}, null, names)
    indep = {k: rng.normal(size=5000) for k in names}
    r2 = DN.calibrate({k: 2.0 for k in names}, indep, names)
    # the duplicated battery must be no HARDER to clear than the independent one
    assert r["max_statistic_p"] <= r2["max_statistic_p"] + 1e-9


def test_grid_wide_max_is_stricter_than_any_single_cell():
    rng = np.random.default_rng(5)
    cells = [np.abs(rng.normal(size=3000)) for _ in range(20)]
    obs = 3.0
    single = (int(np.sum(cells[0] >= obs)) + 1) / 3001.0
    gw = DN.combine_max_statistic(cells, obs)
    assert gw["grid_wide_max_statistic_p"] >= single
    assert gw["n_cells"] == 20


def test_span_sensitivity_reports_the_range():
    rng = np.random.default_rng(2)
    by_span = {d: np.clip(rng.normal(size=20000), -1, 1) for d in (5, 10, 30)}
    rep = DN.span_sensitivity(by_span, BANDS)
    lo, hi = rep["_range"]["TROUGH"]
    assert lo <= hi
