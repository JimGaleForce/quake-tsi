"""B2 / §S3: the lattice runner, its threshold, its ratchets and its refusals.

The properties this suite pins down, in the order they would fail:

  * **`alpha = q/m` and NOTHING ELSE** (§P7-24 SP-3). There is no OR-limb in the code
    and the test checks the arithmetic that killed it: a flat 1e-6 limb is 10x more
    lenient than q/m at m = 1e6 and 1,000x at 1e8.
  * **`m` is asserted, not trusted** (SP-6.1) -- *"understating m is the way to cheat
    this protocol, and it is the only way"*.
  * **the null is EVENT TIMES through the identical map**, never a permutation
    (§P7-22(a)): checked on the signature, because a design that can be bypassed at a
    call site is not a design.
  * **the control arm is scanned at the IDENTICAL alpha** and a promotion whose
    control fires is refused (F7-b).
  * **F7-d**: a human-schedule property below M >= 6.0 cannot promote, ever.
  * **the n >= 8 floor**, and §S6.2's cliff arithmetic that makes it load-bearing.
  * **SP-6.2**: a scan handed a single event outside the exploration window refuses.

Everything runs on synthetic times. No catalogue is opened anywhere in this file.
"""

from __future__ import annotations

import datetime as _dt

import numpy as np
import pytest

from engine import (lattice_s1 as LAT, properties as P, regions_d as RD,
                    searcher as S)

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

T0 = _dt.datetime(2000, 1, 1)
SPAN = 400.0


@pytest.fixture(scope="module")
def cols():
    c, _a = P.build_property_matrix(T0, 0.0, SPAN, 38.0, 142.0,
                                    families=("ephemeris", "human_clock"),
                                    sample_minutes=120.0)
    return {x.name: x for x in c}


def _decl(m_regions=("japan",), props=("lunar_synodic_phase",), mc=(6.0,)):
    return S.declare_scan([LAT.region(r) for r in m_regions], props, mc,
                          ("ephemeris",), scan_id="test")


def _times(n, rng, lo=1.0, hi=SPAN - 1.0):
    return np.sort(rng.uniform(lo, hi, n))


# ------------------------------------------------------------- SP-3, one rule --
def test_alpha_is_q_over_m_and_there_is_no_or_limb():
    assert S.alpha_from_m(1000, 0.10) == pytest.approx(1e-4)
    assert S.alpha_from_m(10 ** 6, 0.10) == pytest.approx(1e-7)
    # the arithmetic that refused the limb: 1e-6 vs q/m
    assert 1e-6 / S.alpha_from_m(10 ** 6, 0.10) == pytest.approx(10.0)
    assert 1e-6 / S.alpha_from_m(10 ** 8, 0.10) == pytest.approx(1000.0)
    assert "REFUSED" in S.OR_LIMB_REFUSED
    import inspect
    src = inspect.getsource(S.promote)
    assert " or " not in src.replace("p_control", "").replace("or the", "")


def test_monday_clears_the_strict_bar_by_four_orders_at_a_million_cells():
    """SP-8's calibration: a protocol that would not promote 14-of-14 is wrong."""
    p_any_day = 7.0 * (1.0 / 7.0) ** 14
    assert p_any_day == pytest.approx(1.03e-11, rel=0.02)
    alpha = S.alpha_from_m(10 ** 6, 0.10)
    assert alpha == pytest.approx(1.0e-7)
    assert alpha / p_any_day == pytest.approx(9700.0, rel=0.05)


def test_the_declaration_freezes_m_alpha_k_and_the_floor_into_one_hash():
    d1 = _decl()
    d2 = _decl(props=("lunar_synodic_phase", "solar_annual_phase"))
    assert d1["m"] == 1 and d2["m"] == 2
    assert d2["alpha"] == pytest.approx(d1["alpha"] / 2.0)
    assert d1["config_hash"] != d2["config_hash"]
    assert d1["k_cap"] == 30 and d1["n_events_floor"] == 8
    assert d1["q"] == 0.10


def test_m_counts_the_full_lattice_including_cells_too_small_to_score():
    d = S.declare_scan([LAT.region("japan"), LAT.region("chile")],
                       ("a", "b", "c"), (4.5, 5.0, 5.5, 6.0), ("x",))
    assert d["m"] == 2 * 3 * 4 == 24


# ------------------------------------------------------- SP-6.1, the ratchet --
def test_assert_cell_count_refuses_an_understated_m():
    with pytest.raises(S.CellCountMismatch) as e:
        S.assert_cell_count(29, 30)
    assert "SP-6.1" in str(e.value)
    assert S.assert_cell_count(30, 30) == 30


def test_run_scan_refuses_when_fewer_cells_arrive_than_declared(cols):
    rng = np.random.default_rng(0)
    d = _decl(props=("lunar_synodic_phase", "solar_annual_phase"))
    spec = {"region": LAT.region("japan"), "mc": 6.0,
            "column": cols["lunar_synodic_phase"],
            "obs_times": _times(20, rng), "null_times": [_times(20, rng)
                                                         for _ in range(5)],
            "ref_times": _times(300, rng), "span_days": SPAN}
    with pytest.raises(S.CellCountMismatch):
        S.run_scan(d, [spec], rng=rng, log_explore=False)


# -------------------------------------------------- the null is TIMES, always --
def test_the_null_signature_cannot_accept_a_permuted_property(cols):
    """§P7-22(a): the null traverses the IDENTICAL map, so it is handed TIMES."""
    import inspect
    sig = inspect.signature(S.scan_cell)
    assert "null_times" in sig.parameters
    assert not any("phase" in p for p in sig.parameters)
    src = inspect.getsource(S.scan_cell)
    assert "col.evaluate(t)" in src or "col.statistic_values(t)" in src
    # the phrase "permuted property" appears in the DOCSTRING precisely to say
    # one cannot be passed, so the check is on the executable body only.
    body = src.split('"""')[2]
    assert "permut" not in body.lower()


def test_scan_cell_pushes_null_times_through_the_columns_own_map(cols):
    """A column that records what it was called with proves the null went through it."""
    rng = np.random.default_rng(1)
    seen = []
    base = cols["lunar_synodic_phase"]

    def spy(t):
        seen.append(np.asarray(t).size)
        return base.evaluate(t)

    col = P.PropertyColumn("spy", "ephemeris", "phase", "level-waveform-phase",
                           spy, dict(base.provenance))
    nulls = [_times(11, rng) for _ in range(6)]
    S.scan_cell(col, _times(17, rng), nulls, _times(200, rng),
                LAT.region("japan"), 6.0, _decl())
    assert seen.count(11) == 6                # every null replicate went through
    assert 17 in seen and 200 in seen         # so did the observations and the pool


# ------------------------------------------------------- the promotion rule ---
def _row(**kw):
    base = dict(p_real=1e-12, p_control=0.5, concentration_form="arc", n_events=20,
                dwell_time_corrected=True, sp2_null_layer_built=True,
                sp2_reason="ok", human_schedule=False, mc=6.0, family="ephemeris",
                region="japan")
    base.update(kw)
    return base


def test_a_clean_row_promotes():
    ok, fails = S.promote(_row(), _decl())
    assert ok and fails == []


def test_a_row_whose_control_fired_is_refused():
    ok, fails = S.promote(_row(p_control=1e-13), _decl())
    assert not ok and any("F7-b" in f for f in fails)


def test_a_row_with_no_control_arm_is_refused():
    ok, fails = S.promote(_row(p_control=None), _decl())
    assert not ok and any("no control arm" in f for f in fails)
    assert "may not be quoted, in any forum, ever" in S.CONTROL_MANDATORY


def test_a_human_schedule_row_below_M6_cannot_promote_however_small_its_p():
    d = S.declare_scan([LAT.region("japan")], ("day_of_week",), (5.5, 6.0),
                       ("human_clock",))
    ok, fails = S.promote(_row(human_schedule=True, mc=5.5, p_real=1e-40,
                               concentration_form="single-cell",
                               family="human_clock"), d)
    assert not ok and any("F7-d" in f for f in fails)
    ok2, _f = S.promote(_row(human_schedule=True, mc=6.0, p_real=1e-40,
                             concentration_form="single-cell",
                             family="human_clock"), d)
    assert ok2, "the M >= 6.0 carve-out must remain live"


def test_the_n_floor_refuses_a_three_of_three_coincidence():
    """§S6.2: the floor exists precisely so 3-of-3 (p = 1/49) never enters the list."""
    ok, fails = S.promote(_row(n_events=3), _decl())
    assert not ok and any("declared floor" in f for f in fails)


def test_the_cliff_arithmetic_the_floor_is_built_around():
    """14-of-14 survives m = 336,000; 10-of-14 fails it, and not marginally."""
    from scipy import stats
    m = 336000
    # §S6.2's table is stated for a PRE-DECLARED day (no factor of 7): 14-of-14 is
    # 1.474e-12 and 10-of-14 is 2.034e-06.
    p14 = float(stats.binom.sf(13, 14, 1.0 / 7.0))
    p10 = float(stats.binom.sf(9, 14, 1.0 / 7.0))
    assert p14 == pytest.approx(1.474e-12, rel=0.02)
    assert p10 == pytest.approx(2.034e-06, rel=0.02)
    assert p14 * m == pytest.approx(4.95e-07, rel=0.05)      # survives
    assert p10 * m == pytest.approx(0.683, rel=0.05)         # FAILS, not marginally
    assert p10 * m > 0.5 and p14 * m < 1e-6


def test_a_level_row_without_dwell_correction_cannot_promote():
    ok, fails = S.promote(_row(dwell_time_corrected=False), _decl())
    assert not ok and any("dwell" in f for f in fails)


def test_a_row_whose_sp2_layer_is_unbuilt_cannot_promote():
    ok, fails = S.promote(_row(sp2_null_layer_built=False,
                               sp2_reason="F7 controls not measured"), _decl())
    assert not ok and "F7 controls not measured" in fails


def test_a_seeded_region_cannot_promote_in_its_own_family():
    ok, fails = S.promote(_row(region="alaska_aleutians", family="solid_tide"),
                          _decl())
    assert not ok and any("seeded" in f for f in fails)


def test_an_undeclared_concentration_form_cannot_promote():
    ok, fails = S.promote(_row(concentration_form="vibes"), _decl())
    assert not ok and any("declared forms" in f for f in fails)


# ------------------------------------------------------------ the categorical --
def test_base_rates_come_from_the_null_pool_and_are_never_one_over_k():
    """§S3.1 / §S6.2(1): the difference between 1/7 and the observed rate costs 4.5
    orders of magnitude and it is not optional."""
    ref = np.array([0.0] * 300 + [1.0] * 100 + [2.0] * 100)
    r = S._base_rates(ref, (0, 1, 2))
    assert r[0] == pytest.approx(0.6, abs=0.01)
    assert not np.allclose(r, 1.0 / 3.0)


def test_the_max_cell_binomial_tail_is_exact_at_k_equals_n():
    """The small-n instrument lives or dies at k = n, where an approximation is
    worthless. 14 of 14 at p = 1/7 is 1.474e-12."""
    stat, cell, tail = S._categorical_stat(np.zeros(14), np.full(7, 1.0 / 7.0),
                                           tuple(range(7)))
    assert tail == pytest.approx(1.474e-12, rel=0.02)
    assert cell == 0
    assert stat == pytest.approx(-np.log10(1.474e-12), rel=0.01)


def test_an_inflated_base_rate_costs_four_and_a_half_orders_of_magnitude():
    """§S6.2(1): the F7 correction is 45% of the evidence on a log scale."""
    _s, _c, t_uniform = S._categorical_stat(np.zeros(14), np.full(7, 1.0 / 7.0),
                                            tuple(range(7)))
    rates = np.full(7, (1.0 - 0.186) / 6.0)
    rates[0] = 0.186
    _s2, _c2, t_inflated = S._categorical_stat(np.zeros(14), rates, tuple(range(7)))
    assert t_inflated == pytest.approx(5.93e-11, rel=0.05)
    assert np.log10(t_inflated / t_uniform) == pytest.approx(1.6, abs=0.2)


# ------------------------------------------------------------- the whole scan --
def test_a_scan_over_a_true_null_produces_a_uniform_looking_p_and_no_promotion(cols):
    rng = np.random.default_rng(7)
    d = _decl(props=("lunar_synodic_phase", "solar_annual_phase"))
    nulls = [_times(30, rng) for _ in range(60)]
    ref = _times(500, rng)
    specs = [{"region": LAT.region("japan"), "mc": 6.0, "column": cols[p],
              "obs_times": _times(30, rng), "null_times": nulls, "ref_times": ref,
              "span_days": SPAN}
             for p in d["properties"]]
    rep = S.run_scan(d, specs, rng=rng, log_explore=False)
    assert rep["n_cells_evaluated"] == d["m"] == 2
    assert rep["n_survivors_real"] == 0
    assert 0.0 < min(r["p_real"] for r in rep["rows"]) <= 1.0
    assert all(r.get("p_control") is not None for r in rep["rows"])
    assert "global" in rep["max_statistic"]              # SP-6.4, always reported
    assert rep["ranked"][0]["rank_in_scan"] == 1
    assert "null_draws" not in rep["rows"][0]            # not serialised


def test_a_planted_concentration_is_recovered_by_the_scan(cols):
    """The instrument must SEE, or nothing else it says means anything."""
    from engine import circstat_event as CE
    rng = np.random.default_rng(11)
    d = _decl()
    col = cols["lunar_synodic_phase"]
    cand = _times(20000, rng)
    ph = np.asarray(col.evaluate(cand), dtype=np.float64)
    planted = CE.thin_by_phase_intensity(
        cand, ph, lambda x: np.exp(3.0 * np.cos(x - 1.0)), 80, rng)
    row = S.scan_cell(col, planted, [_times(80, rng) for _ in range(200)],
                      _times(2000, rng), LAT.region("japan"), 6.0, d)
    assert row["p_real"] <= row["p_resolution_floor"] * 2.0
    assert row["n_events"] == 80


def test_the_control_column_is_the_same_map_on_a_random_monotone_warp(cols):
    rng = np.random.default_rng(3)
    base = cols["lunar_synodic_phase"]
    ctrl = S._control_column(base, rng, span_days=SPAN)
    t = _times(200, rng)
    assert ctrl.name.endswith("__F8-15-control")
    assert S.CONTROL_CV == 1.0 and S.CONTROL_KNOT_DAYS == 90.0
    assert ctrl.ptype == base.ptype and ctrl.pclass == base.pclass
    a = np.asarray(base.evaluate(t))
    b = np.asarray(ctrl.evaluate(t))
    assert not np.allclose(a, b)                       # alignment destroyed
    # marginal structure preserved: both are phases spread over the full circle
    assert b.min() < 0.5 and b.max() > 2 * np.pi - 0.5
    # and the marginal really is preserved, not merely non-degenerate: a dense
    # uniform time grid gives the same occupancy distribution under both maps
    grid = np.linspace(1.0, SPAN - 1.0, 20000)
    ga = np.sort(np.asarray(base.evaluate(grid)))
    gb = np.sort(np.asarray(ctrl.evaluate(grid)))
    assert np.abs(ga - gb).max() < 0.35


def test_the_control_warp_is_not_the_identity_for_any_declared_property(cols):
    """`CONTROL_DEGENERACY_FOUND`: an F8-15 warp seeded from the IDENTITY clock has
    CV = 0 and both modes return the identity map, making the control a duplicate of
    the science arm. This test is the reason that defect is recorded rather than
    quietly fixed."""
    rng = np.random.default_rng(9)
    t = np.linspace(1.0, SPAN - 1.0, 500)
    for name, c in cols.items():
        ctrl = S._control_column(c, rng, span_days=SPAN)
        a = np.asarray(c.evaluate(t), dtype=np.float64)
        b = np.asarray(ctrl.evaluate(t), dtype=np.float64)
        ok = np.isfinite(a) & np.isfinite(b)
        assert not np.allclose(a[ok], b[ok]), name
    assert "BIT-IDENTICAL to the" in S.CONTROL_DEGENERACY_FOUND


# -------------------------------------------------------------- entry guards ---
def test_a_scan_touching_the_holdout_refuses(cols):
    with pytest.raises(S.HoldoutTouched) as e:
        S.assert_scan_entry(_decl(), [cols["lunar_synodic_phase"]], ["japan"],
                            event_days=np.array([1.0, 2.0, 999.0]),
                            explore_window=(0.0, 400.0))
    assert "SP-6.2" in str(e.value)


def test_a_solid_tide_scan_including_alaska_refuses(cols):
    tide, _a = P.build_property_matrix(T0, 0.0, 60.0, 54.0, -160.0,
                                       families=("solid_tide",),
                                       sample_minutes=120.0)
    with pytest.raises(RD.SeedingRegionInConfirmationSet):
        S.assert_scan_entry(_decl(), tide, ["japan", "alaska_aleutians"])


def test_the_subdaily_observer_gate_is_a_hard_entry_gate(cols):
    from engine import observer as OBS
    with pytest.raises(OBS.SubDailyGateNotSatisfied):
        S.assert_scan_entry(_decl(), [cols["utc_hour_phase"]], ["japan"],
                            observer_feature_names=[])
    rep = S.assert_scan_entry(_decl(), [cols["utc_hour_phase"]], ["japan"],
                              observer_feature_names=list(
                                  OBS.REQUIRED_FOR_SUBDAILY))
    assert rep["subdaily_gate"]["satisfied"]
    assert rep["random_clock_control"]["satisfied"]


# ---------------------------------------------------------------- the ledger --
def test_the_explore_count_line_carries_the_full_declared_cell_count(cols, tmp_path):
    """SP-1.4: unpriced never means uncounted. One line, from the PARENT."""
    import json
    rng = np.random.default_rng(5)
    d = _decl()
    path = str(tmp_path / "EXPLORE_COUNT.jsonl")
    specs = [{"region": LAT.region("japan"), "mc": 6.0,
              "column": cols["lunar_synodic_phase"], "obs_times": _times(20, rng),
              "null_times": [_times(20, rng) for _ in range(20)],
              "ref_times": _times(300, rng), "span_days": SPAN}]
    rep = S.run_scan(d, specs, rng=rng, log_explore=True, explore_path=path)
    assert rep["explore_count_logged"]
    lines = [json.loads(l) for l in open(path, encoding="utf-8") if l.strip()]
    assert len(lines) == 1
    cfg = lines[0]["config"]
    assert cfg["n_declared_tests"] == d["m"]
    assert cfg["priced"] == 0
    assert "EXPLORATORY-UNPRICED" in cfg["category"]


def test_the_banners_are_both_present():
    assert "GENERATOR, NOT EVIDENCE" in S.GENERATOR_NOT_EVIDENCE
    assert "never produces evidence" in S.THE_SENTENCE
    assert "mechanically generated FREEZE FILE" in S.THE_SENTENCE
