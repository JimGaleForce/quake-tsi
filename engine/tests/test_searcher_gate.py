"""SP-7: the searcher's own gate -- its arithmetic, its refusals and its non-vacuity.

The gate itself is a several-minute run and is NOT executed here; what is tested is
every piece it is assembled from, plus a MINIATURE end-to-end gate (3 catalogues, a
small lattice) that exercises the identical code path in seconds.

The properties pinned down:

  * **the true-null ETAS simulator really is a Hawkes process** -- it clusters, and its
    branching ratio comes out near the declared K. A null arm without triggering is
    the failure §P7-23(A.3) names, and this is the check that the gate is not making
    it;
  * **the observed arm and the null ensemble come from the SAME generator** -- checked
    on the source, because a gate whose two arms came from different simulators would
    be measuring the difference between two simulators;
  * **the discreteness of the threshold** -- `effective_alpha` is
    `floor(alpha (B+1))/(B+1)`, so the gate prices its expectation against the rate the
    instrument can ATTAIN;
  * **non-vacuity** -- every declared cell must be promotion-eligible apart from its
    p-values, and the planted-signal control must FIRE;
  * **the verdict is computed from the declared criteria and from nothing else.**
"""

from __future__ import annotations

import datetime as _dt
import json
import os

import numpy as np
import pytest

from engine import gate_calibration as GC, searcher as S, searcher_gate as SG

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


# ------------------------------------------------------- the true-null ETAS ----
def test_the_simulator_is_a_hawkes_process_and_actually_clusters():
    """§P7-23(A.3): FULL ETAS WITH TRIGGERING, never background-only. A Poisson draw
    would have inter-event times that are exponential; a Hawkes draw's are far more
    dispersed, because the offspring arrive in bursts."""
    rng = np.random.default_rng(0)
    t = SG.simulate_etas_times(0.0, 20000.0, rng)
    dt = np.diff(t)
    cv = float(dt.std() / dt.mean())
    pois_cv = float(np.mean([
        (lambda d: d.std() / d.mean())(
            np.diff(np.sort(rng.uniform(0.0, 20000.0, t.size))))
        for _ in range(5)]))
    # a Poisson process gives CV = 1 exactly; the Hawkes draw must be clearly over it.
    assert pois_cv == pytest.approx(1.0, abs=0.1), pois_cv
    assert cv > 1.20, cv
    assert cv > 1.20 * pois_cv, (cv, pois_cv)


def test_the_branching_ratio_comes_out_near_the_declared_K():
    """n_total / n_background -> 1/(1-K). Declared K = 0.45 -> ~1.82."""
    rng = np.random.default_rng(2)
    tot, bg, reps = 0, 0, 40
    for _ in range(reps):
        n_bg = int(rng.poisson(SG.ETAS_MU * 20000.0))
        bg += n_bg
        tot += SG.simulate_etas_times(
            0.0, 20000.0, np.random.default_rng(int(rng.integers(1 << 30)))).size
    # each replicate draws its own background, so compare expectations
    expect = SG.ETAS_MU * 20000.0 * reps / (1.0 - SG.ETAS_K)
    assert 0.75 * expect < tot < 1.25 * expect, (tot, expect)


def test_the_background_is_homogeneous_so_the_gate_plants_nothing():
    """A seasonal or diurnal background would plant real structure in
    `solar_annual_phase` / `utc_hour_phase` and the gate would measure the plant."""
    import inspect
    src = inspect.getsource(SG.simulate_etas_times)
    assert "rng.uniform(lo, hi" in src
    for forbidden in ("sin", "cos", "seasonal", "annual", "diurnal"):
        assert forbidden not in src.lower().split('"""')[-1], forbidden


def test_the_observed_arm_and_the_null_ensemble_share_one_generator():
    """A gate whose arms came from different simulators measures the simulators."""
    import inspect
    src = inspect.getsource(SG.run_gate)
    assert src.count("_draw_until") >= 3
    assert "simulate_etas_times" in inspect.getsource(SG._draw_until)


def test_draw_until_conditions_on_N_in_both_arms_by_being_one_function():
    import inspect
    assert "common-mode" in inspect.getdoc(SG._draw_until)


# ------------------------------------------------------------ the arithmetic --
def test_effective_alpha_is_the_attainable_rate_not_the_declared_one():
    assert SG.effective_alpha(0.00333333, 400) == pytest.approx(1.0 / 401.0)
    assert SG.effective_alpha(0.01, 400) == pytest.approx(4.0 / 401.0)
    assert SG.effective_alpha(1.0 / 401.0 - 1e-12, 400) == 0.0
    assert SG.effective_alpha(0.5, 9) == pytest.approx(5.0 / 10.0)


def test_the_expected_promotion_count_is_N_times_q_regardless_of_m():
    """The arithmetic behind the gate's `cap_is_the_mean_note`: at alpha = q/m the
    expectation across N catalogues is N * m * (q/m) = N * q, so SP-7's cap of 3 is
    exactly the MEAN, not an upper bound on it."""
    for m in (30, 1000, 10 ** 6):
        alpha = S.alpha_from_m(m, 0.10)
        assert SG.N_CATALOGS * m * alpha == pytest.approx(3.0)
    assert SG.MAX_TOTAL == 3
    assert SG.EXPECTED_PER_CATALOG == 0.10


def test_the_binomial_tail_is_the_exact_one():
    from scipy import stats
    for k, n, p in ((4, 900, 0.0025), (0, 900, 0.0025), (10, 100, 0.05)):
        assert SG._binom_sf(k, n, p) == pytest.approx(
            float(stats.binom.sf(k - 1, n, p)), rel=1e-9)
    assert SG._binom_sf(0, 10, 0.1) == 1.0
    assert SG._binom_sf(11, 10, 0.1) == 0.0


# --------------------------------------------------------------- the fixture --
def test_the_gate_lattice_is_declared_and_its_m_is_the_full_product():
    d = SG.build_gate_declaration()
    assert d["m"] == len(SG.GATE_REGIONS) * len(SG.GATE_PROPERTIES) * \
        len(SG.GATE_MAG_STRATA)
    assert d["alpha"] == pytest.approx(0.10 / d["m"])
    assert SG.GATE_MAG_STRATA == (6.0,), "human-schedule cells are live only at M>=6"


def test_the_gate_columns_include_a_categorical_and_a_human_schedule_property():
    """Otherwise whole statistic paths and whole SP-2 classes go ungated."""
    per = SG.build_gate_columns(_dt.datetime(2000, 1, 1), span_days=200.0)
    cols = per[SG.GATE_REGIONS[0]]["columns"]
    kinds = {c.ptype for c in cols}
    classes = {c.pclass for c in cols}
    assert kinds == {"phase", "level", "categorical"}
    assert classes == {"level-waveform-phase", "human-schedule"}


def test_every_gate_cell_is_promotion_eligible_apart_from_its_p_values():
    """Non-vacuity: a gate where nothing could ever promote proves nothing."""
    per = SG.build_gate_columns(_dt.datetime(2000, 1, 1), span_days=200.0)
    d = SG.build_gate_declaration()
    for r in SG.GATE_REGIONS:
        for c in per[r]["columns"]:
            ok, why = c.may_promote()
            assert ok, (r, c.name, why)
            row = dict(p_real=1e-30, p_control=0.5,
                       concentration_form=S.concentration_form(c),
                       n_events=40, dwell_time_corrected=True,
                       sp2_null_layer_built=True, sp2_reason=why,
                       human_schedule=(c.pclass == "human-schedule"),
                       mc=6.0, family=c.family, region=r)
            good, fails = S.promote(row, d)
            assert good, (r, c.name, fails)


def test_the_observer_report_travels_beside_every_human_schedule_cell():
    """F7-a. In the gate the catalogues are simulated, so the report says exactly
    that rather than implying a measurement was made."""
    per = SG.build_gate_columns(_dt.datetime(2000, 1, 1), span_days=200.0)
    obs = per[SG.GATE_REGIONS[0]]["observer"]
    assert obs["measured"] is False
    assert "SIMULATED" in obs["why"] and "REAL scan" in obs["why"]


# ----------------------------------------------------------- end to end, small --
def test_a_miniature_gate_runs_the_identical_code_path_and_writes_its_artifact(
        tmp_path, monkeypatch):
    monkeypatch.setattr(SG, "SPAN_DAYS", 400.0)
    monkeypatch.setattr(SG, "GATE_REGIONS", ("japan", "chile"))
    monkeypatch.setattr(SG, "GATE_PROPERTIES",
                        ("lunar_synodic_phase", "day_of_week"))
    monkeypatch.setattr(SG, "ETAS_MU", 0.06)
    out = str(tmp_path / "gate.json")
    art = SG.run_gate(n_catalogs=3, seed=5, out_path=out, verbose=False,
                      n_null=60, n_ref=300)
    assert os.path.exists(out)
    on_disk = json.load(open(out, encoding="utf-8"))
    assert on_disk["artifact_hash"] == art["artifact_hash"]
    assert art["m"] == 4 and art["n_catalogs"] == 3
    assert art["n_trials"] == 12
    assert art["verdict"] in ("PASS", "FAIL", "VACUOUS-FAIL")
    assert art["n_eligible_min"] == art["m"], "the miniature gate must not be vacuous"
    assert art["vacuity_control"]["promoted"] is True
    assert len(art["per_catalog"]) == 3
    assert "what_a_pass_licenses" in art and "right for a real scan" in \
        art["what_a_pass_licenses"]
    assert art["alpha_effective"] <= art["alpha_declared"]


def test_the_verdict_is_a_pure_function_of_the_declared_criteria():
    """Read off the source: the verdict must not depend on anything but the three
    declared criteria, or a future edit could tune it after seeing a result."""
    import inspect
    src = inspect.getsource(SG.run_gate)
    i = src.index("passed = bool(")
    expr = src[i:i + 120]
    assert 'v2["passed"]' in expr, expr
    assert "seed" not in expr and "MAX_TOTAL" not in expr
    # the band itself must contain no written-down threshold
    gc_src = inspect.getsource(GC.calibration_band)
    assert "= 3" not in gc_src and "MAX" not in gc_src


def test_the_cap_diagnostic_is_labelled_a_diagnostic_and_not_a_criterion(
        tmp_path, monkeypatch):
    monkeypatch.setattr(SG, "SPAN_DAYS", 400.0)
    monkeypatch.setattr(SG, "GATE_REGIONS", ("japan",))
    monkeypatch.setattr(SG, "GATE_PROPERTIES", ("lunar_synodic_phase",))
    monkeypatch.setattr(SG, "ETAS_MU", 0.06)
    art = SG.run_gate(n_catalogs=2, seed=1, out_path=str(tmp_path / "g.json"),
                      verbose=False, n_null=60, n_ref=300)
    assert "A DIAGNOSTIC, NOT A CRITERION" in art["cap_is_the_mean_note"]
    assert "DOES NOT RELAX THE CAP" in art["cap_is_the_mean_note"]
    assert art["n_cells_clearing_alpha_before_control"] >= art["n_promotions_total"]


def test_the_gate_reports_a_control_arm_count_beside_the_real_one():
    """§S3.3(1): only the DIFFERENCE is interpretable, so both are always present."""
    import inspect
    src = inspect.getsource(SG.run_gate)
    assert '"n_control_promotions_total"' in src
    assert '"n_promotions_total"' in src


# -------------------------------------------------------------- the standing ---
def test_the_gate_constants_are_sp7s_own_numbers():
    assert SG.N_CATALOGS >= 30
    assert SG.MAX_TOTAL == 3 and SG.EXPECTED_PER_CATALOG == 0.10
    assert "binding" in SG.__doc__.lower()
    assert "No real scan runs before this passes" in SG.__doc__


def test_the_cli_refuses_a_real_scan_without_a_passing_artifact(tmp_path, capsys):
    """§P7-24 SP-7 is binding and no flag bypasses it."""
    from engine import cli
    bad = str(tmp_path / "absent.json")
    rc = cli.main(["search", "--run", "--gate-artifact", bad])
    assert rc == 2
    out = capsys.readouterr().out
    assert "REFUSED" in out and "No flag bypasses this refusal" in out

    failing = str(tmp_path / "failed.json")
    json.dump({"verdict": "FAIL", "passed": False, "n_promotions_total": 9,
               "n_catalogs": 30, "max_total_allowed": 3,
               "binomial_upper_tail_p": 0.0,
               "vacuity_control": {"verdict": "FIRES"}},
              open(failing, "w", encoding="utf-8"))
    rc2 = cli.main(["search", "--run", "--gate-artifact", failing])
    assert rc2 == 2
    assert "REFUSED" in capsys.readouterr().out


def test_v1_cap_is_retained_for_audit_but_is_not_the_verdict():
    """The cap survives as a reported number and cannot decide anything."""
    a_eff = GC.effective_alpha(0.10 / 30.0, 400)
    band = GC.calibration_band(900, a_eff, SG.FALSE_FAIL_RATE)
    assert band["k_fail_high"] > SG.MAX_TOTAL + 1
    v1_ff = GC.binom_sf(SG.MAX_TOTAL + 1, 900, a_eff)
    assert v1_ff > 0.10, "the v1 cap rejects a good instrument this often"
    assert band["achieved_false_fail_total"] <= SG.FALSE_FAIL_RATE + 1e-12


def test_n_catalogues_is_derived_from_the_declared_sensitivity():
    """v1's 30 was undefended. v2 computes N from the detection curve."""
    a_eff = GC.effective_alpha(0.10 / 30.0, SG.N_NULL)
    got = GC.required_catalogues(m=30, alpha_eff=a_eff,
                                 r_min=SG.DETECT_INFLATION_R,
                                 power=SG.DETECT_POWER,
                                 false_fail_rate=SG.FALSE_FAIL_RATE)
    assert got is not None
    assert got["detectable_inflation"] <= SG.DETECT_INFLATION_R + 1e-9
    assert got["n_catalogues"] > 30, "30 catalogues cannot reach the declared 3x"
