"""The G-M1 + R1 gate harness: sandboxing, plant construction, acceptance arithmetic.

What is unit-testable here is exactly the part of `engine/gate_r1.py` that is NOT
the pipeline. The pipeline itself is `mine_session.run` and is already pinned by
test_phase2a_session / test_mine_parallel / test_regsum; re-testing it here would
be duplicate coverage. What IS new, and what fails silently if wrong, is:

  1. SANDBOXING -- a simulated catalog that reaches engine/EXPLORE_COUNT.jsonl or
     engine/out/mine would put sim hashes into the real multiplicity accounting
     and the real §P7-8(c)(5) build-invariant registry. Three separations,
     asserted before any catalog runs.
  2. PLANT CONSTRUCTION ABOVE THE FLOOR -- §P7-8(d). A plant below 2x the
     operative floor produces a miss that reads as an instrument failure and is
     really a power failure. Also that the plant is in the ESTIMATOR'S units:
     A in exp(A cos(theta - phi)) is what the GLM's |beta| estimates.
  3. ACCEPTANCE ARITHMETIC -- the three R1 conditions and the 80% recovery bar.
     These decide the verdict, so a sign error in them is the whole gate.
"""

import json
import math
import os

import numpy as np
import pytest

from engine import floors
from engine import gate_r1 as G
from engine import mine as M


# =============================================================== 1. sandboxing =
def test_sandbox_accepts_a_session_under_the_gate_sandbox():
    info = G.assert_sandboxed(os.path.join(G.SANDBOX_SESSIONS, "session_null_000"),
                              G.SANDBOX_LEDGER)
    assert info["session_dir"].endswith("session_null_000")
    # the build-invariant registry root is dirname(session_dir) -- run() puts it
    # there -- so it must itself be inside the sandbox.
    assert info["artifact_registry_root"].endswith("gate_r1/sessions")


def test_sandbox_refuses_the_real_mine_directory():
    with pytest.raises(AssertionError, match="gate sandbox"):
        G.assert_sandboxed(os.path.join("engine", "out", "mine", "session_x"),
                           G.SANDBOX_LEDGER)


def test_sandbox_refuses_the_real_explore_ledger():
    with pytest.raises(AssertionError, match="REAL exploration ledger"):
        G.assert_sandboxed(os.path.join(G.SANDBOX_SESSIONS, "session_a"),
                           os.path.join("engine", "EXPLORE_COUNT.jsonl"))


def test_sandbox_refuses_a_ledger_outside_the_sandbox(tmp_path):
    with pytest.raises(AssertionError, match="not under the sandbox"):
        G.assert_sandboxed(os.path.join(G.SANDBOX_SESSIONS, "session_a"),
                           str(tmp_path / "somewhere_else.jsonl"))


def test_gate_config_is_the_full_pipeline_and_flat_bh():
    cfg = G.gate_config()
    assert cfg["ladder"]["n_max"] == G.GATE_SURROGATES == cfg["n_surrogates"]
    assert cfg["gpd"]["enabled"] is True
    assert cfg["gpd"]["calibration"]["pass"] is True      # §P6-2(7) licence
    assert cfg["regions"]["enabled"] is True
    assert cfg["regions"]["battery"] is False
    assert "strata" not in cfg                            # flat (default) BH
    assert cfg["fdr_q"] == M.FDR_Q


def test_the_monte_carlo_floor_permits_a_bh_rejection():
    """The reason the budget is 10,000 and not the quick preset's 200.

    A gate in which BH can never reject at any effect size is not a gate. The
    smallest attainable p is 1/(N_max+1); BH's most stringent rung is q/m. At
    the realised m (~157 tests) the 200-surrogate floor of 1/201 sits ABOVE
    q/m and the arithmetic is vacuous -- so this pins the inequality rather
    than the constant.
    """
    m = 157                       # the realised declared count for this config
    floor_at_200 = 1.0 / 201.0
    floor_at_gate = 1.0 / (G.GATE_SURROGATES + 1.0)
    assert floor_at_200 > M.FDR_Q / m           # vacuous
    assert floor_at_gate < M.FDR_Q / m          # can reject


# ========================================================= 2. plant construction
class _Feat:
    """A stand-in Feature: only `values`, `kind`, `name` are read by the planter."""

    def __init__(self, name, values, period_hint=None):
        self.name = name
        self.kind = "phase"
        self.values = np.asarray(values, dtype=np.float64)
        self.period_hint = period_hint
        self.block_days = 60.0


def _prep(n_days=2000, n_cells=8, rate=3.0):
    theta = np.mod(2 * np.pi * np.arange(n_days) / 29.53, 2 * np.pi)
    window = slice(0, n_days)
    rate_w = np.full((n_cells, n_days), rate, dtype=np.float64)
    return {
        "feats": [_Feat("moon_synodic_phase", theta, 29.53),
                  _Feat("annual_phase",
                        np.mod(2 * np.pi * np.arange(n_days) / 365.24,
                               2 * np.pi), 365.24)],
        "window": window, "rate_w": rate_w,
        "offset": rate_w.sum(axis=0),
        "y": np.zeros((n_cells, n_days), dtype=np.int64),
        "counts": rate_w.sum(axis=0),
        "marks": {"mag": np.full(500, 4.7), "depth": np.full(500, 10.0)},
    }


def test_global_plant_sits_at_exactly_two_floors_and_says_so():
    prep = _prep()
    mod, rep = G.global_plant(prep, "moon_synodic_phase")
    n = float(prep["offset"].sum())
    assert rep["N"] == pytest.approx(n)
    assert rep["planted_amplitude"] == pytest.approx(
        floors.PLANT_FACTOR * floors.a_min(rep["vif"], rep["alpha"], n))
    assert rep["amplitude_over_floor"] == pytest.approx(floors.PLANT_FACTOR)
    assert rep["compliant"] is True
    # the report carries the floor next to the plant, which is the §P7-8(d) demand
    assert rep["operative_floor_A_min"] > 0
    assert rep["vif"] == floors.MEASURED_VIF_DF2_PHASE
    assert mod.shape == prep["rate_w"].shape


def test_a_plant_below_the_floor_is_refused_not_silently_run():
    prep = _prep()
    n = float(prep["offset"].sum())
    below = 0.5 * floors.a_min(floors.MEASURED_VIF_DF2_PHASE,
                               floors.ALPHA_TRANCHE_A, n)
    with pytest.raises(floors.PlantBelowFloor, match="below"):
        G.global_plant(prep, "moon_synodic_phase", amplitude=below)


def test_plant_preserves_the_expected_event_count():
    """The I0(A) normalisation: the plant moves PHASE, not the total rate.

    Without it a planted catalog carries more events than a null one and the
    'recovery' would be partly a rate change the offset never saw.
    """
    prep = _prep()
    mod, _rep = G.global_plant(prep, "moon_synodic_phase")
    planted_total = float((prep["rate_w"] * mod).sum())
    null_total = float(prep["rate_w"].sum())
    assert planted_total == pytest.approx(null_total, rel=2e-3)


def test_plant_amplitude_is_in_the_estimators_own_units():
    """A in exp(A cos(theta - phi)) IS what the GLM's |beta| estimates.

    The phase design is the UNSTANDARDISED [sin, cos], so log lambda =
    a + A sin(phi) sin(theta) + A cos(phi) cos(theta) and hypot(beta) = A. This
    is what makes the recovered/planted ratio comparable to G-M1's [0.8, 1.2],
    so it is pinned against the engine's own glm_fit rather than assumed.
    """
    rng = np.random.default_rng(7)
    n = 6000
    theta = np.mod(2 * np.pi * np.arange(n) / 29.53, 2 * np.pi)
    A, phi = 0.30, 0.7
    base = np.full(n, 60.0)
    lam = base * np.exp(A * np.cos(theta - phi)) / float(
        __import__("scipy.special", fromlist=["i0"]).i0(A))
    counts = rng.poisson(lam).astype(np.float64)
    X = np.column_stack([np.sin(theta), np.cos(theta)])
    fit = M.glm_fit(X, counts, base)
    amp = float(np.hypot(*np.asarray(fit["beta"])))
    assert amp == pytest.approx(A, rel=0.06)
    phi_hat = math.atan2(fit["beta"][0], fit["beta"][1])
    assert abs(math.degrees(phi_hat - phi)) < 15.0     # G-M1's own phase bar


def test_regional_plant_uses_each_regions_own_N_and_opposes_phases():
    """Rule 4.7 item 5: re-run PER REGION, because recovery is N-dependent.

    Regions with fewer events must therefore get LARGER plants, and the phases
    must be spread so the domain sum cancels -- that cancellation is the
    §K87-0(d)(i) blind spot the regsum exists to kill, and planting it is the
    only way arm (ii) tests anything the global arm does not.
    """
    prep = _prep(n_cells=6)
    R = 3
    roc = np.array([0, 0, 0, 1, 2, 2])
    # region 1 is a single cell -> smallest N -> largest required plant
    part = {"R": R, "region_of_cell": roc}
    mod, reps = G.regional_plant(prep, part, "moon_synodic_phase")
    assert len(reps) == R
    for r in reps:
        assert r["compliant"] is True
        assert r["amplitude_over_floor"] == pytest.approx(floors.PLANT_FACTOR)
    ns = [r["N"] for r in reps]
    amps = [r["planted_amplitude"] for r in reps]
    # smaller N -> bigger plant, strictly
    assert [a for _, a in sorted(zip(ns, amps))] == sorted(amps, reverse=True)
    phases = sorted(r["phase_rad"] for r in reps)
    assert phases == pytest.approx([0.0, 2 * math.pi / 3, 4 * math.pi / 3])
    # the domain-summed modulation is far flatter than any single region's
    dom = (prep["rate_w"] * mod).sum(axis=0)
    per_region = (prep["rate_w"][roc == 1] * mod[roc == 1]).sum(axis=0)
    assert dom.std() / dom.mean() < 0.5 * (per_region.std() / per_region.mean())


def test_simulated_catalog_is_a_true_null_by_construction():
    """Poisson from the fitted lambda; marks redrawn independently of day."""
    prep = _prep(n_days=1500, n_cells=5, rate=2.0)
    rng = np.random.default_rng(3)
    y, counts, marks = G.simulate_catalog(prep, rng)
    assert y.shape == prep["y"].shape
    assert counts.sum() == marks["day"].size
    assert counts.sum() == pytest.approx(prep["offset"].sum(), rel=0.05)
    # the mark pool is the real one, resampled -- so every drawn mark is a real
    # value and the (day, mark) pairing carries no information
    assert set(np.unique(marks["mag"])).issubset(set(np.unique(prep["marks"]["mag"])))
    assert (np.diff(marks["day"]) >= 0).all()
    assert ((marks["day_float"] - marks["day"] >= 0)
            & (marks["day_float"] - marks["day"] < 1)).all()


# ====================================================== 3. acceptance arithmetic
def _cat(n_tests=157, n_surv=0, p=0.5, gpd=(0, 0), mc=(0, 0)):
    return {"kind": "null", "n_tests": n_tests, "n_bh_survivors": n_surv,
            "max_stat_p": p,
            "p_method_counts": {"GPD_EXTRAPOLATED": {"n": gpd[0],
                                                     "n_survivors": gpd[1]},
                                "MC_RESOLVED": {"n": mc[0],
                                                "n_survivors": mc[1]}}}


def test_r1a_is_mean_survivors_against_q_times_m():
    cats = [_cat(n_surv=s) for s in [0, 0, 1, 0, 2]]
    r = G.r1_survivor_arithmetic(cats)
    assert r["mean_bh_survivors"] == pytest.approx(0.6)
    assert r["bound_q_times_m"] == pytest.approx(15.7)
    assert r["pass"] is True
    # and it FAILS when the mean exceeds the bound, not merely when a max does
    bad = [_cat(n_surv=20) for _ in range(5)]
    assert G.r1_survivor_arithmetic(bad)["pass"] is False


def test_r1b_gates_on_d_plus_only_and_reports_d_minus():
    """§P6-9(a). D- is permitted by super-uniformity and must NOT gate."""
    # a perfectly conservative (super-uniform) sample: every p pushed high.
    conservative = [0.5 + 0.5 * (i + 1) / 30 for i in range(30)]
    r = G.ks_uniform_one_sided(conservative)
    assert r["d_plus"] <= r["d_crit"]
    assert r["pass"] is True                       # conservatism is permitted
    assert r["d_minus_reported_not_gating"] > r["d_crit"]   # and would have failed
    assert r["two_sided_p_not_gating"] < 0.01      # two-sided would reject
    # the critical value is the §P6-9(a) formula, checked at the ledger's own
    # worked point: a = 0.01, n = 20000 -> 0.010730
    assert G.ks_uniform_one_sided([0.5] * 20000, a=0.01)["d_crit"] == \
        pytest.approx(0.010730, abs=5e-6)


def test_r1b_fails_an_anti_conservative_sample():
    anti = list(np.linspace(0.0005, 0.05, 30))
    r = G.ks_uniform_one_sided(anti)
    assert r["pass"] is False
    assert r["d_plus"] > r["d_crit"]


def test_r1c_flags_an_elevated_gpd_survivor_rate_and_only_that():
    flat = [_cat(gpd=(100, 5), mc=(1000, 50))]
    assert G.gpd_vs_mc_survivor_rate(flat)["pass"] is True
    hot = [_cat(gpd=(200, 60), mc=(2000, 40))]
    r = G.gpd_vs_mc_survivor_rate(hot)
    assert r["elevated"] is True and r["pass"] is False
    assert r["fisher_one_sided_p_gpd_greater"] < 0.05
    # a GPD arm that was never exercised is REPORTED, not read as a pass
    none = [_cat(gpd=(0, 0), mc=(1000, 3))]
    z = G.gpd_vs_mc_survivor_rate(none)
    assert z["n_gpd_extrapolated"] == 0
    assert z["gpd_rate"] is None
    assert "never exercised" in z["note"]


def _rec(detected, ratio=None):
    return {"planted_recovery": {"detected": detected, "amplitude_ratio": ratio,
                                 "best_p_raw": 1e-4, "best_p_floor": 1e-4}}


def test_recovery_rate_bar_is_80_percent_and_amplitude_band_is_08_to_12():
    eight = [_rec(True, 1.0)] * 8 + [_rec(False)] * 2
    r = G.recovery_rate(eight)
    assert r["recovery_rate"] == pytest.approx(0.8)
    assert r["detection_pass"] is True and r["pass"] is True
    seven = [_rec(True, 1.0)] * 7 + [_rec(False)] * 3
    assert G.recovery_rate(seven)["pass"] is False
    # detection can pass while the AMPLITUDE fails, and the gate is the AND
    out_of_band = [_rec(True, 1.45)] * 10
    r2 = G.recovery_rate(out_of_band)
    assert r2["detection_pass"] is True
    assert r2["amplitude_pass"] is False and r2["pass"] is False


def test_unresolved_amplitudes_assert_detection_only():
    """§P7-1(d): where the pipeline reports UNRESOLVED by rule, the arm asserts
    the detection and never an amplitude ratio it was never given."""
    recs = [_rec(True, None)] * 10
    r = G.recovery_rate(recs)
    assert r["amplitude_quoted"] is False
    assert r["amplitude_pass"] is None
    assert r["pass"] is True


def test_a_miss_is_labelled_power_or_instrument_never_left_ambiguous():
    plants_ok = [{"amplitude_over_floor": 2.0, "compliant": True}]
    d = G.miss_diagnosis([_rec(False)], plants_ok)
    assert d["n_misses"] == 1 and d["verdict"].startswith("INSTRUMENT")
    plants_bad = [{"amplitude_over_floor": 0.9, "compliant": False}]
    d2 = G.miss_diagnosis([_rec(False)], plants_bad)
    assert d2["verdict"].startswith("POWER")
    assert G.miss_diagnosis([_rec(True)], plants_ok)["n_misses"] == 0


def test_max_stat_T_is_minus_log10_p_not_a_chi2():
    """The scale the max-statistic actually lives on.

    `strata.max_statistic_p` documents its own statistic as "max over tests of
    -log10(empirical p within own null)". Its `t_obs` is therefore ~3.9 where the
    chi2 of the same test is ~2000; comparing a chi2 to `t_obs` would make every
    attribution false, which is why the conversion is a named function.
    """
    row = {"p_circular_shift": 1.0 / 7657.0, "chi2_score": 2030.0}
    assert G.max_stat_T(row) == pytest.approx(math.log10(7657.0), abs=1e-9)
    assert G.max_stat_T(row) == pytest.approx(3.8841, abs=1e-3)
    assert G.max_stat_T({"p_circular_shift": None}) is None


def _grow(feature, p_shift, chi2, passes=False, amp=0.30):
    return {"feature": feature, "test": "glm_poisson_offset_etas", "lag": 0,
            "p_raw": 1e-4, "chi2_score": chi2, "passes_fdr": passes,
            "amplitude_log_rate": amp, "p_method": "MC_RESOLVED",
            "p_floor": 1e-4, "p_circular_shift": p_shift, "bh_eligible": True}


def test_max_statistic_detection_requires_attaining_the_family_maximum():
    """Single-step Westfall-Young, not a family p read as a per-test one."""
    floor = 1.0 / 7657.0
    planted = {"feature": "moon_synodic_phase",
               "test_kind": "glm_poisson_offset_etas", "truth_amplitude": 0.30}
    t_max = math.log10(7657.0)
    tests = [_grow("moon_synodic_phase", floor, 2030.0),
             _grow("other", 0.02, 9.0, amp=0.01)]
    hit = G._planted_recovery(tests, {"max_stat_t_obs": t_max,
                                      "max_stat_p": 0.00013}, planted)
    assert hit["detected"] is True and hit["max_stat_detection"] is True
    assert hit["amplitude_ratio"] == pytest.approx(1.0)

    # a TIE at the resolution ceiling still counts: both tests have the same
    # Westfall-Young adjusted p, so a tie-break would be arbitrary
    tied = [_grow("moon_synodic_phase", floor, 2030.0),
            _grow("sun_moon_elongation", floor, 2008.0, amp=0.21)]
    assert G._planted_recovery(tied, {"max_stat_t_obs": t_max,
                                      "max_stat_p": 0.00013},
                               planted)["max_stat_detection"] is True

    # the family p is significant but the maximum belongs to ANOTHER test
    below = [_grow("moon_synodic_phase", 0.01, 2030.0),
             _grow("other", floor, 9.0, amp=0.01)]
    miss = G._planted_recovery(below, {"max_stat_t_obs": t_max,
                                       "max_stat_p": 0.00013}, planted)
    assert miss["max_stat_detection"] is False and miss["detected"] is False

    # and a significant family p is still required
    weak = G._planted_recovery(tests, {"max_stat_t_obs": t_max,
                                       "max_stat_p": 0.4}, planted)
    assert weak["detected"] is False


def test_regsum_rows_are_not_eligible_for_the_max_statistic_route():
    """`max_statistic_matrix` admits glm rows only; the regional arm is BH."""
    rows = [{"feature": "moon_synodic_phase", "test": "regsum_score_2Rdf",
             "lag": 0, "p_raw": 1e-4, "chi2_score": 2027.0, "passes_fdr": False,
             "p_method": "MC_RESOLVED", "p_floor": 1e-4,
             "p_circular_shift": 1.0 / 7657.0}]
    planted = {"feature": "moon_synodic_phase",
               "test_kind": "regsum_score_2Rdf", "truth_amplitude": None}
    r = G._planted_recovery(rows, {"max_stat_t_obs": math.log10(7657.0),
                                   "max_stat_p": 1e-4}, planted)
    assert r["max_stat_covered"] is False
    assert r["max_stat_detection"] is False
    assert r["detected"] is False            # BH is the only route, and it failed
    rows[0]["passes_fdr"] = True
    r2 = G._planted_recovery(rows, {"max_stat_t_obs": 3.9, "max_stat_p": 1e-4},
                             planted)
    assert r2["detected"] is True and r2["bh_survivor"] is True
    assert r2["amplitude_ratio"] is None     # §P7-1(d): UNRESOLVED by rule


def test_a_bh_ineligible_row_still_reports_its_amplitude():
    """§P6-2(1)/(3): a failed GPD gate suppresses the P, never the EFFECT SIZE.

    This is a real case from the gate run: a plant so strong that the GPD's
    CI-upper landed 26 decades below the one-decade cap, so the row was ruled
    UNRESOLVED and made BH-INELIGIBLE -- while still estimating the planted
    amplitude to within 0.3%. Conflating the two labels would have thrown that
    measurement away and recorded a G-M1 miss that did not happen.
    """
    row = _grow("moon_synodic_phase", 1.0 / 7657.0, 2030.38, amp=0.29969)
    row.update({"p_method": "UNRESOLVED", "bh_eligible": False,
                "p_method_reason": "§P6-2(3) one-decade cap"})
    r = G._planted_recovery([row], {"max_stat_t_obs": math.log10(7657.0),
                                    "max_stat_p": 0.00013},
                            {"feature": "moon_synodic_phase",
                             "test_kind": "glm_poisson_offset_etas",
                             "truth_amplitude": 0.2990})
    assert r["bh_survivor"] is False and r["bh_eligible"] is False
    assert r["detected"] is True                       # via the max-statistic
    assert r["amplitude_unresolved"] is False
    assert r["amplitude_ratio"] == pytest.approx(1.0023, abs=1e-3)


def test_bh_survival_alone_counts_as_recovery():
    tests = [_grow("annual_phase", 0.5, 40.0, passes=True, amp=0.27)]
    r = G._planted_recovery(tests, {"max_stat_t_obs": 99.0, "max_stat_p": 0.9},
                            {"feature": "annual_phase",
                             "test_kind": "glm_poisson_offset_etas",
                             "truth_amplitude": 0.30})
    assert r["bh_survivor"] is True and r["detected"] is True
    assert r["max_stat_detection"] is False
    assert r["amplitude_ratio"] == pytest.approx(0.9)
