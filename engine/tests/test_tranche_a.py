"""Tranche A (§P7-2): the F9-20 control battery, the declared arithmetic, R3/F10-25.

These pin the DECLARATION, not a result. The whole point of §P7-2(a) is that the
control battery is priced, and a battery that silently changed size -- by gaining a
mark axis, by losing a feature, by drifting off 31 lags -- would be a different
declaration than the one recorded in EXPLORE_COUNT.jsonl.
"""

import json
import math
import os

import numpy as np
import pytest

from engine import mine as M, mine_session as ms, strata as strata_mod, tranche_a as TA


def _donors(n_days=2000, n=20):
    """A stand-in real feature list: mixed kinds, mixed period hints, mixed families."""
    rng = np.random.default_rng(7)
    t = np.arange(n_days, dtype=float)
    out = []
    for i in range(n):
        if i % 2:
            out.append(M.Feature(f"donor_phase_{i}", 1, "phase",
                                 np.mod(2 * np.pi * t / (20.0 + i), 2 * np.pi),
                                 period_hint=20.0 + i, lags=(0,)))
        else:
            out.append(M.Feature(f"donor_lin_{i}", 4, "linear",
                                 np.cumsum(rng.normal(size=n_days)) / 30.0
                                 + np.sin(2 * np.pi * t / (50.0 + i)),
                                 lags=(0,)))
    return out


# ---------------------------------------------------- the declared arithmetic --
def test_declared_price_is_713_and_is_23_features_times_31_lags():
    assert M.F9_20_N_NAMED + M.F9_20_N_MATCHED == 23
    assert len(M.F9_20_LAGS) == 31
    assert M.F9_20_N_DECLARED_TESTS == 713
    assert TA.N_PRICED_NEW == 713
    assert TA.N_DECLARED_VECTOR == TA.N_PRICED_NEW + TA.N_ALREADY_DECLARED == 1263


def test_partition_file_declares_the_two_arms_and_satisfies_the_budget_identity():
    part = strata_mod.load_partition(TA.STRATA_FILE, q=0.10)
    names = {s["name"]: s for s in part["strata"]}
    assert names["f9_20_negative_controls"]["m_s"] == M.F9_20_N_DECLARED_TESTS
    assert names["f9_20_negative_controls"]["feature_family"] == M.F9_20_FAMILY
    assert names["f9_20_negative_controls"]["test_kind"] == "glm"
    assert part["catch_all"] == "real_declared"
    assert part["m"] == TA.N_DECLARED_VECTOR
    ident = strata_mod.assert_budget_identity(part["strata"], 0.10)
    assert ident["ok"] and math.isclose(ident["sum_ms_qs"], ident["m_q"])


def test_battery_builds_exactly_713_tests_and_refuses_a_short_donor_list():
    f = M.negative_control_features(None, 2000, _donors(), seed=11)
    assert len(f) == 23
    assert sum(len(x.lags) for x in f) == 713
    assert all(x.control and x.family == M.F9_20_FAMILY for x in f)
    with pytest.raises(ValueError, match="matched controls"):
        M.negative_control_features(None, 2000, _donors(n=5), seed=11)


def test_battery_is_reproducible_from_the_declared_seed_and_moves_with_it():
    a = M.negative_control_features(None, 2000, _donors(), seed=11)
    b = M.negative_control_features(None, 2000, _donors(), seed=11)
    c = M.negative_control_features(None, 2000, _donors(), seed=12)
    assert all(np.array_equal(x.values, y.values) for x, y in zip(a, b))
    assert not all(np.array_equal(x.values, y.values) for x, y in zip(a, c))


def test_matched_controls_carry_the_donor_null_and_block_length():
    donors = _donors()
    ctl = M.negative_control_features(None, 2000, donors, seed=11)
    by_donor = {c.control_of: c for c in ctl if c.control_of in
                {d.name for d in donors}}
    assert len(by_donor) == M.F9_20_N_MATCHED
    for d in donors[:M.F9_20_N_MATCHED]:
        c = by_donor[d.name]
        # same null (periodic decides which null is primary), same df, same
        # bootstrap block length -- otherwise the two arms are not comparable.
        assert c.periodic == d.periodic
        assert c.kind == d.kind and c.df == d.df
        assert math.isclose(c.block_days, d.block_days, rel_tol=1e-9)


def test_phase_randomised_surrogate_preserves_the_power_spectrum_exactly():
    rng = np.random.default_rng(3)
    x = np.cumsum(rng.normal(size=1024)) + np.sin(np.arange(1024) / 7.0)
    y = M.phase_randomised_surrogate(x, rng)
    fx, fy = np.abs(np.fft.rfft(x - x.mean())), np.abs(np.fft.rfft(y - y.mean()))
    assert np.allclose(fx, fy, atol=1e-8)
    assert not np.allclose(x, y)                       # phases really were redrawn


def test_a_control_is_never_the_donor_of_another_control():
    donors = _donors()
    ctl = M.negative_control_features(None, 2000, donors, seed=11)
    assert M.negative_control_features(None, 2000, donors + ctl, seed=11)[3].control_of \
        == donors[0].name


# ------------------------------------------------------------- the read-outs --
def _fake_rows(n_real=40, n_ctl=60, seed=5):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_real):
        rows.append({"test": "glm_poisson_offset_etas", "feature": f"r{i // 4}",
                     "family": 1, "lag": i % 4, "kind": "linear",
                     "amplitude_log_rate": float(abs(rng.normal(0, 0.05))),
                     "p_raw": float(rng.uniform()), "passes_fdr": False,
                     "stratum": "real_declared", "order_key": [0, i, 0, 0, "r"]})
    for i in range(n_ctl):
        rows.append({"test": "glm_poisson_offset_etas", "feature": f"nc{i // 4}",
                     "family": M.F9_20_FAMILY, "lag": i % 4, "kind": "linear",
                     "amplitude_log_rate": float(abs(rng.normal(0, 0.05))),
                     "p_raw": float(rng.uniform()), "passes_fdr": False,
                     "stratum": "f9_20_negative_controls",
                     "order_key": [0, i, 0, 0, "nc"]})
    for t in rows:
        t["p_bh"] = t["p_raw"]
    return rows


def test_control_calibration_counts_both_arms_and_prints_q_times_m():
    rows = _fake_rows()
    part = strata_mod.load_partition(TA.STRATA_FILE, q=0.10)
    cc = ms.control_calibration(rows, part["strata"])
    assert cc["control_arm"]["n_declared"] == 713
    assert cc["real_arm"]["n_declared"] == 550
    assert cc["control_arm"]["expected_false_by_chance_q_times_m"] == pytest.approx(71.3)
    assert cc["F10_25_survivor_ratio"] is None          # zero-denominator, not 0/0


def test_winners_curse_labels_and_deduplicates_the_lag_axis():
    rows = _fake_rows()
    wc = ms.winners_curse_report(rows)
    assert wc["label_applied_to_every_amplitude"] == ms.SELECTION_BIASED
    assert len(wc["real_arm_top_k_features"]) == wc["selection_size_k"]
    # the whole point of the deduplication: one feature cannot occupy two slots
    assert len(set(wc["real_arm_top_k_features"])) == wc["selection_size_k"]
    assert len(set(wc["null_arm_top_k_features"])) == wc["selection_size_k"]
    assert wc["row_level_secondary"]["k"] >= wc["selection_size_k"]


def test_s15_uses_the_stratum_alpha_and_the_per_feature_vif_when_it_exists():
    rows = _fake_rows()
    part = strata_mod.load_partition(TA.STRATA_FILE, q=0.10)
    out = ms.s15_by_stratum(rows, 46585.0, part["strata"],
                            vif_table={"r0": {"vif": 100.0}})
    by = {r["stratum"]: r for r in out["rows"]}
    assert by["f9_20_negative_controls"]["alpha_s"] == pytest.approx(0.10 / 713)
    assert by["real_declared"]["alpha_s"] == pytest.approx(0.10 / 550)
    # the measured VIF really was used, and it really did raise that row's floor
    r0 = [t for t in rows if t["feature"] == "r0"][0]
    other = [t for t in rows if t["feature"] == "r1"][0]
    assert "per-feature" in r0["a_min_vif_source"]
    assert "fallback" in other["a_min_vif_source"]
    assert r0["a_min_formula"] > other["a_min_formula"]
    # §P7-8(d)'s published Tranche A floor, reproduced from the formula
    assert by["f9_20_negative_controls"]["A_min_max"] == pytest.approx(0.149, abs=2e-3)


def test_rank_stability_is_perfect_against_itself_and_labels_its_own_arm():
    rows = _fake_rows()
    same = ms.rank_stability(rows, rows, label="self")
    assert same["spearman_rho"] == pytest.approx(1.0)
    assert same["verdict"] == "STABLE"
    shuffled = [dict(t) for t in rows]
    rng = np.random.default_rng(1)
    for t in shuffled:
        t["p_bh"] = float(rng.uniform())
    noisy = ms.rank_stability(rows, shuffled, label="shuffled")
    assert noisy["spearman_rho"] < 0.5
    assert noisy["verdict"].startswith("UNSTABLE")


# --------------------------------------------------- the F10-24 data resample --
def test_block_resample_preserves_lag_structure_inside_blocks():
    rng = np.random.default_rng(2)
    idx = TA.block_resample_index(1000, 100, rng)
    assert idx.size == 1000
    # contiguity: at most one break per block
    breaks = int((np.diff(idx) != 1).sum())
    assert breaks <= int(np.ceil(1000 / 100))


def test_resample_prepared_keeps_every_series_on_one_index():
    n_days, w0, w1 = 900, 100, 800
    t = np.arange(n_days, dtype=float)
    feats = [M.Feature("f", 1, "linear", np.sin(t / 11.0), lags=(0, 5, 30),
                       period_hint=11.0)]
    counts = np.arange(w1 - w0, dtype=float)
    offset = np.ones(w1 - w0)
    d = np.arange(w0, w1, 3)
    marks = {"day": d, "day_float": d + 0.25,
             "mag": np.full(d.size, 4.6), "depth": np.full(d.size, 10.0)}
    prepared = (None, None, None, slice(w0, w1), counts, offset, marks, feats,
                {}, None)
    (_, _, _, win, c_r, o_r, m_r, f_r, _, _), meta = TA.resample_prepared(
        prepared, 42, block=200, max_lag=30)
    assert c_r.size == o_r.size == counts.size
    assert win.start == 30 and win.stop == counts.size + 30
    assert f_r[0].values.size == counts.size + 30
    # the design at lag 0 lines up with the resampled counts, by construction
    assert f_r[0].design(win, 0).shape[0] == c_r.size
    assert f_r[0].design(win, 30).shape[0] == c_r.size
    # marks stay addressable: run() recovers `day - window.start` as the ordinal
    fe_day = m_r["day"] - win.start
    assert fe_day.min() >= 0 and fe_day.max() < c_r.size
    assert meta["block_boundary_fraction"] == pytest.approx(30 / 200)
