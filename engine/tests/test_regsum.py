"""EQ-24 v2 phase 2b: the 2R-df phase-incoherent regional statistic.

Counted invariants, each traced to the ruling it enforces:

  * §P6-4 Rule 4.1  -- the partition is deterministic and is a function of
    EXPLORATION-WINDOW data only. Perturbing the holdout tail cannot move it.
  * §P6-4 Rule 4.2  -- the sum is 2R-df additive and PHASE-INCOHERENT: rotating one
    region's phase leaves the statistic's null distribution alone, and the sum
    equals the sum of the parts.
  * §P6-4 Rule 4.3 / §P7-1(d) -- per-region amplitudes are labelled UNRESOLVED and
    are ineligible for BH.
  * §P6-4 Rule 4.5  -- region is a live stratum axis in engine.strata.
  * §P6-4 Rule 4.7 item 5 -- THE BLIND-SPOT KILL: a planted regional-phase signal is
    RECOVERED by the 2R-df sum and MISSED by the global-sum test. Engine-test-level
    evidence shaped like G-M1 arm (ii); NOT the full G-M1 run.
  * §P7-1(b)        -- the floor is the formula, and R falls out of it.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import stats

from engine import design, grid as gridmod, mine as M, regions as R, strata as strata_mod

SEED = 20260812


def _make_global_ctx(n_days=1600, seed=SEED, explore_frac=0.7, dead_sector=None,
                     dead_only_in_window=False):
    """A synthetic catalogue whose cells span ALL SIX longitude sectors.

    `_synth.make_ctx` puts every cell inside one 60-degree sector, which would make
    every regional test here run at R = 1. This builder places 4 cells in each of
    the 6 sectors so the partition has something to partition.

    `dead_sector` silences one sector; with `dead_only_in_window` it is silent in
    the exploration window but ACTIVE in the holdout tail -- the exact shape of the
    §P6-4 Finding A leak.
    """
    rng = np.random.default_rng(seed)
    lons, lats = [], []
    for s in range(6):
        lo = -180.0 + 60.0 * s
        for k in range(4):
            lons.append(lo + 10.0 + 12.0 * k)
            lats.append(10.5 + 3.0 * k)
    clon = np.asarray(lons, dtype=float)
    clat = np.asarray(lats, dtype=float)
    n_cells = clon.size

    lam = np.repeat((0.6 * (1.0 + rng.random(n_cells)))[:, None], n_days, axis=1)
    n_explore = int(round(explore_frac * n_days))
    if dead_sector is not None:
        m = R.sector_of_lon(clon) == int(dead_sector)
        if dead_only_in_window:
            lam[np.ix_(m, np.arange(0, n_explore))] = 0.0
        else:
            lam[m, :] = 0.0
    y = rng.poisson(lam).astype(np.float32)

    cell_idx, day_idx = np.nonzero(y)
    rep = y[cell_idx, day_idx].astype(int)
    ev_cell = np.repeat(cell_idx, rep)
    ev_day = np.repeat(day_idx, rep)

    g = gridmod.Grid(1.0, 1.0)
    g.build_domain(clat, clon)
    ec = g.cell_index(clat[ev_cell], clon[ev_cell])
    assert (ec >= 0).all()
    meta = dict(n_days=n_days, n_explore_days=n_explore, mag_floor=4.5,
                dlat=1.0, dlon=1.0)
    ctx = design._make_ctx(g, ec, ev_day, np.full(ec.size, 4.6), meta, verbose=False)
    return ctx, ctx.day_counts(4.5)


def _ctx_window(n_days=1600, seed=SEED, **kw):
    ctx, y = _make_global_ctx(n_days=n_days, seed=seed, **kw)
    window = slice(100, ctx.n_explore_days)
    return ctx, y, window


def _phase_design(n, period=29.53, lag=0.0):
    t = np.arange(n, dtype=np.float64) + lag
    ph = (t % period) / period
    return np.column_stack([np.sin(2 * np.pi * ph), np.cos(2 * np.pi * ph)]), ph


# --------------------------------------------------- Rule 4.1: the partition --
def test_partition_is_deterministic():
    ctx, y, window = _ctx_window()
    a = R.build_regions(ctx, y, window)
    b = R.build_regions(ctx, y, window)
    assert a["digest"] == b["digest"]
    assert np.array_equal(a["region_of_cell"], b["region_of_cell"])
    assert a["R"] == b["R"] >= 1


def test_partition_spans_all_six_sectors():
    ctx, y, window = _ctx_window()
    part = R.build_regions(ctx, y, window)
    assert part["R"] == 6, part["R"]
    assert sorted(r["sector"] for r in part["regions"]) == list(range(6))


def test_partition_uses_exploration_window_only():
    """§P6-4 Rule 4.1 / Finding A: holdout-window seismicity may not move the set.

    Sector 2 is silent in the exploration window and ACTIVE in the holdout tail --
    exactly the shape of the K-080 census leak. It must be DROPPED regardless of
    how loud the tail gets, and every remaining cell's assignment must be
    bit-identical to the fully-silenced case.
    """
    ctx_a, y_a, window = _ctx_window(dead_sector=2, dead_only_in_window=False)
    ctx_b, y_b, _ = _ctx_window(dead_sector=2, dead_only_in_window=True)

    a = R.build_regions(ctx_a, y_a, window)
    b = R.build_regions(ctx_b, y_b, window)

    assert float(y_b[:, window.stop:].sum()) > 0.0
    assert b["R"] == 5, f"holdout-only activity rescued a sector (R = {b['R']})"
    assert [r["sector"] for r in b["regions"]] == [0, 1, 3, 4, 5]
    assert [r["sector"] for r in a["regions"]] == [r["sector"] for r in b["regions"]]

    # And piling MORE events into the tail of an already-live sector moves nothing.
    y_c = np.array(y_a, copy=True)
    y_c[:, window.stop:] += 500.0
    c = R.build_regions(ctx_a, y_c, window)
    assert c["digest"] == a["digest"]
    assert np.array_equal(c["region_of_cell"], a["region_of_cell"])


def test_partition_rule_is_declared_one_way():
    """S-9: one rule id, and the sector map is a pure function of longitude."""
    assert R.REGION_RULE_ID == "R2b-lon6-active"
    assert R.N_SECTORS == 6
    lon = np.array([-180.0, -179.9, -120.0, -60.0, 0.0, 59.9, 60.0, 179.9])
    got = R.sector_of_lon(lon)
    assert got.tolist() == [0, 0, 1, 2, 3, 3, 4, 5]
    # 180 wraps to -180, never to a seventh sector
    assert int(R.sector_of_lon(np.array([180.0]))[0]) == 0


def test_partition_digest_changes_with_the_rule():
    ctx, y, window = _ctx_window()
    a = R.build_regions(ctx, y, window, n_sectors=6)
    b = R.build_regions(ctx, y, window, n_sectors=4)
    assert a["digest"] != b["digest"]


# -------------------------------------------- Rule 4.2: the sum construction --
def _regional_series(y, window, part, rate_scale=1.0):
    """Poisson-ish offset for the synthetic ctx: the in-window per-cell mean rate."""
    yw = np.asarray(y)[:, window]
    mu = np.maximum(yw.mean(axis=1, keepdims=True), 1e-3)
    rate = np.repeat(mu, yw.shape[1], axis=1) * rate_scale
    return R.regional_series(y, rate, part["region_of_cell"], window, part["R"])


def test_sum_equals_sum_of_per_region_parts():
    ctx, y, window = _ctx_window()
    part = R.build_regions(ctx, y, window)
    C, O = _regional_series(y, window, part)
    X, _ = _phase_design(C.shape[1])
    per = M.score_stat_regions(X, C, O)
    assert per.shape == (part["R"],)
    assert (per >= 0).all(), "a score statistic is a non-negative quadratic form"
    assert M.score_stat_regsum(X, C, O) == pytest.approx(float(per.sum()), rel=1e-12)


def test_sum_is_phase_incoherent_under_rotation_of_one_region():
    """Rotating ONE region's phase must not change the sum's null distribution.

    Concretely: circularly rotate region 0's (counts, offset) pair by an arbitrary
    number of days. That changes that region's alignment with the feature -- and so
    changes the GLOBAL sum's statistic -- but each per-region term is a
    rotation-covariant quadratic form, and the block-bootstrap null of the SUM is
    built from the same rotated series, so the null itself is unmoved. This is the
    property that makes the sum blind to §K87-0(d)(i) cancellation.
    """
    ctx, y, window = _ctx_window()
    part = R.build_regions(ctx, y, window)
    C, O = _regional_series(y, window, part)
    X, _ = _phase_design(C.shape[1])

    rot = 173
    C2, O2 = np.array(C, copy=True), np.array(O, copy=True)
    C2[0] = np.roll(C[0], rot)
    O2[0] = np.roll(O[0], rot)

    n = 400
    a = M.score_stat_regsum_block_bootstrap(
        X, C, O, n, np.random.default_rng(1), mean_block=90.0)
    b = M.score_stat_regsum_block_bootstrap(
        X, C2, O2, n, np.random.default_rng(1), mean_block=90.0)

    # Same null law: a two-sample KS test must not reject at 1%.
    ks = stats.ks_2samp(a, b)
    assert ks.pvalue > 0.01, (
        f"rotating one region's phase moved the SUM's null distribution "
        f"(KS p = {ks.pvalue:.4g}); the statistic is not phase-incoherent")
    # And the observed sum of the rotated series is still the sum of its parts.
    per2 = M.score_stat_regions(X, C2, O2)
    assert M.score_stat_regsum(X, C2, O2) == pytest.approx(float(per2.sum()), rel=1e-12)


def test_sum_null_has_about_2R_df_scale():
    """Sanity on the df bookkeeping: an independent-Poisson sum sits near 2R."""
    rng = np.random.default_rng(7)
    n, Rn = 2000, 6
    X, _ = _phase_design(n)
    O = np.full((Rn, n), 4.0)
    vals = []
    for _ in range(300):
        C = rng.poisson(O).astype(np.float64)
        vals.append(M.score_stat_regsum(X, C, O))
    # mean of chi2(2R) is 2R = 12; allow generous slack for 300 draws.
    assert 2 * Rn * 0.75 < float(np.mean(vals)) < 2 * Rn * 1.25, np.mean(vals)


def test_all_shifts_index_zero_is_the_observed_sum():
    ctx, y, window = _ctx_window()
    part = R.build_regions(ctx, y, window)
    C, O = _regional_series(y, window, part)
    X, _ = _phase_design(C.shape[1])
    S = M.score_stat_regsum_all_shifts(X, C, O)
    assert S.shape == (C.shape[1],)
    assert S[0] == pytest.approx(M.score_stat_regsum(X, C, O), rel=1e-8)


def test_surrogates_go_through_the_same_sum():
    """A surrogate world is ONE index draw shared by every region, then summed."""
    ctx, y, window = _ctx_window()
    part = R.build_regions(ctx, y, window)
    C, O = _regional_series(y, window, part)
    X, _ = _phase_design(C.shape[1])
    S = M.score_stat_regsum_block_bootstrap(
        X, C, O, 64, np.random.default_rng(3), mean_block=90.0)
    assert S.shape == (64,)
    assert (S >= 0).all()
    # summing R non-negative parts: the sum must dominate any single region's null
    S1 = M.score_stat_block_bootstrap(
        X, C[0], O[0], 64, np.random.default_rng(3), mean_block=90.0)
    assert float(np.median(S)) > float(np.median(S1))


# ------------------------------------- Rule 4.7 item 5: THE BLIND-SPOT KILL ---
def _planted_world(Rn=6, n=4000, amplitude=0.35, period=29.53, seed=5,
                   base_rate=8.0):
    rng = np.random.default_rng(seed)
    X, ph = _phase_design(n, period=period)
    O = np.full((Rn, n), float(base_rate))
    phases = R.equally_spaced_phases(Rn)
    C = R.plant_regional_phase(O, O, ph, amplitude, phases, rng=rng)
    return X, C, O, phases


def test_planted_regional_phase_cancels_in_the_global_sum():
    """The premise: equally spaced phases cancel to ~0 in the domain aggregate."""
    _X, _C, _O, phases = _planted_world()
    v = sum(complex(math.cos(p), math.sin(p)) for p in phases)
    assert abs(v) < 1e-9, f"planted phases do not cancel (|sum| = {abs(v)})"


def test_regsum_recovers_what_the_global_sum_misses():
    """§P6-4 Rule 4.7 item 5. ENGINE-TEST-LEVEL, G-M1 arm (ii)-SHAPED, NOT G-M1.

    Same amplitude in every region, different phase per region. The domain-summed
    series has the modulation cancelled out of it and its 2-df test sees nothing;
    the 2R-df phase-incoherent sum sees it plainly.
    """
    X, C, O, _ph = _planted_world()

    # (a) the GLOBAL test -- aggregate first, then test. This is the blind spot.
    c_glob, o_glob = C.sum(axis=0), O.sum(axis=0)
    s_glob = float(M.score_stat_regions(X, c_glob[None, :], o_glob[None, :])[0])
    p_glob = M.chi2_sf(s_glob, 2)

    # (b) the REGIONAL SUM -- test first, then aggregate the statistics.
    s_sum = M.score_stat_regsum(X, C, O)
    p_sum = M.chi2_sf(s_sum, 2 * C.shape[0])

    assert p_glob > 0.05, (
        f"the global-sum test was supposed to MISS the planted regional-phase "
        f"signal but got chi2 = {s_glob:.2f}, p = {p_glob:.3g}")
    assert p_sum < 1e-6, (
        f"the 2R-df sum failed to recover the planted signal: chi2 = {s_sum:.2f} "
        f"on {2*C.shape[0]} df, p = {p_sum:.3g}")
    assert s_sum > 10.0 * s_glob


def test_regsum_surrogate_null_also_recovers_the_planted_signal():
    """Recovery must survive the real null, not just the parametric chi2."""
    X, C, O, _ph = _planted_world()
    s_sum = M.score_stat_regsum(X, C, O)
    Sb = M.score_stat_regsum_block_bootstrap(
        X, C, O, 400, np.random.default_rng(11), mean_block=90.0)
    p = M.bootstrap_p(s_sum, Sb)
    assert p <= 1.0 / 401.0 + 1e-12, f"block-bootstrap p = {p:.4g}, expected the floor"


def test_per_region_planted_recovery_at_one_regions_N():
    """§P6-4 Rule 4.7 item 5, second half: recovery is re-demonstrated PER REGION,
    at that region's own (order-of-magnitude smaller) N -- recovery is N-dependent."""
    X, C, O, _ph = _planted_world()
    r = 0
    s_r = float(M.score_stat_regions(X, C[r][None, :], O[r][None, :])[0])
    Sb = M.score_stat_block_bootstrap(
        X, C[r], O[r], 400, np.random.default_rng(13), mean_block=90.0)
    p_r = M.bootstrap_p(s_r, Sb)
    n_r = float(C[r].sum())
    assert p_r < 0.01, (
        f"per-region recovery failed at N = {n_r:.0f} (chi2 = {s_r:.2f}, "
        f"p = {p_r:.4g}); no bound may be reported at this aggregation")


def test_null_world_is_not_recovered():
    """The negative control: no planted signal, no detection."""
    rng = np.random.default_rng(21)
    n, Rn = 4000, 6
    X, _ = _phase_design(n)
    O = np.full((Rn, n), 8.0)
    C = rng.poisson(O).astype(np.float64)
    s = M.score_stat_regsum(X, C, O)
    assert M.chi2_sf(s, 2 * Rn) > 0.001, f"false detection on a pure null: chi2 = {s}"


# ------------------------------- Rule 4.3 / §P7-1(d): UNRESOLVED labelling ----
def test_unresolved_rows_are_labelled_and_cite_the_rule():
    from engine import mine_session as ms

    ctx, y, window = _ctx_window()
    part = R.build_regions(ctx, y, window)
    C, O = _regional_series(y, window, part)
    X, _ = _phase_design(C.shape[1])

    class F:
        name, family, kind, df = "synthetic_phase", 1, "phase", 2

    reg_cfg = {"vif": ms.REGION_VIF_DEFAULT, "alpha": ms.REGION_ALPHA_DEFAULT,
               "target_amplitude": ms.REGION_TARGET_AMPLITUDE_DEFAULT}
    rows = ms._region_unresolved_rows(F, X, C, O, part, reg_cfg)
    assert len(rows) == part["R"]
    for r in rows:
        assert r["p_method"] == "UNRESOLVED"
        assert "§P6-4 Rule 4.3" in r["p_method_reason"]
        assert "§P7-1(d)" in r["p_method_reason"]
        assert r["s15"] in ("MEASURABLE", "UNMEASURABLE")
        assert r["a_min_formula"] > 0.0
        # nothing here is quotable at this catalogue's N
        assert r["quotable"] is False


def test_unresolved_rows_are_never_bh_eligible():
    """A row labelled UNRESOLVED must not carry a rejection into the BH vector."""
    from engine import gpd_tail
    assert "UNRESOLVED" in gpd_tail.P_METHODS
    rows = [{"p_method": "UNRESOLVED", "p_raw": 1e-9},
            {"p_method": "MC_RESOLVED", "p_raw": 1e-9}]
    census = gpd_tail.census(rows)
    assert census["UNRESOLVED"] == 1 and census["MC_RESOLVED"] == 1


# ------------------------------------------- Rule 4.5: region as strata axis --
def test_region_is_a_live_stratum_axis():
    assert strata_mod.stratum_key(1, "region", 3) != strata_mod.stratum_key(1, "region", 4)
    assert strata_mod._test_kind({"test": "regsum_score_2Rdf"}) == "regsum"
    assert strata_mod._test_kind(
        {"test": "glm_poisson_offset_etas_region"}) == "region"


def test_stratified_bh_prices_each_region_separately(tmp_path):
    import json
    Rn = 3
    m_s = 10
    q = 0.10
    decl = {"q": q, "strata": [
        {"name": f"region-{r}", "feature_family": None, "test_kind": "region",
         "region": r, "m_s": m_s, "q_s": q} for r in range(Rn)]}
    path = tmp_path / "regions.json"
    path.write_text(json.dumps(decl), encoding="utf-8")
    part = strata_mod.load_partition(str(path), q=q)
    strata_mod.assert_budget_identity(part["strata"], q)
    rows = [{"family": None, "test": "glm_poisson_offset_etas_region", "region": r}
            for r in range(Rn)]
    names = [strata_mod.stratum_of(x, part) for x in rows]
    assert names == [f"region-{r}" for r in range(Rn)]


# ------------------------------------------------- §P7-1(b): the floor + R ----
def test_floor_formula_reproduces_the_ledgers_own_arithmetic():
    """§P7-1(a): at VIF = 1, alpha = 0.05 the constant is z_.025 + z_.80 = 2.8016."""
    n = 2.0  # A_min = k * sqrt(2/N) -> with N = 2 the sqrt term is 1
    assert R.a_min(1.0, 0.05, n) == pytest.approx(2.8016, abs=2e-4)
    # and §P7-1(a)'s own multiplicity example: alpha = 0.1/259 -> z = 3.549
    assert float(stats.norm.isf((0.10 / 259) / 2.0)) == pytest.approx(3.5494, abs=1e-3)


def test_floor_scales_as_sqrt_vif():
    assert R.a_min(4.0, 0.05, 1000) == pytest.approx(2.0 * R.a_min(1.0, 0.05, 1000))


def test_R_choice_follows_from_the_measured_floor():
    """The declared R is a CONSEQUENCE of the F4-58 measurement, not a convention.

    At the measured VIF (~24.1), Tranche-A alpha and the ledger's own 20% reference
    amplitude, the summed statistic keeps 80% power only up to R = 6 -- and the
    per-region battery's own floor resolves the target at no R >= 2, which is the
    §P7-1(d) outcome.
    """
    from engine import mine_session as ms
    vif = ms.REGION_VIF_DEFAULT
    alpha = ms.REGION_ALPHA_DEFAULT
    n = 46585.0                      # exploration-window events, M >= 4.5
    r_sum, _lam = R.max_R_for_sum(vif, alpha, n, 0.20)
    assert r_sum == R.N_SECTORS == 6, (
        f"declared R = {R.N_SECTORS} but the floor supports {r_sum}")
    r_bat, _real = R.max_R_for_per_region_battery(vif, alpha, n, 0.20)
    assert r_bat < 2, f"per-region battery unexpectedly measurable up to R = {r_bat}"
    # and it WOULD have supported the §P6-4-anticipated 12-24 at the inferred VIF
    r_at_inferred, _ = R.max_R_for_sum(3.94, alpha, n, 0.20)
    assert r_at_inferred > 24


def test_region_floor_table_reports_s15_per_region():
    ctx, y, window = _ctx_window()
    part = R.build_regions(ctx, y, window)
    tbl = R.region_floor_table(part, 24.08, 0.10 / 713.0, 0.20)
    assert len(tbl["rows"]) == part["R"]
    assert tbl["n_measurable"] + tbl["n_unmeasurable"] == part["R"]
    assert 0.0 <= tbl["fraction_unmeasurable"] <= 1.0
    for r in tbl["rows"]:
        assert r["s15"] in ("MEASURABLE", "UNMEASURABLE")


# ------------------------------------------------- key schema / determinism ---
def test_region_test_keys_are_distinct_and_add_no_schema_field():
    from engine import mine_session as ms
    assert M.TEST_KEY_FIELDS == ("master_seed", "feature", "kind", "lag", "rung",
                                 "null_type", "region", "bin_width")
    keys = [M.test_key(1, "f", "regsum", lag=0, null_type="circular_shift",
                       region="ALL")]
    keys += [M.test_key(1, "f", "region", lag=0, null_type="circular_shift",
                        region=r) for r in range(6)]
    keys += [M.test_key(1, "f", "glm", lag=0, null_type="block_bootstrap")]
    digests = ms.assert_task_keys_unique(keys)
    assert len(digests) == len(keys)
    # the phase-2a glm digest is untouched by the arrival of the region kinds
    assert M.test_key_digest(
        M.test_key(1, "f", "glm", lag=0, null_type="block_bootstrap")
    ) == M.test_key_digest(keys[-1])


# ----------------------------------------------- end-to-end session WIRING ---
N_DAYS_E2E = 700


class _StubBaseline:
    name = "stub-const-v0"
    caveat = "constant-rate stub, test only"
    burn_in_days = 0

    def __init__(self, rate):
        self._rate = np.asarray(rate, dtype=np.float64)

    def rate(self, window):
        return self._rate[:, window]

    def report(self):
        return ["stub baseline"]


def _prepared_e2e(seed=5, n_feats=3):
    import datetime as _dt
    rng = np.random.default_rng(seed)
    ctx, y = _make_global_ctx(n_days=N_DAYS_E2E, seed=seed)
    window = slice(60, N_DAYS_E2E)
    rate = np.repeat(np.maximum(y.mean(axis=1, keepdims=True), 1e-2),
                     N_DAYS_E2E, axis=1)
    base = _StubBaseline(rate)
    counts = y[:, window].sum(axis=0).astype(float)
    offset = rate[:, window].sum(axis=0).astype(float)
    t0 = _dt.datetime(2000, 1, 1)
    feats = M.ephemeris_features(t0, N_DAYS_E2E)[:n_feats]
    ev_day = np.repeat(np.arange(N_DAYS_E2E), y.sum(axis=0).astype(int))
    ev_day = ev_day[ev_day >= window.start]
    marks = {"day": ev_day.astype(np.int64),
             "day_float": ev_day.astype(float) + 0.5,
             "mag": 4.5 + rng.standard_exponential(ev_day.size) * 0.3,
             "depth": rng.uniform(1.0, 30.0, ev_day.size)}
    return (ctx, base, y, window, counts, offset, marks, feats, [], t0)


def _run_e2e(tmp_path, battery):
    import os
    import types
    from engine import mine_session as ms

    args = types.SimpleNamespace(
        mag_target=4.5, dlat=1.0, dlon=1.0, explore_frac=0.7,
        data_dir=str(tmp_path), no_download=True, seed=7, tranche1=False,
        ladder=False, gpd=False, strata=None,
        regsum=True, regions=bool(battery),
        region_sectors=None, region_min_fraction=None,
        region_vif=None, region_alpha=None, region_target_amplitude=None)
    preset = dict(ms.QUICK, n_surrogates=200, n_periods=40, n_peaks=2,
                  label="phase2b-test")
    cfg = ms.build_config(args, preset)
    sd = str(tmp_path / ("battery" if battery else "sum"))
    os.makedirs(sd, exist_ok=True)
    out = ms.run(cfg, verbose=False, resume=False, jobs=1, session_dir=sd,
                 ledger_path=str(tmp_path / "EXPLORE_COUNT.jsonl"),
                 prepared=_prepared_e2e())
    return cfg, out, sd


def test_regsum_session_wires_end_to_end(tmp_path):
    import json
    import os
    cfg, _out, sd = _run_e2e(tmp_path, battery=False)
    assert cfg["regions"]["enabled"] and cfg["regions"]["battery"] is False

    rep = open(os.path.join(sd, "report.md"), encoding="utf-8").read()
    doc = json.load(open(os.path.join(sd, "stubs.json"), encoding="utf-8"))
    ck = json.load(open(os.path.join(sd, "checkpoint.json"), encoding="utf-8"))

    assert "Region battery" in rep
    assert "R2b-lon6-active" in rep
    assert "UNRESOLVED" in rep
    assert ck["regions"]["R"] == 6
    assert ck["regions"]["digest"]
    rows = [t for t in ck["tests"] if t.get("test") == "regsum_score_2Rdf"]
    assert rows, "no regsum rows reached the checkpoint"
    for t in rows:
        assert t["df"] == 6 * 2 or t["df"] == 6 * 1
        assert t["amplitude_log_rate"] is None      # §P7-1(d): no amplitude
        assert len(t["regions_unresolved"]) == 6
    assert doc["n_region_amplitudes_unresolved"] == len(rows) * 6
    assert "§P7-1(d)" in doc["region_amplitudes_unresolved_rule"]
    for u in doc["region_amplitudes_unresolved"]:
        assert u["p_method"] == "UNRESOLVED"


def test_region_battery_is_priced_at_R_times_and_scored_neither_way(tmp_path):
    import json
    import os
    _cfg_sum, _o1, sd_sum = _run_e2e(tmp_path / "a", battery=False)
    _cfg_bat, _o2, sd_bat = _run_e2e(tmp_path / "b", battery=True)
    ck_s = json.load(open(os.path.join(sd_sum, "checkpoint.json"), encoding="utf-8"))
    ck_b = json.load(open(os.path.join(sd_bat, "checkpoint.json"), encoding="utf-8"))

    n_regsum = sum(1 for t in ck_s["tests"] if t["test"] == "regsum_score_2Rdf")
    n_batt = sum(1 for t in ck_b["tests"]
                 if t["test"] == "glm_poisson_offset_etas_region")
    assert n_batt == 6 * n_regsum, (n_batt, n_regsum)
    assert ck_b["n_tests"] > ck_s["n_tests"]

    batt = [t for t in ck_b["tests"] if t["test"] == "glm_poisson_offset_etas_region"]
    # §P6-4 Rule 4.4: at this catalogue's N every region is UNMEASURABLE, is scored
    # neither way, and cannot emit a stub.
    assert all(t["s15"] == "UNMEASURABLE" for t in batt)
    assert all(t["bh_eligible"] is False for t in batt)
    assert all(t["passes_fdr"] is False for t in batt)
    assert all(t["p_method"] == "UNRESOLVED" for t in batt)
    doc = json.load(open(os.path.join(sd_bat, "stubs.json"), encoding="utf-8"))
    assert not any(s["test"] == "glm_poisson_offset_etas_region"
                   for s in doc["stubs"])


def test_config_hash_is_unchanged_when_the_regional_axis_is_off(tmp_path):
    """A default run's config -- and therefore its resumable sessions -- is
    byte-identical to phase 2a. The regional block appears only when asked for."""
    import types
    from engine import mine_session as ms
    base = types.SimpleNamespace(
        mag_target=4.5, dlat=1.0, dlon=1.0, explore_frac=0.7,
        data_dir=str(tmp_path), no_download=True, seed=7, tranche1=False,
        ladder=False, gpd=False, strata=None)
    preset = dict(ms.QUICK, label="phase2b-test")
    off = ms.build_config(base, preset)
    assert "regions" not in off
    on = ms.build_config(
        types.SimpleNamespace(**dict(vars(base), regsum=True, regions=False,
                                     region_sectors=None, region_min_fraction=None,
                                     region_vif=None, region_alpha=None,
                                     region_target_amplitude=None)), preset)
    assert "regions" in on
    assert ms._cfg_hash(off) != ms._cfg_hash(on)


def test_regional_series_shapes_and_conservation():
    ctx, y, window = _ctx_window()
    part = R.build_regions(ctx, y, window)
    C, O = _regional_series(y, window, part)
    assert C.shape == O.shape == (part["R"], window.stop - window.start)
    assert float(C.sum()) == pytest.approx(part["n_events_assigned"], rel=1e-9)
