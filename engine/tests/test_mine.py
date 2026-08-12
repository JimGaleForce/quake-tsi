"""`mine` mode: ephemeris accuracy, null calibration, planted-effect recovery.

The headline test is the PLANTED one, mirroring test_planted.py for the ordinary
engine: a synthetic residual series with an injected 29.53-day cycle at a known
amplitude must be found by the miner, the harmonic ladder must pick the right rung,
and a scrambled control must come out clean.
"""

from __future__ import annotations

import datetime as _dt
import math

import numpy as np
import pytest

from engine import ephemeris as eph, mine as M

SEED = 20260811


# ------------------------------------------------------------- ephemeris ---
def test_ephemeris_known_syzygies():
    """New moon 2000-01-06 18:14 UTC and full moon 2000-01-21 04:40 UTC."""
    for when, want in ((_dt.datetime(2000, 1, 6, 18, 14), 0.0),
                       (_dt.datetime(2000, 1, 21, 4, 40), 180.0)):
        jd = eph.julian_day(when)
        d = (eph.moon_position(jd)["lon_deg"] - eph.sun_position(jd)["lon_deg"]) % 360.0
        err = min(abs(d - want), abs(d - want - 360.0), abs(d - want + 360.0))
        assert err < 0.5, f"{when}: elongation {d:.3f} deg, wanted {want} (+/-0.5)"


def test_ephemeris_cycle_lengths_and_ranges():
    t = eph.ephemeris_table(_dt.datetime(1995, 1, 1), 11544)
    for key, want in (("synodic_rad", eph.SYNODIC_MONTH),
                      ("anomalistic_rad", eph.ANOMALISTIC_MONTH),
                      ("draconic_rad", eph.DRACONIC_MONTH)):
        ph = np.unwrap(t[key])
        got = 2 * np.pi * (ph.size - 1) / (ph[-1] - ph[0])
        assert abs(got - want) < 0.01, f"{key}: period {got:.5f} d, wanted {want}"
    assert 355_000 < t["moon_dist_km"].min() < 360_000
    assert 404_000 < t["moon_dist_km"].max() < 408_000
    assert 28.0 < t["moon_dec_deg"].max() < 29.0
    assert abs(t["sun_dec_deg"].max() - 23.44) < 0.1
    assert np.isfinite(np.concatenate([v.ravel() for v in t.values()])).all()


def test_ephemeris_features_are_deterministic_in_t():
    """Causality exemption, demonstrated: the features do not touch the catalogue."""
    a = M.ephemeris_features(_dt.datetime(1995, 1, 1), 4000)
    b = M.ephemeris_features(_dt.datetime(1995, 1, 1), 4000)
    assert [f.name for f in a] == [f.name for f in b]
    for fa, fb in zip(a, b):
        assert np.array_equal(fa.values, fb.values)
        assert fa.causality_exempt and fa.family in (1, 2)


# ------------------------------------------------- statistics calibration ---
def _series(n=4000, amp=0.0, period=29.53, phase=0.7, seed=SEED, base=6.0):
    rng = np.random.default_rng(seed)
    t = np.arange(n, dtype=np.float64)
    off = base * (1.0 + 0.3 * np.sin(2 * np.pi * t / 900.0))
    lam = off * np.exp(amp * np.cos(2 * np.pi * t / period - phase))
    return t, off, rng.poisson(lam).astype(np.float64)


def test_glm_recovers_planted_amplitude():
    t, off, y = _series(amp=0.08)
    X = np.column_stack([np.sin(2 * np.pi * t / 29.53), np.cos(2 * np.pi * t / 29.53)])
    fit = M.glm_fit(X, y, off)
    assert fit["converged"]
    amp = math.hypot(*fit["beta"])
    assert abs(amp - 0.08) < 0.02, fit
    assert fit["bits_per_event"] > 0
    # the score statistic and the likelihood ratio must agree to ~1%
    S = M.score_stat_all_shifts(X, y, off)
    assert abs(S[0] - fit["lr_chi2"]) < 0.02 * fit["lr_chi2"]


def test_block_bootstrap_null_is_calibrated_and_shift_null_is_not():
    """The documented deviation, asserted rather than asserted-in-prose.

    Under H0 the block-bootstrap statistic tracks chi2_2. Under a planted cycle the
    CIRCULAR-SHIFT surrogates inherit the signal (their median stays huge), which is
    exactly why the miner uses the block bootstrap as the null for periodic features.
    """
    rng = np.random.default_rng(SEED)
    t, off, y0 = _series(amp=0.0)
    X = np.column_stack([np.sin(2 * np.pi * t / 29.53), np.cos(2 * np.pi * t / 29.53)])
    Sb = M.score_stat_block_bootstrap(X, y0, off, 400, rng, mean_block=60)
    assert 0.7 < np.median(Sb) / 1.386 < 1.8, np.median(Sb)     # chi2_2 median

    _t, off1, y1 = _series(amp=0.08)
    S = M.score_stat_all_shifts(X, y1, off1)
    assert S[0] > 50, S[0]
    assert np.median(S[30:-30]) > 0.3 * S[0], (
        "circular-shift surrogates no longer inherit the planted signal -- the "
        "documented degeneracy has changed and the null choice must be revisited")
    Sb1 = M.score_stat_block_bootstrap(X, y1, off1, 400, rng, mean_block=60)
    assert M.bootstrap_p(S[0], Sb1) < 0.01


def test_benjamini_hochberg():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.6])
    q, passed = M.benjamini_hochberg(p, 0.1)
    # textbook BH: p_(6) = 0.06 <= 0.1 * 6/10, so the first SIX are rejected
    assert passed[:6].all() and not passed[6:].any(), (q, passed)
    assert (np.diff(q[np.argsort(p)]) >= -1e-12).all()          # monotone
    assert M.benjamini_hochberg(np.array([0.5]), 0.1)[1].tolist() == [False]


def test_fold_lr_and_ladder_pick_the_true_rung():
    t, off, y = _series(n=6000, amp=0.10, period=29.53)
    lam0 = off * (y.sum() / off.sum())
    lad = M.harmonic_ladder(y, lam0, t, 29.53)
    assert abs(lad["winning_period_days"] - 29.53) < 1e-6, lad
    # a scan that landed on the 2nd harmonic must be walked back down to the truth
    lad2 = M.harmonic_ladder(y, lam0, t, 29.53 / 2.0)
    assert abs(lad2["winning_period_days"] - 29.53) < 1e-6, lad2
    assert "WINS over P" in lad2["verdict"]


# ------------------------------------------------------ the planted test ---
def test_planted_cycle_is_found_and_scrambled_control_is_clean():
    """PLANTED-EFFECT TEST (success criterion 3).

    Inject a 29.53 d cycle of known amplitude into a synthetic daily count series,
    run the miner's own three scoring paths over it, and require:
      * the GLM recovers the amplitude and the block-bootstrap p hits its floor;
      * the Lomb-Scargle scan puts a peak within 1% of 29.53 d;
      * the harmonic ladder reports the true period as the winning rung;
      * a scrambled control (same counts, time order destroyed) comes out clean on
        all three.
    """
    rng = np.random.default_rng(SEED)
    n, amp, period = 6000, 0.10, 29.53
    t, off, y = _series(n=n, amp=amp, period=period, seed=SEED)
    X = np.column_stack([np.sin(2 * np.pi * t / period), np.cos(2 * np.pi * t / period)])

    fit = M.glm_fit(X, y, off)
    assert abs(math.hypot(*fit["beta"]) - amp) < 0.02, fit
    S = M.score_stat_all_shifts(X, y, off)
    Sb = M.score_stat_block_bootstrap(X, y, off, 200, rng, mean_block=60)
    p_boot = M.bootstrap_p(S[0], Sb)
    assert p_boot <= 1.0 / 201.0 + 1e-12, p_boot

    lam0 = off * (y.sum() / off.sum())
    resid = y - lam0
    periods = np.exp(np.linspace(math.log(2.0), math.log(400.0), 900))
    peaks, _meta = M.period_scan(t, resid, periods, 200, rng, n_peaks=3, verbose=False)
    assert peaks, "period scan returned no peaks on a planted 10% cycle"
    top = peaks[0]
    assert abs(math.log(top["period_days"] / period)) < 0.01, top
    assert top["p_raw"] <= 0.01, top
    lad = M.harmonic_ladder(y, lam0, t, top["period_days"])
    assert abs(lad["winning_period_days"] - period) / period < 0.02, lad

    # ---- scrambled control: same marginal counts, no time structure ----
    perm = rng.permutation(n)
    y_s = y[perm]
    off_s = off[perm]
    fit_s = M.glm_fit(X, y_s, off_s)
    assert abs(math.hypot(*fit_s["beta"])) < 0.03, fit_s
    S_s = M.score_stat_all_shifts(X, y_s, off_s)
    Sb_s = M.score_stat_block_bootstrap(X, y_s, off_s, 200, rng, mean_block=60)
    assert M.bootstrap_p(S_s[0], Sb_s) > 0.05, "scrambled control passed the GLM test"
    lam0_s = off_s * (y_s.sum() / off_s.sum())
    resid_s = y_s - lam0_s
    peaks_s, meta_s = M.period_scan(t, resid_s, periods, 200, rng, n_peaks=3,
                                    verbose=False)
    # The control must not REDISCOVER THE PLANTED PERIOD. (Asserting that no peak
    # anywhere in a 900-point scan reaches p<0.05 would be a coin flip: the top
    # peak's p is uniform under the null by construction, so such a test flakes at
    # its own alpha. The scientific requirement is about the planted frequency.)
    assert all(abs(math.log(pk["period_days"] / period)) > 0.03 for pk in peaks_s), \
        f"scrambled control rediscovered the planted period: {peaks_s}"
    p_at_planted = float(M.lomb_scargle_power(t, resid_s, np.array([period]))[0])
    assert p_at_planted < meta_s["max_power_perm_p95"], (
        f"scrambled control retains power at the planted period: {p_at_planted:.5f}")


def test_mark_test_finds_planted_association_and_is_clean_on_noise():
    rng = np.random.default_rng(SEED)
    n = 20000
    theta = rng.random(n) * 2 * np.pi
    mark = 4.5 + 0.30 * np.cos(theta) + rng.normal(0, 0.35, n)
    hit = M.mark_test(theta, mark, "phase", 200, rng, block_events=200)
    assert hit["p_raw"] <= 0.02, hit
    assert hit["statistic"] > 0.2, hit
    null = M.mark_test(theta, rng.normal(4.5, 0.35, n), "phase", 200, rng,
                       block_events=200)
    assert null["p_raw"] > 0.05, null


def test_aliasing_audit_flags_a_lattice_artifact():
    """A 2-day alternating pattern is invisible at 2-day binning: LATTICE-SUSPECT."""
    rng = np.random.default_rng(SEED)
    n = 4000
    t = np.arange(n, dtype=np.float64)
    off = np.full(n, 6.0)
    y = rng.poisson(off * np.exp(0.15 * np.cos(np.pi * t))).astype(np.float64)
    f = M.Feature("alternating_2d", 2, "phase", np.mod(np.pi * t, 2 * np.pi),
                  "2-day alternation", period_hint=2.0)
    aud = M.aliasing_audit_glm(f, 0, slice(0, n), y, off, 100, rng)
    assert aud["verdict"] == "LATTICE-SUSPECT", aud
    assert aud["z_ratio_vs_1d"]["2d"] < 0.5, aud


def test_checkpoint_roundtrip_and_resume_detection(tmp_path):
    from engine import mine_session as ms

    cfg = {"preset": "unit-test", "n_surrogates": 1}
    ck_path = tmp_path / "session_x" / "checkpoint.json"
    ck = ms.Checkpoint(str(ck_path), cfg, "deadbeef")
    ck.save()
    ck.put("glm:foo", [{"p_raw": 0.5}])
    assert ck.done("glm:foo") and not ck.done("glm:bar")

    again = ms.Checkpoint.load(str(ck_path))
    assert again.done("glm:foo")
    assert again.state["config_hash"] == "deadbeef"
    assert not again.state["complete"]

    found = ms.find_resumable("deadbeef", root=str(tmp_path))
    assert found and found[0].endswith("session_x")
    again.state["complete"] = True
    again.save()
    assert ms.find_resumable("deadbeef", root=str(tmp_path)) is None


def test_report_banner_and_stub_caveats_are_present(tmp_path):
    from engine import mine_session as ms

    cfg = {"preset": "unit", "n_surrogates": 10, "fdr_q": 0.1, "engine_version": "x",
           "lags": [0]}
    tests = [{"feature": "moon_synodic_phase", "family": 1, "kind": "phase",
              "test": "glm_poisson_offset_etas", "lag": 0, "df": 2,
              "beta": [0.01, 0.02], "se": [0.01, 0.01], "amplitude_log_rate": 0.022,
              "pct_rate_modulation": 2.2, "bits_per_event": 1e-5, "p_raw": 0.0001,
              "bh_q": 0.001, "passes_fdr": True, "n_surrogates": 10}]
    path = ms.write_stubs(str(tmp_path), cfg, tests)
    import json
    with open(path, encoding="utf-8") as fh:
        payload = json.load(fh)
    assert payload["n_stubs"] == 1
    stub = payload["stubs"][0]
    assert any("GENERATOR, NOT EVIDENCE" in c for c in stub["caveats"])
    assert any("EXPECTED-AMPLITUDE FRAMING" in c for c in stub["caveats"])
    assert "Mignan & Broccardo 2019" in payload["standing_warning_eq24"]
    assert "holdout" in stub["next_step"]


# ------------------------------------------------- K-089-R tranche 1 (lags) ---
def test_ephemeris_features_default_lags_are_unchanged():
    """Backward compatibility, asserted: the default is a single lag-0 test."""
    feats = M.ephemeris_features(_dt.datetime(1995, 1, 1), 4000)
    assert len(feats) == 17
    assert all(f.lags == (0,) for f in feats)
    assert sum(1 for f in feats if f.kind == "phase") == 9
    assert sum(1 for f in feats if f.kind == "linear") == 8


def test_tranche1_lag_grid_selects_the_thirteen_and_counts_390_new_tests():
    """§P5-5: 8 linear x 30 = 240 plus 5 non-lag-free phase x 30 = 150."""
    feats = M.ephemeris_features(_dt.datetime(1995, 1, 1), 4000,
                                 lags=M.TRANCHE1_LAGS, lag_features="tranche1")
    assert M.TRANCHE1_LAGS == tuple(range(31))
    scanned = [f for f in feats if len(f.lags) > 1]
    unscanned = [f for f in feats if len(f.lags) == 1]
    assert len(scanned) == 13 and len(unscanned) == 4
    assert {f.name for f in unscanned} == set(M.LAG_FREE_PHASE_FEATURES)
    assert all(f.lags == M.TRANCHE1_LAGS for f in scanned)
    n_linear = sum(1 for f in scanned if f.kind == "linear")
    n_phase = sum(1 for f in scanned if f.kind == "phase")
    assert (n_linear, n_phase) == (8, 5)
    new = sum(len(f.lags) - 1 for f in feats)
    assert new == 13 * 30 == 390
    assert n_linear * 30 == 240 and n_phase * 30 == 150
    # total GLM rows the session will run over the cyclic features
    assert sum(len(f.lags) for f in feats) == 13 * 31 + 4 == 407

    audit = M.scan_axis_audit(feats, tranche1=True)
    assert audit["tranche"] == M.TRANCHE1_LABEL == "K-089-R tranche 1"
    assert audit["n_new_lag_tests_cyclic"] == 390
    assert len(audit["lag_axis_scanned"]) == 13
    assert len(audit["lag_axis_provably_free_not_scanned"]) == 4
    assert audit["lag_axis_neither_scanned_nor_free"] == []
    # without the flag, thirteen features are neither scanned nor provably free
    plain = M.scan_axis_audit(M.ephemeris_features(_dt.datetime(1995, 1, 1), 4000),
                              tranche1=False)
    assert len(plain["lag_axis_neither_scanned_nor_free"]) == 13
    assert plain["n_new_lag_tests_cyclic"] == 0


def test_tranche1_lag_features_are_rejected_when_unknown():
    with pytest.raises(ValueError):
        M.ephemeris_features(_dt.datetime(1995, 1, 1), 500, lags=(0, 1),
                             lag_features=["no_such_feature"])


def test_lag_invariance_audit_reproduces_the_free_scan_exception():
    """Clause 3 as §P5-5 restates it: FREE is proved numerically, per feature."""
    feats = M.ephemeris_features(_dt.datetime(1995, 1, 1), 8081)
    rows = M.lag_invariance_audit(feats, slice(100, 8081))
    assert len(rows) == 9
    by = {r["feature"]: r for r in rows}
    for name in M.LAG_FREE_PHASE_FEATURES:
        r = by[name]
        assert r["declared"] == "FREE"
        assert r["worst_min_r2"] >= 1.0 - M.LAG_FREE_TOL, r
        assert r["measured_free"] and r["ok"], r
    # the five built on true solar/lunar longitudes are NOT rotations under a lag
    for name in ("moon_synodic_phase", "spring_neap_phase", "perigean_spring_beat",
                 "eclipse_year_beat", "annual_synodic_beat"):
        r = by[name]
        assert r["declared"] == "PRICED"
        assert not r["measured_free"], r
        assert r["ok"], r
    assert by["spring_neap_phase"]["worst_min_r2"] < 0.95
    assert all(r["verdict"] == "OK" for r in rows)
