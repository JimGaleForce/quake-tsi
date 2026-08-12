"""End-to-end: a whole mine session with §P6-2 and §P6-3 switched on.

The unit tests in test_gpd_tail.py / test_strata.py pin the estimator and the
procedure. This one pins the WIRING -- that the labels, the census, the budget
table, the max-statistic and the candidate list actually reach `report.md` and
`stubs.json`, which is where §P6-2(5) and §P6-3(4) say they must appear.

It runs on a synthetic ETAS-free substrate with a stub baseline so it costs seconds,
not the production sweep's minutes.
"""

import datetime as _dt
import json
import os
import types

import numpy as np
import pytest

from engine import mine as M
from engine import mine_session as ms
from engine import strata as S

N_DAYS = 800
N_FEATS = 4


class _StubBaseline:
    """Everything `mine_session.run` asks of a baseline, and nothing else."""
    name = "stub-const-v0"
    caveat = "constant-rate stub, test only"
    burn_in_days = 0

    def __init__(self, rate):
        self._rate = np.asarray(rate, dtype=np.float64)

    def rate(self, window):
        return self._rate[None, window]

    def report(self):
        return ["stub baseline"]


def _prepared(seed=5):
    rng = np.random.default_rng(seed)
    day = np.arange(N_DAYS, dtype=float)
    offset = 0.8 * (1.0 + 0.2 * np.sin(2 * np.pi * day / 365.25))
    counts = rng.poisson(offset).astype(float)
    # the window starts inside the record, exactly as the production burn-in does,
    # so the lag-invariance audit's 30-day shifts have somewhere to shift to
    window = slice(60, N_DAYS)
    base = _StubBaseline(offset)
    t0 = _dt.datetime(2000, 1, 1)
    feats = M.ephemeris_features(t0, N_DAYS)[:N_FEATS]
    ev_day = np.repeat(np.arange(N_DAYS), counts.astype(int))
    ev_day = ev_day[ev_day >= window.start]
    marks = {"day": ev_day.astype(np.int64),
             "day_float": ev_day.astype(float) + 0.5,
             "mag": 4.5 + rng.standard_exponential(ev_day.size) * 0.3,
             "depth": rng.uniform(1.0, 30.0, ev_day.size)}
    ctx = types.SimpleNamespace(n_days=N_DAYS)
    y = counts[None, :]
    return (ctx, base, y, window, counts[window], offset[window], marks, feats,
            [], t0)


def _calibration(tmp_path, passing=True):
    """A §P6-2(7) verdict artifact. `--gpd` refuses to run without a passing one."""
    p = tmp_path / ("cal_pass.json" if passing else "cal_reject.json")
    p.write_text(json.dumps({
        "rule": "§P6-2(7)", "pass": bool(passing), "coverage": 0.96,
        "coverage_bar": 0.90, "worst_understatement": 1.4,
        "understatement_bar": 3.0, "n_comparisons": 27, "n_fitted": 27}),
        encoding="utf-8")
    return str(p)


def _cfg(tmp_path, gpd=True, strata=None, calibration=None):
    args = types.SimpleNamespace(
        mag_target=4.5, dlat=1.0, dlon=1.0, explore_frac=0.7,
        data_dir=str(tmp_path), no_download=True, seed=7, tranche1=False,
        ladder=False, gpd=gpd, gpd_confirm_max_n=20000, strata=strata,
        gpd_calibration=(calibration or _calibration(tmp_path)))
    preset = dict(ms.QUICK, n_surrogates=400, n_periods=60, n_peaks=3,
                  label="phase2a-test")
    cfg = ms.build_config(args, preset)
    return cfg


def _run(tmp_path, name, cfg, prepared):
    sd = str(tmp_path / name)
    os.makedirs(sd, exist_ok=True)
    return ms.run(cfg, verbose=False, resume=False, jobs=1, session_dir=sd,
                  ledger_path=str(tmp_path / "EXPLORE_COUNT.jsonl"),
                  prepared=prepared), sd


@pytest.fixture(scope="module")
def flat_session(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("flat")
    out, sd = _run(tmp, "s", _cfg(tmp, gpd=True), _prepared())
    return out, sd


def test_every_row_of_report_and_stubs_carries_a_p_method(flat_session):
    out, sd = flat_session
    rep = open(os.path.join(sd, "report.md"), encoding="utf-8").read()
    doc = json.load(open(os.path.join(sd, "stubs.json"), encoding="utf-8"))
    # §P6-2(5): the census, by name, with all three numbers
    assert "P-METHOD CENSUS" in rep
    for m in ("MC_RESOLVED", "GPD_EXTRAPOLATED", "UNRESOLVED"):
        assert m in rep
        assert m in doc["p_method_census"]
    assert sum(doc["p_method_census"].values()) == doc["n_tests"]
    assert all(s.get("p_method") for s in doc["stubs"])


def test_banner_carries_the_declared_count_and_expected_false_discoveries(
        flat_session):
    """§P6-4 rule 4.7 item 1."""
    out, sd = flat_session
    rep = open(os.path.join(sd, "report.md"), encoding="utf-8").read()
    assert "DECLARED TEST COUNT" in rep
    assert "q x m" in rep
    m = out["n_tests"]
    assert f"{0.10 * m:.1f}" in rep


def test_banner_carries_the_resolution_floor_count(flat_session):
    """§P6-4 rule 4.7 item 4, second half: the K-of-M resolvability count."""
    out, sd = flat_session
    rep = open(os.path.join(sd, "report.md"), encoding="utf-8").read()
    assert "resolution floor ABOVE that" in rep
    assert f"of {out['n_tests']} declared tests" in rep


def test_max_statistic_is_reported_globally_and_never_a_stratum_alone(
        flat_session):
    """§P6-3(4)."""
    out, sd = flat_session
    rep = open(os.path.join(sd, "report.md"), encoding="utf-8").read()
    ck = json.load(open(os.path.join(sd, "checkpoint.json"), encoding="utf-8"))
    ms_rep = ck["max_statistic"]
    assert "GLOBAL max-statistic p" in rep
    assert "GLOBAL p (adjacent" in rep
    assert ms_rep["global"]["n_tests"] == ms_rep["n_covered"]
    for r in ms_rep["per_stratum"]:
        assert r["global_p"] == ms_rep["global"]["p"]


def test_default_is_unstratified_and_the_identity_is_still_asserted(flat_session):
    out, sd = flat_session
    ck = json.load(open(os.path.join(sd, "checkpoint.json"), encoding="utf-8"))
    tbl = ck["bh"]["table"]
    assert len(tbl) == 1 and tbl[0]["stratum"] == S.UNSTRATIFIED
    assert ck["bh"]["identity"]["ok"] is True
    rep = open(os.path.join(sd, "report.md"), encoding="utf-8").read()
    assert "UNSTRATIFIED (flat BH)" in rep


def test_stratified_run_prints_the_budget_table(tmp_path):
    """§P6-3(1): the (stratum, m_s, q_s, threshold_s) table, from a real session."""
    prepared = _prepared()
    out0, sd0 = _run(tmp_path, "pilot", _cfg(tmp_path, gpd=False), prepared)
    ck = json.load(open(os.path.join(sd0, "checkpoint.json"), encoding="utf-8"))
    counts = {}
    for t in ck["tests"]:
        counts[S.stratum_key(t["family"], S._test_kind(t))] = \
            counts.get(S.stratum_key(t["family"], S._test_kind(t)), 0) + 1
    keys = sorted(counts)
    m = sum(counts.values())
    # a REAL reallocation: the first stratum gets 1.5x the flat budget and the
    # remainder is taken from the others, exactly as §P6-3(2) requires it to be
    # -- with the justification written into the file before the run.
    q = 0.10
    first = keys[0]
    q_first = 1.5 * q
    rest_budget = m * q - counts[first] * q_first
    rest_m = m - counts[first]
    q_rest = rest_budget / rest_m
    doc = {"q": q, "note": "phase-2a end-to-end test partition", "strata": []}
    for k in keys:
        fam, kind = k.split("|")
        doc["strata"].append({
            "name": k, "feature_family": int(fam), "test_kind": kind,
            "region": None, "m_s": counts[k],
            "q_s": q_first if k == first else q_rest,
            "note": "declared before the run: test fixture reallocation"})
    pth = tmp_path / "partition.json"
    pth.write_text(json.dumps(doc), encoding="utf-8")

    out1, sd1 = _run(tmp_path, "strat", _cfg(tmp_path, gpd=False, strata=str(pth)),
                     prepared)
    rep = open(os.path.join(sd1, "report.md"), encoding="utf-8").read()
    assert "STRATIFIED (weighted BH)" in rep
    assert "| stratum | m_s | q_s |" in rep
    assert "threshold_s" in rep
    assert "sum_s m_s q_s" in rep
    assert "THE HONEST TRADE" in rep
    for k in keys:
        assert k in rep
    ck1 = json.load(open(os.path.join(sd1, "checkpoint.json"), encoding="utf-8"))
    assert {r["stratum"] for r in ck1["bh"]["table"]} == set(keys)
    assert ck1["bh"]["identity"]["ok"] is True
    # every stratum row in the max-statistic table carries the global figure
    for r in ck1["max_statistic"]["per_stratum"]:
        assert r["global_p"] == ck1["max_statistic"]["global"]["p"]


def test_rule6_confirmation_runs_a_real_targeted_brute_force():
    """§P6-2(6) against the engine's own block bootstrap, not a stub.

    A GPD survivor is marked CANDIDATE-REQUIRES-BRUTE-FORCE, the targeted MC runs at
    N >= 10/p_gpd, the row is relabelled MC_RESOLVED, and the brute-force p is the
    number that survives into the stub.
    """
    (_ctx, _base, _y, window, counts, offset, marks, feats,
     _dl, _t0) = _prepared()
    f = feats[0]
    X = f.design(window, 0)
    s_obs = float(M.score_stat_all_shifts(X, counts, offset)[0])
    row = {"feature": f.name, "family": f.family, "kind": f.kind, "lag": 0,
           "test": "glm_poisson_offset_etas", "chi2_score": s_obs,
           "p_raw": 1 / 401.0, "p_bh": 0.01, "passes_fdr": True,
           "p_method": "GPD_EXTRAPOLATED", "stratum": "S"}
    cfg = {"seed": 7, "gpd": {"enabled": True, "confirm_max_n": 5000}}
    recs = ms.confirm_gpd_candidates([row], feats, window, counts, offset, marks,
                                     marks["day"] - window.start, cfg,
                                     verbose=False)
    assert row["candidate_label"] == "CANDIDATE-REQUIRES-BRUTE-FORCE"
    assert recs[0]["n_required"] == 1000            # ceil(10 / 0.01)
    assert recs[0]["status"].startswith("CONFIRMED")
    assert row["p_method"] == "MC_RESOLVED"
    assert row["p_brute_force"] == recs[0]["p_brute_force"]


def test_rule6_a_candidate_over_the_ceiling_stays_a_candidate():
    (_ctx, _base, _y, window, counts, offset, marks, feats,
     _dl, _t0) = _prepared()
    row = {"feature": feats[0].name, "family": 1, "kind": feats[0].kind, "lag": 0,
           "test": "glm_poisson_offset_etas", "chi2_score": 1.0,
           "p_raw": 1 / 401.0, "p_bh": 1e-6, "passes_fdr": True,
           "p_method": "GPD_EXTRAPOLATED", "stratum": "S"}
    cfg = {"seed": 7, "gpd": {"enabled": True, "confirm_max_n": 5000}}
    recs = ms.confirm_gpd_candidates([row], feats, window, counts, offset, marks,
                                     marks["day"] - window.start, cfg,
                                     verbose=False)
    assert recs[0]["n_required"] == 10_000_000
    assert "NOT CONFIRMED" in recs[0]["status"]
    assert row["p_method"] == "GPD_EXTRAPOLATED"    # still a candidate, no stub
    assert "p_brute_force" not in row


def test_rule7_gpd_refuses_to_run_without_a_passing_calibration(tmp_path):
    """§P6-2(7) made binding: no verdict, or a REJECT verdict, and --gpd refuses."""
    from engine import gpd_tail as G
    missing = str(tmp_path / "does_not_exist.json")
    with pytest.raises(G.NotCalibratedError) as e:
        _cfg(tmp_path, gpd=True, calibration=missing)
    assert "may not report there" in str(e.value)
    with pytest.raises(G.NotCalibratedError) as e2:
        _cfg(tmp_path, gpd=True, calibration=_calibration(tmp_path, passing=False))
    assert "REJECT" in str(e2.value)
    # the verdict is carried in the config, and therefore in the hash
    cfg = _cfg(tmp_path, gpd=True)
    assert cfg["gpd"]["calibration"]["pass"] is True
    assert cfg["gpd"]["calibration"]["coverage"] == 0.96


def test_a_partition_that_breaks_the_identity_refuses_to_run(tmp_path):
    doc = {"q": 0.10, "strata": [
        {"name": "a", "feature_family": 1, "test_kind": "glm", "m_s": 10,
         "q_s": 0.20, "note": "deliberate over-spend"},
        {"name": "b", "feature_family": 2, "test_kind": "glm", "m_s": 10,
         "q_s": 0.20, "note": "deliberate over-spend"}]}
    pth = tmp_path / "bad.json"
    pth.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(S.BudgetIdentityError):
        _cfg(tmp_path, gpd=False, strata=str(pth))
