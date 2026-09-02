"""Tests for prospective_socal.py (EQ-25 Phase 4). NO NETWORK.

Every test monkeypatches prospective_socal.download_catalog to a synthetic SoCal
catalogue, and redirects LOG and CACHE into tmp_path. The frozen parameters are read
from the real results_exp_h.json, which is what the module under test is supposed to
do; nothing is refit anywhere.
"""

import hashlib
import json
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import prospective_socal as ps  # noqa: E402


T0 = pd.Timestamp("2026-06-01T00:00:00Z")


def synth_catalogue(t_end, per_day=4.0, days_back=420.0, week_counts=None, t_week=None):
    """Uniform synthetic catalogue up to t_end, optionally plus a forecast week.

    week_counts, if given, are seven daily counts placed inside [t_week, t_week+7d).
    """
    rows = []
    n = int(round(per_day * days_back))
    step = days_back / max(n, 1)
    for i in range(n):
        t = t_end - pd.Timedelta(days=days_back) + pd.Timedelta(days=step * i)
        rows.append((t, 34.0, -118.0, 2.6 + (i % 7) * 0.1, f"syn{i}", "earthquake"))
    if week_counts is not None:
        for d, k in enumerate(week_counts):
            for j in range(k):
                t = t_week + pd.Timedelta(days=d) + pd.Timedelta(hours=1 + j % 20)
                rows.append((t, 34.0, -118.0, 2.7, f"wk{d}_{j}", "earthquake"))
    df = pd.DataFrame(rows, columns=["time", "latitude", "longitude", "mag",
                                     "id", "type"])
    df = df.sort_values("time").reset_index(drop=True)
    return df


@pytest.fixture
def env(tmp_path, monkeypatch):
    """Redirect the log and the cache; stub the fetch with a synthetic catalogue."""
    log = tmp_path / "results_prospective_socal_log.json"
    cache = tmp_path / "comcat_socal_m25.csv"
    cache.write_text("time,latitude,longitude,mag,id,type\n", encoding="utf-8")
    monkeypatch.setattr(ps, "LOG", str(log))
    monkeypatch.setattr(ps, "CACHE", str(cache))

    state = {"cat": synth_catalogue(T0), "fail": False}

    def fake_download(now, refresh=False):
        if state["fail"]:
            raise SystemExit("ABORT: synthetic ComCat failure")
        return state["cat"].copy(), {"windows": 0, "cache_file": str(cache)}

    monkeypatch.setattr(ps, "download_catalog", fake_download)
    return {"log": log, "cache": cache, "state": state, "tmp": tmp_path}


def emit(now=T0, nsims=20):
    return ps.main(["emit", "--now", now.isoformat(), "--nsims", str(nsims)])


def load(env):
    return json.loads(env["log"].read_text(encoding="utf-8"))


# ------------------------------------------------------------------ hash tests
def test_commitment_hash_round_trips_and_detects_one_byte_change(env):
    assert emit() == 0
    doc = load(env)
    rec = doc["records"][0]
    assert ps.canonical_hash(rec["commitment"]) == rec["commitment_hash"]

    # one-byte change anywhere in the commitment breaks the hash
    tampered = json.loads(json.dumps(rec["commitment"]))
    tampered["T_utc"] = tampered["T_utc"].replace("2026", "2027", 1)
    assert ps.canonical_hash(tampered) != rec["commitment_hash"]

    # and a one-digit change in a committed number breaks it too
    tampered2 = json.loads(json.dumps(rec["commitment"]))
    tampered2["etas_expected_counts_per_day"][0] += 1e-9
    assert ps.canonical_hash(tampered2) != rec["commitment_hash"]


def test_tampered_commitment_is_refused_at_scoring(env):
    assert emit() == 0
    doc = load(env)
    doc["records"][0]["commitment"]["etas_expected_counts_per_day"][3] = 999.0
    env["log"].write_text(json.dumps(doc, indent=1), encoding="utf-8")
    later = T0 + pd.Timedelta(days=11)
    rc = ps.main(["score", "--record", doc["records"][0]["record_id"],
                  "--now", later.isoformat()])
    assert rc == 2


# -------------------------------------------------------------------- due rule
def test_not_due_before_T_plus_10d_and_due_after(env):
    assert emit() == 0
    doc = load(env)
    rid = doc["records"][0]["record_id"]
    c = doc["records"][0]["commitment"]

    before = T0 + pd.Timedelta(days=9, hours=23)
    after = T0 + pd.Timedelta(days=10, seconds=1)
    assert ps.is_due(c, before)[0] is False
    assert ps.is_due(c, after)[0] is True

    # and the CLI refuses with exit code 2 before it is due
    assert ps.main(["score", "--record", rid, "--now", before.isoformat()]) == 2
    before_bytes = env["log"].read_bytes()

    env["state"]["cat"] = synth_catalogue(after)
    assert ps.main(["score", "--record", rid, "--now", after.isoformat()]) == 0
    assert env["log"].read_bytes() != before_bytes
    assert load(env)["records"][0]["scored"] is True


def test_all_due_scores_nothing_when_nothing_is_due(env):
    assert emit() == 0
    before_bytes = env["log"].read_bytes()
    early = T0 + pd.Timedelta(days=5)
    assert ps.main(["score", "--all-due", "--now", early.isoformat()]) == 0
    assert env["log"].read_bytes() == before_bytes


# -------------------------------------------------------------------- scoring
def _hand_record(env, etas_daily, rate_a, rate_b, t0=T0):
    """Append a hand-built, correctly hashed commitment so the scoring maths is
    exercised against numbers chosen by the test rather than by the simulator."""
    commitment = {
        "record_id": "W-TEST",
        "T_utc": t0.isoformat(),
        "horizon_days": 7.0,
        "etas_expected_counts_per_day": list(etas_daily),
        "baselines": {
            ps.BASELINE_A: {"rate_per_day": rate_a, "role": "PRIMARY"},
            ps.BASELINE_B: {"rate_per_day": rate_b, "role": "SECONDARY"},
        },
    }
    doc = ps.empty_log(t0)
    doc["records"].append({
        "record_id": "W-TEST",
        "commitment": commitment,
        "commitment_hash": ps.canonical_hash(commitment),
        "scored": False,
        "scoring": None,
    })
    env["log"].write_text(json.dumps(doc, indent=1), encoding="utf-8")
    return commitment


def test_perfect_etas_beats_poisson_baseline(env):
    truth = [20, 3, 2, 2, 1, 1, 1]           # a strongly clustered week
    flat = sum(truth) / 7.0
    _hand_record(env, [float(x) for x in truth], rate_a=flat, rate_b=3.343)

    after = T0 + pd.Timedelta(days=11)
    env["state"]["cat"] = synth_catalogue(after, per_day=0.0, days_back=1.0,
                                          week_counts=truth, t_week=T0)
    assert ps.main(["score", "--record", "W-TEST", "--now", after.isoformat()]) == 0

    s = load(env)["records"][0]["scoring"]
    assert s["observed_counts_per_day"] == truth
    assert s["n_events"] == sum(truth)
    assert s["ll_etas"] > s["ll_A"]
    assert s["ll_etas"] > s["ll_B"]
    assert s["bits_per_event_vs_A"] > 0.0
    assert s["bits_per_event_vs_B"] > 0.0
    # bits/event is exactly (ll_model - ll_baseline)/(n ln2)
    assert s["bits_per_event_vs_A"] == pytest.approx(
        (s["ll_etas"] - s["ll_A"]) / (s["n_events"] * ps.LN2))


def test_zero_events_does_not_divide_by_zero(env):
    _hand_record(env, [1.0] * 7, rate_a=2.0, rate_b=3.343)
    after = T0 + pd.Timedelta(days=11)
    # catalogue with events only strictly BEFORE T, none in the committed week
    cat = synth_catalogue(T0 - pd.Timedelta(seconds=1), per_day=4.0, days_back=400.0)
    env["state"]["cat"] = cat
    assert ps.main(["score", "--record", "W-TEST", "--now", after.isoformat()]) == 0

    s = load(env)["records"][0]["scoring"]
    assert s["observed_counts_per_day"] == [0] * 7
    assert s["n_events"] == 0
    assert s["bits_per_event_vs_A"] is None
    assert s["bits_per_event_vs_B"] is None
    assert s["delta_ll_nats_vs_A"] == pytest.approx(7 * (2.0 - 1.0))
    assert s["zero_event_note"]
    su = load(env)["summary"]
    assert su["weeks_scored"] == 1
    assert su["total_events"] == 0
    assert su["cumulative_bits_per_event_vs_" + ps.BASELINE_A] is None


# ------------------------------------------------------------- failure-first
def test_refused_emit_leaves_log_byte_identical(env):
    assert emit() == 0
    before_bytes = env["log"].read_bytes()

    # a catalogue that is quiet in the trailing 30 days: below MIN_TRAIL30
    quiet = synth_catalogue(T0 - pd.Timedelta(days=31), per_day=4.0, days_back=400.0)
    env["state"]["cat"] = quiet
    later = T0 + pd.Timedelta(days=1)
    assert ps.main(["emit", "--now", later.isoformat(), "--nsims", "5"]) == 2
    assert env["log"].read_bytes() == before_bytes

    # and a storm above MAX_TRAIL30 is refused the same way
    storm = synth_catalogue(T0, per_day=400.0, days_back=40.0)
    env["state"]["cat"] = storm
    assert ps.main(["emit", "--now", later.isoformat(), "--nsims", "5"]) == 2
    assert env["log"].read_bytes() == before_bytes


def test_fetch_failure_writes_nothing_and_exits_nonzero(env):
    assert emit() == 0
    before_bytes = env["log"].read_bytes()
    env["state"]["fail"] = True
    later = T0 + pd.Timedelta(days=1)
    rc = ps.main(["emit", "--now", later.isoformat(), "--nsims", "5"])
    assert rc != 0
    assert env["log"].read_bytes() == before_bytes


def test_log_gains_exactly_one_record_per_emit(env):
    for i in range(3):
        t = T0 + pd.Timedelta(days=7 * i)
        env["state"]["cat"] = synth_catalogue(t)
        assert emit(now=t) == 0
        assert len(load(env)["records"]) == i + 1
    ids = [r["record_id"] for r in load(env)["records"]]
    assert len(set(ids)) == 3


def test_duplicate_record_id_is_refused(env):
    assert emit() == 0
    before_bytes = env["log"].read_bytes()
    assert emit() == 2
    assert env["log"].read_bytes() == before_bytes


def test_emitter_drift_warns_but_continues(env, capsys):
    assert emit() == 0
    doc = load(env)
    doc["protocol"]["provenance"]["emitter_sha256"] = "0" * 64
    env["log"].write_text(json.dumps(doc, indent=1), encoding="utf-8")
    t = T0 + pd.Timedelta(days=7)
    env["state"]["cat"] = synth_catalogue(t)
    assert emit(now=t) == 0
    out = capsys.readouterr().out
    assert "WARNING: prospective_socal.py has CHANGED" in out
    assert len(load(env)["records"]) == 2


# ------------------------------------------------------------------ invariants
def test_emitted_record_shape_and_no_refit(env):
    assert emit() == 0
    rec = load(env)["records"][0]
    c = rec["commitment"]
    assert c["refit"] is False
    assert len(c["etas_expected_counts_per_day"]) == 7
    assert len(c["analytic_floor_counts_per_day_no_future_triggering"]) == 7
    assert all(x >= 0 for x in c["etas_expected_counts_per_day"])
    assert set(c["baselines"]) == {ps.BASELINE_A, ps.BASELINE_B}
    assert c["baselines"][ps.BASELINE_A]["role"] == "PRIMARY"
    assert c["input_sha256"]["results_exp_h.json"] == ps.sha256_file(ps.EXP_H)
    assert c["input_sha256"]["prospective_socal.py"] == ps.sha256_file(ps.SELF)
    assert rec["scored"] is False
    assert rec["scoring"] is None
    # the ETAS mean must sit at or above its own analytic no-future-triggering floor
    assert sum(c["etas_expected_counts_per_day"]) >= 0.9 * sum(
        c["analytic_floor_counts_per_day_no_future_triggering"])


def test_poisson_ll_matches_scipy_free_formula():
    import math
    assert ps.poisson_ll([0], [2.0]) == pytest.approx(-2.0)
    assert ps.poisson_ll([3], [2.0]) == pytest.approx(
        3 * math.log(2.0) - 2.0 - math.log(6))
    # zero rate is floored, not a crash
    assert math.isfinite(ps.poisson_ll([0, 1], [0.0, 1.0]))


def test_canonical_json_is_sorted_and_whitespace_free():
    b = ps.canonical({"b": 1, "a": {"d": 2, "c": 3}})
    assert b == b'{"a":{"c":3,"d":2},"b":1}'
    assert hashlib.sha256(b).hexdigest() == ps.canonical_hash({"b": 1,
                                                               "a": {"c": 3, "d": 2}})
