"""The p-hacking guard: a config hash may touch the holdout exactly once."""

from __future__ import annotations

import json

import pytest

from engine import splits

CFG = {"covariate": "recent_rate", "params": {"days": 30}, "baseline": "climatology-v1"}


def test_config_hash_is_canonical():
    a = {"b": 1, "a": [1, 2], "c": {"y": 2, "x": 1}}
    b = {"c": {"x": 1, "y": 2}, "a": [1, 2], "b": 1}
    assert splits.config_hash(a) == splits.config_hash(b)
    assert splits.config_hash(a) != splits.config_hash({**a, "b": 2})


def test_holdout_refuses_second_run(tmp_path):
    hlog = str(tmp_path / "HOLDOUT_LOG.jsonl")
    elog = str(tmp_path / "EXPLORE_COUNT.jsonl")

    for _ in range(3):
        splits.log_explore_run(CFG, elog)

    h = splits.assert_holdout_unused(CFG, hlog)          # first contact: allowed
    rec = splits.record_holdout(CFG, {"bits_per_event": 0.42}, hlog, elog)
    assert rec["hash"] == h
    assert rec["n_explore_runs"] == 3, rec

    with pytest.raises(splits.HoldoutRefused) as e:      # second contact: refused
        splits.assert_holdout_unused(CFG, hlog)
    assert "HOLDOUT REFUSED" in str(e.value)
    assert "0.42" in str(e.value)

    # a genuinely different config is still allowed
    splits.assert_holdout_unused({**CFG, "params": {"days": 60}}, hlog)

    lines = [json.loads(x) for x in open(hlog, encoding="utf-8")]
    assert len(lines) == 1, "a refused run must not append to the log"


def test_explore_count_resets_after_a_holdout(tmp_path):
    hlog = str(tmp_path / "H.jsonl")
    elog = str(tmp_path / "E.jsonl")
    splits.log_explore_run(CFG, elog)
    splits.record_holdout(CFG, {}, hlog, elog)
    assert splits.n_explore_runs_since_last_holdout(elog, hlog) == 0
    splits.log_explore_run({**CFG, "params": {"days": 7}}, elog)
    assert splits.n_explore_runs_since_last_holdout(elog, hlog) == 1


def test_temporal_split_is_the_only_split():
    ex, ho = splits.temporal_split(1000, 0.7)
    assert (ex.start, ex.stop) == (0, 700)
    assert (ho.start, ho.stop) == (700, 1000)
    assert not hasattr(splits, "random_split")
    assert not hasattr(splits, "spatial_split")
    with pytest.raises(ValueError):
        splits.temporal_split(10, 1.0)
