"""Temporal split + the holdout gate (the p-hacking guard).

SPEC rule 3: the split is TEMPORAL. Exploration = first 70% of the catalog time
span, holdout = last 30%. Spatial and random splits are not offered, on purpose:
forecast skill only counts forward in time.

SPEC rule 4: --mode holdout requires a config file; the engine hashes it and refuses
to run a hash that already appears in engine/HOLDOUT_LOG.jsonl. No flag bypasses the
refusal. Deleting the log is the human's act, not the engine's.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os

HOLDOUT_LOG = os.path.join("engine", "HOLDOUT_LOG.jsonl")
EXPLORE_COUNT = os.path.join("engine", "EXPLORE_COUNT.jsonl")


class HoldoutRefused(RuntimeError):
    pass


def temporal_split(n_days: int, explore_frac: float = 0.70):
    n_explore = int(round(explore_frac * n_days))
    if not (0 < n_explore < n_days):
        raise ValueError(f"degenerate split: {n_explore}/{n_days}")
    return slice(0, n_explore), slice(n_explore, n_days)


def config_hash(config: dict) -> str:
    canon = json.dumps(config, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canon.encode("utf-8")).hexdigest()


def _read_jsonl(path):
    if not os.path.exists(path):
        return []
    out = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


def _append_jsonl(path, obj):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(obj, sort_keys=True) + "\n")


def log_explore_run(config: dict, log_path: str = EXPLORE_COUNT):
    _append_jsonl(log_path, {
        "ts": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "hash": config_hash(config),
        "config": config,
    })


def n_explore_runs_since_last_holdout(explore_path: str = EXPLORE_COUNT,
                                      holdout_path: str = HOLDOUT_LOG) -> int:
    """Reported multiplicity: exploration runs logged since the last holdout run."""
    holds = _read_jsonl(holdout_path)
    last = holds[-1]["ts"] if holds else ""
    return sum(1 for r in _read_jsonl(explore_path) if r.get("ts", "") > last)


def assert_holdout_unused(config: dict, holdout_path: str = HOLDOUT_LOG) -> str:
    """Raise HoldoutRefused if this exact config has already touched the holdout."""
    h = config_hash(config)
    for rec in _read_jsonl(holdout_path):
        if rec.get("hash") == h:
            raise HoldoutRefused(
                "HOLDOUT REFUSED.\n"
                f"  config hash {h}\n"
                f"  was already run against the holdout at {rec.get('ts')}.\n"
                "  A holdout is spent on first contact. No flag bypasses this refusal.\n"
                "  Change the configuration, or -- knowing exactly what you are giving up --\n"
                f"  remove the line from {holdout_path} yourself. The engine will not.\n"
                f"  Previously recorded result: {json.dumps(rec.get('results'))}"
            )
    return h


def record_holdout(config: dict, results: dict, holdout_path: str = HOLDOUT_LOG,
                   explore_path: str = EXPLORE_COUNT) -> dict:
    rec = {
        "ts": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "hash": config_hash(config),
        "config": config,
        "results": results,
        "n_explore_runs": n_explore_runs_since_last_holdout(explore_path, holdout_path),
    }
    _append_jsonl(holdout_path, rec)
    return rec
