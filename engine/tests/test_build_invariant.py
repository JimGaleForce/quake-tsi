"""§P7-8(c)(5): two sessions whose CONFIGS DIFFER must not produce IDENTICAL artifacts.

The incident this asserts against: `session_20260812T021707` (--mag-target 4.0) is
bitwise identical to `session_20260812T004857` (--mag-target 4.5), because
engine/datasets.py:CATALOG_MAG_FLOOR silently clamped the 4.0 request. An
overnight run was believed to be an independent replicate and was not.

The invariant WARNS and registers; per §P7-8(c)(1) it never edits, reduces or
deletes anything -- a collision is evidence about the build, and the declared
over-count stands because the conservative direction is the safe one.
"""

from __future__ import annotations

import json
import os

from engine import mine_session as MS

TESTS_A = [{"feature": "annual_phase", "chi2_score": 12.5, "p_raw": 0.031, "df": 2},
           {"feature": "b_value_90d", "chi2_score": 3.1, "p_raw": 0.44, "df": 1}]
TESTS_B = [{"feature": "annual_phase", "chi2_score": 99.9, "p_raw": 1e-4, "df": 2}]


def _reg(root):
    path = os.path.join(root, MS.ARTIFACT_REGISTRY)
    if not os.path.exists(path):
        return []
    with open(path, "r", encoding="utf-8") as fh:
        return [json.loads(l) for l in fh if l.strip()]


def test_content_hash_is_stable_and_order_independent():
    h = MS.artifact_content_hash(TESTS_A)
    assert h == MS.artifact_content_hash(list(TESTS_A))
    assert h != MS.artifact_content_hash(TESTS_B)
    # keys reordered inside a row must not change the hash: the hash is about the
    # numbers computed, not about dict insertion order.
    shuffled = [dict(reversed(list(r.items()))) for r in TESTS_A]
    assert MS.artifact_content_hash(shuffled) == h


def test_identical_artifacts_from_DIFFERENT_configs_are_flagged(tmp_path):
    root = str(tmp_path)
    first = MS.check_build_invariant(os.path.join(root, "session_AAA"), "cfghash1",
                                     TESTS_A, root=root)
    assert first["ok"], "the first session in an empty registry cannot collide"
    assert first["message"] is None

    second = MS.check_build_invariant(os.path.join(root, "session_BBB"), "cfghash2",
                                      TESTS_A, root=root)
    assert not second["ok"], (
        "two sessions with different config hashes and identical artifacts must "
        "violate the build invariant")
    assert [c["session"] for c in second["collisions"]] == ["session_AAA"]

    msg = second["message"]
    # the warning must NAME BOTH SESSIONS -- an anonymous warning is unactionable
    assert "session_AAA" in msg and "session_BBB" in msg
    assert "cfghash1" in msg and "cfghash2" in msg
    assert "BUILD INVARIANT VIOLATED" in msg
    assert "P7-8(c)(5)" in msg
    assert "Nothing has been deleted" in msg


def test_identical_artifacts_from_the_SAME_config_are_not_a_violation(tmp_path):
    """A reproduction is the desired behaviour, not a bug."""
    root = str(tmp_path)
    MS.check_build_invariant(os.path.join(root, "session_AAA"), "cfghash1",
                             TESTS_A, root=root)
    same = MS.check_build_invariant(os.path.join(root, "session_CCC"), "cfghash1",
                                    TESTS_A, root=root)
    assert same["ok"]
    assert same["message"] is None


def test_different_artifacts_from_different_configs_are_fine(tmp_path):
    root = str(tmp_path)
    MS.check_build_invariant(os.path.join(root, "session_AAA"), "cfghash1",
                             TESTS_A, root=root)
    other = MS.check_build_invariant(os.path.join(root, "session_DDD"), "cfghash2",
                                     TESTS_B, root=root)
    assert other["ok"]


def test_registry_is_append_only_and_nothing_is_deleted(tmp_path):
    root = str(tmp_path)
    MS.check_build_invariant(os.path.join(root, "session_AAA"), "cfghash1",
                             TESTS_A, root=root)
    MS.check_build_invariant(os.path.join(root, "session_BBB"), "cfghash2",
                             TESTS_A, root=root)
    rows = _reg(root)
    assert [r["session"] for r in rows] == ["session_AAA", "session_BBB"]
    # the colliding session is RECORDED, not withheld, and carries its collision
    assert rows[1]["collides_with"] == ["session_AAA"]
    assert rows[0]["artifact_hash"] == rows[1]["artifact_hash"]
    assert rows[0]["config_hash"] != rows[1]["config_hash"]


def test_re_registering_the_same_session_does_not_duplicate_or_self_collide(tmp_path):
    """A resumed session re-runs this check; it must not collide with itself."""
    root = str(tmp_path)
    a = MS.check_build_invariant(os.path.join(root, "session_AAA"), "cfghash1",
                                 TESTS_A, root=root)
    b = MS.check_build_invariant(os.path.join(root, "session_AAA"), "cfghash1",
                                 TESTS_A, root=root)
    assert a["ok"] and b["ok"]
    assert len(_reg(root)) == 1


def test_register_false_checks_without_writing(tmp_path):
    root = str(tmp_path)
    res = MS.check_build_invariant(os.path.join(root, "session_AAA"), "cfghash1",
                                   TESTS_A, root=root, register=False)
    assert res["ok"]
    assert _reg(root) == []


def test_the_real_incident_reproduces_as_a_violation(tmp_path):
    """The two sessions §P7-8(c) names, replayed through the check."""
    root = str(tmp_path)
    MS.check_build_invariant(os.path.join(root, "session_20260812T004857"),
                             "mag45cfg", TESTS_A, root=root)
    res = MS.check_build_invariant(os.path.join(root, "session_20260812T021707"),
                                   "mag40cfg", TESTS_A, root=root)
    assert not res["ok"]
    assert "session_20260812T004857" in res["message"]
