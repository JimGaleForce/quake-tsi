"""D-6: Kuril and Cascadia as declared constants, and the two guards §P7-22(b) demands.

  * the partition is derived from PUBLISHED tectonic classification, not catalogue data
    (§P7-22(b)(i)) -- checked as the absence of any catalogue import and the presence of
    the citations;
  * `R2b-lon6-active` is FORBIDDEN on this question (§P7-22(b)(iv)), because lon6
    sector 0 carries the seeding region into the confirmation set;
  * Alaska is excluded BY NAME and never by trusting a boundary (§P7-22 Rat. 5).
"""

from __future__ import annotations

import pytest

from engine import regions_d as RD, splits


# --------------------------------------------------------------- the boxes ----
def test_kuril_covers_the_gap_above_the_japan_box():
    """§K92-0: 'the Japan box stops at 46 N'. Kuril must start below it and go north."""
    k = RD.region("kuril")
    assert k["lat_min"] <= 46.0 <= k["lat_max"]
    assert k["lat_max"] >= 55.0                     # into Kamchatka
    assert RD.contains(k, 46.9, 153.3, 40.0)        # central Kuril arc
    assert RD.contains(k, 52.9, 159.7, 60.0)        # Kamchatka
    assert not RD.contains(k, 35.7, 139.7, 30.0)    # Tokyo is not Kuril


def test_cascadia_covers_the_gap_above_the_california_box():
    """§K92-0: 'the California box stops at 42 N'. Cascadia spans Mendocino->Nootka."""
    c = RD.region("cascadia")
    assert c["lat_min"] <= 42.0 <= c["lat_max"]
    assert c["lat_min"] <= 40.3                     # Mendocino triple junction
    assert c["lat_max"] >= 50.0                     # Nootka fault
    assert RD.contains(c, 47.6, -122.3, 50.0)       # Puget Sound
    assert RD.contains(c, 44.0, -125.0, 20.0)       # offshore Oregon
    assert not RD.contains(c, 34.0, -118.2, 10.0)   # Los Angeles is not Cascadia


def test_both_are_classified_SUB_with_a_published_source():
    for name in ("kuril", "cascadia"):
        r = RD.region(name)
        assert r["class"] == "SUB"
        assert "Bird (2003)" in r["sources"]
        assert r["held_out_role"].startswith("CONFIRMATION")


def test_the_boxes_carry_the_caveat_that_a_rectangle_is_not_a_plate_boundary():
    assert "not a tectonic claim" in RD.region("kuril")["box_caveat"]


def test_the_module_reads_no_catalogue():
    """§P7-22(b)(i)/§P6-4 Rule 4.1: nothing in the derivation touches data.

    Checked on the module's IMPORTS and its executable statements, not on its prose:
    the docstring names the K-080 census precisely to say it is forbidden as a
    selector, so a text search over the whole file would fail on the disclaimer.
    """
    import ast
    import inspect
    tree = ast.parse(inspect.getsource(RD))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported.add(node.module or "")
            imported.update(a.name for a in node.names)
    # `splits` supplies config_hash only; nothing that can reach a catalogue.
    assert imported <= {"", "splits", "annotations", "__future__", "json"}, imported
    # and no executable statement mentions a data source
    code = "\n".join(
        ast.unparse(n) for n in ast.walk(tree)
        if isinstance(n, (ast.Call, ast.Attribute, ast.Name)))
    for forbidden in ("datasets", "load_event_marks", "k080", "census",
                      "checkpoint", "loadtxt", "read_csv", "open"):
        assert forbidden not in code, forbidden


# ---------------------------------------------------------------- the guards --
def test_lon6_is_forbidden_on_this_question():
    """§P7-22(b)(iv), and for the substantive reason: sector 0 contains Alaska."""
    with pytest.raises(RD.ForbiddenPartitionRule):
        RD.assert_partition_rule("R2b-lon6-active")
    assert RD.assert_partition_rule(RD.REGION_RULE_ID_D) == RD.REGION_RULE_ID_D
    with pytest.raises(RD.ForbiddenPartitionRule):
        RD.assert_partition_rule("some-other-rule")
    assert "SEEDING REGION" in RD.LON6_FORBIDDEN_ON_THIS_QUESTION


def test_alaska_is_excluded_by_name_in_every_spelling_we_use():
    for spelling in ("alaska", "Alaska", "Alaska-Aleutians", "alaska_aleutians",
                     "ALEUTIANS"):
        with pytest.raises(RD.SeedingRegionInConfirmationSet):
            RD.assert_alaska_excluded(["japan", spelling, "chile"])
    ok = RD.assert_alaska_excluded(["japan", "kuril", "chile", "cascadia"])
    assert "kuril" in ok and "cascadia" in ok


def test_the_exclusion_is_by_name_and_not_by_boundary():
    """A boundary can be got wrong; a name cannot. The check must be on the NAME."""
    import inspect
    src = inspect.getsource(RD.assert_alaska_excluded)
    assert "EXCLUDED_BY_NAME" in src
    assert "lat_min" not in src and "lon_min" not in src


# --------------------------------------------------------------- the record ---
def test_the_declaration_is_hash_affecting():
    cfg = {"partition": RD.declaration_block()}
    h0 = splits.config_hash(cfg)
    blk = RD.declaration_block()
    blk["regions"][0]["lat_max"] = 57.0
    assert splits.config_hash({"partition": blk}) != h0
    assert len(RD.declaration_digest()) == 12
