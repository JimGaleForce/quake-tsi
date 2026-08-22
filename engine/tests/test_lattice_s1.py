"""B6 / §S1: `S1-lattice-v1`, the three-layer region lattice, and its three refusals.

  * the lattice is derived from PUBLISHED classification and GEOMETRY only, and the
    module reads no catalogue in its DEFINITION (§P7-22(b)(i)) -- checked on the
    imports and the executable statements, as `test_regions_d.py` does, and for the
    same reason;
  * the L3 measurability filter is a MEASURABILITY filter and not a SELECTION filter
    (§S1.1(b)): a cell is dropped only for whether a statistic could exist;
  * the per-property-family exclusion refuses Alaska in a solid-tide scan BY NAME and
    DELEGATES to `regions_d.assert_alaska_excluded`, so there is one implementation;
  * the {L1, L2, L3} strata satisfy `sum_s m_s q_s = m q` under the module that
    already owns that identity (`engine/strata.py`), not under a private copy.
"""

from __future__ import annotations

import numpy as np
import pytest

from engine import lattice_s1 as LAT, regions_d as RD, strata as ST


# ------------------------------------------------------------- the three layers --
def test_l1_is_the_13_comcat_boxes_plus_kuril_and_cascadia():
    names = [r["name"] for r in LAT.L1_REGIONS]
    assert len(names) == 15 == LAT.N_L1
    for n in ("japan", "chile", "california", "alaska_aleutians", "greece_aegean"):
        assert n in names
    # the two D-6 boxes arrive VERBATIM from regions_d, not restated here
    for r in RD.REGIONS_D:
        got = LAT.region(r["name"])
        assert (got["lat_min"], got["lat_max"], got["lon_min"], got["lon_max"]) == \
               (r["lat_min"], r["lat_max"], r["lon_min"], r["lon_max"])
        assert got["class"] == r["class"]


def test_every_l1_box_carries_a_bird_2003_class():
    for r in LAT.L1_REGIONS:
        assert r["class"] in LAT.BIRD_CLASSES


def test_l2_is_the_union_of_the_l1_boxes_sharing_a_class():
    for l2 in LAT.L2_REGIONS:
        expect = sorted(r["name"] for r in LAT.L1_REGIONS if r["class"] == l2["class"])
        assert sorted(l2["members"]) == expect
        assert expect, l2["name"]


def test_l3_is_648_cells_of_pure_geometry():
    cells = LAT.l3_grid_cells()
    assert len(cells) == 648                     # (180/10) * (360/10)
    assert all(c["class"] is None for c in cells)
    lat = np.array([c["lat_min"] for c in cells])
    lon = np.array([c["lon_min"] for c in cells])
    assert lat.min() == -90.0 and lat.max() == 80.0
    assert lon.min() == -180.0 and lon.max() == 170.0


def test_l3_cells_tile_without_double_counting():
    """A closed grid double-counts every shared edge; `contains` is half-open above."""
    cells = LAT.l3_grid_cells()
    pts = [(0.0, 0.0), (10.0, -180.0), (-45.3, 122.7), (79.9, 169.9)]
    for la, lo in pts:
        hits = [c for c in cells if LAT.contains(c, la, lo)]
        assert len(hits) == 1, (la, lo, [h["name"] for h in hits])


# ------------------------------------------------- the measurability filter -----
def test_l3_measurability_keeps_only_cells_that_could_carry_a_statistic():
    rng = np.random.default_rng(0)
    # 500 events packed into one cell, 5 scattered elsewhere
    lat = np.concatenate([rng.uniform(30.0, 40.0, 500), rng.uniform(-80, -70, 5)])
    lon = np.concatenate([rng.uniform(130.0, 140.0, 500), rng.uniform(0, 10, 5)])
    kept, audit = LAT.l3_measurable_cells(lat, lon, "M>=4.5")
    assert [c["name"] for c in kept] == ["g+030+0130"]
    assert kept[0]["n_explore_events"] == 500
    assert audit["n_cells_declared"] == 648
    assert audit["n_cells_occupied"] == 2
    assert "EXPLORATION ONLY" in audit["window"]


def test_l3_filter_is_measurability_not_selection():
    """§S1.1(b): the filter sees COUNTS ONLY -- it is handed no statistic at all."""
    import ast
    import inspect
    fn = ast.parse(inspect.getsource(LAT.l3_measurable_cells).lstrip()).body[0]
    body = [n for n in fn.body
            if not (isinstance(n, ast.Expr) and isinstance(n.value, ast.Constant))]
    code = chr(10).join(ast.unparse(n) for n in body).lower()
    assert "count_nonzero" in code
    # the EXECUTABLE body, not the prose: the docstring names a statistic precisely
    # to say the filter never sees one, so a text search over the whole source would
    # fail on the disclaimer -- the same trap `test_regions_d.py` records.
    for forbidden in ("p_value", "pval", "statistic", "kuiper", "concentration"):
        assert forbidden not in code, forbidden
    assert "MEASURABILITY filter and not a SELECTION filter" in \
           LAT.L3_MEASURABILITY_NOTE


# ---------------------------------------------------- catalogue-free definition --
def test_the_module_reads_no_catalogue_in_its_definition():
    """§P7-22(b)(i)/Rule 4.1, checked as `test_regions_d.py` checks it."""
    import ast
    import inspect
    tree = ast.parse(inspect.getsource(LAT))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported.add(node.module or "")
            imported.update(a.name for a in node.names)
    assert imported <= {"", "regions_d", "RD", "splits", "annotations", "__future__",
                        "json", "numpy", "np"}, imported
    code = "\n".join(ast.unparse(n) for n in ast.walk(tree)
                     if isinstance(n, (ast.Call, ast.Attribute, ast.Name)))
    for forbidden in ("datasets", "load_catalog", "read_csv", "loadtxt", "k080",
                      "census", "holdout"):
        assert forbidden not in code, forbidden


# ------------------------------------------------------------- the exclusions ---
def test_alaska_is_refused_in_a_solid_tide_scan_by_name():
    with pytest.raises(RD.SeedingRegionInConfirmationSet):
        LAT.assert_family_exclusions("solid_tide", ["japan", "alaska_aleutians"])
    with pytest.raises(RD.SeedingRegionInConfirmationSet):
        LAT.assert_family_exclusions("solid_tide", ["Alaska-Aleutians"])


def test_alaska_may_appear_for_a_family_it_never_seeded():
    """§S1.1(c): the exclusion is PER PROPERTY FAMILY, not global."""
    got = LAT.assert_family_exclusions("human_clock", ["japan", "alaska_aleutians"])
    assert "alaska_aleutians" in got


def test_the_solid_tide_refusal_delegates_to_regions_d():
    """One implementation of the by-name check, not two that can drift apart."""
    import inspect
    src = inspect.getsource(LAT.assert_family_exclusions)
    assert "RD.assert_alaska_excluded" in src


def test_filter_regions_drops_the_seeded_region_and_says_so():
    regs = [LAT.region(n) for n in ("japan", "alaska_aleutians", "chile")]
    kept, rep = LAT.filter_regions_for_family("solid_tide", regs)
    assert [r["name"] for r in kept] == ["japan", "chile"]
    assert rep["dropped"] == ["alaska_aleutians"]


def test_the_seeded_region_list_contains_the_k092_seed():
    assert "alaska_aleutians" in LAT.excluded_names_for_family("solid_tide")
    assert LAT.excluded_names_for_family("ephemeris") == ()


# ------------------------------------------------------------- the redundancy ---
def test_l2_finding_supplied_by_one_l1_box_is_labelled_redundant():
    l2 = ["e%d" % i for i in range(10)]
    by_l1 = {"japan": l2[:9], "chile": l2[9:]}      # 90% from japan
    r = LAT.redundancy_label(l2, by_l1)
    assert r["redundant"] and r["dominant_l1"] == "japan"
    assert r["label"] == LAT.REDUNDANT_WITH_L1
    assert r["fraction"] == pytest.approx(0.9)


def test_l2_finding_spread_over_boxes_is_not_redundant():
    l2 = ["e%d" % i for i in range(10)]
    by_l1 = {"japan": l2[:5], "chile": l2[5:]}
    r = LAT.redundancy_label(l2, by_l1)
    assert not r["redundant"] and r["label"] is None


# -------------------------------------------------------------- the declaration --
def test_the_strata_satisfy_the_budget_identity_under_strata_py():
    """§S1.1(a): the layer axis is a STRATUM, and the identity is asserted by the
    module that already owns it -- never by a private copy."""
    cells = [{"name": "g+030+0130"}] * 7
    table = LAT.region_strata(cells, q=0.10, n_properties=5, n_mag_strata=4)
    m = sum(t["m_s"] for t in table)
    ident = ST.assert_budget_identity(table, q=0.10)
    assert ident["m"] == m
    assert m == (LAT.N_L1 + LAT.N_L2 + 7) * 5 * 4


def test_the_declaration_embeds_regions_d_verbatim_and_hashes():
    b = LAT.declaration_block()
    assert b["regions_d"] == __import__(
        "engine.regions_d", fromlist=["x"]).declaration_block()
    assert b["lattice_rule_id"] == "S1-lattice-v1"
    d1 = LAT.declaration_digest()
    d2 = LAT.declaration_digest([{"name": "g+030+0130"}])
    assert d1 != d2, "the frozen L3 cell list must be hash-affecting (§S1.1(b))"
    assert len(d1) == 12


def test_region_lookup_is_name_normalised_and_raises_on_unknown():
    assert LAT.region("Alaska-Aleutians")["name"] == "alaska_aleutians"
    assert LAT.region("class_sub")["layer"] == "L2"
    with pytest.raises(KeyError):
        LAT.region("atlantis")


def test_l2_region_mask_is_the_union_of_its_members():
    lat = np.array([35.0, -30.0, 64.0])
    lon = np.array([140.0, -70.0, -19.0])
    sub = LAT.region_mask(LAT.region("class_sub"), lat, lon)
    crb = LAT.region_mask(LAT.region("class_crb"), lat, lon)
    assert sub.tolist() == [True, True, False]      # japan, chile
    assert crb.tolist() == [False, False, True]     # iceland
