"""B6 / §S1 -- `S1-lattice-v1`: THE DECLARED REGION LATTICE, three layers, frozen.

SEARCHER.md §S1 asks for a region lattice with a rule id, derived from published
classification and geometry only, frozen before any scan. §P7-24 governs where the two
differ. This module is that declaration and nothing else: **it reads no catalogue in
its definition**, exactly as `engine/regions_d.py` reads none in its own, and for the
same reason -- §P7-22(b)(i), *"defined from published tectonic classification, not from
catalogue data"*, is what removes the forking hazard entirely.

It EXTENDS `regions_d.declaration_block()` by composition, never by mutation: the
D-6 rule id `D6-tectonic-bird2003-v1` and its two boxes are embedded verbatim inside
`declaration_block()` here, so the Tranche D declaration and the SEARCHER declaration
can never drift apart, and nothing in `regions_d` changes behaviour by this module
existing.

THE THREE LAYERS (§S1.1), and the exact price of their overlap
---------------------------------------------------------------
* **L1 -- named arcs/provinces (15).** The 13 ComCat boxes on disk (`data/comcat_world/
  *.csv`, whose lat/lon extents are the ones `exp_m_world_transfer.py:REGIONS` states
  verbatim from `PATTERN_PROTOCOL.md`) plus KURIL and CASCADIA from `regions_d`. Each
  carries a Bird (2003) plate-boundary class.
* **L2 -- tectonic classes (4).** The union of the L1 boxes sharing a Bird class.
  Overlaps L1 by construction, and that is the point: §P7-23's lesson is that a
  regional claim and a class claim are different hypotheses.
* **L3 -- declared global grid.** 10 deg x 10 deg cells. **The grid is declared here in
  full and unconditionally** (`l3_grid_cells()` enumerates all 648 cells from geometry
  alone); the MEASURABILITY filter `N >= 200 at the stratum's Mc` is applied by
  `l3_measurable_cells()` **on exploration-window events only** and its surviving list
  is frozen and hash-affecting. §S1.1(b) is explicit that this is the honest form: no
  cell is dropped for what its statistic looked like, only for whether a statistic
  could exist at all.

**Overlap is ALLOWED and is priced as a STRATUM, not as a flat multiplicity factor**
(§S1.1(a)): strata = {L1, L2, L3}, budget identity `sum_s m_s q_s = m q` enforced by
`engine/strata.py:assert_budget_identity`, which already exists and is used rather than
re-implemented. An L2 finding whose covered events are >= 80% supplied by a single L1
box is labelled `REDUNDANT-WITH-L1` by `redundancy_label` and is one finding, not two.

THE PER-PROPERTY-FAMILY EXCLUSION LIST, which only ever grows
--------------------------------------------------------------
§S1.1(c): Alaska-Aleutians is excluded BY NAME from every tidal-phase-family scan
(§P7-22 Ratification 5). It may appear as an L1 region for property families unrelated
to K-092, so the exclusion is **per property family**, and `assert_family_exclusions`
is the check at scan entry. `regions_d.assert_alaska_excluded` is CALLED rather than
copied for the tidal families, so there is exactly one implementation of the by-name
refusal in the engine and a future edit cannot fix one and miss the other.

Nothing in this module is evidence. It declares where one may look.
"""

from __future__ import annotations

from . import regions_d as RD, splits

# --------------------------------------------------------------- the rule id --
LATTICE_RULE_ID = "S1-lattice-v1"

LATTICE_RULE = (
    "S1-lattice-v1: a three-layer region lattice. L1 = 15 named arc/province boxes "
    "(13 ComCat boxes on disk + KURIL + CASCADIA from D6-tectonic-bird2003-v1), each "
    "carrying Bird (2003)'s plate-boundary class. L2 = the 4 Bird classes, each the "
    "union of the L1 boxes sharing it. L3 = the global 10x10 degree grid, declared "
    "from geometry alone and filtered for MEASURABILITY (N >= 200 at the stratum's "
    "Mc) on the EXPLORATION WINDOW ONLY. Overlap between layers is allowed and is "
    "priced as a declared STRATUM (§P6-3), never as a flat multiplicity factor.")

OVERLAP_PRICE = (
    "SEARCHER.md §S1.1(a): overlapping regions are CORRELATED tests. Naive Bonferroni "
    "over R = 15 + 4 + n_L3 over-corrects (they are not independent); naive BH "
    "under-controls (BH assumes independence or PRDS, and nested-set statistics are "
    "neither in general). Disposition: the region layer is a declared STRATUM under "
    "§P6-3 with the budget identity sum_s m_s q_s = m q enforced by "
    "engine/strata.py:assert_budget_identity.")

L3_MEASURABILITY_NOTE = (
    "SEARCHER.md §S1.1(b): the L3 count threshold touches catalogue data, so it is "
    "declared as a MEASURABILITY filter and not a SELECTION filter -- no cell is "
    "dropped for what its statistic looked like, only for whether a statistic could "
    "exist. It is applied on the EXPLORATION WINDOW ONLY (Rule 4.1: nothing from the "
    "holdout enters), it is declared before the scan, and the surviving cell list is "
    "frozen and hash-affecting.")

REDUNDANT_WITH_L1 = "REDUNDANT-WITH-L1"
REDUNDANCY_FRACTION = 0.80
REDUNDANCY_RULE = (
    "SEARCHER.md §S1.1(a): an L2 (class) finding whose covered events are >= %.0f%% "
    "supplied by ONE L1 box is labelled %s and is reported as ONE finding, not two. "
    "The L1 hit and the L2 hit are the same events seen twice."
    % (100.0 * REDUNDANCY_FRACTION, REDUNDANT_WITH_L1))

# --------------------------------------------------------------------- L1 -----
# lat_min, lat_max, lon_min, lon_max in degrees, lon in [-180, 180].
#
# The 13 box extents are VERBATIM from `exp_m_world_transfer.py:REGIONS`, which states
# them verbatim from PATTERN_PROTOCOL.md and which is the declaration the CSVs on disk
# were downloaded against. They are repeated here rather than imported because
# `exp_m_world_transfer` is an experiment script that opens catalogues at import-adjacent
# scope, and §P7-22(b)(i) wants this module provably catalogue-free.
#
# `bird_class` maps the protocol's tectonic word to Bird (2003)'s boundary class:
#   subduction -> SUB, transform -> OTF, collision -> CCB, rift -> CRB.
_COMCAT_BOXES = (
    # name,               protocol type,  bird, lat_min, lat_max, lon_min, lon_max
    ("japan",             "subduction", "SUB",   30.0,  46.0,  129.0,  147.0),
    ("chile",             "subduction", "SUB",  -46.0, -17.0,  -76.0,  -66.0),
    ("indonesia",         "subduction", "SUB",  -11.0,   6.0,   95.0,  130.0),
    ("california",        "transform",  "OTF",   31.5,  42.0, -125.0, -113.0),
    ("turkey",            "transform",  "OTF",   35.0,  42.0,   25.0,   45.0),
    ("himalaya",          "collision",  "CCB",   25.0,  38.0,   70.0,   98.0),
    ("iceland",           "rift",       "CRB",   62.0,  67.0,  -25.0,  -13.0),
    ("alaska_aleutians",  "subduction", "SUB",   50.0,  62.0, -180.0, -140.0),
    ("mexico",            "subduction", "SUB",   14.0,  20.0, -105.0,  -92.0),
    ("philippines",       "subduction", "SUB",    5.0,  20.0,  120.0,  128.0),
    ("caribbean",         "transform",  "OTF",   17.0,  20.0,  -75.0,  -68.0),
    ("iran",              "collision",  "CCB",   26.0,  36.0,   44.0,   62.0),
    ("greece_aegean",     "rift",       "CRB",   34.0,  41.0,   19.0,   29.0),
)

BIRD_CLASS_SOURCE = (
    "Bird, P. (2003), 'An updated digital model of plate boundaries', Geochem. "
    "Geophys. Geosyst. 4(3), 1027; Bird & Kagan (2004), BSSA 94(6), 2380. Classes: "
    "SUB subduction, OTF/CTF transform, CRB/OSR rift, CCB continental convergent.")

BOX_SOURCE_COMCAT = (
    "PATTERN_PROTOCOL.md region boxes, verbatim, as stated in "
    "exp_m_world_transfer.py:REGIONS -- the same extents the data/comcat_world CSVs "
    "were downloaded against. A rectangle is a DATA-SELECTION ENVELOPE, not a "
    "tectonic claim (regions_d.BOX_CAVEAT).")


def _l1_from_comcat():
    out = []
    for name, typ, cls, la0, la1, lo0, lo1 in _COMCAT_BOXES:
        out.append({
            "name": name, "layer": "L1", "class": cls, "protocol_type": typ,
            "lat_min": la0, "lat_max": la1, "lon_min": lo0, "lon_max": lo1,
            "depth_max_km": None,
            "source": BOX_SOURCE_COMCAT,
        })
    return out


def _l1_from_regions_d():
    out = []
    for r in RD.REGIONS_D:                      # KURIL, CASCADIA -- verbatim
        out.append({
            "name": r["name"], "layer": "L1", "class": r["class"],
            "protocol_type": "subduction",
            "lat_min": r["lat_min"], "lat_max": r["lat_max"],
            "lon_min": r["lon_min"], "lon_max": r["lon_max"],
            "depth_max_km": r["depth_max_km"],
            "source": "engine/regions_d.py (%s), verbatim" % RD.REGION_RULE_ID_D,
        })
    return out


L1_REGIONS = tuple(_l1_from_comcat() + _l1_from_regions_d())
N_L1 = len(L1_REGIONS)

# --------------------------------------------------------------------- L2 -----
BIRD_CLASSES = ("SUB", "OTF", "CCB", "CRB")


def l2_regions():
    """The 4 Bird classes, each as the declared union of the L1 boxes sharing it."""
    out = []
    for cls in BIRD_CLASSES:
        members = [r["name"] for r in L1_REGIONS if r["class"] == cls]
        out.append({
            "name": "class_%s" % cls.lower(), "layer": "L2", "class": cls,
            "members": members, "source": BIRD_CLASS_SOURCE,
            "overlap_note": OVERLAP_PRICE,
        })
    return out


L2_REGIONS = tuple(l2_regions())
N_L2 = len(L2_REGIONS)

# --------------------------------------------------------------------- L3 -----
L3_CELL_DEG = 10.0
L3_MIN_EVENTS = 200


def l3_grid_cells(cell_deg: float = L3_CELL_DEG):
    """EVERY 10x10 degree cell of the globe, from geometry alone. No catalogue.

    648 cells at 10 deg. This is the DECLARATION; `l3_measurable_cells` is the
    measurability filter, and the two are separate functions precisely so that the
    catalogue-touching half is impossible to confuse with the catalogue-free half.
    """
    d = float(cell_deg)
    n_lat = int(round(180.0 / d))
    n_lon = int(round(360.0 / d))
    out = []
    for i in range(n_lat):
        la0 = -90.0 + i * d
        for j in range(n_lon):
            lo0 = -180.0 + j * d
            out.append({
                "name": "g%+04d%+05d" % (int(la0), int(lo0)),
                "layer": "L3", "class": None,
                "lat_min": la0, "lat_max": la0 + d,
                "lon_min": lo0, "lon_max": lo0 + d,
                "depth_max_km": None,
                "source": "geometry: global %g deg grid" % d,
            })
    return out


def l3_measurable_cells(lat, lon, mc_label, cell_deg: float = L3_CELL_DEG,
                        min_events: int = L3_MIN_EVENTS):
    """The L3 cells that COULD carry a statistic, from exploration-window events only.

    `lat`/`lon` must already be the EXPLORATION-WINDOW events at the stratum's Mc --
    this function does not slice a window and does not know where the holdout is, so
    it cannot spend one. Returns (cells, audit); the audit carries the counts and the
    §S1.1(b) note so the frozen list is self-describing.
    """
    import numpy as np
    la = np.asarray(lat, dtype=np.float64)
    lo = np.asarray(lon, dtype=np.float64)
    if la.size != lo.size:
        raise ValueError("lat and lon must be the same length")
    d = float(cell_deg)
    kept, counts = [], {}
    for c in l3_grid_cells(d):
        n = int(np.count_nonzero(
            (la >= c["lat_min"]) & (la < c["lat_max"])
            & (lo >= c["lon_min"]) & (lo < c["lon_max"])))
        if n:
            counts[c["name"]] = n
        if n >= int(min_events):
            cc = dict(c)
            cc["n_explore_events"] = n
            cc["mc_label"] = str(mc_label)
            kept.append(cc)
    audit = {
        "cell_deg": d, "min_events": int(min_events), "mc_label": str(mc_label),
        "n_cells_declared": len(l3_grid_cells(d)),
        "n_cells_measurable": len(kept),
        "n_cells_occupied": len(counts),
        "note": L3_MEASURABILITY_NOTE,
        "window": "EXPLORATION ONLY (Rule 4.1)",
    }
    return kept, audit


# ------------------------------------------------------------- membership -----
def contains(box: dict, lat, lon, depth_km=None) -> bool:
    """Box membership, half-open in the upper edge so the L3 grid tiles exactly.

    L1 boxes are inclusive at both edges (that is how `regions_d.contains` reads them
    and how the CSVs were downloaded); L3 grid cells are half-open above, because a
    closed grid double-counts every shared edge and the cell count would then not be a
    partition of anything.
    """
    if box.get("layer") == "L3":
        ok = (box["lat_min"] <= lat < box["lat_max"]
              and box["lon_min"] <= lon < box["lon_max"])
    else:
        ok = (box["lat_min"] <= lat <= box["lat_max"]
              and box["lon_min"] <= lon <= box["lon_max"])
    if ok and depth_km is not None and box.get("depth_max_km") is not None:
        ok = depth_km <= box["depth_max_km"]
    return bool(ok)


def region_mask(box: dict, lat, lon, depth_km=None):
    """Vectorised `contains` -- the form the scan actually uses."""
    import numpy as np
    la = np.asarray(lat, dtype=np.float64)
    lo = np.asarray(lon, dtype=np.float64)
    if box.get("layer") == "L2":
        m = np.zeros(la.size, dtype=bool)
        for name in box["members"]:
            m |= region_mask(region(name), la, lo, depth_km)
        return m
    if box.get("layer") == "L3":
        m = ((la >= box["lat_min"]) & (la < box["lat_max"])
             & (lo >= box["lon_min"]) & (lo < box["lon_max"]))
    else:
        m = ((la >= box["lat_min"]) & (la <= box["lat_max"])
             & (lo >= box["lon_min"]) & (lo <= box["lon_max"]))
    if depth_km is not None and box.get("depth_max_km") is not None:
        m &= (np.asarray(depth_km, dtype=np.float64) <= box["depth_max_km"])
    return m


def region(name: str) -> dict:
    """An L1 or L2 region by name, with its lattice rule id attached."""
    key = str(name).strip().lower().replace("-", "_").replace(" ", "_")
    for r in L1_REGIONS + L2_REGIONS:
        if r["name"] == key:
            out = dict(r)
            out["lattice_rule_id"] = LATTICE_RULE_ID
            out["box_caveat"] = RD.BOX_CAVEAT
            return out
    raise KeyError("no S1-lattice-v1 region named %r" % (name,))


def redundancy_label(l2_event_ids, l1_event_ids_by_region,
                     fraction: float = REDUNDANCY_FRACTION):
    """§S1.1(a): is this L2 (class) finding just one L1 box wearing a class hat?

    Returns {'redundant': bool, 'dominant_l1': name|None, 'fraction': f, 'label': ...}.
    Computed on EVENT IDS, not on box geometry, because two boxes can overlap in area
    and share no events, and it is the events that carry the statistic.
    """
    ids = set(l2_event_ids)
    if not ids:
        return {"redundant": False, "dominant_l1": None, "fraction": 0.0,
                "label": None, "rule": REDUNDANCY_RULE}
    best_name, best_f = None, 0.0
    for name, sub in l1_event_ids_by_region.items():
        f = len(ids & set(sub)) / float(len(ids))
        if f > best_f:
            best_name, best_f = name, f
    red = bool(best_f >= float(fraction))
    return {"redundant": red, "dominant_l1": best_name, "fraction": float(best_f),
            "label": (REDUNDANT_WITH_L1 if red else None), "rule": REDUNDANCY_RULE}


# ------------------------------------------- per-property-family exclusions ----
# §S1.1(c): "Any region that has ever seeded a hypothesis is excluded by name from
# that hypothesis family's scan, permanently. This list only grows."
#
# Keys are PROPERTY FAMILIES as `engine/properties.py` names them. A family absent
# from this map has no seeded region YET; adding one is an append, never an edit.
SEEDED_REGIONS = {
    "solid_tide": tuple(RD.EXCLUDED_BY_NAME),        # K-092, §P7-22 Ratification 5
}

SEEDING_EXCLUSION_RULE = (
    "SEARCHER.md §S1.1(c): any region that has ever SEEDED a hypothesis is excluded "
    "BY NAME from that hypothesis family's scan, permanently, and the list only "
    "grows. The exclusion is PER PROPERTY FAMILY -- Alaska-Aleutians is out of every "
    "solid-tide scan (§P7-22 Ratification 5, K-092) and may still appear as an L1 "
    "region for a family unrelated to that seed, e.g. day-of-week. Checked on the "
    "region NAME and never on a boundary, because a boundary can be got wrong.")


class SeededRegionInScan(AssertionError):
    """A region that seeded a hypothesis appeared in that family's scan (§S1.1(c))."""


def excluded_names_for_family(family: str) -> tuple:
    """The by-name exclusion list for one property family. Empty tuple if none."""
    return tuple(SEEDED_REGIONS.get(str(family), ()))


def assert_family_exclusions(family: str, region_names):
    """Scan-entry gate: refuse a region set carrying that family's seeding region.

    For the solid-tide family this DELEGATES to `regions_d.assert_alaska_excluded`
    rather than re-implementing the refusal, so the engine has exactly one by-name
    check and a future edit cannot repair one copy and miss the other.
    """
    names = [str(n).strip().lower().replace("-", "_").replace(" ", "_")
             for n in region_names]
    fam = str(family)
    if fam == "solid_tide":
        RD.assert_alaska_excluded(names)        # one implementation, not two
        return names
    bad = [n for n in names if n in excluded_names_for_family(fam)]
    if bad:
        raise SeededRegionInScan(
            "property family %r declares seeded region(s) %s, and they are present "
            "in the scan's region set. %s" % (fam, ", ".join(bad),
                                              SEEDING_EXCLUSION_RULE))
    return names


def filter_regions_for_family(family: str, regions):
    """Drop the family's seeded regions from a region list, and say which were dropped."""
    ex = set(excluded_names_for_family(family))
    kept = [r for r in regions if r["name"] not in ex]
    dropped = [r["name"] for r in regions if r["name"] in ex]
    return kept, {"family": str(family), "dropped": dropped,
                  "rule": SEEDING_EXCLUSION_RULE}


# --------------------------------------------------------------- the record ---
def declaration_block(l3_cells=None) -> dict:
    """The frozen, hash-affecting lattice declaration.

    `l3_cells` is the FROZEN measurable-cell list from `l3_measurable_cells`. It is
    hash-affecting (§S1.1(b)): a scan run against a different surviving cell list is a
    different scan and gets a different hash and a new EXPLORE_COUNT line (SP-1.5).
    Passing None declares the L3 layer as GRID-ONLY, i.e. the geometry with no
    measurability filter applied yet.
    """
    cells = list(l3_cells or [])
    return {
        "lattice_rule_id": LATTICE_RULE_ID,
        "rule": LATTICE_RULE,
        "overlap_price": OVERLAP_PRICE,
        "redundancy_rule": REDUNDANCY_RULE,
        "layers": {
            "L1": [{k: v for k, v in r.items() if k != "source"}
                   for r in L1_REGIONS],
            "L2": [{k: v for k, v in r.items() if k != "source"}
                   for r in L2_REGIONS],
            "L3": {
                "cell_deg": L3_CELL_DEG,
                "min_events": L3_MIN_EVENTS,
                "n_cells_declared": len(l3_grid_cells()),
                "measurability_note": L3_MEASURABILITY_NOTE,
                "frozen_cells": [c["name"] for c in cells],
                "n_frozen_cells": len(cells),
            },
        },
        "n_L1": N_L1, "n_L2": N_L2, "n_L3_frozen": len(cells),
        "bird_class_source": BIRD_CLASS_SOURCE,
        "box_source_comcat": BOX_SOURCE_COMCAT,
        "seeded_regions_by_family": {k: list(v) for k, v in SEEDED_REGIONS.items()},
        "seeding_exclusion_rule": SEEDING_EXCLUSION_RULE,
        "regions_d": RD.declaration_block(),        # embedded verbatim, not restated
    }


def declaration_digest(l3_cells=None) -> str:
    return splits.config_hash(declaration_block(l3_cells))[:12]


def region_strata(l3_cells=None, q: float = 0.10, n_properties: int = 1,
                  n_mag_strata: int = 1):
    """The {L1, L2, L3} strata table in `engine/strata.py`'s shape, budget-identical.

    `m_s` is the layer's cell count = n_regions_in_layer x n_properties x
    n_mag_strata; `q_s` is set EQUAL to `q` in every stratum, which is the one
    assignment that satisfies `sum_s m_s q_s = m q` identically for any layer sizes
    and is therefore the one that cannot be tuned. A different split of the budget
    across layers is a declarable choice and is NOT offered here (S-9: one value per
    construction choice, no alternatives run).
    """
    n3 = len(list(l3_cells or []))
    per = int(n_properties) * int(n_mag_strata)
    table = [
        {"name": "L1", "m_s": N_L1 * per, "q_s": float(q)},
        {"name": "L2", "m_s": N_L2 * per, "q_s": float(q)},
        {"name": "L3", "m_s": n3 * per, "q_s": float(q)},
    ]
    return [t for t in table if t["m_s"] > 0]


if __name__ == "__main__":          # pragma: no cover - operator convenience
    import json
    b = declaration_block()
    b["digest"] = declaration_digest()
    print(json.dumps(b, indent=2))
