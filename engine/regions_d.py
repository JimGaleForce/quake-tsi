"""D-6 -- KURIL and CASCADIA as DECLARED CONSTANTS, from published tectonic classification.

Two things this module is, and one thing it is emphatically not.

**What it is (1): the two missing boxes.** §K92-0's data audit found that the two
segments the seed names are the two that are not on disk -- *"Kuril (the Japan box
stops at 46 N) and Cascadia (the California box stops at 42 N)"*. Their bounding boxes
are declared here as CONSTANTS so that a later download is a mechanical fetch against a
frozen declaration rather than a boundary chosen while looking at a catalogue.

**What it is (2): the tectonic partition rule id §P7-22(b) permits.** RULED:

  > **A new declared tectonic partition rule id is PERMITTED for K-092**, on four
  > conditions: **(i)** it is defined from **published tectonic classification, not
  > from catalogue data** -- which is what removes the forking hazard entirely, since a
  > physically-defined partition cannot have been tuned to the outcome; **(ii)**
  > declared, frozen and hash-affecting **before any K-092 statistic is computed**;
  > **(iii)** Rule 4.1 respected -- nothing in its derivation touches the holdout
  > window or the K-080 census; **(iv)** `R2b-lon6-active` is **NOT also run on this
  > question.**

  > **And a trap in the existing partition that must be recorded:** `R2b-lon6-active`
  > sector 0 (-180 to -120) **contains Alaska-Aleutians almost exactly** ... **Any
  > regsum arm run under lon6 would therefore include the seeding region in the
  > confirmation set.** The tectonic partition must exclude Alaska explicitly, by
  > name, not by trusting a boundary.

`LON6_FORBIDDEN_ON_THIS_QUESTION` and `assert_partition_rule` are condition (iv) in
code: asking this module for `R2b-lon6-active` raises. `ALASKA_EXCLUDED_BY_NAME` and
`assert_alaska_excluded` are the by-name exclusion, and `EXCLUDED_BY_NAME` is checked
on the region NAME, never on a boundary, because a boundary can be got wrong and a name
cannot.

**What it is NOT: a selector.** Nothing here reads a catalogue, a census, an event
count, a holdout window or `results_k080_census.json`. The boxes and the class labels
come from the published plate-boundary literature cited below and from nothing else --
which is precisely §P7-22(b)(i)'s reason for permitting the new rule at all. Rule 4.1
(§P6-4) forbids the K-080 census as a selector regardless.

PUBLISHED SOURCES FOR THE CLASSIFICATION (cited, per §P7-22(b)(i))
------------------------------------------------------------------
* **Bird, P. (2003), "An updated digital model of plate boundaries", Geochemistry
  Geophysics Geosystems 4(3), 1027, doi:10.1029/2001GC000252.** The standard digital
  plate-boundary model, and -- decisively for this module -- the standard source of the
  **plate-boundary CLASS** labels used here: SUB (subduction zone), OTF/CTF (transform),
  CRB/OSR (rift / spreading), CCB (continental convergent boundary). Every region's
  `class` field below is Bird's class for the boundary that region straddles.
* **Bird, P. & Kagan, Y. Y. (2004), "Plate-tectonic analysis of shallow seismicity",
  Bull. Seismol. Soc. Am. 94(6), 2380-2399.** The seismicity-by-boundary-class
  companion, which is the published precedent for partitioning a global catalogue by
  Bird's classes rather than by geography.
* **Coats, R. R. (1962) / Jicha & Kay (2018)** and the USGS Alaska-Aleutian arc
  descriptions for the Aleutian segment's extent; **Fedotov et al. (1982)** and the
  Kamchatka-Kuril arc literature for the Kuril-Kamchatka trench segment;
  **Heaton & Hartzell (1987), "Earthquake hazards on the Cascadia subduction zone",
  Science 236, 162-168**, and **Wang & Trehu (2016), J. Geodynamics 98, 1-18**, for the
  Cascadia subduction zone's along-strike extent (Mendocino triple junction ~40.3 N to
  the Nootka fault ~50 N) and its downdip width.

The boxes below are generous rectangular hulls of the published arc traces, rounded to
whole or half degrees. **A rectangle is not a plate boundary**, and `BOX_CAVEAT` says
so on every record: the boxes are a data-selection envelope, not a tectonic claim.
"""

from __future__ import annotations

from . import splits

# --------------------------------------------------------------- the rule id --
REGION_RULE_ID_D = "D6-tectonic-bird2003-v1"

PARTITION_RULE = (
    "D6-tectonic-bird2003-v1: regions are named boxes whose CLASS is Bird (2003)'s "
    "plate-boundary class for the boundary the box straddles (SUB / OTF / CRB / CCB). "
    "The partition is PUBLISHED-CLASSIFICATION-DERIVED: no catalogue, no census, no "
    "event count, no holdout window enters its definition (§P7-22(b)(i), §P6-4 Rule "
    "4.1). Alaska-Aleutians is excluded BY NAME, not by boundary (§P7-22(b)).")

LON6_FORBIDDEN_ON_THIS_QUESTION = (
    "§P7-22(b) condition (iv): `R2b-lon6-active` is NOT also run on the K-092 "
    "question. Running both would be precisely the alternative S-9 prices, and it "
    "would cost 12 more tests to buy an ambiguity. It is additionally FORBIDDEN on "
    "this question for a substantive reason: lon6 sector 0 (-180..-120) contains "
    "Alaska-Aleutians almost exactly, so a regsum arm under lon6 would carry the "
    "SEEDING REGION into the confirmation set.")

ALASKA_EXCLUDED_BY_NAME = (
    "§P7-22 Ratification 5: Alaska-Aleutians is exploration-only in BOTH directions "
    "-- the M7.3 sits in the temporal holdout AND the region generated the hypothesis "
    "-- so it is unscoreable twice over and is excluded from the pooled statistic, "
    "from the partition BY NAME, and from any downstream summary. Excluded by NAME "
    "and never by trusting a boundary, because a boundary can be got wrong.")

BOX_CAVEAT = (
    "A rectangular box is a DATA-SELECTION ENVELOPE, not a tectonic claim. The boxes "
    "are generous rectangular hulls of published arc traces (Bird 2003 boundary "
    "geometry), rounded to whole or half degrees; membership in a box is not a "
    "statement that every event in it occurred on the named boundary.")

# Names that may never enter a K-092 confirmation set, checked as NAMES.
EXCLUDED_BY_NAME = ("alaska_aleutians", "alaska", "aleutians")

FORBIDDEN_RULE_IDS = ("R2b-lon6-active",)


# ------------------------------------------------------------- the two boxes --
# lat_min, lat_max, lon_min, lon_max in degrees; lon in [-180, 180].
KURIL = {
    "name": "kuril",
    "long_name": "Kuril-Kamchatka trench segment",
    "class": "SUB",
    "lat_min": 43.0, "lat_max": 56.5,
    "lon_min": 144.0, "lon_max": 167.0,
    "depth_max_km": 200.0,
    "why_this_box": ("the Japan box on disk stops at 46 N (§K92-0), so the "
                     "Kuril-Kamchatka arc from the Hokkaido corner (~43 N, 145 E) to "
                     "the Kamchatka peninsula (~56 N, 163 E) is unrepresented. Extent "
                     "from Bird (2003)'s KU/PA-OK boundary trace and the "
                     "Kamchatka-Kuril arc literature (Fedotov et al. 1982)."),
    "held_out_role": "CONFIRMATION (subduction set, §P7-22 Ratification 5)",
}

CASCADIA = {
    "name": "cascadia",
    "long_name": "Cascadia subduction zone",
    "class": "SUB",
    "lat_min": 39.5, "lat_max": 51.0,
    "lon_min": -131.0, "lon_max": -120.5,
    "depth_max_km": 100.0,
    "why_this_box": ("the California box on disk stops at 42 N (§K92-0). Cascadia "
                     "runs from the Mendocino triple junction (~40.3 N) to the Nootka "
                     "fault (~50 N); the box is extended half a degree beyond each "
                     "end. Extent from Heaton & Hartzell (1987) and Wang & Trehu "
                     "(2016); class SUB from Bird (2003) JF-NA."),
    "held_out_role": "CONFIRMATION (subduction set, §P7-22 Ratification 5)",
}

REGIONS_D = (KURIL, CASCADIA)


def region(name: str) -> dict:
    """The declared box by name, with its caveats attached."""
    for r in REGIONS_D:
        if r["name"] == str(name).lower():
            out = dict(r)
            out["box_caveat"] = BOX_CAVEAT
            out["partition_rule_id"] = REGION_RULE_ID_D
            out["sources"] = ("Bird (2003) G3 4(3) 1027; Bird & Kagan (2004) BSSA "
                              "94(6) 2380; Heaton & Hartzell (1987) Science 236 162; "
                              "Wang & Trehu (2016) J. Geodyn. 98 1")
            return out
    raise KeyError("no declared Tranche D region named %r" % (name,))


def contains(box: dict, lat, lon, depth_km=None) -> bool:
    """Box membership. Longitudes are compared in [-180, 180], no wrap case arises."""
    ok = (box["lat_min"] <= lat <= box["lat_max"]
          and box["lon_min"] <= lon <= box["lon_max"])
    if ok and depth_km is not None and box.get("depth_max_km") is not None:
        ok = depth_km <= box["depth_max_km"]
    return bool(ok)


# ------------------------------------------------------------------- guards ---
class ForbiddenPartitionRule(AssertionError):
    """`R2b-lon6-active` was requested on the K-092 question (§P7-22(b)(iv))."""


class SeedingRegionInConfirmationSet(AssertionError):
    """Alaska-Aleutians was found in a K-092 confirmation set (§P7-22 Rat. 5)."""


def assert_partition_rule(rule_id: str) -> str:
    """Refuse `R2b-lon6-active` on this question. Condition (iv), in code."""
    r = str(rule_id)
    if r in FORBIDDEN_RULE_IDS:
        raise ForbiddenPartitionRule(
            "partition rule %r is FORBIDDEN on the K-092 question. %s"
            % (r, LON6_FORBIDDEN_ON_THIS_QUESTION))
    if r != REGION_RULE_ID_D:
        raise ForbiddenPartitionRule(
            "K-092 declares exactly one partition rule id, %r; got %r. S-9 permits "
            "one rule, written one way." % (REGION_RULE_ID_D, r))
    return r


def assert_alaska_excluded(region_names) -> list:
    """Refuse any confirmation set containing the seeding region, checked by NAME."""
    names = [str(n).strip().lower().replace("-", "_").replace(" ", "_")
             for n in region_names]
    bad = [n for n in names if n in EXCLUDED_BY_NAME]
    if bad:
        raise SeedingRegionInConfirmationSet(
            "the seeding region is present in a K-092 region set: %s. %s"
            % (", ".join(bad), ALASKA_EXCLUDED_BY_NAME))
    return names


# --------------------------------------------------------------- the record ---
def declaration_block() -> dict:
    """The frozen, hash-affecting block a Tranche D config embeds for the partition."""
    return {
        "partition_rule_id": REGION_RULE_ID_D,
        "rule": PARTITION_RULE,
        "regions": [{k: v for k, v in r.items() if k != "why_this_box"}
                    for r in REGIONS_D],
        "excluded_by_name": list(EXCLUDED_BY_NAME),
        "forbidden_rule_ids": list(FORBIDDEN_RULE_IDS),
        "lon6_note": LON6_FORBIDDEN_ON_THIS_QUESTION,
        "alaska_note": ALASKA_EXCLUDED_BY_NAME,
        "box_caveat": BOX_CAVEAT,
    }


def declaration_digest() -> str:
    return splits.config_hash(declaration_block())[:12]


if __name__ == "__main__":          # pragma: no cover - operator convenience
    import json
    b = declaration_block()
    b["digest"] = declaration_digest()
    print(json.dumps(b, indent=2))
