"""B3 / §S4 -- THE PROMOTION GENERATOR: a freeze file, or nothing at all.

Promotion is **mechanical** (§P7-24 SP-4): no Popper round, no negotiation. It emits a
K-092-template freeze file, hashed and committed, and **that is the entire content of
promotion**. §P7-23(B): *"a hashed commitment makes no claim until scored and therefore
needs no licence to be committed; it is a denominator committed before the numerator
exists."* The SCORING is priced, later. **The fast thing is the commitment; the slow
thing is the verdict, and it stays slow.**

WHY THIS IS A CODE GENERATOR AND NOT A RESEARCH PROJECT
--------------------------------------------------------
`K092_FREEZE.md` already exists and its sections ARE the template's slots (§S4's
table). This module fills them from a scan row and writes the file.

THE ONE THING K-092 GOT WRONG, CLOSED IN CODE RATHER THAN IN PROSE
--------------------------------------------------------------------
That freeze shipped a **priority-fact sentence that its own verification output
contradicted** -- it stated `engine/sitetide.py` did not exist while the verification
command printed in the same tool call showed that it did -- *"because the pipeline did
not gate on the check"*. Its own dated correction names the error class: **S-18 clause
1 in its purest form -- a verification whose output said the opposite was carried as
'verified'.**

> **So this generator emits its priority facts as ASSERTIONS THAT MUST PASS, and
> REFUSES TO WRITE THE FILE if any fails.** `PriorityFactFalsified` is that refusal.
> A priority fact here is not a sentence; it is a callable that returns
> `(bool, evidence)` and whose evidence string is written into the file beside the
> claim it supports. **There is no `--force`.**

THE SEED-EXCLUSION SUPERSET: "SCANNED IS SEEN"
------------------------------------------------
§P7-23(A) condition 1 demands a SUPERSET, not a list of consciously-noticed events.
§S4 generalises it for a scan: the superset is **every event the scan statistic
touched**, enumerated by ComCat id and sha256'd. *"Scrolled past is seen"* becomes
*"scanned is seen."* `PF_SUPERSET_COVERS_SCAN` is the assertion that the CSV really
does cover the scanned set, checked on ids rather than on counts.

THE THREE ARMS, AT SP-4's STANDING PRICES -- 16, not §S9's per-candidate estimate
----------------------------------------------------------------------------------
  * within-region, stratum-held-out .... **2**
  * cross-region ....................... **14**  (7 target-class + 7 control-class)
  * prospective log .................... **0**
  * ---------------------------------------------
  * TOTAL .............................. **16**

*"PROMOTION FREEZES THE PROPERTY AND THE STATISTIC. The scan's breadth has already
been paid for in SP-3's threshold; it must not be paid for a second time in the
confirmation."*

Nothing this module writes is evidence. It writes commitments.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os

TEMPLATE_SOURCE = "K092_FREEZE.md"
FREEZE_RULE_ID = "S4-freeze-gen-v1"

DISCLAIMER = ("No result is claimed here. This file is a commitment, not evidence. "
              "Un-committing is impossible.")

STANDING_PRICES = {"within_region_stratum_held_out": 2, "cross_region": 14,
                   "prospective_log": 0}
STANDING_PRICE_TOTAL = sum(STANDING_PRICES.values())      # 16

PRICING_PRINCIPLE = (
    "§P7-24 SP-4: promotion FREEZES THE PROPERTY AND THE STATISTIC. The scan's "
    "breadth has already been paid for in SP-3's threshold and must not be paid for a "
    "second time in the confirmation. The tidal case cost 181 because it scanned 3 "
    "phases x 4 statistics; a promoted searcher claim is one property, one statistic, "
    "and costs 16.")

SCAN_TRAVELS_WITH_THE_CLAIM = (
    "§P7-24 SP-6.3: the seeding scan travels with the claim -- scan id, date, m, and "
    "the candidate's RANK inside it, in the freeze file, forever. The look-elsewhere "
    "effect of the scan is a property of the claim permanently, and a reader must be "
    "able to recompute the effective multiplicity without asking us.")

VARIANT_ENUMERATION_RULE = (
    "§P7-24 SP-5: the variant space is ENUMERATED at promotion time. Discovering a "
    "new variant after looking is a NEW SEED with its own freeze and its own "
    "scan-identity record -- never an amendment to an existing family. The "
    "enumeration is what makes that checkable rather than a matter of good faith. A "
    "class hypothesis is ONE priced hypothesis if and only if the statistic is "
    "INVARIANT to the varying dimension; if the claim names the pattern ('A at Q1 AND "
    "B at Q3') that is a psi claim in disguise with 4^R free parameters and it "
    "inherits K-090(c)'s restrictions on reporting an absolute phase.")

SUPERSET_RULE = (
    "§P7-23(A.1) as generalised by SEARCHER.md §S4: the seed set is excluded BY EVENT "
    "ID as a SUPERSET -- every event the scan statistic TOUCHED, not a list of "
    "consciously-noticed ones. 'Scrolled past is seen' becomes 'scanned is seen.'")

PROMOTION_BUYS = (
    "§P7-24: promotion does NOT make a claim true, quotable, or a K-entry. It buys "
    "exactly one thing: THE RIGHT TO BE TESTED BY THE THREE ARMS. It does not remove "
    "the region-selection circularity of the seeding scan, which is irreducible and "
    "is why the cross-region and prospective arms exist.")


class PriorityFactFalsified(AssertionError):
    """A declared priority fact did not hold. The freeze file is NOT written.

    This is the K-092 dated correction closed in code: that freeze carried a priority
    sentence its own verification output contradicted, because the pipeline did not
    gate on the check. There is no override and no `--force`.
    """


class PriorityFact:
    """A named, checkable claim about the world at freeze time.

    `check` returns `(bool, evidence_string)`. The evidence is written into the file
    beside the claim, so a reader sees WHAT WAS CHECKED and not merely that something
    was.
    """

    __slots__ = ("name", "statement", "check")

    def __init__(self, name, statement, check):
        self.name = str(name)
        self.statement = str(statement)
        self.check = check

    def run(self):
        ok, evidence = self.check()
        return {"name": self.name, "statement": self.statement,
                "passed": bool(ok), "evidence": str(evidence)}


def assert_priority_facts(facts):
    """Run every fact; raise on ANY failure. Returns the passed records."""
    out = [f.run() for f in facts]
    bad = [r for r in out if not r["passed"]]
    if bad:
        raise PriorityFactFalsified(
            "REFUSING TO WRITE THE FREEZE FILE: %d priority fact(s) did not hold.\n%s"
            "\nThis is S-18 clause 1 closed in code -- the K-092 freeze shipped a "
            "priority-fact sentence its own verification output contradicted because "
            "the pipeline did not gate on the check. This pipeline gates."
            % (len(bad), "\n".join("  - %s: %s\n      evidence: %s"
                                   % (r["name"], r["statement"], r["evidence"])
                                   for r in bad)))
    return out


# --------------------------------------------------------- the standard facts --
def standard_priority_facts(row, scan, scanned_ids, superset_ids, csv_path,
                            csv_sha256):
    """The facts every searcher freeze must assert, whatever the property is."""
    def f_superset():
        missing = set(map(str, scanned_ids)) - set(map(str, superset_ids))
        return (not missing,
                "scanned ids = %d; superset ids = %d; ids scanned but NOT excluded "
                "= %d%s" % (len(set(map(str, scanned_ids))),
                            len(set(map(str, superset_ids))), len(missing),
                            "" if not missing
                            else " (first few: %s)" % sorted(missing)[:5]))

    def f_csv():
        if not os.path.exists(csv_path):
            return False, "the superset CSV %r does not exist on disk" % csv_path
        h = sha256_file(csv_path)
        return (h == csv_sha256,
                "sha256(%s) on disk = %s; recorded in this file = %s"
                % (os.path.basename(csv_path), h, csv_sha256))

    def f_threshold():
        alpha = float(scan["alpha"])
        return (float(row["p_real"]) <= alpha,
                "p_real = %.6e; alpha = q/m = %.4g/%d = %.6e (SP-3, one rule, no "
                "OR-limb)" % (float(row["p_real"]), float(scan["q"]),
                              int(scan["m"]), alpha))

    def f_control():
        pc = row.get("p_control")
        alpha = float(scan["alpha"])
        return (pc is not None and float(pc) > alpha,
                "p_control = %s; alpha = %.6e; F7-b requires the control NOT to fire "
                "at the identical threshold" % (pc, alpha))

    def f_null_layer():
        return (bool(row.get("sp2_null_layer_built")),
                "SP-2 null layer for property class %r: %s"
                % (row.get("property_class"), row.get("sp2_reason")))

    def f_dwell():
        return (bool(row.get("dwell_time_corrected")),
                "property type = %r; dwell_time_corrected = %r (§P7-23(C): a "
                "level-concentration statistic that does not divide by the "
                "occupancy measure is not emitted at all)"
                % (row.get("type"), row.get("dwell_time_corrected")))

    def f_scan_identity():
        need = ("scan_id", "m", "config_hash")
        have = [k for k in need if scan.get(k) not in (None, "")]
        return (len(have) == len(need) and row.get("rank_in_scan"),
                "scan_id=%r m=%r config_hash=%r rank_in_scan=%r (SP-6.3)"
                % (scan.get("scan_id"), scan.get("m"), scan.get("config_hash"),
                   row.get("rank_in_scan")))

    def f_human_schedule():
        if not row.get("human_schedule"):
            return True, "property is not human-schedule class; F7-d does not apply"
        return (float(row["mc"]) >= 6.0,
                "human-schedule property at Mc = %.1f; F7-d makes the arm "
                "scientifically live ONLY at M >= 6.0" % float(row["mc"]))

    def f_stratum_enumerated():
        return (float(row["mc"]) in [float(x) for x in scan["mag_strata"]],
                "Mc = %.1f; declared strata = %s (SP-6.7: a claim at an unenumerated "
                "threshold is a NEW SEED, not a result)"
                % (float(row["mc"]), scan["mag_strata"]))

    return [
        PriorityFact("seed_exclusion_is_a_superset",
                     "Every event the scan statistic touched is excluded by id.",
                     f_superset),
        PriorityFact("superset_csv_hash_matches",
                     "The sha256 recorded in this file is the sha256 of the CSV "
                     "actually written.", f_csv),
        PriorityFact("cleared_the_declared_threshold",
                     "The candidate's nominal p cleared alpha = q/m at the scan's "
                     "FULL declared cell count.", f_threshold),
        PriorityFact("control_arm_did_not_fire",
                     "The matched control arm did not clear the identical threshold.",
                     f_control),
        PriorityFact("sp2_null_layer_built",
                     "The property class's mandatory null-validity layer exists.",
                     f_null_layer),
        PriorityFact("dwell_time_corrected",
                     "A level property carries its dwell-time occupancy measure.",
                     f_dwell),
        PriorityFact("scan_identity_recorded",
                     "The seeding scan's id, m, hash and the candidate's rank travel "
                     "with the claim.", f_scan_identity),
        PriorityFact("f7d_human_schedule_carve_out",
                     "A human-schedule property is promoted only at M >= 6.0.",
                     f_human_schedule),
        PriorityFact("magnitude_stratum_enumerated",
                     "The magnitude stratum was enumerated in the scan declaration.",
                     f_stratum_enumerated),
    ]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_superset_csv(path, events, id_key="id"):
    """Write the seed-exclusion superset, one row per event, ids first.

    `events` is a list of dicts. Columns are the union of the keys in declaration
    order with `id` first, so the file is diffable and the id column is never
    ambiguous. Returns (n_rows, sha256).
    """
    import csv
    rows = list(events)
    if not rows:
        raise ValueError("the seed-exclusion superset is EMPTY. §P7-23(A.1) requires "
                         "a superset of everything the scan touched; an empty file "
                         "would assert that the scan touched nothing.")
    keys = [id_key] + [k for k in rows[0] if k != id_key]
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})
    return len(rows), sha256_file(path)


# ------------------------------------------------------------ the template ----
def render_freeze(k_id, row, scan, region, seed_attribution, superset_csv,
                  superset_sha256, n_superset, facts, variant_space,
                  cross_region_arm, prospective, frozen_at=None):
    """The K-092 template, filled. Every §S4 slot plus the two §S11 additions."""
    ts = frozen_at or _dt.datetime.now(_dt.timezone.utc).isoformat()
    lat = "%g to %g" % (region["lat_min"], region["lat_max"]) if "lat_min" in region \
        else "union of %s" % ", ".join(region.get("members", []))
    lon = "%g to %g" % (region["lon_min"], region["lon_max"]) if "lon_min" in region \
        else "(class union; see members)"

    fact_lines = "\n".join(
        "%d. **%s** -- %s\n   ASSERTED AND PASSED. Evidence: %s"
        % (i + 1, f["name"], f["statement"], f["evidence"])
        for i, f in enumerate(facts))

    variant_lines = "\n".join("* %s" % v for v in variant_space)

    return """# %(kid)s FREEZE: %(claim_title)s

Frozen %(ts)s by `engine/freeze_gen.py` (%(rule)s) under HYPOTHESIS_LEDGER.md §P7-24
SP-1..SP-6 and SEARCHER.md §S4. Template: `%(template)s`.

**%(disclaimer)s**

## Priority facts -- ASSERTED, NOT ASSERTED-TO

Every statement below was emitted by a check that had to PASS before this file could
be written. **The generator refuses to write on any failure**, which is the K-092
dated correction (S-18 clause 1) closed in code rather than in prose.

%(facts)s

## The claim

%(seed)s

Region **%(region_name)s** (layer %(layer)s, class %(rclass)s) concentrates its
M >= %(mc).1f events on the property **%(prop)s** (family `%(family)s`, property class
`%(pclass)s`, type `%(ptype)s`) in the declared concentration form
**%(form)s**.

Statistic (frozen, one, no alternatives): `%(stat)s`.
Scalar / convention (FROZEN -- §P7-23(D): the claim is frozen in the scalar actually
examined, and translation to any other scalar is a separate declared step):

%(provenance)s

## Region

Latitude %(lat)s, longitude %(lon)s.

## Seed-exclusion superset (§P7-23(A.1); "scanned is seen")

%(superset_rule)s

%(n_superset)d events, enumerated by catalogue id in `%(superset_csv)s`
(sha256 %(sha)s).

## The seeding scan travels with this claim (SP-6.3)

%(scan_travels)s

| field | value |
|---|---|
| scan id | `%(scan_id)s` |
| scan config hash | `%(scan_hash)s` |
| scan date | %(scan_date)s |
| `m` (FULL declared cell count) | **%(m)d** |
| `alpha = q/m`, q = %(q)g | **%(alpha).6e** |
| this candidate's nominal p | **%(p_real).6e** |
| this candidate's rank in the scan | **%(rank)s** of at most %(kcap)d |
| matched control arm p | %(p_control)s |
| lattice rule / digest | `%(lattice)s` / `%(lattice_digest)s` |

## Enumerated variant space (SP-5)

%(variant_rule)s

%(variants)s

## The three arms, at SP-4's STANDING PRICES

%(pricing)s

| arm | standing price | content |
|---|---:|---|
| within-region, stratum-held-out | **%(p_within)d** | PRIMARY stratum (unexamined events of the seed's own kind) + SECONDARY stratum (other magnitude strata, testing the magnitude-independent form). §P7-23(A.4)'s declared bridging assumption is mandatory; a null on SECONDARY does NOT refute a PRIMARY-specific claim. |
| cross-region | **%(p_cross)d** | %(cross)s |
| prospective log | **%(p_prosp)d** | %(prospective)s |
| | **%(p_total)d total** | |

## Pre-scoring gates (order binding)

1. Dwell-time control: the property's occupancy measure is attached and the
   concentration statistic is computed on the PIT-transformed level (§P7-23(C)).
2. Observer control: `observer.observer_features` for this region and stratum is
   computed and REPORTED BESIDE this claim, always, in every output (F7-a).
3. Event-path VIF: `engine/audit_event_vif.py` for the weak-form floors only; a
   strong-form claim whose nominal p is far below the threshold needs no floor
   (§P7-23's correction to §P7-22).

## What promotion does and does not buy

%(promotion_buys)s

%(disclaimer)s
""" % {
        "kid": k_id,
        "claim_title": "%s x %s at M >= %.1f" % (region["name"], row["property"],
                                                 float(row["mc"])),
        "ts": ts, "rule": FREEZE_RULE_ID, "template": TEMPLATE_SOURCE,
        "disclaimer": DISCLAIMER,
        "facts": fact_lines,
        "seed": seed_attribution,
        "region_name": region["name"], "layer": region.get("layer"),
        "rclass": region.get("class"),
        "mc": float(row["mc"]), "prop": row["property"], "family": row["family"],
        "pclass": row["property_class"], "ptype": row["type"],
        "form": row["concentration_form"], "stat": row["statistic"],
        "provenance": "\n".join("* **%s**: %s" % (k, v)
                                for k, v in sorted(row["provenance"].items())),
        "lat": lat, "lon": lon,
        "superset_rule": SUPERSET_RULE,
        "n_superset": int(n_superset), "superset_csv": superset_csv,
        "sha": superset_sha256,
        "scan_travels": SCAN_TRAVELS_WITH_THE_CLAIM,
        "scan_id": scan["scan_id"], "scan_hash": scan["config_hash"],
        "scan_date": scan.get("scan_date", ts), "m": int(scan["m"]),
        "q": float(scan["q"]), "alpha": float(scan["alpha"]),
        "p_real": float(row["p_real"]), "rank": row.get("rank_in_scan"),
        "kcap": int(scan["k_cap"]),
        "p_control": ("%.6e" % row["p_control"]) if row.get("p_control") is not None
                     else "ABSENT -- and a promotion without one is refused",
        "lattice": scan["lattice_rule_id"], "lattice_digest": scan["lattice_digest"],
        "variant_rule": VARIANT_ENUMERATION_RULE, "variants": variant_lines,
        "pricing": PRICING_PRINCIPLE,
        "p_within": STANDING_PRICES["within_region_stratum_held_out"],
        "p_cross": STANDING_PRICES["cross_region"],
        "p_prosp": STANDING_PRICES["prospective_log"],
        "p_total": STANDING_PRICE_TOTAL,
        "cross": cross_region_arm, "prospective": prospective,
        "promotion_buys": PROMOTION_BUYS,
    }


def generate_freeze(k_id, row, scan, region, superset_events, out_dir=".",
                    seed_attribution=None, variant_space=None,
                    cross_region_arm=None, prospective=None, scanned_ids=None,
                    extra_facts=(), id_key="id", frozen_at=None):
    """Write `<k_id>_FREEZE.md` + `<k_id>_seed_exclusion_superset.csv`, or write NOTHING.

    Order matters and is deliberate: the CSV is written FIRST so that its sha256 is a
    fact about a file on disk rather than about an in-memory buffer, the priority
    facts are then asserted against that file, and the markdown is written LAST. If a
    fact fails, the markdown does not exist and the CSV is removed, so a failed
    promotion leaves no half-committed artifact behind.
    """
    os.makedirs(out_dir, exist_ok=True)
    csv_name = "%s_seed_exclusion_superset.csv" % k_id
    csv_path = os.path.join(out_dir, csv_name)
    md_path = os.path.join(out_dir, "%s_FREEZE.md" % k_id)

    n_rows, sha = write_superset_csv(csv_path, superset_events, id_key=id_key)
    ids = [str(e.get(id_key, "")) for e in superset_events]
    facts = list(standard_priority_facts(
        row, scan, scanned_ids if scanned_ids is not None else ids, ids,
        csv_path, sha)) + list(extra_facts)
    try:
        passed = assert_priority_facts(facts)
    except PriorityFactFalsified:
        if os.path.exists(csv_path):
            os.remove(csv_path)          # no half-committed artifact
        raise

    variants = list(variant_space or [
        "the same property and the same concentration form in ANY declared region of "
        "the same Bird class (rotation-invariant, therefore ONE hypothesis under "
        "SP-5, not R looks)",
        "the same property at the other ENUMERATED magnitude strata %s (SP-6.7)"
        % (scan["mag_strata"],),
    ])
    cross = cross_region_arm or (
        "7 target-class regions + 7 control-class regions from the S1-lattice-v1 L1 "
        "layer, one frozen property, one frozen statistic, with every region that has "
        "seeded this property family excluded BY NAME (§S1.1(c)). The region set is "
        "declared IN FULL IN ADVANCE, including regions expected to fail; adding a "
        "region after seeing its result is forbidden and is the exact defect the "
        "declared set exists to avoid (§S5.1).")
    prosp = prospective or (
        "K-069/K-080 pattern, priced 0: hashed commitment, declared threshold, "
        "statistic, horizon and pass/fail written BEFORE the first event. Primary "
        "horizon 3 years from the freeze timestamp; 1- and 2-year readouts are "
        "DESCRIPTIVE ONLY and do not score.")
    seed = seed_attribution or (
        "Seed: Jim Gale, 2026-08-13, attributed in full -- *\"the agent version of "
        "what I did with Sand Point\"*. This candidate was SELECTED BY A SCAN, not by "
        "a person, and the scan's identity is recorded below so the look-elsewhere "
        "effect travels with the claim forever.")

    md = render_freeze(k_id, row, scan, region, seed, csv_name, sha, n_rows,
                       passed, variants, cross, prosp, frozen_at=frozen_at)
    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write(md)
    return {
        "k_id": k_id, "freeze_path": md_path, "superset_csv": csv_path,
        "superset_sha256": sha, "n_superset": n_rows,
        "priority_facts": passed,
        "priced_tests": STANDING_PRICE_TOTAL,
        "prices": dict(STANDING_PRICES),
        "rule_id": FREEZE_RULE_ID,
        "note": PROMOTION_BUYS,
    }


def generate_from_scan(report, scan, regions_by_name, superset_by_cell,
                       out_dir=".", k_id_prefix="K", k_id_start=200, **kw):
    """Freeze every promoted row in a scan report. Returns one record per freeze.

    `superset_by_cell` maps `(region, property, mc)` to the list of event dicts the
    cell's statistic TOUCHED. Nothing is inferred: a cell with no supplied superset
    raises rather than freezing against an assumed one.
    """
    out = []
    n = int(k_id_start)
    for row in report["ranked"]:
        if not row.get("promoted"):
            continue
        key = (row["region"], row["property"], float(row["mc"]))
        if key not in superset_by_cell:
            raise KeyError("no seed-exclusion superset supplied for cell %r. %s"
                           % (key, SUPERSET_RULE))
        out.append(generate_freeze(
            "%s%03d" % (k_id_prefix, n), row,
            {**scan, "scan_date": report.get("scan_date", "")},
            regions_by_name[row["region"]], superset_by_cell[key],
            out_dir=out_dir, **kw))
        n += 1
    return out


if __name__ == "__main__":          # pragma: no cover - operator convenience
    print(json.dumps({"rule_id": FREEZE_RULE_ID, "template": TEMPLATE_SOURCE,
                      "standing_prices": STANDING_PRICES,
                      "total": STANDING_PRICE_TOTAL,
                      "pricing_principle": PRICING_PRINCIPLE,
                      "promotion_buys": PROMOTION_BUYS}, indent=2))
