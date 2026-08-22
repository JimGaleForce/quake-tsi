"""B3 / §S4: the promotion generator, and the one refusal that is its whole point.

`K092_FREEZE.md` shipped a **priority-fact sentence that its own verification output
contradicted** -- *"because the pipeline did not gate on the check"* -- and its dated
correction names the class: **S-18 clause 1 in its purest form**. This suite exists to
prove that this pipeline gates:

  * every priority fact is a CALLABLE that must return True, and its evidence string
    is written into the file beside the claim it supports;
  * a single failing fact means **NO FILE IS WRITTEN AT ALL**, and the CSV that was
    written first is REMOVED, so a failed promotion leaves no half-committed artifact;
  * there is no `--force` -- checked as the absence of any override parameter.

Plus the slots §P7-24 adds to §S4's table: SP-6.3's scan identity + rank, SP-5's
enumerated variant space, and SP-4's standing price of **16**, not §S9's per-candidate
estimate.
"""

from __future__ import annotations

import hashlib
import os

import pytest

from engine import freeze_gen as FG, lattice_s1 as LAT


def _row(**kw):
    base = dict(region="japan", layer="L1", region_class="SUB",
                property="tide_phase", family="solid_tide", type="phase",
                property_class="level-waveform-phase", mc=6.0,
                statistic="kuiper_V_event_path", value=2.5, p_real=1.0e-9,
                p_control=0.42, concentration_form="arc",
                dwell_time_corrected=True, sp2_null_layer_built=True,
                sp2_reason="dwell attached", human_schedule=False,
                rank_in_scan=1,
                provenance={"source": "engine/sitetide.py",
                            "convention": "D0-tanaka-stressmax-v1",
                            "scalar_definition": "areal_strain"})
    base.update(kw)
    return base


def _scan(**kw):
    base = dict(scan_id="scan_20260813T000000", config_hash="deadbeef", m=1000,
                q=0.10, alpha=1.0e-4, k_cap=30, mag_strata=[4.5, 5.0, 5.5, 6.0],
                lattice_rule_id=LAT.LATTICE_RULE_ID,
                lattice_digest=LAT.declaration_digest())
    base.update(kw)
    return base


def _events(n=12):
    return [{"id": "us700%04d" % i, "time": "2001-0%d-01T00:00:00Z" % (i % 9 + 1),
             "mag": 6.1 + 0.01 * i, "lat": 35.0, "lon": 140.0} for i in range(n)]


# ------------------------------------------------------------- the happy path --
def test_a_clean_promotion_writes_both_files_and_all_facts_pass(tmp_path):
    rec = FG.generate_freeze("K900", _row(), _scan(), LAT.region("japan"),
                             _events(), out_dir=str(tmp_path))
    assert os.path.exists(rec["freeze_path"]) and os.path.exists(rec["superset_csv"])
    assert all(f["passed"] for f in rec["priority_facts"])
    assert len(rec["priority_facts"]) == 9
    md = open(rec["freeze_path"], encoding="utf-8").read()
    assert FG.DISCLAIMER in md
    assert "ASSERTED AND PASSED" in md
    for f in rec["priority_facts"]:
        assert f["name"] in md and f["evidence"] in md


def test_the_recorded_sha256_is_the_sha256_of_the_file_on_disk(tmp_path):
    """The K-092 defect in miniature: a hash recorded from a buffer is not a fact
    about a file. The CSV is written FIRST so the hash is a fact about disk."""
    rec = FG.generate_freeze("K901", _row(), _scan(), LAT.region("japan"),
                             _events(), out_dir=str(tmp_path))
    h = hashlib.sha256(open(rec["superset_csv"], "rb").read()).hexdigest()
    assert h == rec["superset_sha256"]
    assert h in open(rec["freeze_path"], encoding="utf-8").read()


# ------------------------------------------------------------- the refusals ---
@pytest.mark.parametrize("kw,expect", [
    (dict(p_real=1.0e-2), "cleared_the_declared_threshold"),
    (dict(p_control=1.0e-9), "control_arm_did_not_fire"),
    (dict(p_control=None), "control_arm_did_not_fire"),
    (dict(sp2_null_layer_built=False), "sp2_null_layer_built"),
    (dict(dwell_time_corrected=False), "dwell_time_corrected"),
    (dict(rank_in_scan=None), "scan_identity_recorded"),
    (dict(human_schedule=True, mc=5.0), "f7d_human_schedule_carve_out"),
    (dict(mc=6.7), "magnitude_stratum_enumerated"),
])
def test_a_failing_priority_fact_refuses_to_write_anything(tmp_path, kw, expect):
    with pytest.raises(FG.PriorityFactFalsified) as e:
        FG.generate_freeze("K902", _row(**kw), _scan(), LAT.region("japan"),
                           _events(), out_dir=str(tmp_path))
    assert expect in str(e.value)
    assert "REFUSING TO WRITE" in str(e.value)
    assert "S-18 clause 1" in str(e.value)
    assert not os.path.exists(str(tmp_path / "K902_FREEZE.md"))
    assert not os.path.exists(str(tmp_path / "K902_seed_exclusion_superset.csv"))


def test_a_superset_that_does_not_cover_the_scanned_set_refuses(tmp_path):
    """§P7-23(A.1) as generalised by §S4: 'scanned is seen'. The superset must be a
    SUPERSET, and the check is on IDS, not on counts."""
    ev = _events(5)
    scanned = [e["id"] for e in ev] + ["us7000zzzz"]      # one the CSV lacks
    with pytest.raises(FG.PriorityFactFalsified) as e:
        FG.generate_freeze("K903", _row(), _scan(), LAT.region("japan"), ev,
                           out_dir=str(tmp_path), scanned_ids=scanned)
    assert "seed_exclusion_is_a_superset" in str(e.value)
    assert "us7000zzzz" in str(e.value)


def test_an_equal_count_but_different_ids_still_refuses(tmp_path):
    """A count check would pass here. An id check does not, which is the point."""
    ev = _events(5)
    scanned = ["other%d" % i for i in range(5)]
    with pytest.raises(FG.PriorityFactFalsified):
        FG.generate_freeze("K904", _row(), _scan(), LAT.region("japan"), ev,
                           out_dir=str(tmp_path), scanned_ids=scanned)


def test_an_empty_superset_refuses_before_any_fact_runs(tmp_path):
    with pytest.raises(ValueError) as e:
        FG.generate_freeze("K905", _row(), _scan(), LAT.region("japan"), [],
                           out_dir=str(tmp_path))
    assert "EMPTY" in str(e.value)


def test_there_is_no_force_flag(tmp_path):
    """The refusal must not be overridable. There is no parameter that lets it be."""
    import inspect
    for fn in (FG.generate_freeze, FG.assert_priority_facts):
        params = set(inspect.signature(fn).parameters)
        assert not (params & {"force", "override", "ignore_facts", "skip_checks"}), fn
    src = inspect.getsource(FG.assert_priority_facts)
    assert "raise PriorityFactFalsified" in src


def test_a_custom_extra_fact_is_run_and_can_refuse(tmp_path):
    good = FG.PriorityFact("sitetide_committed", "the module is in a commit",
                           lambda: (True, "git log --all on the path is non-empty"))
    rec = FG.generate_freeze("K906", _row(), _scan(), LAT.region("japan"), _events(),
                             out_dir=str(tmp_path), extra_facts=[good])
    assert "sitetide_committed" in [f["name"] for f in rec["priority_facts"]]
    bad = FG.PriorityFact("sitetide_absent", "engine/sitetide.py does not exist",
                          lambda: (os.path.exists("engine/sitetide.py") is False,
                                   "os.path.exists says it DOES exist"))
    with pytest.raises(FG.PriorityFactFalsified) as e:
        FG.generate_freeze("K907", _row(), _scan(), LAT.region("japan"), _events(),
                           out_dir=str(tmp_path), extra_facts=[bad])
    # this is literally the K-092 sentence, and this pipeline refuses to write it
    assert "sitetide_absent" in str(e.value)


# ---------------------------------------------------------------- the slots ---
def test_the_seeding_scan_travels_with_the_claim(tmp_path):
    """SP-6.3: scan id, date, m, and the candidate's rank, in the file, forever."""
    rec = FG.generate_freeze("K908", _row(), _scan(), LAT.region("japan"), _events(),
                             out_dir=str(tmp_path))
    md = open(rec["freeze_path"], encoding="utf-8").read()
    assert "scan_20260813T000000" in md
    assert "deadbeef" in md
    assert "**1000**" in md                       # m
    assert "**1** of at most 30" in md            # rank in scan
    assert "1.000000e-04" in md                   # alpha = q/m
    assert "look-elsewhere effect" in md


def test_the_standing_prices_are_sp4s_sixteen_not_s9s_estimate(tmp_path):
    assert FG.STANDING_PRICES == {"within_region_stratum_held_out": 2,
                                  "cross_region": 14, "prospective_log": 0}
    assert FG.STANDING_PRICE_TOTAL == 16
    rec = FG.generate_freeze("K909", _row(), _scan(), LAT.region("japan"), _events(),
                             out_dir=str(tmp_path))
    assert rec["priced_tests"] == 16
    md = open(rec["freeze_path"], encoding="utf-8").read()
    assert "**16 total**" in md
    assert "7 target-class regions" in md
    assert "FREEZES THE PROPERTY AND THE STATISTIC" in md


def test_the_variant_space_is_enumerated_at_promotion_time(tmp_path):
    """SP-5: discovering a variant after looking is a NEW SEED, never an amendment."""
    rec = FG.generate_freeze("K910", _row(), _scan(), LAT.region("japan"), _events(),
                             out_dir=str(tmp_path),
                             variant_space=["same property in any SUB region",
                                            "same property at M >= 5.5"])
    md = open(rec["freeze_path"], encoding="utf-8").read()
    assert "* same property in any SUB region" in md
    assert "NEW SEED" in md and "never an amendment" in md


def test_the_frozen_scalar_and_convention_are_in_the_file(tmp_path):
    """§P7-23(D): the claim is frozen in the SCALAR ACTUALLY EXAMINED."""
    rec = FG.generate_freeze("K911", _row(), _scan(), LAT.region("japan"), _events(),
                             out_dir=str(tmp_path))
    md = open(rec["freeze_path"], encoding="utf-8").read()
    assert "D0-tanaka-stressmax-v1" in md
    assert "areal_strain" in md
    assert "separate declared step" in md


def test_the_file_never_claims_a_result(tmp_path):
    rec = FG.generate_freeze("K912", _row(), _scan(), LAT.region("japan"), _events(),
                             out_dir=str(tmp_path))
    md = open(rec["freeze_path"], encoding="utf-8").read()
    assert md.count(FG.DISCLAIMER) >= 2           # head and foot
    assert "THE RIGHT TO BE TESTED BY THE THREE ARMS" in md
    assert "irreducible" in md                    # the circularity is not hidden


def test_generate_from_scan_freezes_only_promoted_rows(tmp_path):
    report = {"ranked": [_row(promoted=True, rank_in_scan=1),
                         _row(promoted=False, rank_in_scan=2,
                              property="lunar_synodic_phase", family="ephemeris")]}
    supers = {("japan", "tide_phase", 6.0): _events()}
    out = FG.generate_from_scan(report, _scan(), {"japan": LAT.region("japan")},
                                supers, out_dir=str(tmp_path), k_id_start=800)
    assert len(out) == 1 and out[0]["k_id"] == "K800"
    assert os.path.exists(str(tmp_path / "K800_FREEZE.md"))
    assert not os.path.exists(str(tmp_path / "K801_FREEZE.md"))


def test_generate_from_scan_refuses_a_cell_with_no_supplied_superset(tmp_path):
    """Nothing is INFERRED: a cell without its superset raises, never freezes."""
    report = {"ranked": [_row(promoted=True, rank_in_scan=1)]}
    with pytest.raises(KeyError):
        FG.generate_from_scan(report, _scan(), {"japan": LAT.region("japan")}, {},
                              out_dir=str(tmp_path))


def test_the_csv_has_the_id_column_first_and_one_row_per_event(tmp_path):
    p = str(tmp_path / "s.csv")
    n, sha = FG.write_superset_csv(p, _events(7))
    lines = open(p, encoding="utf-8").read().strip().splitlines()
    assert n == 7 and len(lines) == 8
    assert lines[0].split(",")[0] == "id"
    assert len(sha) == 64
