"""D-0: the phase-convention gate. Literature-derived, declared once, hash-affecting.

§P7-22 Q5 RULED that D-7 does not run and that the convention is set from published
literature, with Jim's provenance *"sought but non-determinative"*. The tests below are
that ruling as invariants: the block cites the literature, the digest moves when the
convention moves, an undeclared or altered convention is REFUSED, and the provenance
slot cannot change the hash.
"""

from __future__ import annotations

import pytest

from engine import conventions_d as D0, splits


def test_the_convention_is_the_tanaka_one_and_cites_it():
    assert D0.CONVENTION_ID == "D0-tanaka-stressmax-v1"
    assert "local MAXIMUM" in D0.CONVENTION_RULE or "MAXIMUM" in D0.CONVENTION_RULE
    joined = " ".join(D0.CITATIONS)
    assert "Tanaka" in joined and "2002" in joined
    assert "Cochran" in joined and "2004" in joined


def test_the_direction_derivation_is_literature_and_says_so():
    """§P7-22 Q5: a look-derived DIRECTION is the exact defect D-0 prevents."""
    b = D0.convention_block()
    assert b["direction_derivation"].startswith("LITERATURE")
    assert "no catalogue" in b["direction_derivation"]
    assert "NOT derived from any look" in b["rule"]


def test_the_scope_flag_forbids_reading_phase_zero_as_failure_promoting():
    """This build has no receiver mechanism, so Cochran's direction cannot transfer."""
    assert "psi remains unreportable" in D0.CONVENTION_SCOPE
    assert "RECEIVER-RESOLVED" in D0.CONVENTION_SCOPE
    assert "may NOT be transferred" in D0.CONVENTION_SCOPE


def test_the_block_is_hash_affecting():
    cfg = {"mag_target": 4.5, "phase_convention": D0.convention_block()}
    h0 = splits.config_hash(cfg)
    tampered = dict(cfg)
    blk = dict(D0.convention_block())
    blk["rule"] = blk["rule"].replace("MAXIMUM", "MINIMUM")
    tampered["phase_convention"] = blk
    assert splits.config_hash(tampered) != h0


def test_an_undeclared_or_altered_convention_is_refused():
    with pytest.raises(D0.ConventionNotDeclared):
        D0.assert_declared({})
    blk = dict(D0.convention_block())
    blk["scalar"] = "something_else"
    with pytest.raises(D0.ConventionNotDeclared):
        D0.assert_declared({"phase_convention": blk})
    ok = D0.assert_declared({"phase_convention": D0.convention_block()})
    assert ok["digest"] == D0.convention_digest()


def test_provenance_is_required_for_the_scalar_and_not_for_the_direction():
    """§P7-23(D) amends §P7-22 Q5, and only one half of it moved."""
    assert D0.PROVENANCE_JIM["determinative_for_scalar"] is True
    assert D0.PROVENANCE_JIM["determinative_for_direction"] is False
    assert D0.PROVENANCE_JIM["status"].startswith("RECEIVED")
    assert D0.PROVENANCE_JIM["answer"] is not None and "earth-tides-globe" in D0.PROVENANCE_JIM["answer"]
    # the direction is fixed by the literature and an answer does not move it
    assert "does not set the direction" in D0.PROVENANCE_JIM["rule"]
    # the answer's prose is not hashed -- the scalar it fixes is
    assert "provenance_jim" not in D0.convention_block()


def test_a_frozen_quadrant_claim_is_refused_until_provenance_answers():
    """§P7-23(D): 'below neutral and falling' is a different quadrant in each scalar.
    The answer was RECEIVED 2026-08-13; the refusal path is preserved by clearing it
    temporarily, so the gate is proven to be a gate and the received state a state."""
    saved = D0.PROVENANCE_JIM["answer"]
    try:
        D0.PROVENANCE_JIM["answer"] = None
        with pytest.raises(D0.ScalarProvenanceRequired):
            D0.assert_scalar_provenance("the D-12 frozen quadrant claim")
    finally:
        D0.PROVENANCE_JIM["answer"] = saved
    assert "earth-tides-globe" in D0.assert_scalar_provenance()


def test_the_default_scalar_is_labelled_provisional_and_says_what_it_does_not_license():
    b = D0.convention_block()
    assert b["scalar_status"].startswith("PROVISIONAL, PENDING PROVENANCE")
    assert "NOT licensed" in b["scalar_status"]
    assert "rotation-invariant" in b["scalar_status"]


def test_the_convention_names_the_sitetide_scalar_it_is_attached_to():
    from engine import sitetide as ST
    assert D0.SCALAR_NAME == "sitetide_" + ST.SCALAR_FOR_PHASE
    assert "EXTENSION POSITIVE" in D0.SCALAR_SIGN
