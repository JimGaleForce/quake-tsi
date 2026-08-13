"""D-0 -- THE TRANCHE D PHASE CONVENTION, DECLARED IN ONE PLACE, LITERATURE-DERIVED.

HYPOTHESIS_LEDGER.md §P7-22 Q5, RULED:

  > **RULED. D-7 does not run.** Its only output is the convention, and **a convention
  > derived from a look at the seeding region is the exact defect D-0 exists to
  > prevent** -- an unscored look is still a look, and it makes the eventual direction
  > *more* look-dependent, not less. **A literature-derived convention is
  > pre-registered by someone who had no stake in our outcome, which is strictly
  > stronger than anything we could declare about ourselves.**
  >
  > **D-0 stands as the gate, with one amendment:** ask Jim which display he was
  > reading -- the provenance is worth having and may explain the seed outright --
  > **but the declared convention is set from published literature independent of his
  > answer. Provenance informs understanding; it does not set the direction.**

So: the convention below is taken from the published tidal-triggering literature. It
is NOT derived from any look at Alaska-Aleutians, at any catalogue, or at any
exploration window. `PROVENANCE_JIM` is a slot that records Jim's answer when he gives
it and is explicitly NON-DETERMINATIVE: nothing in this module reads it.

THE CONVENTION (S-9: one value, written one way)
------------------------------------------------
**Phase 0 = the LOCAL MAXIMUM of the declared tidal scalar at the site.** Phase
increases linearly (in TIME, not in angle of any single constituent) from one local
maximum to the next, so phase 0 and phase 2*pi are successive maxima and phase pi is
the intervening minimum.

This is the **Tanaka convention**, from:

  * Tanaka, S., M. Ohtake and H. Sato (2002), "Evidence for tidal triggering of
    earthquakes as revealed from statistical analysis of global data",
    J. Geophys. Res. 107(B10), 2211, doi:10.1029/2001JB001577. The tidal phase angle
    is defined with **0 deg at the maximum of the tidal stress** and +-180 deg at the
    minimum, with the angle obtained by **linear interpolation in time between
    successive extrema** of the computed tidal series at the epicentre.
  * Cochran, E. S., J. E. Vidale and S. Tanaka (2004), "Earth tides can trigger
    shallow thrust fault earthquakes", Science 306, 1164-1166,
    doi:10.1126/science.1103961. Uses the same phase definition and reports the excess
    of shallow thrust events near the failure-promoting extremum.
  * Tsuruoka, H., M. Ohtake and H. Hayakawa (1995), Geophys. J. Int. 122, 183-194 --
    the earlier statement of the same "phase measured between successive extrema of
    the local tidal stress" construction.

WHAT THIS BUILD IS AND IS NOT ENTITLED TO CALL "STRESS"
-------------------------------------------------------
Tanaka/Cochran place phase 0 at the maximum of a **receiver-resolved Coulomb** tidal
stress. **This build cannot compute that**: MINING_CATALOG F1-30's mechanism piece is
NOT built (no receiver fault geometry, no ocean loading), and F1-30's own Pit says
*"do not quote signs you cannot defend"* (M-006). `engine/sitetide.py` therefore emits
a **mechanism-free** scalar -- the degree-2 body-tide areal strain and the mean normal
stress derived from it, EXTENSION POSITIVE -- and this module declares phase 0 at
**that** scalar's local maximum.

**The consequence, stated before any run.** The convention fixes the phase ORIGIN
identically for signal and for null, which is all a rotation-invariant concentration
statistic needs (§P7-22(a): Kuiper V, Watson U^2 and R1 are rotation-invariant, so the
origin cannot manufacture or destroy concentration). It does **NOT** license reading
"phase near 0" as "failure-promoting": that identification needs the receiver
mechanism this build does not have. **psi stays unreportable** (§P7-22(a)), and the
direction-of-effect language in Cochran 2004 may not be transferred to our scalar
without the mechanism. This is a SCOPE FLAG, not a caveat, and `CONVENTION_SCOPE`
carries it in every record this module emits.

HASH-AFFECTING
--------------
`convention_block()` returns the frozen dict that a Tranche D config must embed. It
goes through `splits.config_hash` exactly like every other construction choice, so a
silent change of convention changes the config hash and cannot be resumed into an
existing session. `assert_declared` is the check that a config carries it unmodified.

Nothing in this module is evidence.
"""

from __future__ import annotations

from . import splits

# ------------------------------------------------------------- the declaration --
CONVENTION_ID = "D0-tanaka-stressmax-v1"

CONVENTION_RULE = (
    "Phase 0 rad = the local MAXIMUM of the declared site tidal scalar; phase "
    "increases LINEARLY IN TIME between successive maxima, so phase pi is the "
    "intervening minimum and phase 2*pi is the next maximum. Tanaka, Ohtake & Sato "
    "(2002, JGR 107(B10) 2211) convention, as used by Cochran, Vidale & Tanaka "
    "(2004, Science 306, 1164). Literature-derived per HYPOTHESIS_LEDGER.md "
    "§P7-22 Q5; NOT derived from any look at any catalogue or region.")

# The scalar the convention is attached to. One name, one sign, one unit (S-9).
SCALAR_NAME = "sitetide_areal_strain"
SCALAR_SIGN = "EXTENSION POSITIVE (areal strain e_tt + e_pp > 0 is dilatation)"
SCALAR_UNITS = "dimensionless strain"

# The literature the convention is taken from, verbatim enough to be checkable.
CITATIONS = (
    "Tanaka, S., Ohtake, M. & Sato, H. (2002), J. Geophys. Res. 107(B10), 2211, "
    "doi:10.1029/2001JB001577 -- phase angle 0 deg at the maximum of the tidal "
    "stress, +-180 deg at the minimum, linearly interpolated in time between "
    "successive extrema.",
    "Cochran, E. S., Vidale, J. E. & Tanaka, S. (2004), Science 306, 1164-1166, "
    "doi:10.1126/science.1103961 -- same phase definition, shallow thrust faults.",
    "Tsuruoka, H., Ohtake, M. & Hayakawa, H. (1995), Geophys. J. Int. 122, 183-194 "
    "-- earlier statement of the between-successive-extrema construction.",
)

CONVENTION_SCOPE = (
    "SCOPE FLAG (F1-30 Pit / M-006 / §P7-22(a)). Tanaka's phase 0 is the maximum of "
    "a RECEIVER-RESOLVED COULOMB tidal stress. This build has no receiver mechanism "
    "and no ocean loading, so phase 0 here is the maximum of a MECHANISM-FREE "
    "body-tide scalar. That fixes the phase ORIGIN identically for signal and null, "
    "which is all a rotation-invariant concentration statistic requires -- and it "
    "does NOT license reading 'phase near 0' as 'failure-promoting'. psi remains "
    "unreportable under §P7-22(a); Cochran 2004's direction-of-effect language may "
    "NOT be transferred to this scalar.")

# --------------------------------- §P7-23(D): THE AMENDMENT, AND WHAT IT MOVES --
# §P7-23(D) landed AFTER this module's first version and amends §P7-22 Q5. The
# amendment is narrow and it is worth stating exactly which half moved:
#
#   > **RULED: the claim is frozen in THE SCALAR JIM'S VIEWER USED**, and any
#   > translation to a stress scalar is a **separate, separately declared step** that
#   > inherits K-090(c)'s restrictions. **This amends §P7-22's D-0 ruling: the
#   > provenance question is no longer "worth having," it is REQUIRED**, because it
#   > identifies the scalar in which the frozen claim is stated. My earlier ruling
#   > treated provenance as non-determinative for the *direction*; that stands. It is
#   > determinative for the *scalar*, and I had that wrong.
#
# So: the DIRECTION stays literature-derived and provenance-independent, exactly as
# above. The SCALAR is now provenance-determined, because vertical displacement, tidal
# potential, volumetric strain and Coulomb stress on a local thrust differ by
# geometry-dependent phase offsets -- "below neutral and falling" is a DIFFERENT
# QUADRANT in each. `DECLARED_SCALAR_STATUS` and `assert_scalar_provenance` are that
# ruling in code: a FROZEN claim may not be stated in this module's default scalar
# while the provenance is unanswered.
PROVENANCE_JIM = {
    "question": ("which tidal display was Jim reading when he made the Sand Point "
                 "observation -- vertical displacement, water level, tidal potential, "
                 "volumetric strain, or a resolved stress -- and from which source?"),
    "answer": None,
    "status": "REQUIRED, NOT RECEIVED (§P7-23(D))",
    "determinative_for_scalar": True,
    "determinative_for_direction": False,
    "rule": ("§P7-23(D) amends §P7-22 Q5. DIRECTION: still literature-derived and "
             "provenance-independent -- 'Provenance informs understanding; it does "
             "not set the direction.' SCALAR: provenance is REQUIRED and "
             "DETERMINATIVE, because the same words name a different quadrant in "
             "each candidate scalar. A frozen claim may not be stated in this "
             "module's default scalar until the question is answered."),
}

DECLARED_SCALAR_STATUS = (
    "PROVISIONAL, PENDING PROVENANCE (§P7-23(D)). `%s` is the scalar this build "
    "EMITS and the one the D-8/D-9 class-level concentration statistics run on -- "
    "which is legitimate, because a rotation-invariant concentration statistic does "
    "not depend on which scalar fixes the origin so long as signal and null share it. "
    "It is NOT licensed as the scalar a FROZEN quadrant claim is stated in: that one "
    "is fixed by Jim's viewer and by nothing else." % SCALAR_NAME)


class ScalarProvenanceRequired(AssertionError):
    """A frozen quadrant claim was stated in a scalar provenance has not fixed."""


def assert_scalar_provenance(what="a frozen quadrant claim"):
    """§P7-23(D): refuse to freeze a quadrant claim in an unprovenanced scalar.

    Returns the answer when it exists. This does NOT gate the concentration
    statistics -- those are rotation-invariant and scalar-agnostic by §P7-22(a) -- it
    gates the one thing §P7-23(D) says it must: stating a claim whose content is a
    named quadrant.
    """
    ans = PROVENANCE_JIM.get("answer")
    if not ans:
        raise ScalarProvenanceRequired(
            "%s cannot be stated: §P7-23(D) makes the provenance question REQUIRED "
            "and DETERMINATIVE FOR THE SCALAR, and it is unanswered. %s %s"
            % (what, PROVENANCE_JIM["question"], DECLARED_SCALAR_STATUS))
    return ans


def convention_block() -> dict:
    """The frozen, hash-affecting dict a Tranche D config must embed.

    Deliberately excludes `PROVENANCE_JIM`'s ANSWER TEXT: the answer fixes which
    scalar a frozen claim is stated in, and that arrives through `scalar` (which IS
    hashed) rather than through free prose whose wording would move the hash without
    moving the meaning.
    """
    return {
        "convention_id": CONVENTION_ID,
        "rule": CONVENTION_RULE,
        "scalar": SCALAR_NAME,
        "scalar_sign": SCALAR_SIGN,
        "scalar_units": SCALAR_UNITS,
        "scalar_status": DECLARED_SCALAR_STATUS,
        "direction_derivation": ("LITERATURE (§P7-22 Q5, unchanged by §P7-23(D)); "
                                 "no catalogue, no region, no window"),
        "scalar_derivation": ("PROVENANCE (§P7-23(D)): determinative for the scalar, "
                              "not for the direction"),
        "citations": list(CITATIONS),
        "scope_flag": CONVENTION_SCOPE,
    }


def convention_digest() -> str:
    """SHA-256 of the frozen block, first 12 hex -- what a report quotes."""
    return splits.config_hash(convention_block())[:12]


class ConventionNotDeclared(AssertionError):
    """A Tranche D statistic was asked for without the D-0 convention in the config."""


def assert_declared(cfg: dict) -> dict:
    """Refuse to proceed unless `cfg['phase_convention']` is the declared block.

    Hash-affecting by consequence, not by exhortation: an altered block has a
    different digest and this raises rather than silently scoring under a convention
    nobody declared.
    """
    got = (cfg or {}).get("phase_convention")
    want = convention_block()
    if got is None:
        raise ConventionNotDeclared(
            "Tranche D requires the D-0 phase convention in the config under key "
            "'phase_convention'. " + CONVENTION_RULE)
    if got != want:
        raise ConventionNotDeclared(
            "the config's phase convention is not the declared D-0 block "
            "(declared digest %s, config digest %s). §P7-22 Q5 permits exactly one, "
            "literature-derived. " % (convention_digest(),
                                      splits.config_hash(got)[:12]))
    return {"convention_id": CONVENTION_ID, "digest": convention_digest()}


def record(extra: dict | None = None) -> dict:
    """The block plus its digest and the non-determinative provenance slot."""
    out = convention_block()
    out["digest"] = convention_digest()
    out["provenance_jim"] = dict(PROVENANCE_JIM)
    if extra:
        out.update(extra)
    return out


if __name__ == "__main__":          # pragma: no cover - operator convenience
    import json
    print(json.dumps(record(), indent=2))
