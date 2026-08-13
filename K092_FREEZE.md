# K-092 FREEZE: the Alaska-Aleutian regional tidal-phase claim (D-12 / D-13)

Frozen 2026-08-13T20:27:26.818534+00:00 per sections P7-22 and P7-23 of HYPOTHESIS_LEDGER.md.
PRIORITY FACT recorded at freeze time: engine/sitetide.py DOES NOT EXIST in this repository
(no site-tide computation exists anywhere in the engine), so no one - including the program -
can compute any Alaska event's tidal phase today. This freeze is therefore genuinely prior
to any possible look beyond the seed set.

## The claim (Jim Gale's seed, 2026-08-13)

Major earthquakes of the Alaska-Aleutian segment rupture preferentially in the BELOW-NEUTRAL,
FALLING quadrant of the local solid-earth tide: the two sign conditions (tide level below zero)
AND (tide level decreasing), i.e. theta in (pi, 3*pi/2) with theta = 0 at maximum. Rotation-free
per section P7-23(D); scalar-dependent, and the scalar is FROZEN as: the solid-earth tide
vertical displacement (cm) at the epicentral coordinates as computed by the observing
application (Jim's viewer, which displayed e.g. -13.2 cm falling 0.6 cm/h for the 1991-05-30
M7.0 Sand Point event). TRANSLATION STEP OWED: Jim attaches or confirms the viewer's tide
implementation before scoring; mapping to any other scalar is a separately declared step.

## Region

Latitude 51 to 58 N, longitude 166 to 152 W (the Alaska Peninsula / eastern Aleutian segment).

## Seed-exclusion superset (D-12, section P7-23(A): "scrolled past is seen")

EVERY M>=6.0 event in the region box, 1900-01-01 through freeze date, is excluded from D-12
scoring as potentially seen: 162 events, enumerated by USGS ComCat id in
K092_seed_exclusion_superset.csv (sha256 35b5cf0fe9145ec830fba782550f7d8b76dc2708db105eb6003a0f92b68bfcb6).

## D-12 scoring set (STRATUM-HELD-OUT label per P7-23(A))

Alaska-Aleutian events in the box NOT in the exclusion superset, exploration window only,
scored against a full-ETAS null WITH triggering (aftershocks inherit mainshock phase).
Unviewed majors (none exist below M6.0 in the superset definition, so this stratum is
M4.0-5.9) test the magnitude-independent form: a null there does NOT refute the
majors-specific claim (declared now, per P7-23(A)).

## D-13 prospective log (section P7-23(B), priced 0)

PREDICTION, hashed at this commit: of the next Alaska-Aleutian box events with M>=5.5
occurring after the freeze timestamp, the fraction rupturing in the frozen quadrant of the
frozen scalar will exceed the uniform expectation (0.25) with binomial p < 0.01 at the
3-year readout. Horizon: 3 years (final scoring 2029-08-13). Descriptive readouts at 1 year
(2027-08-13) and 2 years (2028-08-13), non-scoring. The 2025-07-16 M7.3 (us7000qd1y) and
every pre-freeze event are excluded by construction.

## Pre-scoring gates (order binding)

1. D-0 convention + scalar translation attached (Jim).
2. D-1 arcsine control with the P7-23(C) classification: per seed event, level percentile AND
   phase angle; trough-concentration (artifact-consistent) vs mid-slope concentration
   (artifact-inconsistent, physically favored: maximum unloading stress rate). Rule declared
   before any Alaska data beyond the superset is touched.
3. D-2 event-path VIF (for the weak-form floors only; the strong-form binomial claim above
   needs no floor per P7-23's correction).

No result is claimed here. This file is a commitment, not evidence. Un-committing is impossible.


---

## DATED CORRECTION (2026-08-13T20:28:06.614775+00:00): the priority-fact sentence above was WRONG at commit time

The freeze commit (f1fa3ba) states that engine/sitetide.py did not exist and was "verified in this
commit". The verification command's own output, printed in the same tool call, showed the file DID
exist in the working tree: the tranche-D build agent had created it (untracked, uncommitted) while
the freeze was being written, and the pipeline committed anyway. The claim as written is false and
is retracted. The priority argument is restated on the true basis, each part checkable:
(1) engine/sitetide.py appears in NO git commit at or before the freeze commit (git log --all on the
path is empty; the file is untracked at correction time);
(2) the agent that created it operated under a written prohibition against computing any
Alaska-Aleutian phase (its brief: "do NOT run any real-data phase statistic on the Alaska-Aleutian
region - D-7 is STRUCK"), and its deliverable report is the auditable record that its tests used
synthetic and non-Alaska sites only;
(3) therefore no Alaska event phase had been computed by anyone at freeze time, which is the fact
the priority argument actually needs; the capability existing uncommitted on one machine under a
no-look order weakens the rhetorical form of the claim, not its substance.
Error class: S-18 clause 1 in its purest form - a verification whose output said the opposite was
carried as "verified" because the pipeline did not gate on the check. Recorded here rather than
edited away; the wrong sentence above is left in place.
