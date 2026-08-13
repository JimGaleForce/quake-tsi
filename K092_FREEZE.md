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
