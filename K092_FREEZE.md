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


---

## SECOND DATED CORRECTION (2026-08-22T08:27:21.294232+00:00): the D-13 null constant is WRONG, and the quadrant is named twice with inconsistent labels

Ruled at HYPOTHESIS_LEDGER.md section P7-25 (Popper, 2026-08-22) on the evidence of the D-1 run
(commit 1ca3bf2, results_k092_d1.json + results_k092_d1_null.json). Both wrong passages above are
LEFT IN PLACE and retracted by amendment. Un-editing is impossible; no line is removed.

### (A) "the uniform expectation (0.25)" is the WRONG NULL and is retracted

The D-13 section above hashes the prediction against "the uniform expectation (0.25)". That constant
is exact for a PURE SINUSOID. The frozen scalar is not one. The time-uniform occupancy of the frozen
quadrant IN THE FROZEN SCALAR, measured over 4,504,845 samples at the 162 seed sites and epochs, is
0.38267 pooled, and per site it ranges 0.3167 to 0.4508 with standard deviation 0.0275. It is below
0.25 at NONE of the 162 sites. A prospective test against 0.25 therefore inflates a nominal 1 percent
test to a 4 to 8 percent false-fire probability at the realistic N, rising with N.

THE CORRECTION WAS PRE-AUTHORISED AND THE RIGHT NUMBER WAS ALREADY ON THE RECORD BEFORE THE SEED WAS
SCORED. Section P7-24 SP-2, committed 2026-08-13, the same day as this freeze, states that a
level/waveform-phase statistic against a uniform-phase null is wrong by construction.
engine/audit_arcsine.py (commit 8e7ef6a) measured the real-waveform artifact at 48.9 percent in the
lowest quarter where a sinusoid gives 33.3. K092_SCALAR_PROVENANCE.md section 4 measured the app's own
quadrant duty cycle at 0.4000 over 60 days at Sand Point. This freeze's own attachment contradicted
this freeze's own constant and the two were never read together. Recorded as S-18 clause 1.

REPLACEMENT, ruled as repair option (c), a PER-EVENT NULL. Each prospective event i, at its own
epicentral latitude, longitude and origin time, is scored against ITS OWN quadrant duty cycle q_i.
The test statistic is unchanged: the count k of prospective events in the frozen quadrant out of N.
The null is the POISSON-BINOMIAL with parameters {q_i}. The success rule is unchanged: p < 0.01,
one-sided, at the 2029-08-13 readout. The prediction itself is NOT restated.

THE q_i RECIPE IS FROZEN BY THIS CORRECTION AND IS HASH-AFFECTING. Plus or minus 10.0 days centred on
the event time; 1.0-minute step; plus or minus 10-minute central-difference rate; local-cycle minimum
and maximum taken between the refined bracketing maxima; the app's astro.ts scalar through
exp_k092_d1_bridge.mjs; sign-condition quadrant (level < 0 AND rate < 0). Computed by
exp_k092_d1_null.py's existing machinery. engine/sitetide.py is reported alongside as robustness and
is NEVER the scorer. ANY change to span, step, rate half-width, normalisation, scalar or quadrant
definition VOIDS D-13 OUTRIGHT, with no discretion. This is the ratchet that stops the repair from
becoming a reusable escape hatch.

WHY THIS AMENDMENT IS HONEST DESPITE A RUNNING CLOCK: q_i is a deterministic astronomical function of
(lat, lon, t) with zero free parameters and zero catalogue input, so it cannot be tuned and any third
party can recompute it; and the correction runs AGAINST the claim, since 0.25 -> 0.3827 makes the
prediction strictly harder to satisfy.

### (B) D-13 IS RE-LABELLED: DESCRIPTIVE-PRIMARY, NOT DECISIVE. It was never powered.

The repair exposes a defect the repair cannot fix. The seed file gives 33 M>=6.0 in the box 1990-2025
= 0.917/yr, so Gutenberg-Richter at b = 1 puts M>=5.5 at about 2.9/yr and N at about 9 over the
3-year horizon. Section P7-23(B)'s assumed 5-10/yr was an ARC-WIDE rate applied to a sub-segment box
and is corrected. At N = 9 to 12 the corrected D-13 has power 0.07 to 0.08 against a true quadrant
rate of 0.60, an effect nearly 60 percent larger than the null.

D-13 stands and scores at 2029-08-13 exactly as hashed-and-amended, and its result is reported
whatever it is. But A D-13 NULL MAY NOT BE QUOTED AS A BOUND. The decisive within-region arm is D-12,
which has the events.

### (C) THE QUADRANT IS NAMED TWICE AND THE TWO LABELS DISAGREE. The SIGN CONDITIONS BIND.

The claim section above gives the quadrant both as two sign conditions ("tide level below zero AND
tide level decreasing") and as "theta in (pi, 3*pi/2) with theta = 0 at maximum". Those labels are
inconsistent as written: with theta = 0 at the MAXIMUM the scalar reads A*cos(theta) and the
below-neutral-falling quarter is theta in (pi/2, pi). The interval (pi, 3*pi/2) is that SAME PHYSICAL
QUARTER-CYCLE under the x = A*sin(theta) convention engine/audit_arcsine.py uses. Same quarter, same
exact 1/4 sinusoid null, two mutually inconsistent angular labels.

RULED: THE SIGN-CONDITION FORM BINDS AND IS THE FROZEN DEFINITION. The angle-interval form is
commentary and MAY NOT SCORE ANYTHING. Grounds, all predating D-1: the claim sentence states the sign
conditions first and glosses them with the angle; section P7-23(D) ruled the sign-condition form
rotation-free and made that the whole of the licence; it is the form the application displays and the
form Jim read.

AND THE DIFFERENCE IS NOT COSMETIC ON REAL DATA. D-1 measured the two forms on the 162 seed events:
sign-condition 64/162 = 0.3951, angle-interval 49/162 = 0.3025, with 47 events classified differently,
29 percent of the set. They coincide only for a pure sinusoid. The angle-interval figure has NO
matched null and is therefore UNINTERPRETABLE in either direction; it must not be quoted.

### (D) D-13b: A SECOND PROSPECTIVE ARM, FROZEN AND HASHED AT THIS COMMIT, PRICED 0

D-13 cannot answer the question at its own N. D-13b can. It is legitimate to freeze now precisely
because it is WHOLLY PROSPECTIVE: no data for it exists.

PREDICTION, hashed at this commit: of the Alaska-Aleutian box events (51-58 N, 166-152 W) with
M >= 5.0 occurring after this correction's timestamp, the fraction rupturing in the frozen
sign-condition quadrant of the frozen scalar will exceed the PER-EVENT Poisson-binomial expectation
{q_i} under the recipe frozen in (A), at one-sided p < 0.01.

- Threshold M >= 5.0, declared as testing the MAGNITUDE-INDEPENDENT FORM per section P7-23(A.4): a
  null here does NOT refute a majors-specific claim, and that sentence travels with any D-13b headline.
- Horizon 3 years from this correction; final scoring 2029-08-22. Descriptive readouts at 1 year
  (2027-08-22) and 2 years (2028-08-22), non-scoring.
- DECLUSTERED PRIMARY and FULL-CATALOGUE SECONDARY, both declared, both scored, neither substituting
  for the other. An M7 with six M>=5.0 aftershocks delivers seven correlated successes, and D-13's
  hashed binomial has no declustering clause at all, which is its third defect.
- Declustering rule frozen NOW, before any event: an event is a dependent event if it falls within
  30 days AND 150 km of an equal-or-larger prior event in the scoring set.
- alpha = 0.01, one-sided.

No result is claimed here. This correction is a commitment and a retraction, not evidence.
