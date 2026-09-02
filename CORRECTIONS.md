# Dated corrections to frozen artifacts

Frozen results files are never edited; corrections to their prose fields are recorded here.

## 2026-08-13: results_phase_clock_null.json, gaps[0], ocean-loading sentence

The frozen F-016 results file (committed d5d656d) states that neither real series "carries ocean
loading". That is WRONG for the Xue-Lu Zenodo volumetric series (data/xue_lu_zenodo/Tidal_Vol.txt):
it was computed with SPOTL including TPXO ocean loading plus local west-coast models, per Lu, Xue,
Yue, Zhuang and Zhao (2025) section 2.3 and this repo's own data/xue_lu_zenodo/INVENTORY.md. The
measured instrument-response and sampling-floor numbers in the file are unaffected (they are
properties of the series as given, whatever its physical content); only the gaps-field prose was
wrong. The generating script's gaps text is corrected as of this date for future runs. Surfaced by
the sixth adversarial verification pass on outgoing correspondence, which checked the artifact
against the data provider's own paper.

## 2026-08-13: data/xue_lu_zenodo/INVENTORY.md, dt / duration / units flags

The inventory's "262,800 rows = exactly 30 years hourly" guess (flagged VERIFY-before-use, never
cleared) is superseded: the F-016 calibration used the native sampling dt = 6000 s recorded in its
protocol, under which 262,800 rows span ~50 years, and an independent spectral check places the
principal constituents at their correct periods (O1/P1/K1/N2/M2/S2/K2) under dt = 6000 s and at
non-physical periods under dt = 3600 s. The single-column values (range ~ +/-46, sd ~14) are
consistent with nanostrain rather than Pa; the F-016 protocol's E = 75 GPa / nu = 0.25 conversion
(shear modulus 30 GPa, matching Lu et al. section 2.3) treated them as strain. INVENTORY.md is
annotated accordingly (its original guess left in place, marked superseded).


## 2026-09-02: 9,131-day epoch defect in four catalogue loaders

**What was wrong.** `exp_highn.load_zenodo`, `exp_mass_screen.load_region_full`,
`exp_fluid_driven.load` and `exp_nearcritical.load` returned event times as days since
`exp_world_harmonics.SPAN_START` (1995-01-01). The shared feature builders
(`exp_world_harmonics.features`, `exp_mass_screen.raw_series`) convert with
`jd = t + UNIX_EPOCH_JD`, which assumes days since 1970-01-01. Every tidal and lunar
covariate in the arms fed by those four loaders was therefore evaluated 9,131 days
(25.0 yr) before the event it belonged to. The dwell nulls used the same wrong base, so
the arms remained internally consistent and returned honest-looking nulls against a
scrambled lunar covariate.

**What it does and does not void.** It voids the TIDAL/LUNAR conclusions of P-1.1
(`exp_highn`), P-1.3 (`exp_fluid_driven`), P-1.4 (`exp_learned`,
`exp_learned_sensitivity`), `exp_nearcritical` and the 3,000-cell `exp_mass_screen`
as first computed. It does NOT void the OBSERVER conclusions: 9,131 days is an integer
number of days and 25.0 tropical years to within 0.05 d, so local solar hour and season
are preserved. `exp_diurnal_discriminator`'s OBSERVER verdict is confirmed unchanged on
rerun, and P-1.1's family p = 0.0003 at `sun_hourangle.R1` is byte-for-byte the same
before and after - it was and remains the day/night detection artifact.

**What changed on rerun.** `exp_learned` P-1.4 stays NULL and gets stronger
(dAUC +0.00155, p = 0.34 -> -0.00148, p = 0.83), with the design's sensitivity floor
unchanged at eps = 0.05. `exp_nearcritical` stays null (p 0.80 -> 0.62).
`exp_mass_screen`'s family p moves 0.4034 -> 0.3234 with zero candidates over the null
95th, and **the pre-specified holdout candidate `shallow_lt70km | areal.bot_decile`
(pooled z = -3.70, p = 0.064) VANISHES** - it drops out of the top 15 pooled keys
entirely. It was an artifact of the wrong epoch and the holdout must not be spent on it.
`exp_fluid_driven`'s declared 7/7 test stays null but reverses sign (p 0.814 -> 0.462),
while two of its secondary readings that were null now sit at p = 0.02 (primary) and
p = 0.0007 (full-sequence-null secondary). Those are new exploration findings, not
claims.

**What prevents recurrence.** `tests/test_epoch_consistency.py` states the invariant
once - every loader feeding `raw_series`/`features` returns days since 1970-01-01Z - and
asserts it for every loader against an independently computed Julian date, plus a
cross-loader agreement check. `exp_mass_screen.assert_epoch` prints and hard-checks the
reconstructed first-event date at the top of every affected arm's run.

First noticed by the Phase-2 (arm D) worker while auditing its own epoch handling; confirmed by the supervisor from a code read and a direct check (circular correlation 0.28 between the lunar phase used and the true lunar phase over the 33,293 QTM exploration events); fixed and re-run 2026-09-02.
Superseded outputs preserved as `results_*_epochbug_20260822.json`.
