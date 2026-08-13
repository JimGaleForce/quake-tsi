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
