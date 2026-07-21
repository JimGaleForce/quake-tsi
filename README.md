# Tidal Triggering Replication & Blind Cross-Test

Code, pre-registered protocols, and results for testing the claims of *"Conditional Tidal
Triggering of Microseismicity … and the Tidal Susceptibility Index (TSI)"* (Gale, 2025) —
including a blind cross-test on the independent Southern California dataset of
[Lu, Xue et al. (2025, JGR Solid Earth)](https://doi.org/10.1029/2025JB032249), out-of-sample
tests in three additional regions, and a pre-registered site-level test of the TSI prediction
at Long Valley.

**Headline outcomes** (details in the protocol files and `results_*.json`):

1. The original "8x enrichment of small early normal-fault aftershocks during perigee-syzygy"
   does **not** replicate anywhere out of sample. Root cause identified: a 24-h aftershock
   sequence inherits its mainshock's multi-day tidal-window membership as a block, so the
   effective sample size is the number of mainshocks (~1 in the original result — the 1986
   Chalfant Valley mainshock happened to strike at full moon near perigee). The strike-slip
   negative control reproduces a spurious "6.6x enrichment" by the same construction.
2. Mainshock-level timing shows no perigee-syzygy preference in Southern California, the Central
   Apennines, Taupō, or Iceland (L'Aquila 2009, Amatrice/Visso/Norcia 2016, and the large South
   Iceland events all fell outside windows).
3. A full Southern California 0.2°-bin phase-modulation scan (per-event focal-mechanism stress
   resolution, following Lu et al. §2.5) detects significant modulation in the San Jacinto/Anza
   area independently in both catalogs, sees an elevated-but-underpowered amplitude at Coso, and
   finds **no positive spatial correlation between the TSI ratio and tidal modulation**
   (Spearman ρ = −0.44 QTM / −0.03 SCSN) — retiring TSI as a susceptibility proxy while leaving
   the underlying fluid/effective-stress physics (Lu et al.; tremor; Axial Seamount) intact.
4. Long Valley caldera — the strongest positive TSI anomaly in our geodetic map — shows no
   tidal phase modulation in 19,787 declustered background events (pre-registered prediction
   test; protocol frozen and hashed before data download).

## Methodological discipline

Every confirmatory analysis has a **frozen protocol** (`PROTOCOL.md`, `XUE_LU_PROTOCOL.md`,
`LONG_VALLEY_PROTOCOL.md`) whose SHA-256 hash was recorded in `download_log.md` **before** the
corresponding data was analyzed. Deviations are dated amendments appended before the affected
step ran. Failed passes (a declustering method that collapsed, two parsing bugs) are logged, not
hidden. `KOSMOS_METHODOLOGY.md` documents the original analysis pipeline being tested, with
provenance to the public Kosmos run.

## Reproducing

Python ≥3.11 with `numpy pandas scipy astropy matplotlib`.

Data (not committed; all public):
- Lu/Xue catalogs, tidal stress series, CFM 5.3: https://doi.org/10.5281/zenodo.18491845
  → `data/xue_lu_zenodo/`
- Regional catalogs (INGV, GeoNet, USGS ComCat): `python download_catalogs.py`
- SCEDC focal mechanisms (YSH + annual): https://scedc.caltech.edu/ftp/catalogs/hauksson/Socal_focal/
  → `data/scedc_fm/`
- NGL MIDAS velocities: https://geodesy.unr.edu/velocities/midas.IGS14.txt → `data/ngl/`
- Long Valley catalog: NCEDC FDSN event service → `data/raw/long_valley/`

Order of execution mirrors the audit trail:
`download_catalogs.py` → `xue_lu_crosstest.py` (+ `--h2`) → `coso_positive_control.py` →
`coso_fault_resolved.py` → `coso_fm_test.py` → `tsi_map.py` → `three_region_test.py` →
`long_valley_test.py` → `socal_bin_scan.py`.

## Acknowledgments / data credits

Lu, Xue, Yue, Zhuang & Zhao (2025) for the openly shared dataset (CC-BY 4.0); SCEDC and the
Yang-Hauksson-Shearer focal mechanism catalog; NCEDC/NCSS; INGV; GeoNet; USGS; Nevada Geodetic
Laboratory. The original analysis under test was produced with the Kosmos AI scientist platform
(Edison Scientific); this repository's reanalysis is independent local code.
