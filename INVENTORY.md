# DATA INVENTORY — a living document

**Purpose.** So that no future session re-derives what we hold, and so that nobody again says
"we need to download X" when X has been on disk for a month. Both of those happened on
2026-08-22: the supervisor told Jim that the N regime change required a QTM download and that
ocean loading required FES2014/TPXO registration. **Both were already on disk.** That mistake
cost real time and is the reason this file exists.

**Rules for this file.**
1. Every entry states row counts, span and units as MEASURED, not as remembered.
2. When a dataset is first USED by an arm, add the arm's filename to its "used by" line.
3. When a dataset is found to have a defect, record it here, not only in the ledger.
4. `data/` is gitignored and must never be committed. This file is the map, not the territory.
5. Anything added to `data/` gets an entry here in the same commit that adds it.

Last measured: **2026-08-22**; §1f added and used-by lines refreshed **2026-09-02**. Totals below are raw file contents, before any magnitude
floor, exploration split or declustering.

---

## 1. EARTHQUAKE CATALOGUES

### 1a. The high-N Southern California set — `data/xue_lu_zenodo/`

**This is the program's single largest untapped asset.** Lu, Xue, Yue, Zhuang & Zhao (2025),
JGR Solid Earth, doi:10.1029/2025JB032249; data at doi:10.5281/zenodo.18491845, CC-BY 4.0.
Whitespace-delimited, columns `year mo dy hr mn sec lat lon depth_km mag event_id`.

| file | rows | span | notes |
|---|---|---|---|
| `QTM_12dev.txt` | 898,597 | 2008–2017 | template-matched, M −1.95 to 7.20 |
| `QTM_decluster_m0.1.txt` | 45,069 | 2008–2017 | **ETAS-declustered, M ≥ 0.1** |
| `SCSN_original_catalog.txt` | 634,252 | 1981–2018 | relocated SCSN (Hauksson et al.) |
| `SCSN_decluster_m1.5.txt` | 50,313 | 1981–2018 | **ETAS-declustered, M ≥ 1.5** |

QTM magnitude ladder: M≥0.0 → 508,093; M≥0.5 → 272,746; M≥1.0 → 127,789; M≥1.5 → 55,241;
M≥2.0 → 23,924; M≥2.5 → 9,210.

**Why it matters.** The two declustered sets together are ~95,000 independent events against
the 7,139 every world arm has used. Detection floor moves from ≈5% to ≈1.4% relative
modulation, which is the first time this program would be inside the range where the
literature's claimed effect (~1%) actually lives. Every null we have produced so far is
consistent with a real effect simply smaller than the instrument.

**Ocean-loaded tidal stress, precomputed** — also here, also previously believed missing.
Computed with SPOTL (solid earth + **TPXO ocean loading**) at 117°W, 35°N; the paper argues
principal tidal stresses are near-uniform across the study area.

| file | rows | content |
|---|---|---|
| `Tidal_N_0.txt` | 262,800 | normal stress, orientation 0° |
| `Tidal_N_90.txt` | 262,800 | normal stress, orientation 90° |
| `Tidal_S_0.txt` | 262,800 | shear stress |
| `Tidal_Vol.txt` | 262,800 | volumetric stress |

**DATED CORRECTION (2026-08-13, carried over from the dataset's own INVENTORY):** native
`dt = 6000 s` (≈50 years of coverage), NOT hourly/30 years as first guessed; values are
consistent with **nanostrain, not Pa**. Verify against `calc_stress.py` before use.

**Fault model:** `CFM5.3_traces.lonLat`, `CFM5.3_traceslonLat_fault_geometry.txt` (SCEC
Community Fault Model 5.3; geometry segments headed `> strike dip rake`), `CFM5.3_Metadata.xlsx`.
**Code:** `calc_stress.py` — their rotation onto fault planes (X=N, Y=E, Z=up).

**FROZEN USAGE RULES that travel with this dataset** (from its own inventory, and they are
binding):
1. Freeze the analysis plan in a protocol file BEFORE running anything against these catalogues.
2. Positive control: their reported geothermal-region modulation (Coso). Negative control:
   whole-catalogue null.
3. **Unit of inference: SEQUENCE, not event.**

> Rule 3 was violated on 2026-08-22 by `exp_nearcritical.py` v1, which used a per-event null on
> clustered aftershocks and returned max |z| = 26.27 — pure dependence. With a sequence-level
> null the same arm returns 2.51, p = 0.80. The rule is not bureaucratic; it is the difference
> between a result and an artifact.

> **DEFECT, dated 2026-09-02 (rule 3).** Until this date `exp_highn.load_zenodo` returned days since
> 1995-01-01 while the shared feature builders assumed days since 1970-01-01, so every tidal/lunar
> feature in the arms below was evaluated 9,131 d before its event. Fixed at the source, regression-
> tested (`tests/test_epoch_consistency.py`), all affected arms re-run; see `CORRECTIONS.md`
> 2026-09-02. The observer (solar-hour) conclusions were unaffected; the tidal ones were re-derived.

*Used by:* `exp_highn.py` (P-1.1), `exp_learned.py` + `exp_learned_sensitivity.py` (P-1.4),
`exp_learned_ext.py` (Phase 2), `exp_diurnal_discriminator.py`, `exp_bvalue_skill.py` (K-405/P-2.2),
`exp_neural_tpp.py` (neural TPP vs ETAS), `exp_productivity_assim.py` (K-436), `exp_h_etas.py` (B-2).

### 1b. Global regional catalogues — `data/comcat_world/` (13 files, 71,803 events at M≥4.5)

ComCat, 1995-01-01 to 2026-08-09, M ≥ 4.5. Standard ComCat CSV header.

| region | n | max M | | region | n | max M |
|---|---|---|---|---|---|---|
| Indonesia | 24,437 | 9.1 | | Alaska-Aleutians | 3,714 | 8.2 |
| Japan | 12,765 | 9.1 | | Mexico | 2,267 | 8.2 |
| Chile | 9,499 | 8.8 | | Greece-Aegean | 1,970 | 7.0 |
| Philippines | 8,833 | 7.8 | | Iran | 1,779 | 7.7 |
| Himalaya | 3,937 | 7.8 | | Turkey | 1,433 | 7.8 |
| | | | | Iceland | 475 | 6.5 |
| | | | | California | 441 | 7.2 |
| | | | | Caribbean | 253 | 7.2 |

**Note the headroom nobody has used:** every arm to date sets `MAG_MIN = 5.0` and takes a 70%
exploration split, ending at 7,139 declustered events. The files hold M ≥ 4.5. Dropping to the
catalogue's own floor is a free ~1.6x in N.

*Used by:* `exp_world_harmonics.py`, `exp_world_faultrelative.py`, `exp_world_amplitude.py`,
`exp_mass_screen.py`, `exp_dispersion_gate.py`, `exp_dispersion_dof.py`, `exp_null_calibration.py`.

### 1c. SoCal M≥2.5 — `data/comcat_socal_m25.csv`

18,389 events, 2010-01-01 to 2026, M 2.5–7.2. ComCat CSV.
*Used by:* `exp_nearcritical.py` (15,154 in the exploration split, 972 sequences).

### 1d. Geothermal / volcanic / fluid-driven regions — `data/k034/` (14 files, 137,480 events)

1985–2022, M ≥ 1.5. **These are the low-effective-stress, fluid-driven populations where
tidal sensitivity should be largest if it exists anywhere** — the regular-seismicity analogue
of the tremor/LFE populations that are known to be strongly tidally modulated.

| region | n | | region | n | | region | n |
|---|---|---|---|---|---|---|---|
| geysers | 34,546 | | salton_brawley | 13,998 | | yellowstone | 6,862 |
| long_valley | 29,241 | | san_jacinto | 11,310 | | parkfield | 3,308 |
| coso | 22,320 | | mendocino | 11,118 | | wasatch_slc | 1,016 |
| | | | mono_west_nv_mina | 992 | | lassen | 944 |
| | | | cedar_city_ut | 789 | | smith_valley_nv | 662 |
| | | | little_skull_mtn | 374 | | | |

`manifest.json` and `download.log` present. **Coso is the declared positive control** for the
Lu/Xue cross-test.
*Used by:* `results_coso_*.json` (earlier program phases), `exp_fluid_driven.py` (P-1.3),
`exp_k034_landers_control.py` (K-034), `exp_dyntrig_skill.py` (K-038 arm A, 2026-09-02).

### 1e. Other regional raw catalogues — `data/raw/` (5 regions, 191 year-files)

`apennines` (1990–2026, .txt), `iceland` (1990–2026, .csv), `long_valley` (1984–2026, .txt),
`parkfield_control` (1990–2026, .csv), `taupo` (1990–2026, .txt). Per-year files; formats
differ by region and need a loader each.

### 1f. Global trigger catalogue — `data/global_m55/` (1 file, 18,769 events at M≥5.5)

USGS ComCat FDSN event query, worldwide, **M ≥ 5.5, 1985-01-01 → 2023-01-01**, standard
ComCat CSV header. Downloaded 2026-09-02 by `exp_dyntrig_skill.download_triggers()`.

| file | rows | span | notes |
|---|---|---|---|
| `global_m55.csv` | **18,769** | 1985-01-02T05:32:49.140Z – 2022-12-28T16:34:19.271Z | M 5.5 – 9.1, `orderby=time-asc` |
| `download_log.jsonl` | 38 records | one per calendar year | URL, HTTP status, row count, UTC timestamp |
| `manifest.json` | — | — | minmag, span, page count, pre/post-dedup rows, sha256 |

`sha256(global_m55.csv) = c6ae76d8fc3dfabfbc0c238f2e363c98c384609c64a3135123dd3e2ea441827e`
(first 16: `c6ae76d8fc3dfabf`). 38 pages, one per year; **no page returned 0 rows and no page
hit the 20,000 cap** — both conditions raise rather than pass silently. 18,769 rows before
de-duplication on `id`, 18,769 after (no duplicates across year boundaries).

**Why it exists.** This is the trigger side of K-038: the "dynamic-stress weather" a target
region is exposed to from the whole planet, rather than from one named mainshock. The target
side is §1d (`data/k034/`, the 14 fluid/geothermal cells). Amplitude at a target is computed
with the **frozen K-034 axis** (van der Elst & Brodsky 2010 eq. 6, `K034_SEALED_LITERATURE.md`)
— it is not refit here and must not be refit by anything downstream.

**Defects and cautions, recorded per rule 3.**
1. **No `Ms` column.** ComCat reports a mixed magnitude type (`mw`/`mww`/`mb` …) in `mag`.
   The K-034 amplitude axis is written in `Ms`. Every use of this file therefore rides on
   K-034's already-carried "Mw-substituted axis", not on its primary Ms axis.
2. **Magnitude-scale evolution 1985→2022** is not corrected. Pre-1990 events are more likely
   to carry `mb`/`ms` than `mww`. This is a completeness-style exposure on the trigger side
   and no arm has yet measured it.
3. **18 events lie inside one of the 14 K-034 target boxes** and must be dropped by any arm
   that treats these as *remote* triggers (`exp_dyntrig_skill.py` drops them globally).
4. Sub-M5.5 teleseisms are absent by construction; an arm needing the low-amplitude tail of
   the exposure distribution needs a different download.

**Frozen usage rules that travel with this dataset.**
- Trigger exclusions used by the first arm, and the ones any comparable arm should reuse:
  ≥ 300 km from the target **box** (not the centroid), ≥ 2 rupture lengths, and not inside
  any K-034 box. Rupture length from `log10 L[km] = -2.44 + 0.59 M` (Wells & Coppersmith 1994
  RLD, all slip types) — **UNVERIFIED against the paper in session**; it sets an exclusion
  gate only, never a statistic.
- Surface-wave arrival at the target is `t_trigger + r / 3.5 km s⁻¹`, K-034's phase velocity.

*Used by:* `exp_dyntrig_skill.py` (K-038 arm A, 2026-09-02).

---

## 2. FORCING AND ENVIRONMENTAL SERIES — `data/spaceweather/` (35 files)

Named for its origin, but it holds three quite different things and **two of them are not
space weather at all**.

| file | content | span | cadence |
|---|---|---|---|
| `eopc04_14_IAU2000.62-now.csv` | **IERS Earth orientation: polar motion (x_pole, y_pole), UT1-UTC, and LOD** with per-field sigmas | 1962 (MJD 37665) – now | daily |
| `Kp_ap_since_1932.txt` | Kp and ap geomagnetic indices, GFZ Potsdam, CC-BY 4.0 | 1932 – now | 3-hourly |
| `omni2_YYYY.dat` (32 files) | OMNI2 hourly solar wind: IMF, plasma, indices | 1995 – 2026 | hourly, 8,760 rows/yr |

**LOD and polar motion are the K-428 kill-or-confirm dataset** (Bendick & Bilham decadal
claim) and were sitting in a directory named "spaceweather". Note also that `download_log.jsonl`
records provenance.

**Honest prior:** the physical coupling from geomagnetic activity to fault stress is orders of
magnitude below tidal. The defensible use of Kp/OMNI is as a **matched placebo arm** that
calibrates how often the pipeline manufactures survivors on a property class with real temporal
structure and near-zero physical prior (Kepler K-427). LOD/polar motion is a different case:
it is a genuine (if small) mechanical channel and has a specific published claim to test.

*Used by:* `engine_ext_forcing.py` -> `exp_learned_ext.py` + `exp_learned_ext_sensitivity.py`
(Phase 2 step 1 and the P-2.4 placebo reading, 2026-09-02: all three blocks null; survivor rate 0.052/feature).

---

## 3. GEODESY AND STRAIN

| path | content | size |
|---|---|---|
| `data/ngl/midas.IGS14.txt` | Nevada Geodetic Lab MIDAS velocities, **20,168 stations**, IGS14 frame; per-station east/north/up rates with sigmas and span | 4.8 MB |
| `data/kreemer_young/GPS_table.txt` | GPS velocity table, columns `lon lat ve vn sdve sdvn ve_p vn_p stat source` | 188 KB |
| `data/socal_strain_grid.npz` | 66 × 86 grid; keys `lats`, `lons`, `dilatation_nstrain_yr`, `max_shear_nstrain_yr`, `tsi` | 108 KB |

**Caveat already established as B-5:** the dilatation component of geodetic strain carries ±2×
measurement uncertainty across velocity solutions; **shear is robust**. Any covariate built on
dilatation inherits that uncertainty; build on shear by preference.

*Used by:* `results_strain_comparison.json`, `results_exp_j/k.json` (stress-ledger, B-4).

---

## 4. FAULT GEOMETRY AND MECHANISM

| path | content |
|---|---|
| `data/slab2/Slab2Distribute_Mar2018/` | Slab2 subduction interface geometry (Hayes et al. 2018, doi:10.5066/F7PV6JNV); plain-text `.xyz` grids, plus `_cache/` |
| `data/scedc_fm/` | SCEDC focal mechanisms: `sc2011.focmec` … `sc2017.focmec`, plus `YSH_2010.hash` (25 MB total) |
| `data/xue_lu_zenodo/CFM5.3_*` | SCEC Community Fault Model 5.3 traces and geometry |

**Why `scedc_fm` matters and is unused:** it gives each event its own strike/dip/rake, which
means resolved Coulomb stress can be computed **on the plane that actually slipped** rather
than on a regional guess. Every fault-relative arm so far has used Slab2 interface geometry or
an assumed thrust. Per-event mechanisms remove that approximation entirely for SoCal.

*Used by:* Slab2 → `engine/slab2.py`, `exp_world_faultrelative.py`. `scedc_fm` and CFM5.3 →
nothing yet.

---

## 5. PRECOMPUTED TIDAL / EPHEMERIS GRIDS

| path | keys | shape |
|---|---|---|
| `data/lunar_grid_1980_2027.npz` | `t_unix`, `elong`, `ndist` | 412,009 samples |
| `data/lv_tidal_vol.npz` | `t_unix`, `eps` | 2,182,320 samples (Long Valley volumetric strain) |
| `data/xue_lu_derived/lunar_grid.npz` | — | derived |
| `data/xue_lu_derived/QTM_aftershocks.csv`, `SCSN_aftershocks.csv` | derived aftershock selections | — |

---

## 6. LOGS AND PROVENANCE

`data/comcat_log_20260809.txt`, `data/comcat_world_log.txt`,
`data/retrieval_log_2026-07-20T230712-0700.tsv`, `data/spaceweather/download_log.jsonl`,
`data/k034/manifest.json` + `download.log`, `data/xue_lu_zenodo/INVENTORY.md`.

---

## 7. WHAT WE DO **NOT** HAVE (so nobody assumes we do)

- **Ocean tide gauge records** (UHSLC/GESLA) — needed for the storm-surge step-loading design.
  Small download. *Not present.*
- **Global ocean load model** (FES2014/TPXO) as a computable field. We have the Lu/Xue SoCal
  point series only; there is no global loadtide capability.
- **Repeating-earthquake catalogues** (Parkfield HRSN, Japan/Igarashi). *Not present.*
- **Slow-slip event catalogues** (Cascadia, Guerrero, Boso, Hikurangi). *Not present.*
- **GCMT global focal mechanisms** — needed for outer-rise selection. *Not present.*
- **ERA5 / ERA5-Land / GRACE** hydrology and pressure. *Not present.*
- **Station-noise PSDs** for a measured Mc(x,t) field. *Not present.*
- **Groningen / Oklahoma induced catalogues + pressure histories.** *Not present.*
- **ISC-GEM** global M≥7 back to 1900 (for the LOD decadal test we DO have the LOD half of).
  *Not present.*

---

## 8. STANDING CAPACITY SUMMARY

| what we can test today | best N | floor (4.29/√N) |
|---|---|---|
| SoCal declustered, QTM + SCSN | ~95,000 | **≈1.4%** |
| SoCal declustered, QTM alone M≥0.1 | 45,069 | ≈2.0% |
| Fluid-driven regions (k034) | 137,480 raw | sequence-limited |
| World M≥4.5 declustered | ~11,000 est. | ≈4.1% |
| World M≥5.0 declustered (what we HAVE used) | 7,139 | ≈5.1% |

Effect size claimed in the literature for ordinary seismicity: **~1%**. Only the first row is
close, and only the first row has ocean-loaded stresses to go with it.
