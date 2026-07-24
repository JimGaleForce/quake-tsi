# Kosmos "quakes4" Methodology Playbook (living doc)

Purpose: a self-contained, re-runnable description of the analysis pipeline the Kosmos AI run
(`quakes4`, Edison playground run 7e2d1b15-4b0b-48db-adaa-2a8dfaea4150, Nov 2025) used to produce
the paper's California "tidally sensitive aftershock" result - so it can be repeated on any catalog
without Kosmos. Local reference implementation: `xue_lu_crosstest.py` (+ `download_catalogs.py`).
Update this file whenever a step is clarified from the playground or changed in our reimplementation.

Provenance map (playground subtask IDs):
- r6: Gardner-Knopoff declustering of the global M≥4.5 catalog
- r89 (a.k.a. "Task 88" in Jim's email): definitive identification of the n=14 sensitive events -
  THE canonical selection pipeline, restated below
- r82/r83/r84: spatial logistic-regression / Random-Forest models (source of the max_shear_rate
  vs TSI finding, see EQ-17)
- r85: 3D tidal stress / Coulomb comparison (source of the paper's "100% of 13, Wilcoxon p=0.0002")
- r10/r14/r68: California strain-rate grid & TSI
- r56/r20/r28: Tidal-ETAS variants (global null / no ΔLL gain)

## Stage 0 - Base catalog

- USGS ComCat, global, 1900–2025, M ≥ 4.5, n = 303,702 (`quakes.csv`).
- Known quirks Kosmos recorded: some negative depths; timestamp format inconsistencies
  ("+00:00" vs "T...Z") - robust join key = `time[:19]` + rounded lat/lon; duplicate location
  rows in derived files need dedup.

## Stage 1 - Declustering / event classes (r6)

- Gardner-Knopoff (1974) windows: T(M) = 10^(0.032M+2.7389) days, L(M) = 10^(0.1238M+0.983) km.
- Output classes: mainshock 11.0%, foreshock 28.6%, aftershock 47.8%, background 12.7%.
  ⚠ 28.6% foreshocks is anomalous vs standard GK (few %) - unaudited implementation quirk (EQ-2).
- **Important:** the sensitive-event selection did NOT use GK links directly. r89 re-associated
  aftershocks itself (Stage 3). GK classes only gated which events counted as "aftershock".

## Stage 2 - Lunar parameters (derived artifact `earthquakes_with_lunar_parameters.csv`)

Per event timestamp, astropy built-in ephemeris:
- **Lunar phase angle** 0–180°: geocentric Moon–Sun elongation (0 = new, 180 = full).
  Verified against the n=14 CSV: 1986-07-21 (full moon) → 174.5°; 1999-05-15 (new moon) → 6.0°.
- **Earth–Moon distance** (km) and **normalized lunar distance** = (r − 356,500)/(406,700 − 356,500),
  0 = perigee, 1 = apogee (global normalization, not per-lunation).

## Stage 3 - Sensitive-class selection (r89, "finding f6" criteria)

Applied in order to the classified catalog:
1. California box: lat 32–42, lon −125 to −114.
2. Class = aftershock, M < 5.0.
3. Aftershock→mainshock link: haversine to nearest preceding mainshock within 100 km;
   keep events with Δt < 24 h. (Custom association, not GK's internal linking.)
4. Normal-fault association: USGS Quaternary Fault Database (`Qfaults_2020_WGS84.gdb`), normal-type
   faults, event within 10 km via scipy cKDTree. (10.5% of normal faults lost to strike-notation
   parsing - recorded caveat.)
5. Perigee-syzygy: phase angle < 30° or > 150°, AND normalized distance < 0.25.

Result: n = 14 (12× Chalfant Valley 1986-07-21/22, 1× 1999 Long Valley area, 1× 2013 Canyondam).

## Stage 4 - Tidal volumetric strain rate & percentile (r89)

- εv(t) ∝ Σ_{Moon,Sun} GM/r(t)³ · (3cos²θ(t) − 1), θ = body zenith angle at the event epicenter
  (spherical trig from geocentric RA/Dec + local sidereal time - astropy AltAz was avoided).
- Time series at 1-minute resolution over the 24-hour window CENTERED on the event.
- |dεv/dt| by differencing; **percentile of the event-time value within its own 24-h window**.
  ⚠ Window-relative percentiles are NOT comparable across events for global thresholds (Kosmos's
  own warning); fine for the median-vs-random test only.
- Statistics: one-tailed Wilcoxon of percentiles vs 50% (random-timing null) → median 76.4%,
  p = 0.034; bootstrap 95% CI on median [68.6%, 91.7%].

## Stage 5 - Coulomb-vs-volumetric comparison (r85; paper's p = 0.0002)

For n = 13 events with usable geometry: percentile rank of |dεv/dt| vs percentile rank of
fault-resolved Coulomb stress RATE at origin time; Wilcoxon signed-rank on the paired difference -
volumetric-rate percentile higher in 13/13 cases. ⚠ Needs null calibration (the two percentile
constructions may not be exchangeable - EQ-3); mainshock fault geometries were hand-compiled
(a known scale-limitation Kosmos recorded).

## Stage 6 - Spatial context: strain grid & TSI (r10/r14/r68, r82/r83)

- NGL MIDAS IGS14 velocities, 1,612 stations → regularized inversion → 0.1° grid;
  dilatation rate ε̇vol = ε̇xx + ε̇yy and max shear strain rate; TSI = ε̇vol / max(ε̇shear).
- Join to events via cKDTree nearest grid cell.
- ⚠ r82/r83 logistic regressions found **max_shear_rate alone** (not TSI) is the primary
  significant spatial predictor; Random Forest (r84) did not beat logistic. The paper's TSI-first
  framing is not what the models found (EQ-17).

## Stage 7 - Global tests & ETAS benchmark (r3/r4/r15, r56 etc.)

- Global null: tidal stress distribution at event times vs 10,000 random-time catalogs
  (KS p = 0.187, MWU p = 0.675, Cohen's d = 0.008). "18,265 Pa mean tidal stress" is a full-tensor
  scalar magnitude including the solid-earth body tide - NOT the ~1 kPa fault-resolved cycle
  (source of Bürgmann's puzzlement; must be defined in Paper v2).
- Tidal-ETAS: μ(t) = μ_base(1 + β·I_tidal) with I_tidal = perigee-syzygy indicator; no significant
  ΔLL over ETAS globally; negative on held-out data.

## Our reimplementation deltas (xue_lu_crosstest.py - EQ-18)

| Step | Kosmos | Ours | Why |
|---|---|---|---|
| Aftershock gate | GK class then custom link | Direct link: mainshock M≥5.0, Δt≤24h, ≤100 km, closest | GK class adds nothing after the link; avoids the r6 quirk |
| Fault styles | Qfaults normal type | CFM 5.3 rake ∈ [−135,−45]; strike-slip control set | CFM is the SoCal community standard (matches Lu et al.); adds negative control |
| Null model | (none for enrichment per se) | 10,000 circular ephemeris shifts over the synodic month | preserves clustering exactly; the load-bearing improvement |
| Inference unit | events | + per-mainshock jackknife | Chalfant lesson (EQ-3) |
| Ephemeris | astropy built-in | same | comparability |
| Strain rate | Stage 4 formula | identical formula/resolution | faithfulness |

## Outcome log

- **2026-07-21 - EQ-18 blind cross-test (Lu/Xue SoCal catalogs): NULL.** Normal-fault early
  aftershocks show no perigee-syzygy enrichment (QTM ER=0.02 p=0.50; SCSN ER=0.05 p=0.87 - wide
  nulls, correctly) and no |dεv/dt| elevation (medians 39%/38%). Strike-slip control ER=2.1
  p=0.04 (uncorrected) → construction, not mechanism.
- **Root-cause lesson (add to any future use of Stage 3):** the 24-h aftershock window is ~1/30 of
  the syzygy cycle, so an entire early-aftershock sequence inherits its MAINSHOCK's window
  membership. Effective n = number of mainshocks, not aftershocks. Chalfant's mainshock struck at
  full moon near perigee → all 12 early aftershocks "sensitive" → the original 8.19x / p≈10⁻⁸ was
  ~1 effective observation. Any event-level statistic on early aftershocks vs a slowly-varying
  (multi-day) tidal indicator is pseudo-replicated; use mainshock-level or background-seismicity
  constructions (Lu et al. phase-histogram on declustered catalogs) instead.

## Re-run recipes

- **New region/catalog:** put a whitespace catalog (yr mo dy hr mi sec lat lon depth mag id) under
  `data/`, add an entry to `xue_lu_params.json` `catalogs`, run `python xue_lu_crosstest.py`
  then `--h2`. Freeze a protocol variant FIRST if it's a confirmatory test.
- **Original-flavor run (Qfaults, GK-gated):** swap Stage-3 fault source for Qfaults normal traces
  and gate candidates on GK aftershock class - parameters already in `protocol_params.json`.
- Lunar grid cache: `data/xue_lu_derived/lunar_grid.npz` (hourly 1980–2019); delete to rebuild,
  extend the date range in `build_lunar_grid` for other eras.
