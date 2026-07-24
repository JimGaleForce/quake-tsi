# Frozen Protocol - Long Valley TSI Out-of-Sample Prediction Test (EQ-22)

**Status: FROZEN 2026-07-21 (America/Los_Angeles), before any Long Valley catalog data is
downloaded or inspected.** Hash recorded in download_log.md.

## 1. The pre-registered prediction

Our SoCal geodetic map (tsi_map.py, NGL MIDAS, maps/socal_strain_tsi.png) identifies Long Valley
caldera as the strongest positive TSI anomaly in the mapped region (TSI ≈ +6.3 at the caldera vs
regional median ≈ −0.11), while Coso - the one site with independently established tidal
modulation (Lu et al. 2025; Wang et al. 2022) - has TSI ≈ −0.43.

**H-LV (if the TSI ratio has predictive merit):** declustered background seismicity in the Long
Valley box exhibits statistically significant tidal phase modulation (95% level against
orientation/time-preserving synthetic catalogs), detected by the SAME statistic that detects the
Coso positive control.

**Interpretation commitments (written before data):**
- Modulation at Long Valley → TSI survives as a prospecting variable; proceed to the full
  label-map correlation test and nominate further high-TSI candidates.
- Null at Long Valley (with the statistic demonstrably able to detect Coso) + modulation at
  low-TSI Coso → the TSI ratio is retired as a susceptibility proxy; the susceptibility hunt
  moves to effective-stress proxies (heat flow, b-value, swarm fraction, vp/vs).
- If NO statistic variant detects the Coso positive control, the Long Valley test is NOT
  interpretable as a TSI test and is reported only as methodology work.

## 2. Data

- Catalog: NCEDC (NCSS) via FDSN event service, box lat 37.40–38.10, lon −119.20 to −118.60,
  1984-01-01 to 2026-06-30, M ≥ 1.0, depth ≤ 20 km, chunked per year, raw responses preserved
  under data/raw/long_valley/.
- Era restriction for analysis: 1985+ (post-network buildout); Mc verified by
  Gutenberg-Richter fit before use; analysis floor = max(1.0, Mc).
- Focal mechanisms (if obtainable from NCEDC/NCSS): matched by event ID for the per-event
  variant. If unavailable, the fault-geometry variant uses mapped Long Valley normal faults
  (USGS Qfaults) and caldera ring-fault orientations; this limitation is recorded.

## 3. Declustering

Window-based (Gardner-Knopoff, same windows as protocol_params.json): remove every event inside
the space-time window of any larger event; the remainder is the background set. (Lu et al. used
ETAS declustering; GK-windowing is our pre-declared variant - coarser but deterministic. If the
background set shows modulation, an ETAS-declustered robustness check follows before claiming.)

## 4. Detection statistic

Primary: the statistic variant that SUCCEEDS on the Coso positive control (per results_coso_fm.json
- per-event focal-mechanism resolution if that is what works, else whichever fixed-orientation
variant does), applied to Long Valley UNCHANGED except for location-specific inputs.
Tidal forcing at Long Valley: astropy tidal-potential volumetric strain series computed at the
caldera (37.68N, −118.90E) at 600-s sampling (body tide; ocean loading negligible this far inland
relative to SoCal coast - recorded as an assumption), plus fault-resolved projections via the
calc_stress.py rotation for the mechanism/geometry variants using the same elastic constants.
Secondary (always reported): volumetric-phase modulation, regardless of primary outcome.
Significance: 3,000 synthetic catalogs (random times over the analyzed era, orientations
preserved where applicable), 97.5th percentile threshold, empirical p reported.

## 5. Multiplicity

One primary endpoint (H-LV via the Coso-validated statistic). Everything else exploratory.
This protocol is the SECOND site-level test of the TSI map (Coso, retrospectively, is the first);
any claims about TSI's overall merit require the full label-map correlation (planned separately).

## Amendments

### Amendment 1 - 2026-07-21 (declustering method, documented before any interpreted result)

The §3 GK-window declustering is inapplicable at this catalog's magnitudes: T(M=1.2) ≈ 598 days,
so windows cascade and only 81 of 61,297 events survived - no usable background set (this outcome
is recorded in results_long_valley.json as the failed first pass; its modulation numbers are void
for the primary/secondary endpoints). Replacement, pre-declared here before rerun:
**Zaliapin & Ben-Zion (2013) nearest-neighbor declustering** with standard parameters
b = 1.0, d_f = 1.6, threshold log10(η0) = −5.0 (fixed, not tuned); background = events whose
nearest-neighbor proximity η exceeds η0 (or with no qualifying parent). ETAS-declustered
robustness check still applies per §3 if a positive is claimed. No other parameter changes.
