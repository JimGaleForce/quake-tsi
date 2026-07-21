# Frozen Protocol — Blind Cross-Test on the Lu/Xue Southern California Dataset (EQ-18)

**Status: FROZEN 2026-07-21 (America/Los_Angeles), before any analysis of the catalogs.**
Disclosure at freeze time: the only inspection of the data so far is file head-rows (first 3 lines)
and row counts, done to document formats in `data/xue_lu_zenodo/INVENTORY.md`. No event selection,
filtering, or statistics have been run. Parameters live in `xue_lu_params.json`; analysis code must
read them from there. Changes after this point require a dated amendment appended here BEFORE the
affected step runs.

## 1. Question

Does the Kosmos-r89 "predictability class" — small (M < 5.0), early (< 24 h) aftershocks on normal
faults, occurring during perigee-syzygy — show tidal-window enrichment in the Lu/Xue Southern
California catalogs (QTM, Mc 0.1, 2008–2017; relocated SCSN, Mc 1.5, 1981–2018)?

This is out-of-sample in two ways: (a) the original n=14 came almost entirely from Chalfant Valley
(37.5°N), which lies OUTSIDE the Lu et al. Southern California study region — so nothing here can
re-find the original events; (b) the QTM catalog reaches Mc 0.1, ~3 magnitude units deeper than the
original M≥4.5 catalog. The poroelastic mechanism predicts the signal should be PRESENT and
STRONGER at small magnitudes in fluid-rich (geothermal) zones — Lu et al.'s own finding of
modulation at Coso and other low-effective-stress areas makes this a fair, live test.

## 2. Event classes (faithful to r89's operational definitions, protocol v1.1.0 Amendment 1)

- **Mainshock:** any event with M ≥ 5.0 in the raw catalog (QTM_12dev / SCSN_original).
- **Candidate aftershock:** event with Mc ≤ M < 5.0 (Mc = 0.1 QTM, 1.5 SCSN), occurring
  0 < Δt ≤ 24 h after a mainshock and within 100 km horizontal (haversine); associated to the
  closest-in-distance qualifying mainshock. Each event counted once.
- **Normal-fault subset:** aftershocks within 10 km horizontal of any CFM 5.3 trace point whose
  segment rake ∈ [−135°, −45°].
- **Strike-slip negative-control subset:** same, rake within ±45° of 0° or 180° (and not qualifying
  as normal). The r89 class was mechanism-specific; comparable enrichment in the strike-slip subset
  argues artifact, not mechanism.
- Note: this replaces Gardner-Knopoff classification with r89's actual operational association
  (nearest preceding larger event within 100 km / 24 h) — which is what produced the original n=14.

## 3. Tidal windows (identical to Amendment 1)

- Syzygy: geocentric Moon–Sun elongation < 30° (new) or > 150° (full).
- Perigee: normalized lunar distance < 0.25, normalized = (r − 356,500 km)/(50,200 km), clamped.
- Perigee-syzygy: both true at origin time. Ephemeris: astropy built-in (r89 used the same; DE440
  differences are negligible at these thresholds).

## 4. Statistics

- **H1 (primary):** the normal-fault aftershock class in the QTM catalog is enriched in
  perigee-syzygy windows. Null: 10,000 surrogates re-evaluating each event's tidal mask at
  t + δ, δ ~ Uniform(0, 29.53059 d) (circular ephemeris shift; event times untouched — preserves
  all clustering exactly). ER = observed in-window count / mean surrogate in-window count;
  one-sided p = exceedance fraction (add-one). Decision: supported if ER > 1 with p < 0.05;
  strong if ER ≥ 1.5. SCSN is the secondary readout (longer era, higher Mc).
- **Sequence-level robustness (mandatory):** jackknife by mainshock — drop each mainshock's
  aftershocks and recompute ER; report range. A result driven by one sequence (Chalfant redux)
  is reported as such.
- **H2 (secondary):** for in-window normal-fault events, |dεv/dt| percentile within the 24-h
  window centered on the event (1-min resolution, Moon+Sun, εv ∝ Σ GM/r³·(3cos²θ−1), r89's
  formulation) has median > 50% (one-sided Wilcoxon). Calibration: same statistic on 1,000
  phase-shifted surrogate time sets must center on 50%.
- **Controls:** strike-slip subset through the identical pipeline (expect no comparable enrichment);
  all-mechanism class reported for completeness.
- Multiplicity: one primary endpoint (H1/QTM/normal). Everything else labeled secondary.

## 5. Interpretation commitments (written before running)

- **H1 null in both catalogs:** the r89 class does not generalize beyond Chalfant Valley at any
  magnitude — the California section of the paper is downgraded to a single-swarm case study, and
  we report this back to Bürgmann/Vidale/Xue as the blind-test outcome.
- **H1 supported, control clean:** first out-of-sample support for the predictability class;
  immediately map WHERE the in-window events sit vs Lu et al.'s low-effective-stress regions
  (Coso etc.) and vs TSI (EQ-20), and report both ways.
- **H1 supported but strike-slip control also enriched:** artifact suspected (window construction,
  catalog structure); investigate before reporting any positive.

## Amendments

*(none yet)*
