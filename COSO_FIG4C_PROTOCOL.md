# Coso Figure 4c exact reproduction — frozen protocol

Frozen 2026-08-06, BEFORE any analysis with these parameters was run. SHA-256 hash recorded in
`download_log.md` (Amendment 2) at freeze time.

## Source of parameters

Weifan Lu, email 2026-08-06 (archived in `../email.txt`), answering the panel reply's blocking
question: which exact configuration produces the significant Coso tidal-modulation result
(Lu et al. 2025 Figure 4c).

## Weifan's specification (verbatim intent)

- The ONLY significant Coso result is Figure 4c: tidal modulation computed with focal mechanisms.
- Spatial bin: single 0.4° × 0.4° bin centered at (−117.8, 36.4);
  latitude range [36.2, 36.6], longitude range [−118.0, −117.6].
- Catalog: SCSN, magnitude ≥ 1.5.
- Stress component: **shear stress** (τ), resolved on the fault plane chosen according to the
  focal mechanism.
- Focal mechanism catalog: Yang-Hauksson-Shearer, the version at
  https://scedc.caltech.edu/data/alt-2011-yang-hauksson-shearer.html — one nodal plane per
  solution, so no nodal-plane choice is made.
- Robustness: their Supplementary Figs. 4(a–d) vary the bin size/extent; result stable.

## Divergence from our prior null (results_coso_fm.json)

Our earlier Coso tests used `COSO_BOX = lat (35.60, 36.25), lon (−118.05, −117.60)` — the Coso
geothermal field proper. Weifan's Fig. 4c bin lies almost entirely NORTH of that box (overlap:
lat 36.20–36.25 sliver only). The prior null therefore does not test their claim.

## Test plan (frozen)

1. Script: `coso_fig4c_test.py`, derived from `coso_fm_test.py` (identical stress/phase/synthetic
   machinery — Hooke E=75 GPa ν=0.25, CubicSpline upsampling, per-event FM resolution with
   orientation groups rounded to 15°, phase from find_peaks on the upsampled τ series).
2. Primary test (the reproduction):
   - Catalog: `SCSN_decluster_m1.5.txt` (their Zenodo declustered SCSN, M≥1.5 — the catalog
     their modulation analysis uses).
   - Events: lat ∈ [36.2, 36.6], lon ∈ [−118.0, −117.6].
   - FM match: by SCSN event id against YSH_2010.hash + sc2011–2018 annual files
     (same products the alt-2011 page links; already on disk in `data/scedc_fm/`).
   - Component: τ (shear on the FM plane). σn reported as secondary for completeness only.
   - Statistic: Pm/P0 sinusoid amplitude over 36 phase bins; significance vs 3,000 synthetic
     catalogs that keep the orientation set and randomize event times (Lu et al. §2.5 style,
     no analytic bias correction). Seed 20260806.
3. Secondary (pre-declared sensitivity, not part of the pass/fail call):
   - Same bin on the QTM declustered catalog (M≥0.1 file, cut to M≥1.5 and native Mc).
4. Interpretation rule (frozen): the reproduction "succeeds" iff the primary τ test is
   significant at 95% (empirical p < 0.05 one-sided vs synthetics) with phase φ near their
   reported preferred phase. A null here, with their exact box/catalog/component, would be a
   genuine failure-to-reproduce to raise with Weifan (likely suspects: their 6000-s tidal
   series sampling vs our upsampling, FM-catalog vintage, or declustering differences).

## Outputs

`results_coso_fig4c.json` + console log. Failure or success, the result is committed.
