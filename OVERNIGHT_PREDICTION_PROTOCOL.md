# Overnight prediction experiments — frozen protocol

Frozen 2026-08-09 (early AM), BEFORE any analysis of the test window under these definitions.
SHA-256 recorded in download_log.md at freeze time. Author: Claude session under Jim's
directive ("find a way to predict events... pick a date in the past to separate learn from
test"); Jim asleep — experiments run autonomously, results reported win or lose.

## Common definitions

- Catalog: SCSN declustered M>=1.5 (Zenodo `SCSN_decluster_m1.5.txt`, 1981–2018).
- **Split: train = events before 2010-01-01 UTC; test = 2010-01-01 through catalog end.**
  Chosen from span/power considerations only; no test-window statistic was computed before
  this freeze.
- Phases: per-event focal-mechanism shear-stress (τ) tidal phase via the validated
  `coso_fm_test.py` machinery (YHS/SCEDC FM matched by event id; Hooke E=75 GPa ν=0.25;
  CubicSpline ×10 upsampling; trough–peak–trough phase interpolation). Known limitation
  (6000-s native sampling) applies equally to train and test — cannot create spurious skill.
- Spatial bins: 0.4° × 0.4°, grid aligned to integer multiples of 0.4° from (−122.0, 31.5),
  SoCal box lat [31.5, 38.0], lon [−122.0, −113.5] (matches Lu et al.'s bin scale).
- Null model everywhere: circular time-shift of the tidal-series↔event-time alignment
  (uniform random shifts of 2–2000 days applied to event times before phase lookup),
  preserving all event clustering and series sampling artifacts. 3,000 shifts unless stated.
  Empirical one-sided p.

## EXP-A (confirmatory headline): out-of-sample tidal-phase timing skill

1. TRAIN: per-bin modulation on train events only (FM-matched, n_train >= 100). A bin is
   **selected** if empirical p < 0.05 (orientation-preserving time-randomized synthetics,
   1,000). Record per-bin preferred phase φ_b and amplitude a_b (sinusoid fit as in
   coso_fm_test.fit).
2. TEST (single unblinding): for test-window FM-matched events in selected bins, statistic
   S_b = mean_i cos(φ_i − φ_b), pooled S = Σ_b n_b S_b / Σ_b n_b. Significance vs the
   circular-shift null (φ_b FIXED from train; only test phases re-assigned under shift).
3. Also report: per-bin S_b, n_b; pooled log-likelihood gain (bits/event) of
   rate ∝ 1 + a_b·cos(φ − φ_b) (a_b clipped to [0, 0.9]) vs uniform.
4. **Anti-leak control**: identical scoring on train-null bins (train p > 0.5, n_train >= 100)
   — expected null; a "significant" result there indicates a pipeline artifact and voids A.
5. Success rule (frozen): pooled test p < 0.05 AND control non-significant.

## EXP-B (exploratory, split-clean): spatial transfer of susceptibility

Features per bin computed from TRAIN period only: event rate, b-value (MLE, Mc=1.5), median
depth, swarm fraction (fraction of events with a neighbor within 2 km & 3 days that is not
its own repeat), max shear strain rate at bin center (data/socal_strain_grid.npz — geodetic,
time-independent). Label: train-period modulation amplitude Pm/P0 (n >= 100 bins). Model:
rank correlation of each feature with label + leave-one-longitude-stripe-out logistic
prediction of "significant" bins, scored by AUC. Exploratory: features chosen this session.

## EXP-C (exploratory): time-varying susceptibility as a stress gauge

For the 2 highest-FM-count bins in the full-period scan: sliding 3-yr windows stepped 1 yr
(train+test span allowed — this is exploratory, not a forecast claim): window Pm/P0 and
preferred phase vs (a) next-year event count in bin, (b) occurrence of M >= 4 in bin within
the following year. Report correlations + plot data; no significance claims beyond
descriptive (multiple-comparison burden acknowledged).

## EXP-D (forward statement — only from what survives)

For bins passing EXP-A's success rule: using the bin's dominant train-period FM orientation
group, compute the τ phase series 2026-08-09 → 2026-10-08 and the predicted rate factor
1 + a_b·cos(φ(t) − φ_b). Deliverable: list of the top favorable-phase windows per bin with
the honest gain magnitude (expected ~1.1–1.6× background rate — a timing modulation of small
events, NOT an event prediction). If no bin passes EXP-A, EXP-D reports "no validated basis
for a forward statement" — that is the required output in that case.

## Discipline

- This file is hashed before EXP-A step 2 (the unblinding) runs. Train-side computation may
  iterate freely; the test window is scored ONCE per the rules above.
- All results (including nulls and the control) are committed to the public repo.
- Scripts: exp_a_phase_skill.py, exp_b_transfer.py, exp_c_susceptibility_drift.py,
  exp_d_forward.py; results: results_exp_{a,b,c,d}.json.
