# Temporal-pattern prediction experiments — frozen protocol (round 2)

Frozen 2026-08-09 before any test-window analysis under these definitions. SHA-256 in
download_log.md. Directive (Jim): pure-math temporal structure — specific intervals,
Fibonacci-type ratio sequences, "patterns within patterns and meta patterns outside
patterns" — mined on TRAIN, tested ONCE on TEST.

## Common

- Split: **train < 2010-01-01 UTC, test 2010-01-01 → catalog end** (same as round 1).
- Catalogs: SCSN_original_catalog.txt (full, for clustering/ETAS — NOTE its column order
  differs from declustered files; auto-detect) and SCSN_decluster_m1.5.txt (for periodicity,
  where clustering is a confound). SoCal box lat [31.5, 38.0], lon [−122.0, −113.5].
- All mined structures/parameters are frozen from TRAIN before the single TEST scoring.
- Every experiment reports its result win or lose; all committed.

## EXP-H (confirmatory): temporal ETAS vs Poisson — the bar

Magnitude floor M >= 2.5 (worker verifies train/test Mc stability; raise to 3.0 if unstable
and record). Fit temporal ETAS λ(t) = μ + Σ_{t_i<t} K·10^{α(M_i−M0)}·(t−t_i+c)^{−p} on train
by MLE. TEST: walk-forward daily-rate forecast over the test window (parameters FROZEN from
train; history includes past test events — standard walk-forward). Score: mean log2-likelihood
gain per event vs (a) stationary Poisson with train-period rate, (b) Poisson with TEST-period
rate (rate-oracle control — removes "the 2010s were quieter/busier" trivial skill).
Success (frozen): gain vs BOTH baselines > 0.5 bits/event. Expected: pass by a wide margin —
this is the validated-physics baseline everything else must beat.

## EXP-F: periodicity comb ("specific intervals / geo-influenced intervals")

Declustered catalog. Periods tested: log-spaced comb 0.25 d → 3650 d (60 periods) PLUS named
candidates {0.5 d, 1 d, 7 d, 13.66 d, 14.77 d, 27.32 d, 27.55 d, 29.53 d, 182.62 d,
365.25 d}. Statistic: Rayleigh R of event phases (t mod P). Null: 500 surrogate catalogs per
period — inhomogeneous Poisson with rate = Gaussian-kernel-smoothed observed rate, kernel
σ = 5×P (preserves slow rate structure, destroys phase locking at P). Pre-labeled artifact
periods: 0.5 d, 1 d, 7 d (detection threshold / blasting) — reported separately, never as
prediction. TRAIN: detect (p < 0.01 uncorrected). TEST: each train detection re-scored on
test events (p < 0.05 = confirmed). Prior stated in advance: 365.25 d (hydrologic loading)
is the live physical hypothesis; 1 d and 7 d should fire as artifacts (positive controls for
the method — if they don't fire on the ORIGINAL catalog, the method is broken; run that
method-check on the original catalog, train period only).

## EXP-G: ratio-sequence miner (the Fibonacci question, disciplined)

Sequences: mainshocks M >= 5.0 in train (original catalog), aftershocks = events within
100 d & 50 km linked to the largest preceding event. For each sequence with >= 20 aftershocks:
interevent ratios r_i = Δt_{i+1}/Δt_i. MINE (train): histogram of log(r) pooled across
sequences vs null from per-sequence Omori simulations (fit K,c,p per sequence; simulate 500
synthetic sequences each; same n). Any ratio band (width 0.1 in log10) with observed density
exceeding the null 99th percentile is a candidate "favored ratio" — golden ratio φ=1.618
(log10=0.209) evaluated explicitly whether or not it is a candidate. TEST: candidate bands
re-scored on test-period sequences (>= 20 aftershocks), p < 0.05 vs their Omori nulls.
Honest prior: nothing survives; φ shows nothing special. If something survives, it is
reported with its effect size and flagged for replication on an independent region.

## EXP-I: patterns-within-patterns / meta-patterns

(i) WITHIN — sequence-shape recurrence: fit per-sequence Omori (K, c, p) + b-value for all
train sequences (>= 30 aftershocks). Question: do sequences from the same 0.5°-neighborhood
resemble each other more than the global population (ICC / variance decomposition)? TEST:
for each test-period sequence, predict its (p, b) from its neighborhood's train sequences vs
from the global train mean; success = lower squared error for >= 60% of test sequences
(binomial p < 0.05).
(ii) OUTSIDE — meta-clustering of mainshocks: are M >= 5 mainshock times themselves clustered
beyond Poisson (coefficient of variation + Ripley-K in time, train)? TEST: hazard-doubling
rule mined on train (e.g., "P(M>=5 within X d | M>=5 occurred) vs base rate", X chosen on
train) scored on test mainshocks. This is inter-sequence triggering — ETAS at the top of the
hierarchy.
(iii) DRIFT of the pattern-generator: b-value and interevent-CV in 5-yr train windows —
descriptive only unless a monotone trend appears in train, in which case its sign is frozen
and checked on test.

## Scripts / outputs

exp_h_etas.py / exp_f_periodicity.py / exp_g_ratios.py / exp_i_meta.py →
results_exp_{h,f,g,i}.json. Runtime guards: each worker cuts surrogate counts in half
(recording it) if projected runtime exceeds 90 min.

## EXP-J (added 2026-08-09, frozen before any computation): the forward stress ledger

Directive (Jim): model the system "from the other side" - stress/loading/triggers - and find
where the model says events SHOULD fire but do not, then correlate that negative space.

- J1 LEDGER (train period): per 0.2-deg cell in the SoCal box: loading rate
  Mdot_geo = 2*mu*H*A*max_shear_rate (mu = 30 GPa, seismogenic H = 11 km, A = cell area;
  max shear from data/socal_strain_grid.npz, true nanostrain/yr) vs seismic release
  Mdot_seis = sum of Hanks-Kanamori moments (original catalog, train period, M >= 2.5) / 29 yr.
  Coupling chi = Mdot_seis / Mdot_geo. Classes frozen in advance: SILENT-LOADING
  (loading in top quartile, chi < 0.01), COUPLED (chi in [0.01, 1]), OVERSHOOT (chi > 1,
  possible aftershock transients). Requires >= 20 train events per cell for a chi estimate;
  cells below that with top-quartile loading are SILENT-LOADING by definition.
- J2 PREDICTION (single test unblinding): two rival frozen hypotheses about the negative
  space: (a) CATCH-UP - train silent-loading cells have elevated test-period event rates
  relative to their train rates (vs coupled cells, Mann-Whitney on rate ratios);
  (b) PERSISTENCE - cell character (chi quartile) persists train -> test (Spearman of
  log-chi train vs test on cells measurable in both). Both scored; the data decides.
  Stated prior: persistence wins (creep and coupling are material properties), but a
  catch-up signal in any subset would be the interesting hazard-relevant finding.
- J3 CORRELATES (exploratory): silent-loading cells vs covariates - distance to the SAF
  creeping section (approx trace lat 36.0-36.8 along the fault), distance to geothermal
  fields (tsi_map.GEOTHERMAL), median event depth, b-value, and the round-1 tidal modulation
  amplitude (exp_a_train_bins.csv where overlapping). Rank correlations, exploratory-labeled.
- Script exp_j_stress_ledger.py -> results_exp_j.json + maps/exp_j_ledger.png (loading map,
  release map, log-chi map, class map).

## EXP-K (added 2026-08-09, frozen before computation): style-stratified fault-resolved ledger

Upgrades EXP-J's style-blind ledger. Per 0.2-deg cell:
- Full strain-rate tensor (exx, eyy, exy) recomputed from NGL MIDAS via
  strain_comparison.strain_grid_full (nanostrain/yr), aggregated 2x2 to 0.2 deg.
- STYLE class (frozen): styleness s = dilatation / (2*max_shear); TRANSTENSIONAL s > 0.25,
  CONTRACTIONAL s < -0.25, STRIKE-SLIP otherwise.
- FAULT GEOMETRY: CFM5.3 traces (data/xue_lu_zenodo/CFM5.3_traces.lonLat + the fault-geometry
  file); per cell: distance to nearest trace and that trace's local strike; fault-resolved
  loading = shear strain rate resolved onto the local strike (tensor rotation); ON-FAULT if
  distance < 10 km.
- VARIABLE H (frozen coarse model): H = 6 km within 0.35 deg of a tsi_map.GEOTHERMAL field;
  H = 25 km inside the Transverse Ranges box lat [34.0, 34.8], lon [-119.8, -117.3];
  H = 11 km elsewhere. Sensitivity of the silent list to H reported (flat-11 vs variable).
- Outputs: chi and class per cell as in EXP-J but per (style x on/off-fault) stratum;
  Kruskal-Wallis of log-chi across styles; within-stratum J2-persistence (Spearman train/test
  log-chi, >= 20 events both periods); revised silent-loading list EXCLUDING cells within
  25 km of the SAF creeping segment and cells within 0.35 deg of geothermal fields
  ("unexplained silent" = the hazard-candidate list, reported with fault distance/strike).
- Script exp_k_stratified_ledger.py -> results_exp_k.json + maps/exp_k_stratified.png.

## EXP-L (application, not a hypothesis test): live 7-day forecast from frozen ETAS

Uses results_exp_h.json frozen params (mu, K, alpha, c, p, M0=2.5) and train b=1.0654 -
NO refitting (the pre-registered parameters are the point). Fresh catalog: USGS ComCat FDSN
(earthquake.usgs.gov/fdsnws/event/1/), SoCal box, M >= 2.5, 2010-01-01 -> now (paged);
retrieval logged. Forecast: lambda(now) from full history; next-7-day M >= 5 and M >= 4
probabilities via 1,000 forward ETAS simulations (thinning, GR magnitudes b=1.0654, cap
horizon 7 d); also the no-future-triggering analytic lower bound. Output: results_forecast_
<date>.json + console statement with honest framing (probability forecast, conditional
skill, not an event prediction). Script: etas_forecast.py (reusable daily).
