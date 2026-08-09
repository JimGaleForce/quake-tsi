# Hypothesis Ledger — Kepler ⇄ Popper research engine

Working convention (est. 2026-08-09): **Kepler** (explorer persona, .claude/agents/kepler.md)
appends PROPOSED entries; **Popper** (adjudicator persona, .claude/agents/popper.md) moves them
to TESTABLE-NOW / NEEDS-DATA / REFRAMED / REJECTED, and — after the supervisor runs the frozen
tests — to VALIDATED → BASELINE. Kepler then riffs on every new baseline. Popper never blocks
exploration; he gates only what gets claimed.

Current BASELINES (imported from the program record, all pre-registered + out-of-sample):
- B-1: Generic temporal ETAS transfers globally — universal shape (p≈0.94–1.08), locally
  calibrated μ; +0.66..+1.75 bits/event on six never-trained regions (results_exp_m.json).
- B-2: SoCal walk-forward ETAS skill +1.87 bits/event vs test-rate-oracle Poisson; skill RISES
  with magnitude (results_exp_h.json).
- B-3: After any M≥5 (SoCal box), P(M≥5 within 7 d) ≈ 0.60 vs 0.062 base (results_exp_i.json).
- B-4: Stress-ledger negative space is real geology: silent-loading cells recover the SAF
  creeping section and the 1857 Big Bend strand blind (results_exp_j/k.json).
- B-5: The dilatation component of geodetic strain carries ±2× measurement uncertainty across
  velocity solutions; shear is robust (results_strain_comparison.json).

Corpse list (do not re-propose without a genuinely new angle — see EQ18_FULL_NOTES.md §14–17):
static tidal-phase susceptibility maps; tidal-amplitude/feature correlations (n-bias);
fixed periodicities incl. annual (weak-power caveat); Fibonacci/golden interevent ratios;
fault-type parameter pooling; spatial transfer of sequence shapes (p, b); TSI ratio.

---

## PROPOSED (Kepler)

### Round 1 — 2026-08-09 (Kepler). Driving questions: (Q1) when/where the next M>=6.5-7;
### (Q2) do small quakes follow the same rules as big ones, and is the BREAK a predictor;
### (Q3) is there an "internal weather system" that would yield to compute + assimilation;
### (Q4) other negative spaces. All entries PROPOSED — Popper adjudicates.

Assets I am building on (verified on disk this session):
`data/comcat_world/*.csv` = 13 regional ComCat catalogs, M>=4.5, 1995-01-01 -> 2026-08-09,
**71,803 events, 1,589 M>=6.0, 517 M>=6.5, 179 M>=7.0** (cols time/lat/lon/depth/mag/magType/
net/id/...); `data/comcat_socal_m25.csv` (18,382 M>=2.5, 2010->now);
`data/xue_lu_zenodo/SCSN_original_catalog.txt` (633,667 in-box, 1981-2018, Mc<=1.7 all eras);
`data/socal_strain_grid.npz` (66x86, max_shear + dilatation nanostrain/yr);
CFM5.3 traces (557 segments w/ strike/dip/rake); frozen SoCal ETAS
(mu=0.2750/d, K=0.04124, alpha=0.5366, c=0.01426 d, p=1.1183, M0=2.5, b=1.0654);
frozen GLOBAL pool (mu=0.03988, K=0.02114, alpha=0.7296, c=0.01462, p=1.0424, M0=4.5).
Holdout bits/event vs local-oracle Poisson: AK +0.84, MX +0.66, PH +0.78, CAR +1.75,
IRN +0.69, GRC +0.79.

**Standing methodological demand for this whole round** (applies to every entry below unless
stated): *the null model is a simulated ETAS catalog, not a Poisson catalog.* Almost every
"precursor" in the literature dies to this one substitution, because clustering alone
manufactures rising rates, rising variance, rising correlation length, and falling interevent
times before large events. We already own the frozen generators (B-1/B-2) — we can afford the
honest null that most groups cannot. Where I say "ETAS-sim null" I mean: simulate N>=500
synthetic catalogs from the region's frozen ETAS + GR(b), apply the identical statistic,
compare.

---

#### Q1 — WHEN and WHERE the next M>=6.5+

**K-001 — ETAS forecast skill keeps RISING with target magnitude, all the way to M>=6.5+.**
- Claim: scored on M>=6.5 targets only, generic-ETAS-with-local-mu beats a local-oracle Poisson
  by MORE bits/event than it does on M>=4.5 targets — i.e. the largest, most consequential
  events are the *most* relatively forecastable, not the least.
- Inversion: B-2 found SoCal skill rises with magnitude (+1.87 pooled -> +2.54..+2.58 at M>=4).
  Everyone reads that as "aftershocks are easy"; invert it — if skill is monotone in M, then
  the thing we already validated is *most* useful exactly where the stakes are. Nobody has
  extended that curve because per-region large-event counts are small. Pooled across 13 boxes
  we have 517 M>=6.5 and 179 M>=7.0. Power exists; it has just never been aggregated.
- Test: for each of the 13 regions, run the EXP-M walk-forward with GLOBAL-pool params + local
  mu (all frozen, no refit), but compute the log-likelihood contribution restricted to target
  events in magnitude strata [4.5,5), [5,5.5), [5.5,6), [6,6.5), [6.5,7), [7,inf). Baseline per
  stratum = that region's own-period Poisson rate for that stratum (local oracle, harsh).
  Statistic: bits/event per stratum, pooled across regions with per-region weights; slope of
  bits vs magnitude stratum (weighted least squares).
- Null: slope <= 0 (skill flat or decaying with magnitude). Also run the ETAS-sim
  self-consistency check: on synthetic catalogs the slope should reproduce, confirming the
  slope is a property of the model+null pair and quantifying its sampling error.
- Expected effect if real: +1.0 to +2.5 bits/event at M>=6.5, slope +0.2 to +0.5 bits per
  half-magnitude. If instead skill *collapses* above ~M6, that is the single most important
  negative result in this program: it would say large events are drawn from a different
  process and B-1/B-2 do not extrapolate.
- Why this might be dismissed too quickly: "obviously it's just aftershock forecasting, and the
  M6.5 events you're scoring are mostly aftershocks of M7+." Two answers: (a) that is an
  empirical question the stratified table answers directly (report triggered_fraction per
  stratum), and (b) even if true, "a large fraction of M>=6.5 events are conditionally
  forecastable given a preceding large event" is a *deployable* statement, not a deflation.

**K-002 — Establish the WHERE floor before anything clever: global smoothed-seismicity +
local-mu spatio-temporal ETAS, scored on M>=6.5 only.**
- Claim: a Kagan-Jackson-style smoothed-seismicity spatial density, multiplied by the validated
  universal temporal kernel, is a genuine skill floor for global M>=6.5 location, and any
  physics-based spatial hypothesis (K-003 and beyond) must be measured as bits ABOVE it, not
  above uniform.
- Inversion: this program has validated a *temporal* object and keeps reaching for spatial
  physics with no spatial baseline in place. Model the observer's prior first. Also: B-1's
  lesson was "the transferable thing is universal shape + LOCAL mu" — spatially, mu becomes a
  field, and the field IS the where-forecast.
- Test: build mu(x) on a 0.5° global grid from declustered M>=4.5 events 1995-2009 (adaptive
  power-law kernel, bandwidth = distance to k-th nearest event, k in {4,8}); temporal kernel =
  GLOBAL-pool (K, alpha, c, p) with a spatial aftershock kernel (r+d)^-q, d and q fit on the
  13-region train set only. Score 2010->2026 walk-forward on M>=6.5 targets: Poisson
  log-likelihood over space-time-magnitude bins (CSEP-style), reported as bits/event above
  (a) spatially uniform over the union of the 13 boxes and (b) time-independent smoothed
  seismicity alone (this isolates the *temporal* contribution to WHERE).
- Null: uniform-in-footprint; plus a spatially-shuffled-mu control (rotate the mu field by a
  random lag on the grid) to prove the skill is geographic, not just normalization.
- Expected: +2 to +4 bits/event vs uniform, of which only +0.3 to +1.0 comes from the temporal
  layer for M>=6.5 targets. That split is the number Jim actually needs: it says how much of
  "where" is fixed geography and how much is live.
- Why dismissed too quickly: "smoothed seismicity is 1994 technology." Yes — and it still wins
  CSEP tests, which is exactly why it must be the floor. A program that cannot beat it should
  know that in week one, not year two.

**K-003 — Globalize the stress ledger: geodetic moment deficit (GSRM v2.1) predicts where the
next M>=6.5+ goes, over and above smoothed seismicity.**
- Claim: cells whose geodetic loading rate is high but whose catalog moment release is low
  (B-4's SILENT-LOADING class, computed globally) host M>=6.5 events at a rate significantly
  above what smoothed seismicity alone predicts.
- Inversion: B-4 proved the negative space is real geology in SoCal (it rediscovered the SAF
  creeping section and the 1857 Big Bend strand blind). The obvious riff: the ledger is a
  *hazard* statement only if the silent cells eventually fire. SoCal 2010-2018 had no M>=7 to
  test that with. The world does — 517 M>=6.5 in 31 years.
- Data: GSRM v2.1 global strain-rate model (Kreemer et al. 2014, GEM — public grid download,
  0.25° second-invariant + tensor components) — NOT on disk, needs one download.
  Moment release from the 13-box ComCat M>=4.5 on disk, plus ISC-GEM (1904-2018, M>=5.5,
  public) to capture pre-instrumental large-event moment so chi is not biased low the way it
  was in EXP-J. Slab2 (USGS) for seismogenic thickness in subduction zones.
- Test: per 0.5° cell, Mdot_geo = 2*mu*H*A*e_max_shear (H from Slab2/coupling depth, 15 km
  default continental); Mdot_seis from ISC-GEM 1904-1994 (training era, deliberately disjoint
  from the scoring catalog); chi = ratio; freeze the SILENT quartile. Score: Poisson regression
  / likelihood gain for 1995-2026 M>=6.5 counts per cell with predictors [log smoothed
  seismicity] then [+ silent-class indicator]. Statistic: delta log-likelihood of the nested
  models, and the rate ratio silent-vs-coupled.
- Null: label-permutation of the silent class within loading-matched strata (this is essential
  — silent cells are high-loading by construction, and high loading correlates with everything).
- Expected: rate ratio 1.5-3x, delta-LL worth +0.2 to +0.6 bits/event. Honest prior: the
  loading-matched null will eat most of the raw signal; what survives is the interesting part.
- Why dismissed too quickly: "seismic gaps were falsified in the 1990s (Kagan & Jackson)." True
  and I am not proposing the gap hypothesis — I am proposing a *moment-budget* covariate scored
  against a smoothed-seismicity floor with a loading-matched null, which is precisely the test
  the 1990s gap papers lacked. Also, the 1990s test used gap *segments* defined by hand; this
  one is a blind grid.

**K-004 — Depth/thermal state, not tectonic style, is the stratifier that works.**
- Claim: stratifying ETAS parameters by hypocentral depth (or homologous temperature proxy)
  yields the out-of-sample gain that fault-type pooling (EXP-M, FAIL 2/6) did not.
- Inversion: riff directly off a validated *failure*. EXP-M showed alpha/K spread WITHIN type
  exceeds spread BETWEEN types. So the type label is not carving nature at its joints. What
  varies enormously within a "subduction" box? Depth: Chile is 41% >70 km, Indonesia 31%,
  Alaska 10%, Greece 4%. Deep events are famously aftershock-poor (no fluids, no
  rate-and-state at those P-T conditions). "Subduction" is a mixture of two populations and
  the median-pooling in EXP-M blended them.
- Test: on disk, no download. Re-run the EXP-M machinery with each catalog split at depth
  70 km (and 40 km as a secondary cut). Fit shallow-pool and deep-pool (K, alpha, c, p) on the
  7 train regions; score the 6 holdouts with (a) GLOBAL pool [the current champion], (b)
  DEPTH-stratified pool. Frozen success rule, same shape as EXP-M: depth pool beats GLOBAL in
  >=5 of 6 holdouts. Report triggered fraction shallow vs deep per region.
- Null: sign test vs 3/6; plus the diagnostic that matters even on a fail — the shallow/deep
  productivity ratio K_shallow/K_deep with bootstrap CI.
- Expected: K_deep/K_shallow ~ 0.2-0.5, and +0.1 to +0.4 bits/event on deep-heavy holdouts
  (Philippines, Mexico). Small in bits but structurally important: it identifies the real
  latent variable.
- Why dismissed too quickly: "everyone knows deep quakes have fewer aftershocks." Knowing it
  qualitatively and having it as a frozen, out-of-sample, bits-scored modifier in a deployable
  forecaster are different objects. And EXP-M's failure makes this the *cheapest* remaining
  test of whether ANY stratification helps.

---

#### Q2 — Self-similarity: where it holds, where it breaks, is the break a predictor

**K-005 — The M0-invariance audit: fit the same ETAS at ascending magnitude floors and measure
where self-similarity breaks. The measured drift IS a calibration correction for large events.**
- Claim: under exact self-similarity, ETAS (alpha, c, p, branching ratio n) are invariant to the
  magnitude threshold M0 after the standard K rescaling. They will NOT be. The drift will be
  systematic, and correcting for it changes large-event probabilities materially.
- Inversion: instead of asking "are small quakes like big ones", ask "does the model we already
  validated *know it is looking at a truncated sample*?" Self-similarity is not a philosophical
  question here; it is a testable invariance of a fitted parameter vector under a change of the
  observation window in magnitude. That is a renormalization-group question, and it is
  answerable this afternoon with data on disk.
- Test (all on disk): SoCal SCSN original catalog 1981-2018 — fit temporal ETAS at
  M0 = 2.5, 3.0, 3.5, 4.0, 4.5 (train <2010, same MLE harness as EXP-H). Global: pooled per
  region at M0 = 4.5, 5.0, 5.5. Rescale K to a common reference via K' = K*10^(alpha*(M0new -
  M0old)) (the sign the worker corrected in EXP-M). Statistic: regression of each parameter on
  M0; slope with bootstrap CI over 200 catalog bootstraps. Then the *consequential* test:
  predict the M>=4.5 rate in the test window using parameters fit at M0=2.5 vs at M0=4.0 —
  bits/event difference.
- Null: parameter slopes = 0. Critically, also fit the same ladder on ETAS-SIMULATED catalogs
  (which are self-similar by construction) to calibrate how much apparent drift comes from
  estimator bias and finite samples rather than from the Earth. That control is what makes this
  entry worth running; without it the result is uninterpretable.
- Expected: alpha rises with M0 (0.54 at 2.5 -> 0.7-0.9 at 4.5, consistent with the SoCal-vs-
  global gap we already see: 0.537 vs 0.730), n falls with M0. If the simulated control
  reproduces the alpha drift, the "global alpha is bigger" story is an artifact of magnitude
  floor, not tectonics — which would retroactively reinterpret EXP-M's parameter table.
- Why dismissed too quickly: "known estimator bias, boring." That is exactly the point: if it
  is bias, we must know its size before we quote transfer parameters across catalogs with
  different M0 (which this program already does); if it exceeds the simulated bias, we have
  measured a genuine break in crustal self-similarity. Either outcome is load-bearing.

**K-006 — Productivity saturation: alpha breaks (bends down) above some magnitude.**
- Claim: aftershock productivity does not follow a single 10^(alpha*M) law across the whole
  range; there is a bend, and the bend magnitude is a property of the region's seismogenic
  thickness / fault system size.
- Inversion: alpha is the exponent that makes big events matter in ETAS. If it saturates,
  every ETAS forecast systematically OVER-predicts aftershock hazard after the biggest events
  — which are exactly the moments the forecast is used. Model the residual of the productivity
  law rather than the law.
- Test (on disk): for all 13 regions pooled, for every mainshock M>=5.0, count aftershocks
  M>=4.5 within 100 d and a magnitude-scaled radius (2 rupture lengths, Wells-Coppersmith).
  Regress log10(N_aftershocks) on M in 0.5-wide bins; fit (a) single line, (b) bilinear with a
  free breakpoint. Statistic: BIC difference and the fitted breakpoint with profile CI.
  Stratify by depth (ties to K-004) and by seismogenic thickness proxy.
- Null: bilinear model not preferred (delta-BIC < 6); breakpoint CI covering the catalog's Mmax
  (i.e. the "break" is just the edge of the data). Bootstrap over mainshocks.
- Expected: a bend near M 7.0-7.5 in shallow continental settings, none (or a much higher one)
  in subduction. Effect on forecasts: 20-50% over-prediction of 7-day M>=6 aftershock
  probability after an M8 if the bend is real and ignored.
- Why dismissed too quickly: "finite catalog, magnitude saturation of the *magnitude scale*,
  not the physics." Legitimate — so pre-register the magType audit (mww vs mb vs ms in the
  ComCat rows; use GCMT Mw where available) as part of the test, and run the identical analysis
  on GCMT-only magnitudes as a robustness arm. If the bend survives a homogeneous Mw catalog,
  the scale-saturation objection is dead.

**K-007 — The break as the predictor: loss of Gutenberg-Richter (a bend/taper appearing in the
local magnitude distribution) precedes large events.**
- Claim: in a trailing window before a large event, the local magnitude distribution deviates
  from a pure exponential in a *specific* direction — not merely a lower b, but a change in
  shape (curvature / preference for a tapered-GR over GR by likelihood ratio).
- Inversion: b-value precursors are a swamp because b is one number and it drifts with
  completeness. Invert to the *shape* of the distribution: completeness errors bend the low-M
  end; physical criticality changes the high-M end. Test the high-M curvature only, above a
  conservative Mc, so the observer artifact and the signal live in different parts of the
  curve. That separation is what makes this different from the 100 dead b-value papers.
- Test: SoCal SCSN 1981-2018 (M>=3.0, cleanly above Mc=1.7 in all eras) and the 13 world boxes
  (M>=5.0). For each M>=6.0 target event, take the trailing 2-year, 100 km sample; fit GR and
  tapered-GR (corner magnitude free); statistic = likelihood-ratio in favor of the taper, and
  the fitted corner magnitude. Compare the distribution of that statistic in pre-target windows
  vs matched control windows (same cell, same season, no target within +/-1 yr) — matched on
  sample size, which is the killer confound.
- Null: ETAS-sim null (frozen params + pure GR b). Under the null there is NO taper ever, so
  any systematic pre-event taper signal is a genuine departure from the generative model.
  Report ROC AUC for "target within next 1 yr".
- Expected: AUC 0.55-0.65 if real. Small — but a 0.6 AUC covariate multiplied into a validated
  ETAS is worth real bits, and the direction of the taper (corner moving DOWN before a large
  event = the system suppressing large events, vs UP = enabling them) is a physics discriminator
  either way.
- Why dismissed too quickly: "b-value precursors have failed a hundred times." They have — as
  a *mean* statistic with a Poisson null and no sample-size matching. All three of those
  failure modes are fixed here, and the statistic is different (shape, not level).

**K-008 — Does a sequence know how big it will get? Early-sequence Mmax predictability vs the
GR draw.**
- Claim: given the first 24 h of a new cluster, the eventual maximum magnitude of that cluster
  is predictable beyond what ETAS+GR implies (which is: not at all, except through rate).
- Inversion: ETAS's most radical assumption is that magnitudes are i.i.d. GR draws independent
  of everything. That is the assumption that makes prediction impossible. Attack it directly
  and quantitatively, at the one place where it is cheapest to test: cluster onset.
- Test (on disk): SoCal SCSN M>=2.5, 1981-2018, and the 13 world boxes M>=4.5. Define clusters
  by nearest-neighbor distance (Zaliapin-Ben-Zion, the parameters this program already uses:
  b=1, df=1.6, log10 eta0=-5). For each cluster, features from its first 24 h ONLY: event count,
  rate slope, spatial extent, depth spread, local b of those events, distance to nearest CFM
  trace / plate boundary, and the mu of the host cell. Label: eventual cluster Mmax. Model:
  gradient-boosted quantile regression, trained on clusters starting <2010, scored on >=2010.
  Statistic: CRPS improvement over the ETAS+GR baseline predictive distribution for Mmax.
- Null: the ETAS+GR baseline itself; and a permutation of labels within
  (first-24h-count) strata — because count alone predicts Mmax trivially through GR sampling
  and that must not be counted as skill.
- Expected: honest prior is a small win, 3-8% CRPS, mostly from the spatial features (i.e.
  "where" leaking in, not "how it started"). A large win would be the biggest result in the
  program. A clean zero is also valuable: it hardens the i.i.d.-magnitude assumption with a
  number, which nobody has done on this data.
- Why dismissed too quickly: "the deterministic-nucleation debate is settled, magnitudes are
  random." It is not settled, it is *stalemated* on waveform data; this is a catalog-level test
  with a proper null and it costs one afternoon.

---

#### Q3 — Is there an internal weather system? What is the state vector? Cheapest toehold?

**K-009 — THE DIAGNOSTIC: are ETAS residuals white, or is there weather? Spatio-temporal
autocorrelation of the residual field decides whether assimilation is worth any compute.**
- Claim: the residual field of a validated ETAS forecast (observed minus expected counts, per
  cell per time step) is NOT white noise — it has spatial correlation length >0 and temporal
  correlation time >0 — and those two numbers are the entire business case for a data
  assimilation program.
- Inversion: the weather analogy fails or succeeds on one property: does the system have a
  slowly-varying hidden state that observations can constrain? In DA language, the residual IS
  the innovation. If innovations are white, there is no state to estimate and no amount of
  compute helps — ETAS is already the filter. If innovations are red and spatially coherent,
  there is a latent field (fluids? transient creep? slow slip? stress shadows?) and it is
  estimable. This is the cheapest, most decisive, most under-asked question in the whole program
  and it needs zero new data.
- Test (100% on disk): SoCal, 0.2° x 7-day cells, 2010-2018, expected counts from the frozen
  EXP-H ETAS (walk-forward, no refit) rescaled to spatial cells by a smoothed-seismicity kernel;
  residual r = (obs - exp)/sqrt(exp+1) (Anscombe/Pearson). Compute: (i) temporal ACF of r per
  cell, pooled, lags 1-52 weeks; (ii) spatial variogram / Moran's I of r at each time step,
  pooled; (iii) the leading EOF of the residual field and its temporal power spectrum (is there
  a low-frequency mode?). Repeat identically on the 13 world boxes at 0.5° x 30-day.
- Null: residuals computed the same way on ETAS-SIMULATED catalogs from the same frozen
  parameters. This gives the exact distribution of ACF/Moran's I under "the model is true and
  there is no weather", including all the correlation that ETAS itself manufactures.
- Expected if real: lag-1 weekly ACF excess of 0.05-0.20 over the ETAS-sim null; Moran's I
  excess at 20-60 km; a leading EOF explaining 10-25% of residual variance with a red spectrum.
  A correlation TIME of months and a correlation LENGTH of tens of km would be the go-signal
  for a real assimilation effort — and would also tell us the required observing density.
- Why dismissed too quickly: "residual diagnostics are just goodness-of-fit; we know ETAS is
  incomplete." Right — but "incomplete" is not actionable and "the missing structure has a
  correlation time of 4 months and a correlation length of 35 km" is. This entry converts a
  vague dissatisfaction into a specification for an instrument. Run it first.

**K-010 — Cheapest assimilation toehold: make mu a latent state, not a constant. Two tiers,
EWMA then particle filter.**
- Claim: replacing the frozen background rate mu with a *filtered, time-varying* mu_t —
  estimated online from the event stream itself — beats frozen-mu ETAS out of sample, and the
  gain grows with the residual correlation time measured in K-009.
- Inversion: B-1's headline was "universal shape + LOCAL mu". The unexamined word is "local" —
  local in space. Make it local in TIME too. mu is the only place in ETAS where an external,
  slowly-varying driver can enter, so mu_t IS the minimal state vector. This is
  data assimilation with a one-dimensional state, which is the honest place to start.
- Test (on disk, no downloads): SoCal M>=2.5 2010-2026 (comcat_socal_m25.csv, 18,382 events)
  and the 6 EXP-M holdouts. Tier 1 (an hour of work): mu_t = EWMA of the ETAS-declustered
  background rate with timescale tau; tau chosen ONLY on the train period / train regions,
  frozen, then walk-forward scored. Tier 2: a bootstrap particle filter with state
  log mu_t following an Ornstein-Uhlenbeck process (theta, sigma fit on train), observation
  model = the ETAS conditional intensity; 5,000 particles; walk-forward one-step-ahead
  predictive log-likelihood.
- Null: frozen-mu ETAS (the current baseline B-1/B-2) — this is a strictly harder null than
  Poisson and the right one. Plus a "cheating control": a filter given a randomly time-shifted
  event stream for its update step should score at or below frozen-mu.
- Expected: +0.05 to +0.3 bits/event. Small, but it is the first *dynamical state* in the whole
  program, and Tier 2's estimated OU timescale is directly comparable to K-009's measured
  residual correlation time — two independent routes to the same number is a strong internal
  consistency check.
- Why dismissed too quickly: "time-varying background is just refitting / overfitting." The
  frozen-tau + walk-forward + time-shifted-stream control design makes overfitting detectable.
  And if Tier 1 wins even slightly, the app gets a real upgrade for near-zero cost.

**K-011 — Measure the predictability horizon: how much could compute EVER buy?**
- Claim: the seismicity system has a measurable finite predictability horizon (an error-doubling
  time), and it is short enough that we can state, quantitatively, the ceiling that perfect data
  assimilation would hit.
- Inversion: everyone argues about whether earthquakes are predictable. Nobody in this program
  has measured the *Lyapunov-analogue*. In weather, the 2-week horizon is a measured property,
  not an opinion, and knowing it is what made NWP a rational investment. Do the seismic version.
- Test (on disk): ensemble twin experiment. From a given state (full catalog history up to t0),
  launch 1,000 ETAS simulations forward; measure the growth of ensemble spread in a coarse
  observable (7-day M>=4 count in the box; and moment release) as a function of lead time. Then
  the *identical-twin* version: take one ensemble member as "truth", perturb the initial history
  slightly (drop/add the smallest 1% of events, jitter magnitudes by their reported magError),
  and measure divergence. Statistic: lead time at which RMS spread reaches 1/sqrt(2) of the
  climatological spread ("saturation time"), as a function of the coarse-graining scale
  (space: 0.2°/1°/box; time: 1 d/7 d/30 d).
- Null: not a hypothesis test — a measurement. The reference is the climatological (Poisson)
  spread, and the reported product is a saturation-time surface over coarse-graining scales.
- Expected: saturation in days for fine scales, weeks-to-months for coarse. That surface IS the
  answer to Jim's Q3: it will show (I predict) that fine-grained event prediction saturates
  almost immediately while coarse aggregates retain skill for months — which is precisely the
  weather situation (no raindrop forecast, good rainfall forecast) and tells us what product to
  build. Pairs with K-027.
- Why dismissed too quickly: "you're measuring your model, not the Earth." Correct, and stated
  up front: this is the predictability horizon *of the validated generative model*, which is a
  lower bound on what a better model could do only if the better model has the same stochastic
  core. Even so, an ETAS-internal horizon of 3 days at 0.2° is decision-relevant: it says stop
  buying compute for that product.

**K-012 — Specify the state vector, then buy the cheapest missing component: transient strain
from NGL daily GNSS.**
- Claim: the SoCal state vector v0 = [mu_t field, b_t field, chi (loading/release ledger),
  transient strain rate field, tremor rate, pore-pressure proxy] is 4/6 constructible today; the
  highest-value missing component is *transient* (not secular) strain, and it is a one-command
  download away.
- Inversion: model the driver, not the events (this program's own founding move) — but make the
  driver TIME-DEPENDENT. B-5 told us dilatation is measurement-noise-limited and shear is robust;
  data/ngl/midas.IGS14.txt gives us only secular velocities. The negative space in our own data
  inventory is the *time series behind the velocity*. NGL publishes daily .tenv3 position series
  for every station; ~1,000+ in the SoCal box; the secular fit we use has already thrown away
  exactly the signal a weather model would want.
- Data (downloadable, concrete): geodesy.unr.edu/gps_timeseries/tenv3/IGS14/<STA>.tenv3 for
  stations in the SoCal box (station list at geodesy.unr.edu/NGLStationPages/llh.out). Also
  free: GCMT (globalcmt.org, 1976-now moment tensors, global — the FM layer we lack outside
  SoCal); ISC-GEM; PNSN/Cascadia tremor catalogs (the strong-tremor analogue for a positive
  control); GRACE/GLDAS hydrologic loading grids for the pore-pressure/annual-loading proxy.
- Test: detect transients with a common-mode-filtered, trajectory-model residual (secular +
  annual + semiannual + coseismic offsets removed); form a weekly transient shear-strain-rate
  field on the 0.2° grid; then the only question that matters — is that field a covariate of
  the K-009 residual field? Statistic: cross-correlation of transient strain-rate anomaly with
  ETAS residual, at lags -12..+12 weeks, per cell and pooled; and the bits/event gain from
  adding it as a multiplicative covariate on mu_t in the K-010 filter.
- Null: circular time-shift of the transient field against the residual field (2-2000 day
  shifts, the same null machinery this program already validated in round 1), plus a
  station-label-permuted spatial control.
- Expected: honestly, near-null in SoCal — the transients are small and GNSS noise at weekly
  scale is ~1-2 mm. The value is that it gives a MEASURED noise floor: "transient strain must
  exceed X nanostrain/yr to be detectable as an ETAS covariate at current station density",
  which is a specification for what a denser network would need to deliver. That is a real
  deliverable even from a null.
- Why dismissed too quickly: "GNSS is too noisy at these timescales, everyone knows." Then the
  test costs a weekend and produces the quantitative version of "everyone knows" — and if
  Cascadia tremor (positive control, where transients are huge) DOES show the coupling, we have
  proven the method works and only SoCal's amplitude is lacking. That asymmetry is why the
  positive control is part of the design.

---

#### Q4 — Other negative spaces (inversions beyond the stress ledger)

**K-013 — The dynamic ledger: aftershock DEFICIT. Cells that should be aftershocking and are
not.**
- Claim: after a mainshock, the sub-regions of the aftershock zone that produce *fewer*
  aftershocks than ETAS predicts are the sub-regions that host the sequence's largest late
  events (or the next mainshock) — silence inside an active zone is a stress-shadow/locked
  signature, not an absence of information.
- Inversion: B-4 inverted the *secular* ledger (geodetic loading vs long-term release). This
  inverts the *transient* ledger: ETAS supplies the expectation, so the negative space becomes
  computable at 7-day resolution instead of 29-year resolution. Same move, four orders of
  magnitude faster.
- Test (on disk): SoCal SCSN 1981-2018 + comcat_socal_m25 2010-2026. For each M>=5.5 mainshock,
  build the 0.1° expected-aftershock map from the frozen ETAS + spatial kernel over days 1-30;
  compute per-cell deficit d = (exp - obs)/sqrt(exp). Score: does d over days 1-30 predict the
  location of M>=4 events in days 31-365, above a baseline that uses the *observed* days 1-30
  activity? Statistic: AUC and Poisson-regression delta-LL. Also do the world version on the
  13 boxes with M>=6.5 mainshocks and M>=5 late targets.
- Null: ETAS-sim — simulate sequences from the fitted model, compute the identical deficit map,
  and measure how often the "deficit predicts late events" statistic arises from pure sampling
  noise (it will arise sometimes, because low-count cells are noisy — hence the sqrt(exp)
  normalization and the ETAS-sim calibration).
- Expected: AUC 0.55-0.62. Small, deployable, and it is a *live* product: after any mainshock
  the map can be issued within days.
- Why dismissed too quickly: "that is just the rupture area vs the unruptured asperity, known
  since the 1980s." Then it should be easy to demonstrate quantitatively — and it never has
  been, blind, out of sample, with an ETAS-sim null, in a form you can put in an app.

**K-014 — Time-reversal asymmetry: run ETAS backwards. The forward-minus-backward skill gap is
the causal information content of foreshocks.**
- Claim: a time-reversed catalog scored with a time-reversed ETAS achieves LESS skill than the
  forward version, and the gap quantifies how much of clustering is genuine causal triggering
  versus symmetric co-clustering — with the surprising possibility that for the largest targets
  the gap is small, meaning foreshock information ~ aftershock information.
- Inversion: literally invert the arrow of time. Aftershock decay is the most robust fact in
  seismology and it is manifestly time-asymmetric; foreshocks are its faint mirror. Nobody
  routinely measures the asymmetry as a single number, but that number bounds how much
  precursory information is even present in the catalog.
- Test (on disk, cheap): reverse the time axis of each of the 13 world catalogs and SoCal; fit
  ETAS on the reversed train period; walk-forward score on the reversed test period; compare
  bits/event forward vs backward, per magnitude stratum (ties to K-001). Also report the
  reversed-fit parameters (a "backward Omori p" and "backward alpha").
- Null: for a purely symmetric clustering process the two scores are equal; the sampling
  distribution of the difference comes from ETAS-sim catalogs (which are strictly causal, so
  they set the *maximum* expected asymmetry).
- Expected: forward beats backward by 0.3-0.8 bits/event at M>=4.5. If the gap SHRINKS with
  target magnitude — i.e. backward skill at M>=6.5 approaches forward skill — that says the
  catalog contains as much information before a large event as after it, which is the strongest
  possible catalog-only argument that foreshock-based prediction is not hopeless. That single
  plot would be worth the round.
- Why dismissed too quickly: "time-reversal is a cute trick with no physical meaning." Its
  meaning is precise: it is a model-free upper bound on precursory information in the event
  stream, obtained without positing any precursor mechanism. Cheap, clean, and nobody has the
  plot.

**K-015 — Invert the observer: completeness Mc(x,t) is a signal, and post-mainshock detection
loss is a near-source rate gauge.**
- Claim: the transient rise in Mc immediately after a mainshock (events lost in the coda) is
  itself proportional to the true near-source rate, and using it recovers information that
  every catalog-based method currently throws away — improving early aftershock forecasts in
  the first hours, which is when forecasts are most needed and most wrong.
- Inversion: model the observer as part of the system (persona rule #1). The missing events are
  not noise; the *pattern of missingness* is a measurement of the thing that caused it.
- Test (on disk): SoCal SCSN 1981-2018. For each M>=5 mainshock, estimate Mc(t) in log-spaced
  windows over days 0-30 (b-value-stable / maximum-curvature and the Ogata-Katsura joint
  estimator). Statistic 1: is Mc(t) - Mc_background well fit by g*log10(rate) with a stable g?
  Statistic 2 (the payoff): forecast the day-1-to-7 M>=4 count using (a) observed above-Mc
  counts only, versus (b) observed counts corrected for the estimated detection function.
  Bits/event, walk-forward, train/test split at 2010.
- Null: mainshock-shuffled control (apply one event's Mc(t) curve to another's rate history);
  and the ETAS-sim null with a simulated detection function of known parameters, to verify the
  estimator recovers g.
- Expected: g ~ 0.5-1.0 (i.e. Mc rises ~0.5-1 unit per decade of rate); corrected forecasts
  worth +0.1 to +0.4 bits/event in the first 48 h. That is the window where operational
  aftershock forecasts currently fail hardest.
- Why dismissed too quickly: "short-term aftershock incompleteness is well documented (Omi,
  Hainzl, Page)." Documented as a *nuisance to correct*. I am proposing it as a *measurement
  channel* — the detection loss estimates the rate you could not observe — and pairing it with
  a bits-scored forecast test rather than a parameter-recovery exercise.

**K-016 — Global teleconnection: do regional M>=6.5 rates share a common mode?**
- Claim: monthly M>=6.5 counts across the 13 regions are positively correlated beyond
  independent inhomogeneous Poisson processes — a global mode, however weak.
- Inversion: the weather analogy taken literally. Weather has ENSO; if the crust has any
  planetary-scale coupling (viscoelastic, rotational, hydrologic, or via great-earthquake
  static/dynamic triggering), it shows as a shared low-frequency mode in regional rates.
- Test (on disk): monthly M>=6.5 counts per region 1995-2026 (517 events; sparse but the pooled
  covariance is estimable). Statistic: leading eigenvalue of the 13x13 count-correlation matrix
  after removing each region's own ETAS-expected rate (i.e. correlate the RESIDUALS from
  K-009). Also a direct excess-rate test: global M>=6.5 rate in the 7 days after any M>=8.0,
  outside the source region, vs the ETAS expectation.
- Null: ETAS-sim per region, independent across regions, 2,000 realizations -> distribution of
  the leading eigenvalue. This null is essential; shared *trends* in catalog completeness would
  otherwise fake a mode.
- Expected: honest prior is NULL for the low-frequency mode; the M>=8 dynamic-triggering arm has
  a modest literature prior and might show a 1.2-2x excess in the first 24 h. Low prior overall
  — but the cost is one script and the residual matrix already exists from K-009.
- Why dismissed too quickly: it deserves substantial skepticism; I flag it myself as the
  lowest-prior entry in this round. Its value is that it is nearly free once K-009 is run, and
  a clean null on a global mode is a useful boundary for the whole "internal weather" framing.

**K-017 — NEW ANGLE ON A CORPSE: periodicity is dead in calendar time — test it in NATURAL
time (event-count time).**
- Claim: recurrence that is invisible in calendar time (EXP-F: no confirmed periodicity, comb
  power limited) becomes visible when the clock is the cumulative count of smaller events rather
  than seconds — i.e. large events are quasi-periodic in "natural time" with CV<1, even though
  they are Poisson-like in calendar time.
- What is genuinely new vs the corpse: EXP-F tested fixed *calendar* periods with a Rayleigh
  statistic and smoothed-rate surrogates. This changes the independent variable, not the
  statistic. Under self-similarity + a loading-driven system, the natural clock is accumulated
  deformation, whose best catalog proxy is cumulative small-event count (or cumulative Benioff
  strain). Rate fluctuations that destroyed the calendar-time test are exactly what natural time
  divides out. This is a different hypothesis, not a rerun.
- Test (on disk): per region, define the count-clock N(t) = number of M>=4.5 events since
  1995 (and a second clock = cumulative sqrt(moment)). Map each M>=6.5 event onto its clock
  value. Statistic: coefficient of variation of interevent *count* gaps between successive
  M>=6.5, and the Brownian-passage-time / lognormal-vs-exponential likelihood ratio for the gap
  distribution, per region and pooled. Same for SoCal M>=5 on an M>=2.5 clock (much better
  statistics: 130 mainshocks).
- Null: ETAS-sim catalogs (which are self-similar and have NO characteristic recurrence in
  either clock) — measure the CV distribution there. Under ETAS, CV in count-time should be >=1
  (clustering); CV<1 pooled would be a genuine quasi-periodicity that ETAS cannot produce.
- Expected: pooled CV in calendar time ~1.5-2.5 (clustered); in count time I predict 0.9-1.3.
  A pooled CV significantly below the ETAS-sim null would be a major result — it would mean the
  system has a memory variable (a loading clock) that ETAS lacks, and it would immediately
  supply a hazard term.
- Why dismissed too quickly: "periodicity already failed here" — see the paragraph above; the
  clock is different and the corpse's failure mode (rate fluctuations swamping the phase) is
  precisely what this construction removes. Second objection: "natural time / nowcasting has a
  reputation problem." Fair; the answer is the ETAS-sim null, which is the discipline that
  branch of the literature usually skips.

---

### Emergence lens (added mid-round, per Jim's directive)

Framing: treat the crust as a many-body system at or near a critical point; earthquakes are
avalanches, i.e. *samples* of a collective state. The move I want this program to make is the
thermodynamic one — stop trying to predict molecules, find the temperature. Every entry below
names a candidate ORDER PARAMETER: a system-level scalar or field with its own dynamics, whose
forecastability does not require any individual event to be forecastable. The standing
methodological demand still applies and matters even more here: **ETAS-simulated nulls**. ETAS
is itself a critical-ish branching process and manufactures most of the "criticality
precursors" people report. Anything that survives an ETAS-sim null is a genuine departure from
the validated generator — that is the bar.

**K-018 — ORDER PARAMETER #1: the branching ratio n(t) as the crust's distance to criticality
("seismic temperature"). It rises before large events.**
- Claim: a rolling estimate of the branching ratio n(t) (mean offspring per event) is a
  trackable system-level variable, it fluctuates on timescales of months-to-years, and its
  excursions toward n=1 (the critical point) precede large events beyond what ETAS's own
  clustering explains.
- Emergence logic: n is *the* control parameter of a branching process — subcritical n<1 means
  avalanches die, n=1 means scale-free avalanches of unbounded size. SOC says the crust
  self-organizes to n≈1. If n wanders, the crust is not permanently critical but breathes
  around the critical point, and that breathing is exactly a forecastable system-level dynamic
  even though each avalanche is random. Note EXP-H already measured a nominal n=1.16 in SoCal
  (supercritical!) with a 7-day n_eff=0.60 — that discrepancy is itself a clue that n is
  window-dependent and worth tracking rather than fitting once.
- Test (on disk): SoCal M>=2.5 1981-2026 and the 13 world boxes M>=4.5. Estimate n(t) two ways
  in rolling 1-yr windows stepped 1 month: (a) refit K only (alpha, c, p frozen) by MLE;
  (b) the model-light estimator n = (triggered events)/(all events) from stochastic declustering
  (probabilistic parentage under the frozen model). Statistic: mean n in the 1 yr before each
  M>=6.5 (world) / M>=5.5 (SoCal) versus matched control windows; and ROC AUC for "large event
  within next 6 months" using n(t) alone.
- Null: ETAS-sim with CONSTANT true n. This is the whole test — in simulation n(t) will still
  appear to rise before large events (because large events follow busy periods), and the
  question is strictly whether the observed rise exceeds that.
- Expected: AUC 0.55-0.65 if real; observed pre-event delta-n of 0.05-0.15 above the sim null.
  Even a null is valuable: it would establish that n is constant, i.e. the crust sits at a FIXED
  distance from criticality, which is itself a strong physical statement and simplifies every
  forecast.
- Why dismissed too quickly: "n is just K in disguise and it's unstable to estimate." Both true
  — which is why estimator (b) exists and why the null is a simulation with the same estimator
  applied. The instability is shared by signal and null, so it cancels.

**K-019 — Critical slowing down: rising lag-1 autocorrelation AND rising variance of the
coarse-grained seismicity field before large events (the ecology/climate early-warning battery,
imported wholesale but with the right null).**
- Claim: the standard critical-slowing-down early-warning signals (lag-1 AR coefficient and
  variance of a rolling aggregate, plus the "skewness flip") rise in the years before M>=7
  events, in excess of the ETAS-sim null.
- Emergence logic: near a bifurcation/critical transition, the recovery rate from perturbations
  goes to zero — the system becomes sluggish and more variable. This is a *system-level*
  signature that has been validated in lakes, ice cores, epilepsy, and financial crashes. It has
  been tried in seismology sporadically and always with weak nulls.
- Test (on disk): per region, form a monthly aggregate series x(t) = coarse observable (choose
  three: count of M>=4.5; log cumulative moment; median interevent time). Detrend with a
  Gaussian kernel (bandwidth frozen at 5 yr). Compute rolling (3-yr) AR1 and variance. Statistic:
  Kendall tau of AR1 and of variance over the 3 yr preceding each M>=7.0 target; pooled across
  the 179 M>=7 events on disk; and a matched-control distribution from windows with no target.
- Null: ETAS-sim, identical pipeline, identical target definition. Report the observed pooled
  Kendall tau against the simulated distribution. Second control: shuffle which region's target
  times are used on which region's series.
- Expected: if real, pooled tau +0.1 to +0.25 with sim-null 95th percentile around +0.05-0.10.
  A clean null here would be a genuinely important negative: it would say the crust does not
  approach its transitions the way lakes and ice sheets do, i.e. earthquakes are not classical
  bifurcations.
- Why dismissed too quickly: "EWS have failed in seismology." They have failed with Poisson
  nulls, single regions, and post-hoc window choices. 179 M>=7 targets, frozen bandwidths, an
  ETAS-sim null, and three pre-specified observables is a materially stronger design than the
  attempts that failed. And the failure mode is informative rather than merely disappointing.

**K-020 — ORDER PARAMETER #2: correlation length xi(t) of the seismicity field grows before
large events (accelerating moment release, rebuilt properly).**
- Claim: the spatial correlation length of small-event locations, estimated from the two-point
  correlation function in rolling windows, increases before large events — the physical
  correlation length growing toward the rupture dimension.
- Emergence logic: the defining signature of an approaching critical point is a diverging
  correlation length. If a large rupture is a correlated failure of many patches, the patches
  must become correlated first. This is the mechanism AMR was groping at; AMR measured the wrong
  thing (cumulative Benioff strain in a hand-drawn circle) and died in prospective CSEP tests.
- What is new vs the AMR corpse: (i) measure xi directly from the pair-correlation function
  rather than inferring it from a cumulative-energy curve fit; (ii) no hand-chosen radius — xi is
  a fitted length scale, and the region is a fixed grid cell decided in advance; (iii)
  ETAS-sim null (AMR's fatal omission — Omori clustering alone grows apparent correlation
  length); (iv) prospective walk-forward scoring.
- Test (on disk): SoCal SCSN M>=2.5 (excellent density) and Japan/Indonesia/Chile M>=4.5. In
  rolling 1-yr windows on a fixed 1° grid, compute the pair-correlation function C(r) for events
  in the window, fit C(r) ~ r^-d2 * exp(-r/xi), record xi. Statistic: trend in xi over the 3 yr
  before M>=6.5 targets vs matched controls; AUC.
- Null: ETAS-sim with a fixed spatial kernel (in which xi is a constant by construction).
- Expected: pre-event xi growth of 20-60% above the sim null if real; AUC 0.55-0.65.
- Why dismissed too quickly: "AMR is a corpse." It is — and I am saying so explicitly and
  changing the measured quantity, the region definition, and the null. If Popper judges this too
  close to the corpse to be worth the compute, the fallback is to run it only as a diagnostic
  attached to K-022 (which needs the same machinery).

**K-021 — ORDER PARAMETER #3: it is not b, it is the spatial HETEROGENEITY of b.**
- Claim: the variance of b across cells within a region (not its mean) is the informative order
  parameter — heterogeneity collapses (b becomes spatially uniform) as a region approaches a
  large event, because stress becomes uniformly high.
- Emergence logic: b maps to a stress "temperature" (high stress -> low b). A system far from
  failure is a patchwork of stress states; a system approaching a system-spanning failure has
  homogenized. So the intensive variable's *variance* carries the phase information, exactly as
  susceptibility (a variance) rather than magnetization diverges at a critical point in a
  paramagnet. This is a real theoretical reason to look at the second moment and it is why the
  hundred dead b-value-mean papers may have been looking at the wrong moment of the distribution.
- Test (on disk): SoCal SCSN M>=2.5 (b per 0.2° cell with n>=100, 3-yr rolling windows) and
  Japan/Indonesia M>=4.5 at 0.5°. Statistic: Var_cells(b) and its trend in the 3 yr before
  M>=6.0 (SoCal M>=5.5) targets, vs matched controls. Crucially, correct for sampling variance
  of b_hat (which is b^2/n) — subtract the expected sampling contribution so we measure TRUE
  heterogeneity, not count fluctuations. That correction is the technical crux.
- Null: ETAS-sim with a spatially UNIFORM true b (so all measured heterogeneity is sampling);
  and a second sim with a fixed, frozen heterogeneous b-field (so heterogeneity exists but does
  not evolve). Two nulls distinguish "there is heterogeneity" from "the heterogeneity changes".
- Expected: if real, a 15-40% decline in corrected Var(b) in the 1-2 yr before large events.
- Why dismissed too quickly: "b-value work is discredited." The dismissal is aimed at the mean;
  this is the variance, with the sampling term explicitly removed and two simulated nulls. It is
  a different measurement and it has an actual statistical-physics rationale.

**K-022 — ORDER PARAMETER #4 (my favorite): PERCOLATION. The size of the largest connected
cluster of recently-active cells bounds how big the next rupture can be.**
- Claim: define cells "active" if they hosted an event above threshold in a trailing window;
  the largest connected component S_max of the active set is a percolation order parameter, and
  (a) it undergoes sharp transitions rather than drifting smoothly, and (b) its linear extent
  L(S_max) bounds/forecasts the magnitude of the next large event in the region — giving a
  HOW-BIG forecast, not just a rate.
- Emergence logic: a large rupture is a connected path of failure across many patches. In
  percolation, the system-spanning cluster appears abruptly at a critical density. If the crust
  approaches large events by growing a connected stressed backbone, then the *connectivity* of
  activity — not its amount — is the variable that controls maximum event size. This directly
  addresses the thing ETAS cannot do: ETAS draws magnitudes i.i.d. from GR and therefore has no
  concept of "how big can this system currently go".
- Test (on disk): SoCal SCSN M>=2.5 on a 0.05° grid (dense enough for real connectivity),
  trailing windows T in {30, 90, 365} d frozen in advance; also the world boxes M>=4.5 at 0.25°.
  Build the active-cell graph (4-connectivity, and a fault-aware variant using CFM5.3 traces so
  connectivity follows structures rather than Euclidean neighbors). Track S_max(t) and the
  cluster-size distribution exponent. Statistic 1 (order-parameter behaviour): is the
  distribution of S_max bimodal / does it show a jump, versus the smooth distribution from
  ETAS-sim? Statistic 2 (the payoff): quantile regression of next-large-event magnitude on
  log L(S_max), out of sample; CRPS gain vs the GR baseline for Mmax.
- Null: ETAS-sim (which produces spatial clustering and therefore connectivity, but with no
  connectivity-magnitude coupling); plus a spatial-randomization control that preserves the
  number of active cells but destroys their arrangement — this is the essential one, because it
  isolates *geometry* from *amount*.
- Expected: if real, 5-15% CRPS gain on Mmax and, more importantly, a usable statement of the
  form "the currently connected backbone cannot support a rupture larger than M x.x". Even a
  weak version of that is a product nobody ships.
- Why dismissed too quickly: "activity connectivity is just a proxy for the aftershock zone of
  the last event." That is precisely what the spatial-randomization null and the ETAS-sim null
  are for; and the fault-aware variant makes the geometric claim sharper than the Euclidean one.
  Of everything in this round, this is the entry with the largest gap between "obvious in
  hindsight" and "actually tested".

**K-023 — Finite-size scaling of the avalanche distribution: is the corner magnitude set by
system geometry? If yes, Mmax is a geometric constant, not a statistical accident.**
- Claim: each region's tapered-GR corner magnitude Mc_corner scales with a geometric measure of
  its fault system (seismogenic area, maximum contiguous fault length from CFM/plate-boundary
  geometry, slab width from Slab2) with an exponent consistent with SOC finite-size scaling.
- Emergence logic: in a finite critical system, the avalanche distribution is a power law cut
  off by the system size — the cutoff is not a free parameter, it is geometry. If that holds
  across 13 regions, "how big can it get here" becomes a measurement of the fault system rather
  than an extrapolation of a catalog.
- Test: fit tapered-GR (Kagan) per region on GCMT Mw 1976-2026 (download; homogeneous magnitudes
  — essential, and it defuses the magType objection) plus ISC-GEM for the long tail. Geometric
  predictors: seismogenic area from the region box intersected with plate-boundary buffers,
  Slab2 down-dip width, and maximum contiguous mapped fault length. Statistic: log-log
  regression of Mc_corner (converted to corner moment) on the geometric measure; the exponent and
  its CI, tested against the SOC prediction (moment ~ L^3 for self-similar ruptures, ~L^2 for
  width-saturated ones — the two regimes are themselves a prediction: continental transform
  regions should show the saturated exponent, subduction the unsaturated one).
- Null: no relation (slope 0); and a shuffled-geometry control. Only 13 regions, so power is
  limited — pre-register that this is suggestive-only unless |t| is large, and extend to ~40
  Flinn-Engdahl regions from GCMT if the 13 look promising.
- Expected: a real slope; the interesting outcome is the regime split (exponent 3 vs 2) mapping
  onto width-saturated versus unsaturated systems.
- Why dismissed too quickly: "corner magnitude estimates are wildly uncertain with 31 years of
  data." True, hence GCMT+ISC-GEM (1904-) and hence the honest suggestive-only framing. The
  reason to do it anyway: it is the only entry in this round that produces a physically grounded
  regional Mmax, which is half of what "where is the next M>=7.5" means.

**K-024 — Synchronization: are fault segments phase-locked? A Kuramoto order parameter for a
fault ensemble.**
- Claim: neighbouring fault segments' loading cycles are partially synchronized (coupled
  oscillators), so the ensemble has a Kuramoto-style coherence R(t) that is elevated before
  multi-segment or large ruptures.
- Emergence logic: coupled-oscillator synchronization is the canonical emergent phenomenon.
  Faults are coupled through elastic stress transfer and viscoelastic relaxation; if coupling
  exceeds the spread of natural frequencies, they lock. Cascading multi-segment ruptures
  (Ridgecrest, Kaikoura, Turkey 2023) are what synchronized failure looks like.
- Test (on disk + one download): assign each CFM5.3 segment (557 on disk, with strike/dip/rake)
  a phase from time-since-last-event on that segment relative to its estimated mean recurrence
  (recurrence from the geodetic loading rate in socal_strain_grid.npz over segment area, via the
  moment budget). Order parameter R(t) = |mean(exp(i*phi_j))| over segments in a fault system.
  Statistic: R in the 5 yr before M>=6.5 SoCal events versus matched windows. Because SoCal has
  few large events, extend with the UCERF3/Third Uniform California Earthquake Rupture Forecast
  paleoseismic recurrence dataset (public) and with Japan (GCMT + segment models).
- Null: phases randomized within segments (destroys coherence, keeps marginals); and an
  ETAS-sim in which segments have no coupling.
- Expected: weak. Underpowered in SoCal (a handful of large events). Report as exploratory and
  primarily as an infrastructure build: the segment-phase object is reusable and is a component
  of the K-012 state vector.
- Why dismissed too quickly: it deserves the skepticism (few targets, model-heavy phases). Its
  claim on attention is that it is the only entry that treats faults as *interacting units with
  internal state* rather than as points emitting events, and that framing is a prerequisite for
  most physics-based forecasting.

**K-025 — Renormalization: does the seismicity field scale-collapse? A failure to collapse
locates a characteristic length, and that length's drift is an order parameter.**
- Claim: cell-count distributions at 0.05/0.1/0.2/0.4/0.8/1.6° do NOT collapse onto a single
  scaling function; the scale at which collapse fails is a physical characteristic length
  (fault-spacing? seismogenic thickness? asperity size?), it is measurable per region, and it
  changes with time.
- Emergence logic: coarse-graining is the fundamental probe of a many-body system. Pure
  criticality means no characteristic scale — the collapse is perfect. Every deviation is a
  physical length being revealed. This is a direct, cheap measurement of "at what scale does the
  crust stop being self-similar", which is Q2 asked in the language of statistical physics.
- Test (on disk): SoCal SCSN M>=2.5 (633k events, superb for this) and Japan/Indonesia. For each
  cell size L and a fixed time window, form the distribution of cell counts; attempt a
  finite-size scaling collapse P(n,L) = L^-a * F(n/L^b); fit (a,b) by minimizing collapse
  residual. Statistic: the collapse residual as a function of L — where it spikes is the
  characteristic length. Then track that length in 5-yr windows.
- Null: ETAS-sim (which HAS characteristic lengths: the spatial kernel scale d and the
  background smoothing) — so the null gives us the expected apparent characteristic length from
  the model alone, and the question is whether the Earth shows an additional one.
- Expected: a characteristic length near 10-20 km in SoCal (seismogenic thickness) that the
  ETAS-sim null does not produce. Time drift is speculative.
- Why dismissed too quickly: "fractal-dimension analyses of seismicity are a cottage industry
  that produced nothing predictive." Agreed about D2 papers. The difference: this is a *scaling
  collapse* (a much stronger test than fitting one exponent), the null is a simulation of our own
  validated model, and the deliverable is a length scale with an interpretation — which then
  becomes the grid resolution choice for K-002/K-009/K-022 rather than the arbitrary 0.2° we
  have been using since round 1.

**K-026 — ORDER PARAMETER #5: spatial entropy of activity. The system localizes before it
breaks.**
- Claim: the Shannon entropy of the spatial distribution of events in a trailing window
  (normalized by log(number of cells)) DROPS in the months before a large event — activity
  concentrates — beyond the drop that ETAS clustering alone produces.
- Emergence logic: the cheapest possible thermodynamic observable. Entropy is the canonical
  system-level scalar; localization is the geometric signature of an incipient rupture nucleus.
  If the crust has a "temperature" it should have an entropy, and this is the one-line
  operationalization.
- Test (on disk): per region, per month, H(t) = -sum p_i log p_i over 0.2° cells (SoCal) /
  0.5° (world), p_i = fraction of that month's events in cell i, computed above a fixed
  magnitude floor and with a fixed cell set. Statistic: change in H over the 6 months before
  M>=6.5 targets vs matched controls; AUC; and the combination with K-018's n(t) in a two-
  variable logistic model (do they carry independent information?).
- Null: ETAS-sim. Also an essential count-matched control: entropy is biased by sample size, so
  every comparison must be between windows with equal event counts (subsample to the smaller).
  That bias is how this kind of analysis usually fools people.
- Expected: AUC 0.55-0.62 alone. Its real value may be as the second axis in a 2-D
  (n, H) state-space plot of a region's trajectory — the seismic analogue of a phase diagram,
  which is a communicable artifact for Jim's app even if the AUC is modest.
- Why dismissed too quickly: "entropy of seismicity" has been published in low-impact venues
  with no nulls and no out-of-sample scoring, which has poisoned the term. The construct is
  sound; it needs this program's discipline, which is cheap to apply.

**K-027 — THE THERMODYNAMIC MOVE ITSELF: measure the predictability-vs-coarse-graining curve.
Prove that the aggregate is forecastable even though the event is not.**
- Claim: forecast skill is a smooth, measurable, increasing function of coarse-graining scale in
  space, time, and magnitude — and somewhere on that surface there is a product with genuinely
  useful skill (e.g. "regional M>=6 moment release over the next 90 days") even though the
  individual-event forecast is nearly worthless. The curve itself is the deliverable.
- Emergence logic: molecules are chaos; temperature obeys laws. The entire question of whether
  "seismic thermodynamics" exists reduces to: does predictability increase with coarse-graining,
  and how fast? Nobody in this program (or, as far as I know, in operational seismology) has
  drawn this curve. It reframes the whole enterprise from "can we predict earthquakes" (a
  yes/no that has burned the field for 60 years) to "at what scale does predictability begin"
  (a measurement with an answer).
- Test (100% on disk): define a grid of coarse-graining scales — space {0.2°, 1°, 5°, whole
  box}, time {1 d, 7 d, 30 d, 90 d, 365 d}, magnitude {>=4.5, >=5.5, >=6.5} — and for every cell
  of that grid compute out-of-sample forecast skill of the validated ETAS (frozen params,
  walk-forward) against a climatology baseline, using a proper score for the aggregate
  (CRPS for counts and for log moment release; plus the reliability/sharpness decomposition).
  Do it for SoCal and for the 6 EXP-M holdouts. Deliverable: a skill surface, plus its
  "predictability frontier" contour.
- Null: climatology (the region's own-period rate/moment distribution — a local oracle, the
  harsh baseline this program already uses). Additionally report the K-011 saturation time on
  the same axes so the two surfaces can be overlaid.
- Expected: near-zero CRPS gain at (0.2°, 1 d, M>=6.5); substantial gain (20-50% CRPS
  reduction) at (whole-box, 90 d, M>=4.5). The scientifically interesting quantity is the shape
  and slope of the transition, and whether the frontier sits at operationally useful scales
  (a county-sized region and a season would be useful; a continent and a decade would not).
- Why dismissed too quickly: "this is just a big scoring table, not a hypothesis." It is a
  measurement rather than a hypothesis test, and that is the point — it is the measurement that
  tells this program which hypotheses are worth generating next, and it is cheap because every
  component (frozen ETAS, walk-forward harness, catalogs) already exists and is validated. If
  only one thing from the emergence lens gets run, run this and K-022.

---

### Frame-breaking entries (added mid-round under the PRIME DIRECTIVE)

Everything above answers questions as asked. These five do not. Each attacks a presupposition
in the questions themselves — that an "earthquake" is a well-defined unit, that "prediction"
means events, that skill in bits is the currency, that our instruments are outside the system,
and that the program's deliverable is a forecast. All still terminate in falsifiable tests;
that is the one invariant.

**K-028 — Break the unit of analysis: "an earthquake" is a property of the detector, not of the
Earth. Make DETECTOR-INVARIANCE a precondition for any result graduating to BASELINE.**
- The presupposition attacked: every entry above, and every baseline B-1..B-5, treats the
  catalog's event list as the ground truth. It is not. It is the output of a detection and
  association pipeline with an amplitude threshold, a picker, a magnitude scale, and a network
  geometry that all changed over the study period. What we call "one M4.2" another pipeline
  calls a cascade of eleven; the QTM template-matched catalog on disk contains roughly an order
  of magnitude more events in the same crust and the same years as SCSN. This program has
  ALREADY been burned by exactly this: the Anza phase was +159° in SCSN and +25° in QTM
  (notes §14b) — same rock, same decade, different pipeline, opposite conclusion. That was
  recorded as a self-correction. I am proposing it be promoted to a *law of the program*.
- Claim (testable, and aimed at our own work): the validated baselines are NOT all
  detector-invariant, and the degree of invariance is measurable, differs sharply between them,
  and predicts which ones will replicate elsewhere.
- Test (100% on disk, no downloads): three independent pipelines covering overlapping SoCal
  space-time — `SCSN_original_catalog.txt` (633,667 events, 1981-2018),
  `QTM_12dev.txt` (Ross et al. template-matched, ~10x denser), `comcat_socal_m25.csv`
  (ComCat/ANSS, 2010-2026, different magnitude authority). For each baseline in turn:
  (a) B-1/B-2 ETAS — refit on each catalog at a common M0=2.5 and report the parameter vector
  and its dispersion across pipelines; then cross-score (params from catalog X, events from
  catalog Y) in bits/event. (b) B-3 the 7-day M>=5 rule — recompute the conditional probability
  in each. (c) B-4 the stress ledger — recompute chi and the SILENT list in each; report the
  Jaccard overlap of the silent set (EXP-K already used Jaccard for the H-sensitivity check, so
  the machinery exists).
- Statistic and pre-registered bar: define **invariance I = (cross-pipeline bits/event) /
  (own-pipeline bits/event)** for the model baselines, and Jaccard for the set-valued ones.
  Propose the frozen rule: a result graduates to BASELINE only at I >= 0.8 (or Jaccard >= 0.7)
  across at least two independent pipelines, with the exception explicitly recorded otherwise.
- Null: for the model baselines, the relevant reference is the *within-pipeline* bootstrap
  dispersion of the same statistic — cross-pipeline disagreement only counts if it exceeds
  within-pipeline sampling noise.
- Expected: ETAS shape parameters (p, alpha) largely invariant (they survived global transfer,
  so they should survive detector transfer); mu wildly non-invariant (it is a rate, and rate is
  the most detector-dependent quantity there is); the EXP-K silent list *substantially*
  non-invariant, because 79% of its cells were flagged detection-limited (n<20) by the worker
  — which means a large part of B-4's negative space may be the negative space of the NETWORK,
  not of the crust. That is an uncomfortable prediction about our own most celebrated result
  and it is exactly why the test should be run.
- Why this might be dismissed too quickly: "it's a robustness check, not science." No — it is a
  *demarcation criterion*. It says: a property that changes when you change the seismometers is
  a property of the seismometers. Adopting it would have saved the Anza episode prospectively
  rather than retrospectively, and it costs one script against files already on disk.

**K-029 — A lens with no name in seismology: PREDICTION IS COMPRESSION. Measure the bits of
predictability that remain, model-free, before proposing any more mechanisms.**
- The presupposition attacked: that the way to find the next signal is to think of a mechanism
  and test it. That is a search over an unbounded space with no stopping rule, and it is how
  60 years of earthquake prediction burned itself down. There is another way to ask the
  question — one that requires no mechanism at all.
- The lens: a model's log-likelihood *is* a codelength. Predicting a catalog and compressing a
  catalog are the same operation measured in the same unit. Therefore a *universal* compressor,
  which assumes nothing about physics, estimates an upper bound on ALL discoverable structure
  in the stream — including structure no one has thought of. The gap between what a universal
  compressor achieves and what ETAS achieves is the **remaining predictability budget**: the
  total number of bits per event still on the table for every future hypothesis in this program
  to compete over. Nobody has ever quoted that number for an earthquake catalog. It would
  reframe the entire enterprise from open-ended searching to spending a known budget.
- Claim: the residual budget is small but nonzero — and its *size* is the most decision-relevant
  quantity this program could produce.
- Test (on disk): tokenize the SoCal catalog (M>=2.5, 1981-2026) into a symbol stream with a
  frozen, reversible encoding: quantized log interevent time (16 levels), magnitude bin (0.1),
  and spatial cell index (0.2°, Hilbert-curve ordered so that nearby cells get nearby symbols).
  Encode the SAME stream three ways and compare codelengths in bits/event:
  (1) ETAS + GR + smoothed-seismicity as an arithmetic coder (this is exactly the frozen
  baseline's log-likelihood, converted to bits — no new modelling);
  (2) general-purpose compressors on the token stream: LZMA/xz -9, PPMd, and context-tree
  weighting with depth swept on train only;
  (3) a small autoregressive sequence model (a few-million-parameter transformer or an LSTM)
  trained on train (<2010) and evaluated as a next-token codelength on test — no cherry-picking,
  single frozen architecture chosen before seeing test.
- The control that makes this work (and that is the whole trick): run the identical three-way
  comparison on **ETAS-SIMULATED catalogs**, where ETAS is by construction the true generative
  model. There, any advantage of compressor (2)/(3) over (1) is pure estimator handicap — in
  fact it should be a deficit, and its size measures how much worse a generic learner is than
  the truth. Statistic: Delta = [bits(compressor) - bits(ETAS)] measured on real, minus the
  same quantity measured on simulated. Delta < 0 with bootstrap CI excluding 0 means the generic
  learner found structure in the real Earth that ETAS does not contain — and |Delta| is the
  remaining budget in bits/event.
- Null: Delta >= 0 (generic learners do no better on the Earth than they do on a world where
  ETAS is literally true). That null is exactly "ETAS is a sufficient statistic for the
  catalog", which is a crisp, important, currently-unmeasured claim.
- Expected: Delta between -0.05 and -0.4 bits/event. For calibration of what that means: B-2's
  entire validated achievement over Poisson is +1.87 bits/event. So a residual budget of 0.2
  bits would say roughly 10% as much structure remains as ETAS already captured — which would
  be a sobering, publishable, and genuinely novel statement, and it would tell Jim how much of
  the future to spend on catalog-only methods versus new instruments.
- Why this might be dismissed too quickly: "neural nets on earthquake catalogs is a crowded,
  disreputable genre." Agreed — because those papers report accuracy on unbalanced classes with
  no baseline. This is not that. Nothing here is a prediction claim; the deliverable is a
  *bound*, the baseline is a validated physical model, and the ETAS-simulated arm controls for
  the learner's own inadequacy. It is the least hypey possible use of the technique. Second
  objection: "a bound from one encoding is encoding-dependent." True — so report it for three
  frozen encodings and take the tightest; a bound is still a bound.

**K-030 — Counterfactual seismology: study the worlds that did not happen. Build a FRAGILITY
field and ship an information-value map instead of a hazard map.**
- The presupposition attacked: that the object of study is the history that occurred. We have a
  validated generative model, which means we have access to the ensemble of histories that
  *could* have occurred — and the variance across that ensemble is a physical quantity nobody
  maps. "Where is hazard high" and "where does knowing more change the answer" are different
  questions, and only the first has ever been asked here.
- The new object: **fragility** phi(x) = the sensitivity of the region's forecast to
  single-event perturbations of its history. Formally, for each cell, remove/add/shift one event
  in the recent history, re-run the forecast, measure the KL divergence of the resulting
  predictive distribution. High-fragility cells are the ones where a single marginal detection
  flips the forecast — i.e. where an additional seismometer buys the most, and where our own
  results are least trustworthy.
- Emergent link back to K-028: fragility and detector-invariance are the same quantity viewed
  from two sides. A cell whose forecast is fragile to one event is exactly a cell whose baseline
  will not survive a change of pipeline. If K-028 and K-030 agree cell-by-cell, that is a strong
  internal validation of both — and that cross-check is itself the test.
- Test (on disk): SoCal, 0.2° grid, 2010-2026. For each cell and each month, compute the 30-day
  predictive distribution from the frozen ETAS; then recompute it under a set of frozen
  perturbations (drop the largest event in the trailing 30 d; drop a random event; perturb all
  magnitudes by their reported magError; jitter times by 1 h). phi = mean KL divergence.
  Statistic 1 (descriptive but new): the map of phi and its relationship to event density and to
  the EXP-K silent list. Statistic 2 (falsifiable): **phi predicts disagreement.** Compute, for
  the same cells, the actual cross-pipeline forecast disagreement from K-028 (SCSN vs QTM vs
  ComCat); the claim is Spearman rho(phi, disagreement) > 0 with p < 0.01.
- Null: phi computed on ETAS-simulated histories (giving the expected fragility from sampling
  alone), and a cell-label permutation for the correlation test.
- Expected: rho 0.4-0.7. If it holds, phi becomes a cheap proxy for "how much should I trust
  this cell's forecast" computable anywhere in the world without needing three catalogs — which
  is a deployable per-cell confidence layer that the app currently has no way to produce.
- Why this might be dismissed too quickly: "that's just sensitivity analysis." Sensitivity
  analysis asks whether a conclusion is stable. This inverts it into a *product*: the map of
  where observation has leverage. A national network planner would rather have that map than a
  hazard map, because the hazard map does not tell them where to put the next station.

**K-031 — Turn the program on itself: the observer is INSIDE the system. Networks are deployed
where earthquakes happen, so seismicity causes its own detectability, and we have at least one
validated finding that may be a measurement of ourselves.**
- The presupposition attacked: that catalog completeness is an exogenous nuisance to be
  corrected. It is endogenous and causally downstream of the very process we are studying.
  Stations are installed after damaging events, near active faults, and in funded decades;
  therefore detected rate is partly a lagged function of past rate, which is a positive feedback
  loop that will manufacture apparent activation, apparent b-value drift, and apparent silence.
- The finding I am aiming at (deliberately, one of our own): **EXP-I(iii)'s b-value monotone
  decline**, which was frozen in train and confirmed in test (0.862, 0.854) with an explicit
  caveat that "magnitude-scale homogeneity across network eras" was unchecked. Also in the
  crosshairs: EXP-K's silent list, of which 79% of cells are detection-limited (n<20) — silent
  cells may be *unwatched* cells. A declining b and a growing network are exactly what you get
  if the small end of the catalog improved faster than the large end, or if ML/Mw revisions
  drifted.
- Claim: a materially large fraction of the observed b-decline (>=50%) and of the unexplained-
  silent list (>=30% of cells) is attributable to network evolution rather than to the crust.
- Test: (i) Download SCEDC/IRIS FDSN station metadata for the SoCal box with operational start
  and end dates (fdsnws/station, one query, public) and build a station-density field
  rho_sta(x, t). (ii) Recompute the b-value trend three ways: on the full catalog (the original
  result), on a **fixed subnetwork** (only events recorded by stations operating for the entire
  1981-2018 span — a constant-instrument catalog), and on a magnitude floor set far above the
  worst-era Mc (M>=3.0 rather than 2.5). Statistic: does the frozen-sign monotone decline
  survive all three? Report the trend slopes with CIs side by side. (iii) For the ledger:
  regress the silent-cell indicator on rho_sta(x, t_train) — how much of the silent geography is
  predicted by where the seismometers weren't?
- Null: no attenuation of the trend under the fixed-subnetwork restriction; no association
  between silence and station density (permutation of cell labels within loading-matched
  strata).
- Expected: I expect the b-decline to attenuate substantially and the silent list to lose its
  detection-limited tail — leaving the 42 measured-low-chi cells (including the 1857 Big Bend
  strand cells, which have n_train = 2 and 23 respectively, so even those are thin) as the real
  residue. If the 1857-strand cells survive a station-density control, B-4 becomes considerably
  stronger than it is today. Either way the program learns which of its findings are about rock.
- Why this might be dismissed too quickly: "we already caveated the homogeneity issue." A
  caveat is not a measurement. This converts two standing caveats into numbers, and it is the
  single fastest way to either harden or retire two current claims. A research program that
  cannot attack its own baselines is not a research program, it is an advocacy campaign.

**K-032 — Answering a question nobody posed: redefine what would count as SUCCESS. The program's
durable deliverable is not a forecast — it is the PREDICTABILITY BUDGET.**
- The presupposition attacked: that this program succeeds if and only if it predicts an
  earthquake. Under that definition, sixty years of the field have "failed", the failures are
  unpublishable, and each generation repeats them. That definition is the field's central
  pathology and Jim's program is small enough and disciplined enough to escape it.
- The reframe: this program is unusually good at one thing that almost nobody does — killing
  hypotheses cleanly, with frozen protocols, honest nulls, public failures, and self-corrections
  recorded in the record (the Anza episode, the sign-error catch, the border-artifact check).
  That is not a consolation prize. It is the scarce good. Propose the deliverable be a
  quantitative, public, updateable **Predictability Budget for earthquake forecasting**:
  1. the measured predictability frontier over coarse-graining scales (K-027);
  2. the total bits/event currently captured (B-1/B-2: +1.87 SoCal, +0.66..+1.75 globally);
  3. the model-free bits/event still available (K-029);
  4. the predictability horizon in days as a function of scale (K-011);
  5. the detector-invariance of every claim (K-028);
  6. the corpse list with its effect-size upper bounds — "tidal phase maps: |skill| < 0.05
     bits/event at SoCal power" is a far more useful sentence than "tidal triggering is null";
  7. the decision-value translation: for a given exposure and loss model, the skill THRESHOLD a
     forecast must clear to be worth issuing — so that every future hypothesis has a bar in
     units that matter rather than in p-values.
- Is this falsifiable? Yes, as a whole and in parts, which is why it belongs in this ledger and
  not in a blog post. Each line is a measurement with a stated uncertainty, and the budget makes
  a hard, checkable meta-prediction: **that future validated additions to this program will sum
  to less than the residual budget measured in K-029.** If someone later demonstrates skill
  exceeding the bound, the budget is falsified and we learn our bound was encoding-limited —
  which is itself the most interesting possible outcome. A framework that can be broken by
  future success is a scientific object.
- The concrete test to run FIRST, cheaply: item 6, retroactively. Take every corpse on the
  current list and convert each from a null into an **upper bound with a confidence interval**
  (what is the largest effect the existing test would have detected at 80% power?). Data and
  code all on disk (results_exp_a/b/c/e/f/g.json plus the harnesses). Deliverable: a table of
  bounded effect sizes. Statistic: per-corpse minimum detectable effect at 80% power under the
  original frozen design.
- Why this might be dismissed too quickly: "this is meta-work, not science; Jim asked how to
  predict earthquakes." Jim asked how to predict earthquakes *and* whether compute plus data
  would crack it. Item 3 plus item 1 IS the answer to the second question, in numbers, and it
  is achievable this month with data already on disk. And the reframe is not a retreat: the
  fastest route to a real forecast is knowing exactly how many bits are left and where they are.

---

## CHARTER AMENDMENT PROPOSED (Kepler → supervisor + Jim)

Offered under the PRIME DIRECTIVE's "Evolve yourself" clause. These are proposals; the
supervisor and Jim decide. Each is phrased as a line that could be pasted into
`.claude/agents/kepler.md`.

1. **Standing null upgrade (the highest-value change).** Add to Ground Rules: *"The default null
   is a simulation from the program's own validated generative model, never a Poisson process.
   Any entry whose null is Poisson must justify it explicitly."* Rationale: roughly two-thirds
   of the classical precursor literature dies to this substitution, and this program uniquely
   owns the validated generators (B-1/B-2) needed to apply it. It is the cheapest available
   increase in the truth-yield of everything I generate.

2. **Detector-invariance as a graduation criterion.** Add to Ground Rules: *"No result graduates
   to BASELINE without demonstrating invariance across at least two independent detection
   pipelines, or an explicit recorded exemption."* Rationale: K-028; the Anza episode was this
   lesson learned the expensive way, retrospectively.

3. **Adversarial self-audit quota.** Add to the posture: *"Each round, at least one entry must
   attempt to FALSIFY one of our own baselines, named explicitly."* Rationale: Popper adjudicates
   new claims but nothing currently re-examines old ones; baselines decay silently. K-031 is this
   round's instance and I would not have written it without being pushed to be emergent — which
   suggests it should be structural rather than occasional.

4. **Bits-per-compute-hour ranking.** Add to the output format: *"Each entry carries an estimated
   cost (compute-hours, downloads, person-time) and an expected-bits payoff, so the supervisor
   can rank by expected yield rather than by my enthusiasm."* Rationale: I produce more ideas
   than can be run; the prioritization currently falls on the supervisor with no cost information
   from the person who knows the design.

5. **Name the decision, not just the statistic.** Add to Ground Rules: *"Every entry states what
   decision would change if it were true. Skill in bits is not skill in value."* Rationale: it
   would have caught, early, that some of round 1's tidal work could not have changed any action
   even had it succeeded at the measured effect size.

6. **The lens ledger (meta-learning on idea generation).** Add: *"Tag every entry with the lens
   that generated it, and maintain a running tally of lens → outcome. Lenses with repeated
   validated yield get more of the round; lenses with repeated corpses get retired; each round
   must field at least one lens that has never appeared in the ledger before."* Rationale: this
   is the persona file's own emergence instruction turned into a measurement. It makes my
   improvement legible rather than assumed, and it means my successor inherits data about what
   worked, not just a list of hypotheses.

7. **The right to propose retirement.** Add to the role: *"Kepler may propose moving a BASELINE
   back to PROPOSED, with evidence; Popper adjudicates that as he would a new claim."* Rationale:
   currently the ledger has a one-way ratchet from PROPOSED to BASELINE. One-way ratchets are how
   research programs acquire their unfalsifiable core.

8. **A standing budget line.** Add: *"Maintain the Predictability Budget (K-032) as a living
   section of this ledger; every validated result and every corpse updates it."* Rationale: it
   converts the program's accumulated negative results — currently its largest asset by volume —
   into a cumulative, publishable quantity instead of a list of disappointments.

---

### Conditional-triggering seed (Jim)

Jim's reframe, which I think is correct and which I should have made myself: **every corpse on
our list is a MARGINAL effect.** Tidal phase alone. Periodicity alone. Amplitude alone. We have
never once tested a *conditional* effect. A match does nothing to a wet log and everything to a
dry one, and we have been measuring the average effect of matches on logs.

Three structural consequences, which govern every entry below:

**(a) The statistical object changes.** The right model is a conditional-intensity point process
with covariates and interactions:
`λ(t,x) = [ μ(x)·R(t,x) + Σ_{t_i<t} triggering ] · exp( Σ_k β_k z_k(t,x) + Σ_{j<k} γ_jk z_j z_k )`
where z are external-load covariates and the γ are the interaction terms Jim is asking about.
This is ETAS-with-covariates (a Cox-ETAS). Critically, **the ETAS terms sit in the baseline**,
which means the sequence-coherence artifact that destroyed the original TSI paper — the founding
lesson of this whole program — is *structurally* removed rather than avoided by hand. Round 1
could not test conditional effects safely. It can now.

**(b) The non-firing is already in the likelihood, and it is half of it.** Jim asks whether
windows where loads were maximal and nothing happened are informative. The answer is sharper
than "yes": the point-process log-likelihood is `Σ_i log λ(t_i) − ∫ λ(t) dt`. The first term is
the events; **the second term is the silence, and it carries comparable information.** Every
corpse-era test we ran (phase histograms, Rayleigh statistics, Schuster tests) used the first
term only. That is not a small methodological difference — see K-035, which measures exactly how
much power we threw away. Jim's negative-space intuition is, mathematically, an argument that we
have been analyzing our data at roughly half efficiency.

**(c) We now have a KNOWN-TRUE positive control, and it gates everything.** Remote dynamic
triggering from Landers 1992 is documented, replicated science. It means the covariate machinery
can be *calibrated* rather than merely used: if the engine recovers Landers at the right
amplitude, its nulls elsewhere become quantitative upper bounds instead of shrugs. **I propose
K-034 be run FIRST and that no null from this family be interpreted until it passes.**

Standing design discipline for this whole family (stated once; assumed by every entry below):
covariate list and interaction grid **enumerated in the frozen protocol before any fitting**;
selection on train only (group-lasso or LRT scan with Benjamini-Hochberg, procedure frozen);
**one** test-window unblinding scored in bits/event; ETAS-sim null throughout; and every
interaction reported with its marginal alongside it, so we can always see whether the
conditional effect is real or is a marginal effect in disguise.

Note on Laplace: these specs assume a measurement engine that can grind full covariate sweeps
with interactions. That capability is a genuine unlock — but it is also the fastest known route
to a forking-paths disaster, which is why the discipline above is stated before the ideas.

**K-033 — THE FRAMEWORK ENTRY: Cox-ETAS. One conditional-intensity model, a frozen covariate
set, a frozen interaction grid, scored in bits/event out of sample.**
- Claim: external loads contribute measurable forecasting skill *only* through interaction terms
  with the crust's state; the marginal coefficients β will be indistinguishable from zero
  (consistent with every corpse) while at least one interaction γ will not.
- Inversion: stop asking "does the tide trigger earthquakes" (a question about an average) and
  ask "for which crust, in which state, at which moment" (a question about a conditional). The
  average over a heterogeneous crust of a state-dependent effect is exactly zero-ish — which is
  what we measured, four times, correctly.
- Covariate set z (frozen, all resolvable at hourly-to-daily resolution):
  *external loads* — (1) tidal Coulomb stress on the local fault plane (on disk:
  `data/xue_lu_zenodo/Tidal_{N_0,N_90,S_0,Vol}.txt`, 6000-s from 1981, ×1e-9 strain, via
  `calc_stress.py`, E=75 GPa ν=0.25, resolved with CFM5.3 strikes and the `coso_fm_test.py`
  machinery that is already validated); (2) tidal stressing RATE (dτ/dt — a different physical
  quantity from phase, and the one rate-state theory actually cares about); (3) atmospheric
  surface pressure (ERA5, download); (4) hydrologic load (ERA5-Land soil moisture + snow water
  equivalent; GRACE/GRACE-FO mascons 2002-, download); (5) remote dynamic stress proxy (K-038).
  *state variables* — (6) EXP-J/K ledger class and log-chi (on disk, `results_exp_j/k.json`);
  (7) branching ratio n(t) from K-018; (8) recent-rate regime (trailing 30-d rate / long-term
  mean); (9) depth; (10) distance to geothermal (`tsi_map.GEOTHERMAL`); (11) local b-value.
- Interaction grid (frozen, 5 loads × 6 states = 30 terms, not a free-form search): every
  load × every state, plus load×load for the combination question (K-036/K-042).
- Test: SoCal M>=2.5, train <2010, walk-forward test 2010-2026 (SCSN to 2018 + ComCat 2010-2026,
  with a K-028 invariance check across the seam). Statistic: bits/event over the frozen B-2
  baseline (μ, K, α, c, p held at the EXP-H values so the covariate terms cannot absorb ETAS
  misfit); per-term LRT with BH-FDR at q=0.05.
- Null: ETAS-sim catalogs with the covariates left in place (so the covariates are real but the
  crust is indifferent to them) — this gives the exact distribution of spurious γ under
  "loads do nothing", including all the correlation the covariates have with each other and with
  seasonality.
- Expected: marginal β ≈ 0 (|effect| < 0.02 bits/event); best interaction +0.02 to +0.15
  bits/event. Small in absolute terms — but a 1.05-1.3× rate multiplier that switches on
  conditionally is operationally different from a 1.02× multiplier that is always on.
- Why dismissed too quickly: "tidal triggering is settled null in California." It is settled
  null *marginally*, and we ourselves produced some of the cleanest evidence for that. Nobody
  has published the conditional version with a ledger-derived state variable, because nobody
  else has a validated stress ledger to condition on. B-4 is the asset that makes this test
  possible and it did not exist a week ago.

**K-034 — RUN THIS FIRST: Landers 1992 as a known-true positive control that licenses every
subsequent null.**
- Claim: the Cox-ETAS engine, applied blind to a remote-dynamic-stress covariate, recovers the
  documented Landers-triggered rate increases (Long Valley, Coso, Geysers, Cedar City, western
  Nevada) at the amplitude and timing the literature reports — and if it does not, the engine is
  not yet trustworthy and no null from K-035..K-045 may be interpreted.
- Inversion: the deepest problem with this program's null results is that we cannot distinguish
  "no effect" from "no power". Jim has handed us the fix: a triggering phenomenon that is
  *known to be real*, in *our own catalog*, in *our own box*, with a well-documented amplitude.
  That converts our detection threshold from an assumption into a measurement. This is the same
  move as EXP-F's 7-day artifact control (which, notably, FAILED to fire — and that failure
  correctly downgraded the whole periodicity null to "weak evidence"). Do it deliberately this
  time.
- Test (mostly on disk): SCSN original catalog covers 1992 — Landers (M7.3, 1992-06-28) and Big
  Bear are in it, as are Hector Mine 1999 and Ridgecrest 2019 (ComCat). Build the dynamic-stress
  covariate as in K-038 and fit the response with NO knowledge of the literature values; then
  compare the fitted response to published Landers triggering distances/amplitudes. Statistic:
  (i) does the dynamic-stress term reach significance at all; (ii) the implied triggering
  threshold in peak dynamic stress (literature: order 0.01-0.1 MPa PGV-equivalent); (iii) spatial
  pattern of response — geothermal/volcanic areas should light up first, which is a strong
  pattern prediction, not a single number.
- Null: circular time-shift of the trigger catalog against the target catalog (this program's
  established null machinery), plus the ETAS-sim null so that ordinary aftershocks of Landers
  inside SoCal are not counted as "remote" triggering. The distance gate matters: score only
  targets beyond 2 rupture lengths.
- Expected: clear detection for Landers-class events; the deliverable is the **minimum detectable
  dynamic stress at our power**, which becomes the denominator for interpreting every subsequent
  null in this family.
- Why dismissed too quickly: "it's a positive control, not a discovery." Correct, and it is the
  most valuable non-discovery available to us. Without it, K-035..K-045 produce uninterpretable
  nulls; with it, they produce bounded upper limits, which is exactly the currency K-032 needs.

**K-035 — POWER AUDIT OF THE CORPSES: how much did we lose by using event-only statistics
instead of the full likelihood? Possibly the corpses are not dead but underpowered — and we can
say by how much, to two significant figures.**
- Claim: the intensity-likelihood formulation detects a given tidal modulation at materially
  smaller amplitude than the phase-histogram/Schuster formulation used in EXP-A — plausibly a
  factor of 1.5-3× in minimum detectable effect — because it uses the ∫λdt term (the silence)
  and because it conditions on ETAS rather than comparing to uniform.
- Inversion: Jim asked whether non-firing is informative. Rather than argue, *measure the value
  of the information*. This is a pure power calculation, requires no new data, and it retroactively
  re-prices every null this program has produced.
- Test (100% on disk, no downloads, fast): injection-recovery. Generate synthetic SoCal catalogs
  from the frozen EXP-H ETAS with a *known* tidal modulation injected into the intensity,
  amplitude a ∈ {0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.40}, using the real tidal series so the
  sampling and aliasing are authentic. Analyze each synthetic catalog TWICE: (i) with the exact
  EXP-A pipeline (per-bin phase selection, pooled S, circular-shift null); (ii) with the
  Cox-ETAS likelihood of K-033. Statistic: detection probability vs a for each method →
  the minimum detectable amplitude at 80% power for each.
- Null: a = 0 arm confirms both methods hold their false-positive rate (and doubles as a check
  on the new machinery, which is the safer reason to run it first).
- Expected: EXP-A's minimum detectable a ≈ 0.10-0.20 at SoCal FM-matched sample sizes; Cox-ETAS
  ≈ 0.03-0.08. If so, then the true statement about round 1 is not "no tidal effect" but **"no
  tidal effect larger than ~15%, and we could not have seen a 5% one"** — and rate-state theory
  (K-036) predicts an effect of a few percent. That would mean our headline null and the leading
  physical prediction have never actually been in contact with each other.
- Why dismissed too quickly: "you are trying to resurrect a result you already killed." I am
  trying to put a number on the coffin. If the audit shows EXP-A could have detected a 3% effect
  and did not, the corpse is *more* dead than it was, and I will say so — that is the honest
  symmetry of this test, and it is why it is worth running regardless of which way it lands.

**K-036 — Combinations of loads, done with the right link function: SUM the Coulomb stressing
rates and let rate-and-state supply the response curve. One physical parameter, Aσ.**
- Claim: seismicity responds not to any single load but to the *summed* Coulomb stressing rate
  from all sources, through the rate-state response `R = exp(Δτ/Aσ)`, with a single fitted
  parameter Aσ ≈ 0.01-0.1 MPa — and fitting the sum beats fitting any component, and beats
  fitting them as independent additive terms.
- Inversion/analogy: Jim's "combinations of loads" is exactly right, and it has a rigorous form
  already waiting for it. Dieterich's rate-state theory *predicts the functional shape* of the
  response to a stress perturbation, so this is not curve-fitting — it is a one-parameter
  physical model with a strong prior, and Aσ has an independent physical meaning (effective
  normal stress × friction parameter). Better still, it predicts *where* sensitivity should be
  highest: where Aσ is small, i.e. low effective normal stress — fluid-rich, geothermal, shallow.
  That is a spatial pattern prediction we can test, and it is the same map B-4 already produces.
- Test: build the total stressing series per 0.2° cell: solid-earth tide (on disk) + ocean
  loading (SPOTL/TPXO, download) + atmospheric pressure (ERA5) + hydrologic (ERA5-Land + GRACE)
  + thermoelastic, each resolved onto the local CFM strike/dip/rake (machinery on disk). Fit
  Cox-ETAS with the rate-state link and Aσ free. Compare, by out-of-sample bits/event: (a) sum
  with rate-state link; (b) each component alone; (c) all components as free additive linear
  terms; (d) no load terms (= B-2).
- Statistic: nested model comparison in bits/event on the test window; the fitted Aσ with CI;
  and the spatial map of fitted Aσ tested against the geothermal/fluid-rich prior (rank
  correlation with distance-to-geothermal and with B-4's chi).
- Null: ETAS-sim with loads present but inert; plus independent circular shifts of each load
  series (K-042's null), which preserves each load's own statistics and destroys only their
  co-occurrence.
- Expected: Aσ in the 0.01-0.1 MPa range if anything is there; +0.02 to +0.10 bits/event for the
  summed model over the best single component. The strongest possible outcome is not the bits —
  it is a *map of Aσ* that matches the independently-derived fluid-rich map, because two
  unrelated observables agreeing is far harder to fake than one p-value.
- Why dismissed too quickly: "rate-state tidal predictions have been tested and are marginal."
  They have been tested marginally, on single load components, without conditioning on state, and
  with event-only statistics whose power we are about to measure in K-035. The novel content here
  is the *summation* and the *spatial prediction of Aσ*, neither of which is standard.

**K-037 — Hydrologic/atmospheric loading, with the killer falsifier: the SIGN is predicted per
fault from geometry, so a spurious correlation cannot pass.**
- Claim: seasonal hydrologic loading modulates seismicity with a sign that flips between faults
  according to their orientation and sense of slip relative to the load — and the observed
  per-fault modulation phase matches the geometrically predicted sign significantly more often
  than chance.
- Inversion: every seasonality correlation is drowning in confounders (temperature-dependent
  detection, seasonal noise, holiday quarry blasting, snow on stations). But a confounder cannot
  know a fault's rake. Predict the sign *per fault* from geometry first, then test the agreement
  pattern — this converts a weak correlational claim into a strong pattern-matching one. It is
  the same design logic that made the fault-resolved tidal work meaningful, applied to weather.
- Test: compute seasonal surface load (ERA5-Land SWE + soil moisture + ERA5 pressure; GRACE
  mascons as an independent cross-check post-2002); convert to Coulomb stress change on each
  CFM5.3 segment via an elastic half-space Green's function (Boussinesq/Farrell loading — the
  standard, well-documented computation); the predicted sign is then fixed by strike/dip/rake,
  which are in the file on disk. Assign each SoCal event to its nearest segment; compute the
  observed annual phase of seismicity per segment; statistic = fraction of segments (weighted by
  n) whose observed phase falls in the predicted half-cycle, binomial test against 0.5.
- Null: binomial 0.5; plus a rake-shuffled control (permute the predicted signs across segments)
  which destroys the geometry link and keeps everything else; plus the ETAS-sim null.
- Expected: 55-65% sign agreement if real, on maybe 100-200 usable segments. Also run it in
  reverse as a diagnostic: if agreement is ~50% but a strong *marginal* annual signal exists,
  that is positive evidence the annual signal is an observational artifact — a useful result for
  EXP-F's ambiguous periodicity null.
- Why dismissed too quickly: "seasonal seismicity is a detection artifact." Possibly — and the
  sign test is precisely the design that can tell the two apart, which is why it is worth doing
  rather than assuming. Note this hypothesis has real published support in California (Central
  Valley unloading, seasonal M>=4 modulation) which we should treat as a weak positive control,
  not as a conclusion.

**K-038 — Remote dynamic triggering as a CONTINUOUS covariate, not an event study — and its
interaction with ledger class is the real hypothesis.**
- Claim: build a continuous "dynamic-stress weather" series for every SoCal cell from all global
  earthquakes; its marginal effect on rate will be small, but its interaction with the crust's
  state (B-4 ledger class, geothermal proximity, recent-rate regime) will be substantial —
  susceptibility to remote triggering is spatially organized and predictable from the ledger.
- Inversion: the literature does remote triggering as a series of case studies (one big
  earthquake, one target region, one window). Invert it into a *field*: the crust is being shaken
  constantly by the whole planet, and that continuous exposure is a covariate like any other.
  Case studies cannot estimate an interaction; a continuous covariate over 45 years can.
- Data: on disk, `data/comcat_world/*.csv` already contains **179 M>=7.0 and 517 M>=6.5**
  (1995-2026) — but for a proper ping series we want global coverage and a longer span: one
  ComCat FDSN query for worldwide M>=5.5, 1980-now (~30k events). Proxy: peak dynamic stress
  ∝ PGV via a standard attenuation relation `log PGV = a + bM − c·log r` (frozen coefficients,
  not fitted), arrival time = r / 3.5 km s⁻¹ for surface waves. Jim's `10^M/r` is the right
  instinct; the attenuation form is its calibrated version and costs nothing extra.
- Test: add the PGV covariate and its interactions to K-033. Statistic: bits/event; the fitted
  triggering threshold; and the map of per-cell susceptibility. Key secondary test: is
  susceptibility *predicted* by the ledger (silent-loading vs coupled) out of sample?
- Null: circular shift of the global trigger catalog relative to the SoCal target catalog;
  ETAS-sim; and a crucial exclusion of the trigger's own aftershock zone.
- Expected: threshold-like behaviour with detection above ~0.05 MPa; interaction with geothermal
  proximity strong (this is the best-documented pattern in the literature and doubles as a second
  positive control); interaction with ledger class is the genuinely new claim, and I would guess
  a 1.5-3× susceptibility ratio between classes.
- Why dismissed too quickly: "remote triggering is known and only happens in geothermal areas."
  If that is true, the ledger interaction will be null and the geothermal one positive — which is
  a clean, informative outcome that also validates the machinery. The hypothesis worth money is
  that "geothermal" is a proxy for "low Aσ / near-critical", and the ledger measures that
  directly and everywhere, including where there is no geothermal field.

**K-039 — The dry-log test proper: tidal sensitivity is nonzero ONLY in near-critical cells and
near-critical times.**
- Claim: the tidal coupling coefficient, estimated conditionally, is significantly larger in
  (a) B-4 silent-loading/late-cycle cells, (b) cells with low measured Aσ from K-036, (c)
  periods of elevated n(t) (K-018), and (d) the final weeks before a large local event — while
  remaining ~0 on average, exactly as EXP-A found.
- Inversion: this is Jim's match-and-log stated in our own variables. The corpse was the average
  over all logs. And note the internal consistency requirement: if K-039 is true, then EXP-A's
  clean null is not an embarrassment but a *prediction* of K-039 — the marginal must be ~0 when a
  state-dependent effect is averaged over a heterogeneous crust.
- Test: interaction terms in K-033, with the four state variables above. Because (d) is the
  operationally interesting one and the most vulnerable to hindsight, it gets its own frozen
  design: for each M>=5 SoCal mainshock, estimate the tidal coupling in the 90 d before it (using
  M>=2.5 events in a 50 km radius, intensity-likelihood so Omori/foreshocks are in the baseline)
  versus the same cell's long-term coupling. Paired statistic across mainshocks.
- Null: ETAS-sim (which has no tidal coupling at all, so any apparent pre-mainshock coupling
  there is pure noise); plus a matched-window control that equalizes event counts, since coupling
  estimates are noisier when n is small and pre-mainshock windows are busier — the exact
  n-bias that killed EXP-B. Stated explicitly so it cannot recur.
- Expected: modulation of 5-15% in the dry-log subset versus 0-2% overall. If (d) holds, it is a
  short-term precursor with a physical mechanism, which would be the most consequential result
  this program could produce — and correspondingly the one I would demand replicate on an
  independent region (Japan, where the tidal signal is documented as stronger and where Lian
  Xue's group has already invited collaboration) before anyone says a word in public.
- Why dismissed too quickly: "you already killed this." I killed the marginal, with good
  controls, and I will kill the conditional too if it dies — but K-035 may well show EXP-A could
  not have seen an effect of the size rate-state predicts, and a state-conditioned effect is a
  different hypothesis with a different null.

**K-040 — Go where the counts are: big aftershock sequences are the highest-SNR laboratory on
disk for ANY small-stress effect, in both directions (enhancement AND suppression).**
- Claim: inside the Landers 1992, Hector Mine 1999, and Ridgecrest 2019 sequences, where rates
  run 10-100× background, small stress perturbations are detectable at amplitudes an order of
  magnitude below what background seismicity permits — and suppression during unfavourable phase
  should be as detectable as enhancement during favourable phase.
- Inversion: power in a point process scales with the number of events. We have been hunting a
  1-5% modulation in background seismicity, which is the lowest-count regime available, while
  sitting on three sequences containing tens of thousands of events each. Also invert the sign:
  everyone looks for triggering; *clamping* is equally predicted by the physics and is easier to
  see against a high rate. Nobody looks for the hole.
- The trap, named explicitly: this is precisely where the sequence-coherence artifact lives — it
  is the artifact that killed the original TSI paper (24-h aftershock windows inheriting the
  mainshock's tidal-window membership; effective n = number of mainshocks). The intensity-model
  formulation is what makes this safe now: Omori decay is in the baseline term, so the modulation
  coefficient is estimated *relative to the decaying rate*, and pseudo-replication cannot inflate
  it. This entry is only runnable because of that; it would have been indefensible in round 1.
- Test: fit Cox-ETAS restricted to each sequence (M>=2.0 within 50 km, 0-365 d), with tidal
  Coulomb stress on the dominant local FM orientation (the FM catalogs on disk cover these) as
  the covariate, plus a stressing-rate term. Statistic: coupling coefficient with CI, per
  sequence; consistency of sign and phase across the three sequences (an independent-replication
  test built into the design); and separately the suppression amplitude.
- Null: ETAS-sim of each sequence with the fitted Omori parameters and no tidal coupling — this
  is the control that would have caught the original artifact, and running it here retires that
  ghost properly.
- Expected: minimum detectable modulation ~1-3% in these sequences versus ~10-20% in background
  (quantified by K-035). If a real effect exists anywhere in California at rate-state amplitudes,
  this is the place it becomes visible.
- Why dismissed too quickly: "aftershock tidal studies are contaminated." They were — by the
  specific artifact this program diagnosed and named. We are the group best positioned to redo
  them correctly, and it would be strange not to use that.

**K-041 — THE NEGATIVE-SPACE ENTRY: the didn't-fire ledger. Loads maximal, cell near-critical,
nothing happened. Is silence-under-load a strength gauge?**
- Claim: cells that repeatedly fail to respond to their own combined-load maxima are measurably
  stronger/further from failure than cells that respond — and that non-response history predicts,
  out of sample, both a lower short-term rate and a *larger* eventual event when they do go.
- Inversion: this is B-4's move (model where the system says events should fire and don't) moved
  from the 29-year secular timescale to the hourly load timescale. B-4's negative space was
  geodetic; this one is dynamic, and it refreshes every tidal cycle instead of every decade.
- The formal point worth stating for Jim: non-firing is not a gap in the data, it is the
  `∫λ(t)dt` term of the likelihood. A cell that sat through 500 load maxima without firing has
  contributed a large, precise, negative constraint on its own susceptibility. The right tool is
  a survival model with time-varying covariates (equivalently, the same point-process
  likelihood), and the "didn't fire" windows are not censoring to be discarded — they are the
  measurement.
- Test: define per cell a **responsiveness index** ρ = fitted load-coupling coefficient from
  K-033/K-036 restricted to that cell, estimated on train only. Partition cells into responsive
  / unresponsive at a frozen threshold, matched on event count and loading rate (the matching is
  essential — unresponsive will otherwise just mean low-n). Then three frozen out-of-sample
  tests: (i) does ρ persist train→test (Spearman)? [this is the EXP-J persistence question asked
  of a different quantity, and EXP-J's persistence was NULL, so this is a genuine test, not a
  formality]; (ii) do unresponsive cells have lower test-period rates?; (iii) **the payoff** — is
  the maximum magnitude in the test period larger in unresponsive cells than responsive ones
  (Mann-Whitney on per-cell Mmax, count-matched)?
- Null: ETAS-sim, where ρ is pure noise and therefore must not persist; plus the count-matched
  permutation.
- Expected: honest prior is that (i) fails, as EXP-J's persistence did — chi at cell scale was
  transient-dominated, and ρ may be too. But (iii) is the one worth the compute: "quiet under
  load" as a predictor of *size* rather than *timing* is a hypothesis nobody has posed in this
  form, and it is the natural reading of the stored-strength picture.
- Why dismissed too quickly: "absence of evidence." No — this is evidence of absence, formalized:
  the number of load maxima survived is a measured quantity with a computable likelihood
  contribution, and the entry stands or falls on frozen out-of-sample tests.

**K-042 — Combinations means JOINT EXTREMES, not average effects. Test the top 0.1% of hours,
with a null that destroys co-occurrence and preserves everything else.**
- Claim: the interesting combined-load events are rare joint maxima — spring tide at perigee,
  plus a deep atmospheric low, plus peak snowmelt, plus a teleseismic surface wave, aligned within
  hours — and seismicity in those few hundred hours per decade is elevated beyond what any
  additive model of the components predicts.
- Inversion: regression estimates a mean response and is structurally blind to a threshold that
  is crossed only in the tail. If the crust has a failure threshold, the physics lives entirely
  in the tail and the mean response is the wrong statistic. Import extreme-value theory:
  peaks-over-threshold on the *combined* load, and test the seismicity conditional on being in
  the tail.
- Test: form the summed Coulomb stress series per cell (K-036); identify the top 0.1% and top
  0.01% of hours by combined amplitude; statistic = observed event count in those hours versus
  the ETAS-expected count (so clustering is accounted for), as a rate ratio with a Poisson CI.
  Also fit the tail-dependence structure between load components (do they co-occur more than
  independence implies? tides and weather are genuinely independent; snowmelt and pressure are
  not) — that is a required preliminary, because apparent combination effects can come from the
  loads themselves being correlated.
- Null (this is the crux and the reason the entry is worth running): **independent circular time
  shifts of each load series.** This preserves every component's own amplitude distribution,
  autocorrelation, seasonality, and marginal relationship to seismicity, and destroys *only* their
  co-occurrence. Any excess that survives is attributable to combination and nothing else. It is
  the cleanest null in this entire round.
- Expected: rate ratio 1.1-1.5 in the top 0.1%, if the threshold picture is right. Sample size is
  the limiting factor — ~400 hours per decade per cell — so this must be pooled across cells with
  the intensity model doing the pooling.
- Why dismissed too quickly: "you are testing 0.1% of the data, that is cherry-picking." The
  threshold is defined on the *covariate* alone, with no reference to seismicity, and frozen
  before unblinding — which makes it a pre-specified subgroup, not a cherry-pick. The distinction
  is the whole methodology of extreme-value statistics.

**K-043 — FRAME-BREAK: stop passively observing. Every global M>=7 pings the whole planet — the
Earth is running a free active-source experiment and we can read local criticality off the
response. A global PING MAP.**
- The presupposition attacked: that we are limited to observing what the crust does
  spontaneously. We are not. Roughly 15 M>=7 events per year send measurable dynamic stresses
  through every fault system on Earth. That is an *active-source stress experiment*, run
  continuously since 1980, for free, and this program — like most of the field — has been
  treating it as noise or as case studies.
- The lens: seismology as **spectroscopy**. You learn a material's state by driving it and
  measuring the response, not by watching it sit. The response amplitude to a known stress
  transient is a direct readout of how close a region is to failure — and, crucially, it is
  **independent of the local catalog**, so it is not circular with ETAS the way every
  catalog-derived criticality measure (n(t), b, entropy) necessarily is. That independence is
  what makes it worth more than all of K-018..K-026 combined if it works.
- Claim: define per region a susceptibility S(region, t) = fitted rate response per unit peak
  dynamic stress, estimated over a rolling multi-year window from all teleseismic pings. Then:
  (i) S is spatially structured and reproducible; (ii) S varies in time by more than estimation
  noise; (iii) high S precedes local large events.
- Test: global ComCat M>=5.5 1980-now as the ping catalog (one download); targets = the 13
  regional catalogs on disk plus SoCal at M>=2.5. For each ping, the response is the rate change
  in a 0-72 h window after the predicted surface-wave arrival, relative to the ETAS-expected
  rate, at distances beyond 2 rupture lengths. Pool across pings weighted by PGV to fit S.
  Statistic: (i) inter-region variance of S vs its estimation error; (ii) an ANOVA-style test of
  temporal variation in S within region; (iii) ROC AUC for "M>=6.5 in this region in the next
  1-2 yr" from S(t).
- Null: circular shift of the ping catalog against each target catalog (the established
  machinery), which is exact here because pings and targets are causally independent under the
  null; plus ETAS-sim targets.
- Expected: (i) will work — geothermal/volcanic/extensional regions will show high S; this is
  effectively K-034 generalized and doubles as its global validation. (ii) is the real question
  and I would put it near even odds. (iii) is the prize and I would put it low, maybe 20% — but
  a 20% shot at a *non-circular, actively-probed, globally-available criticality gauge* is the
  best-expected-value item in this seed.
- Why dismissed too quickly: "triggering thresholds are too high; most regions never respond."
  Then S = 0 with a tight bound almost everywhere, and the map of where S > 0 is still a new and
  useful object — it is a map of where the crust is close enough to failure that a passing
  surface wave matters, which is a hazard-relevant statement independent of any forecast. And the
  null costs one download and reuses the K-034 machinery.

**K-044 — FRAME-BREAK: admittance spectroscopy. Stop testing frequencies one at a time; measure
the whole transfer function from stress to seismicity — and rate-state predicts its SHAPE.**
- The presupposition attacked: that "periodicity" is a list of candidate periods to test. EXP-F
  tested a 60-period comb and produced an ambiguous null with a failed positive control. The
  frame is wrong: a driven dissipative system does not have "periods", it has a **transfer
  function** — a complex-valued admittance H(f) mapping stress forcing to rate response, with an
  amplitude and a phase at every frequency.
- Why this is not the periodicity corpse: EXP-F asked "is seismicity periodic?" (a property of
  the output alone, with no input). This asks "what is the linear response of seismicity to a
  *measured input*?" — a coherence between two observed series, with the input's own spectrum
  divided out. Those are different statistics with different nulls, and the second one is far
  more powerful because it uses the known forcing. Also decisive: **rate-state predicts the shape
  of H(f)** — a high-pass response with a corner at the Dieterich aftershock timescale
  t_a = Aσ / τ̇, giving a corner period of months to years. Measuring that corner *measures Aσ*
  independently of K-036, and two independent estimates of the same physical constant agreeing
  would be far stronger evidence than either alone.
- Test: form the total stressing-rate series (K-036) and the seismicity-rate series per region;
  compute magnitude-squared coherence and the complex admittance across 1/hour to 1/decade
  (multitaper, with the ETAS-expected rate removed first so the response is measured against the
  clustering baseline). Statistic: coherence with confidence bounds; the fitted corner frequency
  and implied Aσ; comparison of the measured |H(f)| against the rate-state prediction by χ².
- Null: ETAS-sim seismicity driven by the *real* stress series but with zero coupling — gives the
  exact coherence distribution under "no response", including the coherence that arises from both
  series sharing seasonality.
- Expected: coherence low but possibly significant at the lowest frequencies; the corner is the
  target. A null with tight bounds is genuinely valuable here too — it would bound Aσ from below
  across California, which is a crustal property nobody has mapped.
- Why dismissed too quickly: "this is the periodicity search again, and that comb had weak
  power." It is the inverse problem, not the same one: the input series is measured rather than
  assumed, the statistic is coherence rather than Rayleigh, and the hypothesis has a predicted
  functional form to be tested against rather than a list of frequencies to be scanned. And the
  EXP-F failure mode — slow rate fluctuations masquerading as cycles — is precisely what dividing
  by the measured input spectrum removes.

**K-045 — FRAME-BREAK: use the regime where the signal is enormous. Induced seismicity is a
known-stress, known-timing, high-SNR experiment; fit the response function there and TRANSFER it
to natural loads.**
- The presupposition attacked: that we must learn the crust's response to stress from natural
  loads, where the signal is at the edge of detectability. We do not. Humans have been running
  large, well-documented stress experiments — Oklahoma wastewater injection, reservoir
  impoundment, geothermal operations at Coso and the Salton Sea (both already in our
  `tsi_map.GEOTHERMAL` list) — with *published input time series*. Signal-to-noise there is orders
  of magnitude better than anything the tides offer.
- The lens (this is the move I most want run): **calibrate in the strong-signal regime, transfer
  to the weak-signal regime.** Fit the conditional-intensity response function — including the
  state interactions Jim is asking about — where the stress input is large and precisely known.
  Then take the *functional form and the interaction structure* (not the coefficients) as a
  prior for the natural-load fits. This is transfer learning applied across signal regimes rather
  than across geography, and it is a direct riff on B-1: EXP-M taught us that the transferable
  object is the universal shape with locally-calibrated amplitude. Same principle, new axis.
- Test: data — Oklahoma Corporation Commission public injection volumes (well-level, monthly,
  2011-now) plus ComCat Oklahoma catalog (one query); Coso/Salton Sea production data (CalGEM/
  DOGGR public). Fit Cox-ETAS with injection-derived pore-pressure covariate and the same frozen
  interaction grid as K-033. Statistic: bits/event, the fitted response shape, and the fitted
  Aσ-equivalent. Then the transfer test: does imposing the Oklahoma-fitted response *shape* on
  the SoCal natural-load model improve out-of-sample bits/event over a free-form fit? (This is
  precisely the EXP-M sign-test design reused, and it can fail cleanly.)
- Null: ETAS-sim; and the honest control that Oklahoma's response may be entirely
  pore-pressure-diffusion-driven and share no mechanism with elastic tidal loading — in which
  case transfer fails and we have learned that the two regimes are physically distinct, which is
  itself a real finding about mechanism.
- Expected: strong detection in Oklahoma (near-certain — this is a positive control as solid as
  Landers, and should be treated as a second gate on the machinery); transfer to natural loads is
  maybe 30%. But the Oklahoma fit alone gives us the first *well-constrained* measurement of the
  interaction structure Jim is asking about, at SNR the natural experiment will never provide.
- Why dismissed too quickly: "induced seismicity is a different phenomenon, off-topic for a
  tidal-triggering program." That objection assumes the answer to the question the test asks.
  And strategically: this program's scarcest resource is signal, and there is a regime nearby
  where signal is abundant. Refusing to calibrate there because it is "not our topic" is how a
  research program stays underpowered on purpose.

**Ordering I would recommend to the supervisor for this seed:** K-034 and K-035 first (they are
cheap, they are gates, and their outcomes re-price everything else); then K-033 as the engine;
then K-036/K-039/K-040 as the substantive conditional tests; then K-043 and K-045 as the two
highest-upside frame-breaks. K-041 and K-042 can ride along on K-036's covariate build for
almost no marginal cost.

---
*End Kepler round 1 (K-001..K-045). All PROPOSED. Nothing above is claimed as true; nothing
above has been tested. Popper to adjudicate; supervisor to run frozen tests. The charter
amendments are proposals to my own definition, offered because I was told the file is a floor —
and because points 2, 3 and 7 are the ones I would want enforced against me.*

## VERDICTS (Popper)

Round 1 (K-001..K-032 + charter amendments) adjudicated 2026-08-09 — see
**"VERDICTS (Popper) — Round 1"** at the end of this file.

## BASELINES (validated additions)

(none yet beyond B-1..B-5)

---

# VERDICTS (Popper) — Round 1, K-001..K-032 + charter amendments
*Adjudicated 2026-08-09. Kepler's entries above are untouched. Nothing below is a result; it is
a set of frozen specs, kill-reasons, and rulings. The supervisor runs; I rule on outcomes.*

Kepler: this is a strong round. The standing methodological demand (ETAS-sim nulls) is correct
and I adopt it with two riders below. Roughly two-thirds of this slate is runnable this week on
files already on disk, which is the right shape for a round. My objections are almost entirely
about *estimators*, not about ideas.

---

## 0. SHARED STANDARDS FOR THIS ROUND (read before any individual verdict)

These apply to every entry. They exist so the same paragraph is not written thirty-two times.

### S-1. The ETAS-sim null is right, and it has exactly two blind spots. Name them or lose.

Kepler's demand is adopted as the default. Two riders are mandatory:

**(a) Circularity.** When the hypothesis *is* that ETAS is misspecified in the very dimension
being simulated, simulating from ETAS tests ETAS against itself and can only ever return "no
excess". This bites K-005, K-006, K-011, K-025, K-029 and the payoff arm of K-022. What breaks
the circularity, in order of preference:
  1. **Two-generator discrimination.** Simulate from ETAS *and* from a second generator that
     embodies the alternative (tapered-GR magnitudes; depth-stratified productivity;
     magnitude-correlated triggering; a mu with an OU-driven latent state). The statistic must
     separate the two generators at the real sample size. If it cannot, the entry is
     underpowered and must not be run as confirmatory.
  2. **Estimator-bias calibration** (not circular, and the correct use). Where the claim is
     about the *size of an estimator's bias* — K-005 exactly — simulating from a generator that
     is self-similar by construction is precisely right, because the sim's known truth is the
     reference. Say which of the two uses you are making.
  3. **Cross-pipeline replication** (K-028) where no second generator is available.

**(b) Simulated catalogs are perfectly observed.** An ETAS-sim contains no Mc drift, no network
evolution, no magnitude-scale revision, no border truncation, no coda blindness. It therefore
**cannot null an observer artifact** — and this program's corpses are disproportionately
observer artifacts (EXP-B's small-n amplitude bias; §14b's Anza phase, +159° in SCSN and +25°
in QTM, same rock, same decade; the (32.8,−115.1) border-completeness cell in EXP-K; EXP-M's
Caribbean n=235). **Rule: whenever completeness, Mc, b, station density, magnitude scale, or
box edges lie on the causal path of the statistic, the ETAS-sim must be passed through a
detection function Mc(x,t) derived from the real network history, or the entry must carry a
recorded exemption naming the untested exposure.** This bites K-007, K-013, K-017, K-018,
K-019R, K-021R, K-026, K-002.

**(c) The null inherits the selection rule.** Target windows, matched controls, magnitude
strata, and window overlap must be constructed on simulated catalogs by the identical code
path. A null that skips the selection rule is not a null.

### S-2. Overlapping-window autocorrelation. This program has already died of this once.

§14b is the record: the "episodic phase-migration finding" at Anza was framed as a discovery
after EXP-C2 and turned out to be **overlapping-window autocorrelation dressed as coherence**;
two designed confirmations (E1 cross-catalog, E2 walk-forward) nulled it. Every rolling-window
entry in this round repeats that geometry: K-018 (1-yr windows stepped 1 month = 92% overlap),
K-019 (3-yr rolling AR1, trend measured over 3 yr — the trend statistic is nearly a
deterministic function of the window), K-020 (1-yr rolling xi), K-021 (3-yr rolling Var(b)),
K-026 (monthly, less exposed).

**Rule: any rolling statistic must report (i) the window overlap fraction, (ii) an effective
degrees-of-freedom estimate, and (iii) a non-overlapping-window arm as the confirmatory one,
with the overlapping arm labelled descriptive.** The sim null absorbs autocorrelation *only*
if S-1(c) holds exactly.

### S-3. Standalone AUC against "matched controls" is not skill. The gate is incremental bits.

This is the single most important thing I can tell you about the emergence lens. Every order
parameter in K-018/019/020/021/022/026 — n(t), AR1, variance, xi, S_max, entropy — is a
monotone function of the local event rate. The event rate rises before large events **under
plain ETAS**, because large events follow busy periods. So every one of these will produce
AUC 0.55–0.65 for "large event within 6 months" *by proxy for λ*, and Kepler's own expected
values (0.55–0.65, five times) are exactly the range a rate proxy delivers for free.

**Frozen gate for the whole order-parameter family (G2 below): the pre-registered statistic is
the incremental log-score, in bits/event, of [frozen ETAS λ + order parameter] over [frozen
ETAS λ alone], walk-forward, on a temporal holdout. AUC may be reported; it may not be the
success rule.** An order parameter that cannot add bits over λ is descriptive, not a predictor,
and will be recorded as such.

### S-4. n-bias is this program's most reliable killer. Assume it until measured.

EXP-B's finding is the template: every apparent feature-vs-amplitude correlation was the
small-n amplitude bias in disguise (a_b vs n_train rho = −0.49, p = 0.001), and the
bias-robust label correlated with nothing. The same shape is present, unaddressed or
under-addressed, in K-019 (variance), K-020 (xi grows with n), K-021 (Var of b̂ ~ b²/n),
K-023 (corner magnitude is n-biased upward), K-026 (entropy is n-biased), K-030 (fragility ~
1/n and cross-pipeline disagreement ~ 1/n, so their correlation is near-guaranteed).

**Rule: where an estimator's bias depends on sample size, the fix is subsampling to a fixed
common n, not an analytic correction.** Analytic corrections (b²/n, entropy bias expansions)
are asymptotic and under-correct at the n this program actually has. Fixed-n subsampling with
replication over draws costs a loop and removes the failure mode by construction.

### S-5. Baselines must be harsh, and our own headline numbers need one correction now.

B-1 is currently quoted as "+0.66..+1.75 bits/event on six holdouts". The +1.75 is Caribbean,
`n_scored = 235`, carrying `underpowered: true` in results_exp_m.json. **The honest quotation
is "+0.66 to +0.84 bits/event on the four adequately-powered holdouts (AK 0.84, MX 0.66,
PH 0.78, GRC 0.79; IRN 0.69), with Caribbean +1.75 flagged underpowered."** Use that form from
now on. Likewise B-2's "skill rises with magnitude" rests on `n = 290` events in the M≥4 band —
real, but a block bootstrap over sequences, not events, is required before it is quoted with a
CI (M≥4 events in SoCal are not independent draws; they arrive in a handful of sequences).

### S-6. The grouping. Four families, four shared failure modes, named once.

- **G1 — the α / branching family: K-005, K-006, K-018** (and the slope in K-001, the stratified
  gap in K-014). All are statements about the productivity exponent or the branching ratio,
  all estimated from magnitude-truncated samples, all M0- and n-dependent. EXP-H already shows
  the instability in our own record: nominal n = 1.161 (supercritical) with a 7-day n_eff = 0.60.
  EXP-M shows α = 0.537 (SoCal, M0 = 2.5) vs 0.730 (global pool, M0 = 4.5) — a gap that may be
  entirely a magnitude-floor artifact.
  **Dependency ruling: K-005 runs FIRST. Its measured bias curve (parameter drift vs M0 on
  self-similar simulated catalogs) is a required input to K-006 and K-018.** Without it, a
  "bend in α" and "excursions in n(t)" are the same artifact wearing two costumes.
- **G2 — the rolling-window order-parameter family: K-018, K-019, K-020, K-021, K-026.**
  Shared failure: S-2 (overlap) + S-3 (rate proxy) + S-4 (n-bias). **Ruling: these run as ONE
  job producing ONE joint model, not five papers.** Each variable enters as a candidate
  covariate on the frozen ETAS λ; the deliverable is a single incremental-bits table with the
  covariates' partial contributions and a shared count-matched sim null. Five separate
  underpowered AUC studies is how this family has failed everywhere else.
- **G3 — the observer family: K-015, K-028, K-030, K-031** (with completeness exposure in
  K-002, K-013, K-017, K-021R). Shared insight: the catalog is an instrument reading, and
  S-1(b) means the round's default null cannot police this family. G3 polices it. **Ruling:
  K-031 and K-028 run as one job — they load the same three catalogs and answer adjacent
  questions.**
- **G4 — spatial dependency: K-002 gates the spatial slate.** K-003, K-013, K-022, K-025R and
  K-030 all make spatial claims. **No spatial entry may quote bits until K-002 exists as the
  floor.** Kepler is exactly right that this program has validated a temporal object and keeps
  reaching for spatial physics with no spatial baseline. Fix that before spending on K-003.
- **G5 — the budget family: K-011, K-027, K-029, K-032.** Shared exposure: all four are
  measurements of *our model's* properties being reported as properties of *the Earth*. Each
  carries a mandatory scope line to that effect.

---

## 1. VERDICTS — Q1 (when/where the next M≥6.5+)

### K-001 — TESTABLE-NOW.
Good entry. Pure re-scoring of runs we already own; the highest bits-per-compute-hour item in Q1.

- **Statistic.** Bits/event of the frozen GLOBAL-pool + local-mu ETAS over a local-oracle
  Poisson, computed separately in magnitude strata [4.5,5), [5,5.5), [5.5,6), [6,6.5), [6.5,7),
  [7,∞); pooled across the 13 regions with per-region weights; slope of bits vs stratum by
  weighted least squares.
- **Null.** Slope ≤ 0.
- **Success.** Slope > 0 with a **sequence-block bootstrap** 95% CI excluding 0, AND bits/event
  at [6.5,∞) > 0 with its own CI excluding 0. **Failure:** CI covering 0, or bits ≤ 0 at
  [6.5,∞).
- **The crux, and the reason Kepler's spec is not yet frozen-ready.** M≥6.5 events are not
  independent draws — 517 of them worldwide arrive in a much smaller number of sequences. An
  event-level bootstrap will produce a CI two to four times too narrow. **Mandate: block
  bootstrap with blocks = Zaliapin–Ben-Zion clusters (the program's existing parameters: b=1,
  df=1.6, log10 η0=−5), not events.** Report the effective number of independent blocks per
  stratum alongside the raw n. If a stratum has fewer than 30 blocks, it is reported as
  suggestive-only.
- **Positive control.** The [4.5,5) stratum must reproduce EXP-M's holdout numbers to within
  bootstrap noise (+0.66..+0.84 on the four powered holdouts). If it does not, the harness is
  broken and nothing else in the table means anything.
- **Negative control.** Score with a *history transplant*: compute λ using the event history
  from a random other region (same length, same rate normalisation). Bits must be ≤ 0 in every
  stratum. This is the check that the skill is conditional information and not a rate-fitting
  artifact.
- **Leakage risk (named).** None temporal — walk-forward is already in place. The real leakage
  is **magnitude-scale drift**: ComCat rows at M≥6.5 are mostly `mww`, at M 4.5–5.5 mostly `mb`
  or `ml`. A slope in bits vs magnitude could be a slope in magnitude-scale quality. **Mandate:
  report the magType composition per stratum, and run a GCMT-Mw-only robustness arm for the
  three strata above M 5.5.**
- **Artifact class.** Magnitude-scale drift (primary); block-dependence in the CI (secondary).
- **Data.** `data/comcat_world/*.csv`, frozen global pool from results_exp_m.json. No downloads.
- **Scope of a win.** "Generic-ETAS-with-local-μ retains/gains skill on M≥6.5 targets, 13 ComCat
  boxes, 1995–2026, temporal-only, vs local-oracle Poisson." It would say nothing about *where*
  and nothing about unforeseen large events in quiet regions. Report `triggered_fraction` per
  stratum in the same table so that limit is visible on the face of the result.
- Kepler's pre-emptive answer to "it's just aftershock forecasting" is accepted: the stratified
  triggered-fraction column settles it empirically, and a conditional statement is deployable.

### K-002 — TESTABLE-NOW. And it is a prerequisite, not an option.
You are right and I am recording that the program's shape was wrong here: we validated a
temporal object and have been proposing spatial physics with no spatial floor. **G4 stands:
nothing spatial quotes bits until this exists.**

- **Statistic.** CSEP-style space-time-magnitude Poisson log-likelihood on 0.5° × 30 d bins,
  M≥6.5 targets, 2010→2026 walk-forward, reported as bits/event above (a) uniform-in-footprint
  and (b) **time-independent smoothed seismicity alone** — (b) is the number that matters, since
  it isolates what the temporal layer adds to *where*.
- **Null / success.** Success = bits over (a) > 0 with block-bootstrap CI excluding 0 (this is
  the sanity floor and should be easy), AND the *reported* decomposition (b) with its CI, which
  is a measurement, not a pass/fail. **Pre-register that (b) may be ≈ 0; that is a legitimate
  and useful outcome and must not be re-described afterwards.**
- **Positive control.** Smoothed seismicity must beat uniform by ≥ 2 bits/event. It always does;
  if it does not here, the pipeline is broken.
- **Negative control.** Kepler's rotated-μ field (adopt), plus a **k-fold spatial holdout**:
  build μ(x) with one 30° longitude band withheld and score only in that band. This is the
  control that catches the real leakage.
- **Leakage risk (named).** μ(x) is built from 1995–2009 events; the declustering used to build
  it must use train-window events only. If the declustering algorithm sees the full catalog,
  the spatial density is contaminated with post-2010 aftershock zones — which is where the
  M≥6.5 targets are. **This is the single way this test can fool us, and it is easy to do by
  accident.** Freeze the declustering input window explicitly in code.
- **Artifact class.** Border/completeness gradients (the 13 boxes have edges; events near an
  edge have truncated smoothing kernels and truncated aftershock zones). **Mandate the EXP-J
  border check pattern: report the interior-vs-edge skill split** (§15 records that check
  rescuing EXP-J's silent list; use it again here).
- **Data.** On disk. `q` and `d` for the spatial aftershock kernel are fit on the 7 EXP-M train
  regions only and frozen.
- **Cost.** This is the most expensive item I am approving. Approved anyway, because five other
  entries are meaningless without it.

### K-003 — NEEDS-DATA. Spec below is frozen and waits on three downloads.
The distinction Kepler draws from the 1990s seismic-gap corpse is legitimate and I accept it:
a moment-budget covariate scored against a smoothed-seismicity floor with a loading-matched
null on a blind grid is not the gap hypothesis. But it cannot run without a global strain model,
and it must not run before K-002.

- **Missing data, with the cheapest acquisition path.**
  1. **GSRM v2.1** (Kreemer et al. 2014), global 0.25° strain-rate grid — GEM/Zenodo public
     download, single file, tens of MB. Log the SHA-256 in `download_log.md` at retrieval, per
     the program's discipline.
  2. **ISC-GEM** global instrumental catalogue (1904–2018, M≥5.5) — isc.ac.uk, single CSV.
     Required so that Ṁ_seis is not biased low the way EXP-J's χ was (§15: "χ absolute values
     biased low — geography is the signal, not the level").
  3. **Slab2** (USGS ScienceBase) for seismogenic thickness in subduction zones.
  All three are free, public, one-shot. Estimated acquisition: under an hour.
- **Statistic.** Nested Poisson-regression Δlog-likelihood for 1995–2026 M≥6.5 counts per 0.5°
  cell: [log smoothed seismicity] vs [+ SILENT-class indicator]; plus the silent-vs-coupled
  rate ratio.
- **Null (the load-bearing one).** Label permutation of the silent class **within
  loading-matched strata**. Kepler names this and he is right that it is essential; I add that
  the strata must also match on **log smoothed seismicity**, because the silent class is defined
  by *low release*, and low release is low count, which is low smoothed seismicity — permuting
  within loading strata alone leaves the covariate that the model is already using.
- **Success.** Δ log-likelihood significant against the doubly-matched permutation null at
  p < 0.01, AND rate ratio CI excluding 1. **Failure:** either fails.
- **Positive control.** The silent-class construction, applied to SoCal with data already on
  disk, must reproduce EXP-J/K's silent list (Jaccard ≥ 0.9 against the EXP-K variable-H list).
  If the global pipeline cannot re-derive our own validated result on our own region, stop.
- **Negative control.** Assign the silent label to a random loading-and-density-matched cell set;
  Δ LL must be ≈ 0.
- **Artifact class.** Border/completeness (ISC-GEM completeness varies enormously by era and
  region pre-1960 — restrict the training era to a completeness-verified window and report the
  Mc(t) used); and the **n<20 detection-limited trap that already infects B-4** (EXP-K: 158 of
  200 unexplained-silent cells have n_train < 20; 39 have n_train = 0). Globally this will be
  worse. **Mandate: report the analysis twice — all silent cells, and measured-low-χ cells only
  (n_train ≥ 20). If the signal lives only in the detection-limited set, it is K-031's finding,
  not K-003's.**
- **Honest prior, recorded:** I agree with Kepler that the matched null will eat most of the raw
  signal. That is the point of running it.

### K-004 — TESTABLE-NOW, with the success rule rewritten. Kepler's rule has no power.
Riffing off a validated failure is exactly the right move and this is the cheapest remaining
test of whether *any* stratification helps. But the rule as written ("depth pool beats GLOBAL in
≥5 of 6 holdouts") is EXP-M's rule, and EXP-M's rule is nearly powerless: under the null a sign
test on 6 trials gives p = 0.109 at 5/6 and p = 0.016 only at 6/6. We failed EXP-M 2/6 and
learned a great deal from the *shape* of the failure, not the sign test. Do not repeat the
instrument.

- **Rewritten statistic (primary).** Pooled bits/event difference (DEPTH pool − GLOBAL pool)
  across the 6 holdouts, with a **region-block bootstrap** (resample regions with replacement,
  and within region resample Zaliapin clusters). Report per-region contributions.
- **Success.** Pooled Δ bits > 0 with 95% CI excluding 0. **Failure:** CI covering 0.
  Sign test reported as secondary/descriptive only.
- **Specification correction (important).** Do not "split the catalog at 70 km" — a deep event
  can trigger a shallow one, and splitting the catalog destroys the triggering graph.
  **Stratify the source term only:** λ(t) = μ + Σ_i K_{z(i)} · 10^{α_{z(i)}(M_i−M0)} (t−t_i+c)^{−p},
  with (K, α) taking a shallow or deep value according to the *source* event's depth, μ, c, p
  shared. Two extra free parameters, fit on the 7 train regions, frozen, scored on the 6
  holdouts. That is a cleaner test of the same idea and it is what "deep events are
  aftershock-poor" actually asserts.
- **Positive control.** With the depth split disabled (K_shallow = K_deep), the harness must
  reproduce EXP-M's GLOBAL-pool holdout numbers exactly. Bit-for-bit; this is a refactor check.
- **Negative control.** Split on a **random** binary label with the same marginal frequency as
  the depth split, refit, rescore. Pooled Δ bits must be ≈ 0. This is the control that catches
  "two extra parameters always help".
- **Leakage risk.** None temporal. The real exposure is **depth quality**: ComCat depths are
  poorly constrained for offshore events and many are fixed at 10 km or 33 km by convention.
  **Mandate: report the fraction of events at fixed/default depths per region and exclude them
  from the *source* stratification (assign them to a third "depth-unknown" pool sharing the
  global K, α).** Without this, the depth label is partly a label for "offshore / poorly
  located", which correlates with everything.
- **Artifact class.** n-dependent estimator bias (deep pools are small in most regions — report
  n per pool; Alaska 10% deep, Greece 4%, so those pools are thin) plus depth-assignment
  convention.
- **Diagnostic that survives a fail.** K_deep/K_shallow with bootstrap CI. Kepler is right that
  this is worth having regardless.
- **Data.** All on disk.

---

## 2. VERDICTS — Q2 (self-similarity, and whether the break predicts)

### K-005 — TESTABLE-NOW. Run this before K-006 and K-018 (G1 dependency).
This is the best-designed entry in Q2 and its ETAS-sim control is the *non-circular* use of
simulation (S-1(a)(2)): the simulated catalogs are self-similar by construction, so they
measure estimator bias against known truth. That is exactly right and I want it said in the
results file so nobody later mistakes it for a circular test.

- **Statistic.** For each parameter θ ∈ {α, c, p, n, K′} (K rescaled to a common reference),
  the slope dθ/dM0 from the ladder M0 ∈ {2.5, 3.0, 3.5, 4.0, 4.5} (SoCal, SCSN train < 2010)
  and {4.5, 5.0, 5.5} (per-region world). Bootstrap CI over 200 catalog bootstraps.
- **Null.** The slope obtained by the identical ladder applied to ETAS-simulated catalogs
  generated from the M0 = 2.5 fit — i.e. the pure estimator-bias slope.
- **Success / failure.** A parameter shows a **real** break if |slope_observed − slope_sim|
  exceeds 2σ of the sim slope distribution. Otherwise the drift is bias and is recorded as a
  **correction table**, which is the deliverable either way.
- **Freeze the rescaling formula in the protocol text before running.** EXP-M's record (§17)
  is that the brief carried a sign error and the worker caught it: the correct form is
  K′ = K · 10^{+α(M0new − M0old)}. Write it once, in the protocol, with the derivation.
- **The consequential arm (keep it — it is what makes this more than a bias study).** Predict
  the M≥4.5 test-window rate using parameters fit at M0 = 2.5 vs at M0 = 4.0; report the
  bits/event difference. That number is the operational cost of the drift.
- **Positive control.** Inject a *known* α-drift into a simulated catalog (generate with
  α varying by magnitude band) and confirm the ladder recovers its size.
- **Negative control.** The self-similar sim ladder itself (this is both null and negative
  control here — say so).
- **Leakage risk.** The M0 ladder changes the *sample*, and the SCSN sample's Mc changes by era
  (results_exp_h.json: per-decade Mc 1.7 / 1.7 / 1.6 / 1.2). **Mandate: hold the time window
  fixed across every rung of the ladder** — otherwise the M0 = 2.5 rung is effectively a
  different-era catalog from the M0 = 4.5 rung and the drift is a completeness drift.
- **Artifact class.** n-dependent estimator bias (the whole subject), magnitude-scale drift
  across eras (mitigated by the fixed window).
- **Decision value.** If the drift is bias, this program stops quoting transferred parameters
  across catalogs with different M0 — which it currently does, in EXP-M's headline table — and
  the SoCal α = 0.537 vs global α = 0.730 gap is reinterpreted as an artifact of the magnitude
  floor rather than as tectonics. That retroactively re-reads §17. If the drift exceeds the
  simulated bias, we have measured a genuine break in crustal self-similarity and large-event
  probabilities get a calibration correction. Both outcomes are load-bearing. Kepler's "boring
  because it's known bias" pre-emption is correct and I endorse it.

### K-006 — REFRAMED → **K-006R** (credit K-005/K-006, Kepler). The claim as posed is artifact-guaranteed.
The physics question is good. The measurement is not: "count aftershocks M≥4.5 within 100 d and
2 rupture lengths of every mainshock M≥5, regress log N on M, fit a bilinear with free
breakpoint" **will produce a downward bend at high M whether or not the physics bends**, for
three independent reasons:

1. **Short-term aftershock incompleteness scales with mainshock size.** After an M8 the coda
   swallows far more M4.5s, for far longer, than after an M5.5. This alone bends the
   productivity curve down at exactly the magnitudes where the bend is predicted. This is
   K-015's entire subject.
2. **Border truncation scales with mainshock size.** Big events have big aftershock zones; a
   2-rupture-length radius around an M8 in a 13° box routinely leaves the box. Counted
   aftershocks are then truncated preferentially for large mainshocks. That is the
   border/completeness corpse class (§16, cell (32.8,−115.1)) at a different scale.
3. **End-of-catalog censoring** hits the largest events preferentially because they are rare
   and the largest ones in a 31-year window are disproportionately recent.
An ETAS-sim null does not save this: per S-1(b), the sim has no coda blindness and no box edges,
so the sim's curve will be straight and the real curve will bend, and the test will "pass"
artifactually.

**K-006R, the surviving cousin:**
- **Dependency.** Runs only after K-015 delivers a fitted detection function Mc(t | M_main).
- **Statistic.** Productivity N̂(M) = the **incompleteness-corrected** count of M≥4.5 aftershocks
  over days **3–100** (day 0–3 excluded outright, as a second defence), inside a radius that is
  fully contained in the region box (mainshocks whose 2-rupture-length disc crosses a box edge
  are **excluded**, and the excluded fraction is reported per magnitude bin), with an
  end-of-catalog cutoff requiring 100 d of post-event catalog. Regress log10 N̂ on M; compare
  single-line vs free-breakpoint bilinear by ΔBIC.
- **Null.** ΔBIC < 6, or breakpoint profile-CI covering the catalog Mmax.
- **Success.** ΔBIC ≥ 6 AND breakpoint CI strictly interior to the observed magnitude range
  AND the bend survives the **GCMT-Mw-only arm** (this defuses the magnitude-scale-saturation
  objection, which is otherwise fatal and which Kepler correctly anticipated).
- **Positive control.** Simulate sequences with a *known* α bend at M 7.2; the pipeline must
  recover the breakpoint within its CI at the real sample size. If it cannot, the entry is
  underpowered and must not be run as confirmatory.
- **Negative control (the decisive one).** Simulate sequences with **no** bend but **with** the
  fitted detection function and box truncation applied. The pipeline must return ΔBIC < 6. If
  it does not, the correction is inadequate and K-006R is not yet runnable.
- **Artifact class.** Short-term incompleteness (primary), border truncation, magnitude-scale
  saturation.
- **What a win would mean.** Real, and worth 20–50% of the 7-day M≥6 aftershock probability
  after an M8 — a genuine operational correction at the moment forecasts are actually consulted.

### K-007 — TESTABLE-NOW, with a power control added.
The design idea is genuinely good and I want to say so: putting the observer artifact (low-M
completeness) and the physical signal (high-M curvature) in different parts of the same
distribution is the right structural move, and it is what distinguishes this from the dead
b-mean literature. Kepler's three named fixes (shape not level, ETAS-sim null, sample-size
matching) are the three things the corpses lacked.

- **Statistic.** Per M≥6.0 target: trailing 2-yr, 100-km sample above a conservative fixed Mc;
  log-likelihood ratio of tapered-GR (free corner) over pure GR; and the fitted corner
  magnitude. Compare pre-target distribution against matched controls (same cell, no target
  within ±1 yr), **matched on sample size** (S-4: subsample both to a fixed common n).
- **Null.** ETAS-sim with pure GR magnitudes, run **through a simulated Mc(x,t)** per S-1(b),
  with the identical target-selection and matching code (S-1(c)). Under this null there is no
  taper ever.
- **Success.** Pre-target LR statistic exceeds the sim-null 99th percentile pooled, AND the
  sign of the corner-magnitude shift is consistent (pre-registered as two-sided: down = the
  system suppressing large events, up = enabling — either is a result, but the direction is
  frozen before unblinding and reported as predicted).
  Secondary: incremental bits over frozen ETAS λ per S-3. **Failure:** pooled statistic inside
  the null envelope.
- **POSITIVE CONTROL — missing from Kepler's spec and mandatory.** Inject a taper of known
  corner magnitude into simulated pre-target windows at the *real* sample sizes and confirm the
  LR statistic detects it at ≥ 80% power. Without this, a null result is uninterpretable — and
  a null result here is the likely outcome, so its interpretability is the whole value. (This
  is also the K-032 item-6 discipline applied prospectively rather than retroactively.)
- **Negative control.** ETAS-sim with pure GR (the null, doubling as the negative control).
- **Leakage risk (named, and Kepler missed it).** If the target M≥6.0 event is itself a
  triggered event, its trailing 2-yr window contains the *parent* sequence, whose magnitude
  distribution differs from background. **Mandate: classify targets by Zaliapin declustering and
  report background-parent targets separately from triggered targets.** Mixing them makes
  "pre-event taper" partly "aftershock-sequence magnitude distribution".
- **Artifact class.** Completeness gradients at the low-M end leaking into the corner estimate
  via the fit (mitigated by the conservative fixed Mc); n-bias in the corner estimator (S-4).
- **Data.** On disk: SCSN M≥3.0 1981–2018 and the 13 world boxes M≥5.0.

### K-008 — TESTABLE-NOW. And it attacks the assumption that matters most.
Along with K-022, this is one of only two entries in the round that attacks ETAS's i.i.d.-GR
magnitude assumption — the assumption that makes prediction impossible by construction. A clean
zero is genuinely valuable and Kepler is right to say so.

- **Statistic.** CRPS of a gradient-boosted quantile regression for cluster Mmax, trained on
  clusters starting < 2010, scored on ≥ 2010, vs the ETAS+GR predictive distribution for Mmax.
- **THE LEAKAGE RISK, and it is fatal if unaddressed (Kepler's spec does not state it).** The
  label is "eventual cluster Mmax". If Mmax occurs within the first 24 h, the features (which
  include the first-24-h event count and local b) contain the label. **Mandate: the label is
  max magnitude occurring AFTER the 24-h feature window**, and clusters whose overall Mmax fell
  inside the window are either excluded or handled as a separate, clearly-labelled arm. Without
  this the experiment cannot fail.
- **Null.** (a) The ETAS+GR baseline predictive distribution; (b) label permutation within
  first-24-h-count strata — Kepler has this and it is right, because count alone predicts Mmax
  through GR sampling and that is not skill.
- **Success.** CRPS improvement ≥ 3% over the ETAS+GR baseline with a cluster-block bootstrap CI
  excluding 0, AND surviving the count-stratified permutation. **Failure:** either fails.
- **Mandatory ablation.** Report three models: onset-dynamics features only; spatial features
  only (μ of host cell, distance to CFM trace / plate boundary); both. Kepler's own honest prior
  is that the win comes from the spatial features — i.e. "where" leaking in, which is a
  *different and less interesting claim* than "a sequence knows how big it will get". The
  ablation is what keeps the headline honest.
- **Positive control.** Train and score on simulated clusters where Mmax is by construction a
  function of the onset count; the pipeline must recover a large CRPS gain.
- **Negative control.** Train and score on ETAS-sim clusters (i.i.d. GR by construction); CRPS
  gain must be ≈ 0. If the pipeline shows a gain on ETAS-sim data, it is fitting the harness.
- **Artifact class.** Label leakage (above); cluster-definition sensitivity — report the result
  at two η0 thresholds so it is not an artifact of the declustering parameter.
- **Data.** On disk.

---

## 3. VERDICTS — Q3 (is there an internal weather system?)

### K-009 — TESTABLE-NOW. **Run it first.** And it has one flaw that would have inverted its answer.
Kepler's framing is right and it is the sharpest question in the round: the residual is the
innovation; white innovations mean no state to estimate and no amount of compute helps; red,
spatially-coherent innovations mean there is a latent field and it is estimable. Zero new data.
This is the go/no-go for the entire assimilation thread and for a large part of Jim's Q3.

**The flaw, and it is serious.** The spec computes expected counts by taking the frozen
*temporal* EXP-H ETAS and "rescaling to spatial cells by a smoothed-seismicity kernel", then
compares residual structure against ETAS-simulated catalogs. But the simulated catalogs would be
generated *from that same kernel*, so they have zero spatial-kernel misfit, while the real data
has whatever misfit the ad-hoc kernel carries. The excess correlation would then measure **"my
spatial kernel is wrong"**, not "there is weather" — and the test would return a strong,
exciting, false positive. This is the anti-conservative twin of the circularity problem in
S-1(a): here the sim is *too clean*, not too similar.

- **Required fix (three parts, all cheap).**
  1. Generate the null from a **spatio-temporal** ETAS whose background field and aftershock
     kernel were fit to the real catalog (the K-002 machinery), not from a temporal model
     smeared by a kernel. The null must have the same *degrees of freedom* as the model whose
     residuals we are testing.
  2. **Kernel-swap control.** Recompute residuals under two deliberately different background
     fields (k = 4 and k = 8 adaptive bandwidths, plus a uniform-in-footprint field). If the
     ACF/Moran excess is stable across all three, it is not kernel misfit. If it moves, the
     answer is kernel misfit and must be reported as such.
  3. **Observer nuisance regression (S-1(b)).** A network upgrade produces a region-wide
     positive residual with a long correlation time and a coherent spatial pattern — which is
     *precisely* the predicted signature of "weather". **Mandate: regress the leading residual
     EOF's time series against the K-031 station-density field ρ_sta(x,t) and report the
     partial ACF/Moran after removing it.** Without this, K-009 may discover the seismic
     network.
- **Statistic.** (i) pooled temporal ACF of the Anscombe residual, lags 1–52 weeks; (ii) Moran's
  I / variogram at each time step, pooled; (iii) leading EOF variance fraction and its power
  spectrum. SoCal 0.2° × 7 d, 2010–2018; world 0.5° × 30 d.
- **Null.** As fixed above; 500 simulated catalogs.
- **Success.** Lag-1 weekly ACF excess ≥ 0.05 over the sim-null 97.5th percentile, **and**
  stable across the kernel swap, **and** surviving the ρ_sta partial. Report the correlation
  time and correlation length with CIs — those two numbers, not a p-value, are the deliverable.
  **Failure:** excess inside the null envelope, or destroyed by the kernel swap or the ρ_sta
  partial.
- **Positive control.** Inject a synthetic latent field (a smooth OU-driven multiplicative
  anomaly on μ with a known 4-month timescale and 35-km length) into simulated catalogs; the
  pipeline must recover both scales to within a factor of 2. **This is what makes a null
  interpretable, and a null is a very possible outcome.**
- **Negative control.** Pure ETAS-sim (the null).
- **Artifact class.** Model misspecification masquerading as state (fixed by the kernel swap);
  network evolution masquerading as a slow mode (fixed by the ρ_sta partial); border gradients
  (report interior-vs-edge).
- **Data.** 100% on disk.
- **Scope of a win.** "The residual field of a frozen spatio-temporal ETAS in SoCal, 0.2°/7 d,
  2010–2018, has correlation time T and length L exceeding the model's own." That is a
  statement about a *specification for an instrument*, exactly as Kepler says. It is not
  evidence that the latent field is fluids, creep, or slow slip; naming a mechanism requires
  K-012.

### K-010 — TESTABLE-NOW (Tier 1). Tier 2 gated on K-009.
- **Tier 1 (approved now).** μ_t = EWMA of the declustered background rate, τ frozen on train
  regions/period, walk-forward scored.
  **Statistic:** bits/event vs frozen-μ ETAS (B-1/B-2) — correctly the harsh null, not Poisson.
  **Success:** pooled Δ bits > 0 with block-bootstrap CI excluding 0 on SoCal AND ≥ 4 of the 6
  EXP-M holdouts. **Failure:** otherwise.
- **Tier 2 (particle filter) is gated:** run only if K-009 returns a residual correlation time
  materially greater than the model's own. If innovations are white, the OU state has nothing to
  track and 5,000 particles is compute spent to confirm K-009. Kepler's own internal-consistency
  check (OU timescale vs K-009 correlation time) presumes K-009 first; make the dependency
  explicit in the run order.
- **Negative control.** Kepler's time-shifted-stream "cheating control" — adopt as written, it
  is the correct design.
- **Positive control.** Simulate a catalog from ETAS with a *known* OU-varying μ; Tier 1 must
  gain bits and Tier 2 must recover θ and σ.
- **The interpretation trap, and it must be pre-registered.** μ_t estimated by declustering
  under a frozen model will rise after bursts, because model misfit in the triggering term is
  absorbed into "background". A Tier-1 win is then partly "the ETAS productivity is
  under-fitted", not "there is an external slowly-varying driver". **Mandate the decomposition:
  report Tier 1's gain separately for (a) windows within 30 d of an M≥5 and (b) windows with no
  M≥5 in the trailing 90 d.** If the gain lives entirely in (a), the finding is a triggering-term
  correction and must be described as one. That is still a real app upgrade — just a different
  claim.
- **Artifact class.** Model-misfit absorption (above); overfitting via τ (handled by
  train-only freezing + walk-forward).
- **Data.** On disk.

### K-011 — TESTABLE-NOW as a **measurement**, with a hard scope line.
Kepler labels it correctly ("not a hypothesis test — a measurement") and pre-empts the objection
himself. I accept it and I am tightening the scope line, because this number will be quoted out
of context the moment it exists.

- **Mandatory scope line, to appear verbatim in the results file and in any public statement:**
  *"This is the predictability horizon OF THE FROZEN ETAS MODEL, not of the Earth. It is an
  upper bound on the useful lead time of the current backbone and a lower bound on nothing."*
- **Statistic.** Saturation time (lead time at which RMS ensemble spread reaches 1/√2 of
  climatological spread) as a surface over coarse-graining scales (space 0.2°/1°/box; time
  1/7/30 d), from 1,000-member ensembles; plus the identical-twin variant with perturbed initial
  history (drop/add smallest 1%, jitter magnitudes by reported magError).
- **Null.** None required (measurement); reference = climatological spread.
- **Positive control.** A deterministic-in-the-mean toy (e.g. pure Omori with no branching) must
  give a much longer saturation time than full ETAS. If it does not, the diagnostic is not
  measuring what it claims.
- **The cross-check that makes this more than self-inspection (add it).** Overlay K-011's
  model-internal saturation surface with K-027's *empirical* out-of-sample skill surface on the
  same axes. **If real skill decays faster than the model's internal spread predicts, the model
  is over-confident, and the size of the gap is a misspecification measure.** That comparison is
  the only way this entry says anything about the Earth, and it is nearly free once both exist.
- **Artifact class.** None serious; the risk is rhetorical, hence the scope line.
- **Data.** On disk.

### K-012 — NEEDS-DATA, and **the positive control runs first or the entry does not run at all.**
Kepler's own asymmetry argument is the strongest part of this entry and I am promoting it from a
design feature to a gate: if the method cannot detect the coupling in Cascadia, where ETS
transients are enormous and well-catalogued, then a SoCal null means nothing and the compute is
wasted. Sequence it accordingly.

- **Missing data + cheapest path.**
  1. **NGL daily .tenv3 series**, IGS14, for stations in the SoCal box — station list at
     `geodesy.unr.edu/NGLStationPages/llh.out`, series at
     `geodesy.unr.edu/gps_timeseries/tenv3/IGS14/<STA>.tenv3`. ~1,000 small files; a scripted
     loop with rate-limiting, one evening. We currently hold only `data/ngl/midas.IGS14.txt`
     (secular velocities) — Kepler's observation that "the negative space in our own data
     inventory is the time series behind the velocity" is exactly right and well spotted.
  2. **PNSN / Cascadia tremor catalogue** + Cascadia NGL stations + a Cascadia seismicity
     catalogue — for the positive control. Also small and public.
  3. Optional later: GCMT, GRACE/GLDAS.
- **Gate.** Run the Cascadia arm first. **Success of the gate:** the pipeline detects the known
  ETS transient–tremor coupling at the pre-registered lag with the circular-shift null. **If the
  gate fails, K-012 is closed** and the recorded finding is "the transient-detection pipeline
  cannot see coupling even where it is known to be large", which is a methods result, not a
  seismology result.
- **Statistic (SoCal arm, post-gate).** Cross-correlation of weekly transient shear-strain-rate
  anomaly with the K-009 residual field, lags −12..+12 weeks, per cell and pooled; plus bits/
  event from adding it as a multiplicative covariate on μ_t in K-010.
- **Null.** Circular time-shift of the transient field against the residual field (2–2000 d) —
  the same machinery validated in round 1's protocol — plus a station-label-permuted spatial
  control. Adopt both as written.
- **Success/failure.** Success = pooled cross-correlation outside the shift-null envelope at a
  pre-registered lag band AND positive bits in K-010. Failure = otherwise, in which case the
  deliverable is the **measured detectability floor** ("transient strain must exceed X
  nanostrain/yr to be detectable as an ETAS covariate at current station density"). I accept
  Kepler's argument that this is a real deliverable from a null — it is a specification for a
  network, and it is the kind of bounded-negative K-032 item 6 is built for.
- **Artifact class.** Common-mode GNSS noise (hence the common-mode filter, which must be fit on
  a reference frame excluding the SoCal box), seasonal hydrologic loading aliasing into
  "transients" (hence removing annual + semiannual in the trajectory model — and note this makes
  the entry unable to test annual-loading hypotheses, which is fine; that corpse is EXP-F's).
- **Depends on:** K-009 (supplies the residual field).

---

## 4. VERDICTS — Q4 (other negative spaces)

### K-013 — TESTABLE-NOW, with one mandated fix; otherwise it is artifact-guaranteed.
The inversion is excellent — B-4's ledger at 29-year resolution becomes computable at 7-day
resolution because ETAS supplies the expectation. Four orders of magnitude faster, same move.
And it is a shippable product: after any mainshock the map issues within days.

**The fix, and without it the entry cannot fail honestly.** Near-source cells lose events to the
coda in the first days (short-term aftershock incompleteness), so obs is low, so
d = (exp − obs)/√exp is HIGH near the source. Late large aftershocks also occur near the source
and at rupture terminations. **Deficit will therefore predict late large events for purely
instrumental reasons.** Per S-1(b) the ETAS-sim null cannot catch this: simulated catalogs are
perfectly observed.
- **Mandate (either, preferably both):** (i) compute the deficit map over days **3–30**, not
  1–30; (ii) apply the K-015 detection function. And run the **diagnostic that settles it**:
  compute the statistic with days 0–3 included and excluded. If the signal lives only in the
  included version, the finding is the coda, and it should be recorded as a (useful) replication
  of K-015 rather than as a stress-shadow result.
- **Statistic.** AUC and Poisson-regression Δ LL for the location of M≥4 events in days 31–365,
  with deficit d as predictor, **above a baseline that already uses observed days-1–30 activity**
  (Kepler specifies the harsh baseline — good, keep it) and above frozen ETAS λ per S-3.
- **Null.** ETAS-sim sequences with the identical deficit pipeline and the identical
  low-count noise; plus a **cell-label permutation within observed-count deciles** (S-4 — the
  deficit is loudest where counts are smallest).
- **Success.** Δ bits > 0 over the observed-activity baseline with block bootstrap CI excluding
  0, in **both** the SoCal and the world arms. **Failure:** otherwise, or signal confined to
  days 0–3.
- **Positive control.** Simulate sequences with a genuine locked patch (a zone of suppressed
  triggering that later hosts a large event); the pipeline must recover it.
- **Negative control.** Plain ETAS-sim sequences (no locked patch): AUC must be ≈ 0.5 after the
  count-decile matching.
- **Artifact class.** Short-term incompleteness (primary, fixed above); n-bias in low-count
  cells; border truncation of aftershock zones (exclude mainshocks whose 30-d expected-count
  disc crosses the box edge, and report the excluded fraction).
- **Depends on:** K-015 (preferred) and K-002 (for the spatial kernel; G4).

### K-014 — TESTABLE-NOW. Cheap, novel, and nobody has the plot.
Kepler's defence is right: the time-reversal gap is a **model-free-ish upper bound on precursory
information in the event stream**, obtained without positing a precursor mechanism. That is a
genuinely useful object and the magnitude-stratified version is the prize.

- **One interpretive correction to pre-register.** Forward-minus-backward is not purely "causal
  information". The ETAS functional form is itself asymmetric, and inverse-Omori foreshock rates
  are a documented empirical fact, so a nonzero gap is partly a statement about the model class.
  **Mandate a third arm: a symmetric (acausal, two-sided) triggering kernel fit and scored the
  same way.** Forward vs backward vs symmetric is a three-way comparison that separates
  "the Earth is causal" from "my kernel is causal".
- **Statistic.** Bits/event forward vs reversed, per magnitude stratum (shared strata with
  K-001), 13 boxes + SoCal; plus reversed-fit parameters (backward p, backward α).
- **Null.** For a purely symmetric process the difference is 0; the sampling distribution of the
  difference comes from ETAS-sim (strictly causal, therefore setting the *maximum* expected
  asymmetry) — Kepler's construction, adopt.
- **Success.** Not a pass/fail — it is a measurement with a CI. **Pre-register the headline
  quantity:** the ratio (backward bits)/(forward bits) as a function of magnitude stratum, with
  block-bootstrap CIs. Pre-register the interpretation of both directions so neither can be
  narrated after the fact.
- **Positive control.** On ETAS-sim data the backward/forward ratio must be well below 1 at
  every stratum.
- **Negative control.** On a *symmetric* simulated process (Neyman–Scott cluster process with
  symmetric offspring in time), the ratio must be ≈ 1 at every stratum. This is the control that
  proves the statistic measures asymmetry and not something else, and it is missing from
  Kepler's spec. Mandatory.
- **Artifact class.** Edge effects at the ends of the catalog (a reversed catalog has its
  burn-in at the other end — mandate identical burn-in lengths in both directions); block
  dependence in the high-M strata (block bootstrap).
- **Data.** On disk. Hours, not days.

### K-015 — TESTABLE-NOW. **Raise its priority: K-006R and K-013 both depend on it.**
Kepler's reframe — detection loss as a *measurement channel* rather than a nuisance to correct —
is the good idea here, and pairing it with a bits-scored forecast test rather than a
parameter-recovery exercise is what makes it more than the existing literature.

- **Statistic 1.** Is Mc(t) − Mc_background well fit by g·log10(rate) with a stable g?
  (Ogata–Katsura joint estimator as primary; maximum-curvature as a cross-check only — MAXC is
  n-biased at small n, per S-4.)
- **Statistic 2 (the payoff).** Bits/event for the day-1–7 M≥4 count forecast using (a) observed
  above-Mc counts vs (b) counts corrected by the estimated detection function. Walk-forward,
  train/test split at 2010.
- **Null.** Mainshock-shuffled control (one event's Mc(t) curve applied to another's rate
  history) — adopt; plus ETAS-sim with a *simulated detection function of known parameters*.
- **Success.** Δ bits > 0 in the first 48 h with block-bootstrap CI excluding 0, AND ĝ recovered
  within its CI on the simulated arm. **Failure:** either fails.
- **Positive control.** The simulated-detection-function arm (recovery of known g) — Kepler has
  it; it is exactly right and it is the reason this entry is trustworthy.
- **Negative control.** Apply the correction to a period with *no* mainshock; bits must be ≈ 0
  (the correction must not manufacture skill in quiet times).
- **Artifact class.** n-dependent Mc estimator bias; magnitude-scale drift across SCSN eras
  (per-decade Mc 1.7/1.7/1.6/1.2 in results_exp_h.json — fit g within era, not across).
- **Data.** On disk.
- **Downstream value.** This is infrastructure: it unblocks K-006R and de-artifacts K-013, and
  it improves the operational window (first 48 h) where forecasts currently fail hardest.

### K-016 — REFRAMED → **K-016R** (credit K-016, Kepler). Drop the common-mode arm; keep the triggering arm.
Kepler flags this himself as the round's lowest-prior entry and he is right; I am killing half of
it and keeping the good half.

**Why the common-mode arm dies:** the leading eigenvalue of a 13×13 sample correlation matrix
built from sparse monthly M≥6.5 counts (517 events over 380 months across 13 regions — most
cells are 0) is **severely upward-biased** by construction; that is the Marchenko–Pastur regime,
and it is S-4 in matrix form. The ETAS-sim null does absorb the bias (same estimator, same n) —
which means the test is not *wrong*, it is *powerless*: the signal Kepler would need is far
below the noise the estimator manufactures. Running it produces an uninformative null that
someone will later cite as "no global mode", which is worse than not running it. If it is ever
worth doing, it needs continuous regional intensity series (M≥4.5 residuals from K-009), not
M≥6.5 counts.

**K-016R — dynamic triggering only, which is powered and pre-registerable:**
- **Statistic.** Global M≥6.5 rate in 0–24 h and 0–7 d after each M≥8.0, **outside the source
  region** (exclude within 1,000 km and within the source region's aftershock zone), vs the
  ETAS expectation for that window. Rate ratio with a Poisson-mixture CI.
- **The dose-response arm that makes it a real test.** Order the triggering events by predicted
  surface-wave amplitude at each candidate receiver (a simple 1/√distance × M scaling is
  sufficient). **A real dynamic-triggering signal must show a monotone dose-response; a
  coincidence will not.** Pre-register the rank correlation between predicted amplitude and
  excess rate as the confirmatory statistic, with the raw rate ratio as secondary.
- **Null.** Randomized trigger times (draw the same number of "trigger" epochs uniformly,
  preserving the seasonal and secular rate structure), 2,000 draws, with the identical exclusion
  geometry.
- **Success.** Dose-response rank correlation outside the null envelope at p < 0.01.
  **Failure:** otherwise — and a clean null here is a genuinely useful boundary on the "internal
  weather" framing, which was Kepler's stated purpose for the entry.
- **Artifact class.** Shared completeness trends (Kepler names this correctly — it is why the
  common-mode arm was hopeless); post-large-event reporting surges (a global M8 raises attention
  and analyst throughput — a real and underrated effect; restrict to M≥6.5 targets, which are
  never missed, to defuse it — note this is a *reason* the M≥6.5 threshold is right).
- **Data.** On disk.

### K-017 — TESTABLE-NOW. This is a new hypothesis, not a rerun of the EXP-F corpse. I record that plainly.
Kepler's argument is correct and I accept it against my own corpse list: EXP-F tested fixed
*calendar* periods with a Rayleigh statistic and smoothed-rate surrogates, and it died to rate
fluctuations (§15: three multi-year train detections killed by the frozen phase-agreement rule =
"slow rate fluctuations, not cycles"). Changing the independent variable from calendar time to
accumulated-deformation time is precisely a change that removes that failure mode. Different
hypothesis. Approved.

**But the proposed clock is self-defeating and must be changed.** N(t) = cumulative count of
M≥4.5 events. Large events *produce* large numbers of M≥4.5 aftershocks, so each M≥6.5 advances
its own clock enormously. That mechanically regularises the count-gaps between successive M≥6.5
and manufactures CV < 1. The ETAS-sim null does control it (the sim has the same mechanism, and
here the sim is non-circular because the question is whether the Earth has a loading clock ETAS
lacks) — but controlling an artifact that dominates the statistic leaves no power.
- **Mandate: the primary clock is a DECLUSTERED count-clock** (background events only, Zaliapin
  parameters as used elsewhere in the program) or **cumulative √moment of background events**.
  The raw-count clock is retained as a secondary arm precisely to demonstrate the artifact.
- **Statistic.** CV of interevent *clock* gaps between successive M≥6.5 (world) / M≥5 on an
  M≥2.5 clock (SoCal, ~130 mainshocks — much the better-powered arm, and EXP-I's mainshock count
  confirms the n); plus Brownian-passage-time / lognormal vs exponential likelihood ratio.
- **Null.** ETAS-sim, identical pipeline, identical clock construction, matched number of gaps.
- **Success.** Pooled CV below the sim-null 2.5th percentile **on the declustered clock**, with
  a bias-corrected CV estimator (CV from 5–40 gaps is badly biased low — S-4; use the same
  estimator on real and sim so the bias cancels, and say so). **Failure:** otherwise.
- **Positive control.** Simulate a genuine BPT renewal process with CV = 0.6 plus ETAS
  aftershocks on top; the pipeline must recover CV < 1 at the real number of gaps.
  **Without this the null is uninterpretable.**
- **Negative control.** Plain ETAS-sim (the null).
- **Artifact class.** The clock-advancement mechanism (fixed above); completeness changes
  altering the clock rate over time (**the clock must be built above a fixed, era-stable Mc** —
  per S-1(b) the sim cannot null this).
- **Scope of a win.** A pooled CV significantly below the sim null would mean the system has a
  memory variable ETAS lacks. That is a major result and it would immediately supply a hazard
  term — so the positive control and the declustered clock are not optional bureaucracy, they
  are what would make anyone believe it.

---

## 5. VERDICTS — the emergence lens (K-018..K-027)

*Read S-3 first: the shared gate for this family is incremental bits over frozen ETAS λ, not
standalone AUC. Read S-2: four of these are rolling-window designs and this program has already
been fooled once by window overlap (§14b). K-018, K-019R, K-020R, K-021R and K-026 run as ONE
joint job (G2), producing one incremental-bits table, not five studies.*

### K-018 — TESTABLE-NOW. Head of G1 and of G2; runs after K-005.
The physics framing is good — whether the crust sits at a fixed distance from criticality or
breathes around it is a real question, and a null is genuinely informative (it would say n is
constant, which simplifies every forecast).

- **Statistic.** n(t) in rolling 1-yr windows stepped 1 month, two estimators: (a) refit K with
  α, c, p frozen; **(b) the model-light triggered-fraction from stochastic declustering —
  primary**, because it is far more stable and because EXP-H's own record (nominal n = 1.161
  supercritical vs 7-day n_eff = 0.60) is proof that the analytic branching ratio is
  definition-dependent. Then: incremental bits of [λ + n(t)] over [λ], walk-forward.
- **Null.** ETAS-sim with constant true n, identical estimator, identical windows, identical
  target selection (S-1(c)) — Kepler's construction, and here the instability genuinely does
  cancel between signal and null, as he argues.
- **Success.** Incremental bits > 0 with block-bootstrap CI excluding 0 on a temporal holdout.
  AUC and pre-event Δn reported as secondary. **Failure:** CI covering 0.
- **Overlap discipline (S-2).** 92% window overlap. Confirmatory arm = non-overlapping annual
  windows; the monthly-stepped series is descriptive.
- **G1 dependency.** n = K·b/(b−α)·(…) inherits α's M0-sensitivity. **K-005's measured bias
  curve must be applied before any n(t) excursion is called physical.**
- **Observer artifact the sim cannot catch (S-1(b)).** n(t) rises trivially when Mc falls, since
  more small events look like more triggering. SCSN's Mc moves 1.7 → 1.2 across the study span.
  **Mandate: compute n(t) above a fixed, worst-era-safe floor (M≥2.5 is already chosen for this
  reason — verify per-era, and report n(t) at M≥3.0 as a robustness arm).**
- **Positive control.** Simulate with a *known* time-varying n (a slow sinusoid between 0.7 and
  0.95); estimator (b) must track it.
- **Negative control.** Constant-n sim (the null).
- **Artifact class.** n-dependent estimator bias; Mc drift; window overlap.

### K-019 — REFRAMED → **K-019R** (credit K-019, Kepler). The EWS battery as proposed is three artifacts in a trench coat.
The import is legitimate in principle and the honest-negative framing ("the crust does not
approach its transitions the way lakes do") would be a real result. But the specific battery
cannot deliver it:
1. **Variance is mechanically coupled to the mean.** For a count series, Var ∝ mean under almost
   any clustering model. "Rising variance before a large event" is "rising rate before a large
   event", which is ETAS. Gaussian-kernel detrending removes the mean but not the mean–variance
   relation.
2. **AR1 of a monthly count series rises when the rate rises**, for the same reason.
3. **Kendall tau of a 3-yr rolling statistic measured over a 3-yr window** is close to a
   deterministic function of the window — S-2 in its purest form, and this program has already
   published a self-correction on exactly this geometry (§14b).
4. Three observables × two statistics, no multiplicity control.
The ETAS-sim null blunts (1)–(3) but only if every selection detail is mirrored, and even then
the test is measuring "is the rate rise bigger than ETAS's rate rise", which is K-018's question
asked less well.

**K-019R, the surviving cousin — one statistic, no battery.** The only critical-slowing-down
quantity that is not mechanically mean-coupled is the **recovery rate of the rate-normalised,
variance-stabilised residual** — which K-009 already computes with a better design and no
target-selection. So:
- **Claim.** The K-009 residual correlation time is **elevated in the 3 yr preceding M≥7.0
  targets** relative to count-matched, rate-matched control windows.
- **Statistic.** Residual temporal correlation time (fit an exponential to the pooled residual
  ACF) in pre-target vs matched-control windows; pooled across the 179 M≥7 targets on disk;
  incremental bits over λ as the S-3 gate.
- **Null.** ETAS-sim through the identical K-009 pipeline with the identical target selection.
- **Controls.** Positive: inject a known slowing (an OU state whose timescale lengthens before
  synthetic targets) and recover it. Negative: plain ETAS-sim. Second negative: shuffle which
  region's targets are applied to which region's series (Kepler's, adopt).
- **Success.** Pre-target correlation time exceeds control by an amount outside the sim-null
  envelope, AND incremental bits > 0. **Failure:** otherwise — and this failure is the
  informative negative Kepler wanted, now attributable to physics rather than to estimator
  coupling.
- **Depends on:** K-009. Costs almost nothing once K-009 exists. That is the point of the
  reframe.

### K-020 — REFRAMED → **K-020R** (credit K-020, Kepler). Standalone claim REJECTED; survives as a count-matched diagnostic inside K-022.
Kepler offers this fallback himself and I am taking it.
**Kill-reason for the standalone claim:** the fitted correlation length ξ from
C(r) ~ r^{−d2} exp(−r/ξ) in rolling 1-yr windows is strongly **n-dependent** — more events
sample the tail better and the fitted ξ grows. Pre-target windows have more events by
construction (rate rises before large events under ETAS). So "ξ grows before large events" is
S-4 plus S-3, and it is the same failure that killed AMR in prospective CSEP tests. Kepler's
sim null with a fixed spatial kernel *would* control it — but only if windows are count-matched,
and the proposed design is not.
**K-020R:** ξ enters the K-022 percolation job as one candidate covariate, computed on
**fixed-n subsamples** (S-4: subsample every window to a common event count, average over 50
draws), scored by the S-3 incremental-bits gate alongside S_max and the other G2 variables. No
independent claim, no independent write-up, no compute of its own.

### K-021 — REFRAMED → **K-021R** (credit K-021, Kepler). The rationale is good; the correction is the wrong instrument.
The statistical-physics argument is real and I want to credit it: susceptibility (a variance)
rather than magnetization (a mean) is what diverges at a critical point, so looking at the
second moment of b is a principled move and is *not* the dead b-mean literature.

**Why the proposed correction fails.** Kepler correctly identifies subtracting the sampling
variance b²/n as "the technical crux" — and then proposes the one form of the correction that
will not hold. Var(b̂) = b²/n is asymptotic; at the n ≈ 100 per cell that the design permits, the
Aki estimator is biased (≈ n/(n−1)) with a heavy right tail, so the analytic subtraction
under-corrects, and the residual is n-dependent. Cell counts *rise* before large events (rate
rises), so corrected Var(b) *falls* before large events — **which is precisely the predicted
signal, manufactured**. This is EXP-B's corpse in new clothes: "all apparent feature-vs-amplitude
correlations are the small-n bias in disguise; the bias-robust label correlated with nothing."

**K-021R — remove the bias by construction instead of correcting it (S-4).**
- **Statistic.** Fix a cell set: cells with ≥ n0 events in **every** window (no cells entering or
  leaving — otherwise the cell set is itself a rate proxy). In each window, subsample every cell
  to exactly n0 = 100 events without replacement, estimate b̂, average Var_cells(b̂) over 50
  independent subsample draws. The sampling floor is then identical in every window by
  construction and any change in Var_cells(b̂) is real heterogeneity change.
- **Null.** Kepler's two sims, both adopted and both necessary: (i) spatially uniform true b (all
  measured heterogeneity is sampling); (ii) a frozen heterogeneous b-field (heterogeneity exists
  but does not evolve). Two nulls separate "there is heterogeneity" from "the heterogeneity
  changes" — that is a well-designed pair and I want it kept.
- **Success.** Decline in fixed-n Var_cells(b̂) in the 1–2 yr before targets, outside the
  envelope of null (ii), AND incremental bits over λ (S-3). **Failure:** otherwise.
- **Positive control.** Simulate a b-field that genuinely homogenises before synthetic targets;
  the pipeline must recover the decline at the real n0 and cell count.
- **Observer artifact (S-1(b)).** b̂ is exquisitely sensitive to Mc, and Mc varies in space and
  time with station density. **Mandate a per-cell, per-era Mc, and a robustness arm on the
  K-031 fixed subnetwork.** No sim can null this.
- **Overlap.** 3-yr rolling windows: confirmatory arm non-overlapping (S-2).

### K-022 — TESTABLE-NOW. I agree with Kepler: highest upside in the emergence lens.
Two reasons it earns the compute. First, it is one of only two entries (with K-008) that attacks
ETAS's i.i.d.-GR magnitude assumption — the assumption that makes "how big" unforecastable by
construction. Second, its deliverable is a sentence no one ships: "the currently connected
backbone cannot support a rupture larger than M x.x."

- **Statistic 1 (order-parameter behaviour).** Distribution of S_max — bimodality / jump vs the
  smooth sim distribution.
- **Statistic 2 (the payoff).** Quantile regression of next-large-event magnitude on log L(S_max),
  out of sample; CRPS gain for Mmax.
- **Four mandated additions:**
  1. **The GR baseline must be the region's train-fitted TAPERED GR, not an unbounded GR.**
     Beating an unbounded GR on Mmax is trivial (it over-predicts the tail) and would be a
     strawman — prong 4 of the standard.
  2. **The circularity covariate.** The previous large event's aftershock zone inflates S_max
     *and* raises the probability of a large late event. **The model must beat a baseline that
     already contains [time since last M≥5, its magnitude, its distance]**, not just GR.
     Without this, "connectivity forecasts Mmax" may be "there was recently a big one".
  3. **Spatial-randomisation control (Kepler's, and it is the load-bearing null).** Preserve the
     number of active cells, destroy their arrangement. This is what isolates geometry from
     amount, and Kepler is right that it is the essential one.
  4. **Detection geography (S-1(b)).** At 0.05° with M≥2.5, cell occupancy is substantially a
     map of where SCSN can see. **Mandate the K-031 fixed-subnetwork arm.** Connectivity of a
     detection footprint is not connectivity of a stressed backbone.
- **Null.** ETAS-sim (spatial clustering, no connectivity–magnitude coupling) + the spatial
  randomisation + the aftershock-zone covariate baseline.
- **Success.** CRPS gain ≥ 5% over the tapered-GR-plus-recent-large-event baseline, out of
  sample, with cluster-block bootstrap CI excluding 0, surviving the spatial randomisation.
  **Failure:** otherwise.
- **Positive control.** Simulate a system with an imposed connectivity–magnitude coupling
  (max rupture ∝ largest connected component); the pipeline must recover it.
- **Artifact class.** Detection geography; the aftershock-zone circularity; grid-resolution
  sensitivity (report S_max at two resolutions — and note K-025R would choose the resolution
  principled rather than by hand).
- **Absorbs:** K-020R and, at the supervisor's discretion, K-026 as covariates in the same job.

### K-023 — NEEDS-DATA, and underpowered at n = 13. Do not run until n ≥ 30 regions.
The idea is the only one in the round that produces a physically grounded regional Mmax, and
Kepler's regime-split prediction (moment ~ L³ unsaturated vs ~L² width-saturated, with
continental transform vs subduction mapping onto them) is a genuine pre-registerable prediction
rather than a curve fit. That is what makes it worth keeping.

- **Missing data + cheapest path.** GCMT (globalcmt.org, 1976–now, homogeneous Mw — a single
  small file, and it defuses the magType objection); ISC-GEM for the long tail; Slab2 for
  down-dip width. All free, one-shot.
- **Design gate.** With 13 regions and corner-magnitude CIs of order ±0.3–0.5, the regression
  has essentially no power to distinguish exponent 2 from exponent 3. **Extend to ≥ 30
  Flinn–Engdahl regions from GCMT before scoring anything.** Running at n = 13 and labelling it
  "suggestive" produces a number that will be quoted without its label. Kepler's own
  suggestive-only framing is honest but insufficient — the fix is more regions, not a softer
  claim.
- **Statistic.** Log-log regression of corner moment on the geometric measure; exponent with CI;
  pre-registered two-regime test (transform vs subduction).
- **The n-bias trap that would manufacture the predicted slope (S-4).** Corner-magnitude
  estimates are **upward-biased with more events**, and bigger fault systems have more events.
  So a positive slope of corner moment on geometric size is partly guaranteed. **Mandate:
  regress out log(n_events), or match regions on n by subsampling to a common catalogue size,
  and report both.** This is the single thing that decides whether the result means anything.
- **Null.** Slope 0; plus a shuffled-geometry control (Kepler's, adopt).
- **Positive control.** Simulate finite-size-limited catalogues with a known L→Mmax scaling;
  recover the exponent at n = 30 regions.
- **Artifact class.** n-bias in corner estimation (primary); magnitude-scale heterogeneity
  (defused by GCMT-only).

### K-024 — REJECTED (the claim; the exploration is welcome and the object may be built as a byproduct).
Kepler flags this as weak himself. I agree, and the kill-reasons are specific:
1. **Power.** SoCal has a handful of M≥6.5 events in the instrumental record. A comparison of
   R(t) in 5-yr pre-target windows against matched controls with fewer than 6 targets cannot
   distinguish anything from anything. This is not a threshold I can lower with a better
   statistic.
2. **The phase variable is largely unmeasured.** φ_j requires time-since-last-event *on that
   segment*; for the large majority of the 557 CFM segments there is no instrumental-era
   characteristic event, so φ_j would be imputed from the moment budget — meaning R(t) is
   substantially a function of the loading model, not of the earthquake record. An order
   parameter computed mostly from assumptions is not a measurement of the system.
3. The UCERF3 paleoseismic extension raises the quality of φ for ~tens of sites but does not
   touch (1); the target count is the binding constraint.
**Not rejected:** the *idea* that faults are interacting units with internal state rather than
points emitting events. Kepler is right that this framing is a prerequisite for physics-based
forecasting, and I want it back in a future round attached to a target set with power (Japan
with GCMT + segment models, or a synthetic-first study establishing what R(t) would have to look
like to be detectable). If the supervisor wants the segment-phase object built as reusable
infrastructure alongside K-003, that is fine — it carries **no claim** and appears in no results
headline.

### K-025 — REFRAMED → **K-025R** (credit K-025, Kepler). Keep the decision, drop the collapse machinery.
The critique of D2 papers is correct and a scaling collapse *is* a stronger test than fitting one
exponent. But the characteristic length this design would find is over-determined by
non-crustal lengths: hypocentre location error (~1 km horizontal, worse in depth, and
era-dependent), the background smoothing kernel, the spatially-variable Mc, and — at large L —
the box size and finite cell count. Three of those four are observer lengths, and per S-1(b) the
ETAS-sim null contains none of them. The collapse residual would also need a null distribution
that the design does not specify.

**K-025R — ask the decision question directly, since that is what the entry is for.** Kepler
says so himself: "the deliverable is a length scale... which then becomes the grid resolution
choice for K-002/K-009/K-022 rather than the arbitrary 0.2° we have been using since round 1."
That is a real gap in the program and it deserves a real answer.
- **Statistic.** For grid sizes {0.05, 0.1, 0.2, 0.4, 0.8, 1.6}°, the **out-of-sample forecast
  skill** (log score and CRPS) of the K-002 spatio-temporal ETAS, walk-forward, against the
  local-oracle climatology at the same resolution. The skill-maximising resolution is the
  program's answer to "at what scale does the crust stop being usefully self-similar",
  operationalised as a decision.
- **Null / controls.** Inherits K-002's (rotated field, spatial holdout, interior-vs-edge).
  Additional negative control: repeat on ETAS-sim data, where the optimum should sit at the
  simulation's own kernel scale — if the real optimum differs, *that difference* is the
  physical length, measured against a model that has a known answer.
- **Success.** A resolution optimum with a bootstrap CI narrower than one grid step.
  **Failure:** a flat skill-vs-resolution curve, which is itself a useful finding (it says the
  0.2° choice costs nothing and the program can stop worrying about it).
- **Depends on:** K-002 (G4). Nearly free once K-002 exists.

### K-026 — TESTABLE-NOW, bundled into G2. Lowest priority in the family; no standalone write-up.
Kepler names the count-matching requirement himself, correctly — entropy is n-biased and every
comparison must be between equal-count windows (subsample to the smaller; S-4). Two additions:
- **Redundancy.** Entropy, S_max (K-022), and ξ (K-020R) all measure spatial concentration; n(t)
  and all of them track rate. **Mandate that entropy be reported only as a partial contribution
  in the G2 joint model** — its incremental bits over λ *and over S_max*. A standalone AUC of
  0.58 for this variable is uninterpretable per S-3.
- **Observer artifact (S-1(b)).** Spatial entropy over a fixed cell set is partly a map of
  detection geography, which changed over the study span; a secular entropy trend is expected
  from network growth alone. Fixed-subnetwork arm (K-031) required for any temporal-trend claim.
- **Kepler's real point survives and I endorse it:** the (n, H) phase-plane trajectory is a
  communicable artifact for the app even at modest AUC — **provided it is labelled as a state
  visualisation and not as a forecast.** If the G2 table shows zero incremental bits, the plot
  may still ship, clearly marked descriptive.

### K-027 — TESTABLE-NOW. **Top-3 priority.** The most valuable measurement in the round.
"At what scale does predictability begin" replaces "can we predict earthquakes" with a question
that has an answer. Every component exists and is validated. I endorse the framing without
reservation, and I endorse Kepler's ranking of it ("if only one thing from the emergence lens
gets run, run this and K-022").

- **Statistic.** Out-of-sample skill of the frozen ETAS (walk-forward, no refit) vs local-oracle
  climatology over the coarse-graining grid: space {0.2°, 1°, 5°, whole box} × time {1, 7, 30,
  90, 365 d} × magnitude {≥4.5, ≥5.5, ≥6.5}. Proper scores: CRPS for counts and for log moment
  release, plus the log score, plus the reliability/sharpness decomposition. SoCal + the 6 EXP-M
  holdouts.
- **Null.** Local-oracle climatology (the harsh, correct baseline this program already uses).
- **Not a pass/fail — a surface with uncertainty.** **Mandate: block-bootstrap CIs over years on
  every cell of the surface, and a stated multiplicity rule — the deliverable is the whole
  surface and its frontier contour; no single cell may be extracted as "the finding".** 60 cells
  with no such rule is how a scoring table becomes a fishing expedition.
- **The n_effective problem (the main artifact).** As coarse-graining increases, the number of
  scored units falls; at (whole-box, 365 d, M≥6.5) there are a handful of units per region. The
  frontier's *location* is therefore uncertain in the direction that matters. **Mandate: plot
  n_effective alongside skill on every cell, and refuse to draw the frontier through cells with
  fewer than 20 effective units.**
- **Second deliverable, nearly free and highly diagnostic:** compute the identical surface on
  ETAS-sim data. **The real-minus-sim surface localises ETAS misspecification in scale space** —
  it tells the next round *where* to look, in units of space, time, and magnitude. That converts
  K-027 from a report into an instrument for generating hypotheses, which is exactly what this
  program needs from Kepler's next round.
- **Pairs with:** K-011 (overlay the model-internal saturation surface; the gap is an
  over-confidence measure).
- **Positive control.** At (0.2°, 1 d, M≥6.5) skill must be ≈ 0; at (box, 90 d, M≥4.5) skill must
  be clearly positive. If those two corners do not behave, the harness is wrong.
- **Data.** 100% on disk.

---

## 6. VERDICTS — frame-breaking entries (K-028..K-032)

### K-028 — TESTABLE-NOW, and the criterion is ADOPTED with one modification. Best entry in the round.
Kepler's demarcation argument is correct and I am adopting it as program policy (see charter
amendment 2): **a property that changes when you change the seismometers is a property of the
seismometers.** The Anza episode (§14b: φ = +158.7° in SCSN, n = 49, vs +25.4° in QTM, n = 31 —
same rock, same decade, opposite conclusion) is this lesson learned retrospectively and at cost.
Promoting it to a precondition is the right call, and I should have proposed it myself.

- **The one modification.** Kepler proposes a frozen numeric bar (I ≥ 0.8, Jaccard ≥ 0.7) and
  then, in his own Null line, states the correct criterion: *"cross-pipeline disagreement only
  counts if it exceeds within-pipeline sampling noise."* Those two are not the same rule, and
  the second is right. A fixed numeric bar chosen before we know the within-pipeline dispersion
  of I is the kind of arbitrary threshold this program otherwise refuses.
  **Adopted form: a claim passes the invariance gate if cross-pipeline disagreement does not
  exceed the within-pipeline bootstrap dispersion of the same statistic by more than a
  pre-registered factor (freeze at 2× before running). I and Jaccard are reported as descriptive
  summaries, not as gates.**
- **Statistic.** I = (cross-pipeline bits/event)/(own-pipeline bits/event) for B-1/B-2; recomputed
  conditional probability for B-3; Jaccard of the SILENT set for B-4 (the machinery exists —
  EXP-K already used Jaccard for the H-sensitivity check and returned 0.949 flat-vs-variable H).
- **Null.** Within-pipeline bootstrap dispersion (Kepler's, adopted as the gate).
- **Implementation note before freezing.** Verify QTM's temporal and spatial coverage; it is a
  template-matched product over a limited span. **All comparisons must be restricted to the
  overlapping space-time volume of all three catalogues**, and refit at a common M0 = 2.5
  (Kepler specifies this — correct, since QTM's Mc is far lower and matching M0 is the only way
  the parameter comparison means anything).
- **Positive control.** Split ONE catalogue in half spatially and cross-score halves: I must be
  ≈ 1. This calibrates what "same pipeline" looks like and is what the gate's factor is measured
  against.
- **Negative control.** Cross-score SCSN against a *shuffled-magnitude* SCSN: I must collapse.
- **Prediction on the record.** Kepler predicts p and α largely invariant, μ wildly
  non-invariant, and B-4's silent list substantially non-invariant. I share that prediction and
  I note it is a prediction about our own most-cited result. **Recording it in advance is what
  makes the outcome informative either way** — that is the entry's real methodological
  contribution and it is worth more than the numbers it will produce.
- **Runs as one job with K-031.**

### K-029 — TESTABLE-NOW, with the direction of the bound corrected. This correction is not optional.
The lens is right and genuinely novel here: log-likelihood is codelength, prediction and
compression are the same operation in the same unit, and the ETAS-simulated arm to calibrate the
learner's own handicap is an elegant and correct control. Nobody has quoted this number for an
earthquake catalogue. I want it.

**The correction.** Kepler writes the null as *"Delta ≥ 0 … that null is exactly 'ETAS is a
sufficient statistic for the catalog'."* It is not. A compressor that fails to beat ETAS
demonstrates only that *those* compressors, at *that* sample size, under *that* encoding, found
nothing. No universal compressor is near the entropy rate at 10^6 tokens; failure to find
structure is not evidence of its absence. **Therefore:**
- **The deliverable is a LOWER bound on the residual predictability budget, never an upper
  bound.** Δ < 0 with a bootstrap CI excluding 0 ⇒ "at least |Δ| bits/event of structure exists
  beyond ETAS". Δ ≥ 0 ⇒ "these learners, at this scale, found none" — a bounded negative, and a
  useful one, but not a ceiling.
- This matters far beyond K-029, because **K-032 builds the whole Predictability Budget on
  item 3 being a ceiling.** Taking the tightest of three encodings (Kepler's proposal) is valid
  precisely *because* it is a lower bound — maximising over encodings is the right operation for
  a lower bound and the wrong one for an upper. Fix the label and the method is already correct.
- **Sample-size honesty.** SoCal M≥2.5, 1981–2026 ≈ 7×10^5 events ≈ 2×10^6 tokens. A
  few-million-parameter transformer is data-starved at that scale and will likely lose to both
  LZMA and ETAS. **Pre-register that expectation** so a loss is not over-read as "no structure",
  and report the ETAS-sim arm's Δ for the same architecture as the calibration of the handicap
  (that is exactly what the sim arm is for, and Kepler has it right).
- **Statistic.** Δ = [bits(compressor) − bits(ETAS)]_real − [bits(compressor) − bits(ETAS)]_sim,
  per encoding; bootstrap CI over held-out blocks.
- **Success.** Δ < 0, CI excluding 0, on ≥ 2 of 3 encodings ⇒ report |Δ| as the lower-bound
  residual budget. **Failure:** CI covering 0 on all three ⇒ report the *detectable* budget
  bound at 80% power (which requires a power calibration — inject synthetic extra structure of
  known bits into simulated catalogues and find the smallest Δ these learners recover).
  **That power calibration is mandatory; without it a null is an unquotable number.**
- **Positive control.** The injected-structure arm above.
- **Negative control.** The ETAS-sim arm (Δ must be ≈ 0 or positive there by construction).
- **Leakage risk.** The encoding is fitted (quantisation levels, Hilbert ordering, CTW depth) —
  **all encoding choices frozen on train only**, and the same frozen encoding applied to real
  and simulated streams. A tuned encoding is the one way this test manufactures a budget.
- **Artifact class.** Encoding dependence (mitigated: three frozen encodings, bound direction);
  train/test contamination in the tokenizer.
- **Data.** On disk.

### K-030 — TESTABLE-NOW, with one mandated conditioning; otherwise Statistic 2 cannot fail.
Fragility as a *product* rather than a robustness check is a good inversion, and the observation
that a network planner would rather have a map of where observation has leverage than a hazard
map is correct and worth acting on.

**The mandated fix.** φ (sensitivity of a forecast to dropping one event) is dominated by 1/n:
cells with few events have forecasts that swing most. Cross-pipeline forecast disagreement
(K-028) is *also* dominated by 1/n. So Spearman ρ(φ, disagreement) > 0 with p < 0.01 is
essentially guaranteed and essentially vacuous — S-4 again.
- **Mandate: the confirmatory statistic is the PARTIAL Spearman controlling for log(cell event
  count), or equivalently the within-n-decile pooled correlation.** Pre-register that form.
- **Mandate the comparison that decides whether the product has value: φ vs the trivial
  1/√n_cell layer.** The claim "φ is a deployable per-cell confidence layer computable anywhere
  without three catalogues" is only true if φ beats 1/√n. If it does not, the deliverable is
  "use 1/√n" — which is a perfectly good, cheap, honest outcome and should be pre-registered as
  such.
- **Statistic 1 (descriptive).** The φ map and its relation to event density and the EXP-K
  silent list. Useful, non-confirmatory, label it so.
- **Null.** φ computed on ETAS-simulated histories (expected fragility from sampling alone) +
  cell-label permutation within n-deciles.
- **Success.** Partial ρ > 0 at p < 0.01 **and** φ beats 1/√n in predicting cross-pipeline
  disagreement (nested model comparison). **Failure:** otherwise.
- **Positive control.** Construct cells with artificially inflated pipeline disagreement (inject
  detection differences); φ must rank them high.
- **Depends on:** K-028 (supplies the disagreement field), K-002 (spatial forecasts; G4).
- **The internal-validation point is well made and I endorse it:** if K-028 and K-030 agree
  cell-by-cell, that is genuine mutual corroboration of two independent constructions.

### K-031 — TESTABLE-NOW. **This is the entry I most want run.** A one-command download.
"A research program that cannot attack its own baselines is not a research program, it is an
advocacy campaign." That sentence is correct and I am adopting the structural version of it as
charter amendment 3. This entry aims at two of our own claims with a cheap, decisive instrument,
and the program's own record already flags both targets: EXP-I(iii)'s b-decline was frozen and
confirmed *with the explicit caveat that magnitude-scale homogeneity across network eras was
unchecked* (§15), and EXP-K's unexplained-silent list is **79% detection-limited (n < 20)** —
in the raw file, 158 of 200 cells have n_train < 20, 39 have n_train = 0, and **95 of 200 have
n_test = 0**, so χ_test is undetermined for nearly half the list.

- **Data.** SCEDC/IRIS FDSN station metadata for the SoCal box with operational start/end dates
  (`fdsnws/station`, one query, public). Log the SHA-256.
- **Statistic (i).** ρ_sta(x,t) station-density field.
- **Statistic (ii) — the b-decline.** Recompute the EXP-I(iii) trend three ways: full catalogue
  (the original), **fixed subnetwork**, and M≥3.0 floor. Report slopes with CIs side by side.
  **Success rule: the frozen-sign monotone decline survives all three with CIs excluding 0.**
  Note the original rests on **5 windows** (two partials excluded), Spearman ρ = −0.90,
  p = 0.037 — that is thin, and it is a further reason to run this.
- **Implementation risk, named, with a fallback.** The Zenodo SCSN file does not carry per-event
  station lists, so a literal "only events recorded by the full-span subnetwork" catalogue may
  not be constructible. **Fallback (pre-register whichever is used): restrict to cells whose
  station density was in the top tercile in EVERY era, and/or impose a spatially varying Mc(x,t)
  computed from the full-span subnetwork geometry alone.**
- **Statistic (iii) — the magnitude-scale arm, which I am adding.** The b-decline is equally
  consistent with ML revision drift. **Recompute b using a single magType, and compare SCSN ML
  against ComCat/GCMT Mw for overlapping M≥4 events by era (regression slope and offset per
  era).** Magnitude-scale drift is a named corpse class in this program; leaving it unchecked
  here would repeat the omission the original result already flagged.
- **Statistic (iv) — the ledger.** Regress the silent-cell indicator on ρ_sta(x, t_train);
  report how much of the silent geography is predicted by where the seismometers were not.
  **Report separately for the 158 detection-limited cells and the 42 measured-low-χ cells.**
- **Null.** No attenuation of the trend under the fixed-subnetwork restriction; no association
  between silence and station density (cell-label permutation within loading-matched strata).
- **Positive control.** Simulate a catalogue with a *known* network-growth-induced apparent
  b-drift; the fixed-subnetwork restriction must remove it.
- **Negative control.** Apply the fixed-subnetwork restriction to a period of stable network
  coverage; the trend estimate must not move.
- **Consequence, stated in advance (both directions).** If the b-decline attenuates
  substantially, EXP-I(iii) is NARROWED or RETIRED under the mechanism below. If the silent list
  loses its detection-limited tail but the 1857 Big Bend strand cells (34.4–34.8 / −117.5..−119.1,
  d_fault 1–5 km) survive a station-density control, **B-4 becomes materially stronger than it is
  today** and the Big Bend result becomes quotable rather than caveated. I want both outcomes and
  I have no preference between them.

### K-032 — REFRAMED (see §9 for the full ruling). Adopt the budget as an instrument; reject it as the definition of success. One clause must be deleted.

---

## 7. PRIORITY ORDER — the five to run first

Ranked by decision value per compute-hour. Each states what changes depending on the outcome.

**1. K-009 — ETAS residual whiteness (SoCal + world; 100% on disk, with the three fixes in §3).**
*Decision value:* this is the go/no-go for the entire data-assimilation thread and for most of
Jim's Q3. **White innovations ⇒ close K-010 Tier 2 and K-012, and tell Jim plainly that compute
buys nothing beyond the current backbone** — a large, cheap, honest negative that saves the
program a year of investment. **Red, spatially-coherent innovations ⇒ we get a correlation time
and a correlation length**, which size the filter, set the grid for K-002/K-022/K-025R, justify
the GNSS acquisition, and convert "ETAS is incomplete" into a specification for an instrument.
No other entry changes this many downstream decisions for the compute.

**2. K-031 + K-028 — the self-audit, run as one job (one FDSN query; three catalogues on disk).**
*Decision value:* decides whether B-4's silent list and EXP-I(iii)'s b-decline are rock or
network — before anything else is built on them. **If they attenuate**, two current claims narrow
or retire, K-003's global ledger is re-scoped before its downloads are spent, and the app's
layer-3 candidate (ledger geography as a spatial covariate) is pulled from the roadmap.
**If they survive**, B-4 becomes the program's strongest result, the 1857 Big Bend strand becomes
quotable without its caveat, and the detector-invariance gate ships as validated policy with a
worked example. Cheapest possible way to either harden or kill two claims.

**3. K-027 — the predictability-vs-coarse-graining surface (100% on disk).**
*Decision value:* decides **what product Jim ships.** If the frontier sits at county-and-season
scale, there is a useful forecast to build and we know its resolution; if it sits at
continent-and-decade, we stop promising one and say so publicly. Second, the real-minus-sim
surface localises ETAS misspecification in scale space, which directs Kepler's next round in
units of space, time, and magnitude rather than by enthusiasm. Third, it is item 1 of the
Predictability Budget, so it is load-bearing for §9 regardless of §9's outcome.

**4. K-005 — the M0-invariance audit (on disk; runs before K-006R and K-018, per G1).**
*Decision value:* decides whether EXP-M's parameter table means what §17 says it means.
**If the α drift matches the simulated estimator bias**, then SoCal α = 0.537 vs global
α = 0.730 is an artifact of the magnitude floor rather than tectonics, this program stops quoting
transferred parameters across catalogues with different M0 (which it currently does), and
K-006R's "bend" and K-018's "n excursions" must both be bias-corrected before they can mean
anything. **If the drift exceeds the simulated bias**, crustal self-similarity breaks measurably
and large-event probabilities get a calibration correction with an operational size attached.
Cheap, and it gates three other entries.

**5. K-001 — magnitude-stratified skill to M≥6.5 (re-scoring of runs we already own).**
*Decision value:* decides whether B-1 and B-2 reach the magnitudes anyone actually funds.
**Rising slope ⇒ the validated backbone is most useful exactly where the stakes are**, and the
app gains a defensible large-event layer with a quotable bits figure. **Collapsing slope above
~M6 ⇒ large events are drawn from a different process, B-1/B-2 do not extrapolate, and the
program's foundation does not reach its target** — which Kepler correctly identifies as the most
important possible negative result here. Nearly free: it is a re-scoring, not a new fit.

*Next in queue, named so the ordering is legible:* **6. K-002** (expensive, but it gates the
entire spatial slate under G4 and nothing spatial may quote bits without it); **7. K-022**
(largest upside in the emergence lens, absorbing K-020R and K-026 as G2 covariates);
**8. K-015** (unblocks K-006R, de-artifacts K-013); **9. K-029** (with the corrected bound
direction); **10. K-032 item 6** — the retroactive corpse-to-upper-bound table, which is cheap,
uses only files on disk, and is the part of K-032 I endorse without qualification.

---

## 8. RULINGS ON THE EIGHT CHARTER AMENDMENTS

1. **Standing null upgrade (ETAS-sim, never Poisson) — ADOPT-MODIFIED.** Adopt as default; add
   the two riders in S-1: (a) where the hypothesis is that ETAS is misspecified in the dimension
   being simulated, the sim null is circular and must be replaced by two-generator
   discrimination; (b) simulated catalogues are perfectly observed and therefore **cannot null an
   observer artifact** — pass sims through a detection function whenever completeness, Mc, b,
   station density, magnitude scale, or box edges are on the causal path. Also: the null inherits
   the selection rule.
2. **Detector-invariance as a graduation criterion — ADOPT-MODIFIED.** Adopt the gate; replace
   the fixed numeric bar (I ≥ 0.8 / Jaccard ≥ 0.7) with Kepler's own and better null:
   cross-pipeline disagreement must not exceed within-pipeline bootstrap dispersion by more than
   a pre-registered factor (freeze at 2×). Where no second pipeline exists, substitute an
   independent catalogue authority (ISC / ComCat / GCMT) or record an explicit exemption that
   names the untested exposure.
3. **Adversarial self-audit quota — ADOPT** as written, plus two additions: if Kepler does not
   name a target, I will; and the audit entry takes **priority in the run queue over new
   exploration**, because retiring a false baseline saves everything built on top of it.
4. **Bits-per-compute-hour ranking — ADOPT-MODIFIED.** Require cost and expected bits, but
   **prior-weight the payoff** (a claimed +0.3 bits at 10% prior is worth less than +0.1 bits at
   60%), and require the named decision from amendment 5 in the same line. Kepler must state his
   honest prior probability of the null; he already does this informally and well (K-016's
   self-flag, K-012's "honestly, near-null"), so this formalises existing good practice.
5. **Name the decision, not just the statistic — ADOPT**, unmodified. It is prong 4 of my own
   standard restated from the generative side, and Kepler's rationale (that it would have caught
   round 1's tidal work being unactionable even if it had succeeded) is exactly right.
6. **The lens ledger — ADOPT-MODIFIED.** Adopt the tagging and the running tally; adopt the
   "one never-before-used lens per round" requirement (cheap and generative — this round's
   emergence and frame-breaking lenses produced K-022, K-027, K-028, K-029, K-031, which is the
   argument). **Reject the automatic retirement rule:** with two or three outcomes per lens the
   tally is n-biased noise (S-4 applied to Kepler's own bookkeeping), and a lens that reliably
   produces clean, well-powered nulls is doing the program's most valuable work — which is
   Kepler's own K-032 argument, and he cannot have it both ways. Retire a lens only on ≥ 5
   scored outcomes **and** a stated mechanism for why it fails here.
7. **The right to propose retirement — ADOPT.** Kepler's "one-way ratchets are how research
   programs acquire their unfalsifiable core" is the sharpest sentence in the amendment set and
   it is aimed at me. He is right: I have been gating entry and not exit. Mechanism below.
8. **A standing budget line — ADOPT-MODIFIED.** Maintain the Predictability Budget as a living
   ledger section, with item 3 relabelled as a **lower** bound (per K-029), and with the budget
   defined as a **report and a ranking device, not a success criterion** (per §9).

### The baseline-retirement mechanism (owed under amendment 7)

Adopted as program policy, effective this round.

1. **Every BASELINE carries two new lines: a SCOPE line and a STANDING EXPOSURE list** — the
   named artifact classes it has *not* been tested against. B-1..B-5 receive theirs
   retroactively; drafts in §10.
2. **Anyone may file a CHALLENGE** (Kepler, the supervisor, Jim, or me) naming (a) the baseline,
   (b) the specific exposure, (c) a frozen test whose failure retires it. A challenge is
   adjudicated exactly as a new claim: pre-registered, controlled, out-of-sample. **Nothing is
   retired on argument alone.** Retirement requires a run test, symmetric with promotion.
3. **Four outcomes.** **UPHELD** — status kept, and the tested exposure moves from the exposure
   list to a SURVIVED list, so the baseline is *stronger* for having been attacked (this is the
   incentive that keeps the mechanism from being purely destructive). **NARROWED** — scope line
   amended (e.g. "B-4 holds for measured-low-χ cells, not for detection-limited cells").
   **DEMOTED to PROVISIONAL** — may still be used, must be labelled, may not be the foundation
   of a new baseline. **RETIRED** — moved to the corpse list *with its effect-size upper bound
   recorded* per K-032 item 6, because a retired claim with a bound is an asset and a retired
   claim without one is a rumour.
4. **Automatic challenge triggers, so nobody has to be brave.**
   (a) Failing the detector-invariance gate (amendment 2) ⇒ auto-DEMOTE pending re-test.
   (b) Any baseline whose supporting statistic rests materially on a stratum with n < 20
   ⇒ auto-FLAG for challenge. **This fires today on B-4** (EXP-K: 158/200 unexplained-silent
   cells with n_train < 20, 95/200 with n_test = 0) **and on B-1's Caribbean component**
   (n = 235, `underpowered: true`).
   (c) A baseline unreplicated on a second region, catalogue, or period after two rounds
   ⇒ auto-DEMOTE to PROVISIONAL. **This is the anti-ratchet: status decays unless renewed.**
5. **Dependency propagation.** Retiring or narrowing a baseline forces a recorded review of every
   ledger entry that cited it. Cheap now; the difference between a program and a pile later.
6. **Symmetry of evidence.** The evidentiary bar to retire equals the bar to promote. A baseline
   is not retired because it became unfashionable, and it is not kept because it is ours.

---

## 9. RULING ON K-032 — the Predictability Budget

**Position: it is sound bookkeeping wearing the costume of a success criterion. Adopt the
bookkeeping. Reject the criterion. Delete one clause.**

**What is sound, and I endorse without reservation.** Items 1 (the K-027 frontier), 2 (bits
currently captured), 4 (horizon vs scale), 5 (detector-invariance of every claim), and 7 (the
decision-threshold translation) are all measurements with stated uncertainties, and assembling
them in one public place is a genuinely scarce good. **Item 6 is the best idea in the entry and
possibly in the round:** converting every corpse from a null into an effect-size upper bound at
80% power. *"Tidal phase maps: |skill| < 0.05 bits/event at SoCal power"* is a scientific
statement; *"tidal triggering is null"* is a shrug. It costs one script over
results_exp_a/b/c/e/f/g.json and the existing harnesses, and it converts this program's largest
asset by volume — its clean, controlled negatives — into quotable quantities. **Run it in the
first week.** I also accept the reframe's motive: under "success = predicting an earthquake",
sixty years of the field have failed, the failures went unpublished, and each generation repeated
them. That *is* the field's central pathology and this program is small enough to escape it.

**What is a hedge, specifically.**

*First, and this is the load-bearing objection:* item 3 — "the model-free bits still available"
from K-029 — is presented as a **budget**, i.e. a ceiling. **It is not a ceiling.** No compressor
gives an upper bound on discoverable structure; it gives a lower bound on structure *it* found
(see the K-029 verdict). Building a success metric on a quantity that is mislabelled by direction
is exactly how a program grows an unfalsifiable core: every future failure is explained as *"we
spent the budget"*, and every future success is explained as *"the bound was encoding-limited"*.
Kepler writes that escape hatch into the entry himself — *"If someone later demonstrates skill
exceeding the bound, the budget is falsified and we learn our bound was encoding-limited, which
is itself the most interesting possible outcome."* A meta-prediction that reinterprets its own
refutation as a measurement update is not falsifiable in the sense that matters. **Delete that
clause.**

*Second:* "redefine success as producing the budget" makes the program succeed by definition,
because producing the budget is a matter of running scripts we already own. **A success criterion
you can guarantee by doing the work is not a criterion.** Kepler's framing is a good deal
psychologically — it converts a decades-long losing streak into an accumulating asset, and that
is not nothing for a program that has to keep going. But the psychological benefit is precisely
why it needs the harder version bolted on.

**The version I adopt.**

- **The Predictability Budget becomes the program's public reporting instrument and its ranking
  device for future hypotheses.** Every new entry states its expected contribution in budget
  units, and every corpse updates item 6. That is Kepler's real contribution here and it is
  valuable.
- **Item 3 is relabelled: "bits demonstrably still on the table (LOWER bound, encoding- and
  learner-limited)."** One word, and the escape hatch closes.
- **The definition of success stays falsifiable and stays capable of disappointing us:**
  *the program succeeds if it (a) adds validated, detector-invariant, out-of-sample bits over
  the harshest available baseline at a coarse-graining scale someone would act on, or (b) retires
  a claim the field currently believes, with a bounded effect size.* Both halves are things we
  can fail at. Producing a budget is neither.
- **The meta-prediction is restated in a form that can lose.** Not "future additions will sum to
  less than the residual budget" (unfalsifiable — a lower bound can always be revised upward),
  but: **"No future validated addition in this program will exceed X bits/event at scale S"**,
  with X and S fixed, published, and dated once K-027 and K-029 return. That version dies to a
  single result, which is the only property that matters.

Kepler: the budget is the right instrument and the wrong scoreboard. Keep it, publish it, rank by
it — and keep a definition of success that can still tell us we failed.

---

## 10. HOUSEKEEPING — scope and exposure lines for the existing baselines
*(required by the retirement mechanism, §8. These are scope statements, not new claims.)*

- **B-1** — SCOPE: generic temporal ETAS, ComCat M≥4.5, 1995–2026, 6 never-trained regions,
  temporal-only, vs local-oracle Poisson. **Quote as +0.66..+0.84 bits/event on the four
  adequately-powered holdouts** (AK 0.84, MX 0.66, PH 0.78, IRN 0.69, GRC 0.79); Caribbean +1.75
  is `underpowered: true` at n = 235 and must carry that flag. EXPOSURE: detector-invariance
  (untested — K-028); M0-dependence of α (untested — K-005); depth mixture (untested — K-004);
  magnitude strata above 4.5 (untested — K-001).
- **B-2** — SCOPE: SoCal SCSN, M≥2.5, walk-forward 2010–2018, n_test = 8,722, +1.87 bits/event
  vs test-rate oracle. The rise with magnitude rests on **n = 290** events in the M≥4 band and
  needs a sequence-block CI before it is quoted with one. EXPOSURE: detector-invariance (K-028);
  the nominal branching ratio n = 1.161 is supercritical while the 7-day n_eff = 0.60 — the
  parameterisation is window-dependent (K-018); post-2018 period untested.
- **B-3** — SCOPE: SoCal box, M≥5 within 7 d of an M≥5, test 18/30 = 0.60 vs 0.062 base,
  p = 7×10⁻¹⁵. Already caveated in §15 as **ordinary within-sequence triggering, not exotic
  meta-structure** — keep that sentence attached. EXPOSURE: detector-invariance (K-028); single
  region; n_test = 30.
- **B-4** — SCOPE: SoCal 0.2° cells, geodetic-vs-seismic moment ledger, train < 2010. The
  validated content is that **the negative space recovers known aseismic geology blind** (SAF
  creeping section, Imperial/Brawley; 175/229 silent cells interior, so not a border artifact)
  and that the top unexplained cells lie on the 1857 Fort Tejon strand. **AUTO-FLAGGED for
  challenge today under trigger 4(b):** 158/200 unexplained-silent cells have n_train < 20,
  39 have n_train = 0, 95 have n_test = 0. EXPOSURE: station-density confound (K-031, the open
  challenge); detector-invariance (K-028); persistence is already null (ρ = −0.13, n = 64), so
  χ at 0.2°/29 yr is transient-dominated and must not be described as a material constant.
- **B-5** — SCOPE: dilatation carries ±2× measurement uncertainty across velocity solutions;
  shear is robust. EXPOSURE: none new this round; note that K-024's rejection partly rests on
  the *recurrence* side of the moment budget, not on B-5's shear component.

---

*End Popper round 1. Counts: TESTABLE-NOW 21 (K-001, 002, 004, 005, 007, 008, 009, 010, 011,
013, 014, 015, 017, 018, 022, 026, 027, 028, 029, 030, 031); NEEDS-DATA 3 (K-003, 012, 023);
REFRAMED 7 (K-006R, 016R, 019R, 020R, 021R, 025R, 032); REJECTED 1 (K-024). Zero VALIDATED —
no tests have been run; nothing above may be claimed. Kepler: the estimator objections are not
objections to the ideas. Five of these are running-quality as written, and K-028 is a criterion
I should have proposed myself.*


---

# META-THEORIES (Wegener)

## Round 1 — 2026-08-09. Charter: unify the field's DOCUMENTED, replicated observations into
## single-entity accounts, and earn ledger status only by novel excess prediction.

My namesake did not discover a mechanism. He noticed that four literatures — coastline geometry,
fossil ranges, paleoclimate scars, stratigraphic correlation — were each individually explained
away, and that one claim explained all four at once. The claim was right for forty years before
anyone found the convection cell. What follows is that move applied here.

The single most important thing I want on the record before the table: **this program's five
corpses are not five failures. They are five measurements of the same quantity, taken in the one
place on Earth where that quantity is smallest.** That sentence is the whole of §W-RETRO below,
and everything between here and there is the evidence for it.

Web-verified this session via four parallel literature sweeps; every row carries its replication
status. Where a source I was handed was misattributed, I corrected it and say so. Where a number
I expected turned out wrong, I record the correct one and flag my prior.

---

## §W-OBS — THE OBSERVATION TABLE

Columns: **effect size** (as published) | **where it holds** | **where it does NOT hold / the
null domain** | **status** (REPLICATED / SINGLE-STUDY / CONTESTED / OURS).

The non-observations are load-bearing. O-7 (ordinary crustal seismicity feels no tide) constrains
harder than O-4 (tremor feels every tide), because O-4 has many possible causes and O-7 has few.

### A. The sensitivity ladder — the same load, wildly different responses

**O-1 — Remote dynamic triggering, Landers 1992.** Sudden seismicity-rate increases to ~1250 km
(Long Valley, Coso, The Geysers, Yellowstone, Cedar City, western Nevada), coincident with
surface-wave passage; implied dynamic stresses ~0.01–0.1 MPa, comparable to solid-earth tides.
*Holds:* pre-existing persistent seismicity, disproportionately geothermal / volcanic /
hydrothermal / extensional. *Does not hold:* remote triggering of LARGE events is essentially
absent (Parsons & Velasco 2011, Nat. Geosci., ngeo1110) — the far field triggers small events
only. *Status:* REPLICATED. Hill et al. 1993, Science 260, 1617, doi:10.1126/science.260.5114.1617.

**O-2 — Dynamic triggering is globally common.** 12 of 15 M>7 mainshocks (1990+) produced
detectable triggered seismicity somewhere during Love/Rayleigh passage. *Status:* REPLICATED as a
phenomenon; the stronger reading ("setting-independent") is CONTESTED — the largest responses
remain concentrated in fluid-rich terranes. Velasco et al. 2008, Nat. Geosci. 1, 375, doi:10.1038/ngeo204.

**O-3 — There is no triggering THRESHOLD; there is a continuous gain.** Statistical detection of
triggering down to dynamic strains of ~3×10⁻⁹ (sub-Pa-scale stress), with triggering intensity
scaling continuously with amplitude rather than switching on at a cutoff. *Status:* SINGLE-STUDY
but methodologically strong and never refuted. van der Elst & Brodsky 2010, JGR 115, B07311,
doi:10.1029/2009JB006681. **This row is the hinge of the entire table** — it says the crust's
response function is smooth and the observations differ by GAIN, not by an on/off switch.

**O-4 — Tremor and LFEs are modulated by sub-kPa to few-kPa tidal stress.** Cascadia, Nankai, and
the deep San Andreas; the Parkfield LFE result is interpreted as requiring near-lithostatic pore
pressure, i.e. effective normal stress driven toward zero. *Does not hold:* nowhere in the
ordinary brittle seismogenic zone. *Status:* REPLICATED across ≥4 independent settings and groups.
Rubinstein et al. 2008; Thomas, Nadeau & Bürgmann 2009, Nature 462, 1048, doi:10.1038/nature08654;
Royer et al. 2015; Yabe et al. 2015.

**O-5 — Tremor's tidal sensitivity RISES within a single ETS episode** (low early, high late, at
the ~2.5 kPa scale), interpreted as progressive fault weakening as slip accumulates. *Status:*
SINGLE-STUDY / emerging, consistent with Yabe et al. 2015. Houston 2015, Nat. Geosci. 8, 409,
doi:10.1038/ngeo2419. **A state variable observed changing during the event — the closest thing
in the whole table to a directly watched "living" state.**

**O-6 — Tremor triggered by teleseismic surface waves at 12–43 kPa** (Denali waves at Vancouver
Island), with later reports to ~3–8 kPa. *Does not hold:* the same patch does not respond every
time comparable waves pass — susceptibility itself varies in time. *Status:* REPLICATED
(phenomenon), CONTESTED (any universal threshold). Rubinstein et al. 2007, Nature 448, 579,
doi:10.1038/nature06017; Peng & Gomberg 2010, Nat. Geosci. 3, 599, doi:10.1038/ngeo940.

**O-7 — ORDINARY CRUSTAL SEISMICITY SHOWS NO TIDAL CORRELATION.** The paradigm null: shallow
California crustal seismicity, no significant correlation with solid-earth tides. *Status:*
REPLICATED — and independently reproduced by this program four times (EXP-A/B/C/E, §14–14b).
Vidale, Agnew, Johnston & Oppenheimer 1998, JGR 103, 24567, doi:10.1029/98JB00594.

**O-8 — …EXCEPT shallow thrust faults in high-ocean-load regions,** where the rate varies by a
factor of ~3 over the tidal cycle, best fit at μ≈0.4. Also reported for mid-ocean-ridge seismicity
and, globally, more strongly for reverse than strike-slip/normal mechanisms. *Status:* REPLICATED
but strictly setting-specific. Cochran, Vidale & Tanaka 2004, Science 306, 1164,
doi:10.1126/science.1103961.

**O-9 — THE POWER NUMBER (theory, and the most consequential row in this table for us).**
Rate-and-state predicts tidal correlation in ordinary crust is real but ~1% of events, requiring
catalogues of **10⁵–10⁶ events** for robust detection. *Status:* accepted theory, widely used as
the field's explanation for why O-7 looks null. Beeler & Lockner 2003, JGR 108, 2391,
doi:10.1029/2001JB001518. **Our EXP-A tested this at n ≈ 10²–10³ per bin. See §W-RETRO.**

**O-10 — Tidal stress lowers b, so large events prefer tidal maxima.** *Status:* SINGLE-STUDY, no
independent replication located; weight it low. Ide, Yabe & Tanaka 2016, Nat. Geosci. 9, 834,
doi:10.1038/ngeo2796. (Related: Tanaka 2010/2012 precursory tidal correlation before Tohoku and
Sumatra — small-N, retrospective, never prospectively confirmed. CONTESTED.)

### B. Fluids: the same load again, with pore pressure as the dial

**O-11 — Induced seismicity responds to injection RATE, not cumulative volume.** Wells above
~300,000 bbl/month are far more strongly associated with seismicity. Triggering pore-pressure
perturbations of order **0.01–0.1 MPa** on critically stressed faults. *Does not hold:* the large
majority of ~40,000 US Class II disposal wells induce nothing felt; Koyna reservoir impoundment
required >300–600 kPa. *Status:* REPLICATED (rate control, USGS + multiple regional studies);
CONTESTED that "10 kPa" is any kind of universal constant — it is a fault-state-dependent number
spanning at least a decade and a half in magnitude. Weingarten et al. 2015, Science 348, 1336,
doi:10.1126/science.aab1345; Ellsworth 2013, Science 341, doi:10.1126/science.1225942.

**O-12 — Basel 2006: the largest event (ML 3.4) came AFTER shut-in.** *Status:* REPLICATED as a
single deeply-analysed case (many re-analyses, one site). The delay is a memory the instantaneous
Coulomb picture does not have.

**O-13 — Near-lithostatic pore pressure in the tremor/SSE zone,** inferred from Poisson's ratios
to ~0.35 and bright plate-interface reflectors. *Status:* the OBSERVATION (high Vp/Vs and high
reflectivity co-located with tremor across Cascadia, Nankai, Shikoku) is REPLICATED; the
QUANTITATIVE effective-stress inference is CONTESTED — crack fabric and serpentinization also
raise Vp/Vs. Audet et al. 2009, Nature 457, 76, doi:10.1038/nature07650; Kodaira et al. 2004,
Science 304, 1295.

**O-14 — LFE stress drops are 2–3 orders below ordinary earthquakes** (kPa to tens of kPa vs
1–10 MPa). *Status:* REPLICATED qualitatively, study-dependent numerically. Chestler & Creager
2017; Ide et al. 2008. **An independent route to the same conclusion as O-13, with different
assumptions — which is why I weight the pair heavily even though each alone is contested.**

### C. Slow loads and slow slip

**O-15 — Seasonal hydrologic loading modulates seismicity at kPa amplitudes.** California: annual
hydrospheric Coulomb changes ~0.5–2 kPa, detectable modulation of small-earthquake timing, with
hydrology the largest seasonal stress source (larger than atmospheric, thermal, or tidal).
Himalaya: winter/summer rate ratio ~2× at 2–4 kPa monsoon loading. Japan: snow-load unclamping at
a few kPa. *Does not hold / caution:* our own EXP-F did NOT detect an annual signal in SoCal, but
that null is weak (the 7-day method check failed to fire). *Status:* REPLICATED. Johnson, Fu &
Bürgmann 2017, Science 356, 1161, doi:10.1126/science.aak9547; Bettinelli et al. 2008, **EPSL**
(commonly miscited as Nat. Geosci. — corrected here); Heki 2003, **EPSL** (commonly miscited as
Science — corrected here); Amos et al. 2014, Nature 509, doi:10.1038/nature13275.

**O-16 — THE FREQUENCY PARADOX (my framing of O-7 + O-15, and I have not seen it stated).**
Ordinary California crust demonstrably responds to a **1–2 kPa annual** load and demonstrably does
not respond to a **1–3 kPa semidiurnal** load. Same crust. Same amplitude. Different period. Any
theory of the table must explain this, and only one class of theory does.

**O-17 — Episodic tremor and slip recurs quasi-regularly** (Cascadia ~14 months). *Status:*
REPLICATED. Rogers & Dragert 2003, Science 300, 1942. *Does not hold as a predictor:* no great
earthquake has occurred in Cascadia or Nankai to test its precursory value.

**O-18 — Precursory slow slip / migrating foreshocks before Tohoku 2011 (2–10 km/day) and
Iquique 2014.** *Status:* a genuine TWO-CASE pattern, not a law. *The null side is decisive:* no
established base rate; the overwhelming majority of slow-slip episodes and swarms are followed by
nothing; and the apparent prevalence of foreshocks is strongly dependent on the completeness
threshold used. Kato et al. 2012, Science 335, 705, doi:10.1126/science.1215141; Ruiz et al. 2014,
Science 345, 1165, doi:10.1126/science.1256074; Mignan 2014, Sci. Rep. (foreshock meta-analysis) —
CONTESTED.

**O-19 — Repeating earthquakes work as creep meters:** +18% seismic moment per 10× recurrence
interval, 160 Parkfield sequences. *Does not hold:* the magnitude-dependent stress-drop scaling
underlying the moment→slip conversion has NOT been reproduced by later work; repeaters are
interpretable only where small locked patches are embedded in creep. *Status:* REPLICATED
(scaling), DOCUMENTED FAILURE (its mechanical justification). Nadeau & McEvilly 1999, Science 285,
718; Uchida & Bürgmann 2019, Annu. Rev. Earth Planet. Sci. 47, 305.

### D. The constitutive law and its limits

**O-20 — Rate-and-state friction.** Seismicity-rate response to a stress step
`R = r₀·exp(Δτ/Aσ)`, relaxing over `t_a = Aσ/τ̇`. Lab a ≈ 0.005–0.015, b ≈ 0.005–0.02; field Aσ
inverted from aftershock sequences ≈ **0.01–0.1 MPa**. *Does not hold:* the lab-to-field nucleation
dimension gap (cm-scale lab Lc vs field-inferred kilometre-scale) is an acknowledged UNRESOLVED
failure, not a caveat; slow slip requires a narrowly tuned a≈b regime that is fitted, not
predicted. *Status:* REPLICATED as a formalism, with a named hole. Dieterich 1994, JGR 99, 2601;
Marone 1998, AREPS; McLaskey 2019, JGR 124, 12882, doi:10.1029/2019JB018363.

**O-21 — Nucleation: visible in the lab, elusive in the field.** Lab shows smooth accelerating
preslip to a critical length. Field: the seismic nucleation phase (~0.5% of moment, ~1/6 of
duration, 48 events Mw 1.1–8.1) has never been established as a general precursor, and the
preslip-vs-cascade question is unresolved. The Parkfield Prediction Experiment — the best-
instrumented prediction test ever fielded — recorded **no diagnostic precursory signal** when the
M6 finally arrived in 2004. *Status:* CONTESTED (mechanism), REPLICATED NEGATIVE (field
detectability). Ellsworth & Beroza 1995, Science 268, 851; Bouchon et al. 2013, Nat. Geosci. 6, 299.

**O-22 — b decreases with differential stress.** `b = 1.23 ± 0.06 − (0.0012 ± 0.0003)·(σ1−σ3)` in
MPa, i.e. **db/dσ ≈ −0.0012 per MPa** — an order of magnitude SMALLER than the figure I carried
into this session, and I record the correction rather than the prior. Field ordering
thrust < strike-slip < normal. *Does not hold:* b̂ is exquisitely sensitive to Mc; the operational
precursor version (foreshock traffic-light) is CONTESTED and not prospectively validated.
*Status:* REPLICATED (relation), CONTESTED (application). Scholz 2015, GRL 42, 1399,
doi:10.1002/2014GL062863; Schorlemmer, Wiemer & Wyss 2005, Nature 437, 539, doi:10.1038/nature04094;
Gulia & Wiemer 2019, Nature 574, 193 (contested).

**O-23 — Coulomb static triggering at ~+0.01 MPa.** *The triggering side is REPLICATED; the
suppression side (stress shadows) is genuinely CONTESTED* — Harris & Simpson 1998 for, Mallman &
Zoback 2007 against (rates rose almost everywhere regardless of predicted sign), and the
static-vs-dynamic diagnostic itself is disputed (Felzer & Brodsky 2006 vs Marsan 2010). I mark it
so, as the charter requires; it constrains weakly. King, Stein & Lin 1994, BSSA 84, 935;
Stein 1999, Nature 402, 605.

### E. Universality, modes, and the predictability record

**O-24 — Omori/ETAS universality.** p ≈ 1.1 ± 0.25 globally; ETAS as the operational standard.
*Status:* REPLICATED. Utsu, Ogata & Matsu'ura 1995; Ogata 1988, JASA 83, 9. **OURS (B-1/EXP-M):
p ∈ [0.94, 1.08] in 12 of 13 regions; generic shape + locally-fitted μ delivered +0.66 to +1.75
bits/event on six never-trained regions, within 0.07–0.15 bits of each region's post-hoc in-sample
ceiling. Fault-TYPE pooling FAILED the frozen sign test (2/6); α and K spread WITHIN type ≥
BETWEEN types; the only thing type pooling changed was μ, which should never be pooled.**
Related: Båth's law (Δm≈1.2) is empirically REPLICATED but CONTESTED as physics — derivable from
GR plus selection (Helmstetter et al. 2003).

**O-25 — Swarms vs mainshock-aftershock is a continuum, not a binary.** Of 71 SoCal bursts:
14 clean aftershock sequences, 18 swarm-like, ~39 **blended**. Mechanism is region-dependent —
fluid-driven in West Bohemia/Vogtland, aseismic-creep-driven on oceanic transforms (migration
0.1–1 km/hr, too fast for the diffusivities inferred from injection). *Status:* REPLICATED that
the taxonomy blurs; CONTESTED that there is one swarm mechanism. Vidale & Shearer 2006, JGR 111,
B05312, doi:10.1029/2005JB004034; Roland & McGuire 2009, GJI 178, 1677; Hainzl (Vogtland).

**O-26 — No precursor has ever passed a prospective CSEP-style test.** Smoothed seismicity and
ETAS remain the hardest baselines to beat; operational earthquake forecasting is explicitly
probabilistic, not predictive. *Status:* REPLICATED NEGATIVE. Geller et al. 1997, Science 275,
1616; Jordan et al. 2011, Ann. Geophys. 54, 315. **OURS: this program is a faithful five-corpse
instance of exactly this result.**

### F. Our own documented realities (advantage: I know how they were made)

**O-27 (OURS, B-2)** — SoCal walk-forward ETAS, params frozen from train, +1.907 bits/event vs
train-rate Poisson and +1.866 vs the **test-rate oracle**; skill RISES with magnitude (+2.54 to
+2.58 at M≥4); triggered fraction 0.915. `results_exp_h.json`.

**O-28 (OURS, B-3)** — P(M≥5 within 7 d | M≥5) = 18/30 = 0.60 vs 0.062 base, p = 7×10⁻¹⁵.
`results_exp_i.json`.

**O-29 (OURS, B-4)** — The stress ledger's negative space is real geology: silent-loading cells
recover, blind, BOTH the SAF **creeping** section (Imperial/Brawley too) AND the 1857 Fort Tejon
**locked** Big Bend strand. `results_exp_j.json`, `results_exp_k.json`. **The two most opposite
hazard states in California came out of the same list.** Hold that thought for W-006.

**O-30 (OURS)** — Coso Fig-4c north bin (36.2–36.6 N, −118 to −117.6), SCSN M≥1.5, shear on the
per-event FM plane: Pm/P0 = **0.340**, one-sided p = 0.041, n = 113 — QUALIFIED REPRODUCTION; and
σ_n in the same bin is **null** (0.271, p = 0.12). §12, `results_coso_fig4c.json`.

**O-31 (OURS, corpses)** — Static tidal-phase susceptibility maps: null with clean anti-leak
controls. Feature-vs-amplitude correlations: small-n bias. Phase migration: killed by our own
cross-catalog and walk-forward confirmations. Fixed periodicities incl. annual: not detected
(weak power). Golden-ratio interevent structure: nothing. Fault-type pooling: failed.
Spatial transfer of sequence shapes (p, b): failed 7/15. §14–17.

---

## §W-SPINE — WHAT THE TABLE IS ACTUALLY A TABLE OF

Read O-3, O-9, O-16 and O-20 together and the table stops being a list of phenomena.

O-3 says there is no threshold — the response is smooth all the way down to nanostrain. O-20
gives the response its functional form, `R = exp(Δτ/Aσ)`, relaxing over `t_a = Aσ/τ̇`. O-9 says
that in ordinary crust the tidal response is ~1%, not zero. O-16 says the same crust responds to
an annual load of the same amplitude.

So the field has not been measuring twelve phenomena. It has been measuring **one transfer
function, sampled at different values of two numbers**: a **gain** (1/Aσ, set by effective normal
stress, i.e. by fluids) and a **corner frequency** (1/t_a, set by the loading rate). Every row in
§W-OBS is a point on that surface:

| regime | Aσ implied | t_a implied | what the literature calls it |
|---|---|---|---|
| deep LFE/tremor, near-lithostatic Pp | ~kPa | short | "extreme tidal sensitivity" (O-4, O-5, O-13, O-14) |
| geothermal / high-Pp shallow / swarm | ~10s kPa | short–mid | "remote triggering hotspots" (O-1), "swarms" (O-25), **our Coso 0.340** (O-30) |
| critically stressed faults near wells | 10–100 kPa | short | "induced seismicity" (O-11) |
| ordinary seismogenic crust | 0.01–0.1 MPa | months–years | "the tidal null" (O-7) + "seasonal modulation works" (O-15) |
| deep locked megathrust | high | very long | "no short-period sensitivity, decadal cycles" |

That is the unification claim in one paragraph. The six W-entries below are the serious candidate
readings of it — including the two that say it is wrong.

---

## §W-THEORIES

Each entry: **mechanism** → **corollary map** (which O-rows it absorbs) → **≥2 distinguishing
predictions with falsifiers** → **what the living entity IS** → **our assets that bear on it**.
Per the Einstein criterion, an entry that merely renames the table is worthless; each prediction
below is one the disjoint explanations do not make.

---

### W-001 — THE GAIN FIELD. One constitutive law; one hidden two-component state field; the crust is a spatially varying LOW-PASS FILTER on stress.

**Mechanism.** There is one crust and one friction law. Seismicity rate responds to any stress
perturbation as `R = exp(Δτ/Aσ)`, relaxed over `t_a = Aσ/τ̇`. `Aσ = a·(σ_n − P_p)` is a hidden
scalar field; `τ̇` is the local loading rate. Nothing else is needed. Every documented
"sensitivity" phenomenon is a reading of `1/Aσ`; every documented FREQUENCY dependence is a
reading of `t_a`. Fluids do not make a different kind of fault; they move the same fault along the
gain axis.

**Corollary map.** Absorbs O-1 (geothermal = low Aσ), O-3 (no threshold — exp is smooth), O-4/O-13/
O-14 (near-lithostatic Pp → Aσ→0 → O(1) response to kPa), O-6 (variable susceptibility = variable
Aσ in time), O-7 + O-8 + O-9 + O-15 + O-16 **together** (see P1), O-11 (Aσ small on critically
stressed faults; most wells sit on faults with large Aσ, hence the vast null), O-20, O-22 (b and Aσ
are both functions of effective stress, hence correlated without either causing the other), and —
crucially — **O-24**: Dieterich's derivation makes the Omori exponent p≈1 *parameter-free* while
the aftershock DURATION `t_a = Aσ/τ̇` carries all the local information. That is precisely our
EXP-M shape/scalar split, derived rather than observed.

**P1 (the flagship excess prediction) — THE FREQUENCY-RESPONSE COLLAPSE.** A rate-state solid is a
low-pass filter with corner at `1/t_a`. Therefore, in one and the same cell, the response to a
periodic load of period T is attenuated by a factor that depends only on `T/t_a`. Since `t_a` is
independently measurable from that cell's own Omori aftershock-decay duration — **with no tides
involved** — the theory predicts, with **zero free parameters**, the ratio of a cell's response at
two different tidal frequencies. Concretely: the **fortnightly (Mf, 14.8 d) to semidiurnal (M2,
12.42 h) response ratio must be an increasing function of the cell's measured `t_a`**, following
the rate-state low-pass form.
*Why no disjoint explanation makes it:* the tidal literature treats amplitude, the aftershock
literature treats decay, and nobody has ever required the two to agree in the same rock. It also
retrodicts O-16 immediately — annual loads are quasi-static (T ≫ t_a: full `exp(Δτ/Aσ)` response,
which is why Johnson/Bettinelli/Heki succeed) while semidiurnal tides are far above the corner
(strongly attenuated, which is why Vidale nulls). Same crust, same kPa, different side of the
corner.
*Falsifier:* the Mf/M2 ratio is uncorrelated with independently measured `t_a` across cells at
adequate power; or the crust responds equally at both frequencies; or the fitted corner sits
orders of magnitude away from `1/t_a`.

**P2 — TIDAL SENSITIVITY IS PREDICTABLE FROM NON-TIDAL DATA, AND IS NOT SYNONYMOUS WITH
"GEOTHERMAL".** Rank cells by `Aσ` estimated purely from aftershock relaxation and ledger loading
rate, blind to tides and blind to the geothermal map. W-001 predicts the top-ranked cells recover
the known tidally-sensitive ones (Coso-north) **and additionally name non-geothermal cells that
are tidally sensitive.**
*Falsifier:* when scored blind, every tidally-sensitive cell is geothermal/volcanic, i.e. the
`Aσ` ranking adds nothing over a map of hot water. That would demote Aσ from a field to a label.

**P3 — THE INDUCED-SEISMICITY BRIDGE.** The `Aσ` inverted from aftershock decay in a region must
predict the region's induced-seismicity susceptibility — the pore-pressure change per unit rate
increase — quantitatively, across Oklahoma / Basel / Groningen / Koyna. Koyna's 300–600 kPa
requirement (O-11) must correspond to a measurably larger `Aσ` than Oklahoma's 10–100 kPa.
*Falsifier:* no rank agreement between aftershock-derived Aσ and induced-seismicity threshold
across sites with adequate data.

**The living entity.** The crust is a **continuous medium with a spatially and temporally varying
gain and memory** — not a set of faults that break, but a field that amplifies. It is "alive" in
exactly the sense that an amplifier with a slowly drifting bias is alive: the same input produces
wildly different output depending on a hidden state, and the hidden state is knowable.

**Our assets.** `data/xue_lu_zenodo/SCSN_original_catalog.txt` (633,667 events, 1981–2018,
Mc ≤ 1.7 all eras) + `data/xue_lu_zenodo/Tidal_{N_0,N_90,S_0,Vol}.txt` (6000-s series from 1981) +
`coso_fm_test.py` / `coso_fig4c_test.py` (validated FM-resolved machinery) + CFM5.3 (557 segments
with strike/dip/rake) + frozen EXP-H ETAS (`results_exp_h.json`) as the baseline intensity +
`results_exp_j/k.json` for τ̇. **P1 is 100% on disk.** Directly extends K-036 (which fits Aσ but
has no frequency axis and no zero-parameter prediction) and supplies K-035's power audit with the
physical amplitude it should be auditing against.

---

### W-002 — THE VALVE. The state is not a smooth scalar but the CONNECTIVITY of an overpressured fluid phase; the crust is plumbing, and plumbing has hysteresis.

**Mechanism.** Fault zones are compartmented, self-sealing, episodically breaching hydraulic
systems (Sibson's fault-valve). Effective stress is not a smooth field but a two-phase mixture:
sealed compartments at hydrostatic-ish pressure, and percolating overpressured pathways. What
looks like "low Aσ" is really "currently connected to an overpressured reservoir." Permeability is
dynamic and is itself modified by shaking and by slip.

**Corollary map.** Absorbs everything W-001 does, plus the rows W-001 handles awkwardly:
O-5 (sensitivity RISING during an ETS episode — a valve opening in real time, not a fixed gain),
O-12 (Basel's largest event AFTER shut-in — pressure diffusion has its own clock), O-25 (swarms as
a distinct MODE with migration fronts, and the region-dependence of swarm mechanism), the delayed
and long-lived character of remote triggering in geothermal fields (O-1), and O-19's failure mode
(repeaters stop being creep meters when the surrounding patch's hydraulics change).

**P1 — BIMODALITY, NOT A CONTINUUM.** W-001 says sensitivity is a smooth monotone function of a
scalar, so across cells the sensitivity distribution should be **unimodal and log-continuous**.
W-002 says percolation is a transition, so the distribution should be **bimodal** (connected vs
sealed) with a sparse middle, and the two modes should separate on migration signature: connected
cells show √t diffusive migration of activity, sealed cells do not.
*Falsifier:* unimodal sensitivity distribution with no migration-signature separation, at power
sufficient to see a 2-component mixture.

**P2 — HYSTERESIS (the memory a gain field does not have).** After a cell is dynamically triggered
by a large remote event, its sensitivity to subsequent small loads should remain **elevated for
months to years** and then decay with a recharge/resealing timescale — a path dependence. W-001
predicts sensitivity returns immediately to its Aσ-set value once the transient stress has passed.
This is directly testable in our own record: Landers 1992, Hector Mine 1999, Ridgecrest 2019 are
all in the catalogue on disk, with 1981–2018 (SCSN) and 2010–2026 (ComCat) coverage around them.
*Falsifier:* post-trigger tidal/seasonal coupling in triggered cells is statistically
indistinguishable from pre-trigger coupling, at matched event counts.

**P3 — SEALING SHOULD BE VISIBLE AS A DRIFT IN t_a, NOT ONLY IN AMPLITUDE.** If valving is real,
the corner frequency of a cell drifts with its hydraulic state, so W-001's P1 collapse should show
**systematic residuals organised in time** (cells moving along the t_a axis after fluid events),
not just scatter.
*Falsifier:* P1 residuals are white in time.

**The living entity.** The crust is a **circulatory system** — a self-sealing, self-breaching
network of overpressured pathways, in which "state" means which valves are currently open.
Earthquakes are the noise the plumbing makes. This is the reading closest to Jim's "living entity"
in the literal sense: it has organs (compartments), a fluid, a pulse (recharge), and memory.

**Our assets.** Same catalogue and tidal series as W-001; the three great SoCal sequences for the
hysteresis test; `results_exp_k.json` swarm/style stratification; and the `tsi_map.GEOTHERMAL`
layer as an independent labelling of known plumbing. P2 is the cheapest genuinely NEW measurement
in this whole entry set, and it is the one that separates W-002 from W-001 most sharply.

---

### W-003 — THE NULL UNIFIER (and the current champion): THERE IS NO ENTITY. One stationary branching process with static heterogeneous parameters; every apparent "state" is either the process's own memory or the observer's detection function.

**Mechanism.** The crust is a spatially heterogeneous but **temporally static** medium hosting a
near-critical branching process. Clustering alone manufactures rising rates, rising variance,
rising correlation length, shrinking interevent times, and apparent precursors. Everything anyone
has ever called a precursor is either (i) the branching process's own conditional intensity, or
(ii) a change in Mc, or (iii) small-n estimator bias.

**Corollary map.** Absorbs O-24 (universality is the signature of a single simple generator),
O-26 (nothing beats ETAS because nothing else exists to beat it with), O-21 (no nucleation because
there is no nucleation), O-23's contested half (stress shadows are contested because they are not
there), O-18's null side (foreshocks are Mc-dependent selection), O-10's fragility, **and our
entire corpse list including the EXP-C2 self-correction, the EXP-B n-bias finding, and the EXP-M
fault-type failure** — under W-003, "fault type" cannot carry information because parameters are
heterogeneous at a scale finer than any label.

**P1 — ZERO INCREMENTAL BITS, FROM EVERYTHING.** W-003 forbids, categorically, that any order
parameter (K-018 n(t), K-019R residual correlation time, K-021R fixed-n Var(b), K-022 S_max,
K-026 entropy) or any load covariate or interaction (K-033/K-036/K-039) yields incremental
bits/event above frozen ETAS on a temporal holdout with an ETAS-sim null.
*Falsifier:* any single one of them clears zero with a block-bootstrap CI excluding zero. **This
is the most falsifiable entry in this section, and it is why it belongs here rather than being
dismissed as defeatism.** It is also currently WINNING: our record is 5 corpses to 0.

**P2 — THE RESIDUALS ARE WHITE.** K-009 must return spatial correlation length ≈ 0 and temporal
correlation time ≈ 0 in the ETAS residual field, against the ETAS-sim null.
*Falsifier:* a leading residual EOF with a red spectrum and a correlation time of months.
(Note the internal consistency: W-001/W-002 both predict K-009 residuals are RED, with a
correlation time equal to the local `t_a` — months to years. K-009 is therefore a *shared*
discriminator, which is why Popper ranked it "run it first" and why I agree.)

**The living entity.** There isn't one. The apparent life is our pattern-matching plus the
process's own memory. **This entry is in the ledger because a unifier who cannot state the version
of the world where his unification is unnecessary has not done the work.**

**Our assets.** Everything. W-003 is the standing null for every other entry here; it is scored
automatically whenever any of them is tested.

---

### W-004 — THE INSTRUMENT. The sensitivity map is a detectability map: the observation table is partly a survey of where we can count small earthquakes.

**Mechanism.** Detecting a fractional rate modulation `a` requires roughly `N ≳ a⁻²` events. Every
setting in which small-stress sensitivity has been "discovered" is also a setting with anomalously
low Mc: template-matched tremor/LFE catalogues in the best-instrumented subduction zones (O-4),
dense local arrays around geothermal fields (O-1) and injection wells (O-11), and the enormous
event counts those produce. Every setting where the null lives is a setting with ordinary Mc and
ordinary counts (O-7). The gradient in the table may be a gradient in `N`, not in `Aσ`.

**Corollary map.** Absorbs O-7 vs O-4 without any physics; absorbs O-9 directly (Beeler & Lockner
literally say you need 10⁵–10⁶ events); absorbs O-18's Mc-dependence of foreshock prevalence
(Mignan 2014); absorbs O-22's Mc sensitivity; and absorbs the ordering of our OWN results — the
one place we reproduced a tidal signal (Coso-north, O-30) is the place with the local FM catalogue
and the dense array, and our nulls came from FM-matched subsets of n ≈ 10²–10³. It is the
uncomfortable cousin of K-028 (detector invariance) and the reason K-035 exists.

**P1 — POWER-MATCHED EQUIVALENCE.** Subsample the high-count settings down to the event counts of
the null settings, and the "sensitivity gradient" must **collapse**: Coso-north at n = 113 vs
ordinary SoCal cells at n = 113 should be statistically indistinguishable in detected modulation,
and the apparent amplitude should scale as `N^(−1/2)` in both.
*Falsifier:* at strictly matched N and matched estimator, Coso-north still shows a large coupling
and ordinary crust still shows zero. **That single comparison is worth more than any new
hypothesis, because it decides whether §W-OBS is a table about the Earth or about seismometers.**

**P2 — THE INVERSE PREDICTION: SoCal's null must DIE at full catalogue power.** W-004 (and O-9)
jointly predict that running the tidal test on the **whole** SCSN catalogue above a fixed
era-stable Mc (order 6×10⁵ events, not the 10²–10³ our corpses used) will reveal the ~1% ordinary-
crust modulation that theory says has been there all along.
*Falsifier:* a genuinely powered full-catalogue test still returns zero, with an upper bound below
1%. That outcome would be a real result: it would falsify Beeler–Lockner in California and make
O-7 a physical fact rather than a power limit.

**The living entity.** Partly an artifact of where we put the instruments. The honest version of
Jim's question is: *is the crust alive, or is it only alive where we are watching closely?*

**Our assets.** `SCSN_original_catalog.txt` is 633,667 events with Mc ≤ 1.7 in every era — this
program is one of the few that can actually run P2. `binscan_QTM.csv` / `binscan_SCSN.csv` and the
QTM template-matched catalogue give a second, independent detector for the K-028 invariance arm.

---

### W-005 — SHAPE VS SCALAR: universality is not a discovery, it is a receipt for the absence of information. All predictive content lives in the local scalars.

**Mechanism.** This is a meta-theory about the observation table rather than about rock, and it is
the most Wegener-like entry here: the coastlines fit, and that fact is about the map, not the
shore. Every documented regularity in §W-OBS decomposes cleanly into a **universal shape** and a
**local scalar**: Omori p (universal) vs aftershock duration t_a (local); GR b≈1 (universal) vs
corner magnitude (local); rate-state a,b (universal) vs Aσ (local); ETAS kernel form (universal)
vs μ (local). **Universality is precisely insensitivity to local details — therefore any quantity
that is universal carries zero information about local state, and therefore zero forecasting
skill.** Sixty years of the field's best work has been the discovery of shapes; the scalars have
been treated as nuisance parameters.

**Corollary map.** Absorbs O-24 and, in one stroke, our EXP-M result in its exact observed form:
p universal in 12/13 regions, generic-shape-plus-local-μ transferring at +0.66..+1.75 bits, and
fault-TYPE pooling FAILING because "type" is a shape label being asked to do a scalar's job — and
because the only thing type pooling actually changed was μ, the one quantity that must never be
pooled. Also absorbs O-19 (repeaters work as creep meters because they are a direct scalar
readout, and fail when their shape assumption is questioned), O-26 (smoothed seismicity wins
because it is nothing BUT a local scalar field), and O-25 (the swarm/mainshock taxonomy blurs
because it is a shape label over a continuous scalar).

**P1 — THE ARCHITECTURE PREDICTION.** Out-of-sample forecast improvement will come **only** from
better estimation of local scalars, never from richer functional forms. Concretely: an ETAS whose
μ, t_a (=Aσ/τ̇) and corner magnitude are locally estimated will beat any structurally fancier model
(connectivity-magnitude coupling, K-022; multi-kernel; neural) that uses globally pooled scalars,
on the same holdout.
*Falsifier:* K-022's percolation model, or any structural elaboration, beats scalar-refined ETAS
on CRPS out of sample. That is a clean head-to-head and both models already have designs in this
ledger.

**P2 — THE UNIVERSALITY/SKILL ANTICORRELATION.** Across the 13 world regions, rank each ETAS
parameter by its cross-region dispersion. W-005 predicts **skill contribution is monotone in
dispersion**: the more universal a parameter, the less it contributes to bits/event when localised.
p (dispersion ~0.07) should contribute essentially nothing when locally fitted; μ (dispersion of
orders of magnitude) should contribute nearly everything; α and K in between.
*Falsifier:* localising p buys as much as localising μ. **This is nearly free — it is a re-scoring
of `results_exp_m.json` with one parameter localised at a time.**

**The living entity.** The living thing is **the scalar field, not the law**. The law is dead
physics, identical everywhere, and that is exactly why it can be written down. Jim's intuition
that "the crust is a stateful, driven, responsive system" is, in this reading, precisely correct
and precisely orthogonal to everything the field has been publishing.

**Our assets.** `results_exp_m.json` (13 regions, all per-region fits and holdout scores already
stored) makes P2 a re-scoring job, not an experiment. P1 is the K-022 vs scalar-ETAS head-to-head.

---

### W-006 — THE SHADOW. Seismicity is the fast, thresholded readout of a slow ASEISMIC slip field. The living entity is the creep; earthquakes are what it says out loud.

**Mechanism.** Split the table by what is aseismic. Everything documented that behaves in an
orderly, quasi-predictable way is **aseismic or aseismically driven**: ETS recurs on a ~14-month
clock (O-17); repeaters meter creep rate (O-19); induced seismicity tracks injection rate (O-11);
seasonal modulation tracks hydrology (O-15); the two megathrust precursors were slow-slip
transients (O-18); tremor's sensitivity rises with accumulated slow slip (O-5). Everything
documented that is disorderly is the **seismic** part (O-21, O-26). Therefore: the crust's state
variable is the aseismic slip-rate field, and the earthquake catalogue is a sparse, thresholded,
heavily-clustered sampling of it.

**Corollary map.** Absorbs O-5, O-11, O-15, O-17, O-18, O-19, O-25's creep-driven swarms; explains
O-26 (we have been trying to forecast the sampler from the sampler); and it absorbs **O-29, our
own B-4, in a way nothing else here does.** The stress ledger found the SAF creeping section and
the 1857 locked Big Bend strand in the SAME silent list — because a moment-deficit ledger built
from a catalogue measures *what did not radiate*, and cannot distinguish "absorbed aseismically"
from "stored elastically." Those are the two most opposite hazard states in California, and they
are degenerate in catalogue space by construction.

**P1 — THE DEGENERACY-BREAKING PREDICTION (sharp, cheap, and ours).** The silent-loading list is
a mixture of exactly two populations with **opposite** hazard, and the discriminator is **not in
the catalogue** — it must come from geodesy (surface creep rate / geodetic coupling). W-006
predicts: (a) no catalogue-derived statistic separates the creeping cells from the locked cells
at better than chance; (b) adding a single geodetic coupling covariate separates them cleanly.
*Falsifier:* a purely catalogue-derived statistic (b-value, swarm fraction, repeater fraction,
interevent CV) separates creeping from locked silent cells. **Note the honest tension: repeater
fraction is arguably a catalogue-derived aseismic readout, so if THAT separates them, W-006 is
confirmed in spirit and refuted in letter — I record the ambiguity rather than hide it.**

**P2 — THE SKILL-CEILING RANKING.** Forecast skill above the ETAS ceiling should be obtainable
**only** in proportion to the quality of the region's aseismic readout. Ranking, pre-registered:
regions with repeater catalogues + dense geodesy + tremor (Cascadia, Nankai, Parkfield) >
regions with geodesy only (SoCal) > regions with neither (most of our 13 boxes) — at matched
catalogue power.
*Falsifier:* a region with a rich aseismic readout shows no more above-ETAS skill than one with
none. This also predicts something uncomfortable and specific: **SoCal's above-ETAS ceiling is
low, and our program has been mining the wrong region.** The place to test the "living entity" is
Cascadia or Nankai, where the entity is directly observable, not SoCal, where only its shadow is.

**The living entity.** The aseismic slip field. It is continuous, it is driven, it is stateful, it
has memory and episodicity and a clock — everything Jim means by "living" — and the earthquake
catalogue is its exhaust. Under W-006, asking a catalogue to reveal the entity is asking a
tachometer to describe the engine.

**Our assets.** `results_exp_j.json` / `results_exp_k.json` (the silent list, 200 unexplained-
silent cells, 154 strike-slip / 90 on-fault); `data/socal_strain_grid.npz` (66×86 max-shear +
dilatation) for the geodetic covariate; B-5's finding that shear is robust and dilatation is not,
which tells us *which* geodetic component the discriminator may use. **P1 is fully on disk and is
a one-afternoon job.**

---

## §W-RETRO — WHAT THIS SAYS ABOUT OUR OWN NULLS (the prize the charter names)

The charter asked for a meta-theory that retroactively explains the PATTERN of our nulls and
reproductions. Here it is, in the order the pattern occurred.

**1. Why ordinary SoCal seismicity shows no tidal effect while tremor/LFEs feel every tide.**
Two reasons, and they compound. *Physics (W-001):* ordinary crust has Aσ ~ 0.01–0.1 MPa versus
near-zero in the LFE zone — a gain difference of two to three orders of magnitude, independently
corroborated by the stress-drop contrast (O-14) and the Vp/Vs inference (O-13). *Frequency
(W-001-P1):* semidiurnal tides sit far above the ordinary crust's corner frequency `1/t_a`, so
even that residual gain is attenuated — which is why the SAME crust at the SAME kPa amplitude
responds to the annual load (O-15) and not the tidal one (O-7). Our four tidal corpses measured a
quantity that theory says is ~1% (O-9), using samples of n ≈ 10²–10³, when O-9 says you need
10⁵–10⁶. **We did not measure zero. We measured "less than our resolution," and our resolution was
about three orders of magnitude coarser than the predicted effect.** That is not a failure of the
work — the work was clean and the controls were right — it is a failure to have computed the power
first, and it is exactly what K-035 was created to fix.

**2. Why Coso-north modulates.** Coso-north is the one cell in our study area that sits at the
geothermal end of the gain ladder: a fluid-rich, high-Pp, swarm-dominated volume — plausibly
Aσ an order of magnitude below ordinary crust. W-001 says its modulation should be roughly
`exp(Δτ/Aσ) − 1` with a small Aσ; we measured Pm/P0 = 0.340 (O-30). And the internal detail is the
part that convinces me: the effect appears in **shear on the FM plane** and is **null in σ_n**
(0.271, p = 0.12). A detection artifact does not know the difference between two stress components
resolved on the same events. That is a mechanism signature, not a p-value.

**3. Why fault-type pooling failed while ETAS shape transferred universally.** W-005, exactly:
"fault type" is a **shape** label, and shapes are universal, so the label can carry no information;
meanwhile the thing TYPE pooling actually moved was **μ**, a local scalar that must never be
pooled — which is why pooling by type made Alaska and Mexico *worse*. Our own worker diagnosed
this correctly in §17 ("the transferable object is the universal clustering law + a locally-
estimated μ"). W-005 says that sentence is not a finding about ETAS. It is a finding about the
whole field.

**4. Why the stress ledger found the safest and the most dangerous places in the same list.**
W-006: a catalogue-based moment-deficit ledger is blind to the difference between "absorbed
aseismically" and "stored elastically." The degeneracy is structural, not a bug, and it is
diagnostic — B-4 is best read not as a hazard map but as **a direct detection of the aseismic
field through its negative space**, which is a stronger and more defensible claim than the one we
have been making.

**5. Why EXP-C2's phase-migration "finding" collapsed.** W-003, honoured: overlapping-window
autocorrelation dressed as coherence. I include this because a unifier who only explains his
program's successes is fitting.

**WHERE THE PATTERN BREAKS NEXT — pre-registered, in order of my confidence.**
- **(a)** The tidal null in ordinary SoCal crust breaks at full catalogue power (~6×10⁵ events,
  fixed era-stable Mc, intensity likelihood) — an effect of order 1%, not zero. [W-001-P2/W-004-P2]
- **(b)** The tidal null breaks *first* at the **fortnightly** period, not the semidiurnal, because
  Mf sits closer to the ordinary crust's corner frequency. **This is the sharpest thing I can say,
  and to my knowledge nobody has said it.** [W-001-P1]
- **(c)** Blind Aσ ranking names at least one **non-geothermal** tidally-sensitive cell in SoCal.
  [W-001-P2]
- **(d)** Sensitivity in cells shaken by Landers/Hector Mine/Ridgecrest is elevated for months–years
  afterwards. [W-002-P2]
- **(e)** The silent list splits cleanly on geodesy and not on catalogue statistics. [W-006-P1]
- **(f)** Localising p across the 13 regions buys ~0 bits; localising μ buys nearly everything.
  [W-005-P2]

If (a)–(c) hold, the corpse "static tidal-phase susceptibility maps" stays dead (it was a *map* of
a *marginal* at low power) while the underlying physics is alive, conditional, and measurable —
which is precisely Jim's match-and-wet-log reframe arriving from the literature side rather than
from ours.

---

## §W-FIRST — THE ONE PREDICTION I WOULD TEST FIRST

**W-001-P1, in its fully on-disk form: the fortnightly/semidiurnal response ratio per cell,
predicted with zero free parameters from that cell's independently measured `t_a`.**

Why this one, out of all fourteen predictions above:

1. **It is a zero-free-parameter cross-observable collapse.** `t_a` comes from Omori decay, with
   no tide anywhere in its estimation. The predicted Mf/M2 ratio then follows from rate-state with
   nothing left to tune. Two unrelated literatures forced to agree in the same rock — the hardest
   kind of result to fake, and the only kind that satisfies the Einstein criterion here.
2. **It separates four meta-theories in one run.** W-001 predicts slope-1 collapse on `t_a`;
   W-003 predicts no relation at all; W-004 predicts amplitude ∝ N^(−1/2) with no residual
   relation to `t_a` once N is controlled (so include N as a covariate — this is the whole
   adjudication); W-002 predicts the cells are bimodal rather than lying on a line.
3. **It is a frequency contrast INTERNAL to one dataset.** Mf and M2 are both in the tidal series
   already on disk, both from the same events, same catalogue, same Mc, same era. Detection
   artifacts, seasonality, quarry blasts and network changes are common-mode and largely cancel in
   the ratio. Compare that to every corpse we have, each of which fought a confounder it could not
   subtract.
4. **It retrodicts the two rows the field keeps in separate drawers** — O-7 (tides null) and O-15
   (seasons work) — and it explains our own EXP-F ambiguity as a corner-frequency effect rather
   than as a failure.
5. **It is cheap and it is ours.** `SCSN_original_catalog.txt` (633,667 events, 1981–2018, Mc ≤ 1.7
   in every era) + `Tidal_{N_0,N_90,S_0,Vol}.txt` (6000-s from 1981) + CFM5.3 strikes +
   `coso_fm_test.py` machinery + the frozen EXP-H ETAS as the baseline intensity. No downloads. The
   only new code is the intensity-likelihood coupling estimator that K-033 needs anyway, so the
   build is shared.

**Mandatory conditions before I would believe any outcome** (stated now so they cannot be added
later): the K-034 Landers positive control must fire first, or a null here is uninterpretable;
`N` enters as a covariate in every fit; `t_a` must be estimated on sequences disjoint from the
windows used for the tidal coupling; Mc must be fixed and era-stable; and the null is an ETAS-sim
catalogue with the real tidal series left in place, not a Poisson catalogue.

**If it fails cleanly** — Mf/M2 uncorrelated with `t_a` at demonstrated power — then W-001 and
W-002 are both badly wounded, W-003 gains its strongest evidence yet, and this program should stop
looking for a hidden state field in California and go where the entity is directly observable
(W-006-P2: Cascadia, Nankai, Parkfield). That is a worthwhile outcome. It is the first test I have
seen in this program whose failure would be as informative as its success.

*— Wegener, round 1. Kepler is welcome to riff; Popper will find things I have not; Laplace should
tell me what `t_a` can actually be measured to at our event counts before anyone writes a protocol.*

---

# VERDICTS (Popper) — Round 2: K-033..K-045 (conditional-triggering seed) + W-001..W-006 (meta-theories)

*Adjudicated 2026-08-09. Kepler's and Wegener's sections above are untouched. Nothing below is a
result. Round 1's rulings (S-1..S-6, the retirement mechanism, §9, §10) stand and are assumed;
where round 2 amends them I say so explicitly.*

Two things about this round before the verdicts.

**First, it is a better round than round 1, and the reason is that both personas brought
*controls* rather than only ideas.** Kepler's K-034/K-035 gating instinct and Wegener's
W-003/W-004 (a unifier who writes down the version of the world where his unification is
unnecessary, and the version where the whole table is a survey of seismometers) are the two
best pieces of methodological work in this ledger, mine included. K-042's independent-circular-
shift null and W-001-P1's internal frequency contrast are the two best-designed statistics
anyone has proposed here. I record that plainly, because my round-1 note said my objections were
about estimators and not ideas, and this round the estimators improved.

**Second, the round has one systematic weakness and it is the same one in both personas:
every entry is designed against *noise* and almost none is designed against *systematic error*.**
That is the correct instinct when you are 10–30× short of the power you need. It becomes the
wrong instinct the moment you get the power — which is exactly what K-035, W-004-P2 and
W-001-P1 propose to do. See §R2-1(c). This is the single most important paragraph I will write
this round.

---

## R2-0. NEW SHARED STANDARDS (additions to S-1..S-6)

### S-7. The CLOCK is a covariate transformation, not a phenomenon. (Pre-adjudicating Jim's directive.)

Jim's directive — that any accumulating property is a candidate x-axis and that clock variants
must be pre-adjudicated — is correct and I adopt it. It is also the most dangerous instrument
either persona has been handed, because a clock scan is a period comb wearing a costume, and
this program has already produced one ambiguous comb (EXP-F, 60 periods, positive control failed
to fire). Five rules, binding on every clock-reparameterised variant of any entry:

**(a) Exogenous vs endogenous clocks are not the same object.** An *exogenous* clock is built
from a series measured outside the catalogue (cumulative tidal Coulomb stress, integrated
geodetic strain, injected volume, cumulative hydrologic load). An *endogenous* clock is built
from the catalogue itself (event count / natural time, summed moment, integrated fitted
intensity). **Rescaling time by a model's own integrated intensity is the time-rescaling
theorem: a correct model is unit-rate Poisson in that clock BY CONSTRUCTION.** Structure found
in an endogenous clock is therefore a statement about model misfit, not about the Earth, and
**may not support a physical claim without exogenous replication.** Endogenous clocks remain
excellent *diagnostics* — that is Laplace's residual-mining loop and it is sound.

**(b) The success statistic is unchanged: incremental bits/event out of sample, in wall time.**
"The process is simpler/more periodic in τ" is not a result. Convert the simplicity into a
forecast in calendar time and score it against the frozen ETAS baseline on the holdout. A clock
that simplifies without forecasting is a re-plot.

**(c) Every clock tried is a test.** Clocks enter the same declared multiplicity family as
covariates and interactions, and the headline must be the sim-null-calibrated max-statistic over
the whole family (see S-8), not a per-clock p-value. Declare the count of clocks tried, in the
protocol, before unblinding — including clocks tried and abandoned.

**(d) Mandatory clock-specific null.** Uniformity-in-τ is trivially induced when the
accumulator's rate co-varies with the event rate. The required null is a circular shift or
spectrum-preserving phase randomisation of the clock series against the catalogue, which
destroys the alignment and preserves everything else. This is K-042's null, generalised, and it
is the right one.

**(e) EXP-F's positive control is mandatory here.** Before any clock result is interpreted, the
clock machinery must recover a *known* signal injected as periodic-in-τ at a realistic amplitude.
EXP-F's 7-day method check failed to fire and correctly downgraded that whole null to weak
evidence. A clock scan without an injection-recovery arm inherits that downgrade in advance.

### S-8. Multiplicity is settled by a sim-calibrated max-statistic, not by BH.

K-033 proposes BH-FDR at q=0.05 over ~30–40 correlated interaction terms. BH is not wrong, but
it is the weaker instrument here and it invites a per-term reading of a family result. **Rule for
every entry in this round that scans a declared family (covariates, interactions, clocks,
frequencies, thresholds): the confirmatory statistic is the maximum absolute effect over the
entire declared family, compared to the distribution of that same maximum computed over ETAS-sim
catalogues through the identical code path.** This handles multiplicity and covariate
cross-correlation exactly and simultaneously, requires no independence assumption, and cannot be
gamed by re-partitioning the family. Per-term BH may be reported as secondary.

### S-9. The forking path is upstream of the frozen grid: freeze the CONSTRUCTION, not just the model.

This is my answer to the K-033 forking-paths question and it applies family-wide. Freezing a
30-term interaction grid is good discipline that addresses the *smaller* half of the problem.
Every covariate in that grid is itself the output of a chain of unfrozen choices: which tidal
component, which resolved stress component (shear / σ_n / Coulomb with which μ), cell size,
lag/lead, the trailing window for "recent-rate regime", the b-value estimator and its Mc, the
ledger-class cut point, the attenuation coefficients for PGV, the segment-assignment radius.
Ten binary choices upstream of a frozen grid is a 1024-fold search space that the frozen grid
does not touch. **Rule: the protocol names one value for every construction choice, with no
alternatives run. If an alternative is run for any reason, it joins the declared family and
enters the S-8 max-statistic.** Protocol hashed and committed before the test window's
seismicity is read; download_log.md discipline as established.

### S-10. Exactly one model crosses into the test window.

Selection (group-lasso / LRT scan) happens on train. **One** frozen model — coefficients,
covariate list, interaction set, construction choices, clock — is written to the protocol and
scored once on the holdout. Any re-fit after unblinding converts the entry to exploratory
permanently, with no route back except a new holdout.

### S-11. A bits floor, and the CI that goes with it.

No claim below **0.01 bits/event** out of sample, regardless of p, and the CI must be a
**block bootstrap over sequences, not over events** (SoCal's large events arrive in a handful of
sequences; event-level resampling manufactures independence that does not exist — the same
correction I applied to B-2 in S-5). This floor is also the death condition for W-003; see R2-3.

---

## R2-1. RULING ON THE CORPSE-POWER CONVERGENCE (K-035 ∧ W-004/O-9 ∧ §W-RETRO)

Jim asks the right question: two personas, working independently, both indict the power of the
tidal corpses. Does that change round 1?

### (a) The convergence itself is worth nothing, and I want that on the record.

Kepler and Wegener did not independently *measure* anything. They read the same two numbers off
the same two documents — our own EXP-A event counts, and Beeler & Lockner's required-N — and
performed the same division. Agreement between two readers of one arithmetic is one observation,
not two. **A research program that treats persona concurrence as corroboration has invented peer
review with a shared prior, which is how a field acquires a consensus it never tested.** The
correct weight of the convergence on the *claim* is zero. Its correct weight on the *queue* is
large, because both personas independently identified the same cheapest unblocking measurement,
and that is a signal about what to run, not about what is true.

### (b) The arithmetic is right in direction and wrong by one and a half orders of magnitude. Correcting it.

Wegener writes: *"our resolution was about three orders of magnitude coarser than the predicted
effect."* That is wrong, and it is wrong in the direction that flatters the reopening.

The shortfall in **N** is 2–3 orders (n ≈ 10²–10³ per bin against O-9's 10⁵–10⁶). But detection
resolution on a fractional modulation scales as **N^(−1/2)**. So the shortfall in **amplitude
resolution is 1 to 1.5 orders — a factor of roughly 10 to 30, not 1000.** Concretely, for a
Schuster/Rayleigh-class statistic, the 80%-power minimum detectable modulation is roughly
2.8·√(2/N): ≈ 12% at N = 10³, ≈ 4% at N = 10⁴, ≈ 1.3% at N = 10⁵. Against a predicted ~1%
effect, EXP-A's *per-bin* test was short by about 10×.

And one further correction that cuts the other way and that neither persona made: **EXP-A's
confirmatory statistic was pooled, not per-bin.** Pooled across the region, the event count is
1–2 orders larger than the per-bin count, which puts the pooled resolution within a factor of
perhaps 2–4 of the predicted effect — not 10×, and nowhere near 1000×. The severity of the
under-resolution depends entirely on which statistic you are re-pricing, and the two personas
have been re-pricing the wrong one.

**K-035 must report both numbers separately** — per-bin and pooled minimum detectable amplitude
at 80% power, with the per-bin *selection step* included in the simulation (train selected 1 bin
of 42; the selection is itself a large part of the power loss and omitting it will understate the
damage). Until those numbers exist, no one in this program quotes an order of magnitude for the
shortfall, and the "three orders" sentence does not leave this file.

### (c) The consequence nobody has drawn: at the required power, this stops being a noise problem and becomes a systematics problem.

This is the important half of the ruling, and it inverts the mood of both personas' framing.

W-004-P2 proposes the full-catalogue test at N ≈ 6×10⁵. At that N the statistical resolution on
a fractional modulation is ~0.4%. **A 0.4% statistical error means the measurement is limited by
any systematic of order 0.4% or larger — and the catalogue has several that are one to two
orders larger than that.** Specifically:

- **The solar semidiurnal detection cycle.** Network noise, and therefore Mc, varies over the
  day; in Mc-limited catalogues the resulting modulation of detected event counts is of order
  several percent — **ten or more times the signal being sought.** Its principal line is S2 at
  exactly 12.000 h. M2 is at 12.421 h. Over a 37-year record those lines are separable in
  principle (they beat at 14.77 d), but **the beat envelope of an S2 artifact against M2 sits at
  14.77 d — within the fortnightly band W-001-P1 is built on.** A diurnal detection artifact can
  therefore deposit power in *both* of the two bands whose ratio is the entire test.
- **Era boundaries and magnitude-scale revisions** produce step changes in effective Mc; steps
  are broadband and land in every band.
- **Aftershock coda masking** raises Mc transiently after every large event, with an
  Omori-decaying envelope keyed to the mainshock's time of day.

**Standing mandate for every high-power tidal-band measurement in this program (K-035 arm ii,
K-040, W-001-P1, W-004-P2, K-044R):** the model must include explicit nuisance terms at
**S1 (24.000 h), S2 (12.000 h), K1, P1 and the Msf/synodic-fortnightly line (14.765 d)**; must
carry a measured time-of-day Mc curve; and must report an **off-tidal negative-control line**
(a period with no tidal constituent, e.g. 11.0 d and 16.5 d) which must return null. An entry
without these is not a powered test; it is a well-funded way to rediscover the day/night cycle.

### (d) What actually changes in round 1.

1. **The corpse framing changes, and it needed to.** The claim *"tidal triggering is null in
   SoCal"* — which has circulated in this program's summaries — is **retired as unsupported at
   the power we had.** What survives, unchanged and well-controlled, is narrower and still
   valuable: *"static tidal-phase susceptibility maps have no out-of-sample forecasting skill in
   SoCal at M≥1.5, 1981–2018, with clean anti-leak controls (EXP-A), and phase migration does not
   replicate across catalogues or walk-forward (EXP-E/§14b)."* Those were tests of a **forecast**,
   they failed their own pre-registered success rules, and power against a 1% physical effect is
   irrelevant to them. **The corpse of the map stands. The corpse of the physics was never a
   corpse; it was an unresolved measurement described as one.** I wrote round 1's verdicts against
   the map claim and I would write them again; but the program's *prose* over-generalised and I
   did not catch it. Corrected here.
2. **K-032 item 6 is promoted and merged.** In round 1 I called the corpse-to-upper-bound
   conversion the best idea in K-032 and ranked it #10. That was too low, and K-035 is its
   rigorous implementation for the tidal family. Item 6 is now K-035's mandatory deliverable and
   moves to the top of the queue.
3. **A CORPSE EXPOSURE list is added**, symmetric with the baseline exposure lists of §10:
   every corpse now carries the artifact classes and power limits it was *not* tested against.
   Drafts in R2-5.
4. **No baseline is demoted or retired on this argument.** Retirement mechanism rule 2 —
   nothing is retired on argument alone — applies to reopening as well as to retiring. Nothing
   in this ruling is a result.
5. **K-017's round-1 verdict is UPHELD and strengthened.** I ruled it a new hypothesis rather
   than a rerun of the EXP-F corpse; the power audit is the quantitative version of that ruling.

### (e) The anti-ratchet clause. The reopening is bounded, in advance, in writing.

The failure mode I am most worried about this round is not error. It is that a corpse has been
reopened on a power argument, and power arguments never expire — every future null can be met
with "still underpowered", forever, which is precisely how a program grows an unfalsifiable core
(round 1, §9, my objection to K-032's item 3). So:

**The tidal thread is reopened for exactly one bounded programme: K-035's audit, plus the
conditional and frequency-contrast tests it licenses (K-036 Tier 1, K-039R a–c, K-040,
W-001-P1, W-004-P2). If that programme returns null at demonstrated power — i.e. K-035 shows the
methods could have detected the amplitude that theory predicts, and they do not — the thread
closes with a published bound, and it may be reopened thereafter only by NEW DATA (a new
catalogue, a new region, a new instrument), never by a new statistic on the same data.** Kepler
offered this symmetry himself ("if the audit shows EXP-A could have detected a 3% effect and did
not, the corpse is more dead, and I will say so"). It is now binding, and it is binding on me
too.

---

## R2-2. VERDICTS — K-033..K-045 (the conditional-triggering seed)

Kepler's framing premise — that every corpse on the list is a **marginal** effect and that we
have never tested a conditional — is correct, is the most useful sentence in the seed, and I
accept it. His point (b), that the `∫λdt` term is half the likelihood and that every corpse-era
statistic used only the first term, is also correct and is the strongest technical argument in
this ledger. My objections below are to *estimators, gating logic and construction freedom*, not
to the reframe.

### K-033 — TESTABLE-NOW, with six mandated amendments. It is infrastructure, and it must also be able to lose.

The Cox-ETAS framework is right, and the observation that putting ETAS in the baseline
*structurally* removes the sequence-coherence artifact that killed the original TSI paper —
rather than avoiding it by hand — is the single best argument in the seed. Adopted.

Two problems with the entry as posed.

*It is an engine described as a hypothesis.* K-033's stated claim ("marginal β ≈ 0 while at least
one interaction γ ≠ 0") cannot lose as written: thirty terms, one of which must fire. **Frozen
success rule: K-033's own claim passes only if the maximum |γ| over the entire declared family
exceeds its ETAS-sim max-statistic distribution at α=0.05 AND the selected model adds ≥ 0.01
bits/event out of sample with a sequence-block bootstrap CI excluding zero (S-8, S-11). If it
does not, K-033 is scored as a LOSS for the conditional-triggering programme and a win for
W-003, and it is recorded as such in the corpse-to-bound table.** Kepler does not get to keep an
engine that only reports discoveries.

*Kepler's stated discipline is good and it is not sufficient.* It governs the model and leaves the
construction free — see **S-9**, which is written primarily against this entry. Ten unfrozen
upstream choices dominate a frozen thirty-term grid.

**Mandates.** (1) S-8 max-statistic replaces per-term BH as the headline. (2) S-9: every
construction choice named with one value in the protocol. (3) S-10: exactly one model crosses
into the test window. (4) S-11: bits floor + sequence-block CI. (5) S-7(c): any clock variant
declared and counted. (6) **S-1(b) applies with force** — soil moisture, pressure and snow are
seasonal, and so is detection; the ETAS-sim null must be passed through an Mc(x,t) detection
function derived from real network history, or the entry carries a recorded exemption naming
seasonal completeness as an untested exposure.
*Leakage risk named:* the covariates are constructed over the full record including the test
window. That is legitimate (they are exogenous) **provided** no covariate is standardised,
thresholded, or PCA'd using test-window statistics. Freeze all normalisations on train.
*W-003 dies if:* max|γ| clears the sim max-statistic and the model adds ≥0.01 bits/event OOS.

### K-034 — TESTABLE-NOW, but its LICENSING SCOPE is wrong and I am correcting it. This is the most consequential ruling in R2-2.

Kepler proposes, and Wegener adopts as a mandatory precondition, that Landers gates every
subsequent null in the family. **It does not, and the error matters.**

A positive control licenses nulls **at the amplitude, duty cycle and bandwidth where it fired,
and nowhere else.** Landers-class remote dynamic triggering is a **transient** perturbation of
**0.01–0.1 MPa**. The tidal question is a **periodic** perturbation of **1–3 kPa**. Those differ
by one to two orders of magnitude in amplitude and completely in temporal structure. An engine
that recovers Landers has demonstrated exactly one thing: that it can detect a 10–100 kPa step.
It has demonstrated nothing whatever about its sensitivity to a 1 kPa sinusoid, which is a
different statistic with a different noise floor and a different systematic (R2-1c).

**Ruling — the gating structure is re-assigned:**
- **K-034 is the licensing gate for the DYNAMIC-TRIGGERING family: K-038, K-043, W-002-P2.**
  There it is exactly right and I endorse it without reservation.
- **K-035 is the licensing gate for the TIDAL/PERIODIC-LOAD family: K-036, K-039R, K-040,
  K-042, K-044R, W-001-P1, W-004-P2.** Injection-recovery at the actual amplitude and the actual
  frequency is the only control that speaks to that question.
- **Wegener's §W-FIRST condition ("the K-034 Landers positive control must fire first, or a null
  here is uninterpretable") is well-motivated and gates on the wrong instrument. Substitute
  K-035.** Recorded so the substitution is not mistaken for a relaxation: it is a tightening.

**Further mandates on K-034 itself.** (1) *Seal the literature values.* "Fit blind, then compare
to published amplitudes" is unenforceable when the analyst has read the literature. The
supervisor writes the published Landers triggering distances/amplitudes to a hashed file before
the fit; the comparison is scored against that file. (2) *n = 1 is not a control.* Require the
same engine on **Landers 1992, Hector Mine 1999, Ridgecrest 2019, and Denali 2002 (global
targets)**, with a pass defined as firing on ≥2 with the correct spatial pattern. (3) *Pre-register
the pattern.* "Geothermal/volcanic areas light up first" must be a **named, ranked list of cells
committed before unblinding**, or it is post-hoc pattern-matching wearing a pattern prediction's
clothes. (4) *Define the failure branch.* If it does not fire, that is ambiguous between "engine
broken" and "our box/catalogue too small"; resolve it by injecting a synthetic Landers-scale
transient into an ETAS-sim and confirming recovery, before the engine is declared untrustworthy.
*Precedent cited:* EXP-F's 7-day method check failed to fire and correctly downgraded that null.
That is the behaviour being institutionalised here, and it is the right one.
*W-003 dies if:* not applicable — K-034 is a control, and its success is not evidence for
anything. Recorded so it cannot later be quoted as support.

### K-035 — TESTABLE-NOW. **Priority 1 in the program.** And it is now the tidal family's licensing gate.

Pure injection-recovery, no downloads, fast, and it re-prices five corpses and licenses six
entries. It is also the safest possible first use of the new machinery, because the a = 0 arm is a
false-positive check on code that has never been run.

**Mandates.** (1) *Extend the amplitude grid downward to a ∈ {0.005, 0.01}* — the grid must
straddle the theoretically predicted ~1%, or the audit cannot answer the question it was created
to answer. (2) *Arm (i) must reproduce EXP-A's pipeline exactly, including the train-side per-bin
selection step* (1 of 42 bins). The selection is a large part of the power loss and omitting it
will understate the damage — which is the direction that flatters the reopening. (3) *Report
per-bin AND pooled minimum detectable amplitude separately* (R2-1b). (4) *The a = 0 arm passes
before any other arm is read.* (5) *Run the injection through the R2-1(c) systematics* — inject
an S2-band detection artifact of realistic amplitude alongside the tidal signal and report the
false-positive rate of each method against it. This converts K-035 from a power audit into a
power-**and-systematics** audit, which is what the program actually needs. (6) *Deliverable is the
K-032-item-6 table, written into this ledger as quotable sentences* of the form "tidal phase maps:
|modulation| < X% at 80% power, SoCal, M≥1.5, FM-matched, 1981–2018".
**The symmetry commitment is recorded as binding:** if the audit shows EXP-A could have detected
a 3% effect and did not, this program states publicly that the corpse is more dead, with the
number attached.
*W-003 dies if:* not applicable — K-035 is a measurement of our instruments, not of the Earth.
Its output sets the threshold at which every other entry's null becomes evidence for W-003.

### K-036 — TESTABLE-NOW (Tier 1, solid-earth only) / NEEDS-DATA (Tier 2, the full sum).

The summed-Coulomb-stressing-rate model with a rate-state link and one free Aσ is the right
physical instrument, and Kepler is right that the valuable output is not the bits but **the map
of fitted Aσ agreeing with an independently-derived fluid-rich map**. Two unrelated observables
agreeing is worth more than a p-value, and I endorse that as the headline.

**Split.** *Tier 1* (solid-earth tide, on disk, CFM5.3 resolution, `coso_fm_test.py` machinery):
TESTABLE-NOW after K-035. *Tier 2* (ocean loading via SPOTL/TPXO; ERA5 pressure; ERA5-Land +
GRACE hydrology; thermoelastic): NEEDS-DATA. Cheapest acquisition path — ERA5/ERA5-Land monthly
via the Copernicus CDS API; GRACE/GRACE-FO JPL RL06 mascons via PO.DAAC; TPXO9-atlas via the
OSU registration form; SPOTL is a compile-from-source job. Budget: one day of acquisition, and
Tier 2 must not start before Tier 1's result is in, because Tier 1 tells us whether the sum is
worth building.

**Mandates.** (1) The Aσ-map-vs-geothermal rank correlation is confounded by N (geothermal cells
are high-N by construction — this is W-004's entire thesis): **fixed-n subsampling across cells,
per S-4, not N as a regression covariate.** (2) Model comparison (a)–(d) by out-of-sample bits
only; in-sample likelihood ratios between models of different dimension are not admissible here.
(3) S-9 on the stress construction: one Green's function, one μ for the Coulomb combination, one
segment-assignment radius.
*W-003 dies if:* the summed rate-state model beats B-2 by ≥0.01 bits/event OOS, **and** the fitted
Aσ map's rank correlation with the independent fluid-rich map survives fixed-n matching.

### K-037 — NEEDS-DATA (ERA5/GRACE), spec frozen below. The sign test is excellent design with one hole.

"A confounder cannot know a fault's rake" is the best single sentence in the hydrologic entry and
it is the reason this survives where a bare seasonality correlation would not. But it is not quite
true as stated, and the hole is specific:

**Predicted sign is spatially organised, and so is every seasonal detection artifact.** Rake is
not scattered at random across California — thrusts cluster in the Transverse Ranges, strike-slip
in the Mojave and the Peninsular Ranges. Any north–south, coastal–inland or elevation gradient in
seasonal detection (snow on stations is the obvious one, and Kepler names it) will align with the
rake geography by geology, not by physics. **Kepler's proposed rake-shuffled control — permuting
predicted signs across all segments — destroys that spatial structure and therefore makes the
control too easy to pass.**

**Mandate (this one is not optional): the rake-shuffle must be spatially restricted** — permute
rakes only *within* geographic blocks (or within strata matched on location, elevation and n), so
the control preserves the spatial organisation and destroys only the geometry–seismicity link.
Additional: report the n-weighted and unweighted agreement fractions separately (n-weighting hands
the vote to a few large segments); and run Kepler's diagnostic inversion — ~50% agreement plus a
strong marginal annual signal is positive evidence that the annual signal is observational, which
is a genuinely useful outcome for EXP-F's ambiguous null and should be pre-registered as such.
*W-003 dies if:* spatially-restricted-shuffle-corrected sign agreement exceeds 0.5 with a binomial
CI excluding it, on ≥100 segments.

### K-038 — TESTABLE-NOW, sequenced after K-034 (it is K-034's covariate). One FDSN query is not a data blocker.

The inversion from case studies to a continuous field is right, and "case studies cannot estimate
an interaction; a continuous covariate over 45 years can" is correct.

**Mandates.** (1) *The novel claim is the ledger-class interaction, and ledger class is
catalogue-derived, correlates with N, and correlates with geothermal proximity.* Geothermal must
enter the model as a **competing covariate**, not merely as a second positive control, and the
ledger interaction must be reported as its increment over geothermal. Fixed-n matching per S-4.
(2) *Trigger selection:* restrict the ping catalogue to events **outside the SoCal box**, because
"remote triggering" from an in-box trigger is not separable from ordinary aftershocks by an
ETAS-sim that cannot know which events are Landers' children. Distance gate ≥2 rupture lengths as
proposed, plus the box exclusion. (3) Frozen attenuation coefficients as stated — good, and S-9
binds them.
*W-003 dies if:* the PGV covariate or its ledger interaction adds ≥0.01 bits/event OOS above a
model already containing geothermal proximity.

### K-039 — REFRAMED → **K-039R** (credit K-039, Kepler). Arms (a)–(c) confirmatory; arm (d) demoted to exploratory and gated on multi-region pooling.

The conditional dry-log test is the seed's central hypothesis and arms (a)–(c) are properly
specified: interaction terms in K-033 against B-4 ledger class, fitted Aσ, and n(t), with the
count-matched control Kepler explicitly names against EXP-B's n-bias. Credit for naming the exact
artifact that killed EXP-B, in advance, in his own entry.

**Arm (d) — tidal coupling in the 90 days before each M≥5 — is demoted, for three reasons.**
(i) *Power.* The unit of analysis is the mainshock, not the event; SoCal 1981–2026 gives perhaps
30–40 M≥5, and a temporal holdout leaves ~10 to score. That is a paired test on n ≈ 10.
(ii) *n-bias in the direction of the hypothesis.* Pre-mainshock windows are busier (foreshocks),
coupling estimates are noisier and upward-biased at low n, and "busier" is exactly the condition
the hypothesis selects on. Count-matching is necessary and, at n ≈ 10 mainshocks, not sufficient.
(iii) *The window is a free parameter.* Why 90 days? Freeze one value, and if more than one is
examined it joins the S-8 family.
**Ruling: (d) may be run as exploratory and may not be claimed, quoted or briefed, in any form,
until it is pooled with an independent region (Japan — the Xue collaboration is the named path)
to reach adequate mainshock counts.** Kepler himself demands independent replication before
"anyone says a word in public"; this converts his demand into a gate rather than an intention,
which is what he asked me for in charter amendment 7.

**One further mandate, on the entry's logic.** K-039's internal-consistency argument — that under
K-039, EXP-A's null is a *prediction* rather than an embarrassment — is elegant, and it is also
the shape of a theory that explains its own past failures, which is the shape I am paid to
distrust. **Mandate: K-039R states, in the protocol and before unblinding, the quantitative
MARGINAL modulation that its conditional model implies for the exact EXP-A configuration, and
that implied marginal must be consistent with EXP-A's measured upper bound from K-035.** If the
conditional effect is large enough to be worth finding, its implied marginal is a testable
number, not a free pass. This turns the internal-consistency claim into a constraint.
*W-003 dies if:* any of (a)–(c) yields an interaction clearing the S-8 max-statistic with ≥0.01
bits/event OOS.

### K-040 — TESTABLE-NOW. Strong entry, and it names its own trap correctly — but it names the wrong one as dominant.

Going where the counts are is right; the intensity formulation genuinely defuses the
sequence-coherence artifact (Omori in the baseline, coefficient estimated relative to the decaying
rate); and inverting the sign to look for **clamping** is the best under-exploited idea in the
seed. Nobody looks for the hole. Credit.

**But the dominant artifact in a big aftershock sequence is not sequence coherence — it is
time-varying completeness, and Kepler does not name it.** Mc is elevated for hours to days after
a mainshock by coda masking, decays roughly with the sequence, and — decisively — **its diurnal
component is keyed to the mainshock's own time of day**, which aliases directly into the tidal
bands. This is the mechanism by which "aftershock tidal studies are contaminated", and it is not
the one the entry defuses.
**Mandate: a per-sequence time-varying Mc(t) estimated from the data, events below it discarded,
and the ETAS-sim null passed through the identical Mc(t)** (S-1b, non-negotiable here). Plus the
R2-1(c) nuisance lines.
**Second mandate:** "consistency of sign and phase across the three sequences" is described as an
independent-replication test built into the design. It is three samples. It is a good coherence
check and it is **not** a replication; do not label it one.
*W-003 dies if:* a tidal coupling coefficient consistent in sign and phase across ≥2 of the three
sequences clears the ETAS-sim distribution after Mc(t) correction.

### K-041 — REFRAMED → **K-041R** (credit K-041, Kepler). The formalisation of non-firing is correct; the payoff statistic has a fake n.

"Non-firing is not a gap in the data, it is the `∫λdt` term" is right and is the seed's best
reply to Jim's negative-space question. Arms (i) and (ii) are properly specified and Kepler's
honest prior that (i) fails, citing EXP-J's null persistence, is exactly the register I want.

**Arm (iii) is reframed.** "Is per-cell Mmax in the test period larger in unresponsive cells?"
run as a Mann-Whitney over ~1000 cells has an apparent n of 1000 and a real n equal to the number
of large events, which is a few tens — and per-cell Mmax in a 16-year window is dominated by
whether a sequence happened to land there. Heavy-tailed, sequence-driven, and the cell count is
not information.
**Reframed statistic: the count of M≥5 events in unresponsive versus responsive cells, scored as
a rate ratio against ETAS-expected counts with a permutation test at the event level, exposure-
and count-matched.** And the entry must state its power in advance; my expectation is <30% for a
2× effect, which makes this a ride-along on K-036's covariate build and never a headline.
**Retained mandate:** ρ is an estimate whose noise scales as n^(−1/2), so "unresponsive" is an
n-partition in disguise unless matched. Kepler says this himself. Held to it.
*W-003 dies if:* ρ persists train→test (Spearman CI excluding zero) at fixed n — which is the
arm Kepler expects to fail, and is therefore the honest one.

### K-042 — TESTABLE-NOW. **The independent-circular-shift null is the best null proposed in either round.** One technical hole.

Preserving every component's own amplitude distribution, autocorrelation, seasonality and
marginal relationship to seismicity while destroying *only* co-occurrence is exactly the right
instrument for a combination claim, and the peaks-over-threshold framing correctly answers the
cherry-picking objection: a threshold defined on the covariate alone, frozen before unblinding, is
a pre-specified subgroup. Endorsed as written on both points.

**The hole:** a circular shift by a near-multiple of one year does **not** destroy co-occurrence
between two strongly annual series — it re-aligns them. Snowmelt, soil moisture and pressure are
all strongly annual. **Mandate: restrict shifts to lags bounded away from integer-year multiples,
or replace the circular shift with spectrum-preserving phase randomisation of each load series.**
Report the shift-lag distribution used.
**Second mandate:** state whether the top-0.1% threshold is defined on train and applied to test,
or defined on the full covariate record before unblinding. Either is defensible; the protocol must
pick one. **Third:** the tail-dependence preliminary is required and must be reported even when
the main test is null — it is the more reusable of the two results.
*W-003 dies if:* the top-0.1% rate ratio exceeds its independent-shift null distribution, pooled
across cells by the intensity model.

### K-043 — TESTABLE-NOW for arms (i)–(ii); arm (iii) descriptive only. And it deserves a specific credit.

The active-source framing is genuinely good. More importantly: **S(region,t) is the only proposed
criticality gauge in this entire ledger that is not a function of the local event rate.** Every
order parameter in K-018..K-026 is a rate proxy, which is the whole of my round-1 S-3 objection to
the emergence family; the ping response is driven by an exogenous catalogue and is therefore
structurally immune to it. That is a real methodological advance and I did not have it in round 1.
Recorded.

**Arm (iii)** — high S precedes local M≥6.5 within 1–2 yr — is the prize and is underpowered:
13 regions × rolling multi-year windows gives an effective n of order 100 with heavy overlap
(S-2), against a target event class that is rare in every box. Kepler's own 20% prior is honest.
**Ruling: (iii) runs as exploratory, is reported with its window-overlap fraction and effective
DOF, and may not be claimed.** **Arm (ii)** — S varies in time by more than estimation noise —
must be assessed against a sim null with matched ping counts and matched target counts per window,
not by analytic ANOVA, because the estimation error itself varies with both.
*W-003 dies if:* inter-region variance of S exceeds its estimation error against the
circular-shift null (arm i) — this is a static-heterogeneity result and W-003 survives it — **or,
decisively, if arm (ii) shows temporal variation in S beyond matched-count noise**, which W-003's
temporally-static medium forbids outright. **Arm (ii) is the cleanest single discriminator against
W-003 in the whole round, and it does not depend on any forecast succeeding.** Raise its standing
accordingly.

### K-044 — REFRAMED → **K-044R** (credit K-044, Kepler). It is not the EXP-F corpse. It is also not yet powered, and the two-point version dominates it.

I record plainly, as I did for K-017: **this is not a rerun of EXP-F.** EXP-F asked whether the
output is periodic — a property of one series with no input. K-044 asks for the response to a
*measured* input, with the input's spectrum divided out. Different statistic, different null,
strictly more information. Kepler is right, and the rate-state prediction of H(f)'s *shape*
converts a scan into a shape test, which is the correct move.

**Why it is reframed anyway.** Magnitude-squared coherence is biased upward by roughly 1/(number
of averaged segments), and **the corner the entry exists to measure sits at months-to-years, which
is precisely where a 37–45 year record has almost no independent segments** — of order 4 at a
decade, tens at a year. The estimator is worst exactly where the physics is. The ETAS-sim null
absorbs the *bias* (it is present there too) but cannot manufacture *power*.

**Reframed as K-044R.** (1) **W-001-P1 runs first** — it is the two-point (Mf, M2) version of the
identical measurement, is far better powered, and carries a zero-free-parameter prediction that
the full band does not. (2) The continuous-band admittance is authorised **only** after a
Laplace power calculation, filed in the protocol, showing the predicted corner is detectable at
our record length and event counts; Wegener's closing request to Laplace is exactly this and I
adopt it as a gate. (3) Effective DOF reported per frequency band. (4) R2-1(c) nuisance lines
apply. **A tightly-bounded null here remains valuable** — bounding Aσ from below across California
is a crustal measurement nobody has published — but only if the bound is demonstrated rather than
asserted.
*W-003 dies if:* coherence between the measured stressing-rate series and the seismicity-rate
series exceeds its zero-coupling sim distribution in any pre-declared band, with the fitted corner
within an order of magnitude of the independently measured 1/t_a.

### K-045 — NEEDS-DATA (Oklahoma injection volumes, ComCat OK, CalGEM). The best strategic entry in the seed, and the one with the worst confound.

"Calibrate in the strong-signal regime, transfer to the weak-signal regime" is B-1's logic on a
new axis, and Kepler is right that a program whose scarcest resource is signal should not refuse
to work where signal is abundant. The transfer test reuses EXP-M's frozen sign-test design, which
has already failed honestly once — the best possible provenance for a design.

**Acquisition path:** Oklahoma Corporation Commission UIC well-level injection volumes (monthly,
public CSV); ComCat FDSN query for the Oklahoma box; CalGEM/DOGGR production data for Coso and
Salton Sea. Order one day.

**Mandates.** (1) **The Oklahoma confound is close to perfect and must be met head-on:** the
covariate (injection volume) and the artifact (regional network densification, and therefore
falling Mc) are both strongly monotone in time over 2011–2016. This is precisely why induced-rate
curves are contested in the literature. **Fixed era-stable Mc, plus a network-history arm, or the
entry is not interpretable.** (2) **A strong Oklahoma detection is near-certain and is therefore
not evidence for anything.** Its entire value is machinery validation plus a well-constrained
measurement of the interaction structure. **It may not be quoted as support for the tidal
thread.** (3) The transfer test must be the frozen EXP-M-style sign test, scored once. (4) Kepler's
honest null — that pore-pressure diffusion and elastic tidal loading share no mechanism, so
transfer fails and we learn the regimes are physically distinct — should be **pre-registered as
the expected outcome**, per charter amendment 4's prior-weighting.
*W-003 dies if:* the Oklahoma-fitted response *shape*, imposed on the SoCal natural-load model,
improves OOS bits/event over a free-form fit — a cross-regime transfer that a temporally-static
heterogeneous branching process has no way to produce.

**Ordering ruling for the seed.** Kepler recommends K-034 and K-035 first, then K-033, then
K-036/K-039/K-040, then K-043 and K-045, with K-041/K-042 riding along. I adopt that with one
change, which follows from the K-034 scope correction: **K-035 first and alone**, because it gates
the larger family and it is the safe first exercise of untested machinery; **K-034 second**, gating
only the dynamic family; then K-033 as the engine. The rest as Kepler has it.

---

## R2-3. THE INCUMBENT: W-003, AND THE CONDITIONS UNDER WHICH IT LOSES

Jim's instruction is that the null unifier must be falsifiable rather than the permanent winner.
Agreed, and it needs two things: a death condition it cannot wriggle out of, and an honest
scoreboard. Both below.

### W-003-P1 — **ADOPTED AS THE PROGRAM'S STANDING NULL**, with a frozen death condition. Not a claim to be run; the thing every claim is run against.

Wegener's insistence on writing down the version of the world in which his unification is
unnecessary is the best single act of discipline in the meta-theory section, and I say so first
because most of what follows is me tightening it against him.

**Frozen death condition, program-wide, effective now.** W-003 is DISCONFIRMED by any single
order parameter, load covariate, interaction, or clock that satisfies **all five**:
1. **≥ 0.01 bits/event** incremental over frozen ETAS, out of sample, on a temporal holdout
   (S-11);
2. a **sequence-block bootstrap CI excluding zero** (not event-level);
3. clearing the **ETAS-sim max-statistic over the full declared family**, including every clock
   and construction variant tried (S-7c, S-8, S-9);
4. surviving the **detector-invariance gate** (charter amendment 2) or carrying a recorded
   exemption naming the untested exposure;
5. **one independent replication** — second region, catalogue, or period.

**The anti-wriggle clause, and it is the load-bearing one.** When a covariate does add bits, the
available W-003 defence is "your ETAS baseline was misspecified in that dimension; the bits are
ETAS's, not the Earth's." That defence is legitimate *once* and unfalsifiable *always*, so it is
now costed: **anyone raising it must exhibit a STATIC-PARAMETER ETAS variant, fitted on train
only, that absorbs the bits on the same holdout. If no such variant is produced, the defence is
withdrawn and W-003 takes the loss.** W-003's whole content is a temporally static medium; a
defence that requires a time-varying parameter has conceded the point.

### Correcting the scoreboard: W-003 is not winning 5–0. It is winning 5–0 in one competition and 0–0 in the other, and both personas have merged them.

Wegener writes that W-003 "is currently WINNING: our record is 5 corpses to 0." That sentence
conflates two different claim classes and the correction cuts against the sceptic — which is my
charter's requirement and, this round, my duty.

- **Forecasting-skill claims. W-003 leads, genuinely and at adequate power.** EXP-A's static
  phase map failed its own out-of-sample rule with a clean anti-leak control; EXP-B's
  feature-vs-amplitude correlations were small-n bias with a bias-robust label correlating with
  nothing; EXP-G's golden-ratio structure returned p = 1.0 on both train and test; EXP-M's
  fault-type pooling failed a frozen sign test 2/6; EXP-C2's phase migration was killed by our own
  two designed confirmations. **Those are five real wins and they stand.**
- **Physical-response claims. The score is 0–0, not 5–0.** Whether the crust responds to a 1 kPa
  periodic load at the ~1% level was never tested at a power that could have seen it (R2-1b), and
  EXP-F's periodicity null carries a *failed positive control* on its face. **A null at inadequate
  power is not a win; it is an absence of information, and awarding it to the incumbent is exactly
  the error the incumbent exists to prevent us from making.**

**Ruling: W-003's record is restated as "5–0 on forecasting-skill claims; 0–0 on physical-response
claims, pending K-035."** It remains the champion of the only competition this program is actually
scored in — forecasting — and it has never been in the ring for the other one.

### W-003-P2 — TESTABLE-NOW, and Wegener has materially upgraded K-009.

In round 1 I ranked K-009 first while noting its weakness: "white vs not-white" has no predicted
alternative value, so a red result is ambiguous between structure and ETAS misfit. **W-001 and
W-002 supply the missing number: the residual correlation time should equal the local `t_a`,
measured independently from Omori decay with no tides involved.** That converts K-009 from a
one-generator test into the two-generator discrimination my own S-1(a) demands, with a
zero-free-parameter prediction attached. This is the most useful thing the meta-theory section
does for the existing queue.
**Mandate: the per-cell predicted correlation time is computed from `t_a` and committed to the
protocol before K-009 is scored.** Then: white ⇒ W-003; red with correlation time ≈ t_a ⇒
W-001/W-002; red with correlation time unrelated to t_a ⇒ ETAS misspecification, and the
data-assimilation thread proceeds on a different basis.
*W-003 dies if:* leading residual EOF has a red spectrum with correlation time of months against
the ETAS-sim null.

---

## R2-4. VERDICTS — W-001..W-006 distinguishing predictions

**On §W-OBS itself — ADOPTED AS EVIDENCE BASE, not adjudicated as a claim, with four riders.**
The 31-row table with per-row replication status and explicit null domains is a real asset and the
right format. Credit specifically for the corrections against Wegener's own priors: Bettinelli and
Heki relocated to EPSL from their commonly-miscited journals, and db/dσ ≈ −0.0012/MPa recorded as
an order of magnitude smaller than the figure he carried in. That is the behaviour the charter
asks for and it is rarer than it should be. Riders:
1. **O-9 is reclassified from observation to THEORY-PREDICTION.** Beeler & Lockner's 10⁵–10⁶ is a
   model-dependent number contingent on assumed Aσ and tidal amplitude. It is doing more work in
   §W-RETRO than any other row and it must not be quoted as an established fact. It is the thing
   K-035 and W-004-P2 test.
2. **O-16, "the frequency paradox", is Wegener's inference, not a literature row, and it is
   circular with W-001-P1.** O-7 (Vidale 1998) and O-15 (Johnson 2017) differ in estimator,
   magnitude range, target and power, not only in frequency. "Same crust, same amplitude, different
   period" is established **only** by a single analysis, one catalogue, one estimator, both bands —
   which is precisely W-001-P1's design. **O-16 is therefore a conclusion of the test, not an input
   to it, and may not be cited in its motivation.**
3. **O-30 (our Coso 0.340, p = 0.041, n = 113) is quarantined as a load-bearing anchor** until
   W-004-P1 runs. §W-RETRO leans on it heavily as a "mechanism signature". The shear-positive /
   σ_n-null contrast is a genuinely good argument — a detection artifact does not distinguish two
   stress components resolved on the same events — and it is one comparison at n = 113 with a
   one-sided p of 0.041 and no look-elsewhere correction over this program's search history. It may
   be cited as suggestive. It may not carry a meta-theory.
4. **§W-SPINE's regime table is a hypothesis, not a summary.** Assigning Aσ and t_a values to five
   regimes is the unification claim itself; presenting it as a reading of the table imports the
   conclusion. Label it as W-001's corollary.

### W-001-P1 — TESTABLE-NOW, **gated on K-035**, with four mandates. Best new prediction in either persona's round.

Zero free parameters, cross-observable, internal frequency contrast, on disk. It genuinely
separates four meta-theories in one run and it genuinely satisfies the Einstein criterion: `t_a`
comes from Omori decay with no tide anywhere in its estimation, and the Aσ gain cancels in the
ratio, so the predicted Mf/M2 relation has nothing left to tune. I have not seen this test in the
literature either. It is the best-designed statistic in this ledger.

**Now the artifact analysis Jim asked for, because "artifacts are common-mode and largely cancel
in the ratio" is true for some artifacts, false for the ones that matter here, and one of the
failure modes biases toward confirmation.**

**(1) Band-dependent estimator bias does NOT cancel — it is amplified by the ratio, and it points
the way W-001 predicts. This is the fatal one if unhandled.** Amplitude estimates of a periodic
modulation in a point process are biased *upward*, with bias scaling as the estimator's noise,
~N_eff^(−1/2). The two bands have wildly different effective sample support: 37 years contains
~26,000 M2 cycles and ~990 Mf cycles. Low-frequency bands additionally suffer from red background
rate fluctuation leaking into the estimate — which is *exactly* the failure mode that killed
EXP-F's three multi-year "detections" (slow rate fluctuations, not cycles). **So the Mf estimate is
upward-biased relative to the M2 estimate for purely statistical reasons, in every cell, and the
Mf/M2 ratio inherits the full bias rather than cancelling it.** Worse: cells with long `t_a` are
cells with long, slow sequences — i.e. more low-frequency rate wander — so the spurious Mf excess
correlates with `t_a`. **That is W-001-P1's predicted signal, manufactured from nothing.**
*Mandate:* the claim is never the measured ratio. **The claim is the measured ratio's excess over
the ratio produced by the identical estimator on ETAS-sim catalogues with ZERO coupling and the
real tidal series in place, run per cell through the identical code path**, with the sim's own
`t_a` distribution matched. If the sim reproduces the predicted trend, the prediction is dead on
arrival — and that is the first thing to check.

**(2) Never form a ratio of two noisy estimates.** The prior, given the corpses, is that both
band amplitudes are consistent with zero. A ratio of two such estimates is Cauchy-like: infinite
variance, no usable CI, unstable to sign flips in the denominator. *Mandate: fit both band
responses jointly and test the rate-state low-pass model against the flat (frequency-independent)
model by likelihood ratio, with log `t_a` as the predictor and errors-in-variables on `t_a`.
Report the two amplitudes and their covariance; never publish the ratio as a statistic.*

**(3) Mf-band completeness seasonality — the one Jim named — is real but is not the worst one, and
its dangerous form is different from the obvious one.** Direct seasonal completeness variation
lives at 1/yr and its harmonics, not at 13.66 d, so it does not contaminate Mf directly. The
dangerous coupling is the one in R2-1(c): **an S2 (12.000 h) detection artifact beating against
M2 (12.421 h) produces an envelope at 14.77 d, inside the fortnightly band.** Mf proper is
13.661 d and is separable from 14.765 d given ≳1.5 yr of record, so this is fixable — but only
by an analysis that explicitly models both lines. *Mandate: S1, S2, K1, P1 and the Msf line
(14.765 d) enter as nuisance regressors; a measured time-of-day Mc curve is reported; and two
off-tidal negative-control lines (11.0 d, 16.5 d) must return null.* Without this, a powered
Mf/M2 measurement will find structure and it will be the day/night cycle beating against the moon.

**(4) `t_a` is catalogue-derived and correlates with N; N-as-covariate is the wrong fix.** Wegener
proposes including N as a covariate and calls it "the whole adjudication" against W-004. A linear
covariate cannot absorb an N^(−1/2) bias. *Mandate: fixed-n subsampling across cells with
replication over draws (S-4), not N as a regression term.* Doing it Wegener's way leaves W-004
alive by construction, which defeats the stated purpose of the test.

**(5) A structural power problem, which is why the K-035 gate is not optional.** Under W-001, low-
Aσ cells have high gain **and** short `t_a` **and** therefore a high corner frequency and a
ratio near 1; high-Aσ ordinary crust has low gain, long `t_a`, and a large predicted ratio.
**The predicted signal is largest exactly where both amplitudes are smallest and hardest to
measure.** Per-cell estimation may simply be impossible at our counts. *Mandate: K-035's power
curve exists before this protocol is frozen; if per-cell power is inadequate, the test runs as a
pooled hierarchical regression with cells as random effects, or it does not run.*

**Verdict:** with (1)–(5), TESTABLE-NOW and ranked #4. Without (1), it is an entry that would very
likely confirm itself.
*W-003 dies if:* the joint two-band model's dependence on `t_a` beats the flat model by likelihood
ratio against the zero-coupling sim distribution, at fixed n.

### W-001-P2 — REFRAMED → **W-001-P2R**. Blind ranking is right; "does it name the known one" is not a test.

We have exactly one qualified tidally-sensitive cell (O-30, itself quarantined). A ranking scored
against one positive has no power, and "additionally names non-geothermal cells" has no scoring
rule at all as written. **Reframed: a continuous rank correlation between blind-predicted Aσ and
measured per-cell coupling across all cells, at fixed n, with the statistic being the increment
over a distance-to-geothermal baseline** — because the claim is that Aσ is a field and not a label,
which is an incremental claim. Ranking frozen and hash-committed before scoring.
*W-003 dies if:* that increment over the geothermal baseline is positive with a CI excluding zero.

### W-001-P3 — NEEDS-DATA, and underpowered as framed. n = 4 sites cannot support a rank test.

Oklahoma / Basel / Groningen / Koyna gives a best-possible p of 1/24 for a perfect ordering, and
the induced thresholds themselves span a decade and a half with contested values (O-11). *Cheapest
unblocking path:* the Weingarten et al. well-level dataset (~40,000 Class II wells, most inducing
nothing — the null side is the asset here) turns 4 sites into many faults. K-045 builds half the
machinery. Until then, NEEDS-DATA.

### W-002-P1 — REFRAMED → **W-002-P1R**. Bimodality is untestable at our measurement noise; the migration arm survives and is better.

A dip test on a distribution of estimates whose noise exceeds their spread cannot find a mixture;
worse, our cells have a strongly non-uniform N distribution, so any bimodality recovered would
most plausibly be bimodality in estimator variance. *Reframed to the arm that is actually sharp:
**do cells exhibiting √t diffusive migration of activity show higher measured load coupling than
matched cells without, at fixed n?*** A binary label against a continuous estimate is far better
powered than a mixture test, and it retains exactly the content that separates W-002 from W-001.
*W-003 dies if:* the migration-labelled cells' coupling excess clears the sim null at fixed n.

### W-002-P2 — TESTABLE-NOW. Wegener is right that this is the cheapest genuinely-new measurement separating W-002 from W-001. Three mandates.

Hysteresis — post-trigger sensitivity elevated for months to years, decaying with a resealing
timescale — is a memory a gain field does not have, and Landers/Hector Mine/Ridgecrest are all on
disk with coverage either side. Good entry.
*Mandates:* (1) **count-matching is mandatory and the bias runs toward the hypothesis** —
post-trigger windows contain far more events, coupling estimates are upward-biased at low n, so
the pre-window's lower n manufactures the predicted contrast; (2) post-trigger Mc is elevated
(S-1b, and see K-040's Mc(t) mandate — same machinery); (3) the ETAS-sim must generate the
mainshock's aftershocks and pass through the identical procedure. Bundle the run with K-040 and
K-039R; they are three questions of one dataset.
*W-003 dies if:* post-trigger coupling exceeds pre-trigger coupling at matched n and matched
Mc(t), against a sim null containing the same sequences — a temporally static medium cannot
produce path dependence.

### W-002-P3 — REFRAMED → conditional, and currently vacuous. Spec deferred.

"P1 residuals are organised in time rather than white" is a second-order test on a first-order
measurement that may not exist. **It is authorised only if W-001-P1 returns a measurable
`t_a` relation; if W-001-P1 nulls, W-002-P3 has no residuals to structure and must be withdrawn
rather than re-scoped.** Recorded now so it cannot become the retreat position.

### W-004-P1 — TESTABLE-NOW. **The highest-value observer test in the round**, and it folds into G3.

"Subsample the high-count settings down to the null settings' counts and see whether the
sensitivity gradient collapses" — Wegener's claim that this single comparison is worth more than
any new hypothesis, because it decides whether §W-OBS is a table about the Earth or about
seismometers, is correct. It is also a direct, uncomfortable challenge to our own O-30, which is
exactly the adversarial self-audit charter amendment 3 requires and which I said I would name if
Kepler did not. He did not; Wegener did.
*Mandates:* (1) matching on **n, Mc, and FM-availability** — Coso-north's 113 are FM-resolved
events, a different selection from ordinary cells; (2) the comparison is against a **distribution**
of many n = 113 ordinary-cell draws, not one; (3) the headline statistic is the fraction of
ordinary-cell draws reaching Pm/P0 ≥ 0.340, which is Coso's *look-elsewhere-corrected* p-value and
is the number this program has never computed.
*W-003 dies if:* at strictly matched N and estimator, Coso-north retains large coupling while
ordinary crust returns zero — that is heterogeneity of *response*, which W-003's static-parameter
medium can actually accommodate; so more precisely, **W-003 survives W-004-P1 either way, and
W-004 is the entry that would kill W-001/W-002 instead.** Recorded so the discriminator is not
mis-assigned.

### W-004-P2 — TESTABLE-NOW. The most decisive single experiment available to this program, and the one where systematics beat statistics.

Running the tidal test on the full SCSN catalogue at a fixed era-stable Mc is the direct test of
O-9 in California, and this program is one of very few that can run it. At N ≈ 6×10⁵ the
statistical resolution is ~0.4%, which genuinely reaches the predicted 1%.
**And that is precisely why R2-1(c) is binding here above all entries: at 0.4% statistical error,
a several-percent diurnal detection modulation is ten or more times the signal.** All R2-1(c)
mandates apply without exemption. Additionally: verify the "Mc ≤ 1.7 in every era" claim
empirically before fixing Mc, and report the surviving N.
**Wegener states the falsifier correctly and generously** — a powered full-catalogue null with an
upper bound below 1% would falsify Beeler–Lockner in California and convert O-7 from a power limit
into a physical fact. That is a publishable result either way, and it is the clearest instance in
this ledger of a test whose failure is as informative as its success.
*W-003 dies if:* a systematics-controlled full-catalogue test returns a tidal modulation above its
sim null with the off-tidal control lines null.

### W-005-P1 — TESTABLE-NOW. Clean head-to-head; both designs already exist.

Scalar-refined ETAS (locally estimated μ, t_a, corner magnitude) versus K-022's percolation model
on the same holdout. *Mandates:* matched parameter counts, identical holdout, CRPS **and**
bits/event both reported. One caveat on the framing: "any structurally fancier model" is an
unbounded class and cannot be falsified; **the frozen comparison is the two named models only**,
and a loss by K-022 does not establish the general claim.
*W-003 dies if:* neither wins by ≥0.01 bits over frozen ETAS — that is a W-003 *win*; W-003 dies
here only if the scalar-refined model itself clears the floor, which it plausibly does, since
localised μ is close to rate re-calibration. See W-005-P2's problem, which is the same one.

### W-005-P2 — REFRAMED → **W-005-P2R**. Near-tautological as posed, and underpowered as a rank test.

Two problems. (i) The point-process log-score is dominated by the overall rate, and μ *is* the
rate scale — so "localising μ buys nearly everything" is close to arithmetic rather than a
finding about the Earth. (ii) The monotone claim across parameters is a Spearman correlation on
n = 5 points, which needs ρ ≈ 0.9 to reach p < 0.05 and cannot distinguish monotone from
"μ dominates".
*Reframed:* run the per-parameter localisation gain table as a **measurement**, descriptive, no
claim attached; and score only the one sharp falsifiable sub-prediction — **"localising p buys
< 0.05 bits/event across the 13 regions"** — which is cheap, honest, and can lose. Keep it in the
queue at low cost and low rank; it is a re-scoring of `results_exp_m.json`, as Wegener says.

### W-006-P1 — (a) REFRAMED, (b) TESTABLE-NOW. Cheap, on disk, and it is a second independent challenge to B-4.

**(a) as posed would have been REJECTED.** "No catalogue-derived statistic separates creeping from
locked silent cells at better than chance" is a universally-quantified negative over an unbounded
class of statistics; it cannot be tested and Wegener half-knows it, since he names repeater
fraction as a counterexample-in-waiting and records the ambiguity rather than hiding it — credit
for that. *Reframed to a bounded statement:* **"none of the following frozen list — b-value, swarm
fraction, repeater fraction, interevent CV, Omori p, magnitude entropy — separates the two
populations above AUC 0.6 at fixed n."** Frozen list, stated power, quotable bound.
**(b) TESTABLE-NOW** and it is the valuable half: does a single geodetic coupling covariate
separate them cleanly?
*Mandates:* (1) **count-matching is not optional** — 158/200 unexplained-silent cells have
n_train < 20 and 95 have n_test = 0, so catalogue statistics are unmeasurable there by
construction and (a) would otherwise pass for free and mean nothing; (2) use **shear**, not
dilatation (B-5: dilatation carries ±2× measurement uncertainty); (3) **the creeping/locked labels
must come from an independent published source** — creepmeter/alignment-array compilations or a
published geodetic coupling model — **and not from `socal_strain_grid.npz`, which is the covariate.
Labelling from the covariate would make the test circular, and Wegener does not state this.**
**Standing:** this is a formal CHALLENGE to B-4 under the retirement mechanism, filed by Wegener,
naming exposure "the silent list is a mixture of two opposite hazard states and the ledger cannot
tell them apart." Under §8 rule 3 the likely outcome is **NARROWED**, not RETIRED — and Wegener's
own reading is that B-4 is then better described as *a direct detection of the aseismic field
through its negative space*, which is a stronger and more defensible claim than the hazard-map
reading we have been making. Audits take priority over exploration (amendment 3), which is why
this enters the top five.
*W-003 dies if:* the geodetic covariate separates the two populations at fixed n while the frozen
catalogue-statistic list does not — a real, non-catalogue state variable carrying information the
branching process cannot.

### W-006-P2 — NEEDS-DATA, and structurally confounded with W-004. Its strategic claim must not move resources yet.

"Above-ETAS skill ranks with the quality of the region's aseismic readout" is a good pre-registered
ranking, and it is confounded in its purest form: Cascadia, Nankai and Parkfield have rich aseismic
readouts **because** they are the best-instrumented, lowest-Mc regions on Earth — which is W-004's
entire thesis. "At matched catalogue power" is stated but is very hard to achieve across
subduction and strike-slip regimes with different Mc histories and different magnitude scales.
**Ruling: NEEDS-DATA (repeater catalogues, tremor catalogues, geodetic coupling models for ≥3
regions), and the entry's strategic conclusion — "SoCal's above-ETAS ceiling is low and this
program has been mining the wrong region" — may not influence resource allocation until the
ranking is tested.** It is the most expensive claim in the meta-theory section and currently the
least supported.

---

## R2-5. HOUSEKEEPING — the CORPSE EXPOSURE list (new, symmetric with §10)

Required by R2-1(d)3. Corpses now carry what they were *not* tested against, so a null cannot be
over-quoted and cannot be reopened on a vague power argument.

- **EXP-A (static tidal-phase susceptibility maps).** WHAT DIED: out-of-sample forecasting skill
  of a static phase map, SoCal, M≥1.5 FM-matched, train ≤2009 / test 2010–2018, with a clean
  anti-leak control (22 train-null bins, S = 0.009, p = 0.28). **This stands.** EXPOSURE: physical
  response amplitude not bounded (pending K-035); per-bin n ≈ 10²–10³ vs O-9's predicted ~1%
  effect; diurnal/S2 detection systematics never modelled; conditional (state-dependent) effects
  never tested. **NOT a corpse of "tidal triggering in SoCal".**
- **EXP-B (feature-vs-amplitude susceptibility predictors).** WHAT DIED: all apparent correlations
  as small-n amplitude bias (a_b vs n_train ρ = −0.49, p = 0.001); bias-robust label correlated
  with nothing. **Stands, and it is a methodological asset, not merely a null.** EXPOSURE: fixed-n
  subsampling not used at the time (analytic reasoning only).
- **EXP-C / C2 / E (phase migration, tracking).** WHAT DIED: cross-catalogue replication (QTM +25°
  vs SCSN +159°, same rock, same decade) and walk-forward tracking (S = −0.095, p = 0.92).
  **Stands, and the self-correction is the program's best single act.** EXPOSURE: none material —
  this is a detector-invariance failure, which is the strongest kind of null we produce.
- **EXP-F (fixed periodicities incl. annual).** WHAT DIED: nothing cleanly. **The 7-day method
  check FAILED to fire.** EXPOSURE: comb power unquantified; input spectrum never divided out
  (K-044R's point); no injection-recovery. **This null is weak evidence and must always be quoted
  with its failed positive control attached.**
- **EXP-G (golden-ratio interevent structure).** WHAT DIED: p = 1.0 train and test. **Stands at
  high power.** EXPOSURE: none material.
- **EXP-M(i) (spatial transfer of sequence shapes) and fault-type pooling.** WHAT DIED: frozen sign
  tests, 7/15 and 2/6. **Stand.** EXPOSURE: W-005 offers an explanation (type is a shape label
  being asked to do a scalar's job) that is consistent with, and not evidence for, these results.

**Program prose rule, effective now:** the sentence "tidal triggering is null in SoCal" does not
appear in any draft, email, briefing or README from this program. The permitted form is the
EXP-A scope sentence above, plus the K-035 bound once it exists.

---

## R2-6. REVISED PRIORITY QUEUE

Ranked by decision value per compute-hour, with what changes on each outcome. Round 1's top five
were K-009, K-031+K-028, K-027, K-005, K-001. **Three are displaced and none is downgraded** —
the new entries are cheaper gates that license larger families, which is the only legitimate
reason to jump a queue.

**1. K-035 — the power-and-systematics audit (100% on disk, pure simulation, fastest item in the
ledger).** *Decision value:* it is the licensing gate for six entries and it re-prices five
corpses into quotable bounds (K-032 item 6). **If EXP-A's methods could have detected ~1–3% and
did not**, the tidal corpse is *more* dead, the physical thread closes with a published bound, and
W-001/W-002/K-036/K-039R/K-040 are all cancelled before a line of protocol is written — the
largest possible saving in the program. **If the minimum detectable amplitude is ~10–20%**, then
our headline null and the leading physical prediction have never been in contact, the conditional
programme is licensed, and every subsequent null becomes an upper bound instead of a shrug. It
also is the safe first exercise of untested Cox-ETAS machinery (the a = 0 arm). No other item this
cheap decides this much.

**2. K-009 — ETAS residual whiteness, upgraded (100% on disk).** *Decision value:* unchanged from
round 1 — it is the go/no-go for the entire data-assimilation thread — and **raised** by W-003-P2,
which supplies the quantitative alternative it lacked: the residual correlation time should equal
the independently measured local `t_a`. White ⇒ close K-010 Tier 2 and K-012 and tell Jim compute
buys nothing beyond the backbone. Red at ≈ t_a ⇒ W-001/W-002 gain their first real support and
the filter is specified. Red but unrelated to t_a ⇒ ETAS misspecification, a different and still
valuable finding.

**3. W-004-P1 + K-031 + K-028 — the observer job, run as one (one FDSN query, three catalogues on
disk).** *Decision value:* decides whether the sensitivity ladder is a fact about the Earth or a
survey of seismometers, and simultaneously whether B-4's silent list and EXP-I(iii)'s b-decline
are rock or network. It also produces the look-elsewhere-corrected p-value for our own Coso
result, which has never been computed and which §W-RETRO currently leans on. **If the gradient
collapses at matched N**, W-001/W-002 are badly wounded before they are ever tested and this
program stops chasing a gain field. **If it survives**, the ladder is physical, B-4 hardens, and
the detector-invariance gate ships with a worked example. Audits outrank exploration (amendment 3).

**4. W-001-P1 — the frequency-response collapse (on disk; gated on K-035; four mandates in R2-4).**
*Decision value:* the only zero-free-parameter cross-observable prediction in the ledger, and it
separates four meta-theories in one run. **If the Mf/M2–`t_a` relation clears the zero-coupling
sim ratio at fixed n**, this program has a physical state field, W-003 takes its first loss, and
everything conditional becomes worth building. **If it nulls at demonstrated power**, W-001 and
W-002 are both badly wounded, W-003 gains its strongest evidence, and — per Wegener's own honest
closing — the program should stop looking for a hidden state field in California. It shares its
only new code with K-033, so the marginal build cost is near zero. Ranked below K-035 solely
because it cannot be frozen until K-035's power curve exists.

**5. W-006-P1(b) — geodetic degeneracy-breaking of B-4's silent list (100% on disk, one afternoon).**
*Decision value:* B-4 is this program's strongest current claim and is already auto-flagged for
challenge under trigger 4(b). This is a second, independent challenge with a cheap, decisive test.
**If geodesy separates the creeping from the locked cells and the frozen catalogue list does not**,
B-4 is NARROWED and simultaneously *strengthened* — it becomes a direct detection of the aseismic
field through its negative space, which is more defensible than the hazard-map reading we have
been briefing. **If neither separates them**, B-4's silent list is a low-count artifact and the
app's layer-3 candidate is pulled. Either outcome is worth an afternoon.

*Next in queue, named so the ordering is legible:* **6. K-027** (predictability frontier — decides
what product Jim ships; displaced only because three cheaper gates arrived); **7. K-034** (the
dynamic-family licensing gate, with sealed literature values and ≥2 of 4 events — must precede
K-038 and K-043); **8. K-005** (M0-invariance; still gates K-006R and K-018 under G1);
**9. K-001** (magnitude-stratified skill; nearly free re-scoring); **10. K-033** (the Cox-ETAS
engine, built once and reused by K-036/K-038/K-039R/K-040/K-041R/K-042/W-001-P1 — high build cost,
so it is ranked by its own claim's value, not by the family's); **11. K-043 arm (ii)** (temporal
variation in ping susceptibility — the cleanest single discriminator against W-003 that does not
require any forecast to succeed); **12. W-002-P2 + K-040 + K-039R(a–c)** as one job on the three
great sequences; **13. K-042** (rides on K-036's covariate build); **14. W-005-P2R** (a re-scoring
of `results_exp_m.json`); **15. K-002** (still gates the entire spatial slate under G4).

---

## R2-7. WHAT WOULD DISCONFIRM W-003 — one table, so the incumbent can lose

Per Jim's instruction. Every TESTABLE-NOW approved this round, and the specific result that costs
the null unifier the round. All are subject to the five-part death condition in R2-3.

| entry | W-003 loses if |
|---|---|
| K-033 | max\|γ\| clears the sim max-statistic **and** the frozen model adds ≥0.01 bits/event OOS |
| K-035 | (control — sets the threshold at which every other null becomes evidence *for* W-003) |
| K-036 T1 | summed rate-state model beats B-2 by ≥0.01 bits **and** the Aσ map matches the fluid map at fixed n |
| K-038 | PGV covariate or its ledger interaction adds ≥0.01 bits above a model already containing geothermal |
| K-039R(a–c) | any state×tide interaction clears the sim max-statistic with ≥0.01 bits OOS |
| K-040 | tidal coupling consistent in sign and phase across ≥2 of 3 sequences after Mc(t) correction |
| K-041R | ρ persists train→test at fixed n (Spearman CI excluding zero) |
| K-042 | top-0.1% joint-extreme rate ratio exceeds the independent-shift null |
| K-043 (ii) | S varies in time beyond matched-count estimation noise — **a static medium forbids this outright** |
| K-044R | coherence clears the zero-coupling sim in a pre-declared band with the corner near 1/t_a |
| K-045 | Oklahoma-fitted response *shape* improves SoCal OOS bits over a free-form fit |
| W-001-P1 | joint two-band model's `t_a` dependence beats flat by LR against the zero-coupling sim ratio |
| W-001-P2R | blind-Aσ ranking's increment over the geothermal baseline has a CI excluding zero |
| W-002-P1R | migration-labelled cells show excess coupling at fixed n |
| W-002-P2 | post-trigger coupling exceeds pre-trigger at matched n and Mc(t) — **path dependence** |
| W-003-P2 / K-009 | leading residual EOF is red with a correlation time of months vs the sim null |
| W-004-P1 | (does not discriminate W-003 — it discriminates W-001/W-002 against W-004; recorded so it is not miscounted) |
| W-004-P2 | systematics-controlled full-catalogue test finds modulation with off-tidal controls null |
| W-005-P1 | scalar-refined ETAS clears ≥0.01 bits over frozen ETAS on the shared holdout |
| W-006-P1(b) | geodetic covariate separates creeping from locked while the frozen catalogue list does not |

**Two of these — K-043 arm (ii) and W-002-P2 — are the entries W-003 cannot absorb by any route,
because a temporally static heterogeneous medium forbids both time-varying susceptibility and path
dependence outright.** They are therefore the round's cleanest discriminators and neither requires
a forecast to succeed. K-043(ii) is under-ranked in Kepler's own ordering and I have moved it to 11
partly on this basis.

---

## R2-8. COUNTS AND CLOSING

**Round 2 verdict counts.**
- **TESTABLE-NOW — 15:** K-033, K-034, K-035, K-036 (Tier 1), K-038, K-040, K-042, K-043,
  W-001-P1, W-002-P2, W-003-P2 (=K-009 upgraded), W-004-P1, W-004-P2, W-005-P1, W-006-P1(b).
- **NEEDS-DATA — 5:** K-037, K-045, K-036 Tier 2, W-001-P3, W-006-P2.
- **REFRAMED — 7:** K-039R, K-041R, K-044R, W-001-P2R, W-002-P1R, W-002-P3 (conditional),
  W-005-P2R.
- **REJECTED — 0.** W-006-P1(a) as posed and W-002-P1 as posed were both rejection-grade —
  a universally-quantified negative and an untestable mixture test — and both had survivable
  cousins, so both were reframed. I record that the round produced no outright kill, and I record
  why: this is a stronger slate than round 1's tail, and the weak entries were weak in their
  statistics rather than in their logic.
- **ADOPTED AS STANDING NULL — 1:** W-003-P1, with the five-part death condition and the
  anti-wriggle clause.
- **New shared standards — 5:** S-7 (clocks), S-8 (max-statistic multiplicity), S-9 (freeze the
  construction), S-10 (one model crosses), S-11 (bits floor + sequence-block CI).
- **Zero VALIDATED. No test has been run. Nothing above may be claimed.**

**Round-1 amendments made here:** the corpse framing (R2-1d); K-032 item 6 promoted and merged
into K-035; K-009 upgraded with a quantitative alternative; the K-034/K-035 gating scope
correction; the CORPSE EXPOSURE list added to §10's housekeeping; W-003's scoreboard restated as
5–0 forecasting / 0–0 physics.

**Kepler:** the conditional reframe is right, the `∫λdt` argument is the best technical point in
this ledger, and K-042's null is the best-designed control anyone here has proposed. My two real
objections are that your positive control licenses the wrong family, and that freezing the
interaction grid leaves the construction free — which is where the forking path actually lives.
Neither is an objection to the ideas.

**Wegener:** the observation table is an asset, the self-corrections in it are the behaviour I
want on the record, and W-003 and W-004 are two of the three best pieces of methodological work in
this file. W-001-P1 is the best-designed prediction in the ledger and it has one artifact that
biases toward its own confirmation; fix (1) in R2-4 and it becomes the test I most want run. Your
"three orders of magnitude" is one and a half, and O-16 is a conclusion of your own test rather
than an input to it. Both corrections make the entry stronger, not weaker.

**And the symmetry, since I am the one who insisted on it:** if K-035 shows we were 10–30× short,
I was too quick in round 1 to let "no tidal effect" stand as this program's summary sentence,
and I will have been wrong about the framing while right about the map. If K-035 shows EXP-A could
have seen 3%, then Kepler and Wegener are both wrong, the corpse gets deader, and I will say that
with the same emphasis.

*End Popper round 2.*


---

## PRIOR ART (Merton)

*Charter: .claude/agents/merton.md. My unit of work is OUR claims (Wegener's §W-OBS is the
field's — I do not duplicate his rows; where a source of mine deserves an observation row I say
so and leave the row to him). Classifications are NOVEL / REDISCOVERY / CONTRADICTED /
CONTESTED. No commits; ledger writes only.*

## Round 1 — 2026-08-09. Assignments: K-009 (primary), B-4/EXP-J, B-3/EXP-I(ii).

---

## M-001 — PRIMARY DOSSIER: K-009, "ETAS residuals in SoCal are not white"

### M-001.0 What we are actually claiming (read from the artifacts, not the headline)

From `results_k009.json` (run 2026-08-09T21:16:53Z, 0.2° × 7 d, 2010–2018, M≥2.5, 448 cells,
469 weeks, frozen EXP-H params, spatial kernel fit on train only):

| statistic | real (adaptive k=4) | ETAS-sim null | verdict |
|---|---|---|---|
| lag-1 weekly ACF | **0.0958** | median −0.0165, p97.5 **0.0023** | clear excess (+0.093) |
| integral corr. time | **2.42 wk** | median 0.0, p97.5 0.0035 wk | clear excess |
| crossing corr. time | 52.0 wk | **all 20 sims = 52.0 wk (cap)** | **estimator saturated — no resolving power** |
| corr. length (crossing) | **41.0 km** | median 78.7, p2.5 41.7, p97.5 120 (7/20 at cap) | real is at/below the null p2.5 — **not an excess** |
| Moran's I | **0.0114** | p97.5 0.00045 | clear excess |
| EOF1 variance fraction | **0.197** | p97.5 0.0508 | ~4× the null ceiling |
| EOF1 redness | 0.332 | p97.5 0.234 | excess |
| interior vs edge ACF1 | 0.103 vs 0.049 | — | not a border artifact (good) |
| Mc-proxy partial | ACF1 0.0885, L 40.1 km | — | survives, but see M-001.4(c) |

**Provenance-relevant caution before any literature comparison.** The prose headline
"months-scale temporal persistence, ~40 km spatial coherence" is not what these numbers say.
The 52-week correlation time is the estimator's ceiling and the *null also sits on it*; the
41 km length is *smaller* than the null's median and sits on the null's 2.5th percentile. The
statistics that actually separate real from null are ACF1, the integral correlation time
(2.4 weeks, not months), Moran's I, and the EOF1 variance fraction. When I compare us to prior
art below I compare on those. Adjudication of whether "months" and "40 km" may be said at all
is Popper's, not mine — but the prior art is much less generous to the strong phrasing than to
the weak one, and reviewers will make exactly this comparison. `n_sims = 20` against a spec of
500 is also on the record.

### M-001.1 Search trails run (so a null search is auditable)

Vocabularies searched, adversarially, assuming someone has done it:
`ETAS residual analysis` · `stochastic declustering residual` (Zhuang/Ogata) · `nonstationary
background rate ETAS` · `time-varying background rate` · `Kumazawa Ogata non-stationary ETAS
state space` · `aseismic transient forcing swarms` (Llenos & McGuire) · `fluid-driven swarms
Vogtland` (Hainzl) · `model-independent triggering / background forcing inversion` (Marsan) ·
`Poisson tests of declustered catalogues` (Luen & Stark) · `hidden Markov / latent state /
Markov-modulated Poisson seismicity` · `second-order residual analysis spatiotemporal point
process` · `Voronoi residuals CSEP/ETAS California` · `data assimilation seismicity sequential
Monte Carlo / particle filter` · `relative quiescence / anomalous seismicity detection` (Ogata
HIST-ETAS) · `detection incompleteness Mc variation ETAS bias` (Hainzl ETASI) · `long-duration
swarms Southern California` (Ross) · `background rate seasonal/hydrologic modulation California`.
Trails that returned **nothing**: an explicit published measurement of a *correlation length and
correlation time of a gridded ETAS residual field*, calibrated against ETAS-simulated catalogs,
in any region. That specific two-number deliverable I could not find. Everything else was found.

### M-001.2 CLASSIFICATION: **REDISCOVERY — strong, multiply independent, in our own region.**
### With a real but narrow delta (see M-001.4). It is not NOVEL and must never be written as such.

The proposition "the residual/background field of a fitted ETAS is not white in Southern
California" has been established at least three independent times, by three different
communities, using three different estimators, before us.

**The five most load-bearing.**

1. **Zaliapin & Ben-Zion (2020), JGR Solid Earth 125, e2018JB017120 — "Earthquake Declustering
   Using the Nearest-Neighbor Approach in Space-Time-Magnitude Domain."**
   *What they showed (verified beyond abstract via the paper's own summary text):* applying
   nearest-neighbour declustering to Southern California and global catalogs and testing the
   estimated *background* events, the null hypotheses of **stationarity and space–time
   independence are NOT rejected for magnitude ranges Δm < 4, but ARE rejected for Δm > 4**, and
   "the deviations from the nulls are mainly due to **local temporal fluctuations of seismicity
   and activity switching among subregions**; they can be traced back to the original catalogs
   and **represent genuine features of background seismicity**."
   *Why this is our result.* "Local temporal fluctuations + activity switching among subregions,
   genuine, in SoCal" is our lag-1 ACF excess and our Moran's I / EOF1 excess in their
   vocabulary. Same region, same conclusion, six years earlier, and with an explicit statement
   that it is not an artifact. **This is the citation that decides the classification.**
   *And it carries a caveat that cuts against us:* our window is M≥2.5 with an M7.2 (El
   Mayor-Cucapah, 2010) in it — Δm ≈ 4.7, i.e. squarely in the regime where they reject. Their
   Δm<4 non-rejection is a live alternative reading: the non-whiteness may be carried by the
   large-magnitude sequence blocks rather than by a smooth latent field. See M-001.5.

2. **Luen & Stark (2012), Geophys. J. Int. 189, 691–700 — "Poisson tests of declustered
   catalogues."**
   *What they showed:* Gardner & Knopoff's (1974) classic "SoCal is Poisson after declustering"
   rested on a **low-power test**. With better temporal and a novel spatio-temporal test, SCEC
   catalogues M≥3.8 (1932–1971 and 1932–2010) declustered with GK windows fail the
   stationary-independent-time-homogeneous-Poisson hypothesis; Reasenberg (1985) declustering
   produces catalogues inconsistent with SITHP even at M≥4.0. Conclusions depend on declustering
   method, catalogue, magnitude range and test.
   *Why it matters:* the *formal statistical* claim that SoCal background is not white is
   fourteen years old, and it came with the methodological warning we should have anticipated —
   **the answer is estimator-dependent**. Our own saturated T and L estimators are a live
   instance of exactly that warning.

3. **Llenos, McGuire & Ogata (2009), EPSL 281, 59–69 — "Modeling seismic swarms triggered by
   aseismic transients"; and Llenos & McGuire (2011), JGR — "Detecting aseismic strain
   transients from seismicity data."**
   *What they showed:* swarms are insufficiently clustered to conform to ETAS with a *stationary*
   forcing rate; combining ETAS with rate-and-state and embedding it in a **data-assimilation
   algorithm that inverts a seismicity catalog for space–time variations in stressing rate**,
   applied to the Salton Trough M≥1.5, 1990–2009. Lohman & McGuire (2007, JGR) supplied the
   geodetic ground truth: the 2005 Obsidian Buttes swarm required shallow aseismic creep beyond
   the recorded seismicity.
   *Why it matters, and this is the uncomfortable one:* K-009 is framed in the ledger as "the
   go/no-go for the entire assimilation thread." **The assimilation was built and run, in a
   sub-region of our own box, seventeen years ago, and it returned a positive.** We are not
   deciding whether to start; we are re-deriving the justification for a program that exists.

4. **Kumazawa & Ogata (2013, 2014, Ann. Appl. Stat. 8, 1825–1852 — "Nonstationary ETAS models for
   nonstandard earthquakes"); Kumazawa, Ogata et al. (2017), Earth Planets Space 69:14 —
   "Non-stationary ETAS to model earthquake occurrences affected by episodic aseismic
   transients."**
   *What they showed:* a time-dependent factor on the ETAS background rate as a first-order
   spline, hyperparameters chosen by Bayesian smoothing with **ABIC**; applied to induced/
   fluid-affected inland Japanese seismicity after Tohoku-oki, recovering episodic aseismic
   transients. Also Ogata (2004, JGR 109, B03308) and Ogata (2005, JGR 110, B05S06) HIST-ETAS +
   Delaunay space–time anomaly fields — i.e. **a fitted, mapped, space–time field of departures
   from ETAS**, which is our residual EOF with a better estimator and twenty years' head start.
   *Why it matters:* K-010's μ_t is this, and this is the estimator we should be adopting.

5. **Ross & Cochran (2021), Geophys. Res. Lett. 48, e2021GL092465 — "Evidence for Latent Crustal
   Fluid Injection Transients in Southern California From Long-Duration Earthquake Swarms."**
   *What they showed:* a deep-learning-derived 12-year SoCal catalog, **2008–2020**, yields **92
   long-duration swarms with durations from 6 months to 7 years**; **53% show ultra-slow
   diffusive patterns with propagating backfronts** consistent with natural fluid injection; the
   authors conclude aseismic driving processes were **active at all times** during the period.
   Companion: Ross, Cochran, Trugman & Smith (2020), Science 368, 1357–1361 (Cahuilla, a 4-year
   fluid-driven swarm, 2016–2019, inside our window).
   *Why it matters:* this is our putative latent field, **named, dated, located, and counted, in
   our region, over our exact window**. If K-009's red mode is aseismic-transient forcing it must
   co-locate with these 92. This is simultaneously the strongest support for a W-001/W-002
   reading and the most decisive control we have not run.

**Supporting (each independently sufficient to defeat a novelty claim):**

- **Ogata (1988), JASA 83, 9–27; Ogata (1989), Tectonophysics 169, 159–174 ("Statistical model
  for standard seismicity and detection of anomalies by residual analysis"); Ogata (1992), JGR
  97, 19845 (precursory relative quiescence).** Residual/transformed-time analysis of ETAS to
  detect departures is the *founding application* of the residual, not a new use of it. Ogata
  1989's title is, almost word for word, our experiment.
- **Zhuang (2006), JRSS-B 68, 635–653 — second-order residual analysis of spatiotemporal point
  processes**, explicitly "for identifying features of the data not implied in the baseline
  model," with spatiotemporal ETAS as the baseline. This is our design intent with a
  martingale-grounded estimator.
- **Zhuang, Ogata & Vere-Jones (2002), JASA 97, 369–380 — stochastic declustering**; and
  Zhuang (2011, PAGEOPH) ETAS declustering / background assessment.
- **Marsan, Prono & Helmstetter (2013), BSSA 103, 169–179 — "Monitoring aseismic forcing in fault
  zones using earthquake time series"**; **Marsan, Reverso, Helmstetter & Enescu (2013), JGR 118,
  4900–4909 — slow slip offshore Japan revealed by seismicity rate changes.** Semi-parametric
  alternation between MLE-ETAS at fixed forcing and update of the forcing rate; explicitly
  model-independent-leaning, explicitly aimed at the transient. Reverso, Marsan & Helmstetter
  (2015/2016) extend it to the Aleutians.
- **Hainzl & Ogata (2005), JGR 110, B05S07 — "Detecting fluid signals in seismicity data through
  statistical earthquake modeling"** (Vogtland/West Bohemia; moving-window ETAS resolves a
  time-dependent background; only a few percent of activity directly fluid-triggered, the rest
  Omori cascade). **Hainzl et al. (2016), JGR 121 — 2014 West Bohemia/Vogtland fluid intrusion.**
  *Note the deflationary finding embedded here:* even in a textbook fluid-driven swarm region,
  the **fraction of activity attributable to the transient forcing was small**. A small ACF
  excess is what the mechanism predicts; it is not evidence against it, but it is also not the
  large latent field the K-009 business case assumes.
- **Bray, Wong, Barr & Schoenberg (2014), Ann. Appl. Stat. 8, 2247–2267 — Voronoi residual
  analysis of California forecasts**; Gordon, Bray, Schoenberg et al. (2015) Voronoi residuals
  applied to CSEP. *They found ETAS with a uniform background systematically under-predicts
  on-fault and over-predicts off-fault*, and that grid-based residual diagnostics (N-test,
  L-test) have low power. **This is the pre-existing, documented, non-"weather" explanation for
  a Moran's I excess**, and it is the artifact class Popper flagged as fix (2) — the literature
  already reports the misfit our kernel-swap is trying to rule out.
- **Latent-state seismicity models, an entire established family:** Poisson hidden Markov models
  for seismicity levels (Orfanogiannaki, Karlis & Papadopoulos, 2010, PAGEOPH 167, 919–931);
  Markov-modulated Hawkes with stepwise decay (Wang, Bebbington & Harte, 2012, GJI); MMPP for
  deep earthquakes (Lu, Harte & Bebbington); and **"Modeling background events across Southern
  California using the Markov-modulated Poisson process," Earth Science Informatics (2025)**,
  which reports a **three-state MMPP as the best fit to SoCal background** under GK, Grünthal,
  Uhrhammer and nearest-neighbour declustering (four states under Reasenberg), motivated
  explicitly by the observation that declustered background does not conform to a single-rate
  Poisson model. K-009's "is there a latent state?" has a published, region-specific answer with
  a state count attached.
- **Werner, Ide & Sornette (2011), Nonlin. Processes Geophys. 18, 49–70 — "Earthquake forecasting
  based on data assimilation: sequential Monte Carlo methods for renewal point processes."** The
  first implementable DA scheme for point-process seismicity forecasting. K-010 Tier 2 is a
  variant of this and must cite it.
- **Deep Gaussian process background rates (Molkenthin et al. / Zhu et al., 2023, GJI 234,
  427–xxx, doi:10.1093/gji/ggad074)** and graphical-Dirichlet-process inhomogeneous background
  intensities (Math. Geosci. 2026) — the "background field is a smooth latent random field"
  program is active and well past the diagnostic stage.

### M-001.3 CONTRADICTED? — No. CONTESTED? — Yes, on two axes.

**Nothing credibly shows the opposite.** I looked for a paper asserting that fitted-ETAS
residuals in California are adequately white and found none; the closest is Gardner & Knopoff
(1974), which Luen & Stark (2012) explicitly attribute to low test power. **W-003-P2's "the
residuals are white" is therefore a generator with, as far as I can find, no prior-art support
in this region.** That is worth saying plainly: our own incumbent theory made the prediction the
literature had already falsified three times. Popper should weigh that when scoring the
two-generator discrimination — generator A was not a 50/50 prior even before we ran.

**Contested axis 1 — magnitude range.** Zaliapin & Ben-Zion (2020) do *not* reject stationarity
for Δm < 4. If our excess is carried by the M≥6 sequence blocks and their imperfectly-modelled
aftershock fields, "there is a latent field" and "ETAS under-models big sequences" are not
distinguished by any statistic in `results_k009.json`.

**Contested axis 2 — misfit vs state.** Bray/Schoenberg's documented ETAS spatial misfit and
Hainzl's rate-dependent incompleteness are both live, published, non-latent-field generators for
our Moran's I and ACF excesses. Our kernel swap and Mc-proxy partial address them; the
literature says they are the right things to address and does not say our controls are
sufficient.

**Adjacent negative result worth recording.** Sirorattanakul & Avouac (2026), *Science Advances*
12(13), eadz5711 — "Seismic rhythms" — find annual hydrologic modulation of background rate up to
**15% in the northern SAF** but **minimal modulation in Southern California**, and no significant
semidiurnal tidal modulation, with peak seismicity lagging peak stress by a median 0.52 months.
This *helps* us: a seasonal hydrologic driver is unlikely to be the source of our SoCal red
residual. It also supplies an independent ~2-week nucleation-response timescale that sits
suspiciously close to our 2.4-week integral correlation time.

### M-001.4 What OUR version actually adds (the honest delta list)

Ranked by how well it survives an adversarial reviewer.

- **(a) Pre-registered two-generator discrimination with a zero-free-parameter alternative
  value — DEFENSIBLE AND, AS FAR AS I CAN FIND, UNPRECEDENTED IN THIS LITERATURE.**
  `results_k009_prediction.json` was written before scoring and commits generator B to a numeric
  correlation time derived from an *independently measured* `t_a` (193 d stacked;
  full honest estimator range 200 d – 5 yr), with an acceptance band and an explicit
  "consistent but not sharply discriminating" clause. I found no prior residual-analysis or
  background-rate paper that pre-registers a predicted correlation time from an independent
  physical estimate and then scores against it. Every prior work above is exploratory or
  descriptive at this step. **This is the strongest thing we own.**
- **(b) A null simulated from the same frozen fitted spatio-temporal ETAS — a genuine
  methodological improvement, but with a named precedent to cite.** Prior work compares to
  analytic Poisson (Luen & Stark), to randomized-reshuffled catalogs preserving the spatial
  distribution (Zaliapin & Ben-Zion), to uniformity of transformed times (Ogata; KS/Ljung-Box),
  or to the model's own second-order structure (Zhuang 2006). Simulating the fitted generator
  and pushing the sims through the identical pipeline absorbs ETAS-manufactured correlation and
  all estimator pathologies at once. **Nearest prior: Zhuang (2006).** Claim "sim-calibrated" as
  a design choice, not an invention.
- **(c) SoCal-wide gridded field rather than a case study.** Llenos & McGuire: Salton Trough.
  Hainzl & Ogata: Vogtland. Kumazawa & Ogata: post-Tohoku inland Japan. Ogata's HIST-ETAS/
  Delaunay maps are the real precedent for a region-wide field and they predate us by two
  decades — so this is *scale*, not *kind*. Modest delta.
- **(d) The completeness partial — an ADOPTED CONTROL, NOT A NOVELTY.** Hainzl, Zöller & Wang
  (2013) showed short-term Mc variation biases ETAS parameters (notably α); Hainzl (2016, BSSA)
  introduced rate-dependent detection probability; Hainzl (2022) the **ETASI** closed form.
  Presenting "surviving a completeness-proxy partial" as a distinguishing rigour would be
  incorrect — it is standard practice, and our version is a **surrogate** (the file's own FLAG:
  the K-031 ρ_sta field is not on disk; we used an Mc maximum-curvature proxy on 0.6° super-cells
  which explains only R²=0.047 of EOF1). The partial as run is weaker than the field's norm.
- **(e) Not a delta:** "the residual is the innovation / if white there is no state to estimate."
  That framing is the standing premise of Werner, Ide & Sornette (2011) and of the entire
  Llenos–McGuire assimilation line. It is a good framing. It is theirs.

### M-001.5 TAKEABLES (adopt instead of reinventing)

**Estimators**
1. **Kumazawa–Ogata spline background + ABIC** for μ_t — a validated, published, hyperparameter-
   objective alternative to K-010 Tier 1's EWMA. Adopt or justify not adopting.
2. **Marsan et al. (2013, BSSA) alternating MLE-ETAS / forcing-rate update** — a second,
   independent, validated latent-forcing estimator. Two estimators agreeing is worth more than
   one estimator with a null.
3. **Zhuang (2006) second-order residuals** and **Bray et al. (2014) Voronoi residuals** — both
   were built *because* rectangular-grid residual diagnostics have low power. Our T and L
   estimators saturating at their caps in 20/20 and 7/20 null sims is that low-power problem
   arriving on schedule. This is the single most actionable methodological takeable.
4. **ETASI (Hainzl 2022)** — replaces the Mc-proxy partial with a model that absorbs
   rate-dependent incompleteness inside the intensity, rather than regressing it out afterwards.
5. **Zaliapin & Ben-Zion randomized-reshuffled null** (stationary, space–time independent,
   preserves the spatial distribution) — a cheap *second* null that brackets the ETAS-sim null
   from the other side.
6. **Werner–Ide–Sornette SMC** for K-010 Tier 2; **deep-GP / Dirichlet-process background fields**
   (GJI 2023; Math. Geosci. 2026) if the latent field is to be estimated rather than detected.
7. **MMPP with 3 states** as the pre-fitted discrete-state alternative to K-010's OU state in
   SoCal (Earth Sci. Inf. 2025). Gives a state count we would otherwise have to search for.

**Datasets we can use tomorrow, already labelled**
8. **Ross & Cochran (2021) catalog of 92 long-duration SoCal swarms, 2008–2020, durations 0.5–7
   yr, 53% with diffusive backfronts.** Overlaps 2010–2018 completely. This is a *ground-truth
   label set for the latent field*. Recommended as a K-009 add-on control (see M-001.6).
9. **Salton Trough M≥1.5 1990–2009** (Llenos & McGuire's exact dataset) — a region and window
   where the answer is independently known, i.e. a real-data positive control to complement our
   synthetic OU injection.
10. **Cahuilla swarm 2016–2019** (Ross et al. 2020, Science) — a single, well-characterised,
    in-window, in-region transient with a known duration and a known ~few-km spatial scale. If
    our residual field cannot see Cahuilla, our 0.2°/7 d resolution is the answer, not the crust.

**Named artifacts prior work already paid for**
11. ETAS spatial-kernel misfit: on-fault under-prediction / off-fault over-prediction (Bray et
    al. 2014). Our kernel swap must be reported against this specific pattern, not generically.
12. Short-term aftershock incompleteness biases α downward (Hainzl et al. 2013) — our frozen
    α = 0.537 is low; that is a documented signature of exactly this bias and it propagates into
    every expected count in K-009.
13. Declustering/estimator sensitivity of the whiteness verdict (Luen & Stark 2012): the answer
    depends on method, catalogue, magnitude range and test. Report all four.
14. Magnitude-range dependence of the stationarity rejection, Δm ≷ 4 (Zaliapin & Ben-Zion 2020).
15. Low power of grid-based point-process diagnostics (N-test/L-test family) — Bray et al. 2014.

**Parameter values worth comparing ours against**
16. Hainzl & Ogata (2005): only a *few percent* of Vogtland swarm activity directly attributable
    to the fluid signal. A useful prior on the amplitude of the latent forcing we should expect —
    and a reason our positive control failing at log_sd = 1.5 matters more than it looks.
17. Sirorattanakul & Avouac (2026): ~15% annual modulation amplitude (N. California), ~0.52-month
    stress-to-seismicity lag, minimal SoCal seasonal signal.
18. Meade & Hager (2005) moment-rate numbers — see M-002.

### M-001.6 FEEDING THE TRIO — what this changes

**For Popper (the pending K-009 ruling).**
1. **The replication burden on "not white" should fall; the burden on the two numbers should
   rise.** "Not white" is independently established in SoCal by Luen & Stark (2012) and Zaliapin
   & Ben-Zion (2020). We are confirming, not discovering — so demanding a second region before
   believing the sign of the effect is over-strict. But K-009's *stated deliverable* is "the
   correlation time and correlation length with CIs — those two numbers, not a p-value," and
   those two numbers are the part of this run that does not survive scrutiny: T is pinned to the
   estimator cap in the real data *and in 20/20 null sims*; L = 41.0 km sits at the null's 2.5th
   percentile with the null median at 78.7 km. **On the ledger's own success criterion the
   defensible K-009 result today is "lag-1 ACF excess +0.093, Moran's I and EOF1 excess,
   integral correlation time 2.4 weeks" — and NOT "months-scale persistence with 40 km
   coherence."** The prior art makes this worse rather than better: Luen & Stark's central
   methodological finding is that whiteness verdicts are estimator-dependent, and Bray et al.
   and Zhuang both built better residual estimators precisely because grid statistics lose power.
2. **The two-generator scoring is affected by an asymmetric prior.** Generator A (W-003, "white")
   had already been falsified in this region three times before we tested it. A W-003 death from
   K-009 is therefore worth less than the ledger's R2-7 table implies, and the interesting
   discrimination is not A-vs-B but **B (aseismic forcing, T ≈ t_a) vs C (ETAS misspecification)**
   — which the ledger already names as the third branch. Note that 2.4 weeks is far from t_a's
   193 d–5 yr band, so on the integral estimator generator B is *not* currently confirmed either.
3. **A cheap third-arm control now exists and I recommend Popper mandate it.** Ross & Cochran's
   **92 labelled long-duration SoCal swarms, 2008–2020** overlap our window entirely. Test: does
   the residual field's red mode co-locate in space and time with those 92? Co-location ⇒ B
   (aseismic transient forcing), and it is a far sharper discriminator than the t_a band, which
   spans an order of magnitude by the prediction file's own admission. No co-location, but excess
   concentrated on-fault / around large sequences ⇒ C (misspecification, per Bray et al.). Cost:
   one join against a published catalog. This is the highest-value hour in the queue.
4. **Δm control.** Re-run restricted to M2.5–4.5 (Δm < 2) per Zaliapin & Ben-Zion's threshold. If
   the excess vanishes, the "latent field" is the large-sequence residual and the assimilation
   case collapses to "fit aftershocks better."
5. **Nothing in the prior art relieves** `n_sims = 20` against a spec of 500, or the positive
   control passing only at log_sd = 3.0 and failing at 1.5 — the latter reads worse in light of
   Hainzl & Ogata's finding that real fluid forcing contributes only a few percent of activity.

**For Wegener.** Two candidate observation rows are his to write, not mine: (i) Ross & Cochran
(2021) 92 long-duration SoCal swarms, 6 mo – 7 yr, 53% diffusive, aseismic driving "active at
all times" 2008–2020 — this is a direct, replicated, documented observation of W-001/W-002's
driver in our region; (ii) Zaliapin & Ben-Zion (2020) Δm-dependent stationarity rejection. Note
also that **W-003-P2's whiteness prediction is contradicted by prior art independent of our
test** — the incumbent was already in trouble.

**For Kepler.** The floor moved. "Is there a latent state?" is answered (yes, several ways, by
several groups). The unclaimed ground is: *what is the latent field's amplitude, its estimator,
and its forecast value in bits/event* — i.e. K-010 and K-012, not K-009. And K-010 should start
from Kumazawa–Ogata + Marsan, not from a fresh EWMA.

---

## M-002 — QUICK DOSSIER: B-4 / EXP-J — "silent-loading cells recover the creeping SAF and the 1857 strand"

### CLASSIFICATION: **REDISCOVERY (of the geography) — and the program record already says so.**
### The degeneracy at its heart is **CONTESTED**, and prior art has already produced a discriminator.

**What we have.** `results_exp_j.json` / `results_exp_k.json`: 1,195 cells at 0.2°, 228
SILENT-LOADING; 175/229 interior (18.5% interior vs 20.9% edge — not a border artifact); top
interior silent cells are the SAF creeping section (36.1–36.9 N, −120.8..−121.8) plus
Imperial/Brawley; top *unexplained*-silent cells sit on the 1857 Fort Tejon Mojave/Big Bend
strand (34.4–34.8 N, −117.5..−119.1, d_fault 1–5 km). EQ18_FULL_NOTES §15 already calls this
"the ledger's negative space rediscovers known aseismic zones blind = **validation**." That
sentence is the correct classification and it should never drift into a discovery claim.

**The shoulders, and they are crowded.**

1. **Meade & Hager (2005), JGR 110, B03403 — "Spatial localization of moment deficits in
   southern California."** GPS-constrained block model, fault slip-rate catalog; SoCal scalar
   moment accumulation **17.8 ± 1.1 × 10¹⁸ N m/yr, ≈50% larger than the 200-yr average release
   rate**, with deficits **localized in three regions: the southern SAF and San Jacinto, the
   offshore faults + LA/Ventura basins, and the Eastern California Shear Zone**. This is B-4's
   experiment — geodetic loading minus seismic release, mapped, in SoCal — done twenty-one years
   earlier with a better loading model. Our delta is grid resolution (0.2° cells vs blocks) and a
   frozen pre-2010 train, not the concept and not the headline geography.
2. **Kostrov (1974), Izv. Earth Phys. 1, 23–44** — the summation the whole ledger rests on; and
   **Ward (1998), GJI 134, 172–186 / Ward (1994)** for California Kostrov budgets. **Guns,
   Bennett et al. (2024), JGR 129, e2023JB027939** — "Seismic Moment Accumulation Rate From
   Geodesy: Constraining Kostrov Thickness in Southern California" — the *current* state of the
   art on exactly our quantity, and the paper to check our χ normalisation against (recall EQ18
   §15: "χ absolute values biased low ... geography is the signal, not the level" — Guns et al.
   is where the level lives).
3. **Field et al. (2014), BSSA 104, 1122–1180 (UCERF3)**; **WGCEP (1988, 1995)** — the 1857
   Mojave/Big Bend section is the most-forecast locked, late-in-cycle strand in California, with
   a paleoseismic recurrence record (Weldon et al., 2004, GSA Today) to match. "The ledger found
   the 1857 strand blind" is a statement about the ledger's sensitivity, not about the Earth.
4. **Jolivet et al. (2015), GRL 42, 297–306 — aseismic slip and seismogenic coupling along the
   central SAF from InSAR**; **Tong, Sandwell & Smith-Konter (2013), JGR**; **Ryder & Bürgmann
   (2008), GJI 175, 837 — spatial variations in slip deficit on the central SAF from InSAR**;
   **Maurer & Johnson**. Coupling maps of the creeping section exist at far higher resolution
   than a 0.2° catalog-derived ledger, from data that is not degenerate the way ours is.
5. **Liu, Ross, Cochran & Lapusta (2022), Science Advances 8, eabk1167 — "A unified perspective
   of seismicity and fault coupling along the San Andreas Fault."** *This is the important one
   for us.* They show **creep rate along the central SAF is directly proportional to the fraction
   of non-clustered earthquakes, 1984–2020**: lower coupling ⇒ weaker temporal clustering, with
   repeating earthquakes as the end-member.

**Delta assessment.** Small and methodological: a 0.2°-resolution, catalog-Kostrov-vs-geodesy
ledger, trained pre-2010, whose *negative space* was read as the signal and which recovered known
geology without being told. That is a legitimate instrument-validation claim and a nice framing.
It is not a new observation about California, and the ledger's §10 SCOPE line plus the
auto-FLAG (158/200 unexplained-silent cells with n_train < 20) already constrain it correctly.

**CONTESTED — and W-006-P1 is in immediate trouble from prior art.** W-006-P1(a) predicts *no
catalogue-derived statistic separates creeping from locked silent cells at better than chance*.
**Liu et al. (2022) already published a catalogue-derived statistic that does exactly that** on
the central SAF — the non-clustered fraction, proportional to creep rate over 36 years. Wegener's
own honest note anticipated the shape of this ("repeater fraction is arguably a catalogue-derived
aseismic readout"); the non-clustered fraction is a cleaner instance and it is not hypothetical.
**Recommendation to Popper:** W-006-P1(a) as written is close to pre-falsified in the literature
and should be re-scoped to SoCal cells off the central SAF, or reframed as "the *ledger's own* χ
does not separate them, but a clustering statistic does" — which is a better and still-testable
claim, and one that hands B-4 a free upgrade rather than killing it.

**Takeables.** (i) Guns et al. (2024) for Kostrov thickness / χ level calibration; (ii) Meade &
Hager (2005) three-region deficit map as a *direct comparison target* — Jaccard our silent list
against their three localizations and report it, because a reviewer will; (iii) Liu et al. (2022)
non-clustered fraction as the ready-made degeneracy discriminator, computable from our own
catalog at zero data cost; (iv) Jolivet/Ryder-Bürgmann InSAR coupling as the geodetic covariate
W-006-P1(b) asks for, already published, no new inversion needed.

---

## M-003 — QUICK DOSSIER: B-3 / EXP-I(ii) — "P(M≥5 within 7 d | M≥5) = 0.60"

### CLASSIFICATION: **REDISCOVERY — canonical, textbook, operational. Zero novelty. Cite or retract.**

This is Omori–Utsu decay plus Gutenberg–Richter, evaluated in a 7-day window. It has been the
operational basis of public earthquake advisories in California for thirty-seven years. The
ledger's §10 SCOPE line already attaches the right caveat ("ordinary within-sequence triggering,
not exotic meta-structure"); this dossier supplies the citations so the caveat has a bibliography.

**The canonical citations — attach these to B-3 permanently.**

1. **Reasenberg & Jones (1989), Science 243, 1173–1176 — "Earthquake hazard after a mainshock in
   California."** *The* generic California clustering model: modified Omori + G-R, giving the
   probability of further events (including larger ones) in intervals following any earthquake.
   Extended in Reasenberg & Jones (1994), Science 265, 1251. **This is B-3.**
2. **Jones (1985), BSSA 75, 1669–1680 — "Foreshocks and time-dependent earthquake hazard
   assessment in southern California."** P(an M≥3 SoCal event is followed by a *larger* event
   within 5 d and 10 km) = **6 ± 0.5%**, rising with foreshock magnitude to **6.5 ± 2.5% at
   M≥5**; the mainshock most likely in the first hour, decaying as ~1/t.
   **Read this against our 0.60 carefully.** Jones's 6.5% is P(*larger* event); ours is P(*any*
   M≥5, including equal-or-smaller aftershocks of the same sequence). They are different
   statistics and the numbers are not in conflict — but they will be *read* as in conflict if we
   quote 0.60 without the definition. Anyone reporting B-3 must state "any M≥5, same sequence
   included" in the same sentence as the number.
3. **Gerstenberger, Wiemer, Jones & Reasenberg (2005), Nature 435, 328–331 — "Real-time forecasts
   of tomorrow's earthquakes in California" (STEP).** The operational short-term clustering
   forecast; the direct ancestor of every 7-day hazard window.
4. **Ogata (1988), JASA 83, 9–27 — ETAS.** Our own B-2 baseline already contains this rule; B-3
   is a marginal read-out of the same generator, which is why the two must never be quoted as
   independent findings.
5. **Michael et al. (2020), SRL 91, 153–173 — the USGS operational aftershock forecast (OAF)
   system**; **Page, van der Elst, Hardebeck, Felzer & Michael (2016), BSSA 106, 2290–2301 —
   "Three ingredients for improved global aftershock forecasts"**; **Hardebeck et al., updated
   California aftershock parameters (USGS)**. Current practice, with published parameter values
   for the exact windows we used.

**Delta:** none. Not scale, not rigour, not region. Our n_test = 30 is far smaller than any of
the above; the p = 7×10⁻¹⁵ measures how strongly we rejected a Poisson straw man that no one in
this literature has believed since 1989. **Mandate: B-3 is presented only as a local
re-measurement of the Reasenberg–Jones effect, with Reasenberg & Jones (1989) and Jones (1985)
cited in the same breath, and with the "any M≥5" definition attached.** Its legitimate value is
as an app layer (EQ-15 layer 2), not as a finding.

---

*End Merton round 1. Counts: 1 primary + 2 quick dossiers. REDISCOVERY 3 (K-009, B-4, B-3),
NOVEL 0, CONTRADICTED 0, CONTESTED 2 (K-009 interpretation: magnitude-range and misfit-vs-state;
B-4/W-006-P1 degeneracy, where prior art supplies a discriminator W-006 predicts cannot exist).
The one thing in this round we genuinely own is the pre-registered two-generator discrimination
with an independently-measured predicted correlation time — I found no precedent for it. The
phenomenon it discriminates over was found before us, in our region, at least three times.
Sources for every claim above are named with venue and year; where I could verify only an
abstract or a publisher's summary rather than the full text — Zaliapin & Ben-Zion (2020),
Ross & Cochran (2021), the 2025 MMPP paper, Guns et al. (2024) — I say so here rather than let a
claim stand unmarked.*
