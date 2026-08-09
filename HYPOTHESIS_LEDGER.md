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
*End Kepler round 1 (K-001..K-032). All PROPOSED. Nothing above is claimed as true; nothing
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
