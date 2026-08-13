# MINING CATALOG (Kepler, PROPOSED) - candidate features, properties and combinations for EQ-24 v2

**Status: PROPOSED. Pure generation. No result is claimed anywhere in this file.** Nothing here is
evidence, nothing here is a finding, and no entry may be entered for or against any ledger entry until
Popper prices it and G-M1 clears at the relevant band and aggregation level (§P6-4 banner item 5).
This file is Kepler's alone this round; a K-entry registering it comes after Popper sees it. I generate,
Popper prices. Every price quoted below is my own arithmetic offered as an input to his pricing, never
as a substitute for it.

Written 2026-08-12 against: §P5 blind-spot rulings, §K87-0(a)-(d), §P6-4 (Rules 4.1-4.7), S-13, S-14,
S-15, S-16, S-17 (candidate), the F-012 bounds paper, K-001..K-091, `engine/mine.py`
(`ephemeris_features`, `catalog_features`, `download_features`), `engine/SPEC.md`, `engine/README.md`.

---

## 0. HOW TO READ AN ENTRY

```
**ID name** | TRANCHE-READINESS
Def:   definition precise enough to implement without asking me a question
Src:   where the data is (on disk / free download with URL / infrastructure to build)
Hunch: the physical or statistical reason, one sentence, honest about how speculative
Pit:   pitfalls: aliasing, sinc zeros, catalog-composition drift, collinearity, corpses to dodge
S15:   measurability at M>=4.5 global daily binning, per the arithmetic in §0.3
Price: my estimate of the priced test count (features x lags x marks x rungs)
```

**Tranche-readiness values.** `READY` = data on disk or closed-form, implementable this week.
`DOWNLOAD` = free source named with URL, one fetch-and-parse away. `INFRASTRUCTURE` = needs a build,
and the build is named.

### 0.1 The five hard constraints every entry in this file respects

1. **Exploration window only** (`[365, 8081)`, 7716 days). No holdout hash may be spent from mine mode.
   Any region or cell selector is derived from exploration-window data only (Rule 4.1). **The K-080
   census cell list is never a selector here** - it is holdout-contaminated and it is on a 0.5 deg grid
   the engine does not use. Where an entry needs a spatial partition, the partition rule is stated and
   is computable from `[365, 8081)` alone.
2. **No per-cell battery at M>=4.5.** §P6-4 Finding B is not negotiable: 4,329 active cells at a mean
   of ~16 events is power-absent, not power-limited. Every spatially-resolved entry below routes
   through the **2R-df phase-incoherent regional sum** (Rule 4.2) or is explicitly marked
   EXPLORATORY-UNPRICED per Rule 4.4.
3. **The tidal bound is spent.** Globally-coherent fortnightly-to-annual rate modulation is bounded at
   ~5.6% (miner) and the F-012 paper bounds tidal phase modulation at ~6.6%. **No entry below
   re-proposes a static tidal-phase map or a plain fortnightly phase test.** Tidal entries appear only
   where the bound does not reach: (a) sub-daily bands, which the day-binning notches out exactly and
   which therefore need the unbinned path (INFRASTRUCTURE); (b) transient-conditioned response, which
   S-14(c) brackets; (c) regionally-incoherent phase, which the domain sum cancels before any test runs.
4. **Second moments are a separate instrument.** §K87-0(c): every statistic currently in the session is
   a first-moment single-phase statistic. Entries in family F9 are the ones that address this, and they
   are the ones I would fight for.
5. **Prefer the boring model.** Mignan & Broccardo stands: any feature whose skill a 2-parameter GLM
   cannot reproduce after conditioning is treated as leakage until proven otherwise.

### 0.2 Corpses this catalog explicitly dodges (checked entry by entry)

static tidal-phase maps; calendar-time periodicity scans as a standalone claim; Fibonacci and other
ratio numerology; fault-type parameter pools; sequence-shape spatial transfer; the TSI ratio;
per-cell M>=4.5 batteries; classical quiescence as a standalone precursor observable (M-007: owned and
prospectively decided - quiescence appears here only as a covariate or as a second-moment object);
isotropic-vs-anisotropic ETAS kernel comparison as a novelty claim (M-008: the mechanism-free arm is
occupied); California seasonal hydrological loading as a standalone claim (Sirorattanakul & Avouac own
it - appears here only as a global-extension or as a conditioning variable, and says so); site
dilatation quoted to better than a factor 2 (M-006); `exp(dtau/Asigma)` used as a power calculator
outside 7-347 kPa (§P5-8, VOID).

### 0.3 THE MEASURABILITY ARITHMETIC I AM PROPOSING (and its one honest gap)

S-15 requires a power floor in the units of the statistic, declared before the run. Rather than
re-derive it 300 times, I propose one rule and use it everywhere below, flagged as PROPOSED:

- **Raw-count floor (Popper, §P6-4 Re-derivation 1):** `A_min = 3.96 / sqrt(N)` for a sinusoidal rate
  modulation on `N` events.
- **Clustering-deflated floor (mine, PROPOSED):** the miner's own reported global bound is ~5.6% at
  N = 46,585, and `3.96/sqrt(46585) = 1.8%`. The factor of ~3 in amplitude is a factor of ~9 in
  variance, which is exactly the size of an aftershock-clustering overdispersion. So I propose
  `N_eff ~ N/9` and `A_min ~ 11.9/sqrt(N)` as the working floor for any daily-binned count statistic
  on the global M>=4.5 catalog.
- **The gap, stated rather than hidden:** I have not verified that the factor is clustering rather
  than 2-df-vs-1-df bookkeeping, BH overhead, or the surrogate floor. **This is a measurement someone
  should make and it is cheap** - it is entry F4-58 below. Until it is made, every S15 line here quotes
  both floors, and Popper should treat the deflated one as the operative number.

Reference points under the deflated rule (`A_min ~ 11.9/sqrt(N)`):

| N events | A_min deflated | A_min raw | verdict at a plausible ~5-10% effect |
|---|---|---|---|
| 46,585 (global, all) | 5.5% | 1.8% | MEASURABLE |
| 20,000 | 8.4% | 2.8% | MARGINAL |
| 10,000 | 11.9% | 4.0% | MARGINAL to UNMEASURABLE |
| 3,500 (R=20 region) | 20.1% | 6.7% | UNMEASURABLE deflated / MARGINAL raw |
| 1,000 | 37.6% | 12.5% | UNMEASURABLE |
| 16 (mean cell) | absurd | 99% | POWER-ABSENT (§P6-4 Finding B) |

**Consequence I want read out loud:** under the deflated rule, even the R=20 regional battery is at the
edge. That is not an argument against the 2R-df instrument (which keeps the full event count and is a
detection statistic, not a per-region estimator) - it is an argument that **per-region amplitude
quotes are UNRESOLVED and only the summed statistic may be quoted**, which is a stricter reading of
Rule 4.3 than Rule 4.3 itself states. I flag it because if I am right, a per-region table of
amplitudes is the K-076 Colombia cell all over again.

### 0.4 Price conventions used below

`phase` features cost 2 df, `linear` features cost 1 df. Lag scans are 0..30 d inclusive = 31 values
unless stated. Lag is free only where theta is linear in t (the four `LAG_FREE_PHASE_FEATURES`);
assume PRICED unless an entry claims otherwise and names the invariance. Mark tests cost
`n_features x n_marks`. Ladder rungs are never free. Where an entry says `Price: 31`, read
"31 priced tests at the S-8 max-statistic, before any BH stratification".

---

# FAMILY 1 (extended) - SOLAR-SYSTEM GEOMETRY BEYOND THE CURRENT NINE

The engine currently carries 9 family-1 features, all geocentric, all evaluated at 12:00 UTC.
§K87-0(d)(iv) is the standing hole: **no local tidal stress exists anywhere in the engine**. Entries
F1-30..F1-40 are the ones that fill it, and they are INFRASTRUCTURE.

**F1-01 moon_libration_lon** | READY
Def: optical libration in selenographic longitude (deg), Meeus ch. 53, daily at 12:00 UTC.
Src: closed form, extend `engine/ephemeris.py`. No download.
Hunch: SPECULATIVE and mechanism-free - libration has no plausible stress channel at Earth; it is
entered deliberately as a **negative-control cycle** with the right period family, so a detection here
diagnoses the pipeline rather than the Earth.
Pit: 27.55 d / 27.21 d beat against `moon_anomalistic_phase` and `moon_draconic_phase`; report the
collinearity (min R^2 against the existing lag-0 column space) or the "new" test is a relabelled old one.
S15: N = 46,585; A_min 5.5% deflated / 1.8% raw. MEASURABLE for a bound.
Price: 31 (1 linear x 31 lags).

**F1-02 moon_libration_lat** | READY
Def: optical libration in selenographic latitude (deg), Meeus ch. 53.
Src: closed form in `ephemeris.py`.
Hunch: SPECULATIVE, mechanism-free, same negative-control role as F1-01, different period (draconic).
Pit: near-degenerate with `moon_declination` at the 27.2 d line; must be entered as a residual against it.
S15: as F1-01. MEASURABLE for a bound.
Price: 31.

**F1-03 moon_parallax** | READY
Def: equatorial horizontal parallax of the Moon (arcsec) = asin(a_E / r).
Src: closed form; monotone in `moon_distance` already present.
Hunch: NOT independent - this is a reparameterisation of distance and I say so; its only value is as a
**nonlinearity probe**, since the tidal potential goes as r^-3 and a linear-in-parallax response versus
a linear-in-r^-3 response are distinguishable.
Pit: pure corpse risk if entered as a new feature; enter it only as the pair difference against
`tidal_potential_proxy`, else it is a duplicate test paying full multiplicity.
S15: MEASURABLE for a bound.
Price: 31, and I recommend 0 unless the nonlinearity framing is adopted.

**F1-04 moon_ecliptic_latitude** | READY
Def: geocentric lunar ecliptic latitude beta (deg), Meeus ch. 47.
Src: closed form.
Hunch: the driver behind draconic phase in an unwrapped form; a response that saturates at |beta| max
looks different in beta than in draconic phase, and the current set only carries the phase.
Pit: |beta| and beta must both be entered (a symmetric response is invisible in the signed one);
collinear with `moon_declination` at the 27.2 d line.
S15: MEASURABLE for a bound.
Price: 62 (2 linear x 31 lags).

**F1-05 lunar_nodal_phase_18_6y** | READY
Def: phase of the regression of the lunar ascending node, period 6798.4 d (18.61 y), as `kind='phase'`.
Src: closed form.
Hunch: the nodal cycle modulates the declination-tide envelope by ~4% and is the longest deterministic
tidal envelope; if any secular tidal modulation exists, this is where it lives.
Pit: **the exploration window is 7716 d = 1.13 nodal cycles.** This is a one-cycle claim and the block
bootstrap will clip at 800 d giving ~10 blocks. I would price it but I would refuse to quote a bound
from it; it is the clearest UNRESOLVED-by-window entry in the catalog.
S15: UNMEASURABLE by window length, not by count. Declare it as such and score it neither way.
Price: 1 (lag scan meaningless at this period).
Corpse-dodge: this is not the dead calendar-periodicity scan; the period is fixed a priori, not searched.

**F1-06 nodal_envelope_amplitude** | READY
Def: the slowly-varying amplitude multiplier `1 + 0.037*cos(nodal_phase)` applied to
`moon_abs_declination`, entered as a linear covariate.
Src: closed form.
Hunch: a **product** feature, which is the honest form of the nodal hypothesis: the nodal cycle should
not modulate rate directly, it should modulate the depth of the 13.6 d declination modulation.
Pit: this is a second-order effect on a first-order effect that is already bounded at 5.6%; the implied
effect is <0.25% and is **below any floor in the table in §0.3**. Say so before running it.
S15: UNMEASURABLE at M>=4.5 global. I include it so nobody proposes it later thinking it is live.
Price: 31, and my recommendation is 0.

**F1-07 apsidal_precession_phase_8_85y** | READY
Def: phase of lunar perigee precession, period 3232.6 d (8.85 y), `kind='phase'`.
Src: closed form.
Hunch: mechanism-free at first order, but it is the beat partner that generates the 411.8 d
perigean-spring cycle already in family 2, so its absence is an asymmetry in the current set.
Pit: 2.4 cycles in the window; block length clips at 800 d; heavily collinear with
`perigean_spring_beat` by construction.
S15: MARGINAL by window. Bound only.
Price: 1.

**F1-08 sun_earth_distance** | READY
Def: heliocentric distance of Earth (AU), Meeus ch. 25.
Src: closed form.
Hunch: the solar tidal amplitude term; the engine carries `sun_declination` but not the solar r^-3
term, so the solar tidal potential is currently incomplete in a way the lunar one is not.
Pit: annual period exactly - collinear with `annual_phase` and with every seasonal environmental
feature in family F6. **Aliases with the seasonal hydrology corpse-adjacent claim** and must be entered
jointly with F6 features or its interpretation is unavailable.
S15: MEASURABLE for a bound; interpretation UNRESOLVED without F6 conditioning.
Price: 31.

**F1-09 solar_tidal_potential_proxy** | READY
Def: `(r0/r_sun)^3 * (0.5 - 1.5*sin^2(sun_dec))`, the exact solar analogue of the existing
`tidal_potential_proxy`.
Src: closed form.
Hunch: the current set has a zonal lunar potential proxy and no solar one; the solar term is ~46% of
the lunar and its omission is a plain gap, not a hypothesis.
Pit: annual collinearity as F1-08.
S15: MEASURABLE for a bound.
Price: 31.

**F1-10 total_zonal_potential** | READY
Def: sum of lunar and solar zonal potential proxies (F1-09 plus the existing feature).
Src: closed form.
Hunch: the physically correct combined quantity; the parts are already tested, the sum is not, and a
linear-response system responds to the sum.
Pit: exactly a linear combination of two existing columns - **free-ish in information, full price in
multiplicity**. Popper should decide whether to charge for it; I would enter it and drop F1-09.
S15: MEASURABLE.
Price: 31.

**F1-11 tidal_potential_rate** | READY
Def: analytic d/dt of F1-10, per day.
Src: closed form.
Hunch: K-087's rate-versus-level axis applied to the combined potential; a rate-regime response and a
level-regime response are different physics and the 2-df omnibus refuses to name which it saw.
Pit: the derivative of a near-sinusoid is a 90 deg rotation, so for `kind='phase'` features this is
**free and already covered**; it is only informative for `kind='linear'` features. Do not pay for the
phase versions.
S15: MEASURABLE.
Price: 31.

**F1-12 tidal_potential_accel** | READY
Def: second analytic derivative of F1-10.
Src: closed form.
Hunch: if the response is a threshold-crossing rather than a level, the crossing rate depends on the
second derivative near the extremum; this is the cheapest probe of a threshold response.
Pit: differentiation amplifies the highest frequency present, which for a day-sampled series is the
band the sinc notch already killed; expect this to be dominated by numerical noise unless computed
analytically.
S15: MEASURABLE for a bound only.
Price: 31.

**F1-13 syzygy_proximity** | READY
Def: `min(|elongation|, 180 - |elongation|)` in deg, the unsigned distance to the nearest syzygy.
Src: closed form from the existing `elongation_deg`.
Hunch: folds new and full moon together, which the signed elongation does not; if the response is
symmetric in syzygy type the current feature halves the effective count.
Pit: this is a **known corpse-adjacent shape** (the classic new/full-moon claim); it must be entered
knowing the 6.6% bound already applies and framed as a bound-tightener, not a discovery.
S15: MEASURABLE.
Price: 31.

**F1-14 perigee_proximity_days** | READY
Def: signed days to the nearest lunar perigee.
Src: closed form.
Hunch: a **time-to-event clock** rather than a phase; if the response is a sharp window around perigee
rather than a sinusoid, a sinusoidal basis loses most of its power, and this feature does not.
Pit: strictly a monotone reparameterisation of anomalistic phase within a half-cycle; the gain is
entirely in the basis shape, so the honest test is a harmonic-ladder comparison, not a new p-value.
S15: MEASURABLE.
Price: 31.

**F1-15 syzygy_perigee_coincidence_window** | READY
Def: binary/soft indicator that a perigee and a syzygy fall within +/- 1.5 d of each other, smoothed.
Src: closed form.
Hunch: the perigean spring tide is the largest routine tidal excursion; a rectangular-window feature
tests the "only the big ones matter" hypothesis that a sinusoid cannot express.
Pit: rare-event feature - roughly 8-10 coincidences per year, so the effective N supporting it is a
small fraction of 46,585; this is the entry most likely to look significant for count reasons.
S15: effective N ~ 3,000-5,000 event-days. A_min ~17-22% deflated. **MARGINAL to UNMEASURABLE.**
Price: 31.

**F1-16 draconic_syzygy_eclipse_proximity** | READY
Def: days to nearest solar or lunar eclipse (true eclipse ephemeris, not the 173.3 d beat proxy).
Src: closed form (eclipse conditions from the existing three phases) or a static table.
Hunch: mechanism-free relative to the existing `eclipse_year_beat`, and I say so; its only distinct
content is that true eclipses are a **conjunction of three conditions**, which is a sparser and
sharper marker than the smooth beat.
Pit: ~4-7 eclipses/yr means very low effective N; and a positive here with a null on the smooth beat
would be a strong hint of a counting artifact rather than physics.
S15: effective N ~1,500. UNMEASURABLE except as a bound.
Price: 31.

**F1-17 moon_azimuth_at_domain_centroid** | INFRASTRUCTURE
Def: lunar azimuth as seen from the seismicity-weighted domain centroid, deg.
Src: needs the site-projection build (see F1-30).
Hunch: azimuth carries the horizontal-shear tidal component, which the zonal potential proxy discards;
strike-slip regions respond to shear, not to the zonal term.
Pit: a single centroid for a global domain is close to meaningless; this only makes sense inside the
2R-df regional instrument with a per-region centroid.
S15: as regional; per-region amplitude UNRESOLVED per §0.3.
Price: 31 per region form, but 1 test per (feature, lag) under the 2R-df sum.

**F1-18 moon_hour_angle** | INFRASTRUCTURE
Def: local hour angle of the Moon at a site, deg.
Src: site-projection build.
Hunch: the diurnal tidal clock, which is exactly the band the day-binning notches to zero.
Pit: **the sinc zeros are exact here**: K1 at 23.93 h suppresses by 366x, S1 by exactly 0. Any
day-binned test on this feature is a null by construction and must not be run as one.
S15: UNMEASURABLE under day binning. Requires the unbinned path (F9-15) before it is priced at all.
Price: 0 under current binning; 31 after the unbinned path exists.

**F1-19 lunar_semidiurnal_phase_M2** | INFRASTRUCTURE
Def: phase of the M2 constituent (12.4206 h) at a site.
Src: site projection.
Hunch: M2 is the largest tidal constituent and is the single most-cited tidal trigger candidate.
Pit: sinc factor 0.0348 under day binning (28.7x suppression) - this is the flagship example of the
notch. It is ALSO the band where the S2/S1 detection-cycle systematic lives, and the miner's current
immunity to that systematic is a genuine asset that the unbinned path gives up. Say so when building it.
S15: UNMEASURABLE day-binned. Post-unbinning, N=46,585 with sub-day timestamps: MEASURABLE.
Price: 0 now; 31 after F9-15.

**F1-20 solar_semidiurnal_phase_S2** | INFRASTRUCTURE
Def: S2 phase (exactly 12.000 h).
Src: site projection.
Hunch: entered only as the **artifact control** for F1-19 - S2 is where a detection-cycle artifact
would sit, so it is the calibration line, not a physics claim.
Pit: sinc exactly 0 day-binned. Post-unbinning it is confounded with the diurnal detection cycle and
must be run jointly with the G-M1 solar-cycle positive control (K-091).
S15: as F1-19.
Price: 0 now; 31 after F9-15, and it should be free-ridden on K-091's existing control.

**F1-21..F1-29 (constituent block)** | INFRASTRUCTURE
Def: the standard tidal constituents as site-local phases: O1 (25.819 h), K1 (23.934 h), P1, Q1, N2,
S2, K2, Mf (13.661 d), Mm (27.555 d), Ssa (182.6 d).
Src: site projection plus a constituent table (SPOTL or GOTIC2 conventions).
Hunch: the fortnightly and longer constituents are already effectively bounded by the miner's null; the
diurnal and semidiurnal block is the unexplored half of the spectrum and is unexplored **for an
instrumental reason, not a scientific one**.
Pit: constituents are mutually near-degenerate over a 21-year window at the sidelines (K1/P1 separate
only over 182.6 d; K2/S2 over 182.6 d); declare the resolution and do not claim to separate what the
window cannot.
S15: the sub-daily block is UNMEASURABLE until F9-15 exists; the long block is already bounded.
Price: 10 constituents x 31 lags = 310, contingent on F9-15, and I would run only the 4 diurnal/
semidiurnal ones first (124).

**F1-30 site_local_tidal_stress_tensor** | INFRASTRUCTURE
Def: the full solid-earth tidal stress tensor at a site, projected onto a receiver fault geometry, to
give Coulomb tidal stress in kPa as a daily (or hourly) series.
Src: **the missing build**. SPOTL (github.com/crustal/spotl) or a direct Love-number computation, plus
ocean loading from FES2014 (free, AVISO registration) and a receiver mechanism from `data/scedc_fm`.
Hunch: this is the single largest gap in the engine - §K87-0(d)(iv). Everything the program has
concluded about tides globally has been concluded from geocentric proxies, and the sign of the
Coulomb tidal stress is a function of fault geometry that a geocentric proxy cannot carry.
Pit: M-006 stands - site dilatation is unreliable to worse than a factor 2 including sign flips across
velocity solutions; the same warning applies to any modelled stress tensor. Ocean loading dominates
near coasts and is the largest error term. Do not quote signs you cannot defend.
S15: MEASURABLE once built, at the full global N; the whole point of the build is that it restores the
regional coherence the domain sum destroys.
Price: 1 feature x 31 lags per mechanism class; if entered with sign, 31 x (number of mechanism strata)
under S-13.
**This is the highest-leverage INFRASTRUCTURE item in the catalog.**

**F1-31 ocean_tidal_loading_amplitude** | INFRASTRUCTURE
Def: site-local ocean tidal loading stress amplitude (kPa) from FES2014 Green's function convolution.
Src: FES2014 via AVISO (free registration), or precomputed loading via the free onsite loading
calculators (holt.oso.chalmers.se/loading).
Hunch: near coasts and subduction zones the ocean load exceeds the solid tide, and most of the global
M>=4.5 catalog is in subduction zones - so the largest tidal stresses in the catalog are the ones the
engine does not model at all.
Pit: it is coastal, so it is spatially structured in a way that correlates with tectonic setting; any
apparent effect is confounded with subduction-versus-continental composition. S-13 stratification is
mandatory, not optional.
S15: MEASURABLE at full N once built.
Price: 31 per stratum.

**F1-32 tidal_stress_rate_site** | INFRASTRUCTURE
Def: d/dt of F1-30, kPa/day.
Src: as F1-30.
Hunch: rate-and-state predicts a rate-regime response below the nucleation time t_a; K-087 makes this
the estimable phase angle, and site-local stress is the only place where the phase angle has a physical
meaning rather than a proxy meaning.
Pit: §P5-8: `exp(dtau/Asigma)` is VOID as a power calculator at transient amplitudes; the link
function used to convert this to expected rate must be one demonstrated in range (S-17 candidate).
S15: MEASURABLE post-build.
Price: 31.

**F1-33 tidal_stress_second_moment_site** | INFRASTRUCTURE
Def: trailing 30 d variance of F1-30.
Src: as F1-30.
Hunch: K-064's superadditivity in its cheapest form - if the response is convex, variance converts to
rate, and a variance covariate detects that with no phase alignment required at all.
Pit: variance of a deterministic near-periodic series is itself near-periodic (the envelope), so this
is close to `perigee_syzygy` in disguise; report the collinearity.
S15: MEASURABLE post-build.
Price: 31.

**F1-34 earth_orientation_tidal_residual** | DOWNLOAD
Def: the tidal (short-period) component of LOD, isolated by removing the IERS-modelled atmospheric and
seasonal terms, in ms.
Src: already fetched: IERS finals.all.iau2000 (`download_features`), plus the IERS tidal model tables.
Hunch: LOD's short-period tidal variation is a direct measurement of the whole-Earth tidal deformation,
which makes it a **measured** rather than modelled tidal amplitude - and measured beats modelled.
Pit: LOD's tidal terms are themselves modelled and removed in the published series; naively using
finals.all gives you the residual after the model, which is the opposite of what you want. Read the
column definitions before implementing.
S15: MEASURABLE for a bound.
Price: 31.

**F1-35 planetary_tide_jupiter** | READY
Def: Jupiter's tidal potential at Earth, `(m_J/M_sun) * (a_sun/r_J)^3` scaled to the solar term.
Src: closed form (low-precision VSOP is plenty).
Hunch: **MECHANISM-FREE, and I am labelling it that explicitly.** Jupiter's tidal potential at Earth is
~1e-5 of the Moon's. There is no mechanism. It is included for exactly one reason: it is the perfect
**scale-calibration null** - a feature with a real astronomical period and a physically certain
zero effect, against which any claimed planetary or long-period detection can be sized.
Pit: proposing this as physics is numerology and is a corpse. Proposing it as a null-calibration ruler
is not. The distinction must survive into the report text or it will be misread.
S15: MEASURABLE for a bound, which is the whole point.
Price: 31.

**F1-36 planetary_synodic_jupiter_saturn** | READY
Def: phase of the Jupiter-Saturn synodic cycle (7253 d).
Src: closed form.
Hunch: **MECHANISM-FREE. No mechanism exists at any amplitude I can defend.** Entered only as a
long-period null-calibration line matched to the nodal cycle's period band, so that F1-05's result has
a same-band reference.
Pit: 1.06 cycles in the window; this is the numerology corpse's home address and must carry the
mechanism-free label in the output, not just here.
S15: UNMEASURABLE by window. Score neither way.
Price: 1.

**F1-37 planetary_synodic_venus** | READY
Def: Earth-Venus synodic phase (583.9 d).
Src: closed form.
Hunch: MECHANISM-FREE; a mid-period null-calibration line (13 cycles in window, so unlike F1-36 it is
actually estimable).
Pit: as F1-36; and 583.9 d beats against the 411.8 d perigean cycle in a way that could manufacture
apparent structure. Report the beat.
S15: MEASURABLE for a bound.
Price: 31.

**F1-38 barycentric_solar_wobble** | READY
Def: distance of the Sun from the solar-system barycentre (solar radii).
Src: closed form from planetary positions.
Hunch: MECHANISM-FREE at Earth (free-fall means the wobble is not felt as a tide); included because it
is the most-cited claim in the amateur literature and having a **priced, pre-registered null** on it is
worth more than not having one.
Pit: pure numerology if framed as physics. Frame it as adjudication of an existing external claim.
S15: MEASURABLE for a bound; window covers ~1 Jupiter orbit so long-period content is UNRESOLVED.
Price: 31.

**F1-39 solar_rotation_phase_27d** | DOWNLOAD
Def: Carrington rotation phase (27.2753 d).
Src: closed form / NOAA Carrington tables.
Hunch: the solar rotation modulates recurrent geomagnetic storms (F3 family), so this is the **clock**
for family 3, not a driver itself; if any geomagnetic effect exists, it should show at this period.
Pit: **27.2753 d is within 0.3% of the draconic month (27.212 d).** Over 7716 d these two lines are
separable but their windowed sidelobes overlap. This is the single most dangerous aliasing coincidence
in the whole catalog, and a detection on either one must report the other's power.
S15: MEASURABLE.
Price: 31.

**F1-40 solar_cycle_phase_11y** | DOWNLOAD
Def: phase of the sunspot cycle, from SILSO monthly numbers.
Src: sidc.be/SILSO (free).
Hunch: K-091 already uses the solar detection cycle as a free always-on positive control (G-M1); this
entry is that control **as a mining feature**, so its recovery calibrates every other family-3 result.
Pit: 1.9 cycles in the window; and the F107 flux feature already in the engine hit the best raw p in
the last session (p = 0.00114 at lag 3 d, BH q = 0.161) - which per §K87-0 is almost certainly the
detection cycle, i.e. an observer effect. Treat a hit here as an observer measurement until proven
otherwise; that is K-015 and K-031 territory, and it is a real result if it holds.
S15: MEASURABLE.
Price: 31.

---

# FAMILY 2 (extended) - BEATS, ENVELOPES, MOIRE, AND ENVELOPE-OF-ENVELOPE

The existing family 2 carries 8 features. The generative rule behind them (the "moire rule") has been
applied once, at depth 1. Depth 2 is unexplored and is cheap.

**F2-01 spring_neap_envelope** | READY
Def: |cos(synodic)| smoothed over 3 d: the amplitude envelope of the semidiurnal tide, not its phase.
Src: closed form.
Hunch: if the response is to tidal stress AMPLITUDE rather than to phase, the envelope is the correct
feature and the phase is the wrong one; the whole tidal corpse is a phase corpse, and an envelope is a
different object.
Pit: this is the closest live entry to the tidal corpse and must state the ~6.6% bound in its own text.
Envelope features are non-negative and skewed; standardisation is not enough, and the null must be the
block bootstrap at 2x the 14.765 d period.
S15: MEASURABLE, A_min 5.5% deflated. The bound already nearly excludes the interesting range.
Price: 31.

**F2-02 declination_envelope** | READY
Def: 13.66 d envelope of the declination tide; |sin(2*draconic)| smoothed.
Src: closed form.
Hunch: same logic as F2-01 for the declination band, which is the band where the day-binning sinc is
~0.99 and the instrument is at full sensitivity.
Pit: near-degenerate with existing `moon_abs_declination`; report min R^2 against it.
S15: MEASURABLE.
Price: 31.

**F2-03 envelope_of_envelope_perigean** | READY
Def: the 411.8 d modulation depth of the 14.765 d envelope, i.e. amplitude of the amplitude.
Src: closed form.
Hunch: the depth-2 moire term; if a response is convex in stress, the second-order envelope carries
signal that the first-order one averages away, and nobody has computed it.
Pit: second-order on an already-bounded first-order effect; the honest expectation is a bound. Report
the implied ceiling BEFORE running.
S15: MARGINAL. Bound-producing at best.
Price: 31.

**F2-04 beat_draconic_anomalistic** | READY
Def: `mod(draconic - anomalistic, 2pi)`, period ~205.9 d.
Src: closed form.
Hunch: the one two-cycle beat missing from the existing family-2 set (which has synodic-anomalistic,
synodic-draconic, annual-synodic); its absence is an asymmetry rather than a hypothesis.
Pit: 205.9 d is close to the 182.6 d semiannual line and to the 173.3 d eclipse half-year already
present; sidelobe overlap must be reported.
S15: MEASURABLE.
Price: 31.

**F2-05 beat_annual_draconic** | READY
Def: `mod(annual - draconic, 2pi)`, ~27.9 d.
Src: closed form.
Hunch: completes the beat lattice; mechanism-free individually, but a lattice with a hole in it cannot
support the claim "we scanned the beat space".
Pit: near the 27.2 d and 27.55 d lines; aliasing declaration mandatory.
S15: MEASURABLE.
Price: 31.

**F2-06 beat_annual_anomalistic** | READY
Def: `mod(annual - anomalistic, 2pi)`, ~27.7 d.
Src: closed form.
Hunch: as F2-05, completing the lattice.
Pit: the three ~27 d lines together are a resolution trap: over 7716 d they separate, but their
sidelobes do not, so a max-statistic over the three is the only honest summary.
S15: MEASURABLE.
Price: 31.

**F2-07 triple_beat_syn_ano_dra** | READY
Def: `mod(synodic - anomalistic - draconic, 2pi)`.
Src: closed form.
Hunch: SPECULATIVE. A three-cycle coincidence marker; the only defence is that the largest tidal
excursions are three-condition coincidences (perigee plus syzygy plus node), and no existing feature
expresses a triple condition.
Pit: full numerology risk. Must carry the mechanism-free label unless framed as the smooth proxy for
the triple-coincidence indicator, in which case say that.
S15: MEASURABLE.
Price: 31.

**F2-08 saros_phase** | READY
Def: phase of the 6585.32 d saros cycle.
Src: closed form.
Hunch: MECHANISM-FREE beyond being the eclipse-repeat period; entered as a long-period null line
matched in period to F1-05 (nodal), giving that entry a same-band control.
Pit: 1.17 cycles in window. Numerology-adjacent; label accordingly.
S15: UNMEASURABLE by window.
Price: 1.

**F2-09 metonic_phase** | READY
Def: phase of the 19-year Metonic cycle (6939.6 d).
Src: closed form.
Hunch: MECHANISM-FREE, calendrical rather than physical. Included so the catalog distinguishes physical
from calendrical cycles, and to give the nodal claim a calendrical control.
Pit: pure numerology as physics; its value is entirely as a control line.
S15: UNMEASURABLE by window.
Price: 1.

**F2-10 potential_times_declination_envelope** | READY
Def: product of `tidal_potential_proxy` with `moon_abs_declination`.
Src: closed form.
Hunch: interaction terms are how a convex response shows up in a linear model; the existing set has
three products and this is the missing one between the two strongest tidal proxies.
Pit: products of correlated features are dominated by the squared term; centre both first or you are
testing the square of one feature under a misleading name.
S15: MEASURABLE.
Price: 31.

**F2-11 potential_squared** | READY
Def: `tidal_potential_proxy^2`, centred.
Src: closed form.
Hunch: K-064's convexity claim in its cheapest form: a convex response to a zero-mean forcing produces
a mean-rate offset proportional to the variance, and the square is the feature that sees it.
Pit: the square of a near-sinusoid is a sinusoid at twice the frequency plus a constant, so this is
largely the ladder rung 2P in disguise. Report the identity or pay twice for one test.
S15: MEASURABLE.
Price: 31.

**F2-12 potential_cubed** | READY
Def: cubed, centred proxy.
Src: closed form.
Hunch: distinguishes a threshold response (odd, sharp) from a convex response (even); the odd/even
distinction is exactly what the 2-df omnibus refuses to report.
Pit: as F2-11 with the 3P rung.
S15: MEASURABLE for a bound.
Price: 31.

**F2-13 exceedance_fraction_above_quantile** | READY
Def: fraction of the day for which the tidal potential proxy exceeds its own 90th percentile.
Src: closed form (needs sub-daily evaluation of the geocentric proxy, but not site projection).
Hunch: a threshold-dwell feature: if triggering is threshold-crossing, the relevant quantity is time
above threshold, which no sinusoidal feature expresses.
Pit: the quantile is a free parameter and therefore a multiplicity axis; declare ONE value per S-9 with
no alternatives run, or price the scan.
S15: MEASURABLE.
Price: 31 per declared quantile; I would declare 0.90 and stop.

**F2-14 time_since_last_extremum** | READY
Def: days since the last local maximum of the combined zonal potential.
Src: closed form.
Hunch: if the response has a fixed delay after peak stress rather than a fixed phase, this feature sees
it and phase features do not.
Pit: strongly sawtooth, so its own autocorrelation sets the block length; do not use a flat 90 d block
(the F107 incident in the `mine.py` docstring is exactly this failure mode).
S15: MEASURABLE.
Price: 31.

**F2-15 stress_cycle_count_natural_time** | READY
Def: cumulative count of tidal-potential zero-upcrossings since window start, entered as an axis rather
than a covariate (see family F8).
Src: closed form.
Hunch: K-017's natural-time move applied to the forcing rather than to the catalog: if the system
integrates cycles rather than days, this is the correct clock.
Pit: monotone increasing, so collinear with time itself and with every secular catalog trend (network
growth, magnitude re-homogenisation). Detrending is mandatory and changes what is tested.
S15: MEASURABLE but interpretation UNRESOLVED without a completeness control (family F7).
Price: 31.

**F2-16 amplitude_modulation_index** | READY
Def: (envelope max minus min) / (max plus min) over a trailing 30 d window of the combined proxy.
Src: closed form.
Hunch: a dimensionless modulation depth; if the system responds to relative rather than absolute stress
variation, this is the invariant quantity and it is not in the set.
Pit: ratio features are unstable when the denominator is small; clip and declare the clip. Also: the
TSI ratio corpse was a ratio feature, and the failure mode there was exactly denominator instability.
S15: MEASURABLE.
Price: 31.

**F2-17 near_degenerate_control** | READY
Def: d/dt of the perigean-spring beat phase, a near-constant series.
Src: closed form.
Hunch: entered as a degenerate-feature control to verify that the engine's zero-variance guard behaves
and that near-degenerate features do not manufacture significance.
Pit: it will nearly trip `design()`'s sd <= 0 guard; that is intentional and must be documented, not
fixed away.
S15: not applicable, it is a pipeline test.
Price: 1.

**F2-18..F2-25 harmonic ladder block (8 entries)** | READY
Def: the {P/3, P/2, P, 2P, 3P} ladder applied to the 8 kind='linear' family-1/2 features, which
currently receive no ladder at all.
Src: closed form.
Hunch: a harmonic is a new feature, never a free axis (§P5-5); a response that is sharp in phase has
most of its power in harmonics, not in the fundamental, and the linear features have never been given
the ladder.
Pit: rungs are never free and the price is real. Rung 2P of a squared feature is the same test as
F2-11; deduplicate before pricing.
S15: MEASURABLE.
Price: 8 features x 4 new rungs = 32 (I would not lag-scan the ladder in tranche 1).

**F2-26 cross_family_moire_with_solar** | READY
Def: `mod(synodic - carrington, 2pi)` where carrington is the 27.2753 d solar rotation.
Src: closed form.
Hunch: SPECULATIVE, and its real value is diagnostic: a lunar-solar-rotation beat has no mechanism, so
a detection here is evidence of a **detection-cycle artifact** interacting with lunar sampling, which is
an observer finding (K-015/K-031) and a real one.
Pit: the 0.3% draconic/Carrington coincidence (F1-39) makes this beat extremely long-period and poorly
resolved; declare the resolution.
S15: MARGINAL.
Price: 31.

---

# FAMILY 3 (extended) - DOWNLOADED GEOPHYSICAL AND SPACE-ENVIRONMENT SERIES

The engine currently downloads two sources giving three features (Ap, F10.7, LOD). This family is where
the cheapest breadth lives: every entry is a free daily series, a parser, and a cache line.

**A standing warning that applies to this whole family.** Solar and geomagnetic indices are correlated
with the **solar detection cycle**, which is an observer effect on catalog completeness (K-015, K-031,
K-091). The last session's best raw p was `F107_solar_flux` at lag 3 d. Under this catalog's rules, a
family-3 hit is an **observer measurement until proven otherwise**, and proving otherwise requires the
completeness controls in family F7. That is not a reason to skip the family; it is the reason the family
is interesting.

## 3A. Geomagnetic and solar-wind

**F3-01 Kp_daily_max** | DOWNLOAD
Def: daily maximum Kp (quasi-log 0-9).
Src: already fetched by `download_features` (GFZ Kp_ap_Ap_SN_F107_since_1932.txt); only Ap is currently
extracted, Kp is in the same file.
Hunch: Kp is the standard storm index; the max is a different statistic from the Ap daily mean, and
extreme-value features often beat mean features for threshold-like physics.
Pit: Kp is quasi-logarithmic and quantised to thirds, so its distribution is lumpy; standardisation is
misleading and a rank transform is safer.
S15: MEASURABLE, A_min 5.5% deflated.
Price: 31.

**F3-02 Ap_daily_max_3h** | DOWNLOAD
Def: maximum of the eight 3-hour ap values in a day.
Src: same file.
Hunch: sub-daily peak amplitude, which the daily Ap average washes out; the physical hypothesis (if
any) is a transient-induced-current effect, and transients are peaks, not means.
Pit: correlated with F3-01 at r > 0.9; the pair must be entered as one plus a residual, not as two.
S15: MEASURABLE.
Price: 31.

**F3-03 Ap_storm_onset_indicator** | DOWNLOAD
Def: 1 on the first day where Ap crosses 40 after >= 5 quiet days, decaying over 3 d.
Src: same file.
Hunch: **storm SEQUENCE rather than storm LEVEL** - the Jim directive to model the driver's structure,
not its amplitude. Onsets are sharp, sudden-commencement events; a level covariate cannot represent one.
Pit: threshold (40) and quiet-run length (5) are free parameters; declare one value each per S-9 or pay
for the grid. Onsets are rare (~50-80 in the window): effective N is small.
S15: effective N ~2,000-4,000 event-days. A_min ~19-27% deflated. MARGINAL to UNMEASURABLE.
Price: 31.

**F3-04 Ap_recovery_phase** | DOWNLOAD
Def: days since last storm onset, capped at 10, 0 when no recent storm.
Src: same file.
Hunch: the negative-space partner of F3-03; if there is a delayed response it lives in the recovery
phase, and the existing lag scan can find it only if the response is a fixed lag rather than a
storm-relative one.
Pit: fully redundant with the lag scan on F3-01 IF the response is linear; not redundant if the response
is storm-conditional. State that distinction or the entry is unpriceable.
S15: MEASURABLE for a bound.
Price: 31.

**F3-05 Ap_storm_count_27d** | DOWNLOAD
Def: trailing 27 d count of storm onsets.
Src: same file.
Hunch: cumulative-dose rather than instantaneous forcing; if any mechanism exists it is more likely to
be integrative (thermal, electrochemical) than instantaneous.
Pit: trailing sums of rare events are lumpy and heavily autocorrelated; block length must be set from
the measured ACF, not defaulted.
S15: MEASURABLE.
Price: 31.

**F3-06 Dst_index** | DOWNLOAD
Def: hourly Dst, taken as daily minimum (nT).
Src: wdc.kugi.kyoto-u.ac.jp/dstdir (free).
Hunch: Dst measures ring-current strength and is the standard magnetospheric storm scalar; it is a
cleaner physical quantity than Kp and it is the index most of the claimed EQ-geomagnetic literature uses.
Pit: provisional versus final Dst differ; the series has a documented baseline drift. Use the final
series only for the exploration window and declare the version.
S15: MEASURABLE.
Price: 31.

**F3-07 SYM_H** | DOWNLOAD
Def: 1-min SYM-H, daily minimum.
Src: OMNIWeb (omniweb.gsfc.nasa.gov), free.
Hunch: the high-resolution Dst analogue; its daily minimum captures sharp storm mains that hourly Dst
smooths.
Pit: near-duplicate of F3-06; enter one and the residual of the other.
S15: MEASURABLE.
Price: 31.

**F3-08 AE_index** | DOWNLOAD
Def: daily mean auroral electrojet index.
Src: WDC Kyoto, free.
Hunch: AE is a high-latitude index; if a geomagnetic effect exists it should be latitude-structured, so
AE versus Dst is a **discriminator between high-latitude and equatorial mechanisms**, not two copies of
one hypothesis.
Pit: AE has known coverage gaps and station-count changes over the record: this is catalog-composition
drift in the covariate, which is the same disease as in the seismic catalog. Declare it.
S15: MEASURABLE.
Price: 31.

**F3-09 solar_wind_speed** | DOWNLOAD
Def: daily mean solar wind bulk speed, km/s.
Src: OMNIWeb.
Hunch: the driver behind recurrent storms; if the 27 d Carrington line matters, speed is the physical
carrier and Ap is the response.
Pit: recurrent high-speed streams put a strong 27 d line into this feature, which collides with the
draconic month at 0.3% (F1-39). This is the aliasing trap of the catalog.
S15: MEASURABLE.
Price: 31.

**F3-10 IMF_Bz** | DOWNLOAD
Def: daily mean southward IMF component, nT.
Src: OMNIWeb.
Hunch: Bz south is the coupling switch for storms; entering it separates coupling from wind energy.
Pit: near-zero mean with heavy tails; the daily mean destroys most of its information. Consider the
daily minimum instead and declare which.
S15: MEASURABLE.
Price: 31.

**F3-11 solar_wind_dynamic_pressure** | DOWNLOAD
Def: daily mean `n * v^2`, nPa.
Src: OMNIWeb.
Hunch: the only solar-wind quantity that exerts an actual mechanical stress on the magnetosphere; it is
the one with an arguable, if tiny, mechanical channel to the solid Earth.
Pit: the implied crustal stress is many orders below the tidal stress that is already bounded at ~6%.
The AMPLITUDE_HONESTY banner in `mine.py` applies verbatim: read any result as an upper bound.
S15: MEASURABLE for a bound only.
Price: 31.

**F3-12 F107_flux** (existing) plus **F3-13 F107_81d_centred** | DOWNLOAD
Def: existing daily flux, plus the 81-day running mean (the standard solar-activity smoother).
Src: same GFZ file.
Hunch: separates the rotation-modulated flux from the slow cycle; the last session's hit was at lag 3 d
on the raw flux, and whether that lives at the 27 d or 11 y timescale is decidable with this pair.
Pit: **the 81 d smoother is a low-pass filter, so its block-bootstrap length must be set from the
smoothed ACF, not the raw.** The `mine.py` docstring records exactly this failure (90 d blocks gave
F10.7 a z of 32).
S15: MEASURABLE.
Price: 31 for the new smoothed feature.

**F3-14 sunspot_number** | DOWNLOAD
Def: daily international sunspot number.
Src: SILSO, sidc.be/SILSO (free).
Hunch: the canonical solar-activity index and the one most likely to correlate with the network's own
duty cycle; it is a **completeness proxy in disguise**, and I want it entered as such.
Pit: version 2 recalibration (2015) introduces a step change inside the exploration window. That step is
a composition artifact and will manufacture apparent structure at the 11 y band.
S15: MEASURABLE.
Price: 31.

**F3-15 cosmic_ray_neutron_count** | DOWNLOAD
Def: daily neutron monitor count rate (Oulu or Moscow station).
Src: cosmicrays.oulu.fi (free), NMDB (nmdb.eu).
Hunch: anticorrelated with solar activity, so it is a **sign-flipped duplicate** of the solar features
and therefore an excellent internal consistency check: a real solar effect must show opposite signs
here and in F10.7, and a detection-cycle artifact need not.
Pit: station-specific cutoff rigidity and pressure corrections; use the corrected series and declare the
station.
S15: MEASURABLE.
Price: 31.

## 3B. Earth rotation and orientation

**F3-16 LOD** (existing) plus **F3-17 dLOD/dt** | DOWNLOAD
Def: existing length-of-day, plus its first difference (ms/d).
Src: IERS finals.all, already fetched.
Hunch: the decadal LOD-seismicity claim in the literature is about **rate of change**, not level; the
engine currently carries only the level, so the published version of the hypothesis is not in the engine.
Pit: differencing a series with a strong annual and semiannual term amplifies those terms; the result
will be dominated by the seasonal cycle unless it is removed first, and removing it makes the entry a
different hypothesis. State which one you are running.
S15: MEASURABLE.
Price: 31.

**F3-18 LOD_decadal_residual** | DOWNLOAD
Def: LOD with tidal, annual, semiannual and 90-day terms removed (the "core" residual).
Src: IERS plus the IERS tidal model.
Hunch: the decadal LOD signal is core-mantle coupling, which is the only geophysically plausible channel
in this sub-family; the literature claim (Bendick and Bilham) is specifically about this component.
Pit: 21 y of data holds ~2 decadal cycles; and the removal model is itself a choice with free
parameters. **The known claim is about M>=7 global rates, which is a much smaller N than 46,585.**
S15: at M>=7 the global count over the window is roughly 300-350 events, giving A_min ~65% deflated:
**UNMEASURABLE at M>=7**, MEASURABLE at M>=4.5 where the claim was never made. Say both.
Price: 31.

**F3-19 LOD_annual_amplitude** | DOWNLOAD
Def: amplitude of the annual LOD term over a trailing 400 d window.
Src: IERS.
Hunch: LOD's annual term is atmospheric angular momentum; this is an **atmospheric loading index in
disguise** and is the cheapest global proxy for family F6.
Pit: interpretation is entirely atmospheric, so it must be conditioned jointly with F6 or it will be
misattributed to rotation.
S15: MEASURABLE.
Price: 31.

**F3-20 polar_motion_x**, **F3-21 polar_motion_y** | DOWNLOAD
Def: IERS pole coordinates x, y (arcsec), daily.
Src: IERS finals.all, already downloaded; columns present, currently unparsed.
Hunch: polar motion changes the centrifugal potential, which produces a genuine (if small) global stress
field with a computable pattern; unlike most of family 3 this one has an actual stress channel.
Pit: the pole path is a spiral, so x and y are strongly phase-shifted copies; enter the pair as a phase
feature (F3-22), not as two linears, or pay twice for one hypothesis.
S15: MEASURABLE.
Price: 62 as two linears; 31 as one phase feature, which is what I recommend.

**F3-22 chandler_wobble_phase** | DOWNLOAD
Def: phase of the ~433 d Chandler component, isolated by band-pass on the pole path.
Src: IERS.
Hunch: the Chandler wobble is a free oscillation of the Earth with a real, computable centrifugal stress
field of order 1-2 kPa - **comparable in size to the tidal stresses the program has spent years on, and
never once tested.** This is my favourite entry in family 3.
Pit: 433 d beats against the annual to give the well-known ~6.4 y envelope, and the annual is entangled
with every seasonal artifact in the catalog. The band-pass has free parameters; declare one.
S15: MEASURABLE, A_min 5.5% deflated. ~17.8 cycles in the window, which is a properly estimable period.
Price: 31.

**F3-23 chandler_annual_beat_6_4y** | DOWNLOAD
Def: phase of the ~2350 d beat between Chandler and annual polar motion.
Src: IERS.
Hunch: the polar-motion amplitude envelope; if the centrifugal stress hypothesis has an amplitude
dependence, this is where it shows.
Pit: 3.3 cycles in window. Marginal resolution.
S15: MARGINAL.
Price: 31.

**F3-24 centrifugal_stress_proxy** | DOWNLOAD
Def: the degree-2 order-1 centrifugal potential change from polar motion, converted to a surface stress
proxy (kPa) with a fixed Love-number set.
Src: IERS plus closed-form Love numbers.
Hunch: the physically correct scalar for F3-20..F3-23; it turns three correlated coordinate features
into one stress feature with units, which is what S-17 wants.
Pit: it is latitude and longitude dependent, so the global domain sum cancels it exactly the way it
cancels the diurnal tide (§K87-0(d)(i)). **This feature is only meaningful through the 2R-df regional
instrument.** Do not run it globally and call it a null.
S15: global sum UNMEASURABLE by cancellation; regional 2R-df MEASURABLE.
Price: 31 under the 2R-df sum (1 test per feature-lag).

**F3-25 UT1_minus_UTC_rate** | DOWNLOAD
Def: rate of change of UT1-UTC.
Src: IERS.
Hunch: equivalent to LOD by construction; entered as a **redundancy control** to check that the pipeline
gives the same answer to the same physical question posed two ways.
Pit: exact redundancy; if it produces a different answer from F3-17 there is a bug, and that is the
value of running it.
S15: MEASURABLE.
Price: 31, and it is a pipeline test, not a hypothesis.

## 3C. Gravity, geodesy, and global loading fields

**F3-26 GRACE_equivalent_water_thickness** | DOWNLOAD
Def: monthly GRACE/GRACE-FO mass anomaly (cm equivalent water) interpolated to daily, per region.
Src: podaac.jpl.nasa.gov GRACE mascon solutions (free), or GFZ.
Hunch: continental water storage produces surface loading stresses of 1-10 kPa, which is larger than the
solid-earth tide - the biggest routinely-varying stress on the crust that is not a tide.
Pit: monthly sampling means everything below 60 d is aliased; the 161 d S2 tidal alias in GRACE is a
notorious artifact. **Interpolating monthly to daily manufactures smoothness that is not in the data**,
and the block bootstrap will be fooled by it. K-037 already flagged this family as NEEDS-DATA.
S15: the effective independent sample is ~250 months, not 7716 days: A_min ~75% deflated on the monthly
axis. **UNMEASURABLE as a daily feature; MEASURABLE only as a slow regional covariate inside a
conditional (family F10).**
Price: 31 per region form, 1 under the 2R-df sum.

**F3-27 GRACE_loading_rate** | DOWNLOAD
Def: month-to-month difference of F3-26.
Src: as above.
Hunch: rate of loading, not level - the rate-and-state relevant quantity, and the form in which the
seasonal-hydrology literature states its result.
Pit: differencing monthly data doubles the noise; and **Sirorattanakul and Avouac own the California
seasonal hydrology result** (M-010). Any use here must be global-extension or conditioning, and must
cite them. Not a novelty claim.
S15: as F3-26.
Price: 31 per region form.

**F3-28 GLDAS_soil_moisture** | DOWNLOAD
Def: daily 0.25 deg soil-moisture and snow-water-equivalent from GLDAS/NOAH.
Src: disc.gsfc.nasa.gov (free, Earthdata login).
Hunch: the daily-resolution version of F3-26; it fixes the aliasing problem GRACE has, at the cost of
being a model output rather than a measurement.
Pit: it is a model, so its variability is the model's, and model reanalysis versions change inside the
window. Composition drift in the covariate.
S15: MEASURABLE at daily resolution regionally; global sum cancels the seasonal phase between
hemispheres exactly (northern and southern seasons are antiphase), so **the global test is
self-cancelling by construction** and only the 2R-df regional instrument can see it. This is the cleanest
example in the catalog of why Rule 4.2 exists.
Price: 1 per (feature, lag) under the 2R-df sum; 31 for the lag scan.

**F3-29 ERA5_surface_pressure** | DOWNLOAD
Def: daily mean surface pressure anomaly, hPa, regionally averaged.
Src: ERA5 via Copernicus CDS (free registration).
Hunch: 10 hPa is 1 kPa of vertical load, applied over synoptic scales; the pressure field is the fastest
large-amplitude surface loading there is, and it has a sign that is computable per fault.
Pit: K-037's killer-falsifier framing requires the predicted SIGN per fault, which requires mechanisms
(`data/scedc_fm` regionally, or the global CMT catalog). Without sign it is a weak two-sided test.
S15: regional MEASURABLE; global sum near-cancelling.
Price: 31 under 2R-df.

**F3-30 ERA5_pressure_rate** | DOWNLOAD
Def: 24 h change in surface pressure.
Src: ERA5.
Hunch: storm passage is a rate phenomenon; the rate-versus-level axis (K-087) applied to the atmosphere.
Pit: as F3-29.
Price: 31 under 2R-df.

**F3-31 typhoon_hurricane_proximity** | DOWNLOAD
Def: minimum distance from a region centroid to any active tropical cyclone centre, and the associated
central pressure deficit.
Src: IBTrACS (ncei.noaa.gov/products/international-best-track-archive), free.
Hunch: the largest, fastest surface unloading events on Earth (50+ hPa in hours), and there is a
published slow-slip-triggering literature; this is a **transient**, which is exactly the class the tidal
bound does not cover (S-14(c)).
Pit: cyclones are seasonal and tropical, so they correlate with subduction-zone seismicity through pure
geography. Stratify by region per S-13 or the result is a map of the tropics.
S15: rare-event feature; effective N per region is small. MARGINAL, and honest only in the tropical
strata.
Price: 31 under 2R-df.

**F3-32 sea_level_anomaly** | DOWNLOAD
Def: regional daily altimetric sea-level anomaly, cm.
Src: AVISO/Copernicus (free).
Hunch: ocean mass loading independent of the tide; a slow, real, measurable load on subduction zones,
where most of the M>=4.5 catalog lives.
Pit: strongly seasonal and ENSO-driven, so it is collinear with F6 entries; and altimeter missions change
inside the window (composition drift).
S15: regional MEASURABLE.
Price: 31 under 2R-df.

## 3D. Space-weather sequence features (the driver's structure, not its level)

**F3-33 storm_interarrival_time** | DOWNLOAD
Def: days since previous storm onset (Ap > 40).
Src: GFZ.
Hunch: the driver's own clock; if the system responds to storm sequences rather than storms, the
interarrival is the state variable, and no level feature carries it.
Pit: as F3-03, small effective N.
S15: MARGINAL.
Price: 31.

**F3-34 storm_sequence_order** | DOWNLOAD
Def: index of the storm within a recurrent-stream sequence (1st, 2nd, 3rd rotation of the same stream).
Src: GFZ plus OMNIWeb speed.
Hunch: **habituation versus sensitisation** - if a system responds to the first shock of a sequence and
not the third, that is a signature no level covariate can produce, and it is the same shape as the
"wetness meter" idea in K-073.
Pit: the sequence assignment is a judgement call and therefore a hidden parameter; write the rule down
per S-9 before running.
S15: MARGINAL.
Price: 31.

**F3-35 quiet_run_length** | DOWNLOAD
Def: consecutive days with Ap < 10.
Src: GFZ.
Hunch: the negative space of family 3 - model the silence, not the activity. If storms matter, the
absence of storms is the control condition, and it is a longer-lived state than the storms are.
Pit: heavily autocorrelated by construction; the block bootstrap must use the measured ACF.
S15: MEASURABLE.
Price: 31.

**F3-36 geomagnetic_jerk_indicator** | DOWNLOAD
Def: second time derivative of the core field's secular variation at a global observatory average.
Src: INTERMAGNET definitive data (free), or the CHAOS field model.
Hunch: geomagnetic jerks are core-driven and are the same physical family as the decadal LOD claim; if
core-mantle coupling has a surface expression, jerks and LOD should agree, and their agreement or
disagreement is more informative than either alone.
Pit: jerks are ~5-10 events per decade: effective N is tiny. This is a bound-producer.
S15: UNMEASURABLE except as a bound.
Price: 31.

**F3-37 schumann_resonance_amplitude** | INFRASTRUCTURE
Def: daily amplitude of the 7.83 Hz Schumann fundamental.
Src: no free continuous global archive covering 2005-2026 that I can name; would need a station
partnership. **Labelled INFRASTRUCTURE honestly rather than DOWNLOAD.**
Hunch: MECHANISM-FREE for triggering; it is a heavily-claimed precursor in the fringe literature and a
priced, pre-registered null on it would be a public good.
Pit: no defensible mechanism, no free data, high prior of instrument artifacts. Low priority and I say so.
S15: unknown until data exists.
Price: 31 if it ever exists.

**F3-38 ionospheric_TEC** | DOWNLOAD
Def: daily global-mean and regional total electron content, TECU.
Src: IGS/CDDIS global ionosphere maps (free).
Hunch: the pre-earthquake TEC anomaly literature is large, contested, and almost entirely
retrospective; **entering TEC as a prospective covariate in a pre-registered miner is the cheapest way
this program could contribute to a real controversy.**
Pit: TEC is dominated by solar activity and the diurnal cycle, so it is a solar feature wearing a
different hat; the residual after solar regression is the only honest form. Also, the published claims
are about the days BEFORE an event, which makes this a **precursor** test, and precursor tests need
the anti-leakage discipline of family F7 plus a strictly causal construction.
S15: MEASURABLE.
Price: 31 for the solar-residualised form.

**F3-39 groundwater_well_levels** | DOWNLOAD
Def: daily water level at long-record USGS wells in seismically active western US basins.
Src: waterdata.usgs.gov (free API).
Hunch: a direct in-situ measurement of pore pressure, which is the quantity rate-and-state actually
cares about; every other loading feature in this catalog is a proxy for it.
Pit: wells are pumped; anthropogenic drawdown dominates natural signal in most basins. Site selection is
a hidden multiplicity axis and the selection rule must be declared and exploration-window-only.
S15: regional (western US) only; N in that region over the window is a small fraction of global.
MARGINAL.
Price: 31 per declared basin, and I would declare at most 3.

**F3-40 reservoir_level_time_series** | DOWNLOAD
Def: daily impoundment level at large reservoirs near active faults.
Src: national agency portals (free, varies by country); US: USBR/USACE.
Hunch: reservoir-induced seismicity is a **known-true positive control** for pore-pressure triggering
(the K-034 logic applied to a different forcing), which makes it a calibration asset rather than a
discovery target.
Pit: known-true means owned; this is not a novelty claim. Its role is to calibrate the response curve
that the tidal entries need and cannot get (K-045's logic, K-036 Tier 2).
S15: local, high-SNR, small N. MEASURABLE within the induced regime only.
Price: 31 per declared site.

---

# FAMILY 4 (extended) - CATALOG-ENDOGENOUS FEATURES

The engine currently carries three: trailing 90 d b-value, 30 d deep fraction, 30 d mean depth. This is
the family with the most headroom, the least download cost, and the most severe leakage risk.

**Three rules that bind this whole family.**
1. **Strict causality.** `z[:, t]` uses only events strictly before day t. The engine's causality test
   (shuffle future events, z must not change) must pass for every entry here.
2. **The baseline absorbs clustering or the feature is measuring clustering.** Under
   `--baseline climatology` every trailing-count feature scores fake skill; that is SPEC decision 2's
   declared gap. Every entry here is priced against `--baseline etas` and says so.
3. **Composition drift is the enemy.** Network growth, magnitude-scale re-homogenisation, and reporting
   changes all produce secular trends in every trailing statistic. Any entry whose signal is secular is
   an observer measurement (family F7) until proven otherwise.

## 4A. b-value and magnitude-distribution dynamics

**F4-01 b_value_30d** | READY
Def: trailing 30 d global Aki b, same estimator as the existing 90 d feature.
Src: on disk.
Hunch: the existing 90 d window is a choice, not a measurement; the timescale on which b carries
information is unknown and is itself the interesting quantity.
Pit: **window length is a multiplicity axis**, and adding 30/180/365 versions is a scan that must be
priced. Short windows at M>=4.5 give ~380 events per 30 d globally, so b is noisy at sigma_b ~ b/sqrt(N).
S15: b's own error at N=380 is ~5%; a rate response to a 5%-noise covariate is attenuated accordingly.
MEASURABLE but attenuated; state the attenuation.
Price: 31 per window length.

**F4-02 b_value_365d** | READY
Def: trailing 365 d global Aki b.
Src: on disk.
Hunch: the slow-drift version; separates a real stress signal (fast) from composition drift (slow), and
the difference of the two is more informative than either.
Pit: at 365 d the feature is nearly a secular trend and will correlate with network growth.
S15: MEASURABLE.
Price: 31.

**F4-03 b_value_ratio_30_365** | READY
Def: b_30d / b_365d.
Src: on disk.
Hunch: **the composition-drift-cancelling form.** A ratio of a fast to a slow estimate of the same
quantity removes the slow secular artifact by construction, which is the honest way to ask "is b
anomalous now" without asking "has the network changed".
Pit: ratios are unstable; and the TSI-ratio corpse died partly on ratio instability. Clip and declare.
S15: MEASURABLE.
Price: 31.

**F4-04 b_value_uncertainty** | READY
Def: Aki-Utsu standard error of the trailing 90 d b, `b/sqrt(N)`.
Src: on disk.
Hunch: it is a **count feature in disguise** and therefore a completeness proxy; entered explicitly so
that any b-value result can be checked against the possibility that it is a count result.
Pit: nearly a deterministic function of the trailing count; report the collinearity with F4-20.
S15: MEASURABLE.
Price: 31.

**F4-05 b_value_spatial_dispersion** | READY
Def: standard deviation of per-region b across the R pre-declared regions, trailing 365 d.
Src: on disk.
Hunch: **K-021's claim that heterogeneity of b, not mean b, is informative** - the heterogeneity is a
system-level order parameter and the mean is a spatial average that destroys it.
Pit: per-region b at M>=4.5 over 365 d has ~120 events per region at R=20, so sigma_b ~9% per region;
the dispersion of a noisy estimator is dominated by estimation noise unless it is bias-corrected.
**Bias correction is mandatory here, not optional.**
S15: MARGINAL. The bias-corrected version may be MEASURABLE; declare the correction before running.
Price: 31.

**F4-06 b_value_depth_contrast** | READY
Def: b(shallow, z<70 km) minus b(deep, z>70 km), trailing 365 d.
Src: on disk.
Hunch: b is stress-dependent and depth is a proxy for the stress and thermal regime; the contrast is a
differential measurement that cancels global reporting artifacts affecting both populations.
Pit: the depth split at 70 km is the existing `deep_fraction` threshold; keep it identical so the two
features are interpretable together. Depth errors in global catalogs are large and depth is often fixed
at 10 or 33 km by convention, which puts a spurious spike in the shallow population.
S15: MEASURABLE.
Price: 31.

**F4-07 magnitude_gap_max_minus_second** | READY
Def: (largest magnitude in trailing 30 d) minus (second largest).
Src: on disk.
Hunch: Bath's-law-like gap statistic as a **state variable of the whole catalog**; a large gap means one
event dominates, which is a different dynamical state from many comparable events.
Pit: extreme-order statistics have heavy-tailed sampling distributions; the block bootstrap handles this
but the standardisation does not. Rank transform.
S15: MEASURABLE.
Price: 31.

**F4-08 max_magnitude_30d** | READY
Def: largest magnitude in trailing 30 d.
Src: on disk.
Hunch: the crudest possible state variable and therefore the necessary control for everything else in
this sub-family: if F4-07 works and this does not, the gap is real; if both work, it is just size.
Pit: trivially correlated with subsequent aftershock rate, so under a climatology baseline it scores fake
skill by construction. **ETAS baseline mandatory.**
S15: MEASURABLE.
Price: 31.

**F4-09 magnitude_distribution_KS_vs_GR** | READY
Def: KS distance between the trailing 365 d magnitude distribution and a fitted Gutenberg-Richter.
Src: on disk.
Hunch: **K-007's loss-of-GR claim as a continuous covariate** rather than a precursor test; a bend or
taper is a departure from GR and the KS distance is agnostic about which.
Pit: KS is insensitive in the tail, which is exactly where the taper lives; an Anderson-Darling or an
upper-tail-weighted statistic is better and should be declared instead. Also, incompleteness produces a
GR departure at the low end that has nothing to do with physics.
S15: MEASURABLE.
Price: 31.

**F4-10 GR_taper_corner_magnitude** | READY
Def: fitted corner magnitude of a tapered-GR (Kagan) fit over trailing 730 d.
Src: on disk.
Hunch: K-023's finite-size question in covariate form; if the corner moves, the system's effective size
is moving, which is an order-parameter statement.
Pit: corner magnitude is estimated from the largest few events, so its error is enormous at any window
length the miner can use. This is the entry most likely to produce a number with no information in it.
S15: **UNMEASURABLE at global M>=4.5 over 730 d windows.** Include as a declared UNMEASURABLE so nobody
proposes it fresh later.
Price: 31, recommendation 0.

**F4-11 fraction_above_M5** | READY
Def: trailing 90 d share of events with M >= 5.0.
Src: on disk.
Hunch: a nonparametric b-proxy immune to Mc estimation, which is the largest error source in the b
estimator; if it agrees with b it validates b, if it disagrees it indicts Mc.
Pit: near-perfectly collinear with b by construction; the value is entirely as a robustness check.
S15: MEASURABLE.
Price: 31.

**F4-12 fraction_above_M6** | READY
Def: trailing 180 d share of events with M >= 6.0.
Src: on disk.
Hunch: the tail version; the tail is where the hazard is and where the count is smallest.
Pit: ~90-110 M>=6 events per 180 d globally is not true; the real number is roughly 6-8, so this feature
is mostly zeros. A mostly-zero feature under a Gaussian standardisation is a spike detector.
S15: UNMEASURABLE as a smooth covariate; MEASURABLE only as an indicator (F4-13).
Price: 31, recommendation: run F4-13 instead.

**F4-13 days_since_last_global_M7** | READY
Def: days since the last global M>=7.0, capped at 365.
Src: on disk.
Hunch: **K-043's frame-break in its cheapest form** - every global M7 is a free active-source experiment,
and "time since the last probe" is the state variable that says how recently the planet was interrogated.
Pit: strongly aliased with global aftershock sequences (an M7 is followed by its own aftershocks at
M>=4.5, which is the target). **ETAS baseline mandatory and even then the residual aftershock leakage is
the dominant concern.** This is the single most leakage-prone entry in the catalog and I flag it as such.
S15: MEASURABLE.
Price: 31.

**F4-14 days_since_last_global_M8** | READY
Def: as F4-13 at M>=8.0, capped at 1825.
Src: on disk.
Hunch: K-074's planetary-probe logic; M8 events are the only sources whose surface waves are large
enough that a global dynamic-triggering response is documented rather than hypothesised.
Pit: ~4-6 M8 events in the window: effective N is tiny and the feature is nearly a step function.
S15: UNMEASURABLE as a rate covariate. Usable only inside the conditional designs of family F10.
Price: 31, recommendation 0 standalone.

## 4B. Inter-event time statistics

**F4-15 mean_interevent_time_30d** | READY
Def: mean of the global inter-event times over the trailing 30 d, in hours.
Src: on disk.
Hunch: the reciprocal of rate, entered so the family has an explicit rate control; every other entry here
must be shown not to be this.
Pit: exactly 1/rate, so it is the definitional control and not a hypothesis.
S15: MEASURABLE.
Price: 31.

**F4-16 cv_interevent_time_30d** | READY
Def: coefficient of variation of trailing 30 d inter-event times.
Src: on disk.
Hunch: **the classic burstiness statistic.** CV = 1 is Poisson, CV > 1 is clustered, CV < 1 is
quasi-periodic. This is a shape statistic that is invariant to the rate, which is exactly the property
that makes it immune to the completeness drift that ruins rate features.
Pit: CV of a mixture of independent regional Poisson processes tends to 1 by superposition, so the global
CV is a weak instrument by construction. **This argues for the regional (2R-df) form, not the global.**
S15: global MEASURABLE but weak; regional preferred.
Price: 31.

**F4-17 burstiness_parameter** | READY
Def: `B = (sigma - mu)/(sigma + mu)` of trailing inter-event times (Goh and Barabasi).
Src: on disk.
Hunch: a bounded, better-behaved reparameterisation of CV; the bounded range fixes the standardisation
problem CV has.
Pit: monotone in CV, so it is one hypothesis in two coordinates; enter one.
S15: MEASURABLE.
Price: 31.

**F4-18 memory_coefficient** | READY
Def: lag-1 correlation of consecutive inter-event times (Goh and Barabasi's M).
Src: on disk.
Hunch: **burstiness and memory are independent axes**, and the program has never computed the memory
axis at all; a process can be bursty without memory (Poisson clusters) or memoryful without burstiness.
Pit: contaminated by the superposition of regions again; and it is a second-order statistic, so its own
sampling variance is large at any window the miner can use.
S15: MARGINAL at 30 d; MEASURABLE at 365 d.
Price: 31 per window.

**F4-19 interevent_time_shape_parameter** | READY
Def: fitted gamma shape parameter of trailing 90 d inter-event times.
Src: on disk.
Hunch: the unified-scaling-law form; the shape parameter is the standard summary of the departure from
exponentiality and it is the natural-time-adjacent quantity.
Pit: near-degenerate with CV (shape = 1/CV^2 for a gamma). One of F4-16, F4-17, F4-19 should be run, not
three; running three is paying triple for one idea and it is exactly the kind of thing this catalog
exists to prevent.
S15: MEASURABLE.
Price: 31, and I recommend running only one of the three.

**F4-20 event_count_30d** | READY
Def: trailing 30 d global event count.
Src: on disk.
Hunch: the master control feature for the whole family; **almost everything else in family 4 is
correlated with it**, and any result that does not survive conditioning on it is a count result.
Pit: under a climatology baseline this scores enormous fake skill (SPEC decision 2); under ETAS it is
mostly absorbed. Its purpose is to be conditioned on, not to be tested.
S15: MEASURABLE.
Price: 31.

**F4-21 count_acceleration** | READY
Def: (count in trailing 15 d) minus (count in the 15 d before that), normalised.
Src: on disk.
Hunch: **accelerating moment release rebuilt as a rate feature** (K-020's live remnant); AMR as a
precursor claim is heavily contested, but AMR as a covariate in a properly-baselined miner is a
different and cheaper question.
Pit: AMR's literature is a graveyard of retrospective fits; and the differencing amplifies the aftershock
signal. ETAS baseline mandatory; and the honest framing is "does it add bits over ETAS", nothing more.
S15: MEASURABLE.
Price: 31.

**F4-22 count_jerk** | READY
Def: second difference of the 15 d count series.
Src: on disk.
Hunch: the curvature version; entered because the acceleration hypothesis has never been separated from
the curvature hypothesis in this program.
Pit: differencing twice makes it noise-dominated at 15 d; use a smoothed derivative or declare it a
noise control.
S15: MARGINAL.
Price: 31.

## 4C. Depth structure and migration

**F4-23 depth_percentile_10** and **F4-24 depth_percentile_90** | READY
Def: 10th and 90th percentiles of trailing 30 d hypocentral depths.
Src: on disk.
Hunch: the existing mean depth conflates a shift of the whole distribution with a change in its shape;
percentiles separate them, and the tails are where subduction-versus-crustal composition shows.
Pit: global depth is bimodal (crustal plus slab); percentiles of a bimodal distribution move for
composition reasons. And catalog depths are frequently fixed at 10/33/35 km, which puts atoms of
probability mass at fixed values and breaks percentile continuity.
S15: MEASURABLE.
Price: 62.

**F4-25 depth_migration_rate** | READY
Def: slope of a robust regression of event depth on time over the trailing 60 d, km/d.
Src: on disk.
Hunch: **fluid migration signatures are depth-migration signatures**; a global aggregate of this is
nearly meaningless, but a regional 2R-df sum of it is a fluid-transport instrument.
Pit: dominated by which regions were active, i.e. composition again. Global form is UNINTERPRETABLE and
I would run only the regional form.
S15: regional MEASURABLE under 2R-df.
Price: 31 under 2R-df.

**F4-26 depth_dispersion_30d** | READY
Def: standard deviation of trailing 30 d depths.
Src: on disk.
Hunch: a localisation order parameter in the vertical: if activity concentrates in depth before something
happens, dispersion falls, and that is a one-dimensional version of K-026's spatial entropy.
Pit: composition again; and the fixed-depth convention inflates the concentration artificially.
S15: MEASURABLE.
Price: 31.

**F4-27 shallow_fraction_below_10km** | READY
Def: trailing 30 d share of events shallower than 10 km.
Src: on disk.
Hunch: the very shallow population is where fluids, induced seismicity, and the tidal response are all
strongest; separating it from the existing deep fraction gives the depth axis three bins instead of two.
Pit: the 10 km fixed-depth convention lands exactly on this threshold and will dominate the feature.
**Use 8 km or 12 km instead and say why**, or the feature is a reporting-convention detector.
S15: MEASURABLE.
Price: 31.

## 4D. Spatial structure, correlation length, entropy, percolation

Everything in this sub-block is a **system-level order parameter**: a variable the individual events are
samples of, with its own dynamics. This is the sub-block I would spend the tranche on.

**F4-28 spatial_correlation_length_30d** | READY
Def: the e-folding distance of the two-point correlation function of trailing 30 d epicentres, computed
as a single global scalar (km).
Src: on disk.
Hunch: **K-020's xi(t), rebuilt as a covariate rather than a precursor.** A growing correlation length is
the canonical approach-to-criticality signature in every other many-body system, and it is a scalar with
its own dynamics whether or not it predicts anything.
Pit: the two-point function of a global catalog is dominated by plate-boundary geometry, which is fixed;
the informative quantity is the **deviation from the long-run geometry**, not the raw value. Compute
against a 5-year reference and declare the reference window (exploration-only).
S15: MEASURABLE as a scalar; its own estimation noise at 30 d is the binding constraint and should be
reported by bootstrap before the run.
Price: 31.

**F4-29 correlation_length_growth_rate** | READY
Def: d/dt of F4-28 over 15 d.
Src: on disk.
Hunch: criticality theory predicts divergence, i.e. a positive growth rate, not a large value; the rate
is the theory's actual observable and the level is not.
Pit: derivative of a noisy estimator.
S15: MARGINAL.
Price: 31.

**F4-30 spatial_entropy_of_activity** | READY
Def: Shannon entropy of the trailing 30 d event counts across the 1 deg active cells, normalised by
log(n_active).
Src: on disk.
Hunch: **K-026's localisation order parameter.** If the system localises before it breaks, entropy falls;
this is measurable globally with the full event count and needs no per-cell power at all, which is
exactly the property §P6-4 Finding B says we need.
Pit: entropy of a count field is biased by the total count (a Miller-Madow or Chao-Shen correction is
mandatory), and the count is the target variable. **Without bias correction this feature IS the rate
feature.** That is the whole game for this entry.
S15: MEASURABLE, and it is one of the very few spatially-resolved ideas that is measurable at M>=4.5
because it aggregates rather than tests per cell.
Price: 31.

**F4-31 spatial_entropy_change** | READY
Def: 30 d change in F4-30.
Src: on disk.
Hunch: as F4-29, the theory's observable is the direction of change.
Pit: as F4-30, plus differencing noise.
S15: MEASURABLE.
Price: 31.

**F4-32 participation_ratio** | READY
Def: `(sum_i n_i)^2 / sum_i n_i^2` over cells, trailing 30 d: the effective number of active cells.
Src: on disk.
Hunch: an entropy alternative with much better small-sample behaviour and a direct interpretation as
"how many places are active"; in condensed matter this is the standard localisation diagnostic.
Pit: same count-bias problem as entropy but smaller; still needs the control on F4-20.
S15: MEASURABLE.
Price: 31.

**F4-33 largest_connected_cluster_size** | READY
Def: size (in cells) of the largest connected component of cells active in the trailing 30 d, 8-connectivity
on the 1 deg grid.
Src: on disk.
Hunch: **K-022's percolation order parameter**, which the ledger records as the highest-upside entry in
its family. Percolation of stressed patches bounding rupture size is a real prediction from a real theory
and it has a single scalar observable.
Pit: connectivity on a 1 deg grid at M>=4.5 is dominated by plate-boundary linearity: subduction zones
are already connected strings, so the largest cluster is a geography measurement first and a physics
measurement second. **The deviation-from-reference form is mandatory**, and the reference must be
exploration-window-only.
S15: MEASURABLE as a global scalar.
Price: 31.

**F4-34 cluster_size_distribution_exponent** | READY
Def: fitted power-law exponent of the connected-cluster-size distribution, trailing 90 d.
Src: on disk.
Hunch: at a percolation threshold the cluster-size exponent takes a universal value; measuring where the
crust sits relative to that value is a **thermodynamic statement about the system**, independent of any
forecast.
Pit: power-law fitting on ~20-50 clusters is unreliable; use the Clauset-Shalizi-Newman MLE with the
declared xmin and report the fit's own uncertainty. Do not eyeball it.
S15: MARGINAL.
Price: 31.

**F4-35 percolation_distance_threshold** | READY
Def: the linkage distance at which the trailing 30 d epicentre set percolates (single-linkage clustering
threshold at which the giant component appears), km.
Src: on disk.
Hunch: the **continuum** version of F4-33 that does not depend on the arbitrary 1 deg grid; it is the
natural length scale of the current activity field, measured rather than assumed.
Pit: computationally heavier (single-linkage on ~400 points per 30 d window is fine); still
geography-dominated, so deviation-from-reference again.
S15: MEASURABLE.
Price: 31.

**F4-36 fractal_dimension_D2** | READY
Def: correlation dimension of trailing 90 d epicentres (Grassberger-Procaccia).
Src: on disk.
Hunch: D2 is the standard measure of how the seismicity field fills space; changes in D2 are the
renormalisation question in K-025 asked as a scalar.
Pit: D2 estimation needs a declared scaling range, and the scaling range is a hidden parameter that
notoriously drives the answer. Declare it once per S-9.
S15: MEASURABLE at 90 d.
Price: 31.

**F4-37 nearest_neighbour_distance_median** | READY
Def: median of the space-time nearest-neighbour distance (Zaliapin-Ben-Zion eta), trailing 90 d.
Src: on disk.
Hunch: the Zaliapin metric is the field's best-validated declustering statistic and its **bimodality is a
state variable**: the ratio of the clustered mode to the background mode is a direct measurement of how
much of the current catalog is aftershocks.
Pit: the metric embeds b and the fractal dimension as parameters; those must be fixed a priori (from the
exploration window) or the feature is circular with F4-01 and F4-36.
S15: MEASURABLE.
Price: 31.

**F4-38 clustered_fraction** | READY
Def: fraction of trailing 90 d events assigned to the clustered mode by the Zaliapin threshold.
Src: on disk.
Hunch: **the single most interpretable state variable in this family**: what fraction of what is
happening is aftershocks. If the background fraction moves, mu moves, and mu is the quantity K-010 wants
to make a latent state.
Pit: threshold is fixed a priori (see F4-37); and under an ETAS baseline much of this is absorbed, which
is the correct outcome and should be reported as such rather than treated as a loss.
S15: MEASURABLE.
Price: 31.

## 4E. Productivity, branching, and ETAS-residual features

**F4-39 branching_ratio_n_90d** | READY
Def: trailing 90 d estimate of the ETAS branching ratio n (fraction of events that are triggered).
Src: on disk, computed from the existing ETAS fit machinery.
Hunch: **K-018's distance-to-criticality order parameter, and the head of the whole order-parameter
programme.** n = 1 is the critical branching point; the distance 1 - n is the system's distance to
runaway, and it is estimable from data the program already has.
Pit: n and the background rate mu trade off almost exactly in a short-window fit, so a moving n may be a
moving mu. **The pair must be estimated jointly and reported as a pair**; reporting n alone is the
classic error here.
S15: MEASURABLE at 90 d globally; per region MARGINAL.
Price: 31.

**F4-40 background_rate_mu_90d** | READY
Def: the joint partner of F4-39.
Src: on disk.
Hunch: K-010's latent-state mu made into an observable; if mu has dynamics, the "weather" K-009 found in
the residuals has a name.
Pit: as F4-39.
S15: MEASURABLE.
Price: 31.

**F4-41 aftershock_productivity_alpha_drift** | READY
Def: rolling 365 d estimate of the ETAS productivity exponent alpha.
Src: on disk.
Hunch: **drift of the generator's parameters is a meta-pattern**: if alpha moves, the universal shape
that the program validated is only locally universal, and that is a finding about the validated baseline
rather than about any new signal.
Pit: alpha is poorly constrained without large events in the window; its rolling estimate will be
dominated by whether an M7+ happened to occur. **Report the estimate's own CI and treat any excursion
inside it as UNRESOLVED**, per the S-17 candidate.
S15: MARGINAL.
Price: 31.

**F4-42 omori_p_drift** | READY
Def: rolling 365 d estimate of the Omori p exponent.
Src: on disk.
Hunch: as F4-41, on the decay side; p and alpha drifting in opposite directions would be a
composition artifact, drifting together would be physics.
Pit: as F4-41.
S15: MARGINAL.
Price: 31.

**F4-43 etas_residual_mean_30d** | READY
Def: trailing 30 d mean of the daily Poisson residual `(N_obs - Lambda)/sqrt(Lambda)`.
Src: on disk (the engine computes Lambda already).
Hunch: **K-009's "weather" as a covariate**: if the residuals are not white, then the recent residual
predicts the next residual, and that is the cheapest possible assimilation toehold.
Pit: this is close to circular - it is a lagged version of the target normalised by the model. Under the
ETAS baseline it tests exactly the residual autocorrelation K-009 already measured, which means **the
honest framing is "does the mining pipeline reproduce K-009", i.e. a G-M1-class pipeline control**, and
that is genuinely useful.
S15: MEASURABLE.
Price: 31.

**F4-44 etas_residual_variance_30d** | READY
Def: trailing 30 d variance of the same residual.
Src: on disk.
Hunch: **critical slowing down predicts rising variance before a transition** (K-019); variance is the
second moment and no first-moment statistic in the current session can see it.
Pit: residual variance is inflated by any model misspecification, and ETAS is misspecified everywhere;
the entry measures misspecification unless it is differenced against a long-run baseline.
S15: MEASURABLE.
Price: 31.

**F4-45 etas_residual_lag1_autocorr_90d** | READY
Def: trailing 90 d lag-1 autocorrelation of the daily residual.
Src: on disk.
Hunch: the other half of critical slowing down: rising autocorrelation. This is the canonical
early-warning-signal pair from ecology and climate (Scheffer et al.), imported wholesale.
Pit: the early-warning-signals literature has a serious record of false positives under non-stationary
noise, and our noise is very non-stationary. The known failure mode must be stated in the entry.
S15: MEASURABLE at 90 d.
Price: 31.

**F4-46 etas_residual_skewness_90d** | READY
Def: trailing 90 d skewness of the residual.
Src: on disk.
Hunch: flickering between states produces skewness before a transition; the third moment is unexplored.
Pit: third moments need large samples; 90 daily values is not large.
S15: MARGINAL.
Price: 31.

**F4-47 release_deficit_D_global** | READY
Def: `D = Lambda_ETAS - N_obs` summed over the domain, trailing 90 d.
Src: on disk.
Hunch: **K-076's deficit field, aggregated to the level where it is measurable.** The per-cell version is
the thing §P6-4 Finding B kills; the global sum keeps the full event count and is a legitimate scalar.
Pit: the classical quiescence observable is owned (M-007) and prospectively decided; this must be framed
as an instrument and a covariate, never as a rediscovery of quiescence. And it is the negative of F4-43
up to normalisation, so run one.
S15: MEASURABLE globally.
Price: 31.

**F4-48 release_deficit_regional_2Rdf** | READY
Def: the per-region deficit vector entering the 2R-df phase-incoherent sum.
Src: on disk, exploration-window region definition.
Hunch: the deficit's whole point is that it is regional, and the 2R-df instrument is precisely the way to
test a regional quantity without paying per-region multiplicity or quoting per-region amplitudes.
Pit: **the region set must be exploration-window-only (Rule 4.1)** and must not descend from the K-080
census list. Per-region amplitudes are UNRESOLVED per §0.3 and must not be tabulated.
S15: the summed statistic is MEASURABLE; per-region values are not.
Price: 1 per (feature, lag) = 31 for the lag scan.

## 4F. Quiescence variants that respect the spent hash

**Standing constraint.** The 70/30 quiescence holdout hash `1afa6cdc...` was spent 2026-08-10. Nothing
below spends a holdout hash; everything is exploration-window mining, and every entry is a covariate or
a second-moment object rather than a re-run of the classical precursor claim that M-007 shows is owned
and prospectively decided.

**F4-49 quiescence_duration_in_natural_time** | READY
Def: for the global aggregate, the number of **expected events** (integrated Lambda) accumulated since
the observed rate last exceeded its expectation.
Src: on disk.
Hunch: **K-079's claim that quiescence must be measured in accumulated strain, not wall time**, in the
cheapest available currency: expected-event-count is the catalog's own natural clock.
Pit: the clock change is the hypothesis; if it does nothing the correct conclusion is that wall time was
fine, which is a real result. Do not present the clock change as a given.
S15: MEASURABLE.
Price: 31.

**F4-50 quiescence_area_fraction** | READY
Def: fraction of the MEASURABLE cell set (per S-15, declared floor) currently in deficit.
Src: on disk.
Hunch: **the aggregate of the per-cell battery without the per-cell battery**: a fraction over 1,000
cells has the counting power the individual cells lack, exactly the aggregation move that Rule 4.2 makes
for phase.
Pit: the MEASURABLE set must be defined on the exploration window (Rule 4.1); and this feature is
strongly anticorrelated with the global count, so F4-20 conditioning is mandatory.
S15: MEASURABLE as an aggregate. **This entry is the honest survivor of the dead per-cell battery and
I would put it in the first tranche.**
Price: 31.

**F4-51 quiescence_persistence** | READY
Def: mean run-length of consecutive deficit days in the global aggregate.
Src: on disk.
Hunch: a second-order property of the deficit; persistence distinguishes a slow state change from noise
in a way the level cannot.
Pit: run-length statistics are extremely sensitive to the threshold; declare one.
S15: MEASURABLE.
Price: 31.

**F4-52 quiescence_spatial_clumping** | READY
Def: Moran's I of the deficit field over the exploration-defined active cells.
Src: on disk.
Hunch: a spatially coherent deficit is a physical statement; a scattered one is noise. Moran's I
separates them with one number and no per-cell claims.
Pit: spatial autocorrelation statistics need a declared weight matrix; and the plate-boundary geometry
gives a large baseline I that must be differenced away.
S15: MEASURABLE as an aggregate.
Price: 31.

**F4-53 anti_quiescence_surplus** | READY
Def: the positive half of the deficit, i.e. `max(N_obs - Lambda, 0)` aggregated.
Src: on disk.
Hunch: the **negative space of the negative space**: quiescence research always looks at the low tail;
the high tail is the same statistic and nobody reports it. If only the low tail carries signal, that
asymmetry is itself the finding.
Pit: the high tail is where aftershock leakage lives; ETAS baseline mandatory.
S15: MEASURABLE.
Price: 31.

## 4G. Composition, observer, and self-audit features

**F4-54 catalog_completeness_Mc_proxy** | READY
Def: trailing 90 d Mc estimated by maximum-curvature on the global magnitude histogram.
Src: on disk.
Hunch: **K-015's invert-the-observer move**: Mc is a measurement of the network, and network changes are
the largest secular signal in any catalog statistic. Entering Mc explicitly lets every other feature be
conditioned on it.
Pit: max-curvature Mc is biased low and jumpy; a goodness-of-fit or b-stability estimator is better.
And at a global aggregate Mc is a mixture of regional Mc's, which is not a physical quantity.
S15: MEASURABLE as a control, not as a hypothesis.
Price: 31.

**F4-55 reporting_latency** | READY
Def: median (catalog insert time minus origin time) where available.
Src: ComCat metadata; **partially on disk**, may need re-download with the extra fields.
Hunch: a direct measurement of the observer's duty cycle; it should show the weekly and diurnal
operational cycles, and if any seismic feature correlates with it, that feature is an observer artifact.
Pit: not present in the current CSVs; needs a re-fetch. Historical latency is not always retrievable.
S15: MEASURABLE if fetched.
Price: 31.

**F4-56 day_of_week_indicator** | READY
Def: 7-level categorical, or a 7 d phase feature.
Src: closed form.
Hunch: **there is no physical 7 d cycle**, so this is the purest available observer control: any weekly
signal is human. Given that the last session's best hit was a solar feature at lag 3 d, having a
calibrated human-cycle line in the same run is worth its price several times over.
Pit: 7.0 d has a sinc factor of ~1 under day binning, so unlike the diurnal artifacts this one is fully
visible. That is the point.
S15: MEASURABLE, A_min 5.5%.
Price: 31 (or 1 as a phase feature at fixed period).

**F4-57 time_of_day_UTC_distribution_skew** | READY
Def: circular mean and resultant length of event times-of-day, trailing 90 d.
Src: on disk (`day_float` already exists in `load_event_marks`).
Hunch: the diurnal detection cycle is real and measurable, and its **drift over the record** is a
completeness time series that costs nothing to compute; this is K-091's positive control turned into a
covariate.
Pit: it lives exactly in the notched band, which is why it must be computed from `day_float` and not from
daily bins. It is currently computed nowhere.
S15: MEASURABLE.
Price: 31.

**F4-58 clustering_deflation_measurement** | READY
Def: not a feature: the direct measurement of the factor between `3.96/sqrt(N)` and the miner's reported
~5.6% bound, by planting a sinusoid of known amplitude at the fortnightly band and measuring recovery,
per G-M1.
Src: on disk, existing planted-signal harness.
Hunch: **§0.3's honest gap, closed.** Every S-15 declaration in this catalog depends on this factor, and
right now it is my inference from two numbers rather than a measurement.
Pit: none that I can see; this is a pipeline measurement with a known answer shape. If the factor turns
out to be BH overhead rather than clustering, every measurability line in this file shifts and should be
rewritten before the tranche runs.
S15: not applicable; it IS the S-15 calibration.
Price: 0 priced tests (it is a control, not a hypothesis). **I would run this before anything else in
this catalog.**

**F4-59 magnitude_rounding_indicator** | READY
Def: trailing 90 d fraction of magnitudes reported at exactly one decimal ending in 0 or 5.
Src: on disk.
Hunch: a pure reporting-convention detector; magnitude rounding changes when agencies change practice,
and such changes are step functions that will manufacture apparent structure in b and in every
magnitude-derived feature.
Pit: none; it is a diagnostic. Its only risk is being mistaken for physics.
S15: MEASURABLE.
Price: 31, and its real use is as a covariate to condition on.

**F4-60 network_contributor_mix** | READY
Def: trailing 90 d Herfindahl index of the contributing-network field in ComCat.
Src: on disk (`net` column).
Hunch: composition drift made explicit; if the mix of reporting agencies changes, every global statistic
changes with it, and this is the one-number summary of that.
Pit: as F4-59.
S15: MEASURABLE.
Price: 31.

**F4-61 moment_release_rate_30d** | READY
Def: trailing 30 d summed seismic moment (from magnitude), in N m.
Src: on disk.
Hunch: moment, not count, is the physically conserved-ish quantity; a count feature and a moment feature
answer different questions and only the count version exists today.
Pit: dominated by the single largest event in the window, so it is nearly `max_magnitude` (F4-08) in
disguise; report the correlation. Heavy-tailed, so use log-moment.
S15: MEASURABLE.
Price: 31.

**F4-62 moment_to_count_ratio** | READY
Def: log(moment_30d) minus a b-implied expectation from count_30d.
Src: on disk.
Hunch: the residual of moment given count is a **magnitude-distribution anomaly detector** that does not
require estimating b at all.
Pit: ratio instability again (the TSI-ratio corpse's failure mode); use the residual form, not the ratio
form, and say so.
S15: MEASURABLE.
Price: 31.

**F4-63 largest_event_share_of_moment** | READY
Def: share of trailing 90 d moment contributed by the single largest event.
Src: on disk.
Hunch: a concentration statistic; high share means one event dominates, which is a distinct dynamical
regime from distributed release, and this is the moment-space analogue of the magnitude gap (F4-07).
Pit: heavy-tailed and near-binary in practice.
S15: MEASURABLE.
Price: 31.

**F4-64 foreshock_fraction_proxy** | READY
Def: trailing 365 d fraction of M>=5.5 events preceded within 3 d and 50 km by an M>=4.5.
Src: on disk.
Hunch: **K-014's time-reversal asymmetry as a slowly-varying state variable**: if the foreshock rate
moves, the amount of causal information available to any forecaster moves with it.
Pit: near-tautologically related to the clustering fraction (F4-38); and the 3 d / 50 km choice is a
hidden 2-parameter grid. Declare one pair.
S15: MARGINAL (M>=5.5 counts over 365 d are modest).
Price: 31.

**F4-65 swarm_index** | READY
Def: trailing 30 d fraction of clusters whose largest event occurs in the middle third of the cluster's
time span (a swarm-versus-mainshock-aftershock discriminator).
Src: on disk.
Hunch: swarms are fluid-driven and mainshock sequences are stress-driven; the mix of the two is a
statement about which process currently dominates, and it is a global scalar.
Pit: cluster identification is the Zaliapin machinery again, with its parameters fixed a priori.
S15: MEASURABLE.
Price: 31.

---

# FAMILY 5 - CROSS-REGIONAL: THE 2R-df INSTRUMENT'S NATURAL EXTENSIONS

Rule 4.2 gives the program a phase-incoherent regional sum: `R` pre-declared regions, per-region 2-df
score statistics, summed to a 2R-df quadratic form, one test per (feature, lag). It was invented to kill
the phase-cancellation blind spot. **It generalises far past that**, and this family is the generalisation.

**Standing constraints for the whole family.** Region partitions are exploration-window-only (Rule 4.1),
declared once with no alternatives run (S-9), frozen in the config hash (Rule 4.5). Per-region amplitudes
are UNRESOLVED and may not be tabulated (my §0.3 reading of Rule 4.3). The K-080 census list may appear
only as a cross-check overlap fraction.

## 5A. Partitions (the choice IS a hypothesis)

**F5-01 partition_longitude_sectors_R18** | READY
Def: 18 sectors of 20 deg longitude.
Src: on disk.
Hunch: the correct partition for anything whose phase is a function of local time (all diurnal and
semidiurnal tides, TEC, the detection cycle); longitude is exactly the variable the global sum cancels.
Pit: sectors have wildly unequal event counts (the Pacific rim dominates); the 2R-df sum handles this
correctly but the per-region interpretation does not exist.
S15: summed statistic MEASURABLE.
Price: 1 partition x per-feature cost.

**F5-02 partition_tectonic_province_R12** | READY
Def: 12 provinces by Flinn-Engdahl or by a declared plate-boundary-type map (subduction, continental
transform, oceanic ridge, intraplate).
Src: on disk (region files in `data/comcat_world` are a natural starting partition: Alaska-Aleutians,
California, Caribbean, Chile, Greece-Aegean, Himalaya, Iceland, Indonesia, Iran, Japan, Mexico,
Philippines, Turkey = 13 regions already downloaded and separated).
Hunch: the correct partition for anything whose response depends on fault mechanism or thermal state.
**And the 13 region files already on disk are a pre-existing, non-post-hoc partition**, which is a
genuine asset: they were defined before this question was asked.
Pit: **fault-type parameter pools are a corpse** (fault-type stratification did not help ETAS
parameters). This partition is not that: it is a phase-coherence partition for a covariate test, not a
parameter pool. The distinction must be stated in the entry or a reviewer will conflate them.
S15: summed MEASURABLE.
Price: 1 partition.

**F5-03 partition_latitude_bands_R8** | READY
Def: 8 latitude bands.
Src: on disk.
Hunch: the correct partition for anything seasonal or hemispheric (hydrology, atmospheric loading,
polar motion), where the northern and southern signals are antiphase and cancel exactly in the global sum.
Pit: latitude correlates with tectonic setting; confounded with F5-02 by geography.
S15: summed MEASURABLE.
Price: 1 partition.

**F5-04 partition_by_depth_class_R3** | READY
Def: shallow (<35 km), intermediate (35-300 km), deep (>300 km) as three "regions" in the 2R-df sum.
Src: on disk.
Hunch: **the partition axis need not be spatial at all.** Depth classes have different physics and
plausibly different response phases; the 2R-df machinery does not care whether the partition is
geographic.
Pit: K-078's depth-class entry is already ADMIT-RESCOPED and bound-producing unless n >= 40; this is a
different statistic on the same axis and must not be presented as independent evidence of the same claim.
S15: summed MEASURABLE; deep class alone is thin.
Price: 1 partition.

**F5-05 partition_by_local_rate_tercile** | READY
Def: three "regions" defined by whether a cell's exploration-window mean rate is in the low, middle, or
high tercile.
Src: on disk, exploration-window only.
Hunch: a **state-based rather than place-based partition** - the "regions" are dynamical classes. If the
response depends on how loaded or how active a place is, this partition sees it and a geographic one
averages it away.
Pit: **this is the closest thing in the catalog to the holdout-contamination failure Rule 4.1 names.**
The tercile assignment must be computed on `[365, 8081)` and frozen; computing it on the full record is
exactly the K-080 mistake. Say so in the config comment, not just here.
S15: summed MEASURABLE.
Price: 1 partition.

**F5-06 partition_by_ledger_class** | READY
Def: regions defined by the stress-ledger class (loading-and-firing versus loading-and-silent) from the
existing B-4 work.
Src: on disk (SoCal only).
Hunch: K-071/K-077's interaction claim in partition form: if the silent cells respond differently, the
2R-df sum with this partition is the cheapest possible test of it.
Pit: SoCal-only, so the event count is a small fraction of global; and the ledger classification has its
own provenance and its own degeneracy (M-002 records that the degeneracy is CONTESTED).
S15: SoCal M>=4.5 counts are far too small; **this entry is UNMEASURABLE at the global mining floor and
is honest only at SoCal's low Mc**, i.e. in a different instrument.
Price: 1 partition, contingent on a SoCal low-Mc mining mode existing.

## 5B. Lead-lag and coupling between regions

**F5-07 region_pair_lagged_cross_correlation** | READY
Def: for each ordered region pair (i, j) and lag L in 1..30 d, the correlation between region i's daily
residual and region j's residual L days later.
Src: on disk.
Hunch: **teleconnection (K-016) with the direction of the arrow retained.** If any region systematically
leads another, that is either a physical coupling or a shared driver, and either one is a finding.
Pit: **the multiplicity price is brutal and must be stated up front: R(R-1) ordered pairs x 30 lags. At
R = 13 that is 156 x 30 = 4,680 tests for one feature.** Under BH at q = 0.1 you expect 468 false
survivors if nothing is there. The honest design is not the full battery: it is a **single summary
statistic** (see F5-08).
S15: each pair's residual series has ~7716 daily values but the effective count is set by the smaller
region; MARGINAL per pair, and the per-pair table should not be published.
Price: 4,680 as a battery; 1 as the summary. Run the summary.

**F5-08 global_coupling_summary_statistic** | READY
Def: the sum over ordered pairs of squared lagged cross-correlations at each lag, i.e. a single
phase-incoherent coupling statistic per lag, with a null from region-wise circular shifts that preserve
each region's own autocorrelation.
Src: on disk.
Hunch: **the 2R-df move applied to coupling instead of to phase.** The blind spot Rule 4.2 kills for
phase (cancellation in the sum) has an exact analogue for coupling (dilution in the battery), and the
same cure works: sum the squares, pay one test.
Pit: the null must destroy cross-region coupling while preserving within-region autocorrelation, which
means independent circular shifts per region, and the shift must exceed the longest within-region
correlation time. Declare that time from data.
S15: MEASURABLE as a single statistic.
Price: 31 (one statistic x 31 lags).
**This is the entry I would defend hardest in family 5.**

**F5-09 coupling_asymmetry** | READY
Def: the antisymmetric part of the lagged coupling matrix summed over pairs: does information flow
preferentially in one direction?
Src: on disk.
Hunch: a symmetric coupling means a shared driver; an antisymmetric one means propagation. **The
antisymmetry separates "common cause" from "cause" without needing a mechanism**, which is a rare and
valuable property.
Pit: an apparent asymmetry can be produced purely by different autocorrelation times in the two regions;
the null must match autocorrelation per region, and a whitened version should be run alongside.
S15: MEASURABLE as a summary.
Price: 31.

**F5-10 transfer_entropy_summary** | READY
Def: sum of pairwise transfer entropies between regional daily residual series, binned to 3 states.
Src: on disk.
Hunch: the nonlinear version of F5-08; if coupling exists but is not linear, correlation misses it.
Pit: transfer entropy needs far more data than we have to be unbiased at 3 states x 30 lags; the bias
correction is essential and the estimator has a known positive bias that manufactures signal. Under
Mignan-Broccardo discipline, **if transfer entropy finds something that the linear statistic cannot see
after conditioning, treat it as estimator bias until proven otherwise.**
S15: MARGINAL to UNMEASURABLE; the bias is the binding constraint, not the count.
Price: 31.

**F5-11 region_synchrony_kuramoto** | READY
Def: Kuramoto order parameter over the per-region phases of the dominant fortnightly component of each
region's residual series.
Src: on disk.
Hunch: **K-024's synchronisation question, made measurable.** The ledger records the standalone claim as
REJECTED but permits the object to be built; as a mining covariate the object is legitimate.
Pit: the ledger's rejection was of the CLAIM, and any output here must not be quoted as reviving it. And
the phase of a noisy residual's fortnightly component is mostly noise at regional counts.
S15: MARGINAL.
Price: 31.

**F5-12 regional_rate_common_mode** | READY
Def: first principal component of the R regional daily rate series (fitted on exploration only), entered
as a covariate.
Src: on disk.
Hunch: **K-016R's common mode**, and it is the single most likely place for a global observer artifact
to be visible as an actual signal: a worldwide reporting or processing change moves all regions together.
Pit: PC1 of rate series will be dominated by the global count, which is the target. Orthogonalise against
F4-20 first or the feature is the target in a hat.
S15: MEASURABLE.
Price: 31.

**F5-13 common_mode_residual_dispersion** | READY
Def: cross-region dispersion after removing PC1: how differently are regions behaving right now?
Src: on disk.
Hunch: the complement of F5-12, and the more interesting half: dispersion rising means the regions are
decoupling, which is a system-state statement.
Pit: dispersion of Poisson counts is driven by the counts.
S15: MEASURABLE.
Price: 31.

**F5-14 antipodal_pairing** | READY
Def: for each region, the daily residual of its geographic antipode, entered as a covariate at lags 0..3 d.
Src: on disk.
Hunch: **K-062's antipodal focusing** in the cheapest available form; surface-wave energy refocuses near
the antipode and O'Malley 2018 reports it. This is an explicit replication, and the entry says so.
Pit: the antipode of most seismic regions is ocean, so the pairing is sparse; and this is a replication,
not a discovery, per the K-062 ruling.
S15: sparse pairs give small N. MARGINAL.
Price: 4 lags x 1 statistic.

**F5-15 great_circle_distance_from_last_M7_5** | READY
Def: for each region, the great-circle distance to the most recent global M>=7.5, and days since.
Src: on disk.
Hunch: **K-070/K-074's planetary-probe design as a mining covariate**: every large event interrogates
every region, and the region-by-distance-by-delay response surface is the observable.
Pit: dynamic triggering is documented to last minutes to days, which is inside the day bin for the
short end (§K87-0's sinc issue in a different guise); and the distance dependence is confounded with
aftershock-zone proximity. Exclude the source region by a declared radius.
S15: ~50-70 M>=7.5 events in the window across R regions gives a decent stacked count. MEASURABLE as a
stacked (epoch-folded) statistic, not as a per-event one.
Price: 31 (distance-binned stack x lags).

**F5-16 region_readiness_rank** | READY
Def: the ORDER in which regions show a rate increase after a shared global probe (F5-15), as a rank
vector, with the statistic being rank stability across probes (Kendall's W).
Src: on disk.
Hunch: **K-075's readiness leaderboard, promoted in the ledger to a family gate.** The claim is that the
ORDER is a state variable: if the same regions always answer first, that ordering is a physical
property of the crust and not of the probe.
Pit: rank stability with ~50 probes and 13-20 regions is estimable but the per-probe ranks are extremely
noisy; Kendall's W has a known upward bias with tied and missing ranks.
S15: MARGINAL, and I would call it MEASURABLE only after a planted-order recovery test (S-17 candidate).
Price: 1 statistic.

**F5-17 region_specific_lag_to_common_mode** | READY
Def: per-region cross-correlation lag with PC1 (F5-12), summarised by its dispersion.
Src: on disk.
Hunch: if there is a global driver, the regional lags to it map the propagation; if there is no driver,
the lags are noise and their dispersion is a null distribution we can measure directly.
Pit: circular reasoning risk with F5-12; orthogonalise.
S15: MARGINAL.
Price: 31.

**F5-18 pairwise_b_value_coupling** | READY
Def: F5-08's summary statistic computed on regional b-value series instead of residual series.
Src: on disk.
Hunch: **coupling need not be in the rate.** If regions share a stress state, their b-values may move
together while their rates do not, and nobody has looked at coupling in any variable except rate.
Pit: b series are heavily smoothed by their 90 d window, so their cross-correlations are inflated and
their effective sample size is 7716/90 ~ 86, not 7716. **The block length must be 4x the b window
minimum**, and the effective N is the binding constraint.
S15: MARGINAL at 90 d windows; MEASURABLE at 30 d with noisier b.
Price: 31.

**F5-19 pairwise_depth_coupling** | READY
Def: as F5-18 on regional mean-depth series.
Src: on disk.
Hunch: as F5-18, in the depth variable; a shared depth migration across unconnected regions would be
either a fluid-scale impossibility or an observer artifact, and both are worth knowing.
Pit: as F5-18, plus the depth-convention artifacts of F4-27.
S15: MARGINAL.
Price: 31.

**F5-20 leave_one_region_out_stability** | READY
Def: not a feature: the mandated stability audit in which every 2R-df result is recomputed R times, each
time dropping one region, reporting the max shift in the statistic.
Src: on disk.
Hunch: **the 2R-df sum's one real vulnerability is that a single dominant region carries it.** This audit
is the cheapest possible defence and it converts "the sum is significant" into "the sum is significant and
not because of Japan".
Pit: none; it is an audit. It multiplies compute by R, not multiplicity by R (it is a robustness check on
a single declared statistic, not R new tests) - but Popper should confirm that reading, since if it is
treated as R tests the price changes completely.
S15: not applicable.
Price: 0 new tests if ruled an audit; R x if ruled a battery. **Popper's call, flagged.**

---

# FAMILY 6 - ENVIRONMENTAL AND HYDROLOGICAL LOADING

**The ownership line, stated before the entries.** Sirorattanakul and Avouac own California seasonal
hydrological loading (M-010). Nothing in this family is a novelty claim about California seasonal
loading. The live space is: (a) global extension, where the phase varies by region and the global sum
cancels it, so the 2R-df instrument is the natural and possibly the only honest instrument; (b) the
non-seasonal, transient part of hydrology, which the seasonal literature does not cover; (c) hydrology
as a **conditioning variable** for other tests rather than as a driver in its own right.

**F6-01 seasonal_phase_by_hemisphere** | READY
Def: `annual_phase` sign-flipped in the southern hemisphere, entered per region.
Src: closed form plus the region partition.
Hunch: **the single cheapest fix to a known exact cancellation.** Northern and southern seasonal loading
are antiphase; the global sum kills the seasonal signal exactly, the way it kills the diurnal tide. One
sign flip recovers it.
Pit: it is not a new hypothesis, it is a bug fix to the instrument, and it must be described that way.
The seasonal result itself is owned.
S15: MEASURABLE under the 2R-df sum.
Price: 31.

**F6-02 regional_precipitation_anomaly** | DOWNLOAD
Def: daily regional precipitation anomaly.
Src: GPCP or CHIRPS or ERA5 (all free).
Hunch: precipitation is the input to the hydrological load and leads the load by weeks; the lag between
them is estimable and the lag structure is what distinguishes direct poroelastic response from
surface-load response.
Pit: the lag is region-specific (snow versus rain, basin size), so a global lag scan is a smear; and
precipitation is intermittent and heavy-tailed.
S15: MEASURABLE regionally.
Price: 31 under 2R-df.

**F6-03 snow_water_equivalent** | DOWNLOAD
Def: regional SWE anomaly, daily.
Src: GLDAS/NOAH, free.
Hunch: snow is the slowest and largest single-season load, and its melt is a **transient unloading with
a sharp date** - transients are the class the tidal bound does not cover.
Pit: mountainous, so it is confounded with orography and therefore with specific fault systems; and this
is the closest entry to the owned California result.
S15: MEASURABLE regionally where snow exists.
Price: 31 under 2R-df.

**F6-04 snowmelt_onset_transient** | DOWNLOAD
Def: the date of maximum negative SWE derivative per region per year, as an epoch-folding anchor.
Src: GLDAS.
Hunch: **an epoch-folded transient test rather than a periodic one.** A sharp annual transient with a
date that moves year to year is invisible to a fixed-phase annual feature and visible to an epoch fold.
Pit: ~21 events per region: tiny N. This is a stacked test or it is nothing.
S15: UNMEASURABLE per region; MARGINAL stacked across regions.
Price: 1 stacked statistic.

**F6-05 ENSO_ONI_index** | DOWNLOAD
Def: monthly Oceanic Nino Index, interpolated to daily.
Src: NOAA CPC (free).
Hunch: ENSO redistributes water mass, sea level, and atmospheric pressure on interannual timescales, and
there is an existing (contested) literature linking it to seismicity rates. It is the largest interannual
mass-redistribution signal on the planet.
Pit: **monthly resolution over 21 years is ~250 independent samples, and ENSO's own autocorrelation
time is ~6-12 months, so the effective sample size is 20-40.** That is the binding constraint and it makes
this feature a bound-producer at best. The block bootstrap must use blocks of 2+ years, giving ~10 blocks.
S15: **UNMEASURABLE for anything but a large effect.** Declare it as such, run it once, publish the bound.
Price: 31, but the honest price is 1 (the lag scan is meaningless at ENSO's autocorrelation time).

**F6-06 ENSO_rate_of_change** | DOWNLOAD
Def: 3-month difference of ONI.
Src: NOAA CPC.
Hunch: transitions rather than states; the ENSO transition is when mass actually moves.
Pit: as F6-05.
S15: UNMEASURABLE except as a bound.
Price: 1.

**F6-07 PDO_index**, **F6-08 NAO_index**, **F6-09 AMO_index** | DOWNLOAD
Def: standard monthly climate indices.
Src: NOAA/NCAR, free.
Hunch: MECHANISM-FREE beyond the general mass-redistribution argument, and I label them so. Their real
value is as a **set of long-period null lines** with realistic autocorrelation structure, against which
any claimed long-period seismic signal can be sized.
Pit: the multi-decadal ones (PDO, AMO) have 1-2 cycles in the window; claiming anything from them is the
long-period numerology trap.
S15: UNMEASURABLE by autocorrelation time.
Price: 3 (one each, no lag scan).

**F6-10 regional_sea_level_seasonal** | DOWNLOAD
Def: the seasonal (annual plus semiannual) fit to regional sea level, as an ocean-load feature.
Src: AVISO.
Hunch: subduction zones are underwater and the annual sea-level cycle is 10-20 cm there, i.e. 1-2 kPa of
load with a regionally-varying phase.
Pit: exactly antiphase between hemispheres again; global sum cancels.
S15: MEASURABLE under 2R-df.
Price: 31.

**F6-11 river_discharge_major_basins** | DOWNLOAD
Def: daily discharge at gauges on the largest rivers near active margins.
Src: GRDC (free, registration) or national agencies.
Hunch: an in-situ measurement of basin water flux, at daily resolution, which GRACE cannot give.
Pit: gauge selection is a multiplicity axis; damming changes the series inside the window (composition
drift in the covariate).
S15: MARGINAL regionally.
Price: 31 per declared basin, max 3 declared.

**F6-12 atmospheric_tide_S1_S2_pressure** | DOWNLOAD
Def: the amplitude of the S1 and S2 atmospheric pressure tides, regionally.
Src: ERA5.
Hunch: **the atmospheric tide sits at exactly 24.000 and 12.000 h, which is exactly where the day-binning
sinc is exactly zero and where the detection-cycle artifact lives.** It is therefore both the perfect
confound and the perfect calibration line for the sub-daily infrastructure build.
Pit: entering it without the unbinned path (F9-15) is meaningless; entering it with the unbinned path
requires the S1/S2 detection control to be run simultaneously or the two are inseparable.
S15: UNMEASURABLE now; MEASURABLE after F9-15, jointly with K-091's control.
Price: 0 now; 31 after.

**F6-13 permafrost_freeze_thaw_index** | DOWNLOAD
Def: regional accumulated freezing/thawing degree-days.
Src: ERA5 2 m temperature.
Hunch: SPECULATIVE. Freeze-thaw changes near-surface elastic properties and pore pressure in high-latitude
regions (Alaska-Aleutians is a downloaded region on disk); the mechanism is real but the depth reach is
shallow and the amplitude at seismogenic depth is likely negligible.
Pit: perfectly seasonal, so it is collinear with every other seasonal feature; its only distinctive
content is its regional phase, which is why it belongs under 2R-df or nowhere.
S15: MEASURABLE under 2R-df in high-latitude strata only.
Price: 31.

**F6-14 glacial_mass_loss_rate** | DOWNLOAD
Def: regional ice-mass change rate from GRACE mascons (Iceland and Alaska are both downloaded regions).
Src: PO.DAAC.
Hunch: glacial unloading is a documented driver of seismicity in Iceland and Alaska; this is a
**known-effect calibration region**, not a discovery target, and its value is that it gives the loading
response curve a real anchor (the K-045 logic).
Pit: known and owned; monthly resolution; and it is a secular trend more than a variable, so it is nearly
collinear with time in those regions.
S15: regionally MARGINAL.
Price: 31 in declared strata.

**F6-15 hydro_conditioning_stratum** | READY
Def: not a feature: a two-level stratum (wet season / dry season, declared per region from exploration
data) used to condition every other test in the catalog.
Src: on disk plus ERA5.
Hunch: **the conditioning use is the live one.** If tidal or any other small-stress response exists only
in the wet state (K-039's dry-log logic, K-077's interaction shape), then unconditional tests are diluted
by construction and the conditional is the correct instrument.
Pit: conditioning doubles the test count for every feature it touches; the price must be paid explicitly
(see F10 for the class price). And the stratum definition is exploration-only.
S15: splits N in two: A_min rises by sqrt(2) to ~7.8% deflated. MEASURABLE but noticeably weaker, and
that cost must be quoted whenever the conditional is proposed.
Price: multiplies the conditioned feature set by 2.

---

# FAMILY 7 - THE OBSERVER AS PART OF THE SYSTEM

K-015, K-028, K-031, K-053 and K-091 all say the same thing from different angles: the catalog is a
measurement made by an instrument that is itself inside the system, and the instrument's state is a
covariate we have never entered. Family 3's standing warning (a solar hit is an observer hit until proven
otherwise) is unresolvable without this family.

**F7-01 diurnal_detection_amplitude** | READY
Def: trailing 365 d amplitude of the 24 h cycle in event counts computed from `day_float` (sub-day
timestamps already returned by `load_event_marks`), at a fixed magnitude band just above global Mc.
Src: on disk. **The data path exists and is used at exactly one call site.**
Hunch: **K-091's always-on positive control, made continuous.** The diurnal detection cycle is a known,
large, physically-certain artifact; its amplitude over time is a direct measurement of the network's
noise floor, and it is free.
Pit: it lives in the notched band, so it must never be computed from daily bins. And it is magnitude-band
dependent by construction: at high magnitude it should vanish, and that vanishing is the calibration.
S15: MEASURABLE, and it is the reference against which every other measurability claim here should be
checked (G-M1).
Price: 31 as a covariate; 0 as a control.

**F7-02 diurnal_amplitude_by_band** | READY
Def: F7-01 computed in magnitude bands 4.5-5.0, 5.0-5.5, 5.5-6.0, >6.0.
Src: on disk.
Hunch: the band at which the diurnal amplitude reaches zero **is** the global completeness magnitude,
measured by a method entirely independent of the magnitude-frequency distribution. That is a real
instrument and I do not believe this program has one.
Pit: high bands have few events; the zero is approached from a noisy direction.
S15: MEASURABLE in the low bands, MARGINAL above 5.5.
Price: 4 controls.

**F7-03 weekly_detection_amplitude** | READY
Def: as F7-01 at the 7 d period.
Src: on disk.
Hunch: the human-operations cycle; a second, independent artifact ruler with a period in the band the
miner is sensitive to (sinc ~1), unlike the diurnal one.
Pit: analysts review events on weekdays, so this can affect the catalog at the review stage rather than
the detection stage, which has a different magnitude dependence. That difference is diagnostic.
S15: MEASURABLE.
Price: 31.

**F7-04 station_count_proxy** | DOWNLOAD
Def: number of contributing stations per event, averaged over trailing 90 d.
Src: ComCat `nst` field (may need re-fetch; sparsely populated for older global events).
Hunch: the most direct available measure of instrument state.
Pit: sparsely populated and inconsistently defined across contributing networks; may be unusable globally
and usable regionally (California).
S15: MEASURABLE only where populated. Declare coverage before use.
Price: 31.

**F7-05 azimuthal_gap_median** | DOWNLOAD
Def: trailing 90 d median azimuthal gap.
Src: ComCat `gap` field.
Hunch: a location-quality index that responds to network changes faster than station count does.
Pit: as F7-04.
S15: as F7-04.
Price: 31.

**F7-06 location_uncertainty_median** | DOWNLOAD
Def: trailing 90 d median horizontal error.
Src: ComCat `horizontalError`.
Hunch: as F7-05; and location uncertainty directly limits every spatial statistic in family 4, so it
belongs as a conditioning variable for those.
Pit: reported inconsistently; many nulls.
S15: as F7-04.
Price: 31.

**F7-07 magnitude_type_mix** | READY
Def: trailing 90 d fractions of mb, Mw, ml, mwc etc.
Src: on disk (`magType` column).
Hunch: **magnitude type changes are magnitude scale changes**, and a shift in the mix shifts b, Mc, and
every magnitude-derived feature at once. This is the most likely single cause of any long-period trend in
family 4.
Pit: none as a diagnostic; the risk is only in mistaking it for physics.
S15: MEASURABLE.
Price: 31, and its best use is conditioning.

**F7-08 mb_minus_Mw_offset** | READY
Def: trailing 365 d mean difference for events reporting both.
Src: on disk.
Hunch: a direct measurement of scale drift between magnitude systems; if it moves, b moves for free.
Pit: only a subset of events report both; selection by size.
S15: MARGINAL.
Price: 31.

**F7-09 detection_completeness_by_region** | READY
Def: per-region Mc (max-curvature or b-stability), trailing 365 d, entered through the 2R-df instrument.
Src: on disk.
Hunch: the regional version of F4-54; regional Mc drift is the mechanism by which a global instrument
change looks like a regional physical signal.
Pit: Mc estimators are biased and noisy at regional counts.
S15: MARGINAL per region; the 2R-df sum is MEASURABLE.
Price: 31.

**F7-10 catalog_revision_indicator** | DOWNLOAD
Def: indicator of days whose events were subsequently revised (magnitude or location changed after first
publication).
Src: ComCat versioned queries; **needs a build to reconstruct historical revisions.**
Hunch: revision is the observer changing its mind, and the revision rate is a measure of how uncertain
the catalog was at the time. If a "signal" lives preferentially in revised events, it is a processing
artifact.
Pit: historical revision history is not fully retrievable for older events; the reconstruction is partial
by construction and its incompleteness is itself time-varying.
S15: unknown until built.
Price: 31 if built.

**F7-11 observer_conditioned_null** | READY
Def: not a feature: a surrogate null in which event times are resampled with the empirical
time-of-day and day-of-week detection kernel imposed, rather than uniformly.
Src: on disk.
Hunch: **the correct null for any feature that could be an observer artifact.** A circular-shift null
preserves the detection cycle trivially; a null that resamples must reimpose it or it is anti-conservative.
Pit: needs care not to destroy the very structure being tested; only applies to the mark and unbinned
tests.
S15: not applicable.
Price: 0 (it is a null, not a test), but it changes the p of every sub-daily test.

**F7-12 detection_cycle_cross_test** | READY
Def: the mandated cross-check that any surviving family-3 feature has its skill destroyed by conditioning
on F7-01/F7-03, and survives conditioning on F4-20.
Src: on disk.
Hunch: this is the operational form of "an observer hit until proven otherwise"; without it, the family-3
tranche cannot produce an interpretable result of any kind.
Pit: none; it is the interpretation gate.
S15: not applicable.
Price: 0 new tests, mandatory for family 3.

**F7-13 network_growth_step_detector** | READY
Def: change-point detection (PELT or binary segmentation) on the global daily count series, with the
detected change points entered as indicator features.
Src: on disk.
Hunch: composition drift is not smooth, it is stepwise (a network comes online on a date). Steps are
detectable and, once detected, can be conditioned away.
Pit: change-point detection on a series with real dynamics will find real dynamics and call them steps.
The detected set must be cross-checked against known network deployment dates (documentable) before any
step is treated as an artifact.
S15: MEASURABLE.
Price: 1 detection run plus k indicators, k declared in advance.

**F7-14 magnitude_of_completeness_free_test** | READY
Def: re-run the entire tranche at M>=5.0 and M>=5.5 as a robustness ladder.
Src: on disk.
Hunch: **K-005's M0-invariance audit applied to the miner.** A real physical signal should survive raising
the floor (with reduced power); an incompleteness artifact should not.
Pit: **this multiplies the whole tranche's test count by 3 if priced as new tests.** It should be priced
as a robustness ladder on declared survivors only, not as a parallel battery, and Popper should rule on
that. At M>=5.5 the global count drops to roughly 9,000 (A_min ~12.5% deflated), so the ladder loses
power fast and a non-survival at M>=5.5 is weak evidence of artifact.
S15: M5.0 MEASURABLE; M5.5 MARGINAL.
Price: 0 if a post-hoc ladder on survivors; 2x tranche if a parallel battery. **Flagged for Popper.**

**F7-15 observer_feature_family_null_calibration** | READY
Def: run the entire observer family as if it were a physics family and report how many survive BH; the
expected number is the family's own false-positive calibration.
Src: on disk.
Hunch: **K-049/K-052's self-audit logic**: the observer family has a known answer (these features affect
the catalog), so the number that survive is a measurement of the pipeline's sensitivity rather than of
the Earth.
Pit: it is only a calibration if the true effects are known in size, and they are not; treat as
qualitative.
S15: not applicable.
Price: counted within the family.

---

# FAMILY 8 - CLOCKS AND MEASURES (WHICH x-AXIS MAKES THE SYSTEM SIMPLE?)

Jim's 2026-08-09 directive: any accumulating or state property is a candidate x-axis, and the clock that
simplifies the system IS a discovery about the system. This family does not add features; it **replays
existing features on a different time axis**, which makes it the highest-leverage-per-line-of-code family
in the catalog and also the one whose multiplicity must be priced most carefully.

**The pricing question I am handing Popper explicitly.** Re-running a feature set on a new clock is a new
test per feature. But a clock is not a feature, and re-running 100 features on 6 clocks is 600 tests, not
106. **I propose that clocks be priced as a declared axis with its own stratum in the BH structure**
(§P6-3 stratified BH), so that the multiplicity is honest but the discovery of a good clock is not made
impossible by its own price. That is a proposal, not a ruling.

**F8-01 clock_natural_time_event_count** | READY
Def: replace the day index with cumulative event count; bin into equal-count bins.
Src: on disk.
Hunch: **K-017's live remnant.** Calendar periodicity is a corpse; natural-time periodicity is not the
same hypothesis, and the transformation is exactly the time-rescaling move that makes ETAS residuals
uniform.
Pit: natural time destroys the alignment with any external clock, so **every ephemeris feature becomes
meaningless in it**; natural time is the right clock for catalog-endogenous features and the wrong clock
for astronomical ones. Do not run family 1 in natural time; that is a category error and it would look
like a result.
S15: equal-count binning at 46,585 events gives whatever bin count you declare; the count per bin is
constant by construction, which removes the Poisson heteroscedasticity that daily binning has. That is
a genuine statistical gain.
Price: 1 clock x (family 4 feature count).

**F8-02 clock_expected_events_ETAS** | READY
Def: rescale time by the integrated ETAS intensity `tau = integral Lambda dt` (the residual-analysis
transformation).
Src: on disk.
Hunch: this is **the** canonical time rescaling, and under a correct model the process is homogeneous
Poisson in tau. Any structure surviving in tau is structure ETAS does not explain, which is exactly what
K-009 measured, generalised to every feature in the catalog.
Pit: the rescaling depends on the ETAS fit, so a bad fit manufactures structure; and the fit is
exploration-only.
S15: MEASURABLE.
Price: 1 clock x (feature count).

**F8-03 clock_cumulative_moment** | READY
Def: rescale by cumulative seismic moment released.
Src: on disk.
Hunch: if the system is driven by energy rather than by events, moment is the clock; and moment time is
dominated by large events, so it stretches the quiet periods, which is where quiescence claims live.
Pit: moment time is nearly a step function (one M8 dominates a decade); the transformation is extremely
non-uniform and most bins will hold almost no wall-clock time.
S15: effective bin count is small. MARGINAL.
Price: 1 clock.

**F8-04 clock_cumulative_tidal_stress** | READY
Def: rescale by the cumulative absolute tidal potential (or, once F1-30 exists, cumulative |Coulomb tidal
stress|).
Src: closed form now; better after F1-30.
Hunch: **Jim's directive, applied literally.** If the system integrates tidal work, the correct clock is
tidal-stress time, and a signal invisible in wall time can be exact in stress time.
Pit: the cumulative absolute tidal potential is nearly linear in wall time (it is a near-stationary
oscillation), so this clock is close to the identity and will produce a near-duplicate of the wall-clock
result. **Say that up front**: the informative version is cumulative stress **above a threshold**, which
is not near-linear.
S15: MEASURABLE.
Price: 1 clock, and I recommend the thresholded variant.

**F8-05 clock_cumulative_ledger_loading** | READY
Def: rescale by the stress-ledger's accumulated loading (SoCal only).
Src: on disk.
Hunch: K-079's "measure quiescence in accumulated strain" as a global axis rather than a per-entry fix.
Pit: SoCal-only; the ledger's degeneracy is CONTESTED (M-002).
S15: SoCal counts too small at M>=4.5 (see F5-06). UNMEASURABLE in the global miner.
Price: 1 clock, contingent on a SoCal instrument.

**F8-06 clock_cumulative_geodetic_strain** | INFRASTRUCTURE
Def: rescale by integrated geodetic strain rate from the GNSS field.
Src: `data/ngl/midas.IGS14.txt` and `data/kreemer_young` are on disk for the velocity field, but a
time-varying global strain series needs a build.
Hunch: the physically correct loading clock; everything else in this family is a proxy for it.
Pit: M-006 stands: site dilatation is unreliable to worse than a factor 2 including sign. A strain clock
built from unreliable dilatation is an unreliable clock, and its unreliability is multiplicative through
every test that uses it.
S15: unknown until built.
Price: 1 clock, after a build.

**F8-07 clock_phase_of_dominant_cycle** | READY
Def: rescale by cumulative synodic phase (i.e. lunar-month time).
Src: closed form.
Hunch: exactly linear in wall time, so this is the **identity-clock control**: it must reproduce the
wall-clock result exactly, and if it does not there is a bug in the clock machinery.
Pit: none; it is a control.
S15: not applicable.
Price: 1 control.

**F8-08 clock_aftershock_free_time** | READY
Def: delete all declustered-aftershock days (Zaliapin) and concatenate the remainder.
Src: on disk.
Hunch: **background time.** If the interesting physics is in the background rate and the aftershocks are
noise, this clock removes the noise entirely rather than modelling it, and it is a completely different
kind of control from an ETAS baseline.
Pit: declustering is not a solved problem and the deleted set is a choice; and concatenation creates
artificial discontinuities that will show up at the periods matching the deleted-interval statistics.
The null must be built on the same concatenation.
S15: removes ~60-80% of events depending on the threshold, so N drops to ~10,000-19,000: A_min rises to
~9-12% deflated. **A real power cost that must be quoted whenever this clock is used.**
Price: 1 clock.

**F8-09 clock_distance_to_criticality** | READY
Def: rescale by cumulative `1/(1-n(t))` using the branching-ratio series (F4-39).
Src: on disk.
Hunch: **the boldest clock in the catalog**: if the system's own distance to criticality sets its
internal tempo, then time measured in units of criticality-approach is the clock in which the dynamics
are simple. Nobody has tried this and it costs one line once F4-39 exists.
Pit: n is estimated with large error, so the clock is noisy, and a noisy clock smears everything.
Estimation error in the clock is a fundamentally different problem from estimation error in a feature and
I do not know how to bound it cleanly. Flagged honestly as the least-defensible entry in the family.
S15: MEASURABLE only if F4-39's error is small; probably MARGINAL.
Price: 1 clock.

**F8-10 clock_completeness_corrected_time** | READY
Def: rescale by the cumulative estimated detection capability (from F7-01's diurnal amplitude, inverted).
Src: on disk.
Hunch: **an observer clock**: time measured in units of "how much the network could see". Any secular
trend produced by instrument growth is exactly flattened by this clock, which makes it the cleanest
possible artifact control for long-period claims.
Pit: circular if the completeness estimate depends on the counts. Use the diurnal-amplitude estimator
(F7-01), which is count-normalised and therefore mostly free of this.
S15: MEASURABLE.
Price: 1 clock.

**F8-11 measure_change_magnitude_weighting** | READY
Def: not a time change but a measure change: weight each event by `10^(1.5 M)` (moment) or by
`10^(-b M)` (GR-flattening) in every count statistic.
Src: on disk.
Hunch: **the y-axis is a choice too.** Counting events treats an M4.5 and an M7 identically; every other
weighting is a different question, and only one has been asked.
Pit: moment weighting makes the statistic a heavy-tailed sum whose CLT convergence is slow; the block
bootstrap handles it but the normal approximation in the score statistic does not. Use the empirical
null exclusively for weighted statistics.
S15: moment weighting reduces the effective N drastically (a few events dominate): A_min likely > 30%.
**MARGINAL to UNMEASURABLE for moment weighting; MEASURABLE for GR-flattened weighting**, which is the
better-behaved choice and the one I recommend.
Price: 1 measure x (feature count).

**F8-12 measure_change_per_region_normalisation** | READY
Def: normalise each region's contribution by its own long-run rate before summing.
Src: on disk.
Hunch: the global sum is currently dominated by the most active regions; equal-weighting regions asks a
different question ("do regions respond?") than count-weighting ("do events respond?"), and both are
legitimate and only one is asked.
Pit: equal-weighting amplifies the noisiest regions, which lowers power. Quantify the loss before running.
S15: reduces effective N; MARGINAL.
Price: 1 measure.

**F8-13 clock_time_since_last_local_event** | READY
Def: per cell, time since the last event, aggregated to a global distribution summary.
Src: on disk.
Hunch: a **hazard-clock**: the natural x-axis for a renewal process, and the axis on which
quasi-periodicity (if any exists) would appear.
Pit: the recurrence-interval literature at M>=4.5 in 1 deg cells is dominated by aftershocks, so this is
mostly an Omori measurement.
S15: MEASURABLE as an aggregate.
Price: 1 clock.

**F8-14 clock_reverse_time** | READY
Def: run the entire tranche on the time-reversed catalog.
Src: on disk.
Hunch: **K-014's time-reversal asymmetry as a blanket control.** Any feature that scores the same
forwards and backwards is measuring a symmetric property (correlation, not causation); features that
score only forwards are the causal candidates. This is a one-line change that doubles the interpretive
content of every result in the tranche.
Pit: ETAS is strongly time-asymmetric, so the baseline must be refitted on the reversed catalog or the
comparison is meaningless. That refit is the only real work.
S15: same as forward.
Price: 1 clock x feature count, and it is the best-value clock in the family.

**F8-15 clock_random_reference** | READY
Def: a random monotone time warp with the same coarse statistics as the real clocks.
Src: closed form.
Hunch: the **null clock**. If six clocks are tried, the improvement of the best over a random warp is the
only honest measure of whether the clock discovery is real.
Pit: none; it is the control that makes the whole family interpretable, and without it a clock scan is a
p-hacking machine.
S15: not applicable.
Price: 1 control per clock scan. **Mandatory if any clock scan is run.**

---

# FAMILY 9 - STATISTICS AND INSTRUMENTS (INCLUDING THE SECOND MOMENT AND THE UNBINNED PATH)

§K87-0(c) is the sharpest limitation on the record: **every statistic in the session is a first-moment,
single-phase statistic, and a cross-cycle two-stage process is invisible to all of them.** No number of
surrogates fixes that. This family is the fix, plus the instruments that unlock the notched bands.

These are not features; they are **new observables computed on existing features**, so their price
multiplies the feature set rather than adding to it. That is stated per entry.

**F9-01 second_phase_moment** | READY
Def: for a candidate cycle, compute the **second circular moment** of the event-phase distribution,
`|mean(exp(2 i theta))|`, in addition to the first.
Src: on disk (`day_float` plus closed-form phases).
Hunch: **K-088's core insight.** If a patch unlocks at phase X in cycle N and releases at a dispersed
phase in cycle N+k, the first moment decays geometrically toward zero while a bimodal or axial structure
survives in the second. The current instrument is measuring the wrong moment.
Pit: the second moment is also driven by any two-lobed structure, including a semidiurnal or half-cycle
response, so a hit is ambiguous between "axial response" and "harmonic response" and must be reported as
such. And §P5-3 mandated two repairs to K-088's permutation null under S-13; those repairs bind here.
S15: the second moment's null variance is 1/N as for the first, so A_min is comparable: MEASURABLE at
the same order, ~5.5% deflated. **This is the key point: the missing moment costs nothing in power.**
Price: doubles the phase-feature test count (17 cyclic features x 1 new statistic = 17 in tranche 1).

**F9-02 phase_phase_bicoherence** | READY
Def: the normalised bispectrum between two cycles: coherence at the sum and difference frequencies of
(synodic, draconic), (synodic, anomalistic), (annual, synodic), etc.
Src: on disk.
Hunch: **phase-phase coupling is the natural language for "the response depends on where you are in TWO
cycles at once"**, which is exactly the perigean-spring hypothesis stated properly, and it is a
second-order statistic that no product feature fully captures.
Pit: bicoherence estimation is data-hungry and biased upward at low counts; the bias is the classic
failure mode of this statistic. Needs a matched surrogate null (phase-randomised, amplitude-preserving),
not a block bootstrap.
S15: MARGINAL at 46,585 events over a 21 y record for the low-frequency pairs; better for the
fortnightly pairs.
Price: `C(k,2)` pairs; at k = 6 primary cycles that is 15 tests.

**F9-03 circular_variance_trailing** | READY
Def: the trailing 90 d circular variance of event phases relative to a given cycle, entered as a
covariate (i.e. a feature built from a statistic).
Src: on disk.
Hunch: **a meta-feature: the concentration of the phase distribution has its own dynamics.** If phase
concentration rises before large events, that is a second-order precursor that no phase test can see,
because phase tests average concentration over the whole record.
Pit: strongly count-dependent (circular variance is biased by N); the bias correction is mandatory.
S15: MEASURABLE at 90 d.
Price: 31 per cycle, and I would run it for 3 declared cycles, not 17.

**F9-04 kuiper_and_watson_omnibus** | READY
Def: replace or supplement the 2-df Rayleigh-type score statistic with Kuiper's V and Watson's U^2 on the
event-phase distribution.
Src: on disk.
Hunch: **the 2-df quadratic form is optimal against a sinusoid and weak against everything else.** A
sharp response confined to 10% of the cycle is nearly invisible to it and obvious to Kuiper. The program
has bounded sinusoidal modulation at 5.6% and has said nothing at all about non-sinusoidal modulation.
**This may be the single largest unexamined gap in the tidal null.**
Pit: Kuiper and Watson have no closed-form null under a non-uniform baseline intensity, so the null must
be simulated from the ETAS Lambda; that is a real build but a small one. And they are less powerful than
Rayleigh when the truth IS sinusoidal, so both must be reported.
S15: MEASURABLE; for a narrow response the effective A_min in "fraction of events displaced" terms is
considerably better than the sinusoidal bound.
Price: 2 statistics x 17 cyclic features = 34 in tranche 1.

**F9-05 harmonic_ladder_on_all_cyclic** | READY
Def: extend the existing {P/3, P/2, P, 2P, 3P} ladder to every cyclic feature and every partition.
Src: on disk.
Hunch: a sharp response has its power in the harmonics; the ladder is the poor man's Kuiper and it is
already implemented.
Pit: rungs are never free (§P5-5). The declared tranche-3 count in the code comments is 68; that number
should be honoured rather than re-derived.
S15: MEASURABLE.
Price: 68 (as already declared in the engine).

**F9-06 epoch_folding_on_transient_anchors** | READY
Def: stack the daily residual around declared transient anchors (storm onsets, snowmelt dates, M7
occurrences, eclipse dates) and test the stacked profile.
Src: on disk.
Hunch: **transients are the class the bounds paper does not cover (S-14(c) bracketed).** A periodic test
cannot see an event-anchored response; an epoch fold can, and it is the standard instrument in every
other field that has this problem.
Pit: anchor counts are small (10-100), so the stack is noisy; and anchors correlate with each other and
with season. The null must shift anchors, not events, and must preserve anchor spacing.
S15: with 50 anchors and a 30 d window, the stack holds ~1,500 event-days of signal: **A_min ~30%
deflated. MARGINAL to UNMEASURABLE for anything but a large transient response.** State this before
proposing any transient test.
Price: 1 per declared anchor set; ~6 anchor sets in this catalog.

**F9-07 wavelet_coherence_feature_target** | READY
Def: time-resolved coherence between a feature and the residual series in the time-frequency plane.
Src: on disk.
Hunch: **the response may be intermittent.** A stationary test of an intermittent coupling is diluted by
exactly the duty cycle; wavelet coherence sees intermittency directly.
Pit: wavelet coherence is notorious for producing visually compelling patches under pure noise; the null
must be the same surrogate machinery and the statistic must be a scalar summary (e.g. total
above-threshold area), not a picture. **A picture is not a test and must not be shipped as one.**
S15: MARGINAL; the multiple-comparison structure across the time-frequency plane is severe.
Price: 1 scalar summary per feature; 17 in tranche 1.

**F9-08 phase_locking_value** | READY
Def: Hilbert-transform phase-locking value between a feature's analytic phase and the residual's analytic
phase, in a declared band.
Src: on disk.
Hunch: PLV is the neuroscience-standard coupling statistic and is amplitude-independent, so it is immune
to the amplitude drift (network growth) that contaminates correlation.
Pit: the Hilbert transform of a count series is ill-defined without pre-filtering, and the filter choice
is a hidden parameter. Declare one band.
S15: MARGINAL.
Price: 1 per feature per declared band.

**F9-09 conditional_intensity_phase_residual** | READY
Def: fit the phase response **jointly** with the ETAS parameters rather than as a covariate on fixed
residuals.
Src: on disk.
Hunch: fixing the baseline and then testing a covariate underestimates the covariate's effect when the
baseline has absorbed part of it; joint fitting is the correct and more powerful procedure.
Pit: joint fitting can also **manufacture** an effect by trading against mu; the parameter correlation
must be reported. And it is much slower.
S15: MEASURABLE, with a modest power gain that should be quantified by injection before it is claimed.
Price: replaces rather than adds; same count, more compute.

**F9-10 magnitude_marked_phase_tests** | READY
Def: extend the existing mark tests (currently 2 marks) to: magnitude, depth, log-moment, distance to
nearest prior event, time since prior event, cluster membership, region.
Src: on disk.
Hunch: **the mark axis is nearly unexplored** (23 features x 2 marks = 46 tests in the last session).
A forcing that does not change the rate may still change WHICH events happen, and mark tests see that
while count tests cannot.
Pit: marks are correlated with each other and with the rate; and the circular-linear correlation
statistic is weak against non-monotone dependence.
S15: mark tests use the full event set (46,585): MEASURABLE, and they escape the day-binning sinc
entirely because they are computed at event times (DATED CORRECTION 2026-08-13, K-092 flag: this holds ONLY if the feature itself is re-derived at event time day_float, as marks_ext.py states in code; a daily-sampled feature looked up per event does NOT escape the sinc). **That last property is underexploited and is worth
a tranche on its own.**
Price: 7 marks x 23+ features = 161+ tests.

**F9-11 peaks_over_threshold_on_daily_counts** | READY
Def: model the upper tail of daily counts with a GPD and test whether feature values differ on exceedance
days.
Src: on disk.
Hunch: **K-042's joint-extremes logic**: if the response is a threshold effect, it lives in the top
percentile of days and the mean-based statistics dilute it by a factor of 100.
Pit: exceedance days are aftershock days, so this is largely an aftershock detector under a climatology
baseline; ETAS baseline mandatory. And §P6-2's severe conditions on GPD extrapolation bind here.
S15: top 1% of 7716 days = 77 days: **MARGINAL.** Top 5% = 386 days: MEASURABLE for a large effect.
Price: 31 per declared threshold; declare one.

**F9-12 quantile_regression_on_counts** | READY
Def: regress the 0.1, 0.5, 0.9 quantiles of the daily count on each feature.
Src: on disk.
Hunch: a forcing might sharpen the distribution without moving its mean; the mean-only GLM cannot see
that and quantile regression can.
Pit: three quantiles is three tests per feature; and count quantiles are integers, which makes quantile
regression awkward at low rates. Global daily counts are ~6/day, so this is workable but lumpy.
S15: MEASURABLE at the median, MARGINAL at the tails.
Price: 3 x feature count.

**F9-13 recurrence_quantification** | READY
Def: recurrence-plot determinism and laminarity measures on the residual series, trailing 365 d.
Src: on disk.
Hunch: RQA detects deterministic structure that autocorrelation misses, and laminarity specifically
detects **intermittency and chaos-to-chaos transitions**, which is the shape a regime change would take.
Pit: RQA has three free parameters (embedding, delay, threshold) and its literature is full of results
that do not survive parameter changes. Declare all three per S-9 and run a sensitivity as an audit, not
as a scan.
S15: MARGINAL.
Price: 31 as a covariate.

**F9-14 permutation_entropy** | READY
Def: ordinal permutation entropy of the daily residual series, trailing 365 d, order 4.
Src: on disk.
Hunch: a robust, parameter-light complexity measure; **complexity dropping before a transition is a
generic critical-slowing-down signature** and it is cheap.
Pit: order and window are parameters; and permutation entropy of a near-Poisson series with many ties is
degenerate. Tie-breaking rule must be declared.
S15: MEASURABLE at 365 d.
Price: 31.

**F9-15 unbinned_event_time_scoring** | INFRASTRUCTURE
Def: score features at true event times using `day_float`, with the baseline intensity evaluated
continuously, instead of on daily bins. **K-090(a), already ADMIT-PROMOTED in the ledger.**
Src: on disk (`load_event_marks` already returns `day_float`; it is consumed at exactly one call site).
The build is the continuous baseline evaluation, not the data.
Hunch: **this single build removes the exact zeros.** Day binning multiplies a period-P sinusoid by
|sinc(pi Delta t/P)|: M2 by 0.0348, K1 by 0.0027, S1 and S2 by exactly zero. The entire diurnal and
semidiurnal band is currently unobservable **for an instrumental reason**, and one build fixes it.
Pit: it gives up the miner's current **immunity** to the S1/S2 detection-cycle systematic, which
§K87-0(a) correctly identifies as a real asset. The build must therefore ship with the observer controls
(F7-01, F7-03, F7-11) as mandatory co-runs, not optional ones. Also: sub-daily tests need site-local
phase (F1-30) to be meaningful, because sub-daily tidal phase is a function of longitude and the domain
sum cancels it (the second exact zero). **One build is not enough; it is two, and they are F9-15 and
F1-30.**
S15: post-build, MEASURABLE at the full N in the sub-daily band, per region.
Price: unlocks ~124-310 tests in family F1 that are currently priced at zero because they are
unmeasurable.
**Together with F1-30 this is the highest-leverage build in the catalog.**

**F9-16 sub_daily_site_local_pipeline** | INFRASTRUCTURE
Def: the composition of F9-15 and F1-30: per-event, site-local, phase-resolved tidal Coulomb stress at
the event's own hypocentre and mechanism.
Src: builds named above.
Hunch: this is what the tidal literature actually tests and what this program has never been able to
compute; every tidal result the program owns is a geocentric-proxy result.
Pit: it is the most confound-dense measurement in the catalog: mechanism uncertainty, ocean loading,
depth error, and the detection cycle all land on it at once. It should ship with a planted-signal
recovery per band **and per mechanism class** before any number is quoted (S-17 candidate, §P6-4 item 5).
S15: MEASURABLE post-build.
Price: 31 x (mechanism strata).

**F9-17 selection_debiased_effect_estimator** | READY
Def: report, for every surviving test, a selection-debiased effect size (conditional-likelihood or
bootstrap-shrinkage), alongside the median effect among survivors of an equal-sized null run.
Src: on disk.
Hunch: **§P6-4 banner item 3, which nobody asked for and which is the biggest new hazard v2 introduces.**
At 10^6 tests every quoted amplitude is the maximum of a huge search and is upward-biased by selection;
without this, v2's amplitudes will exceed v1's for purely combinatorial reasons and someone will read
that as a finding.
Pit: none; it is a mandatory correction. The only risk is not doing it.
S15: not applicable.
Price: 0 tests, mandatory.

**F9-18 GPD_p_extrapolation_with_coverage_check** | READY
Def: extrapolated p-values from the surrogate tail, with the §P6-2(7) coverage demonstration.
Src: on disk.
Hunch: at 50,000 surrogates the p floor is 2e-5 and BH over 10^6 tests cannot reject at all; tail
extrapolation is the only way a large tranche can produce a rejection.
Pit: §P6-2's conditions are severe and they bind: the CI must be shown to cover a brute-force p and never
understate it by more than 3x, at the argument value used (S-17 candidate).
S15: not applicable.
Price: 0 tests, enabling.

**F9-19 planted_signal_recovery_per_band_and_aggregation** | READY
Def: G-M1 extended: plant a known sinusoid at each declared band, at each aggregation level (global,
2R-df, per-region), measure recovered/planted amplitude and phase error.
Src: on disk.
Hunch: §P6-4 item 5 makes this mandatory before any bound is quoted at any band or aggregation; it is
also the only thing that will tell us whether my §0.3 deflation factor is right (see F4-58).
Pit: none; it is the gate.
S15: it IS the S-15 machinery.
Price: 0 priced tests. **Run before the tranche, not after.**

**F9-20 negative_control_feature_battery** | READY
Def: a declared set of features with certain-zero effect (F1-35 Jupiter tide, F1-36 Jupiter-Saturn
synodic, F2-09 Metonic, plus 20 random-phase synthetic cycles at matched periods) run inside the same
tranche.
Src: closed form.
Hunch: **the tranche's own false-positive rate, measured rather than assumed.** BH controls FDR under
assumptions; a matched negative-control battery measures it under the actual dependence structure of the
actual data, which is worth more than the theorem.
Pit: synthetic cycles must be matched in period, autocorrelation and amplitude distribution to the real
features or the calibration is not a calibration. And the controls must be counted in the multiplicity,
not exempted.
S15: not applicable.
Price: ~23 features x 31 lags = 713 tests, which is a real price and worth it. **I would rather spend
700 tests on knowing my false-positive rate than on 700 more speculative features.**

---

# FAMILY 10 - COMBINATIONS: PRODUCTS, CONDITIONALS, AND LAG STRUCTURES

**The multiplicity price is the whole subject of this family**, so it is stated per class before the
entries, as Jim asked.

### The combination-class price table (PROPOSED to Popper)

| Class | Construction | Price formula | At k = 40 features, L = 31 lags |
|---|---|---|---|
| C1 pairwise product | `z_i * z_j`, both centred | `C(k,2) x L` | 24,180 |
| C2 pairwise product, declared subset | only pairs with a named mechanism | `p x L`, p declared | ~10-30 pairs: 310-930 |
| C3 phase-phase coupling | bicoherence over cycle pairs | `C(k_cyc,2)` | 15 at k_cyc = 6 |
| C4 two-stage conditional | feature tested only within a declared state stratum | `k x L x S` strata | S = 2: 2,480 |
| C5 lag-difference | `z(t-L1) - z(t-L2)` | `C(L,2) x k` | 18,600 |
| C6 interaction with region | feature x region under 2R-df | `k x L` (unchanged) | 1,240 |
| C7 triple product | `z_i z_j z_k` | `C(k,3) x L` | 304,720 |
| C8 clock x feature | family 8 replay | `k x L x n_clocks` | 6 clocks: 7,440 |
| C9 statistic x feature | family 9 replay | `k x L x n_stats` | 5 stats: 6,200 |

**My recommendation, stated as a recommendation and not a ruling.** C7 is where a brute-force machine
goes to die: 300,000 triple products at q = 0.1 yields 30,000 expected chance survivors, and the
winner's-curse bias on the top of that pile is severe (F9-17). **C1 is nearly as bad and its products
are mostly collinear with their parents.** The good value is in C2, C3, C4, C6, C8 and C9, all of which
are structured rather than combinatorial. I would spend the machine's new capacity on **structure**
(clocks, statistics, partitions, conditionals) and not on **combinatorics** (all-pairs, all-triples),
and if the machine must run a combinatorial class I would run C1 once, as a declared
EXPLORATORY-UNPRICED sweep whose only published output is the negative-control-calibrated survivor count.

**F10-01 tidal_x_ledger_class** | READY
Def: tidal feature tested only within the high-loading stratum (C4).
Hunch: **K-039R/K-071/K-077's shared shape**: small-stress response should exist only in near-critical
material, so the unconditional test is diluted by the fraction of the crust that is not near-critical.
Pit: the stratum must be defined from exploration data; SoCal-only ledger limits it (see F5-06). And the
sqrt(2) power loss of splitting must be quoted.
S15: MARGINAL at global M>=4.5.
Price: C4 with S = 2.

**F10-02 tidal_x_recent_large_event** | READY
Def: tidal features conditioned on "within 30 d and 1000 km of a global M>=7" (C4).
Hunch: K-043/K-073's wetness-meter logic: a stress-transient-primed crust should be the most tidally
responsive state available, and it is a state we can identify without any geodesy.
Pit: it is also the state with the most aftershocks, i.e. the most model misspecification; ETAS baseline
plus a residual-based target both mandatory.
S15: the conditioned subset is a small fraction of event-days: MARGINAL.
Price: C4, S = 2.

**F10-03 tidal_x_swarm_state** | READY
Def: tidal features conditioned on the swarm index (F4-65) being in its top tercile (C4).
Hunch: swarms are fluid-driven and fluid-driven systems are the ones where documented tidal sensitivity
is strongest (LFEs, tremor); this is the best-motivated conditional in the catalog.
Pit: swarm identification parameters must be fixed a priori.
S15: MARGINAL.
Price: C4, S = 3 terciles (or 2 if declared as top-versus-rest, which I prefer).

**F10-04 tidal_x_depth_class** | READY
Def: tidal features x shallow/intermediate/deep (C6 with the depth partition).
Hunch: tidal stress amplitude and the fluid content both vary strongly with depth; a response averaged
over depth classes with opposite signs cancels.
Pit: K-078 already covers depth-class sign predictions and is bound-producing unless n >= 40.
S15: shallow MEASURABLE, deep thin.
Price: C6.

**F10-05 geomagnetic_x_completeness** | READY
Def: family-3 features conditioned on the observer state (C4 with F7-01's diurnal amplitude tercile).
Hunch: **the discriminator between physics and observer for the whole of family 3**: a real effect should
be independent of the network's state; an artifact should track it.
Pit: the conditioning variable and the confound are the same thing, so the design must be a formal
mediation analysis, not a stratification, if it is to be interpretable. Flagged as a design question.
S15: MARGINAL after splitting.
Price: C4, S = 2.

**F10-06 seasonal_x_hemisphere** | READY
Def: annual features x hemisphere sign (C6).
Hunch: the exact-cancellation fix of F6-01, expressed as a combination class.
Pit: as F6-01; it is a bug fix, not a discovery.
S15: MEASURABLE under 2R-df.
Price: C6.

**F10-07 storm_x_quiet_run** | READY
Def: storm-onset response conditioned on the preceding quiet-run length (C4).
Hunch: **habituation**: the response to a stimulus may depend on how long since the last stimulus, which
is a shape no additive model can express.
Pit: small N on both axes.
S15: UNMEASURABLE except as a bound.
Price: C4, S = 2.

**F10-08 perigee_x_syzygy_explicit_product** | READY
Def: the explicit product of the anomalistic and synodic phase designs, i.e. the full 4-column
interaction rather than the existing single `perigee_syzygy` proxy column (C3-adjacent).
Hunch: the existing proxy is one particular scalar combination; the full interaction spans the space of
all bilinear phase-phase couplings and the proxy is a 1-dimensional slice of a 4-dimensional object.
**We have tested one direction in a four-dimensional space and reported the space as null.**
Pit: 4 df instead of 1, so the power per df is lower; and the interpretation of a hit is harder.
S15: MEASURABLE.
Price: 1 feature x 31 lags at 4 df.

**F10-09 all_cycle_pairs_bilinear** | READY
Def: the F10-08 construction for all `C(9,2) = 36` family-1 phase pairs.
Hunch: the systematic version; this is the honest form of "we scanned the beat space", which the current
five hand-picked beats do not constitute.
Pit: 36 x 4 df x 31 lags = 4,464 df-heavy tests; and adjacent pairs are strongly collinear. This is a C1
-class combinatorial expense and I flag it as the borderline case: it is structured (it spans a specific
space) but it is also large.
S15: MEASURABLE per test.
Price: 4,464.

**F10-10 lag_difference_features** | READY
Def: `z(t-L1) - z(t-L2)` for declared lag pairs (C5).
Hunch: a difference of lags is a **crude band-pass**: it isolates the component of the feature that
changed over the interval, which is the rate hypothesis expressed nonparametrically.
Pit: `C(31,2) = 465` pairs per feature is a combinatorial blow-up for very little new information, since
the analytic derivative (F1-11) already captures the limit. **Declare 3 pairs, not 465.**
S15: MEASURABLE.
Price: C5 restricted: 3 pairs x k features.

**F10-11 cumulative_dose_features** | READY
Def: trailing sums of any family-1/2/3 feature over 7, 30, 90 d.
Hunch: **integrative rather than instantaneous response.** Everything in family 1 is currently entered as
an instantaneous value; if the crust integrates, the integral is the covariate and the instantaneous
value is a poor proxy for it.
Pit: trailing sums of near-periodic features are the same feature attenuated and phase-shifted (a boxcar
filter is a sinc in frequency), so most of these are **exactly the lag scan in disguise** and would be
paid for twice. Compute the boxcar transfer function per feature-period and drop the redundant ones
before pricing. This is a genuine saving of hundreds of tests.
S15: MEASURABLE.
Price: 3 windows x k, minus the analytically-redundant majority.

**F10-12 ratio_to_long_run_mean** | READY
Def: any feature divided by its own trailing 365 d mean.
Hunch: a normalisation that removes secular drift; for family-4 features it is the composition-drift fix.
Pit: **ratios killed the TSI feature.** Use log-differences instead and say so. Where the denominator can
approach zero, this class is unusable.
S15: MEASURABLE.
Price: 1 x k, as a transform not a new hypothesis (Popper should rule whether a transform is a new test;
I think it is).

**F10-13 rank_transformed_versions** | READY
Def: every linear feature replaced by its within-window rank.
Hunch: robustness to heavy tails and to the standardisation problems flagged repeatedly above; a rank
version that agrees with the raw version is evidence the result is not tail-driven.
Pit: doubles the count if priced as new tests; it is better run as a **post-hoc robustness check on
survivors only**, which costs nothing.
S15: MEASURABLE.
Price: 0 as a survivor-only check; 1 x k as a battery.

**F10-14 two_stage_unlock_release** | READY
Def: the explicit K-088 model: an unlock hazard driven by phase in cycle N, a release with a dispersed
delay, fitted jointly.
Hunch: **the only entry in the catalog that directly models the process §K87-0(c) says is invisible.**
Everything else in family 9 detects its symptoms; this one fits it.
Pit: the model has a dispersion parameter that is weakly identified, and §P5-3's two mandated repairs to
K-088's permutation null bind. The ledger records arm (3) as DEFER; this entry is arms (1) and (4).
S15: MARGINAL; identification, not count, is the binding constraint.
Price: 1 model x 17 cyclic features.

**F10-15 conditional_on_own_history** | READY
Def: feature tested only on days where the feature's own trailing 30 d variance is in the top tercile.
Hunch: **an intermittency conditional that requires no external variable**: if the coupling only operates
when the forcing is strong, this finds it using nothing but the forcing itself.
Pit: conditioning on a function of the covariate (not of the target) is legitimate and does not leak, but
it does change the null; the surrogates must be conditioned identically.
S15: splits N by 3: A_min ~9.5% deflated. MARGINAL.
Price: C4, S = 3 (or 2 if top-versus-rest).

**F10-16 feature_x_clock_matrix** | READY
Def: the full family-8 replay: every declared feature on every declared clock (C8).
Hunch: **this is what I would actually spend the 50-100x throughput on.** The clock is the discovery, and
a clock scan is only possible with a machine that can afford it. The wall-clock result is one cell of a
matrix nobody has ever computed.
Pit: 6 clocks x k features x L lags is large but structured; and the random-clock control (F8-15) is
mandatory, without which this is a p-hacking engine. Ephemeris features must be excluded from
catalog-derived clocks (F8-01's category error).
S15: as each clock's own entry.
Price: C8 = 7,440 at k = 40, L = 31, 6 clocks. **With the control and the stratified BH, this is my
top tranche recommendation.**

**F10-17 feature_x_statistic_matrix** | READY
Def: every feature scored by every family-9 statistic (C9).
Hunch: the second-moment gap (§K87-0(c)) and the non-sinusoidal gap (F9-04) are both statistic gaps, not
feature gaps; **the machine's new capacity is better spent on new statistics over old features than on
new features under the old statistic.**
Pit: statistics are correlated with each other (Rayleigh and Kuiper on the same data), so the effective
multiplicity is lower than the nominal, and BH will be conservative. That conservatism is acceptable; the
S-8 max-statistic handles the dependence exactly.
S15: as each statistic's own entry.
Price: C9 = 6,200 at 5 statistics.

**F10-18 partition_x_feature_matrix** | READY
Def: every feature under every declared partition through the 2R-df sum (C6).
Hunch: the partition is a hypothesis about **where the phase is coherent**, and testing 4 partitions is
testing 4 hypotheses about the system's spatial organisation, at a cost of 4x rather than 1000x.
Pit: partitions are nested and correlated (latitude and tectonic province overlap); report the overlap.
S15: summed statistic MEASURABLE per partition.
Price: 4 partitions x k x L = 4,960 at k = 40.

**F10-19 triple_conditional_state_vector** | READY
Def: condition simultaneously on (loading state, fluid state, recent-transient state), 8 cells.
Hunch: K-012's state-vector logic; the response may require a conjunction of conditions.
Pit: 8 cells splits N by 8, giving A_min ~16% deflated per cell. **This is the point at which
conditioning becomes self-defeating**, and I include it to mark exactly where that boundary is.
S15: UNMEASURABLE per cell at M>=4.5.
Price: C4 with S = 8; recommendation 0 at global M>=4.5.

**F10-20 negative_space_conjunction** | READY
Def: test the days that are extreme in a feature and yet had NO events, against the model's expectation
(the didn't-fire ledger as a daily statistic).
Hunch: **K-041's negative space at the day level.** The zero-count days carry information that
count-weighted statistics discard, and a forcing that suppresses would show here first.
Pit: zero-count days are rare globally (the global rate is ~6/day, so a zero day is a 0.2% event); this
is a rare-event statistic at the global level and only becomes usable regionally.
S15: **UNMEASURABLE globally** (too few zero days); MARGINAL regionally.
Price: 31 regionally.

**F10-21 cross_family_residual_correlation** | READY
Def: correlate the residuals of the best family-1 model with the residuals of the best family-4 model.
Hunch: **the meta-pattern move: interesting structure lives in the correlation between the residuals of
different models.** If two unrelated model families fail in the same places, something is there that
neither contains.
Pit: both residual sets are dominated by the same aftershock misfit, so a large correlation is expected
and uninformative; the test must be against a null that shares the aftershock structure.
S15: MEASURABLE.
Price: 1 statistic.

**F10-22 model_disagreement_as_feature** | READY
Def: the daily difference between the climatology and ETAS baselines' predicted intensities, entered as a
covariate.
Hunch: **where the two baselines disagree is where clustering matters most**; a feature that only works
where they agree is a background-rate feature, and one that only works where they disagree is a
clustering feature. This decomposes every result for free.
Pit: it is a deterministic function of the catalog, so it is a clustering proxy; interpret accordingly.
S15: MEASURABLE.
Price: 31.

**F10-23 stacked_ensemble_gate** | READY
Def: a gradient-boosted tree over all declared features, with the mandated comparison against a
2-parameter GLM after conditioning.
Hunch: included because the machine can afford it and because refusing to run it is not the same as
knowing what it would say.
Pit: **Mignan and Broccardo, quoted verbatim in the README as the standing warning**: a boosted-tree find
that a 2-parameter GLM cannot reproduce after conditioning is treated as leakage until proven otherwise.
The tree's honest output is a **variable-importance ranking used as a generator**, never as evidence.
S15: not applicable in the S-15 sense; declare it EXPLORATORY-UNPRICED per Rule 4.4's spirit.
Price: 0 priced tests, generator only.

**F10-24 feature_selection_stability** | READY
Def: bootstrap the exploration window and record how often each feature appears in the top-k; report
stability rather than rank.
Hunch: **K-052's winner's-curse ledger made operational**: an unstable top-10 is a top-10 of noise, and
stability is measurable before any holdout is spent.
Pit: bootstrap of a time series must be a block bootstrap at the feature's own timescale (the `mine.py`
F107 lesson).
S15: not applicable.
Price: 0 priced tests, mandatory reporting.

**F10-25 tranche_level_negative_control_ratio** | READY
Def: the ratio of survivors among real features to survivors among the F9-20 negative controls.
Hunch: **the single number I would put at the top of the v2 report.** It is the tranche's own
signal-to-noise, measured under the actual dependence structure, and it is interpretable by a reader who
does not trust anything else in the report.
Pit: requires the negative-control battery to be run in the same tranche, at the same lags, under the
same nulls. If it is run separately it means nothing.
S15: not applicable.
Price: 0 beyond F9-20's cost.

---

# 11. FOUR THINGS NOBODY ASKED FOR

Kepler's charter requires that each round produce at least one item answering a question nobody posed.
These are PROPOSED, brief, and testable.

**N-1. The catalog's own price is a measurable quantity, and we should measure it before we spend it.**
This file declares roughly 281 candidate entries. The distribution of survivors among the negative
controls (F9-20) versus among the real features (F10-25) is a **direct measurement of how much of the
program's hypothesis space is noise**. K-056 proposed randomised sampling of the hypothesis space as an
object of study; this catalog makes that concrete, because a catalog is a sampling frame. Test: run a
random 10% subsample of the catalog and a matched control battery; the survivor ratio estimates the
whole catalog's yield before the whole catalog is run. Pass: the ratio exceeds 1 with a bootstrap CI
excluding 1. Fail: it does not, and the correct action is to stop generating features and start building
instruments (families 7, 8, 9), which is a decision this program can now make quantitatively.

**N-2. The right unit of the v2 report is the bound, not the survivor.** With 10^6 tests and BH, the
expected output is zero survivors and a very large collection of **upper bounds on effects that do not
exist**. That collection is the actual product. Nobody has proposed shipping it. Test: for every declared
test, emit the S-15 floor and the fitted amplitude with its CI, so the tranche produces a **bounds
atlas** rather than a ranked list. Pass: the atlas is quotable per band and aggregation (G-M1 cleared).
Fail: recovery is not demonstrated at a band, in which case that row reads UNRESOLVED and not a number,
per the S-17 candidate.

**N-3. The measurability floor is a feature of the network, not of nature, and it is improving.** A_min
scales as N^-1/2 and N grows with the network. That means **every null in this program has a known
expiry date**, computable now: the year in which the accumulated catalog would resolve the effect that
today's bound excludes. Test: project N(t) from the observed growth in `data/comcat_world` and publish,
for each bound, the year in which it becomes a detection at a given assumed effect size. Pass: the
projection is monotone and the numbers are decades not centuries for the interesting effects. Fail: the
numbers are centuries, which is itself the most decision-relevant result this program could produce.

**N-4. The exploration window is a resource with a budget, and nobody is accounting for it.** Every
explore run increments `EXPLORE_COUNT.jsonl` and every increment raises the multiplicity that the eventual
holdout test must carry. A 50-100x throughput increase therefore **spends the holdout's power**, not just
compute. Test: compute the current `n_explore_runs` and the implied multiplicity penalty on the next
holdout test, and declare a budget before the tranche runs. Pass: the budget is stated and the tranche
fits inside it. Fail: the tranche's mining alone exhausts the holdout's power, in which case the correct
design is a **pre-declared holdout hierarchy** (several independent holdout windows, one per tranche
generation), which does not currently exist and which I am proposing here.

---

# 12. PROPOSED TOP-3 TRANCHE GROUPINGS BY EXPECTED INFORMATION PER PRICED TEST

**These are PROPOSALS. I generate; Popper prices.** The ordering below is my expected-information
ranking, not a licence, and every one of these is subject to G-M1 clearing at the relevant band and
aggregation before any output may be entered for or against a ledger entry.

## TRANCHE A (first, and it is nearly free): THE INSTRUMENT TRANCHE

**Contents:** F4-58 (measure the clustering-deflation factor), F9-19 (planted-signal recovery per band
AND per aggregation level), F9-20 (negative-control battery), F9-17 (selection-debiased effect
estimator), F7-01/F7-02/F7-03 (diurnal and weekly detection amplitude as continuous observer controls),
F8-15 (random-clock control), F10-24 (feature-selection stability), F10-25 (survivor ratio).

**Priced tests: ~713** (essentially all of it the negative-control battery; the rest are controls and
audits that cost compute, not multiplicity).

**Why first.** Every measurability line in this file, and therefore every S-15 declaration in the whole
v2 programme, currently rests on my inference of a deflation factor from two published numbers. That is
one measurement away from being known. And a 50-100x throughput increase without a measured
false-positive rate under the actual dependence structure is not a scaled-up instrument, it is a
scaled-up generator of ranked noise. **Expected information per priced test here is the highest in the
catalog because the price is nearly zero and the output conditions everything else.**

**What a pass and a fail mean.** Pass: recovery ratios in [0.8, 1.2] with phase error under 15 deg at
each declared band and aggregation; the survivor ratio (real versus control) is computable; the deflation
factor is measured, and §0.3's table is rewritten with it. Fail: recovery is not demonstrated at some
band or aggregation, in which case **no bound may be quoted there at all**, and the honest output of the
whole v2 programme shrinks to the bands where it was demonstrated. Either outcome is worth more than any
speculative feature in this file.

## TRANCHE B: THE STATISTIC TRANCHE (new statistics on old features)

**Contents:** F9-01 (second circular moment), F9-04 (Kuiper and Watson omnibus), F9-10 (the mark axis
extended to 7 marks), F2-18..F2-25 (the harmonic ladder on the 8 linear cyclic features), F9-05 (the
declared 68-test ladder), F10-14 (the two-stage unlock-release model), F10-08 (the full bilinear
perigee-syzygy interaction).

**Priced tests: ~1,000** (17 second-moment plus 34 omnibus plus ~161 mark plus 68 ladder plus 32 linear
ladder plus 17 two-stage plus 31 bilinear, plus overhead).

**Why second.** §K87-0(c) says in the program's own words that the two-stage process **could not have
been tested** by any statistic in the session, and §K87-0's verdict line says the second phase moment was
never computed. Meanwhile the entire tidal bound is a bound on **sinusoidal** modulation, and F9-04 says
plainly that a response confined to 10% of a cycle is nearly invisible to the statistic that produced it.
**These are gaps in what was measured, not gaps in what was proposed**, which makes their expected
information far higher per test than any new covariate: the events are already there, the features are
already there, and only the observable is missing. And the mark tests are computed at event times, which
means they **escape the day-binning sinc entirely** (DATED CORRECTION 2026-08-13, K-092 flag: this holds ONLY if the feature itself is re-derived at event time day_float, as marks_ext.py states in code; a daily-sampled feature looked up per event does NOT escape the sinc) - the most underexploited property of the existing
instrument.

**What a pass and a fail mean.** Pass: any survivor is a genuinely new class of finding, because no
first-moment statistic could have produced it. Fail: the tidal null strengthens substantially, because it
stops being a null about sinusoids and starts being a null about phase structure of any shape and about
two-stage processes. **This tranche cannot waste its budget: a null here is a much better null than the
one we have.**

## TRANCHE C: THE CLOCK AND PARTITION TRANCHE (structure, not combinatorics)

**Contents:** F10-16 (feature x clock matrix) over F8-01 natural time, F8-02 ETAS-rescaled time, F8-08
aftershock-free time, F8-10 completeness-corrected time, F8-14 reverse time, with F8-15 the mandatory
random-clock control; plus F10-18 (feature x partition under the 2R-df sum) over F5-01 longitude, F5-02
tectonic province, F5-03 latitude, F5-04 depth class; plus F5-08 (the global coupling summary statistic)
and F4-50 (quiescence area fraction, the honest survivor of the dead per-cell battery).

**Priced tests: ~10,000-12,000** (C8 ~7,440 plus C6 ~4,960, reduced by excluding ephemeris features from
catalog-derived clocks per F8-01's category-error note).

**Why third and why not larger.** This is the tranche that a 50-100x machine makes possible and that a v1
machine could not attempt, and it is the one that answers Jim's clock directive directly: **the clock that
simplifies the system IS a discovery about the system.** It is also the tranche where the 2R-df instrument
gets used for what it is worth - four partitions is four hypotheses about the crust's spatial
organisation at 4x cost rather than 1000x. I am deliberately **not** proposing the combinatorial classes
(C1 all-pairs at 24,180, C7 all-triples at 304,720) even though the machine could run them: their
products are mostly collinear with their parents, their winner's-curse bias is severe, and their expected
information per priced test is the lowest in the catalog. If Jim wants the machine saturated, saturate it
with clocks, partitions, statistics and conditionals, and run C1 exactly once as declared
EXPLORATORY-UNPRICED whose only published output is the control-calibrated survivor count.

**What a pass and a fail mean.** Pass: one clock or one partition shows a materially larger statistic than
the random-clock control, which localises where the system is simple. Fail: no clock beats the random
warp, which is a real and publishable statement - **that the system's structure is not hiding in the
choice of x-axis** - and it retires an entire class of "we were measuring in the wrong units" objections
that currently has no evidence either way.

**Runner-up, named so it is not lost:** the two INFRASTRUCTURE builds F9-15 (unbinned event-time scoring,
already ADMIT-PROMOTED as K-090(a)) and F1-30 (site-local tidal stress tensor). Together they unlock the
diurnal and semidiurnal band, which is currently unobservable **for an instrumental reason and not a
scientific one** - the sinc factors are 0.0348 for M2, 0.0027 for K1, and exactly zero for S1 and S2.
That is the only place in this catalog where the tidal question is genuinely open, and it costs a build
rather than a tranche. I rank it fourth only because it must ship with the observer controls from Tranche
A, which do not exist yet.

---

# 13. AUDIT LINE

**Entries by family (281 total):**

| Family | Subject | Entries |
|---|---|---|
| F1 | solar-system geometry beyond the current nine | 40 |
| F2 | beats, envelopes, moire, harmonic ladder | 26 |
| F3 | downloaded geophysical and space-environment series | 40 |
| F4 | catalog-endogenous | 65 |
| F5 | cross-regional and the 2R-df instrument's extensions | 20 |
| F6 | environmental and hydrological loading | 15 |
| F7 | the observer as part of the system | 15 |
| F8 | clocks and measures | 15 |
| F9 | statistics and instruments, including the second moment | 20 |
| F10 | combinations, with the class price table | 25 |
| | **total** | **281** |

*Counting convention, so a recount agrees: the file contains 260 entry headers, five of which cover ID
ranges or pairs (F1-21..F1-29 = 9 constituents, F2-18..F2-25 = 8 ladder entries, F3-12/F3-13,
F3-16/F3-17, F3-20/F3-21, F4-23/F4-24, F6-07/F6-08/F6-09). Counted by ID the total is 281.*

**Tranche-readiness counts:** READY 201, DOWNLOAD 59, INFRASTRUCTURE 21. (Sum 281. By header rather than
by ID: 193 / 54 / 13 = 260.)

**INFRASTRUCTURE builds named, not hand-waved:** site-local tidal stress tensor via SPOTL or direct Love
numbers plus FES2014 ocean loading (F1-30, and F1-17..F1-29, F1-31..F1-33 depend on it); unbinned
event-time scoring with a continuously evaluated baseline (F9-15, and F9-16 depends on it); a global
time-varying geodetic strain series (F8-06); a ComCat revision-history reconstruction (F7-10); a Schumann
resonance archive, which I could not name a free source for and labelled INFRASTRUCTURE rather than
pretend otherwise (F3-37).

**Free download sources named with URLs or portals:** GFZ Kp/ap/F10.7 (already wired), IERS finals.all
(already wired, polar-motion columns unparsed), WDC Kyoto Dst and AE, NASA OMNIWeb, SILSO sunspots, Oulu
and NMDB neutron monitors, NASA PO.DAAC GRACE mascons, NASA GES DISC GLDAS, Copernicus CDS ERA5, AVISO
altimetry and FES2014, NOAA IBTrACS, NOAA CPC climate indices, USGS waterdata, GRDC discharge, IGS/CDDIS
ionosphere maps, INTERMAGNET.

**Corpses explicitly dodged, each checked entry by entry:** static tidal-phase maps (no entry proposes
one; tidal entries appear only in the sub-daily, transient-conditioned, or regionally-incoherent slots
where the ~5.6% miner bound and the ~6.6% F-012 bound do not reach); calendar-time periodicity scans as a
standalone claim (periods here are fixed a priori or entered as declared clocks, never searched);
Fibonacci and ratio numerology (three entries are explicitly labelled MECHANISM-FREE null-calibration
lines and say so in their own text: F1-35, F1-36, F1-38, plus F2-08, F2-09, F6-07..F6-09); fault-type
parameter pools (F5-02 is a phase-coherence partition for a covariate test and states the distinction);
sequence-shape spatial transfer (not proposed); the TSI ratio (three entries flag ratio instability as
the failure mode that killed it and prescribe residual or log-difference forms instead: F2-16, F4-03,
F4-62, F10-12); the per-cell M>=4.5 battery (§P6-4 Finding B respected throughout; every spatially
resolved idea routes through the 2R-df sum or an aggregate, and F4-50 is offered as the honest survivor
of the dead design); the holdout-contaminated K-080 census cell list as a selector (Rule 4.1 restated in
§0.1 and in F4-48, F5-05); classical quiescence as a standalone precursor claim (M-007: family 4F entries
are covariates, aggregates and second-moment objects, and the spent hash `1afa6cdc...` is not touched);
anisotropic-versus-isotropic ETAS kernels as novelty (M-008: not proposed as a claim anywhere here);
California seasonal hydrological loading as novelty (M-010: family F6 opens by naming Sirorattanakul and
Avouac as the owners and restricts itself to global extension, transients, and conditioning);
`exp(dtau/Asigma)` as a power calculator outside 7-347 kPa (§P5-8: F1-32 names the VOID ruling and
requires an in-range link function per the S-17 candidate); site dilatation quoted better than a factor 2
(M-006: F1-30 and F8-06 both carry the warning); ML-found effects a boring model cannot reproduce
(Mignan and Broccardo: F10-23 is generator-only and says so).

**Standards respected and named in-file:** S-8 (max statistic), S-9 (one declared value, no alternatives
run), S-13 (structure-aware and magnitude-matched), S-14(c) (transient bracketing), S-15 (the
UNMEASURABLE class is mandatory and is scored neither way - every entry carries an S15 line, and 14
entries are declared UNMEASURABLE with a recommended price of zero rather than quietly dropped), S-16
(conditional offsets), S-17 candidate (estimators only in demonstrated range), G-M1 (band-matched and,
per §P6-4 item 5, aggregation-matched planted-signal recovery), Rules 4.1 through 4.7.

**Open questions I am handing to Popper rather than answering myself:** (1) is the deflation factor in
§0.3 clustering or bookkeeping, and F4-58 measures it; (2) is a clock a new test per feature or a declared
stratum in the BH structure (family F8 preamble); (3) is the leave-one-region-out check an audit or R new
tests (F5-20); (4) is the magnitude-floor ladder a post-hoc robustness check on survivors or a parallel
battery (F7-14); (5) is a transform of an existing feature (rank, ratio, cumulative dose) a new test, and
how much of F10-11 is analytically redundant with the existing lag scan.

**Status of this file: PROPOSED. Generated by Kepler, 2026-08-12. No result is claimed. Nothing here is
evidence. Popper prices; the supervisor runs the frozen tests; a K-entry registering this catalog comes
after Popper sees it.**
