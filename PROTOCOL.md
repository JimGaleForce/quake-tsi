# Frozen Replication Protocol — Out-of-Sample Test of the TSI "Predictability Class"

**Linear issue:** EQ-1 · **Status:** FROZEN as of 2026-07-20 (America/Los_Angeles)
**Rule:** No parameter below may change after regional event data has been downloaded and inspected,
except via a written amendment section appended to this file *before* the corresponding analysis is run.
Analysis code must read parameters from this file's companion `protocol_params.json`, not hardcode them.

## 1. Hypothesis under test

The paper "Conditional Tidal Triggering of Microseismicity … (TSI)" (Gale, 2025-11-25) reports, for
California: small (M < 5.0), early (t < 24 h) aftershocks on normal faults are enriched ~8x during
perigee-syzygy tidal windows, with timing governed by tidal volumetric strain rate |dεv/dt|, and
sensitive events concentrated in high-TSI (fluid-rich, extensional) crust.

**H1 (primary, confirmatory):** In extensional, fluid-rich regions outside California, aftershocks in
the class {M < 5.0, t < 24 h since mainshock, normal-faulting} occur during perigee-syzygy windows at
a rate exceeding the tidal-phase-randomized null (one-sided).

**H2 (secondary):** For events qualifying under H1's class, the percentile rank of |dεv/dt| at origin
time exceeds the percentile rank of fault-resolved Coulomb stress rate (paired, one-sided), i.e. the
rate/poroelastic mechanism, not Coulomb, governs timing.

**H3 (exploratory, not confirmatory):** Enrichment scales with local TSI (events in top-quartile TSI
cells show higher ER than bottom-quartile).

## 2. Test regions (chosen a priori — all predicted "susceptible" by the TSI model)

| Region | Box (lat, lon) | Rationale | Primary catalog | Secondary |
|---|---|---|---|---|
| Central Apennines, Italy | 41.5–44.5N, 11.5–14.5E | Normal faults + CO2-rich fluids; 2009 L'Aquila, 2016–17 Amatrice–Norcia–Campotosto | INGV FDSN | USGS ComCat |
| Taupō Volcanic Zone + adjacent extensional NZ | 40.0–37.0S, 175.0–177.5E | Rifting + geothermal | GeoNet FDSN | USGS ComCat |
| Iceland | 63.0–67.0N, 25.0–13.0W | Plate boundary rift + geothermal | USGS ComCat (IMO/SIL access is a follow-up amendment if obtainable) | — |

Negative-control region (H1 should FAIL here; strike-slip, low-TSI): **central San Andreas /
Parkfield-Cholame corridor** 35.0–37.0N, 121.5–119.5W, strike-slip events only, same era. A "hit" in
the control at comparable ER weakens the mechanism claim regardless of test-region outcomes.

## 3. Data windows and floors

- Time span: 1990-01-01 to 2026-06-30 (regional catalogs); completeness era begins per region:
  Italy 2005-04-16 (post-INGV network upgrade), NZ 2003-01-01, Iceland 1995-01-01 (ComCat, M≥4.0 only).
- **Primary magnitude floor: M ≥ 2.5** (regional catalogs, Italy & NZ). This *extends* the paper
  (whose floor was 4.5) into true small-event territory where the poroelastic mechanism predicts a
  *stronger* signal; it is frozen here as the primary test because power at M≥4.5 alone is inadequate
  outside California.
- **Comparability floor: M ≥ 4.5** run identically as a secondary readout (matches the paper).
- Depth ≤ 30 km. Events lacking depth or magnitude are dropped.

## 4. Declustering and event classes (identical to paper)

Gardner-Knopoff (1974) windows, exactly:

- T(M) = 10^(0.032·M + 2.7389) days
- L(M) = 10^(0.1238·M + 0.983) km

Implementation: pyCSEP's `gardner_knopoff` as reference; a second independent implementation must
agree on ≥99% of labels (this doubles as the EQ-2 audit on the new regions). Classes: mainshock,
aftershock (linked to a larger preceding event within its window), foreshock, background.
**Analysis class:** aftershocks with (t − t_mainshock) < 24 h and M < 5.0.

## 5. Faulting-style assignment (in priority order; first available wins)

1. Regional/global moment tensor (INGV TDMT, GeoNet MT, GCMT): normal if BOTH nodal-plane rakes are
   in [−135°, −45°].
2. No tensor: associate to nearest mapped active fault within 5 km horizontal (Italy: DISS/ITHACA;
   NZ: GNS Active Faults DB; Iceland: rift-segment normal faults) and take its style.
3. Neither: excluded from fault-specific tests (still counted in all-mechanism totals).

## 6. Tidal windows and forcing computations

- Ephemeris: JPL DE440 (as in paper) via `skyfield` or `astropy`.
- **Syzygy window:** |t − t_newmoon| ≤ 36 h or |t − t_fullmoon| ≤ 36 h.
- **Perigee window:** |t − t_perigee| ≤ 48 h.
- **Perigee-syzygy window:** syzygy AND perigee windows overlap the event time.
- ⚠️ *Open item recorded at freeze time:* the paper text does not state the original window
  half-widths. If the original Kosmos run's definitions are recovered, an amendment may substitute
  them **before unblinding** any regional result; otherwise the above defaults stand.
- Solid tide stress/strain tensors at hypocenter: PyGTide/SPOTL, as in paper; |dεv/dt| by central
  difference at 10-min sampling.

## 7. Null model and test statistics (the load-bearing part)

Aftershock rates decay (Omori); a uniform-time null is invalid. The null is built by **tidal-phase
randomization**: recompute the perigee-syzygy mask after circularly shifting the ephemeris time axis
by δ ~ Uniform(0, 29.53 d), independently 10,000 times, leaving all event times untouched. This
preserves the sequence's temporal structure exactly and randomizes only tidal phase.

- **Enrichment Ratio:** ER = (observed in-window count) / (mean in-window count over the 10,000 shifts).
- **p-value:** one-sided exceedance fraction across shifts (add-one corrected).
- **Primary decision rule:** H1 is supported if, in the pooled test regions (Italy + NZ + Iceland,
  event-weighted), ER > 1 with p < 0.05, AND at least one individual region shows ER ≥ 1.5 with
  p < 0.05. Report all per-region ERs with 95% CIs regardless of outcome.
- H2: Wilcoxon signed-rank, one-sided, on percentile(|dεv/dt|) − percentile(Coulomb-rate) with a
  **calibration check**: the same statistic computed at the 10,000 phase-shifted surrogate times must
  center on zero, else the percentile construction is biased and H2 is reported as inconclusive.
- Multiplicity: exactly ONE primary endpoint (above). Everything else is labeled secondary/exploratory.
  No subgroup fishing: the class definition in §4 is fixed.

## 8. TSI computation

TSI = ε̇vol / max(ε̇shear) on a 0.1° grid from NGL MIDAS IGS14 velocities, interpolated with the same
regularized inversion as the paper. Iceland/NZ/Italy all have adequate GNSS density. Grid built ONCE,
before event analysis, and hashed into the download log.

## 9. Blinding & audit trail

1. This file + `protocol_params.json` are hashed (SHA-256) and the hash recorded in
   `download_log.md` BEFORE the first regional event query is issued.
2. Catalog downloads are scripted (`download_catalogs.py`); raw responses are stored unmodified in
   `data/raw/` with retrieval timestamps.
3. The analysis runs blind per region; pooled unblinding happens only after all three regions complete.
4. Every deviation gets an amendment entry below, dated, before the affected step runs.

## 10. Interpretation commitments (written before seeing data)

- **All regions null:** the California signal was most likely a forking-paths artifact; the paper's
  California section is downgraded to hypothesis-generating and the global-null + TSI-mapping
  contributions stand on their own. This outcome is publishable and we commit to publishing it.
- **≥1 region replicates (per §7 rule):** TSI predictability class is real; proceed to EQ-5
  (true microseismicity) and EQ-7 (unified forcing) with the same frozen machinery.
- **Control region "replicates":** mechanism claim is unsupported even if test regions hit;
  investigate window/null construction for artifacts before any positive claim.

## Amendments

### Amendment 1 — 2026-07-20 (pre-unblinding; original Kosmos definitions recovered)

Source: Kosmos run `quakes4` (playground.edisonscientific.com/kosmos/7e2d1b15-…), subtask **r89**
"Identification of Tidally Sensitive California Aftershocks" — the task that produced the paper's
sensitive-event dataset. Exercising the §6 reserved amendment: the original run did NOT use time
half-widths; it used phase-angle/distance thresholds. To replicate the original class definition
faithfully, the following **replace** the corresponding frozen defaults. No regional event data has
been analyzed at amendment time (raw catalogs downloaded but untouched by any analysis code).

1. **Syzygy criterion (replaces §6 ±36 h):** lunar phase angle < 30° or > 150° (astropy convention,
   0°=new, 180°=full).
2. **Perigee criterion (replaces §6 ±48 h):** normalized lunar distance < 0.25, where
   normalized = (r − 356,500 km) / (406,700 − 356,500) clamped to [0,1] (0=perigee, 1=apogee;
   global normalization, matching the Kosmos `earthquakes_with_lunar_parameters.csv` construction).
3. **Perigee-syzygy window:** both criteria simultaneously true at event origin time.
4. **Mapped-fault association distance (replaces §5 item 2, 5 km):** ≤ 10 km (matches r89's
   cKDTree join against USGS Qfaults normal faults).
5. **Aftershock–mainshock association (supplements §4):** in addition to GK class labels, r89
   linked each aftershock to its nearest *preceding larger* event within 100 km and took
   t − t_mainshock from that link; we adopt the same rule for the latency cut.
6. **Strain-rate percentile (H2):** percentile of |dεv/dt| computed within the 24-h window centered
   on the event at 1-minute resolution (astropy tidal potential gradient), matching r89.

The tidal-phase-randomization null (§7) is unaffected: surrogates now recompute the phase-angle /
normalized-distance mask under circular ephemeris shifts, identically.

Also recorded from r89 for the audit trail (affects interpretation, not parameters):
- The original "sensitive" set is **n=14**, of which **12 are from the 1986 Chalfant Valley
  sequence** and 13/14 from the Long Valley–Chalfant corridor → effective sample ≈ 1–3 sequences.
- r89's own Wilcoxon (median percentile > 50%) was p = 0.034 — the paper's stronger p = 0.0002 /
  "100% of 13" figure comes from a different subtask comparing vs Coulomb-rate percentiles.
- Kosmos's internal dataset notes state that in multivariate spatial models (r82/r83),
  **max_shear_rate alone — not the TSI ratio — was the primary significant spatial predictor** of
  tidal sensitivity. H3 will therefore be evaluated BOTH ways (TSI ratio and max-shear-rate),
  pre-declared here.
