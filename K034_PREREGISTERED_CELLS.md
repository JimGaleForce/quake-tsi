# K-034 — PRE-REGISTERED TARGET CELLS, CLASSES, RANKING RULE, SOURCES, WINDOWS, STATISTICS

Written and hashed **before** any K-034 statistic was computed, per Popper R2-2 mandate (3)
("'Geothermal/volcanic areas light up first' must be a named, ranked list of cells committed
before unblinding, or it is post-hoc pattern-matching wearing a pattern prediction's clothes").

Same analyst-authored caveat as `K034_SEALED_LITERATURE.md` §Preamble applies.

## §1. Sources (R2-2 mandate 2: >= 2 of these four must fire)

| id | event | origin (UTC) | lat | lon | Ms (frozen) | Mw | rupture L (km) | gate = 2L (km) |
|---|---|---|---|---|---|---|---|---|
| landers | Landers, CA | 1992-06-28 11:57:34 | 34.200 | -116.436 | 7.3 | 7.3 | 70 | 140 |
| hectormine | Hector Mine, CA | 1999-10-16 09:46:44 | 34.594 | -116.271 | 7.4 | 7.1 | 48 | 96 |
| ridgecrest | Ridgecrest, CA | 2019-07-06 03:19:53 | 35.766 | -117.605 | 7.1 | 7.1 | 50 | 100 |
| denali | Denali, AK | 2002-11-03 22:12:41 | 63.517 | -147.444 | 8.5 | 7.9 | 340 | 680 |

Magnitudes and rupture lengths other than "Landers Ms 7.3" (VERIFIED, Hill et al. 1993 title) are
**PARTIALLY VERIFIED** and are used only to set the amplitude axis and the distance gate. A
sensitivity using Mw in place of Ms on all four sources is reported.

## §2. Target cells

Class **A** = caldera / producing geothermal field; class **B** = volcanic or geothermal-adjacent;
class **C** = non-geothermal tectonic (the negative-pattern controls). Class is assigned here,
before unblinding, from the geology, not from any rate statistic.

| cell | class | lat_min | lat_max | lon_min | lon_max | documented Landers-triggered? |
|---|---|---|---|---|---|---|
| long_valley | A | 37.40 | 37.90 | -119.10 | -118.60 | YES (VERIFIED) |
| coso | A | 35.80 | 36.20 | -118.00 | -117.60 | YES (PARTIALLY VERIFIED) |
| geysers | A | 38.70 | 38.92 | -122.95 | -122.70 | YES (VERIFIED) |
| yellowstone | A | 44.30 | 44.90 | -111.20 | -110.40 | YES (VERIFIED) |
| salton_brawley | A | 32.90 | 33.30 | -115.80 | -115.40 | not named in the sealed set |
| lassen | B | 40.30 | 40.70 | -121.70 | -121.30 | YES (PARTIALLY VERIFIED) |
| mono_west_nv_mina | B | 38.20 | 38.60 | -118.40 | -117.90 | YES (VERIFIED, "Mina, 500 km") |
| little_skull_mtn | B | 36.60 | 37.00 | -116.50 | -116.00 | YES (VERIFIED, "280 km, M5.6") |
| cedar_city_ut | B | 37.40 | 37.90 | -113.40 | -112.80 | YES (PARTIALLY VERIFIED) |
| smith_valley_nv | C | 38.60 | 39.00 | -119.60 | -119.20 | YES (VERIFIED, "590 km, M3.4") |
| parkfield | C | 35.70 | 36.10 | -120.70 | -120.30 | not named |
| mendocino | C | 40.20 | 40.70 | -124.70 | -124.10 | not named |
| wasatch_slc | C | 40.40 | 41.00 | -112.20 | -111.60 | not named |
| san_jacinto | C | 33.20 | 33.70 | -116.90 | -116.40 | not named (near-field for landers/hectormine; gate will exclude) |

## §3. RANKING RULE (the pattern prediction, committed)

For each source, cells passing the distance gate are ranked by the key
`(class_rank, -sigma_dyn_primary)` with `class_rank(A)=0, class_rank(B)=1, class_rank(C)=2`.
The prediction P3 is that the **observed** rate-response ranking correlates positively with this
pre-registered ranking (Spearman rho > 0, one-sided), and specifically that class A+B outranks
class C.

## §4. WINDOWS AND STATISTIC (frozen)

- **Post window:** `[t0, t0 + 5 d]`. Secondary window reported but not headline: `[t0, t0 + 1 d]`.
- **Background:** the 90 days ending at `t0` (i.e. `[t0 - 90 d, t0)`), same cell, same catalogue,
  same magnitude cut.
- **Magnitude cut:** M >= 1.5 for all cells (below the M2.5 quarry-contamination concern, which is
  a SoCal daytime-blast issue; a M >= 2.5 sensitivity arm is reported).
- **Statistic:** `RR = N_post / (lambda_bg * 5 d)` where `lambda_bg = N_bg / 90 d`.
- **Distance gate:** cells with epicentral distance < 2 rupture lengths from the source are
  **excluded, not scored**, per Kepler's entry text and R2-2 mandate on K-038(2).
- **Null (primary):** circular time-shift. 999 pseudo-origin times drawn on a regular circular
  grid over `t0 +/- 3 yr` in the same cell's catalogue, excluding `|dt| < 30 d`; the identical
  statistic is recomputed at each. `p = (1 + #{RR_null >= RR_obs}) / (1 + n_null)`. This null
  inherits the cell's true burstiness/overdispersion, so it is strictly more conservative than
  Poisson.
- **Multiplicity (S-8):** the declared family is {sources x gated cells x 2 windows x 2 magnitude
  cuts x 3 A*sigma points}. Family-wise significance is by the **max-statistic** of the same
  circular-shift null, computed over the declared family. Both per-cell p and family-wise p
  are reported; **the PASS rule uses the family-wise number.**
- **Aftershock / ETAS contamination:** handled by (i) the distance gate, (ii) the shift null
  being drawn from the same cell so that the cell's own background clustering is in the null,
  and (iii) reporting, for each fired cell, whether it lies inside 2L of the source.

## §5. POWER / DETECTION-THRESHOLD CURVE (the deliverable)

For each gated cell, from its observed `lambda_bg`:
1. Monte-Carlo the shift-null distribution of `RR` (as above).
2. For a grid of true rate multipliers `R` in [1.0, 20], simulate `N_post ~ Poisson(R * lambda_bg *
   5 d)` and compute detection probability at alpha = 0.05 against the cell's own null.
3. `R_min` = smallest `R` with power >= 0.80.
4. Convert to stress via rate-state: `dTau_min = A*sigma * ln(R_min)`, at
   **A*sigma in {0.03, 0.10, 0.15} MPa** (S-14). The reported **amplitude floor** is
   `dTau_min` at the **adverse end, A*sigma = 0.15 MPa**, per S-14's rule that the adverse end
   sets the flag.
5. Power for each *claim* is reported at all three A*sigma, evaluated at the cell's own
   predicted `sigma_dyn` from the frozen amplitude model.

## §6. PASS / FAIL RULE (committed)

**PASS** iff the engine fires (family-wise p < 0.05) in at least one documented cell for at least
**2 of the 4 sources**, AND the class A+B vs class C pattern is in the predicted direction.
**FAIL** otherwise, in which case the injection-recovery arm (§5) adjudicates
"engine broken" vs "catalogue too small".
