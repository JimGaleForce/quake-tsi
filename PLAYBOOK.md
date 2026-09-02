# THE PLAYBOOK — an ordered queue, not a wish list

**What this is.** The Kepler rounds generate freely; this file decides what actually runs and
in what order. It exists because a ranked list of 31 ideas is not an executable plan, and
because the program has repeatedly spent effort on the arm that was easiest to reach rather
than the arm with the most information in it.

**Read with `INVENTORY.md`.** Every entry names the data it needs and whether we hold it.
Nothing enters the queue without that line.

---

## THE ORDERING PRINCIPLE

Rank by **information per unit of holdout risk**, which decomposes into three multipliers this
program can actually compute:

| multiplier | what it buys | where we stand |
|---|---|---|
| **N** | floor ∝ 1/√N | 7,139 used vs **95,000 available**. Free 3.6× in floor. |
| **Amplitude** | threshold crossing, not linear gain | body tide ~1 kPa used; **ocean-loaded series on disk**, 5–20 kPa |
| **Dilution (f)** | a whole-catalogue test measures s·f, not s | f ≈ 0.05–0.15 whole-catalogue; ≈1 in near-critical populations |

**These multiply.** The reason every arm so far returned null is not that the hypothesis died;
it is that 5% floor × 0.1 dilution means we were bounding true susceptibility at ~50%, and
nobody claims 50%. **An arm that does not move at least one multiplier is not worth running.**
That single rule retires most of the ranked list.

---

## STANDING RULES (violating one voids the arm, not just the claim)

1. **Unit of inference is the SEQUENCE, not the event**, wherever clustered events are kept.
   Cost of ignoring it, measured 2026-08-22: max |z| went 26.27 → 2.51 when the null was fixed.
2. **Waveform-matched dwell null, never uniform-phase**, for every level/phase property.
3. **Null window ≥ one full period of the slowest feature in the arm.** The ±10 d window is
   1.35 fortnightly cycles and is fatal to any amplitude statistic. Use 14.765294 d half-window.
4. **Multiplicity paid over EVERY cell**, one max-statistic, no sub-family quoted as if declared.
   Report the **pooled-across-regions** max as well: the cell max is provably blind to a weak
   effect present everywhere at once.
5. **Sensitivity measured before a null is quoted.** A statistic that cannot move is a broken
   instrument, not a bound. Measure at the MINIMUM over nuisance phase.
6. **Holdout is spent once.** Exploration split for everything below unless the entry says
   otherwise.
7. **Name the artifact that could fake the result BEFORE the run**, in the module docstring.

---

## TIER 0 — INFRASTRUCTURE. Nothing above it can run without these.

| # | build | why | data | status |
|---|---|---|---|---|
| **P-0.1** | **Unified per-event covariate layer** over QTM/SCSN + all forcing series, on the existing `engine/properties.py` join (which already enforces `property_class` and refuses a property with no null layer attached) | every entry below consumes it; building it once is the difference between a playbook and a pile of scripts | §1a, §2, §3, §4 | **NEXT** |
| **P-0.2** | **Sequence/cluster ID layer** for every catalogue, cached | rule 1 is unenforceable without it | all catalogues | with P-0.1 |
| **P-0.3** | Load and validate the **Lu/Xue ocean-loaded stress series** against `calc_stress.py`; resolve the nanostrain-vs-Pa and dt=6000 s corrections | the only ocean-loading capability we have; unusable until units are pinned | §1a | with P-0.1 |
| **P-0.4** | **Per-event focal mechanism join** (SCEDC `.focmec` → strike/dip/rake) | resolved Coulomb on the plane that ACTUALLY SLIPPED instead of a regional guess; removes the largest approximation in every fault-relative arm | §4 | after P-0.1 |

---

## TIER 1 — RUN THESE. Each moves at least one multiplier, and all data is on disk.

### P-1.1 — THE HIGH-N REPLAY. *Moves N by 13×.*
Re-run the frozen battery on `QTM_decluster_m0.1` + `SCSN_decluster_m1.5` (~95,000 declustered
events) instead of 7,139. **This is not a new hypothesis** — it is every hypothesis we have
already tested, given an instrument that can finally see the claimed effect size. Floor
1.4% vs a literature claim of ~1%.
*Gate:* Mc(x,t) is not measured, so a tidal modulation of DETECTION at M≥0.1 is a live artifact.
Run the arm restricted to M≥1.5 (well above completeness) alongside as the control; if the two
disagree, the difference IS the detection artifact and must be reported as such.
*Either outcome is a result:* a measured effect, or the tightest tidal bound in the literature.

### P-1.2 — OCEAN-LOADED STRESS AT SOCAL. *Moves amplitude by 5–20×.*
Repeat P-1.1 using the SPOTL/TPXO **ocean-loaded** series rather than the body tide alone.
If the response is threshold-like, this is the difference between sitting under a threshold and
crossing it. Requires P-0.3.
*Positive control declared in advance:* Coso geothermal modulation, per the dataset's own rules.

### P-1.3 — FLUID-DRIVEN REGIONS. *Moves dilution: f → ~1.*
The `k034` set (Geysers 34,546; Long Valley 29,241; Coso 22,320; Salton 13,998 …) is the
low-effective-stress population. **Tidal modulation of tremor and LFEs is established and
order-unity**; these are the regular-seismicity populations physically closest to that regime.
Sequence-level nulls mandatory. Coso doubles as the positive control.

### P-1.4 — THE MACHINE-LEARNED ARM. *Escapes multiplicity entirely.*
Train on the exploration split over the FULL covariate layer; score **one number** on the
holdout: out-of-sample log-likelihood gain versus a baseline ETAS **with clustering intact**.
All searching happens in training, where nothing is claimed; the test is a single look, so
there is no multiplicity to pay however vast the space searched.
*Two constraints, both binding:* (a) baseline must be full ETAS or the model just rediscovers
aftershocks and looks brilliant; (b) a win says *something* in the covariate set carries signal
without saying what, so a positive requires a separate attribution stage.
*This is the direct answer to "failure of N": binned statistics compress thousands of events
into one number and throw the rest away; a per-event model uses the whole vector and can
represent interactions nobody declared.*

### P-1.5 — LOD / POLAR MOTION KILL-OR-CONFIRM. *Cheap, one test, contested claim.*
We hold IERS EOP with LOD since 1962. Lag frozen at the published value, **no lag search**,
block-bootstrap null preserving multi-year rate autocorrelation.
*Blocked on:* ISC-GEM global M≥7 back to 1900 (small download). Our world catalogue starts
1995, which is ~2 decadal cycles — underpowered but not useless as a first look.

---

## TIER 2 — WORTH RUNNING, AFTER TIER 1 REPORTS.

| # | entry | multiplier moved | data |
|---|---|---|---|
| P-2.1 | **Per-event-mechanism resolved Coulomb** (P-0.4) rerun of the fault-relative arm | removes geometry approximation | on disk |
| P-2.2 | **Magnitude as the responder** — b(phase), not rate(phase). Every statistic we own is an OCCURRENCE statistic; E[rate]=0 says nothing about the size distribution | new axis entirely | on disk |
| P-2.3 | **World M≥4.5 instead of M≥5.0** — a free 1.6× in N that nobody took | N | on disk |
| P-2.4 | **Matched placebo arm** (Kp/OMNI) — publish the rate at which the pipeline manufactures survivors on a class with real structure and no physical prior | calibrates everything else | on disk |
| P-2.5 | **Nonlinear harmonics** — a threshold response MUST generate overtones at frequencies the forcing does not contain, where the dwell null is empty by construction | artifact-immune by design | on disk |
| P-2.6 | **The negative-space catalogue** — thin the fitted ETAS intensity to draw phantom events, then two-sample test real vs phantom over the whole covariate block; one AUC, multiplicity paid once, and it captures SUPPRESSION which occurrence statistics under-weight | omnibus | on disk |

---

## TIER 3 — BLOCKED ON DOWNLOADS. Listed so the downloads get done in the right order.

| # | entry | needs |
|---|---|---|
| P-3.1 | Repeaters tested at their due date (f ≈ 1 by construction; the cleanest possible dated per-fault prediction) | Parkfield HRSN, Igarashi Japan repeater catalogues |
| P-3.2 | Slow-slip-conditioned regular seismicity (the best exogenous conditioning variable on Earth; Merton's named gap) | Cascadia/Guerrero/Boso/Hikurangi SSE catalogues |
| P-3.3 | Storm surges as step-function loading (event-study design, needs no dwell null at all) | UHSLC/GESLA tide gauges |
| P-3.4 | Outer-rise events (max amplitude, max criticality, known geometry) | GCMT + global load model |
| P-3.5 | Hydrologic loading, heavy-vs-light years | ERA5-Land SWE, GRACE |
| P-3.6 | Groningen Aσ ruler — measures the ONE unknown constant, converting every bound into a physics statement | KNMI Groningen catalogue + NAM pressure |
| P-3.7 | Mc(x,t) field from station noise — gates P-1.1 and five other entries | IRIS/EarthScope PSDs |

---

## RETIRED — do not re-propose without moving a multiplier

- Any further re-slicing of the 7,139-event world M≥5.0 declustered set with new angular or
  amplitude statistics. **3,000 cells have now been run against it** (`exp_mass_screen.py`,
  family p = 0.40, zero candidates). That well is dry, and dry for a computable reason.
- Per-region factor census / "why is Chile special". The dispersion gate found τ² = 0: there is
  no regional scatter to explain. *(Caveat: the gate's own p was overstated ~1000× by assuming
  independent cells; corrected effective d.o.f. 86 not 176, p = 0.032 not 6e-5. The conclusion
  stands, the alarm does not.)*
- Uniform-phase nulls of any kind.

---

## OPEN LOOSE ENDS

**Status refresh 2026-09-02 (supervisor; details in HYPOTHESIS_LEDGER.md "SESSION 2026-09-02").**

- **The 2026-08-22 near-miss is GONE.** `shallow_lt70km | areal.bot_decile` (pooled z -3.70, p 0.064)
  was an artifact of the 9,131-day epoch defect (CORRECTIONS.md 2026-09-02). After the fix it is
  outside the top 15 pooled keys. **Do not spend the holdout on it.**
- **Ran and null tonight:** P-2.2 (`exp_bvalue_skill.py`, bounds local delta-b < 0.2 only); P-2.4
  placebo (`exp_learned_ext.py`, survivor rate 0.052/feature, family-wise P(any) 0.53); Phase 2
  step 1 non-tidal blocks (same file, null); K-038 as forecast skill (`exp_dyntrig_skill.py`,
  +0.002 bits/event, bounds beta < 1.5 only; positive control reproduces); K-436 productivity
  assimilation (`exp_productivity_assim.py`, null; the naive first-hour count is incompleteness).
- **Re-run after the epoch fix (tidal readings re-derived, observer readings unchanged):** P-1.1,
  P-1.3, P-1.4 (+ sensitivity), `exp_nearcritical`, `exp_mass_screen`, `exp_diurnal_discriminator`.
  P-1.4 stays null (dAUC -0.0015, p 0.83, MDE eps 0.05).
- **Two EXPLORATORY leads, neither claimed:**
  1. `exp_neural_tpp.py`: a GRU point process beats ETAS+diurnal by **+0.040 bits/event on QTM
     M>=1.0, block CI [+0.029, +0.054]** (under its frozen 0.05 bar), +0.066 on SCSN M>=2.5;
     leakage calibration passes (loses to ETAS on ETAS sims). Merton: EarthquakeNPP found no NPP
     beating ETAS, so CONTESTED. Next step, if Popper admits it: spatial ETAS baseline, untouched
     window, attribution stage (what does the GRU know that ETAS does not: productivity scatter?
     Omori p variability? magnitude-dependent c?).
  2. `exp_fluid_driven.py` (corrected times): pooled over the 7 fluid-driven regions,
     `moon_sinel_abs.mean` z 3.50, p 0.022 against its own max-null; the 6 tectonic regions pool
     at p 0.77. Declared 7/7 test null (p 0.46). Not claimed; one of several secondary families.
     **Withdrawn by Popper Round 7 (§P7-26(3)):** the sentence "this is the Lu et al. direction".
     `moon_sinel_abs` is |sin(lunar elevation)| (`exp_mass_screen.py:177`), an ephemeris scalar, not
     the fault-plane shear stress Lu et al. test. The pool is also STRONGER under the full-sequence
     null (p 0.0007) than declustered (p 0.022), which under SP-2 reads as clustering-derived. The
     EQ-18 cross-test is still worth pre-registering, for a different reason (the Coso Fig 4c
     detector works); its freeze list is in §P7-26(3).
- **New dataset:** `data/global_m55/` (INVENTORY §1f), 18,769 worldwide M>=5.5 triggers 1985-2022.
- Unexplained at 2.15 sigma (declustered per-region under-dispersion): unchanged, logged, not built on.
- Fable audit findings 3, 4, 6, 7 still open. SP-7 v2 ratification and the K-108 non-independence
  ruling: DEFERRED by Popper Round 7 (§P7-26(7)), top of the Round 8 docket ahead of new arms.
  Popper's Round 7 rulings on every 2026-09-02 arm are summarised below.

**Popper Round 7 (2026-09-02, ledger §P7-26; the first ruling on any of the above).**

- Correction ACCEPTED, one residual exposure: `exp_world_faultrelative.load_region` was a fifth
  loader outside `tests/test_epoch_consistency.py`; the mass-screen re-run needed an
  `engine/EXPLORE_COUNT.jsonl` declaration under SP-1.5. Both being closed this session.
- Arm A NULL at poor sensitivity (bounds beta < 1.5 at D0 = 10 kPa only). Arm B primary is a
  NULL-OF-A-BROKEN-RULE (zero power at every injected delta; not quotable as a bound), but arm B
  **control (c) is PROMOTED to PROVISIONAL (negative)**: a local b DESTROYS 0.48 to 0.88 bits/event
  at M>=2.5 on both catalogues, pre-registered, same sign, far outside CI. Arm C admitted as a
  **PROVISIONAL-LEAD, not reportable**. Arm D null at weak sensitivity; its placebo rate
  (0.052/feature) is accepted as a calibration. Arm G RETIRED AS POSED; its ledger prose claiming
  the alpha-refit baseline "absorbs" the gain is WITHDRAWN (the JSON says A2 is worse than A by
  0.013 bits). P-1.1 observer-confirmed, retired as a tidal arm. P-1.4 is the program's best tidal
  bound (MDE eps 0.05 over the whole covariate space). Holdout candidate from 2026-08-22 RETIRED.
- **Arm C proof path, frozen at §P7-26(2), three stages, nothing skips.** Stage 1: the ETAS-repair
  ladder R0 (alpha free, lb 0.05) -> R1 (hierarchical per-sequence productivity) -> R2 (two-component
  Omori) -> R3 (magnitude-dependent c) -> R4 (fitted Mc(t) incompleteness layer, MANDATORY), on the
  exploration split only, sequence-block CI, positive control (injected productivity scatter must be
  recovered by R1) and negative control (Poisson through every rung) BEFORE the primary is read.
  Stage 2: one scored look on arm C's OWN untouched holdouts (24,532 QTM, 10,604 SCSN), hash into
  `engine/HOLDOUT_LOG.jsonl` first. Stage 3: 2019+ ComCat (8,256 events at M>=2.5, not ~3,000 as
  first estimated), only after a sensitivity run: arm G already measured that a sequence-block CI on
  a SoCal window of this shape gives 0/12 CIs excluding zero at an injected +0.037. The spatial
  `EtasV1` leg is REMOVED (cell-day grid vs per-event continuous-time likelihoods are not
  comparable, and a spatial baseline changes the target). 2019+ reservation by §P3-6 is not
  exclusivity; every 2019+ report must state the cumulative count of independent scored looks.
- K-431..K-448: **K-448, K-440, K-434 elevated** (the last two are the round's named ETAS repairs).
  K-441 decisive but not one evaluation (O(2N+2) fits) and its Mc(t) rider is binding. K-444
  reframed (K-444R; the bias control is the experiment). K-442 contradicted by the same night's
  diurnal reading and its covariates are post-detection (leakage). Full table at §P7-26(4).
- **Priority (Popper):** 1. Phase 4 prospective log NOW, v1 = frozen ETAS vs Poisson only ("every
  week it is not running is unrecoverable evidence"); 2. Phase 2 stage 1 in parallel; 3. K-441 (with
  rider), K-443, K-446; 4. arm C stage 2 only if stage 1 passes; 5. EQ-18 cross-test when the bin
  list arrives; 6. K-440 + K-434; 7. arm C stage 3 after its sensitivity run.

**Retired tonight (do not re-propose without moving a multiplier):** K-436 as posed with naive
first-hour counts (incompleteness); K-038's marginal-effect form (its ledger-class interaction,
the actual claim, is untested and remains open).

<details><summary>Loose ends as of 2026-08-22 (superseded above)</summary>


- **Near-miss on the books:** mass screen pooled `shallow_lt70km | areal.bot_decile`,
  z = −3.70 across 11 regions, p = 0.064. Deficit of shallow events in the most compressive
  decile; sign is what a mechanism predicts. **Not claimed.** It is a single pre-specified
  statistic with a declared sign, which makes it the cleanest possible holdout candidate if
  Jim chooses to spend one.
- **Unexplained at 2.15σ:** declustered per-region z's under-dispersed (variance 0.641). Pipeline
  cleared as the cause — an independent second null ensemble returns 0.9759 ± 0.0099 where 1.000
  is expected — so it is a property of the data. Logged, not built on.
- Fable audit findings 3, 4, 6, 7 still open.
- Popper has not ratified the SP-7 v2 amendment; K-108 non-independence ruling outstanding.
</details>
