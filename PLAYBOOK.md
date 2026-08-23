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
