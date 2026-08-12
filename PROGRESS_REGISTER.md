# Progress Register — the program's shareable record

*Maintained by **Faraday** (`.claude/agents/faraday.md`), the representative/publisher persona.
Opened 2026-08-09. **Round 2 revision 2026-08-11** (cycle-2 step 4). One entry per candidate claim.
Nothing here is a claim until it carries both a Popper standing and a Merton attribution block;
everything else is visible but not claimable.*

**WHAT CHANGED THIS ROUND, IN ONE PARAGRAPH.** Five new entries: **F-012** (K-035, the
corpse-to-bounds conversion — the single most publication-relevant artifact the program owns),
**F-013** (K-034, the Landers positive-control gate, PASS-qualified, certified floor ~34 kPa),
**F-014** (L-1, the S-14(c) two-branch transient link bracket), **F-015** (the EQ-23/EQ-24 engine as
an open methods artifact, superseding the aspirational half of F-010), **F-016** (the phase-clock
instrument response — the newest and, per unit of work, the most interesting thing here). Merton's
**M-009.12(a)** correction is folded into F-005 and everywhere forecast-skill context is quoted: our
+1.907 bits/event is **57% of the closest matched temporal-only comparator (SCEDC_25, 3.319
bits/event)**, and the old "ours is lower because theirs are spatio-temporal" defence is **withdrawn**.
Popper's **§P4-2** ruling on K-046 is folded into F-009. A **PUBLICATION TRIAGE** section is added at
the end — the ranked answer to "what can we publish", with the single next step for each candidate,
and the one-paper recommendation.

**How to read an entry.** CLAIM is the sentence a skeptical seismologist would accept as scoped —
region, magnitude floor, period, method, number. FRAMING (why it matters) is kept separate from
CLAIM on purpose, and no enthusiasm is allowed to cross that line. PROOF CHAIN gives real hashes
and real commits, verified by me against `download_log.md` and `git log` this session. NOT SHOWN
is as load-bearing as the claim. TIER is INTERNAL → PANEL-READY → PUBLISHABLE.

**Verification standard.** Every number below was re-read by me from the primary JSON, not from a
summary. Where the program record and the primary JSON disagree, the JSON wins and the
disagreement is recorded in §OVER-QUOTE AUDIT at the end. Round 1 found five discrepancies, two
material. **Round 2 (2026-08-11) re-read `results_k034.json`, `results_l1.json`, `results_k035.json`
and `weifan_phase_null.py` directly, re-executed the phase-null calibration, and found three further
discrepancies — one of them in an email that has already been drafted for an outside collaborator.**

---

## Register index

| id | claim (short) | class | Popper standing | Merton standing | tier |
|---|---|---|---|---|---|
| F-001 | Coso Fig 4c independent reproduction | REPRODUCTION | validated (frozen rule met) | attributed (direct) | **PANEL-READY** |
| F-002 | Strain-field two-velocity-set comparison; dilatation fragility | EXTENSION (candidate) → **REDISCOVERY** | provisional | attributed — REDISCOVERY (M-006) | INTERNAL |
| F-003 | Tidal-prediction null suite (forecast skill), SoCal | NEW APPLICATION (falsification) | validated + re-scoped R2-1(d) | attributed (own prior work) | **PANEL-READY** |
| F-004 | Pure-math temporal-pattern null suite (periodicity / ratios / shape transfer) | NEW APPLICATION (falsification) | mixed: G stands, F is weak | **PENDING-ATTRIBUTION** | INTERNAL |
| F-005 | B-2 SoCal walk-forward ETAS skill | REPRODUCTION (confirmed) | validated | attributed — REDISCOVERY (M-004) | INTERNAL |
| F-006 | B-1 generic-ETAS cross-region transfer | REPRODUCTION + candidate EXTENSION | validated (with S-5 correction) | attributed — REDISCOVERY (positive) + CONTESTED (negative) (M-005) | INTERNAL |
| F-007 | B-3 seven-day M≥5 conditional rule | REPRODUCTION | validated | attributed (M-003) | **PANEL-READY** |
| F-008 | B-4 stress-ledger negative space recovers known aseismic geology | REDISCOVERY + methods delta | validated, AUTO-FLAGGED | attributed (M-002) | INTERNAL |
| F-009 | K-009 ETAS residual field is not white in SoCal | REPRODUCTION + 1 candidate methods delta | success rule met, partial vs spec | attributed (M-001) | INTERNAL |
| F-010 | The research engine itself (frozen protocols, personas, ledger) | METHODS/PROCESS (candidate) | n/a — not a hypothesis | **PENDING-ATTRIBUTION** | INTERNAL |
| F-011 | Live ETAS forecaster + globe layer | ENGINEERING ARTIFACT (not a claim) | n/a | n/a | INTERNAL |
| F-012 | K-035: the tidal corpses re-priced as powered upper bounds | **EXTENSION** (methods delta on our own null) | validated (frozen mandates 1–6 met; a0 arm PASS) | attributed via F-003 + Beeler–Lockner theory line; **PENDING-ATTRIBUTION on the bounds-reporting norm** | **PANEL-READY** |
| F-013 | K-034: Landers-class remote dynamic triggering detected blind at documented sites; certified floor ~34 kPa | REPRODUCTION (positive control) + **licensing instrument** | validated PASS (qualified — one of four family readings fails) | attributed (sealed literature, Hill et al. 1993; vdE&B 2010) | **PANEL-READY** |
| F-014 | L-1: the standard transient link function is void; two-branch bracket | **EXTENSION** (methods) | provisional — arithmetic validated, link unclaimed | attributed (Dieterich rate-state; K-034 own data) | INTERNAL |
| F-015 | EQ-23/EQ-24 engine: holdout hash-gate, generator-not-evidence, planted-signal tests | METHODS/SOFTWARE ARTIFACT | n/a — not a hypothesis | **PENDING-ATTRIBUTION** | INTERNAL |
| F-016 | Anchor-based tidal phase clocks have a measurable non-uniform null response | **NEW APPLICATION** (candidate NEW — instrument characterisation) | **PENDING-VALIDATION** (no frozen protocol, no results JSON) | **PENDING-ATTRIBUTION** | INTERNAL |

Six of eleven entries were PENDING-ATTRIBUTION or attribution-thin at open. **Merton round 2
(2026-08-09) cleared three of them — F-005 (M-004), F-006 (M-005), F-002 (M-006).** Merton has now
produced **nine** dossiers (M-001 K-009, M-002 B-4, M-003 B-3, M-004 B-2, M-005 B-1, M-006 B-5,
M-007 quiescence, M-008 angular transfer, M-009 the seventeen-entry wave umbrella). **Four entries
remain PENDING-ATTRIBUTION: F-004, F-010, F-015, F-016** — and F-016 is the one where the missing
search is most likely to change the class, because it is the only entry in this register that could
plausibly come back **NEW**.

**Tier counts: 5 PANEL-READY · 11 INTERNAL · 0 PUBLISHABLE · 0 Merton-certified NOVEL.** The honest
headline for a supervisor asking "can we publish": *nothing in this register is submission-ready
today, two things are four-to-six weeks of writing and one verification away, and the register knows
exactly which two.* See §PUBLICATION TRIAGE.

---

## F-001 — Independent reproduction of Lu et al. (2025) Figure 4c at Coso

**CLAIM.** In the single 0.4° × 0.4° bin lat [36.2, 36.6], lon [−118.0, −117.6] (the bin used in
Figure 4c of Lu et al. 2025), using the SCSN catalog at M ≥ 1.5 with per-event YHS/SCEDC
focal-mechanism resolution of the tidal shear stress τ, we independently measure a phase
modulation of Pm/P0 = 0.340 at φ = −133° on n = 113 FM-matched events, one-sided empirical
p = 0.0407 against orientation-preserving time-randomized synthetics. Lu et al. report Pm/P0 =
0.46 on n = 112 in the same bin. In the same bin and the same events, the normal-stress component
σ_n shows Pm/P0 = 0.271, p = 0.12 — not significant. Under the pre-registered rule (one-sided
empirical p < 0.05) this is a **qualified reproduction**: same sample to within one event, signal
present in the claimed stress component and absent in the control component, amplitude 26% lower
than theirs.

**CLASS — REPRODUCTION.** This validates our machinery against a collaborator's published result
computed from their own parameters. It is not a discovery and must never be written as one. Its
value to this program is exactly that our pipeline, pointed at their bin with their recipe,
returns their sign, their component selectivity, and 74% of their amplitude.

**ATTRIBUTION.** Lu, W., Xue, L., Yue, H., Zhuang, J. & Zhao, D. (2025), *JGR Solid Earth*,
doi:10.1029/2025JB032249 — the result being reproduced, the method (§2.5 per-event
focal-mechanism stress resolution), the catalogs (Zenodo doi:10.5281/zenodo.18491845, CC-BY 4.0),
and the Fig 4c parameters, supplied directly by Weifan Lu by email 2026-08-06. Supporting: SCEDC
and the Yang–Hauksson–Shearer focal mechanism catalog. *Attribution for a reproduction is
complete by construction — we name the work we are reproducing.* Merton has **not** searched
whether anyone else has independently reproduced Fig 4c or measured phase modulation at Coso;
that search is not required for this claim but would strengthen the framing.

**PROOF CHAIN.**
- Protocol: `COSO_FIG4C_PROTOCOL.md`, SHA-256
  `aa97685b1000e7f5da0affc05ff81483a0d168a01191da385bb7faccd053bae0`, recorded in
  `download_log.md` **before** any analysis in this bin. **I verified this hash twice this
  session:** the working copy matches, and the file as committed at `5b4997a` matches. Clean.
- Public commit: `5b4997a` "Coso Fig 4c exact reproduction per Weifan Lu's parameters
  (2026-08-06)", pushed to github.com/JimGaleForce/quake-tsi.
- Code: `coso_fig4c_test.py` on `coso_fm_test.py` machinery, seed 20260806.
- Result: `results_coso_fig4c.json`.
- Controls: negative/component control — σ_n in the identical bin and events returns null
  (p = 0.12). Artifact class named: phase-assignment fidelity (their 6,000-s native series,
  ~7.5 samples/semidiurnal cycle, vs our CubicSpline ×10 upsampling) is the leading suspect for
  the amplitude gap, together with FM-catalog vintage (113 vs 112 matches) and declustering
  details. Cross-catalog arm attempted and correctly abandoned: QTM has 29 FM matches (20 at
  M ≥ 1.5) in this bin — too few to test.
- Replication status (Popper ladder): **PROVISIONAL** — one bin, one catalog, one run.

**THE CORRECTION THAT MUST TRAVEL WITH THIS ENTRY.** Three earlier "Coso null" results in this
program tested a *different, mostly non-overlapping* box (COSO_BOX, lat 35.60–36.25). Figure 4c's
bin is almost entirely north of it. **Those nulls do not bear on Lu et al.'s claim, and the
sentence "Coso is null in our hands" is retracted.** The wrong-box error was found by us, before
any reply was sent, and corrected by freezing a new protocol and re-running. Reporting this is not
optional; it is the reason the reproduction is credible.

**WHAT IT DOES NOT SHOW.** It does not show that tidal modulation at Coso is established at the
95% amplitude threshold: our observed 0.340 sits *below* the synthetic 97.5th-percentile amplitude
(0.362), and `results_coso_fig4c.json` carries `significant_95: false`. The reproduction rests on
the frozen one-sided empirical p, and any restatement must say so. It does not show the amplitude
gap is a sampling fluctuation rather than a method difference. It does not show anything about the
*cause* of the modulation, and it does not resurrect the TSI proxy (F-003). Note also that the
Fig 4c bin is centered between the Coso geothermal field and Owens Lake, so a
"geothermal-fluid" reading of this bin is not supported by the bin's location alone.

**HOSTILE-REVIEWER OBJECTION, PRE-ANSWERED.** *"Your amplitude is 26% below theirs and under your
own 97.5th-percentile threshold — this is a failed reproduction dressed up."* Answer: the
pre-registered success rule was the one-sided empirical p, frozen and hashed before the bin was
touched, and it was met at p = 0.041; the amplitude threshold is a second, stricter statistic we
report *against ourselves* rather than omit. The component selectivity (τ positive, σ_n null on
identical events) is the part that is hard to get by chance and it reproduces cleanly.

**TIER: PANEL-READY.** Shareable with Lu/Xue/Vidale/Bürgmann as-is, with the wrong-box correction
and the `significant_95: false` caveat in the same paragraph as the number.

**POWER, ADDED BY K-035 (2026-08-10) — and it is unflattering, so it goes in the entry.** From
`results_k035.json :: corpse_to_bound_table_K032_item6`, re-read by me: the Coso Fig 4c configuration
at n = 113 has a **minimum detectable modulation of 41% at 80% power** (`mda80 = 0.40875`,
bracketed in [0.4, 0.8] where power runs 0.77 → 1.00). Our observed 0.340 and Lu et al.'s 0.46 both
sit **at or below the 80%-power floor of the design**. That does not retract the reproduction — the
frozen rule was the one-sided empirical p and it was met — but it means **this bin cannot
discriminate 0.34 from 0.46, and no one in this program may treat the 26% amplitude gap as a
measured difference.** It is a difference between two numbers drawn from a design that cannot
resolve them.

**SINGLE MISSING VERIFICATION TO REACH PUBLISHABLE.** Re-run the phase assignment at the native
6,000-s sampling (no spline upsampling) and report the amplitude both ways — **and now, per F-016,
report the anchor-clock null response for this bin's own stress series alongside it.** If the gap
closes, the reproduction becomes quantitative and the method difference is measured rather than
suspected. Realistically this entry's ceiling is PANEL-READY as a standalone; its publication route
is as the machinery-validation section of the F-012 bounds paper or the F-016 methods paper.

---

## F-002 — Two-velocity-set strain comparison and the fragility of the dilatation component

**CLAIM (restated by me; see below for why).** Recomputing the same strain-rate estimator on two
independent GNSS velocity solutions for Southern California — NGL MIDAS (1,084 stations) and the
Kreemer & Young (2022) velocity compilation (941 stations after QC) — over 4,679 common grid
nodes: maximum shear strain rate agrees closely (Pearson r = 0.934, Spearman ρ = 0.928, median
|difference| 5.3 nstrain/yr), the second invariant likewise (r = 0.932), and the high-strain mask
(ε_min ≥ 47 nstrain/yr, the Kreemer & Young Fig 8 criterion) has Jaccard 0.830 between solutions.
**Dilatation does not** (r = 0.782, ρ = 0.741, median |difference| 8.1 nstrain/yr), and at
individual sites the two solutions disagree by factors from 2.0× to 83× and can flip sign
(Long Valley +146.0 vs +71.5; Coso −17.2 vs −6.7; Salton Sea +27.7 vs +4.3; Brawley +24.1 vs
+0.29; Cerro Prieto **−52.9 vs +10.7**, nstrain/yr). Consequence for this program: the sign of the
TSI ratio at the two sites that matter is robust across velocity solutions (Coso −0.43/−0.21,
Long Valley +6.3/+2.0), so the TSI retirement does not depend on the velocity field; but any
claim resting on dilatation *magnitude* at a single site does.

**Why restated.** Baseline B-5 is currently written as "the dilatation component carries **±2×**
measurement uncertainty across velocity solutions." **That number is not in
`results_strain_comparison.json` and is not derivable from it.** It matches Long Valley (2.0×) and
Coso (2.6×) and badly understates Salton Sea (6.4×), Brawley (83×), and Cerro Prieto (sign flip).
The defensible statements are the correlation coefficients, the median absolute difference, and
the site table. I recommend B-5's headline be amended to the form above. **This is over-quote
finding #4 in the audit below.**

**CLASS — EXTENSION (candidate).** Kreemer & Young themselves warn that GPS-only models produce
spurious dilatational artifacts along the San Andreas system; our delta would be a *quantified*
version of that caution — a same-estimator, two-velocity-set contrast with a numeric agreement
table on a fixed 4,679-node grid, plus the site-level fragility list. Whether that is a real delta
or a rediscovery of a known result is **Merton's to rule**, and he has not.

**ATTRIBUTION — RESOLVED: REDISCOVERY. Merton dossier M-006 (`HYPOTHESIS_LEDGER.md`).** My prior
was right and Merton confirms it. The corrected claim (§18 wording: "site dilatation unreliable to
worse than 2×, up to sign flips") is owned outright by **Maurer & Materna (2023), *GJI* 234(3),
2128–2142, doi:10.1093/gji/ggad191** — five strain methods on one SoCal GNSS field: *"The standard
deviation of the dilatation rate is as large or larger than the signal in many places"*, with
dilatation masked to "indistinguishable from zero" over much of the map, and **Cerro Prieto named
as the worst-conditioned area** — the exact site of our sign flip. Corroborated by Hearn et al.
(2010) (~100% between-method variability, SoCal), Xu et al. (2021) (cross-model strain correlation
→ 0 at ~30 km), Hackl et al. (2009), Baxter et al. (2011), Titus et al. (2011) (alternating-sign
dilatation as a known network artifact), Pagani et al. (2021), and by Kreemer & Young's own
spurious-dilatation caution. **The ε_min ≥ 47 nstrain/yr criterion is theirs — question closed,
cite every time.** The only unowned element is the *axis*: all prior work varies the estimator on a
fixed velocity field, we varied the velocity field on a fixed estimator — an axis Maurer & Materna
explicitly flag as an unquantified contributor. Merton rules that a methods footnote, not a
finding, and notes our two velocity sets share NGL/UNR processing (so our contrast is a **lower
bound**) and that the 83× Brawley ratio should be dropped (near-zero denominator).

**PROOF CHAIN.**
- Code: `strain_comparison.py`; result `results_strain_comparison.json`; figure
  `maps/strain_comparison_kreemer_young.png`.
- Public commit: `741a178` "Strain-rate comparison vs Kreemer & Young (2022) velocities; fix
  strain unit scale".
- **No frozen protocol.** This was a responsive analysis (Roland Bürgmann's 2026-07-29 ask), not a
  pre-registered confirmatory test. Recorded as such; it is a measurement, not a hypothesis test,
  and does not need pre-registration — but it also may not borrow the credibility of the entries
  that have it.
- Data provenance caveat, on the record: Kreemer & Young's gridded MELD/Haines–Holt strain model
  is **not deposited** (Dataverse doi:10.7910/DVN/BICMWB is the velocity table only; NGL,
  ScienceBase and the Maurer GJI supplement were checked). We therefore compare *velocity
  solutions through our estimator*, not our model against their model. Any reader who thinks this
  is a model-vs-model comparison has been misled and the entry must prevent that.
- Self-caught error in the same work: a ×1000 unit bug in `tsi_map.py` (gradients in
  microstrain/yr labeled nanostrain/yr) was found and fixed; npz/json regenerated. The figure
  previously shared with Bürgmann predates the bug and was correct. Logged, not hidden.

**WHAT IT DOES NOT SHOW.** It does not compare our strain field to Kreemer & Young's *model*
(theirs is unavailable). It does not establish which solution is right — only that they disagree,
and where. It does not bound dilatation error generally; the 4,679 nodes are one region, one
estimator, two solutions. It says nothing about time-dependent strain.

**TIER: INTERNAL.** Blocked on: (a) ~~Merton's dossier~~ **cleared, M-006**; (b) the B-5 headline
restatement above; (c) **a class decision that is now the supervisor's**: M-006 rules F-002 is
promotable to PANEL-READY **as a REDISCOVERY** citing Maurer & Materna (2023) and Kreemer & Young
(2022), *or* it waits for the `Strain_2D` five-estimator × two-velocity-solution run and is
promoted as an EXTENSION. It may not be promoted as an EXTENSION on the current evidence.

**SINGLE MISSING VERIFICATION TO REACH PANEL-READY.** ~~Merton's prior-art dossier~~ — delivered
(M-006). Remaining: the claim-text fix (mine, drafted above) plus the prior-art citation block,
after which this ships as a rediscovery. To ship as an EXTENSION instead, run the two velocity sets
through `Strain_2D` (Materna, Maurer & Sandoe 2021) to separate velocity-solution variance from
method variance — Merton's takeable #1, about a day of work.

---

## F-003 — A multi-angle, frozen-protocol null on tidal-phase *forecasting skill* in Southern California

**CLAIM (the permitted form — see the prose rule below).** Static tidal-phase susceptibility maps
have **no out-of-sample forecasting skill** in Southern California at M ≥ 1.5, 1981–2018, under a
protocol frozen and hashed before the test window was touched: train-period per-bin modulation on
0.4° bins selected exactly one bin (Anza, p_train = 0.048 of 42 bins tested); on the 2010–2018
test window that bin's frozen phase scored pooled S = −0.100 (anti-aligned), p = 0.94,
−0.052 bits/event, while the **anti-leak control on 22 train-null bins was clean**
(S = 0.009, p = 0.28). Phase migration does **not** replicate across catalogs (Anza 2012–2017:
QTM φ = +25.4°, n = 31 vs SCSN φ = +158.7°, n = 49 — same rock, same decade) and walk-forward
phase tracking has no skill (S = −0.095, p = 0.92 against a per-year-shift null with p95 = 0.129).
No credible bin-level susceptibility predictor exists at 0.4° scale: every apparent
feature-vs-amplitude correlation is the small-n amplitude bias (a_b vs n_train ρ = −0.49,
p = 0.001), and the bias-robust label correlates with nothing. Independently, the TSI ratio has no
positive spatial correlation with measured tidal modulation (Spearman ρ = −0.44 QTM, −0.03 SCSN),
and Long Valley — the strongest positive TSI anomaly in our geodetic map, pre-registered before
data download — shows no phase modulation in 19,787 declustered background events.

**MANDATORY PROSE RULE (Popper, R2-5, binding).** The sentence *"tidal triggering is null in
SoCal"* does not appear in any draft, email, briefing, README or panel material from this program.
The permitted form is the scoped sentence above. **The corpse of the map stands; the corpse of the
physics was never a corpse — it was an unresolved measurement described as one.** Any entry-derived
text I produce enforces this.

**CLASS — NEW APPLICATION (falsification).** A pre-registered, multi-angle falsification of a
prior claim of our own (Gale 2025) on an independent dataset, with clean positive/negative and
anti-leak controls, including a self-correction that killed our own most attractive intermediate
finding. The methodological asset — the split / walk-forward / anti-leak / per-year-shift harness
— is reusable and is arguably worth more than the null.

**ATTRIBUTION.** The claim being falsified is our own (Gale 2025, TSI; Kosmos r89 pipeline,
documented in `KOSMOS_METHODOLOGY.md`). Data and method credits: Lu, Xue, Yue, Zhuang & Zhao
(2025) for the openly shared SoCal catalogs and the §2.5 FM stress-resolution method; SCEDC and
the Yang–Hauksson–Shearer FM catalog; NCEDC/NCSS for Long Valley; NGL MIDAS. The self-correction
at §14b (overlapping-window autocorrelation dressed as coherence) is the program's own.
**PENDING-ATTRIBUTION, narrow:** Merton has not searched the tidal-triggering-forecast-skill
literature (Vidale et al. 1998; Cochran, Vidale & Tanaka 2004; Beeler & Lockner 2003; Tanaka;
Métivier; Ide; Scholz; van der Elst) to establish whether "static phase maps have no forecast
skill" has been shown before. My prior: the *physics* literature is enormous and the *forecast
skill* framing may well be less crowded — but I will not say "first" without his trail.

**PROOF CHAIN.**
- Protocols, all hashed in `download_log.md` before the corresponding analysis:
  `PROTOCOL.md` `76404c44df02665645a26f343114f478f8f8944b2bd03723f3c36e30d305c35d` (Amendment 1)
  — **verified by me** against the file as committed at `78f2227`;
  `XUE_LU_PROTOCOL.md` `3e36181a112b758ab7858ce244f2c5baf753fd63bf731d026e32a50a60011e09`
  — **verified** at `78f2227`;
  `protocol_params.json` `ef4d0a37e25ed67aadc245141836276a0b064f289416ee4da563d3ed9f4271cd`
  — **verified** against the working copy;
  `OVERNIGHT_PREDICTION_PROTOCOL.md`
  `62b259286d568a8ef0a446cf9749a8b1850f4ece0a3040fb22b8b02908755adf` — **verified** against the
  working copy.
  **Defect found: `LONG_VALLEY_PROTOCOL.md` does not verify.** `download_log.md` records
  `826ddc072275789e7f85276560f5ad9592033cdc6ea23bcb048e83a74920cffc` at freeze
  (2026-07-21 09:30:40 PDT); the file as committed at `78f2227` hashes
  `fd5b29789768df1f879a32f06570909eb1be29f4ab5a454c1d278af353433e50` and the working copy hashes
  `d3141cd19d638d6cf09cf9efaca0c451360472b0e1fef8607c3e0c70220cc2d4`. No amendment is logged.
  See audit finding #5 — the Long Valley sub-result must not be described as hash-verified until
  this is reconciled.
- Public commits: `78f2227` (tidal replication, blind cross-test, bin scan, TSI evaluation),
  `c8470e1` (overnight protocol freeze), `2f5e3e0` (overnight results, complete null suite with
  clean controls).
- Results: `results_xue_lu.json`, `results_coso_control.json`,
  `results_coso_fault_resolved.json`, `results_coso_fm.json`, `results_three_region.json`,
  `results_long_valley.json`, `results_bin_scan.json`, `binscan_{QTM,SCSN}.csv`,
  `results_exp_{a,b,c,c2,d,e}.json`, `results_tsi_map.json`.
- Controls: positive control (the strike-slip negative control reproduces a spurious "6.6×
  enrichment" by the same block-membership construction that produced the original 8.19×);
  anti-leak control (22 train-null bins, clean); cross-catalog detector-invariance arm (QTM vs
  SCSN — this is the arm that killed our own phase-migration finding); circular time-shift and
  per-year-shift nulls. Failed passes (a declustering method that collapsed, two parsing bugs) are
  logged in `download_log.md`, not hidden.
- Replication status: **VALIDATED** as a forecast-skill null within its stated scope; Popper
  R2-5 confirms "This stands."

**WHAT IT DOES NOT SHOW — and this is the half of the entry that matters most.**
1. It does **not** show that tidal triggering is absent in Southern California. It bounds nothing
   about the physical response amplitude. Popper's R2-1(b) correction is binding: the per-bin
   shortfall against a predicted ~1% effect is roughly **10×** in amplitude resolution (not the
   "three orders of magnitude" that circulated), and the *pooled* confirmatory statistic was
   within a factor of perhaps 2–4 of the predicted effect. **No one in this program quotes an
   order of magnitude for the shortfall until K-035 produces per-bin and pooled minimum-detectable
   amplitudes at 80% power, with the selection step simulated.**
2. It does not test conditional (state-dependent) triggering — never attempted.
3. Diurnal/S2 detection systematics were never modelled. At the N where this becomes a powered
   measurement it becomes a *systematics* problem, not a noise problem.
4. The bin scan does find significant modulation in the San Jacinto/Anza area independently in
   both catalogs, and elevated-but-underpowered amplitude at Coso. The null is about the *map's
   forecast skill*, not about the absence of modulation anywhere.
5. Long Valley: 19,787 declustered events, one caldera, one window — and see the hash defect
   above.

**TIER: PANEL-READY.** This is, essentially, the panel reply (`DRAFT_PANEL_EMAIL.*`,
`DRAFT_FIG4C_REPLY.*`). It ships with the prose rule and the "corpse of the map, not of the
physics" sentence attached.

**SINGLE MISSING VERIFICATION TO REACH PUBLISHABLE — DELIVERED, 2026-08-10. See F-012.** K-035 ran:
per-bin and pooled minimum detectable modulation at 80% power with the selection step simulated,
plus the R2-1(c) systematics arm. **The nulls are now upper bounds, and F-003 no longer stands
alone — it stands as the falsification half of a bounds paper.** Two consequences travel back into
this entry and are binding on any restatement of it:
1. **The R2-1(b) shortfall number is now measured, not estimated.** The best bound this design
   reaches is **6.3% at 80% power** (pooled over all 42 eligible bins, n = 3,920) against a
   Beeler–Lockner theory line of **1%** — a factor of **6.3**, not "three orders of magnitude" and
   not "10×". The full-catalogue intensity-likelihood configuration reaches **2.8%** (n = 23,465),
   a factor of 2.8 — **the closest this program has come to the theory line, and it still does not
   contact it.** `contacts_theory: false` on all six rows.
2. **A method limit was found that F-003's text did not know about.** At full-catalogue n the
   off-tidal negative-control line **fails** (reject rate 0.34 against a nominal 0.05), because ETAS
   clustering deposits power at the 11-day control line faster than the statistical error shrinks.
   It holds at n = 1,906 (0.06). **This bounds the METHOD, not the tidal result** — the tidal-band
   lines stayed calibrated in the a = 0 arm — but it means a powered full-catalogue tidal test must
   condition on an ETAS baseline, not a stationary λ₀. That sentence must appear in any paper built
   on this entry.

**REMAINING TO PUBLISHABLE (now a writing task plus one search).** Merton has still not searched
the tidal-triggering **forecast-skill** literature (the PENDING-ATTRIBUTION above), and the word
"first" may not appear until he has. Everything else is drafting.

---

## F-004 — Pure-math temporal-pattern nulls: periodicity, favored ratios, sequence-shape transfer

**CLAIM.** Under a protocol frozen and hashed before any test-window analysis
(`PATTERN_PROTOCOL.md`, 2026-08-09), on SCSN Southern California 1981–2018 with train < 2010 and a
single test scoring: (i) **no favored interevent-time ratio exists** — pooled log-ratio structure
in M ≥ 5.0 aftershock sequences (≥ 20 aftershocks) shows nothing at the golden ratio, train
p = 1.0 and test p = 1.0; the 41 apparent "candidate bands" are misspecification of a bare-Omori
(unbranched) null and 33 of 41 reproduce on test as a *shape* effect, i.e. burstiness beyond
unbranched Omori, not discrete preferred ratios. (ii) **No periodicity was confirmed**, including
the pre-declared live physical hypothesis of annual hydrologic loading; three multi-year train
detections (441 d, 719 d, 3,102 d) were killed by the frozen phase-agreement rule (~180° flips) and
read as slow rate fluctuation, not cycles. (iii) **Sequence shapes do not transfer spatially**:
predicting a test sequence's (p, b) from its 0.5°-neighborhood's train sequences beat the global
train mean in 7 of 15 cases, against a frozen success rule of ≥ 60%.

**CLASS — NEW APPLICATION (falsification).** Pre-registered, single-scoring, honest-prior-declared
nulls on a family of claims that recur constantly in the informal literature.

**ATTRIBUTION — PENDING-ATTRIBUTION.** What Merton must search: earthquake periodicity /
Schuster-test literature (Ader & Avouac; Beaucé; the seasonal-modulation line —
Sirorattanakul & Avouac 2026 is already in the ledger and reports ~15% annual modulation in
northern California with minimal SoCal signal, which is adjacent and helpful); prior tests of
"golden ratio"/log-periodic structure in seismicity (Sornette's log-periodic corrections are the
obvious neighbor and must be distinguished from what we tested); and the ETAS
parameter-transferability literature for (iii).

**PROOF CHAIN.**
- Protocol: `PATTERN_PROTOCOL.md` at round-2 freeze, SHA-256
  `165527d14b28bd1a0ea1cf5340e6b0252548d0e5ea8c2d33e6f1edca0ad2aa16`, recorded in
  `download_log.md`. **I verified this hash against the file as committed at `675c095`.** Clean.
- Public commits: `675c095` (protocol freeze), `4719739` (results).
- Results: `results_exp_f.json`, `results_exp_g.json`, `results_exp_i.json`.
- Controls: Omori-simulation nulls per sequence (EXP-G); inhomogeneous-Poisson surrogates with
  σ = 5P kernel smoothing (EXP-F); pre-labeled artifact periods (0.5 d, 1 d, 7 d) declared in
  advance and reported separately.
- Replication status: EXP-G **VALIDATED at high power** (Popper R2-5: "Stands at high power",
  no material exposure). EXP-F **provisional and weak** — see below. EXP-I(i) validated.

**WHAT IT DOES NOT SHOW.** **The EXP-F periodicity null is weak evidence and must never be quoted
without its failed positive control attached: the 7-day method check did not fire.** A comb that
cannot detect the artifact it was built to detect cannot be trusted to report the absence of a
real period. Popper's corpse-exposure line: comb power unquantified, input spectrum never divided
out, no injection–recovery. EXP-G's null does not exclude continuous scaling structure — only
discrete favored ratios. EXP-I(i) does not show that sequence shape is unpredictable, only that
0.5°-neighborhood pooling does not beat the global mean.

**TIER: INTERNAL** — held here entirely by EXP-F's failed positive control and the missing
attribution. EXP-G alone would be PANEL-READY.

**SINGLE MISSING VERIFICATION TO REACH PANEL-READY.** An injection–recovery curve for the
periodicity comb (inject a known-amplitude period, measure recovery) plus a repaired 7-day
positive control. Either the method fires and the null becomes quotable with a detectable-amplitude
bound, or it does not and EXP-F is withdrawn from the register entirely. I would rather withdraw it
than ship it.

---

## F-005 — B-2: SoCal walk-forward temporal ETAS skill against a rate-oracle Poisson

**CLAIM.** On SCSN Southern California (lat 31.5–38.0, lon −122.0 to −113.5), M ≥ 2.5, with ETAS
parameters fit by MLE on 1982-01-02 → 2010-01-01 and **frozen** (μ = 0.27502/d, K = 0.041238,
α = 0.53656, c = 0.014264 d, p = 1.11829, M0 = 2.5), walk-forward daily-rate scoring on the
2010-01-01 → 2018-12-31 test window (n_test = 8,722 events) yields **+1.907 bits/event** against a
stationary Poisson at the train-period rate and **+1.866 bits/event** against a Poisson at the
*test-period* rate (an oracle baseline that removes any "the 2010s were quieter" trivial skill).
The frozen success rule was > 0.5 bits/event against both; it passed by 3.7×. Completeness was
checked per decade and per split before the floor was chosen (worst-case Mc + 0.2 = 1.7,
M0 = 2.5 stable). Mean triggered fraction over test events = 0.915.

**Numbers re-verified by me against `results_exp_h.json` this session: 1.9072792871648883 and
1.8656787636987344. Exact.**

**SUB-CLAIM HELD BACK (not claimable).** "Skill rises with magnitude" (M ≥ 4 band:
+2.578 / +2.536 bits/event with the uniform integral correction) rests on **n = 290** events that
are not independent draws — M ≥ 4 events in SoCal arrive in a handful of sequences. Popper S-5:
this requires a **block bootstrap over sequences, not events**, before it is quoted with a CI. It
is therefore listed here and excluded from the claim sentence.

**CLASS — REPRODUCTION (presumed).** Temporal ETAS beating Poisson out-of-sample in California is
Ogata's result and CSEP's daily bread. Our contribution is a clean, frozen-parameter,
harsh-baseline instance of it that establishes our machinery works — and it is the backbone every
other claim in this program is scored against. That is exactly what a reproduction is for.

**ATTRIBUTION — RESOLVED: REDISCOVERY (canonical). Merton dossier M-004 (`HYPOTHESIS_LEDGER.md`).**
Zero novelty in the finding — Ogata (1988, *JASA* 83, 9–27) is the generator and this is the
premise, not a result, of every CSEP short-term experiment. **The question Merton was asked is
answered: +1.907 bits/event = a probability gain of 3.75× per event (1.322 nats; 0.574 log10);
+1.866 = 3.65×. The published California comparators are ≈6× (Werner, Helmstetter, Jackson &
Kagan 2011, *BSSA* 101(4), 1630–1648: "a probability gain per earthquake of about 6" = 2.585
bits/event, all-CA, m≥3.95) and >10× (Helmstetter, Kagan & Jackson 2006, *BSSA* 96(1), 90–106,
SoCal, m≥2 = >3.32 bits/event). Ours sits BELOW both — and must, because theirs are
spatio-temporal gains over a spatially smoothed time-independent model while ours is temporal-only
over a stationary rate with no spatial term. Verdict: ORDINARY, on the conservative side. Not
high, not suspicious.** Two mandates follow. (i) **Units discipline:** the field reports this in
four incompatible conventions (probability gain factor; IGPE in nats; Kagan score in bits; per-bin
log-likelihood). "1.907 bits" must never be set beside a RELM "0.3" — that is a category error in
both units and forecast class. (ii) **The temporal-only qualifier is load-bearing in the opposite
direction from what we assumed:** without it a seismologist reads our number as spatio-temporal and
judges it implausibly *low*. Merton also notes the published gain is not a constant of nature —
Werner et al. document it moving by a factor of two with region size, grid size and target
threshold — so every bits/event quotation in this program must carry those four qualifiers.

**PROOF CHAIN.**
- Protocol: `PATTERN_PROTOCOL.md` §EXP-H, SHA-256
  `165527d14b28bd1a0ea1cf5340e6b0252548d0e5ea8c2d33e6f1edca0ad2aa16` — **verified at `675c095`**.
- Public commit: `4719739`. Result: `results_exp_h.json`. Code: `exp_h_etas.py`.
- Controls: two baselines, the harsher of which (test-rate oracle) is the one quoted;
  four MLE starts, all converging to the same optimum (max LL spread < 1e-5);
  truncation error measured, not assumed (ΔLL/event = −0.0093 for the 1,000-day fit truncation;
  test scoring is untruncated).
- Recorded flags, carried forward honestly: the final untruncated polish **hit its iteration cap**
  (nit = 25 = maxiter) and the frozen parameters come from that polish; α = 0.537 sits near its
  bound (expected for temporal-only, and Hainzl et al. 2013 note short-term aftershock
  incompleteness biases α *downward* — so this value is not to be interpreted physically);
  nominal branching ratio n = 1.161 is **supercritical**, while the 7-day effective n = 0.60 —
  i.e. the parameterisation is window-dependent and the nominal n must not be quoted as a physical
  criticality statement.
- Replication status: **VALIDATED** within scope (one region, one catalog, one nine-year window).

**WHAT IT DOES NOT SHOW.** Temporal only — no spatial forecast is claimed and no spatial baseline
exists (Popper G4: no spatial entry may quote bits until K-002 exists). Not detector-invariant
(K-028 untested — SCSN vs QTM vs ComCat magnitude-scale consistency is **unverified**). Post-2018
untested. The M ≥ 4 rise is held back as above. It does not show ETAS is *correct*, only that it
is much better than Poisson: K-009 (F-009) shows its residuals are not white.

**MERTON M-009.12(a) — THE CORRECTION THAT SUPERSEDES M-004.3, AND IT GOES AGAINST US.** M-004.3
placed +1.907 below every published California comparator and defended that as *expected, because
ours is temporal-only and theirs are spatio-temporal*. **That defence is withdrawn.** Merton
obtained the fitted per-event log-likelihoods from Stockman, Lawson & Werner's own public repository
(`github.com/ss15859/EarthquakeNPP`, `Experiments/ETAS/output_data_SCEDC_25/ll_scores.json`,
retrieved 2026-08-10) — **temporal-only, same metric decomposition, four California datasets:**

| dataset | region / floor | temporal gain, bits/event | prob. gain |
|---|---|---|---|
| `ComCat_25` | all California, M ≥ 2.5 | 1.330 | 2.51× |
| **B-2 (ours)** | **S. California, M ≥ 2.5, temporal-only** | **+1.907** | **3.75×** |
| `SCEDC_20` | S. California, Mw ≥ 2.0 | 2.738 | 6.67× |
| **`SCEDC_25`** | **S. California, Mw ≥ 2.5 — the matched comparator** | **3.319** | **9.98×** |
| `SCEDC_30` | S. California, Mw ≥ 3.0 | 4.191 | 18.26× |

**THE PERMITTED SENTENCE, AND IT IS THE ONLY ONE.** *"Our temporal-only walk-forward ETAS scores
+1.907 bits/event against a test-rate-oracle Poisson on SCSN SoCal M ≥ 2.5, 2010–2018 — inside the
published temporal-only range of 1.33–4.19 bits/event for California, above the statewide
`ComCat_25` comparator and at **57% of the closest matched comparator** (`SCEDC_25`, 3.319
bits/event, Mw ≥ 2.5 SoCal, test 2014–2020)."* **The words "typical", "ordinary", "on the
conservative side" and any framing that reads as reassurance are retired from this entry.** Two
innocent explanations exist and neither is verified: test-window composition (theirs 2014–2020,
essentially Ridgecrest; ours 2010–2018, El Mayor's decay plus quiet stretches) and catalog clipping
(theirs is a shape-filed SCEDC product, ours is not identically clipped). **Until one of those is
measured, "we are 1.4 bits below the matched published number" is the fact and the explanations are
hypotheses.** *Provenance flag carried from Merton: these are the authors' committed model outputs,
not values read off their Figure 2, and no sentence of the form "Figure 2 shows X" may be written
on this basis.*

**TIER: INTERNAL.** *Downgraded from "PANEL-READY on my next pass" by M-009.12(a).* The attribution
block is cleared, but the entry can no longer be promoted on an editorial pass, because the
literature placement it would carry into a room with Werner in it has changed sign in tone.

**SINGLE MISSING VERIFICATION TO REACH PANEL-READY — now mandatory, not optional.** Score B-2 inside
the `SCEDC_25` protocol of Stockman, Lawson & Werner (2026, *TMLR*, arXiv:2410.08226): same SCEDC
catalog, same Mw ≥ 2.5, same train 1985→2014 / test 2014→2020 split, same temporal/spatial
likelihood decomposition, published ETAS and homogeneous-Poisson baselines, public code. Merton has
upgraded this from *desirable* to **required before B-2 is quoted against the literature again**,
and I adopt that. It is the difference between "we are below the published number" and "we are below
the published number and we know why."

---

## F-006 — B-1: generic temporal ETAS transfers to never-trained regions

**CLAIM (Popper's mandated quotation form, S-5 / §10).** Temporal ETAS parameters pooled as a
single **global** median across seven training regions (Japan, Chile, Indonesia, California,
Turkey, Himalaya, Iceland; ComCat M ≥ 4.5, 1995–2026, burn-in 1995–2000) and applied without
refitting to **six regions never used in fitting**, scored walk-forward 2000–2026 against each
region's own-period Poisson rate (a deliberately harsh *local oracle* baseline), yields
**+0.66 to +0.84 bits/event on the four-to-five adequately-powered holdouts** — Alaska-Aleutians
+0.843 (n = 3,061), Philippines +0.780 (n = 7,429), Greece-Aegean +0.792 (n = 1,765), Iran +0.687
(n = 1,521), Mexico +0.660 (n = 2,002) — with **Caribbean +1.753 flagged `underpowered: true`
at n_scored = 235** and quoted only with that flag. **Fault-type parameter pooling FAILED its
frozen sign test**: the TYPE pool beat the GLOBAL pool in 2 of 6 holdouts against a required 5 of
6. The informative shape of the failure: TYPE pools differed from GLOBAL mainly through μ, a
per-region background *rate* that should never be pooled (the subduction pool's μ = 0.310/d,
Chile/Indonesia-dominated, wrecked Alaska and Mexico).

**CANDIDATE EXTENSION (not yet claimable).** "The transferable object is the universal clustering
law plus a locally estimated μ" is a quantified six-holdout demonstration of a proposition the
field broadly assumes. Whether the *quantification* is a delta over known ETAS universality is
**Merton's to rule**.

**ATTRIBUTION — RESOLVED, SPLIT. Merton dossier M-005 (`HYPOTHESIS_LEDGER.md`).**
*The positive is a* **REDISCOVERY**; *the negative is* **CONTESTED, leaning CONTRADICTED on the
adjacent quantity.** Page et al. (2016) does **not** own both halves — **it owns the opposite of
one of them**, and that had to be found before a reviewer found it.

- **Positive (generic parameters transfer) — REDISCOVERY, and the owning paper was missing from
  this entry: Bayona, Savran, Iturrieta, Gerstenberger, Graham, Marzocchi, Schorlemmer & Werner
  (2023), *The Seismic Record* 3(2), 86–95, doi:10.1785/0320230006**, "Are regionally calibrated
  seismicity models more informative than global models?" — global GEAR1 vs **19** regional models,
  **prospectively**, 2014–2021, CSEP metrics, three regions: GEAR1 ranks 1st in New Zealand, 2nd in
  California, 3rd in Italy. Same thesis, larger scale, already in print (time-independent M≥4.95
  class, not temporal ETAS — which is our niche). Add **Chu, Schoenberg, Bird, Jackson & Kagan
  (2011), *BSSA* 101(5), 2323–2339**: *"the ETAS model with few parameters and with the same
  functional form seems to fit reasonably well to the seismicity in each zone."*
- **Negative (fault-type pooling adds little) — CONTESTED.** Page, van der Elst, Hardebeck,
  Felzer & Michael (2016), *BSSA* 106(5), 2290–2301, doi:10.1785/0120160073 make **tectonic region
  their first ingredient** and find **regional variation in mean aftershock productivity of almost
  a factor of 10**; Chu et al. (2011) find productivity varying **>5×** across plate-boundary
  zones; Hardebeck, Llenos, Michael, Page & van der Elst (2019), *SRL* 90(1), 262–270 find
  productivity varying *within* California (SoCal > NoCal, Mendocino low, Long Valley/Coso/Salton
  Sea high). Read naively, we are contradicted.
- **Merton's reconciliation, which must be written into the claim before promotion.** They
  regionalize **productivity**; we pooled a **background rate**. Chu et al. measured that same
  quantity and found background rates **"range by a factor of nearly 500"** — a hundredfold larger
  spread than the productivity spread in the same study. Our TYPE-pool failure is therefore a
  rediscovery of Chu et al.'s 500× μ spread arriving through a forecast score. Further, the
  Reasenberg–Jones/Page framework **has no μ to pool at all** (it forecasts aftershocks of an
  identified mainshock), so the failure mode that killed our TYPE pools is structurally impossible
  in the model where "tectonic region helps" was demonstrated. The two results are not in contact.
- **Permitted claim shrinks accordingly.** "Regional/type-specific tuning adds little" is **not
  supportable** and would be refuted by three papers in one paragraph. The supportable sentence is:
  *"pooling the full temporal-ETAS parameter vector by fault type — including μ — does not transfer
  better than a single global pool, because μ is not a poolable quantity."* Our CANDIDATE EXTENSION
  line below is the right sentence and is now attributable: it is Chu et al. (2011)'s conclusion,
  quantified by a forecast score in six never-trained regions instead of by a parameter table in
  eight tectonic zones. Real delta, small delta.

**PROOF CHAIN.**
- Protocol: `PATTERN_PROTOCOL.md` §EXP-M addendum, SHA-256
  `aca4b729277762fe1ca9f9fdf561291e3527710f5240394b7a7a12e06d6995b2`, recorded in `download_log.md`
  as frozen **before any global data download**. **I verified this hash against both the working
  copy and the file as committed at `a45ca8d`.** Clean, and the four-step PATTERN_PROTOCOL hash
  chain (`165527d1` → `1e126abc` → `4b347599` → `aca4b729`) verifies at commits `675c095`,
  `a40fc9a`, `0ba8c5e`, `a45ca8d` respectively. This is the cleanest freeze record in the program.
- Public commits: `0f4b7ec` (protocol freeze, pre-download), `a45ca8d` (results).
- Result: `results_exp_m.json`; retrieval log `data/comcat_world_log.txt`. Code:
  `exp_m_world_transfer.py`.
- Controls: local-oracle Poisson per region (harsh by design); each holdout's own post-hoc
  in-sample fit reported explicitly as a *ceiling, not a forecast*; SoCal-parameter column
  reported and caveated as a construction artifact (M2.5-fitted sources applied to an M4.5
  catalog, −7 to −10 bits).
- Worker integrity catches, on the record: a **sign error in the supervisor's K-rescaling brief**
  (K′ = K·10^(+α(M0_new − M0_old)), not minus) was caught by the executing worker, both derivations
  documented, the correct one used; one spurious sub-floor ComCat row in Indonesia dropped and
  logged.
- Replication status: **VALIDATED** (six never-trained regions is itself the replication).

**WHAT IT DOES NOT SHOW.** Temporal only, M ≥ 4.5, ComCat only — no detector-invariance test
across catalogs (K-028). **The universality of p is narrower than the program has been saying:**
of the 13 regional fits on record (7 train + 6 holdout post-hoc), p lies in [0.94, 1.08] in
**11**, not 12 or 13 — Iceland (train, p = 1.191) and Caribbean (holdout own-fit, p = 1.208) sit
outside. See audit finding #2. It does not show μ cannot be regionalized by some *other* covariate
than fault type. It does not show TYPE pooling fails at finer magnitude or spatial granularity —
only at temporal-ETAS/M4.5 granularity, with n = 1 pools for collision and rift, which is
acknowledged in the protocol and is a real weakness of the frozen design.

**TIER: INTERNAL.** ~~Blocked on attribution~~ — attribution cleared (M-005), but promotion is
**not** automatic: M-005 raises the bar rather than lowering it. Two things now stand between this
entry and PANEL-READY: (a) the negative must be restated in the narrowed form above (mine,
editorial, no compute), and (b) Merton requires a **positive control** before the word "little" is
used about regional tuning anywhere — see below.

**SINGLE MISSING VERIFICATION TO REACH PANEL-READY.** ~~Merton's B-1 dossier~~ — delivered (M-005),
and it answers the Page-et-al. question in the negative: the published global-parameter transfer
literature owns our positive (Bayona et al. 2023; Chu et al. 2011) but reports the **opposite** of
our negative on the neighbouring quantity (productivity). Remaining verification: **the
Hardebeck et al. (2019) within-California positive control** — reproduce the published SoCal /
NoCal / Mendocino / hydrothermal productivity differences with our pipeline. If we cannot recover a
regional difference that is known to exist, our fault-type null is uninformative rather than
negative. Merton's other takeables (adopt García et al. 2012 or Bird 2003 regionalization instead
of our home-made fault types, whose n = 1 collision and rift pools are an artifact of the taxonomy;
pool (K, α, c, p) by type and estimate μ locally *by construction*) apply to the next transfer
experiment, not to this one.

---

## F-007 — B-3: the seven-day M≥5 conditional rule is a local re-measurement of Reasenberg–Jones

**CLAIM (with its definitional footnote, which is not optional).** In the Southern California box,
with the window length X* = 7 days **frozen from the training period** (train: 37/99 = 0.374),
the probability that **any** M ≥ 5.0 event — *including equal-or-smaller aftershocks of the same
sequence* — occurs within 7 days of an M ≥ 5.0 event is **18/30 = 0.60** on the 2010–2018 test
window, against a Poisson base rate of 0.062 (binomial p = 7.4 × 10⁻¹⁵). This is **ordinary
within-sequence triggering, not exotic meta-structure**, and it is a local re-measurement of an
effect that has been the operational basis of California earthquake advisories since 1989.

**THE FOOTNOTE THAT MUST TRAVEL WITH THE NUMBER (Merton, M-003, mandatory).** Jones (1985) reports
P(*a larger* event within 5 d and 10 km | M ≥ 3 SoCal event) = 6 ± 0.5%, rising to 6.5 ± 2.5% at
M ≥ 5. **Ours is P(any M ≥ 5, same sequence included) — a different statistic.** The numbers are
not in conflict, but they *will be read* as in conflict if 0.60 is quoted without the definition.
Anyone reporting B-3 states "any M ≥ 5, same sequence included" in the same sentence as the number.

**CLASS — REPRODUCTION.** Merton: "REDISCOVERY — canonical, textbook, operational. Zero novelty.
Cite or retract." I adopt that verbatim. Delta: none — not scale, not rigour, not region. Our
n_test = 30 is far smaller than any of the prior work, and the p = 7 × 10⁻¹⁵ measures how strongly
we rejected a Poisson straw man that nobody in this literature has believed since 1989. Its
legitimate value is (a) machinery validation and (b) an app layer.

**ATTRIBUTION (complete — Merton M-001…M-003, dossier M-003).**
1. **Reasenberg & Jones (1989), *Science* 243, 1173–1176** — "Earthquake hazard after a mainshock
   in California": the generic California clustering model (modified Omori + Gutenberg–Richter)
   giving the probability of further events, including larger ones, in intervals after any
   earthquake. **This is B-3.** Extended in Reasenberg & Jones (1994), *Science* 265, 1251.
2. **Jones (1985), BSSA 75, 1669–1680** — foreshocks and time-dependent hazard in southern
   California; the 6–6.5% P(larger) numbers above.
3. **Gerstenberger, Wiemer, Jones & Reasenberg (2005), *Nature* 435, 328–331** — STEP, the
   operational short-term clustering forecast; direct ancestor of every 7-day hazard window.
4. **Ogata (1988), JASA 83, 9–27** — ETAS. B-3 is a marginal read-out of the same generator that
   produces F-005/B-2; **B-2 and B-3 must never be quoted as independent findings.**
5. **Michael et al. (2020), SRL 91, 153–173** (USGS operational aftershock forecast);
   **Page, van der Elst, Hardebeck, Felzer & Michael (2016), BSSA 106, 2290–2301**; Hardebeck
   et al., updated California aftershock parameters (USGS) — current practice, with published
   parameter values for the exact windows we used.

**PROOF CHAIN.**
- Protocol: `PATTERN_PROTOCOL.md` §EXP-I(ii), SHA-256
  `165527d14b28bd1a0ea1cf5340e6b0252548d0e5ea8c2d33e6f1edca0ad2aa16` — **verified at `675c095`**.
  X* was mined on train and frozen before the single test scoring, exactly as specified.
- Public commit: `4719739`. Result: `results_exp_i.json` §ii. Code: `exp_i_meta.py`.
- Controls: 1,000 Poisson simulations for the train-side mainshock-clustering diagnostics
  (interevent CV = 1.473 vs null mean 0.996, p = 0.001; Ripley-K in time exceeds Poisson at 30 d,
  100 d and 365 d, p = 0.001 each). Frozen X* selection from a four-value grid {7, 30, 90, 365 d}.
- Replication status: **VALIDATED** in scope; and independently corroborated by 37 years of
  operational practice, which is the strongest replication status any entry here has.

**WHAT IT DOES NOT SHOW.** n_test = 30 eligible mainshocks, and they are not independent — a
single productive sequence contributes several. The test value (0.60) is substantially above the
train value (0.374); with n = 30 that gap is not evidence of anything and must not be presented as
a trend. Single region, single catalog, no detector-invariance test (K-028). It says nothing about
the probability of a *larger* event, which is the quantity operational forecasting actually cares
about.

**TIER: PANEL-READY** — as a reproduction, with the citation block and the definitional footnote
attached. It is honest, it is useful as an app layer (EQ-15 layer 2), and presenting it as anything
more would be a retraction waiting to happen.

**SINGLE MISSING VERIFICATION TO REACH PUBLISHABLE.** Recompute the *same* statistic in the
Reasenberg–Jones parameterisation — P(larger event | M ≥ 5), 7 d — alongside ours, so the two
numbers sit in one table and the comparison to Jones (1985) is made by us rather than by a
reviewer. One afternoon on data already on disk. (Realistically this entry's ceiling is
PANEL-READY: with zero delta there is nothing here to publish, and saying so is the point.)

---

## F-008 — B-4: the stress ledger's negative space recovers known aseismic geology blind

**CLAIM.** A 0.2°-cell ledger over Southern California comparing geodetic loading
(Ṁ_geo = 2μHA·max-shear from NGL MIDAS, μ = 30 GPa, H = 11 km flat / a frozen variable-H model)
against catalog seismic release (Hanks–Kanamori sum, SCSN M ≥ 2.5, train 1981–2009), computed
**without any fault or creep information as input**, classifies 228 of 1,195 cells as
SILENT-LOADING (top-quartile loading with χ = Ṁ_seis/Ṁ_geo < 0.01 or fewer than 20 train events).
175 of 229 silent cells are grid-**interior** (18.5% silent among interior vs 20.9% among edge
cells — so not a border artifact), and the top interior silent cells are the **San Andreas
creeping section** (36.1–36.9 N, −120.8 to −121.8) and **Imperial/Brawley**. After excluding cells
within 25 km of the creeping segment and within 0.35° of a geothermal field, the 200 remaining
"unexplained-silent" cells are topped by cells lying **1–5 km from the 1857 Fort Tejon
Mojave/Big Bend rupture strand** (34.4–34.8 N, −117.5 to −119.1), nearest-fault
`SAFS-SAFZ-1857-San_Andreas_fault_rupture-CFM5`. **The correct reading is instrument validation:
the ledger's negative space rediscovered known aseismic and known late-cycle geology without being
told where they were.** The silent set is robust to the seismogenic-thickness model (Jaccard 0.949
between flat-11 km and variable-H); style does not significantly stratify coupling
(Kruskal–Wallis p = 0.34, power-limited); and χ does **not** persist train → test
(Spearman ρ = −0.131, n = 64, p = 0.30).

**CLASS — REDISCOVERY (of the geography) + a small methods delta.** Merton, M-002, and the program
record already said so (EQ18 §15: "the ledger's negative space rediscovers known aseismic zones
blind = **validation**"). That sentence is the correct classification and **it must never drift
into a discovery claim.** Our delta is methodological and modest: a 0.2°-resolution
catalog-Kostrov-vs-geodesy ledger, trained pre-2010, whose *negative space* was read as the signal,
plus the interior/edge check, the style stratification, the H-sensitivity Jaccard, and the
explicitly published unexplained-silent list.

**ATTRIBUTION (complete — Merton M-002).**
1. **Meade & Hager (2005), JGR 110, B03403 — "Spatial localization of moment deficits in southern
   California."** GPS-constrained block model: SoCal scalar moment accumulation
   17.8 ± 1.1 × 10¹⁸ N m/yr, ≈50% larger than the 200-year average release rate, with deficits
   localized in three regions — the southern SAF/San Jacinto, the offshore faults + LA/Ventura
   basins, and the Eastern California Shear Zone. **This is B-4's experiment, twenty-one years
   earlier, with a better loading model.**
2. **Kostrov (1974), Izv. Earth Phys. 1, 23–44** — the moment summation the ledger rests on;
   **Ward (1994, 1998, GJI 134, 172)** for California Kostrov budgets.
3. **Guns, Bennett et al. (2024), JGR 129, e2023JB027939** — "Seismic Moment Accumulation Rate
   From Geodesy: Constraining Kostrov Thickness in Southern California" — the current state of the
   art on exactly our quantity, and where our χ *level* should be calibrated.
4. **Field et al. (2014), BSSA 104, 1122–1180 (UCERF3)**; WGCEP (1988, 1995); Weldon et al. (2004)
   — the 1857 Mojave/Big Bend section is the most-forecast locked, late-in-cycle strand in
   California. "The ledger found the 1857 strand blind" is a statement about the **ledger's
   sensitivity**, not about the Earth.
5. **Jolivet et al. (2015), GRL 42, 297**; **Ryder & Bürgmann (2008), GJI 175, 837**;
   **Tong, Sandwell & Smith-Konter (2013), JGR** — InSAR coupling maps of the creeping section at
   far higher resolution than a 0.2° catalog-derived ledger, from data that is not degenerate the
   way ours is.
6. **Liu, Ross, Cochran & Lapusta (2022), *Science Advances* 8, eabk1167** — creep rate along the
   central SAF is directly proportional to the fraction of non-clustered earthquakes, 1984–2020.

**W-006-P1(a) PRE-FALSIFICATION NOTE (Merton, mandatory).** W-006-P1(a) predicts that *no
catalogue-derived statistic separates creeping from locked silent cells at better than chance.*
**Liu et al. (2022) already published a catalogue-derived statistic that does exactly that** on
the central SAF. W-006-P1(a) as written is close to pre-falsified in the literature and must be
re-scoped (to SoCal cells off the central SAF) or reframed as "the *ledger's own* χ does not
separate them, but a clustering statistic does" — which is a better claim and hands B-4 a free
upgrade rather than killing it.

**PROOF CHAIN.**
- Protocols: `PATTERN_PROTOCOL.md` §EXP-J addendum, SHA-256
  `1e126abc8a42a1e6ca7e0a5b9a874d83a9084ff9b5496c2bb5d57047ed119a6d`, **verified at `a40fc9a`**,
  frozen before computation; §EXP-K/L addendum
  `4b347599113aca2dd7ae6313c178f2e142d3fc31e0632e8e82b81a33fa581e54`, **verified at `0ba8c5e`**.
- Public commits: `a40fc9a` (EXP-J freeze), `4719739` (EXP-J results), `0ba8c5e` (EXP-K/L freeze),
  `ca54d56` (EXP-K results). Results: `results_exp_j.json`, `results_exp_k.json`. Maps:
  `maps/exp_j_ledger.png`, `maps/exp_k_stratified.png`.
- Controls: interior-vs-edge border check (the one that matters, and it passes);
  H-model sensitivity (Jaccard 0.949 — H moves the *ranking*, not the *membership*);
  CFM5.3 geometry cross-check (557 segments parsed, 0 bad headers, CFM header strike agrees with
  local trace strike to within 20° for 88.6% of cells); strain-grid cross-check against EXP-J
  (max abs difference 2.8 × 10⁻¹⁴ nstrain/yr — i.e. bit-identical).
- Replication status: **VALIDATED as instrument validation**, and **AUTO-FLAGGED FOR CHALLENGE**
  under the retirement mechanism trigger 4(b).

**WHAT IT DOES NOT SHOW — and this entry carries the heaviest load of caveats in the register.**
- **158 of the 200 unexplained-silent cells have n_train < 20; 39 have n_train = 0; 95 have
  n_test = 0.** Four-fifths of the hazard-candidate list is detection-limited, not measured-silent.
  Only 42 are measured-low-χ.
- χ does not persist train → test (ρ = −0.131, n = 64). **χ at 0.2°/29 yr is transient-dominated
  and must not be described as a material constant.**
- χ absolute values are biased low (M ≥ 2.5 over 29 yr misses rare large-event moment). **The
  geography is the signal; the level is not.** Guns et al. (2024) is where the level lives.
- Cell (32.8, −115.1), n = 0, is a flagged border-completeness artifact.
- Station-density confound untested (K-031, the open challenge); detector-invariance untested
  (K-028).
- **A live over-quote sits in the primary artifact:** `results_exp_j.json` carries
  `"verdict": "CATCH-UP wins (prior overturned)"`. The supervisor's own §15 note says the catch-up
  headline is mostly coupled-cell aftershock-era decay plus the low-n rate regularizer. See audit
  finding #3. **No catch-up claim is made in this register.**
- It is not a hazard map. Wegener's reading — that B-4 is better described as *a direct detection
  of the aseismic slip field through its negative space* — is more defensible than the hazard
  framing and is the direction I expect this entry to move.

**TIER: INTERNAL.** It is the program's strongest current claim and it is simultaneously
auto-flagged for challenge; those are consistent, and publishing before the challenge resolves is
exactly the mistake this register exists to prevent.

**SINGLE MISSING VERIFICATION TO REACH PANEL-READY.** W-006-P1(b): break the degeneracy with
geodesy. Test whether an independent geodetic/InSAR coupling field (Jolivet 2015 /
Ryder & Bürgmann 2008 — published, no new inversion) separates creeping from locked silent cells
where the catalog-derived χ does not. 100% on disk, one afternoon. If it separates them, B-4
narrows *and* hardens into an aseismic-field detection; if nothing separates them, the silent list
is a low-count artifact and we say so. (Free companion, and a reviewer will ask for it anyway:
Jaccard our silent list against Meade & Hager's three localizations.)

---

## F-009 — K-009: the fitted-ETAS residual field in Southern California is not white

**CLAIM.** On a 0.2° × 7-day grid over Southern California, 2010–2018, M ≥ 2.5 (448 cells,
469 weeks, 8,720 events), with temporal ETAS parameters frozen from the pre-2010 training window
and the spatial kernel fit on train only, the residual field departs from the field simulated by
the *same fitted generator*: pooled **lag-1 weekly autocorrelation = 0.0958** against a
sim-null median of −0.0165 and 97.5th percentile of 0.0023 (excess **+0.0935**);
**Moran's I = 0.0114** against a null 97.5th percentile of 0.00045; **leading EOF variance
fraction = 0.197** against a null 97.5th percentile of 0.0508; **integral correlation time
= 2.42 weeks** against a null 97.5th percentile of 0.0035 weeks. The excess is larger in the grid
interior than at the edges (ACF1 0.103 vs 0.049), is stable across three background-kernel
variants (0.076–0.096), and survives an Mc-completeness-proxy partial (0.0885). **Popper's frozen
success rule — lag-1 ACF excess ≥ 0.05 over the sim-null 97.5th percentile, stable across the
kernel swap, surviving the ρ_sta partial — is met.** The two-generator verdict registered before
scoring: **W-003 ("the residuals are white") is falsified; W-001/W-002 ("red at the local t_a")
is consistent but NOT sharply confirmed** — the measured correlation time is a censored lower
bound (> 52 weeks on the crossing estimator) that *overlaps* rather than lands inside the
factor-2 acceptance band of 13.8–55.3 weeks around the independently measured t_a = 193 d.

**WHAT MAY NOT BE SAID (Merton M-001.0, and I am reinforcing it).** The prose headline
"months-scale temporal persistence with ~40 km spatial coherence" **is not what these numbers
say.** The 52-week crossing correlation time is the *estimator's ceiling* and **all 20 null sims
also sit on it** — the statistic has no resolving power there. The 41.0 km e-folding length is
*smaller* than the null's median (78.7 km) and sits at the null's 2.5th percentile — that is not
an excess. The statistics that actually separate real from null are ACF1, the integral correlation
time (**2.4 weeks, not months**), Moran's I, and the EOF1 variance fraction. Those are the only
ones quoted above, and they are the only ones this program may quote.

**AND ONE MORE, WHICH I FOUND IN THE PRIMARY JSON AND WHICH IS NOT IN MERTON'S TABLE.** Excluding
the 52 weeks of 2010 — the year of the M7.2 El Mayor-Cucapah sequence — the temporal excess largely
collapses: **ACF1 falls 0.0958 → 0.0382** (below Popper's own ≥ 0.05 excess bar) and the
**EOF1 variance fraction falls 0.197 → 0.0494, i.e. below the null's 97.5th percentile of
0.0508 — the EOF excess disappears entirely.** Moran's I is comparatively robust
(0.0114 → 0.0099). The leading residual EOF's top five weeks carry 90.6% of its variance and it
peaks in the week of 2010-04-02. **So the defensible statement is: the *spatial* coherence is
robust; the *temporal* persistence is substantially carried by one great sequence.** The commit
message for `c371442` says exactly this and is correct; Merton's M-001.0 table does not carry it,
and anything built from that table alone will over-quote. See audit finding #1 — this is the
most material discrepancy in the register.

**CLASS — REPRODUCTION (Merton: REDISCOVERY, strong, multiply independent, in our own region) with
exactly one candidate methods delta.** "The residual/background field of a fitted ETAS is not white
in Southern California" has been established at least three independent times, by three
communities, with three estimators, before us. **It is not NOVEL and must never be written as
such.**

**ATTRIBUTION (complete — Merton M-001; the five load-bearing, then the rest).**
1. **Zaliapin & Ben-Zion (2020), JGR Solid Earth 125, e2018JB017120.** Nearest-neighbour
   declustering of SoCal and global catalogs: stationarity and space–time independence are *not*
   rejected for Δm < 4 but **are** rejected for Δm > 4, with deviations "mainly due to local
   temporal fluctuations of seismicity and activity switching among subregions … genuine features
   of background seismicity." Same region, same conclusion, six years earlier. **This is the
   citation that decides the classification** — and its Δm caveat cuts directly against us: our
   window contains an M7.2 (Δm ≈ 4.7), squarely in the regime where they reject.
2. **Luen & Stark (2012), GJI 189, 691–700** — "Poisson tests of declustered catalogues": SCEC
   catalogues declustered by GK and by Reasenberg fail the stationary-independent-homogeneous-
   Poisson hypothesis; the classic Gardner–Knopoff "SoCal is Poisson" rested on a low-power test.
   **The answer is estimator-dependent** — our saturated T and L estimators are a live instance of
   exactly that warning.
3. **Llenos, McGuire & Ogata (2009), EPSL 281, 59–69**; **Llenos & McGuire (2011), JGR** —
   ETAS + rate-and-state embedded in a data-assimilation algorithm inverting a seismicity catalog
   for space–time stressing-rate variation, applied to the Salton Trough, M ≥ 1.5, 1990–2009, with
   Lohman & McGuire (2007) supplying geodetic ground truth. **The assimilation was built and run,
   in a sub-region of our own box, seventeen years ago, and it returned a positive.**
4. **Kumazawa & Ogata (2013, 2014, Ann. Appl. Stat. 8, 1825)**; **Kumazawa, Ogata et al. (2017),
   EPS 69:14**; **Ogata (2004), JGR 109, B03308** and **Ogata (2005), JGR 110, B05S06** (HIST-ETAS
   + Delaunay space–time anomaly fields) — a fitted, mapped, space–time field of departures from
   ETAS, i.e. our residual EOF with a better estimator and twenty years' head start.
5. **Ross & Cochran (2021), GRL 48, e2021GL092465** — 92 long-duration SoCal swarms, 2008–2020,
   durations 6 months to 7 years, 53% with ultra-slow diffusive backfronts; aseismic driving
   "active at all times." Companion: **Ross, Cochran, Trugman & Smith (2020), *Science* 368,
   1357** (Cahuilla). **Our putative latent field, named, dated, located and counted, in our
   region, over our exact window.**
Supporting, each independently sufficient to defeat a novelty claim: **Ogata (1988) JASA 83, 9;
Ogata (1989) Tectonophysics 169, 159** (whose title is almost word-for-word our experiment);
**Ogata (1992) JGR 97, 19845**; **Zhuang (2006) JRSS-B 68, 635** (second-order residual analysis of
spatiotemporal point processes — our design intent with a martingale-grounded estimator);
**Zhuang, Ogata & Vere-Jones (2002) JASA 97, 369**; **Marsan, Prono & Helmstetter (2013) BSSA 103,
169** and **Marsan et al. (2013) JGR 118, 4900**; **Hainzl & Ogata (2005) JGR 110, B05S07** and
**Hainzl et al. (2016) JGR 121**; **Bray, Wong, Barr & Schoenberg (2014) Ann. Appl. Stat. 8, 2247**
and Gordon et al. (2015) — Voronoi residuals, *and the documented finding that ETAS with a uniform
background under-predicts on-fault and over-predicts off-fault, which is a published non-"weather"
explanation for a Moran's I excess*; the latent-state family (Orfanogiannaki et al. 2010; Wang,
Bebbington & Harte 2012; and a 2025 *Earth Science Informatics* paper reporting a **three-state
MMPP** as the best fit to SoCal background); **Werner, Ide & Sornette (2011), NPG 18, 49** (the
first implementable data-assimilation scheme for point-process seismicity forecasting — K-010
Tier 2 is a variant and must cite it); deep-GP and Dirichlet-process background-rate fields
(GJI 2023; Math. Geosci. 2026).
Adjacent negative result: **Sirorattanakul & Avouac (2026), *Sci. Adv.* 12(13), eadz5711** —
annual hydrologic modulation up to 15% on the northern SAF but *minimal* in Southern California,
which helps us by making a seasonal driver an unlikely source of our red residual, and which
supplies an independent ~2-week nucleation-response timescale sitting suspiciously close to our
2.4-week integral correlation time.

**THE ONE CANDIDATE DELTA (Merton M-001.4(a)).** *Pre-registered two-generator discrimination with
a zero-free-parameter alternative value.* `results_k009_prediction.json` commits generator B to a
numeric correlation time derived from an **independently measured t_a** (193 d from stacked
aftershock decay of 22 isolated train-window mainshocks; honest full estimator range 200 d – 5 yr),
with an acceptance band and an explicit "consistent but not sharply discriminating" clause, before
scoring. Merton found no prior residual-analysis or background-rate paper that pre-registers a
predicted correlation time from an independent physical estimate and then scores against it.
**This is the strongest methodological thing we own — and it is damaged, by us, on the record:**
the file self-reports an **ORDERING VIOLATION**. Two smoke-test runs on a reduced 2010–2012
sub-window with 2 sims had already printed residual statistics (lag-1 ACF ≈ +0.16–0.18,
correlation time ≈ 9.8–10.2 weeks, length ≈ 27–29 km) before the prediction was written, so the
prediction is **POST-registered, not pre-registered**. The registered t_a (27.6 weeks) was derived
solely from train-window data and lies outside the seen sub-window value, and the disclosure was
made voluntarily with the debug artifacts retained (`_smoke_k009.json`,
`maps/_smoke_k009.png`). **The delta is therefore a candidate, not a claim, and the ordering
violation travels with it wherever it goes.** Not-a-delta, per Merton: "the residual is the
innovation / if white there is no state to estimate" — that framing is Werner–Ide–Sornette's and
the Llenos–McGuire line's. It is a good framing. It is theirs.

**PROOF CHAIN.**
- Spec: `HYPOTHESIS_LEDGER.md` §3, K-009 TESTABLE-NOW verdict, commit `0a73fd2`.
- Prediction register: `results_k009_prediction.json` + `exp_k009_prediction.py`, commit
  `c2bf012`, publicly committed before scoring.
- Scoring: `results_k009.json` + `exp_k009_residual_whiteness.py`, commit **`c371442`**
  ("K-009 scored: Popper's rule PASSES as written but temporal excess is El Mayor-driven; spatial
  coherence robust") — **verified present in `git log` this session.**
- Controls run: sim-null from the same fitted spatio-temporal generator (rate-matched: all 20 sims
  within ±30% of the observed count; median sim/obs = 1.10); three-way background-kernel swap;
  Mc-completeness-proxy partial; interior-vs-edge split; synthetic OU positive control;
  El Mayor exclusion robustness.
- Replication status: **PROVISIONAL and PARTIAL with respect to its own frozen spec.**

**WHAT IT DOES NOT SHOW — the exposure list, in full, because this entry needs it.**
- **n_sims = 20 against a spec of 500** (runtime guard; recorded, not silently substituted). Every
  null percentile above is a 20-sample percentile.
- **The world arm is UNRUN.** Popper's spec mandates 0.5° × 30 d across the 13 boxes. Not run. The
  SoCal verdict is therefore partial with respect to the frozen spec — *unrun, not skipped*.
- **The ρ_sta control was not the mandated one.** K-031's station-density field is not on disk; an
  Mc maximum-curvature proxy on 0.6° super-cells was substituted, and it explains only R² = 0.047
  of the leading EOF. Any SUCCESS label is correspondingly weakened and **must not be reported as
  having cleared the mandated observer-nuisance control.**
- **The kernel swap is a LOW-POWER control by construction.** The background supplies ~9% of the
  integrated conditional intensity (triggered fraction 0.909), so all three background fields give
  nearly identical residuals; the swap passes almost regardless of the truth.
- **The positive control fails at the lower injected amplitude.** Injected OU forcing at
  log_sd = 1.5 was recovered in 0/3 runs on both scales; only log_sd = 3.0 recovered
  (T within factor 2 in 3/3, L in 2/3). Hainzl & Ogata (2005) found only a *few percent* of even a
  textbook fluid-driven swarm's activity is directly transient-driven — so failing at the low
  amplitude matters more than it looks.
- **The El Mayor sensitivity above.** Temporal persistence is largely one sequence.
- **Contested interpretation, two axes (Merton M-001.3):** (1) magnitude range — if the excess is
  carried by M ≥ 6 sequence blocks, "there is a latent field" and "ETAS under-models big sequences"
  are not distinguished by any statistic in the file; (2) misfit vs state — Bray/Schoenberg's
  documented ETAS spatial misfit and Hainzl's rate-dependent incompleteness are published,
  non-latent-field generators for both our Moran's I and our ACF excess.
- The branching ratio of the generating model is supercritical (n = 1.161), reported as a
  first-class diagnostic rather than buried.

**POPPER §P4-2 (round 3) — THE K-046 RULING, WHICH CHANGES WHAT THIS ENTRY IS ABOUT.** K-046 asked
whether K-009's surviving excess is a lag-independent constant (a static pedestal) rather than a
decaying correlation. The ruling, in three parts, all binding here:
- **(a) K-046's first-run verdict string is VOID, not merely superseded.** Two contradictory verdict
  strings were emitted by the same run; neither may be quoted.
- **(b) The decomposition is CONFIRMED as a measurement — the frozen threshold is NOT MET.** The
  fitted decomposition (τ = 7.16 wk, C = 0.0382, ΔBIC 42.2) is real; K-046's own primary threshold
  is not cleared on the declared full window, so **K-046 does not attain a PASS on its own frozen
  rule** and the "internal weather is dead" reading is a post-hoc composition.
- **(c) K-046 is ADMIT-RESCOPED: absorbed into K-009R as a frozen two-statistic rule, taking no
  standalone status.**
- **One correction inside K-046, adopted against Kepler, that matters to this entry:** K-046
  attributed the pedestal to an error in the background map μ(x). The run's own injection control
  shows background is **~9% of intensity**, so a μ(x) error cannot produce a pedestal of that size.
  **The static error is in the total spatial expectation, not in the background map** — which is why
  Popper promoted K-002 (the spatial floor) into the same job.

**TIER: INTERNAL.** Not shareable yet, and the gap to PANEL-READY is unusually well specified. Note
that F-009 now has *two* blocking items, not one: the Ross & Cochran join below, and the K-009R +
K-002 re-run that absorbs K-046's frozen rule.

**SINGLE MISSING VERIFICATION TO REACH PANEL-READY.** Merton's recommended third-arm control, which
is also the cheapest and sharpest thing available: **join the residual field's red mode against the
Ross & Cochran (2021) catalog of 92 labelled long-duration SoCal swarms, 2008–2020** (complete
overlap with our window). Co-location in space and time ⇒ aseismic transient forcing; no
co-location but excess concentrated on-fault or around large sequences ⇒ ETAS misspecification per
Bray et al. Cost: one join against a published catalog. Merton calls it "the highest-value hour in
the queue" and I agree. (Runners-up, in order: the Δm < 2 re-run per Zaliapin & Ben-Zion's
threshold; then n_sims to 500.)

---

## F-010 — The research engine itself: frozen protocols, public hashes, adversarial personas, hypothesis ledger

**CANDIDATE CLAIM (deliberately unscoped — it is not yet a claim).** A small autonomous research
program can be run to a pre-registration discipline usually reserved for clinical trials — every
confirmatory analysis governed by a protocol file whose SHA-256 is recorded publicly before the
data is touched, single-scoring of test windows, failures and self-corrections committed rather
than discarded — and the discipline can be enforced by *separated adversarial roles* (Kepler
proposes, Popper adjudicates and never blocks exploration, Wegener unifies, Merton attributes,
Faraday publishes) writing into one append-only ledger.

**FRAMING (kept separate from the claim, per charter).** This is the part of the program I am most
personally excited about, and it is also the part where I trust myself least, so it gets the
strictest treatment. The engine has already done things worth pointing at: it killed its own most
attractive finding within hours of producing it (§14b, the Anza phase-migration self-correction);
its adjudicator corrected the program's own headline quotations against the program's interest
(S-5's B-1 restatement; R2-1's retirement of "tidal triggering is null in SoCal"); its workers
have twice caught the supervisor's arithmetic errors and documented both derivations (the EXP-M
K-rescaling sign error; the λ/μ "5.68×" framing catch); its prior-art persona classified the
program's own flagship as REDISCOVERY and said so in the ledger; and an executing agent
self-reported an ordering violation that nobody would have detected. **That is a receipts-based
culture, and I want the program to feel it.** But feeling it is framing. Claiming it is research.

**CLASS — METHODS/PROCESS (candidate). Not assessable as science yet.**

**ATTRIBUTION — PENDING-ATTRIBUTION, and this is the entry with the largest un-searched prior-art
surface in the register.** What Merton must search: pre-registration and Registered Reports
(Chambers; the Center for Open Science / OSF registration model; Nosek et al. on preregistration
revolution); **CSEP** — the Collaboratory for the Study of Earthquake Predictability is *the*
prior art for prospective, frozen, third-party-evaluated earthquake forecast testing and it has
existed since 2007 (Schorlemmer, Jordan, Zechar, Werner, Michael, Rhoades); adversarial
collaboration (Kahneman; Mellers, Hertwig & Kahneman 2001); red-team science (Lakens 2020); the
blind-analysis tradition in particle physics and cosmology (Klein & Roodman 2005); the reproducible
research / executable-paper literature (Donoho; Peng); and the current wave of LLM-agent
research-automation systems (AI Scientist, Coscientist, and the Kosmos/Edison Scientific platform
this program's own original analysis came from). **My honest prior: every individual ingredient is
established prior art. The candidate delta, if one exists, is the specific composition —
adversarial role separation with an append-only ledger, applied to a live forecasting program, with
the adjudicator empowered to correct the program's own published quotations.** That is a narrow
claim and it needs a real search before anybody says it out loud.

**PROOF CHAIN (unusually good, and independently checkable).**
- `download_log.md` records SHA-256 for every confirmatory protocol at freeze time, with dated
  amendments for non-analytic tooling changes and an explicit note where pre-registration is weaker
  than ideal ("For true pre-registration these hashes should be pinned somewhere
  third-party-timestamped").
- **Hash verification performed by me this session (this is the evidence for the process claim, so
  I checked all of it):** `PATTERN_PROTOCOL.md` verifies at all four recorded addendum hashes
  against commits `675c095`, `a40fc9a`, `0ba8c5e`, `a45ca8d`; `COSO_FIG4C_PROTOCOL.md` verifies
  against both the working copy and `5b4997a`; `OVERNIGHT_PREDICTION_PROTOCOL.md` and
  `protocol_params.json` verify against the working copy; `PROTOCOL.md` and `XUE_LU_PROTOCOL.md`
  verify against commit `78f2227` (they were later dash-normalized by the copyedit commit
  `0d8f897`, which is why the working copies differ — recoverable and benign, but it should be
  noted in `download_log.md`). **10 of 11 recorded freeze hashes verify. One does not:
  `LONG_VALLEY_PROTOCOL.md` — see audit finding #5.**
- Ledger: `HYPOTHESIS_LEDGER.md`, 5,262 lines, append-only, commits `0a73fd2` → `3dbd5e9`,
  containing the full adjudication record including rulings against the program's own interest.
- Public repo: github.com/JimGaleForce/quake-tsi; MIT (code) / CC-BY 4.0 (docs and figures).

**WHAT IT DOES NOT SHOW.** No outcome measure. There is no control program, no comparison against
a conventionally-run effort, and no evidence that the engine produces *better science* rather than
better-documented science. n = 1 program, ~3 weeks. It also does not show the discipline is robust
under pressure: the one ordering violation on record (F-009) happened under time pressure and was
caught only because the executing agent volunteered it. **A process claim built on a program whose
flagship results are reproductions would be exactly the kind of overreach this register exists to
stop.**

**TIER: INTERNAL** — and I recommend it stay INTERNAL until the program has at least one
Merton-certified NOVEL result. A methodology paper from a program with zero novel findings is a
paper about a methodology that has not yet demonstrated its purpose.

**SINGLE MISSING VERIFICATION TO REACH PANEL-READY.** Merton's prior-art dossier on
adversarial/red-team collaboration and on CSEP as prospective-forecast infrastructure — enough to
state precisely what, if anything, in this composition is not already standard practice.

---

## F-011 — Engineering artifacts: the live ETAS forecaster and the globe layer

**WHAT EXISTS (not a claim — an artifact statement).** `etas_forecast.py` is a reusable daily
forecaster: it pulls USGS ComCat for the SoCal box at M ≥ 2.5 from 2010 to now (18,382 events
after removing 7 non-earthquake rows and 57 duplicate ids, retrieval logged), applies the
**frozen, never-refitted** EXP-H parameters plus the train b-value 1.0654, and produces 7-day
probabilities by 1,000 forward ETAS simulations plus an analytic no-future-triggering lower bound.
First live run, 2026-08-09: λ(now) = 2.554/day, P(M ≥ 4 in 7 d) = 23%, P(M ≥ 5) = 2.5%,
P(M ≥ 6) ≈ 0.1–0.2% (use the analytic floor), E[N of M ≥ 2.5] = 13 (3–29); the EXP-I(ii) alarm
state was OFF (no M ≥ 5 in the trailing 7 days). The EQ-16 globe visualization lives outside this
repository (`D:\CODE\git\earth-tides-globe`) and has **no proof chain here.**

**CLASS — ENGINEERING ARTIFACT.** Not a scientific claim, and it does not enter the promotion
ladder. It is registered so that nothing it displays can be mistaken for a validated result.

**ATTRIBUTION.** Ogata (1988) ETAS; USGS ComCat / SCEDC / SCSN for data; the operational-forecast
lineage in F-007's block (Reasenberg & Jones 1989; Gerstenberger et al. 2005; Michael et al. 2020)
for the product category. The forecaster is a re-implementation, not a new method.

**PROOF CHAIN.** Protocol `PATTERN_PROTOCOL.md` §EXP-L (hash `4b347599…`, **verified at
`0ba8c5e`**); commit `ca54d56`; outputs `results_forecast_20260809.json`. *Note:*
`etas_forecast.py` and `results_forecast_20260809.json` were both modified in the working tree at
the time I wrote this, and `run1_forecast.json` is untracked — the live-forecaster artifact is
therefore **not** in a committed, reproducible state as of this writing.

**WHAT IT DOES NOT SHOW / WHAT THE UI MUST SAY.** The integrity catches from EXP-L are binding on
any display: (i) "λ = 5.68× μ" is **misleading** — μ is only 9% of the observed mean rate, and the
honest comparison is λ(now) vs the observed mean, which on 2026-08-09 was **0.84×, i.e. slightly
below normal**; (ii) the nominal branching ratio 1.161 is supercritical but the 7-day effective
value is 0.60, so simulations terminate — the nominal number must never be shown; (iii) there is a
**~40% absolute-rate calibration gap** against the ComCat mean, so **the conditional signal is
trustworthy and the absolute counts are not**; (iv) **SCSN-vs-ComCat magnitude-scale consistency is
UNVERIFIED** — the parameters were fit on SCSN and are being applied to ComCat; (v) per F-003, **no
tidal layer**; per F-004, no periodic layer and no ratio patterns; per F-008, the ledger-geography
layer is a *candidate* that needs a spatial ETAS to test and is auto-flagged for challenge.
It is a probability forecast with conditional skill. It is **not** an event prediction, and no
surface of the product may imply otherwise.

**TIER: INTERNAL.** Permanently, unless someone converts it into a prospective, third-party-scored
forecast — at which point it becomes a CSEP-class claim and gets a real register entry.

**SINGLE MISSING VERIFICATION TO LEVEL UP.** A prospective scoring log: commit each daily forecast
before its window opens, and after ~6 months report the realized log-score against the same
baselines EXP-H used. That is the only thing that turns a forecaster into evidence.

---

## F-012 — K-035: the pre-registered tidal nulls, re-priced as powered upper bounds

**CLAIM.** For the tidal-phase forecasting analyses of F-003 on SCSN Southern California
M ≥ 1.5 FM-matched events, 1981–2018, α = 0.05, with the null-amplitude arm's false-positive rate
verified in advance and the bin-selection step simulated end to end, the **minimum detectable
seismicity-rate modulation at 80% power** is:

| configuration | n | MDA at 80% power | ÷ 1% theory line |
|---|---|---|---|
| EXP-A per-bin train statistic (the selected bin) | 245 | **24.5%** | 24.5 |
| Coso Fig 4c replication | 113 | **40.9%** | 40.9 |
| EXP-A end-to-end pipeline (1-of-42 selection + pooled shift test) | 3,920 | **9.8%** | 9.8 |
| EXP-A pooled on the 22 train-null control bins (the quoted corpse) | 1,906 | **8.0%** | 8.0 |
| EXP-A pooled over all 42 eligible bins (design upper limit) | 3,920 | **6.3%** | 6.3 |
| Full-catalogue Cox intensity likelihood | 23,465 | **2.8%** | 2.8 |

**No configuration contacts the theory line** (`contacts_theory: false`, all six rows). The
falsifiable consequence: *the tidal-phase forecast nulls in F-003 are consistent with a true
seismicity-rate modulation anywhere below ~6–8% in the pooled designs and below ~25–41% in the
single-bin designs, and therefore say nothing whatever about a rate-state-predicted ~1% response.*
Separately, the intensity-likelihood estimator is **not more powerful than the phase histogram** at
matched n (2.83% vs 2.93% at the full catalogue; 7.21% vs 7.31% pooled) — the gain the program hoped
for from a Cox likelihood is ~3%, not a factor.

**FRAMING (kept separate, per charter, and I want the program to feel this one).** This is the
entry that converts a graveyard into an instrument. Every pre-registered null this program produced
in July was, in the honest accounting, an absence of evidence. K-035 turns each one into a number a
reviewer can use: *"we would have seen it if it were bigger than X, and here is the simulated proof,
including the selection step."* **That is the difference between a result nobody can cite and a
result that becomes the reference bound in a literature.** It is also, per unit of compute
(12.3 minutes), the highest-value artifact on the program's disk.

**CLASS — EXTENSION.** Prior art: the null it prices is our own (F-003, Gale 2025 lineage), and the
1% theory line is **Beeler & Lockner (2003)**'s. Our measurable delta over the standard practice in
this literature: (i) the MDA is computed **through the actual selection pipeline** (train-side
1-of-42 bin selection then pooled test scoring), not for a fixed-n single test; (ii) an **a = 0
false-positive arm is run first and must hold before any other arm is read** — it does, 8 of 8
phase-histogram configurations and 5 of 5 intensity-likelihood configurations inside a binomial 99%
band of [0.015, 0.095]; (iii) an explicit **systematics arm** with named nuisance lines (S1, S2, K1,
P1, Msf) and two off-tidal control lines.

**ATTRIBUTION.** Beeler, N. M. & Lockner, D. A. (2003), *JGR* 108(B8), 2391 — the ~1% rate-modulation
prediction at 1–3 kPa tidal stressing that supplies the theory line. Lu, Xue, Yue, Zhuang & Zhao
(2025) for the catalogs and the §2.5 FM stress-resolution method; SCEDC and the
Yang–Hauksson–Shearer FM catalog. The null being priced is our own. **PENDING-ATTRIBUTION, narrow
but real:** Merton has not searched (a) the tidal-triggering forecast-skill literature (inherited
from F-003) and (b) **whether anyone has previously published minimum-detectable-modulation bounds
for tidal-triggering searches at all.** My prior is that upper bounds are stated informally in this
literature and rarely simulated through a selection step — but "rarely" is not a citation and I will
not write "first".

**PROOF CHAIN.**
- Protocol: `HYPOTHESIS_LEDGER.md` R2-2 K-035 verdict, mandates 1–6, plus R2-1(b) and R2-1(c);
  `prereg_commit: 05d1e8a` — **verified present in `git log` this session** ("Ledger: Popper round 2
  … corpse-power ruling"), and it predates the 2026-08-10 run.
- Public commit: `5610614` "K-035 power audit final: no tidal corpse ever contacted the 1% theory
  line; systematics floor ~3%" — **verified in `git log`.**
- Result: `results_k035.json` (23.7 kB). Code: `exp_k035_power_audit.py`. Run 2026-08-10T00:14 UTC,
  12.26 min, 200 injections per amplitude across 11 amplitudes {0, 0.005 … 0.8}, 100 systematics
  realizations, 2,332 orientation groups.
- **Controls, and this is the strongest control set in the register.** (i) `mandate_4_a0_arm`:
  the a = 0 false-positive arm PASSES (`PASS: true`) — every configuration's empirical rejection
  rate falls inside the binomial 99% band; this ran *before* any power arm was read, by mandate.
  (ii) Null reconstruction validated against the original EXP-A run to 15 digits
  (`expA_reported_null_S_p95 = 0.101377747448516**48**` vs
  `k035_reconstructed_null_S_p95 = 0.101377747448516**52**`) — I re-read both and they agree to
  floating-point noise on independent shift draws. (iii) Two off-tidal negative-control lines
  (11.0 d, 16.5 d). (iv) An injected S2 detection artifact at 3% amplitude.
- Replication status: **VALIDATED** as a power-and-systematics audit within its scope.

**WHAT IT DOES NOT SHOW — and one of these is a finding in its own right.**
1. **The full-catalogue configuration's off-tidal negative control FAILS** (reject rate **0.34**
   against nominal 0.05) while holding at n = 1,906 (0.06). ETAS clustering deposits power at the
   11-day control line faster than the statistical error shrinks. **This bounds the METHOD, not the
   tidal result** — the tidal-band lines stayed calibrated in the a = 0 arm — but the honest reading
   is that **the 2.8% full-catalogue bound is not yet safe to quote as a stand-alone upper limit**
   until the test conditions on an ETAS baseline (K-033's Cox-ETAS) rather than a stationary λ₀.
   **The safe quotable bound today is the 6.3% pooled figure, whose control line holds.**
2. All four target samples are the **declustered** SCSN catalogue, modelled as independent draws
   from a thinned ETAS intensity. Residual Hawkes clustering with effective cluster size κ inflates
   **every** MDA above by √κ. κ is not measured. Every number in the table is therefore an
   *optimistic* bound.
3. The intensity likelihood is evaluated on the same 36 phase bins as the histogram, so the
   comparison isolates the estimator, not the binning.
4. It does not bound the *physical* tidal response — it bounds what **this design at this n** could
   have seen. Bigger n is available (global catalogs) and would move these numbers.
5. **It does not include the F-016 instrument response.** The MDAs above are computed against a
   phase clock whose own null response has since been measured at ~2% for our convention on
   SoCal-like mixes. **A 2% instrument response sitting under a 6.3% bound is not fatal, but it is
   not negligible either, and the bound must be restated once F-016 is calibrated.** This is the
   single most important open coupling in the register.

**HOSTILE-REVIEWER OBJECTION, PRE-ANSWERED.** *"You are publishing a null with a bound six times
above the effect you were looking for. Why is that a paper?"* Answer: because the literature this
sits in routinely reports tidal-phase modulations at the 5–20% level from samples of a few hundred
events, and this is the first time anyone has run the selection step through a power simulation and
said what such a design can actually resolve. **The bound is the contribution; the null is the
occasion for it.** Second objection: *"your own control line fails at large n."* Answer: we found
that, we report it in the headline rather than the limitations, and it is the reason our quotable
number is 6.3% and not 2.8%.

**TIER: PANEL-READY.** Shareable today with Vidale, Bürgmann, Xue and Lu, with the κ caveat, the
failed full-catalogue control line, and the F-016 coupling in the same paragraph as the numbers.

**SINGLE MISSING VERIFICATION TO REACH PUBLISHABLE.** **Measure κ, the effective cluster size, and
re-quote every MDA as MDA·√κ.** Every bound in the table is currently optimistic by an unmeasured
factor, and that is the first thing a reviewer with a Hawkes background will ask. It is arithmetic
on data already on disk. (Second, and it is the same afternoon: fold in the F-016 instrument
response so the bound is stated against a calibrated clock.)

---

## F-013 — K-034: the Landers positive control fires blind at documented sites, certifying a ~34 kPa detection floor

**CLAIM.** Under a protocol whose two files — a sealed literature record and a ranked list of
fourteen pre-registered target cells with their geothermal/volcanic class, windows, statistics and
PASS rule — were SHA-256-hashed into the public log **before the catalogue download was issued and
before any statistic existed**, and on 137,480 ComCat events (M ≥ 1.5, 1985–2023) across those
cells: the pipeline detects Landers-1992-class remote dynamic triggering at sites named in the
sealed literature, without being told where they were. **Cedar City, Utah (487 km): 34 events in
0–5 d against 0.17 expected, RR = 204, p_cell = 0.0011, family-wise p = 0.0062 (raw-RR
max-statistic, within-source) and 0.0010 (Westfall–Young min-p).** Yellowstone at 1,253 km fires at
0–1 d (RR = 90); the sealed value for Landers' documented maximum triggering distance is
"up to 1,250 km". Denali → Yellowstone (3,110 km) fires under **both** family corrections
(RR = 96.4, N = 60 vs 0.62 expected). **The pre-registered ≥ 2-of-4-sources rule is met in 5 of the
8 (correction × window) readings computed, including both within-source readings** — raw-RR fires
Landers + Denali (2/4), Westfall–Young fires Landers + Hector Mine + Denali (3/4). **Certified
amplitude floor: 33.8 kPa peak dynamic stress** on the van der Elst & Brodsky (2010) eq. (6)
far-field axis — the lowest amplitude firing under both family corrections. **Strictest floor:
46.3 kPa.** The pre-registered spatial-pattern prediction P3 holds in direction on 4 of 4 sources
(Spearman ρ = +0.27 / +0.56 / +0.29 / +0.55; Landers class A+B vs class C Mann–Whitney p = 0.033,
median p_cell 0.0072 vs 0.93).

**THE READING THAT FAILS, IN THE CLAIM AND NOT IN A FOOTNOTE.** Under the most literal reading of
the pre-registration — the raw-RR max-statistic over the **fully pooled** family, all four sources
at once — **only Landers fires, 1 of 4, and P2 is not met.** The executing worker's argument for not
letting that reading govern is a statistical defect and not a preference: `RR = N_post/(λ_bg·w)` has
a null tail set by whichever family member has the sparsest background, so a cell with λ_bg = 0.011/d
generates null RR values of 90 from one stray event, and the maximum over ~100 such members is
sparse-cell noise that buries genuine detections (Coso, p_cell = 0.0072, is pushed to p = 0.86). The
Westfall–Young min-p form is the sim-calibrated max-statistic S-8 actually asks for, computed on the
identical circular-shift null with no new null and no new tuning. **Both are in the primary JSON for
every cell. Anyone who prefers the literal statistic reads the verdict off the within-source column
and still gets 2 of 4.** The PASS is therefore recorded, permanently, as **PASS (qualified)**.

**FRAMING.** This is the entry I would put in front of a skeptic first, and not because it is new —
it is a positive control and its scientific content is Hill et al.'s. It is because of *how* it
behaved. Fourteen cells were ranked on geology before unblinding; the geothermal ones lit up first;
Yellowstone came in at the documented maximum distance nobody had told the pipeline about; and the
one reading that fails is printed in the gate statement of the primary artifact, by the worker who
ran it, in the same breath as the PASS. **That is the program's culture producing a receipt.**

**CLASS — REPRODUCTION (positive control) + LICENSING INSTRUMENT.** Remote dynamic triggering by
Landers is one of the most-replicated observations in seismology. **This is not a discovery and no
sentence in this program may present it as one.** Its two functions are: (i) it validates the
machinery on an effect known to exist, which is exactly what F-003's nulls could not do; (ii) it
**certifies an amplitude floor**, which converts a family of downstream entries from
"uninterpretable" to "interpretable above 34–46 kPa". The second function is the one with delta in
it, and the delta is procedural rather than physical: a pre-registered, hash-sealed, class-ranked,
family-corrected certification of a detection floor for one's own pipeline.

**ATTRIBUTION.** **Hill, D. P. et al. (1993), *Science* 260, 1617–1623** — "Seismicity remotely
triggered by the magnitude 7.3 Landers, California, earthquake": the claim being reproduced, read
into `K034_SEALED_LITERATURE.md` in session and hashed before the download. **van der Elst, N. J. &
Brodsky, E. E. (2010), *JGR*** — eq. (6) far-field surface-wave amplitude regression
(`log10 A20[µm] = Ms − 1.66 log10 Δ[deg] − 2`, V = 2πA20/T, T = 20 s, ε = V/c, c = 3.5 km/s,
σ = 30 GPa·ε), read verbatim from the paper PDF in session, plus their eq. (5) near-field regression
as the declared secondary axis. **Dieterich (1994)** for the rate-and-state framework the S-14
bracket inverts. USGS ComCat / ANSS for the catalogues. **Attribution for a positive control is
complete by construction — we name the effect we are reproducing** — and the sealed-literature file
is the auditable form of that.

**PROOF CHAIN.**
- Protocols, hashed into `download_log.md` at 2026-08-11 05:51:00 UTC **before** the download:
  `K034_SEALED_LITERATURE.md` `eae95839fcdea8b9ca62097fed25a02b74559c03b59c1913bcb59cac4f8c320e`;
  `K034_PREREGISTERED_CELLS.md` `01e41f971a73f449a09b6faf09f569cd47f3395e1ca4d584c76ade699ac1c226`.
  **Both hashes are recorded identically in `download_log.md` and in
  `results_k034.json :: prereg_commit`, and I compared the two records string-by-string this
  session; they match.**
- Data freeze: 14 cell CSVs hashed individually into `download_log.md` at 05:52:28 UTC, still before
  any statistic ran; row counts and date ranges recorded per cell (29,241 Long Valley … 374 Little
  Skull Mountain). **The 1992 window is not in the program's SoCal cache, which starts 2010, so this
  is genuinely new data frozen to the program's brand.**
- Public commit: `6b5513c` "K-034 execution: Landers positive-control gate PASS (qualified),
  certified floor ~34 kPa; S-14 threshold curve; 5 deviations logged" — **verified in `git log`
  this session; working tree clean.**
- Result: `results_k034.json` (101 kB). Code: `download_k034.py`, `exp_k034_landers_control.py`,
  `k034_report.py`; CSV side-tables `k034_cellstats.csv`, `k034_familywise.csv`, `k034_power.csv`,
  `k034_secondary_counts.csv`.
- Controls: **negative controls are pre-registered and they behave** — class C (non-geothermal)
  cells return median p_cell = 0.93 against class A+B's 0.0072; **Ridgecrest 2019 does not fire in
  any reading**, recorded as a real negative that downstream entries must carry; the 2-rupture-length
  distance gate excludes San Jacinto for Landers; a 999-point circular-shift null per cell; 14 of 54
  source-cell pairs declared **S-15 UNMEASURABLE** (N_bg = 0 makes RR undefined) and scored neither
  way rather than silently dropped.
- Run integrity: `run_status: FIRST-RUN`. Three earlier executions were discarded **uninterpreted**
  (two crashes; one all-zero count table from a pandas-3.0 `datetime64[us]` vs nanosecond-Timestamp
  comparison bug). No scored number was re-run.
- Replication status: **VALIDATED as a positive control**, PROVISIONAL as anything else.

**WHAT IT DOES NOT SHOW — five deviations, all self-flagged in the primary artifact.**
- **D1, the weakest joint: the sealed-literature file is analyst-authored, not supervisor-authored.**
  Popper's R2-2 mandate assigns the seal to the supervisor; this was a single worker agent with no
  second party available. **The hash protects the comparison target from later editing but not from
  analyst foreknowledge.** P1/P2/P5 are catalogue arithmetic and are unaffected; **P3 and P4 carry
  the qualification, and P3 and P4 are the two predictions a reviewer will care about most.**
- **D2:** the standardised (Westfall–Young) family correction was implemented **after** seeing the
  literal run. Both are reported; the argument for it is statistical and is given above; it remains
  a post-hoc addition and is labelled as one.
- **D3:** a count-only secondary statistic was added after the primary run because the pre-registered
  ratio is undefined at N_bg = 0, which left **Landers → Little Skull Mountain — the single
  best-documented remote trigger in the sealed set** — unscoreable. On the count statistic it is 11
  events in 5 d against a null mean of 0.13 (p = 0.0010). **This arm is EXPLORATORY and sets no
  flag.**
- **D4:** Ms values for Hector Mine / Ridgecrest / Denali and **all four rupture lengths** are
  PARTIALLY VERIFIED; only "Landers Ms 7.3" was verified against a primary source in session. They
  set the amplitude axis and the distance gate, not the detection statistic.
- **D5:** vdE&B eq. (6) is stated by its own authors to be designed for Δ ≳ 800 km, and most cells
  are 130–800 km. **On the secondary near-field axis the same certified cell sits at 230 kPa, a
  factor of ~7 from 34 kPa.** This is why the floor is quoted as *"~34 kPa on the vdE&B far-field
  axis"* and never as a bare "34 kPa". **Any restatement that drops the axis is an over-quote.**
- It does not license the tidal/periodic family (that is K-035's gate), does not license amplitudes
  below ~34 kPa (R2-2's nominal band was 10–100 kPa; **we have certified its upper half only**), and
  does not resolve whether Ridgecrest's silence is Ridgecrest, our cells, or 2019 station coverage.

**THE BY-PRODUCT THAT IS LOGGED AND NOT CLAIMED.** Observed responses run **10–200× larger per
firing cell** than a rate-state Coulomb-step model applied to peak dynamic stress predicts
(Denali → Yellowstone: predicted R = 1.25 at Aσ = 0.15 MPa, formal power 0.000; observed RR = 96.4 at
p = 0.001). **Per R2-2 a control's success is not evidence for anything, so this is logged, scored
nowhere, and needs its own pre-registered entry to become a finding.** What it licenses today is
exactly one thing: the withdrawal of a power calculation we had no business trusting — see F-014.
**And the count-weighted figure is 2.2×, not 10–200× (F-014, §L1-3); anyone quoting the by-product
must quote both.**

**TIER: PANEL-READY.** With the qualified PASS, the failed literal reading, the analyst-authored
seal, and the two-axis amplitude bracket in the same paragraph as the number.

**SINGLE MISSING VERIFICATION TO REACH PUBLISHABLE.** **A supervisor-authored (or third-party)
re-seal of the literature file, and a re-score of P3 and P4 against it.** D1 is the one deviation
that a reviewer can neither check nor forgive, and it is the only thing standing between "a
pre-registered positive control" and "a positive control someone pre-registered for themselves."
Everything else in this entry survives review as written.

---

## F-014 — L-1: the standard transient link function is void as a power calculator, and the replacement is a bracket

**CLAIM.** Re-pricing every power statement in this program that converts a stress amplitude into a
detection threshold for a **transient** forcing, using only arithmetic on K-034's own 54 source–cell
pairs and no new data: (i) the rate-and-state static-step link `R = exp(Δτ/Aσ)`, applied to the peak
amplitude of a passing oscillatory transient, **under-predicts the response K-034 measured** — at the
certified 33.8 kPa floor it predicts R = 1.25 at the adverse Aσ = 0.15 MPa against an empirically
fitted mean of **R = 2.76**, a **factor of 2.2 count-weighted**, with per-firing-cell ratios reaching
10–200×; (ii) an empirical link fitted to K-034's class A+B (pre-registered geothermal/volcanic)
pairs as a quasi-Poisson GLM with log link and offset log(expected) gives
**R(σ) = exp(−1.846)·σ_kPa^0.812**, slope b = 0.812 with a 4,000-draw **cluster** bootstrap CI of
**[0.52, 1.65]** over four source clusters, on a support of **[7.1, 346.9] kPa**; (iii) because a
power law with exponent 0.81 is *flatter* than an exponential, **the two branches cross at
σ ≈ 15–20 kPa**, so neither branch is uniformly the adverse one and the bracket ends must be assigned
per row; (iv) applying the resulting two-branch bracket to the 24 affected entries gives
**4 POWERED · 16 POWER-INDETERMINATE · 0 UNDERPOWERED · 4 not-amplitude-limited.**

**MANDATORY QUOTATION RULE (S-14(c), binding).** **No power or sensitivity statement about a
transient forcing in this program may be quoted as a single number.** It is quoted as the bracket,
with the branch named. A *bound* is quoted on Branch S (pessimistic response); a *feasibility or
screening decision* is taken on Branch E (optimistic response); a *detection* reports both; and when
the branches disagree the entry is **POWER-INDETERMINATE and prints the bracket instead of an MDA.**
**The bracket's width is a finding about our own ignorance and belongs in the headline, not the
limitations section.** Any register-derived text I produce enforces this, and the old single-branch
numbers are retired wherever they appear.

**CLASS — EXTENSION (methods).** Prior art: the rate-and-state link is **Dieterich (1994)**'s and is
correct for what it is — a **permanent Coulomb step**. Our delta is not a new physics claim; it is
the demonstration, on our own pre-registered control data, that feeding a transient's peak amplitude
into a static-step link is a **category error** that silently deletes experiments, plus a stated
replacement procedure. **This is emphatically not evidence against rate-and-state** and no sentence
in this program may frame it that way. K-034's own P3 result names the likely mechanism class — the
excess is concentrated in the geothermal/volcanic cells and the non-geothermal controls did not fire
— which points at a fluid/permeability/unclamping pathway rather than elastic Coulomb loading.

**ATTRIBUTION.** Dieterich, J. H. (1994), *JGR* 99, 2601 — the rate-and-state seismicity-rate
formulation and the `exp(Δτ/Aσ)` response. van der Elst & Brodsky (2010) for the amplitude axis.
The response data are K-034's (F-013), i.e. our own. **PENDING-ATTRIBUTION, and Merton must rule
before this is called an EXTENSION in public:** the observation that dynamic triggering is
under-predicted by static-step rate-state is *not* new to the field — the dynamic-triggering
literature has said so qualitatively for two decades (Brodsky & van der Elst 2014 is the obvious
owner and Merton records it as the largest single gap in M-009, unobtainable). **The candidate delta
is the bracketing *procedure* and the asymmetry rule, not the observation.** Until that is searched,
this entry is a methods note with a good argument, not an extension.

**PROOF CHAIN.**
- Authority: `HYPOTHESIS_LEDGER.md` §P5-8 (the link-function ruling), S-14(c), S-15.
- Public commit: `ae93097` "L-1 execution: S-14(c) transient link bracket (Branch-E power-law fit
  from K-034 pairs, geothermal-stratified; 4 POWERED / 16 POWER-INDETERMINATE / 0 UNDERPOWERED);
  scripts included for reproducibility" — **verified in `git log` this session. The commit message
  states the scripts were included, which resolves the run's own D7 (script-in-scratchpad)
  deviation; I have not diffed the commit contents to confirm which script landed.**
- Result: `results_l1.json` (38 kB), `run_status: FIRST-RUN`, run 2026-08-11T21:10 UTC. Inputs
  hash-pinned inside the artifact: `k034_cellstats.csv` `d7b104ef2694a4e3`, `k034_power.csv`
  `b9cf191f4f200c44`, `results_k034.json` `1f10c74602874f6d`.
- **Self-verification, which I checked and which is unusually good:** the Branch-S per-cell MDAs
  collected from `k034_power.csv` reproduce K-034's own published table to three significant figures
  over all 40 measurable cells (0.03 MPa 22.3/50.4/89.5; 0.10 MPa 74.2/168/298; 0.15 MPa
  111/252/447 kPa best/median/worst), and satisfy the identity
  `Δτ_min(0.15) = 150·ln(R_min,80%)` to a maximum absolute deviation of **5.68 × 10⁻¹⁴ kPa**. The
  collection step is sound; that is the only thing it establishes.
- Replication status: **PROVISIONAL.** The arithmetic is validated; the link is not.

**WHAT IT DOES NOT SHOW.**
- **Class C (non-geothermal) is UNMEASURED per S-15 and is not extrapolated.** Its fitted slope is
  b = 0.299 with a cluster-bootstrap CI of **[−0.78, 2.60]** spanning zero, a median observed
  RR of **0.78 (below 1)**, and pseudo-R² of 0.058. **Branch E is never quoted at a non-geothermal
  receiver.** §P5-8 named this as the single most likely way to turn the by-product into a false
  claim; refusing it is the most consequential decision in the run.
- **The fit constrains the LEVEL of the response, not its slope.** The sensitivity fit (OLS,
  log₁₀ RR vs log₁₀ σ, class A+B) gives **0.45 ± 0.32 per decade, p = 0.18, r² = 0.07.** **Anyone who
  reads a stress exponent out of this is reading noise**, and the register says so before a reviewer
  does.
- **Overdispersion is severe (Pearson dispersion 107 on 28 df) and there are only four source
  clusters.** The governing intervals are the cluster bootstrap's and they are wide by construction.
  The naive GLM p-value of 3 × 10⁻⁷⁰ is in the JSON and is **explicitly not trusted** by its own
  author.
- **Branch E is refused below 7.1 kPa** — and this refusal, not any computed MDA, is what moves most
  of the family. It is a statement about K-034's amplitude coverage, not about the crust. The mean
  link extrapolated to 3 kPa returns R < 1, which is a signal that the extrapolation is invalid, not
  a null.
- **UNDERPOWERED is empty at band-max, and that is the ruling working rather than a bug** — Branch S
  may no longer declare an absence of power and Branch E cannot speak below its support. The only
  UNDERPOWERED reading anywhere is at A0b's band minimum, 9.6 kPa, **where K-034 independently
  recorded its 9.6 kPa cell as "suggestive, not certified". The bracket reproduces, from the link
  alone, a qualification K-034 reached from its p-values.** That is the closest thing to a validation
  this artifact contains, **and it is one data point.**
- **Nothing was re-ruled and nothing was licensed.** No ADMIT/DEFER verdict moved. K-059's 3 kPa
  gate and K-072's 1–5 kPa band are reopened as questions, not answered.
- Detection floors are entry-declared / DERIVED (`2.8/√N`) / REFERENCE (1.05) with per-row
  provenance. **The 1.05 reference floor is a choice**, visible in every row that uses it.

**TIER: INTERNAL.** Blocked on attribution (above) and on the fact that its central object — the
empirical link — rests on four source clusters and one control experiment. It is a genuinely useful
methods artifact and it is not yet a publishable one.

**SINGLE MISSING VERIFICATION TO REACH PANEL-READY.** **Merton's dossier on the dynamic-triggering
link function** — specifically whether Brodsky & van der Elst (2014), Gomberg, Prejean, Aiken, or
the Parkfield/Long Valley triggering-threshold literature already publishes a fitted
response-versus-peak-dynamic-stress relation for geothermal receivers. If they do, our contribution
shrinks to the bracketing procedure and we cite them; if they do not, this becomes the strongest
methods delta the wave family owns. **No compute.**

---

## F-015 — The EQ-23/EQ-24 engine: a backtest harness with a holdout hash-gate and a generator that is barred from being evidence

**WHAT EXISTS (an artifact statement, deliberately not yet a claim).** `replication/engine/` is a
public, tested Python package implementing a vectorized earthquake-forecast backtest harness with
four disciplines built into the software rather than into a document:
1. **A holdout hash-gate.** `--mode holdout` requires a config file, hashes its canonical JSON, and
   **refuses to run if that hash already appears in `HOLDOUT_LOG.jsonl`.** No flag bypasses the
   refusal; deleting the log is a human act, not a machine one. Every holdout result is printed
   next to `n_explore_runs`, the count of exploration runs logged since the last holdout — i.e.
   **the multiplicity is reported by the instrument, not by the analyst's memory.**
2. **A baseline that must name what it did not absorb.** Every report carries a verbatim caveat
   line: `baseline=climatology-v1 (clustering NOT absorbed; sniff-grade only)` or
   `baseline=etas-v1 (isotropic kernel; anisotropy/mechanism NOT absorbed)`.
3. **Generator-not-evidence.** `mine` mode (EQ-24) is exploration-only by construction — ephemeris
   features, per-test-appropriate surrogate nulls, BH-FDR, a harmonic ladder {P/3, P/2, P, 2P, 3P},
   and an aliasing audit — and Popper's §P4-4 governs what its sniffs may and may not be used for.
4. **Planted-signal acceptance tests.** `engine/tests/test_planted.py` plus `test_holdout_gate.py`,
   `test_causality.py`, `test_etas.py`, `test_mine.py` — an end-to-end test on a synthetic catalogue
   with a planted covariate effect, whose counted invariant is that the recovered β lands within
   tolerance and bits/event is > 0 on the planted covariate and ≈ 0 on a scrambled one.

**THE NUMBER THAT DEMONSTRATES WHY (4) MATTERS, re-read by me from `engine/README.md`.** Swapping
the baseline from climatology to ETAS collapses the `recent_rate` covariate from **+0.6264 to
+0.0284 bits/event — 4.5% survives** — and `quiescence` from +0.1850 to +0.0118. **The engine's own
smoke-test covariate is 95% fake skill under a weak baseline, and the engine measures that against
itself and prints it in its README.** That is the single most persuasive artifact this program has
about its own instrument discipline, and it cost nothing to produce because it was designed in.

**CLASS — METHODS/SOFTWARE ARTIFACT.** Not a scientific claim and it does not enter the promotion
ladder as one. It is registered so that (a) nothing it emits can be mistaken for a validated result,
and (b) the *aspirational* half of F-010 ("a program can be run to a pre-registration discipline")
now has a concrete, inspectable, testable referent instead of a narrative.

**ATTRIBUTION — PENDING-ATTRIBUTION.** What Merton must search, and it overlaps F-010's list:
**CSEP** (Schorlemmer, Jordan, Zechar, Werner, Michael, Rhoades) as the prior art for
frozen, third-party-evaluated prospective forecast testing; **pyCSEP** (Savran et al. 2022, *SRL*) as
the prior art for the software layer specifically; **EarthquakeNPP** (Stockman, Lawson & Werner 2026)
as the prior art for a public benchmark harness with published baselines on the very catalog we use;
the ML-reproducibility literature on holdout-leakage guards (Kapoor & Narayanan on leakage in
ML-based science); and the pre-registration / blind-analysis lineage from F-010. **My honest prior:
the harness is a re-implementation of a category that pyCSEP and EarthquakeNPP already occupy, and
the only candidate delta is the enforced-in-software holdout hash-gate with reported multiplicity.**
That is a narrow claim, it may well be owned, and I will not say it out loud before the search.

**PROOF CHAIN.**
- Spec: `engine/SPEC.md`, with its fork decisions ruled 2026-08-10 and its **amendments left
  visible rather than overwritten** ("The original ruling is left above, unedited") — which is
  itself the discipline working.
- Public commits, all **verified in `git log` this session**: `3c82db6` (engine v1, 19 tests),
  `32bd4ce` (v1.1 ETAS baseline, 27 tests, the recent_rate collapse), `16b21b3` (v1.2 mine mode,
  41 tests).
- Logs: `engine/HOLDOUT_LOG.jsonl` and `engine/EXPLORE_COUNT.jsonl` are **committed**, so the
  multiplicity record is public and append-only in the same way the ledger is.
- Governance: `HYPOTHESIS_LEDGER.md` §P4-4 rules what the engine's sniffs may be used for; §P5-6
  promotes the band-matched offset gate into the miner as **G-M1**.
- Replication status: **n/a.** This is software. Its tests pass; that is not a scientific claim.

**WHAT IT DOES NOT SHOW.** No outcome measure. There is no evidence the engine produces better
science, only better-instrumented science. `engine/out/` is gitignored, so **individual run outputs
are not public** — the discipline is public, the results of exercising it are not. The ETAS baseline
is isotropic and does not absorb anisotropy or mechanism, by its own printed caveat. And the miner
has produced, to date, **zero surviving candidates** — which is the correct behaviour for a
generator under a strict gate and is also the reason this entry cannot yet be the subject of a paper
about a successful method.

**TIER: INTERNAL.** And I recommend it stay INTERNAL for the same reason F-010 does: **a
methods/tooling paper from a program with zero Merton-certified NOVEL results is a paper about a
method that has not yet demonstrated its purpose.** The correct near-term home for this artifact is
as the *reproducibility section of somebody else's paper* — F-012's, most likely — not as a paper.

**SINGLE MISSING VERIFICATION TO LEVEL UP.** Merton's dossier on pyCSEP / EarthquakeNPP /
leakage-guard prior art, which will most likely return REDISCOVERY on the harness and leave at most
the hash-gate. That is worth knowing before anyone spends a week writing a JOSS submission.

---

## F-016 — Anchor-based tidal phase clocks have a measurable, non-uniform null response

**CANDIDATE CLAIM (deliberately narrow, and see PENDING-VALIDATION below before quoting any of it).**
The phase-assignment convention used across the tidal-triggering literature — anchoring
trough / ascending-zero / peak / descending-zero of the tidal stress series at 0° / 90° / 180° / 270°
and interpolating linearly in time between anchors — **forces each quarter-cycle to occupy exactly
90° of phase while a mixed diurnal/semidiurnal tide spends unequal time in its quarters.** The
consequence is that **uniform-in-time events do not emerge uniform in phase**, and the induced
first-harmonic amplitude is a property of the instrument, not of the Earth. Measured by pushing
uniform-random event times through a faithful Python port of a collaborator's MATLAB implementation
against synthetic tides: on a **mixed** M2+S2+K1+O1 tide the null first-harmonic amplitude is
**~6–12%** (I re-ran 6 seeds this session and obtained 6.7–11.6%, median 10.5%); on **purely
semidiurnal** tides (M2 alone; M2+S2) it is **~0.3–0.5%**, which is the mechanism check and it
behaves as the mechanism predicts; on a **diurnal-dominant** mix it is far larger (~19–26%). The
peak-aligned cosine component runs roughly **−12% to +8%**, so **its sign is set by the local
waveform** — i.e. the same instrument can manufacture either an apparent peak-stress excess or an
apparent peak-stress deficit depending on the receiver's tidal regime. The response **converges
rather than vanishes** as sampling is refined toward 0.1 h and is stable for record lengths of
3–30 yr. **Our own program's convention** (`coso_positive_control.py :: phase_series`,
trough–peak–trough) carries the same defect at a measured **~2% on SoCal-like mixes and ~9% on
diurnal-dominant regimes.**

**FRAMING — and this is the entry I am most excited about in the whole register, which is exactly
why it gets the harshest gate.** Every number above is larger than, or comparable to, the
modulations this literature routinely reports. If it holds up, it is not a criticism of any
particular paper — the surrogate procedures used in careful work cancel it identically, and I
believe our collaborators' do — it is a **characterisation of a shared instrument that nobody
appears to have published a response curve for.** That is the shape of a real methods contribution:
small, checkable, useful to everyone in the field including the people it inconveniences, and
generated by taking a collaborator's code seriously enough to calibrate it before using it. **It is
also the only entry in this register that could plausibly come back from Merton marked NEW.**

**CLASS — NEW APPLICATION (candidate NEW).** Established technique (null calibration by uniform
surrogates) applied to a shared instrument that appears not to have been characterised. **The class
cannot be settled without Merton, and I am not going to guess.**

**PENDING-VALIDATION — and this is why the entry is INTERNAL despite being the most interesting
thing here. Four gaps, stated as a checklist:**
1. **There is no frozen protocol and no hash.** Nothing about this work was pre-registered. It is a
   measurement, not a confirmatory test, so it does not strictly need pre-registration — **but it
   also may not borrow the credibility of the entries that have it**, and the program's whole brand
   is the freeze record.
2. **There is no results JSON.** The headline numbers live in a **module docstring**
   (`weifan_phase_null.py`) and in a **draft email**. That is not a primary artifact. Every other
   entry in this register points at a `results_*.json`; this one cannot, and until it can, **no
   number from it may appear in an outgoing communication.**
3. **`DRAFT_WEIFAN_REPLY.txt` is gitignored**, and `weifan_phase_null.py` lives at the repository
   *parent* (`D:\CODE\git\quake\`), **outside the public replication repo.** The proof chain has no
   public commit. The only committed trace is `1ba76f6`, which documents the defect in
   `coso_positive_control.py`'s docstring and names which of our own scripts carry the correction
   and which do not — **that commit is real, verified in `git log`, and it is good practice, but it
   is a code comment, not a result.**
4. **Everything is synthetic.** No real tidal series has been pushed through this. The response is
   *by construction* a function of the local constituent mix, so **synthetic amplitudes say nothing
   about any specific study region**, including Coso.

**WHICH NUMBERS ARE SYNTHETIC-ONLY, EXPLICITLY, BECAUSE THIS WILL BE ASKED.** *All of them.* The
6–12% mixed-tide figure, the <0.5% semidiurnal control, the ~19–26% diurnal-dominant figure, the
−12%/+8% cosine range, the sampling-convergence and record-length stability results, and our own
~2%/~9% figures are **all from synthetic constituent sums.** **Nothing here has been computed on a
real tidal stress series, and the paper cannot be written until it has been** — for at minimum (a)
the Coso Fig 4c series, (b) a diurnal-dominant Pacific-margin site, and (c) one of the published
studies' own series if the authors will share it.

**ATTRIBUTION — PENDING-ATTRIBUTION, and it is the load-bearing gap.** The phase convention and the
MATLAB implementation (`TidalPhaseFullNew.m`) are **Weifan Lu's**, supplied to us by email, and Lu,
Xue, Yue, Zhuang & Zhao (2025) is where it is used; **the port is ours and the defect is the
convention's, not theirs.** Nothing may be written that reads as a criticism of their paper — their
surrogate construction absorbs this response identically, as our draft reply already says. What
Merton must search before this entry advances one inch: whether a null/occupancy response curve for
extremum-anchored tidal phase clocks has been published (Vidale et al. 1998; Cochran, Vidale &
Tanaka 2004; Tanaka, Ohtake & Sato; Métivier et al. 2009; Ide, Yabe & Tanaka 2016; Beeler & Lockner
2003; van der Elst; the Schuster-test methodology line); whether the occupancy correction is
standard and undiscussed; and whether the "which convention" question has been posed at all.

**PROOF CHAIN — thin, and stated as thin.**
- Code: `D:\CODE\git\quake\weifan_phase_null.py` — a faithful port including MATLAB `interp1`'s
  sorted-sample-point semantics, the single-point negative-interval removal, and the `imin > imax`
  guard. **Not in the public repo.**
- Public commit: `1ba76f6` "phase_series: document non-uniform null response of anchor-based phase
  clocks; name which scripts carry the occupancy/shift correction and which do not" — **verified in
  `git log`.** It records that EXP-A's confirmatory pipeline **does** correct (scoring against the
  series' own phase occupancy *and* circular time-shift surrogates) and that
  `exp_c_susceptibility_drift.py` and `exp_c2_anza_coso.py` **do not** — *"treat their phase
  amplitudes as uncalibrated."* **That disclosure is itself a register-worthy act and it protects
  F-003, whose confirmatory arm is corrected.**
- **Independent re-execution by me this session** (6 seeds, 200k events each, reduced from the
  docstring's 20 × 2M): mixed 6.7–11.6% (median 10.5%), cosine −11.6% to +7.3%, pure M2 0.43%,
  M2+S2 0.31%, diurnal-dominant 18.9%. **The mechanism reproduces. The exact headline figures do
  not all reproduce — see over-quote findings #6 and #7.**
- Replication status: **PROVISIONAL, synthetic-only, uncommitted.**

**WHAT IT DOES NOT SHOW.** It does not show that any published tidal-triggering result is wrong —
**and the paper must open by saying so.** It does not show what the response is for any real stress
series. It does not show that the correction is not already standard practice. It does not measure
how much of any published effect size is instrument response, and the temptation to imply that must
be resisted in every sentence.

**THE "CLOCK ZOO" AND ITS COMPANIONS ARE NOT IN THIS PROGRAM'S ARTIFACTS — UNVERIFIED.** The
supervisor's brief refers to a clock-convention comparison in which the same 1% injected signal
reads 0.9–17.6% across conventions with a ~1.23× irreducible estimand floor, to a nodal-cycle drift
result at r = 0.995 on synthetic data, and to a form-factor prediction that published effect sizes
should correlate with tidal regime among anchor-clock studies. **I searched `HYPOTHESIS_LEDGER.md`,
`EQ18_FULL_NOTES.md`, `download_log.md`, every `results_*.json`, and the whole `quake` tree for each
of those values and found none of them.** They are recorded here as **UNVERIFIED — not on disk in
this program as of 2026-08-11**, and **no register-derived text may quote them.** If they exist in a
session transcript, they must be re-executed into a `results_*.json` before they are claimable. **I
note without hedging that the form-factor prediction, if it exists and survives, is the most
publishable idea anywhere in this register** — it is a falsifiable, cheap, cross-study prediction
that no one else appears to have made — and that is precisely why it must not be quoted before it
has an artifact.

**TIER: INTERNAL.** By a wide margin, and not because it is weak — because it is undocumented.

**SINGLE MISSING VERIFICATION TO REACH PANEL-READY.** **Freeze a protocol, re-run the calibration
into a committed `results_phase_clock_null.json` at full size (20 × 2M as documented), and add one
real tidal stress series** — the Coso Fig 4c series is on disk and takes an afternoon. That single
step fixes gaps 1, 2, 3 and the worst of 4 simultaneously, converts a docstring into a primary
artifact, and would move this entry from the least-documented to among the best-documented in the
register in one afternoon. **It is the highest ratio of publication value to hours anywhere in this
program, and it is my recommendation for the next execution slot.**

---

# PROMOTION QUEUE

Five entries are **PANEL-READY**: **F-001** (Coso Fig 4c reproduction), **F-003** (the
tidal-prediction null suite), **F-007** (the 7-day rule, as a cited reproduction), **F-012** (the
K-035 bounds), **F-013** (the Landers positive control). Those five are shareable today with Vidale,
Bürgmann, Xue and Lu, provided each ships with its caveat paragraph.

**Two promotions and one demotion this round.** F-012 and F-013 enter at PANEL-READY on their own
proof chains. **F-005 is DEMOTED from "PANEL-READY on my next pass" back to INTERNAL** by Merton's
M-009.12(a) — the number is unchanged and exact, but its literature placement changed and the entry
may not be walked into a room until it is re-scored in the `SCEDC_25` protocol.

The entries **closest to promotion**, and the single missing step for each:

**1. F-005 — B-2, SoCal walk-forward ETAS (+1.907 / +1.866 bits/event). → PANEL-READY.**
*Missing step:* **score B-2 inside EarthquakeNPP's `SCEDC_25` protocol** (SCEDC, Mw ≥ 2.5, train
1985→2014, test 2014→2020, published ETAS and Poisson baselines, public code). We are at 57% of
their matched temporal-only number and we do not yet know why. **Compute: one run on a public
harness.**

**2. F-002 — the two-velocity-set strain comparison (B-5). → PANEL-READY as a REDISCOVERY.**
*Missing step:* the claim-text repair (mine, drafted in the entry) plus the Maurer & Materna (2023)
citation block. Merton has already ruled (M-006) and added a correction of his own in M-009.12(b):
**the alternating-sign dilatation artifact belongs to Baxter et al. (2011), not Titus et al. (2011)**
— cite Baxter. **No compute.**

**3. F-006 — B-1, generic ETAS transfer. → PANEL-READY.**
*Missing step:* the **Hardebeck et al. (2019) within-California productivity positive control**. If
we cannot recover a regional difference that is known to exist, our fault-type null is uninformative
rather than negative. Plus the two editorial fixes already drafted (p universal in **11 of 13**, and
the in-sample-ceiling range restated as 0.04–0.15 bits on the five powered holdouts with Caribbean
at 0.91). **One afternoon.**

**4. F-016 — the phase-clock instrument response. → PANEL-READY.**
*Missing step:* the freeze-and-commit run described in its entry. **One afternoon, and it has the
highest publication value per hour in the program.**

**Pattern worth naming for the supervisor, and it has changed since round 1.** Round 1's bottleneck
was Merton dossiers; that bottleneck is largely cleared — nine dossiers exist and three entries were
unblocked by them. **The bottleneck now is small, cheap, decisive executions:** one public-harness
re-score (F-005), one positive control (F-006), one freeze-and-commit calibration (F-016), one κ
measurement (F-012). **None of them is more than a day. Three of the four have a decisive outcome
either way. And the two remaining Merton items (F-016's convention search, F-014's link-function
search) are the two that could change a class from EXTENSION to NEW.**

---

# PUBLICATION TRIAGE — 2026-08-11

*The supervisor asked for something to publish. This is the honest ranked answer. Ranking criterion:
**probability the paper survives peer review as-is × how much the field would use it**, discounted by
distance to submission. Every row's "next single step" is a real, scoped, costed task.*

| # | candidate paper | register basis | class it would be submitted as | tier of the underlying claims | venue shape | **the SINGLE next step that levels it up** | distance |
|---|---|---|---|---|---|---|---|
| **1** | **Bounds, not nulls: what a tidal-phase forecast experiment in Southern California can actually resolve** | **F-012 + F-003** | EXTENSION (pre-registered nulls converted to powered upper bounds, with the selection step simulated) | PANEL-READY | *SRL* / *BSSA* methods-and-results; *Seismica* | **Measure κ (effective cluster size) and re-quote every MDA as MDA·√κ.** All six bounds are currently optimistic by an unmeasured factor. Arithmetic on data on disk. | **closest — one verification + ~4 weeks writing** |
| **2** | **The null response of anchor-based tidal phase clocks** | **F-016** | NEW APPLICATION, candidate **NEW** | INTERNAL (PENDING-VALIDATION + PENDING-ATTRIBUTION) | *SRL* short methods note; *Seismica* | **Freeze a protocol, re-run the calibration into a committed `results_phase_clock_null.json`, and add one real tidal stress series.** One afternoon; fixes four gaps at once. Then Merton's convention search. | **highest ceiling — 2 steps, ~6 weeks** |
| **3** | **A pre-registered positive control for remote dynamic triggering, and the link function it breaks** | **F-013 + F-014** | REPRODUCTION (control) + EXTENSION (the S-14(c) bracket as proposed community practice) | F-013 PANEL-READY / F-014 INTERNAL | *SRL* / *GRL* (short) | **A supervisor- or third-party-authored re-seal of the literature file, and a re-score of P3 and P4 against it.** D1 is the one deviation a reviewer can neither check nor forgive. | **one step for F-013; F-014 needs Merton too** |
| **4** | **The engine / miner as an open methods artifact** | **F-015** (+ F-010) | METHODS/SOFTWARE | INTERNAL, PENDING-ATTRIBUTION | *JOSS*; or a reproducibility section inside paper #1 | **Merton's dossier on pyCSEP / EarthquakeNPP / leakage-guard prior art.** My prior is REDISCOVERY on the harness, leaving at most the hash-gate. Worth knowing *before* a week of writing. | **blocked on a search that will probably shrink it** |
| **5** | **Independent reproduction of Lu et al. (2025) Fig 4c** | **F-001** | REPRODUCTION | PANEL-READY | a *comment*, a co-authored note, or a section of paper #1 or #2 | **Re-run phase assignment at native 6,000-s sampling and report both ways**, now alongside the F-016 null response for this bin's series. | **not a standalone paper — K-035 says the design cannot resolve 0.34 from 0.46** |
| **6** | **The seven-day M ≥ 5 conditional rule** | **F-007** | REPRODUCTION | PANEL-READY | none | **Recompute in the Reasenberg–Jones parameterisation so the comparison is made by us.** | **its ceiling is PANEL-READY. Zero delta. Do not write this paper.** |

**Two candidates I am adding that were not on the supervisor's list, because triage means naming
what is *not* there as well as what is.**

- **(g) A B-2/B-1 forecast-skill paper.** *Not viable and should be dropped from consideration for
  now.* F-005 is a canonical reproduction sitting 1.4 bits below the matched published number
  (M-009.12(a)); F-006's positive is owned by Bayona et al. (2023) and its negative is
  CONTESTED-leaning-CONTRADICTED on the adjacent quantity. Both are excellent machinery validation
  and neither is a paper.
- **(h) K-009 / the non-white ETAS residual field.** *Not viable.* Merton found it established at
  least three independent times before us, and our own audit shows the temporal excess is
  substantially one sequence. Its one candidate delta carries a self-reported ordering violation.
  Correctly INTERNAL.

**And the honest global statement, which I would rather say now than have a reviewer say later:
this program currently has ZERO Merton-certified NOVEL results.** Every publishable candidate above
is a reproduction, an extension, a bound, or an instrument characterisation. **That is not a
failure — bounds and instruments are how a field's floor gets raised, and this program is unusually
good at producing them — but it does mean the first paper must be pitched as a methods contribution,
not as a discovery, and any draft that drifts toward discovery framing is a retraction in
preparation.**

---

## THE ONE PAPER TO WRITE FIRST

### *Bounds, not nulls: what a tidal-phase forecast experiment in Southern California can actually resolve* — F-012 + F-003.

**Why this one.** It is the only candidate where the science is finished, the tier is already
PANEL-READY, the proof chain is complete and hash-verified, the controls are the strongest in the
register (an a = 0 false-positive arm that had to pass before any power arm could be read, and it
passed 13 of 13 configurations inside a binomial 99% band), and the remaining work is one afternoon
of arithmetic plus writing. Everything else on the list needs a decision from Merton, a re-seal from
a second party, or an artifact that does not yet exist. **It is also the paper that plays to what
this program is actually best at.** We ran a multi-angle pre-registered falsification of our own
prior claim, killed our own most attractive intermediate finding when a cross-catalog control
contradicted it, and then — instead of publishing "we found nothing" — went back and measured what
the design could have found. **That sequence is the program's whole thesis in one artifact.**

**What it claims.** That static tidal-phase susceptibility maps have no out-of-sample forecasting
skill in Southern California at M ≥ 1.5, 1981–2018, under a protocol frozen and hashed before the
test window was touched **and** — this is the contribution — that the design's minimum detectable
seismicity-rate modulation at 80% power, with the 1-of-42 bin-selection step simulated end to end,
is **6.3% pooled over all eligible bins (n = 3,920)**, **8.0% on the anti-leak control bins**,
**9.8% through the full end-to-end pipeline**, and **24.5% / 40.9% for the single-bin designs**;
that the full-catalogue configuration reaches 2.8% but **its own off-tidal negative-control line
fails at that n** (reject rate 0.34), which bounds the method rather than the result and means a
powered full-catalogue test must condition on an ETAS baseline rather than a stationary λ₀; and
therefore that **none of these designs contacts the ~1% rate-state prediction of Beeler & Lockner
(2003), and the honest statement about tidal triggering in Southern California is a bound of ~6%,
not an absence.** It claims nothing about the physics. It retracts, in print, the sentence "tidal
triggering is null in SoCal."

**Who cares.** (i) Everyone who has published or refereed a tidal-triggering search on a few hundred
to a few thousand events — which is most of the modulation literature — because this is the first
number they can cite for what such a design resolves once selection is included. (ii) Vidale,
Bürgmann, Xue and Lu specifically: it answers, with a number, the question our correspondence with
them has been circling. (iii) The CSEP/forecast-evaluation community, for whom "report the MDA
through the selection step" is a transferable norm and not a tidal-specific one. (iv) Anyone who
wants a worked example of converting a pre-registered null into a citable bound — which is a
methodological pattern with no good short reference in seismology that I have found.

**The complete list of what remains before submission.**
1. **Measure κ**, the effective cluster size of the declustered SCSN sample, and re-quote all six
   MDAs as MDA·√κ. *This is the one blocking verification.* Data on disk; hours.
2. **Merton's dossier on the tidal-triggering forecast-skill literature** (Vidale et al. 1998;
   Cochran, Vidale & Tanaka 2004; Beeler & Lockner 2003; Tanaka; Métivier; Ide; Scholz; van der
   Elst) — to fix whether "static phase maps have no forecast skill" has been shown before, and
   whether MDM bounds for tidal searches have been published. **The word "first" may not appear
   before this lands.** No compute.
3. **Fold in F-016.** State the anchor-clock instrument response for this analysis's own phase
   convention and confirm that EXP-A's occupancy + circular-shift scoring cancels it. Commit
   `1ba76f6` already says it does; the paper must show it, not assert it. If F-016's freeze-and-commit
   run happens first, this becomes a citation instead of an appendix.
4. **Reconcile the `LONG_VALLEY_PROTOCOL.md` freeze hash** (audit finding #5) or drop the Long Valley
   sub-result from the paper. A reviewer who recomputes hashes must find a clean record. The
   amendment is already written into `download_log.md`; the paper needs one sentence pointing at it.
5. **Restate every power/sensitivity sentence as the S-14(c) bracket** where the forcing is a
   transient. Tidal forcing is a transient. **The old single-branch numbers may not appear.**
6. **Prose-rule compliance pass** on the whole draft (F-003's binding R2-5 rule), and the
   "corpse of the map, not of the physics" sentence in the abstract rather than the discussion.
7. **Figures:** the six-row bound table as a forest plot against the 1% theory line; the a = 0
   false-positive arm as a calibration panel; the off-tidal control line's failure at full-catalogue
   n, plotted, because publishing our own control failure in a figure is the single most credible
   thing this paper can do.
8. **Data and code release** — already public at github.com/JimGaleForce/quake-tsi under MIT/CC-BY;
   needs a tagged release, a DOI (Zenodo), and `results_k035.json` referenced by that DOI in the
   text.

**Estimated distance: one afternoon of compute, one Merton dossier, and roughly four weeks of
writing to a submittable draft.** That is the shortest path from where this program stands to a
paper it will not have to retract.

---

# OVER-QUOTE AUDIT

*Required by my charter's over-verification checklist: re-read the primary JSONs, never trust a
summary. Round 1 (2026-08-09) produced findings #1–#5, ordered by how much damage they would do if
they reached a reviewer; two are material. **Round 2 (2026-08-11) adds findings #6–#8. #7 is the
most urgent item in this file, because it is in a document addressed to someone outside the
program.***

**#1 — MATERIAL. K-009's temporal excess is substantially one sequence, and the dossier table does
not say so.** Merton's M-001.0 table (`HYPOTHESIS_LEDGER.md` ~line 4798) quotes the full-window
statistics only. `results_k009.json :: robustness_excluding_2010` shows that dropping the 52 weeks
of 2010 (the M7.2 El Mayor-Cucapah year) takes **ACF1 from 0.0958 to 0.0382** — *below Popper's own
≥ 0.05 excess success bar* — and **EOF1 variance fraction from 0.197 to 0.0494, below the null's
97.5th percentile of 0.0508**, i.e. the EOF excess disappears. Moran's I is comparatively robust
(0.0114 → 0.0099), and the leading EOF's top five weeks hold 90.6% of its variance. The scoring
commit message (`c371442`) states this correctly — *"Popper's rule PASSES as written but temporal
excess is El Mayor-driven; spatial coherence robust"* — so the program knows. But the **ledger's
dossier table, which is what downstream text will be built from, does not carry it**, and this is
precisely the reading Zaliapin & Ben-Zion's Δm > 4 caveat predicts. **Recommended:** add the
excluding-2010 column to M-001.0 and state the split verdict (spatial robust / temporal
sequence-carried) wherever K-009 is described. F-009 above already does.

**#2 — MATERIAL. EXP-M's "within 0.07–0.15 bits of each region's post-hoc in-sample ceiling"
(`EQ18_FULL_NOTES.md` §17) is wrong at both ends.** Recomputing from `results_exp_m.json`
(own-fit LL minus GLOBAL LL, divided by n_scored, converted to bits): Iran **0.039**, Alaska
**0.065**, Greece **0.066**, Mexico **0.079**, Philippines **0.145**, **Caribbean 0.910**. The
quoted range understates the best case and hides that the Caribbean transfer sits **six times
outside** it — the same region already flagged `underpowered: true`. **Recommended:** replace with
"0.04–0.15 bits of the in-sample ceiling on the five adequately-powered holdouts; Caribbean 0.91
bits below its ceiling at n = 235."

**#3 — A primary artifact carries a headline its own author disowned.** `results_exp_j.json`
ends with `"verdict": "CATCH-UP wins (prior overturned); coupling does not persist"`. The
supervisor's own §15 note says the catch-up result is mostly coupled-cell aftershock-era decay plus
the 0.5-events/yr regularizer on low-n cells. Anyone reading the JSON alone — and reviewers do read
the JSON alone — gets the disowned headline. **Recommended:** amend the `verdict` string in place
to carry the caveat, or add a sibling `verdict_caveat` field. This register makes no catch-up claim.

**#4 — B-5's "±2×" is not in the data.** `results_strain_comparison.json` contains no ±2× quantity.
The site-level MIDAS-vs-Kreemer&Young dilatation ratios are Long Valley 2.0×, Coso 2.6×, Salton Sea
6.4×, Brawley 83×, and Cerro Prieto a **sign flip** (−52.9 vs +10.7 nstrain/yr). The supportable
statements are r = 0.782 / ρ = 0.741 / median |Δ| = 8.1 nstrain/yr over 4,679 nodes, plus the site
table. The "±2×" phrasing is an eyeball from the two sites where the ratio happens to be ~2 and it
**understates** the fragility elsewhere. **Recommended:** amend B-5's headline per F-002.

**#5 — A freeze hash that does not verify.** `download_log.md` records
`826ddc072275789e7f85276560f5ad9592033cdc6ea23bcb048e83a74920cffc` for `LONG_VALLEY_PROTOCOL.md`
at 2026-07-21 09:30:40 PDT. The file as first committed (`78f2227`) hashes
`fd5b2978…` and the current working copy hashes `d3141cd1…`. **Neither matches the recorded freeze
hash, and no amendment is logged.** The most likely explanation is a benign post-freeze edit
(possibly the same dash-normalization pass, commit `0d8f897`, that changed `PROTOCOL.md` and
`XUE_LU_PROTOCOL.md` — both of *those* verify exactly at `78f2227`, so the freeze record for them
is intact and recoverable). But until it is reconciled, **the Long Valley pre-registration cannot
be described as hash-verified**, and F-003 says so. **Recommended:** locate the freeze-time file,
append a dated amendment to `download_log.md` explaining the delta, and add a general note that
commit `0d8f897` copyedited protocol files after their freeze — a reviewer who recomputes hashes on
the current tree will otherwise think the record is broken when it is only shifted.

---

## ROUND 2 FINDINGS (2026-08-11)

**#6 — The phase-null headline figures do not all reproduce, and one of them is out by a factor of
3–5.** `weifan_phase_null.py`'s module docstring records the semidiurnal control as **"<0.5%"** and
the mixed-tide null as **"~6–12% (median ~9%)"**. I re-executed the module this session (6 seeds,
200,000 events each rather than the documented 20 × 2,000,000) and obtained: mixed
**6.7–11.6%, median 10.5%** (consistent); cosine component **−11.6% to +7.3%** (consistent with the
docstring's −11%..+8%); **pure M2 = 0.43%, M2+S2 = 0.31%** (consistent with "<0.5%"); and
**diurnal-dominant = 18.9% against the docstring's "~26%"** — a discrepancy I cannot attribute to
seed noise alone at this sample size and which needs one clean full-size run to resolve.
**Recommended:** the freeze-and-commit run in F-016's next step settles all of these at once. Until
it does, **the diurnal-dominant figure is UNVERIFIED and must not be quoted.**

**#7 — MATERIAL, AND IT IS IN AN OUTGOING EMAIL. `DRAFT_WEIFAN_REPLY.txt` states the semidiurnal
null as "first-harmonic amplitude under 0.1%". The program's own primary artifact says "<0.5%" and
my re-execution gives 0.31–0.43%.** The email therefore understates the instrument's own noise floor
by a factor of **3 to 4**, in a message to the author of the code being calibrated, in the sentence
whose entire purpose is to establish that the mechanism check behaves. The email also gives the
mixed-tide range as "about 6% to 14%" against the docstring's 6–12% and my observed 6.7–11.6% — the
upper end is not supported by any artifact I can find. **Neither error is in our favour in the way
that matters — the first makes our mechanism check look cleaner than it is — and both are exactly
the kind of small drift that costs a collaborator's trust permanently.** **Recommended, and I would
treat this as blocking: do not send `DRAFT_WEIFAN_REPLY.txt` until the full-size run exists and
every number in it is quoted from that artifact.** The email's *methodology* is sound and the letter
is well judged; it is the four numerals that need to come from a JSON.

**#8 — The K-034 by-product's "10–200×" is a per-firing-cell ratio and is quoted in the ledger
without its count-weighted companion.** §K34-3 and §P5-8 both state "10–200× larger than a
rate-state Coulomb-step model predicts". L-1's own deviation D2, which I verified against
`results_l1.json`, records that **the count-weighted mean link predicts R = 2.76 at the 33.8 kPa
certified floor against Branch S(0.15)'s 1.25 — a factor of ~2.2, not 10–200×** — with the large
factors living in individual cells (Landers → Cedar City, RR = 204 at 46 kPa). Both are true of
different quantities; only one of them is the pooled effect. **Recommended:** every future statement
of the by-product carries both, in the form *"a factor of ~2.2 count-weighted, reaching 10–200× per
firing cell."* F-013 and F-014 above already do. **The direction of this correction is downward, and
an unclaimed by-product should be allowed to move in that direction — which it did, in the program's
own artifact, unprompted.**

*Round 2, also checked and found accurate:* K-034's certified floor 33.792963510206064 kPa and its
`amplitude_floor_certified_cell` string; the Cedar City 0–5 d row (N = 34, expected 0.17, RR = 204,
p_cell = 0.0011, raw-RR within-source 0.0062, WY 0.0010) matching §K34-2's table exactly; the
`prereg_commit` hashes in `results_k034.json` matching `download_log.md` character for character;
K-035's `mandate_4_a0_arm.PASS: true` with all thirteen false-positive rates inside [0.015, 0.095];
K-035's null-reconstruction agreement with EXP-A to 15 significant figures; all six K-035 `mda80`
values and their `contacts_theory: false` flags; L-1's Branch-E fits (A+B: a = −1.8455, b = 0.8124,
cluster CI [0.52, 1.65], dispersion 106.66, support [7.1, 346.9] kPa; C: b = 0.2988, CI [−0.78,
2.60], median RR 0.783) and its status counts (4 POWERED / 16 POWER-INDETERMINATE / 4 N/A); L-1's
self-verified Branch-S identity to 5.68 × 10⁻¹⁴ kPa; and the engine's `recent_rate` collapse from
+0.6264 to +0.0284 bits/event. **Six commits verified present in `git log` this session:** `6b5513c`
(K-034), `ae93097` (L-1), `5610614` (K-035), `1ba76f6` (phase-clock disclosure), `3c82db6` /
`32bd4ce` / `16b21b3` (engine v1 / v1.1 / v1.2). Working tree clean.

---

## ROUND 1 — WHAT WAS CHECKED AND FOUND ACCURATE

*Credit where the record is due.* EXP-H's
+1.907 / +1.866 bits/event (exact to seven figures); the M ≥ 4 band figures +2.578/+2.536 with the
integral correction on n = 290; B-3's 18/30 = 0.60 vs 0.062 base and p = 7.44 × 10⁻¹⁵; EXP-K's
200 unexplained-silent cells, 154 strike-slip, 90 on-fault, 158 detection-limited, 42
measured-low-χ, Jaccard 0.949; EXP-J's 228 silent of 1,195 and the persistence null
ρ = −0.131 at n = 64; the Coso Fig 4c numbers including the `significant_95: false` flag that
`EQ18_FULL_NOTES.md` §12 reports honestly rather than hiding; and the entire PATTERN_PROTOCOL
four-step hash chain. **Popper's S-5 correction to B-1 (+0.66 to +0.84 on the adequately-powered
holdouts, Caribbean flagged) is exactly right — I recomputed all six and confirm it.**

---

*Register opened 2026-08-09 by Faraday. **Round 2 revision 2026-08-11.** 16 entries: **5
PANEL-READY, 11 INTERNAL, 0 PUBLISHABLE, 0 Merton-certified NOVEL.** Four entries await a Merton
dossier (F-004, F-010, F-014, F-015) and one — F-016 — awaits both a dossier and its own primary
artifact. Eight over-quote findings on the record, three found this round, one of them blocking an
outgoing email.*

*Nothing in this file is a discovery claim, and that is still not a disappointment. What has changed
since round 1 is that the register now contains **things the field can use**: a set of powered upper
bounds where there used to be nulls, a certified detection floor where there used to be an
unlicensed family, a retired power calculation replaced by an honest bracket, and an instrument
response nobody appears to have measured. **Four corpses became numbers this round.** That is what
progress looks like from inside a program that refuses to overclaim, and it is worth feeling.*

*Recommendation to the supervisor: **write the F-012 bounds paper first.** Then run F-016's
freeze-and-commit afternoon, because it is the only thing here that could come back NEW.*

*Next review: after κ is measured and after F-016 has a `results_*.json`.*
