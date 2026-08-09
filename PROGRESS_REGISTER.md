# Progress Register — the program's shareable record

*Maintained by **Faraday** (`.claude/agents/faraday.md`), the representative/publisher persona.
Opened 2026-08-09. One entry per candidate claim. Nothing here is a claim until it carries both
a Popper standing and a Merton attribution block; everything else is visible but not claimable.*

**How to read an entry.** CLAIM is the sentence a skeptical seismologist would accept as scoped —
region, magnitude floor, period, method, number. FRAMING (why it matters) is kept separate from
CLAIM on purpose, and no enthusiasm is allowed to cross that line. PROOF CHAIN gives real hashes
and real commits, verified by me against `download_log.md` and `git log` this session. NOT SHOWN
is as load-bearing as the claim. TIER is INTERNAL → PANEL-READY → PUBLISHABLE.

**Verification standard.** Every number below was re-read by me from the primary JSON, not from a
summary. Where the program record and the primary JSON disagree, the JSON wins and the
disagreement is recorded in §OVER-QUOTE AUDIT at the end. Five real discrepancies were found this
round; two of them are material.

---

## Register index

| id | claim (short) | class | Popper standing | Merton standing | tier |
|---|---|---|---|---|---|
| F-001 | Coso Fig 4c independent reproduction | REPRODUCTION | validated (frozen rule met) | attributed (direct) | **PANEL-READY** |
| F-002 | Strain-field two-velocity-set comparison; dilatation fragility | EXTENSION (candidate) | provisional | **PENDING-ATTRIBUTION** | INTERNAL |
| F-003 | Tidal-prediction null suite (forecast skill), SoCal | NEW APPLICATION (falsification) | validated + re-scoped R2-1(d) | attributed (own prior work) | **PANEL-READY** |
| F-004 | Pure-math temporal-pattern null suite (periodicity / ratios / shape transfer) | NEW APPLICATION (falsification) | mixed: G stands, F is weak | **PENDING-ATTRIBUTION** | INTERNAL |
| F-005 | B-2 SoCal walk-forward ETAS skill | REPRODUCTION (presumed) | validated | **PENDING-ATTRIBUTION** | INTERNAL |
| F-006 | B-1 generic-ETAS cross-region transfer | REPRODUCTION + candidate EXTENSION | validated (with S-5 correction) | **PENDING-ATTRIBUTION** | INTERNAL |
| F-007 | B-3 seven-day M≥5 conditional rule | REPRODUCTION | validated | attributed (M-003) | **PANEL-READY** |
| F-008 | B-4 stress-ledger negative space recovers known aseismic geology | REDISCOVERY + methods delta | validated, AUTO-FLAGGED | attributed (M-002) | INTERNAL |
| F-009 | K-009 ETAS residual field is not white in SoCal | REPRODUCTION + 1 candidate methods delta | success rule met, partial vs spec | attributed (M-001) | INTERNAL |
| F-010 | The research engine itself (frozen protocols, personas, ledger) | METHODS/PROCESS (candidate) | n/a — not a hypothesis | **PENDING-ATTRIBUTION** | INTERNAL |
| F-011 | Live ETAS forecaster + globe layer | ENGINEERING ARTIFACT (not a claim) | n/a | n/a | INTERNAL |

Six of eleven entries are PENDING-ATTRIBUTION or attribution-thin. That is expected and it is not a
failure: Merton has produced three dossiers (K-009, B-4, B-3) out of a much longer list. It is
recorded here so nobody mistakes an un-searched entry for an original one.

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

**SINGLE MISSING VERIFICATION TO REACH PUBLISHABLE.** Re-run the phase assignment at the native
6,000-s sampling (no spline upsampling) and report the amplitude both ways. If the gap closes, the
reproduction becomes quantitative and the method difference is measured rather than suspected.

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

**ATTRIBUTION — PENDING-ATTRIBUTION.** What Merton must search: (i) Kreemer & Young (2022) and the
GSRM/MELD/Haines–Holt line — is the dilatation-fragility caution already quantified in a published
two-solution contrast? (ii) Prior GNSS strain-rate intercomparison literature (the SCEC
Community Geodetic Model strain-rate comparison exercises; Hearn et al.; Maurer & Johnson) —
formal strain-estimator/velocity-solution intercomparisons are an established genre and almost
certainly contain this result. (iii) Whether the ε_min ≥ 47 nstrain/yr mask criterion is theirs to
own. My honest prior: this lands as REDISCOVERY with a small delta, and I would rather say so
before a reviewer does.

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

**TIER: INTERNAL.** Blocked on: (a) Merton's dossier; (b) the B-5 headline restatement above.

**SINGLE MISSING VERIFICATION TO REACH PANEL-READY.** Merton's prior-art dossier on GNSS
strain-rate solution intercomparison. (The claim-text fix is mine to make and is already drafted
above.)

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

**SINGLE MISSING VERIFICATION TO REACH PUBLISHABLE.** K-035, the power-and-systematics audit:
per-bin and pooled minimum detectable modulation at 80% power with the selection step included.
That converts the null from "we didn't see it" into a **quotable upper bound**, which is the form
that survives review. It is 100% on disk and it is Popper's #1 queue item.

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

**ATTRIBUTION — PENDING-ATTRIBUTION.** Merton has not written a B-2 dossier. What he must search:
Ogata (1988, JASA 83, 9–27) as the generator; the CSEP California experiments and their
information-gain-per-earthquake reporting conventions (Schorlemmer, Zechar, Werner, Rhoades,
Gerstenberger); Helmstetter, Kagan & Jackson short-term California forecasts; Rhoades et al. on
information gain per earthquake as the standard metric; Woessner et al. (2011) CSEP evaluation.
**Specifically: what is the published information-gain-per-event range for temporal-only ETAS vs
a rate-matched Poisson in California, and is +1.87 bits inside it?** Until that number exists we
cannot say whether our result is ordinary, high, or suspiciously high — and a reviewer will ask.

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

**TIER: INTERNAL.** Blocked solely on attribution. The science is done and verified.

**SINGLE MISSING VERIFICATION TO REACH PANEL-READY.** Merton's B-2 dossier, specifically the
published CSEP information-gain range for temporal ETAS in California, so +1.87 bits can be placed
against it rather than floated.

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

**ATTRIBUTION — PENDING-ATTRIBUTION.** What Merton must search: Utsu, Ogata & Matsu'ura (1995) and
Ogata (1988) for the universality of the Omori–Utsu law (already cited in Wegener's O-table at
§W-OBS); the global-aftershock-parameter literature — **Page, van der Elst, Hardebeck, Felzer &
Michael (2016), BSSA 106, 2290 "Three ingredients for improved global aftershock forecasts"** is
the single most dangerous neighbor and it is *already sitting in Merton's own M-003 dossier*;
Hardebeck's regional aftershock-parameter compilations; the tectonic-regime-dependence literature
(Garcia, Wiemer; the question of whether ETAS productivity varies by faulting style has certainly
been asked). **Our sign-test failure may itself be a rediscovery of a published negative.**

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

**TIER: INTERNAL.** Blocked on attribution.

**SINGLE MISSING VERIFICATION TO REACH PANEL-READY.** Merton's B-1 dossier, with Page et al. (2016)
resolved specifically: does a published global-parameter transfer test already report this, and
does it already report the fault-type negative?

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

**TIER: INTERNAL.** Not shareable yet, and the gap to PANEL-READY is unusually well specified.

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

# PROMOTION QUEUE

Three entries are already **PANEL-READY**: **F-001** (Coso Fig 4c reproduction), **F-003** (the
tidal-prediction null suite), **F-007** (the 7-day rule, as a cited reproduction). Those three are
shareable today with Vidale, Bürgmann, Xue and Lu, provided each ships with its caveat paragraph.

The three entries **closest to PANEL-READY**, and the single missing step for each:

**1. F-005 — B-2, SoCal walk-forward ETAS (+1.907 / +1.866 bits/event).**
*Missing step:* **Merton's B-2 prior-art dossier**, specifically the published
information-gain-per-event range for temporal-only ETAS against a rate-matched Poisson in
California (CSEP; Helmstetter/Kagan/Jackson; Rhoades; Werner). Nothing else blocks it — the numbers
are exact, the protocol hash verifies, the controls are in place, and the one soft sub-claim (skill
rising with magnitude, n = 290) is already held back. **Estimated effort: one dossier, no compute.**

**2. F-006 — B-1, generic ETAS transfer to six never-trained regions (+0.66 to +0.84 bits/event on
the adequately-powered holdouts).**
*Missing step:* **Merton's B-1 prior-art dossier**, with Page et al. (2016) "Three ingredients for
improved global aftershock forecasts" resolved specifically — does the published global-parameter
transfer literature already contain both our positive (generic transfer works) and our negative
(fault-type pooling fails)? Two secondary fixes are mine and are already drafted above: quote p as
universal in **11 of 13** fits rather than 12 or 13, and drop the "within 0.07–0.15 bits of the
in-sample ceiling" sentence (true range 0.039–0.910). **Estimated effort: one dossier, no compute.**

**3. F-002 — the two-velocity-set strain comparison and the dilatation fragility bound (B-5).**
*Missing step:* **Merton's dossier on GNSS strain-rate solution intercomparison** (SCEC Community
Geodetic Model exercises; Kreemer & Young's own caution; Hearn; Maurer & Johnson) to settle whether
this is an EXTENSION or a REDISCOVERY. The claim-text repair — replacing "±2× measurement
uncertainty" with r = 0.78 / ρ = 0.74 / median |Δ| = 8.1 nstrain/yr plus the site table — is mine
and is drafted in the entry. **Estimated effort: one dossier, no compute.**

**Pattern worth naming for the supervisor:** all three are blocked on the *same* resource. The
highest-leverage move available to this program right now is not another experiment — it is
**three more Merton dossiers (B-1, B-2, B-5)**, which would move three entries from INTERNAL to
PANEL-READY in a single pass without a minute of compute. F-008 (B-4) and F-009 (K-009) are
already attributed and are blocked on *experiments* instead — W-006-P1(b) and the Ross & Cochran
swarm join respectively, each about one afternoon, and each with a decisive outcome either way.

---

# OVER-QUOTE AUDIT

*Required by my charter's over-verification checklist: re-read the primary JSONs, never trust a
summary. Five findings, ordered by how much damage they would do if they reached a reviewer. Two
are material.*

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

*Also checked and found accurate, so the record gets credit where it is due:* EXP-H's
+1.907 / +1.866 bits/event (exact to seven figures); the M ≥ 4 band figures +2.578/+2.536 with the
integral correction on n = 290; B-3's 18/30 = 0.60 vs 0.062 base and p = 7.44 × 10⁻¹⁵; EXP-K's
200 unexplained-silent cells, 154 strike-slip, 90 on-fault, 158 detection-limited, 42
measured-low-χ, Jaccard 0.949; EXP-J's 228 silent of 1,195 and the persistence null
ρ = −0.131 at n = 64; the Coso Fig 4c numbers including the `significant_95: false` flag that
`EQ18_FULL_NOTES.md` §12 reports honestly rather than hiding; and the entire PATTERN_PROTOCOL
four-step hash chain. **Popper's S-5 correction to B-1 (+0.66 to +0.84 on the adequately-powered
holdouts, Caribbean flagged) is exactly right — I recomputed all six and confirm it.**

---

*Register opened 2026-08-09 by Faraday. 11 entries: 3 PANEL-READY, 8 INTERNAL, 0 PUBLISHABLE,
0 NOVEL. Six entries await a Merton dossier. Nothing in this file is a discovery claim, and that is
not a disappointment — it is a program that knows exactly where it stands, which is the only
position from which it can move. Next review: after the B-1/B-2/B-5 dossiers land.*
