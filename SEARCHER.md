# SEARCHER — the miner's younger brother

**Status: PROPOSED (Kepler). Nothing here is licensed, nothing here is evidence, no claim is made
anywhere in this file. Popper prices and adjudicates; the supervisor runs frozen tests.**

*Seed — Jim Gale, 2026-08-13, compressed from his words and attributed to him in full:* "the agent
version of what I did with Sand Point. Look at a few earthquakes in a region, notice coincidences
across MANY properties at once — solid tide, aurora/geomagnetic storms, time of day, day of week,
anything. Then at scale, ML-statistically test ALL properties x ALL regions. Anything outside null
gets QUICKLY promoted to a properly tested theory, regionally and globally. Chase cross-region
variants — another region matching in a different quadrant. And honor small n with strong
concentration: in THIS odd region all 14 major quakes in 100 years are on a Monday — that is still
significant." He added, unprompted and correctly: *"I understand that often most ARE outside of
null."*

That last sentence is the whole design problem, and Jim named it before I did. This file exists to
make the SEARCHER a **generator with an accountant attached**, not a p-hacking machine with a
dashboard.

---

## §S0. WHAT THE SEARCHER IS — and the single sentence that keeps it honest

The MINER (`engine/mine.py`, `mine_session.py`) asks: *does this declared feature predict the
ETAS-residualised daily count?* It is a regression instrument on a target, priced test by test, and
it is deliberately narrow.

The SEARCHER asks a different question, the one Jim actually performed with his eyes:
**given the events themselves, is there any property of the world on which the events of some region
are CONCENTRATED?** No target. No regression. Unbinned, per-event, circular-and-categorical. It is a
**concentration detector over a (region x property x magnitude-stratum) lattice**, and it is the
formalisation of noticing.

> **THE SENTENCE. The SEARCHER never produces evidence. It produces a RANKING and a
> CONTROL-CALIBRATED SURVIVOR COUNT, and its only licensed output is a mechanically generated FREEZE
> FILE. Every claim this program ever makes from a SEARCHER hit is made by the frozen test that the
> freeze file specifies, run afterwards, on data the scan did not touch.**

This is not a new discipline. It is exactly the shape Popper granted C1 at **§P7-4** under Rule 4.4
(EXPLORATORY-UNPRICED, five conditions), with condition 5 — *"the follow-up rule is declared BEFORE
the run… without this condition, 'unpriced' is a selection stage with no accounting, which is how
p-hacking enters through a door marked exploratory"* — promoted from a condition to the SEARCHER's
architecture. The SEARCHER **is** condition 5, industrialised.

**One difference from C1 that must be stated plainly, because it is a relaxation and I am asking for
it explicitly.** §P7-4 condition 1 said the only outputs are two integers, no feature names, no
ranked list. **The SEARCHER cannot work under that condition** — a scan whose output is an integer
cannot promote anything. What I propose instead, and what Popper must rule on:

> **The ranked list may leave the run IF AND ONLY IF the ranking rule, the promotion threshold, the
> freeze-file template, and the confirmation design are all declared and hashed BEFORE the scan
> runs, and the ONLY thing a scan hit is permitted to do is instantiate that template.** The scan
> chooses *which* hypothesis gets frozen. It never chooses *how* it gets tested, *what* threshold it
> faces, or *what* would count as a pass. Those are fixed in advance, identically for every
> candidate. A selection stage with a pre-declared, candidate-independent follow-up is not a forking
> path; it is a sampling rule.

If Popper does not grant that, the SEARCHER is dead as designed and should be told so now rather
than after it is built.

---

## §S1. STAGE 0 — REGIONALIZE

**Deliverable: a declared region lattice with a rule id, frozen before Stage 2, derived from
published classification and geometry only.**

### S1.1 The rule

Extend `engine/regions_d.py`'s existing `D6-tectonic-bird2003-v1` — which already satisfies
§P7-22(b)(i) (*"defined from published tectonic classification, not from catalogue data"*) — to a
full lattice. Proposed rule id: **`S1-lattice-v1`**, composed of exactly three declared layers:

| layer | definition | source | count (proposed) |
|---|---|---|---|
| **L1 — named arcs/provinces** | The 13 boxes on disk (`data/comcat_world/*.csv`) plus KURIL and CASCADIA from `regions_d.py`, each carrying Bird (2003) class SUB/OTF/CRB/CCB. | published boundary traces | 15 |
| **L2 — tectonic classes** | The union of L1 boxes sharing a Bird class. Overlaps L1 by construction. | Bird (2003) | 4 |
| **L3 — declared grid** | Global 10° x 10° cells with N ≥ 200 events at the stratum's Mc. Overlaps L1 and L2 by construction. | geometry, then a **count threshold — see the hazard below** | ~40–70 |

**Overlap is ALLOWED and is the point** — Jim's "another region matching in a different quadrant"
requires nested and crossing regions. §P7-23's whole lesson is that a regional claim and a class
claim are different hypotheses, and a lattice that cannot express both cannot express Jim's question.

> **But overlap is NOT free, and here is the exact price and the exact hazard.**
>
> **(a) Overlapping regions are correlated tests, not independent ones.** An L1 hit inside an L2
> class will usually reproduce as an L2 hit from the same events. Naive Bonferroni over
> R = 15+4+60 = 79 regions therefore **over-corrects** (they are not independent) while naive BH
> **under-controls** (BH assumes independence or PRDS, and nested-set statistics are neither in
> general). **Disposition: the region axis is a declared STRATUM under §P6-3, not a flat multiplicity
> factor.** Strata = {L1, L2, L3}, budget identity `Σ m_s q_s = m q` enforced by
> `engine/strata.py:assert_budget_identity`, which already exists. An L2 hit whose covered events are
> ≥ 80% supplied by one L1 box is labelled **REDUNDANT-WITH-L1** and reported as one finding, not two.
> **(b) The L3 count threshold touches catalogue data and therefore violates the spirit of
> §P7-22(b)(i).** Honest fix: the threshold is applied **on the exploration window only**, is declared
> before the scan, and the surviving cell list is frozen and hash-affecting. It is a *measurability*
> filter, not a *selection* filter — no cell is dropped for what its statistic looked like, only for
> whether a statistic could exist. Rule 4.1 respected: nothing from the holdout window enters.
> **(c) Alaska-Aleutians is EXCLUDED BY NAME from every tidal-phase-family scan**, per §P7-22
> Ratification 5 and `regions_d.EXCLUDED_BY_NAME`. It may appear as an L1 region for property families
> unrelated to K-092 (e.g. day-of-week), and the exclusion list is therefore **per-property-family**,
> declared, and checked by `assert_alaska_excluded` at scan entry. Any region that has ever seeded a
> hypothesis is excluded by name from that hypothesis family's scan, permanently. This list only grows.

---

## §S2. STAGE 1 — THE PROPERTY LATTICE

**Deliverable: a per-event property vector `P(e)` for every event in every region, built once,
cached, and hash-affecting.**

This is the object the program does not have. The miner works on **daily binned counts** with
covariate columns; the SEARCHER needs **per-event property values** — the `circstat_event.py` path,
which exists, but with nothing yet joined to it at scale.

### S2.1 The families, instantiated from MINING_CATALOG.md

| # | property family | example properties | type | status post-tranche-B |
|---|---|---|---|---|
| **A** | **Solid tide** (K-092 family) | site tide level, level percentile within local range, phase angle θ, quadrant, dθ/dt sign, M2/O1/K1 constituent phases | phase + level | **EXISTS.** `engine/sitetide.py` landed (`site_scalar_at`, `tanaka_phase`, `constituent_phase`, `constituent_spectrum`). Untracked at K-092 freeze time per the K092 dated correction; now available. |
| **B** | **Lunar/solar ephemeris** (Family 1) | lunar phase, lunar declination, perigee/apogee distance, draconic phase, solar declination, Sun–Moon elongation, syzygy proximity | phase | **EXISTS.** `engine/ephemeris.py` + `data/lunar_grid_1980_2027.npz`. Currently day-gridded (`day_grid(hour=12.0)`) — **needs a per-event evaluation path**, `julian_day_at` already supports it. |
| **C** | **Human clock** (Family 7) | day-of-week, UTC hour, local-solar hour, month, day-of-year, holiday flag | categorical + phase | **PARTIAL.** `marks_ext.py:utc_hour_phase`, `local_solar_hour_phase` exist. Day-of-week exists only as a **control feature** on the count path (`observer.py:obs_day_of_week_phase`). **Needs promotion to a scannable per-event property under F7 controls — see §S3.4.** |
| **D** | **Geomagnetic / space weather** (Family 3A, 3D) | Kp and Ap at event time, Dst, SYM-H, solar-wind speed and Bz, storm-phase label (quiet/initial/main/recovery), days-since-storm-onset, F10.7 | level + categorical | **DOES NOT EXIST.** Needs download: Kp/Ap (GFZ Potsdam, 3-hourly, 1932–present, small), Dst (WDC Kyoto), OMNI solar wind (NASA, 1963–present). All public, all scriptable, total < 500 MB. **This is the largest single build in the design and it is the one Jim named first after tide.** |
| **E** | **Earth rotation** (Family 3B) | LOD, polar motion x/y, Chandler-wobble phase, annual-wobble phase | level + phase | **DOES NOT EXIST.** IERS EOP C04 download, tiny. |
| **F** | **Season / annual** (Family 6) | day-of-year phase, local hydrological-season index | phase | **PARTIAL** — day-of-year is free; hydrological loading is a build. |
| **G** | **Event's own marks** (Family 4) | depth, magnitude, dt-to-prior-event, distance-to-prior, cluster-member flag, local b-value at event time, focal mechanism class where available | level + categorical | **EXISTS.** `marks_ext.py:build_marks`, `_nearest_prior`. **Carries the tranche-B `VIF_mark` warning** — `dt_prior` and `cluster_member` measured above the 4.575 pooled fallback. |
| **H** | **Clocks as properties** (Family 8) | position in natural time, ETAS-rescaled time, cumulative-moment coordinate | level | **PARTIAL.** `engine/clocks.py` has the machinery + `random_clock_control`; needs a per-event evaluation. **Priced as a new test per feature per §P7-5(2)** — clocks are never free. |

### S2.2 The three rules that govern the vector

1. **Every property is evaluated at the event's own timestamp and coordinates**, never at a daily
   bin centre. The K-092 lesson: a daily bin has exactly one hour phase, which is why
   `obs_utc_hour_phase` is zero on the count path by construction. **The SEARCHER lives entirely on
   the sub-daily path** and therefore inherits `observer.assert_subdaily_gate` as a hard entry gate.
2. **Level properties carry their DWELL-TIME density with them.** For any property derived from a
   quasi-sinusoid (tide level, lunar distance, Dst during a storm), the null density of the *level*
   is **not uniform**. §P7-23(C) verified it exactly: for `x = A·sin θ` under uniform phase the level
   density is arcsine, `f(x) = 1/(π√(1−x²))` — **0.318 at mid-level, 1.019 at |x| = 0.95, and exactly
   1/3 of all events in the lowest quarter of range.** Any level property in the lattice ships with
   its analytic or simulated dwell-time density attached, and **a level-concentration statistic that
   does not divide by it is not emitted at all.** This is the single most likely way the SEARCHER
   would manufacture Jim's own Sand Point observation out of nothing, and it is disarmed at the data
   layer rather than at the review layer.
3. **Provenance is a field, not a memory.** Every property carries `{source, version, url_or_module,
   convention, scalar_definition}`. §P7-23(D)'s ruling — the quadrant is convention-free but only
   *once the scalar is fixed* — generalises: **a property without its convention attached is an
   S-18 clause 1 defect waiting to happen**, and the SEARCHER's whole output is properties.

---

## §S3. STAGE 2 — THE COINCIDENCE SCAN (the noticing)

**Deliverable: for every cell of the (region x property x magnitude-stratum) lattice, one
concentration statistic, one null-calibrated p, and one control-calibrated flag. Priced
EXPLORATORY-UNPRICED under Rule 4.4 and §P7-4/C1.**

### S3.1 The statistics — one per property TYPE, declared, no alternatives

Per **S-9** (*"the protocol names one value for every construction choice, with no alternatives
run"*), the SEARCHER declares **exactly one primary statistic per property type**, fixed for all
cells:

| property type | primary statistic | why | secondary (reported, labelled) |
|---|---|---|---|
| **phase** (circular) | **Kuiper V**, event-path, via `circstat_event.event_kuiper_watson` | rotation-invariant (§P7-22(a) — exactly, not approximately, which is what survives an unknown ocean-loading or admittance offset); omnibus, so it detects arcs a sinusoid test misses; demonstrated 100% recovery on the three-arc construction | Watson U², R₁, second moment — all four already priced and built |
| **categorical, k levels** (day-of-week, storm phase, quadrant) | **exact multinomial concentration**: the pre-declared max-cell binomial tail, scored against the **catalogue's own empirical base rates**, not against 1/k | see §S6 — the difference between 1/7 and the observed weekday rate costs 4.5 orders of magnitude and it is not optional | χ² omnibus |
| **level** (Kp, depth, tide level) | **dwell-time-corrected concentration**: the event-level distribution against the property's own occupancy measure, Kuiper on the PIT-transformed level | the arcsine lesson (§P7-23(C)), applied as a transform rather than as a caveat | quantile-band binomials |

**The null is ALWAYS ETAS-simulated event times pushed through the IDENTICAL property code path.**
Never a permutation of the property. §P7-22(a) is explicit about why: a deterministic warp of the
property map (ocean loading, admittance, ephemeris approximation) is **common-mode to signal and
null and cancels exactly** only if the null traverses the same code. `circstat_event.event_omnibus`
already takes `phase_fn` and `null_times` in exactly this shape. **Full ETAS with triggering, not
background-only** (§P7-23(A) condition 3): aftershocks inherit their mainshock's property
neighbourhood, and a background-only null scores that inheritance as signal.

### S3.2 Magnitude strata

Declared set **{M ≥ 4.5, M ≥ 5.0, M ≥ 5.5, M ≥ 6.0}**, and per **§P7-5(4)** the headline effect is
the **minimum over the declared Mc set** — a worst-case reading, demote-only, priced 0. A survivor at
one Mc that fails at a lower one is **labelled**, never claimed, unless it is separately declared as a
magnitude-dependent hypothesis and separately priced (which multiplies the tranche by 4, and the
SEARCHER does not do it by default).

### S3.3 What actually leaves the run

Per §P7-4 conditions 1 and 4, as amended by §S0:

1. **The control-calibrated survivor count**: `(n_survivors_real, n_survivors_control)` where the
   control arm is the matched F9-20-style scan over synthetic properties with identical marginal
   structure. **Only the difference is interpretable.** A raw survivor count with no control arm may
   not be quoted, in any forum, ever.
2. **The ranked list**, capped at a **pre-declared K** (proposed K = 30), sorted by a **pre-declared
   ranking function**, for the sole purpose of instantiating freeze files.
3. **Zero stubs. Zero K-entries. Zero findings.** Per §P7-16(4), unpriced rows cannot produce
   findings, and per §P7-22 Ratification the same rows cannot drive a confirmatory statistic either.
4. **One `EXPLORE_COUNT.jsonl` line** with the scan's `n_declared_tests`, config hash, and session
   id. The line is append-only and is never retroactively lowered (§P7-16(7) precedent).

### S3.4 THE F7 OBSERVER GATE — MANDATORY, AND THIS IS WHERE THE SCAN WILL DIE FIRST

**Day-of-week and hour-of-day are exactly where catalogue-composition artifacts live, and this
program measured them THIS WEEK.** Tranche B's F7 result is on the record: *"F7 was built to detect
diurnal and weekly detection structure and it detected it, at the resolution floor… a measurement of
the observer, not of the Earth."* The global max-statistic p = 1.31e-4 in that run was driven
**entirely by two unpriced control rows**.

> **RULES, all four binding, proposed as non-negotiable conditions of the SEARCHER existing at all.**
>
> **(F7-a)** Any scan cell whose property is a **human-schedule property** (day-of-week, UTC hour,
> local hour, month-boundary, holiday) runs **only after** `observer.observer_features` has been
> computed for that region and stratum, and its result is reported **beside** the cell, always, in
> every output.
> **(F7-b)** The **control arm is scanned identically and its survivors are counted in the same
> multiplicity** — Kepler's own F9-20 correction to Popper, ratified at §P7-6. A control channel that
> faces a different threshold than the science channel is not a control.
> **(F7-c)** A human-schedule survivor at a magnitude stratum where `obs_diurnal_amplitude_b*` is
> **non-zero** is **presumed artifact** and is not eligible for promotion. The band where that
> amplitude reaches zero **is** the completeness magnitude, measured independently of the
> magnitude-frequency distribution — `observer.py` says so in its own docstring. Use it.
> **(F7-d)** **THE COMPLETENESS CARVE-OUT, and it is the one place a human-schedule property can be
> real.** At **M ≥ 6.0 globally, the catalogue is complete**: no M6 goes undetected on a Sunday.
> A weekday concentration at M ≥ 6.0 therefore **cannot be a detection artifact**, because there is
> no detection to miss. **This is the honest scope of the carve-out, and three residual channels
> survive it and must be named rather than waved past:** (i) **magnitude assignment** can still drift
> with analyst schedule near the M6 boundary, so an M ≥ 6.0 result must be re-run at M ≥ 6.5 as a
> demote-only audit; (ii) **origin-time** is waveform-derived and schedule-free, so the *time* is
> clean even where the *magnitude* is not — which is why the day-of-week label is safer than the
> magnitude threshold that selects it; (iii) **historical completeness** — "100 years" of M ≥ 6
> global completeness is a claim about the ISC-GEM era, roughly 1904 onward for M ≥ 7 and later for
> M ≥ 6, and the declared start year is a hash-affecting choice that must be set from published
> completeness studies, not from where our CSV happens to begin.
>
> **The clean consequence, stated because it is the design's sharpest single point: the SEARCHER's
> human-schedule arm is scientifically live ONLY at M ≥ 6.0 and above, and is a pure observer
> measurement below it. Both are worth running. Only one is worth promoting.**

### S3.5 Cross-property redundancy

Lunar phase, tidal level, and syzygy proximity are not independent properties; nor are Kp, Ap and
Dst. Before the scan, **`engine/marks_ext.py:redundancy_audit` is run over the full property matrix**
and properties with pairwise |r| > 0.9 (or circular equivalent) are **collapsed into one declared
representative**, chosen by a pre-declared rule (earliest in the family table). This reduces the
declared multiplicity honestly rather than inflating it for appearance.

---

## §S4. STAGE 3 — FAST PROMOTION (the "QUICKLY")

**Deliverable: for every candidate above the pre-declared promotion score, a mechanically generated,
hashed, committed FREEZE FILE. Nothing else. Promotion is not a claim; it is a commitment to be
tested.**

**Promotion is cheap because the template already exists.** `K092_FREEZE.md` is the exemplar and it
is the reason this stage is a code generator rather than a research project. Its structure, section
by section, becomes the template's slots:

| K092_FREEZE.md section | template slot | filled from |
|---|---|---|
| "The claim", with seed attribution | `claim` | region + property + concentration form, rendered from the cell |
| "Region", explicit lat/lon box | `region` | Stage 0 lattice entry, verbatim |
| "Seed-exclusion superset… every M ≥ 6.0 event in the box, by ComCat id, sha256'd" | `seed_exclusion` | **§P7-23(A) condition 1 — exclude a SUPERSET, not a list.** For a scan, the superset is *every event the scan statistic touched*, enumerated by id and hashed. "Scrolled past is seen" becomes "scanned is seen." |
| "D-12 scoring set — STRATUM-HELD-OUT label" | `within_region_arm` | the region's events outside the exclusion superset, full-ETAS null |
| "D-13 prospective log, priced 0" | `prospective_log` | threshold, horizon (3 y primary, 1 y and 2 y descriptive-only), pass/fail written before the first event |
| "scalar/property provenance… TRANSLATION STEP OWED" | `provenance` | the property's provenance field from §S2.2(3) |
| "Pre-scoring gates (order binding)" | `gates` | dwell-time control, observer control, event-path VIF, in that order |
| "No result is claimed here. This file is a commitment, not evidence." | `disclaimer` | verbatim, always |
| — *new slot, not in K-092* — | `cross_region_arm` | §S5 |

**And the K-092 dated correction is itself a template requirement.** That freeze shipped a
priority-fact sentence that its own verification output contradicted, because the pipeline did not
gate on the check. **The generator therefore emits its priority facts as ASSERTIONS THAT MUST PASS,
and refuses to write the file if any fails** — S-18 clause 1 closed in code rather than in prose.

### S4.1 The promotion score, declared BEFORE the scan

```
promote(cell) := (p_real ≤ τ_p)
             AND (p_control > τ_p)                       # F7-b: control must NOT fire
             AND (concentration_form ∈ DECLARED_FORMS)   # quadrant / single-cell / arc
             AND (n_events ≥ 8)                          # §S6's floor, see below
             AND (dwell_time_corrected == True)
             AND NOT (human_schedule AND Mc < 6.0)       # F7-d
             AND (region NOT IN seeded_regions[family])
```
with **τ_p, the declared-forms list, and the n floor all frozen in the config hash before the scan**.
Ranked by `-log10(p_real)` within stratum. Cap K = 30. **These constants are set once and are not
tuned after seeing a scan; changing any of them is a new declaration and a new `EXPLORE_COUNT` line.**

### S4.2 What promotion costs

**A freeze file is priced 0** and needs no licence — §P7-23(B): *"a hashed commitment makes no claim
until scored and therefore needs no licence to be committed; it is a denominator committed before
the numerator exists."* The **scoring** of it is priced, later, when Popper prices it. This is what
makes "QUICKLY" both possible and safe: **the fast thing is the commitment; the slow thing is the
verdict, and it stays slow.**

---

## §S5. STAGE 4 — VARIANT FAMILIES (cross-region pattern algebra)

**Jim's move:** region A concentrates in quadrant Q1, region B concentrates in quadrant Q3 — not the
same phase, but the same *kind* of structure. That is not two findings and it is not one failed
finding. It is **a third, different, and stronger hypothesis**:

> **"Regional coherence with region-specific phase": each region concentrates its events at SOME
> preferred value of property P, with the preferred value varying by region.**

This is **exactly the construction §P7-22(a) ruled measurable** — concentration without ψ. Kuiper V,
Watson U² and R₁ are rotation-invariant *by construction*, so a per-region concentration statistic is
**blind to which quadrant each region picked** and remains valid even though the phase offsets are
unknown, unreportable (K-090(c) stands), and region-specific. The regional-sum-of-concentrations arm
already priced at §P7-22 (the 12-test regsum) is the same instrument.

### S5.1 The rule that keeps this from being free looks

> **A variant family is ONE declared hypothesis with ONE freeze file, not R separate looks.** The
> family declares: the property, the concentration form, the region set (all of it, in advance,
> including regions expected to fail), and the pooling statistic. It is scored **once**, on the
> pooled rotation-invariant statistic across the whole declared region set. **Adding a region after
> seeing its quadrant is forbidden and is the exact defect the family exists to avoid.**

Two sub-forms, both declarable, priced separately:

| form | statistic | what a pass means |
|---|---|---|
| **V-coherent** | pooled Kuiper/Watson over regions, phase-free | regions concentrate; the *offsets* mean nothing and are not reported |
| **V-structured** | pooled concentration **plus** a declared prediction of *which* offset each region takes, from an independent physical variable (fault dip, latitude, ocean-loading admittance) | far stronger — this is the arm that would turn a pattern into a mechanism, and it is the one that can be got wrong in an interesting way |

**V-structured is the entry nobody asked for, and it is the one I would run.** V-coherent can pass
on a whole-Earth artifact affecting all regions equally. V-structured cannot: it requires the offsets
to be *predicted*, from a variable that no catalogue statistic supplied. **A pass on V-structured is
a mechanism; a pass on V-coherent is a phenomenon; a fail on V-structured with a pass on V-coherent
is the most informative outcome of the three** and says the effect is real and our physical model of
it is wrong.

---

## §S6. THE SMALL-n INSTRUMENT — working Jim's Monday honestly

**Jim's principle is correct and Popper has already ratified it at §P7-23:** *"10 events all in one
pre-declared quarter gives p = 9.5e-7; 20 gives 9.1e-13. He is right."* And the reason the program
initially got this wrong is worth restating, because the SEARCHER will meet it constantly:
**§P7-22 evaluated a strong claim with a weak claim's instrument.** VIF floors and minimum
detectable amplitudes answer *"what is the smallest effect this region can resolve"* — they say
nothing about whether a **large** effect is detectable at small n. A quarter-cycle concentration at
R₁ near 1 is decisive at n = 10.

### S6.1 The arithmetic, exactly, for 14-of-14 on a Monday

| quantity | value |
|---|---|
| P(all 14 on Monday \| uniform, p = 1/7) | **1.474e-12** |
| × 7, for "some weekday", not a pre-named one | **1.032e-11** |
| × multiplicity, R = 20 regions, P = 50 properties, S = 3 strata (m = 21,000) | **3.10e-08** |
| × multiplicity, R = 60, P = 200, S = 4 (m = 336,000) | **4.95e-07** |
| × multiplicity, R = 200, P = 500, S = 5 (m = 3,500,000) | **5.16e-06** |

**14-of-14 survives full Bonferroni pricing at every plausible lattice size, by four to five orders
of magnitude.** It would take m ≈ 3.4e10 declared tests to break it. Jim is right, and it is not
close.

### S6.2 The two corrections that make it honest — and one of them nearly kills it

**(1) The null is not 1/7.** It is the region-and-stratum's **own empirical day-of-week rate**,
which F7 measured to be non-uniform. Suppose analyst-schedule structure lifts the observed Monday
rate to 0.186 (a 30% excess, well within what F7-03 sits at):

| under an inflated null p = 0.186 | value |
|---|---|
| P(all 14) | **5.93e-11** |
| × m = 336,000 | **1.99e-05** |

**Still survives — but the artifact cost 4.5 orders of magnitude.** That single row is the argument
for F7 controls being mandatory rather than advisory: **the correction is not decorative; it is 45%
of the evidence, on a logarithmic scale.**

**(2) Near-total concentration is NOT total concentration, and the cliff is brutal.**

| | p | × m = 336,000 |
|---|---|---|
| 14 of 14 on Monday | 1.474e-12 | 4.95e-07 — **survives** |
| **10 of 14** on Monday | **2.034e-06** | **0.683 — FAILS, and not marginally** |

**This is the honest headline of the whole small-n instrument.** Jim's example works *because it is
perfect*. Drop four events and the same lattice, the same regions, the same properties reduce it to
noise. **The SEARCHER must therefore report the exact tail, never a "strong concentration" adjective,
and the n ≥ 8 floor in §S4.1 exists precisely so that a 3-of-3 coincidence — p = 1/49 before
multiplicity — never enters the ranked list at all.**

### S6.3 The dependence correction, which is the real risk at small n

Fourteen events in one region across a century are **not** 14 independent draws if any of them are
aftershocks of each other. §K92-0(3) makes the point for Alaska: five M ≥ 7.2 events in ~150 km since
2020 are **one interacting sequence**. **Disposition: the ETAS-simulated null handles this
automatically and analytically-computed tails do not.** So:

> **The reported p is ALWAYS the full-ETAS-null p. The analytic (1/7)ⁿ is a sanity check printed
> beside it and is never the headline** — §P7-23(A)'s own instruction: *"report the ETAS-null p, not
> the analytic (1/4)ⁿ; the analytic figure assumes independence and uniformity that the null
> measures directly."* At M ≥ 6 with 14 events across 100 years in one box, declustering will
> usually leave 14 ≈ 14 — but *"usually"* is not a number and the null is.

---

## §S7. THE TRANSFORMER ANALOGY — accepted as metaphor and as a ranker, refused as evidence

Jim's framing: attention over the property x region lattice, with promotion as the head.

**As an architecture metaphor it is apt and useful.** The scan computes a compatibility score over a
product space; the promotion stage selects a sparse top-K; the freeze file is the emitted token. It
names the shape correctly.

**As a future component it is genuinely valuable in exactly one place: SCAN ORDERING.** A learned
ranker — gradient-boosted or attention-based — trained on (cell features → did the frozen test later
pass) would let the SEARCHER spend its compute on the cells most likely to survive confirmation.
Given a lattice of 10⁵–10⁶ cells and an ETAS null that costs simulation per cell, **ordering is worth
real money**, and this is a legitimate ML application inside a Popperian program.

> **AND IT CHANGES NOTHING ABOUT WHAT COUNTS AS EVIDENCE. Stated flatly so it can never be
> misremembered: a learned ranker changes what gets TESTED FIRST. It never changes what counts as a
> PASS.** The evidential layer stays the ruled classical machinery — declared statistic, ETAS null,
> declared multiplicity, strata budget identity, frozen threshold, held-out confirmation. Three
> consequences, all binding:
>
> 1. **The ranker may never touch the threshold.** τ_p is frozen in the config hash and is identical
>    for every cell regardless of rank.
> 2. **The ranker is trained on CONFIRMED OUTCOMES only** — cells whose freeze files have been
>    scored — and therefore cannot exist until the SEARCHER has produced a scored corpus. It is a
>    **phase-2 component**, and proposing it now is proposing a dependency, not a feature.
> 3. **The ranker's training data is a holdout spend if it includes holdout-window outcomes.** Rule
>    4.1. The ranker trains on exploration-window confirmations, or it does not train.

**A learned ranker in a program with a frozen evidential layer is safe. In a program without one it
is the most efficient p-hacking engine ever built.** That is the entire content of this section.

---

## §S8. ENGINE REQUIREMENTS

### S8.1 What exists (verified by reading the modules, 2026-08-13)

| capability | module | note |
|---|---|---|
| event-path circular statistics | `engine/circstat_event.py` | `event_kuiper_watson`, `event_omnibus(theta_obs, null_times, phase_fn, …)`, `resample_event_times`, `event_r1_stat`, `event_second_moment_stat`. **Inherits the §P7-21(c) Kuiper plateau-crossing fix — binding.** |
| solid-earth site tide | `engine/sitetide.py` | `site_scalar_at`, `tanaka_phase`, `tidal_maxima`, `constituent_spectrum` |
| ephemeris | `engine/ephemeris.py` | `julian_day_at` supports arbitrary event times |
| observer controls | `engine/observer.py` | F7-01/02/03, Mc drift, network era, day-of-week, UTC hour; `assert_subdaily_gate` |
| declared strata + budget identity | `engine/strata.py` | `stratified_bh`, `assert_budget_identity`, `max_statistic_p`, `max_statistic_report` with covered-set reporting |
| region lattice primitives | `engine/regions_d.py` | rule id, `assert_partition_rule`, `assert_alaska_excluded`, `declaration_digest` |
| clocks + random-clock control | `engine/clocks.py` | `assert_random_clock_control` |
| per-event marks | `engine/marks_ext.py` | `build_marks`, `_nearest_prior`, `redundancy_audit` |
| declaration logging | `engine/EXPLORE_COUNT.jsonl`, `HOLDOUT_LOG.jsonl` | append-only |
| freeze template exemplar | `K092_FREEZE.md` + `K092_seed_exclusion_superset.csv` | the thing that makes Stage 3 cheap |

### S8.2 What must be built

| # | build | size | depends on |
|---|---|---|---|
| **B1** | **`engine/properties.py` — the per-event property join.** Takes `(events, property_family)` → `(N_events, N_properties)` matrix with provenance and dwell-time density attached per column. **This is the missing object.** | medium, ~400 lines | S2 families A/B/C/F/G exist; D/E need B4 |
| **B2** | **`engine/searcher.py` — the region-lattice runner.** Iterates the declared lattice, computes the declared statistic per cell, runs the matched control arm, emits `(n_real, n_control)` and the capped ranked list, writes one `EXPLORE_COUNT` line. | medium | B1, `circstat_event`, `strata` |
| **B3** | **`engine/freeze_gen.py` — the promotion generator.** Renders the K-092 template from a cell; enumerates and sha256s the seed-exclusion superset; **asserts every priority fact and refuses to write on failure**; emits `K###_FREEZE.md` + `.csv`. | small, ~250 lines | `K092_FREEZE.md` as template |
| **B4** | **Geomagnetic / space-weather download** (Kp, Ap, Dst, OMNI) + **IERS EOP**. All public, scriptable, small. | small-medium | none |
| **B5** | **`engine/cli.py search` mode**, alongside `explore` / `holdout` / `mine`, with the same GENERATOR_NOT_EVIDENCE banner `mine.py:84` already prints — **and a second banner naming §S0's sentence.** | small | B2 |
| **B6** | Region lattice L3 grid-cell definition + declaration digest, extending `regions_d.declaration_block()`. | small | none |

**Sequencing: B1 → B2 → B3 gets a complete SEARCHER over the properties that already exist (tide,
ephemeris, human clock, marks, season). B4 adds the family Jim named second (geomagnetic) and is
parallelisable.** Nothing here waits on tranche D.

---

## §S9. TRANCHE SKETCH, FOR POPPER TO PRICE

**Proposed as TRANCHE S. Counts are rough and are offered for pricing, not asserted.**

| arm | what | proposed class | rough count |
|---|---|---|---|
| **S-0** | Region lattice `S1-lattice-v1` declared, frozen, hash-affecting | GATE | 0 |
| **S-1** | Property lattice B1 built; provenance + dwell-time density per column; `redundancy_audit` run and collapse applied | GATE | 0 |
| **S-2** | **Instrument recovery on synthetic**: planted quadrant concentration, planted single-cell categorical, planted arc — recovery curves per statistic per n ∈ {8, 14, 30, 100}; plus **zero-recovery on the negative controls** (day-lattice, random clock) | GATE, and mandatory before any real scan — §P7-6's "demonstrated recovery" demand | 0 |
| **S-3** | **F7 observer baseline** per region per stratum, reported and stored | GATE | 0 |
| **S-4** | **THE SCAN**, real arm | **EXPLORATORY-UNPRICED** (Rule 4.4, §P7-4 five conditions as amended by §S0) | ~10⁴–10⁵ cells, **0 priced**; declared in `EXPLORE_COUNT` |
| **S-5** | **THE SCAN**, matched control arm, identical threshold | **EXPLORATORY-UNPRICED**, counted in the same multiplicity | same count, 0 priced |
| **S-6** | **Freeze-file generation** for up to K = 30 promoted candidates | 0 — a commitment is not a claim (§P7-23(B)) | 0 |
| **S-7** | **Within-region held-out scoring** of promoted candidates, STRATUM-HELD-OUT label | **PRICED** | ~2 per candidate → **up to 60** |
| **S-8** | **Cross-region confirmation** of promoted candidates, declared region set, Alaska-class exclusions by name | **PRICED** | ~4–7 per candidate; realistically only the top 3–5 candidates go here → **~20–35** |
| **S-9** | **Variant families** (V-coherent and V-structured), one declared hypothesis each | **PRICED** | 2 per family, ~3 families → **6** |
| **S-10** | **Prospective logs** for every promoted candidate, hashed, 3-year horizon | 0 (K-080/D-13 pattern) | 0 |
| | **TOTAL PRICED** | | **~90–100**, plus one unpriced scan pair |

**The shape is deliberate: the scan is free and enormous; the confirmation is priced and small; the
prospective arm is free and is the one that actually settles things.** That is the same shape §P7-23
found for K-092 and it is not a coincidence — it is the only shape a scan can have in this program.

---

## §S10. HOW THIS DESIGN FAILS, WRITTEN BEFORE IT RUNS

Four ways, named in advance so that meeting one is a result rather than a surprise:

1. **The control arm matches the real arm.** `n_control ≈ n_real`. This is the **expected** outcome
   and Jim predicted it (*"often most ARE outside of null"*). It is a real and publishable
   measurement: *the rate at which a property x region lattice produces apparent concentration under
   a matched null* is a number nobody in this field has published, and it is the number every future
   coincidence claim should be measured against. **The SEARCHER's most likely contribution is a
   calibration constant, not a discovery.** I would take that trade.
2. **Every survivor is a human-schedule property below M6.** Then the SEARCHER has re-measured F7 at
   greater expense, and F7-c/F7-d catch it before promotion. Cost: the build. Value: a
   region-resolved observer map, which the program does not have.
3. **The dwell-time correction eats everything.** Every level-property survivor turns out to be
   arcsine. Then §P7-23(C) generalises from tide to the whole property lattice and that is worth
   knowing precisely.
4. **Promoted candidates fail confirmation at exactly the base rate.** Then the promotion score is
   uninformative and the ranker of §S7 has nothing to learn — which is itself the measurement that
   says scan-based discovery does not work in this domain, and it is the strongest negative result
   available here.

**None of those four is a wasted tranche. Three of them are publishable as bounds, which is the
property §P7-4 said Tranche C had and is the reason I am proposing this at all.**

---

## AUDIT LINE

*Kepler, 2026-08-13. Seed: Jim Gale, attributed in full at the head of this file — the SEARCHER is
his construct, formalised, not mine. This file is NEW: created at
`D:\CODE\git\quake\replication\SEARCHER.md`, zero existing lines modified anywhere. One K-entry
(K-094) appended to `HYPOTHESIS_LEDGER.md` under a new `# PROPOSED (Kepler)` heading; append-only,
zero existing ledger lines modified. Nothing committed. No claim is made in this file. Every stage
above is grounded in an existing ruling — §P7-4/C1 (five conditions, Rule 4.4), §P7-5(2)/(3)/(4),
§P7-6, §P7-16(4)/(7), §P7-19(c), §P7-22(a)/(b)/Q1/Q4/Q5 and its five ratifications, §P7-23(A)/(B)/
(C)/(D), §P6-3, §P6-5, Rule 4.1, S-8, S-9, S-13, S-15(c), S-17, S-18 — and the one place I am asking
for a RELAXATION of an existing ruling (§P7-4 condition 1, the ranked-list output) is flagged as such
in §S0 rather than assumed. Statistics reported in §S6 are arithmetic on declared constants, computed
this session, not measurements on any catalogue; no phase, concentration, or day-of-week statistic
was run on any real catalogue in the writing of this file.*

---

## §S11. RECONCILIATION WITH §P7-24 (appended after the fact; nothing above modified)

**§P7-24 — "THE SEARCHER — STANDING PROTOCOL (SP-1..SP-8)" — was appended to `HYPOTHESIS_LEDGER.md`
by the Popper seat during the same round in which this file was written, and I read it only after
this file and the K-094 entry were already on disk.** It adjudicates this design directly. Recorded
here rather than edited in, per the program's append-only convention; **where §P7-24 and §S0–§S10
differ, §P7-24 governs.**

**What §P7-24 grants that §S0 asked for.** SP-1 condition 2 makes the **ranked list** publishable
output alongside the control-calibrated survivor count. **The relaxation of §P7-4 condition 1 that
§S0 flagged as make-or-break is GRANTED.** The design lives.

**Where §P7-24 is stricter than this file, and is right.**

1. **SP-3 pins the threshold this file left declarable.** §S4.1 wrote `τ_p` as a frozen-but-free
   constant. SP-3 fixes it at **`alpha = q/m`, q = 0.10, m = the scan's FULL declared cell count**,
   and **refuses an OR-limb**. That is stricter and it removes the one dial a future operator could
   have turned. Adopt as written; `τ_p` in §S4.1 is superseded.
2. **SP-7 adds a gate I did not have, and it is the one that matters.** §S9's S-2 gate demanded
   recovery on planted signals. SP-7 demands the complement: **the full searcher run over ≥ 30
   true-null ETAS-simulated catalogues, expecting ≤ 0.1 promotions each.** A recovery gate proves the
   instrument can see; SP-7 proves it does not hallucinate — and it catches an invalid SP-2 null
   layer, which no multiplicity arithmetic can. **No real scan before it passes.**
3. **SP-6.3 — "the seeding scan travels with the claim"** (scan id, date, `m`, and the candidate's
   rank inside the freeze file, forever). Absent from §S4's template table. **Add it as a required
   slot.** It lets any reader recompute the look-elsewhere effect without asking us.
4. **SP-2 organises the null-validity layer by PROPERTY CLASS, built once and reused forever.** This
   file scattered the same requirements across §S2.2(2), §S3.1 and §S3.4. SP-2's table is the better
   object: it makes "a property whose class has no null layer cannot promote" a checkable precondition
   rather than a review-time judgement. **`engine/properties.py` (B1) should carry `property_class`
   as a required column and refuse to emit a property without one.**
5. **SP-6.7 — magnitude strata enumerated in the scan declaration**, and a claim at an unenumerated
   threshold is a **new seed, not a result.** §S3.2 declared the set but did not say what happens off
   it. SP-6.7 does.

**Where §P7-24 supersedes this file's arithmetic.** §S9's priced sketch (~90–100) is **wrong in
structure, not just in number**: SP-4 sets a **standing price of 16 per promoted claim** (2
within-region + 14 cross-region + 0 prospective), on the principle that *promotion freezes the
property and the statistic, so the scan's breadth is paid for once in SP-3's threshold and never
again in the confirmation.* **Use SP-4's 16, not §S9's per-candidate estimate.**

**The one thing in this file that §P7-24 does not contain, offered as a complement to SP-8.** SP-8
works Jim's Monday to a promotion. §S6.2 works the **cliff beside it**: on the same lattice,
**10-of-14 on Monday gives 2.03e-06, which fails m = 336,000 at 0.68** — while 14-of-14 clears by
four to five orders. The two together say the instrument is a detector of near-total concentration at
small n and of nothing weaker, which is why §S4.1's **n ≥ 8 floor** and §S6.3's insistence on the
**full-ETAS-null p over the analytic tail** are load-bearing rather than decorative.

*Kepler, 2026-08-13, appended after reading §P7-24. Nothing above §S11 modified; the original audit
line stands as written. §P7-24 governs wherever it and this file differ.*
