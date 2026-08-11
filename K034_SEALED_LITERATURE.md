# K-034 — SEALED LITERATURE VALUES AND FROZEN AMPLITUDE MODEL

Written and hashed **before** any K-034 statistic was computed, per Popper R2-2 mandate (1)
("Seal the literature values ... the supervisor writes the published Landers triggering
distances/amplitudes to a hashed file before the fit; the comparison is scored against
this file").

**DEVIATION, FLAGGED:** R2-2 assigns the authoring of this file to *the supervisor*, not the
analyst. This program is being executed by a single worker agent; there is no second party
available to author the seal. The seal is therefore **analyst-authored**, and its protection is
weaker than R2-2 intends: it prevents post-hoc *editing* of the comparison target (hash-pinned in
`download_log.md` before the fit) but not analyst foreknowledge. Recorded so the licence carries
this qualification.

## §1. Verification status of every value below

Each line is tagged:
- **VERIFIED** — retrieved from a source in this session, source named.
- **PARTIALLY VERIFIED** — recalled/secondary, consistent with a source seen this session but
  the primary number was not read. May not be quoted in a protocol without a primary read (S-14
  precedent on Merton's Aσ figures).

## §2. The triggering claim being scored against (Hill et al. 1993, Science 260, 1617-1623)

- **VERIFIED** (science.org abstract + search retrieval, this session): "The magnitude 7.3 Landers
  earthquake of 28 June 1992 triggered a remarkably sudden and widespread increase in earthquake
  activity across much of the western United States. The triggered earthquakes ... occurred at
  distances up to **1250 kilometers (17 source dimensions)** from the Landers mainshock, were
  **confined to areas of persistent seismicity and strike-slip to normal faulting**. Many of the
  triggered areas also are **sites of geothermal and recent volcanic activity**. Static stress
  changes ... appear to be too small to have caused the triggering."
- **VERIFIED** (search retrieval, this session): **Little Skull Mountain, NV — 280 km from
  Landers, M5.6, 22.3 hours after Landers.**
- **VERIFIED** (search retrieval, this session): **Mina, NV — 500 km, M4.0, 36 minutes after.**
- **VERIFIED** (search retrieval, this session): **Smith Valley, NV — 590 km, M3.4, 56 minutes after.**
- **VERIFIED** (search retrieval, this session): triggering confirmed at **Long Valley caldera**,
  **The Geysers geothermal field**, and **Yellowstone**.
- **PARTIALLY VERIFIED**: the Hill et al. (1993) site set also includes **Coso**, **Cedar City,
  Utah**, **Lassen Peak / Burney, CA**, **western Idaho**, and the **Nevada Test Site region**.

**The scored predictions (committed here, before the fit):**

- **P1 (detection).** The dynamic-stress / post-window rate statistic fires (family-wise
  significant) in at least one documented site for Landers.
- **P2 (n>1 control, R2-2 mandate 2).** The same engine, unchanged, fires on **at least 2 of
  {Landers 1992, Hector Mine 1999, Ridgecrest 2019, Denali 2002}**. This is the PASS rule.
- **P3 (spatial pattern, R2-2 mandate 3).** Geothermal/volcanic cells rank above non-geothermal
  cells of comparable predicted amplitude. Scored against the ranked list in
  `K034_PREREGISTERED_CELLS.md`, by rank correlation, not by eyeball.
- **P4 (amplitude).** The implied triggering threshold in peak dynamic stress lies in the
  literature band **0.01-0.1 MPa** (Kepler's entry text) / R2-2's certified band **10-100 kPa**.
- **P5 (distance).** Detections occur beyond 2 rupture lengths and out to ~10^3 km.

## §3. Frozen amplitude (dynamic stress) model — the covariate, not the response

**Primary, far field — van der Elst & Brodsky (2010), JGR 115 B07311, their eq. (6), read
verbatim from the paper PDF in this session (VERIFIED):**

```
log10 A20 = Ms - 1.66 * log10(D) - 2
```
where `A20` is in micrometres and `D` is epicentral distance in **degrees**; the paper states
this is the surface-wave magnitude relation of Lay & Wallace (1995), used in reverse (catalog
magnitude -> amplitude), taking the T = 20 s long-period waves as the indicator of peak dynamic
strain.

**VERIFIED, same paper, same paragraph:** displacement is converted to velocity by
`V ~= 2*pi*A20 / T` with `T = 20 s` (Aki & Richards, 2002).

**VERIFIED, same paper:** dynamic stress is obtained from dynamic strain with a **crustal shear
modulus of 30 GPa** ("For a crustal shear modulus of 30 GPa, this corresponds to a dynamic stress
of 0.1 kPa").

Strain is taken as `eps = V / c` with **phase velocity c = 3.5 km/s** frozen (20 s Rayleigh
phase velocity, continental crust). Hence

```
sigma_dyn [Pa] = 30e9 * V[m/s] / 3500
```

**VERIFIED caveat, same paper:** "The surface wave magnitude equation is designed for distances
on the order of at least 800 km." Most K-034 target cells are 150-900 km from their source, i.e.
**inside the extrapolation zone**. This is a declared limitation of the amplitude axis, not of the
detection statistic, and it is why a second amplitude axis is carried:

**Secondary / near-to-regional-field bracket — same paper, their eq. (5) with the unconstrained
Table 1 constants, read verbatim in this session (VERIFIED):**

```
log10 PGV[cm/s] = c1 + c2*M - c3*log10( sqrt(r^2 + c4^2) ),  c1=-2.29, c2=0.85, c3=1.29, c4=0
```
with `r` = hypocentral distance in km.

Both axes are reported for every cell. The certified amplitude floor is reported as the interval
spanned by the two axes; no single-axis floor is claimed.

**Sanity anchor, computed here before the fit and recorded so the arithmetic is auditable:**
Landers Ms 7.3 at D = 415 km (3.732 deg) gives, on the primary axis, A20 = 2.24e4 um,
V = 7.04e-3 m/s, eps = 2.0e-6, sigma_dyn = 60 kPa. This lands inside the sealed 10-100 kPa band
without any tuning to it. On the secondary axis the same source/receiver gives ~296 kPa.

## §4. Rate-state response and the S-14 bracket

Per K-036 / S-14, the response to a stress perturbation is `R = exp(dTau / A*sigma)`, and power is
reported at **A*sigma in {0.03, 0.10, 0.15} MPa**, with the **status flag set by the adverse end,
A*sigma = 0.15 MPa**. The three figures are Merton's, adopted by Popper in S-14, and are
**PARTIALLY VERIFIED** there; they are used here for the power arithmetic only and none is quoted
as a measurement.

## §5. What a FAIL means (R2-2 mandate 4)

If the engine does not fire on >= 2 sources, the failure is ambiguous between "engine broken" and
"box/catalogue too small". It is resolved by injecting a synthetic Landers-scale transient into a
simulated catalogue at the observed cell background rates and confirming recovery. That injection
is run **unconditionally** in this execution, because it is also the deliverable
(the detection-threshold curve), not only the failure branch.
