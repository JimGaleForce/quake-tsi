# K-092 SCALAR PROVENANCE ATTACHMENT (§P7-23(D))

Status: DRAFT ATTACHMENT, NOT COMMITTED. Prepared 2026-08-13 to satisfy the gate that
`engine/conventions_d.py::assert_scalar_provenance` enforces and that `K092_FREEZE.md`
lists as pre-scoring gate 1 ("D-0 convention + scalar translation attached").

Nothing here is evidence. No Alaska-Aleutian phase statistic was computed. The only
Alaska datum touched is the ONE displayed number already public in the freeze file
(the 1991-05-30 Sand Point readout), reproduced solely to prove the provenance chain.

---

## 1. THE ANSWER TO `PROVENANCE_JIM["question"]`

> *"which tidal display was Jim reading ... and from which source?"*

**VERTICAL (RADIAL) SOLID-EARTH BODY-TIDE DISPLACEMENT, IN CENTIMETRES, AT THE SURFACE
POINT**, computed in-browser by Jim's earthquake viewer
`D:\CODE\git\earth-tides-globe`, function

    solidTideDisplacementCm(timeMs, latDeg, lonDeg)   -- src/utils/astro.ts:189-212

with the "rising/falling X cm/h" trend from

    solidTideRateCmPerHour(timeMs, latDeg, lonDeg)    -- src/utils/astro.ts:215-227

and rendered by `src/components/TideInspector.tsx` (the panel headed
`SOLID-EARTH TIDE`, value line at :170-182, `falling` string at :179). The same
expression is duplicated in the GPU fragment shader `src/scene/tides.ts:34-67` for the
globe overlay and in the quake callout at `src/components/Globe.tsx:735,814`; the three
are the same formula and the repo's `scripts/verify-tides.mjs` is the parity gate.

**It is NOT** water level, NOT tidal potential, NOT volumetric strain, NOT a resolved
Coulomb stress, and NOT ocean tide. It is body-tide vertical displacement only.

## 2. THE ALGORITHM, EXACTLY

Per body (Moon, Sun), with `n` the site unit vector and `d` the sub-body unit vector:

    disp_cm = moonAmpCm * ( P2(cos z_moon) + (R/d_moon) * P3(cos z_moon) )
            + sunAmpCm  *   P2(cos z_sun)

    ampCm(body) = 100 * h2 * ( G * M_body * R^2 / d^3 ) / g

so it is exactly `h2 * W2 / g` per body, plus a **lunar degree-3 term** (`P3`, scaled by
R/d_moon ~ 1/60, i.e. ~1.7 % of the lunar degree-2 amplitude).

| item | app value (`src/utils/astro.ts`) |
|---|---|
| Love number | `LOVE_H2 = 0.612` (h2; used for displacement). `LOVE_K2 = 0.303` exists but is used ONLY in the separate Longman gravimetric routine, not in the displayed scalar. |
| Shida number l2 | **none** — no horizontal displacement, no strain |
| harmonics | degree 2 (Moon + Sun) + degree 3 (Moon only) |
| ocean loading | **absent** |
| body constants | `G=6.674e-11`, `M_MOON=7.342e22`, `M_SUN=1.989e30`, `R_EARTH=6.371e6` (constant, spherical), `g=9.81`, `AU=1.496e11` |
| latitude | **geographic latitude used directly as geocentric** (no WGS-84 flattening correction) |
| depth | not modelled — surface point only |

**Ephemeris.** Low-precision Meeus, coded inline (`sunPosition` :63-73, `moonPosition`
:75-131): 13-term lunar longitude series, 7-term latitude series, 8-term distance
series; 2-term solar equation of centre. GMST from the IAU-1982 series (:51-61). Both
bodies are reduced to a sub-body geographic point (declination, RA minus GMST) and the
site zenith angle is the geocentric spherical-trig one.

**Time handling.** `timeMs` is a JavaScript epoch (UTC). `dateToJulianDay` converts
straight to JD with no ΔT: TT is not distinguished from UTC, GMST is evaluated on the
same JD, positions are **geocentric** (not topocentric), and there is **no light-time
correction**. All UTC in, UTC out.

**Sign convention.** From the source comment at `astro.ts:187`: *"Positive = crust
lifted away from the Earth's center."* So the displayed **-13.2 cm means the ground is
displaced DOWNWARD by 13.2 cm**, and "falling 0.6 cm/h" means the displacement is
becoming more negative at 0.6 cm/h. The "neutral" the UI draws (the zero line in
`TideInspector.tsx:92-96`) is **literal zero of this scalar**, not a running mean.

Rate estimator: symmetric central difference over **±10 minutes**, converted to cm/h,
with a UI dead-band — `|rate| < 0.05 cm/h` prints "at turning point" rather than
rising/falling (`TideInspector.tsx:176-181`).

## 3. REPRODUCTION OF THE FROZEN 1991 READOUT — **EXACT MATCH**

Run against the app's own TypeScript, unmodified, via Node 25.2.1 native type
stripping:

    import('D:/CODE/git/earth-tides-globe/src/utils/astro.ts')
    t = Date.UTC(1991,4,30,13,17,41)          // 1991-05-30T13:17:41.000Z
    solidTideDisplacementCm(t, 54.57, -161.61)  ->  -13.206888677138735
    solidTideRateCmPerHour (t, 54.57, -161.61)  ->   -0.5757627956078135

The UI formats `disp.toFixed(1)` and, for a negative rate, `falling ${(-rate).toFixed(1)} cm/h`:

    -13.206888... -> "-13.2 cm"        -0.575762... -> "falling 0.6 cm/h"

which is **character-for-character the string recorded in `K092_FREEZE.md`**:
"-13.2 cm falling 0.6 cm/h". The provenance chain is closed: the frozen scalar is
`solidTideDisplacementCm` at the epicentre in UTC, and this event sits in the
below-neutral/falling quadrant in it.

## 4. MAPPING TO `replication/engine/sitetide.py`

**Is the app's scalar our `radial_disp_m`?** Structurally yes, up to constants:

| | app | `sitetide.py` |
|---|---|---|
| formula | `h2 * W2 / g` (+ lunar P3) | `radial_disp_m = H2 * W2 / g` (:286) |
| h2 | 0.612 | 0.6078 (IERS 2010 TN36 Table 6.3) |
| degree 3 | lunar P3 included | **absent** |
| site geometry | sphere R=6.371e6, geographic lat as geocentric | WGS-84 ellipsoid radius, geodetic→geocentric lat, depth |
| GM_moon | 4.9002e12 (G·M) | 4.9028001e12 |
| g | 9.81 | 9.80665 |
| time scale | UTC-as-TT, geocentric, no light-time | **same** (`ephemeris.py:19`) |
| primary emitted scalar | this one | `areal_strain` (`SCALAR_FOR_PHASE`) |

**The key algebraic fact.** `sitetide.areal_strain = AREAL_FACTOR * W2 / (g·a)` and
`sitetide.radial_disp_m = H2 * W2 / g`, both **positive multiples of the same W2**.
They are therefore exactly proportional, so **they induce identical sign, identical
slope sign, and identical quadrant**. Verified numerically: over the 60-day series
below, the quadrant masks from `areal_strain` and `radial_disp_m` are bit-identical.
So the translation question is *not* "displacement vs areal strain" — those two are the
same quadrant by construction — it is only "app implementation vs our implementation".

### Measured agreement, Sand Point coordinates (54.57 N, 161.61 W), surface

Both series on the same 5-minute grid; app values dumped from its own TypeScript,
sitetide values from `site_tide(jd, 54.57, -161.61, 0.0)`.

| span | Pearson r (app disp vs `radial_disp_m`) | RMS amplitude ratio | max abs diff | RMS diff |
|---|---|---|---|---|
| 7 days (2020-03-01 → 03-08) | **0.99995128** | 1.00231 | 0.293 cm | 0.086 cm |
| 60 days (2020-01-01 → 03-01) | **0.99991063** | 1.00522 | 0.376 cm | 0.124 cm |

Pearson r of app displacement against `areal_strain` is identical to eight digits
(0.99995128 on the week), confirming the exact proportionality.

**Phase lag: none detectable.** Best cross-correlation lag on the week = **0 steps
(0.0 min)** on a 5-minute grid. Extremum-by-extremum (parabolic sub-sample refinement,
nearest-neighbour matched over the 60 days, 103 matched pairs each):

* maxima: app minus sitetide = **+0.035 min mean, sd 1.08 min, max 3.9 min**
* minima: app minus sitetide = **+0.161 min mean, sd 4.23 min, max 18.1 min**

i.e. sub-minute mean offset against a ~12.4 h cycle (≈ 5e-4 rad). No systematic lag.
The minima spread is larger because the diurnal-modulated minima are flatter.

### Quadrant transfer — the actual §P7-23(D) requirement

Quadrant = (level < 0) AND (rate < 0), each series judged in its own scalar with the
app's own ±10-min central-difference rate estimator applied to both:

| span | app quadrant duty cycle | sitetide quadrant duty cycle | agreement | disagreeing samples |
|---|---|---|---|---|
| 7 days | 0.3699 | 0.3659 | **99.405 %** | 12 / 2017 |
| 60 days | 0.4000 | 0.3998 | **99.369 %** | 109 / 17281 |

Excluding the app's own "at turning point" dead-band (`|rate| ≥ 0.05 cm/h`) and
near-zero levels (`|level| ≥ 0.5 cm`), 60-day agreement rises to **99.868 %**
(22 / 16659).

**Every one of the 109 disagreements is a boundary case.** Of the 109: 85 are
slope-sign flips, 24 are level-sign flips, none are both. Bounding them:

* every slope flip has `min(|rate_app|, |rate_sitetide|) ≤ **0.038 cm/h**` — i.e. one of
  the two series is within a few minutes of a turning point;
* every level flip has `min(|level_app|, |level_sitetide|) ≤ **0.095 cm**` — i.e. within
  a hair of the zero crossing.

There is **no region of the cycle where the two implementations disagree away from a
quadrant boundary.** The 5 % of time spent inside those two boundary bands carries
100 % of the disagreement.

## 5. TRANSLATION STATEMENT (what §P7-23(D) is owed)

**The frozen scalar is: solid-earth body-tide VERTICAL DISPLACEMENT in cm at the
epicentral surface point, UP POSITIVE, as computed by
`earth-tides-globe/src/utils/astro.ts::solidTideDisplacementCm` (degree-2 lunar+solar
plus lunar degree-3, h2 = 0.612, low-precision Meeus geocentric ephemeris, UTC, no
ocean loading), with "falling" = the ±10-minute central-difference time derivative
negative.** The frozen quadrant is (displacement < 0) AND (d displacement/dt < 0).

**Does quadrant classification transfer to `engine/sitetide.py`?** — **APPROXIMATELY,
AND THE ERROR IS BOUNDED AND BOUNDARY-CONFINED, NOT SYSTEMATIC.**

1. **Within our engine the translation is EXACT.** `radial_disp_m` and the declared
   `areal_strain` are positive multiples of the same W2, so the quadrant is *identical*
   in the two. The scalar choice inside `sitetide.py` costs nothing. (This retires the
   §P7-23(D) worry for this particular pair: displacement and areal strain do **not**
   differ by a phase offset. Coulomb-on-a-thrust still would; that translation remains
   un-attempted and separately declared.)
2. **Across implementations the disagreement rate is 0.63 % of time (60 days at the Sand
   Point site), all of it within `|rate| < 0.04 cm/h` of a turning point or
   `|level| < 0.1 cm` of a zero crossing.** For an event set sampled effectively
   uniformly in time, that is the expected misclassification rate.
3. **Recommendation, if the program wants an exactness guarantee rather than a
   0.63 % rate:** score D-12/D-13 with the app's code verbatim (the Node import used in
   §3 above runs `astro.ts` unmodified and is reproducible in one line), and report the
   `sitetide.py` classification alongside as a robustness check. If instead
   `sitetide.py` is used as the scorer, the pre-registered protocol should declare the
   boundary bands above as the known translation tolerance BEFORE scoring, since
   whether a marginal event lands in-quadrant is exactly what those bands decide.
4. **Unchanged by any of this:** the DIRECTION stays literature-derived (§P7-22 Q5 as
   preserved by §P7-23(D)), ocean loading is absent from **both** implementations, and
   `SCOPE_FLAGS` / `CONVENTION_SCOPE` apply verbatim — a subduction-zone coastal site is
   exactly where the missing ocean load is largest, so no amplitude and no absolute psi
   from either implementation is a measurement.

## 6. WHAT IS UNVERIFIED

* **UNVERIFIED:** that the specific on-screen readout Jim saw came from `TideInspector`
  rather than from the `Globe.tsx` quake callout. Both call the same
  `solidTideDisplacementCm` / `solidTideRateCmPerHour` pair and format identically, so
  the scalar is the same either way; only the UI path is unattested.
* **UNVERIFIED:** that the app's source at the git revision Jim was viewing is the
  revision read here (working tree of `D:\CODE\git\earth-tides-globe`, 2026-08-13). The
  exact reproduction of "-13.2 cm falling 0.6 cm/h" is strong circumstantial evidence
  that it is, but it is not a commit-hash attestation. If a hash is wanted, pin it.
* **NOT ATTEMPTED:** any phase, quadrant, or statistic for any Alaska-Aleutian event
  other than the single already-public 1991-05-30 readout reproduced in §3. D-7 remains
  STRUCK.
* The 7-day and 60-day comparison spans (2020) were chosen to be seismically
  uninformative for this claim; they exercise the two implementations, not the
  catalogue.

## 7. REPRODUCTION RECIPE

    # app scalar, its own code, unmodified (Node >= 23 for TS type stripping)
    node --input-type=module -e "
      const m = await import('file:///D:/CODE/git/earth-tides-globe/src/utils/astro.ts');
      const t = Date.UTC(1991,4,30,13,17,41);
      console.log(m.solidTideDisplacementCm(t,54.57,-161.61),
                  m.solidTideRateCmPerHour(t,54.57,-161.61));"
    # -> -13.206888677138735 -0.5757627956078135

    # engine scalar at the same instant
    python -c "
    import numpy as np, datetime as dt, sys; sys.path.insert(0,r'D:\CODE\git\quake\replication')
    from engine import sitetide as S
    jd = dt.datetime(1991,5,30,13,17,41,tzinfo=dt.timezone.utc).timestamp()/86400.0 + 2440587.5
    r = S.site_tide(np.array([jd-10/1440., jd, jd+10/1440.]), 54.57, -161.61, 0.0)
    c = r['radial_disp_m']*100
    print(c[1], (c[2]-c[0])/(1/3.))"
    # -> -12.9850 cm, -0.5963 cm/h   (same sign, same quadrant, 0.22 cm apart)

Both place the event in the below-neutral, falling quadrant.

---

*Prepared as the §P7-23(D) attachment. Not committed. `PROVENANCE_JIM["answer"]` in
`engine/conventions_d.py` is left unset — filling it is Jim's act, not the pipeline's.*
