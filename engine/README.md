# EQ-23 predictability engine (v1)

Vectorized earthquake-forecast backtesting harness: covariate plugins in, bits/event
and green/red cells out. Design is ruled in [SPEC.md](SPEC.md); this file is usage.

Run everything from `replication/`, always with `python -u`.

## Usage

```bash
# 1. build (or refresh) the design cache from data/comcat_world/*.csv
python -u -m engine.cli cache
python -u -m engine.cli cache --rebuild --dlat 2 --dlon 2     # coarser grid

# 2. see what's registered
python -u -m engine.cli list

# 3. fit the ETAS baseline (once; ~6 min, cached to out/cache/etas_params.json)
python -u -m engine.cli fit-etas
python -u -m engine.cli fit-etas --refit          # discard the cache and refit

# 4. explore: all sliders allowed, exploration period only
python -u -m engine.cli run --covariate recent_rate --mode explore
python -u -m engine.cli run --covariate recent_rate --mode explore --baseline etas
python -u -m engine.cli run --covariate quiescence  --mode explore \
       --params '{"window_days": 180, "radius_km": 150}'

# 5. holdout: config must be frozen in a file first, and it runs ONCE
python -u -m engine.cli run --covariate recent_rate --mode holdout \
       --config engine/configs/holdout_recent_rate.json

# tests
python -u -m pytest engine/tests -q
```

Outputs land in `engine/out/` (gitignored): `cache/design_<key>.npz`,
`run_<runid>.json` (full result incl. Molchan trajectory), and
`cells_<runid>.json` (per-cell green/red export for later visualization; there is no
plotting in v1).

Logs that are *not* disposable: `engine/HOLDOUT_LOG.jsonl` and
`engine/EXPLORE_COUNT.jsonl`.

## The two baselines — which one absorbs what

Every report prints exactly one caveat line, and it switches with `--baseline`. It is
the single most important line the engine emits, because it says what the number you
are reading has already been charged for:

| `--baseline` | printed caveat | absorbed by the baseline | still available to a covariate |
| --- | --- | --- | --- |
| `climatology` (default, `climatology-v1`) | `baseline=climatology-v1 (clustering NOT absorbed; sniff-grade only)` | *where* earthquakes happen (time-stationary per-cell rate `mu(x)`) | everything about *when* — including all aftershock clustering |
| `etas` (`etas-v1`) | `baseline=etas-v1 (isotropic kernel; anisotropy/mechanism NOT absorbed)` | Omori-Utsu clustering: background `nu*mu(x)` plus an isotropic space-time triggering kernel | anisotropy of the aftershock cloud, focal mechanism / fault geometry, depth, rupture length, secular rate change, completeness drift |

### climatology-v1

Time-stationary per-cell Poisson rate `mu(x)`, fitted on the exploration period. It
knows *where* earthquakes happen and nothing about *when* — specifically, it does not
know that earthquakes come in aftershock cascades. So any covariate carrying the
information "there were recent nearby earthquakes" scores large, highly significant,
entirely real information gain against it, and essentially none of that gain is
evidence of predictability beyond Omori-Utsu clustering.

That is why `recent_rate` exists and why it must score positive under this baseline:
it is the engine's smoke test (if it comes out flat, the plumbing is broken) **and**
the calibration of how much apparent skill this baseline hands out for free.

### etas-v1

```
lambda(c,t) = nu*mu(c) + K * sum_{i: 0 < t-t_i <= 365 d} 10^(alpha*(m_i-4.5))
                             * (t-t_i+c0)^(-p) * G_sigma(dist(c, cell_i))
```

Six parameters `(nu, K, alpha, c, p, sigma)` fitted by maximising the Poisson
log-likelihood **on the exploration window only** (L-BFGS-B on log-parameters within
declared bounds, then a Nelder-Mead polish), cached with the fit window in
`engine/out/cache/etas_params.json`, and printed in the run header. `--baseline etas`
loads that cache, or fits if it is absent.

Fitted on the shipped `comcat_world` catalogue (1628 cells, exploration days
365–8081, 46,585 M≥4.5 events, ~6 min, 1100 objective evaluations):

```
nu = 0.5444   K = 0.0602   alpha = 0.5000*   c = 0.0010* d   p = 1.1599   sigma = 45.6 km
                                    * pinned at the lower bound (see below)
background share 0.546 of total intensity / triggered 0.454
```

`alpha` and `c` both sit on their lower bounds, and the run header prints a
`WARNING: parameters pinned at a bound` line when they do. `c` is expected there —
at day resolution nothing distinguishes `c = 0.001 d` from `c = 0.05 d`. `alpha = 0.5`
is *not* expected (published ETAS fits usually land near 0.8–1.0) and is a real
finding about this configuration: with a 1° cell, a 1-day step and an M4.5 floor, the
productivity contrast between a M5 and a M7 mainshock is being compressed, most
plausibly because large events' aftershocks spill outside a single cell and the
isotropic 45 km kernel under-collects them. Treat the magnitude scaling of this
baseline as weakly identified.

Declared approximations, all of them cheap-on-purpose:

- **Day resolution.** `t - t_i` is a whole number of days and is `>= 1` by
  construction — day `t` sees only events strictly before `t`, the same causality
  contract the covariates sign. A consequence: `c` is essentially unidentifiable
  below one day and runs to its lower bound.
- **Omori truncation at 365 d.** This keeps only ~64% of the kernel's formal mass at
  the fitted `p = 1.16` — an Omori tail is heavy, and a third of it lives beyond a
  year. It costs little in likelihood (measured on this catalogue: `--omori-trunc 20`
  scores +0.911 bits/event over climatology, `--omori-trunc 365` scores +0.962),
  because the far tail is nearly flat in time and gets absorbed by the `nu*mu`
  background term instead. So the truncation mostly *reassigns* decade-long
  sequences from triggering to background: do not read `nu*mu` as a declustered
  background rate. The first 365 days of the record are dropped from every ETAS run
  as un-warmed burn-in.
- **Isotropic Gaussian kernel on cell-centre ground distance**, column-normalised
  over the restricted domain (`sum_c G(c,s) = 1`), so no triggered mass leaks out of
  the domain and `K` reads as "triggered events per unit productivity".
- **Background = climatology's spatial shape**, rescaled by `nu`. The background is
  not re-estimated by stochastic declustering; `mu(x)` is the same
  aftershock-contaminated climatology, so `nu*mu` is a smoothed-background
  approximation, not a declustered one.

What it buys, and what it then takes away (exploration window, days 365–8081, 46,585
events, all **in-sample**):

| | climatology-v1 | etas-v1 |
| --- | --- | --- |
| baseline log-likelihood | −255,210.0 | −224,138.7 (**+0.962 bits/event**) |
| `recent_rate` bits/event | +0.6264 | **+0.0284** (4.5% survives) |
| `quiescence` bits/event | +0.1850 (beta −0.105) | **+0.0118** (beta −0.027) |
| `anisotropy` bits/event | +0.000001 (beta +0.0014 ± 0.0057) | +0.000009 (beta −0.0044 ± 0.0058) |

That collapse *is* the point: nearly all of what `recent_rate` and `quiescence`
appeared to know was aftershock clustering. `anisotropy` was flat before and stays
flat — as built (azimuths to cell centres, no mechanisms) it never had anything for
either baseline to eat.

`etas-v1` is still not a discovery machine. It absorbs clustering; it does not absorb
anything on the right-hand column of the table above, and `--mode explore` numbers
remain in-sample upper bounds.

Two further limitations, stated because they change how numbers should be read:

- **Domain restriction.** A global 1° grid is ~98% zero-rate ocean and craton. The
  forecast domain is restricted to cells with ≥1 event during the *exploration*
  period (1628 cells for the shipped data). Events that land in never-before-active
  cells during the holdout are counted and printed (284 for the shipped data), not
  silently dropped — they are a real cost charged to this domain definition.
- **`--mode explore` numbers are IN-SAMPLE.** Beta is fitted and evaluated on the
  exploration period, and the run prints a note saying so. They are upper bounds. Only
  `--mode holdout` produces a genuinely out-of-sample bits/event.

## The holdout gate

`--mode holdout` requires `--config <file.json>`. The engine takes the sha256 of the
canonical JSON, refuses to run if that hash already appears in `HOLDOUT_LOG.jsonl`,
and on completion appends `{hash, config, results, n_explore_runs}`. `n_explore_runs`
is the number of exploration runs logged since the last holdout run — your reported
multiplicity — and it is printed next to the holdout result.

No flag bypasses the refusal. Deleting a line from the log is a human act, not
something the engine will do for you.

## Adding a covariate

Drop a file in `engine/covariates/`; it is auto-discovered.

```python
# engine/covariates/my_thing.py
import numpy as np
from . import register

DEFAULTS = {"m_min": 5.0, "radius_km": 200.0, "days": 90}

@register("my_thing", defaults=DEFAULTS, burn_in=90, describe="one-line description")
def compute(ctx, p):
    p.setdefault("burn_in", int(p["days"]))          # days of history you need
    n = ctx.past_counts(p["m_min"], p["radius_km"], p["days"])
    return np.log1p(n)                                # -> (n_cells, n_days) float32
```

Rules the harness enforces for you:

- **Return shape** `(ctx.n_cells, ctx.n_days)`, float32, all finite. Violations raise.
- **Causality.** `z[:, t]` may use only events strictly *before* day `t`. Build from
  `ctx.past_counts(...)` / `ctx.trailing_sum(...)`, which are exclusive-of-today by
  construction. `engine/tests/test_causality.py` scrambles all events after a cut day
  and asserts your covariate does not move before it; add your covariate to the
  `COVARIATES` list there. The same shuffle is applied to the ETAS baseline's own
  intensity in `engine/tests/test_etas.py::test_etas_intensity_is_causal` — the
  baseline signs the identical contract.
- **Baseline offset.** Under `--baseline climatology` the offset is a 1-D per-cell
  rate; under `--baseline etas` it is the 2-D `lambda_etas(cell, day)`. Nothing in a
  covariate changes: `score.py` fits `lambda = offset * exp(a + beta*z)` with the
  same intercept-refit convention either way, and the ETAS baseline's 365-day
  burn-in is unioned with your covariate's.
- **Never touch raw CSVs.** Everything comes through `ctx`.
- **Standardisation is automatic** and uses exploration-window mean/sd only, so a
  holdout run applies exactly the transform that was fitted in training.
- **Burn-in** is dropped from both the fit and the evaluation window.

Context helpers available to plugins:

| helper | returns |
| --- | --- |
| `ctx.day_counts(m_min)` | `(n_cells, n_days)` daily event counts |
| `ctx.trailing_sum(counts, days)` | strictly-causal trailing window sum |
| `ctx.neighbour_matrix(radius_km, r_inner_km=0)` | `(n_cells, n_cells)` annulus mask |
| `ctx.past_counts(m_min, radius_km, days, r_inner_km=0)` | counts in a space-time neighbourhood |
| `ctx.past_azimuth_resultant(m_min, r0, r1, days)` | `(R, N)` circular resultant of azimuths |
| `ctx.dist_km`, `ctx.az_rad`, `ctx.grid` | cell-pair geometry |

Every neighbourhood helper is one float32 GEMM against a trailing cumsum, so a
covariate over 1628 cells × 11544 days computes in well under a second.
