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

# 3. explore: all sliders allowed, exploration period only
python -u -m engine.cli run --covariate recent_rate --mode explore
python -u -m engine.cli run --covariate quiescence  --mode explore \
       --params '{"window_days": 180, "radius_km": 150}'

# 4. holdout: config must be frozen in a file first, and it runs ONCE
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

## What v1 actually measures — and the gap

**`baseline=climatology-v1 (clustering NOT absorbed; sniff-grade only)`** is printed on
every report, and it is the single most important line the engine emits.

The v1 baseline is a time-stationary per-cell Poisson rate `mu(x)`, fitted on the
exploration period. It knows *where* earthquakes happen. It knows nothing about *when*
— specifically, it does not know that earthquakes come in aftershock cascades. So any
covariate that carries the information "there were recent nearby earthquakes" scores
large, highly significant, entirely real information gain against it, and essentially
none of that gain is evidence of predictability beyond Omori-Utsu clustering.

That is why `recent_rate` exists and why it must score positive: it is the engine's
smoke test (if it comes out flat, the plumbing is broken) **and** the calibration of
how much apparent skill this baseline hands out for free. Read every other covariate's
bits/event against `recent_rate`'s, not against zero.

Closing the gap means an ETAS baseline. The slot exists —
`engine/baseline.py:EtasBaseline` — and raises `NotImplementedError` with this
explanation rather than quietly doing something weaker. Nothing this engine prints in
v1 is a discovery.

Two further v1 limitations, stated because they change how numbers should be read:

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
  `COVARIATES` list there.
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
