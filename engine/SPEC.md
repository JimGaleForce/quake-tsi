# EQ-23 Predictability Engine — v1 spec

Vectorized backtest harness: covariate plugins in, bits/event and green/red cells out.
Tracker: Linear EQ-23. Design ruled 2026-08-10 by supervisor session; do not re-litigate
the fork decisions below — flag disagreements, don't silently change them.

## Fork decisions (ruled)

1. **Location:** this directory (`replication/engine/`), a plain python package run as
   `python -m engine.cli ...` from `replication/`. Outputs to `engine/out/` (gitignored).
2. **Baseline (v1):** per-cell Poisson climatology `mu(x)` fitted on the EXPLORATION
   period only. NOT ETAS. This is v1's declared gap: clustering leaks into covariates, so
   any covariate that proxies "recent aftershocks" will score fake skill. Every report
   the engine prints must carry the line: `baseline=climatology-v1 (clustering NOT
   absorbed; sniff-grade only)`. An `etas` baseline slot exists in the interface,
   unimplemented, raising NotImplementedError with that explanation.

   **AMENDED 2026-08-10 (ledger P4-4, ruled GATING). The original ruling is left
   above, unedited.** The slot is now filled: `engine/baseline.py:EtasV1` is a
   space-time ETAS fitted by Poisson ML on the exploration window only. Decision 2
   stands as the description of `--baseline climatology`, which remains the default
   and still prints the line above verbatim. Under `--baseline etas` the printed
   line becomes `baseline=etas-v1 (isotropic kernel; anisotropy/mechanism NOT
   absorbed)`. The requirement is unchanged in substance: every report carries a
   caveat naming exactly what its baseline did and did not absorb. Model, declared
   approximations, and the measured collapse of `recent_rate`/`quiescence` are in
   README.md, "The two baselines".
3. **Split:** temporal. Exploration = first 70% of the catalog time span, holdout = last
   30%. Forecast skill only counts forward in time; spatial or random splits are wrong
   here and must not be offered.
4. **Holdout discipline (the p-hacking guard, non-negotiable):**
   - `--mode explore` (default): all sliders/params allowed, runs only on exploration.
   - `--mode holdout`: requires `--config <file.json>`; the engine hashes the config
     (sha256 of canonical JSON), refuses to run if that hash already appears in
     `engine/HOLDOUT_LOG.jsonl`, and appends `{hash, config, results, n_explore_runs}`
     on completion. `n_explore_runs` = count of exploration runs logged since the last
     holdout run (from `engine/EXPLORE_COUNT.jsonl`, one line per explore run) —
     reported multiplicity, always printed next to the holdout result.
   - No flag bypasses the refusal. Deleting the log is the human's act, not the engine's.
5. **Scoring:** per space-time cell (1° × 1 day default, configurable), target = count of
   M≥m events. Model: `lambda = mu(x) * exp(beta * z(x,t))`, beta fitted by Poisson
   regression on exploration. Metrics: information gain in bits/event vs baseline
   (primary), Molchan trajectory (tau vs miss rate), and a green/red JSON export
   (`out/cells_<runid>.json`: per-cell sign of realized log-likelihood gain) for later
   visualization. No plotting in v1.

## Layout

```
engine/
  SPEC.md            this file
  __init__.py
  datasets.py        load replication/data/comcat_world/*.csv -> one deduped catalog
  grid.py            grid spec; (lat,lon,t) -> cell indexing
  design.py          precompute + cache per-cell/per-day arrays (out/cache/*.npz)
  baseline.py        climatology mu(x); etas slot (NotImplementedError)
  splits.py          temporal split; holdout gate + logs
  score.py           Poisson regression, bits/event, Molchan, green/red export
  covariates/
    __init__.py      registry + plugin interface
    recent_rate.py   log(1+count of M>=4.5 within R km, past D days) — sanity covariate
                     (MUST come out positive-skill under climatology baseline; that is
                     the engine's own smoke test AND the demonstration of gap #2)
    quiescence.py    K-076-style: -(observed - expected)/sqrt(expected) over trailing
                     window W vs exploration-period climatology
    anisotropy.py    K-081-sniff-style: resultant length R of azimuths of prior M>=5.5
                     in an annulus [r0,r1] within past D days (structure-aware caveat
                     printed, not solved, in v1)
  cli.py             run/list/cache subcommands
  tests/             pytest; unit tests + one end-to-end on a synthetic catalog with a
                     PLANTED covariate effect (counted invariant: recovered beta within
                     tolerance; bits/event > 0 on the planted covariate, ~0 on a
                     scrambled one)
```

**AMENDED 2026-08-11 (supervisor; Linear EQ-24).** `mine` mode added — an
exploration-only pattern miner (GENERATOR, never evidence): `ephemeris.py`
(zero-download Meeus-style solar/lunar features), `mine.py` (features, GLM/period/mark
statistics, per-test-appropriate surrogate nulls, BH-FDR, harmonic ladder {P/3, P/2, P,
2P, 3P}, aliasing audit), `mine_session.py` (checkpointed/resumable session driver,
report.md + stubs.json). CLI: `mine --quick|--overnight [--surrogates N]`. Mine sessions
append to EXPLORE_COUNT.jsonl so holdout multiplicity includes mining. Runs only on the
exploration window; no holdout hash may be spent from mine mode. Design rulings and the
two mandatory warnings (generator-not-evidence; Mignan & Broccardo caution) live in
Linear EQ-24 and README.md.

## Plugin interface

```python
# covariates/__init__.py
@register("name")
def compute(ctx: EngineContext, params: dict) -> np.ndarray:  # (n_cells, n_days) float32
```
`EngineContext` carries: catalog dataframe, grid, day index, cached design arrays, and
`ctx.past_counts(m_min, radius_km, days)` style helpers so plugins never touch raw CSVs.
Plugins must be causal: `z[:, t]` may use only events strictly before day `t`. The
end-to-end test includes a causality check (shuffle future events; z must not change).

## Invariants (assert, loudly)

- Catalog load: expected ~72k events from comcat_world; assert 50k < n < 150k, assert
  dedup by event id actually removed duplicates when boxes overlap.
- Design matrix: assert no NaN/Inf; assert >0 events land in >100 distinct cells.
- Any bulk transform reports counts in/out; "produced 0" raises.
- Poisson fit: assert convergence flag; report beta with SE.

## Environment gotchas (known, do not rediscover)

- miniconda base; numpy 2.5.1 / pandas 3.0.3 / scipy 1.18.0; always `python -u`.
- pandas 3.0 datetime parsing quirks; parse ComCat `time` with `utc=True`.
- ComCat CSVs: column order varies by source; select by NAME. Sub-M2.5 spurious rows
  exist in some files; filter to the file's stated floor.
- Windows paths; run tests from `replication/` with `python -u -m pytest engine/tests`.
- Do NOT commit anything to git.
