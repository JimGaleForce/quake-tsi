# Prospective SoCal weekly forecast: frozen protocol

EQ-25 Phase 4. Popper, `HYPOTHESIS_LEDGER.md` P7-26(5) rank 1: "v1 = frozen ETAS vs Poisson and nothing else." B-2's scope line (L2989) records "post-2018 period untested"; this log closes that gap one week at a time. Emitter and scorer: `prospective_socal.py`. Log: `results_prospective_socal_log.json` (append only). Tests: `tests/test_prospective_socal.py`.

## What is committed, and when

At an emission time T (UTC, to the second), `emit` writes exactly one record:

- Seven per-day expected M>=2.5 counts for the bins [T+d, T+d+1), d = 0..6, as the mean over 1000 ETAS cluster simulations (`etas_forecast.simulate`, seed 20260813).
- The per-day analytic no-future-triggering lower bound (recorded, not the forecast), and the simulation count.
- Two Poisson baseline rates, both from data strictly before T.
- sha256 of the ComCat cache as used, of `results_exp_h.json`, of `prospective_socal.py` and of `etas_forecast.py`; the frozen parameters echoed; the trailing-30-day count.
- `commitment_hash`: sha256 of the canonical JSON (sorted keys, no whitespace) of the committed block. Scoring adds a sibling block and never edits the commitment.

Parameters are TAKEN, not fit: `results_exp_h.json -> train_fit.frozen_params` (mu, K, alpha, c, p, M0 = 2.5), fit on SCSN 1982-2010 and scored once on 2010-2018. Nothing is ever refit. Domain: ComCat, type == earthquake, lat [31.5, 38.0], lon [-122.0, -113.5], M >= 2.5, whole box, temporal only.

## The two baselines

- **PRIMARY, `poisson_trailing_365d`**: observed M>=2.5 count in [T-365d, T) / 365, constant across the seven bins. The yardstick a working forecaster actually has.
- **SECONDARY, `poisson_train_rate`**: 3.3430565924513234 per day, the frozen B-2 training rate (`results_exp_h.json -> test.baselines.poisson_train_rate`), so the prospective number is directly comparable to B-2's retrospective +1.87 bits/event.

## The scoring rule

A record is DUE when now >= T + 7 d + 3 d. The **3-day latency is frozen**: ComCat revises recent solutions, so a week is scored only once its catalogue has had three days to settle, and the rule is frozen so it cannot be tuned per week. Scoring before due is refused (exit 2).

At scoring: refresh the catalogue, count observed M>=2.5 events in the seven one-day bins, and take the Poisson log-likelihood of those counts under the committed ETAS expectations and under each constant baseline (rate x 1 day). Information gain is `(ll_model - ll_baseline) / (n_events * ln2)` bits per event, the `engine/score.py` `information_gain` definition. If n_events == 0, bits per event are recorded as null and the total log-likelihood difference is recorded in nats instead; no division by zero occurs. The log carries a cumulative summary (weeks scored, total events, total bits over each baseline).

Every emitted record is scored when due. Records are never withdrawn and never rescored, the summary is over all of them, and `commitment_hash` must verify before any score is written.

## What a positive result would and would not license

Exactly one statement: **a frozen temporal ETAS beats a Poisson rate on unseen, post-commitment SoCal M>=2.5 seismicity, whole box.** That is B-2 proven prospectively and nothing else. It licenses nothing spatial (the model has no spatial component), nothing about magnitude beyond the M>=2.5 counting threshold, no beyond-ETAS mechanism, and it is not earthquake prediction. A negative is equally publishable and is pre-committed to be reported.

Standing risk, recorded before week 1 was scored: the parameters were fit on SCSN while the driving catalogue is ComCat, and from a quiet state the frozen model's absolute 7-day count runs below the observed box rate. Count calibration and timing skill can therefore move in opposite directions. Both are in the log; neither is adjustable.

## Scheduling and publication

Running this weekly, and publishing the curve, are **separate decisions** not made here. This file freezes the instrument only.
