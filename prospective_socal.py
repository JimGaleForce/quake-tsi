"""Phase 4 / EQ-25: WEEKLY, HASH-COMMITTED, PROSPECTIVE SoCal forecast.

Popper, HYPOTHESIS_LEDGER.md section P7-26(5) rank 1: "v1 = frozen ETAS vs Poisson
and nothing else." B-2's own scope line (L2989) records "post-2018 period untested".
This module closes that gap by committing, in advance and with a hash, seven per-day
expected M>=2.5 counts for southern California from the FROZEN EXP-H ETAS, and by
scoring them as bits per event over a Poisson baseline once the week has elapsed.

NOTHING here is refit. All ETAS parameters come from
results_exp_h.json -> train_fit.frozen_params via etas_forecast.load_frozen(), which
carries its own no-refit assertions. The catalogue fetch, the ETAS cluster simulation
and the analytic no-future-triggering bound are IMPORTED from etas_forecast.py, not
copied.

Two subcommands:

    python -u prospective_socal.py emit  [--now ISO]
    python -u prospective_socal.py score [--record ID | --all-due]

emit  writes exactly one new record to the append-only log
      results_prospective_socal_log.json, containing the seven ETAS daily
      expectations, two Poisson baseline rates, the input hashes, and a
      commitment_hash = sha256 of the canonical JSON of the committed block.

score is refused until now >= T + 7 d + LATENCY_DAYS. The 3-day latency is FROZEN:
      ComCat revises recent solutions, so a week is only scored once its catalogue
      has had three days to settle.

Failure-first (all refusals exit 2; a fetch failure exits 1 and writes nothing):
  * ComCat fetch fails            -> abort, no partial record is ever written.
  * trailing-30-day count outside
    [MIN_TRAIL30, MAX_TRAIL30]    -> refuse to emit, print why.
  * record scored while not due   -> refuse.
  * emitter hash in the log's
    protocol block != this file   -> LOUD WARNING, continue (each record carries the
                                     emitter sha256 of the version that emitted it).

Deterministic: the simulation seed is etas_forecast.SEED, declared and frozen.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import sys

import numpy as np
import pandas as pd

from etas_forecast import (
    CACHE,
    EXP_H,
    HORIZON_DAYS,
    LAT_MAX,
    LAT_MIN,
    LON_MAX,
    LON_MIN,
    MINMAG,
    N_SIMS,
    SEED,
    TRUNC_DAYS,
    download_catalog,
    expected_from_history,
    load_frozen,
    simulate,
)

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(HERE, "results_prospective_socal_log.json")
SELF = os.path.join(HERE, "prospective_socal.py")
ETAS_FORECAST = os.path.join(HERE, "etas_forecast.py")

# ---- frozen protocol constants (one value each, declared before any emission) ----
N_DAYS = int(round(HORIZON_DAYS))     # 7 committed daily bins
LATENCY_DAYS = 3.0                    # FROZEN: ComCat settling time before scoring
TRAIL_DAYS = 365.0                    # PRIMARY baseline window
TRAIL30_DAYS = 30.0                   # sanity-count window
MIN_TRAIL30 = 5                       # refuse to emit below this
MAX_TRAIL30 = 5000                    # refuse to emit above this
EPOCH = pd.Timestamp("2010-01-01", tz="UTC")
LN2 = math.log(2.0)
RATE_FLOOR = 1e-12                    # guards log(0) in the Poisson log-likelihood

BASELINE_A = "poisson_trailing_365d"  # PRIMARY
BASELINE_B = "poisson_train_rate"     # SECONDARY (frozen B-2 training rate)


# ================================================================ small helpers
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical(obj):
    """Canonical JSON bytes: sorted keys, no whitespace, UTF-8."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True).encode("utf-8")


def canonical_hash(obj):
    return hashlib.sha256(canonical(obj)).hexdigest()


def poisson_ll(counts, rates):
    """sum_d [ k_d log(lam_d) - lam_d - log(k_d!) ], one day per bin."""
    ll = 0.0
    for k, lam in zip(counts, rates):
        lam = max(float(lam), RATE_FLOOR)
        ll += float(k) * math.log(lam) - lam - math.lgamma(float(k) + 1.0)
    return ll


def parse_utc(s):
    t = pd.Timestamp(s)
    t = t.tz_localize("UTC") if t.tz is None else t.tz_convert("UTC")
    return t


def floor_second(t):
    return t.floor("s")


# ============================================================== catalogue access
def load_catalogue(now, refresh=False):
    """Refresh the ComCat cache through `now` and return the earthquake rows.

    Raises whatever etas_forecast.download_catalog raises. On failure NOTHING in
    this module has been written: the log is only touched at the very end of emit,
    and download_catalog itself writes its cache only after every window succeeded.
    """
    cat, dl = download_catalog(now, refresh=refresh)
    cat = cat[cat["type"] == "earthquake"].reset_index(drop=True)
    cat = cat[cat["time"] <= now].reset_index(drop=True)
    return cat, dl


def daily_counts(cat, t0, n_days=N_DAYS):
    """Observed events in the n_days one-day bins [t0+d, t0+d+1)."""
    out = []
    for d in range(n_days):
        a = t0 + pd.Timedelta(days=d)
        b = t0 + pd.Timedelta(days=d + 1)
        out.append(int(((cat["time"] >= a) & (cat["time"] < b)).sum()))
    return out


# ==================================================================== protocol
def protocol_block():
    return {
        "entry": "EQ-25 Phase 4 (Popper P7-26(5) rank 1)",
        "artifact": ("weekly hash-committed PROSPECTIVE SoCal forecast log "
                     "(commits AND scores; scoring never edits a commitment)"),
        "version": "v1: frozen ETAS vs Poisson and nothing else",
        "opened_utc": None,          # filled at creation
        "append_only": ("Records are APPENDED to `records`. A record's `commitment` "
                        "block is never edited. Scoring adds a sibling `scoring` "
                        "block and flips `scored`; `commitment_hash` must still "
                        "verify before any score is written."),
        "model": {
            "name": "frozen temporal ETAS (EXP-H / B-2)",
            "refit": False,
            "params_source": ("results_exp_h.json :: train_fit.frozen_params "
                              "(+ b_value_train_aki); fit on SCSN 1982-2010, scored "
                              "once walk-forward on 2010-2018"),
            "intensity": ("lambda(t) = mu + sum_{t_i<t} K*10^(alpha*(M_i-M0))"
                          "*(t-t_i+c)^(-p), t in days"),
            "forecast": ("N_SIMS ETAS cluster simulations over [T, T+7d), the seven "
                         "committed numbers are the per-day MEAN simulated count of "
                         "M>=2.5 events (etas_forecast.simulate)"),
            "analytic_floor": ("per-day no-future-triggering lower bound from "
                               "etas_forecast.expected_from_history; recorded but NOT "
                               "the committed forecast"),
            "history_truncation_days": TRUNC_DAYS,
            "seed": SEED,
        },
        "domain": {
            "catalog": "USGS ComCat FDSN event/1/query, format=csv, type==earthquake",
            "box": {"lat": [LAT_MIN, LAT_MAX], "lon": [LON_MIN, LON_MAX]},
            "minmagnitude": MINMAG,
            "cache_file": os.path.basename(str(CACHE)),
            "coverage_note": ("temporal only, whole box. No spatial resolution of any "
                              "kind is committed or claimed."),
        },
        "baselines": {
            "primary": BASELINE_A,
            "primary_rule": ("observed M>=2.5 count in [T-365d, T) divided by 365, "
                             "per day, constant across the seven bins. Uses only data "
                             "strictly before T."),
            "secondary": BASELINE_B,
            "secondary_rule": ("results_exp_h.json :: "
                               "test.baselines.poisson_train_rate.rate_per_day, the "
                               "frozen B-2 training rate. Uses no data after 2010."),
        },
        "scoring": {
            "due_rule": (f"now >= T + {HORIZON_DAYS:.0f} d + {LATENCY_DAYS:.0f} d"),
            "latency_days": LATENCY_DAYS,
            "latency_reason": ("FROZEN 3-day ComCat settling latency. Recent solutions "
                               "are preliminary and are revised; a week is scored only "
                               "after its catalogue has had three days to settle. The "
                               "rule is frozen so it cannot be tuned per week."),
            "statistic": ("Poisson log-likelihood of the seven observed daily counts "
                          "under the committed ETAS daily expectations vs under each "
                          "constant baseline rate; information gain in bits per event, "
                          "(ll_model - ll_baseline)/(n_events*ln2), the engine/score.py "
                          "`information_gain` definition."),
            "zero_event_rule": ("if n_events == 0 the bits/event are recorded as null "
                                "and the total log-likelihood difference is recorded in "
                                "NATS instead. No division by zero occurs."),
            "cherry_picking": ("every emitted record is scored when due. Records are "
                               "never withdrawn, and the cumulative summary is over ALL "
                               "scored records."),
        },
        "refusals": {
            "exit_2": ["trailing-30-day count outside "
                       f"[{MIN_TRAIL30}, {MAX_TRAIL30}]",
                       "score requested on a record that is not due",
                       "commitment_hash does not verify"],
            "exit_nonzero_no_write": ["ComCat fetch failure (no partial record is "
                                      "ever written)"],
            "warn_and_continue": ["the log's protocol emitter sha256 differs from the "
                                  "current prospective_socal.py (each record carries "
                                  "the emitter hash of the version that emitted it)"],
        },
        "licence": {
            "what_a_positive_licenses": ("B-2 prospectively: a frozen temporal ETAS "
                                         "beats a Poisson rate on unseen SoCal M>=2.5 "
                                         "seismicity, whole box. Nothing else."),
            "what_it_does_not_license": ["anything spatial", "anything about magnitude "
                                         "beyond the M>=2.5 counting threshold",
                                         "any claim of earthquake prediction",
                                         "any beyond-ETAS mechanism"],
        },
        "provenance": {
            "emitter": "prospective_socal.py",
            "emitter_sha256": sha256_file(SELF),
            "etas_forecast_sha256": sha256_file(ETAS_FORECAST),
            "results_exp_h_sha256": sha256_file(EXP_H),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
    }


def empty_log(now):
    proto = protocol_block()
    proto["opened_utc"] = now.isoformat()
    return {"protocol": proto, "records": [],
            "summary": {"weeks_scored": 0, "total_events": 0,
                        f"total_bits_vs_{BASELINE_A}": 0.0,
                        f"total_bits_vs_{BASELINE_B}": 0.0,
                        "note": "cumulative over ALL scored records; no cherry-picking"}}


def read_log():
    if not os.path.exists(LOG):
        return None
    with open(LOG, encoding="utf-8") as fh:
        return json.load(fh)


def write_log(doc):
    with open(LOG, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1)
        fh.write("\n")


def check_emitter_hash(doc):
    """LOUD warning if the log was opened by a different emitter version."""
    logged = doc.get("protocol", {}).get("provenance", {}).get("emitter_sha256")
    current = sha256_file(SELF)
    if logged and logged != current:
        print("!" * 78)
        print("WARNING: prospective_socal.py has CHANGED since this log was opened.")
        print(f"  log protocol emitter_sha256 = {logged}")
        print(f"  current file      sha256    = {current}")
        print("  Continuing. Each record carries the emitter sha256 of the version")
        print("  that emitted it; the record hash is what is authoritative.")
        print("!" * 78)
        return False
    return True


def recompute_summary(doc):
    weeks = 0
    n_ev = 0
    bits_a = 0.0
    bits_b = 0.0
    for r in doc["records"]:
        if not r.get("scored"):
            continue
        s = r["scoring"]
        weeks += 1
        n_ev += int(s["n_events"])
        bits_a += float(s["bits_total_vs_A"])
        bits_b += float(s["bits_total_vs_B"])
    doc["summary"] = {
        "weeks_scored": weeks,
        "total_events": n_ev,
        f"total_bits_vs_{BASELINE_A}": bits_a,
        f"total_bits_vs_{BASELINE_B}": bits_b,
        f"cumulative_bits_per_event_vs_{BASELINE_A}": (bits_a / n_ev) if n_ev else None,
        f"cumulative_bits_per_event_vs_{BASELINE_B}": (bits_b / n_ev) if n_ev else None,
        "note": "cumulative over ALL scored records; no cherry-picking",
    }


# ======================================================================== emit
def do_emit(args):
    now = floor_second(parse_utc(args.now) if args.now else pd.Timestamp.now(tz="UTC"))

    par, exp_h = load_frozen()
    rate_b = float(exp_h["test"]["baselines"]["poisson_train_rate"]["rate_per_day"])

    print("=" * 78)
    print("EQ-25 Phase 4  PROSPECTIVE SoCal weekly commitment (frozen ETAS vs Poisson)")
    print(f"  T (emission / forecast origin) = {now.isoformat()}")
    print("=" * 78)

    # ---- failure site 1: the ComCat fetch. Nothing is written on failure.
    try:
        cat, dl = load_catalogue(now, refresh=args.refresh)
    except BaseException as exc:                       # noqa: BLE001
        print(f"ABORT: ComCat fetch failed: {exc}")
        print("Nothing was written. The log on disk is unchanged.")
        return 1

    # ---- failure site 2: the trailing-30-day sanity gate.
    t30 = now - pd.Timedelta(days=TRAIL30_DAYS)
    n30 = int(((cat["time"] >= t30) & (cat["time"] < now)).sum())
    if not (MIN_TRAIL30 <= n30 <= MAX_TRAIL30):
        print(f"REFUSE: trailing-{TRAIL30_DAYS:.0f}-day M>={MINMAG} count = {n30}, "
              f"outside the frozen admissible band [{MIN_TRAIL30}, {MAX_TRAIL30}].")
        print("  SoCal M>=2.5 runs at roughly 20/week outside sequences; a big sequence")
        print("  can push it far higher, hence the wide upper bound. A count outside")
        print("  the band means the catalogue is broken or the region is doing")
        print("  something the frozen parameters were never scored on. Not emitting.")
        return 2

    t365 = now - pd.Timedelta(days=TRAIL_DAYS)
    n365 = int(((cat["time"] >= t365) & (cat["time"] < now)).sum())
    rate_a = n365 / TRAIL_DAYS

    # ---- the forecast itself: frozen ETAS, imported machinery, no refit.
    t_hist = (cat["time"] - EPOCH).dt.total_seconds().to_numpy() / 86400.0
    m_hist = cat["mag"].to_numpy(dtype=float)
    t_now = (now - EPOCH).total_seconds() / 86400.0

    counts, mmax, capped, binned = simulate(
        t_now, t_hist, m_hist, par, n_sims=args.nsims, T=float(N_DAYS), seed=SEED)
    etas_daily = [float(binned["day_count"][:, d].mean()) for d in range(N_DAYS)]

    cum = [0.0]
    for d in range(1, N_DAYS + 1):
        lam, _, _ = expected_from_history(t_now, t_hist, m_hist, par, T=float(d))
        cum.append(float(lam))
    floor_daily = [cum[d + 1] - cum[d] for d in range(N_DAYS)]

    cat_sha = sha256_file(str(CACHE))
    commitment = {
        "record_id": "W-" + now.strftime("%Y%m%dT%H%M%SZ"),
        "T_utc": now.isoformat(),
        "horizon_days": float(N_DAYS),
        "bins": "seven one-day bins [T+d, T+d+1), d = 0..6",
        "etas_expected_counts_per_day": etas_daily,
        "etas_expected_total": float(sum(etas_daily)),
        "analytic_floor_counts_per_day_no_future_triggering": floor_daily,
        "n_sims": int(args.nsims),
        "seed": SEED,
        "baselines": {
            BASELINE_A: {"rate_per_day": rate_a, "role": "PRIMARY",
                         "n_events_trailing_365d": n365,
                         "window": [t365.isoformat(), now.isoformat()]},
            BASELINE_B: {"rate_per_day": rate_b, "role": "SECONDARY",
                         "source": ("results_exp_h.json :: "
                                    "test.baselines.poisson_train_rate.rate_per_day")},
        },
        "frozen_params": par,
        "refit": False,
        "catalog": {
            "n_events_used": int(len(cat)),
            "t_first": cat["time"].min().isoformat() if len(cat) else None,
            "t_last": cat["time"].max().isoformat() if len(cat) else None,
            "n_events_trailing_30d": n30,
            "box": {"lat": [LAT_MIN, LAT_MAX], "lon": [LON_MIN, LON_MAX]},
            "minmagnitude": MINMAG,
        },
        "input_sha256": {
            os.path.basename(str(CACHE)): cat_sha,
            "results_exp_h.json": sha256_file(EXP_H),
            "prospective_socal.py": sha256_file(SELF),
            "etas_forecast.py": sha256_file(ETAS_FORECAST),
        },
        "emitted_at_utc": pd.Timestamp.now(tz="UTC").isoformat(),
    }
    record = {
        "record_id": commitment["record_id"],
        "commitment": commitment,
        "commitment_hash": canonical_hash(commitment),
        "scored": False,
        "scoring": None,
    }

    doc = read_log()
    if doc is None:
        doc = empty_log(now)
    else:
        # ---- failure site 4: emitter drift. Loud warning, then continue.
        check_emitter_hash(doc)
    n_before = len(doc["records"])
    if any(r["record_id"] == record["record_id"] for r in doc["records"]):
        print(f"REFUSE: record_id {record['record_id']} already exists in the log. "
              "Records are append-only and never overwritten.")
        return 2
    doc["records"].append(record)
    recompute_summary(doc)
    write_log(doc)

    # ---- bulk invariant
    check = read_log()
    assert len(check["records"]) == n_before + 1, (
        f"BULK INVARIANT VIOLATED: {n_before} records before, "
        f"{len(check['records'])} after one emit")

    print(f"\n[emit] record_id       = {record['record_id']}")
    print(f"[emit] T               = {commitment['T_utc']}")
    print("[emit] ETAS expected M>=2.5 counts per day (committed):")
    for d in range(N_DAYS):
        print(f"         day {d + 1}  E[N] = {etas_daily[d]:8.3f}   "
              f"(analytic floor {floor_daily[d]:7.3f})")
    print(f"[emit] ETAS 7-day total = {sum(etas_daily):.3f}")
    print(f"[emit] baseline A {BASELINE_A} (PRIMARY)   = "
          f"{rate_a:.6f} /day  (n365 = {n365})")
    print(f"[emit] baseline B {BASELINE_B} (SECONDARY) = {rate_b:.6f} /day")
    print(f"[emit] trailing-30-day count = {n30}")
    print(f"[emit] commitment_hash = {record['commitment_hash']}")
    print(f"[emit] records in log  = {n_before} -> {len(check['records'])}")
    print(f"[emit] wrote {LOG}")
    return 0


# ======================================================================= score
def is_due(commitment, now):
    t0 = parse_utc(commitment["T_utc"])
    due_at = t0 + pd.Timedelta(days=commitment["horizon_days"] + LATENCY_DAYS)
    return now >= due_at, due_at


def score_one(record, cat, now, cat_sha):
    """Score one record in place. Returns the scoring block."""
    c = record["commitment"]
    t0 = parse_utc(c["T_utc"])
    obs = daily_counts(cat, t0, N_DAYS)
    n_ev = int(sum(obs))

    etas = c["etas_expected_counts_per_day"]
    rate_a = float(c["baselines"][BASELINE_A]["rate_per_day"])
    rate_b = float(c["baselines"][BASELINE_B]["rate_per_day"])

    ll_etas = poisson_ll(obs, etas)
    ll_a = poisson_ll(obs, [rate_a] * N_DAYS)
    ll_b = poisson_ll(obs, [rate_b] * N_DAYS)

    d_a = ll_etas - ll_a
    d_b = ll_etas - ll_b
    scoring = {
        "scored_at_utc": pd.Timestamp.now(tz="UTC").isoformat(),
        "observed_counts_per_day": obs,
        "n_events": n_ev,
        "ll_etas": ll_etas,
        "ll_A": ll_a,
        "ll_B": ll_b,
        "baseline_A": BASELINE_A,
        "baseline_B": BASELINE_B,
        "delta_ll_nats_vs_A": d_a,
        "delta_ll_nats_vs_B": d_b,
        "bits_total_vs_A": d_a / LN2,
        "bits_total_vs_B": d_b / LN2,
        "bits_per_event_vs_A": (d_a / (n_ev * LN2)) if n_ev > 0 else None,
        "bits_per_event_vs_B": (d_b / (n_ev * LN2)) if n_ev > 0 else None,
        "zero_event_note": (None if n_ev > 0 else
                            "n_events == 0: bits/event undefined and recorded as null; "
                            "the total log-likelihood difference in NATS "
                            "(delta_ll_nats_vs_A / _B) is the whole result for this "
                            "week. No division by zero was performed."),
        "catalog_sha256_at_scoring": cat_sha,
        "scored_by_sha256": sha256_file(SELF),
    }
    record["scoring"] = scoring
    record["scored"] = True
    return scoring


def do_score(args):
    now = floor_second(parse_utc(args.now) if args.now else pd.Timestamp.now(tz="UTC"))
    doc = read_log()
    if doc is None:
        print(f"REFUSE: no log at {LOG}. Emit something first.")
        return 2
    check_emitter_hash(doc)

    if args.record:
        targets = [r for r in doc["records"] if r["record_id"] == args.record]
        if not targets:
            print(f"REFUSE: no record with id {args.record}")
            return 2
    elif args.all_due:
        targets = [r for r in doc["records"] if not r.get("scored")
                   and is_due(r["commitment"], now)[0]]
        if not targets:
            print("[score] nothing due. Log unchanged.")
            return 0
    else:
        print("REFUSE: score requires --record ID or --all-due")
        return 2

    # ---- failure site 3: not due yet.
    for r in targets:
        due, due_at = is_due(r["commitment"], now)
        if not due:
            print(f"REFUSE: record {r['record_id']} is not due. "
                  f"T = {r['commitment']['T_utc']}, "
                  f"due at {due_at.isoformat()} "
                  f"(T + {r['commitment']['horizon_days']:.0f} d + "
                  f"{LATENCY_DAYS:.0f} d frozen latency), now = {now.isoformat()}.")
            return 2
        if r.get("scored"):
            print(f"REFUSE: record {r['record_id']} is already scored. "
                  "Scored records are never rescored or edited.")
            return 2
        # ---- commitment hash must still verify BEFORE anything is scored.
        recomputed = canonical_hash(r["commitment"])
        if recomputed != r["commitment_hash"]:
            print(f"REFUSE: commitment_hash does not verify for {r['record_id']}.")
            print(f"  stored     = {r['commitment_hash']}")
            print(f"  recomputed = {recomputed}")
            print("  The commitment has been altered. Refusing to score it.")
            return 2

    # ---- failure site 1 again: the fetch. Nothing is written on failure.
    try:
        cat, dl = load_catalogue(now, refresh=args.refresh)
    except BaseException as exc:                       # noqa: BLE001
        print(f"ABORT: ComCat fetch failed: {exc}")
        print("Nothing was written. The log on disk is unchanged.")
        return 1
    cat_sha = sha256_file(str(CACHE))

    for r in targets:
        s = score_one(r, cat, now, cat_sha)
        print("-" * 78)
        print(f"[score] {r['record_id']}  T = {r['commitment']['T_utc']}")
        print(f"[score] observed per day = {s['observed_counts_per_day']}  "
              f"n_events = {s['n_events']}")
        print(f"[score] committed  ETAS  = "
              f"{[round(x, 2) for x in r['commitment']['etas_expected_counts_per_day']]}")
        print(f"[score] ll_etas = {s['ll_etas']:.4f}  ll_A = {s['ll_A']:.4f}  "
              f"ll_B = {s['ll_B']:.4f}")
        if s["n_events"] > 0:
            print(f"[score] bits/event vs A ({BASELINE_A}, PRIMARY)   = "
                  f"{s['bits_per_event_vs_A']:+.4f}")
            print(f"[score] bits/event vs B ({BASELINE_B}, SECONDARY) = "
                  f"{s['bits_per_event_vs_B']:+.4f}")
        else:
            print("[score] n_events = 0: bits/event null; "
                  f"delta LL nats vs A = {s['delta_ll_nats_vs_A']:+.4f}, "
                  f"vs B = {s['delta_ll_nats_vs_B']:+.4f}")

    recompute_summary(doc)
    write_log(doc)
    su = doc["summary"]
    print("-" * 78)
    print(f"[cumulative] weeks_scored = {su['weeks_scored']}  "
          f"total_events = {su['total_events']}")
    print(f"[cumulative] total bits vs A = {su['total_bits_vs_' + BASELINE_A]:+.4f}  "
          f"vs B = {su['total_bits_vs_' + BASELINE_B]:+.4f}")
    print(f"[cumulative] bits/event vs A = "
          f"{su['cumulative_bits_per_event_vs_' + BASELINE_A]}")
    print(f"[cumulative] bits/event vs B = "
          f"{su['cumulative_bits_per_event_vs_' + BASELINE_B]}")
    print(f"[score] wrote {LOG}")
    return 0


# ======================================================================== main
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="EQ-25 Phase 4: weekly hash-committed prospective SoCal forecast")
    sub = ap.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("emit", help="commit one week")
    e.add_argument("--now", default=None, help="ISO UTC override for T")
    e.add_argument("--refresh", action="store_true")
    e.add_argument("--nsims", type=int, default=N_SIMS)

    s = sub.add_parser("score", help="score a due week")
    s.add_argument("--record", default=None)
    s.add_argument("--all-due", dest="all_due", action="store_true")
    s.add_argument("--now", default=None, help="ISO UTC override for the clock")
    s.add_argument("--refresh", action="store_true")

    args = ap.parse_args(argv)
    if args.cmd == "emit":
        return do_emit(args)
    return do_score(args)


if __name__ == "__main__":
    sys.exit(main())
