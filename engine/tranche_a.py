"""TRANCHE A driver -- the first REAL-DATA v2 exploration run (HYPOTHESIS_LEDGER.md §P7-2).

WHAT THIS RUNS, and what it deliberately does NOT re-run
--------------------------------------------------------
§P7-2(b)'s overlap reconciliation exists so that nothing runs twice under two
names. Read against it, Tranche A on 2026-08-12 has already discharged most of
itself, and this driver runs only the remainder:

  ALREADY DONE, NOT RE-RUN HERE
    F4-58  the deflation measurement -- committed (results_f4_58_vif.json), and its
           per-feature output is CONSUMED here as the operative floor input.
    F9-19  planted recovery per band AND aggregation -- this IS G-M1 arm (ii)
           extended, not a new gate; it ran inside the v2 acceptance gate
           (results_gate_r1.json). One execution, one artifact, one ledger object.
    R1     null DATA through the real pipeline -- 30 true-null ETAS catalogues,
           also in the gate.

  WHAT IS LEFT, AND WHAT THIS DRIVER DOES
    F9-20 arm 1   null FEATURES on the REAL data. R1 and F9-20 are TWO ARMS OF ONE
                  GATE and neither discharges the other: R1 catches pipeline and
                  dependence miscalibration, F9-20 catches feature-construction
                  leakage and residual real-data structure. This is the whole of
                  Tranche A's price: 23 controls x 31 lags = 713 declared tests.
    F9-17 / R3    winner's-curse handling on every reported amplitude.
    F10-25        the survivor ratio, real arm against control arm.
    F10-24        rank stability under DATA resampling.
    R4            rank stability under RESEEDING. Adjacent to F10-24, not identical:
                  one asks whether the DATA moves the ranking, the other whether the
                  RNG does. Both required, both reported, never conflated.

PRICE, AND WHERE IT IS DECLARED
-------------------------------
The NEW price of this tranche is 713 -- Kepler's integer, confirmed by Popper at
§P7-2(a) -- and every other Tranche A item is priced at 0 because it reports a
property of the run rather than making a rejection. The session's BH DENOMINATOR is
larger than 713 and that is not a contradiction: the real arm must sit in the same
declared vector or F10-25's ratio compares two arms facing two thresholds, which
§P7-2(a) says makes it meaningless. Its 550 already-declared tests therefore
re-occupy 550 slots in this session's denominator -- already-declared is not free --
but they are NOT new scope and are not presented as any. m = 713 + 550 = 1263.

THE AUDITS ARE PRICED AT 0 AND ARE SANDBOXED
--------------------------------------------
The four stability re-runs make no rejection and enter no BH vector, so they write
to their own session directories under a scratch root and append to their own
ledger. Only the ONE priced session touches engine/out/mine and
engine/EXPLORE_COUNT.jsonl, and only the parent ever appends to either.

NOTHING HERE IS EVIDENCE. This is a generator. No output may be entered for or
against any ledger entry.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import types

import numpy as np

from . import mine as M, mine_session as ms, splits

TRANCHE_A_ID = "TRANCHE-A-F9-20-arm1"
STRATA_FILE = os.path.join("engine", "configs", "strata_tranche_a.json")
AUDIT_ROOT = os.path.join("engine", "out", "tranche_a")
RESULTS_JSON = "results_tranche_a.json"

# The declared arithmetic, asserted in code before anything runs.
N_PRICED_NEW = M.F9_20_N_DECLARED_TESTS            # 713 -- Tranche A's whole price
N_ALREADY_DECLARED = 550                           # K-089-R tranche 1, 20260812T004857
N_DECLARED_VECTOR = N_PRICED_NEW + N_ALREADY_DECLARED

# Audit seeds, DECLARED HERE and not chosen after seeing a result (S-9).
RESEED_SEEDS = (20260812, 20260813)
RESAMPLE_SEEDS = (90260812, 90260813)
# Moving-block length for the F10-24 data resample. Long against the 30 d maximum
# lag so that the within-block lag structure survives the resample; the fraction of
# positions whose lag history crosses a block boundary is stated in the output
# rather than hoped away.
RESAMPLE_BLOCK_DAYS = 365
MAX_LAG = 30


def build_args(seed=20260811, jobs=8):
    """The frozen Tranche A config inputs. One declared value each, no alternatives."""
    return types.SimpleNamespace(
        mag_target=4.5, dlat=1.0, dlon=1.0, explore_frac=0.70,
        data_dir="data/comcat_world",
        # Downloads OFF: the real arm this tranche calibrates is the already-declared
        # K-089-R tranche-1 sweep, which ran with downloads off. Turning them on here
        # would change the real arm and the comparison would not be the one priced.
        no_download=True,
        seed=int(seed), tranche1=True, controls=True, strata=STRATA_FILE,
        ladder=False, gpd=False, regsum=False, regions=False, jobs=int(jobs))


def frozen_config(seed=20260811):
    cfg = ms.build_config(build_args(seed=seed), ms.OVERNIGHT)
    assert cfg["controls"]["n_declared_tests"] == N_PRICED_NEW
    assert cfg["strata"]["m"] == N_DECLARED_VECTOR, (
        f"partition declares m = {cfg['strata']['m']}, expected "
        f"{N_DECLARED_VECTOR} = {N_PRICED_NEW} priced + {N_ALREADY_DECLARED} "
        f"already-declared")
    return cfg


# ------------------------------------------------------- F10-24 data resample --
def block_resample_index(n, block, rng):
    """Moving-block bootstrap index over [0, n), blocks of `block` days.

    Contiguous blocks, so within a block every lag relation in the original record
    survives untouched; only the joins are artificial. `block_boundary_fraction`
    below reports exactly how much of the record that costs.
    """
    block = int(min(block, n))
    starts = rng.integers(0, n - block + 1, size=int(np.ceil(n / block)))
    idx = np.concatenate([np.arange(s, s + block) for s in starts])[:n]
    return idx.astype(np.int64)


def resample_prepared(prepared, seed, block=RESAMPLE_BLOCK_DAYS, max_lag=MAX_LAG):
    """A block-bootstrap resample of the MINING WINDOW, jointly across every series.

    Target, ETAS offset, every feature and every event mark are re-indexed by the
    SAME day index, so the alignment that the miner is testing is preserved inside
    each block and destroyed only at the joins. That is the point: F10-24 asks
    whether the DATA moves the ranking, which requires a resample that could still
    carry the effect, not a null.
    """
    (ctx, base, y, window, counts, offset, marks, feats, dl_log, t0) = prepared
    rng = np.random.default_rng(
        M.seed_sequence_for(M.test_key(int(seed), TRANCHE_A_ID, "f10_24_resample")))
    n_win = int(counts.size)
    idx = block_resample_index(n_win, block, rng)
    # An extended index so that lags 0..max_lag never run off the head of the record.
    idx_ext = np.concatenate([idx[:max_lag], idx])
    n_ext = idx_ext.size
    new_window = slice(max_lag, n_ext)

    counts_r = np.concatenate([counts[idx[:max_lag]], counts[idx]])[max_lag:]
    offset_r = np.concatenate([offset[idx[:max_lag]], offset[idx]])[max_lag:]

    feats_r = []
    for f in feats:
        vals_win = f.values[window]
        g = M.Feature(f.name, f.family, f.kind, vals_win[idx_ext], f.describe,
                      lags=f.lags, causality_exempt=f.causality_exempt,
                      period_hint=f.period_hint,
                      control=getattr(f, "control", False),
                      periodic_override=bool(f.periodic),
                      control_of=getattr(f, "control_of", None),
                      # the resample must not silently re-tune the null: the
                      # bootstrap block is the ORIGINAL feature's, carried.
                      block_days_override=f.block_days)
        # Attributes set AFTER construction do not survive a rebuild. `subdaily_only`
        # is one: lose it and the resampled session hands a zero-variance diurnal
        # column to the count path, which raises. Carried explicitly rather than
        # hoped for.
        if getattr(f, "subdaily_only", False):
            g.subdaily_only = True
        feats_r.append(g)

    # Marks follow their day. Each resampled day j carries the events of original
    # window-day idx[j]; an original day drawn twice contributes its events twice,
    # which is what a bootstrap of a point process on a day lattice means.
    day_local = marks["day"] - window.start
    order = np.argsort(day_local, kind="stable")
    day_sorted = day_local[order]
    lo = np.searchsorted(day_sorted, idx, side="left")
    hi = np.searchsorted(day_sorted, idx, side="right")
    take, newday = [], []
    for j in range(idx.size):
        if hi[j] > lo[j]:
            take.append(order[lo[j]:hi[j]])
            newday.append(np.full(hi[j] - lo[j], j, dtype=np.int64))
    take = (np.concatenate(take) if take else np.zeros(0, dtype=np.int64))
    newday = (np.concatenate(newday) if newday else np.zeros(0, dtype=np.int64))
    # `day` is expressed on the NEW window's own axis: run() recovers the
    # in-window ordinal as `marks["day"] - window.start`, and the new window
    # starts at `max_lag` because the extended head carries the lag history.
    marks_r = {}
    for k, v in marks.items():
        v = np.asarray(v)
        if k == "day":
            marks_r[k] = newday + new_window.start
        elif k == "day_float":
            frac = v[take] - np.floor(v[take])
            marks_r[k] = newday.astype(np.float64) + frac + new_window.start
        else:
            marks_r[k] = v[take]

    meta = {
        "block_days": int(block),
        "n_blocks": int(np.ceil(n_win / min(block, n_win))),
        "n_window_days": n_win,
        "block_boundary_fraction": float(min(max_lag, block) / float(block)),
        "block_boundary_note": (
            "the share of positions whose 0..%d d lag history crosses a block join "
            "and is therefore artificial. Stated, not hidden: it is the price of "
            "resampling a series that is being tested at a lag." % max_lag),
        "n_distinct_days_drawn": int(np.unique(idx).size),
        "n_events_after_resample": int(marks_r["day"].size),
        "n_events_before": int(marks["day"].size),
    }
    return (ctx, base, y, new_window, counts_r, offset_r, marks_r, feats_r,
            dl_log, t0), meta


# --------------------------------------------------------------------- driver --
def _load_tests(session_dir):
    with open(os.path.join(session_dir, "checkpoint.json"), "r",
              encoding="utf-8") as fh:
        return json.load(fh)["tests"]


def run_audit(kind, seed, prepared_main, jobs, root=AUDIT_ROOT, verbose=True):
    """One stability re-run. PRICED AT 0: sandboxed session dir, sandboxed ledger."""
    os.makedirs(root, exist_ok=True)
    tag = f"{kind}_{seed}"
    sess = os.path.join(root, "session_audit_" + tag)
    ledger = os.path.join(root, "AUDIT_EXPLORE_COUNT.jsonl")
    meta = None
    if kind == "reseed":
        cfg = frozen_config(seed=seed)
        prepared = prepared_main            # same data, different master seed
    elif kind == "resample":
        cfg = frozen_config(seed=20260811)  # same seed, different data
        prepared, meta = resample_prepared(prepared_main, seed)
    else:
        raise ValueError(kind)
    print(f"\n=== AUDIT ({kind}, seed {seed}) -> {sess} "
          f"[PRICED 0, sandboxed ledger {ledger}] ===", flush=True)
    out = ms.run(cfg, verbose=verbose, resume=False, session_dir=sess, jobs=jobs,
                 ledger_path=ledger, prepared=prepared)
    return {"kind": kind, "seed": int(seed), "session_dir": sess,
            "resample_meta": meta, "n_tests": out["n_tests"],
            "n_pass": out["n_pass"], "elapsed": out["elapsed"]}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--session-dir", default=None)
    ap.add_argument("--skip-audits", action="store_true",
                    help="run only the ONE priced session")
    ap.add_argument("--audit-jobs", type=int, default=None)
    a = ap.parse_args(argv)

    t_open = _dt.datetime.now(_dt.timezone.utc).isoformat()
    cfg = frozen_config()
    print("=" * 78)
    print(f"TRANCHE A -- {M.F9_20_LABEL}")
    print("=" * 78)
    print(f"declared price, NEW      : {N_PRICED_NEW} "
          f"(23 controls x {len(M.F9_20_LAGS)} lags; §P7-2(a))")
    print(f"already declared, re-run : {N_ALREADY_DECLARED} "
          f"(K-089-R tranche 1; occupies its slots, is NOT new scope)")
    print(f"BH denominator m         : {N_DECLARED_VECTOR}")
    print(f"strata file sha256       : {cfg['strata']['sha256']}")
    print(f"config hash              : {splits.config_hash(cfg)}")
    print("=" * 78, flush=True)

    # THE ONE PRICED SESSION. Real ledger, real session root.
    prepared = ms.prepare(cfg, verbose=True)
    main_out = ms.run(cfg, verbose=True, resume=False,
                      session_dir=a.session_dir, jobs=a.jobs, prepared=prepared)
    main_tests = _load_tests(main_out["session_dir"])

    audits, stability = [], {}
    if not a.skip_audits:
        jobs = a.audit_jobs or a.jobs
        for s in RESEED_SEEDS:
            audits.append(run_audit("reseed", s, prepared, jobs))
        for s in RESAMPLE_SEEDS:
            audits.append(run_audit("resample", s, prepared, jobs))

        rs = [x for x in audits if x["kind"] == "reseed"]
        rr = [x for x in audits if x["kind"] == "resample"]
        stability["R4_reseeding"] = [
            ms.rank_stability(main_tests, _load_tests(rs[0]["session_dir"]),
                              label=f"R4: master seed {cfg['seed']} vs "
                                    f"{rs[0]['seed']} (same data)"),
            ms.rank_stability(_load_tests(rs[0]["session_dir"]),
                              _load_tests(rs[1]["session_dir"]),
                              label=f"R4: master seed {rs[0]['seed']} vs "
                                    f"{rs[1]['seed']} (same data)"),
        ]
        stability["F10_24_data_resampling"] = [
            ms.rank_stability(main_tests, _load_tests(rr[0]["session_dir"]),
                              label=f"F10-24: original window vs block resample "
                                    f"{rr[0]['seed']} (same seed)"),
            ms.rank_stability(_load_tests(rr[0]["session_dir"]),
                              _load_tests(rr[1]["session_dir"]),
                              label=f"F10-24: block resample {rr[0]['seed']} vs "
                                    f"{rr[1]['seed']} (same seed)"),
        ]
        stability["never_conflate"] = (
            "§P7-2(b): R4 and F10-24 are ADJACENT AND BOTH REQUIRED. R4 asks "
            "whether the RNG moves the ranking; F10-24 asks whether the DATA "
            "does. Different failure modes, labelled so they are never conflated.")

    with open(os.path.join(main_out["session_dir"], "checkpoint.json"),
              encoding="utf-8") as fh:
        st = json.load(fh)
    payload = {
        "id": TRANCHE_A_ID,
        "title": "EQ-24 v2 Tranche A -- F9-20 arm 1 on real data, with R3, "
                 "F10-25, F10-24 and R4",
        "banner": M.GENERATOR_NOT_EVIDENCE,
        "standing_warning_eq24": M.MIGNAN_BROCCARDO,
        "opened_utc": t_open,
        "closed_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "count_reconciliation": {
            "priced_new_this_tranche": N_PRICED_NEW,
            "priced_new_source": ("MINING_CATALOG.md F9-20 (~23 features x 31 lags "
                                  "= 713), CONFIRMED at HYPOTHESIS_LEDGER.md "
                                  "§P7-2(a) as the declared integer"),
            "already_declared_reoccupying_slots": N_ALREADY_DECLARED,
            "bh_denominator_m": N_DECLARED_VECTOR,
            "priced_at_zero": {
                "F4-58": "ALREADY DONE and committed; consumed here as floor input",
                "F9-19": "IS G-M1 arm (ii)-extended; ALREADY RUN in the v2 gate",
                "R1": "null DATA arm; ALREADY RUN in the v2 gate",
                "F9-17/R3": "reports a property of the run; no rejection",
                "F10-24": "audit; demote-only; no rejection",
                "R4": "audit; demote-only; no rejection",
                "F10-25": "a ratio of two counts already in the vector",
            },
        },
        "config": cfg,
        "config_hash": splits.config_hash(cfg),
        "strata_file": STRATA_FILE,
        "strata_sha256": cfg["strata"]["sha256"],
        "priced_session": {k: v for k, v in main_out.items() if k != "report"},
        "bh": st.get("bh"),
        "max_statistic": st.get("max_statistic"),
        "tranche_a": st.get("tranche_a"),
        "p_method_census": st.get("p_method_census"),
        "audits_priced_zero": audits,
        "rank_stability": stability,
        "what_none_of_this_licenses": (
            "nothing. This is a GENERATOR on the exploration window under a v1 ETAS "
            "baseline, blind through the diurnal and semidiurnal band by two exact "
            "zeros, blind to the second phase moment, and forbidden from being "
            "entered for or against any ledger entry. A survivor in the control "
            "stratum is a measured false positive; a survivor in the real stratum "
            "is a hypothesis stub and nothing more."),
    }
    with open(RESULTS_JSON, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1, default=float)
    print(f"\nwrote {RESULTS_JSON}")

    if stability:
        rep = os.path.join(main_out["session_dir"], "report.md")
        with open(rep, "a", encoding="utf-8") as fh:
            fh.write("\n\n## R4 (reseeding) and F10-24 (data resampling) -- "
                     "rank stability of the top-100\n\n")
            fh.write("*Both are AUDITS, priced at 0, run in sandboxed session "
                     "directories against a sandboxed ledger. " +
                     stability["never_conflate"] + "*\n\n")
            fh.write("| audit | comparison | top-100 Spearman rho | set overlap | "
                     "verdict |\n| --- | --- | ---: | ---: | --- |\n")
            for grp in ("R4_reseeding", "F10_24_data_resampling"):
                for r in stability[grp]:
                    rho = ("n/a" if r["spearman_rho"] is None
                           else f"{r['spearman_rho']:.4f}")
                    fh.write(f"| `{grp}` | {r['label']} | {rho} | "
                             f"{r.get('top_k_set_overlap')}/{r['n_top']} | "
                             f"{r['verdict']} |\n")
            fh.write("\n")
        print(f"appended the stability banner to {rep}")
    return payload


if __name__ == "__main__":
    main()
