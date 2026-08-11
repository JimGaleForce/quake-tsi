"""Session driver for `mine` mode: checkpointed sweep -> report.md + stubs.json.

Split out of engine/mine.py (which holds the statistics) so the orchestration --
resume logic, task ordering, ledger accounting and report writing -- reads in one
sitting. Everything this module writes carries the generator-not-evidence banner.
"""

from __future__ import annotations

import datetime as _dt
import json
import math
import os
import time

import numpy as np

from . import (baseline as bl, design, mine as M, splits, __version__)

QUICK = {
    "n_surrogates": 200, "n_periods": 800, "n_peaks": 5,
    "lags": (0, 1, 3, 7, 14, 30), "label": "quick",
}
OVERNIGHT = {
    "n_surrogates": 10000, "n_periods": 3000, "n_peaks": 10,
    "lags": tuple(range(0, 31)), "label": "overnight",
}
DEFAULT = {
    "n_surrogates": 1000, "n_periods": 1500, "n_peaks": 8,
    "lags": (0, 1, 2, 3, 5, 7, 10, 14, 21, 30), "label": "default",
}

PERIOD_MIN, PERIOD_MAX = 2.0, 4000.0


# ------------------------------------------------------------ checkpointing ---
def _cfg_hash(cfg):
    return splits.config_hash(cfg)[:12]


def find_resumable(cfg_hash, root=M.MINE_DIR):
    """The newest incomplete session directory for this exact configuration."""
    if not os.path.isdir(root):
        return None
    best = None
    for name in sorted(os.listdir(root)):
        path = os.path.join(root, name, "checkpoint.json")
        if not (name.startswith("session_") and os.path.exists(path)):
            continue
        try:
            with open(path, "r", encoding="utf-8") as fh:
                ck = json.load(fh)
        except (OSError, json.JSONDecodeError):
            continue
        if ck.get("config_hash") == cfg_hash and not ck.get("complete"):
            best = (os.path.join(root, name), ck)
    return best


class Checkpoint:
    def __init__(self, path, cfg, cfg_hash):
        self.path = path
        self.state = {
            "engine_version": __version__, "kind": "mine",
            "config": cfg, "config_hash": cfg_hash,
            "created": _dt.datetime.now(_dt.timezone.utc).isoformat(),
            "complete": False, "completed_tasks": [], "results": {},
            "banner": M.GENERATOR_NOT_EVIDENCE,
        }

    @classmethod
    def load(cls, path):
        with open(path, "r", encoding="utf-8") as fh:
            state = json.load(fh)
        obj = cls.__new__(cls)
        obj.path = path
        obj.state = state
        return obj

    def done(self, key):
        return key in self.state["results"]

    def put(self, key, value):
        self.state["results"][key] = value
        if key not in self.state["completed_tasks"]:
            self.state["completed_tasks"].append(key)
        self.save()

    def save(self):
        os.makedirs(os.path.dirname(self.path) or ".", exist_ok=True)
        tmp = self.path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(self.state, fh, indent=1)
        os.replace(tmp, self.path)


# ------------------------------------------------------------------- driver ---
def build_config(args, preset):
    return {
        "engine_version": __version__,
        "mode": "mine",
        "preset": preset["label"],
        "n_surrogates": int(preset["n_surrogates"]),
        "n_periods": int(preset["n_periods"]),
        "n_peaks": int(preset["n_peaks"]),
        "lags": list(preset["lags"]),
        "fdr_q": M.FDR_Q,
        "baseline": "etas",
        "mag_target": float(args.mag_target),
        "grid": {"dlat": float(args.dlat), "dlon": float(args.dlon)},
        "explore_frac": float(args.explore_frac),
        "data_dir": args.data_dir.replace("\\", "/"),
        "downloads": bool(not args.no_download),
        "seed": int(args.seed),
    }


def prepare(cfg, verbose=True):
    """Design + ETAS baseline + target series + the full feature list."""
    ctx = design.build_design(
        data_dir=cfg["data_dir"], dlat=cfg["grid"]["dlat"], dlon=cfg["grid"]["dlon"],
        explore_frac=cfg["explore_frac"], verbose=verbose)
    explore, _hold = splits.temporal_split(ctx.n_days, cfg["explore_frac"])
    y = ctx.day_counts(cfg["mag_target"])

    base = bl.EtasV1(verbose=verbose, mag_target=cfg["mag_target"])
    base.fit(ctx, y, explore)
    if verbose:
        for line in base.report():
            print(line)

    burn = int(base.burn_in_days)
    window = slice(burn, explore.stop)
    counts, offset = M.build_target(ctx, base, y, window)
    if verbose:
        print(f"mining window (days) = [{window.start}, {window.stop}) "
              f"= {counts.size} days")
        print(f"  observed events    = {counts.sum():.0f}")
        print(f"  ETAS expectation   = {offset.sum():.1f}")
        print(f"  residual mean/sd   = {(counts-offset).mean():+.4f} / "
              f"{(counts-offset).std():.4f} events/day")

    t0 = _dt.datetime.fromisoformat(str(ctx.meta["t0"]))
    all_marks = M.load_event_marks(ctx, cfg["data_dir"], ctx.meta["mag_floor"])
    in_win = (all_marks["day"] >= window.start) & (all_marks["day"] < window.stop)
    marks = {k: v[in_win] for k, v in all_marks.items()}
    if verbose:
        print(f"  marks in window    = {marks['day'].size} events "
              f"(magnitude + depth)")

    lags = tuple(cfg["lags"])
    feats = M.ephemeris_features(t0, ctx.n_days)
    dl_feats, dl_log = M.download_features(t0, ctx.n_days, lags,
                                           enabled=cfg["downloads"], verbose=verbose)
    feats += dl_feats
    feats += M.catalog_features(ctx, all_marks, lags)
    if verbose:
        by_fam = {}
        for f in feats:
            by_fam.setdefault(f.family, []).append(f.name)
        for fam in sorted(by_fam):
            print(f"  family {fam}: {len(by_fam[fam])} features -> "
                  f"{', '.join(by_fam[fam])}")
    return ctx, base, y, window, counts, offset, marks, feats, dl_log, t0


def run(cfg, verbose=True, resume=True, session_dir=None):
    os.makedirs(M.MINE_DIR, exist_ok=True)
    ch = _cfg_hash(cfg)
    ckpt = None
    if session_dir is None and resume:
        found = find_resumable(ch)
        if found:
            session_dir, _state = found
            ckpt = Checkpoint.load(os.path.join(session_dir, "checkpoint.json"))
            if verbose:
                print(f"RESUMING session {session_dir} "
                      f"({len(ckpt.state['results'])} tasks already complete)")
    if session_dir is None:
        stamp = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
        session_dir = os.path.join(M.MINE_DIR, f"session_{stamp}")
        os.makedirs(session_dir, exist_ok=True)
    if ckpt is None:
        ckpt = Checkpoint(os.path.join(session_dir, "checkpoint.json"), cfg, ch)
        ckpt.save()

    t_start = time.time()
    (ctx, base, y, window, counts, offset, marks, feats,
     dl_log, t0) = prepare(cfg, verbose=verbose)
    ckpt.state["data_log"] = dl_log
    ckpt.state["window"] = [window.start, window.stop]
    ckpt.state["session_dir"] = session_dir.replace("\\", "/")
    ckpt.save()

    rng = np.random.default_rng(cfg["seed"])
    n_surr = int(cfg["n_surrogates"])
    resid = counts - offset
    days = np.arange(counts.size, dtype=np.float64)

    # ---------------------------------------------------- (a) GLM sweep ----
    for f in feats:
        key = f"glm:{f.name}"
        if ckpt.done(key):
            if verbose:
                print(f"  [skip, checkpointed] {key}")
            continue
        t_f = time.time()
        rows = []
        for lag in f.lags:
            X = f.design(window, lag)
            fit = M.glm_fit(X, counts, offset)
            S = M.score_stat_all_shifts(X, counts, offset)
            p_shift, n_used = M.empirical_p(S, n_surr, rng)
            Sb = M.score_stat_block_bootstrap(X, counts, offset, n_surr, rng,
                                              mean_block=f.block_days)
            p_boot = M.bootstrap_p(S[0], Sb)
            # circular shifts are provably powerless for a deterministic cycle, so
            # for periodic features the block bootstrap IS the null; elsewhere both
            # are valid and the more conservative one is reported.
            p = p_boot if f.periodic else max(p_shift, p_boot)
            b = np.asarray(fit["beta"])
            amp = float(np.hypot(*b)) if f.kind == "phase" else float(abs(b[0]))
            rows.append({
                "feature": f.name, "family": f.family, "kind": f.kind, "lag": int(lag),
                "test": "glm_poisson_offset_etas", "df": f.df,
                "beta": fit["beta"], "se": fit["se"], "amplitude_log_rate": amp,
                "pct_rate_modulation": 100.0 * (math.exp(amp) - 1.0),
                "bits_per_event": fit["bits_per_event"],
                "chi2_score": float(S[0]), "p_parametric": M.chi2_sf(S[0], f.df),
                "p_circular_shift": p_shift, "p_block_bootstrap": p_boot,
                "block_days": f.block_days,
                "null": ("block-bootstrap (periodic feature)" if f.periodic
                         else "max(circular-shift, block-bootstrap)"),
                "p_raw": p, "n_surrogates": min(n_used, n_surr),
                "converged": fit["converged"],
            })
        ckpt.put(key, rows)
        if verbose:
            best = min(rows, key=lambda r: r["p_raw"])
            print(f"  glm {f.name:<24s} {len(rows):>2d} lags in "
                  f"{time.time()-t_f:5.1f}s   best p={best['p_raw']:.4g} "
                  f"(lag {best['lag']}, {best['pct_rate_modulation']:+.2f}%/sd)")

    # -------------------------------------------------- (c) mark tests -----
    fe_day = marks["day"] - window.start
    for f in feats:
        key = f"marks:{f.name}"
        if ckpt.done(key):
            if verbose:
                print(f"  [skip, checkpointed] {key}")
            continue
        vals = f.values[window][fe_day]
        rows = []
        for mk in ("mag", "depth"):
            r = M.mark_test(vals, marks[mk], f.kind, n_surr, rng)
            r.update({"feature": f.name, "family": f.family, "kind": f.kind,
                      "mark": mk, "n_events": int(marks[mk].size)})
            rows.append(r)
        ckpt.put(key, rows)
        if verbose:
            print(f"  marks {f.name:<22s} mag p={rows[0]['p_raw']:.4g} "
                  f"depth p={rows[1]['p_raw']:.4g}")

    # ------------------------------------------------ (b) period scan ------
    if not ckpt.done("period_scan"):
        periods = np.exp(np.linspace(math.log(PERIOD_MIN), math.log(PERIOD_MAX),
                                     int(cfg["n_periods"])))
        peaks, meta = M.period_scan(days, resid, periods, n_surr, rng,
                                    n_peaks=int(cfg["n_peaks"]), verbose=verbose)
        lam0 = offset * (counts.sum() / offset.sum())
        for pk in peaks:
            pk["ladder"] = M.harmonic_ladder(counts, lam0, days, pk["period_days"],
                                             max_period=counts.size / 3.0)
            pk["feature"] = f"period_{pk['period_days']:.4g}d"
            pk["family"] = 0
            pk["test"] = "lomb_scargle_peak"
            wp = pk["ladder"]["winning_period_days"]
            ph = np.mod(days, wp) / wp
            b = np.minimum((ph * 8).astype(int), 7)
            Y = np.bincount(b, weights=counts, minlength=8)
            L = np.bincount(b, weights=lam0, minlength=8)
            ratio = Y / np.maximum(L * (Y.sum() / L.sum()), 1e-9)
            pk["pct_rate_modulation"] = float(100.0 * (ratio.max() - ratio.min()) / 2.0)
            pk["fold_ratio_by_phase_bin"] = [round(float(v), 4) for v in ratio]
        ckpt.put("period_scan", {"peaks": peaks, "scan": meta,
                                 "n_trial_periods": int(cfg["n_periods"]),
                                 "period_range_days": [PERIOD_MIN, PERIOD_MAX]})
    elif verbose:
        print("  [skip, checkpointed] period_scan")

    # --------------------------------------------------- multiplicity ------
    tests = []
    for f in feats:
        tests += ckpt.state["results"][f"glm:{f.name}"]
        tests += ckpt.state["results"][f"marks:{f.name}"]
    tests += ckpt.state["results"]["period_scan"]["peaks"]
    p = np.array([t["p_raw"] for t in tests])
    q, passed = M.benjamini_hochberg(p, M.FDR_Q)
    for t, qq, pa in zip(tests, q, passed):
        t["bh_q"] = float(qq)
        t["passes_fdr"] = bool(pa)
    n_tests = len(tests)
    if verbose:
        print(f"multiplicity: {n_tests} tests, BH-FDR at q={M.FDR_Q} -> "
              f"{int(passed.sum())} survive")
        floor = 1.0 / (n_surr + 1.0)
        k_min = int(math.ceil(floor * max(n_tests, 1) / M.FDR_Q))
        n_at_floor = sum(1 for t in tests if t["p_raw"] <= floor + 1e-12)
        if k_min > 1:
            print(f"  WARNING: {n_surr} surrogates censor every p at the floor "
                  f"{floor:.5f}. BH can reject only if >= {k_min} tests tie at that "
                  f"floor ({n_at_floor} did). Survivors are provisional; rerun with "
                  f"more surrogates before trusting the ordering.")

    # ------------------------------------------------- aliasing audit ------
    fmap = {f.name: f for f in feats}
    for t in tests:
        if not t["passes_fdr"]:
            continue
        key = ("audit:" + t["test"] + ":" + t["feature"] + ":"
               + str(t.get("lag", t.get("mark", ""))))
        if ckpt.done(key):
            t["aliasing"] = ckpt.state["results"][key]
            continue
        if t["test"] == "glm_poisson_offset_etas":
            aud = M.aliasing_audit_glm(fmap[t["feature"]], t["lag"], window,
                                       counts, offset, n_surr, rng)
        elif t["test"] == "lomb_scargle_peak":
            aud = M.aliasing_audit_period(t["ladder"]["winning_period_days"], counts,
                                          offset, marks["day_float"], n_surr, rng,
                                          float(window.start))
        else:
            aud = {"verdict": "N/A (mark test: no time lattice claim)"}
        ckpt.put(key, aud)
        t["aliasing"] = aud
        if verbose:
            print(f"  aliasing audit {t['feature']} ({t['test']}): {aud['verdict']}")

    # ------------------------------------------------------- outputs -------
    ckpt.state["tests"] = tests
    ckpt.state["n_tests"] = n_tests
    ckpt.state["elapsed_seconds"] = round(time.time() - t_start, 1)
    ckpt.state["complete"] = True
    ckpt.save()

    rep = write_report(session_dir, cfg, ckpt.state, tests, feats, counts, offset,
                       marks, base, window)
    stubs = write_stubs(session_dir, cfg, tests)

    # One ledger line per SESSION, not per invocation: a resumed session is the same
    # sweep continued, so counting it twice would inflate the reported multiplicity.
    if not ckpt.state.get("ledger_logged"):
        splits._append_jsonl(splits.EXPLORE_COUNT, {
            "ts": _dt.datetime.now(_dt.timezone.utc).isoformat(),
            "kind": "mine",
            "n_tests": n_tests,
            "hash": splits.config_hash(cfg),
            "config": cfg,
            "session_dir": session_dir.replace("\\", "/"),
        })
        ckpt.state["ledger_logged"] = True
        ckpt.save()
        if verbose:
            print(f"session ledger       -> {splits.EXPLORE_COUNT} "
                  f"(kind=mine, n_tests={n_tests})")
    elif verbose:
        print(f"session ledger       -> already logged for this session "
              f"(resumed run; multiplicity is not double-counted)")
    return {"session_dir": session_dir, "report": rep, "stubs": stubs,
            "n_tests": n_tests, "n_pass": int(sum(t["passes_fdr"] for t in tests)),
            "elapsed": ckpt.state["elapsed_seconds"]}


# ------------------------------------------------------------------ report ---
def _row_effect(t):
    if t["test"] == "glm_poisson_offset_etas":
        if t["kind"] == "phase":
            return f"amp {t['amplitude_log_rate']:.4f} log-rate"
        return f"beta {t['beta'][0]:+.4f} +/- {t['se'][0]:.4f} /sd"
    if t["test"] == "lomb_scargle_peak":
        return f"LS power {t['power']:.4f}"
    return f"{t['test']} {t['effect']:+.4f}"


def _row_where(t):
    if t["test"] == "glm_poisson_offset_etas":
        return f"lag {t['lag']} d"
    if t["test"] == "lomb_scargle_peak":
        return f"P = {t['period_days']:.4g} d"
    return f"mark {t['mark']}"


def _amp_note(t):
    if t["test"] == "lomb_scargle_peak":
        return f"+/-{t['pct_rate_modulation']:.1f}% folded rate"
    if t["test"] == "glm_poisson_offset_etas":
        return f"+/-{t['pct_rate_modulation']:.2f}% rate"
    return "rank correlation (no rate amplitude)"


def write_report(session_dir, cfg, state, tests, feats, counts, offset, marks,
                 base, window):
    path = os.path.join(session_dir, "report.md")
    order = sorted(tests, key=lambda t: (t["bh_q"], t["p_raw"]))
    n_pass = sum(t["passes_fdr"] for t in tests)
    L = []
    A = L.append

    A(f"# mine session report -- {os.path.basename(session_dir)}")
    A("")
    A("> **" + M.GENERATOR_NOT_EVIDENCE.split(".")[0] + ".**")
    A(">")
    for line in _wrap(M.GENERATOR_NOT_EVIDENCE):
        A("> " + line)
    A("")
    A("> **Standing warning (EQ-24, verbatim):** " + M.MIGNAN_BROCCARDO)
    A("")
    A("## Configuration")
    A("")
    A("```json")
    A(json.dumps(cfg, indent=1))
    A("```")
    A("")
    A(f"- engine v{state['engine_version']}, preset `{cfg['preset']}`, "
      f"elapsed {state.get('elapsed_seconds')} s")
    A(f"- baseline: **{base.name}** -- {base.caveat}")
    A(f"- mining window: exploration days [{window.start}, {window.stop}) = "
      f"{counts.size} days (365 d ETAS burn-in dropped)")
    A(f"- target: daily domain-wide counts vs sum of lambda_etas; "
      f"{counts.sum():.0f} observed events, {offset.sum():.1f} expected")
    A(f"- marks: {marks['mag'].size} events (magnitude, depth)")
    A(f"- features: {len(feats)} "
      + ", ".join(f"family {k}: {v}" for k, v in sorted(
          _count_by(feats, lambda f: f.family).items())))
    A(f"- **{state['n_tests']} tests** in this session; BH-FDR at q = {cfg['fdr_q']} "
      f"-> **{n_pass} survive**")
    floor = 1.0 / (cfg["n_surrogates"] + 1.0)
    m = max(state["n_tests"], 1)
    k_min = int(math.ceil(floor * m / cfg["fdr_q"]))
    n_at_floor = sum(1 for t in tests if t["p_raw"] <= floor + 1e-12)
    A(f"- **surrogate resolution (read this before the table):** with "
      f"{cfg['n_surrogates']} surrogates the smallest attainable empirical p is "
      f"1/(N+1) = {floor:.5f}. Every p at that value is CENSORED -- the true p is "
      f"somewhere at or below it -- so the BH q attached to it is an upper bound "
      f"computed from a tie, not a measurement. "
      + (f"**No test in this session can pass FDR at all** ({k_min} > "
         f"{m} tests would have to tie at the floor)."
         if k_min > m else
         f"At this resolution BH can only reject if at least {k_min} "
         f"{'test ties' if k_min == 1 else 'tests tie'} at the floor "
         f"simultaneously; {n_at_floor} did."
         + ("" if k_min == 1 else
            " A survivor list produced from a tie at the floor is provisional --"
            " rerun with more surrogates before believing the ordering.")))
    A(f"- the period scan carries its own, coarser floor: its max-power Monte Carlo "
      f"is capped at {state['results']['period_scan']['scan']['n_mc']} draws "
      f"(each draw is a full periodogram), so no period peak can report p below "
      f"{1.0 / (state['results']['period_scan']['scan']['n_mc'] + 1):.5f}.")
    A(f"- ledger: one line appended to `engine/EXPLORE_COUNT.jsonl` "
      f"(kind=mine, n_tests={state['n_tests']}), so holdout multiplicity reporting "
      f"includes this sweep. No holdout hash was spent.")
    A("")
    dl = state.get("data_log", [])
    if dl:
        A("### Optional downloads (family 3)")
        A("")
        A("| source | status | detail |")
        A("| --- | --- | --- |")
        for r in dl:
            det = (f"{r.get('n_bytes','?')} B, sha256 `{r.get('sha256','')[:16]}`, "
                   f"coverage {r.get('coverage')}" if r.get("status") == "ok"
                   else str(r.get("reason", "")))
            A(f"| `{r['key']}` | {r['status']} | {det} |")
        A("")
        A("Frozen copies and hashes are in `engine/out/mine/data_log.jsonl` "
          "(sniff-grade hygiene). This is deliberately NOT `download_log.md`, which "
          "belongs to frozen-protocol runs only.")
        A("")

    A("## Ranked candidates (top 25 by BH q, then raw p)")
    A("")
    A("| # | feature | lag/period | effect | raw p | BH q | surrogates | ladder | "
      "aliasing | amplitude honesty |")
    A("| ---: | --- | --- | --- | ---: | ---: | ---: | --- | --- | --- |")
    for i, t in enumerate(order[:25], 1):
        lad = (t["ladder"]["verdict"] if t["test"] == "lomb_scargle_peak" else "n/a")
        ali = t.get("aliasing", {}).get("verdict", "-- (did not pass FDR)")
        A(f"| {i} | `{t['feature']}` | {_row_where(t)} | {_row_effect(t)} | "
          f"{t['p_raw']:.4g} | {t['bh_q']:.4g} | {t.get('n_surrogates', '')} | "
          f"{lad} | {ali} | {_amp_note(t)} |")
    A("")

    A("## Survivors (post-FDR) in full")
    A("")
    if n_pass == 0:
        A("**Nothing survived BH-FDR at q = %.2f.**" % cfg["fdr_q"])
        A("")
        A("That is a clean result for the instrument and is reported as such: with "
          "the ETAS baseline absorbing Omori-Utsu clustering, the global daily "
          "residual carries no ephemeris-, space-weather- or catalogue-mark "
          "association that survives this session's multiplicity. It is NOT evidence "
          "of absence -- see the power note below.")
    for t in order:
        if not t["passes_fdr"]:
            continue
        A(f"### `{t['feature']}` -- {_row_where(t)}")
        A("")
        A(f"- test: `{t['test']}`, {t.get('df', '?')} df")
        A(f"- effect: {_row_effect(t)}; {_amp_note(t)}")
        A(f"- raw p = {t['p_raw']:.4g} (empirical, {t.get('n_surrogates')} "
          f"circular-shift surrogates); BH q = {t['bh_q']:.4g}")
        if t["test"] == "glm_poisson_offset_etas":
            A(f"- in-sample bits/event vs etas-v1: {t['bits_per_event']:+.6f}")
        if t["test"] == "lomb_scargle_peak":
            lad = t["ladder"]
            A(f"- harmonic ladder ({lad['n_bins']} phase bins, epoch-folding "
              f"likelihood ratio): **{lad['verdict']}**")
            A("")
            A("  | rung | period (d) | fold LR |")
            A("  | --- | ---: | ---: |")
            for r in lad["rungs"]:
                A(f"  | {r['multiplier']:.4g} x P | {r['period_days']:.4g} | "
                  f"{r['fold_lr']:.2f} |")
            A("")
        ali = t.get("aliasing")
        if ali:
            A(f"- aliasing audit: **{ali.get('verdict')}** -- {ali.get('rule','')}")
            if "z_ratio_vs_1d" in ali:
                A(f"  - surrogate-calibrated effect ratio vs 1 d binning: "
                  f"{ali['z_ratio_vs_1d']}")
            if "fold_lr_ratio_vs_1d" in ali:
                A(f"  - folded-concentration ratio vs 1 d binning: "
                  f"{ali['fold_lr_ratio_vs_1d']}")
                u = ali["unbinned_rayleigh"]
                A(f"  - unbinned event-time Rayleigh (n = {u['n_events']} events, "
                  f"ETAS-simulated null): R = {u['R']:.3f}, p = {u['p']:.4g}")
        A("")
        A("- " + M.AMPLITUDE_HONESTY)
        A("")

    A("## Instrument validation -- what the miner rediscovered")
    A("")
    A("A miner run against a clean residual should find NOTHING at strong "
      "significance. Anything it does find that is already known physics or a known "
      "catalogue artifact is **validation of the instrument, not discovery**. The "
      "candidates below are pre-labelled as such wherever they match a known "
      "structure:")
    A("")
    A("| candidate | raw p | known-structure label |")
    A("| --- | ---: | --- |")
    for t in order[:25]:
        A(f"| `{t['feature']}` {_row_where(t)} | {t['p_raw']:.4g} | "
          f"{_known_label(t)} |")
    A("")

    A("## Method notes, including one honest deviation")
    A("")
    A("1. **Surrogate nulls.** Every GLM and mark test is calibrated against "
      "CIRCULAR TIME SHIFTS of the target (the counts/offset pair for GLM tests, "
      "the mark series along the time-ordered event sequence for mark tests). "
      "Because the score statistic and the rank correlations are exact "
      "cross-correlations, ALL admissible shifts are evaluated in closed form "
      "rather than sampled; when the requested surrogate count exceeds the number "
      "of admissible shifts the run uses every one of them and reports the actual "
      "count in the `surrogates` column. Shifts within 30 days of zero are excluded.")
    A(f"2. **DEVIATION, flagged not hidden.** A circular shift of an evenly sampled "
      f"series does not change its Lomb-Scargle power at all -- a shift is a pure "
      f"phase rotation -- so circular-shift surrogates are mathematically vacuous "
      f"as a null for the PERIOD SCAN. The period scan therefore uses two Monte "
      f"Carlo nulls instead: an AR(1) red-noise null matched to the residual's own "
      f"lag-1 autocorrelation (phi = "
      f"{state['results']['period_scan']['scan']['ar1_phi']:.3f}) and a permutation "
      f"(white) null, and reports the MORE CONSERVATIVE of the two p-values. "
      f"Circular shifts remain the null everywhere else.")
    A("2b. **Block length is measured, not guessed.** A block bootstrap is a valid "
      "null only if its blocks are long compared to the structure being tested; too "
      "short and the surrogates cannot represent the feature's own timescale, which "
      "makes the null too narrow and the test ANTI-CONSERVATIVE. Cyclic features use "
      "2x their known period; everything else uses 4x the e-folding time of its own "
      "autocorrelation, clipped to [30, 800] days. The per-test value is recorded as "
      "`block_days`. Residual risk, stated: a feature whose autocorrelation is fast "
      "but whose ENVELOPE is slow (Ap is the example -- storm-scale ACF, solar-cycle "
      "envelope) gets a short block and its p is therefore optimistic.")
    A("3. **Harmonic ladder.** Every candidate period P is scored by epoch-folding "
      "likelihood ratio at {P/3, P/2, P, 2P, 3P}, extended at the winning edge "
      "while the score keeps improving. The reported period is the winning rung and "
      "all rung scores are printed. Rungs longer than a third of the record "
      f"({counts.size / 3.0:.0f} d here) are not scored at all: with fewer than "
      "three observed cycles an epoch fold measures the record length, not a "
      "period.")
    A("3b. **The two period-scan nulls disagree on purpose, and the report takes "
      "the loser.** The global daily residual is red (AR(1) phi = "
      f"{state['results']['period_scan']['scan']['ar1_phi']:.3f}), so a permutation "
      "(white) null calls essentially every peak significant. The AR(1) null is the "
      "one that knows the residual is autocorrelated. Reporting the more "
      "conservative of the two is what keeps the period scan from manufacturing a "
      "candidate list out of red noise -- and the gap between the two columns is "
      "the size of the mistake that would have been made.")
    A("4. **Aliasing audit.** Every post-FDR survivor is re-tested at 2-day and "
      "7-day binning, and period claims are additionally re-tested on UNBINNED "
      "event times (real catalogue timestamps, sub-day precision) against an "
      "ETAS-simulated inhomogeneous-Poisson null. A pattern whose effect halves "
      "under re-binning, or fails unbinned, is flagged LATTICE-SUSPECT: a pattern "
      "that moves when the lattice moves is the lattice.")
    A("5. **Causality.** Family-1 and family-2 features are deterministic functions "
      "of t (ephemeris) and are therefore EXEMPT from the causality-shuffle test by "
      "construction: no rearrangement of the catalogue can change them. Family-3 "
      "(downloaded indices) are lagged one day so day t uses only values published "
      "strictly before t. Family-4 (catalogue-derived) are strictly trailing "
      "windows, exclusive of today, and ARE included in "
      "`engine/tests/test_causality.py`.")
    A("6. **In-sample.** Every effect size here is fitted and evaluated on the "
      "exploration window. They are upper bounds, exactly as `--mode explore` is.")
    A("")
    A("## Power / bounds note")
    A("")
    sd = float((counts - offset).std())
    mean = float(counts.mean())
    A(f"The global daily residual has mean count {mean:.2f} events/day and residual "
      f"sd {sd:.2f} over {counts.size} days. A sinusoidal rate modulation of "
      f"amplitude A is detectable at this session's thresholds only above roughly "
      f"{100 * 2.5 * sd / (mean * math.sqrt(counts.size / 2.0)):.2f}% "
      f"(a 2.5-sigma-equivalent rule of thumb, not a formal power calculation). "
      f"Absence of a survivor is a BOUND at about that amplitude, not an absence of "
      f"effect.")
    A("")
    A("That rule of thumb assumes independent days, so it is OPTIMISTIC -- and most "
      "optimistic exactly where the block bootstrap is longest. A feature tested "
      "with 700-day blocks has an effective sample of ~11 independent pieces, not "
      "7716, and its true detection threshold is several times the number above. "
      "Read the per-test `block_days` before quoting any bound: the longer the "
      "block, the weaker the bound.")
    A("")
    A("---")
    A("")
    A(M.AMPLITUDE_HONESTY)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return path


KNOWN = [
    (("annual_phase", "sun_declination"), "KNOWN: annual/seasonal cycle -- seasonal "
     "catalogue completeness, hydrological loading and reporting cadence all live "
     "here. Validation of the instrument, not discovery."),
    (("b_value_90d", "deep_fraction_30d", "mean_depth_30d"),
     "KNOWN: catalogue-composition drift (network eras, Mc drift, deep-slab "
     "sequences). Expected to be found; validation, not discovery."),
    (("moon_synodic_phase", "moon_anomalistic_phase", "moon_draconic_phase",
      "spring_neap_phase", "half_draconic_phase", "perigean_spring_beat",
      "eclipse_year_beat", "annual_synodic_beat", "perigee_syzygy",
      "tidal_potential_proxy", "moon_distance", "sun_moon_elongation",
      "moon_declination", "moon_abs_declination", "declination_product"),
     "KNOWN-CLASS: solid-earth tidal forcing. The tidal corpse says expect a bound "
     "at the sub-percent level, not a detection."),
    (("F107_solar_flux", "Ap_geomagnetic", "length_of_day"),
     "OPEN-CLASS: space weather / rotation. Contested literature; no accepted "
     "mechanism at this amplitude."),
]


def _known_label(t):
    name = t["feature"]
    for names, label in KNOWN:
        if name in names:
            return label
    if t["test"] == "lomb_scargle_peak":
        p = t["period_days"]
        for ref, lab in ((365.25, "annual"), (182.6, "semi-annual"),
                         (29.53, "synodic month"), (14.77, "spring-neap"),
                         (27.55, "anomalistic month"), (13.66, "declination tide"),
                         (7.0, "WEEKLY -- catalogue/reporting cadence, an artifact"),
                         (1.0, "diurnal")):
            if abs(math.log(p / ref)) < 0.03:
                return f"KNOWN period: {lab}. Validation of the instrument."
        return "UNLABELLED period -- no known structure matched; treat as candidate."
    return "unlabelled"


def _count_by(items, keyfn):
    out = {}
    for i in items:
        out[keyfn(i)] = out.get(keyfn(i), 0) + 1
    return out


def _wrap(text, width=88):
    words, line, out = text.split(), "", []
    for w in words:
        if len(line) + len(w) + 1 > width:
            out.append(line)
            line = w
        else:
            line = (line + " " + w).strip()
    if line:
        out.append(line)
    return out


# ------------------------------------------------------------------- stubs ---
def write_stubs(session_dir, cfg, tests):
    path = os.path.join(session_dir, "stubs.json")
    entries = []
    for t in sorted(tests, key=lambda t: (t["bh_q"], t["p_raw"])):
        if not t["passes_fdr"]:
            continue
        entries.append({
            "status": "DRAFT STUB -- generator output, not a K-entry",
            "observable": _observable(t),
            "feature": t["feature"],
            "family": t["family"],
            "test": t["test"],
            "where": _row_where(t),
            "effect_size": _row_effect(t),
            "rate_modulation": _amp_note(t),
            "p_raw": t["p_raw"],
            "bh_q": t["bh_q"],
            "n_surrogates": t.get("n_surrogates"),
            "ladder": t.get("ladder", {}).get("verdict"),
            "aliasing_verdict": t.get("aliasing", {}).get("verdict"),
            "known_structure_label": _known_label(t),
            "caveats": [M.GENERATOR_NOT_EVIDENCE, M.AMPLITUDE_HONESTY],
            "next_step": ("Write up as a K-entry with a Popper ruling, then -- and "
                          "only then -- freeze a config and spend ONE holdout hash "
                          "via `python -u -m engine.cli run --mode holdout`."),
        })
    payload = {
        "banner": M.GENERATOR_NOT_EVIDENCE,
        "standing_warning_eq24": M.MIGNAN_BROCCARDO,
        "session": os.path.basename(session_dir),
        "config": cfg,
        "n_tests": len(tests),
        "n_stubs": len(entries),
        "stubs": entries,
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    return path


def _observable(t):
    if t["test"] == "lomb_scargle_peak":
        return (f"Global daily M>={4.5} occurrence residual against etas-v1 carries a "
                f"periodicity at {t['ladder']['winning_period_days']:.4g} d.")
    if t["test"] == "glm_poisson_offset_etas":
        return (f"Global daily M>={4.5} occurrence rate depends on `{t['feature']}` "
                f"at lag {t['lag']} d, after ETAS residualization.")
    return (f"Event {t['mark']} is associated with `{t['feature']}` at the time of "
            f"the event, after ETAS residualization of occurrence.")
