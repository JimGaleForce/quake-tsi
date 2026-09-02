"""PHASE 2 STEP 1 -- THE LEARNED ARM EXTENDED TO NON-TIDAL FORCING. Priced 1 (block 3).

A NEW FILE ON PURPOSE. `exp_learned.py` and `results_learned.json` are the committed P-1.4
result and are not touched, so that number stays reproducible. This module copies its
design exactly -- same catalogue, same exploration split, same matched case-control
strata, same K, same seeds, same hyperparameters, same within-stratum permutation null --
and adds one thing: a THIRD covariate block of non-tidal forcing, built by
`engine_ext_forcing.py`.

THE NESTED LADDER. Model A is the nuisance/observer block and is identical to the
committed arm's model A, which is what makes the reproduction check meaningful.

    A   nuisance      local solar hour, day of week, day of year
    B   A + tidal     EXACTLY the committed model B; reproduces the committed P-1.4 dAUC
    C   A + block3    LOD / polar motion / Kp / ap / OMNI at event time
    D   B + block3    the full stack

    DECLARED: dAUC(C vs A) and dAUC(D vs B), each against the same permutation null.

WHY BLOCK 3 IS TWO DIFFERENT THINGS AND IS REPORTED AS TWO DIFFERENT THINGS. LOD and polar
motion are a genuine, small mechanical channel with a specific published claim attached.
Kp / ap / OMNI are a MATCHED PLACEBO with a near-zero physical prior. A fire on the placebo
half is a false positive by construction, which is exactly why the block is worth running:
playbook P-2.4 asks for the rate at which this pipeline manufactures survivors on a
property class that has real temporal structure and no physics. That rate is measured here
(`placebo`), over the same features, under the same permutation null, at nominal alpha.

===========================================================================
THE ARTIFACTS THAT COULD FAKE THIS RESULT -- NAMED BEFORE THE RUN (rule 7)
===========================================================================

(1) HUMAN-SCHEDULE LEAKAGE IN Kp. Kp is a UT-gridded index built from magnetometer
    observatories on human operating schedules. If it carries a diurnal or weekly cycle of
    its own it can proxy the catalogue's cultural-noise day/night artifact, and SoCal sits
    at a fixed UT offset so a UT cycle aliases straight onto local solar hour. Model A is
    handed solar hour, day of week and day of year for free, which is the primary defence.
    `engine_ext_forcing.schedule_audit` MEASURES the residual cycle rather than asserting
    it away, and the measurement is in the output.

(2) COUNT-PATH PERIODICITY / QUANTISATION. Kp is quantised to thirds and ap to a discrete
    ladder; their within-stratum ties are plateaus, not structure. The audit counts
    distinct values. The permutation null is immune to this by construction -- it permutes
    labels, not features -- which is the reason the null is done this way and not
    analytically.

(3) LONG-PERIOD LOD TREND CONFOUNDED WITH CATALOGUE COMPLETENESS DRIFT. LOD has decadal
    and interannual structure; SoCal network completeness also drifts across 2008-2015.
    RESIDUAL EXPOSURE, stated plainly: the matching removes most of this and NOT all of
    it. Each stratum's four phantom times are drawn from the event's OWN +/- fortnight, so
    epoch is matched to within a fortnight and the decadal channel is differenced away.
    That is also why this arm CANNOT test the Bendick & Bilham decadal claim at all -- it
    only sees sub-fortnight structure. What survives the matching is any LOD variation
    steep at a two-week scale (the seasonal and fortnightly-tidal terms in LOD itself),
    and those are correlated with the tidal block by construction, which is precisely why
    D vs B, not D vs A, is the declared comparison for the full stack.

(4) THE DAY/NIGHT DETECTION ARTIFACT. Roughly 2 percent of catalogued events are displaced
    into local night by cultural noise (`exp_diurnal_discriminator.py`). It is handed to
    model A for free, so C and D can only win beyond it. Model A's absolute AUC of 0.5038
    is how big that artifact looks inside a matched stratum, and it is the scale against
    which every dAUC here should be read.

(5) MISSINGNESS AS A COVARIATE. OMNI2 has real data gaps. They are carried as NaN plus an
    explicit `miss_*` indicator column, never zero-filled. A gap is time-correlated, so in
    principle the indicator could discriminate a case from its own phantoms. The fill
    counts are reported per feature, and the permutation null sees the identical
    indicators, so anything the indicators manufacture appears in the null too.

2026-09-02: this arm runs on the corrected loaders (days since 1970-01-01Z; see
CORRECTIONS.md 2026-09-02 and tests/test_epoch_consistency.py). Model B reproduces the
corrected exp_learned.py P-1.4 result.

NOTHING IS EDITED IN `exp_learned.py`, `exp_highn.py` or `exp_mass_screen.py`. The
reproduction check below reads the committed dAUC(B vs A) from `results_learned.json` at
runtime rather than a hard-coded number, so it stays correct across any future re-run of
the committed arm.

STATE CLASS: first-run, exploration split. The holdout is not touched.
"""

from __future__ import annotations

import datetime as _dt
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import engine_ext_forcing as F
import exp_highn as HN
import exp_learned as L
import exp_mass_screen as MS

OUT_JSON = HERE / "results_learned_ext.json"
RESULTS_LEARNED_JSON = HERE / "results_learned.json"
ZEN = HERE / "data" / "xue_lu_zenodo"

# every knob inherited from the committed arm, by reference not by copy
CONTROLS_PER_CASE = L.CONTROLS_PER_CASE
NULL_HALF_WINDOW_DAYS = L.NULL_HALF_WINDOW_DAYS
TRAIN_FRAC = L.TRAIN_FRAC
RNG_SEED = L.RNG_SEED
NUISANCE = L.NUISANCE
TIDAL = L.TIDAL

N_PERMUTATIONS = 200            # the first 40 use the committed arm's own seeds
COMMITTED_N_PERM = L.N_PERMUTATIONS
BLOCK3 = list(F.BLOCK3)

# committed P-1.4 numbers, READ AT RUNTIME from results_learned.json -- not hard-coded, so
# this stays correct across any future re-run of the committed arm (e.g. the 2026-09-02
# epoch fix, which changed these numbers from their pre-fix values dAUC +0.001552,
# p = 0.3415).
_committed = json.loads(RESULTS_LEARNED_JSON.read_text(encoding="utf-8"))

# counted invariants, from INVENTORY.md as MEASURED there
EXPECT = {
    "qtm_rows": 45069,          # INVENTORY 1a
    "n_cases_explore": 33293,   # results_learned.json
    "n_holdout": 11776,         # results_learned.json
    "committed_auc_a": _committed["model_A"]["auc"],
    "committed_auc_b": _committed["model_B"]["auc"],
    "committed_dauc": _committed["declared_statistic_delta_auc"],
    "committed_p_dauc": _committed["p_delta_auc"],
    "eop_rows": 23381, "kp_rows": 276536, "omni_rows": 280512, "omni_files": 32,
}

# P-2.4 placebo split of block 3: which half has a physical prior and which does not
BLOCK3_MECHANICAL = list(F.EOP_FEATURES)
BLOCK3_PLACEBO = [c for c in F.BLOCK3_VALUES if c not in BLOCK3_MECHANICAL]

ALPHA = 0.05
Z_CRIT = 1.959963985                       # two-sided nominal 0.05
N_PLACEBO_DRAWS = 4000
PLACEBO_CHUNK = 250


def stratum_z(values, member):
    """Within-stratum standardisation across the K+1 candidates. [n, K+1] -> [n, K+1].

    Strata holding a non-finite candidate, or with zero spread (Kp plateaus), are marked
    unusable and returned as all-zero with a False in the second return value.
    """
    v = values[member]
    ok = np.all(np.isfinite(v), axis=1)
    mu = np.where(ok, v.mean(axis=1), 0.0)[:, None]
    sd = np.where(ok, v.std(axis=1), 0.0)[:, None]
    usable = ok & (sd[:, 0] > 0)
    z = np.where(usable[:, None], (v - mu) / np.where(sd > 0, sd, 1.0), 0.0)
    return z, usable


def placebo_Z(cols, names, member, seed=770000):
    """P-2.4 machinery: observed and null matched-set Z for every named feature.

    For each feature the statistic is the standard conditional one for a matched set: the
    sum over strata of the CASE's within-stratum z, divided by sqrt(n usable strata).
    Under H0 the case is uniform over the K+1 candidates, so E[z_case] = 0 and
    Var[z_case] = 1 by construction and Z is unit normal.

    The SAME random re-pick is used across features within a draw, so the cross-feature
    correlation of block 3 is preserved and "survivors per draw" is a family-wise count
    rather than a product of independent Bernoullis.

    Returns (obs[n_features], null[D, n_features], n_usable[n_features]).
    """
    n, Kp1 = member.shape
    Zs, Us = [], []
    for nm in names:
        z, u = stratum_z(cols[nm], member)
        Zs.append(z)
        Us.append(u)
    obs = np.array([float(Zs[j][Us[j], 0].sum() / np.sqrt(max(Us[j].sum(), 1)))
                    for j in range(len(names))])
    denom = np.array([np.sqrt(max(int(u.sum()), 1)) for u in Us])
    null = np.empty((N_PLACEBO_DRAWS, len(names)), dtype=np.float64)
    rng = np.random.default_rng(seed)
    for a in range(0, N_PLACEBO_DRAWS, PLACEBO_CHUNK):
        b = min(a + PLACEBO_CHUNK, N_PLACEBO_DRAWS)
        pick = rng.integers(0, Kp1, size=(n, b - a))
        for j, u in enumerate(Us):
            sel = np.take_along_axis(Zs[j], pick, axis=1)      # [n, chunk]
            null[a:b, j] = sel[u].sum(axis=0) / denom[j]
    return obs, null, np.array([int(u.sum()) for u in Us])


def survivor_report(names, idx, obs, null, nusable):
    """Survivor counts at nominal alpha for the subset of features in `idx`."""
    o = obs[idx]
    nl = null[:, idx]
    hit = np.abs(nl) > Z_CRIT
    per_draw = hit.sum(axis=1).astype(np.float64)
    return {
        "features_tested": [names[i] for i in idx],
        "n_features": int(len(idx)),
        "alpha_nominal_two_sided": ALPHA,
        "n_null_draws": int(N_PLACEBO_DRAWS),
        "n_usable_strata": {names[i]: int(nusable[i]) for i in idx},
        "observed_Z": {names[i]: float(obs[i]) for i in idx},
        "observed_n_survivors": int((np.abs(o) > Z_CRIT).sum()),
        "observed_survivor_names": [names[i] for i in idx if abs(obs[i]) > Z_CRIT],
        "null_survivor_rate_per_feature": {names[i]: float(hit[:, j].mean())
                                           for j, i in enumerate(idx)},
        "null_mean_survivor_rate_over_features": float(hit.mean()),
        "null_mean_survivors_per_draw": float(per_draw.mean()),
        "null_expected_if_independent": float(ALPHA * len(idx)),
        "null_familywise_any_survivor_rate": float((per_draw > 0).mean()),
        "null_p95_survivors_per_draw": float(np.quantile(per_draw, 0.95)),
        "null_max_survivors_per_draw": int(per_draw.max()),
    }


def main():
    t_start = time.time()
    rng = np.random.default_rng(RNG_SEED)

    # ---------------------------------------------------------- catalogue, counted
    t_raw, la_raw, lo_raw, dp_raw, mg_raw = HN.load_zenodo(ZEN / "QTM_decluster_m0.1.txt")
    MS.assert_epoch(t_raw, 2008, "QTM_declustered")
    inv = {"qtm_rows_parsed": int(t_raw.size),
           "qtm_rows_expected_INVENTORY": EXPECT["qtm_rows"],
           "qtm_rows_match": bool(t_raw.size == EXPECT["qtm_rows"])}
    t, la, lo, dp, mg = HN.split(t_raw, la_raw, lo_raw, dp_raw, mg_raw)[:5]
    n_hold = HN.split(t_raw, la_raw, lo_raw, dp_raw, mg_raw)[-1]
    n = t.size
    inv.update({"n_cases_explore": int(n),
                "n_cases_expected_committed": EXPECT["n_cases_explore"],
                "n_cases_match": bool(n == EXPECT["n_cases_explore"]),
                "n_holdout": int(n_hold),
                "n_holdout_match": bool(n_hold == EXPECT["n_holdout"])})
    print("QTM declustered: %d rows parsed (INVENTORY %d), exploration %d, holdout %d"
          % (t_raw.size, EXPECT["qtm_rows"], n, n_hold), flush=True)

    # ------------------------------------------------------------ the matched design
    K = CONTROLS_PER_CASE
    off = rng.uniform(-NULL_HALF_WINDOW_DAYS, NULL_HALF_WINDOW_DAYS, size=(n, K))
    t_all = np.concatenate([t, (t[:, None] + off).ravel()])
    la_all = np.concatenate([la, np.repeat(la, K)])
    lo_all = np.concatenate([lo, np.repeat(lo, K)])
    y = np.concatenate([np.ones(n), np.zeros(n * K)])
    stratum = np.concatenate([np.arange(n), np.repeat(np.arange(n), K)])
    member = np.empty((n, K + 1), dtype=np.int64)
    member[:, 0] = np.arange(n)
    member[:, 1:] = n + np.arange(n)[:, None] * K + np.arange(K)[None, :]

    # EPOCH: HN.load_zenodo returns days since 1970-01-01Z (fixed at the source,
    # 2026-09-02; see CORRECTIONS.md and tests/test_epoch_consistency.py), so the
    # explicit epoch here is the Unix epoch.
    jd_true = F.jd_from_days_since(t_all, _dt.datetime(1970, 1, 1,
                                                       tzinfo=_dt.timezone.utc))

    # ------------------------------------------------------------------ the features
    print("building features for %d rows ..." % t_all.size, flush=True)
    cols = L.build_features(t_all, la_all, lo_all)            # committed, unmodified

    series, load_audit = F.load_forcing()
    for k, v in (("iers_eop", "eop_rows"), ("kp_ap", "kp_rows"), ("omni2", "omni_rows")):
        got = load_audit["sources"][k]["n_rows"]
        load_audit["sources"][k]["n_rows_expected"] = EXPECT[v]
        load_audit["sources"][k]["n_rows_match"] = bool(got == EXPECT[v])
        print("  %-9s %7d rows parsed (expected %7d)  %s"
              % (k, got, EXPECT[v], "OK" if got == EXPECT[v] else "MISMATCH"), flush=True)
    b3, fill = F.event_block(jd_true, series)
    print("  block3 fill: %d/%d rows have at least one missing source"
          % (fill["n_rows_any_source_missing"], fill["n_rows"]), flush=True)

    sched = F.schedule_audit(series, float(jd_true.min()), float(jd_true.max()))

    allcols = dict(cols)
    allcols.update(b3)

    FEAT = {
        "A": NUISANCE,
        "B": NUISANCE + TIDAL,
        "C": NUISANCE + BLOCK3,
        "D": NUISANCE + TIDAL + BLOCK3,
    }
    X = {k: np.column_stack([allcols[c] for c in v]) for k, v in FEAT.items()}

    # ------------------------------------------------------- the split (unchanged)
    order = np.argsort(t)
    cut = order[int(TRAIN_FRAC * n)]
    train_strata = np.zeros(n, dtype=bool)
    train_strata[t < t[cut]] = True
    tr = train_strata[stratum]
    te = ~tr
    inv.update({"train_rows": int(tr.sum()), "test_rows": int(te.sum()),
                "test_cases": int(y[te].sum())})
    print("  train rows %d   test rows %d   (test cases %d)"
          % (tr.sum(), te.sum(), int(y[te].sum())), flush=True)

    # ------------------------------------------------------------------- observed
    obs = {}
    for k in ("A", "B", "C", "D"):
        a, l = L.fit_score(X[k][tr], y[tr], X[k][te], y[te], 0)
        obs[k] = {"auc": a, "log_loss": l, "n_features": len(FEAT[k])}
        print("  MODEL %-2s  AUC %.6f  logloss %.6f  (%d features)"
              % (k, a, l, len(FEAT[k])), flush=True)

    d_obs = {
        "B_vs_A": obs["B"]["auc"] - obs["A"]["auc"],
        "C_vs_A": obs["C"]["auc"] - obs["A"]["auc"],
        "D_vs_B": obs["D"]["auc"] - obs["B"]["auc"],
        "D_vs_A": obs["D"]["auc"] - obs["A"]["auc"],
    }
    repro = {
        "committed_auc_A": EXPECT["committed_auc_a"], "reproduced_auc_A": obs["A"]["auc"],
        "committed_auc_B": EXPECT["committed_auc_b"], "reproduced_auc_B": obs["B"]["auc"],
        "committed_dAUC_B_vs_A": EXPECT["committed_dauc"],
        "reproduced_dAUC_B_vs_A": d_obs["B_vs_A"],
        "abs_difference": abs(d_obs["B_vs_A"] - EXPECT["committed_dauc"]),
        "exact_match": bool(abs(d_obs["B_vs_A"] - EXPECT["committed_dauc"]) < 1e-12),
    }
    print("  REPRODUCTION: committed dAUC(B-A) %+.6f, reproduced %+.6f  -> %s"
          % (EXPECT["committed_dauc"], d_obs["B_vs_A"],
             "EXACT" if repro["exact_match"] else "DIFFERS"), flush=True)

    # ---------------------------------------------------------------- the null
    print("\npermutation null: %d refits of ALL FOUR models "
          "(first %d seeds are the committed arm's own) ..."
          % (N_PERMUTATIONS, COMMITTED_N_PERM), flush=True)
    null = {k: [] for k in ("A", "B", "C", "D")}
    for i in range(N_PERMUTATIONS):
        r2 = np.random.default_rng(RNG_SEED + 1000 + i)
        yp = np.zeros_like(y)
        pick = r2.integers(0, K + 1, size=n)
        yp[np.arange(n)[pick == 0]] = 1.0
        nz = pick > 0
        yp[n + np.arange(n)[nz] * K + (pick[nz] - 1)] = 1.0
        for k in null:
            a, _ = L.fit_score(X[k][tr], yp[tr], X[k][te], yp[te], i + 1)
            null[k].append(a)
        if (i + 1) % 25 == 0:
            print("    %d/%d  (%.0f s elapsed)"
                  % (i + 1, N_PERMUTATIONS, time.time() - t_start), flush=True)
    null = {k: np.asarray(v) for k, v in null.items()}

    def pval(nd, d):
        return float((np.sum(nd >= d) + 1) / (nd.size + 1))

    nd = {
        "B_vs_A": null["B"] - null["A"],
        "C_vs_A": null["C"] - null["A"], "D_vs_B": null["D"] - null["B"],
        "D_vs_A": null["D"] - null["A"],
    }
    contrasts = {}
    for k, v in nd.items():
        contrasts[k] = {
            "observed_dAUC": d_obs[k],
            "null_mean": float(v.mean()), "null_sd": float(v.std(ddof=1)),
            "null_p95": float(np.quantile(v, 0.95)),
            "p_one_sided": pval(v, d_obs[k]),
            "n_permutations": int(v.size),
        }
    # the committed arm used exactly 40 permutations with these seeds; quote both
    v40 = nd["B_vs_A"][:COMMITTED_N_PERM]
    repro.update({
        "committed_p_dAUC_B_vs_A": EXPECT["committed_p_dauc"],
        "reproduced_p_first_%d_permutations" % COMMITTED_N_PERM:
            pval(v40, d_obs["B_vs_A"]),
        "reproduced_null_mean_first_%d" % COMMITTED_N_PERM: float(v40.mean()),
        "reproduced_null_sd_first_%d" % COMMITTED_N_PERM: float(v40.std(ddof=1)),
    })

    for k in ("B_vs_A", "C_vs_A", "D_vs_B"):
        c = contrasts[k]
        print("  %-8s dAUC %+.6f   null %+.6f +/- %.6f   p = %.4f"
              % (k, c["observed_dAUC"], c["null_mean"], c["null_sd"],
                 c["p_one_sided"]), flush=True)

    # ------------------------------------------------------------- P-2.4 placebo
    print("\nP-2.4 placebo survivor rate over the block-3 features ...", flush=True)
    pnames = list(F.BLOCK3_VALUES) + list(TIDAL)
    pobs, pnull, pnu = placebo_Z(allcols, pnames, member)
    ix = {nm: i for i, nm in enumerate(pnames)}
    plac_all = survivor_report(pnames, [ix[c] for c in F.BLOCK3_VALUES], pobs, pnull, pnu)
    plac_pl = survivor_report(pnames, [ix[c] for c in BLOCK3_PLACEBO], pobs, pnull, pnu)
    plac_me = survivor_report(pnames, [ix[c] for c in BLOCK3_MECHANICAL], pobs, pnull, pnu)
    plac_tid = survivor_report(pnames, [ix[c] for c in TIDAL], pobs, pnull, pnu)
    print("  block3 (all %d features): null survivor rate %.4f per feature, "
          "%.3f survivors per draw, family-wise any-survivor %.3f"
          % (plac_all["n_features"],
             plac_all["null_mean_survivor_rate_over_features"],
             plac_all["null_mean_survivors_per_draw"],
             plac_all["null_familywise_any_survivor_rate"]), flush=True)
    print("  observed survivors on real labels: block3 %d/%d, placebo-half %d/%d, "
          "mechanical-half %d/%d, tidal %d/%d"
          % (plac_all["observed_n_survivors"], plac_all["n_features"],
             plac_pl["observed_n_survivors"], plac_pl["n_features"],
             plac_me["observed_n_survivors"], plac_me["n_features"],
             plac_tid["observed_n_survivors"], plac_tid["n_features"]), flush=True)

    # ------------------------------------------------------------------- output
    runtime = time.time() - t_start
    out = {
        "arm": ("Phase 2 step 1: the learned arm extended with a non-tidal forcing block "
                "(LOD / polar motion / Kp / ap / OMNI2), plus the P-2.4 placebo reading"),
        "state_class": "first-run, exploration split",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 1,
        "priced_what": ("one declared statistic per new block: dAUC(C vs A) and "
                        "dAUC(D vs B). B vs A is a REPRODUCTION of the committed P-1.4 "
                        "number, not a new test."),
        "holdout_touched": False,
        "catalogue": "QTM_decluster_m0.1 exploration split",
        "controls_per_case": int(K),
        "null_half_window_days": float(NULL_HALF_WINDOW_DAYS),
        "train_frac": TRAIN_FRAC,
        "rng_seed": RNG_SEED,
        "counted_invariants": inv,
        "forcing_load_audit": load_audit,
        "block3_fill_report": fill,
        "schedule_audit": sched,
        "schedule_audit_note": (
            "ARTIFACT 1 measured: the UT-diurnal and weekday amplitudes of each index over "
            "the arm's own span, in units of that index's own sd. A large UT-diurnal "
            "amplitude is the channel by which a UT-gridded index could proxy the "
            "catalogue's local-solar-hour detection artifact; model A holds solar hour, "
            "day of week and day of year, so that channel is already priced into the "
            "nested comparison."),
        "feature_sets": FEAT,
        "models": obs,
        "observed_delta_auc": d_obs,
        "contrasts": contrasts,
        "reproduction_check": repro,
        "placebo": {
            "question": ("P-2.4: how many survivors at nominal alpha does this pipeline "
                         "manufacture on a property class with real temporal structure and "
                         "a near-zero physical prior?"),
            "statistic": ("per feature, the matched-set conditional z: sum over strata of "
                          "the case's within-stratum standardised value, over sqrt(n "
                          "usable strata). Unit normal under H0 by construction."),
            "block3_all": plac_all,
            "block3_placebo_half_kp_ap_omni": plac_pl,
            "block3_mechanical_half_lod_polar": plac_me,
            "tidal_block_for_contrast": plac_tid,
        },
        "runtime_seconds": runtime,
        "n_permutations": N_PERMUTATIONS,
    }

    pC, pD, pDA = (contrasts["C_vs_A"]["p_one_sided"], contrasts["D_vs_B"]["p_one_sided"],
                   contrasts["D_vs_A"]["p_one_sided"])
    out["verdict"] = (
        "C vs A: dAUC = %+.5f, p = %.4f. D vs B: dAUC = %+.5f, p = %.4f. "
        "D vs A: dAUC = %+.5f, p = %.4f. %s "
        "P-2.4 placebo: nominal-alpha survivor rate %.4f per feature over the block-3 "
        "features, family-wise P(any survivor) %.3f over %d features -- the calibration "
        "for what nominal alpha buys on a property class with real temporal structure and "
        "a near-zero physical prior."
        % (d_obs["C_vs_A"], pC, d_obs["D_vs_B"], pD, d_obs["D_vs_A"], pDA,
           ("BOTH DECLARED CONTRASTS NULL. Adding non-tidal forcing to the observer block "
            "does not improve out-of-sample discrimination, on a design whose sensitivity "
            "is measured separately in exp_learned_ext_sensitivity.py."
            if (pC > 0.05 and pD > 0.05) else
            "AT LEAST ONE CONTRAST FIRED. Read it against the placebo survivor rate below "
            "before treating it as anything but a false positive: the placebo half of this "
            "block has a near-zero physical prior by construction."),
           plac_all["null_mean_survivor_rate_over_features"],
           plac_all["null_familywise_any_survivor_rate"], plac_all["n_features"]))

    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print("\n" + "=" * 78)
    print(out["verdict"])
    print("\nruntime %.1f s   wrote %s" % (runtime, OUT_JSON.name))


if __name__ == "__main__":
    main()
