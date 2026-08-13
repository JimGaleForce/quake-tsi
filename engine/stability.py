"""§P7-14(b) -- the PER-CANDIDATE resampling-stability line. LABEL, DO NOT KILL.

WHERE THIS COMES FROM, and why an audit turned into a condition
---------------------------------------------------------------
Tranche A measured two stability numbers that look contradictory and are not: the
ranking was robust to the RNG (R4 reseeding rho 0.926 / 0.922) and mobile under the
data (F10-24 moving-block resampling rho **0.487 / 0.044**). §P7-14(b) reads it
correctly and then does something with it:

  > *"With zero survivors, the tranche is ranking noise, and noise is SUPPOSED to
  > reorder under resampling. A high F10-24 here would have been the anomaly... This
  > converts F10-24 from an audit into a prediction with teeth for B."*
  >
  > **CONDITION ON TRANCHE B.** Every candidate that surfaces in B carries a
  > MANDATORY resampling-stability line: (i) R4 reseeding stability, (ii) F10-24
  > data-resampling selection frequency in the top-k, (iii) the explicit statement
  > that (ii) is a **LOWER** bound. **Label, do not kill.** ... **And the sharp
  > version:** Tranche A has supplied the NULL BASELINE for this statistic
  > (0.487 / 0.044). A genuine candidate in B should show resampling stability
  > MARKEDLY ABOVE that baseline; one that does not is indistinguishable from the
  > noise Tranche A just characterised. **The R4 pair alone is NOT sufficient** -- it
  > measures the RNG, and the RNG was never the threat.

THREE THINGS THIS MODULE REFUSES TO DO, each for a stated reason
----------------------------------------------------------------
1. **It never deletes a candidate.** The metric is biased LOW by construction (the
   moving-block resample perturbs the dependence structure, so measured stability is
   a lower bound on true stability), and a hard kill on a biased-low statistic
   discards true positives. `UNSTABLE` is a TAG that travels with the row; §P7-14(b)
   forbids writing such a row up as a K-entry WITHOUT the tag, not writing it up.
2. **It never reports (ii) without the lower-bound sentence.** That sentence is not
   decoration: "the ranking moves with the data" is the claim the evidence supports,
   and "the ranking is noise" is not, and the difference is exactly the bound.
3. **It never accepts the R4 pair as sufficient.** `stability_line` requires at least
   one data-resampling arm and raises without one, because the failure mode is a
   report that shows two reassuring reseeding numbers and no data arm at all.

THE THRESHOLDS ARE DECLARED HERE, BEFORE ANY RUN (S-9)
------------------------------------------------------
§P7-14(b) says "markedly above" and leaves the number to the build. With exactly two
moving-block resamples the selection frequency lives in {0.0, 0.5, 1.0}, so the only
honest reading of "markedly above 0.044" is BOTH: a candidate must re-enter the top-k
under BOTH resamples. And rank stability must clear 0.70 against the 0.487 null. Both
are frozen in `THRESHOLDS` and both are stated in the output of every line, so a
reader never has to ask what bar was applied or when it was set.
"""

from __future__ import annotations

import numpy as np

# §P7-14(b): Tranche A's measured null baseline for this statistic. These are the
# numbers a Tranche B candidate must be MARKEDLY ABOVE to be distinguishable from
# the noise Tranche A characterised. Source: results_tranche_a.json,
# rank_stability.F10_24_data_resampling, commit 497b196.
NULL_BASELINE = {
    "f10_24_rho_pair": (0.487, 0.044),
    "f10_24_rho_mean": 0.2655,
    "r4_rho_pair": (0.926, 0.922),
    "source": ("HYPOTHESIS_LEDGER.md §P7-14(b) / results_tranche_a.json "
               "(session_20260813T001742, commit 497b196)"),
}

# DECLARED BEFORE THE RUN (S-9). One value each, no alternatives.
THRESHOLDS = {
    "min_selection_frequency": 1.0,     # must re-enter the top-k under BOTH resamples
    "min_resample_rho": 0.70,           # markedly above the 0.487 null arm
    "top_k": 100,
    "declared": ("frozen in engine/stability.py before Tranche B's declaration; "
                 "§P7-14(b) asks for 'markedly above' the 0.487/0.044 baseline and "
                 "these are that phrase made into numbers, in advance"),
}

LOWER_BOUND_SENTENCE = (
    "(iii) THIS IS A LOWER BOUND. The moving-block resample perturbs the dependence "
    "structure itself, so injected structure-noise pushes the measured stability "
    "DOWN: the true stability is at least this. 'The ranking moves with the data' is "
    "what this supports; 'the ranking is noise' is not, and the difference is "
    "exactly this bound (§P7-14(b), the measurer's own bounding, ratified).")

LABEL_NOT_KILL = (
    "LABEL, DO NOT KILL (§P7-14(b)). A candidate below the declared stability "
    "threshold is TAGGED UNSTABLE and may not be written up as a K-entry without the "
    "tag. It is NOT deleted: the metric is biased low by construction and a hard "
    "kill on a biased-low statistic would discard true positives.")

R4_NOT_SUFFICIENT = (
    "The R4 reseeding pair alone is NOT sufficient (§P7-14(b)): it measures whether "
    "the RNG moves the ranking, and the RNG was never the threat. A stability line "
    "without a DATA-resampling arm is not a stability line.")

TAG_STABLE = "STABLE"
TAG_UNSTABLE = "UNSTABLE"
TAG_BORDERLINE = "BORDERLINE"


class NoDataResamplingArm(AssertionError):
    """A stability line was asked for with reseeding arms only (§P7-14(b))."""


def candidate_key(row):
    """The identity a candidate keeps across resamples: (test, feature, lag, mark)."""
    return (row["test"], row["feature"], row.get("lag"), row.get("mark"))


def _ranked(tests):
    """Rows ordered as the session orders them, mapped to rank. Ties broken as in
    `mine_session.rank_stability`, so the two never disagree about a ranking."""
    s = sorted(tests, key=lambda t: (float(t.get("p_bh", t["p_raw"])),
                                     t.get("order_key", [])))
    return {candidate_key(t): i for i, t in enumerate(s)}


def _spearman(a, b):
    from scipy import stats
    return float(stats.spearmanr(np.asarray(a, float),
                                 np.asarray(b, float)).statistic)


def stability_line(candidate, ref_tests, reseed_runs, resample_runs,
                   top_k=None, thresholds=None):
    """The mandatory §P7-14(b) line for ONE candidate. Returns the report dict.

    `reseed_runs` / `resample_runs` are lists of test-row lists -- one per audit run.
    §P7-14(b) asks for the moving-block resample x2, so `resample_runs` normally has
    length 2; the code takes whatever it is given and states how many it used.
    """
    th = dict(THRESHOLDS if thresholds is None else thresholds)
    k = int(top_k or th["top_k"])
    if not resample_runs:
        raise NoDataResamplingArm(R4_NOT_SUFFICIENT)

    key = candidate_key(candidate) if isinstance(candidate, dict) else tuple(candidate)
    ref = _ranked(ref_tests)
    if key not in ref:
        raise KeyError("candidate %r is not in the reference test set" % (key,))
    ref_rank = ref[key]
    ref_top = [x for x, v in sorted(ref.items(), key=lambda kv: kv[1])[:k]]

    def _arm(runs, label):
        rows, in_top = [], 0
        for i, ts in enumerate(runs):
            r = _ranked(ts)
            rank = r.get(key)
            common = [x for x in ref_top if x in r]
            rho = (_spearman([ref[x] for x in common], [r[x] for x in common])
                   if len(common) >= 3 else None)
            present = rank is not None and rank < k
            in_top += int(present)
            rows.append({"run": i, "label": "%s[%d]" % (label, i),
                         "candidate_rank": (None if rank is None else int(rank)),
                         "in_top_k": bool(present),
                         "n_common_top_k": len(common),
                         "top_k_spearman_rho": rho})
        freq = in_top / float(max(len(runs), 1))
        rhos = [r["top_k_spearman_rho"] for r in rows
                if r["top_k_spearman_rho"] is not None]
        return {"runs": rows, "n_runs": len(runs), "selection_frequency": freq,
                "mean_top_k_rho": (float(np.mean(rhos)) if rhos else None)}

    r4 = _arm(reseed_runs, "R4_reseed") if reseed_runs else {
        "runs": [], "n_runs": 0, "selection_frequency": None, "mean_top_k_rho": None}
    f1024 = _arm(resample_runs, "F10-24_resample")

    freq = f1024["selection_frequency"]
    rho = f1024["mean_top_k_rho"]
    pass_freq = freq is not None and freq >= float(th["min_selection_frequency"])
    pass_rho = rho is not None and rho >= float(th["min_resample_rho"])
    tag = (TAG_STABLE if (pass_freq and pass_rho)
           else TAG_BORDERLINE if (pass_freq or pass_rho) else TAG_UNSTABLE)

    return {
        "candidate": {"test": key[0], "feature": key[1], "lag": key[2],
                      "mark": key[3], "reference_rank": int(ref_rank)},
        "top_k": k,
        "i_R4_reseeding": r4,
        "ii_F10_24_data_resampling": f1024,
        "iii_lower_bound": LOWER_BOUND_SENTENCE,
        "null_baseline": dict(NULL_BASELINE),
        "thresholds": th,
        "vs_null_baseline": {
            "selection_frequency_vs_0.044": (
                None if freq is None else float(freq) - 0.044),
            "resample_rho_vs_0.487": (None if rho is None else float(rho) - 0.487),
            "reading": ("MARKEDLY ABOVE the Tranche A null baseline"
                        if (pass_freq and pass_rho) else
                        "NOT markedly above the Tranche A null baseline: "
                        "indistinguishable from the noise Tranche A characterised"),
        },
        "tag": tag,
        "label_not_kill": LABEL_NOT_KILL,
        "r4_not_sufficient": R4_NOT_SUFFICIENT,
        "headline": ("%s | R4 rho %s | F10-24 selection freq %s, rho %s "
                     "(LOWER BOUND) | null baseline 0.487/0.044"
                     % (tag,
                        "n/a" if r4["mean_top_k_rho"] is None
                        else "%.3f" % r4["mean_top_k_rho"],
                        "n/a" if freq is None else "%.2f" % freq,
                        "n/a" if rho is None else "%.3f" % rho)),
    }


def stability_block(candidates, ref_tests, reseed_runs, resample_runs,
                    top_k=None, thresholds=None):
    """The stability line for EVERY candidate that must carry one.

    §P7-14(b) says "every candidate that surfaces in B", which is every BH survivor
    and the max-statistic detection. `require_lines_for` derives that set from a
    session's rows so no caller has to remember which rows qualify.
    """
    return [stability_line(c, ref_tests, reseed_runs, resample_runs, top_k,
                           thresholds) for c in candidates]


def require_lines_for(tests, max_stat_key=None):
    """The rows that MUST carry a stability line: BH survivors + the max-stat row."""
    out = [t for t in tests if bool(t.get("pass_bh") or t.get("bh_pass"))]
    keys = {candidate_key(t) for t in out}
    if max_stat_key is not None and tuple(max_stat_key) not in keys:
        for t in tests:
            if candidate_key(t) == tuple(max_stat_key):
                out.append(t)
                break
    return out


def assert_lines_present(tests, lines, max_stat_key=None):
    """Refuse a report whose candidates are not all carrying their line."""
    need = {candidate_key(t) for t in require_lines_for(tests, max_stat_key)}
    have = {(l["candidate"]["test"], l["candidate"]["feature"],
             l["candidate"]["lag"], l["candidate"]["mark"]) for l in lines}
    missing = need - have
    if missing:
        raise AssertionError(
            "§P7-14(b): %d candidate(s) would be reported without a "
            "resampling-stability line: %s. The line is MANDATORY on every BH "
            "survivor and on the max-statistic detection."
            % (len(missing), sorted(str(m) for m in missing)))
    return {"n_candidates": len(need), "n_lines": len(lines), "satisfied": True}
