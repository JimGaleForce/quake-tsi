"""P-1.1 THE HIGH-N REPLAY: the same hypotheses, an instrument that can finally see. Priced 1.

This is not a new idea. It is EVERY idea this program has already tested, run for the first
time at a sample size where the claimed effect is detectable.

THE ARITHMETIC THAT MAKES IT THE TOP OF THE QUEUE. Every arm to date has used 7,139
declustered events, giving a detection floor near 5% relative modulation. The effect
claimed in the literature for ordinary seismicity is around 1%. **This program has never
once run a test powered at the size the thing is claimed to be**, so every null it has
produced is consistent with the effect being real and simply smaller than the instrument.
That is not a finding about the world; it is a finding about N.

The Lu/Xue release has been on disk since 2026-07-21 and holds two ETAS-declustered
catalogues: QTM at M >= 0.1 (45,069 events, 2008-2017) and SCSN at M >= 1.5 (50,313
events, 1981-2018). Together, after the 70% exploration split, roughly 66,000 independent
events. Floor about 1.7%. That is the first time this program crosses into the range where
the literature's own claim lives, and it required no download at all.

WHY A PER-EVENT NULL IS LEGITIMATE HERE, having just been shown to be catastrophic
elsewhere. `exp_nearcritical.py` returned max |z| = 26.27 from a per-event null because it
deliberately KEPT clustered aftershocks, and co-temporal events are not independent samples
of a slowly-varying global quantity. These catalogues are ETAS-DECLUSTERED BY THEIR
AUTHORS: one event per sequence, which is exactly the condition the sequence-level null was
constructed to restore. The dataset's own frozen rule -- "unit of inference: sequence, not
event" -- is satisfied because here the two coincide.

THE ARTIFACT THAT COULD FAKE THIS, NAMED BEFORE THE RUN. At M >= 0.1 the QTM catalogue sits
near its completeness threshold, and if DETECTION is tidally modulated -- through noise,
through station uptime, through coda blindness -- then a tidal signal in the catalogue is a
signal in the instrument, not in the Earth. We have no measured Mc(x,t) field. So the arm
runs a magnitude ladder as its own control:

    QTM M >= 0.1   near completeness, maximum N, maximum artifact exposure
    QTM M >= 1.5   well above completeness, ~9x fewer events
    QTM M >= 2.5   far above completeness
    SCSN M >= 1.5  a different network, different era, different processing

**A real geophysical effect should persist up the ladder and across the two networks. A
detection artifact should die as the magnitude floor rises, and need not agree between
networks.** That contrast is the actual test, and it is declared here rather than invented
afterwards. If M >= 0.1 fires and M >= 2.5 does not, the honest reading is the artifact.

BOTH MULTIPLICITY CORRECTIONS ARE PAID: cell max over every catalogue-subpopulation-feature
cell, and the pooled max over features across subpopulations, because a cell max is blind to
a weak effect present in every subpopulation at once.

PRICED 1: one declared family, one max-statistic, one p. Exploration split only.
"""

from __future__ import annotations

import datetime as _dt
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import exp_mass_screen as MS
import exp_world_harmonics as W

OUT_JSON = HERE / "results_highn.json"
ZEN = HERE / "data" / "xue_lu_zenodo"

EXPLORE_FRAC = 0.70
NULL_HALF_WINDOW_DAYS = MS.FORTNIGHT_DAYS
S_NULL_PER_EVENT = 120
N_REPLICATES = 3000
REP_CHUNK = 500
RNG_SEED = 20260822
MIN_SUBPOP = 200

SCALARS = MS.SCALARS
ANGLES = MS.ANGLES
ORDERS = MS.ORDERS


def load_zenodo(path):
    """`year mo dy hr mn sec lat lon depth_km mag event_id`, whitespace-delimited."""
    t, la, lo, dp, mg = [], [], [], [], []
    with open(path, encoding="utf-8", errors="replace") as f:
        for ln in f:
            p = ln.split()
            if len(p) < 10:
                continue
            try:
                sec = float(p[5])
                base = _dt.datetime(int(p[0]), int(p[1]), int(p[2]),
                                    int(p[3]), int(p[4]), 0,
                                    tzinfo=_dt.timezone.utc)
                ts = base + _dt.timedelta(seconds=sec)
                t.append((ts - W.SPAN_START).total_seconds() / 86400.0)
                la.append(float(p[6])); lo.append(float(p[7]))
                dp.append(float(p[8])); mg.append(float(p[9]))
            except (ValueError, OverflowError):
                continue
    o = np.argsort(np.asarray(t))
    return (np.asarray(t)[o], np.asarray(la)[o], np.asarray(lo)[o],
            np.asarray(dp)[o], np.asarray(mg)[o])


def split(t, *arrs):
    cut = t.min() + (t.max() - t.min()) * EXPLORE_FRAC
    m = t < cut
    return (t[m],) + tuple(a[m] for a in arrs) + (int((~m).sum()),)


def null_mean_chunked(pool_arr, idx_all, mask, denom):
    """Mean over masked events of the selected null draw, evaluated in replicate chunks
    so that an [n_events x n_replicates] gather is never materialised at full size."""
    sub = pool_arr[mask]
    out = np.empty(idx_all.shape[1])
    for a in range(0, idx_all.shape[1], REP_CHUNK):
        b = min(a + REP_CHUNK, idx_all.shape[1])
        out[a:b] = np.take_along_axis(sub, idx_all[mask, a:b], axis=1).sum(axis=0) / denom
    return out


def main():
    rng = np.random.default_rng(RNG_SEED)
    running_max = np.zeros(N_REPLICATES)
    records = []
    pooled = {}
    catalogues = []

    specs = [
        ("QTM_declustered", ZEN / "QTM_decluster_m0.1.txt",
         [("M_ge0.1_ALL", 0.1), ("M_ge1.5_control", 1.5), ("M_ge2.5_control", 2.5)]),
        ("SCSN_declustered", ZEN / "SCSN_decluster_m1.5.txt",
         [("M_ge1.5_ALL", 1.5), ("M_ge2.5_control", 2.5)]),
    ]

    for cat_name, path, bands in specs:
        t, la, lo, dp, mg = load_zenodo(path)
        t, la, lo, dp, mg, n_hold = split(t, la, lo, dp, mg)
        print("%s: %d exploration events (%d holdout reserved), M %.2f..%.2f"
              % (cat_name, t.size, n_hold, mg.min(), mg.max()), flush=True)
        catalogues.append({"name": cat_name, "n_explore": int(t.size),
                           "n_holdout_reserved": int(n_hold)})

        obs_raw = MS.raw_series(t, la, lo)
        off = rng.uniform(-NULL_HALF_WINDOW_DAYS, NULL_HALF_WINDOW_DAYS,
                          size=(t.size, S_NULL_PER_EVENT))
        pool = MS.raw_series((t[:, None] + off).ravel(),
                             np.repeat(la, S_NULL_PER_EVENT),
                             np.repeat(lo, S_NULL_PER_EVENT))
        pool = {k: v.reshape(t.size, S_NULL_PER_EVENT) for k, v in pool.items()}
        del off
        idx = rng.integers(0, S_NULL_PER_EVENT, size=(t.size, N_REPLICATES))

        pops = []
        for bname, mlo in bands:
            m = mg >= mlo
            if int(m.sum()) >= MIN_SUBPOP:
                pops.append(("%s|%s" % (cat_name, bname), m))
                print("    %-34s n=%6d" % (bname, int(m.sum())), flush=True)

        def emit(feature, pop, mask, obs_val, null_col):
            mu, sd = float(null_col.mean()), float(null_col.std(ddof=1))
            if not np.isfinite(sd) or sd <= 0:
                return
            z = (float(obs_val) - mu) / sd
            zn = (null_col - mu) / sd
            np.maximum(running_max, np.abs(zn), out=running_max)
            records.append({"cell": pop, "feature": feature,
                            "n": int(mask.sum()), "z": z})
            pooled.setdefault(feature, {"z": [], "zn": [], "w": []})
            pooled[feature]["z"].append(z)
            pooled[feature]["zn"].append(zn)
            pooled[feature]["w"].append(np.sqrt(float(mask.sum())))

        for s in SCALARS:
            thr_hi = float(np.quantile(pool[s], 0.90))
            thr_lo = float(np.quantile(pool[s], 0.10))
            hi_p = (pool[s] > thr_hi).astype(np.float64)
            lo_p = (pool[s] < thr_lo).astype(np.float64)
            o = obs_raw[s]
            for pn, m in pops:
                d = float(m.sum())
                emit("%s.mean" % s, pn, m, o[m].mean(),
                     null_mean_chunked(pool[s], idx, m, d))
                emit("%s.top_decile" % s, pn, m, (o[m] > thr_hi).mean(),
                     null_mean_chunked(hi_p, idx, m, d))
                emit("%s.bot_decile" % s, pn, m, (o[m] < thr_lo).mean(),
                     null_mean_chunked(lo_p, idx, m, d))
            del hi_p, lo_p
        print("    scalars done (%d cells)" % len(records), flush=True)

        for a in ANGLES:
            th_p, th_o = pool[a], obs_raw[a]
            for mm in ORDERS:
                cp = np.cos(mm * th_p)
                sp = np.sin(mm * th_p)
                for pn, m in pops:
                    d = float(m.sum())
                    gc = null_mean_chunked(cp, idx, m, d)
                    gs = null_mean_chunked(sp, idx, m, d)
                    emit("%s.R%d" % (a, mm), pn, m,
                         np.hypot(np.cos(mm * th_o[m]).mean(),
                                  np.sin(mm * th_o[m]).mean()),
                         np.hypot(gc, gs))
                del cp, sp
        print("    angles done (%d cells total)" % len(records), flush=True)
        del pool, obs_raw, idx

    obs_max = max(abs(r["z"]) for r in records)
    fam_p = float((np.sum(running_max >= obs_max) + 1) / (N_REPLICATES + 1))
    crit95 = float(np.quantile(running_max, 0.95))
    records.sort(key=lambda r: -abs(r["z"]))

    pool_running = np.zeros(N_REPLICATES)
    prows = []
    for feat, d in pooled.items():
        if len(d["z"]) < 3:
            continue
        w = np.asarray(d["w"]); den = np.sqrt((w * w).sum())
        po = float((w * np.asarray(d["z"])).sum() / den)
        pn = (w[:, None] * np.stack(d["zn"])).sum(axis=0) / den
        np.maximum(pool_running, np.abs(pn), out=pool_running)
        prows.append({"feature": feat, "pooled_z": po, "n_cells": len(d["z"])})
    prows.sort(key=lambda r: -abs(r["pooled_z"]))
    pool_obs = abs(prows[0]["pooled_z"]) if prows else float("nan")
    pool_p = float((np.sum(pool_running >= pool_obs) + 1) / (N_REPLICATES + 1))

    # --- the declared ladder read: does the top feature survive up the magnitude floor?
    top_feat = records[0]["feature"]
    ladder = {r["cell"]: r["z"] for r in records if r["feature"] == top_feat}

    n_tot = sum(c["n_explore"] for c in catalogues)
    out = {
        "arm": "P-1.1 high-N replay on the Lu/Xue declustered catalogues",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 1, "holdout_touched": False,
        "catalogues": catalogues,
        "n_events_total_explore": int(n_tot),
        "detection_floor_estimate_relative": float(4.29 / np.sqrt(n_tot)),
        "prior_arms_used": 7139,
        "prior_floor_estimate_relative": float(4.29 / np.sqrt(7139)),
        "null": {"kind": "waveform-matched dwell-time, per-event",
                 "why_per_event_is_legitimate_here": (
                     "these catalogues are ETAS-declustered by their authors, so one "
                     "event per sequence and the dataset's 'unit of inference: sequence' "
                     "rule is satisfied because sequence and event coincide"),
                 "half_window_days": NULL_HALF_WINDOW_DAYS,
                 "draws_per_event": S_NULL_PER_EVENT, "replicates": N_REPLICATES},
        "n_cells": len(records),
        "cell_max_family": {"observed_max_abs_z": obs_max, "null_max_p95": crit95,
                            "p": fam_p,
                            "where": "%s / %s" % (records[0]["cell"],
                                                  records[0]["feature"])},
        "pooled_max_family": {"observed_max_abs_pooled_z": pool_obs,
                              "null_max_p95": float(np.quantile(pool_running, 0.95)),
                              "p": pool_p,
                              "where": prows[0]["feature"] if prows else None},
        "magnitude_ladder_for_top_feature": {"feature": top_feat, "z_by_cell": ladder,
                                             "read": ("a real effect persists as the "
                                                      "magnitude floor rises and agrees "
                                                      "across the two networks; a "
                                                      "detection artifact dies")},
        "top_30_cells": records[:30],
        "top_15_pooled": prows[:15],
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    print("\n" + "=" * 78)
    print("HIGH-N REPLAY: %d cells, %d exploration events (was 7,139)"
          % (len(records), n_tot))
    print("  detection floor  %.2f%%   (previous arms %.2f%%)"
          % (100 * out["detection_floor_estimate_relative"],
             100 * out["prior_floor_estimate_relative"]))
    print("  CELL   max |z| = %.3f at %s ; null 95th %.3f  ->  FAMILY p = %.4f"
          % (obs_max, out["cell_max_family"]["where"], crit95, fam_p))
    print("  POOLED max |z| = %.3f at %s ; null 95th %.3f  ->  p = %.4f"
          % (pool_obs, out["pooled_max_family"]["where"],
             out["pooled_max_family"]["null_max_p95"], pool_p))
    print("\n  magnitude ladder for the top feature (%s):" % top_feat)
    for k, v in sorted(ladder.items()):
        print("    %-40s z = %+.3f" % (k, v))
    print("\n  top 12 cells:")
    for r in records[:12]:
        print("    %+6.2f  %-34s %-24s n=%d"
              % (r["z"], r["cell"], r["feature"], r["n"]))
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
