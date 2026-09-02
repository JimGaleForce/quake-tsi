"""THE NEAR-CRITICAL ARM: test the faults we KNOW are close to failure. Priced 1.

The program has spent its whole life measuring a diluted quantity and calling the
dilution a null.

THE DILUTION ARGUMENT, WHICH IS ARITHMETIC AND NOT A HOPE. A main-effect test on a whole
catalogue does not measure the tidal susceptibility s. It measures s * f, where f is the
fraction of faults sampled that are close enough to failure for a kilopascal to matter. If
f is 0.05 to 0.15 -- and every physical estimate says something like that -- then a true
susceptibility of 20 percent shows up as 1 to 3 percent, under every detection floor this
program has ever achieved. **Every null we have produced is consistent with a real and
large effect seen through a factor of ten of dilution.**

So stop averaging over faults that are nowhere near failure. AFTERSHOCKS ARE THE ONE
POPULATION WE KNOW IS NEAR FAILURE, because they just failed. The program has spent its
life DISCARDING them on a reflex that treats the clustered signal as false by assumption.
That reflex was never tested. It is the target here.

=============================================================================
THE NULL HAD TO BE REBUILT, AND THE FIRST VERSION OF THIS FILE WAS WRONG.
=============================================================================

Run with the program's standard PER-EVENT dwell null, this arm returned max |z| = 26.27
at `aftershock_lt1d / elongation.R1`, with a declared family p of 0.0002. That is not a
discovery. It is the single largest artifact this program has produced, and it is
recorded here in full rather than quietly deleted, because the mechanism is instructive:

  * Every one of the top cells was `elongation` or `elong_cos_abs` -- the Moon-Sun angle,
    a GLOBAL quantity with a 29.53 day period and no site dependence at all.
  * The 9,028 events in `aftershock_lt1d` are not 9,028 independent samples of it. They
    come from a small number of sequences, and every event in a sequence shares
    essentially the same elongation, because a 29.53 day cycle does not move during the
    day after a mainshock.
  * The per-event null draws each event's comparison times INDEPENDENTLY, so it believes
    it is looking at 9,028 independent samples. Its standard deviation is therefore too
    small by roughly the square root of the ratio of events to sequences, and z inflates
    by that same factor. z = 26 is that ratio, not a signal.

This is exactly the dependence the program declusters to avoid, met head-on because this
arm deliberately KEEPS the dependent events. The fix is not to decluster -- that would
throw away the whole point -- but to make the null respect the dependence:

**ONE OFFSET PER SEQUENCE, APPLIED RIGIDLY TO EVERY EVENT IN IT.** The null shifts whole
sequences in time rather than shuffling events within them. Within-sequence structure --
Omori decay, spatial footprint, magnitude distribution, and above all the near-identical
tidal phase shared by co-temporal events -- is preserved exactly. The effective sample
size becomes the NUMBER OF SEQUENCES, which is what it always physically was.

The honest cost is stated rather than hidden: effective N falls from ~14,000 events to
~1,000 sequences, so the detection floor RISES to roughly 14 percent. That is still the
best-placed test this program has run, because it is 14 percent of an UNDILUTED
susceptibility rather than 5 percent of one diluted by a factor of ten.

WHAT IS DECLARED, AND THE GRADIENT IS THE REAL TEST. A single significant cell is a
lottery ticket. Aftershocks are split by time since mainshock (under 1 day, 1 to 7 days,
7 to 30 days); stress is most elevated and the population most critical immediately after
and decays thereafter, so a mechanism predicts a MONOTONE DECLINE across those bands with
the mainshock/background population flat at the bottom. The differential
(aftershock minus mainshock, same catalogue, same region, same epoch, same network, same
code path) is declared alongside, because a detection artifact, an ephemeris error or a
dwell-null defect hits both populations alike and cannot produce a contrast.

THE ARTIFACT THAT COULD STILL FAKE THIS, NAMED IN ADVANCE. Aftershock catalogues are
incomplete right after the mainshock: coda blindness swallows small events for minutes to
hours. It does not alias onto tidal phase, being keyed to time-since-mainshock which is
uniform with respect to lunar time, but a magnitude-restricted arm (M >= 3.0, above
completeness even in coda) is run alongside so a survivor is tested rather than argued.

Exploration split only; the SoCal holdout is not touched.
PRICED 1: one declared family, one max-statistic over every cell, one p.
"""

from __future__ import annotations

import csv
import datetime as _dt
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import exp_mass_screen as MS
import exp_world_harmonics as W

OUT_JSON = HERE / "results_nearcritical.json"
CATALOG = HERE / "data" / "comcat_socal_m25.csv"

MAG_MIN = 2.5
EXPLORE_FRAC = 0.70
DECLUSTER_DAYS = 30.0
DECLUSTER_KM = 150.0
NULL_HALF_WINDOW_DAYS = MS.FORTNIGHT_DAYS
S_NULL_PER_SEQUENCE = 200
N_REPLICATES = 4000
RNG_SEED = 20260822
MIN_SUBPOP = 50
MIN_SEQUENCES = 30

SCALARS = ("areal", "areal_abs", "shear", "rate", "rate_abs", "elong_cos_abs")
ANGLES = ("bearing", "moon_hourangle", "moon_az", "elongation")
ORDERS = (1, 2, 3, 4)


def load():
    rows = []
    for r in csv.DictReader(open(CATALOG, newline="", encoding="utf-8")):
        try:
            m = float(r["mag"])
            if m < MAG_MIN:
                continue
            ts = _dt.datetime.fromisoformat(r["time"].replace("Z", "+00:00"))
            rows.append((ts, float(r["latitude"]), float(r["longitude"]),
                         float(r["depth"] or 0.0), m))
        except (ValueError, KeyError, TypeError):
            continue
    rows.sort(key=lambda x: x[0])
    t0, t1 = rows[0][0], rows[-1][0]
    cut = t0 + (t1 - t0) * EXPLORE_FRAC
    ex = [r for r in rows if r[0] < cut]
    # EPOCH INVARIANT: days since 1970-01-01Z, as MS.raw_series requires.
    t = np.array([r[0].timestamp() / 86400.0 for r in ex])
    return (t, np.array([r[1] for r in ex]), np.array([r[2] for r in ex]),
            np.array([r[3] for r in ex]), np.array([r[4] for r in ex]),
            len(rows) - len(ex), cut)


def classify(t, lat, lon, mag):
    """Aftershock flag, days since the triggering mainshock, and a SEQUENCE ID.

    The sequence id is the root of the triggering chain: an aftershock inherits the
    sequence of the event that triggered it, so a whole cascade shares one id and will
    therefore share one null offset. Background events are each their own sequence.
    """
    n = t.size
    is_after = np.zeros(n, dtype=bool)
    since = np.full(n, np.nan)
    seq = np.arange(n)
    order = np.argsort(t)
    recent = []
    for ii in order:
        while recent and t[ii] - t[recent[0]] > DECLUSTER_DAYS:
            recent.pop(0)
        for jj in reversed(recent):
            if mag[jj] >= mag[ii] and MS.W._great_circle_km(
                    lat[jj], lon[jj], lat[ii], lon[ii]) <= DECLUSTER_KM:
                is_after[ii] = True
                since[ii] = t[ii] - t[jj]
                seq[ii] = seq[jj]
                break
        recent.append(ii)
    return is_after, since, seq


def cluster_sums(F, seq_idx, n_seq, mask):
    """Sum F over the masked events of each sequence -> [n_seq, S].

    Every statistic used here is a SUM over events, so aggregating to sequence level
    first is exact, not an approximation, and it makes the sequence-level null cheap.
    """
    out = np.zeros((n_seq, F.shape[1]))
    np.add.at(out, seq_idx[mask], F[mask])
    return out


def main():
    rng = np.random.default_rng(RNG_SEED)
    t, la, lo, dp, mg, n_holdout, cut = load()
    MS.assert_epoch(t, 2010, "comcat_socal_m25")
    print("SoCal exploration split: %d events (holdout %d reserved, cutoff %s)"
          % (t.size, n_holdout, cut.date()), flush=True)
    is_after, since, seq_raw = classify(t, la, lo, mg)
    uniq, seq_idx = np.unique(seq_raw, return_inverse=True)
    n_seq = uniq.size
    print("  aftershocks %d  background/mainshock %d  SEQUENCES %d"
          % (int(is_after.sum()), int((~is_after).sum()), n_seq), flush=True)

    pops = [
        ("aftershock_lt1d", is_after & (since < 1.0)),
        ("aftershock_1to7d", is_after & (since >= 1.0) & (since < 7.0)),
        ("aftershock_7to30d", is_after & (since >= 7.0)),
        ("mainshock_background", ~is_after),
        ("aftershock_all", is_after),
        ("aftershock_M3plus", is_after & (mg >= 3.0)),
        ("mainshock_M3plus", (~is_after) & (mg >= 3.0)),
    ]
    kept = []
    for nm, m in pops:
        ns = np.unique(seq_idx[m]).size
        if int(m.sum()) >= MIN_SUBPOP and ns >= MIN_SEQUENCES:
            kept.append((nm, m, ns))
            print("    %-24s n=%6d  sequences=%5d" % (nm, int(m.sum()), ns), flush=True)
        else:
            print("    %-24s n=%6d  sequences=%5d   SKIPPED"
                  % (nm, int(m.sum()), ns), flush=True)
    pops = kept

    # ---- ONE offset per SEQUENCE, applied rigidly to every event in it
    off_seq = rng.uniform(-NULL_HALF_WINDOW_DAYS, NULL_HALF_WINDOW_DAYS,
                          size=(n_seq, S_NULL_PER_SEQUENCE))
    off = off_seq[seq_idx]                                   # [n_events, S]
    obs_raw = MS.raw_series(t, la, lo)
    pool = MS.raw_series((t[:, None] + off).ravel(),
                         np.repeat(la, S_NULL_PER_SEQUENCE),
                         np.repeat(lo, S_NULL_PER_SEQUENCE))
    pool = {k: v.reshape(t.size, S_NULL_PER_SEQUENCE) for k, v in pool.items()}
    idx_seq = rng.integers(0, S_NULL_PER_SEQUENCE, size=(n_seq, N_REPLICATES))

    running_max = np.zeros(N_REPLICATES)
    records, zmap = [], {}

    def emit(feature, pop, mask, nev, obs_val, num_seq_sums, denom):
        """num_seq_sums: [n_seq, S] cluster sums of the numerator; denom: event count."""
        g = np.take_along_axis(num_seq_sums, idx_seq, axis=1)      # [n_seq, R]
        null_col = g.sum(axis=0) / denom
        mu, sd = float(null_col.mean()), float(null_col.std(ddof=1))
        if not np.isfinite(sd) or sd <= 0:
            return
        z = (float(obs_val) - mu) / sd
        zn = (null_col - mu) / sd
        np.maximum(running_max, np.abs(zn), out=running_max)
        records.append({"population": pop, "feature": feature, "n": int(nev), "z": z})
        zmap[(pop, feature)] = (z, zn)

    for s in SCALARS:
        thr_hi = float(np.quantile(pool[s], 0.90))
        thr_lo = float(np.quantile(pool[s], 0.10))
        hi_p = (pool[s] > thr_hi).astype(np.float64)
        lo_p = (pool[s] < thr_lo).astype(np.float64)
        o = obs_raw[s]
        for nm, m, _ns in pops:
            d = float(m.sum())
            emit("%s.mean" % s, nm, m, d, o[m].mean(),
                 cluster_sums(pool[s], seq_idx, n_seq, m), d)
            emit("%s.top_decile" % s, nm, m, d, (o[m] > thr_hi).mean(),
                 cluster_sums(hi_p, seq_idx, n_seq, m), d)
            emit("%s.bot_decile" % s, nm, m, d, (o[m] < thr_lo).mean(),
                 cluster_sums(lo_p, seq_idx, n_seq, m), d)
        del hi_p, lo_p

    for a in ANGLES:
        th_p, th_o = pool[a], obs_raw[a]
        for mm in ORDERS:
            cp = np.cos(mm * th_p)
            sp = np.sin(mm * th_p)
            for nm, m, _ns in pops:
                d = float(m.sum())
                gc = np.take_along_axis(cluster_sums(cp, seq_idx, n_seq, m),
                                        idx_seq, axis=1).sum(axis=0) / d
                gs = np.take_along_axis(cluster_sums(sp, seq_idx, n_seq, m),
                                        idx_seq, axis=1).sum(axis=0) / d
                null_col = np.hypot(gc, gs)
                mu, sd = float(null_col.mean()), float(null_col.std(ddof=1))
                if not np.isfinite(sd) or sd <= 0:
                    continue
                obs_val = np.hypot(np.cos(mm * th_o[m]).mean(),
                                   np.sin(mm * th_o[m]).mean())
                z = (float(obs_val) - mu) / sd
                zn = (null_col - mu) / sd
                np.maximum(running_max, np.abs(zn), out=running_max)
                records.append({"population": nm, "feature": "%s.R%d" % (a, mm),
                                "n": int(d), "z": z})
                zmap[(nm, "%s.R%d" % (a, mm))] = (z, zn)
            del cp, sp

    obs_max = max(abs(r["z"]) for r in records)
    fam_p = float((np.sum(running_max >= obs_max) + 1) / (N_REPLICATES + 1))
    crit95 = float(np.quantile(running_max, 0.95))
    records.sort(key=lambda r: -abs(r["z"]))

    diff_running = np.zeros(N_REPLICATES)
    diffs = []
    for (pop, feat), (z, zn) in list(zmap.items()):
        if pop != "aftershock_all":
            continue
        other = zmap.get(("mainshock_background", feat))
        if other is None:
            continue
        dn = zn - other[1]
        sd = float(dn.std(ddof=1))
        if not np.isfinite(sd) or sd <= 0:
            continue
        dz = ((z - other[0]) - float(dn.mean())) / sd
        np.maximum(diff_running, np.abs((dn - dn.mean()) / sd), out=diff_running)
        diffs.append({"feature": feat, "z_aftershock": z, "z_mainshock": other[0],
                      "diff_z": dz})
    diffs.sort(key=lambda r: -abs(r["diff_z"]))
    diff_obs = max(abs(d["diff_z"]) for d in diffs) if diffs else float("nan")
    diff_p = (float((np.sum(diff_running >= diff_obs) + 1) / (N_REPLICATES + 1))
              if diffs else float("nan"))

    bands = ("aftershock_lt1d", "aftershock_1to7d", "aftershock_7to30d")
    grad = {}
    have = [b for b in bands if any(k[0] == b for k in zmap)]
    if len(have) == 3:
        feats = sorted({k[1] for k in zmap if k[0] == have[0]})
        mono = sum(1 for f in feats
                   if all((b, f) in zmap for b in have)
                   and abs(zmap[(have[0], f)][0]) > abs(zmap[(have[1], f)][0])
                   > abs(zmap[(have[2], f)][0]))
        grad = {"features_tested": len(feats), "monotone_declining": mono,
                "expected_by_chance": len(feats) / 6.0,
                "p_binomial_greater": float(sum(
                    math.comb(len(feats), i) * (1 / 6) ** i * (5 / 6) ** (len(feats) - i)
                    for i in range(mono, len(feats) + 1))),
                "caveat": ("the three bands share sequences, so this binomial is "
                           "indicative only and is NOT a calibrated family p")}

    out = {
        "arm": "near-critical: aftershocks as the target, SEQUENCE-level null",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 1, "holdout_touched": False,
        "catalogue": "comcat_socal_m25.csv, M>=2.5, exploration split",
        "n_events": int(t.size), "n_holdout_reserved": int(n_holdout),
        "n_sequences": int(n_seq),
        "n_aftershock": int(is_after.sum()),
        "n_mainshock_background": int((~is_after).sum()),
        "null": {"kind": "waveform-matched dwell-time, ONE OFFSET PER SEQUENCE",
                 "why": ("a per-event null treats co-temporal aftershocks as independent "
                         "samples of a slowly-varying global quantity and inflates z by "
                         "sqrt(events/sequences); the first run of this arm returned "
                         "max |z| = 26.27 at aftershock_lt1d / elongation.R1 for exactly "
                         "that reason and it was an artifact, not a signal"),
                 "half_window_days": NULL_HALF_WINDOW_DAYS,
                 "draws_per_sequence": S_NULL_PER_SEQUENCE,
                 "replicates": N_REPLICATES},
        "populations": [{"name": nm, "n": int(m.sum()), "sequences": int(ns)}
                        for nm, m, ns in pops],
        "n_cells": len(records),
        "cell_max_family": {"observed_max_abs_z": obs_max, "null_max_p95": crit95,
                            "p": fam_p,
                            "where": "%s / %s" % (records[0]["population"],
                                                  records[0]["feature"])},
        "differential_aftershock_minus_mainshock": {
            "observed_max_abs_diff_z": diff_obs,
            "null_max_p95": float(np.quantile(diff_running, 0.95)),
            "p": diff_p, "top": diffs[:10]},
        "gradient_with_time_since_mainshock": grad,
        "top_20_cells": records[:20],
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    print("\n" + "=" * 78)
    print("NEAR-CRITICAL ARM (sequence-level null): %d cells, %d events, %d sequences"
          % (len(records), t.size, n_seq))
    print("  CELL max |z| = %.3f at %s ; null 95th %.3f  ->  FAMILY p = %.4f"
          % (obs_max, out["cell_max_family"]["where"], crit95, fam_p))
    print("  DIFFERENTIAL (aftershock - mainshock) max |z| = %.3f ; null 95th %.3f "
          "->  p = %.4f"
          % (diff_obs, out["differential_aftershock_minus_mainshock"]["null_max_p95"],
             diff_p))
    if grad:
        print("  GRADIENT (indicative): %d of %d features decline monotonically "
              "(chance %.1f)" % (grad["monotone_declining"], grad["features_tested"],
                                 grad["expected_by_chance"]))
    print("\n  top 10 cells:")
    for r in records[:10]:
        print("    %+6.2f  %-24s %-24s n=%d"
              % (r["z"], r["population"], r["feature"], r["n"]))
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
