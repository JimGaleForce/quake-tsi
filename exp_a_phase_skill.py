"""EXP-A (confirmatory headline) - out-of-sample tidal-phase timing skill.

Implements OVERNIGHT_PREDICTION_PROTOCOL.md EXP-A exactly, reusing the validated per-event
focal-mechanism machinery of coso_fm_test.py (load_fm, base_series, combo_coeffs) and
coso_positive_control.py (load_declustered, phase_series, T0/DT/UPSAMPLE) via importlib.

TRAIN (events before 2010-01-01 UTC): per 0.4 deg bin with >= 100 FM-matched train events,
tau (shear) phase per event on its own mechanism (orientations rounded to 15 deg), 36-bin
sinusoid fit -> a_b = Pm/P0 and preferred phase phi_b, empirical p vs 1,000 orientation-
preserving time-randomized synthetics.  SELECTED: p < 0.05.  CONTROL: p > 0.5.

TEST (events from 2010-01-01 UTC, single unblinding): S_b = mean_i cos(phi_i - phi_b) with
phi_b FIXED from train, pooled S = sum n_b S_b / sum n_b; null = 3,000 circular time shifts
(2-2000 d) added to test event times before phase lookup.  Anti-leak control: identical
scoring on the train-null (CONTROL) bins with an independent shift set.

Nothing in the repository is modified by this script; it only writes exp_a_train_bins.csv
and results_exp_a.json.
"""
import json
import time
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, HERE / fname)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cc = _load("cc", "coso_positive_control.py")
fmt = _load("fmt", "coso_fm_test.py")

# ---------------------------------------------------------------- frozen constants
CATALOG = "SCSN_decluster_m1.5.txt"
SPLIT = pd.Timestamp("2010-01-01", tz="UTC")
BIN_DEG = 0.4
LON_ORIGIN, LAT_ORIGIN = -122.0, 31.5
BOX_LAT, BOX_LON = (31.5, 38.0), (-122.0, -113.5)
MIN_TRAIN = 100
COMPONENT = "tau"
ROUND_DEG = fmt.ROUND_DEG          # 15 deg orientation grouping
N_PHASE_BINS = fmt.N_BINS          # 36
N_SYNTH = 1000
SEED_TRAIN = 20260809
N_SHIFT = 3000
SEED_TEST = 20260810
SEED_CONTROL = 20260811
SHIFT_MIN_DAYS, SHIFT_MAX_DAYS = 2.0, 2000.0
A_CLIP = 0.9
RUNTIME_BUDGET_S = 90 * 60

PHASE_EDGES = np.linspace(-180, 180, N_PHASE_BINS + 1)


def fit(h):
    """Sinusoid fit of a 36-bin phase histogram -> (Pm/P0, phi_deg).

    Numerically identical to the `fit` closure inside coso_fm_test.run; it lives inside that
    function so it cannot be imported, hence this minimal re-statement (no other internals of
    coso_fm_test are duplicated - load_fm/base_series/combo_coeffs are imported and used).
    """
    d = h / h.sum() * N_PHASE_BINS
    th = np.radians(PHASE_EDGES[:-1] + 180 / N_PHASE_BINS)
    A = np.c_[np.ones_like(th), np.cos(th), np.sin(th)]
    coef, *_ = np.linalg.lstsq(A, d, rcond=None)
    return np.hypot(coef[1], coef[2]) / coef[0], np.degrees(np.arctan2(coef[2], coef[1]))


def fit_many(rows):
    """Vectorised Pm/P0 for a stack of histograms (synthetics)."""
    d = rows / rows.sum(axis=1, keepdims=True) * N_PHASE_BINS
    th = np.radians(PHASE_EDGES[:-1] + 180 / N_PHASE_BINS)
    A = np.c_[np.ones_like(th), np.cos(th), np.sin(th)]
    coef = np.linalg.lstsq(A, d.T, rcond=None)[0]
    return np.hypot(coef[1], coef[2]) / coef[0]


def main():
    t_start = time.time()
    global N_SYNTH
    synth_reduced = False

    # ---------------------------------------------------------------- catalog + FM
    fm = fmt.load_fm()
    cat = cc.load_declustered(CATALOG)
    m = cat.merge(fm[["eid", "strike", "dip", "rake"]], on="eid", how="inner").reset_index(drop=True)
    print(f"catalog {len(cat)} declustered events, {len(m)} FM-matched")

    m = m[(m.lat >= BOX_LAT[0]) & (m.lat < BOX_LAT[1]) &
          (m.lon >= BOX_LON[0]) & (m.lon < BOX_LON[1])].reset_index(drop=True)
    m["bi"] = np.floor((m.lon - LON_ORIGIN) / BIN_DEG).astype(int)
    m["bj"] = np.floor((m.lat - LAT_ORIGIN) / BIN_DEG).astype(int)
    m["bin_lon0"] = LON_ORIGIN + BIN_DEG * m.bi
    m["bin_lat0"] = LAT_ORIGIN + BIN_DEG * m.bj
    m["is_train"] = m.t_unix.to_numpy() < SPLIT.timestamp()
    print(f"in SoCal box: {len(m)} ({int(m.is_train.sum())} train / {int((~m.is_train).sum())} test)")

    n_train_bin = m[m.is_train].groupby(["bj", "bi"]).size()
    eligible = [k for k, v in n_train_bin.items() if v >= MIN_TRAIN]
    print(f"bins with >= {MIN_TRAIN} FM-matched train events: {len(eligible)} "
          f"(of {len(n_train_bin)} occupied)")

    w = m[m.set_index(["bj", "bi"]).index.isin(eligible)].reset_index(drop=True)
    for c in ["strike", "dip", "rake"]:
        w[f"{c}_r"] = (np.round(w[c] / ROUND_DEG) * ROUND_DEG).astype(int)

    # ---------------------------------------------------------------- tidal series
    SXX, SYY, SXY, SZZ = fmt.base_series()
    nfine = len(SXX)
    span_s = nfine * cc.DT
    wrap_mod = span_s - 86400.0        # protocol: wrap modulo series duration minus 1 day
    print(f"tidal series: {nfine} upsampled samples, dt={cc.DT:.0f}s, span={span_s/86400:.0f} d")

    t_rel = w.t_unix.to_numpy() - cc.T0.timestamp()
    idx_ev = np.clip((t_rel / cc.DT).astype(np.int64), 0, nfine - 1)

    # circular time shifts (test + independent control set), drawn up front
    shifts = np.concatenate([
        np.random.default_rng(SEED_TEST).uniform(SHIFT_MIN_DAYS, SHIFT_MAX_DAYS, N_SHIFT),
        np.random.default_rng(SEED_CONTROL).uniform(SHIFT_MIN_DAYS, SHIFT_MAX_DAYS, N_SHIFT),
    ]) * 86400.0

    is_test = ~w.is_train.to_numpy()
    test_pos = np.flatnonzero(is_test)
    test_row = np.full(len(w), -1, dtype=np.int64)
    test_row[test_pos] = np.arange(len(test_pos))
    shift_phase = np.full((len(test_pos), 2 * N_SHIFT), np.nan, dtype=np.float32)
    print(f"working set: {len(w)} events in eligible bins ({len(test_pos)} test)")

    # ---------------------------------------------------------------- orientation pass
    ev_phase = np.full(len(w), np.nan, dtype=np.float32)
    ev_group = np.full(len(w), -1, dtype=np.int32)
    groups = sorted(w.groupby(["strike_r", "dip_r", "rake_r"]).groups.items(), key=lambda kv: kv[0])
    ngrp = len(groups)
    print(f"orientation groups (15 deg, {COMPONENT}): {ngrp}")
    grp_prob = np.zeros((ngrp, N_PHASE_BINS))

    t_pass = time.time()
    for gi, ((st, di, ra), gidx) in enumerate(groups):
        cn, ct = fmt.combo_coeffs(st, di, ra)
        c = ct if COMPONENT == "tau" else cn
        series = c[0] * SXX + c[1] * SYY + c[2] * SXY + c[3] * SZZ
        phase = cc.phase_series(series)
        gpos = w.index.get_indexer(gidx)
        ev_phase[gpos] = phase[idx_ev[gpos]]
        ev_group[gpos] = gi
        h, _ = np.histogram(phase[~np.isnan(phase)], bins=PHASE_EDGES)
        grp_prob[gi] = h / h.sum()
        gt = gpos[is_test[gpos]]
        if len(gt):
            rel = (t_rel[gt][:, None] + shifts[None, :]) % wrap_mod
            sidx = np.clip((rel / cc.DT).astype(np.int64), 0, nfine - 1)
            shift_phase[test_row[gt]] = phase[sidx]
        if gi == 29:
            proj = (time.time() - t_pass) / 30 * ngrp
            print(f"  [runtime guard] projected orientation pass {proj/60:.1f} min")
            if proj > RUNTIME_BUDGET_S:
                N_SYNTH, synth_reduced = 500, True
                print("  [runtime guard] synthetics reduced to 500")
        if gi % 250 == 0:
            print(f"  group {gi}/{ngrp} ({time.time()-t_pass:.0f}s)")
    t_orient = time.time() - t_pass
    print(f"orientation pass done in {t_orient/60:.1f} min "
          f"(phase NaN rate {float(np.isnan(ev_phase).mean()):.4f})")

    # ---------------------------------------------------------------- TRAIN fits
    rng = np.random.default_rng(SEED_TRAIN)
    tr_mask = w.is_train.to_numpy() & ~np.isnan(ev_phase)
    rows = []
    for (bj, bi) in sorted(eligible):
        sel = np.flatnonzero(tr_mask & (w.bj.to_numpy() == bj) & (w.bi.to_numpy() == bi))
        obs_h, _ = np.histogram(ev_phase[sel], bins=PHASE_EDGES)
        a_b, phi_b = fit(obs_h.astype(float))
        syn_h = np.zeros((N_SYNTH, N_PHASE_BINS))
        gg, counts = np.unique(ev_group[sel], return_counts=True)
        for g, n_g in zip(gg, counts):
            # identical construction to coso_fm_test.run: n_g phases drawn uniformly from the
            # group's valid phase series then binned into 36 bins == multinomial(n_g, p_g)
            syn_h += rng.multinomial(int(n_g), grp_prob[g], size=N_SYNTH)
        syn_a = fit_many(syn_h)
        p_train = float((np.sum(syn_a >= a_b) + 1) / (N_SYNTH + 1))
        rows.append({"bin_lat0": round(LAT_ORIGIN + BIN_DEG * bj, 4),
                     "bin_lon0": round(LON_ORIGIN + BIN_DEG * bi, 4),
                     "n_train": int(len(sel)), "a_b": float(a_b), "phi_b": float(phi_b),
                     "p_train": p_train, "bj": int(bj), "bi": int(bi)})
        print(f"  bin lat0={rows[-1]['bin_lat0']:.1f} lon0={rows[-1]['bin_lon0']:.1f} "
              f"n={len(sel)} a={a_b:.3f} phi={phi_b:7.1f} p={p_train:.4f}"
              f"{'  <-- SELECTED' if p_train < 0.05 else ''}")

    bins_df = pd.DataFrame(rows)
    bins_df.drop(columns=["bj", "bi"]).to_csv(HERE / "exp_a_train_bins.csv", index=False)
    selected = bins_df[bins_df.p_train < 0.05].reset_index(drop=True)
    control = bins_df[bins_df.p_train > 0.5].reset_index(drop=True)
    print(f"\nSELECTED bins (p<0.05): {len(selected)}   CONTROL bins (p>0.5): {len(control)}"
          f"   of {len(bins_df)} eligible")

    # ---------------------------------------------------------------- TEST scoring
    te_mask = is_test & ~np.isnan(ev_phase)
    bj_arr, bi_arr = w.bj.to_numpy(), w.bi.to_numpy()

    def score(table, shift_slice):
        per_bin, num, den = [], 0.0, 0.0
        bits_num, bits_den = 0.0, 0.0
        shift_num = np.zeros(N_SHIFT)
        shift_den = np.zeros(N_SHIFT)
        for r in table.itertuples():
            sel = np.flatnonzero(te_mask & (bj_arr == r.bj) & (bi_arr == r.bi))
            n_b = len(sel)
            entry = {"bin_lat0": r.bin_lat0, "bin_lon0": r.bin_lon0, "n_b": int(n_b),
                     "a_b": float(r.a_b), "phi_b": float(r.phi_b), "p_train": float(r.p_train)}
            if n_b:
                d = np.radians(ev_phase[sel].astype(float) - r.phi_b)
                S_b = float(np.mean(np.cos(d)))
                a_use = float(np.clip(r.a_b, 0.0, A_CLIP))
                bits_b = float(np.mean(np.log2(1.0 + a_use * np.cos(d))))
                entry.update({"S_b": S_b, "bits_per_event": bits_b})
                num += n_b * S_b
                den += n_b
                bits_num += n_b * bits_b
                bits_den += n_b
                sp = shift_phase[test_row[sel]][:, shift_slice].astype(np.float64)
                cs = np.cos(np.radians(sp - r.phi_b))
                ok = ~np.isnan(sp)
                shift_num += np.nansum(np.where(ok, cs, 0.0), axis=0)
                shift_den += ok.sum(axis=0)
            else:
                entry.update({"S_b": None, "bits_per_event": None})
            per_bin.append(entry)
        if den == 0:
            return {"per_bin": per_bin, "n_events": 0, "pooled_S": None, "p": None,
                    "bits_per_event": None}
        S_obs = num / den
        S_shift = shift_num / np.maximum(shift_den, 1)
        p = float((np.sum(S_shift >= S_obs) + 1) / (N_SHIFT + 1))
        return {"per_bin": per_bin, "n_events": int(den), "pooled_S": float(S_obs), "p": p,
                "bits_per_event": float(bits_num / bits_den),
                "null_S_mean": float(S_shift.mean()), "null_S_p95": float(np.percentile(S_shift, 95))}

    print("\n=== TEST (single unblinding) ===")
    test_res = score(selected, slice(0, N_SHIFT))
    print("\n=== ANTI-LEAK CONTROL (train-null bins) ===")
    ctrl_res = score(control, slice(N_SHIFT, 2 * N_SHIFT))

    runtime = time.time() - t_start
    out = {
        "protocol": "OVERNIGHT_PREDICTION_PROTOCOL.md EXP-A",
        "catalog": CATALOG, "component": COMPONENT,
        "split_utc": str(SPLIT), "bin_deg": BIN_DEG,
        "grid_origin": {"lon": LON_ORIGIN, "lat": LAT_ORIGIN},
        "box": {"lat": list(BOX_LAT), "lon": list(BOX_LON)},
        "n_synth_train": N_SYNTH, "synth_reduced": synth_reduced,
        "n_shift": N_SHIFT, "shift_days": [SHIFT_MIN_DAYS, SHIFT_MAX_DAYS],
        "seeds": {"train_synth": SEED_TRAIN, "test_shift": SEED_TEST, "control_shift": SEED_CONTROL},
        "runtime_s": round(runtime, 1),
        "train": {
            "n_fm_matched": int(len(m)), "n_train_events": int(m.is_train.sum()),
            "n_test_events": int((~m.is_train).sum()),
            "n_bins_eligible": int(len(bins_df)),
            "n_selected": int(len(selected)), "n_control": int(len(control)),
            "bins": bins_df.drop(columns=["bj", "bi"]).to_dict(orient="records"),
        },
        "test": {k: v for k, v in test_res.items()},
        "control": {"pooled_S": ctrl_res["pooled_S"], "p": ctrl_res["p"],
                    "n_events": ctrl_res["n_events"],
                    "bits_per_event": ctrl_res["bits_per_event"],
                    "per_bin": ctrl_res["per_bin"]},
    }
    (HERE / "results_exp_a.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

    print("\n================ EXP-A SUMMARY ================")
    print(f"eligible bins (n_train>={MIN_TRAIN}): {len(bins_df)}   selected: {len(selected)}   "
          f"control: {len(control)}")
    for e in test_res["per_bin"]:
        print(f"  SEL lat0={e['bin_lat0']:.1f} lon0={e['bin_lon0']:.1f} a_b={e['a_b']:.3f} "
              f"phi_b={e['phi_b']:7.1f} n_b={e['n_b']:5d} "
              f"S_b={'  n/a' if e['S_b'] is None else format(e['S_b'], '+.4f')}")
    print(f"TEST     pooled S = {test_res['pooled_S']}  p = {test_res['p']}  "
          f"n = {test_res['n_events']}  bits/event = {test_res['bits_per_event']}")
    print(f"CONTROL  pooled S = {ctrl_res['pooled_S']}  p = {ctrl_res['p']}  "
          f"n = {ctrl_res['n_events']}")
    if test_res["p"] is not None and ctrl_res["p"] is not None:
        verdict = "PASS" if (test_res["p"] < 0.05 and ctrl_res["p"] >= 0.05) else "FAIL"
        print(f"frozen success rule (test p<0.05 AND control non-significant): {verdict}")
    print(f"runtime {runtime/60:.1f} min  -> results_exp_a.json, exp_a_train_bins.csv")


if __name__ == "__main__":
    main()
