"""F-012 BLOCKING VERIFICATION - MEASURE kappa, THE CLUSTERING VARIANCE-INFLATION FACTOR.

Faraday's definition (PROGRESS_REGISTER.md F-012 "WHAT IT DOES NOT SHOW" item 2, and
results_k035.json :: limitations[0]):

    "All four target samples come from the DECLUSTERED SCSN M>=1.5 catalogue, so events are
     modelled as independent draws from the thinned ETAS intensity.  Residual Hawkes
     clustering with effective cluster size kappa inflates every MDA by sqrt(kappa)."

WHAT THE K-035 MDA ARITHMETIC ASSUMED (read off exp_k035_power_audit.py, not inferred):
  * injected catalogues:      H[cname] += rng.multinomial(n_g, P)
  * per-bin null (A1, C, D):  syn += nrng.multinomial(n_g, grp_prob[g])   (cached_null_a)
  n_g INDEPENDENT draws over the 36 phase bins per orientation group g.  The variance of the
  resultant is therefore the multinomial (independent-events) variance.
  * pooled null (A2, A3, A3b): the EXP-A circular-shift null, which DOES carry the real
  inter-event time structure -- those three paid part of the penalty inside their threshold.

TARGET QUANTITY.  With u_i = (cos, sin) of the 36-bin-centre tidal phase of event i (exactly
what fit_hist and the pooled S test consume), and (Cbar, Sbar) = mean_i u_i,
    kappa = E[Var(Cbar) + Var(Sbar)] / [Var(Cbar) + Var(Sbar)]_independent
          = 1 + ( sum_{i != j} Cov(u_i, u_j) ) / ( sum_i Var(u_i) )
          = 1 + 2 * sum_{i<j} R_ij / sum_i V_i .
The expectation is over the point process; Cov is taken over the tidal-phase alignment, i.e.
over the EXP-A circular-shift ensemble.  MDA* = MDA * sqrt(kappa).

TWO ESTIMATORS WERE TRIED AND REJECTED THIS SESSION; BOTH ARE RECORDED IN THE OUTPUT.
  v1  kappa = V_shift(real times) / V_indep_analytic.  On a Poisson-time control (truth 1)
      it returned 0.463.  DIAGNOSED, not guessed: over 40 independent Poisson catalogues the
      ratio has mean 1.010 and sd 0.482 at n = 1906 (0.885 +/- 0.406 at n = 200;
      1.102 +/- 0.780 at n = 20000).  The estimator is UNBIASED but carries ~50-80% relative
      scatter at every n, because a rigid shift gives all events ONE shared random variable
      and the tidal phase field's power sits in a few coherent lines, so sum_{i!=j} R(t_i-t_j)
      fluctuates at the same order as the diagonal.  Unusable at the +/-10% target.
  v2  kappa = V_shift(real) / mean_j V_shift(jittered_j).  Same defect: jitter re-randomises
      the numerator's chance term instead of cancelling it.

ESTIMATOR v3 -- TRUNCATED PAIR-CORRELATION WITH A MATCHED-JITTER REFERENCE (reported):
    kappa(T) = 1 + 2 * [ P_real(T) - mean_j P_jit,j(T) ] / sum_i V_i
    P(T)     = sum over ordered pairs with 0 < t_j - t_i <= T of  R_ij
    R_ij     = shift-ensemble sample covariance of u_i and u_j (unbiased, ddof = 1)
    V_i      = shift-ensemble sample variance of u_i (both components summed)
  Truncating at T excludes the far-pair terms, whose expectation is ~0 but whose realisation
  noise is what destroyed v1/v2.  The clustering that survives declustering lives at lags of
  hours to days, so the truncation bias is auditable: kappa is reported over a lag ladder
  T in {0.25 ... 30 d} and must plateau.  PRIMARY T = 10 d, fixed before the run.
  The jitter reference (u_i ~ U(-W/2, W/2), W = 30 d, J realisations) supplies the
  unclustered expectation of P(T) at the same local rate; it is precise because the truncated
  sum has few terms.

BIAS AUDITS (in the same pass, reported before the configurations are read):
  B1  Poisson-time control.  TRUTH kappa = 1.
  B2  Planted overdispersion: every event replicated m = 3 at the identical instant.  Then
      sum_i V_i triples while sum_{i<j} R_ij gains 3n exact-duplicate pairs of covariance V:
      TRUTH kappa_planted = 3 * kappa_base, EXACTLY.
  B3  UNDECLUSTERED branching-ETAS at the frozen EXP-H parameters.  Not a truth case -- the
      upper bracket on what the Xue/Lu declustering removed.

Deterministic: every RNG is np.random.default_rng(<fixed literal>).  No time-based seeding.
Writes results_kappa.json.  Modifies nothing else.
"""
import json
import math
import time
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
T_START = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, HERE / fname)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cc = _load("cc", "coso_positive_control.py")
fmt = _load("fmt", "coso_fm_test.py")
K35 = _load("K35", "exp_k035_power_audit.py")

CATALOG = K35.CATALOG
SPLIT = K35.SPLIT
BIN_DEG = K35.BIN_DEG
LON_ORIGIN, LAT_ORIGIN = K35.LON_ORIGIN, K35.LAT_ORIGIN
BOX_LAT, BOX_LON = K35.BOX_LAT, K35.BOX_LON
MIN_TRAIN = K35.MIN_TRAIN
ROUND_DEG = K35.ROUND_DEG
N_BINS = K35.N_BINS
PHASE_EDGES = K35.PHASE_EDGES
COS_CTR, SIN_CTR = K35.COS_CTR, K35.SIN_CTR
FIG4C_BOX = K35.FIG4C_BOX
ETAS = K35.ETAS
SHIFT_MIN_DAYS, SHIFT_MAX_DAYS = K35.SHIFT_MIN_DAYS, K35.SHIFT_MAX_DAYS

N_SHIFT = 600                  # shift ensemble over which R_ij and V_i are averaged
J_JIT = 6                      # jitter realisations supplying the unclustered reference
W_JIT_D = 30.0                 # jitter window, days (>> the clustering timescale)
T_LADDER_D = [0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0]
T_PRIMARY_D = 10.0             # PRE-DECLARED primary truncation
SEED_SHIFT = 20260811
SEED_JITTER = 20260815
SEED_POISSON = 20260812
SEED_ETAS = 20260813
SEED_GRPMIX = 20260814
PLANT_M = 3
N_ETAS_REAL = 1
J_ETAS = 4
PAIR_CHUNK = 8_000

V1_REJECTED = {"estimator": "kappa = V_shift(real)/V_indep_analytic",
               "B1_poisson_control_kappa": 0.4630345839207085, "truth": 1.0,
               "diagnosis_40_poisson_catalogs": {
                   "n200": {"mean": 0.885, "sd": 0.406},
                   "n1906": {"mean": 1.010, "sd": 0.482},
                   "n20000": {"mean": 1.102, "sd": 0.780}},
               "verdict": "REJECTED: unbiased but ~50-80% relative scatter at every n. The "
                          "rigid-shift variance of ONE catalogue is not an estimator of the "
                          "unconditional variance. Not used for any reported number."}


def main():
    import os
    global N_SHIFT, J_JIT, J_ETAS, N_ETAS_REAL, MIN_TRAIN
    SMOKE = os.environ.get("KAPPA_SMOKE") == "1"   # development path only; never recorded
    if SMOKE:
        N_SHIFT, J_JIT, J_ETAS, N_ETAS_REAL, MIN_TRAIN = 120, 3, 2, 1, 15
    out = {"experiment": "F-012 kappa (clustering variance inflation of the K-035 MDA bounds)",
           "state_class": "first-run for the reported estimator (v3). Two earlier estimators "
                          "(v1, v2) were run and rejected on their bias audits in the same "
                          "session; both are recorded below.",
           "run_utc": pd.Timestamp.now("UTC").isoformat(),
           "prices": "results_k035.json :: corpse_to_bound_table_K032_item6",
           "estimator": {
               "reported": "v3 truncated pair-correlation with matched-jitter reference",
               "target": "kappa = 1 + 2*sum_{i<j} Cov(u_i,u_j) / sum_i Var(u_i), u = "
                         "(cos,sin) of the 36-bin-centre tidal phase; expectation over the "
                         "EXP-A circular-shift ensemble.",
               "formula": "kappa(T) = 1 + 2*[P_real(T) - mean_j P_jit,j(T)] / sum_i V_i",
               "P": "sum of R_ij over ordered pairs with 0 < t_j - t_i <= T",
               "R_ij": "unbiased shift-ensemble sample covariance of u_i and u_j (ddof=1)",
               "reference": "J = %d jitter realisations, u_i ~ U(-%g/2, +%g/2) d" % (
                   J_JIT, W_JIT_D, W_JIT_D),
               "n_shift": N_SHIFT, "T_primary_days": T_PRIMARY_D,
               "T_ladder_days": T_LADDER_D,
               "correction": "MDA* = MDA * sqrt(kappa)"},
           "rejected_estimators": {"v1": V1_REJECTED,
                                   "v2": {"estimator": "V_shift(real)/mean_j V_shift(jitter_j)",
                                          "verdict": "REJECTED: jitter re-randomises the "
                                                     "numerator's chance term rather than "
                                                     "cancelling it; same 50-80% scatter."}},
           "seeds": {"shift": SEED_SHIFT, "jitter": SEED_JITTER, "poisson_control": SEED_POISSON,
                     "etas_bracket": SEED_ETAS, "group_mix": SEED_GRPMIX,
                     "policy": "fixed integer literals only; no time-based seeding"},
           "multiplicity": "ONE pre-declared measurement (F-012's single missing verification), "
                           "not a scan. Configurations are exactly the six rows of the K-035 "
                           "bounds table; the primary truncation T = 10 d and the jitter window "
                           "W = 30 d were fixed before the run, and the full T ladder is "
                           "reported rather than selected from. No p-value is reported."}

    # ------------------------------------------------------------ catalogue (K-035 identical)
    fm = fmt.load_fm()
    cat = cc.load_declustered(CATALOG)
    m = cat.merge(fm[["eid", "strike", "dip", "rake"]], on="eid", how="inner").reset_index(drop=True)
    m = m[(m.lat >= BOX_LAT[0]) & (m.lat < BOX_LAT[1]) &
          (m.lon >= BOX_LON[0]) & (m.lon < BOX_LON[1])].reset_index(drop=True)
    if SMOKE:
        m = m.iloc[::17].reset_index(drop=True)
    m["bi"] = np.floor((m.lon - LON_ORIGIN) / BIN_DEG).astype(int)
    m["bj"] = np.floor((m.lat - LAT_ORIGIN) / BIN_DEG).astype(int)
    m["is_train"] = m.t_unix.to_numpy() < SPLIT.timestamp()
    for c in ["strike", "dip", "rake"]:
        m[f"{c}_r"] = (np.round(m[c] / ROUND_DEG) * ROUND_DEG).astype(int)
    n_train_bin = m[m.is_train].groupby(["bj", "bi"]).size()
    eligible = sorted([k for k, v in n_train_bin.items() if v >= MIN_TRAIN])
    ebidx = {k: i for i, k in enumerate(eligible)}
    bin_of = np.array([ebidx.get(k, -1) for k in zip(m.bj, m.bi)], dtype=np.int64)
    in_fig4c = ((m.lat.between(*FIG4C_BOX["lat"])) & (m.lon.between(*FIG4C_BOX["lon"]))).to_numpy()
    gkeys, gcode = np.unique(m[["strike_r", "dip_r", "rake_r"]].to_numpy(), axis=0,
                             return_inverse=True)
    NG = len(gkeys)
    tr = m.is_train.to_numpy()
    elig = bin_of >= 0
    sel_key = (int(round((33.5 - LAT_ORIGIN) / BIN_DEG)),
               int(round((-116.4 - LON_ORIGIN) / BIN_DEG)))
    SEL_B = ebidx[sel_key] if sel_key in ebidx else int(np.argmax(np.bincount(bin_of[bin_of >= 0])))
    ctrl_bins = []
    tb = pd.read_csv(HERE / "exp_a_train_bins.csv")
    for r in tb[tb.p_train > 0.5].itertuples():
        k = (int(round((r.bin_lat0 - LAT_ORIGIN) / BIN_DEG)),
             int(round((r.bin_lon0 - LON_ORIGIN) / BIN_DEG)))
        if k in ebidx:
            ctrl_bins.append(ebidx[k])
    ctrl_bins = np.array(sorted(ctrl_bins), dtype=np.int64)
    ctrl_mask = np.isin(bin_of, ctrl_bins) if len(ctrl_bins) else np.zeros(len(m), bool)

    masks = {"A1_train_bin245": tr & (bin_of == SEL_B),
             "A2_end_to_end": (~tr) & elig,
             "A3_pooled_test_all_eligible": (~tr) & elig,
             "A3b_pooled_control_bins": (~tr) & ctrl_mask,
             "C_coso_fig4c": in_fig4c,
             "D_full_catalog": np.ones(len(m), bool)}
    print(f"FM-matched in box: {len(m)}; eligible bins {len(eligible)}; groups {NG}; "
          f"control bins {len(ctrl_bins)}")
    for k, v in masks.items():
        print(f"  {k}: n = {int(v.sum())}")
    out["catalog"] = {"file": CATALOG, "n_fm_matched_in_box": int(len(m)),
                      "n_orientation_groups": int(NG), "n_eligible_bins": len(eligible),
                      "n_control_bins": int(len(ctrl_bins)),
                      "config_n": {k: int(v.sum()) for k, v in masks.items()}}

    SXX, SYY, SXY, SZZ = fmt.base_series()
    nfine = len(SXX)
    span_s = nfine * cc.DT
    wrap_mod = span_s - 86400.0
    t_rel_ev = m.t_unix.to_numpy() - cc.T0.timestamp()
    t_lo, t_hi = float(t_rel_ev.min()), float(t_rel_ev.max())

    # ------------------------------------------------------------ event sets
    # jitter of a subset == the subset of a jittered full set (jitter is i.i.d. per event),
    # so the six configurations share ONE real set and ONE jitter ensemble.
    jrng = np.random.default_rng(SEED_JITTER)
    prng = np.random.default_rng(SEED_POISSON)
    mixrng = np.random.default_rng(SEED_GRPMIX)

    def jit(t, J):
        return [t + jrng.uniform(-W_JIT_D / 2, W_JIT_D / 2, len(t)) * 86400.0 for _ in range(J)]

    universes = {}       # uname -> dict(t=..., g=..., jits=[...])
    universes["MAIN"] = {"t": t_rel_ev, "g": gcode, "jits": jit(t_rel_ev, J_JIT)}

    # Poisson controls are DENSITY-MATCHED to the configurations they audit: same n, same
    # time window, same orientation mix. B1a matches D (largest n -> tightest audit);
    # B1b matches A3b (the pooled corpse, on its own 2010-2018 test window).
    t_a = t_rel_ev[masks["A3b_pooled_control_bins"]]
    g_a = gcode[masks["A3b_pooled_control_bins"]]
    t_pois = np.sort(prng.uniform(t_lo, t_hi, len(t_rel_ev)))
    universes["B1a_poisson_D"] = {"t": t_pois, "g": gcode.copy(),
                                  "jits": jit(t_pois, J_JIT)}
    t_pois_b = np.sort(prng.uniform(t_a.min(), t_a.max(), len(t_a)))
    universes["B1b_poisson_A3b"] = {"t": t_pois_b, "g": g_a.copy(),
                                    "jits": jit(t_pois_b, J_JIT)}
    t_pl = np.concatenate([t_a] * PLANT_M)
    g_pl = np.tile(g_a, PLANT_M)
    o = np.argsort(t_pl, kind="stable")
    t_pl, g_pl = t_pl[o], g_pl[o]
    universes["B2_plant_m3"] = {"t": t_pl, "g": g_pl, "jits": jit(t_pl, J_JIT)}

    span_days = (m.t_unix.max() - m.t_unix.min()) / 86400.0
    mags = cat.mag.to_numpy()
    bval = float(1.0 / (math.log(10.0) * (mags[mags >= 2.5].mean() - 2.45)))
    etas_meta = []
    for r in range(N_ETAS_REAL):
        tt = K35.simulate_etas(span_days, ETAS["mu"], ETAS["K"], ETAS["alpha"], ETAS["c"],
                               ETAS["p"], ETAS["M0"], ETAS["tmax"], bval,
                               np.random.default_rng(SEED_ETAS + r))
        t_sim = (m.t_unix.min() - cc.T0.timestamp()) + tt * 86400.0
        g_sim = mixrng.choice(g_a, len(t_sim), replace=True)
        universes[f"B3_etas_r{r}"] = {"t": t_sim, "g": g_sim, "jits": jit(t_sim, J_ETAS)}
        etas_meta.append({"realization": r, "n_events": int(len(tt))})
    print(f"B3 undeclustered ETAS: {[e['n_events'] for e in etas_meta]} events")

    shifts = np.random.default_rng(SEED_SHIFT).uniform(
        SHIFT_MIN_DAYS, SHIFT_MAX_DAYS, N_SHIFT) * 86400.0

    # flat list of (universe, replicate index) -> arrays to fill with binned phase indices
    reps = []
    for un, u in universes.items():
        reps.append((un, -1, u["t"], u["g"]))
        for j, tj in enumerate(u["jits"]):
            reps.append((un, j, tj, u["g"]))
    store = {}
    for un, j, t, g in reps:
        store[(un, j)] = np.full((len(t), N_SHIFT), -1, np.int8)
    tot = sum(v.size for v in store.values())
    print(f"replicates: {len(reps)}; phase-index store: {tot/1e6:.0f} M cells "
          f"({tot/1e6:.0f} MB)")

    # ------------------------------------------------------------ single group pass
    idx_by = {(un, j): {int(gg): np.flatnonzero(g == gg) for gg in np.unique(g)}
              for un, j, t, g in reps}
    t_by = {(un, j): t for un, j, t, g in reps}
    groups_needed = sorted({int(gg) for un, j, t, g in reps for gg in np.unique(g)})
    step = max(1, 4_000_000 // N_SHIFT)
    t_gp = time.time()
    for jj, gi in enumerate(groups_needed):
        st, di, ra = gkeys[gi]
        cn, ct = fmt.combo_coeffs(st, di, ra)
        ser = ct[0] * SXX + ct[1] * SYY + ct[2] * SXY + ct[3] * SZZ
        ph = K35.phase_series_fast(ser)
        good_s = ~np.isnan(ph)
        b = np.full(nfine, -1, np.int8)
        b[good_s] = np.clip(((ph[good_s] + 180.0) / (360.0 / N_BINS)).astype(np.int64),
                            0, N_BINS - 1).astype(np.int8)
        for key in store:
            ii = idx_by[key].get(gi)
            if ii is None or not len(ii):
                continue
            tt_ = t_by[key]
            arr = store[key]
            for a0 in range(0, len(ii), step):
                sub = ii[a0:a0 + step]
                rel = (tt_[sub][:, None] + shifts[None, :]) % wrap_mod
                arr[sub] = b[np.clip((rel / cc.DT).astype(np.int64), 0, nfine - 1)]
        if jj % 300 == 0:
            print(f"  group {jj}/{len(groups_needed)}  ({time.time()-t_gp:.0f}s, "
                  f"total {(time.time()-T_START)/60:.1f} min)", flush=True)
    print(f"group pass done in {(time.time()-t_gp)/60:.1f} min")

    # ------------------------------------------------------------ pair machinery
    T_MAX_S = max(T_LADDER_D) * 86400.0

    def uvec(key, sel, order):
        b = store[key][sel][order]
        ok = b >= 0
        bi = np.where(ok, b, 0)
        cu = np.where(ok, COS_CTR[bi], 0.0).astype(np.float32)
        su = np.where(ok, SIN_CTR[bi], 0.0).astype(np.float32)
        return cu, su

    def var_sum(cu, su):
        """sum_i [Var_s(cos_i) + Var_s(sin_i)] over the shift ensemble (ddof=1)."""
        return float((cu.var(axis=1, ddof=1) + su.var(axis=1, ddof=1)).sum())

    def pair_profile(key, sel):
        """Cumulative sum of R_ij over pairs with lag <= each T in the ladder, for events
        selected by boolean mask `sel` within the replicate's own event list."""
        t = t_by[key][sel]
        order = np.argsort(t, kind="stable")
        t = t[order]
        cu, su = uvec(key, sel, order)
        n = len(t)
        hi = np.searchsorted(t, t + T_MAX_S, side="right")
        cnt = hi - np.arange(n) - 1
        tot = int(cnt.sum())
        prof = np.zeros(len(T_LADDER_D))
        if tot == 0:
            return prof, 0, var_sum(cu, su), n
        ii = np.repeat(np.arange(n), cnt)
        off = np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        jj = ii + 1 + off
        lag = t[jj] - t[ii]
        S = N_SHIFT
        for a0 in range(0, tot, PAIR_CHUNK):
            sl = slice(a0, a0 + PAIR_CHUNK)
            a, bidx = ii[sl], jj[sl]
            ca, cb = cu[a], cu[bidx]
            sa, sb = su[a], su[bidx]
            r = ((ca * cb).sum(1) - ca.sum(1) * cb.sum(1) / S +
                 (sa * sb).sum(1) - sa.sum(1) * sb.sum(1) / S) / (S - 1)
            lg = lag[sl]
            for ti, T in enumerate(T_LADDER_D):
                prof[ti] += float(r[lg <= T * 86400.0].sum())
        return prof, tot, var_sum(cu, su), n

    def kappa_ladder(uname, sel, J):
        key0 = (uname, -1)
        p_real, np_real, vsum, n = pair_profile(key0, sel)
        pj = np.array([pair_profile((uname, j), sel)[0] for j in range(J)])
        ref = pj.mean(axis=0)
        sd = pj.std(axis=0, ddof=1)
        kap = 1.0 + 2.0 * (p_real - ref) / vsum
        se = 2.0 * sd * math.sqrt(1.0 + 1.0 / J) / vsum
        return {"n": int(n), "n_pairs_le_Tmax": int(np_real), "sum_var": vsum,
                "T_ladder_days": T_LADDER_D,
                "kappa_ladder": [float(x) for x in kap],
                "kappa_ladder_se": [float(x) for x in se],
                "kappa": float(kap[T_LADDER_D.index(T_PRIMARY_D)]),
                "kappa_se": float(se[T_LADDER_D.index(T_PRIMARY_D)]),
                "ci95": [float(kap[T_LADDER_D.index(T_PRIMARY_D)]
                               - 1.96 * se[T_LADDER_D.index(T_PRIMARY_D)]),
                         float(kap[T_LADDER_D.index(T_PRIMARY_D)]
                               + 1.96 * se[T_LADDER_D.index(T_PRIMARY_D)])]}

    def finish(d):
        k = d["kappa"]
        lo, hi = d["ci95"]
        d["sqrt_kappa"] = float(math.sqrt(max(k, 1e-9)))
        d["sqrt_ci95"] = [float(math.sqrt(max(lo, 1e-9))), float(math.sqrt(max(hi, 1e-9)))]
        d["n_eff"] = float(d["n"] / k) if k > 0 else None
        return d

    # ------------------------------------------------------------ BIAS AUDITS
    print("\n=== bias audits ===")
    audits = {}
    base = finish(kappa_ladder("MAIN", masks["A3b_pooled_control_bins"], J_JIT))
    audits["base_configuration"] = "A3b_pooled_control_bins"
    audits["kappa_base"] = base["kappa"]

    b1all = {}
    for tag, un, nn in (("B1a_matched_to_D", "B1a_poisson_D", len(t_pois)),
                        ("B1b_matched_to_A3b", "B1b_poisson_A3b", len(t_pois_b))):
        d = finish(kappa_ladder(un, np.ones(nn, bool), J_JIT))
        d["truth_kappa"] = 1.0
        d["bias"] = float(d["kappa"] - 1.0)
        # pre-declared: consistent with the truth within 2 SE, and no worse than 5% off
        d["pass"] = bool(abs(d["kappa"] - 1.0) <= max(0.05, 2.0 * d["kappa_se"]))
        b1all[tag] = d
        print(f"[{tag}] Poisson n = {d['n']}: kappa = {d['kappa']:.4f} "
              f"+/- {d['kappa_se']:.4f}  truth 1.0  bias {d['bias']:+.4f}  PASS={d['pass']}")
        print(f"      ladder {[round(x,4) for x in d['kappa_ladder']]}")
    audits["B1_poisson_controls"] = b1all
    b1 = b1all["B1a_matched_to_D"]

    b2 = finish(kappa_ladder("B2_plant_m3", np.ones(len(t_pl), bool), J_JIT))
    truth = PLANT_M * base["kappa"]
    b2["truth_kappa"] = float(truth)
    b2["ratio_to_truth"] = float(b2["kappa"] / truth)
    b2["pass"] = bool(abs(b2["kappa"] / truth - 1.0) <= 0.05)
    audits["B2_planted_replication_m3"] = b2
    print(f"[B2] planted m=3: kappa = {b2['kappa']:.4f} +/- {b2['kappa_se']:.4f}  "
          f"truth {truth:.4f}  ratio {b2['ratio_to_truth']:.4f}  PASS={b2['pass']}")

    b3 = {}
    for r in range(N_ETAS_REAL):
        un = f"B3_etas_r{r}"
        d = finish(kappa_ladder(un, np.ones(len(universes[un]["t"]), bool), J_ETAS))
        b3[f"r{r}"] = d
        print(f"[B3] undeclustered ETAS r{r}: n = {d['n']}, kappa = {d['kappa']:.3f} "
              f"+/- {d['kappa_se']:.3f}")
    audits["B3_etas_undeclustered_bracket"] = {
        "realizations": etas_meta, "per_realization": b3,
        "note": "NOT a truth case: UNDECLUSTERED branching-ETAS at the frozen EXP-H parameters "
                "over the same span, full simulated catalogue, orientation mix resampled from "
                "the A3b sample. Upper bracket on what the Xue/Lu declustering removed."}
    out["bias_audits"] = audits
    out["audit_verdict"] = ("PASS" if (all(d["pass"] for d in b1all.values()) and b2["pass"])
                            else "FAIL")
    out["audit_pass_rule"] = ("pre-declared: |kappa_Poisson - 1| <= max(0.05, 2 SE) on BOTH "
                              "density-matched Poisson controls, and "
                              "|kappa_planted/(3*kappa_base) - 1| <= 0.05")
    print(f"AUDIT VERDICT: {out['audit_verdict']}")

    # ------------------------------------------------------------ kappa per configuration
    print("\n=== kappa per K-035 configuration (T = 10 d primary) ===")
    per_cfg = {}
    for cname in [c for c in masks if c != "A2_end_to_end"] + ["A2_end_to_end"]:
        if cname == "A2_end_to_end":
            kk = dict(per_cfg["A3_pooled_test_all_eligible"])
            kk["note"] = "identical test-event set to A3_pooled_test_all_eligible"
        elif cname == "A3b_pooled_control_bins":
            kk = dict(base)          # already computed for the bias-audit baseline
        else:
            kk = finish(kappa_ladder("MAIN", masks[cname], J_JIT))
        per_cfg[cname] = kk
        print(f"  {cname}: n = {kk['n']:6d}  pairs<=30d {kk['n_pairs_le_Tmax']:8d}  "
              f"kappa = {kk['kappa']:.3f} +/- {kk['kappa_se']:.3f}  "
              f"sqrt(kappa) = {kk['sqrt_kappa']:.3f}  n_eff = {kk['n_eff']:.0f}")
        print(f"      ladder {[round(x,3) for x in kk['kappa_ladder']]}")
    out["kappa_by_configuration"] = per_cfg

    # ------------------------------------------------------------ re-quoted bounds
    k35 = json.loads((HERE / "results_k035.json").read_text())
    row_cfg = {
        "EXP-A per-bin train statistic (the selected bin, n=245)": "A1_train_bin245",
        "EXP-A end-to-end pipeline (1-of-42 selection + pooled circular-shift test)":
            "A2_end_to_end",
        "EXP-A pooled statistic on the 22 train-null CONTROL bins (the quoted corpse: "
        "n=1906, S=0.009, p=0.28)": "A3b_pooled_control_bins",
        "EXP-A pooled statistic over all 42 eligible bins (upper limit of the design)":
            "A3_pooled_test_all_eligible",
        "Coso Fig 4c replication (n=113)": "C_coso_fig4c",
        "Full-catalogue intensity likelihood (n=23465)": "D_full_catalog"}
    shift_null_cfg = {"A2_end_to_end", "A3_pooled_test_all_eligible", "A3b_pooled_control_bins"}
    z95, z80 = 1.6448536269514722, 0.8416212335729143

    table = []
    for row in k35["corpse_to_bound_table_K032_item6"]:
        cfg = row_cfg[row["corpse"]]
        kk = per_cfg[cfg]
        sk = kk["sqrt_kappa"]
        mda = row["mda80"]
        new = mda * sk
        new_lo, new_hi = mda * kk["sqrt_ci95"][0], mda * kk["sqrt_ci95"][1]
        f_ta = sk * (z95 + z80) / (z95 * sk + z80) if cfg in shift_null_cfg else sk
        table.append({
            "corpse": row["corpse"], "config": cfg, "n": row["n"],
            "kappa": kk["kappa"], "kappa_ci95": kk["ci95"], "sqrt_kappa": sk,
            "mda80_old": mda, "mda80_old_pct": 100 * mda,
            "mda80_kappa_corrected": new, "mda80_kappa_corrected_pct": 100 * new,
            "mda80_kappa_corrected_pct_ci95": [100 * new_lo, 100 * new_hi],
            "null_carried_clustering": cfg in shift_null_cfg,
            "threshold_aware_factor": float(f_ta),
            "mda80_threshold_aware_pct": float(100 * mda * f_ta),
            "mda_over_theory_old": row["mda_over_theory"],
            "mda_over_theory_new": float(100 * new),
            "contacts_theory_after_kappa": bool(100 * new <= 1.0)})
        print(f"  {cfg}: MDA {100*mda:.2f}% -> {100*new:.2f}% "
              f"[{100*new_lo:.2f}, {100*new_hi:.2f}]  "
              f"(threshold-aware {100*mda*f_ta:.2f}%)")
    out["requoted_bounds"] = table
    pooled = [r for r in table if r["config"] == "A3_pooled_test_all_eligible"][0]
    out["headline"] = {
        "pooled_all_eligible_old_pct": pooled["mda80_old_pct"],
        "pooled_all_eligible_new_pct": pooled["mda80_kappa_corrected_pct"],
        "pooled_all_eligible_new_pct_ci95": pooled["mda80_kappa_corrected_pct_ci95"],
        "pooled_stays_single_digit": bool(pooled["mda80_kappa_corrected_pct"] < 10.0),
        "max_kappa_over_configs": float(max(r["kappa"] for r in table)),
        "any_bound_contacts_theory_after_kappa":
            bool(any(r["contacts_theory_after_kappa"] for r in table))}
    out["runtime_minutes"] = (time.time() - T_START) / 60.0
    (HERE / "results_kappa.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"\n-> results_kappa.json  ({out['runtime_minutes']:.1f} min)")


if __name__ == "__main__":
    main()
