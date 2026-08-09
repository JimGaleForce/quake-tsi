"""EXP-B (exploratory, split-clean) - spatial transfer of susceptibility.

Implements OVERNIGHT_PREDICTION_PROTOCOL.md EXP-B: for each of the 42 bins in
exp_a_train_bins.csv (the label: train-period modulation amplitude a_b and its
significance p_train), compute five TRAIN-period-only features -- event rate,
Aki (1965) MLE b-value, median depth, swarm fraction, and max shear strain rate
at the bin center -- and report their rank correlation with the label plus a
leave-one-longitude-stripe-out logistic AUC for the (degenerate, single-positive)
"significant" label.

Nothing in the repository is modified by this script; it only writes
results_exp_b.json.
"""
import json
import time
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).parent


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, HERE / fname)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


cc = _load("cc", "coso_positive_control.py")

# ---------------------------------------------------------------- frozen constants
CATALOG = "SCSN_decluster_m1.5.txt"
SPLIT = pd.Timestamp("2010-01-01", tz="UTC")
TRAIN_YEARS = 29.0          # 1981-01-01 -> 2010-01-01
BIN_DEG = 0.4
LON_ORIGIN, LAT_ORIGIN = -122.0, 31.5
MC = 1.5                    # completeness magnitude for b-value
HALF_BIN = 0.05             # Aki (1965) half-bin correction
MIN_N_BVAL = 50
SWARM_KM = 2.0
SWARM_DAYS = 3.0
STRAIN_GRID = "data/socal_strain_grid.npz"
BINS_CSV = "exp_a_train_bins.csv"
SEED = 20260809


def km_projection(lat_deg, lon_deg, lat0):
    """Equirectangular projection to km, centered on lat0 (fine at 0.4 deg scale)."""
    km_per_deg_lat = 111.32
    km_per_deg_lon = 111.32 * np.cos(np.radians(lat0))
    return lon_deg * km_per_deg_lon, lat_deg * km_per_deg_lat


def swarm_fraction(t_unix, lat, lon):
    """Fraction of events with >=1 OTHER event within 2 km epicentral distance
    AND 3 days (either direction). Sort-by-time windowing: O(n log n + n*w)."""
    n = len(t_unix)
    if n == 0:
        return None
    if n == 1:
        return 0.0
    order = np.argsort(t_unix)
    t_s = t_unix[order]
    lat0 = float(np.mean(lat))
    x, y = km_projection(lat[order], lon[order], lat0)
    window_s = SWARM_DAYS * 86400.0
    has_neighbor = np.zeros(n, dtype=bool)
    lo = 0
    hi = 0
    for i in range(n):
        while lo < n and t_s[lo] < t_s[i] - window_s:
            lo += 1
        if hi < i:
            hi = i
        while hi < n and t_s[hi] <= t_s[i] + window_s:
            hi += 1
        cand = np.arange(lo, hi)
        cand = cand[cand != i]
        if len(cand):
            d = np.hypot(x[cand] - x[i], y[cand] - y[i])
            if np.any(d <= SWARM_KM):
                has_neighbor[i] = True
    return float(has_neighbor.mean())


def b_value_mle(mags):
    """Aki (1965) MLE b-value with half-bin correction, M >= MC events only."""
    m = mags[mags >= MC]
    n = len(m)
    if n < MIN_N_BVAL:
        return None, n
    mbar = float(np.mean(m))
    denom = mbar - (MC - HALF_BIN)
    if denom <= 0:
        return None, n
    b = np.log10(np.e) / denom
    return float(b), n


def nearest_strain(lat_c, lon_c, lats, lons, grid):
    j = int(np.argmin(np.abs(lats - lat_c)))
    i = int(np.argmin(np.abs(lons - lon_c)))
    v = grid[j, i]
    return None if np.isnan(v) else float(v)


def main():
    t_start = time.time()

    bins_df = pd.read_csv(HERE / BINS_CSV)
    print(f"loaded {len(bins_df)} bins from {BINS_CSV}")

    cat = cc.load_declustered(CATALOG)
    train = cat[cat.t_unix < SPLIT.timestamp()].reset_index(drop=True)
    print(f"catalog {len(cat)} events, {len(train)} in TRAIN period (< {SPLIT})")

    strain = np.load(HERE / STRAIN_GRID)
    s_lats, s_lons, s_grid = strain["lats"], strain["lons"], strain["max_shear_nstrain_yr"]
    print(f"strain grid: {s_grid.shape} nodes, "
          f"lat[{s_lats.min()},{s_lats.max()}] lon[{s_lons.min()},{s_lons.max()}]")

    t_lat = train.lat.to_numpy()
    t_lon = train.lon.to_numpy()
    t_mag = train.mag.to_numpy()
    t_depth = train.depth.to_numpy()
    t_unix = train.t_unix.to_numpy()

    rows = []
    for r in bins_df.itertuples():
        lat0, lon0 = float(r.bin_lat0), float(r.bin_lon0)
        m = ((t_lat >= lat0) & (t_lat < lat0 + BIN_DEG) &
             (t_lon >= lon0) & (t_lon < lon0 + BIN_DEG))
        n_bin = int(m.sum())
        rate = n_bin / TRAIN_YEARS
        b_val, n_bval = b_value_mle(t_mag[m])
        med_depth = float(np.median(t_depth[m])) if n_bin else None
        swarm_frac = swarm_fraction(t_unix[m], t_lat[m], t_lon[m])
        lat_c, lon_c = lat0 + BIN_DEG / 2, lon0 + BIN_DEG / 2
        strain_val = nearest_strain(lat_c, lon_c, s_lats, s_lons, s_grid)

        rows.append({
            "bin_lat0": lat0, "bin_lon0": lon0,
            "n_train_all": n_bin, "n_train_fm": int(r.n_train),
            "a_b": float(r.a_b), "phi_b": float(r.phi_b), "p_train": float(r.p_train),
            "significant": bool(r.p_train < 0.05),
            "rate_per_yr": rate, "b_value": b_val, "n_for_bvalue": n_bval,
            "median_depth_km": med_depth, "swarm_fraction": swarm_frac,
            "max_shear_nstrain_yr": strain_val,
        })
        print(f"  lat0={lat0:.1f} lon0={lon0:.1f} n={n_bin:5d} rate={rate:6.2f}/yr "
              f"b={('n/a' if b_val is None else format(b_val, '.3f'))} "
              f"depth={med_depth:.1f}km swarm={swarm_frac:.3f} "
              f"strain={('n/a' if strain_val is None else format(strain_val, '.1f'))}"
              f"{'  <-- SIGNIFICANT' if r.p_train < 0.05 else ''}")

    feat_df = pd.DataFrame(rows)
    print(f"\nfeature table built for {len(feat_df)} bins in {time.time()-t_start:.1f}s")

    # ---------------------------------------------------------------- Spearman correlations
    FEATURES = ["rate_per_yr", "b_value", "median_depth_km", "swarm_fraction", "max_shear_nstrain_yr"]

    def spearman_block(target_col, invert=False):
        out = {}
        y_full = feat_df[target_col].to_numpy(dtype=float)
        if invert:
            y_full = -y_full
        for f in FEATURES:
            x = feat_df[f].to_numpy(dtype=float)
            mask = ~np.isnan(x) & ~np.isnan(y_full)
            n = int(mask.sum())
            if n < 3:
                out[f] = {"rho": None, "p": None, "n": n}
                continue
            rho, p = spearmanr(x[mask], y_full[mask])
            out[f] = {"rho": float(rho), "p": float(p), "n": n}
        return out

    spearman_vs_a_b = spearman_block("a_b")
    spearman_vs_p_train = spearman_block("p_train", invert=True)

    print("\n=== Spearman rank correlation: feature vs a_b (label) ===")
    print(f"{'feature':22s} {'rho':>8s} {'p':>10s} {'n':>4s}")
    for f in FEATURES:
        d = spearman_vs_a_b[f]
        rho_s = "n/a" if d["rho"] is None else f"{d['rho']:+.3f}"
        p_s = "n/a" if d["p"] is None else f"{d['p']:.4f}"
        print(f"{f:22s} {rho_s:>8s} {p_s:>10s} {d['n']:4d}")

    print("\n=== Spearman rank correlation: feature vs -p_train (robustness cut; "
          "positive rho = feature co-varies with significance) ===")
    print(f"{'feature':22s} {'rho':>8s} {'p':>10s} {'n':>4s}")
    for f in FEATURES:
        d = spearman_vs_p_train[f]
        rho_s = "n/a" if d["rho"] is None else f"{d['rho']:+.3f}"
        p_s = "n/a" if d["p"] is None else f"{d['p']:.4f}"
        print(f"{f:22s} {rho_s:>8s} {p_s:>10s} {d['n']:4d}")

    # ---------------------------------------------------------------- pairwise feature intercorrelations
    intercorr = {}
    for i, fi in enumerate(FEATURES):
        for fj in FEATURES[i + 1:]:
            xi = feat_df[fi].to_numpy(dtype=float)
            xj = feat_df[fj].to_numpy(dtype=float)
            mask = ~np.isnan(xi) & ~np.isnan(xj)
            n = int(mask.sum())
            if n < 3:
                intercorr[f"{fi}__vs__{fj}"] = {"rho": None, "p": None, "n": n}
                continue
            rho, p = spearmanr(xi[mask], xj[mask])
            intercorr[f"{fi}__vs__{fj}"] = {"rho": float(rho), "p": float(p), "n": n}

    print("\n=== Pairwise feature intercorrelations (Spearman) ===")
    for k, d in intercorr.items():
        rho_s = "n/a" if d["rho"] is None else f"{d['rho']:+.3f}"
        print(f"  {k:45s} rho={rho_s:>8s}  n={d['n']}")

    # ---------------------------------------------------------------- LOO-longitude-stripe AUC
    n_pos = int(feat_df["significant"].sum())
    degenerate = n_pos <= 1
    print(f"\n=== Leave-one-longitude-stripe-out logistic AUC (significant = p_train<0.05) ===")
    print(f"n positive = {n_pos} / {len(feat_df)}  "
          f"{'(DEGENERATE: AUC undefined/unstable with a single positive)' if degenerate else ''}")

    X_full = feat_df[FEATURES].to_numpy(dtype=float, copy=True)
    # median-impute missing (b_value can be None for n<50); track which columns had gaps
    col_medians = np.nanmedian(X_full, axis=0)
    n_missing = {f: int(np.isnan(X_full[:, i]).sum()) for i, f in enumerate(FEATURES)}
    for i in range(X_full.shape[1]):
        X_full[np.isnan(X_full[:, i]), i] = col_medians[i]
    y_full = feat_df["significant"].to_numpy(dtype=int)
    stripes = feat_df["bin_lon0"].to_numpy()
    uniq_stripes = np.unique(stripes)

    oof_pred = np.full(len(feat_df), np.nan)
    fold_status = []
    for lon0 in uniq_stripes:
        test_mask = stripes == lon0
        train_mask = ~test_mask
        y_tr = y_full[train_mask]
        if len(np.unique(y_tr)) < 2:
            # training fold has a single class -> logistic regression degenerates;
            # predict the (constant) training class's empirical rate.
            fold_status.append({"lon0": float(lon0), "status": "single_class_train_fold",
                                 "n_test": int(test_mask.sum())})
            oof_pred[test_mask] = float(y_tr.mean())
            continue
        scaler = StandardScaler()
        X_tr = scaler.fit_transform(X_full[train_mask])
        X_te = scaler.transform(X_full[test_mask])
        clf = LogisticRegression(max_iter=1000, random_state=SEED)
        clf.fit(X_tr, y_tr)
        oof_pred[test_mask] = clf.predict_proba(X_te)[:, 1]
        fold_status.append({"lon0": float(lon0), "status": "fit", "n_test": int(test_mask.sum())})

    try:
        auc = float(roc_auc_score(y_full, oof_pred))
    except ValueError as e:
        auc = None
        print(f"  roc_auc_score failed: {e}")

    print(f"LOO-stripe AUC = {auc}  (stripes = {len(uniq_stripes)}, "
          f"features median-imputed where missing: {n_missing})")

    # ---------------------------------------------------------------- write results
    out = {
        "exploratory": True,
        "protocol": "OVERNIGHT_PREDICTION_PROTOCOL.md EXP-B",
        "catalog": CATALOG, "split_utc": str(SPLIT), "train_years": TRAIN_YEARS,
        "bin_deg": BIN_DEG, "grid_origin": {"lon": LON_ORIGIN, "lat": LAT_ORIGIN},
        "mc": MC, "half_bin_correction": HALF_BIN, "min_n_bvalue": MIN_N_BVAL,
        "swarm_km": SWARM_KM, "swarm_days": SWARM_DAYS,
        "n_bins": len(feat_df),
        "features_per_bin": feat_df.to_dict(orient="records"),
        "spearman_vs_a_b": spearman_vs_a_b,
        "spearman_vs_p_train": {
            "note": "computed against -p_train so positive rho = feature co-varies with "
                    "significance (lower p_train); this is the same ranking as spearman "
                    "vs p_train with the sign of rho flipped.",
            **spearman_vs_p_train,
        },
        "feature_intercorrelations": intercorr,
        "auc": {
            "label": "significant (p_train < 0.05)",
            "n_positive": n_pos,
            "n_total": len(feat_df),
            "degenerate_single_positive": degenerate,
            "method": "leave-one-longitude-stripe-out logistic regression, "
                      "median-imputed missing features, standardized within fold",
            "n_stripes": len(uniq_stripes),
            "features": FEATURES,
            "n_missing_per_feature": n_missing,
            "fold_status": fold_status,
            "auc": auc,
        },
        "runtime_s": round(time.time() - t_start, 1),
    }
    (HERE / "results_exp_b.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

    print("\n================ EXP-B SUMMARY ================")
    print(f"{len(feat_df)} bins, {n_pos} significant (p_train<0.05)")
    print(f"{'feature':22s} {'rho vs a_b':>10s} {'p':>8s} {'rho vs -p_train':>16s} {'p':>8s}")
    for f in FEATURES:
        da = spearman_vs_a_b[f]
        dp = spearman_vs_p_train[f]
        ra = "n/a" if da["rho"] is None else f"{da['rho']:+.3f}"
        pa = "n/a" if da["p"] is None else f"{da['p']:.4f}"
        rp = "n/a" if dp["rho"] is None else f"{dp['rho']:+.3f}"
        pp = "n/a" if dp["p"] is None else f"{dp['p']:.4f}"
        print(f"{f:22s} {ra:>10s} {pa:>8s} {rp:>16s} {pp:>8s}")
    print(f"LOO-stripe AUC (degenerate single positive): {auc}")
    print(f"runtime {time.time()-t_start:.1f}s -> results_exp_b.json")


if __name__ == "__main__":
    main()
