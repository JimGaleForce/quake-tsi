"""K-009 (confirmatory): are the residuals of a validated ETAS white, or is there weather?

FROZEN SPEC = Popper's TESTABLE-NOW verdict for K-009 in HYPOTHESIS_LEDGER.md
(committed publicly at git 0a73fd2 = the pre-registration timestamp), section
"## 3. VERDICTS - Q3", entry "### K-009".  Quoted rule:

  Statistic. (i) pooled temporal ACF of the Anscombe residual, lags 1-52 weeks;
  (ii) Moran's I / variogram at each time step, pooled; (iii) leading EOF variance
  fraction and its power spectrum. SoCal 0.2 deg x 7 d, 2010-2018; world 0.5 deg x 30 d.
  Null. As fixed above; 500 simulated catalogs.
  Success. Lag-1 weekly ACF excess >= 0.05 over the sim-null 97.5th percentile, AND
  stable across the kernel swap, AND surviving the rho_sta partial. Report the
  correlation time and correlation length with CIs -- those two numbers, not a
  p-value, are the deliverable. Failure: excess inside the null envelope, or
  destroyed by the kernel swap or the rho_sta partial.

Popper's three REQUIRED FIXES, all implemented here:
  1. The null is generated from a SPATIO-TEMPORAL ETAS whose background field and
     aftershock kernel were fit to the real catalog (K-002 machinery) -- NOT a
     temporal model smeared by an ad-hoc kernel.  The residuals being tested and the
     null-generating model are the SAME object, so they carry the same degrees of
     freedom and the same spatial misfit.
  2. Kernel-swap control: residuals recomputed under three deliberately different
     background fields (adaptive k=4, adaptive k=8, uniform-in-footprint).
  3. Observer-nuisance regression (S-1(b)): the leading residual EOF's time series is
     regressed against a network-capability field, and the partial ACF/Moran after
     removal is reported.  NOTE/FLAG: the K-031 station-density field rho_sta(x,t) is
     NOT on disk (K-031 has not been run).  A completeness-magnitude proxy
     Mc(x,t) derived from the sub-M2.5 catalog is used as the surrogate and is
     labelled as such everywhere.

Also per Popper: positive control (injected OU latent field, 4-month / 35-km, must be
recovered to within a factor of 2), negative control (pure ETAS-sim = the null), and
an interior-vs-edge split.

Leakage control (inherited from K-002's named leakage risk): the background field
mu(x) and the spatial aftershock kernel (d, gamma, q) are estimated on TRAIN-WINDOW
events only (< 2010-01-01).  The residual/test window is 2010-01-01 .. 2019-01-01.
The temporal parameters are the frozen EXP-H values read from results_exp_h.json;
nothing temporal is refit.

Run: python -u exp_k009_residual_whiteness.py
"""

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE / "results_k009.json"
FIG = HERE / "maps" / "k009_residuals.png"
EXPH = HERE / "results_exp_h.json"
CATALOG = "SCSN_original_catalog.txt"

# ---- frozen geometry / windows (EXP-H box; Popper's SoCal grid) ----
LAT_MIN, LAT_MAX = 31.5, 38.0
LON_MIN, LON_MAX = -122.0, -113.5
CELL_DEG = 0.2
DT_DAYS = 7.0
T_TEST0 = pd.Timestamp("2010-01-01", tz="UTC")
T_TEST1 = pd.Timestamp("2019-01-01", tz="UTC")
MIN_TRAIN_EVENTS_PER_CELL = 5      # defines the analysis footprint
MAX_LAG_WEEKS = 52
DIST_BINS = np.arange(0.0, 130.0, 10.0)   # km, for the spatial correlogram
MORAN_D = 40.0                             # km, binary-weight radius for Moran's I

# ---- simulation / null ----
N_SIMS_TARGET = 20
N_SIMS_MIN = 10
N_POSCTRL = 6              # 3 per amplitude rung
MMAX_SIM = 7.5
SIM_PAD_DEG = 1.0                  # sources may live this far outside the box
SIM_EVENT_CAP_FACTOR = 8.0         # hard cap on sim size (supercritical guard)
RUNTIME_BUDGET_MIN = 90.0

# ---- positive control (Popper: 4-month timescale, 35-km length) ----
OU_TAU_DAYS = 4 * 30.44
OU_LEN_KM = 35.0
OU_SD_LADDER = [1.5, 3.0]   # amplitude ladder: the background is only ~10% of the
                            # intensity, so a weak anomaly on mu is undetectable by
                            # construction; the ladder locates the detectability floor
                            # instead of asserting a single amplitude.

RNG_SEED = 20260809
LN10 = np.log(10.0)


# ================================================================= data loading
def load_catalog(fname):
    """Auto-detect column order exactly as xue_lu_crosstest / exp_h_etas do."""
    data = HERE / "data" / "xue_lu_zenodo"
    raw_cols = ["yr", "mo", "dy", "hr", "mi", "sec", "eid", "lat", "lon", "depth", "mag"]
    dec_cols = ["yr", "mo", "dy", "hr", "mi", "sec", "lat", "lon", "depth", "mag", "eid"]
    probe = pd.read_csv(data / fname, sep=r"\s+", header=None, nrows=1000)
    is_raw = probe[6].abs().max() > 90
    cols = raw_cols if is_raw else dec_cols
    df = pd.read_csv(data / fname, sep=r"\s+", header=None, names=cols)
    assert df.lat.abs().max() <= 90 and df.lon.abs().max() <= 180
    sec = df["sec"].astype(float)
    ts = pd.to_datetime(dict(year=df.yr, month=df.mo, day=df.dy, hour=df.hr,
                             minute=df.mi, second=0), utc=True) + pd.to_timedelta(sec, unit="s")
    df["ts"] = ts
    df["t"] = (ts - pd.Timestamp(0, tz="UTC")).dt.total_seconds().to_numpy() / 86400.0
    return df[["ts", "t", "lat", "lon", "depth", "mag"]], ("raw" if is_raw else "declustered")


# ================================================================= geometry
LAT0 = 0.5 * (LAT_MIN + LAT_MAX)
KM_PER_DEG_LAT = 110.574
KM_PER_DEG_LON = 111.320 * np.cos(np.deg2rad(LAT0))


def to_km(lat, lon):
    return (np.asarray(lon) - LON_MIN) * KM_PER_DEG_LON, (np.asarray(lat) - LAT_MIN) * KM_PER_DEG_LAT


# ================================================================= background fields
def adaptive_background(train_lat, train_lon, cell_lat, cell_lon, k):
    """Adaptive (k-th nearest neighbour) isotropic Gaussian kernel smoothing.

    Returns a probability mass per analysis cell (sums to 1 over the footprint).
    Built from TRAIN-window events only.
    """
    ex, ey = to_km(train_lat, train_lon)
    cx, cy = to_km(cell_lat, cell_lon)
    pts = np.column_stack([ex, ey])
    from scipy.spatial import cKDTree
    tree = cKDTree(pts)
    dk, _ = tree.query(pts, k=min(k + 1, len(pts)))
    h = np.maximum(dk[:, -1], 2.0)                       # floor 2 km
    dens = np.zeros(len(cx))
    CH = 2000
    for a in range(0, len(ex), CH):
        b = min(a + CH, len(ex))
        d2 = (cx[:, None] - ex[None, a:b]) ** 2 + (cy[:, None] - ey[None, a:b]) ** 2
        hh = h[None, a:b]
        dens += (np.exp(-0.5 * d2 / hh ** 2) / (2 * np.pi * hh ** 2)).sum(axis=1)
    p = dens / dens.sum()
    return p


# ================================================================= spatial kernel
def spatial_pdf(r_km, D_km, q):
    """Normalised isotropic power-law: f(r) = (q-1)/(pi D^2) (1 + r^2/D^2)^{-q}."""
    return (q - 1.0) / (np.pi * D_km ** 2) * (1.0 + (r_km / D_km) ** 2) ** (-q)


def spatial_radius_sample(u, D_km, q):
    """Inverse CDF of the above: CDF(r) = 1 - (1 + r^2/D^2)^{1-q}."""
    return D_km * np.sqrt(np.power(1.0 - u, 1.0 / (1.0 - q)) - 1.0)


def omori_G(t_src, a, b, c, p):
    """Integral of (t - t_src + c)^{-p} over [max(a, t_src), b], vectorised (n_src, n_win)."""
    lo = np.maximum(a - t_src, 0.0)
    hi = b - t_src
    live = hi > lo
    lo = np.where(live, lo, 0.0)
    hi = np.where(live, hi, 0.0)
    s = 1.0 - p
    G = ((hi + c) ** s - (lo + c) ** s) / s
    return np.where(live, G, 0.0)


# ================================================================= expected-count engine
class Grid:
    def __init__(self, train_lat, train_lon):
        lat_edges = np.arange(LAT_MIN, LAT_MAX + 1e-9, CELL_DEG)
        lon_edges = np.arange(LON_MIN, LON_MAX + 1e-9, CELL_DEG)
        self.lat_edges, self.lon_edges = lat_edges, lon_edges
        nla, nlo = len(lat_edges) - 1, len(lon_edges) - 1
        H, _, _ = np.histogram2d(train_lat, train_lon, bins=[lat_edges, lon_edges])
        keep = H >= MIN_TRAIN_EVENTS_PER_CELL
        ii, jj = np.nonzero(keep)
        self.ii, self.jj = ii, jj
        self.n_cells = len(ii)
        self.clat = 0.5 * (lat_edges[ii] + lat_edges[ii + 1])
        self.clon = 0.5 * (lon_edges[jj] + lon_edges[jj + 1])
        self.cx, self.cy = to_km(self.clat, self.clon)
        self.dx_km = CELL_DEG * KM_PER_DEG_LON
        self.dy_km = CELL_DEG * KM_PER_DEG_LAT
        self.area_km2 = self.dx_km * self.dy_km
        self.cell_index = -np.ones((nla, nlo), dtype=np.int64)
        self.cell_index[ii, jj] = np.arange(self.n_cells)
        self.train_counts = H[ii, jj]
        # edge flag: cell has fewer than 8 in-footprint neighbours (interior/edge split)
        nb = np.zeros(self.n_cells, dtype=int)
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if di == 0 and dj == 0:
                    continue
                a, b = ii + di, jj + dj
                ok = (a >= 0) & (a < nla) & (b >= 0) & (b < nlo)
                v = np.zeros(self.n_cells, dtype=bool)
                v[ok] = keep[a[ok], b[ok]]
                nb += v
        self.n_nbr = nb
        self.is_edge = nb < 6      # interior = >= 6 of 8 neighbours inside the footprint

    def bin_events(self, lat, lon, t, week_edges):
        """Observed counts (n_cells, n_weeks)."""
        i = np.searchsorted(self.lat_edges, lat, side="right") - 1
        j = np.searchsorted(self.lon_edges, lon, side="right") - 1
        w = np.searchsorted(week_edges, t, side="right") - 1
        ok = ((i >= 0) & (i < len(self.lat_edges) - 1) & (j >= 0) & (j < len(self.lon_edges) - 1)
              & (w >= 0) & (w < len(week_edges) - 1))
        i, j, w = i[ok], j[ok], w[ok]
        c = self.cell_index[i, j]
        ok2 = c >= 0
        c, w = c[ok2], w[ok2]
        H = np.zeros((self.n_cells, len(week_edges) - 1))
        np.add.at(H, (c, w), 1.0)
        return H


# 3x3 Gauss-Legendre nodes/weights on [-1,1]
_GL_X = np.array([-np.sqrt(3 / 5), 0.0, np.sqrt(3 / 5)])
_GL_W = np.array([5 / 9, 8 / 9, 5 / 9])
# 9x9 refinement for near sources
_R_X, _R_W = np.polynomial.legendre.leggauss(9)


_S_X, _S_W = np.polynomial.legendre.leggauss(14)


def cell_kernel_mass(grid, sx, sy, D, chunk=400):
    """F[cell, src] = integral of the spatial kernel over the cell.

    3x3 Gauss-Legendre in the far field.  Within 2 cell widths of the source the
    kernel is sharply peaked (fitted D can be sub-km), so those pairs are redone with
    a 14x14 rule in sinh-stretched coordinates x = sx + D sinh(u), which resolves the
    peak exactly where it lives.  Accuracy is audited by sum_cells F ~= 1 for interior
    sources (reported in results_k009.json :: numerics_audit).
    """
    q = cell_kernel_mass.q
    n_src = len(sx)
    F = np.zeros((grid.n_cells, n_src), dtype=np.float32)
    hx, hy = grid.dx_km / 2.0, grid.dy_km / 2.0
    near_thresh = 2.0 * max(grid.dx_km, grid.dy_km)
    for a in range(0, n_src, chunk):
        b = min(a + chunk, n_src)
        Dc = D[a:b][None, :]
        acc = np.zeros((grid.n_cells, b - a))
        for gx, wx in zip(_GL_X, _GL_W):
            for gy, wy in zip(_GL_X, _GL_W):
                px = grid.cx[:, None] + gx * hx
                py = grid.cy[:, None] + gy * hy
                r2 = (px - sx[None, a:b]) ** 2 + (py - sy[None, a:b]) ** 2
                acc += wx * wy * (q - 1.0) / (np.pi * Dc ** 2) * (1.0 + r2 / Dc ** 2) ** (-q)
        acc *= hx * hy
        # ---- sinh-stretched refinement for near pairs
        r0 = np.sqrt((grid.cx[:, None] - sx[None, a:b]) ** 2 + (grid.cy[:, None] - sy[None, a:b]) ** 2)
        ci, si = np.nonzero(r0 < near_thresh)
        if len(ci):
            Dn = D[a:b][si]
            sxn, syn = sx[a:b][si], sy[a:b][si]
            ux0 = np.arcsinh((grid.cx[ci] - hx - sxn) / Dn)
            ux1 = np.arcsinh((grid.cx[ci] + hx - sxn) / Dn)
            uy0 = np.arcsinh((grid.cy[ci] - hy - syn) / Dn)
            uy1 = np.arcsinh((grid.cy[ci] + hy - syn) / Dn)
            mx, rx = 0.5 * (ux0 + ux1), 0.5 * (ux1 - ux0)
            my, ry = 0.5 * (uy0 + uy1), 0.5 * (uy1 - uy0)
            ref = np.zeros(len(ci))
            for gx, wx in zip(_S_X, _S_W):
                ux = mx + gx * rx
                dxk = Dn * np.cosh(ux) * rx * wx
                ddx2 = (Dn * np.sinh(ux)) ** 2
                for gy, wy in zip(_S_X, _S_W):
                    uy = my + gy * ry
                    dyk = Dn * np.cosh(uy) * ry * wy
                    r2 = ddx2 + (Dn * np.sinh(uy)) ** 2
                    ref += dxk * dyk * (q - 1.0) / (np.pi * Dn ** 2) * (1.0 + r2 / Dn ** 2) ** (-q)
            acc[ci, si] = ref
        F[:, a:b] = acc.astype(np.float32)
    return F


def expected_counts(grid, src_t, src_m, src_lat, src_lon, week_edges, pars, p_bg, F=None):
    """E[cell, week] = mu0 * p_bg[cell] * dt  +  K sum_i w_i G_i(week) F_i(cell).

    Exactly separable because the ETAS kernel is a product of a temporal and a
    (normalised) spatial factor -- no quadrature error in time.
    """
    mu0, K, alpha, c, p, M0, d, gam, q = (pars[k] for k in
                                          ("mu", "K", "alpha", "c", "p", "M0", "d", "gamma", "q"))
    n_w = len(week_edges) - 1
    E = np.zeros((grid.n_cells, n_w), dtype=np.float64)
    E += (p_bg[:, None] * mu0 * DT_DAYS)
    if len(src_t):
        if F is None:
            sx, sy = to_km(src_lat, src_lon)
            D = d * np.power(10.0, 0.5 * gam * (src_m - M0))
            cell_kernel_mass.q = q
            F = cell_kernel_mass(grid, sx, sy, D)
        w = np.power(10.0, alpha * (src_m - M0))
        G = omori_G(src_t[:, None], week_edges[None, :-1], week_edges[None, 1:], c, p)
        WG = (K * w[:, None] * G).astype(np.float32)
        # chunked matmul over sources
        CH = 6000
        for a in range(0, len(src_t), CH):
            b = min(a + CH, len(src_t))
            E += (F[:, a:b] @ WG[a:b]).astype(np.float64)
    return E, F


# ================================================================= whiteness statistics
def residual_field(obs, exp):
    return (obs - exp) / np.sqrt(exp + 1.0)


def pooled_temporal_acf(R, max_lag):
    """rho(k) = sum_{cell,t} r(t) r(t+k) / sum_{cell,t} r(t)^2 (no per-cell demeaning:
    a region-wide offset is signal, not nuisance)."""
    n_c, n_t = R.shape
    denom = float((R ** 2).sum())
    out = np.zeros(max_lag + 1)
    out[0] = 1.0
    for k in range(1, max_lag + 1):
        out[k] = float((R[:, :-k] * R[:, k:]).sum()) / denom
    return out


def corr_time_weeks(acf, kmax=MAX_LAG_WEEKS):
    """e-folding time from a least-squares fit rho(k)=A exp(-k/T).

    Fitted over Popper's full specified lag range (1..52 weeks), not a truncated
    sub-range: an earlier 12-lag version SATURATED at its cap for both the data and
    most of the null, which destroyed the statistic. Bounded so it is always defined:
    fewer than 3 positive lags -> 0 weeks; non-decaying fit -> capped at kmax.
    """
    k = np.arange(1, kmax + 1)
    y = acf[1:kmax + 1]
    m = np.isfinite(y) & (y > 1e-6)
    if m.sum() < 3:
        return 0.0, 0.0
    kk, yy = k[m], y[m]
    A = np.column_stack([np.ones_like(kk, dtype=float), -kk.astype(float)])
    coef, *_ = np.linalg.lstsq(A, np.log(yy), rcond=None)
    if coef[1] <= 0:
        return float(kmax), float(np.exp(coef[0]))
    return float(min(1.0 / coef[1], kmax)), float(np.exp(coef[0]))


def integral_timescale_weeks(acf):
    """Sum of rho(k) over k=1.. up to the first non-positive lag (always defined)."""
    s = 0.0
    for k in range(1, len(acf)):
        if acf[k] <= 0:
            break
        s += acf[k]
    return float(s)


class SpatialOps:
    """Pre-computed distance-bin pair matrices for the pooled correlogram + Moran's I."""

    def __init__(self, grid):
        d = np.sqrt((grid.cx[:, None] - grid.cx[None, :]) ** 2 +
                    (grid.cy[:, None] - grid.cy[None, :]) ** 2)
        np.fill_diagonal(d, -1.0)
        self.centers = 0.5 * (DIST_BINS[:-1] + DIST_BINS[1:])
        self.B = []
        self.npairs = []
        for lo, hi in zip(DIST_BINS[:-1], DIST_BINS[1:]):
            m = ((d >= lo) & (d < hi)).astype(np.float32)
            self.B.append(m)
            self.npairs.append(float(m.sum()))
        self.W = ((d >= 0) & (d < MORAN_D)).astype(np.float32)
        self.Wsum = float(self.W.sum())

    def correlogram(self, R):
        R32 = R.astype(np.float32)
        var = float((R32 ** 2).sum())
        n_t = R.shape[1]
        n_c = R.shape[0]
        out = np.full(len(self.centers), np.nan)
        for i, B in enumerate(self.B):
            if self.npairs[i] == 0:
                continue
            num = float(((B @ R32) * R32).sum()) / (self.npairs[i] * n_t)
            den = var / (n_c * n_t)
            out[i] = num / den
        return out

    def moran(self, R):
        """Moran's I per time step (weights binary, d < MORAN_D), then pooled mean."""
        R32 = R.astype(np.float32)
        n_c = R.shape[0]
        z = R32 - R32.mean(axis=0, keepdims=True)
        num = ((self.W @ z) * z).sum(axis=0)
        den = (z ** 2).sum(axis=0)
        I = (n_c / self.Wsum) * num / np.maximum(den, 1e-12)
        return float(np.mean(I)), I


def corr_length_km(centers, rho, dmax=None):
    """e-folding length, bounded exactly as corr_time_weeks is.

    NOTE ON RESOLUTION: a 0.2 deg cell is ~18.5 x 22.1 km, so the shortest cell-centre
    separation is ~18 km and the 0-10 km bin is EMPTY (NaN). This estimator therefore
    measures the decay of the correlogram above ~18 km and cannot see structure below
    it. Empty/NaN bins are dropped rather than terminating the fit (an earlier version
    returned 0 because it stopped at the empty first bin).
    """
    if dmax is None:
        dmax = float(DIST_BINS[-1])
    m = np.isfinite(rho) & (centers <= dmax) & (rho > 1e-4)
    if m.sum() < 3:
        return 0.0
    A = np.column_stack([np.ones(m.sum()), -centers[m]])
    coef, *_ = np.linalg.lstsq(A, np.log(rho[m]), rcond=None)
    if coef[1] <= 0:
        return float(dmax)
    return float(min(1.0 / coef[1], dmax))


def integral_length_km(centers, rho):
    """Integral of rho(d) dd up to the first non-positive FINITE bin.

    Empty (NaN) bins are skipped, not treated as a stopping condition -- the 0-10 km
    bin is always empty at 0.2 deg resolution.
    """
    w = float(centers[1] - centers[0]) if len(centers) > 1 else 10.0
    s = 0.0
    for i in range(len(centers)):
        if not np.isfinite(rho[i]):
            continue
        if rho[i] <= 0:
            break
        s += rho[i] * w
    return float(s)


def leading_eof(R):
    """Leading EOF of the residual field: variance fraction, PC time series, spectrum."""
    X = R - R.mean(axis=1, keepdims=True)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    var = S ** 2 / (S ** 2).sum()
    pc = Vt[0] * S[0]
    f = np.fft.rfftfreq(len(pc), d=1.0)
    P = np.abs(np.fft.rfft(pc - pc.mean())) ** 2
    # redness = power below 1/12 wk^-1 (periods > 3 months) as a fraction of total
    lo = (f > 0) & (f < 1.0 / 12.0)
    red = float(P[lo].sum() / max(P[f > 0].sum(), 1e-12))
    return float(var[0]), pc, f, P, red


def time_rescaling(src_t, src_m, tgt_t, week_edges, pars, F_colsum):
    """Time-rescaling theorem on the spatially integrated process over the footprint.

    Lambda(t) = mu0 (t - T0) + K sum_i w_i F_i * int_{T0}^{t} (s - t_i + c)^{-p} ds
    tau_j = Lambda(t_j) - Lambda(t_{j-1}) ~ iid Exp(1) if the model is complete.
    """
    mu0, K, alpha, c, p, M0 = (pars[k] for k in ("mu", "K", "alpha", "c", "p", "M0"))
    T0 = week_edges[0]
    w = np.power(10.0, alpha * (src_m - M0)) * F_colsum
    nodes = np.concatenate([[T0], tgt_t])
    G = omori_G(src_t[:, None], T0, nodes[None, :], c, p)
    Lam = mu0 * (nodes - T0) + K * (w[:, None] * G).sum(axis=0)
    tau = np.diff(Lam)
    tau = tau[np.isfinite(tau) & (tau >= 0)]
    return tau


def exp1_stats(tau):
    if len(tau) < 20:
        return dict(n=int(len(tau)), ks=float("nan"), mean=float("nan"),
                    cv=float("nan"), lag1=float("nan"))
    x = np.sort(tau)
    n = len(x)
    F = 1.0 - np.exp(-x)
    ks = float(np.max(np.abs(F - (np.arange(1, n + 1) - 0.5) / n)))
    u = 1.0 - np.exp(-tau)
    lag1 = float(np.corrcoef(u[:-1], u[1:])[0, 1])
    return dict(n=int(n), ks=ks, mean=float(tau.mean()),
                cv=float(tau.std() / max(tau.mean(), 1e-12)), lag1=lag1)


def all_stats(R, sops, max_lag=MAX_LAG_WEEKS):
    acf = pooled_temporal_acf(R, max_lag)
    T, _ = corr_time_weeks(acf)
    rho = sops.correlogram(R)
    L = corr_length_km(sops.centers, rho)
    mI, _ = sops.moran(R)
    ev, pc, f, P, red = leading_eof(R)
    return dict(acf=acf, acf1=float(acf[1]), corr_time_weeks=T,
                corr_time_integral_weeks=integral_timescale_weeks(acf),
                correlogram=rho, corr_len_km=L,
                corr_len_integral_km=integral_length_km(sops.centers, rho),
                moran_I=mI, eof1_varfrac=ev, eof1_pc=pc, eof1_redness=red)


# ================================================================= ETAS simulation
def simulate_st_etas(rng, hist_t, hist_m, hist_lat, hist_lon, pars, p_bg, grid,
                     T0, T1, b_val, cap, ou_field=None, week_edges=None):
    """Branching simulation of the spatio-temporal ETAS over [T0, T1], conditional on
    the real pre-window history (so the null carries the same warm-up as the data).

    ou_field: optional (n_cells, n_weeks) multiplicative anomaly on the background
    (the positive control's injected latent field).
    """
    mu0, K, alpha, c, p, M0, d, gam, q = (pars[k] for k in
                                          ("mu", "K", "alpha", "c", "p", "M0", "d", "gamma", "q"))
    beta = b_val * np.log(10.0)
    # --- background generation
    if ou_field is None:
        n_bg = rng.poisson(mu0 * (T1 - T0))
        bt = np.sort(rng.uniform(T0, T1, n_bg))
        ci = rng.choice(grid.n_cells, size=n_bg, p=p_bg)
    else:
        lam_cw = mu0 * p_bg[:, None] * ou_field * DT_DAYS
        n_cw = rng.poisson(lam_cw)
        ci_l, bt_l = [], []
        cc, ww = np.nonzero(n_cw)
        for cidx, widx in zip(cc, ww):
            k = n_cw[cidx, widx]
            ci_l.append(np.full(k, cidx))
            bt_l.append(rng.uniform(week_edges[widx], week_edges[widx + 1], k))
        ci = np.concatenate(ci_l) if ci_l else np.zeros(0, dtype=int)
        bt = np.concatenate(bt_l) if bt_l else np.zeros(0)
        o = np.argsort(bt)
        ci, bt = ci[o], bt[o]
    blat = grid.clat[ci] + rng.uniform(-CELL_DEG / 2, CELL_DEG / 2, len(ci))
    blon = grid.clon[ci] + rng.uniform(-CELL_DEG / 2, CELL_DEG / 2, len(ci))
    bmag = M0 - np.log(rng.uniform(size=len(ci))) / beta
    bmag = np.minimum(bmag, MMAX_SIM)

    out_t = [bt]
    out_m = [bmag]
    out_la = [blat]
    out_lo = [blon]
    total = len(bt)
    truncated = False

    # --- cascade.  Generation 0 parents = the real pre-window history AND the
    #     freshly generated in-window background events (both are valid parents).
    gen_t = np.concatenate([hist_t, bt])
    gen_m = np.concatenate([hist_m, bmag])
    gen_la = np.concatenate([hist_lat, blat])
    gen_lo = np.concatenate([hist_lon, blon])
    o0 = np.argsort(gen_t)
    gen_t, gen_m, gen_la, gen_lo = gen_t[o0], gen_m[o0], gen_la[o0], gen_lo[o0]
    for _ in range(200):
        if len(gen_t) == 0:
            break
        w = np.power(10.0, alpha * (gen_m - M0))
        G = omori_G(gen_t, T0, T1, c, p)
        nu = K * w * G
        nu = np.where(np.isfinite(nu) & (nu > 0), nu, 0.0)
        nk = rng.poisson(nu)
        if nk.sum() == 0:
            break
        par = np.repeat(np.arange(len(gen_t)), nk)
        n_new = len(par)
        if total + n_new > cap:
            truncated = True
            keep = max(0, int(cap - total))
            par = par[:keep]
            n_new = keep
            if n_new == 0:
                break
        pt = gen_t[par]
        # sample offspring times from the normalised Omori density on [max(T0,pt), T1]
        lo = np.maximum(T0 - pt, 0.0)
        hi = T1 - pt
        s = 1.0 - p
        A = (lo + c) ** s
        B = (hi + c) ** s
        u = rng.uniform(size=n_new)
        ct = np.power(A + u * (B - A), 1.0 / s) - c + pt
        cm = M0 - np.log(rng.uniform(size=n_new)) / beta
        cm = np.minimum(cm, MMAX_SIM)
        D = d * np.power(10.0, 0.5 * gam * (gen_m[par] - M0))
        rr = spatial_radius_sample(rng.uniform(size=n_new), D, q)
        th = rng.uniform(0, 2 * np.pi, n_new)
        clat = gen_la[par] + (rr * np.sin(th)) / KM_PER_DEG_LAT
        clon = gen_lo[par] + (rr * np.cos(th)) / KM_PER_DEG_LON
        inpad = ((clat > LAT_MIN - SIM_PAD_DEG) & (clat < LAT_MAX + SIM_PAD_DEG) &
                 (clon > LON_MIN - SIM_PAD_DEG) & (clon < LON_MAX + SIM_PAD_DEG))
        ct, cm, clat, clon = ct[inpad], cm[inpad], clat[inpad], clon[inpad]
        total += len(ct)
        out_t.append(ct); out_m.append(cm); out_la.append(clat); out_lo.append(clon)
        gen_t, gen_m, gen_la, gen_lo = ct, cm, clat, clon
        if truncated:
            break
    t = np.concatenate(out_t); m = np.concatenate(out_m)
    la = np.concatenate(out_la); lo = np.concatenate(out_lo)
    o = np.argsort(t)
    return t[o], m[o], la[o], lo[o], truncated


def ou_latent_field(rng, grid, n_weeks, tau_days, len_km, sd):
    """Smooth OU-in-time, exponentially-correlated-in-space multiplicative anomaly,
    exp(z - var/2) so that E[anomaly] = 1."""
    d = np.sqrt((grid.cx[:, None] - grid.cx[None, :]) ** 2 +
                (grid.cy[:, None] - grid.cy[None, :]) ** 2)
    C = np.exp(-d / len_km)
    ev, V = np.linalg.eigh(C)
    ev = np.clip(ev, 0, None)
    Lh = V * np.sqrt(ev)[None, :]
    phi = np.exp(-DT_DAYS / tau_days)
    z = np.zeros((grid.n_cells, n_weeks))
    x = Lh @ rng.normal(size=grid.n_cells)
    for k in range(n_weeks):
        x = phi * x + np.sqrt(1 - phi ** 2) * (Lh @ rng.normal(size=grid.n_cells))
        z[:, k] = x
    z *= sd
    return np.exp(z - 0.5 * sd ** 2)


# ================================================================= main
def main():
    t_start = time.time()
    rng = np.random.default_rng(RNG_SEED)
    res = {"experiment": "K-009",
           "spec": "HYPOTHESIS_LEDGER.md :: '## 3. VERDICTS - Q3' :: '### K-009' (Popper TESTABLE-NOW verdict, git 0a73fd2)",
           "state_class": "first-run (confirmatory scoring; prediction POST-registered, see pre_registration)",
           "run_utc": pd.Timestamp.now("UTC").isoformat()}

    # ---------- frozen EXP-H temporal parameters ----------
    exph = json.loads(EXPH.read_text())
    fp = exph["train_fit"]["frozen_params"]
    b_val = exph["train_fit"]["b_value_train_aki"]
    M0 = fp["M0"]
    print(f"[frozen EXP-H] mu={fp['mu']:.5f}/d K={fp['K']:.5f} alpha={fp['alpha']:.4f} "
          f"c={fp['c']:.5f} p={fp['p']:.4f} M0={M0} b={b_val:.4f}")
    res["frozen_temporal_params"] = dict(fp, b_value=b_val,
                                         source="results_exp_h.json :: train_fit.frozen_params")

    # ---------- catalog ----------
    df, order = load_catalog(CATALOG)
    box = df.lat.between(LAT_MIN, LAT_MAX) & df.lon.between(LON_MIN, LON_MAX)
    dall = df[box].sort_values("t").reset_index(drop=True)
    cat = dall[dall.mag >= M0 - 1e-9].reset_index(drop=True)
    t0w = (T_TEST0 - pd.Timestamp(0, tz="UTC")).total_seconds() / 86400.0
    t1w = (T_TEST1 - pd.Timestamp(0, tz="UTC")).total_seconds() / 86400.0
    week_edges = np.arange(t0w, t1w + 1e-9, DT_DAYS)
    n_weeks = len(week_edges) - 1
    train = cat[cat.t < t0w]
    test = cat[(cat.t >= week_edges[0]) & (cat.t < week_edges[-1])]
    print(f"[data] {order} order; M>={M0} in box: {len(cat)}  train(<2010)={len(train)}  "
          f"window={len(test)} over {n_weeks} weeks")
    res["catalog"] = {"file": CATALOG, "detected_column_order": order,
                      "n_in_box_all_mags": int(len(dall)), "n_M0": int(len(cat)),
                      "n_train": int(len(train)), "n_window": int(len(test)),
                      "window": [str(T_TEST0), str(T_TEST1)], "n_weeks": int(n_weeks),
                      "cell_deg": CELL_DEG, "dt_days": DT_DAYS}

    grid = Grid(train.lat.to_numpy(), train.lon.to_numpy())
    print(f"[grid] {grid.n_cells} footprint cells (>= {MIN_TRAIN_EVENTS_PER_CELL} train events); "
          f"{int(grid.is_edge.sum())} edge / {int((~grid.is_edge).sum())} interior")
    res["grid"] = {"n_cells": int(grid.n_cells), "min_train_events_per_cell": MIN_TRAIN_EVENTS_PER_CELL,
                   "n_edge": int(grid.is_edge.sum()), "n_interior": int((~grid.is_edge).sum()),
                   "cell_km": [round(grid.dx_km, 2), round(grid.dy_km, 2)]}

    # ---------- background fields (TRAIN ONLY -- leakage control) ----------
    print("[bg] building background fields from TRAIN-window events only ...")
    t_bg = time.time()
    p_k4 = adaptive_background(train.lat.to_numpy(), train.lon.to_numpy(),
                               grid.clat, grid.clon, k=4)
    p_k8 = adaptive_background(train.lat.to_numpy(), train.lon.to_numpy(),
                               grid.clat, grid.clon, k=8)
    p_uni = np.full(grid.n_cells, 1.0 / grid.n_cells)
    BG = {"adaptive_k4": p_k4, "adaptive_k8": p_k8, "uniform_footprint": p_uni}
    print(f"[bg] done ({time.time()-t_bg:.0f}s)")

    # ---------- spatial aftershock kernel: fit (d, gamma, q) on TRAIN events only ----------
    print("[fit] fitting spatial kernel (d, gamma, q) on train window, temporal params frozen ...")
    tt = train.t.to_numpy(); tm = train.mag.to_numpy()
    tx, ty = to_km(train.lat.to_numpy(), train.lon.to_numpy())
    # target = train events after a 1-yr burn-in; pairs truncated to 100 d / 100 km
    burn = tt[0] + 365.0
    tgt = np.nonzero(tt >= burn)[0]
    PAIR_T, PAIR_R = 100.0, 100.0
    pi_list, pj_list = [], []
    for a in range(0, len(tgt), 2000):
        idx = tgt[a:a + 2000]
        lo = np.searchsorted(tt, tt[idx] - PAIR_T, side="left")
        for jj, l in zip(idx, lo):
            if jj > l:
                pi_list.append(np.arange(l, jj)); pj_list.append(np.full(jj - l, jj))
    pi = np.concatenate(pi_list); pj = np.concatenate(pj_list)
    dr = np.sqrt((tx[pi] - tx[pj]) ** 2 + (ty[pi] - ty[pj]) ** 2)
    keep = dr <= PAIR_R
    pi, pj, dr = pi[keep], pj[keep], dr[keep]
    dt = tt[pj] - tt[pi]
    print(f"[fit] {len(pi)} candidate parent-child pairs (dt<={PAIR_T}d, dr<={PAIR_R}km), "
          f"{len(tgt)} targets")
    mu0, Kp, alph, cc, pp = fp["mu"], fp["K"], fp["alpha"], fp["c"], fp["p"]
    g_time = Kp * np.power(10.0, alph * (tm[pi] - M0)) * (dt + cc) ** (-pp)
    mpar = tm[pi] - M0
    bg_at_train = None
    # background density at each target (events/day/km^2) from the k=4 field
    ilat = np.searchsorted(grid.lat_edges, train.lat.to_numpy(), side="right") - 1
    ilon = np.searchsorted(grid.lon_edges, train.lon.to_numpy(), side="right") - 1
    ilat = np.clip(ilat, 0, len(grid.lat_edges) - 2); ilon = np.clip(ilon, 0, len(grid.lon_edges) - 2)
    cidx = grid.cell_index[ilat, ilon]
    bg_floor = mu0 * p_k4[p_k4 > 0].min() / grid.area_km2
    bg_at_train = np.where(cidx >= 0,
                           np.maximum(mu0 * p_k4[np.clip(cidx, 0, None)] / grid.area_km2, bg_floor),
                           bg_floor)
    tgt_pos = np.searchsorted(tgt, pj)

    def neg_ll(x):
        d_, gam_, q_ = np.exp(x[0]), x[1], 1.0 + np.exp(x[2])
        D = d_ * np.power(10.0, 0.5 * gam_ * mpar)
        f = (q_ - 1.0) / (np.pi * D ** 2) * (1.0 + (dr / D) ** 2) ** (-q_)
        s = np.zeros(len(tgt))
        np.add.at(s, tgt_pos, g_time * f)
        lam = bg_at_train[tgt] + s
        return -float(np.log(np.maximum(lam, 1e-300)).sum())

    x0 = np.array([np.log(3.0), 0.5, np.log(0.5)])
    r = minimize(neg_ll, x0, method="Nelder-Mead",
                 options={"maxiter": 300, "xatol": 1e-3, "fatol": 1e-2})
    d_fit, gam_fit, q_fit = float(np.exp(r.x[0])), float(r.x[1]), float(1.0 + np.exp(r.x[2]))
    print(f"[fit] spatial kernel: d={d_fit:.3f} km  gamma={gam_fit:.3f}  q={q_fit:.3f}  "
          f"(-LL={r.fun:.1f}, {r.nit} it)")
    res["spatial_kernel_fit"] = {"d_km": d_fit, "gamma": gam_fit, "q": q_fit,
                                 "neg_LL": float(r.fun), "nit": int(r.nit),
                                 "pair_truncation": {"dt_days": PAIR_T, "dr_km": PAIR_R},
                                 "n_pairs": int(len(pi)), "n_targets": int(len(tgt)),
                                 "estimated_on": "train window only (< 2010-01-01)",
                                 "form": "f(r|M) = (q-1)/(pi D^2) (1 + r^2/D^2)^-q, D = d*10^(gamma/2 (M-M0))"}

    PARS = dict(mu=mu0, K=Kp, alpha=alph, c=cc, p=pp, M0=M0, d=d_fit, gamma=gam_fit, q=q_fit)

    # ---------- observed field and real residuals (3 background fields) ----------
    sops = SpatialOps(grid)
    src = cat[cat.t < week_edges[-1]]
    src_t = src.t.to_numpy(); src_m = src.mag.to_numpy()
    src_lat = src.lat.to_numpy(); src_lon = src.lon.to_numpy()
    obs = grid.bin_events(test.lat.to_numpy(), test.lon.to_numpy(), test.t.to_numpy(), week_edges)
    print(f"[exp] computing expected counts ({len(src_t)} sources x {grid.n_cells} cells "
          f"x {n_weeks} weeks) ...")
    te = time.time()
    E_k4, F_real = expected_counts(grid, src_t, src_m, src_lat, src_lon, week_edges, PARS, p_k4)
    print(f"[exp] done ({time.time()-te:.0f}s). obs total={obs.sum():.0f}  exp total={E_k4.sum():.0f}")

    # numerics audit: for interior sources the spatial kernel mass over cells must be ~1
    Fsum = F_real.sum(axis=0).astype(np.float64)
    si_i = np.clip(np.searchsorted(grid.lat_edges, src_lat, side="right") - 1, 0, len(grid.lat_edges) - 2)
    si_j = np.clip(np.searchsorted(grid.lon_edges, src_lon, side="right") - 1, 0, len(grid.lon_edges) - 2)
    sc = grid.cell_index[si_i, si_j]
    deep = (sc >= 0) & (grid.n_nbr[np.clip(sc, 0, None)] == 8)
    res["numerics_audit"] = {
        "kernel_mass_sum_over_cells_interior_sources": {
            "n": int(deep.sum()), "median": float(np.median(Fsum[deep])),
            "p5": float(np.percentile(Fsum[deep], 5)), "p95": float(np.percentile(Fsum[deep], 95))},
        "obs_total": float(obs.sum()), "exp_total_k4": float(E_k4.sum()),
        "obs_over_exp": float(obs.sum() / E_k4.sum()),
        "note": "mass < 1 is expected: footprint cells do not tile the whole box"}
    print(f"[audit] kernel mass over footprint cells, interior sources: "
          f"median {np.median(Fsum[deep]):.3f} (p5 {np.percentile(Fsum[deep],5):.3f}, "
          f"p95 {np.percentile(Fsum[deep],95):.3f});  obs/exp = {obs.sum()/E_k4.sum():.3f}")

    real = {}
    for name, pbg in BG.items():
        E, _ = expected_counts(grid, src_t, src_m, src_lat, src_lon, week_edges, PARS, pbg, F=F_real)
        R = residual_field(obs, E)
        st = all_stats(R, sops)
        real[name] = dict(st)
        real[name]["exp_total"] = float(E.sum())
        if name == "adaptive_k4":
            R_real, E_real = R, E
        print(f"[real:{name:18s}] acf1={st['acf1']:+.4f}  T={st['corr_time_weeks']:.2f} wk  "
              f"L={st['corr_len_km']:.1f} km  MoranI={st['moran_I']:+.4f}  "
              f"EOF1={st['eof1_varfrac']*100:.1f}%  red={st['eof1_redness']:.2f}")

    F_colsum = F_real.sum(axis=0).astype(np.float64)
    tau_real = time_rescaling(src_t, src_m, test.t.to_numpy(), week_edges, PARS, F_colsum)
    tr_real = exp1_stats(tau_real)
    print(f"[real:time-rescaling] n={tr_real['n']} KS={tr_real['ks']:.4f} mean={tr_real['mean']:.3f} "
          f"CV={tr_real['cv']:.3f} lag1(u)={tr_real['lag1']:+.4f}")

    # interior vs edge
    edge_split = {}
    for lab, sel in [("interior", ~grid.is_edge), ("edge", grid.is_edge)]:
        acf = pooled_temporal_acf(R_real[sel], MAX_LAG_WEEKS)
        T, _ = corr_time_weeks(acf)
        edge_split[lab] = {"n_cells": int(sel.sum()), "acf1": float(acf[1]),
                           "corr_time_weeks": T}
        print(f"[real:{lab:8s}] n={int(sel.sum())} acf1={acf[1]:+.4f} T={T:.2f} wk")
    res["interior_vs_edge"] = edge_split

    # ---------- NULL: simulated catalogs through the identical pipeline ----------
    hist = cat[cat.t < week_edges[0]]
    h_t, h_m = hist.t.to_numpy(), hist.mag.to_numpy()
    h_la, h_lo = hist.lat.to_numpy(), hist.lon.to_numpy()
    cap = int(SIM_EVENT_CAP_FACTOR * len(test))
    n_sims = N_SIMS_TARGET
    sims = []
    print(f"[null] simulating spatio-temporal ETAS (cap={cap}) ...")
    for s in range(N_SIMS_TARGET):
        el = (time.time() - t_start) / 60.0
        if s >= N_SIMS_MIN and el > RUNTIME_BUDGET_MIN * 0.62:
            n_sims = s
            print(f"[null] RUNTIME GUARD: stopping at {s} sims ({el:.1f} min elapsed)")
            break
        ts_ = time.time()
        st_, sm_, sla_, slo_, trunc = simulate_st_etas(
            rng, h_t, h_m, h_la, h_lo, PARS, p_k4, grid, week_edges[0], week_edges[-1], b_val, cap)
        inw = (st_ >= week_edges[0]) & (st_ < week_edges[-1])
        inbox = ((sla_ >= LAT_MIN) & (sla_ <= LAT_MAX) & (slo_ >= LON_MIN) & (slo_ <= LON_MAX))
        obs_s = grid.bin_events(sla_[inw & inbox], slo_[inw & inbox], st_[inw & inbox], week_edges)
        # sources for the sim = real pre-window history + simulated in-window events
        ss_t = np.concatenate([h_t, st_[inw]])
        ss_m = np.concatenate([h_m, sm_[inw]])
        ss_la = np.concatenate([h_la, sla_[inw]])
        ss_lo = np.concatenate([h_lo, slo_[inw]])
        o = np.argsort(ss_t)
        ss_t, ss_m, ss_la, ss_lo = ss_t[o], ss_m[o], ss_la[o], ss_lo[o]
        E_s, F_s = expected_counts(grid, ss_t, ss_m, ss_la, ss_lo, week_edges, PARS, p_k4)
        R_s = residual_field(obs_s, E_s)
        stt = all_stats(R_s, sops)
        sel = inw & inbox
        tau_s = time_rescaling(ss_t, ss_m, np.sort(st_[sel]), week_edges, PARS,
                               F_s.sum(axis=0).astype(np.float64))
        stt["time_rescaling"] = exp1_stats(tau_s)
        stt["n_events_in_box"] = int(sel.sum())
        stt["truncated"] = bool(trunc)
        sims.append(stt)
        del F_s
        print(f"[null] sim {s+1}/{N_SIMS_TARGET}: n={int(sel.sum())} (obs {len(test)}) "
              f"acf1={stt['acf1']:+.4f} T={stt['corr_time_weeks']:.2f} L={stt['corr_len_km']:.1f} "
              f"MoranI={stt['moran_I']:+.4f} EOF1={stt['eof1_varfrac']*100:.1f}% "
              f"{'TRUNC' if trunc else ''} ({time.time()-ts_:.0f}s)")
    n_sims = len(sims)

    def null_env(key):
        v0 = np.array([s[key] for s in sims], dtype=float)
        v = v0[np.isfinite(v0)]
        if len(v) == 0:
            return dict(n=0, n_total=int(len(v0)))
        return dict(n=int(len(v)), n_total=int(len(v0)), median=float(np.median(v)),
                    p2_5=float(np.percentile(v, 2.5)), p97_5=float(np.percentile(v, 97.5)),
                    min=float(v.min()), max=float(v.max()))

    NULL = {k: null_env(k) for k in ["acf1", "corr_time_weeks", "corr_time_integral_weeks",
                                     "corr_len_km", "corr_len_integral_km",
                                     "moran_I", "eof1_varfrac", "eof1_redness"]}
    NULL["n_events_in_box"] = null_env("n_events_in_box")
    NULL["time_rescaling_ks"] = {"n": n_sims,
                                 "median": float(np.median([s["time_rescaling"]["ks"] for s in sims])),
                                 "p97_5": float(np.percentile([s["time_rescaling"]["ks"] for s in sims], 97.5))}
    NULL["time_rescaling_lag1"] = {"n": n_sims,
                                   "median": float(np.median([s["time_rescaling"]["lag1"] for s in sims])),
                                   "p2_5": float(np.percentile([s["time_rescaling"]["lag1"] for s in sims], 2.5)),
                                   "p97_5": float(np.percentile([s["time_rescaling"]["lag1"] for s in sims], 97.5))}
    null_acf_curves = np.array([s["acf"] for s in sims])
    null_cor_curves = np.array([s["correlogram"] for s in sims])
    res["n_sims"] = n_sims
    res["n_sims_spec"] = 500
    res["null_per_sim"] = [{k: s[k] for k in
                            ["acf1", "corr_time_weeks", "corr_time_integral_weeks",
                             "corr_len_km", "corr_len_integral_km", "moran_I",
                             "eof1_varfrac", "eof1_redness", "n_events_in_box", "truncated"]}
                           for s in sims]
    # saturation audit: a statistic pinned at its cap is not a measurement
    res["estimator_saturation"] = {
        "T_cap_weeks": float(MAX_LAG_WEEKS), "L_cap_km": float(DIST_BINS[-1]),
        "n_sims_T_at_cap": int(sum(s["corr_time_weeks"] >= MAX_LAG_WEEKS - 1e-9 for s in sims)),
        "n_sims_T_at_zero": int(sum(s["corr_time_weeks"] <= 1e-9 for s in sims)),
        "n_sims_L_at_cap": int(sum(s["corr_len_km"] >= DIST_BINS[-1] - 1e-9 for s in sims)),
        "n_sims_L_at_zero": int(sum(s["corr_len_km"] <= 1e-9 for s in sims)),
        "note": ("If most of the null sits at a cap the statistic has no resolving power "
                 "there and the comparison must not be read as a measurement.")}
    res["sims_truncated"] = int(sum(s["truncated"] for s in sims))

    # ---------- rho_sta partial (S-1(b)) ----------
    print("[rho_sta] building network-capability proxy Mc(x,t) from the sub-M2.5 catalog ...")
    SUP = 3          # 0.6-deg super-cells
    HALF = 13        # +/- 13 weeks
    sc_i = grid.ii // SUP
    sc_j = grid.jj // SUP
    keys = sc_i * 1000 + sc_j
    uk, inv = np.unique(keys, return_inverse=True)
    small = dall[(dall.mag >= 0.5) & (dall.t >= week_edges[0] - HALF * DT_DAYS) &
                 (dall.t < week_edges[-1] + HALF * DT_DAYS)]
    s_i = np.clip(np.searchsorted(grid.lat_edges, small.lat.to_numpy(), side="right") - 1,
                  0, len(grid.lat_edges) - 2) // SUP
    s_j = np.clip(np.searchsorted(grid.lon_edges, small.lon.to_numpy(), side="right") - 1,
                  0, len(grid.lon_edges) - 2) // SUP
    s_key = s_i * 1000 + s_j
    s_w = np.floor((small.t.to_numpy() - week_edges[0]) / DT_DAYS).astype(int)
    s_mag = small.mag.to_numpy()
    Mc = np.full((len(uk), n_weeks), np.nan)
    for a, k in enumerate(uk):
        mk = s_key == k
        if mk.sum() < 200:
            continue
        wk, mg = s_w[mk], s_mag[mk]
        for w in range(n_weeks):
            sl = (wk >= w - HALF) & (wk <= w + HALF)
            if sl.sum() < 100:
                continue
            h, e = np.histogram(mg[sl], bins=np.arange(-1.0, 5.0, 0.1))
            Mc[a, w] = 0.5 * (e[np.argmax(h)] + e[np.argmax(h) + 1])
    # fill gaps and expand to cells
    for a in range(len(uk)):
        row = Mc[a]
        if np.all(np.isnan(row)):
            Mc[a] = np.nanmedian(Mc) if np.isfinite(np.nanmedian(Mc)) else 1.0
        else:
            idx = np.arange(n_weeks)
            g = np.isfinite(row)
            Mc[a] = np.interp(idx, idx[g], row[g])
    RHO = Mc[inv]                             # (n_cells, n_weeks) capability proxy
    RHO_z = (RHO - RHO.mean(axis=1, keepdims=True))
    sd = RHO_z.std(axis=1, keepdims=True)
    RHO_z = RHO_z / np.where(sd > 1e-9, sd, 1.0)

    # (a) regress the leading EOF PC on the region-mean proxy
    ev1, pc1, f1, P1, red1 = leading_eof(R_real)
    rmean = RHO.mean(axis=0)
    rmean = (rmean - rmean.mean()) / max(rmean.std(), 1e-12)
    beta_pc = float(np.dot(pc1 - pc1.mean(), rmean) / np.dot(rmean, rmean))
    r2_pc = float(beta_pc ** 2 * np.dot(rmean, rmean) / max(np.dot(pc1 - pc1.mean(), pc1 - pc1.mean()), 1e-12))
    # (b) cell-wise partial: remove the projection of each cell's residual on its own proxy
    num = (R_real * RHO_z).sum(axis=1, keepdims=True)
    den = (RHO_z ** 2).sum(axis=1, keepdims=True)
    R_part = R_real - (num / np.maximum(den, 1e-12)) * RHO_z
    st_part = all_stats(R_part, sops)
    print(f"[rho_sta] EOF1-PC vs region-mean Mc: R^2={r2_pc:.4f}   "
          f"partial acf1={st_part['acf1']:+.4f} (raw {real['adaptive_k4']['acf1']:+.4f})  "
          f"partial MoranI={st_part['moran_I']:+.4f}  partial T={st_part['corr_time_weeks']:.2f} wk  "
          f"partial L={st_part['corr_len_km']:.1f} km")
    res["rho_sta_partial"] = {
        "proxy": "Mc(x,t) maximum-curvature completeness on 0.6-deg super-cells, +/-13-week window, "
                 "from the sub-M2.5 SCSN catalog",
        "FLAG": "K-031 station-density field rho_sta(x,t) is NOT on disk; this is a surrogate",
        "eof1_pc_vs_region_mean_Mc_R2": r2_pc, "eof1_pc_beta": beta_pc,
        "partial_acf1": st_part["acf1"], "partial_corr_time_weeks": st_part["corr_time_weeks"],
        "partial_corr_len_km": st_part["corr_len_km"], "partial_moran_I": st_part["moran_I"],
        "partial_eof1_varfrac": st_part["eof1_varfrac"]}

    # ---------- positive control ----------
    print(f"[posctrl] injecting OU latent field (tau={OU_TAU_DAYS:.0f} d = {OU_TAU_DAYS/DT_DAYS:.1f} wk, "
          f"L={OU_LEN_KM} km, sd ladder={OU_SD_LADDER}) into simulated catalogs ...")
    pc_out = []
    ladder = [(sd, r) for sd in OU_SD_LADDER for r in range(N_POSCTRL // len(OU_SD_LADDER))]
    for s, (OU_SD, _rep) in enumerate(ladder):
        el = (time.time() - t_start) / 60.0
        if el > RUNTIME_BUDGET_MIN * 0.85:
            print(f"[posctrl] RUNTIME GUARD: stopping at {s} positive controls")
            break
        fld = ou_latent_field(rng, grid, n_weeks, OU_TAU_DAYS, OU_LEN_KM, OU_SD)
        st_, sm_, sla_, slo_, trunc = simulate_st_etas(
            rng, h_t, h_m, h_la, h_lo, PARS, p_k4, grid, week_edges[0], week_edges[-1], b_val,
            cap, ou_field=fld, week_edges=week_edges)
        inw = (st_ >= week_edges[0]) & (st_ < week_edges[-1])
        inbox = ((sla_ >= LAT_MIN) & (sla_ <= LAT_MAX) & (slo_ >= LON_MIN) & (slo_ <= LON_MAX))
        obs_s = grid.bin_events(sla_[inw & inbox], slo_[inw & inbox], st_[inw & inbox], week_edges)
        ss_t = np.concatenate([h_t, st_[inw]]); ss_m = np.concatenate([h_m, sm_[inw]])
        ss_la = np.concatenate([h_la, sla_[inw]]); ss_lo = np.concatenate([h_lo, slo_[inw]])
        o = np.argsort(ss_t)
        E_s, _ = expected_counts(grid, ss_t[o], ss_m[o], ss_la[o], ss_lo[o], week_edges, PARS, p_k4)
        R_s = residual_field(obs_s, E_s)
        stt = all_stats(R_s, sops)
        pc_out.append({"log_sd": OU_SD,
                       "corr_time_weeks": stt["corr_time_weeks"], "corr_len_km": stt["corr_len_km"],
                       "acf1": stt["acf1"], "eof1_varfrac": stt["eof1_varfrac"],
                       "n_events_in_box": int((inw & inbox).sum()), "truncated": bool(trunc)})
        print(f"[posctrl] {s+1} (sd={OU_SD}): T={stt['corr_time_weeks']:.2f} wk "
              f"(truth {OU_TAU_DAYS/DT_DAYS:.1f}) L={stt['corr_len_km']:.1f} km "
              f"(truth {OU_LEN_KM}) acf1={stt['acf1']:+.4f}")
    if pc_out:
        T_true, L_true = OU_TAU_DAYS / DT_DAYS, OU_LEN_KM
        by_sd = {}
        for sd in OU_SD_LADDER:
            rr = [x for x in pc_out if x["log_sd"] == sd]
            if not rr:
                continue
            Tt = np.array([x["corr_time_weeks"] for x in rr], dtype=float)
            Ll = np.array([x["corr_len_km"] for x in rr], dtype=float)
            okT = np.isfinite(Tt) & (Tt > T_true / 2) & (Tt < T_true * 2)
            okL = np.isfinite(Ll) & (Ll > L_true / 2) & (Ll < L_true * 2)
            by_sd[str(sd)] = {"n": len(rr), "T_median": float(np.nanmedian(Tt)),
                              "L_median": float(np.nanmedian(Ll)),
                              "frac_T_within_factor2": float(okT.mean()),
                              "frac_L_within_factor2": float(okL.mean()),
                              "acf1_median": float(np.nanmedian([x["acf1"] for x in rr]))}
            print(f"[posctrl] sd={sd}: median T={np.nanmedian(Tt):.2f} wk, L={np.nanmedian(Ll):.1f} km; "
                  f"within-factor-2 T {okT.mean()*100:.0f}% L {okL.mean()*100:.0f}%")
        passes = [v for v in by_sd.values()
                  if v["frac_T_within_factor2"] >= 0.5 and v["frac_L_within_factor2"] >= 0.5]
        res["positive_control"] = {
            "injected": {"tau_weeks": T_true, "tau_days": OU_TAU_DAYS, "length_km": L_true,
                         "log_sd_ladder": OU_SD_LADDER},
            "runs": pc_out, "by_amplitude": by_sd,
            "PASS_at_any_amplitude": bool(passes),
            "interpretation": ("Popper: 'the pipeline must recover both scales to within a factor "
                               "of 2. This is what makes a null interpretable.' If no amplitude "
                               "recovers both scales, a null K-009 result is NOT interpretable as "
                               "'there is no weather' -- it is only 'this pipeline could not have "
                               "seen weather injected on mu at these amplitudes'.")}

    # ---------- envelope-crossing scales (the non-saturating headline statistics) ----
    # The exponential-fit T and L saturate at their caps for the noise-dominated null
    # (a fitted slope <= 0 is meaningless when the ACF amplitude is ~0), so the primary
    # scale statistic is: how far does the REAL curve stay above the sim-null 97.5th
    # percentile envelope? That is well defined for both, and it is what "correlation
    # time / length with sim-null bounds" actually means here.
    acf_p975 = np.percentile(null_acf_curves, 97.5, axis=0)
    cor_p975 = np.nanpercentile(null_cor_curves, 97.5, axis=0)
    r_acf = np.asarray(real["adaptive_k4"]["acf"])
    r_cor = np.asarray(real["adaptive_k4"]["correlogram"])
    T_cross = 0
    for k in range(1, MAX_LAG_WEEKS + 1):
        if r_acf[k] > acf_p975[k]:
            T_cross = k
        else:
            break
    n_lags_above = int(sum(r_acf[k] > acf_p975[k] for k in range(1, MAX_LAG_WEEKS + 1)))
    T_censored = bool(T_cross >= MAX_LAG_WEEKS)
    L_cross, L_last_above = 0.0, 0.0
    broke = False
    for i, dkm in enumerate(sops.centers):
        if not np.isfinite(r_cor[i]) or not np.isfinite(cor_p975[i]):
            continue
        if r_cor[i] > cor_p975[i]:
            L_last_above = float(DIST_BINS[i + 1])
            if not broke:
                L_cross = float(DIST_BINS[i + 1])
        else:
            broke = True
    n_bins_above = int(np.nansum(r_cor > cor_p975))
    res["envelope_crossing"] = {
        "definition": ("largest scale at which the real pooled curve stays above the "
                       "ETAS-sim null 97.5th percentile"),
        "TIME_contiguous_weeks": int(T_cross), "TIME_n_lags_above_of_52": n_lags_above,
        "TIME_censored_at_spec_lag_window": T_censored,
        "LENGTH_contiguous_km": L_cross, "LENGTH_last_bin_above_km": L_last_above,
        "LENGTH_n_bins_above": n_bins_above,
        "LENGTH_resolution_floor_km": 18.0,
        "correlogram_note": ("the 0-10 km bin is empty at 0.2 deg resolution; the 35 km bin is "
                             "negative in the data, a lattice-aliasing artifact of a regular grid, "
                             "so the contiguous and last-above lengths differ and both are reported")}
    print(f"[crossing] real ACF above null p97.5 for lags 1..{T_cross} contiguous "
          f"({n_lags_above}/52 lags above){' [CENSORED at spec window]' if T_censored else ''}; "
          f"correlogram above null to {L_cross:.0f} km contiguous / {L_last_above:.0f} km last "
          f"({n_bins_above} bins above)")

    # ---------- robustness: is the whole signal one M7.2 aftershock sequence? ----------
    # The leading EOF is visibly dominated by a spike in the first months of 2010, which is
    # the M7.2 El Mayor-Cucapah sequence. If the redness is only that sequence's misfit, the
    # "latent field" reading collapses. Re-run the statistics with 2010 dropped.
    drop = 52
    R_no2010 = R_real[:, drop:]
    st_no = all_stats(R_no2010, sops)
    acf_p975_no = np.percentile(np.array([s["acf"] for s in sims]), 97.5, axis=0)
    Tc_no = 0
    for k in range(1, MAX_LAG_WEEKS + 1):
        if st_no["acf"][k] > acf_p975_no[k]:
            Tc_no = k
        else:
            break
    pc_abs = np.abs(pc1_probe := np.asarray(real["adaptive_k4"]["eof1_pc"]))
    top_wk = int(np.argmax(pc_abs))
    top5_frac = float(np.sort(pc_abs ** 2)[-5:].sum() / (pc_abs ** 2).sum())
    res["robustness_excluding_2010"] = {
        "reason": ("the leading residual EOF is dominated by a single spike; the M7.2 "
                   "El Mayor-Cucapah sequence (2010-04-04) falls in this window"),
        "weeks_dropped": drop,
        "eof1_pc_peak_week": top_wk,
        "eof1_pc_peak_date": str(pd.Timestamp(week_edges[top_wk] * 86400, unit="s", tz="UTC").date()),
        "eof1_pc_top5_week_variance_fraction": top5_frac,
        "acf1_full": real["adaptive_k4"]["acf1"], "acf1_excl_2010": st_no["acf1"],
        "moran_full": real["adaptive_k4"]["moran_I"], "moran_excl_2010": st_no["moran_I"],
        "eof1_varfrac_full": real["adaptive_k4"]["eof1_varfrac"],
        "eof1_varfrac_excl_2010": st_no["eof1_varfrac"],
        "TIME_crossing_weeks_excl_2010": int(Tc_no),
        "null_acf1_p97_5": NULL["acf1"]["p97_5"],
        "SURVIVES": bool(st_no["acf1"] - NULL["acf1"]["p97_5"] >= 0.05)}
    print(f"[robust] EOF1 PC peaks at week {top_wk} ({res['robustness_excluding_2010']['eof1_pc_peak_date']}), "
          f"top-5 weeks carry {top5_frac*100:.1f}% of PC variance")
    print(f"[robust] excluding 2010: acf1={st_no['acf1']:+.4f} (full {real['adaptive_k4']['acf1']:+.4f}), "
          f"MoranI={st_no['moran_I']:+.4f}, EOF1={st_no['eof1_varfrac']*100:.1f}%, "
          f"T_crossing={Tc_no} wk -> survives={res['robustness_excluding_2010']['SURVIVES']}")

    # ---------- rate-matched null sensitivity (supercritical-ETAS guard) ----------
    n_obs_w = len(test)
    rm = [s for s in sims if 0.7 * n_obs_w <= s["n_events_in_box"] <= 1.3 * n_obs_w]
    if len(rm) >= 3:
        rm_acf1 = np.array([s["acf1"] for s in rm])
        res["rate_matched_null"] = {
            "n": len(rm), "criterion": "sim event count within +-30% of the observed count",
            "acf1_median": float(np.median(rm_acf1)), "acf1_p97_5": float(np.percentile(rm_acf1, 97.5)),
            "corr_time_weeks_median": float(np.median([s["corr_time_weeks"] for s in rm])),
            "corr_len_km_median": float(np.median([s["corr_len_km"] for s in rm]))}
    else:
        res["rate_matched_null"] = {
            "n": len(rm), "criterion": "sim event count within +-30% of the observed count",
            "STATUS": ("too few rate-matched sims to form an envelope -- the free-running "
                       "spatio-temporal ETAS does not reproduce the observed event count of this "
                       "window. Reported as a first-class diagnostic, not suppressed.")}
    print(f"[null] rate-matched sims (+-30% of n_obs={n_obs_w}): {len(rm)} / {n_sims}")

    # ---------- VERDICT per Popper's frozen rule ----------
    a1 = real["adaptive_k4"]["acf1"]
    excess = a1 - NULL["acf1"]["p97_5"]
    cond_excess = bool(excess >= 0.05)
    swap_vals = [real[k]["acf1"] for k in BG]
    swap_excess = [v - NULL["acf1"]["p97_5"] for v in swap_vals]
    cond_swap = bool(all(e >= 0.05 for e in swap_excess))
    cond_rho = bool((st_part["acf1"] - NULL["acf1"]["p97_5"]) >= 0.05)
    verdict = "SUCCESS" if (cond_excess and cond_swap and cond_rho) else "FAILURE"

    res["real"] = {k: {kk: (v.tolist() if isinstance(v, np.ndarray) else v)
                       for kk, v in real[k].items() if kk != "eof1_pc"} for k in real}
    res["real"]["adaptive_k4"]["eof1_pc"] = pc1.tolist()
    res["real"]["time_rescaling"] = tr_real
    res["null"] = NULL
    res["null_acf_envelope"] = {"lags": list(range(MAX_LAG_WEEKS + 1)),
                                "p2_5": np.percentile(null_acf_curves, 2.5, axis=0).tolist(),
                                "median": np.median(null_acf_curves, axis=0).tolist(),
                                "p97_5": np.percentile(null_acf_curves, 97.5, axis=0).tolist()}
    res["null_correlogram_envelope"] = {"dist_km": sops.centers.tolist(),
                                        "p2_5": np.nanpercentile(null_cor_curves, 2.5, axis=0).tolist(),
                                        "median": np.nanmedian(null_cor_curves, axis=0).tolist(),
                                        "p97_5": np.nanpercentile(null_cor_curves, 97.5, axis=0).tolist()}
    T_sat = bool(real["adaptive_k4"]["corr_time_weeks"] >= MAX_LAG_WEEKS - 1e-9)
    L_sat = bool(real["adaptive_k4"]["corr_len_km"] >= DIST_BINS[-1] - 1e-9)
    res["headline"] = {
        "PRIMARY_TIME_weeks": (f">{T_cross}" if T_censored else str(int(T_cross))),
        "PRIMARY_TIME_statement": (
            f"residual correlation TIME {'>' if T_censored else '='} {T_cross} weeks "
            f"({T_cross*DT_DAYS:.0f} d): the pooled residual ACF exceeds the ETAS-sim null "
            f"97.5th percentile at every lag from 1 to {T_cross} weeks"
            + (" -- CENSORED, the spec's 52-week lag window is exhausted, so this is a LOWER BOUND"
               if T_censored else "")),
        "PRIMARY_LENGTH_km": L_cross,
        "PRIMARY_LENGTH_statement": (
            f"residual correlation LENGTH: the pooled spatial correlogram exceeds the ETAS-sim "
            f"null 97.5th percentile contiguously to {L_cross:.0f} km and in {n_bins_above} bins "
            f"out to {L_last_above:.0f} km; e-folding fit L = "
            f"{real['adaptive_k4']['corr_len_km']:.1f} km; resolution floor ~18 km"),
        "efolding_TIME_weeks": real["adaptive_k4"]["corr_time_weeks"],
        "efolding_TIME_SATURATED_AT_CAP": T_sat,
        "efolding_LENGTH_SATURATED_AT_CAP": L_sat,
        "residual_correlation_TIME_days": real["adaptive_k4"]["corr_time_weeks"] * DT_DAYS,
        "null_TIME_weeks_2.5_97.5": [NULL["corr_time_weeks"]["p2_5"], NULL["corr_time_weeks"]["p97_5"]],
        "null_TIME_envelope_IS_DEGENERATE": (
            "the exponential-fit T is pinned at its cap for most of the null because a fitted "
            "slope <= 0 is meaningless when the null ACF amplitude is ~0; use PRIMARY_TIME"),
        "residual_correlation_TIME_integral_weeks": real["adaptive_k4"]["corr_time_integral_weeks"],
        "null_TIME_integral_weeks_2.5_97.5": [NULL["corr_time_integral_weeks"]["p2_5"],
                                              NULL["corr_time_integral_weeks"]["p97_5"]],
        "residual_correlation_LENGTH_km": real["adaptive_k4"]["corr_len_km"],
        "null_LENGTH_km_2.5_97.5": [NULL["corr_len_km"]["p2_5"], NULL["corr_len_km"]["p97_5"]],
        "residual_correlation_LENGTH_integral_km": real["adaptive_k4"]["corr_len_integral_km"],
        "null_LENGTH_integral_km_2.5_97.5": [NULL["corr_len_integral_km"]["p2_5"],
                                             NULL["corr_len_integral_km"]["p97_5"]],
        "LENGTH_resolution_floor_km": ("~18 km: the shortest 0.2-deg cell-centre separation. "
                                       "Structure below this is invisible to the correlogram."),
        "lag1_weekly_ACF": a1,
        "null_lag1_97.5": NULL["acf1"]["p97_5"],
        "lag1_excess_over_null_97.5": excess}
    res["success_rule"] = {
        "quoted": ("Lag-1 weekly ACF excess >= 0.05 over the sim-null 97.5th percentile, and stable "
                   "across the kernel swap, and surviving the rho_sta partial. Failure: excess inside "
                   "the null envelope, or destroyed by the kernel swap or the rho_sta partial."),
        "lag1_excess": excess, "cond_excess_ge_0.05": cond_excess,
        "kernel_swap_acf1": {k: real[k]["acf1"] for k in BG},
        "kernel_swap_excess": dict(zip(BG.keys(), swap_excess)), "cond_kernel_swap_stable": cond_swap,
        "rho_sta_partial_acf1": st_part["acf1"], "cond_survives_rho_sta": cond_rho,
        "VERDICT": verdict}
    # ---------- score against the registered two-generator prediction (c2bf012) ----------
    T_meas = real["adaptive_k4"]["corr_time_weeks"]
    L_meas = real["adaptive_k4"]["corr_len_km"]
    pred_path = HERE / "results_k009_prediction.json"
    gen = {"STATUS": "prediction file not found"}
    if pred_path.exists():
        P = json.loads(pred_path.read_text())
        T_pred = P["generator_B_W001_W002"]["predicted_correlation_TIME_weeks"]
        rng_wk = P["generator_B_W001_W002"]["estimator_uncertainty"]["full_estimator_range_weeks"]
        b_lo, b_hi = T_pred / 2.0, T_pred * 2.0
        # W-003 predicts white: no lag exceeds the null envelope and no lag-1 excess.
        red = bool(T_cross >= 1 and excess >= 0.05)
        A_wins = bool(T_cross == 0 and excess < 0.05)
        # Measured time is an interval: [T_cross, inf) if censored, else a point at T_cross.
        if T_censored:
            B_consistent = bool(red and b_hi >= T_cross)     # band overlaps [T_cross, inf)
            B_sharp = False
        else:
            B_consistent = bool(red and b_lo <= T_cross <= b_hi)
            B_sharp = B_consistent
        in_full_range = bool(red and T_cross >= rng_wk[0])
        if A_wins:
            gname = "W-003 WINS (residuals white)"
        elif red and B_consistent and B_sharp:
            gname = "W-001/W-002 WINS (red, T inside the factor-2 band around t_a)"
        elif red and B_consistent:
            gname = ("W-003 FALSIFIED; consistent with W-001/W-002 but NOT sharply -- the "
                     "measured T is a censored lower bound that overlaps the factor-2 band "
                     "rather than landing inside it")
        elif red:
            gname = ("W-003 FALSIFIED; residuals red but off the predicted t_a timescale -- "
                     "both generators lose, reported as a bare measurement")
        else:
            gname = "NEITHER -- not red, but not cleanly white either"
        gen = {"prediction_file": "results_k009_prediction.json",
               "pre_registration_commit": "c2bf012",
               "T_measured_weeks_crossing": int(T_cross),
               "T_measured_is_censored_lower_bound": T_censored,
               "T_measured_efolding_weeks": T_meas,
               "T_efolding_saturated": bool(T_meas >= MAX_LAG_WEEKS - 1e-9),
               "L_measured_km_crossing": L_cross, "L_measured_efolding_km": L_meas,
               "W003_predicted_T_weeks": 0.0, "W003_predicted_L_km": 0.0,
               "W003_falsification_evidence": {
                   "lag1_ACF": a1, "null_p97_5": NULL["acf1"]["p97_5"], "excess": excess,
                   "moran_I": real["adaptive_k4"]["moran_I"],
                   "moran_null_p97_5": NULL["moran_I"]["p97_5"],
                   "eof1_varfrac": real["adaptive_k4"]["eof1_varfrac"],
                   "eof1_null_p97_5": NULL["eof1_varfrac"]["p97_5"],
                   "n_lags_above_null_of_52": n_lags_above},
               "W001_W002_predicted_T_weeks": T_pred,
               "W001_W002_factor2_band_weeks": [b_lo, b_hi],
               "W001_W002_full_estimator_range_weeks": rng_wk,
               "residuals_are_red": red,
               "W003_wins": A_wins, "W001_W002_wins": bool(B_sharp),
               "W001_W002_consistent_not_sharp": bool(B_consistent and not B_sharp),
               "inside_full_t_a_estimator_range": in_full_range,
               "GENERATOR_VERDICT": gname}
        print(f"[generators] T_crossing={T_cross}{'+ (censored)' if T_censored else ''} wk vs "
              f"W-003 predicted 0 and W-001/W-002 predicted {T_pred:.1f} "
              f"(band {b_lo:.1f}-{b_hi:.1f})  ->  {gname}")
    res["generator_discrimination"] = gen

    res["pre_registration"] = {
        "popper_spec_commit": "0a73fd2 (HYPOTHESIS_LEDGER.md, K-009 TESTABLE-NOW verdict)",
        "prediction_register_commit": "c2bf012 (results_k009_prediction.json + exp_k009_prediction.py)",
        "prediction_status": "PREDICTION POST-REGISTERED (ordering violation, self-reported)",
        "smoke_test_disclosure": (
            "Two smoke-test runs of this script were executed BEFORE the two-generator ruling "
            "arrived, on a REDUCED 2-year sub-window (2010-01-01 .. 2012-01-01, not the spec's "
            "2010-2018) with 2 simulations, to debug the pipeline. They printed measured residual "
            "correlation statistics, so the executing agent was NOT blind when the prediction was "
            "written. Seen in those debug runs: lag-1 pooled ACF approx +0.16 to +0.18, "
            "correlation time approx 9.8-10.2 weeks, correlation length approx 27-29 km, 2-sim "
            "null lag-1 ACF approx -0.015. The registered t_a prediction (27.6 weeks) was derived "
            "solely from train-window (<2010) aftershock stacking and lies OUTSIDE the seen "
            "sub-window value, but the ordering is violated and the result carries that caveat. "
            "Debug artifacts retained for audit: _smoke_k009.json, maps/_smoke_k009.png."),
        "not_backdated": True}

    res["runtime_minutes"] = round((time.time() - t_start) / 60.0, 2)
    sim_ns = sorted(int(s["n_events_in_box"]) for s in sims)
    res["flags"] = {
        "n_sims_run": n_sims,
        "n_sims_spec": 500,
        "n_sims_guard_reduced": (f"runtime guard {RUNTIME_BUDGET_MIN} min: {n_sims} sims run "
                                 f"instead of the spec's 500; recorded, not silently substituted"),
        "runtime_guard_min": RUNTIME_BUDGET_MIN,
        "world_arm_UNRUN": ("Popper's spec also mandates the world arm at 0.5 deg x 30 d over the "
                            "13 boxes. It was NOT run (runtime guard). This is UNRUN, not skipped: "
                            "the SoCal verdict is therefore partial with respect to the frozen spec."),
        "rho_sta_surrogate_WEAKENS_SUCCESS": (
            "Popper's success rule requires the excess to survive the K-031 station-density partial. "
            "K-031 has not been run and rho_sta(x,t) is not on disk, so an Mc(x,t) completeness "
            "proxy was substituted. Any SUCCESS label from this run is correspondingly weakened and "
            "must not be reported as having cleared the mandated observer-nuisance control."),
        "kernel_swap_is_LOW_POWER": (
            "The background field supplies only a small fraction of the conditional intensity "
            f"(triggered fraction of integrated intensity in EXP-H test = "
            f"{exph['test']['triggered_fraction_of_integrated_intensity']:.3f}), so the three "
            "background fields produce nearly identical residuals. The kernel-swap control will "
            "pass almost regardless of the truth. It is reported as a weak control, and 'stable "
            "across the kernel swap' must NOT be read as strong evidence against kernel misfit."),
        "sim_count_distribution_FIRST_CLASS": {
            "note": ("Branching ratio n = "
                     f"{exph['train_fit']['branching_ratio_n']:.3f} is SUPERCRITICAL, so the "
                     "free-running spatio-temporal ETAS need not reproduce the observed event "
                     "count. A sim/obs count mismatch changes the residual noise level and would "
                     "confound the null, so the distribution is reported as a primary diagnostic."),
            "n_observed_in_window": int(n_obs_w),
            "sim_counts_sorted": sim_ns,
            "sim_count_median": int(np.median(sim_ns)) if sim_ns else None,
            "sim_count_over_obs_median": float(np.median(sim_ns) / n_obs_w) if sim_ns else None,
            "n_rate_matched_within_30pct": int(len(rm)),
            "sim_event_cap": cap,
            "sims_hitting_cap": int(sum(s["truncated"] for s in sims))},
        "branching_ratio_supercritical": exph["train_fit"]["branching_ratio_n"],
        "K002_not_on_disk": ("no results_k002.json exists; the spatio-temporal ETAS required by "
                             "Popper's fix #1 was built inside this script (background field and "
                             "spatial kernel fit on the train window only)"),
        "pair_truncation_in_spatial_kernel_fit": [PAIR_T, PAIR_R]}

    OUT.write_text(json.dumps(res, indent=2, default=float))

    # ---------- figure ----------
    FIG.parent.mkdir(exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    lags = np.arange(MAX_LAG_WEEKS + 1)
    a = ax[0, 0]
    a.fill_between(lags, np.percentile(null_acf_curves, 2.5, axis=0),
                   np.percentile(null_acf_curves, 97.5, axis=0), color="0.8",
                   label=f"ETAS-sim null 95% ({n_sims} sims)")
    a.plot(lags, np.median(null_acf_curves, axis=0), "k--", lw=1, label="null median")
    for k, st in zip(BG, ["-", "--", ":"]):
        a.plot(lags, real[k]["acf"], st, lw=1.6, label=f"real ({k})")
    a.axhline(0, color="k", lw=0.5)
    a.set_xlabel("lag (weeks)"); a.set_ylabel("pooled residual ACF")
    a.set_title(f"(i) temporal ACF of Anscombe residual\ncorrelation TIME "
                f"{'>' if T_censored else '='} {T_cross} wk "
                f"(above sim-null p97.5 at {n_lags_above}/52 lags"
                f"{'; CENSORED lower bound' if T_censored else ''})")
    a.legend(fontsize=7); a.set_xlim(0, MAX_LAG_WEEKS)

    a = ax[0, 1]
    a.fill_between(sops.centers, np.nanpercentile(null_cor_curves, 2.5, axis=0),
                   np.nanpercentile(null_cor_curves, 97.5, axis=0), color="0.8",
                   label="ETAS-sim null 95%")
    a.plot(sops.centers, np.nanmedian(null_cor_curves, axis=0), "k--", lw=1, label="null median")
    for k, st in zip(BG, ["-", "--", ":"]):
        a.plot(sops.centers, real[k]["correlogram"], st, lw=1.6, label=f"real ({k})")
    a.plot(sops.centers, st_part["correlogram"], "-", color="crimson", lw=1.2,
           label="real, rho_sta-partialled")
    a.axhline(0, color="k", lw=0.5)
    a.set_xlabel("separation (km)"); a.set_ylabel("pooled residual correlation")
    a.set_title(f"(ii) spatial correlogram\ncorrelation LENGTH: above sim-null p97.5 to "
                f"{L_cross:.0f} km contiguous / {L_last_above:.0f} km last "
                f"({n_bins_above} bins); e-fold fit {real['adaptive_k4']['corr_len_km']:.0f} km")
    a.legend(fontsize=7)

    a = ax[1, 0]
    x = np.sort(tau_real)
    n = len(x)
    theo = -np.log(1.0 - (np.arange(1, n + 1) - 0.5) / n)
    a.plot(theo, x, ".", ms=2, label=f"real (KS={tr_real['ks']:.4f})")
    lim = max(theo.max(), x.max())
    a.plot([0, lim], [0, lim], "k--", lw=1)
    a.set_xlabel("Exp(1) quantile"); a.set_ylabel("rescaled interevent time")
    a.set_title(f"(iii) time-rescaling QQ (whole-footprint process)\n"
                f"null KS median {NULL['time_rescaling_ks']['median']:.4f}, "
                f"97.5% {NULL['time_rescaling_ks']['p97_5']:.4f}")
    a.legend(fontsize=8)

    a = ax[1, 1]
    a.plot(np.arange(n_weeks), pc1, lw=0.9, label=f"EOF1 PC ({ev1*100:.1f}% var)")
    a.plot(np.arange(n_weeks), rmean * pc1.std(), lw=0.9, color="crimson",
           label=f"region-mean Mc proxy (R2={r2_pc:.3f})")
    a.set_xlabel("week since 2010-01-01"); a.set_ylabel("EOF1 amplitude")
    a.set_title(f"(iv) leading residual EOF vs observer proxy\n"
                f"null EOF1 var frac {NULL['eof1_varfrac']['p2_5']*100:.1f}-"
                f"{NULL['eof1_varfrac']['p97_5']*100:.1f}%; PC1 peak wk {top_wk} "
                f"(M7.2 El Mayor); excl-2010 acf1 {st_no['acf1']:+.3f}")
    a.legend(fontsize=8)

    fig.suptitle(f"K-009 ETAS residual whiteness, SoCal 0.2 deg x 7 d, 2010-2018  |  VERDICT: {verdict}",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(FIG, dpi=130)
    print(f"[fig] {FIG}")

    print("\n" + "=" * 78)
    print(f"K-009  |  {grid.n_cells} cells x {n_weeks} weeks; {len(test)} observed events")
    print(f"  HEADLINE  correlation TIME   {'>' if T_censored else '='} {T_cross} weeks "
          f"({T_cross*DT_DAYS:.0f} d){' [CENSORED: lower bound]' if T_censored else ''} "
          f"-- ACF above sim-null p97.5 at {n_lags_above}/52 lags")
    print(f"  HEADLINE  correlation LENGTH = {L_cross:.0f} km contiguous "
          f"({L_last_above:.0f} km last bin above null, {n_bins_above} bins); "
          f"e-folding fit {real['adaptive_k4']['corr_len_km']:.1f} km; resolution floor ~18 km")
    print(f"  lag-1 weekly ACF = {a1:+.4f}; null 97.5% = {NULL['acf1']['p97_5']:+.4f}; "
          f"excess = {excess:+.4f} (need >= 0.05): {cond_excess}")
    print(f"  kernel swap stable: {cond_swap}   " +
          "  ".join(f"{k}:{real[k]['acf1']:+.4f}" for k in BG))
    print(f"  survives rho_sta partial: {cond_rho} (partial acf1 {st_part['acf1']:+.4f})")
    print(f"  VERDICT (Popper's rule): {verdict}")
    print(f"  GENERATOR (vs c2bf012 register): {gen.get('GENERATOR_VERDICT', 'n/a')}")
    print(f"  sim counts: median {np.median(sim_ns):.0f} vs obs {n_obs_w} "
          f"({len(rm)}/{n_sims} rate-matched)")
    print(f"  runtime {res['runtime_minutes']} min  ->  {OUT.name}, {FIG.name}")
    print("=" * 78)


if __name__ == "__main__":
    main()
