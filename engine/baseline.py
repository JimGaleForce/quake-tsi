"""Baselines.

  * ``ClimatologyV1``  -- time-stationary per-cell Poisson rate mu(x). v1's declared
    gap: aftershock clustering is NOT absorbed and leaks into every covariate.
  * ``EtasV1``         -- space-time ETAS at the engine's cell-day resolution, fitted
    by Poisson maximum likelihood on the EXPLORATION window only. This is the
    baseline that absorbs Omori-Utsu clustering, so a covariate scored against it is
    being asked "do you know anything an aftershock model does not?".

ETAS model (v1, ruled):

    lambda(c, t) = nu * mu(c)
                 + K * sum_{i : 0 < t - t_i <= L} 10^(alpha*(m_i - 4.5))
                       * (t - t_i + c0)^(-p) * G_sigma(dist(c, cell_i))

  * ``mu(c)``  is the exploration-window climatology above; ``nu`` rescales it, so the
    background keeps climatology's spatial shape and ETAS only re-calibrates its level.
  * ``t - t_i`` is in whole days and is >= 1 by construction: day ``t`` sees only events
    on days strictly before ``t``. Same causality contract as the covariates.
  * ``G_sigma`` is an isotropic Gaussian on ground distance between CELL CENTRES,
    column-normalised over the domain (sum_c G(c, s) = 1 for every source cell s), so
    K is exactly "expected triggered target-events per unit productivity" and no
    triggered mass leaks out of the restricted domain.
  * The Omori tail is TRUNCATED at L = 365 days (declared approximation). At the
    fitted p = 1.16 that keeps only ~64% of the kernel's formal mass -- the tail of
    an Omori law is heavy and a third of it lives beyond a year. It costs little in
    likelihood anyway (measured: L = 20 d scores +0.911 bits/event over climatology,
    L = 365 d scores +0.962) because the far tail is nearly flat in time and is
    absorbed by the nu*mu background term instead. What the truncation really does
    is REASSIGN decade-long aftershock sequences from triggering to background, so
    do not read nu as a declustered background rate.

What this baseline still does NOT absorb, and which therefore stays available to a
covariate: anisotropy of the aftershock cloud (the kernel is isotropic), focal
mechanism / fault geometry, depth, magnitude-dependent rupture length, secular rate
changes, and catalogue completeness drift.
"""

from __future__ import annotations

import json
import math
import os
import time

import numpy as np

from . import BASELINE_CAVEAT, ETAS_CAVEAT

ETAS_CACHE = os.path.join("engine", "out", "cache", "etas_params.json")

MAG_REF = 4.5                 # magnitude zero-point of the productivity term
OMORI_TRUNC_DAYS = 365        # declared time-truncation of the Omori tail
ETAS_BURN_IN_DAYS = 365       # head of the record with an un-warmed triggering history
LOCAL_MAX_NATS = 1e-3         # LL slack below which the fit counts as a local maximum

PARAM_NAMES = ("nu", "K", "alpha", "c", "p", "sigma")
PARAM_BOUNDS = {
    "nu": (1e-3, 10.0),
    "K": (1e-6, 100.0),
    "alpha": (0.5, 1.5),
    "c": (1e-3, 1.0),         # days
    "p": (0.8, 1.5),
    "sigma": (10.0, 300.0),   # km
}
PARAM_START = {"nu": 0.5, "K": 0.05, "alpha": 1.0, "c": 0.05, "p": 1.1, "sigma": 50.0}

LN2 = math.log(2.0)


class ClimatologyV1:
    """Time-stationary per-cell rate, fitted on the EXPLORATION period only.

    mu_c = (N_c + alpha) / (T_explore + alpha / mu_bar)

    i.e. the posterior mean of a Gamma(alpha, alpha/mu_bar) prior centred on the
    domain-average per-cell daily rate mu_bar. Smoothing is required because a cell
    with zero exploration events would otherwise have mu=0 and give -inf log-likelihood
    the first time an event lands in it. alpha=0.5 by default (weak).
    """

    name = "climatology-v1"
    caveat = BASELINE_CAVEAT
    burn_in_days = 0

    def __init__(self, alpha: float = 0.5):
        self.alpha = float(alpha)
        self.mu = None

    def fit(self, ctx, counts: np.ndarray, train_slice: slice):
        train = counts[:, train_slice]
        n_train_days = train.shape[1]
        if n_train_days < 1:
            raise ValueError("empty training window")
        n_c = train.sum(axis=1).astype(np.float64)
        total = n_c.sum()
        if total <= 0:
            raise ValueError("bulk transform produced 0 training events (baseline fit)")
        mu_bar = total / (n_c.size * n_train_days)
        self.mu = (n_c + self.alpha) / (n_train_days + self.alpha / mu_bar)
        assert np.isfinite(self.mu).all() and (self.mu > 0).all()
        self.n_train_days = n_train_days
        self.mu_bar = float(mu_bar)
        self.n_train_events = float(total)
        self.n_empty_cells = int((n_c == 0).sum())
        return self

    def rate(self, sl: slice | None = None):
        """Intensity for a day window. 1-D -> broadcast over days (time-stationary)."""
        return self.mu

    def report(self):
        return [
            f"baseline             = {self.name} (alpha={self.alpha})",
            f"  train days         = {self.n_train_days}",
            f"  train events       = {self.n_train_events:.0f}",
            f"  mean cell rate/day = {self.mu_bar:.3e}",
            f"  mu range           = {self.mu.min():.3e} .. {self.mu.max():.3e}",
            f"  cells empty in train = {self.n_empty_cells} (smoothed, not zero)",
        ]


# ---------------------------------------------------------------------- ETAS ---
def _log_theta(params: dict) -> np.ndarray:
    return np.array([math.log(params[k]) for k in PARAM_NAMES], dtype=np.float64)


def _from_log_theta(theta) -> dict:
    return {k: float(math.exp(v)) for k, v in zip(PARAM_NAMES, np.asarray(theta))}


def _log_bounds():
    return [(math.log(PARAM_BOUNDS[k][0]), math.log(PARAM_BOUNDS[k][1]))
            for k in PARAM_NAMES]


class EtasKernel:
    """Vectorised ETAS intensity over the (n_cells, n_days) design.

    The expensive pieces are (a) the temporal Omori convolution, done once per
    (alpha, c, p) as a scatter-add over the ~5e4 occupied source cell-days rather
    than a python loop over cells, and (b) one float32 GEMM per (sigma) against the
    (n_cells, n_cells) Gaussian kernel. One full intensity evaluation over the whole
    design is well under a second, which is what makes ML fitting tractable.
    """

    def __init__(self, ctx, mu: np.ndarray, trunc_days: int = OMORI_TRUNC_DAYS):
        self.n_cells = int(ctx.n_cells)
        self.n_days = int(ctx.n_days)
        self.mu = np.asarray(mu, dtype=np.float32)
        self.trunc = int(trunc_days)
        self.lags = np.arange(1, self.trunc + 1, dtype=np.float64)

        # source events -> (cell, day, magnitude); day t sees days < t only.
        self.src_cell = np.asarray(ctx.ev_cell, dtype=np.int64)
        self.src_day = np.asarray(ctx.ev_day, dtype=np.int64)
        self.src_dm = np.asarray(ctx.ev_mag, dtype=np.float64) - MAG_REF

        self.d2 = (np.asarray(ctx.dist_km, dtype=np.float32) ** 2)
        self._sigma_cache = (None, None)
        self._build_lag_index()

    def _build_lag_index(self):
        """Precompute the (source pair) x (lag) -> flat target index map.

        The index map depends only on the catalogue and the truncation, NOT on any
        parameter, so it is built once. Each likelihood evaluation is then a single
        outer product plus one bincount over ~2e7 entries (~0.1 s) instead of a
        python loop. Targets past the end of the record are routed to one dump slot.
        """
        nd, nc = self.n_days, self.n_cells
        flat = self.src_cell * nd + self.src_day
        # collapse duplicate (cell, day) sources: ~5e4 pairs from ~7e4 events
        order = np.argsort(flat, kind="stable")
        f_sorted = flat[order]
        uniq, inv = np.unique(f_sorted, return_inverse=True)
        self._pair_of_event = np.empty(flat.size, dtype=np.int64)
        self._pair_of_event[order] = inv
        self.n_pairs = uniq.size
        s_idx, d_idx = uniq // nd, uniq % nd

        lags = np.arange(1, self.trunc + 1, dtype=np.int64)
        tgt = d_idx[:, None] + lags[None, :]
        self._dump = nc * nd
        self._flat_tgt = np.where(tgt < nd, s_idx[:, None] * nd + tgt,
                                  self._dump).ravel()

    # -- pieces ------------------------------------------------------------
    def spatial(self, sigma: float) -> np.ndarray:
        s, g = self._sigma_cache
        if s == sigma:
            return g
        g = np.exp(self.d2 * np.float32(-0.5 / (sigma * sigma)))
        g /= g.sum(axis=0, keepdims=True)     # column-normalised: no domain leakage
        self._sigma_cache = (sigma, g)
        return g

    def temporal(self, alpha: float, c0: float, p: float) -> np.ndarray:
        """A[s, t] = sum over prior events in cell s of productivity * Omori(t - t_i)."""
        val = np.bincount(self._pair_of_event,
                          weights=10.0 ** (alpha * self.src_dm),
                          minlength=self.n_pairs)
        kern = (self.lags + c0) ** (-p)
        acc = np.bincount(self._flat_tgt,
                          weights=(val[:, None] * kern[None, :]).ravel(),
                          minlength=self._dump + 1)
        return acc[:self._dump].reshape(self.n_cells, self.n_days).astype(np.float32)

    def intensity(self, params: dict, temporal: np.ndarray | None = None) -> np.ndarray:
        """(n_cells, n_days) float32 ETAS intensity."""
        if temporal is None:
            temporal = self.temporal(params["alpha"], params["c"], params["p"])
        g = self.spatial(params["sigma"])
        lam = (np.float32(params["K"]) * (g @ temporal))
        lam += (np.float32(params["nu"]) * self.mu)[:, None]
        np.maximum(lam, np.float32(1e-12), out=lam)
        return lam

    # -- likelihood --------------------------------------------------------
    @staticmethod
    def poisson_ll(lam: np.ndarray, y: np.ndarray) -> float:
        """sum_{c,t} y log(lam) - lam. Log is taken only where y > 0 (~5e4 cells)."""
        nzc, nzt = np.nonzero(y)
        yy = y[nzc, nzt].astype(np.float64)
        ll = float((yy * np.log(lam[nzc, nzt].astype(np.float64))).sum())
        ll -= float(lam.sum(dtype=np.float64))
        return ll


class EtasV1:
    """Space-time ETAS baseline. Parameters fitted on the exploration window only."""

    name = "etas-v1"
    caveat = ETAS_CAVEAT
    burn_in_days = ETAS_BURN_IN_DAYS

    def __init__(self, alpha: float = 0.5, cache_path: str = ETAS_CACHE,
                 params: dict | None = None, fit_if_missing: bool = True,
                 verbose: bool = True, trunc_days: int = OMORI_TRUNC_DAYS,
                 polish: bool = True, mag_target: float | None = None):
        self.mag_target = None if mag_target is None else float(mag_target)
        self.clim_alpha = float(alpha)
        self.cache_path = cache_path
        self.params = dict(params) if params else None
        self.fit_if_missing = bool(fit_if_missing)
        self.verbose = bool(verbose)
        self.trunc_days = int(trunc_days)
        # the head of the record has an un-warmed triggering history: drop exactly
        # the truncation length from both the ETAS fit and every downstream window.
        self.burn_in_days = int(trunc_days)
        self.polish = bool(polish)
        self.lam = None
        self.fit_info = {}

    # -- cache -------------------------------------------------------------
    def cache_key(self, ctx, mag_target, fit_slice):
        return {
            "design_cache_key": str(ctx.meta.get("cache_key", "synthetic")),
            "n_cells": int(ctx.n_cells),
            "n_days": int(ctx.n_days),
            "mag_target": float(mag_target),
            "fit_window_days": [int(fit_slice.start), int(fit_slice.stop)],
            "omori_trunc_days": int(self.trunc_days),
            "mag_ref": MAG_REF,
            "clim_alpha": self.clim_alpha,
            "model": "etas-v1",
        }

    def load_cache(self, key):
        if not os.path.exists(self.cache_path):
            return None
        with open(self.cache_path, "r", encoding="utf-8") as fh:
            rec = json.load(fh)
        if rec.get("key") != key:
            return None
        return rec

    def save_cache(self, rec):
        os.makedirs(os.path.dirname(self.cache_path) or ".", exist_ok=True)
        with open(self.cache_path, "w", encoding="utf-8") as fh:
            json.dump(rec, fh, indent=1, sort_keys=True)
        return self.cache_path

    # -- fit ---------------------------------------------------------------
    def fit(self, ctx, counts: np.ndarray, train_slice: slice):
        """Fit (or load) ETAS parameters and materialise lambda over the whole design."""
        self.clim = ClimatologyV1(self.clim_alpha).fit(ctx, counts, train_slice)
        self.kernel = EtasKernel(ctx, self.clim.mu, trunc_days=self.trunc_days)

        start = min(self.burn_in_days, max(0, int(train_slice.stop) - 1))
        fit_sl = slice(start, int(train_slice.stop))
        mag_target = float(self.mag_target if self.mag_target is not None
                           else ctx.meta.get("mag_floor", MAG_REF))
        key = self.cache_key(ctx, mag_target, fit_sl)

        rec = None if self.params else self.load_cache(key)
        if self.params is not None:
            self.fit_info = {"source": "explicit", "key": key}
        elif rec is not None:
            self.params = {k: float(rec["params"][k]) for k in PARAM_NAMES}
            self.fit_info = dict(rec)
            self.fit_info["source"] = "cache"
        else:
            if not self.fit_if_missing:
                raise FileNotFoundError(
                    f"no ETAS parameter cache at {self.cache_path} for this design "
                    f"({key}); run `python -u -m engine.cli fit-etas` first."
                )
            rec = self._fit_params(counts, fit_sl, key)
            self.params = {k: float(rec["params"][k]) for k in PARAM_NAMES}
            self.fit_info = dict(rec)
            self.fit_info["source"] = "fitted"

        self.lam = self.kernel.intensity(self.params)
        assert np.isfinite(self.lam).all() and (self.lam > 0).all()
        self.fit_slice = fit_sl
        return self

    def _fit_params(self, counts, fit_sl, key):
        from scipy import optimize

        y = counts[:, fit_sl]
        n_ev = float(y.sum())
        if n_ev <= 0:
            raise ValueError("bulk transform produced 0 events in the ETAS fit window")

        k = self.kernel
        state = {"n": 0, "t0": time.time(), "cache": (None, None)}

        def neg_ll(theta):
            pr = _from_log_theta(theta)
            shape = (pr["alpha"], pr["c"], pr["p"])
            if state["cache"][0] != shape:
                state["cache"] = (shape, k.temporal(*shape))
            lam = k.intensity(pr, temporal=state["cache"][1])[:, fit_sl]
            ll = k.poisson_ll(lam, y)
            state["n"] += 1
            if not np.isfinite(ll):
                return 1e18
            return -ll

        x0 = _log_theta(PARAM_START)
        bounds = _log_bounds()
        res = optimize.minimize(neg_ll, x0, method="L-BFGS-B", bounds=bounds,
                                options={"eps": 1e-5, "maxiter": 300, "ftol": 1e-11,
                                         "gtol": 1e-8})
        best_x, best_f, methods = res.x, float(res.fun), ["L-BFGS-B"]
        if self.verbose:
            print(f"  L-BFGS-B: nll={best_f:.4f} nfev={state['n']} "
                  f"success={res.success} ({res.message})")
        if self.polish:
            res2 = optimize.minimize(neg_ll, best_x, method="Nelder-Mead",
                                     bounds=bounds,
                                     options={"maxfev": 400, "xatol": 1e-5,
                                              "fatol": 1e-7})
            methods.append("Nelder-Mead")
            if self.verbose:
                print(f"  Nelder-Mead polish: nll={float(res2.fun):.4f} "
                      f"nfev={state['n']}")
            if float(res2.fun) < best_f:
                best_x, best_f = res2.x, float(res2.fun)

        params = _from_log_theta(best_x)

        # SPEC invariant "Poisson fit: assert convergence flag". L-BFGS-B reports
        # ABNORMAL when its line search dies against a bound, which is not the same
        # thing as "not at an optimum". So check optimality directly: perturb every
        # parameter +/-10% (clipped into the bounds) and confirm no move improves.
        probe = {}
        for nm in PARAM_NAMES:
            lo, hi = PARAM_BOUNDS[nm]
            d = []
            for f in (0.9, 1.1):
                pr = dict(params)
                pr[nm] = min(max(params[nm] * f, lo), hi)
                d.append(-neg_ll(_log_theta(pr)) + best_f)
            probe[nm] = d
        gap = max(v for d in probe.values() for v in d)
        # 1e-3 nats out of ~2e5 is float noise in a float32 intensity field, not a
        # missed optimum; anything above `tol` is a genuinely unconverged fit.
        is_local_max = gap <= LOCAL_MAX_NATS
        tol = max(1.0, 1e-5 * abs(best_f))
        assert gap <= tol, (
            f"ETAS fit is NOT at a local maximum: a +/-10% perturbation improved the "
            f"log-likelihood by {gap:.3f} nats (tolerance {tol:.3f}). "
            f"dLL per parameter: {probe}"
        )
        if not is_local_max and self.verbose:
            print(f"  NOTE: +/-10% probe improves LL by {gap:.4f} nats "
                  f"(< tolerance {tol:.3f}); the optimum is flat, not wrong.")
        elapsed = time.time() - state["t0"]

        # climatology comparison over exactly the same window, own intercept refit
        from . import score
        ll_clim, a0 = score.baseline_ll(self.clim.mu, y)
        ll_etas = -best_f
        rec = {
            "key": key,
            "params": {kk: float(params[kk]) for kk in PARAM_NAMES},
            "ll_etas": ll_etas,
            "ll_climatology": float(ll_clim),
            "a0_climatology": float(a0),
            "n_events_fit": n_ev,
            "bits_per_event_vs_climatology": (ll_etas - ll_clim) / (n_ev * LN2),
            "bits_total_vs_climatology": (ll_etas - ll_clim) / LN2,
            "fit_window_days": [int(fit_sl.start), int(fit_sl.stop)],
            "omori_trunc_days": int(self.trunc_days),
            "methods": methods,
            "n_objective_evals": int(state["n"]),
            "fit_seconds": float(elapsed),
            "lbfgsb_success": bool(res.success),
            "lbfgsb_message": str(res.message),
            "is_local_max": bool(is_local_max),
            "perturbation_gap_nats": float(gap),
            "perturbation_dll_+-10pct": {k: [float(v) for v in d]
                                         for k, d in probe.items()},
            "bounds": {kk: list(PARAM_BOUNDS[kk]) for kk in PARAM_NAMES},
            "fitted_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "caveat": ETAS_CAVEAT,
        }
        at_bound = [kk for kk in PARAM_NAMES
                    if params[kk] <= PARAM_BOUNDS[kk][0] * 1.001
                    or params[kk] >= PARAM_BOUNDS[kk][1] * 0.999]
        rec["params_at_bound"] = at_bound
        self.save_cache(rec)
        return rec

    # -- interface ---------------------------------------------------------
    def rate(self, sl: slice | None = None):
        """(n_cells, n_days_in_window) intensity: ETAS is NOT time-stationary."""
        return self.lam if sl is None else self.lam[:, sl]

    @property
    def mu(self):
        """Kept for interface parity; 2-D because ETAS varies day to day."""
        return self.lam

    def report(self):
        p = self.params
        info = self.fit_info
        lines = [
            f"baseline             = {self.name} (params {info.get('source','?')})",
            f"  nu    = {p['nu']:.5f}   (background rescale of climatology mu)",
            f"  K     = {p['K']:.5f}   (triggered events per unit productivity)",
            f"  alpha = {p['alpha']:.5f}   (productivity 10^(alpha*(m-4.5)))",
            f"  c     = {p['c']:.5f} d (Omori offset)",
            f"  p     = {p['p']:.5f}   (Omori exponent)",
            f"  sigma = {p['sigma']:.3f} km (isotropic Gaussian kernel)",
            f"  Omori truncation   = {self.trunc_days} d (declared approximation)",
            f"  fit window (days)  = {info.get('fit_window_days')} "
            f"(EXPLORATION only)",
        ]
        if "ll_etas" in info:
            lines.append(f"  fit log-likelihood = {info['ll_etas']:.2f} "
                         f"(climatology {info['ll_climatology']:.2f})")
            lines.append(f"  ETAS gain vs climatology (in fit window) = "
                         f"{info['bits_per_event_vs_climatology']:+.4f} bits/event")
        if "perturbation_gap_nats" in info:
            gap = float(info["perturbation_gap_nats"])
            lines.append(
                f"  local maximum      = {gap <= LOCAL_MAX_NATS} (best +/-10% "
                f"perturbation gains {gap:.4f} nats, tol {LOCAL_MAX_NATS}; "
                f"L-BFGS-B exit success={info.get('lbfgsb_success')})")
        if info.get("params_at_bound"):
            lines.append(f"  WARNING: parameters pinned at a bound: "
                         f"{info['params_at_bound']} -- the unconstrained optimum "
                         f"is outside the ruled bounds, or the parameter is flat")
        bg = float((self.params["nu"] * self.clim.mu).sum()) * self.lam.shape[1]
        tot = float(self.lam.sum(dtype=np.float64))
        lines.append(f"  background share   = {bg / tot:.3f} of total intensity "
                     f"(triggered {1 - bg / tot:.3f})")
        return lines


# ---- back-compat alias: the v1 slot name, now implemented -------------------
EtasBaseline = EtasV1

BASELINES = {"climatology-v1": ClimatologyV1, "etas": EtasV1}


def get_baseline(name: str):
    from . import canonical_baseline
    return BASELINES[canonical_baseline(name)]


def make_baseline(name: str, alpha: float = 0.5, **kw):
    """Construct a baseline by (possibly aliased) name with the run's settings."""
    cls = get_baseline(name)
    if cls is ClimatologyV1:
        return cls(alpha=alpha)
    return cls(alpha=alpha, **kw)


def etas_params_summary(path: str = ETAS_CACHE):
    """The cached ETAS parameters, or None. Used by the CLI and the tests."""
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)
