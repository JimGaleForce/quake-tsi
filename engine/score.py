"""Poisson regression, information gain in bits/event, Molchan, green/red export.

Model (SPEC rule 5):   lambda(x,t) = mu(x) * exp(a + beta * z(x,t))

`a` is a nuisance re-normalisation fitted jointly with beta. The BASELINE is refit
with its own intercept (closed form) over the same evaluation window, so bits/event
compares two equally-renormalised models and cannot be manufactured by baseline
miscalibration alone.

`mu` here is whatever the chosen baseline offers as its OFFSET, and may be either:
  * 1-D, shape (n_cells,)          -- a time-stationary rate (climatology-v1), or
  * 2-D, shape (n_cells, n_window) -- a per-cell-per-day rate (etas-v1),
already sliced to the same day window as `z` and `y`. The intercept-refit convention
is identical in both cases: under `--baseline etas`, lambda_etas simply takes the
place of mu as the offset.
"""

from __future__ import annotations

import json
import math
import os

import numpy as np

LN2 = math.log(2.0)


# ------------------------------------------------------- offset handling ---
def _mu_chunk(mu, s, e):
    """Offset for day-window columns [s, e): (n_cells, 1) if 1-D else (n_cells, e-s)."""
    mu = np.asarray(mu)
    if mu.ndim == 1:
        return mu[:, None].astype(np.float64)
    return mu[:, s:e].astype(np.float64)


def _mu_total(mu, n_days):
    """Total baseline expectation over a window of n_days days."""
    mu = np.asarray(mu)
    if mu.ndim == 1:
        return float(mu.sum()) * n_days
    return float(mu.sum(dtype=np.float64))


def check_offset(mu, z):
    mu = np.asarray(mu)
    if mu.ndim == 1:
        if mu.shape[0] != z.shape[0]:
            raise ValueError(f"offset {mu.shape} does not match design {z.shape}")
    elif mu.shape != z.shape:
        raise ValueError(f"offset {mu.shape} does not match design window {z.shape}")
    if not np.isfinite(mu).all() or (mu <= 0).any():
        raise ValueError("baseline offset must be finite and strictly positive")
    return mu


# ---------------------------------------------------------------- fitting ---
def _accumulate(mu, z, y, a, b, chunk=512):
    """Return LL, grad(2,), hess(2,2) at (a, b), chunked over days."""
    n_days = z.shape[1]
    ll = 0.0
    g = np.zeros(2)
    h = np.zeros((2, 2))
    for s in range(0, n_days, chunk):
        e = min(s + chunk, n_days)
        zz = z[:, s:e].astype(np.float64)
        yy = y[:, s:e].astype(np.float64)
        lam = _mu_chunk(mu, s, e) * np.exp(a + b * zz)
        ll += float((yy * np.log(lam) - lam).sum())
        r = yy - lam
        g[0] += float(r.sum())
        g[1] += float((zz * r).sum())
        h[0, 0] += float(lam.sum())
        h[0, 1] += float((zz * lam).sum())
        h[1, 1] += float((zz * zz * lam).sum())
    h[1, 0] = h[0, 1]
    return ll, g, h


def fit_poisson(mu, z, y, max_iter=50, tol=1e-9, verbose=False):
    """Newton fit of (a, beta). Returns dict with beta, se_beta, converged, ll."""
    mu = check_offset(mu, z)
    a, b = 0.0, 0.0
    # sensible start for the intercept: match total counts at beta=0
    tot_y = float(y.sum())
    tot_mu = _mu_total(mu, z.shape[1])
    if tot_y <= 0:
        raise ValueError("bulk transform produced 0 events in the fitting window")
    a = math.log(tot_y / tot_mu)

    ll, g, h = _accumulate(mu, z, y, a, b)
    converged = False
    n_iter = 0
    d = float("nan")
    for n_iter in range(1, max_iter + 1):
        try:
            step = np.linalg.solve(h, g)     # H_negative = h, Newton step = h^-1 g
        except np.linalg.LinAlgError:
            break
        t = 1.0
        for _ in range(30):                  # step halving
            a2, b2 = a + t * step[0], b + t * step[1]
            ll2, g2, h2 = _accumulate(mu, z, y, a2, b2)
            if np.isfinite(ll2) and ll2 >= ll:
                break
            t *= 0.5
        else:
            break
        d = abs(ll2 - ll)
        a, b, ll, g, h = a2, b2, ll2, g2, h2
        if verbose:
            print(f"    newton {n_iter:2d}  a={a:+.6f} beta={b:+.6f} ll={ll:.6f}")
        if d < tol * max(1.0, abs(ll)):
            converged = True
            break

    cov = np.linalg.inv(h)
    se = float(np.sqrt(max(cov[1, 1], 0.0)))
    assert converged, (
        f"Poisson fit did NOT converge in {n_iter} iterations (last |dLL|={d})"
    )
    return {
        "a": float(a), "beta": float(b), "se_beta": se,
        "z_stat": float(b / se) if se > 0 else float("nan"),
        "converged": bool(converged), "n_iter": int(n_iter), "ll": float(ll),
    }


def baseline_ll(mu, y, chunk=512):
    """Log-likelihood of mu with its own fitted intercept over this window."""
    n_days = y.shape[1]
    tot_y = float(y.sum())
    tot_mu = _mu_total(mu, n_days)
    a0 = math.log(tot_y / tot_mu)
    ll = 0.0
    for s in range(0, n_days, chunk):
        e = min(s + chunk, n_days)
        yy = y[:, s:e].astype(np.float64)
        lam = (_mu_chunk(mu, s, e) * math.exp(a0)) * np.ones((1, e - s))
        ll += float((yy * np.log(lam) - lam).sum())
    return ll, a0


def model_ll(mu, z, y, a, b, chunk=512):
    ll, _, _ = _accumulate(mu, z, y, a, b, chunk=chunk)
    return ll


# ---------------------------------------------------------------- metrics ---
def information_gain(mu, z, y, a, b):
    """bits/event of the fitted model over the (re-normalised) baseline."""
    ll_m = model_ll(mu, z, y, a, b)
    ll_b, a0 = baseline_ll(mu, y)
    n_ev = float(y.sum())
    if n_ev <= 0:
        raise ValueError("no events in the evaluation window")
    return {
        "ll_model": ll_m, "ll_baseline": ll_b, "a0_baseline": a0,
        "n_events_eval": n_ev,
        "bits_total": (ll_m - ll_b) / LN2,
        "bits_per_event": (ll_m - ll_b) / (n_ev * LN2),
    }


def molchan(mu, z, y, a, b, n_bins=400, chunk=512):
    """Molchan trajectory: tau (alarmed fraction of baseline expectation) vs miss rate.

    Cells are ordered by the alarm function exp(beta*z) (highest first). Binned into
    `n_bins` quantile-free log-ratio bins so this stays O(n) in memory.
    """
    lo, hi = None, None
    n_days = z.shape[1]
    for s in range(0, n_days, chunk):
        e = min(s + chunk, n_days)
        v = (b * z[:, s:e]).astype(np.float64)
        lo = v.min() if lo is None else min(lo, v.min())
        hi = v.max() if hi is None else max(hi, v.max())
    if not np.isfinite(lo) or hi <= lo:
        hi = lo + 1e-6
    edges = np.linspace(lo, hi, n_bins + 1)

    w_mu = np.zeros(n_bins)     # baseline expectation mass per bin
    w_y = np.zeros(n_bins)      # observed events per bin
    for s in range(0, n_days, chunk):
        e = min(s + chunk, n_days)
        v = (b * z[:, s:e]).astype(np.float64)
        idx = np.clip(np.digitize(v, edges) - 1, 0, n_bins - 1)
        m = np.broadcast_to(_mu_chunk(mu, s, e), (z.shape[0], e - s))
        np.add.at(w_mu, idx.ravel(), np.ascontiguousarray(m).ravel())
        np.add.at(w_y, idx.ravel(), y[:, s:e].astype(np.float64).ravel())

    order = np.arange(n_bins)[::-1]           # highest alarm first
    tau = np.cumsum(w_mu[order]) / max(w_mu.sum(), 1e-30)
    hit = np.cumsum(w_y[order]) / max(w_y.sum(), 1e-30)
    nu = 1.0 - hit
    tau = np.concatenate([[0.0], tau])
    nu = np.concatenate([[1.0], nu])
    # area skill score: 1 - 2*area under (tau, nu); 0.5 == random, higher is better
    area = float(np.trapezoid(nu, tau))
    return {
        "tau": tau.tolist(), "nu": nu.tolist(),
        "area_under_molchan": area,
        "molchan_skill": 1.0 - 2.0 * area,
        "nu_at_tau_0.05": float(np.interp(0.05, tau, nu)),
        "nu_at_tau_0.10": float(np.interp(0.10, tau, nu)),
        "nu_at_tau_0.20": float(np.interp(0.20, tau, nu)),
    }


def green_red_cells(ctx, mu, z, y, a, b, a0, chunk=512):
    """Per-cell realized log-likelihood gain (model - baseline). Sign = green/red."""
    n_cells = z.shape[0]
    gain = np.zeros(n_cells)
    nev = np.zeros(n_cells)
    n_days = z.shape[1]
    for s in range(0, n_days, chunk):
        e = min(s + chunk, n_days)
        zz = z[:, s:e].astype(np.float64)
        yy = y[:, s:e].astype(np.float64)
        mc = _mu_chunk(mu, s, e)
        lam_m = mc * np.exp(a + b * zz)
        lam_b = mc * math.exp(a0)
        g = yy * (np.log(lam_m) - np.log(lam_b)) - (lam_m - lam_b)
        gain += g.sum(axis=1)
        nev += yy.sum(axis=1)
    return gain, nev


def export_cells(path, ctx, gain, nev, runid, header):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    cells = []
    for c in np.nonzero((nev > 0) | (np.abs(gain) > 1e-9))[0]:
        cells.append({
            "cell": int(c),
            "lat": float(ctx.grid.cell_lat[c]),
            "lon": float(ctx.grid.cell_lon[c]),
            "n_events": float(nev[c]),
            "bits": float(gain[c] / LN2),
            "color": "green" if gain[c] > 0 else ("red" if gain[c] < 0 else "grey"),
        })
    payload = {
        "runid": runid,
        "header": header,
        "n_green": sum(1 for c in cells if c["color"] == "green"),
        "n_red": sum(1 for c in cells if c["color"] == "red"),
        "cells": cells,
    }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    return payload["n_green"], payload["n_red"]
