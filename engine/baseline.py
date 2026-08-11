"""Baselines. v1 = per-cell Poisson climatology; ETAS slot is deliberately empty."""

from __future__ import annotations

import numpy as np

from . import BASELINE_CAVEAT

ETAS_MESSAGE = (
    "baseline='etas' is not implemented in v1. This is the engine's declared gap: the "
    "v1 baseline is a time-stationary per-cell Poisson climatology, so aftershock "
    "clustering is NOT absorbed by the baseline and leaks into every covariate. Any "
    "covariate that proxies 'there were recent nearby earthquakes' will therefore score "
    "apparent skill against climatology-v1 that an ETAS baseline would eat. Implement "
    "engine/baseline.py:EtasBaseline before treating any result here as a discovery."
)


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

    def report(self):
        return [
            f"baseline             = {self.name} (alpha={self.alpha})",
            f"  train days         = {self.n_train_days}",
            f"  train events       = {self.n_train_events:.0f}",
            f"  mean cell rate/day = {self.mu_bar:.3e}",
            f"  mu range           = {self.mu.min():.3e} .. {self.mu.max():.3e}",
            f"  cells empty in train = {self.n_empty_cells} (smoothed, not zero)",
        ]


class EtasBaseline:
    name = "etas"

    def fit(self, *a, **k):
        raise NotImplementedError(ETAS_MESSAGE)


BASELINES = {"climatology-v1": ClimatologyV1, "etas": EtasBaseline}


def get_baseline(name: str):
    if name not in BASELINES:
        raise KeyError(f"unknown baseline {name!r}; have {sorted(BASELINES)}")
    return BASELINES[name]
