"""Covariate plugin registry.

    @register("name")
    def compute(ctx: EngineContext, params: dict) -> np.ndarray  # (n_cells, n_days) f32

Contract every plugin must honour:
  * CAUSAL: z[:, t] may use only events strictly before day t. Use ctx.past_counts /
    ctx.trailing_sum, which are exclusive-of-today by construction. The test suite
    shuffles future events and asserts z does not change.
  * FINITE: no NaN/Inf anywhere in the returned array.
  * Never touch raw CSVs; go through ctx.
Register a `burn_in` (days of history the covariate needs) via the DEFAULTS dict so
the runner can drop the un-warmed head of the record.
"""

from __future__ import annotations

import importlib
import pkgutil

import numpy as np

_REGISTRY = {}


def register(name, defaults=None, burn_in=0, describe=""):
    def deco(fn):
        fn.cov_name = name
        fn.defaults = dict(defaults or {})
        fn.burn_in = burn_in
        fn.describe = describe
        _REGISTRY[name] = fn
        return fn
    return deco


def _discover():
    for m in pkgutil.iter_modules(__path__):
        if not m.name.startswith("_"):
            importlib.import_module(f"{__name__}.{m.name}")


def available():
    _discover()
    return dict(sorted(_REGISTRY.items()))


def get(name):
    _discover()
    if name not in _REGISTRY:
        raise KeyError(f"unknown covariate {name!r}; have {sorted(_REGISTRY)}")
    return _REGISTRY[name]


def compute(name, ctx, params=None):
    """Compute, validate and standardise a covariate.

    Standardisation uses EXPLORATION-period statistics only (mean/sd over the training
    window), so a holdout run sees exactly the transform that was fitted in training.
    """
    fn = get(name)
    p = dict(fn.defaults)
    p.update(params or {})
    z = fn(ctx, p)
    z = np.asarray(z, dtype=np.float32)
    if z.shape != (ctx.n_cells, ctx.n_days):
        raise ValueError(f"{name}: expected {(ctx.n_cells, ctx.n_days)}, got {z.shape}")
    if not np.isfinite(z).all():
        raise ValueError(f"{name}: covariate contains NaN/Inf")
    burn = int(p.get("burn_in", fn.burn_in))
    train = z[:, burn:ctx.n_explore_days]
    mean = float(train.mean())
    sd = float(train.std())
    if sd <= 0:
        raise ValueError(f"{name}: covariate is constant over the training window "
                         f"(produced 0 variance)")
    zs = ((z - mean) / sd).astype(np.float32)
    stats = {"raw_mean": mean, "raw_sd": sd, "burn_in": burn, "params": p,
             "min": float(zs.min()), "max": float(zs.max())}
    return zs, stats
