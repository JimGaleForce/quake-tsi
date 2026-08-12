"""The cached-basis Lomb-Scargle kernel (EQ-24 v2 phase 1c). Its whole licence.

Phase 1c replaced `scipy.signal.lombscargle` with a kernel that precomputes the
x-independent half of the periodogram once per (time grid, period grid) pair and
then costs one BLAS pass per series. That is a DELIBERATE KERNEL CHANGE under
§P6-5, and a deliberate kernel change is licensed by evidence, not by intent. This
file is the evidence:

  1. NUMERICAL EQUIVALENCE. Max relative deviation from the library it replaced,
     below 1e-10, measured on THREE REAL SERIES and THREE SYNTHETIC ones across the
     FULL production period grid -- not at a few frequencies, and not on toy data.
  2. ONE KERNEL, STRUCTURALLY. The observed statistic and every surrogate of every
     null are scored by the same function because there is no other function to
     reach and no selector to set. A rank is only meaningful if both sides of the
     comparison were computed the same way, and the way to guarantee that is to make
     the alternative unreachable rather than to remember not to use it.
  3. THE MEMORY GUARD. Each worker process holds its own basis. The footprint is
     computed, printed, and refused above a per-worker ceiling with a message that
     names the fix.
  4. NO HIDDEN BATCH STATE. One series in, one periodogram out, whatever was
     evaluated before it -- which is what keeps the surrogate chunk size a pure
     execution knob (see test_mine_parallel.py::test_period_chunk_size_invariance).
"""

from __future__ import annotations

import inspect
import math

import numpy as np
import pytest

from engine import mine as M, mine_session as ms

# The production grid at the quick preset: exactly what a real scan evaluates.
N_PERIODS = int(ms.QUICK["n_periods"])
TOL = 1e-10


def _production_periods(n=N_PERIODS):
    return np.exp(np.linspace(math.log(ms.PERIOD_MIN), math.log(ms.PERIOD_MAX), n))


def _real_series():
    """Three REAL series on the real mining time grid, from the cached design.

    Real means the actual global earthquake catalogue through the actual design
    pipeline -- daily domain counts at two magnitude targets and the trailing
    global b-value -- not a synthetic stand-in. The ETAS fit is deliberately NOT
    run here (it costs minutes and the kernel does not know what the numbers mean);
    the series below have the real record's gaps, clustering and heavy tail, which
    is what a periodogram kernel can actually be wrong about.
    """
    from engine import datasets, design, splits
    ctx = design.build_design(data_dir=datasets.DEFAULT_DATA_DIR, dlat=0.5,
                              dlon=0.5, explore_frac=0.7, verbose=False)
    explore, _hold = splits.temporal_split(ctx.n_days, 0.7)
    window = slice(365, explore.stop)          # the production 365 d ETAS burn-in
    out = {}
    for mag in (4.5, 5.0):
        out[f"real_daily_counts_M{mag}"] = ctx.day_counts(mag)[:, window].sum(
            axis=0).astype(np.float64)
    n_d = np.zeros((1, ctx.n_days), dtype=np.float64)
    np.add.at(n_d[0], ctx.ev_day, 1.0)
    m_d = np.zeros((1, ctx.n_days), dtype=np.float64)
    np.add.at(m_d[0], ctx.ev_day, ctx.ev_mag.astype(np.float64))
    tr = design.EngineContext.trailing_sum
    n90 = tr(n_d.astype(np.float32), 90)[0].astype(np.float64)
    m90 = tr(m_d.astype(np.float32), 90)[0].astype(np.float64)
    mean_m = np.where(n90 > 0, m90 / np.maximum(n90, 1e-9), 4.8)
    b90 = np.clip(1.0 / (math.log(10.0) * np.maximum(mean_m - 4.45, 1e-3)), 0.2, 3.0)
    out["real_trailing_b_value_90d"] = b90[window]
    return out


def _synthetic_series(n):
    rng = np.random.default_rng(20260812)
    t = np.arange(n, dtype=np.float64)
    return {
        "synthetic_white_noise": rng.normal(size=n),
        "synthetic_ar1_phi0.6": M.ar1_draw(n, 0.6, 1.0, np.random.default_rng(11)),
        "synthetic_planted_29.53d": (0.3 * np.cos(2 * np.pi * t / 29.53 - 0.7)
                                     + rng.normal(size=n)),
    }


def _max_rel_dev(t, x, periods):
    got = M.lomb_scargle_power(t, x, periods)
    ref = M.lomb_scargle_power_scipy(t, x, periods)
    assert got.shape == ref.shape
    return float(np.max(np.abs(got - ref) / np.maximum(np.abs(ref), 1e-300)))


# ============================================================================
# 1. NUMERICAL EQUIVALENCE -- the §P6-5 licence for the kernel swap
# ============================================================================
def test_cached_kernel_matches_scipy_on_three_real_and_three_synthetic_series():
    real = _real_series()
    assert len(real) == 3
    n = next(iter(real.values())).size
    t = np.arange(n, dtype=np.float64)
    periods = _production_periods()
    assert periods.size == N_PERIODS

    worst, report = 0.0, []
    for name, x in list(real.items()) + list(_synthetic_series(n).items()):
        assert x.size == n
        d = _max_rel_dev(t, x, periods)
        report.append(f"    {name:<28s} max|rel dev| = {d:.3e}")
        worst = max(worst, d)
    print(f"\n  cached-basis LS kernel vs scipy.signal.lombscargle, "
          f"{n} days x {periods.size} trial periods:")
    print("\n".join(report))
    print(f"    {'WORST':<28s} max|rel dev| = {worst:.3e}  (bar {TOL:.0e})")
    assert worst < TOL, f"kernel deviates from scipy by {worst:.3e} >= {TOL:.0e}"


def test_peak_ordering_and_argmax_are_identical_to_scipy():
    """Deviation at the 1e-12 level is only harmless if it moves no RANK.

    The scan uses exactly two things from a periodogram: the argmax (every
    surrogate) and the ordered local maxima (the observed series). Both are
    compared here directly, because a max relative deviation is not by itself an
    argument that a near-tie did not flip.
    """
    n = 4000
    t = np.arange(n, dtype=np.float64)
    periods = _production_periods()
    for name, x in _synthetic_series(n).items():
        got = M.lomb_scargle_power(t, x, periods)
        ref = M.lomb_scargle_power_scipy(t, x, periods)
        assert int(got.argmax()) == int(ref.argmax()), name
        assert np.array_equal(np.argsort(got), np.argsort(ref)), \
            f"{name}: the full power ORDERING moved, so a rank could flip"
        pk_got = M.period_observed(t, x, periods, n_peaks=8)[1]
        assert pk_got, name


# ============================================================================
# 2. ONE KERNEL, STRUCTURALLY -- not one kernel by convention
# ============================================================================
def test_there_is_exactly_one_kernel_and_no_way_to_select_another():
    """No flag, no kwarg, no environment switch: nothing to get wrong per call."""
    sig = inspect.signature(M.lomb_scargle_power)
    assert list(sig.parameters) == ["t", "x", "periods"], (
        "lomb_scargle_power grew a parameter. If it is a kernel selector, the "
        "observed statistic and its surrogates can be scored differently and every "
        "rank in the scan becomes unsafe -- see the header of this file.")
    for p in sig.parameters.values():
        assert p.default is inspect.Parameter.empty

    # the scipy reference must be reachable from tests and from NOTHING else
    src = inspect.getsource(M)
    body_uses = [ln.strip() for ln in src.splitlines()
                 if "lomb_scargle_power_scipy" in ln
                 and not ln.strip().startswith("#")]
    assert body_uses == ["def lomb_scargle_power_scipy(t, x, periods):"], (
        f"the pre-1c scipy kernel is referenced inside engine.mine: {body_uses}. "
        f"It exists only as a validation reference and must never be reachable "
        f"from a scan.")
    assert "lomb_scargle_power_scipy" not in inspect.getsource(ms)

    # and both consumers go through the one kernel
    for fn in (M.period_observed, M.period_surrogate_maxima):
        s = inspect.getsource(fn)
        assert "lomb_scargle_power(" in s
        assert "scipy" not in s and "signal." not in s


def test_observed_and_surrogates_share_the_basis_object_bit_for_bit():
    """Same grids => literally the same cached arrays on both sides of the rank."""
    n = 600
    t = np.arange(n, dtype=np.float64)
    periods = _production_periods(200)
    b1 = M.ls_basis(t, periods)
    M.period_observed(t, np.random.default_rng(1).normal(size=n), periods, n_peaks=3)
    M.period_surrogate_maxima(t, np.random.default_rng(2).normal(size=n), periods,
                              "permutation", 0, 3, M.period_seed_seqs(5)["permutation"])
    assert M.ls_basis(t, periods) is b1


# ============================================================================
# 3. THE MEMORY GUARD -- each worker pays for its own copy
# ============================================================================
def test_quick_preset_basis_footprint_is_stated_and_under_the_ceiling():
    n_time, n_freq = 7716, N_PERIODS          # the real mining window, quick preset
    nbytes = M.ls_basis_footprint_bytes(n_time, n_freq)
    assert nbytes == 2 * n_time * n_freq * 8
    print(f"\n  per-worker basis at {n_time} x {n_freq} = {nbytes / 2**20:.1f} MiB "
          f"(ceiling {M.LS_BASIS_MAX_BYTES / 2**20:.0f} MiB); at --jobs 30 the "
          f"total is {30 * nbytes / 2**30:.2f} GiB")
    assert nbytes < M.LS_BASIS_MAX_BYTES


def test_a_grid_over_the_ceiling_fails_loudly_and_names_the_fix():
    n_time = 7716
    n_freq = int(M.LS_BASIS_MAX_BYTES // (2 * n_time * 8)) + 1000
    assert M.ls_basis_footprint_bytes(n_time, n_freq) > M.LS_BASIS_MAX_BYTES
    t = np.arange(n_time, dtype=np.float64)
    periods = np.exp(np.linspace(math.log(2.0), math.log(4000.0), n_freq))
    with pytest.raises(MemoryError) as exc:
        M.build_ls_basis(t, periods, announce=False)
    msg = str(exc.value)
    assert "--jobs" in msg and "MiB" in msg and "n_periods" in msg


def test_every_declared_preset_footprint_is_stated_for_the_reader():
    """The per-worker ceiling is not the whole bill: --jobs multiplies it.

    Printed rather than asserted for the larger presets, because the aggregate is a
    property of --jobs and the worker cannot see --jobs. What IS asserted is that
    the number is computable in advance, which is what lets `mine_session.run`
    print the total before the pool is created rather than after the machine swaps.
    """
    rows = []
    for preset in (ms.QUICK, ms.DEFAULT, ms.OVERNIGHT):
        nb = M.ls_basis_footprint_bytes(7716, int(preset["n_periods"]))
        rows.append(f"    {preset['label']:<10s} n_periods={preset['n_periods']:>5d} "
                    f"-> {nb / 2**20:7.1f} MiB/worker, "
                    f"{30 * nb / 2**30:5.2f} GiB at --jobs 30, per-worker ceiling "
                    f"{'OK' if nb <= M.LS_BASIS_MAX_BYTES else 'REFUSED'}")
    print("\n  cached-basis footprint by preset (7716-day mining window):")
    print("\n".join(rows))
    assert M.ls_basis_footprint_bytes(7716, int(ms.QUICK["n_periods"])) \
        <= M.LS_BASIS_MAX_BYTES


# ============================================================================
# 4. CACHING IS AN OPTIMISATION, NOT A STATE MACHINE
# ============================================================================
def test_basis_is_keyed_on_content_and_reused_within_a_process():
    n = 400
    t = np.arange(n, dtype=np.float64)
    periods = _production_periods(120)
    a = M.ls_basis(t, periods)
    assert M.ls_basis(t, periods) is a                    # identity fast path
    assert M.ls_basis(t.copy(), periods.copy()) is a       # content digest path
    assert M.ls_basis_digest(t, periods) == M.ls_basis_digest(t.copy(), periods)
    other = M.ls_basis(t + 1.0, periods)
    assert other is not a
    assert M.ls_basis_digest(t + 1.0, periods) != M.ls_basis_digest(t, periods)
    assert M.ls_basis(t, _production_periods(121)) is not a


def test_kernel_has_no_batch_state_so_results_cannot_depend_on_call_order():
    n = 500
    t = np.arange(n, dtype=np.float64)
    periods = _production_periods(150)
    rng = np.random.default_rng(4)
    xs = [rng.normal(size=n) for _ in range(4)]
    alone = [M.lomb_scargle_power(t, x, periods) for x in xs]
    interleaved = []
    for x in xs:
        M.lomb_scargle_power(t, rng.normal(size=n), periods)   # noise between calls
        M.lomb_scargle_power(t + 3.0, x, periods)              # and another grid
        interleaved.append(M.lomb_scargle_power(t, x, periods))
    for a, b in zip(alone, interleaved):
        assert np.array_equal(a, b), "kernel output depended on what ran before it"


def test_kernel_rejects_a_series_of_the_wrong_length():
    t = np.arange(200, dtype=np.float64)
    periods = _production_periods(50)
    with pytest.raises(ValueError, match="samples"):
        M.lomb_scargle_power(t, np.zeros(199), periods)
