"""Tests for exp_bvalue_skill.py (P-2.2 / K-405).

Three things must be true before any number from this arm is worth reading:
  1. the Aki-Utsu estimator recovers a known b from a synthetic GR sample;
  2. the bin-aware truncated-GR likelihood is a genuine normalised pmf over the bins
     it claims as its support, and the local-vs-global bit gain is exactly the
     log2 likelihood ratio;
  3. the whole vectorised local pipeline recovers an INJECTED b(x) field on a small
     synthetic catalogue, and returns ~0 bits when there is no b structure to find.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import exp_bvalue_skill as E  # noqa: E402


def _gr_sample(rng, n, b, mc, binw=0.1, mmax=8.0):
    """Discrete truncated-GR magnitudes on the `binw` grid, floor mc, top mmax."""
    beta = b * E.LN10
    J = int(round((mmax - mc) / binw))
    z = rng.random(n)
    span = -np.expm1(-beta * binw * (J + 1))
    j = np.clip(np.floor(-np.log1p(-z * span) / (beta * binw)).astype(int), 0, J)
    return mc + j * binw


def _gr_sample_vec(rng, b, mc, binw=0.1, mmax=8.0):
    """Per-event draw from a truncated GR with a per-event b."""
    b = np.asarray(b, dtype=float)
    beta = b * E.LN10
    J = int(round((mmax - mc) / binw))
    z = rng.random(b.size)
    span = -np.expm1(-beta * binw * (J + 1))
    j = np.clip(np.floor(-np.log1p(-z * span) / (beta * binw)).astype(int), 0, J)
    return mc + j * binw


# --------------------------------------------------------------- 1. Aki-Utsu
@pytest.mark.parametrize("b_true", [0.7, 1.0, 1.3])
def test_aki_recovers_known_b(b_true):
    rng = np.random.default_rng(7)
    mags = _gr_sample(rng, 200000, b_true, mc=1.0)
    b_hat, n = E.aki_b(mags, 1.0)
    assert n == 200000
    assert abs(b_hat - b_true) < 0.02, (b_hat, b_true)


def test_aki_bin_correction_matters():
    """Without the binw/2 correction the estimator is biased; with it, it is not."""
    rng = np.random.default_rng(11)
    mags = _gr_sample(rng, 200000, 1.0, mc=1.0)
    corrected, _ = E.aki_b(mags, 1.0, binw=0.1)
    uncorrected = 1.0 / (E.LN10 * (mags.mean() - 1.0))
    assert abs(corrected - 1.0) < 0.02
    assert abs(uncorrected - 1.0) > 0.02
    assert abs(corrected - 1.0) < abs(uncorrected - 1.0)


def test_maxc_mc_finds_the_modal_bin():
    rng = np.random.default_rng(3)
    mags = _gr_sample(rng, 50000, 1.0, mc=1.4)
    assert abs(E.maxc_mc(mags) - (1.4 + E.MC_OFFSET)) < 1e-9


# ------------------------------------------------- 2. the bin-aware likelihood
def test_gr_bin_logprob_is_a_normalised_pmf():
    for b in (0.6, 1.0, 1.5):
        for mc in (0.5, 1.7):
            bins = np.round(np.arange(mc, E.MMAX_TRUNC + 1e-9, E.BIN_W), 6)
            p = np.exp(E.gr_bin_logprob(bins, mc, b))
            assert p.min() > 0
            assert abs(p.sum() - 1.0) < 1e-9, (b, mc, p.sum())


def test_gr_bin_logprob_matches_the_analytic_bin_mass():
    b, mc, binw = 1.1, 1.0, 0.1
    beta = b * E.LN10
    lo = mc - binw / 2.0
    top = E.MMAX_TRUNC + binw / 2.0 - lo
    for m in (1.0, 1.5, 3.2):
        a = m - binw / 2.0 - lo
        want = ((np.exp(-beta * a) - np.exp(-beta * (a + binw)))
                / (1.0 - np.exp(-beta * top)))
        got = np.exp(E.gr_bin_logprob(m, mc, b, binw=binw))
        assert abs(got - want) < 1e-12


def test_bit_gain_is_the_log2_likelihood_ratio_and_is_zero_when_b_is_equal():
    m, mc = 2.3, 1.0
    same = (E.gr_bin_logprob(m, mc, 1.0) - E.gr_bin_logprob(m, mc, 1.0)) / E.LN2
    assert abs(same) < 1e-15
    diff = (E.gr_bin_logprob(m, mc, 0.8) - E.gr_bin_logprob(m, mc, 1.2)) / E.LN2
    ref = np.log2(np.exp(E.gr_bin_logprob(m, mc, 0.8))
                  / np.exp(E.gr_bin_logprob(m, mc, 1.2)))
    assert abs(diff - ref) < 1e-10


def test_true_b_beats_a_wrong_b_in_expectation():
    """A correctly specified model must earn positive bits against a wrong one."""
    rng = np.random.default_rng(19)
    mags = _gr_sample(rng, 100000, 1.2, mc=1.0)
    good = E.gr_bin_logprob(mags, 1.0, 1.2)
    bad = E.gr_bin_logprob(mags, 1.0, 0.9)
    assert (good - bad).mean() / E.LN2 > 0.05


# -------------------------------------------- 3. the vectorised local estimator
def test_local_mc_b_matches_the_scalar_estimator_row_by_row():
    rng = np.random.default_rng(23)
    rows = []
    for _ in range(40):
        mc = rng.choice([0.4, 0.9, 1.5])
        b = rng.uniform(0.7, 1.4)
        m = _gr_sample(rng, E.K_LOCAL, b, mc=mc)
        rows.append(np.round(m / E.BIN_W).astype(np.int32))
    qn = np.stack(rows)
    mc_bin, b_vec, n_sub, ok = E.local_mc_b(qn)
    assert ok.sum() >= 35
    for i in range(len(rows)):
        if not ok[i]:
            continue
        mags = rows[i] * E.BIN_W
        mc_scalar = E.maxc_mc(mags)
        assert abs(mc_bin[i] * E.BIN_W - mc_scalar) < 1e-9
        b_scalar, ns = E.aki_b(mags, mc_scalar)
        assert ns == n_sub[i]
        assert abs(b_vec[i] - b_scalar) < 1e-9


# ------------------------------------------------------- injection recovery
def _synthetic_catalogue(n, rng, b_field=None, mc=0.5):
    """A small clustered SoCal-shaped synthetic: uniform in a 2 deg box, 4 years."""
    lat = 33.0 + rng.random(n) * 1.0
    lon = -118.0 + rng.random(n) * 1.0
    t = np.sort(rng.random(n) * 1460.0)
    if b_field is None:
        b = np.full(n, 1.0)
    else:
        b = b_field(lon)
    mags = _gr_sample_vec(rng, b, mc=mc)
    return t, lat, lon, mags


def _run_small(t, lat, lon, mags):
    q_all = np.round(mags / E.BIN_W).astype(np.int32)
    x, y, _, _ = E.project_km(lat, lon)
    is_train = t < 0.6 * t.max()
    pool = np.flatnonzero(~is_train)
    nbr, ok, _ = E.neighbour_sets(pool, t, x, y, chunk=400, kq=400, verbose=False)
    assert ok.sum() > 100, "expected many local estimates, got %d" % ok.sum()
    res = E.score_pass(q_all, nbr, ok, pool, q_all[is_train])
    assert res is not None
    return res


def test_injection_recovery_positive_and_null_flat():
    rng = np.random.default_rng(101)
    n = 12000
    t, lat, lon, _ = _synthetic_catalogue(n, rng)

    # (i) NO b structure -> the local estimator must not manufacture bits
    flat = _gr_sample_vec(rng, np.full(n, 1.0), mc=0.5)
    gain_flat, _, _, _ = _run_small(t, lat, lon, flat)
    assert abs(gain_flat.mean()) < 0.02, gain_flat.mean()

    # (ii) a strong injected b(x) -> the identical pipeline must earn bits
    b = 1.0 + 0.45 * (2.0 * (lon - lon.min()) / (lon.max() - lon.min()) - 1.0)
    inj = _gr_sample_vec(rng, b, mc=0.5)
    gain_inj, _, _, _ = _run_small(t, lat, lon, inj)
    assert gain_inj.mean() > 0.03, gain_inj.mean()
    assert gain_inj.mean() > gain_flat.mean() + 0.02


# --------------------------------------------------------------- housekeeping
def test_neighbour_sets_are_strictly_prior_and_within_the_window():
    rng = np.random.default_rng(5)
    n = 20000
    t, lat, lon, mags = _synthetic_catalogue(n, rng)
    x, y, _, _ = E.project_km(lat, lon)
    pool = np.flatnonzero(t > 0.6 * t.max())
    nbr, ok, _ = E.neighbour_sets(pool, t, x, y, chunk=300, kq=300, verbose=False)
    assert ok.any()
    sub = nbr[ok]
    tgt = pool[ok]
    assert (sub < tgt[:, None]).all()
    d = np.hypot(x[sub] - x[tgt][:, None], y[sub] - y[tgt][:, None])
    assert (d <= E.RADIUS_KM + 1e-9).all()
    assert (t[sub] >= (t[tgt] - E.LOOKBACK_DAYS)[:, None]).all()
    assert sub.shape[1] == E.K_LOCAL


def test_block_bootstrap_ci_brackets_the_mean():
    rng = np.random.default_rng(13)
    vals = rng.normal(0.5, 1.0, 4000)
    blocks = np.repeat(np.arange(400), 10)
    lo, hi, nb = E.block_bootstrap_mean(vals, blocks, n_boot=400)
    assert nb == 400
    assert lo < vals.mean() < hi


def test_measured_rounding_of_a_known_grid():
    rng = np.random.default_rng(2)
    assert E.measure_rounding(np.round(rng.random(5000) * 5, 2)) == 0.01
    assert E.measure_rounding(np.round(rng.random(5000) * 5, 1)) == 0.1


def test_shuffle_preserves_the_per_group_magnitude_multiset():
    rng = np.random.default_rng(17)
    n = 5000
    lat = 33.0 + rng.random(n) * 2
    lon = -118.0 + rng.random(n) * 2
    year = rng.integers(2008, 2012, n)
    q = rng.integers(-10, 50, n).astype(np.int32)
    order, starts = E.shuffle_groups(lat, lon, year)
    qs = E.shuffle_groups and E.shuffled_mag(q, order, starts, rng)
    key_lat = np.floor(lat / E.NULL_CELL_DEG).astype(np.int64)
    key_lon = np.floor(lon / E.NULL_CELL_DEG).astype(np.int64)
    key = key_lat * 1000000 + key_lon * 10 + 0
    key = key * 10000 + year
    for k in np.unique(key)[:25]:
        m = key == k
        assert sorted(q[m].tolist()) == sorted(qs[m].tolist())
