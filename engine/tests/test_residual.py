"""Acceptance tests for engine/residual.py.

Identity tests where possible: a residual must be exactly orthogonal to what was
projected out, a series projected against itself must leave nothing, and a battery of
copies must have effective rank 1. Those cannot be satisfied by an implementation that
is merely approximately right.
"""

from __future__ import annotations

import datetime as _dt

import numpy as np
import pytest

from engine import ephemeris as E
from engine import residual as R
from engine import sitetide as ST
from engine import tidal_tensor as TT

LAT, LON = 55.35, -160.5
T0 = _dt.datetime(2020, 6, 1)
T = np.arange(0.0, 20.0, 1.0 / 1440.0)
JD = E.julian_day_at(T0, T)


def test_residual_is_exactly_orthogonal_to_the_basis():
    rng = np.random.default_rng(0)
    b = rng.normal(size=5000)
    y = 3.0 * b + rng.normal(size=5000)
    r, shared = R.project_out(y, [b])
    assert abs(float(np.dot(r - r.mean(), b - b.mean()))) < 1e-8 * np.dot(b, b)
    assert 0.0 < shared < 1.0


def test_projecting_a_series_against_itself_leaves_nothing():
    rng = np.random.default_rng(1)
    b = rng.normal(size=2000)
    r, shared = R.project_out(b, [b])
    assert shared == pytest.approx(1.0, abs=1e-12)
    assert np.max(np.abs(r)) < 1e-9 * np.std(b)


def test_projecting_out_an_unrelated_series_changes_nothing():
    rng = np.random.default_rng(2)
    y = rng.normal(size=20000)
    b = rng.normal(size=20000)
    r, shared = R.project_out(y, [b])
    assert shared < 0.01
    assert R.new_information_fraction(y, [b]) > 0.99


def test_shared_variance_equals_r_squared_for_one_basis():
    rng = np.random.default_rng(3)
    b = rng.normal(size=8000)
    y = 2.0 * b + 0.5 * rng.normal(size=8000)
    r2 = float(np.corrcoef(y, b)[0, 1]) ** 2
    assert R.shared_variance(y, [b]) == pytest.approx(r2, rel=1e-9)


def test_empty_basis_is_a_no_op():
    y = np.arange(10.0)
    r, shared = R.project_out(y, [])
    assert shared == 0.0
    assert np.allclose(r, y - y.mean())


# ------------------------------------------------------------- effective rank --
def test_effective_rank_of_identical_copies_is_one():
    rng = np.random.default_rng(4)
    base = rng.normal(size=3000)
    rep = R.effective_rank([base + 1e-12 * rng.normal(size=3000) for _ in range(10)])
    assert rep["n_nominal"] == 10
    assert rep["effective_rank"] == 1
    assert rep["variance_explained_by_first"] > 0.99


def test_effective_rank_of_independent_columns_is_full():
    rng = np.random.default_rng(5)
    cols = [rng.normal(size=20000) for _ in range(6)]
    rep = R.effective_rank(cols, tol=0.01)
    assert rep["effective_rank"] >= 5


def test_effective_rank_handles_constant_columns():
    rep = R.effective_rank([np.ones(50), np.ones(50)])
    assert rep["effective_rank"] == 0


# ----------------------------------------------------- the real tidal geometry --
def test_resolved_coulomb_is_mostly_explained_by_the_bounded_scalar():
    """K-108's number, as a standing test rather than a one-off measurement."""
    scal = ST.site_scalar(JD, LAT, LON, 0.0)
    st = TT.stress_tensor(JD, LAT, LON, 0.0)
    cfs = TT.coulomb(st, 250.0, 20.0, 90.0, 0.4)["coulomb_pa"]
    sh = R.shared_variance(cfs, [scal])
    assert sh > 0.85, sh
    assert R.new_information_fraction(cfs, [scal]) < 0.15


def test_stressing_rate_is_almost_entirely_new_information():
    """The orthogonal complement, by identity rather than by luck."""
    scal = ST.site_scalar(JD, LAT, LON, 0.0)
    rate = np.gradient(scal, T * 24.0)
    assert R.new_information_fraction(rate, [scal]) > 0.99


def test_horizontal_shear_is_almost_entirely_new_information():
    scal = ST.site_scalar(JD, LAT, LON, 0.0)
    st = TT.stress_tensor(JD, LAT, LON, 0.0)
    tau = TT.shear_along_bearing(st["s_NN"], st["s_EE"], st["s_NE"], 340.0)
    assert R.new_information_fraction(tau, [scal]) > 0.90


def test_axial_shared_variance_runs_on_the_doubled_angle():
    scal = ST.site_scalar(JD, LAT, LON, 0.0)
    st = TT.stress_tensor(JD, LAT, LON, 0.0)
    b = TT.principal_bearing(st["e_NN"], st["e_EE"], st["e_NE"])
    rep = R.axial_shared_variance(b, [scal])
    assert 0.0 <= rep["shared_variance"] <= 1.0
    assert len(rep["residual_cos"]) == len(b)


def test_orthogonality_report_ranks_by_novelty():
    scal = ST.site_scalar(JD, LAT, LON, 0.0)
    st = TT.stress_tensor(JD, LAT, LON, 0.0)
    cands = {
        "cfs_thrust": TT.coulomb(st, 250.0, 20.0, 90.0, 0.4)["coulomb_pa"],
        "rate": np.gradient(scal, T * 24.0),
        "shear": TT.shear_along_bearing(st["s_NN"], st["s_EE"], st["s_NE"], 340.0),
    }
    rep = R.orthogonality_report(cands, [scal])
    assert rep["ranked_by_novelty"][-1] == "cfs_thrust", rep["ranked_by_novelty"]
    assert rep["ranked_by_novelty"][0] in ("rate", "shear")
    assert rep["battery_effective_rank"]["n_nominal"] == 3
