import numpy as np

from engine.grid import Grid, day_index


def test_cell_index_roundtrip_and_outside_domain():
    g = Grid(1.0, 1.0)
    lat = np.array([0.2, 10.7, -5.4, 33.9])
    lon = np.array([0.1, 20.3, -100.9, -118.2])
    g.build_domain(lat, lon)
    assert g.n_cells == 4
    idx = g.cell_index(lat, lon)
    assert sorted(idx.tolist()) == [0, 1, 2, 3]
    # same cell, different point inside it
    assert g.cell_index(np.array([0.9]), np.array([0.9]))[0] == g.cell_index(
        np.array([0.2]), np.array([0.1]))[0]
    # a point in an inactive cell is outside the domain
    assert g.cell_index(np.array([70.0]), np.array([70.0]))[0] == -1


def test_coarser_grid_merges_cells():
    lat = np.array([0.2, 0.9, 1.4])
    lon = np.array([0.2, 0.9, 1.4])
    assert Grid(1.0, 1.0).build_domain(lat, lon) == 2
    assert Grid(2.0, 2.0).build_domain(lat, lon) == 1


def test_distance_and_azimuth():
    g = Grid(1.0, 1.0)
    g.build_domain(np.array([0.5, 0.5, 1.5]), np.array([0.5, 1.5, 0.5]))
    d = g.pair_distance_km()
    assert np.allclose(np.diag(d), 0.0, atol=1e-3)
    assert np.allclose(d, d.T, atol=1e-2)
    # one degree of latitude near the equator is ~111 km
    off = d[d > 0]
    assert 100 < off.min() < 125
    az = g.pair_azimuth_rad()
    # cells sorted by (lat, lon): 0=(0.5,0.5) 1=(0.5,1.5) 2=(1.5,0.5)
    assert abs(np.rad2deg(az[0, 1]) - 90.0) < 1.0    # due east
    assert abs(np.rad2deg(az[0, 2]) - 0.0) < 1.0     # due north


def test_day_index_clips():
    assert day_index(np.array([-3.0, 0.4, 9.9, 40.0]), 10).tolist() == [0, 0, 9, 9]
