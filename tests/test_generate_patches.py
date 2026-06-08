"""Tests for the pure-geometry helpers in generate_patches (no pycurv needed)."""
import numpy as np
from scipy.spatial.distance import pdist

from surface_morphometrics import generate_patches as gp


def test_nearest_triangle():
    tri = np.array([[float(x), 0, 0] for x in range(11)])
    part = np.array([[0.2, 0, 0], [9.6, 0, 0]])
    d, i = gp.nearest_triangle(tri, part)
    assert i.tolist() == [0, 10]
    assert np.allclose(d, [0.2, 0.4])


def test_assign_patches_overlap_nearest_center():
    # Points along a line; two centers at the ends, radius covers the middle.
    tri = np.array([[float(x), 0, 0] for x in range(11)])
    res = gp.assign_patches(tri, [0, 10], [1, 2], patch_radius=8)
    num = res["number"]
    assert num[4] == 1          # closer to center 0
    assert num[6] == 2          # closer to center 1
    assert num[5] == 2          # tie -> last-added wins
    assert num[0] == 1 and num[10] == 2
    assert abs(res["center_distance"][4] - 4.0) < 1e-9


def test_assign_patches_protein_distance():
    tri = np.array([[float(x), 0, 0] for x in range(11)])
    pxyz = np.array([[0, 0, 2.0], [10, 0, 2.0]])  # particles 2 nm off membrane
    res = gp.assign_patches(tri, [0, 10], [1, 2], 8, particle_xyz=pxyz)
    assert abs(res["protein_distance"][0] - 2.0) < 1e-9
    # triangles in no patch have NaN center distance
    assert np.isnan(res["center_distance"][res["number"] == 0]).all()


def test_choose_random_centers_min_distance_and_seed():
    grid = np.array([[float(x), float(y), 0] for x in range(10) for y in range(10)])
    chosen = gp.choose_random_centers(grid, 5, min_distance=3.0, rng=np.random.default_rng(0))
    assert len(chosen) == 5
    assert pdist(grid[chosen]).min() >= 3.0
    # reproducible for the same seed
    assert chosen == gp.choose_random_centers(grid, 5, 3.0, np.random.default_rng(0))


def test_tomo_name_matches():
    m = gp.tomo_name_matches
    assert m("/data/YTC042_2_lam10_ts_004.mrc_6.65Apx.mrc", "YTC042_2_lam10_ts_004")
    assert m("TS_004.mrc", "TS_004_bin4")          # reverse substring
    assert m("TS_004_bin4.mrc", "TS_004")          # forward substring
    assert not m("/data/TS_005.mrc", "TS_004")     # different tomogram
    assert not m("", "TS_004") and not m("TS_004", "")  # empty guards
