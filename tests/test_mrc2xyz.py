"""Tests for mrc2xyz.mrc_to_xyz (mrcfile/numpy/pandas; no pymeshlab)."""
import os
import tempfile

import numpy as np
import mrcfile
import pandas as pd

from surface_morphometrics import mrc2xyz


def _write_labeled_mrc(path, voxel_size=10.0):
    # (z, y, x) volume; label 2 at three known voxels, label 1 elsewhere-ish
    data = np.zeros((4, 5, 6), dtype=np.int8)
    coords_zyx = [(1, 2, 3), (0, 0, 0), (3, 4, 5)]
    for z, y, x in coords_zyx:
        data[z, y, x] = 2
    with mrcfile.new(path, overwrite=True) as m:
        m.set_data(data)
        m.voxel_size = voxel_size
    return coords_zyx


def test_mrc_to_xyz_writes_scaled_points():
    with tempfile.TemporaryDirectory() as tmp:
        mrc = os.path.join(tmp, "seg.mrc")
        xyz = os.path.join(tmp, "out.xyz")
        coords_zyx = _write_labeled_mrc(mrc, voxel_size=10.0)  # -> 1 nm voxels

        ret = mrc2xyz.mrc_to_xyz(mrc, xyz, label=2, angstrom=False)
        assert ret == 0 and os.path.exists(xyz)

        pts = pd.read_csv(xyz, sep=" ", header=None).to_numpy()
        assert pts.shape == (3, 3)                      # 3 labeled voxels, x y z
        # voxel (z,y,x) -> row (x,y,z) * 1 nm; origin defaults to 0
        expected = {(x, y, z) for (z, y, x) in coords_zyx}
        got = {tuple(row) for row in pts.astype(int)}
        assert got == expected


def test_mrc_to_xyz_missing_label_returns_error():
    with tempfile.TemporaryDirectory() as tmp:
        mrc = os.path.join(tmp, "seg.mrc")
        xyz = os.path.join(tmp, "out.xyz")
        _write_labeled_mrc(mrc)
        assert mrc2xyz.mrc_to_xyz(mrc, xyz, label=99, angstrom=False) == 1
        assert not os.path.exists(xyz)
