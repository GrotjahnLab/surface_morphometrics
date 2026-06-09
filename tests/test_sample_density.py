"""Tests for the pure parts of sample_density (numpy/scipy/mrcfile; no graph-tool)."""
import os
import tempfile

import numpy as np
import mrcfile

from surface_morphometrics import sample_density as sd


def test_interpolate_samples_along_normals():
    # Volume whose value equals its first-axis index: data[i,j,k] = i
    n = 10
    data = np.broadcast_to(np.arange(n)[:, None, None], (n, n, n)).astype(float).copy()
    data_matrix = (np.arange(n), np.arange(n), np.arange(n))

    # One triangle at (5,5,5) with unit normal along axis 0 (already in voxel steps)
    xyz = np.array([[5.0], [5.0], [5.0]])
    n_v = np.array([[1.0], [0.0], [0.0]])

    out = sd.interpolate(data, data_matrix, xyz, n_v, sample_spacing=1.0, scan_range=2)
    # samples = [-2,-1,0,1,2] along axis 0 -> values 3..7
    assert out.shape == (1, 5)
    assert np.allclose(out[0], [3, 4, 5, 6, 7])


def test_load_mrc_units_and_axis_order():
    data = np.arange(2 * 3 * 4, dtype=np.float32).reshape(2, 3, 4)  # (z, y, x)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "vol.mrc")
        with mrcfile.new(path, overwrite=True) as m:
            m.set_data(data)
            m.voxel_size = 10.0  # Angstrom -> 1.0 nm
        arr, data_matrix, voxsize, origin = sd.load_mrc(path)   # nm by default
        assert np.isclose(voxsize, 1.0)                          # 10 A -> 1 nm
        assert arr.shape == (4, 3, 2)                            # swapaxes(0, 2)
        assert len(data_matrix) == 3 and len(origin) == 3
        # angstrom mode keeps raw voxel size
        _, _, voxsize_a, _ = sd.load_mrc(path, angstroms=True)
        assert np.isclose(voxsize_a, 10.0)
