"""Tests for the OBJ/MTL/colormap writer (numpy + matplotlib; no vtk needed)."""
import os
import tempfile

import numpy as np

from surface_morphometrics import export_obj as eo


def _mesh():
    pts = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0]], float)
    faces = np.array([[0, 1, 2], [0, 2, 3], [1, 4, 2]])
    return pts, faces


def test_write_obj_mtl_creates_files_and_uvs():
    pts, faces = _mesh()
    vals = np.array([3.0, 4.0, 5.0])
    out = os.path.join(tempfile.mkdtemp(), "s_thickness")
    lo, hi = eo.write_obj_mtl(out, pts, faces, vals, "thickness", "viridis",
                              vmin=3.0, vmax=5.0, nan_color="lightgrey")
    for ext in (".obj", ".mtl", ".png"):
        assert os.path.exists(out + ext)
    obj = open(out + ".obj").read()
    assert "mtllib s_thickness.mtl" in obj and "usemtl quant_thickness" in obj
    assert obj.count("\nv ") == len(pts)
    assert obj.count("\nf ") == len(faces)
    assert "map_Kd s_thickness.png" in open(out + ".mtl").read()


def test_write_obj_mtl_nan_swatch():
    import matplotlib.image as mpimg
    from matplotlib.colors import to_rgb
    pts, faces = _mesh()
    vals = np.array([3.0, np.nan, 5.0])
    out = os.path.join(tempfile.mkdtemp(), "s_thickness")
    eo.write_obj_mtl(out, pts, faces, vals, "thickness", "viridis",
                     vmin=3.0, vmax=5.0, nan_color="lightgrey")
    img = mpimg.imread(out + ".png")
    assert img.shape[1] == eo._RAMP_W + eo._NAN_W          # ramp + swatch
    assert np.allclose(img[0, eo._RAMP_W + 4, :3], to_rgb("lightgrey"), atol=0.02)
    # the NaN face's u samples the swatch
    u = [float(l.split()[1]) for l in open(out + ".obj") if l.startswith("vt ")]
    assert abs(u[1] - (eo._RAMP_W + eo._NAN_W / 2.0) / (eo._RAMP_W + eo._NAN_W)) < 1e-6


def test_write_obj_mtl_nan_disabled():
    import matplotlib.image as mpimg
    pts, faces = _mesh()
    vals = np.array([3.0, np.nan, 5.0])
    out = os.path.join(tempfile.mkdtemp(), "s_thickness")
    eo.write_obj_mtl(out, pts, faces, vals, "thickness", "viridis",
                     vmin=3.0, vmax=5.0, nan_color=None)
    img = mpimg.imread(out + ".png")
    assert img.shape[1] == eo._RAMP_W                       # no swatch
    u = [float(l.split()[1]) for l in open(out + ".obj") if l.startswith("vt ")]
    assert abs(u[1] - 0.5 / eo._RAMP_W) < 1e-6              # NaN -> low end
