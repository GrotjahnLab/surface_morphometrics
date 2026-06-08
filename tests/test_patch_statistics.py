"""Tests for the area-weighted per-region aggregator (pandas only)."""
import numpy as np
import pandas as pd

from surface_morphometrics import patch_statistics as ps


def _df():
    return pd.DataFrame({
        "patch_number": [1, 1, 2, 2, 0],
        "curvedness_VV": [0.1, 0.3, 0.2, 0.0, 0.9],
        "thickness": [3.0, 4.0, 3.5, 3.5, 9.0],
        "area": [1.0, 1.0, 2.0, 2.0, 1.0],
    })


def test_aggregate_regions_area_weighted():
    df = pd.DataFrame({"component_number": [1, 1], "thickness": [2.0, 4.0], "area": [3.0, 1.0]})
    r = ps.aggregate_regions(df, "component_number", ["thickness"])
    # weighted mean = (2*3 + 4*1) / 4 = 2.5
    assert abs(r.iloc[0]["thickness_mean"] - 2.5) < 1e-9
    assert r.iloc[0]["total_area"] == 4.0


def test_aggregate_regions_drops_zero_and_nan():
    df = pd.DataFrame({
        "patch_number": [1, 1, 1],
        "curvedness_VV": [0.1, 0.3, 0.0],   # the 0 is dropped
        "thickness": [3.0, 4.0, np.nan],    # the NaN is dropped
        "area": [1.0, 1.0, 1.0],
    })
    r = ps.aggregate_regions(df, "patch_number", ["curvedness_VV", "thickness"]).iloc[0]
    assert abs(r["curvedness_VV_mean"] - 0.2) < 1e-9
    assert abs(r["thickness_mean"] - 3.5) < 1e-9
    assert r["n_triangles"] == 3            # n_triangles counts all region rows


def test_aggregate_regions_excludes_region_zero():
    r = ps.aggregate_regions(_df(), "patch_number", ["thickness"])
    assert set(r["region_id"]) == {1, 2}    # region 0 excluded
