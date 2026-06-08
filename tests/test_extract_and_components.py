"""Tests for the extractor mask and connected-component numbering (numpy only)."""
import numpy as np

from surface_morphometrics import extract_patches as ep
from surface_morphometrics import label_connected_components as lcc


def test_mask_from_property():
    v = np.array([1.0, 5.0, np.nan, 30.0, 31.0, 0.0])
    assert ep.mask_from_property(v, vmax=30).tolist() == [True, True, False, True, False, True]
    assert ep.mask_from_property(v, vmin=5).tolist() == [False, True, False, True, True, False]
    assert ep.mask_from_property(v, vmin=1, vmax=30).tolist() == [True, True, False, True, False, False]
    # NaN always excluded; open bounds keep everything else
    assert ep.mask_from_property(v).tolist() == [True, True, False, True, True, True]


def test_label_component_numbers_ordering_and_min_size():
    ci = np.array([0, 0, 0, 1, 1, 2, 3, 3, 3, 3])  # sizes: c0=3, c1=2, c2=1, c3=4
    # 1-based by descending size: c3->1, c0->2, c1->3, c2->4
    assert lcc.label_component_numbers(ci, min_size=0).tolist() == [2, 2, 2, 3, 3, 4, 1, 1, 1, 1]
    # min_size drops the size-1 component to 0
    assert lcc.label_component_numbers(ci, min_size=2).tolist() == [2, 2, 2, 3, 3, 0, 1, 1, 1, 1]
