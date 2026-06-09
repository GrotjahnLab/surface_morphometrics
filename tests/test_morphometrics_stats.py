"""Tests for the pure statistics helpers in morphometrics_stats (numpy/pandas)."""
import numpy as np

from surface_morphometrics import morphometrics_stats as ms


def test_weighted_avg_and_std_uniform_matches_numpy():
    v = np.array([1.0, 2.0, 3.0, 4.0])
    w = np.ones_like(v)
    avg, std = ms.weighted_avg_and_std(v, w)
    assert np.isclose(avg, v.mean())
    assert np.isclose(std, v.std())            # population std


def test_weighted_avg_and_std_weighting():
    v = np.array([0.0, 10.0])
    avg, std = ms.weighted_avg_and_std(v, np.array([3.0, 1.0]))
    assert np.isclose(avg, 2.5)                # (0*3 + 10*1)/4
    assert np.isclose(std, np.sqrt(np.average((v - 2.5) ** 2, weights=[3, 1])))


def test_weighted_median_example_from_docstring():
    # values [1,3,0], weights [0.1,0.3,0.6] -> weighted median 0
    assert ms.weighted_median([1, 3, 0], [0.1, 0.3, 0.6]) == 0
    # uniform weights -> behaves like an upper median on the sorted values
    assert ms.weighted_median([1, 2, 3], [1, 1, 1]) == 2
    # split exactly at 0.5 -> average of the two central values
    assert np.isclose(ms.weighted_median([0, 10], [1, 1]), 5.0)


def test_weighted_histogram_peak():
    # most weight near 8 -> peak bin centered there
    values = np.array([1.0, 8.0, 8.2, 8.4, 2.0])
    weights = np.array([1.0, 5.0, 5.0, 5.0, 1.0])
    peak = ms.weighted_histogram_peak(values, weights, bins=10, bin_range=(0, 10))
    assert 7.0 <= peak <= 9.0
