"""Tests for the pure fitting/alignment helpers in _thickness_worker (numpy/scipy)."""
import numpy as np

from surface_morphometrics import _thickness_worker as tw


def test_monogaussian_peak_and_symmetry():
    x = np.linspace(-10, 10, 401)
    y = tw._monogaussian(x, h=2.0, c=1.0, w=1.5)
    assert np.isclose(y.max(), 2.0, atol=1e-3)        # peak height
    assert np.isclose(x[np.argmax(y)], 1.0, atol=0.1)  # peak position
    # symmetric about the center
    assert np.isclose(tw._monogaussian(0.0, 2, 1, 1.5), tw._monogaussian(2.0, 2, 1, 1.5))


def test_dual_gaussian_is_sum_plus_offset():
    x = np.linspace(-5, 5, 50)
    a = tw._monogaussian(x, 1, -1, 1)
    b = tw._monogaussian(x, 1, 1, 1)
    assert np.allclose(tw._dual_gaussian(x, 1, -1, 1, 1, 1, 1, 0.3), a + b + 0.3)
    assert np.allclose(tw._monogaussian_with_offset(x, 1, -1, 1, 0.3), a + 0.3)


def test_compute_r_squared():
    y = np.array([1.0, 2.0, 3.0, 4.0])
    assert tw._compute_r_squared(y, y) == 1.0          # perfect fit
    # predicting the mean everywhere -> R^2 == 0
    assert np.isclose(tw._compute_r_squared(y, np.full_like(y, y.mean())), 0.0)
    # constant truth -> defined as 0
    assert tw._compute_r_squared(np.ones(4), np.ones(4)) == 0.0


def test_weighted_centroid_recovers_peak_position():
    x = np.linspace(-10, 10, 201)
    y = tw._monogaussian(x, 1.0, 2.0, 1.0)             # peak at +2 (within center_range)
    assert np.isclose(tw._weighted_centroid(x, y, center_range=5.0), 2.0, atol=0.2)
    # no samples in range -> 0.0
    assert tw._weighted_centroid(np.array([20.0, 30.0]), np.array([1.0, 1.0])) == 0.0


def test_xcorr_offset_recovers_known_shift():
    x = np.linspace(-10, 10, 201)
    ref = tw._monogaussian(x, 1.0, 0.0, 1.5)
    shifted = tw._monogaussian(x, 1.0, 1.0, 1.5)       # shifted +1 nm
    off = tw._xcorr_offset(shifted, ref, x)
    assert np.isclose(off, 1.0, atol=0.15)
    # flat profile -> zero offset (no information)
    assert tw._xcorr_offset(np.ones_like(x), ref, x) == 0.0


def test_xcorr_offsets_batch_matches_known_shifts():
    x = np.linspace(-10, 10, 201)
    ref = tw._monogaussian(x, 1.0, 0.0, 1.5)
    shifts = [-2.0, 0.0, 1.5]
    profiles = np.vstack([tw._monogaussian(x, 1.0, s, 1.5) for s in shifts])
    out = tw._xcorr_offsets_batch(profiles, ref, x)
    assert np.allclose(out, shifts, atol=0.2)
