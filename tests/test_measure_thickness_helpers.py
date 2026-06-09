"""Tests for the pure fitting helpers in measure_thickness (no pycurv/graph-tool)."""
import numpy as np

from surface_morphometrics import measure_thickness as mt


def test_find_mins_locates_left_and_right_minima():
    x = np.linspace(-10, 10, 41)
    # two wells: minima near -5 and +5
    y = -np.exp(-(x + 5) ** 2) - np.exp(-(x - 5) ** 2) + 1.0
    left, right, _ = mt.find_mins(y)
    mid = len(y) // 2
    assert left < mid <= right
    assert abs(x[left] - (-5)) < 1.0
    assert abs(x[right] - 5) < 1.0


def test_gauss_is_unit_area_normal():
    x = np.linspace(-20, 20, 4001)
    y = mt.gauss(x, [0.0, 2.0])
    assert np.isclose(np.trapz(y, x), 1.0, atol=1e-3)   # normalized to area 1
    assert np.isclose(x[np.argmax(y)], 0.0, atol=0.05)


def test_monogaussian_and_dual_gaussian():
    x = np.linspace(-5, 5, 100)
    assert np.isclose(mt.monogaussian(1.0, 2.0, 1.0, 1.5), 2.0)  # value at center == height
    expected = mt.monogaussian(x, 1, -1, 1) + mt.monogaussian(x, 1, 1, 1) + 0.5
    assert np.allclose(mt.dual_gaussian(x, 1, -1, 1, 1, 1, 1, 0.5), expected)


def test_sinc_and_dual_sinc():
    # squared sinc peaks at mu (limit -> A)
    assert np.isclose(mt.sinc(1.0 + 1e-6, 3.0, 1.0, 0.5), 3.0, atol=1e-3)
    x = np.array([0.3, 0.7])
    p = (3.0, 1.0, 0.5, 2.0, -1.0, 0.5, 0.25)
    assert np.allclose(mt.dual_sinc(x, p), mt.sinc(x, *p[0:3]) + mt.sinc(x, *p[3:6]) + p[6])


def test_fit_gaussian_recovers_fwhm():
    x = np.linspace(-10, 10, 41)
    sigma = 2.0
    y = mt.gauss(x, [0.0, sigma])
    _, fwhm = mt.fit_gaussian(x, y, skipedge=3)
    assert np.isclose(fwhm, 2 * np.sqrt(2 * np.log(2)) * sigma, rtol=0.05)


def test_peak_fit_and_find_two_peaks():
    x = np.linspace(-10, 10, 201)
    # single peak -> peak_fit returns a positive width
    y1 = np.exp(-(x ** 2))
    width, height, h0, h1 = mt.peak_fit(x, y1)
    assert width > 0 and h0 < h1

    # two resolved peaks -> find_two_peaks returns their separation
    y2 = np.exp(-((x - 3) ** 2)) + np.exp(-((x + 3) ** 2))
    sep, p1, p2 = mt.find_two_peaks(x, y2)
    assert np.isclose(sep, 6.0, atol=0.3)
    # fewer than two peaks -> zeros
    assert mt.find_two_peaks(x, y1) == (0, 0, 0)


def test_func_sum_of_exponentials():
    x = np.array([0.0, 1.0])
    # a1=2,b1=1, a2=3,b2=0  -> 2 e^-x + 3
    out = mt.func(x, 2.0, 1.0, 3.0, 0.0)
    assert np.allclose(out, [2.0 + 3.0, 2.0 * np.exp(-1) + 3.0])
