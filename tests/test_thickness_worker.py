"""Tests for the pure fitting/alignment helpers in _thickness_worker (numpy/scipy)."""
import numpy as np
import scipy.optimize as opt

from surface_morphometrics import _thickness_worker as tw


def _fit_bilayer(dat, x):
    """Mirror the dual-Gaussian offset fit used in fit_triangle_chunk_offsets."""
    mid = len(dat) // 2
    lm = np.argmin(dat[:mid])
    rm = np.argmin(dat[mid:]) + mid
    a, b = x[lm + 2:rm - 2], dat[lm + 2:rm - 2]
    center_seed, half_seed, _ = tw._seed_bilayer_center(a, b)
    window = tw.MAX_THICKNESS / 2.0 + 2.0
    m = (a >= center_seed - window) & (a <= center_seed + window)
    af, bf = (a[m], b[m]) if m.sum() >= 5 else (a, b)
    p0 = [0.02, 1.5, 0.02, 1.5, center_seed, half_seed, 0.0]
    bounds = ([0.005, 0.8, 0.005, 0.8, center_seed - 3.0, tw.MIN_THICKNESS / 2.0, -1],
              [0.04, 2.2, 0.04, 2.2, center_seed + 3.0, tw.MAX_THICKNESS / 2.0, 1])
    p, _ = opt.curve_fit(tw._dual_gaussian_centered, af, bf, p0, bounds=bounds)
    return p[4], p[5]  # center, half_sep


def test_dual_gaussian_centered_is_two_separated_peaks():
    x = np.linspace(-8, 8, 200)
    y = tw._dual_gaussian_centered(x, 1.0, 1.0, 1.0, 1.0, center=0.5, half_sep=2.0, o=0.1)
    # equals two Gaussians at center -/+ half_sep, plus offset
    expected = (tw._monogaussian(x, 1, 0.5 - 2.0, 1)
                + tw._monogaussian(x, 1, 0.5 + 2.0, 1) + 0.1)
    assert np.allclose(y, expected)


def test_seed_bilayer_center_finds_two_leaflets():
    x = np.linspace(-10, 10, 201)
    # leaflets at -0.6 and +3.0 -> center +1.2, half 1.8
    b = tw._monogaussian(x, 1.0, -0.6, 1.0) + tw._monogaussian(x, 0.8, 3.0, 1.0)
    center, half, n_peaks = tw._seed_bilayer_center(x, b)
    assert n_peaks == 2
    assert abs(center - 1.2) < 0.2 and abs(half - 1.8) < 0.2


def test_seed_bilayer_center_reports_single_unresolved_peak():
    x = np.linspace(-10, 10, 201)
    b = tw._monogaussian(x, 1.0, 1.0, 2.5)   # one broad, unresolved peak
    center, half, n_peaks = tw._seed_bilayer_center(x, b)
    assert n_peaks == 1                       # only one leaflet resolved
    assert abs(center - 1.0) < 0.2            # center seeded at the single peak


def test_seed_rejects_small_shoulders_as_single_peak():
    # A single central membrane peak flanked by two smaller protein-density
    # shoulders. The shoulders are local maxima but have tiny prominence, so they
    # must NOT be taken as leaflets -- the region is one (unresolved) membrane.
    x = np.linspace(-10, 10, 161)
    b = (tw._monogaussian(x, 0.020, 0.0, 1.5)      # central membrane
         + tw._monogaussian(x, 0.008, -4.0, 1.0)   # protein shoulder
         + tw._monogaussian(x, 0.008, 4.0, 1.0))   # protein shoulder
    center, half, n_resolved = tw._seed_bilayer_center(x, b)
    assert n_resolved == 1            # shoulders rejected, only the membrane counts
    assert abs(center) < 0.3          # centered on the membrane, not a shoulder/midpoint


def test_seed_keeps_asymmetric_bilayer_as_two_leaflets():
    # Genuinely asymmetric bilayer (one leaflet taller) must still resolve as two.
    x = np.linspace(-10, 10, 161)
    b = tw._monogaussian(x, 0.018, -0.6, 1.0) + tw._monogaussian(x, 0.015, 3.0, 1.0)
    center, half, n_resolved = tw._seed_bilayer_center(x, b)
    assert n_resolved == 2
    assert abs(center - 1.2) < 0.3


def test_single_resolved_peak_falls_back_to_single_gaussian():
    # A poorly-resolved region (one broad peak): the dual fit must NOT be forced;
    # n_resolved < 2 routes it to the single-Gaussian step, which still centers on
    # the peak (no spurious ~2 nm bilayer reported).
    x = np.linspace(-10, 10, 161)
    dat = 0.011 + tw._monogaussian(x, 0.02, 1.0, 2.5)
    mid = len(dat) // 2
    lm = np.argmin(dat[:mid]); rm = np.argmin(dat[mid:]) + mid
    a, b = x[lm + 2:rm - 2], dat[lm + 2:rm - 2]
    _, _, n_resolved = tw._seed_bilayer_center(a, b)
    assert n_resolved < 2                      # dual Gaussian is skipped
    # the single-Gaussian fallback recovers the peak center
    p_mono, _ = opt.curve_fit(tw._monogaussian_with_offset, a, b,
                              [0.03, 1.0, 2.5, 0],
                              bounds=([0.005, -6, 1.0, -1], [0.06, 6, 5.0, 1]))
    assert abs(p_mono[1] - 1.0) < 0.3


def test_fit_recovers_center_of_asymmetric_bilayer_with_confounder():
    # The reported failure: taller leaflet + nearby density made the old fit put
    # both Gaussians in one leaflet. The centered+windowed fit must recover the
    # true bilayer center (+1.2 nm), not a single leaflet position.
    x = np.linspace(-10, 10, 161)
    dat = (0.011
           + tw._monogaussian(x, 0.018, 1.2 - 1.8, 1.0)   # taller left leaflet
           + tw._monogaussian(x, 0.015, 1.2 + 1.8, 1.0)   # shorter right leaflet
           + tw._monogaussian(x, 0.014, 7.0, 1.0))        # confounding nearby density
    center, half = _fit_bilayer(dat, x)
    assert abs(center - 1.2) < 0.3                          # bilayer center, not a leaflet
    assert tw.MIN_THICKNESS <= 2 * half <= tw.MAX_THICKNESS  # peaks stay separated
    assert half >= tw.MIN_THICKNESS / 2.0


def test_fit_keeps_well_centered_bilayer_centered():
    x = np.linspace(-10, 10, 161)
    dat = 0.011 + tw._monogaussian(x, 0.017, -2.0, 1.1) + tw._monogaussian(x, 0.017, 2.0, 1.1)
    center, half = _fit_bilayer(dat, x)
    assert abs(center) < 0.2 and abs(2 * half - 4.0) < 0.5  # no spurious move


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
