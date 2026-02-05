"""
Worker functions for parallel thickness measurement.

This module is intentionally minimal to avoid importing heavy dependencies
(pycurv, graph_tool) in multiprocessing worker processes.
"""

import numpy as np
import scipy.optimize as opt


def _monogaussian(x, h, c, w):
    """Single gaussian function."""
    return h * np.exp(-(x - c)**2 / (2 * w**2))


def _dual_gaussian(x, h1, c1, w1, h2, c2, w2, o):
    """Dual gaussian function for curve fitting."""
    return _monogaussian(x, h1, c1, w1) + _monogaussian(x, h2, c2, w2) + o


def _monogaussian_with_offset(x, h, c, w, o):
    """Single Gaussian with offset for curve_fit."""
    return _monogaussian(x, h, c, w) + o


def _compute_r_squared(y_true, y_pred):
    """Compute R² (coefficient of determination)."""
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    if ss_tot == 0:
        return 0.0
    return 1 - (ss_res / ss_tot)


def _weighted_centroid(x, y, center_range=5.0):
    """
    Compute weighted centroid of the profile within center_range of x=0.

    More robust fallback when Gaussian fitting fails.
    """
    central_mask = (x > -center_range) & (x < center_range)
    if not np.any(central_mask):
        return 0.0

    x_central = x[central_mask]
    y_central = y[central_mask]

    # Shift weights to be positive
    weights = y_central - y_central.min()
    weight_sum = weights.sum()

    if weight_sum == 0:
        return 0.0

    return np.average(x_central, weights=weights)


def _xcorr_offset(profile, reference, x_positions):
    """
    Find the offset that best aligns profile with reference using cross-correlation.

    Uses normalized cross-correlation with sub-pixel interpolation for precision.

    Parameters
    ----------
    profile : np.ndarray
        The local density profile to align
    reference : np.ndarray
        The reference profile (e.g., global average)
    x_positions : np.ndarray
        The x positions corresponding to profile samples

    Returns
    -------
    float
        The optimal shift in the same units as x_positions
    """
    # Normalize both signals (zero mean, unit variance)
    prof_norm = profile - np.mean(profile)
    ref_norm = reference - np.mean(reference)

    prof_std = np.std(prof_norm)
    ref_std = np.std(ref_norm)

    if prof_std == 0 or ref_std == 0:
        return 0.0

    prof_norm = prof_norm / prof_std
    ref_norm = ref_norm / ref_std

    # Cross-correlation
    correlation = np.correlate(prof_norm, ref_norm, mode='full')

    # Find peak
    peak_idx = np.argmax(correlation)
    center_idx = len(reference) - 1

    # Sub-pixel interpolation using parabolic fit around peak
    if 1 <= peak_idx < len(correlation) - 1:
        y0 = correlation[peak_idx - 1]
        y1 = correlation[peak_idx]
        y2 = correlation[peak_idx + 1]
        # Parabolic interpolation
        denom = 2 * (2 * y1 - y0 - y2)
        if abs(denom) > 1e-10:
            sub_pixel = (y0 - y2) / denom
        else:
            sub_pixel = 0.0
    else:
        sub_pixel = 0.0

    # Convert to position offset
    sample_spacing = x_positions[1] - x_positions[0] if len(x_positions) > 1 else 1.0
    shift_samples = (peak_idx + sub_pixel) - center_idx
    offset = shift_samples * sample_spacing

    return offset


# Global variables for multiprocessing workers (set by initializer)
_worker_thickness_arr = None
_worker_distances = None
_worker_neighbor_indices = None
_worker_x = None
_worker_monolayer = False
_worker_reference = None  # Reference profile for cross-correlation
_worker_use_xcorr = False  # Whether to use cross-correlation method


def init_worker(thickness_arr, distances, neighbor_indices, x, monolayer=False,
                reference=None, use_xcorr=False):
    """Initialize worker process with shared data."""
    global _worker_thickness_arr, _worker_distances, _worker_neighbor_indices
    global _worker_x, _worker_monolayer, _worker_reference, _worker_use_xcorr
    _worker_thickness_arr = thickness_arr
    _worker_distances = distances
    _worker_neighbor_indices = neighbor_indices
    _worker_x = x
    _worker_monolayer = monolayer
    _worker_reference = reference
    _worker_use_xcorr = use_xcorr


def fit_triangle_chunk(indices):
    """
    Fit dual gaussian to a chunk of triangles' density profiles.

    Processes multiple triangles per call to reduce multiprocessing overhead.

    Returns
    -------
    list of tuples
        [(thickness, offset), ...] for each triangle in chunk
    """
    results = []
    for i in indices:
        # Get valid neighbors (those within distance bound)
        valid_mask = _worker_distances[i] != np.inf
        l = _worker_distances[i][valid_mask]
        neighbors = _worker_neighbor_indices[i][valid_mask]

        if len(neighbors) == 0:
            results.append((np.nan, 0))
            continue

        # Vectorized weight calculation
        weights = 1.0 / (1.0 + l)

        # Compute weighted average of density profiles
        dat = np.average(_worker_thickness_arr[neighbors], weights=weights, axis=0) * -1
        dat = dat - dat.min()
        dat_sum = dat.sum()
        if dat_sum == 0:
            results.append((np.nan, 0))
            continue
        dat = dat / (80/81 * dat_sum)

        # Find initial peak and set up fitting bounds
        ipk = _worker_x[np.argmax(dat)]
        p0 = [0.02, ipk-0.5, 1.5, 0.02, ipk+0.5, 1.5, 0]
        bounds = ([0.005, ipk-6, 0.8, 0.005, ipk-1.5, 0.8, -1],
                  [0.04, ipk+1.5, 2.2, 0.04, ipk+6, 2.2, 1])

        # Find minima to constrain fitting region
        mid = len(dat) // 2
        left_min = np.argmin(dat[:mid])
        right_min = np.argmin(dat[mid:]) + mid

        a = _worker_x[left_min+2:right_min-2]
        b = dat[left_min+2:right_min-2]

        if len(a) < 7:  # Need enough points for 7-parameter fit
            results.append((np.nan, 0))
            continue

        try:
            p3, _ = opt.curve_fit(_dual_gaussian, a, b, p0, bounds=bounds)
            thickness = np.abs(p3[4] - p3[1])
            offset = (p3[4] + p3[1]) / 2
            results.append((thickness, offset))
        except Exception:
            results.append((np.nan, 0))

    return results


def fit_triangle_chunk_offsets(indices):
    """
    Fit gaussian(s) to density profiles and return center offsets and sigmas.

    Two modes available (controlled by _worker_use_xcorr):

    Cross-correlation mode (_worker_use_xcorr=True):
        Uses normalized cross-correlation with reference profile to find
        optimal shift. Very robust to noise and shape variations.

    Gaussian fitting mode (_worker_use_xcorr=False):
        Uses adaptive model selection:
        1. Try dual Gaussian (bilayer model)
        2. If R² > 0.6 and thickness in 2-7 nm range, use it
        3. Otherwise try single Gaussian
        4. If that also fails, fall back to weighted centroid

    When _worker_monolayer is True in Gaussian mode, skips dual Gaussian
    and goes straight to single Gaussian fitting.

    Returns
    -------
    list of tuples
        [(offset, sigma1, sigma2, fit_method), ...] for each triangle in chunk.
        fit_method: 0=failed, 1=dual_gaussian, 2=single_gaussian, 3=centroid, 4=xcorr
    """
    # Fit quality thresholds
    R2_THRESHOLD = 0.6
    MIN_THICKNESS = 2.0  # nm
    MAX_THICKNESS = 7.0  # nm

    results = []
    for i in indices:
        # Get valid neighbors (those within distance bound)
        valid_mask = _worker_distances[i] != np.inf
        l = _worker_distances[i][valid_mask]
        neighbors = _worker_neighbor_indices[i][valid_mask]

        if len(neighbors) == 0:
            results.append((0, np.nan, np.nan, 0))
            continue

        # Vectorized weight calculation
        weights = 1.0 / (1.0 + l)

        # Compute weighted average of density profiles
        dat = np.average(_worker_thickness_arr[neighbors], weights=weights, axis=0) * -1
        dat = dat - dat.min()
        dat_sum = dat.sum()
        if dat_sum == 0:
            results.append((0, np.nan, np.nan, 0))
            continue
        dat = dat / (80/81 * dat_sum)

        # Cross-correlation mode - use template matching with reference profile
        if _worker_use_xcorr and _worker_reference is not None:
            offset = _xcorr_offset(dat, _worker_reference, _worker_x)
            results.append((offset, np.nan, np.nan, 4))
            continue

        # Find initial peak and set up fitting bounds
        ipk = _worker_x[np.argmax(dat)]

        # Find minima to constrain fitting region
        mid = len(dat) // 2
        left_min = np.argmin(dat[:mid])
        right_min = np.argmin(dat[mid:]) + mid

        a = _worker_x[left_min+2:right_min-2]
        b = dat[left_min+2:right_min-2]

        if len(a) < 5:  # Need enough points for fitting
            # Fall back to centroid
            offset = _weighted_centroid(_worker_x, dat)
            results.append((offset, np.nan, np.nan, 3))
            continue

        # Track best result
        best_offset = None
        best_sigma1 = np.nan
        best_sigma2 = np.nan
        best_method = 0

        # Step 1: Try dual Gaussian (unless monolayer mode)
        if not _worker_monolayer:
            try:
                p0 = [0.02, ipk-0.5, 1.5, 0.02, ipk+0.5, 1.5, 0]
                bounds = ([0.005, ipk-6, 0.8, 0.005, ipk-1.5, 0.8, -1],
                          [0.04, ipk+1.5, 2.2, 0.04, ipk+6, 2.2, 1])
                p3, _ = opt.curve_fit(_dual_gaussian, a, b, p0, bounds=bounds)

                # Validate fit
                y_pred = _dual_gaussian(a, *p3)
                r2 = _compute_r_squared(b, y_pred)
                thickness = np.abs(p3[4] - p3[1])

                if r2 > R2_THRESHOLD and MIN_THICKNESS <= thickness <= MAX_THICKNESS:
                    best_offset = (p3[1] + p3[4]) / 2
                    best_sigma1 = p3[2]
                    best_sigma2 = p3[5]
                    best_method = 1
            except Exception:
                pass

        # Step 2: Try single Gaussian if dual didn't work
        if best_method == 0:
            try:
                p0_mono = [0.03, ipk, 2.5, 0]
                bounds_mono = ([0.005, ipk-6, 1.0, -1],
                               [0.06, ipk+6, 5.0, 1])
                p_mono, _ = opt.curve_fit(_monogaussian_with_offset, a, b, p0_mono, bounds=bounds_mono)

                # Validate fit
                y_pred = _monogaussian_with_offset(a, *p_mono)
                r2 = _compute_r_squared(b, y_pred)

                if r2 > R2_THRESHOLD:
                    best_offset = p_mono[1]
                    best_sigma1 = p_mono[2]
                    best_sigma2 = np.nan
                    best_method = 2
            except Exception:
                pass

        # Step 3: Fall back to weighted centroid
        if best_method == 0:
            best_offset = _weighted_centroid(_worker_x, dat)
            best_sigma1 = np.nan
            best_sigma2 = np.nan
            best_method = 3

        results.append((best_offset, best_sigma1, best_sigma2, best_method))

    return results
