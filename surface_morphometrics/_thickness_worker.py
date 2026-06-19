"""
Worker functions for parallel thickness measurement.

This module is intentionally minimal to avoid importing heavy dependencies
(pycurv, graph_tool) in multiprocessing worker processes.
"""

import numpy as np
import scipy.optimize as opt

# Bilayer fit-quality thresholds, shared by the per-triangle and global fits.
R2_THRESHOLD = 0.6
MIN_THICKNESS = 2.0   # nm (minimum accepted peak-to-peak leaflet separation)
MAX_THICKNESS = 7.0   # nm (maximum accepted peak-to-peak leaflet separation)


def _monogaussian(x, h, c, w):
    """Single gaussian function."""
    return h * np.exp(-(x - c)**2 / (2 * w**2))


def _dual_gaussian(x, h1, c1, w1, h2, c2, w2, o):
    """Dual gaussian function for curve fitting."""
    return _monogaussian(x, h1, c1, w1) + _monogaussian(x, h2, c2, w2) + o


def _monogaussian_with_offset(x, h, c, w, o):
    """Single Gaussian with offset for curve_fit."""
    return _monogaussian(x, h, c, w) + o


def _dual_gaussian_centered(x, h1, w1, h2, w2, center, half_sep, o):
    """Dual Gaussian parameterized by the bilayer center and half-separation.

    The two leaflet peaks sit at ``center - half_sep`` and ``center + half_sep``,
    so bounding ``half_sep`` to a physical range forces one Gaussian into each
    leaflet -- the fit can never collapse both peaks into the same leaflet (the
    failure mode of seeding both Gaussians at the global-max peak). ``center`` is
    the membrane center (the refinement offset); ``2 * half_sep`` is the thickness.
    """
    return (_monogaussian(x, h1, center - half_sep, w1)
            + _monogaussian(x, h2, center + half_sep, w2) + o)


def _seed_bilayer_center(a, b, default_center=0.0):
    """Estimate (center, half_sep) of a bilayer profile region for fit seeding.

    Locates the two tallest peaks in the profile (the two leaflets) and seeds the
    center at their midpoint and half_sep at half their separation. Falls back to
    ``default_center`` if fewer than two peaks are found. ``half_sep`` is clamped
    to the physical thickness range.
    """
    from scipy.signal import find_peaks

    a = np.asarray(a)
    b = np.asarray(b)
    peaks, _ = find_peaks(b)
    if len(peaks) >= 2:
        # the two tallest peaks are the leaflets
        tallest = peaks[np.argsort(b[peaks])[-2:]]
        p_left, p_right = sorted(a[tallest])
        center = 0.5 * (p_left + p_right)
        half = 0.5 * (p_right - p_left)
    elif len(peaks) == 1:
        center, half = float(a[peaks[0]]), MIN_THICKNESS / 2.0
    else:
        center, half = default_center, MIN_THICKNESS / 2.0
    half = float(np.clip(half, MIN_THICKNESS / 2.0, MAX_THICKNESS / 2.0))
    return float(center), half


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


def _xcorr_offsets_batch(profiles, reference, x_positions):
    """
    Vectorized cross-correlation of N profiles against a single reference.

    Uses FFT-based batch correlation instead of N separate np.correlate calls,
    which is ~10x faster for typical neighborhood sizes (~50 profiles).

    Returns shifts in the same units as x_positions (one value per profile).
    """
    N, L = profiles.shape
    # Normalize each profile (zero mean, unit variance)
    means = profiles.mean(axis=1, keepdims=True)
    prof_norm = profiles - means
    stds = prof_norm.std(axis=1, keepdims=True)
    nonzero = stds.flatten() > 0
    prof_norm[nonzero] /= stds[nonzero]

    ref_norm = reference - reference.mean()
    ref_std = reference.std()
    if ref_std == 0:
        return np.zeros(N)
    ref_norm = ref_norm / ref_std

    # Batch cross-correlation via FFT: correlate(a, b) = IFFT(FFT(a) * conj(FFT(b)))
    # Use next power-of-2 for FFT efficiency
    n_full = 2 * L - 1
    n_fft = 1 << (n_full - 1).bit_length()
    P = np.fft.rfft(prof_norm, n=n_fft, axis=1)   # (N, n_fft//2+1)
    R = np.fft.rfft(ref_norm, n=n_fft)             # (n_fft//2+1,)
    corr = np.fft.irfft(P * np.conj(R)[np.newaxis, :], n=n_fft, axis=1)
    # Full-length slice (same as mode='full' output)
    corr = np.concatenate([corr[:, -(L - 1):], corr[:, :L]], axis=1)  # (N, 2L-1)

    # Find peak and do parabolic sub-pixel refinement (vectorized)
    peak_idx = np.argmax(corr, axis=1)  # (N,)
    center_idx = L - 1
    sub_pixel = np.zeros(N)
    valid_pk = (peak_idx >= 1) & (peak_idx < n_full - 1)
    if np.any(valid_pk):
        rows = np.where(valid_pk)[0]
        pi = peak_idx[rows]
        y0 = corr[rows, pi - 1]
        y1 = corr[rows, pi]
        y2 = corr[rows, pi + 1]
        denom = 2 * (2 * y1 - y0 - y2)
        good = np.abs(denom) > 1e-10
        sub_pixel[rows[good]] = (y0[good] - y2[good]) / denom[good]

    sample_spacing = x_positions[1] - x_positions[0] if len(x_positions) > 1 else 1.0
    shift_samples = (peak_idx + sub_pixel) - center_idx
    return shift_samples * sample_spacing


# Global variables for local thickness computation workers
_lt_value_array = None
_lt_distances = None
_lt_neighbors = None
_lt_x_positions = None


def init_local_thickness_worker(value_array, distances, neighbors, x_positions):
    """Initialize worker process with shared data for local thickness computation."""
    global _lt_value_array, _lt_distances, _lt_neighbors, _lt_x_positions
    _lt_value_array = value_array
    _lt_distances = distances
    _lt_neighbors = neighbors
    _lt_x_positions = x_positions


def compute_thickness_chunk(indices):
    """
    Compute local thickness for a chunk of triangle indices.

    Returns list of thickness values (np.nan for failed fits).
    """
    results = []
    for i in indices:
        valid_mask = _lt_distances[i] != np.inf
        if not np.any(valid_mask):
            results.append(np.nan)
            continue

        l = _lt_distances[i][valid_mask]
        neighbors = _lt_neighbors[i][valid_mask]
        weights = 1.0 / (1.0 + l)

        dat = np.average(_lt_value_array[neighbors], weights=weights, axis=0) * -1
        dat = dat - dat.min()
        dat_sum = dat.sum()
        if dat_sum == 0:
            results.append(np.nan)
            continue
        dat = dat / (80/81 * dat_sum)

        try:
            ipk = _lt_x_positions[np.argmax(dat)]
            mid = len(dat) // 2
            left_min = np.argmin(dat[:mid])
            right_min = np.argmin(dat[mid:]) + mid

            a = _lt_x_positions[left_min+2:right_min-2]
            b = dat[left_min+2:right_min-2]

            if len(a) >= 7:
                p0 = [0.02, ipk-0.5, 1.5, 0.02, ipk+0.5, 1.5, 0]
                bounds = ([0.005, ipk-6, 0.8, 0.005, ipk-1.5, 0.8, -1],
                          [0.04, ipk+1.5, 2.2, 0.04, ipk+6, 2.2, 1])
                popt, _ = opt.curve_fit(_dual_gaussian, a, b, p0, bounds=bounds)
                thickness = np.abs(popt[4] - popt[1])
                if 2.0 <= thickness <= 7.0:
                    results.append(thickness)
                else:
                    results.append(np.nan)
            else:
                results.append(np.nan)
        except Exception:
            results.append(np.nan)

    return results


# Global variables for offset fitting workers (set by initializer)
_worker_thickness_arr = None
_worker_distances = None
_worker_neighbor_indices = None
_worker_x = None
_worker_monolayer = False
_worker_reference = None  # Reference profile for cross-correlation
_worker_use_xcorr = False  # Whether to use cross-correlation method
_worker_global_fit_params = None  # (c1, w1, c2, w2) from global dual Gaussian fit
_worker_rolling_xcorr = False  # Whether fit_triangle_chunk aligns neighbors before averaging
_worker_raw_average = True    # If True: average raw profiles then normalize (original behavior)


def init_worker(thickness_arr, distances, neighbor_indices, x, monolayer=False,
                reference=None, use_xcorr=False, global_fit_params=None,
                rolling_xcorr=False, raw_average=True):
    """Initialize worker process with shared data."""
    global _worker_thickness_arr, _worker_distances, _worker_neighbor_indices
    global _worker_x, _worker_monolayer, _worker_reference, _worker_use_xcorr
    global _worker_global_fit_params, _worker_rolling_xcorr, _worker_raw_average
    _worker_thickness_arr = thickness_arr
    _worker_distances = distances
    _worker_neighbor_indices = neighbor_indices
    _worker_x = x
    _worker_monolayer = monolayer
    _worker_reference = reference
    _worker_use_xcorr = use_xcorr
    _worker_global_fit_params = global_fit_params
    _worker_rolling_xcorr = rolling_xcorr
    _worker_raw_average = raw_average


def fit_triangle_chunk(indices):
    """
    Fit dual gaussian to a chunk of triangles' density profiles.

    Processes multiple triangles per call to reduce multiprocessing overhead.

    Returns
    -------
    list of tuples
        [(thickness, offset), ...] for each triangle in chunk
    """
    sample_spacing = _worker_x[1] - _worker_x[0] if len(_worker_x) > 1 else 1.0
    max_shift_samples = max(1, int(round(2.0 / sample_spacing)))

    results = []
    for i in indices:
        # Get valid neighbors (those within distance bound)
        valid_mask = _worker_distances[i] != np.inf
        l = _worker_distances[i][valid_mask]
        neighbors = _worker_neighbor_indices[i][valid_mask]

        if len(neighbors) == 0:
            results.append((np.nan, 0))
            continue

        weights = 1.0 / (1.0 + l)

        if _worker_raw_average:
            # Original measure_thickness.py behavior: weighted-average raw profiles,
            # then invert and normalize the composite.
            dat = np.average(_worker_thickness_arr[neighbors], weights=weights, axis=0) * -1
            dat = dat - dat.min()
            dat_sum = dat.sum()
            if dat_sum == 0:
                results.append((np.nan, 0))
                continue
            dat = dat / (80 / 81 * dat_sum)
        else:
            # Invert and normalize each neighbor's profile individually before averaging.
            indiv = _worker_thickness_arr[neighbors] * -1
            mins = indiv.min(axis=1, keepdims=True)
            indiv = indiv - mins
            row_sums = indiv.sum(axis=1, keepdims=True)
            valid_rows = row_sums.flatten() > 0
            if not np.any(valid_rows):
                results.append((np.nan, 0))
                continue
            indiv[valid_rows] /= (80 / 81 * row_sums[valid_rows])

            valid_profs = indiv[valid_rows]
            valid_weights = weights[valid_rows]

            dat0 = np.average(valid_profs, weights=valid_weights, axis=0)
            if _worker_rolling_xcorr:
                # xcorr-align each neighbor to the initial average before re-averaging,
                # sharpening bilayer peaks and improving dual Gaussian success rate.
                shifts_nm = _xcorr_offsets_batch(valid_profs, dat0, _worker_x)
                shift_samples = np.clip(
                    np.round(shifts_nm / sample_spacing).astype(int),
                    -max_shift_samples, max_shift_samples
                )
                n_cols = valid_profs.shape[1]
                col_idx = (np.arange(n_cols)[np.newaxis, :] + shift_samples[:, np.newaxis]) % n_cols
                aligned_profs = valid_profs[np.arange(len(shift_samples))[:, np.newaxis], col_idx]
                dat = np.average(aligned_profs, weights=valid_weights, axis=0)
            else:
                dat = dat0

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
    # Fit-quality thresholds are module-level constants (R2_THRESHOLD,
    # MIN_THICKNESS, MAX_THICKNESS).

    # Sample spacing — constant across all triangles, computed once per chunk
    sample_spacing = _worker_x[1] - _worker_x[0] if len(_worker_x) > 1 else 1.0
    max_shift_samples = max(1, int(round(2.0 / sample_spacing)))

    results = []
    for i in indices:
        # Get valid neighbors (those within distance bound)
        valid_mask = _worker_distances[i] != np.inf
        l = _worker_distances[i][valid_mask]
        neighbors = _worker_neighbor_indices[i][valid_mask]

        if len(neighbors) == 0:
            results.append((0, np.nan, np.nan, 0))
            continue

        weights = 1.0 / (1.0 + l)

        # Invert and normalize each neighbor's profile individually
        indiv = _worker_thickness_arr[neighbors] * -1
        mins = indiv.min(axis=1, keepdims=True)
        indiv = indiv - mins
        row_sums = indiv.sum(axis=1, keepdims=True)
        valid_rows = row_sums.flatten() > 0
        if not np.any(valid_rows):
            results.append((0, np.nan, np.nan, 0))
            continue
        indiv[valid_rows] /= (80 / 81 * row_sums[valid_rows])

        valid_profs = indiv[valid_rows]
        valid_weights = weights[valid_rows]

        if _worker_use_xcorr:
            # xcorr mode: simple weighted average, then match against global reference
            dat = np.average(valid_profs, weights=valid_weights, axis=0)
        else:
            # Gaussian mode: xcorr-align each neighbor profile to the initial average
            # before re-averaging.  This removes per-triangle positional variation so the
            # re-averaged profile has sharper bilayer peaks — improving dual Gaussian
            # success rate in low-SNR regions.
            dat0 = np.average(valid_profs, weights=valid_weights, axis=0)
            # Vectorized batch xcorr + roll (no Python loop over neighbors)
            shifts_nm = _xcorr_offsets_batch(valid_profs, dat0, _worker_x)
            shift_samples = np.clip(
                np.round(shifts_nm / sample_spacing).astype(int),
                -max_shift_samples, max_shift_samples
            )
            n_cols = valid_profs.shape[1]
            col_idx = (np.arange(n_cols)[np.newaxis, :] + shift_samples[:, np.newaxis]) % n_cols
            aligned_profs = valid_profs[np.arange(len(shift_samples))[:, np.newaxis], col_idx]
            dat = np.average(aligned_profs, weights=valid_weights, axis=0)

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

        # Step 1: Try dual Gaussian (unless monolayer mode). Parameterized by the
        # bilayer center + half-separation so the two leaflet peaks are forced to
        # opposite sides and cannot collapse into the same leaflet.
        if not _worker_monolayer:
            try:
                if _worker_global_fit_params is not None:
                    # Seed from the (more reliable) global average fit.
                    c1_g, w1_g, c2_g, w2_g = _worker_global_fit_params
                    center_seed = 0.5 * (c1_g + c2_g)
                    half_seed = float(np.clip(0.5 * (c2_g - c1_g),
                                              MIN_THICKNESS / 2.0, MAX_THICKNESS / 2.0))
                    w1_seed, w2_seed = w1_g, w2_g
                else:
                    center_seed, half_seed = _seed_bilayer_center(a, b)
                    w1_seed = w2_seed = 1.5
                # Restrict the fit to a window around the seeded bilayer so distant
                # density (adjacent membranes, CTF ringing) cannot pull the fit out.
                fit_window = MAX_THICKNESS / 2.0 + 2.0
                wmask = (a >= center_seed - fit_window) & (a <= center_seed + fit_window)
                af, bf = (a[wmask], b[wmask]) if int(wmask.sum()) >= 5 else (a, b)

                # params: h1, w1, h2, w2, center, half_sep, offset
                p0 = [0.02, w1_seed, 0.02, w2_seed, center_seed, half_seed, 0.0]
                bounds = ([0.005, 0.8, 0.005, 0.8, center_seed - 3.0, MIN_THICKNESS / 2.0, -1],
                          [0.04,  2.2, 0.04,  2.2, center_seed + 3.0, MAX_THICKNESS / 2.0,  1])
                p3, _ = opt.curve_fit(_dual_gaussian_centered, af, bf, p0, bounds=bounds)

                # Validate fit
                y_pred = _dual_gaussian_centered(af, *p3)
                r2 = _compute_r_squared(bf, y_pred)
                thickness = 2.0 * p3[5]

                if r2 > R2_THRESHOLD and MIN_THICKNESS <= thickness <= MAX_THICKNESS:
                    best_offset = p3[4]          # bilayer center
                    best_sigma1 = p3[1]
                    best_sigma2 = p3[3]
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
