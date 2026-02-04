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


# Global variables for multiprocessing workers (set by initializer)
_worker_thickness_arr = None
_worker_distances = None
_worker_neighbor_indices = None
_worker_x = None
_worker_monolayer = False


def init_worker(thickness_arr, distances, neighbor_indices, x, monolayer=False):
    """Initialize worker process with shared data."""
    global _worker_thickness_arr, _worker_distances, _worker_neighbor_indices, _worker_x, _worker_monolayer
    _worker_thickness_arr = thickness_arr
    _worker_distances = distances
    _worker_neighbor_indices = neighbor_indices
    _worker_x = x
    _worker_monolayer = monolayer


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

    Used by refine_mesh.py for mesh refinement. Supports both bilayer (dual
    Gaussian) and monolayer (single Gaussian) fitting modes via _worker_monolayer.

    Returns
    -------
    list of tuples
        [(offset, sigma1, sigma2), ...] for each triangle in chunk.
        For monolayer mode, sigma2 will be np.nan.
    """
    results = []
    for i in indices:
        # Get valid neighbors (those within distance bound)
        valid_mask = _worker_distances[i] != np.inf
        l = _worker_distances[i][valid_mask]
        neighbors = _worker_neighbor_indices[i][valid_mask]

        if len(neighbors) == 0:
            results.append((0, np.nan, np.nan))
            continue

        # Vectorized weight calculation
        weights = 1.0 / (1.0 + l)

        # Compute weighted average of density profiles
        dat = np.average(_worker_thickness_arr[neighbors], weights=weights, axis=0) * -1
        dat = dat - dat.min()
        dat_sum = dat.sum()
        if dat_sum == 0:
            results.append((0, np.nan, np.nan))
            continue
        dat = dat / (80/81 * dat_sum)

        # Find initial peak and set up fitting bounds
        ipk = _worker_x[np.argmax(dat)]

        # Find minima to constrain fitting region
        mid = len(dat) // 2
        left_min = np.argmin(dat[:mid])
        right_min = np.argmin(dat[mid:]) + mid

        a = _worker_x[left_min+2:right_min-2]
        b = dat[left_min+2:right_min-2]

        if len(a) < 5:  # Need enough points for fitting
            results.append((0, np.nan, np.nan))
            continue

        try:
            if _worker_monolayer:
                # Single Gaussian fit for monolayer/high defocus
                p0_mono = [0.03, ipk, 2.5, 0]
                bounds_mono = ([0.005, ipk-6, 1.0, -1],
                               [0.06, ipk+6, 5.0, 1])
                p_mono, _ = opt.curve_fit(_monogaussian_with_offset, a, b, p0_mono, bounds=bounds_mono)
                offset = p_mono[1]
                sigma1 = p_mono[2]
                sigma2 = np.nan
            else:
                # Dual Gaussian fit for bilayer
                p0 = [0.02, ipk-0.5, 1.5, 0.02, ipk+0.5, 1.5, 0]
                bounds = ([0.005, ipk-6, 0.8, 0.005, ipk-1.5, 0.8, -1],
                          [0.04, ipk+1.5, 2.2, 0.04, ipk+6, 2.2, 1])
                p3, _ = opt.curve_fit(_dual_gaussian, a, b, p0, bounds=bounds)
                offset = (p3[1] + p3[4]) / 2
                sigma1 = p3[2]
                sigma2 = p3[5]
            results.append((offset, sigma1, sigma2))
        except Exception:
            results.append((0, np.nan, np.nan))

    return results
