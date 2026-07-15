"""Tapering windows for slope-field stacks and time series."""

import numpy as np
from scipy.signal.windows import tukey
from scipy.stats import t as _student_t

__all__ = ["circular_tukey", "blackman_harris", "taper_2d"]


def _radial_tukey(s1, s2, taper_width):
    """Unnormalized circular Tukey window on an (s1, s2) grid."""
    min_dim = min(s1, s2)
    y, x = np.mgrid[1:s1 + 1, 1:s2 + 1].astype(float)
    x -= x.mean()
    y -= y.mean()
    r = np.sqrt(x**2 + y**2) / (min_dim / 2)

    # Cosine taper between r = 1 - taper_width and r = 1
    cos_period = taper_width * 2
    rshift = 1 - taper_width
    w2d = (np.cos(2 * np.pi / cos_period * (r - rshift)) + 1) / 2
    w2d[r < (1 - taper_width)] = 1.0
    w2d[r > 1] = 0.0
    return w2d


def circular_tukey(in_array, taper_width=0.2, normalization="power",
                   temporal_alpha=0.0, dtype=None):
    """Circular Tukey (tapered cosine) window applied to a 2D or 3D array.

    A 3D (ny, nx, nt) array is tapered radially in (y, x), optionally also
    along t with a 1D Tukey window.

    normalization:
      "power"    : scale so mean(w^2) = 1; the windowed variance is an
                   unbiased estimate of the field variance (standard PSD
                   window-power correction)
      "spectral" : scale the spectral norm of the spatial window to that
                   of a uniform window (MATLAB-library parity; biases the
                   variance high by ~7% at taper_width = 0.2)
      "none"     : raw window, max 1

    Args:
        in_array       : (ny, nx) or (ny, nx, nt) array
        taper_width    : cosine taper fraction of the inscribed radius
        temporal_alpha : Tukey cosine fraction along t (0 disables)
        dtype          : output/working precision (default float64); the
                         window itself is also cast, so a real-valued dtype
                         such as float32 halves the memory footprint of the
                         windowed array for large stacks

    Returns:
        Windowed array, same shape as input, cast to dtype.
    """
    dtype = np.float64 if dtype is None else dtype
    in_array = np.asarray(in_array, dtype=dtype)
    if in_array.ndim == 2:
        s1, s2 = in_array.shape
        s3 = 1
    else:
        s1, s2, s3 = in_array.shape

    w2d = _radial_tukey(s1, s2, taper_width)
    w_t = tukey(s3, temporal_alpha) if (temporal_alpha > 0 and s3 > 1) \
        else np.ones(s3)

    if normalization == "power":
        C = 1.0 / np.sqrt(np.mean(w2d**2) * np.mean(w_t**2))
    elif normalization == "spectral":
        C = np.sqrt(s1 * s2) / np.linalg.norm(w2d, 2)
    elif normalization == "none":
        C = 1.0
    else:
        raise ValueError("normalization must be 'power', 'spectral' or "
                         "'none'")

    w2d = (C * w2d).astype(dtype)
    w_t = w_t.astype(dtype)
    if in_array.ndim == 3:
        return w2d[:, :, None] * w_t[None, None, :] * in_array
    return w2d * in_array


def blackman_harris(n):
    """N-point 4-term Blackman-Harris window (column-vector analogue)."""
    m = np.arange(n) * (2 * np.pi / (n - 1))
    return (0.35875 - 0.48829 * np.cos(m) + 0.14128 * np.cos(2 * m)
            - 0.01168 * np.cos(3 * m))


def taper_2d(X, kind="cosine", w=0.2, s=3):
    """Border taper of a 2D array (or each plane of a 3D array).

    kind='cosine'  : separable Tukey window, w = cosine fraction
    kind='t'       : student's-t CDF border ramp, w = integer half-height
                     position, s = degrees of freedom

    Returns:
        (tapered_array, window)
    """
    X = np.asarray(X, dtype=float)
    if X.ndim == 2:
        n, m = X.shape
        p = 1
    else:
        n, m, p = X.shape

    if kind == "cosine":
        window = (np.outer(tukey(n, w), np.ones(m))
                  * np.rot90(np.outer(tukey(m, w), np.ones(n))))
    elif kind == "t":
        x = np.arange(-w, w + 1, dtype=float)
        ramp = _student_t.cdf(x, s)
        ramp[0] = 0.0
        rng = len(x)
        vert = np.ones(n)
        vert[:rng] = ramp
        vert[n - rng:] = ramp[::-1]
        horz = np.ones(m)
        horz[:rng] = ramp
        horz[m - rng:] = ramp[::-1]
        window = np.outer(vert, np.ones(m)) * np.rot90(np.outer(horz, np.ones(n)))
    else:
        raise ValueError("kind must be 'cosine' or 't'")

    if X.ndim == 3:
        return X * window[:, :, None], window
    return X * window, window
