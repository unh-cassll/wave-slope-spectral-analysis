"""Directional wavenumber-frequency spectra from slope-field stacks."""

import os
from dataclasses import dataclass, field

import numpy as np
from scipy import fft as _scipy_fft

__all__ = ["compute_slope_spectrum", "compute_segment_power", "finish_spectrum",
           "azimuthal_integral", "wavenumber_grids", "DirectionalSpectrum"]

# Threads for the 3-D FFTs (scipy pocketfft; numerically identical to numpy).
# Default 1 preserves single-threaded behavior; set SLOPESPECTRA_FFT_WORKERS
# to the desired thread count for large stacks.
FFT_WORKERS = int(os.environ.get("SLOPESPECTRA_FFT_WORKERS", "1"))


def _fftn(a, s, axes):
    return _scipy_fft.fftn(a, s=s, axes=axes, workers=FFT_WORKERS)


def wavenumber_grids(framesize, dx_m):
    """Wavenumber label arrays for an FFT of size (framesize, framesize).

    Returns (kx, ky): kx ascending along columns, ky descending along rows
    (first row = +ky), matching the image-row convention in which row 0 is
    the +y edge of the field.
    """
    kmin = 2 * np.pi / (framesize * dx_m)
    m = framesize
    Kx = kmin * np.arange(-np.ceil((m - 1) / 2), np.floor((m - 1) / 2) + 1)
    Ky = kmin * np.arange(np.ceil((m - 1) / 2), -np.floor((m - 1) / 2) - 1, -1)
    kx = np.tile(Kx, (m, 1))
    ky = np.tile(Ky[:, None], (1, m))
    return kx, ky


def azimuthal_integral(Skxky, dx_m):
    """Omnidirectional wavenumber spectrum S(k) by annular averaging.

    Each ring's density is the azimuthal MEAN of the Cartesian cells it
    contains, multiplied by the analytic ring area 2*pi*k*dk:

        S(k_n) = mean_{annulus} S(kx, ky) * 2*pi*k_n .

    Using the analytic ring area rather than the discrete cell count
    removes the bin-to-bin sawtooth that arises because the number of grid
    cells per annulus fluctuates (the Gauss circle problem). The result is
    rescaled so total variance is conserved exactly:
    sum(S(k)) dk = sum(S) dk^2.

    Args:
        Skxky : (m, m) Cartesian spectral density, NaN beyond Nyquist
        dx_m  : grid spacing [m]

    Returns:
        (k, S_k) 1D arrays.
    """
    Skxky = np.asarray(Skxky, dtype=float)
    m = Skxky.shape[0]
    kx, ky = wavenumber_grids(m, dx_m)
    dk = 2 * np.pi / (m * dx_m)
    kmat = np.sqrt(kx**2 + ky**2)
    num_bins = m // 2
    k = dk * np.arange(1, num_bins + 1)

    flat = Skxky.ravel()
    bin_idx = np.round(kmat / dk).astype(int).ravel()
    finite = np.isfinite(flat) & (bin_idx >= 1) & (bin_idx <= num_bins)
    sums = np.bincount(bin_idx[finite], weights=flat[finite],
                       minlength=num_bins + 1)[1:num_bins + 1]
    counts = np.bincount(bin_idx[finite],
                         minlength=num_bins + 1)[1:num_bins + 1]
    mean_density = np.divide(sums, counts, out=np.zeros_like(sums),
                             where=counts > 0)
    S_k = mean_density * 2 * np.pi * k

    # Conserve total variance: sum(S_k) dk == sum(S) dk^2
    total = np.nansum(flat[finite]) * dk * dk
    norm = S_k.sum() * dk
    if norm > 0:
        S_k *= total / norm
    return k, S_k


@dataclass
class DirectionalSpectrum:
    """Directional spectrum products of a slope-field stack.

    Skf is the directional wavenumber-frequency slope spectrum with axes
    (ky descending, kx ascending, f ascending positive); Skxky is its
    frequency integral. Units: slope^2 / (rad/m)^2 / Hz and
    slope^2 / (rad/m)^2.
    """
    kx: np.ndarray
    ky: np.ndarray
    dk: float
    Skxky: np.ndarray
    omni_k: np.ndarray
    omni_S: np.ndarray
    f: np.ndarray | None = None
    df: float | None = None
    Skf: np.ndarray | None = None
    attrs: dict = field(default_factory=dict)


def compute_segment_power(sx, sy, dx_m, fs_hz=None, framesize=None,
                          dtype=np.float64):
    """Unshifted directional wavenumber(-frequency) power of one slope-field
    segment, before the Nyquist-circle mask, physical normalization, and
    spatial fftshift.

    To average many segments (e.g. successive records from a long time
    series): sum or average the returned `power` across segments, then call
    `finish_spectrum` once on the result. The mask, normalization, and
    spatial shift are independent of the data, so applying them after
    averaging is exact and cheaper than once per segment.

    Args:
        sx, sy    : (ny, nx) or (ny, nx, nt) slope-field stacks; the time
                    dimension, if present, must be of even length
        dx_m      : pixel size [m]
        fs_hz     : frame rate [Hz]; required when nt > 1
        framesize : FFT size in the spatial dimensions (fields are
                    zero-padded up to this size); defaults to ny
        dtype     : real dtype for the FFT input (default float64); a
                    32-bit dtype roughly halves memory bandwidth for the
                    dominant FFT stage at ~1e-7 relative precision cost

    Returns:
        (power, framesize, s3): power has shape (framesize, framesize,
        s3//2) holding ascending positive frequencies f = df*(1..s3//2), in
        the unshifted spatial-frequency layout (index 0 = DC on axes 0-1);
        for a single frame, shape (framesize, framesize). s3 is the frame
        count, needed by `finish_spectrum` to recover df.
    """
    sx = np.asarray(sx, dtype=dtype)
    sy = np.asarray(sy, dtype=dtype)
    if sx.shape != sy.shape:
        raise ValueError("sx and sy must have the same shape")
    if sx.ndim == 2:
        sx = sx[:, :, None]
        sy = sy[:, :, None]
    s1, s2, s3 = sx.shape
    if framesize is None:
        framesize = s1
    if s1 > framesize or s2 > framesize:
        raise ValueError("spatial dimensions exceed framesize")
    if s3 > 1:
        if fs_hz is None:
            raise ValueError("fs_hz is required for a multi-frame stack")
        if s3 % 2:
            raise ValueError("the time dimension must be of even length")

    Ax = _fftn(sx, s=(framesize, framesize, s3), axes=(0, 1, 2))
    Ay = _fftn(sy, s=(framesize, framesize, s3), axes=(0, 1, 2))

    if s3 > 1:
        # Natural FFT order stores negative frequencies in the back half of
        # axis 2; selecting and reversing it gives ascending positive-f
        # labels without ever shifting the temporal axis.
        Ax = Ax[:, :, :s3 // 2 - 1:-1]
        Ay = Ay[:, :, :s3 // 2 - 1:-1]
        power = (Ax * np.conj(Ax) + Ay * np.conj(Ay)).real
    else:
        power = (Ax[:, :, 0] * np.conj(Ax[:, :, 0])
                + Ay[:, :, 0] * np.conj(Ay[:, :, 0])).real
    return power, framesize, s3


def finish_spectrum(power, framesize, dx_m, s3, fs_hz=None):
    """Shift, mask, and normalize a (possibly segment-averaged) power cube
    from `compute_segment_power` into a DirectionalSpectrum.

    Args:
        power     : power cube from `compute_segment_power` (or the mean of
                    several such cubes sharing the same framesize/s3)
        framesize, dx_m, s3, fs_hz : as passed to `compute_segment_power`

    Returns:
        DirectionalSpectrum, as `compute_slope_spectrum`.
    """
    kmin = 2 * np.pi / (framesize * dx_m)
    kmax = np.pi / dx_m
    kx, ky = wavenumber_grids(framesize, dx_m)

    k = np.sqrt(kx**2 + ky**2)
    mask = np.ones_like(k)
    mask[k > kmax] = np.nan

    dk = kmin
    df = fs_hz / s3 if fs_hz is not None else None
    power = np.fft.fftshift(np.asarray(power, dtype=float), axes=(0, 1))

    if s3 > 1:
        f = df * np.arange(1, s3 // 2 + 1)
        N = framesize * framesize * s3
        C = 2.0 / (N**2 * dk * dk * df)
        Skf = mask[:, :, None] * power * C
        # Temporal Nyquist plane has no conjugate partner: undo the
        # one-sided doubling there
        Skf[:, :, -1] *= 0.5
        Skxky = np.nansum(Skf, axis=2) * df
        Skxky[np.isnan(mask)] = np.nan
    else:
        f = None
        Skf = None
        N = framesize * framesize
        C = 1.0 / (N**2 * dk * dk)
        Skxky = mask * power * C

    omni_k, omni_S = azimuthal_integral(Skxky, dx_m)

    return DirectionalSpectrum(
        kx=kx, ky=ky, dk=dk, Skxky=Skxky, omni_k=omni_k, omni_S=omni_S,
        f=f, df=df, Skf=Skf,
        attrs={"dx_m": dx_m, "fs_hz": fs_hz, "framesize": framesize},
    )


def compute_slope_spectrum(sx, sy, dx_m, fs_hz=None, framesize=None,
                           dtype=np.float64):
    """Directional wavenumber(-frequency) spectrum of slope-field stacks.

    Args:
        sx, sy    : (ny, nx) or (ny, nx, nt) slope-field stacks; the time
                    dimension, if present, must be of even length
        dx_m      : pixel size [m]
        fs_hz     : frame rate [Hz]; required when nt > 1
        framesize : FFT size in the spatial dimensions (fields are
                    zero-padded up to this size); defaults to ny
        dtype     : real dtype for the FFT input (default float64)

    Returns:
        DirectionalSpectrum. For nt > 1, Skf holds the spectral density on
        positive frequencies f = df * (1..nt/2) and Skxky its integral
        over f; for a single frame only Skxky is computed. The f = 0 plane
        is excluded, so static (per-pixel temporal mean) variance is
        absent; remove the temporal mean beforehand for exact variance
        accounting.
    """
    power, framesize, s3 = compute_segment_power(sx, sy, dx_m, fs_hz,
                                                 framesize, dtype)
    return finish_spectrum(power, framesize, dx_m, s3, fs_hz)
