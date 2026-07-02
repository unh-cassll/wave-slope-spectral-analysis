"""Fourier-domain filters for image stacks and time series."""

import numpy as np
from scipy.signal import detrend

__all__ = ["filter_dispersive_bandstop", "notch_filter"]

GRAV = 9.81  # m/s^2


def filter_dispersive_bandstop(stack, dx_m, dt_s, ux=0.0, uy=0.0,
                               halfwidth_hz=0.1):
    """Remove energy along the (Doppler-shifted) gravity wave dispersion
    shell from an image stack; the retained coefficients are rescaled so a
    spectrally flat field keeps its total energy.

    Both conjugate branches of the shell are stopped, so the filtered
    stack is real-valued by construction.

    Args:
        stack        : (ny, nx, nt) image stack, rows with the +y edge first
        dx_m         : pixel size [m]
        dt_s         : frame interval [s]
        ux, uy       : advection components Doppler-shifting the shell [m/s]
        halfwidth_hz : half-width of the stopped band about the shell [Hz]

    Returns:
        Filtered stack, same shape.
    """
    stack = np.asarray(stack, dtype=float)
    ny, nx, nt = stack.shape

    # Label arrays matching the fftshifted FFT layout (ky descending)
    kx = 2 * np.pi * np.fft.fftshift(np.fft.fftfreq(nx, dx_m))
    ky = -2 * np.pi * np.fft.fftshift(np.fft.fftfreq(ny, dx_m))
    f = np.fft.fftshift(np.fft.fftfreq(nt, dt_s))
    df = 1.0 / (nt * dt_s)

    Kx, Ky, F = np.meshgrid(kx, ky, f)
    K = np.sqrt(Kx**2 + Ky**2)

    Mf = np.fft.fftshift(np.fft.fftn(stack, s=(ny, nx, nt), axes=(0, 1, 2)))

    # Stop both branches of the Doppler-shifted dispersion shell
    f_intr = np.sqrt(GRAV * K) / (2 * np.pi)
    f_dopp = (ux * Kx + uy * Ky) / (2 * np.pi)
    band = halfwidth_hz + df / 2
    mask = np.ones_like(stack)
    mask[np.abs(F - (f_intr - f_dopp)) < band] = 0.0
    mask[np.abs(F + (f_intr + f_dopp)) < band] = 0.0

    # Rescale retained coefficients: preserves a flat spectrum's energy
    coeff = np.sqrt(stack.size / mask.sum())

    return np.real(np.fft.ifftn(np.fft.ifftshift(coeff * mask * Mf)))


def notch_filter(signal_in, dt_s, low_period_s, high_period_s):
    """Remove energy between periods low_period_s and high_period_s from a
    time series; mean and trend are preserved.

    Args:
        signal_in     : 1D input signal
        dt_s          : sampling interval [s]
        low_period_s  : low period cutoff [s]
        high_period_s : high period cutoff [s]

    Returns:
        (real_part, imag_part) of the filtered signal.
    """
    si = np.asarray(signal_in, dtype=float).ravel()
    lsi = len(si)
    odd = lsi % 2 == 1
    if odd:
        si = np.concatenate([si, si[-1:]])
    n = len(si)

    trend = si - detrend(si, type="linear")
    cls_ = int(np.ceil(lsi / 2))
    f = np.arange(-cls_, cls_) / lsi / dt_s
    f = f[:n] if len(f) > n else np.pad(f, (0, n - len(f)), mode="edge")

    SI = np.fft.fftshift(np.fft.fft(detrend(si, type="linear")))
    SI[cls_] = 0.0

    band = (((f >= -1 / low_period_s) & (f <= -1 / high_period_s))
            | ((f <= 1 / low_period_s) & (f >= 1 / high_period_s)))
    SI[band] = 0.0

    out = np.fft.ifft(np.fft.ifftshift(SI))
    sireal = np.real(out) + trend
    siimag = np.imag(out)
    if odd:
        sireal = sireal[:lsi]
        siimag = siimag[:lsi]
    return sireal, siimag
