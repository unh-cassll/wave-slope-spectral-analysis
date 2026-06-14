"""Frequency-averaged 1D auto- and cross-spectra (geophysical normalization).

Port of K. K. Kahma's SPECTF (1990): Blackman-Harris window, linear detrend,
band-averaged one-sided spectra in unit^2 * s (i.e. unit^2 / Hz).
"""

import numpy as np
from scipy.signal import detrend

from .windows import blackman_harris

__all__ = ["spectf"]


def _band_average(P, nfa, a0, maxb, C, df):
    """Average elementary bands a0..maxb-1 in groups of nfa."""
    if nfa == 1:
        f = np.arange(a0, maxb) * df
        return f, P[a0:maxb] * C
    m = (maxb - a0) // nfa
    idx0 = a0
    fout = np.empty(0)
    Pout = np.empty(0)
    if m > 0:
        block = P[idx0:idx0 + m * nfa].reshape(m, nfa)
        Pout = block.sum(axis=1) * C
        fout = (np.arange(1, m + 1) * nfa + (a0 - 0.5 - nfa / 2)) * df
    a = a0 + m * nfa
    if a < maxb:
        c = maxb - a
        Pout = np.append(Pout, P[a:maxb].sum() * C * nfa / c)
        fout = np.append(fout, df * (a + maxb - 1) / 2)
    return fout, Pout


def spectf(x, dt_s, nfa=31, a0=0, y=None):
    """Frequency-averaged power (or cross) spectrum estimate.

    Args:
        x    : data vector
        dt_s : sampling interval [s]
        nfa  : number of elementary frequency bands averaged per output band
        a0   : number of low elementary bands to discard (0 keeps the mean)
        y    : optional second data vector for the cross spectrum

    Returns:
        Auto-spectrum: (nbands, 2) array with columns [f_Hz, Sxx].
        Cross-spectrum: (nbands, 6) array with columns
        [f_Hz, Sxx, Syy, Sxy, phase_deg, coherence]; positive phase means
        x leads y.
    """
    x = np.asarray(x, dtype=float).ravel()
    if x.size % 2:
        x = x[:-1]
    N = x.size
    window = blackman_harris(N)

    Xx = np.fft.fft(window * detrend(x, type="linear"))
    maxb = N // 2 + 1
    Xx = Xx[:maxb].copy()
    Xx[-1] /= 2

    C = dt_s / (nfa * np.pi * np.linalg.norm(window)**2)
    df = 2 * np.pi / (dt_s * N)

    if y is None:
        f, Pxx = _band_average(np.abs(Xx)**2, nfa, a0, maxb, C, df)
        return np.column_stack([f / (2 * np.pi), 2 * np.pi * Pxx])

    y = np.asarray(y, dtype=float).ravel()
    if y.size % 2:
        y = y[:-1]
    if y.size != N:
        raise ValueError("x and y are not of the same size")
    Yy = np.fft.fft(window * detrend(y, type="linear"))
    Yy = Yy[:maxb].copy()
    Yy[-1] /= 2

    f, Pxx = _band_average(np.abs(Xx)**2, nfa, a0, maxb, C, df)
    _, Pyy = _band_average(np.abs(Yy)**2, nfa, a0, maxb, C, df)
    _, Pxy = _band_average(np.conj(Xx) * Yy, nfa, a0, maxb, C, df)

    phase = np.rad2deg(np.arctan2(-np.imag(Pxy), np.real(Pxy)))
    coh = np.abs(Pxy / np.sqrt(Pxx * Pyy))
    return np.column_stack([f / (2 * np.pi), 2 * np.pi * Pxx, 2 * np.pi * Pyy,
                            2 * np.pi * Pxy, phase, coh])
