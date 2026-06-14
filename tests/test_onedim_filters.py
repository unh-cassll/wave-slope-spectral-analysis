import numpy as np
import pytest

from slopespectra import spectf, notch_filter, filter_dispersive_bandstop


def test_spectf_recovers_sine_variance():
    fs = 10.0
    t = np.arange(4096) / fs
    a, f0 = 1.5, 1.0
    x = a * np.sin(2 * np.pi * f0 * t)
    S = spectf(x, 1 / fs, nfa=1)
    f, Sxx = S[:, 0], S[:, 1]
    var_spec = np.trapezoid(Sxx, f)
    assert np.isclose(var_spec, a**2 / 2, rtol=0.05)
    assert np.isclose(f[np.argmax(Sxx)], f0, atol=2 * np.diff(f).mean())


def test_spectf_band_averaging_shapes():
    rng = np.random.default_rng(3)
    x = rng.standard_normal(2048)
    S1 = spectf(x, 0.1, nfa=1)
    S31 = spectf(x, 0.1, nfa=31)
    assert S31.shape[0] < S1.shape[0]
    # Band averaging conserves total variance
    v1 = np.trapezoid(S1[:, 1], S1[:, 0])
    v31 = np.trapezoid(S31[:, 1], S31[:, 0])
    assert np.isclose(v1, v31, rtol=0.05)


def test_spectf_cross_spectrum_phase():
    fs = 10.0
    t = np.arange(4096) / fs
    f0 = 0.5
    x = np.sin(2 * np.pi * f0 * t)
    y = np.sin(2 * np.pi * f0 * (t - 0.25))   # y lags x by 45 deg
    S = spectf(x, 1 / fs, nfa=1, y=y)
    f = np.real(S[:, 0])
    i0 = np.argmin(np.abs(f - f0))
    assert np.isclose(np.real(S[i0, 4]), 45.0, atol=3.0)
    assert np.real(S[i0, 5]) > 0.99


def test_notch_filter_removes_band_preserves_mean():
    fs = 5.0
    t = np.arange(2000) / fs
    x = 3.0 + np.sin(2 * np.pi * 0.5 * t) + 0.5 * np.sin(2 * np.pi * 0.05 * t)
    out, _ = notch_filter(x, 1 / fs, low_period_s=1.0, high_period_s=4.0)
    # The 2 s period component is in the stopped band
    S = spectf(out - out.mean(), 1 / fs, nfa=1)
    f, Sxx = S[:, 0], S[:, 1]
    in_band = (f > 0.45) & (f < 0.55)
    out_band = (f > 0.04) & (f < 0.06)
    assert Sxx[in_band].max() < 1e-3 * Sxx[out_band].max()
    assert np.isclose(out.mean(), x.mean(), atol=1e-6)


def test_notch_filter_odd_length():
    x = np.sin(np.arange(501) * 0.1)
    out, _ = notch_filter(x, 0.5, 1.0, 10.0)
    assert out.shape == x.shape


def test_dispersive_bandstop_removes_shell_energy():
    # Deep-water wave on the dispersion shell
    g = 9.81
    m, nt, dx, fs = 32, 64, 0.2, 4.0
    k0 = 2 * np.pi / (m * dx) * 4
    f0 = np.sqrt(g * k0) / (2 * np.pi)
    x = np.arange(m) * dx
    t = np.arange(nt) / fs
    X, _, T = np.meshgrid(x, x, t, indexing="xy")
    stack = np.cos(k0 * X - 2 * np.pi * f0 * T)
    rng = np.random.default_rng(5)
    noise = 0.1 * rng.standard_normal(stack.shape)
    filtered = filter_dispersive_bandstop(stack + noise, dx, 1 / fs,
                                          halfwidth_hz=0.2)
    assert filtered.shape == stack.shape
    # Shell energy is strongly suppressed
    assert filtered.var() < 0.5 * (stack + noise).var()
