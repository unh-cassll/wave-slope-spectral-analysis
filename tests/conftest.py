"""Shared fixtures: synthetic propagating-wave slope fields."""

import numpy as np
import pytest

GRAV = 9.81


def synthetic_wave_stack(m=64, nt=64, dx=0.05, fs=8.0, n_bin=8,
                         direction_deg=90.0, amplitude=0.05, phase=0.3):
    """Slope-field stack of a monochromatic deep-water wave.

    direction_deg is the compass going-to direction (0 = +y, 90 = +x).
    Returns (sx, sy, meta) with meta = dict(k0, f0, dx, fs).
    """
    k0 = 2 * np.pi / (m * dx) * n_bin
    f0 = np.sqrt(GRAV * k0) / (2 * np.pi)
    az = np.deg2rad(direction_deg)
    kx0 = k0 * np.sin(az)
    ky0 = k0 * np.cos(az)

    x = np.arange(m) * dx
    y = np.arange(m) * dx
    t = np.arange(nt) / fs
    # rows = y ascending toward +y; compass +y is the -row direction in the
    # image convention, hence the sign flip on ky0
    X, Y, T = np.meshgrid(x, y, t, indexing="xy")
    arg = kx0 * X - ky0 * Y - 2 * np.pi * f0 * T + phase
    a = amplitude
    sx = -a * kx0 * np.sin(arg)
    sy = a * ky0 * np.sin(arg)
    return sx, sy, {"k0": k0, "f0": f0, "dx": dx, "fs": fs,
                    "kx0": kx0, "ky0": ky0}


@pytest.fixture
def eastward_wave():
    return synthetic_wave_stack(direction_deg=90.0)


# Off-bin wavenumbers so spectral leakage spreads each component across
# several k bins (the dispersion-shell median filter rejects single-bin
# spikes by design)
BROADBAND_COMPONENTS = [(5.3, 80.0, 0.04), (8.4, 90.0, 0.05),
                        (11.6, 100.0, 0.03), (14.2, 95.0, 0.02)]


def broadband_wave_stack():
    """Superposition of dispersive waves spread about east."""
    rng = np.random.default_rng(7)
    sx = sy = None
    for n_bin, az, amp in BROADBAND_COMPONENTS:
        sx_i, sy_i, meta = synthetic_wave_stack(
            n_bin=n_bin, direction_deg=az, amplitude=amp,
            phase=rng.uniform(0, 2 * np.pi))
        sx = sx_i if sx is None else sx + sx_i
        sy = sy_i if sy is None else sy + sy_i
    return sx, sy, meta


@pytest.fixture
def broadband_wave():
    return broadband_wave_stack()
