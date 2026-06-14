import numpy as np

from slopespectra import polar_from_cartesian, wavenumber_grids


def _gaussian_blob(m, dx, kx0, ky0, width):
    kx, ky = wavenumber_grids(m, dx)
    S = np.exp(-((kx - kx0)**2 + (ky - ky0)**2) / (2 * width**2))
    kmax = np.pi / dx
    S[np.sqrt(kx**2 + ky**2) > kmax] = np.nan
    return S, kx, ky


def test_variance_conservation():
    m, dx = 64, 0.05
    dk = 2 * np.pi / (m * dx)
    S, kx, ky = _gaussian_blob(m, dx, 10 * dk, 0.0, 3 * dk)
    Sk, k, theta = polar_from_cartesian(S, kx, ky)
    var_cart = np.nansum(S) * dk * dk
    var_polar = np.nansum(k[:, None] * Sk) * dk * np.median(np.diff(theta))
    assert np.isclose(var_polar, var_cart, rtol=1e-6)


def test_directional_peak_east():
    m, dx = 64, 0.05
    dk = 2 * np.pi / (m * dx)
    S, kx, ky = _gaussian_blob(m, dx, 10 * dk, 0.0, 2 * dk)
    Sk, k, theta = polar_from_cartesian(S, kx, ky)
    D = np.nansum(k[:, None] * Sk, axis=0)
    az = np.rad2deg(theta[np.argmax(D)])
    assert abs((az - 90 + 180) % 360 - 180) <= 10


def test_heading_offset_rotates_labels():
    m, dx = 64, 0.05
    dk = 2 * np.pi / (m * dx)
    S, kx, ky = _gaussian_blob(m, dx, 10 * dk, 0.0, 2 * dk)
    Sk0, k, theta = polar_from_cartesian(S, kx, ky, heading_deg=0.0)
    Sk45, _, _ = polar_from_cartesian(S, kx, ky, heading_deg=45.0)
    D0 = np.nansum(k[:, None] * Sk0, axis=0)
    D45 = np.nansum(k[:, None] * Sk45, axis=0)
    az0 = np.rad2deg(theta[np.argmax(D0)])
    az45 = np.rad2deg(theta[np.argmax(D45)])
    assert abs(((az45 - az0) - 45 + 180) % 360 - 180) <= 10


def test_radial_peak_location():
    m, dx = 64, 0.05
    dk = 2 * np.pi / (m * dx)
    S, kx, ky = _gaussian_blob(m, dx, 12 * dk, 0.0, 2 * dk)
    Sk, k, theta = polar_from_cartesian(S, kx, ky)
    dtheta = np.median(np.diff(theta))
    Sk_omni = np.nansum(k[:, None] * Sk, axis=1) * dtheta
    assert np.isclose(k[np.nanargmax(Sk_omni)], 12 * dk, atol=2 * dk)


def test_near_nyquist_rings_populated_from_valid_cells():
    # Rings overlapping the masked corner are reconstructed from their
    # valid cells (azimuthal mean x analytic ring area), not dropped
    m, dx = 64, 0.05
    dk = 2 * np.pi / (m * dx)
    S, kx, ky = _gaussian_blob(m, dx, 25 * dk, 0.0, 2 * dk)
    Sk, k, theta = polar_from_cartesian(S, kx, ky)
    # The outermost annulus has valid cells, so it is finite, not NaN
    assert np.isfinite(Sk[-1, :]).any()


def test_radial_smoothness_isotropic():
    # Constant Cartesian density -> S(k) proportional to k, monotone and
    # smooth (no Gauss-circle sawtooth)
    m, dx = 128, 0.05
    dk = 2 * np.pi / (m * dx)
    kx, ky = wavenumber_grids(m, dx)
    S = np.ones_like(kx)
    S[np.sqrt(kx**2 + ky**2) > np.pi / dx] = np.nan
    Sk, k, theta = polar_from_cartesian(S, kx, ky)
    dtheta = np.median(np.diff(theta))
    S_k = np.nansum(k[:, None] * Sk, axis=1) * dtheta
    ratio = (S_k / (2 * np.pi * k))[2:m // 4]
    # Flat to within 1% (was a ~25% sawtooth before the area normalization)
    assert np.std(ratio) / np.mean(ratio) < 0.01
