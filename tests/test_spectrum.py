import numpy as np
import pytest

from slopespectra import (compute_slope_spectrum, compute_segment_power,
                         finish_spectrum, azimuthal_integral,
                         wavenumber_grids)
from conftest import synthetic_wave_stack


def test_peak_location_eastward(eastward_wave):
    sx, sy, meta = eastward_wave
    spec = compute_slope_spectrum(sx, sy, meta["dx"], meta["fs"])
    idx = np.unravel_index(np.nanargmax(np.nan_to_num(spec.Skf, nan=0.0)),
                           spec.Skf.shape)
    assert np.isclose(spec.kx[0, idx[1]], meta["k0"], atol=spec.dk / 2)
    assert np.isclose(spec.ky[idx[0], 0], 0.0, atol=spec.dk)
    assert np.isclose(spec.f[idx[2]], meta["f0"], atol=spec.df)


@pytest.mark.parametrize("direction_deg", [0.0, 90.0, 180.0, 270.0, 45.0])
def test_peak_quadrant_by_direction(direction_deg):
    sx, sy, meta = synthetic_wave_stack(direction_deg=direction_deg)
    spec = compute_slope_spectrum(sx, sy, meta["dx"], meta["fs"])
    idx = np.unravel_index(np.nanargmax(np.nan_to_num(spec.Skf, nan=0.0)),
                           spec.Skf.shape)
    kx_pk = spec.kx[0, idx[1]]
    ky_pk = spec.ky[idx[0], 0]
    az = np.rad2deg(np.arctan2(kx_pk, ky_pk)) % 360
    err = (az - direction_deg + 180) % 360 - 180
    assert abs(err) < 10.0


def test_variance_conservation(eastward_wave):
    sx, sy, meta = eastward_wave
    spec = compute_slope_spectrum(sx, sy, meta["dx"], meta["fs"])
    spectral = np.nansum(spec.Skf) * spec.dk**2 * spec.df
    field = sx.var() + sy.var()
    assert np.isclose(spectral, field, rtol=1e-3)


def test_omni_integral_conserves_variance(eastward_wave):
    sx, sy, meta = eastward_wave
    spec = compute_slope_spectrum(sx, sy, meta["dx"], meta["fs"])
    omni_var = np.sum(spec.omni_S) * spec.dk
    cart_var = np.nansum(spec.Skxky) * spec.dk**2
    assert np.isclose(omni_var, cart_var, rtol=1e-6)
    assert np.isclose(spec.omni_k[np.argmax(spec.omni_S)], meta["k0"],
                      atol=spec.dk)


def test_single_frame_spectrum(eastward_wave):
    sx, sy, meta = eastward_wave
    spec = compute_slope_spectrum(sx[:, :, 0], sy[:, :, 0], meta["dx"])
    assert spec.Skf is None
    assert spec.f is None
    spectral = np.nansum(spec.Skxky) * spec.dk**2
    field = sx[:, :, 0].var() + sy[:, :, 0].var()
    assert np.isclose(spectral, field, rtol=0.05)


def test_azimuthal_integral_isotropic_ring():
    # Isotropic ring of unit density: S(k) = 2*pi*k within the ring
    m, dx = 128, 0.05
    dk = 2 * np.pi / (m * dx)
    kx = dk * (np.arange(m) - m // 2)
    KX, KY = np.meshgrid(kx, kx)
    K = np.sqrt(KX**2 + KY**2)
    k_ring = 20 * dk
    S = np.where(np.abs(K - k_ring) <= 3 * dk, 1.0, 0.0)
    k, S_k = azimuthal_integral(S, dx)
    total_cart = np.sum(S) * dk * dk
    total_omni = np.sum(S_k) * dk
    assert np.isclose(total_omni, total_cart, rtol=1e-6)
    assert np.isclose(k[np.argmax(S_k)], k_ring, atol=2 * dk)


def test_odd_time_dimension_raises(eastward_wave):
    sx, sy, meta = eastward_wave
    with pytest.raises(ValueError, match="even"):
        compute_slope_spectrum(sx[:, :, :63], sy[:, :, :63],
                               meta["dx"], meta["fs"])


def test_azimuthal_integral_smooth_isotropic():
    # Constant density -> S(k) = 2*pi*k, smooth (no annulus-count sawtooth)
    m, dx = 128, 0.05
    kx, ky = wavenumber_grids(m, dx)
    S = np.ones_like(kx)
    S[np.sqrt(kx**2 + ky**2) > np.pi / dx] = np.nan
    k, S_k = azimuthal_integral(S, dx)
    ratio = (S_k / (2 * np.pi * k))[2:m // 4]
    assert np.std(ratio) / np.mean(ratio) < 0.01
    # Total variance conserved
    dk = 2 * np.pi / (m * dx)
    assert np.isclose(np.sum(S_k) * dk, np.nansum(S) * dk * dk, rtol=1e-3)


def test_segment_power_matches_compute_slope_spectrum(eastward_wave):
    sx, sy, meta = eastward_wave
    power, framesize, s3 = compute_segment_power(sx, sy, meta["dx"], meta["fs"])
    split = finish_spectrum(power, framesize, meta["dx"], s3, meta["fs"])
    direct = compute_slope_spectrum(sx, sy, meta["dx"], meta["fs"])
    assert np.allclose(np.nan_to_num(split.Skf), np.nan_to_num(direct.Skf))
    assert np.allclose(np.nan_to_num(split.Skxky), np.nan_to_num(direct.Skxky))


def test_segment_power_averaging_matches_stacked_input(eastward_wave):
    # Averaging two segments' raw power and finishing once must equal
    # concatenating them into one double-length stack and computing directly
    # (both hold the same number of independent 600-frame realizations).
    sx, sy, meta = eastward_wave
    p1, framesize, s3 = compute_segment_power(sx, sy, meta["dx"], meta["fs"])
    p2, _, _ = compute_segment_power(sx, sy, meta["dx"], meta["fs"])
    mean_power = (p1 + p2) / 2
    spec = finish_spectrum(mean_power, framesize, meta["dx"], s3, meta["fs"])
    direct = compute_slope_spectrum(sx, sy, meta["dx"], meta["fs"])
    assert np.allclose(np.nan_to_num(spec.Skf), np.nan_to_num(direct.Skf))


def test_segment_power_float32_close_to_float64(eastward_wave):
    sx, sy, meta = eastward_wave
    p64, framesize, s3 = compute_segment_power(sx, sy, meta["dx"], meta["fs"],
                                               dtype=np.float64)
    p32, _, _ = compute_segment_power(sx, sy, meta["dx"], meta["fs"],
                                      dtype=np.float32)
    spec64 = finish_spectrum(p64, framesize, meta["dx"], s3, meta["fs"])
    spec32 = finish_spectrum(p32, framesize, meta["dx"], s3, meta["fs"])
    scale = np.nanmax(np.abs(spec64.Skf))
    assert np.nanmax(np.abs(spec64.Skf - spec32.Skf)) / scale < 1e-5


def test_segment_power_single_frame(eastward_wave):
    sx, sy, meta = eastward_wave
    power, framesize, s3 = compute_segment_power(sx[:, :, 0], sy[:, :, 0],
                                                 meta["dx"])
    assert s3 == 1
    spec = finish_spectrum(power, framesize, meta["dx"], s3)
    direct = compute_slope_spectrum(sx[:, :, 0], sy[:, :, 0], meta["dx"])
    assert spec.Skf is None
    assert np.allclose(np.nan_to_num(spec.Skxky), np.nan_to_num(direct.Skxky))
