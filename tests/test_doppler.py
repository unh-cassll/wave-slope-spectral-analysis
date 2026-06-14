import numpy as np
import pytest

from slopespectra import (compute_slope_spectrum, moving_window_fourier,
                          upwind_downwind_wave_frequencies)
from conftest import synthetic_wave_stack


@pytest.fixture(scope="module")
def northward_wave():
    # Toward +y ("downwind" half of the spectrum), well inside the
    # 50 rad/m analysis cutoff
    return synthetic_wave_stack(m=64, nt=64, dx=0.2, fs=4.0, n_bin=8,
                                direction_deg=0.0)


def test_downwind_shell_frequency(northward_wave):
    sx, sy, meta = northward_wave
    spec = compute_slope_spectrum(sx, sy, meta["dx"], meta["fs"])
    ds = upwind_downwind_wave_frequencies(spec.Skf, meta["dx"], meta["fs"])
    i_k = np.argmin(np.abs(ds["k"].values - meta["k0"]))
    f_down = float(ds["f_down"][i_k])
    assert np.isfinite(f_down)
    assert np.isclose(f_down, meta["f0"], atol=2 * spec.df)
    # Downwind ridge direction near 0 deg
    assert abs(float(ds["dir_down"][i_k])) < 20


def test_moving_window_fourier_runs(northward_wave):
    sx, sy, meta = northward_wave
    ds = moving_window_fourier(sx, sy, meta["dx"], meta["fs"],
                               window_length_s=8.0, window_step_s=4.0)
    assert "t" in ds.dims
    assert ds.sizes["t"] >= 2
    # Each window resolves the same dominant downwind frequency
    i_k = np.argmin(np.abs(ds["k"].values - meta["k0"]))
    f_down = ds["f_down"][:, i_k].values
    assert np.all(np.isfinite(f_down))
    assert np.allclose(f_down, meta["f0"], atol=0.3)
    assert "S_omni" in ds
