import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest
import xarray as xr

from slopespectra import (compute_wave_spectra, plot_kf_slice,
                          plot_omnidirectional, plot_sub_spectra)
from conftest import synthetic_wave_stack


@pytest.fixture(scope="module")
def pipeline_result():
    sx, sy, meta = synthetic_wave_stack(m=70, nt=64)   # non-power-of-two
    ds, dirspect = compute_wave_spectra(sx, sy, dx_m=meta["dx"],
                                        fs_hz=meta["fs"],
                                        return_dirspect=True)
    return ds, dirspect, meta


def test_pipeline_crops_and_runs(pipeline_result):
    ds, dirspect, meta = pipeline_result
    assert ds.attrs["framesize"] == 64
    assert ds.attrs["mss"] > 0
    # Windowed spectrum still localizes the input wave
    k_pk = float(ds["k"][int(np.argmax(ds["S_k"].values))])
    assert np.isclose(k_pk, meta["k0"], atol=2 * ds.attrs["dk"])


def test_pipeline_accepts_dataset(pipeline_result):
    sx, sy, meta = synthetic_wave_stack(m=64, nt=32)
    src = xr.Dataset(
        {"slope_x": (("y", "x", "time"), sx),
         "slope_y": (("y", "x", "time"), sy),
         "fs": meta["fs"]},
        attrs={"dx_m": meta["dx"]},
    )
    ds = compute_wave_spectra(src)
    assert ds.attrs["fs_hz"] == meta["fs"]


def test_plotting_functions(pipeline_result):
    ds, dirspect, _ = pipeline_result
    fig = plot_sub_spectra(ds)
    assert len(fig.axes) >= 6
    ax = plot_kf_slice(dirspect, axis="kx")
    assert ax.get_xlabel().startswith("$k_x")
    ax = plot_kf_slice(dirspect, axis="ky")
    assert ax.get_xlabel().startswith("$k_y")
    ax = plot_omnidirectional(ds)
    assert ax.get_ylabel().startswith("$F(f)$")
    matplotlib.pyplot.close("all")


def test_spectral_total_unbiased():
    from conftest import broadband_wave_stack
    sx, sy, meta = broadband_wave_stack()
    ds = compute_wave_spectra(sx, sy, dx_m=meta["dx"], fs_hz=meta["fs"],
                              temporal_alpha=0.1)
    # Power-normalized window: spectral total tracks the analyzed variance
    ratio = ds.attrs["total_slope_variance"] / ds.attrs["mss_crop"]
    assert 0.9 < ratio < 1.1


def test_spectral_normalization_matlab_parity():
    from conftest import broadband_wave_stack
    sx, sy, meta = broadband_wave_stack()
    ds = compute_wave_spectra(sx, sy, dx_m=meta["dx"], fs_hz=meta["fs"],
                              window_normalization="spectral")
    assert ds.attrs["window_normalization"] == "spectral"
    # Spectral-norm scaling biases the total high relative to 'power'
    ds_p = compute_wave_spectra(sx, sy, dx_m=meta["dx"], fs_hz=meta["fs"])
    assert (ds.attrs["total_slope_variance"]
            > ds_p.attrs["total_slope_variance"])
