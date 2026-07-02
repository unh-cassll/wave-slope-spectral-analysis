"""Integration test against the bundled georectified sample stack,
mirroring sample_spectra_calculations.m."""

from pathlib import Path

import numpy as np
import pytest

DATA_DIR = Path(__file__).resolve().parents[1] / "_data"

pytestmark = pytest.mark.sample_data

h5py = pytest.importorskip("h5py")

FPS = 15.0
M_PER_PX = 2 * 0.01074


@pytest.fixture(scope="module")
def sample_stacks():
    fx = DATA_DIR / "angle_stack_x.mat"
    fy = DATA_DIR / "angle_stack_y.mat"
    if not (fx.exists() and fy.exists()):
        pytest.skip("sample _data stacks not present")
    # MATLAB v7.3 files are HDF5 with reversed (column-major) axis order
    with h5py.File(fx) as h:
        ax = np.transpose(np.asarray(h["ax_int"]), (2, 1, 0))
    with h5py.File(fy) as h:
        ay = np.transpose(np.asarray(h["ay_int"]), (2, 1, 0))
    sx = np.tan(np.deg2rad(ax.astype(float) * 1e-3))
    sy = np.tan(np.deg2rad(ay.astype(float) * 1e-3))
    return sx, sy


def test_sample_spectra_calculations(sample_stacks):
    from slopespectra import compute_wave_spectra

    sx, sy = sample_stacks
    assert sx.shape == (349, 373, 300)

    mss_field = np.nanvar(sx) + np.nanvar(sy)

    ds, dirspect = compute_wave_spectra(sx, sy, dx_m=M_PER_PX, fs_hz=FPS,
                                        return_dirspect=True)

    # Power-of-two center crop
    assert ds.attrs["framesize"] == 256
    assert ds.attrs["mss"] == pytest.approx(mss_field, rel=1e-6)

    # Sane, finite spectra
    for var in ("S_k", "F_k", "B_k", "S_f", "F_f"):
        vals = ds[var].values
        assert np.isfinite(vals).any()
        assert np.nanmin(vals) >= 0

    # Windowed total slope variance within a factor of two of the field mss
    assert ds.attrs["total_slope_variance"] == pytest.approx(mss_field,
                                                             rel=1.0)

    # Short gravity waves dominate: spectral peak below 5 Hz and within
    # the resolved wavenumber band
    f_pk = float(ds["f"][int(np.argmax(ds["S_f"].values))])
    assert 0.1 < f_pk < 5.0
    assert ds.attrs["k_cutoff"] > float(ds["k"].min())

    # Directional spectrum integrates back to the omnidirectional one
    dtheta = ds.attrs["dtheta"]
    S_k_re = (ds["S_k_theta"] * ds["k"]).sum("theta") * dtheta
    ratio = S_k_re.values / ds["S_k"].values
    good = np.isfinite(ratio) & (ds["S_k"].values > 0)
    assert np.allclose(ratio[good], 1.0, rtol=1e-6)
