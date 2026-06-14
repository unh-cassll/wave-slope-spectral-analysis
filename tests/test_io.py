import numpy as np
import pytest
import xarray as xr

from slopespectra import SlopeStack, load_slope_stack


@pytest.fixture
def arrays():
    rng = np.random.default_rng(0)
    sx = rng.standard_normal((32, 40, 10))
    sy = rng.standard_normal((32, 40, 10))
    return sx, sy


def test_plain_arrays(arrays):
    sx, sy = arrays
    stack = load_slope_stack(sx, sy, dx_m=0.02, fs_hz=10.0)
    assert stack.shape == (32, 40, 10)
    assert stack.dx_m == 0.02
    assert stack.fs_hz == 10.0
    assert stack.heading_deg == 0.0


def test_plain_arrays_leading_time_axis(arrays):
    sx, sy = arrays
    stack = load_slope_stack(np.moveaxis(sx, -1, 0), np.moveaxis(sy, -1, 0),
                             dx_m=0.02, fs_hz=10.0, time_axis=0)
    assert stack.shape == (32, 40, 10)
    assert np.allclose(stack.sx, sx)


def test_2d_arrays(arrays):
    sx, sy = arrays
    stack = load_slope_stack(sx[:, :, 0], sy[:, :, 0], dx_m=0.02)
    assert stack.shape == (32, 40, 1)
    assert stack.fs_hz is None


def test_missing_metadata_raises(arrays):
    sx, sy = arrays
    with pytest.raises(ValueError, match="dx_m"):
        load_slope_stack(sx, sy)
    with pytest.raises(ValueError, match="fs_hz"):
        load_slope_stack(sx, sy, dx_m=0.02)
    with pytest.raises(ValueError, match="sy is required"):
        load_slope_stack(sx, dx_m=0.02, fs_hz=10.0)


def test_xarray_dataset_with_metadata(arrays):
    sx, sy = arrays
    ds = xr.Dataset(
        {
            "slope_x": (("time", "y", "x"), np.moveaxis(sx, -1, 0)),
            "slope_y": (("time", "y", "x"), np.moveaxis(sy, -1, 0)),
            "fs": 10.0,
            "camera_azimuth": 247.0,
        },
        attrs={"dx_m": 0.02},
    )
    stack = load_slope_stack(ds)
    assert stack.shape == (32, 40, 10)
    assert np.allclose(stack.sx, sx)
    assert stack.fs_hz == 10.0
    assert stack.dx_m == 0.02
    assert stack.heading_deg == 247.0


def test_xarray_fs_from_time_coordinate(arrays):
    sx, sy = arrays
    ds = xr.Dataset(
        {
            "slope_x": (("y", "x", "time"), sx),
            "slope_y": (("y", "x", "time"), sy),
        },
        coords={"time": np.arange(10) / 5.0,
                "x": np.arange(40) * 0.03,
                "y": np.arange(32) * 0.03},
    )
    stack = load_slope_stack(ds)
    assert stack.fs_hz == pytest.approx(5.0)
    assert stack.dx_m == pytest.approx(0.03)


def test_explicit_kwargs_override_dataset(arrays):
    sx, sy = arrays
    ds = xr.Dataset(
        {
            "slope_x": (("y", "x", "time"), sx),
            "slope_y": (("y", "x", "time"), sy),
            "fs": 10.0,
        },
        attrs={"dx_m": 0.02},
    )
    stack = load_slope_stack(ds, dx_m=0.05, fs_hz=2.0, heading_deg=90.0)
    assert stack.dx_m == 0.05
    assert stack.fs_hz == 2.0
    assert stack.heading_deg == 90.0


def test_netcdf_roundtrip(tmp_path, arrays):
    pytest.importorskip("netCDF4")
    sx, sy = arrays
    original = SlopeStack(sx, sy, dx_m=0.02, fs_hz=10.0, heading_deg=12.0)
    path = tmp_path / "slopes.nc"
    original.to_dataset().to_netcdf(path)
    stack = load_slope_stack(path)
    assert np.allclose(stack.sx, sx)
    assert stack.dx_m == pytest.approx(0.02)
    assert stack.fs_hz == pytest.approx(10.0)
    assert stack.heading_deg == pytest.approx(12.0)


def test_epss_result_ducktype(arrays):
    sx, sy = arrays

    class FakeEpssResult:
        slope_x = np.moveaxis(sx, -1, 0)
        slope_y = np.moveaxis(sy, -1, 0)
        slope_dx_m = 0.0228

    stack = load_slope_stack(FakeEpssResult(), fs_hz=30.0)
    assert stack.shape == (32, 40, 10)
    assert np.allclose(stack.sx, sx)
    assert stack.dx_m == pytest.approx(0.0228)


def test_center_crop_and_fill():
    rng = np.random.default_rng(1)
    sx = rng.standard_normal((70, 90, 4))
    sx[0, 0, 0] = np.nan
    stack = SlopeStack(sx, sx.copy(), dx_m=0.01, fs_hz=1.0)
    cropped = stack.center_crop()
    assert cropped.shape == (64, 64, 4)
    filled = stack.fill_nan()
    assert np.isfinite(filled.sx).all()
    assert "mss_before_fill" in filled.attrs


def test_mss():
    sx = np.full((8, 8, 2), 2.0)
    sx[..., 1] = -2.0
    stack = SlopeStack(sx, np.zeros_like(sx), dx_m=0.01, fs_hz=1.0)
    assert stack.mss == pytest.approx(4.0)
