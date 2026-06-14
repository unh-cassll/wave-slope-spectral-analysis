"""Slope-field ingestion: plain arrays, xarray Datasets, NetCDF records, or
E-PSS pipeline results, normalized to a common (ny, nx, nt) container."""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import xarray as xr

__all__ = ["SlopeStack", "load_slope_stack"]

# Variable/dimension/attribute name candidates in E-PSS-style datasets
_SX_NAMES = ("slope_x", "sx", "slopefield_stack_x", "slope_east", "sE")
_SY_NAMES = ("slope_y", "sy", "slopefield_stack_y", "slope_north", "sN")
_TIME_DIMS = ("time", "t", "frame")
_FS_NAMES = ("fs", "fs_hz", "framerate", "framerate_hz", "sampling_rate")
_DX_NAMES = ("dx", "dx_m", "slope_dx_m", "ground_dx_m", "m_per_px")
_HEADING_NAMES = ("camera_azimuth", "camera_heading_deg", "heading_deg")


@dataclass
class SlopeStack:
    """Two-component slope-field stack with acquisition metadata.

    sx, sy are (ny, nx, nt) with rows as the +y edge first (image
    convention); dx_m is the ground pixel size and fs_hz the frame rate.
    """
    sx: np.ndarray
    sy: np.ndarray
    dx_m: float
    fs_hz: float | None = None
    heading_deg: float = 0.0
    attrs: dict = field(default_factory=dict)

    @property
    def shape(self):
        return self.sx.shape

    @property
    def nt(self):
        return self.sx.shape[2]

    @property
    def mss(self):
        """Total mean square slope (x-variance + y-variance)."""
        return float(np.nanvar(self.sx) + np.nanvar(self.sy))

    def center_crop(self, n=None):
        """Center-cropped copy with square spatial extent n (default: the
        largest power of two fitting both spatial dimensions)."""
        s1, s2, _ = self.sx.shape
        if n is None:
            n = 2 ** int(np.floor(np.log2(min(s1, s2))))
        n2 = n // 2
        rows = slice(s1 // 2 - n2, s1 // 2 + n - n2)
        cols = slice(s2 // 2 - n2, s2 // 2 + n - n2)
        return SlopeStack(self.sx[rows, cols, :], self.sy[rows, cols, :],
                          self.dx_m, self.fs_hz, self.heading_deg,
                          dict(self.attrs))

    def fill_nan(self, value=0.0):
        """Copy with NaNs replaced by `value` (after recording mss)."""
        attrs = dict(self.attrs)
        attrs.setdefault("mss_before_fill", self.mss)
        return SlopeStack(np.nan_to_num(self.sx, nan=value),
                          np.nan_to_num(self.sy, nan=value),
                          self.dx_m, self.fs_hz, self.heading_deg, attrs)

    def to_dataset(self):
        """xarray.Dataset representation with acquisition metadata attrs."""
        ny, nx, nt = self.sx.shape
        coords = {
            "y": ("y", np.arange(ny) * self.dx_m, {"units": "m"}),
            "x": ("x", np.arange(nx) * self.dx_m, {"units": "m"}),
        }
        if self.fs_hz:
            coords["time"] = ("time", np.arange(nt) / self.fs_hz,
                              {"units": "s"})
        else:
            coords["time"] = ("time", np.arange(nt))
        return xr.Dataset(
            {
                "slope_x": (("y", "x", "time"), self.sx),
                "slope_y": (("y", "x", "time"), self.sy),
            },
            coords=coords,
            attrs={"dx_m": self.dx_m, "fs_hz": self.fs_hz,
                   "camera_heading_deg": self.heading_deg, **self.attrs},
        )


def _first_match(candidates, available):
    for name in candidates:
        for avail in available:
            if str(avail).lower() == name.lower():
                return avail
    return None


def _scalar_from_dataset(ds, names):
    name = _first_match(names, list(ds.variables))
    if name is not None:
        return float(np.asarray(ds[name]).ravel()[0])
    name = _first_match(names, list(ds.attrs))
    if name is not None:
        return float(np.asarray(ds.attrs[name]).ravel()[0])
    return None


def _stack_from_dataarray(da):
    """(ny, nx, nt) array from a DataArray in any of the common layouts."""
    dims = [str(d).lower() for d in da.dims]
    tdim = _first_match(_TIME_DIMS, da.dims)
    if da.ndim == 2:
        return np.asarray(da, dtype=float)[:, :, None]
    if da.ndim != 3:
        raise ValueError(f"slope variable must be 2D or 3D, got {da.dims}")
    if tdim is None:
        # Assume trailing time axis when no dimension is named
        return np.asarray(da, dtype=float)
    other = [d for d in da.dims if d != tdim]
    return np.asarray(da.transpose(*other, tdim), dtype=float)


def _from_dataset(ds, dx_m, fs_hz, heading_deg):
    sx_name = _first_match(_SX_NAMES, list(ds.data_vars))
    sy_name = _first_match(_SY_NAMES, list(ds.data_vars))
    if sx_name is None or sy_name is None:
        raise ValueError(
            "could not locate slope variables in the dataset; expected "
            f"names like {_SX_NAMES[0]}/{_SY_NAMES[0]}")
    sx = _stack_from_dataarray(ds[sx_name])
    sy = _stack_from_dataarray(ds[sy_name])

    if dx_m is None:
        dx_m = _scalar_from_dataset(ds, _DX_NAMES)
    if dx_m is None:
        # Fall back to the spacing of a metric spatial coordinate
        sdim = _first_match(("x", "y"), list(ds.coords))
        if sdim is not None and ds[sdim].size > 1:
            dx_m = float(np.median(np.diff(np.asarray(ds[sdim], dtype=float))))
    if fs_hz is None:
        fs_hz = _scalar_from_dataset(ds, _FS_NAMES)
    if fs_hz is None:
        tdim = _first_match(_TIME_DIMS, list(ds.coords))
        if tdim is not None and ds[tdim].size > 1:
            tvals = np.asarray(ds[tdim], dtype=float)
            dt = float(np.median(np.diff(tvals)))
            if dt > 0:
                fs_hz = 1.0 / dt
    if heading_deg is None:
        heading_deg = _scalar_from_dataset(ds, _HEADING_NAMES)

    return sx, sy, dx_m, fs_hz, heading_deg, dict(ds.attrs)


def load_slope_stack(source, sy=None, *, dx_m=None, fs_hz=None,
                     heading_deg=None, time_axis=-1):
    """Normalize a slope-field source to a SlopeStack.

    Accepted sources:
      * a pair of numpy arrays (source=sx, sy=sy), (ny, nx) or 3D with the
        time axis given by `time_axis` (default trailing, MATLAB layout;
        pass time_axis=0 for E-PSS (T, Ny, Nx) stacks)
      * an xarray.Dataset with slope_x/slope_y-style variables and
        metadata (fs, dx, camera_azimuth) as variables or attrs
      * a path to a NetCDF file containing such a dataset
      * an E-PSS result object exposing slope_x/slope_y (T, Ny, Nx) arrays

    Explicit keyword arguments override any metadata found in the source.

    Returns:
        SlopeStack with sx, sy as float64 (ny, nx, nt).
    """
    attrs = {}

    if isinstance(source, (str, Path)):
        with xr.open_dataset(source) as ds:
            ds = ds.load()
        source = ds

    if isinstance(source, xr.Dataset):
        sx, sy_arr, dx_m, fs_hz, heading_deg, attrs = _from_dataset(
            source, dx_m, fs_hz, heading_deg)
    elif hasattr(source, "slope_x") and hasattr(source, "slope_y"):
        # E-PSS pipeline result (duck-typed): (T, Ny, Nx) layout
        sx = np.moveaxis(np.asarray(source.slope_x, dtype=float), 0, -1)
        sy_arr = np.moveaxis(np.asarray(source.slope_y, dtype=float), 0, -1)
        if dx_m is None:
            for name in ("slope_dx_m", "dx_m"):
                val = getattr(source, name, None)
                if val is not None:
                    dx_m = float(val)
                    break
    else:
        if sy is None:
            raise ValueError("sy is required when passing plain arrays")
        sx = np.asarray(source, dtype=float)
        sy_arr = np.asarray(sy, dtype=float)
        if sx.ndim == 2:
            sx = sx[:, :, None]
            sy_arr = sy_arr[:, :, None]
        elif time_axis != -1 and time_axis != 2:
            sx = np.moveaxis(sx, time_axis, -1)
            sy_arr = np.moveaxis(sy_arr, time_axis, -1)

    if sx.shape != sy_arr.shape:
        raise ValueError("slope component shapes differ: "
                         f"{sx.shape} vs {sy_arr.shape}")
    if dx_m is None:
        raise ValueError("dx_m could not be determined; pass dx_m=...")
    if fs_hz is None and sx.shape[-1] > 1:
        raise ValueError("fs_hz could not be determined; pass fs_hz=...")

    return SlopeStack(sx, sy_arr, float(dx_m),
                      None if fs_hz is None else float(fs_hz),
                      0.0 if heading_deg is None else float(heading_deg),
                      attrs)
