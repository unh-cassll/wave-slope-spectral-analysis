"""Reduce the ASIT2019 raw DoFP record to an orthorectified slope-field
stack with the E-PSS pipeline, and save it as a metadata-tagged NetCDF
ready for slopespectra ingestion.

Run with the polarimetric-slope-sensing venv, from that repo's root:

    cd ../polarimetric-slope-sensing
    ./.venv/bin/python ../wave-slope-spectral-analysis/examples/run_epss_reduction.py
"""

import sys
from pathlib import Path

import numpy as np
import xarray as xr

PSS_REPO = Path(__file__).resolve().parents[2] / "polarimetric-slope-sensing"
sys.path.insert(0, str(PSS_REPO))

from eta_field_recon.eta_pipeline import reconstruct_eta_from_record  # noqa: E402

RECORD = PSS_REPO / "_data" / "asit_2019_raw_pol_stack.nc"
GAIN_REF = PSS_REPO / "_data" / "asit_2019_raw_pol_median.nc"
OUT = Path(__file__).resolve().parents[1] / "_data" / "asit2019_slope_stack_ds8.nc"

REDUCE_DOWNSAMPLE = 8
FS_HZ = 30.0
CAMERA_AZIMUTH_DEG = 190.0


def main():
    result = reconstruct_eta_from_record(
        RECORD,
        orthorectify=True,
        reduce_downsample=REDUCE_DOWNSAMPLE,
        gain_mode="empirical",
        gain_reference_path=GAIN_REF,
        return_slopes=True,
        short_wave=False,
        force_long_wave=False,
        verbose=True,
    )

    sx = np.asarray(result.slope_x, dtype=np.float32)   # (T, Ny, Nx)
    sy = np.asarray(result.slope_y, dtype=np.float32)
    dx = float(result.slope_dx_m)
    nt, ny, nx = sx.shape
    print(f"slope stack {sx.shape}, dx = {dx*1e3:.2f} mm")

    ds = xr.Dataset(
        {
            "slope_x": (("time", "y", "x"), sx),
            "slope_y": (("time", "y", "x"), sy),
            "fs": FS_HZ,
            "camera_azimuth": CAMERA_AZIMUTH_DEG,
        },
        coords={
            "time": np.arange(nt) / FS_HZ,
            "y": np.arange(ny) * dx,
            "x": np.arange(nx) * dx,
        },
        attrs={
            "title": "ASIT2019 orthorectified wave slope stack (E-PSS)",
            "source_file": RECORD.name,
            "dx_m": dx,
            "processing_gain_mode": "empirical",
            "reduce_downsample": REDUCE_DOWNSAMPLE,
            "orthorectified": "yes (static geometry)",
        },
    )
    OUT.parent.mkdir(exist_ok=True)
    encoding = {v: {"zlib": True, "complevel": 4} for v in ("slope_x", "slope_y")}
    ds.to_netcdf(OUT, encoding=encoding)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
