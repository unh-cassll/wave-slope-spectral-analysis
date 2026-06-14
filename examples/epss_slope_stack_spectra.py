"""Spectral analysis of an E-PSS-produced slope-field stack.

Consumes the metadata-tagged NetCDF written by run_epss_reduction.py
(ASIT2019 record reduced with the polarimetric-slope-sensing pipeline) and
produces the standard slopespectra figures.

    python examples/epss_slope_stack_spectra.py
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))   # allow running without installation

from slopespectra import (compute_wave_spectra, figure_style,  # noqa: E402
                          plot_kf_slice, plot_omnidirectional,
                          plot_sub_spectra)
STACK = REPO / "_data" / "asit2019_slope_stack_ds8.nc"
FIGS = REPO / "_figures"


def main():
    FIGS.mkdir(exist_ok=True)

    # All metadata (dx, fs, camera azimuth) rides along in the dataset.
    # Band-average 5 frequency bins (df: 1/60 -> 1/12 Hz) and use 32
    # log-spaced nu bins to reduce bin-to-bin variance.
    ds, dirspect = compute_wave_spectra(STACK, n_freq_avg=5, n_nu_log=32,
                                        temporal_alpha=0.1,
                                        return_dirspect=True)
    print(f"mss full field = {ds.attrs['mss']:.4f}; analyzed crop = "
          f"{ds.attrs['mss_crop']:.4f}; spectral total = "
          f"{ds.attrs['total_slope_variance']:.4f}")
    print(f"framesize = {ds.attrs['framesize']} px, "
          f"dx = {ds.attrs['dx_m']*1e3:.2f} mm, fs = {ds.attrs['fs_hz']} Hz")
    print(f"k_cutoff = {ds.attrs['k_cutoff']:.1f} rad/m")

    figure_style()

    fig = plot_sub_spectra(ds)
    fig.savefig(FIGS / "asit2019_sub_spectra.png", dpi=200)

    # Skip the two lowest k bins: waves longer than the 2.6 m FOV pile
    # ~95% of the elevation variance into them
    ax = plot_omnidirectional(ds, min_k_bins=2)
    ax.figure.savefig(FIGS / "asit2019_Ff_comparison.png", dpi=200)

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
    plot_kf_slice(dirspect, axis="ky", ax=axs[0])
    plot_kf_slice(dirspect, axis="kx", ax=axs[1])
    fig.tight_layout()
    fig.savefig(FIGS / "asit2019_kf_slices.png", dpi=200)

    # Directional distribution at the dominant wavenumber vs wind direction
    # (wind from 181.9 deg -> waves going toward ~2 deg true)
    theta_deg = np.rad2deg(ds["theta"].values)
    D = (ds["S_k_theta"] * ds["k"]).sum("k").values
    _, ax = plt.subplots(figsize=(6, 4))
    ax.plot(theta_deg, D / D.max(), linewidth=2)
    ax.axvline(2.0, color="k", linestyle="--", linewidth=1,
               label="wind going-to direction")
    ax.set_xlabel(r"$\theta$ [$^\circ$ true]")
    ax.set_ylabel("normalized directional slope variance")
    ax.legend()
    ax.figure.tight_layout()
    ax.figure.savefig(FIGS / "asit2019_directional_distribution.png", dpi=200)

    print(f"figures written to {FIGS}")


if __name__ == "__main__":
    main()
