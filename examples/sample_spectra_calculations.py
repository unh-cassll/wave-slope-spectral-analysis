"""Spectral analysis applied to the georectified sample wave image stack.

Python port of sample_spectra_calculations.m: loads the bundled slope-angle
stacks, computes the wavenumber-frequency directional spectrum and all
sub-spectra, and reproduces the driver's figures.

    python examples/sample_spectra_calculations.py
"""

import sys
from pathlib import Path

import h5py
import numpy as np
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))   # allow running without installation

from slopespectra import (compute_wave_spectra, figure_style,  # noqa: E402
                          plot_kf_slice, plot_omnidirectional,
                          plot_sub_spectra)
DATA = REPO / "_data"
FIGS = REPO / "_figures"

FPS = 15.0
M_PER_PX = 2 * 0.01074   # from the Motion-Correction-Georeferencing Library
GRAV = 9.81


def load_angle_stack(path, var):
    """(ny, nx, nt) slope stack from a MATLAB v7.3 milli-degree angle file."""
    with h5py.File(path) as h:
        ang = np.transpose(np.asarray(h[var]), (2, 1, 0))
    return np.tan(np.deg2rad(ang.astype(float) * 1e-3))


def main():
    FIGS.mkdir(exist_ok=True)
    sx = load_angle_stack(DATA / "angle_stack_x.mat", "ax_int")
    sy = load_angle_stack(DATA / "angle_stack_y.mat", "ay_int")

    # Heading is zero: absolute heading is already handled in rectification
    ds, dirspect = compute_wave_spectra(sx, sy, dx_m=M_PER_PX, fs_hz=FPS,
                                        heading_deg=0.0, n_freq_avg=3, temporal_alpha=0.1,
                                        n_nu_log=48, return_dirspect=True)
    print(f"mss full field = {ds.attrs['mss']:.4f}; analyzed crop = "
          f"{ds.attrs['mss_crop']:.4f}; spectral total = "
          f"{ds.attrs['total_slope_variance']:.4f}")
    print(f"k_cutoff = {ds.attrs['k_cutoff']:.1f} rad/m")

    figure_style()

    fig = plot_sub_spectra(ds)
    fig.savefig(FIGS / "sample_sub_spectra.png", dpi=200)

    # Compare F(f) with F(f) mapped from F(k) through linear dispersion
    ax = plot_omnidirectional(ds)
    ax.figure.savefig(FIGS / "sample_Ff_comparison.png", dpi=200)

    # k-f slices along both axes with the gravity-capillary dispersion curve
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))
    plot_kf_slice(dirspect, axis="ky", klim=100, ax=axs[0])
    plot_kf_slice(dirspect, axis="kx", klim=100, ax=axs[1])
    fig.tight_layout()
    fig.savefig(FIGS / "sample_kf_slices.png", dpi=200)

    print(f"figures written to {FIGS}")


if __name__ == "__main__":
    main()
