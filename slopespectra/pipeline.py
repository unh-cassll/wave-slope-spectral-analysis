"""End-to-end driver: slope-field source to sub-spectra Dataset."""

import numpy as np

from .io import load_slope_stack
from .spectrum import compute_slope_spectrum
from .subspectra import compute_sub_spectra
from .windows import circular_tukey

__all__ = ["compute_wave_spectra"]


def compute_wave_spectra(source, sy=None, *, dx_m=None, fs_hz=None,
                         heading_deg=None, time_axis=-1, framesize=None,
                         taper_width=0.2, temporal_alpha=0.0,
                         window_normalization="power",
                         dtheta=np.deg2rad(5.0),
                         n_freq_avg=1, dnu_base=0.01, n_nu_log=None,
                         return_dirspect=False):
    """Wavenumber-frequency directional spectrum and all sub-spectra from a
    slope-field source.

    Mirrors the MATLAB sample driver: center-crop to a power-of-two square,
    circular Tukey taper, 3D FFT, polar projection, sub-spectra.

    Args:
        source, sy   : any source accepted by load_slope_stack
        dx_m, fs_hz, heading_deg, time_axis : metadata overrides
        framesize    : square center-crop / FFT size, at most the spatial
                       extent; default largest power of two fitting it
        taper_width  : circular Tukey taper fraction
        temporal_alpha : Tukey cosine fraction along t (0 disables)
        window_normalization : 'power' (unbiased variance; default),
                       'spectral' (MATLAB parity) or 'none'
        dtheta       : direction bin width [rad]
        n_freq_avg   : adjacent frequency bins averaged before the polar
                       projection (band averaging)
        dnu_base     : base inverse-phase-speed bin width [s/m]
        n_nu_log     : contiguous log-spaced nu bin count (overrides
                       dnu_base when set)
        return_dirspect : also return the intermediate DirectionalSpectrum

    Returns:
        xarray.Dataset from compute_sub_spectra, with mean-square-slope and
        acquisition metadata in attrs (and the DirectionalSpectrum when
        return_dirspect=True).
    """
    stack = load_slope_stack(source, sy, dx_m=dx_m, fs_hz=fs_hz,
                             heading_deg=heading_deg, time_axis=time_axis)
    if stack.fs_hz is None:
        raise ValueError("a multi-frame record with fs_hz is required")

    mss_x = float(np.nanvar(stack.sx))
    mss_y = float(np.nanvar(stack.sy))
    stack = stack.fill_nan(0.0)
    stack = stack.center_crop(framesize)
    n_spect = stack.sx.shape[0]
    nt = stack.nt
    nt -= nt % 2

    # Remove the static (per-pixel temporal-mean) slope pattern: it carries
    # no wave signal on the retained f > 0 planes, and a temporal taper
    # would otherwise modulate it into them
    sx = stack.sx[:, :, :nt]
    sy_w = stack.sy[:, :, :nt]
    sx = sx - sx.mean(axis=2, keepdims=True)
    sy_w = sy_w - sy_w.mean(axis=2, keepdims=True)

    # Wave variance of the analyzed (cropped, demeaned) region: the
    # reference against which the spectral totals should be checked
    mss_crop = float(sx.var() + sy_w.var())

    sx = circular_tukey(sx, taper_width,
                        normalization=window_normalization,
                        temporal_alpha=temporal_alpha)
    sy_w = circular_tukey(sy_w, taper_width,
                          normalization=window_normalization,
                          temporal_alpha=temporal_alpha)

    dirspect = compute_slope_spectrum(sx, sy_w, stack.dx_m, stack.fs_hz,
                                      framesize=n_spect)
    ds = compute_sub_spectra(dirspect.Skf, dirspect.dk, dirspect.df,
                             heading_deg=stack.heading_deg, dtheta=dtheta,
                             n_freq_avg=n_freq_avg, dnu_base=dnu_base,
                             n_nu_log=n_nu_log)
    ds.attrs.update({
        "mss_x": mss_x,
        "mss_y": mss_y,
        "mss": mss_x + mss_y,
        "mss_crop": mss_crop,
        "dx_m": stack.dx_m,
        "fs_hz": stack.fs_hz,
        "framesize": n_spect,
        "taper_width": taper_width,
        "temporal_alpha": temporal_alpha,
        "window_normalization": window_normalization,
    })
    if return_dirspect:
        return ds, dirspect
    return ds
