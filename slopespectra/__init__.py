"""slopespectra: Fourier analysis of ocean wave slope fields.

Python port of the wave-slope-spectral-analysis MATLAB library: given
two-component water surface slope-field stacks (e.g. from polarimetric
slope sensing / E-PSS), computes wavenumber-frequency directional spectra
and associated products.

Typical use:

    from slopespectra import compute_wave_spectra
    ds = compute_wave_spectra(sx, sy, dx_m=0.02148, fs_hz=15.0)
"""

from .io import SlopeStack, load_slope_stack
from .windows import circular_tukey, blackman_harris, taper_2d
from .spectrum import (compute_slope_spectrum, compute_segment_power,
                       finish_spectrum, azimuthal_integral,
                       wavenumber_grids, DirectionalSpectrum)
from .polar import polar_from_cartesian, PolarBinner
from .subspectra import compute_sub_spectra
from .pipeline import compute_wave_spectra
from .doppler import upwind_downwind_wave_frequencies, moving_window_fourier
from .filters import filter_dispersive_bandstop, notch_filter
from .onedim import spectf
from .gram_charlier import gram_charlier_pdf, fit_gram_charlier_slope_pdf
from .plotting import (figure_style, plot_sub_spectra, plot_kf_slice,
                       plot_omnidirectional)

__version__ = "0.1.0"

__all__ = [
    "SlopeStack", "load_slope_stack",
    "circular_tukey", "blackman_harris", "taper_2d",
    "compute_slope_spectrum", "compute_segment_power", "finish_spectrum",
    "azimuthal_integral", "wavenumber_grids",
    "DirectionalSpectrum",
    "polar_from_cartesian", "PolarBinner",
    "compute_sub_spectra", "compute_wave_spectra",
    "upwind_downwind_wave_frequencies", "moving_window_fourier",
    "filter_dispersive_bandstop", "notch_filter",
    "spectf",
    "gram_charlier_pdf", "fit_gram_charlier_slope_pdf",
    "figure_style", "plot_sub_spectra", "plot_kf_slice",
    "plot_omnidirectional",
]
