# wave-slope-spectral-analysis

Fourier analysis of short ocean wave slope data in order to obtain wavenumber-frequency directional spectra and associated products.

The repository contains two parallel implementations:

- **`slopespectra/`** — a pip/uv-installable Python package (the maintained implementation)
- **`matlab/`** — the original MATLAB implementation (`sample_spectra_calculations.m` driver + `_functions/` library)

## Installation

```bash
pip install .            # or: uv pip install .
pip install '.[all]'     # + NetCDF ingestion, sample-data readers, pytest
```

## What it computes

Given a time series of two-component water surface slope fields — as produced by polarimetric slope sensing ([polarimetric-slope-sensing](https://github.com/unh-cassll/polarimetric-slope-sensing)) or any other slope-imaging modality — `slopespectra` computes:

- the directional wavenumber-frequency slope spectrum S(kx, ky, f) (3D FFT; power-normalized circular Tukey taper with optional temporal taper, so the spectral total is an unbiased estimate of the analyzed-region slope variance)
- polar-projected spectra S(k, θ, f) with per-bin and global variance conservation; each wavenumber ring is normalized by its analytic area (2πk·dk) rather than the fluctuating grid-cell count, so the radial spectra are free of the discrete-annulus (Gauss circle) sawtooth
- wavenumber/frequency directional and omnidirectional spectra for slope S, elevation F = S/k², and saturation B = k²·S
- inverse-phase-speed spectra Q(ν, θ)
- dispersion-shell diagnostics: ridge frequency f(k, θ), upwind/downwind phase speed difference Δc (Plant & Wright [1980]), and the resolved-wavenumber cutoff
- moving-window Doppler analysis, 1D band-averaged auto/cross spectra (`spectf`), dispersive bandstop and notch filters, and Gram-Charlier slope-PDF fits (Cox & Munk [1954])

Results are returned as metadata-tagged `xarray.Dataset`s.

## Quick start

```python
from slopespectra import compute_wave_spectra

# Dumb arrays: (ny, nx, nt) stacks plus explicit metadata
# (n_freq_avg band-averages adjacent frequency bins; n_nu_log selects
# log-spaced inverse-phase-speed bins with cloud-in-cell scatter)
ds = compute_wave_spectra(sx, sy, dx_m=0.0215, fs_hz=15.0,
                          n_freq_avg=3, n_nu_log=48)

# Or an E-PSS-style metadata-tagged dataset / NetCDF path: dx, fs and
# camera azimuth are discovered from the file
ds = compute_wave_spectra("slope_stack.nc")

ds["S_k"]                          # omnidirectional slope spectrum
ds["F_f_theta"]                    # frequency-directional elevation spectrum
ds.attrs["k_cutoff"]               # resolved-wavenumber cutoff [rad/m]
```

Slope sources accepted by `compute_wave_spectra` / `load_slope_stack`:

1. plain numpy array pairs (trailing or leading time axis)
2. `xarray.Dataset` with `slope_x`/`slope_y`-style variables and metadata as variables or attrs
3. a NetCDF path to such a dataset
4. an E-PSS pipeline result object (`run_epss` / `reconstruct_eta_from_record` with `return_slopes=True`)

## Plotting

`slopespectra.plotting` carries the E-PSS paper figure style:

```python
from slopespectra import figure_style, plot_sub_spectra, plot_kf_slice

fig = plot_sub_spectra(ds)                 # six-panel overview
ax = plot_kf_slice(dirspect, axis="kx")    # k-f slice + dispersion curve
```

## Examples

- `examples/sample_spectra_calculations.py` — port of the MATLAB demo driver; runs on the bundled `_data/angle_stack_*.mat` sample stacks
- `examples/run_epss_reduction.py` — reduces the ASIT2019 raw DoFP record with the E-PSS pipeline into a metadata-tagged slope-stack NetCDF
- `examples/epss_slope_stack_spectra.py` — spectral analysis and figures from that E-PSS product

## Tests

```bash
python -m pytest                      # unit + sample-data integration tests
python -m pytest -m "not sample_data" # skip the bundled-data test
```

Conventions: stacks are `(ny, nx, nt)` with row 0 at the +y edge (image convention); θ is the compass-style going-to wave direction (0° = +y, 90° = +x) plus the camera heading; spectra are one-sided in frequency.
