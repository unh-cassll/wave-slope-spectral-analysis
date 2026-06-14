"""Upwind/downwind dispersion-shell analysis of wavenumber-frequency spectra.

Technique following Plant & Wright [1980], JGR: Oceans: the difference in
phase speed between upwind- and downwind-propagating waves of the same
wavenumber carries the Doppler shift of the near-surface current.
"""

import warnings

import numpy as np
import xarray as xr
from scipy.ndimage import rotate as _rotate

from .spectrum import compute_slope_spectrum
from .windows import circular_tukey

__all__ = ["upwind_downwind_wave_frequencies", "moving_window_fourier"]

# Analysis constants (Plant & Wright shell extraction)
HIGH_WAVENUMBER_CUTOFF = 50.0   # rad/m
LOW_WAVENUMBER_CUTOFF = 0.0     # rad/m
MAX_PHASE_SPEED = 5.0           # m/s, above f = 1 Hz
MAX_FREQUENCY = 10.0            # Hz
HALF_SECTOR_DEG = 45.0          # downwind/upwind sector half-width
PERCENTILE_CUTOFF = 99.9        # spectral energy percentile per ring


def upwind_downwind_wave_frequencies(Skf, dx_m, fs_hz):
    """Energy-weighted frequency and direction of the dominant up/downwind
    spectral ridges, per wavenumber ring.

    Args:
        Skf   : (m, m, nf) directional spectrum on positive frequencies
                f = (1..nf) * fs / (2 nf), axes (ky descending, kx, f)
        dx_m  : pixel size [m]
        fs_hz : frame rate [Hz]

    Returns:
        xarray.Dataset with coord k and variables f_down, f_up [Hz],
        dir_down, dir_up [deg], S_down, S_up.
    """
    Skf = np.asarray(Skf, dtype=float)
    s1 = Skf.shape[0]
    nf = Skf.shape[2]

    f_obs = np.arange(1, nf + 1) * fs_hz / (2 * nf)
    kmin = 2 * np.pi / (s1 * dx_m)
    kmax = np.pi / dx_m

    # Wavenumber label arrays (rows = ky descending)
    edges = np.arange(-kmax, kmax + kmin / 2, kmin)[1:]
    kx = np.tile(edges, (s1, 1))
    ky = np.flipud(np.tile(edges[:, None], (1, s1)))

    # Trim to the high-wavenumber cutoff
    keep = np.abs(edges) < HIGH_WAVENUMBER_CUTOFF
    kx = kx[np.ix_(keep, keep)]
    ky = ky[np.ix_(keep, keep)]
    S = Skf[np.ix_(keep, keep, np.arange(nf))]
    k = np.sqrt(kx**2 + ky**2)

    fmat = np.broadcast_to(f_obs[None, None, :], S.shape)
    with np.errstate(divide="ignore", invalid="ignore"):
        Cp = 2 * np.pi * fmat / k[:, :, None]

    # Mask outside the wavenumber/phase-speed/frequency/direction limits
    mask = np.where((k < HIGH_WAVENUMBER_CUTOFF)
                    & (k >= LOW_WAVENUMBER_CUTOFF), 1.0, np.nan)
    mask3 = np.repeat(mask[:, :, None], nf, axis=2)
    mask3[(Cp > MAX_PHASE_SPEED) & (fmat > 1.0)] = np.nan
    mask3[fmat > MAX_FREQUENCY] = np.nan

    angle_down = np.rad2deg(np.arctan2(kx, ky))
    angle_up = np.flipud(angle_down)
    bad_sector = (np.abs(angle_down) > HALF_SECTOR_DEG) \
        & (np.abs(angle_up) > HALF_SECTOR_DEG)
    mask3[np.repeat(bad_sector[:, :, None], nf, axis=2)] = np.nan

    S_masked = mask3 * S
    rows_down = ky[:, 0] > 0
    rows_up = ky[:, 0] < 0
    S_down = S_masked.copy()
    S_up = S_masked.copy()
    S_down[rows_up, :, :] = np.nan
    S_up[rows_down, :, :] = np.nan

    k_vec = edges[edges > 0]
    kL = len(k_vec)
    r = np.floor(k / np.nanmax(k_vec) * kL).astype(int)

    out = {name: np.full(kL, np.nan)
           for name in ("f_down", "f_up", "dir_down", "dir_up",
                        "S_down", "S_up")}
    angle_down3 = np.repeat(angle_down[:, :, None], nf, axis=2)
    angle_up3 = np.repeat(angle_up[:, :, None], nf, axis=2)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for i in range(1, kL + 1):
            ring = r == i
            for half, S_half, ang3, suffix in (
                    ("down", S_down, angle_down3, "_down"),
                    ("up", S_up, angle_up3, "_up")):
                vals = S_half[ring, :]
                finite = vals[np.isfinite(vals)]
                if finite.size == 0:
                    continue
                cutoff = np.sort(finite)[
                    int(np.floor(PERCENTILE_CUTOFF / 100 * finite.size)) - 1]
                sel = vals > cutoff
                Ssel = vals[sel]
                if Ssel.size == 0:
                    continue
                fsel = np.broadcast_to(f_obs, vals.shape)[sel]
                tsel = ang3[ring, :][sel]
                out["S" + suffix][i - 1] = np.nanmean(Ssel)
                out["f" + suffix][i - 1] = (np.nanmean(fsel * Ssel)
                                            / np.nanmean(Ssel))
                out["dir" + suffix][i - 1] = (np.nanmean(tsel * Ssel)
                                              / np.nanmean(Ssel))

    return xr.Dataset(
        {name: (("k",), arr) for name, arr in out.items()},
        coords={"k": ("k", k_vec, {"units": "rad m-1"})},
        attrs={"dx_m": dx_m, "fs_hz": fs_hz},
    )


def moving_window_fourier(sx, sy, dx_m, fs_hz, window_length_s, window_step_s,
                          ccw_rot_angle_deg=0.0, taper_width=0.1,
                          verbose=False):
    """Moving-window Doppler-shell analysis of a slope-field stack.

    For each temporal window, computes the directional wavenumber-frequency
    spectrum, rotates it so the wind direction is along +ky, and extracts
    the upwind/downwind ridge frequencies and the omnidirectional spectrum.

    Args:
        sx, sy            : (ny, nx, nt) slope-field stacks
        dx_m              : pixel size [m]
        fs_hz             : frame rate [Hz]
        window_length_s   : window duration [s]
        window_step_s     : window step [s]
        ccw_rot_angle_deg : CCW spectrum rotation aligning wind with +ky [deg]
        taper_width       : circular Tukey taper fraction

    Returns:
        xarray.Dataset over (t, k) with f_down, f_up, dir_down, dir_up,
        S_down, S_up and the per-window omnidirectional spectrum S_omni.
    """
    sx = circular_tukey(np.asarray(sx, dtype=float), taper_width)
    sy = circular_tukey(np.asarray(sy, dtype=float), taper_width)
    s1, _, s3 = sx.shape

    frames_per_bit = int(fs_hz * window_length_s)
    frames_per_bit -= frames_per_bit % 2
    frames_per_shift = max(int(np.floor(fs_hz * window_step_s)), 1)
    starts = np.arange(0, s3 - frames_per_bit + 1, frames_per_shift)

    records = []
    for n, istart in enumerate(starts):
        xbit = sx[:, :, istart:istart + frames_per_bit]
        ybit = sy[:, :, istart:istart + frames_per_bit]
        spect = compute_slope_spectrum(xbit, ybit, dx_m, fs_hz, framesize=s1)
        Skf = spect.Skf
        if ccw_rot_angle_deg:
            Skf = _rotate(np.nan_to_num(Skf, nan=0.0), ccw_rot_angle_deg,
                          axes=(1, 0), reshape=False, order=1)
        ds = upwind_downwind_wave_frequencies(Skf, dx_m, fs_hz)
        ds["dir_down"] = np.mod(ds["dir_down"] + ccw_rot_angle_deg, 360.0)
        ds["dir_up"] = np.mod(ds["dir_up"] + ccw_rot_angle_deg + 180.0, 360.0)
        ds = ds.assign(S_omni=(("k_omni",), spect.omni_S))
        ds = ds.assign_coords(k_omni=("k_omni", spect.omni_k))
        ds = ds.expand_dims(t=[istart / fs_hz])
        records.append(ds)
        if verbose:
            print(f"window {n + 1}/{len(starts)} done")

    return xr.concat(records, dim="t")
