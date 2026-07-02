"""Sub-spectra derived from the wavenumber-frequency directional spectrum.

Given S(kx, ky, f), computes:

* wavenumber directional spectra: slope S(k, theta), elevation F(k, theta),
  saturation B(k, theta), and their omnidirectional integrals
* frequency directional spectra and omnidirectional integrals
* inverse phase speed directional spectra Q(nu, theta) and integrals
* dispersion-shell phase speed difference about +/-180 degrees
  (Plant & Wright [1980])
"""

import warnings

import numpy as np
import xarray as xr
from scipy.ndimage import median_filter

from .polar import PolarBinner
from .spectrum import wavenumber_grids

__all__ = ["compute_sub_spectra"]

# Frequency band retained for the dispersion-shell analysis [Hz]
F_LOW_HZ = 0.0
F_HIGH_HZ = 15.0


def _inverse_phase_speed_grid(dnu_base=0.01, n_nu_log=None,
                              nu_min=0.05, nu_max=4.0):
    """Inverse phase speed bin centers, widths and edges [s/m].

    Default: piecewise-uniform grid with width doubling each octave above
    nu = 0.5. With n_nu_log set: contiguous log-spaced bins on
    [nu_min, nu_max], whose nu-proportional widths suppress empty-bin
    noise where spectral cells are sparse.
    """
    if n_nu_log is not None:
        edges = np.geomspace(nu_min, nu_max, n_nu_log + 1)
        nu = np.sqrt(edges[:-1] * edges[1:])
        dnu = np.diff(edges)
        return nu, dnu, edges[:-1], edges[1:]
    nu = np.concatenate([
        np.arange(0.05, 0.5 + dnu_base / 2, dnu_base),
        np.arange(0.5 + dnu_base * 2, 1 + dnu_base, dnu_base * 2),
        np.arange(1 + dnu_base * 4, 2 + dnu_base * 2, dnu_base * 4),
        np.arange(2 + dnu_base * 8, 4 + dnu_base * 4, dnu_base * 8),
    ])
    dnu = np.concatenate([[dnu_base], np.diff(nu)])
    return nu, dnu, nu - dnu / 2, nu + dnu / 2


def compute_sub_spectra(Skf, dk, df, heading_deg=0.0, dtheta=np.deg2rad(5.0),
                        n_freq_avg=1, dnu_base=0.01, n_nu_log=None):
    """Sub-spectra from the wavenumber-frequency directional slope spectrum.

    Args:
        Skf         : (m, m, nf) directional slope spectral density on
                      positive frequencies f = df * (1..nf), axes
                      (ky descending, kx ascending, f)
        dk          : wavenumber resolution [rad/m]
        df          : frequency resolution [Hz]
        heading_deg : camera heading [deg]; zero for georectified stacks
        dtheta      : direction bin width [rad]
        n_freq_avg  : number of adjacent frequency bins averaged before the
                      polar projection (band averaging; reduces bin-to-bin
                      variance at the cost of frequency resolution)
        dnu_base    : base inverse-phase-speed bin width [s/m] (the grid
                      width doubles per octave above nu = 0.5)
        n_nu_log    : if set, use this many contiguous log-spaced nu bins
                      on [0.05, 4] s/m instead of the dnu_base grid
                      (band averaging in nu; suppresses empty-bin noise)

    Returns:
        xarray.Dataset with coords k [rad/m], f [Hz], nu [s/m], theta [rad]
        and the slope (S), elevation (F), saturation (B) and inverse phase
        speed (Q) spectra, plus dispersion-shell diagnostics.
    """
    Skf = np.asarray(Skf, dtype=float)
    Nx, Ny, Nt = Skf.shape
    m = int(np.mean([Nx, Ny]))
    kx, ky = wavenumber_grids(m, 2 * np.pi / (m * dk))

    if n_freq_avg > 1:
        nf_keep = (Nt // n_freq_avg) * n_freq_avg
        # Band centers of the original f = df * (1..nf) grid
        f = df * (np.arange(nf_keep // n_freq_avg) * n_freq_avg
                  + (n_freq_avg + 1) / 2)
        Skf = Skf[:, :, :nf_keep].reshape(
            Nx, Ny, nf_keep // n_freq_avg, n_freq_avg).mean(axis=3)
        df = df * n_freq_avg
        Nt = Skf.shape[2]
    else:
        f = df * np.arange(1, Nt + 1)
    k = dk * np.arange(1, m // 2 + 1)

    # Cartesian -> polar per frequency plane
    binner = PolarBinner(kx, ky, dtheta=dtheta, heading_deg=heading_deg)
    theta = binner.theta_vec
    Skftheta = np.full((m // 2, len(theta), Nt), np.nan)
    for fnum in range(Nt):
        Skftheta[:, :, fnum] = binner(Skf[:, :, fnum])

    # Enforce variance conservation after the polar conversion
    kcol = k[:, None, None]
    total_cart = np.nansum(Skf) * dk**2 * df
    total_polar = np.nansum(kcol * Skftheta) * dk * dtheta * df
    Skftheta_norm = Skftheta * total_cart / total_polar
    Fkftheta_norm = kcol**-2 * Skftheta_norm
    Bkftheta_norm = kcol**2 * Skftheta_norm

    # Trim the noise-afflicted frequency range
    inds_f = np.flatnonzero((f > F_LOW_HZ) & (f < F_HIGH_HZ))

    # Locate the dispersion shell in (f, theta) space
    Skftheta_filt = median_filter(Skftheta_norm, size=(3, 1, 1), mode="reflect")
    ind_p = np.argmax(np.nan_to_num(Skftheta_filt, nan=-np.inf), axis=0)
    Smax = np.take_along_axis(Skftheta_filt, ind_p[None, :, :], axis=0)[0]
    kp = k[ind_p]                      # (theta, f)
    Sp = Smax * kp * dk                # slope^2/df/dtheta

    # Convert the shell to (k, theta) space
    f_k_theta = np.full((len(k), len(theta)), np.nan)
    Smax_k_theta = np.full((len(k), len(theta)), np.nan)
    f_keep = f < F_HIGH_HZ
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for i in range(1, len(k)):
            sel = (np.abs(kp - k[i]) < dk) & f_keep[None, :]   # (theta, f)
            fbit = np.where(sel, f[None, :], np.nan)
            f_k_theta[i, :] = np.nanmean(fbit, axis=1)
            Smax_k_theta[i, :] = (np.nansum(np.where(sel, Sp, 0.0), axis=1)
                                  * df / (k[i] * dk))

    # Phase speed difference between a wave and its 180-degree-flipped pair
    half = len(theta) // 2
    f_k_theta_opposed = np.concatenate(
        [f_k_theta[:, half:], f_k_theta[:, :half]], axis=1)
    delta_c = 2 * np.pi * (f_k_theta - f_k_theta_opposed) / k[:, None]

    # Fraction of good phase speed difference sets the cutoff wavenumber
    frac_good = np.sum(~np.isnan(delta_c), axis=1) / len(theta)
    good = np.flatnonzero(frac_good > 0.1)
    if good.size:
        inds_k = np.arange(min(good[-1] + 2, len(k)))
    else:
        warnings.warn("no wavenumber bin passed the phase-speed-difference "
                      "quality check; retaining all wavenumbers")
        inds_k = np.arange(len(k))

    # Integrate with respect to frequency
    S_k_theta = np.nansum(Skftheta_norm, axis=2) * df
    F_k_theta = np.nansum(Fkftheta_norm, axis=2) * df
    B_k_theta = np.nansum(Bkftheta_norm, axis=2) * df

    # Integrate with respect to wavenumber
    S_f_theta = np.nansum(kcol * Skftheta_norm, axis=0) * dk   # (theta, f)
    F_f_theta = np.nansum(kcol * Fkftheta_norm, axis=0) * dk
    B_f_theta = np.nansum(kcol * Bkftheta_norm, axis=0) * dk

    # Integrate with respect to direction
    S_k = np.sum(k[:, None] * S_k_theta, axis=1) * dtheta
    F_k = np.sum(k[:, None] * F_k_theta, axis=1) * dtheta
    B_k = np.sum(k[:, None] * B_k_theta, axis=1) * dtheta
    S_f = np.sum(S_f_theta, axis=0) * dtheta
    F_f = np.sum(F_f_theta, axis=0) * dtheta
    B_f = np.sum(B_f_theta, axis=0) * dtheta

    # Slope and elevation variances
    slope_var = np.nansum(k[:, None] * S_k_theta) * dk * dtheta
    elev_var = np.nansum(k[:, None] * F_k_theta) * dk * dtheta

    # Cut to the resolved wavenumber/frequency ranges
    Skf_cut = Skftheta_norm[np.ix_(inds_k, np.arange(len(theta)), inds_f)]
    Fkf_cut = Fkftheta_norm[np.ix_(inds_k, np.arange(len(theta)), inds_f)]
    Bkf_cut = Bkftheta_norm[np.ix_(inds_k, np.arange(len(theta)), inds_f)]

    # Inverse phase speed spectra
    nu, dnu, nu_lo, nu_hi = _inverse_phase_speed_grid(dnu_base, n_nu_log)
    k_mat = k[inds_k][:, None]
    f_mat = f[inds_f][None, :]
    nu_mat = k_mat / (2 * np.pi * f_mat)

    Qs_nu_theta = np.full((len(theta), len(nu)), np.nan)
    Qeta_nu_theta = np.full_like(Qs_nu_theta, np.nan)
    Qb_nu_theta = np.full_like(Qs_nu_theta, np.nan)
    if n_nu_log is not None:
        # Cloud-in-cell scatter in log(nu): each cell's variance is split
        # linearly between the two nearest bin centers, suppressing the
        # quantization noise of hard binning on sparse cells
        h = np.log(nu[1]) - np.log(nu[0])
        pos = (np.log(np.maximum(nu_mat, 1e-300)) - np.log(nu[0])) / h
        j = np.floor(pos).astype(int)
        frac = pos - j
        cic = [(j, 1.0 - frac), (j + 1, frac)]
    else:
        cic = None
        nu_masks = [(nu_mat >= nu_lo[n]) & (nu_mat < nu_hi[n])
                    for n in range(len(nu))]
    for tn in range(len(theta)):
        # Spectral variance density per direction
        Sslice = k_mat * Skf_cut[:, tn, :] * dk * df    # slope^2/rad
        Fslice = k_mat * Fkf_cut[:, tn, :] * dk * df    # m^2/rad
        Bslice = k_mat * Bkf_cut[:, tn, :] * dk * df    # m^-2/rad
        if cic is not None:
            for arr, out in ((Sslice, Qs_nu_theta), (Fslice, Qeta_nu_theta),
                             (Bslice, Qb_nu_theta)):
                vals = np.nan_to_num(arr, nan=0.0).ravel()
                acc = np.zeros(len(nu))
                for jj, w in cic:
                    idx = np.clip(jj.ravel(), 0, len(nu) - 1)
                    inside = (jj.ravel() >= 0) & (jj.ravel() <= len(nu) - 1)
                    acc += np.bincount(idx, weights=np.where(
                        inside, w.ravel() * vals, 0.0), minlength=len(nu))
                out[tn, :] = acc / (nu * dnu)
        else:
            for n, sel in enumerate(nu_masks):
                norm = nu[n] * dnu[n]
                Qs_nu_theta[tn, n] = np.nansum(
                    np.where(sel, Sslice, 0.0)) / norm
                Qeta_nu_theta[tn, n] = np.nansum(
                    np.where(sel, Fslice, 0.0)) / norm
                Qb_nu_theta[tn, n] = np.nansum(
                    np.where(sel, Bslice, 0.0)) / norm

    Qs = np.sum(nu[None, :] * Qs_nu_theta, axis=0) * dtheta      # s/m
    Qeta = np.sum(nu[None, :] * Qeta_nu_theta, axis=0) * dtheta  # m^3/s
    Qb = np.sum(nu[None, :] * Qb_nu_theta, axis=0) * dtheta      # 1/(s*m)

    coords = {
        "k": ("k", k, {"units": "rad m-1", "long_name": "wavenumber"}),
        "f": ("f", f, {"units": "Hz", "long_name": "frequency"}),
        "nu": ("nu", nu, {"units": "s m-1", "long_name": "inverse phase speed"}),
        "theta": ("theta", theta,
                  {"units": "rad", "long_name": "wave going-to direction"}),
    }
    ds = xr.Dataset(
        data_vars={
            # wavenumber spectra
            "S_k": (("k",), S_k, {"units": "slope^2 / (rad m-1)"}),
            "S_k_theta": (("k", "theta"), S_k_theta,
                          {"units": "slope^2 / (rad m-1)^2 / rad"}),
            "F_k": (("k",), F_k, {"units": "m^2 / (rad m-1)"}),
            "F_k_theta": (("k", "theta"), F_k_theta,
                          {"units": "m^2 / (rad m-1)^2 / rad"}),
            "B_k": (("k",), B_k, {"units": "rad m-1"}),
            "B_k_theta": (("k", "theta"), B_k_theta, {"units": "rad-1"}),
            # frequency spectra
            "S_f": (("f",), S_f, {"units": "slope^2 / Hz"}),
            "S_f_theta": (("theta", "f"), S_f_theta,
                          {"units": "slope^2 / Hz / rad"}),
            "F_f": (("f",), F_f, {"units": "m^2 / Hz"}),
            "F_f_theta": (("theta", "f"), F_f_theta, {"units": "m^2 / Hz / rad"}),
            "B_f": (("f",), B_f, {"units": "m^-2 / Hz"}),
            "B_f_theta": (("theta", "f"), B_f_theta, {"units": "m^-2 / Hz / rad"}),
            # inverse phase speed spectra
            "Qs": (("nu",), Qs, {"units": "m s-1"}),
            "Qs_nu_theta": (("nu", "theta"), Qs_nu_theta.T,
                            {"units": "m^2 s^-2 / rad"}),
            "Qeta": (("nu",), Qeta, {"units": "m^3 s-1"}),
            "Qeta_nu_theta": (("nu", "theta"), Qeta_nu_theta.T,
                              {"units": "m^4 s^-2 / rad"}),
            "Qb": (("nu",), Qb, {"units": "s-1 m-1"}),
            "Qb_nu_theta": (("nu", "theta"), Qb_nu_theta.T,
                            {"units": "s^-2 / rad"}),
            # dispersion shell diagnostics
            "f_down_k_theta": (("k", "theta"), f_k_theta, {"units": "Hz"}),
            "f_up_k_theta": (("k", "theta"), f_k_theta_opposed, {"units": "Hz"}),
            "delta_c_k_theta": (("k", "theta"), delta_c, {"units": "m s-1"}),
            "Sp_k_theta": (("k", "theta"), Smax_k_theta),
            "Sp_f_theta": (("f", "theta"), Sp.T),
            "kp_f_theta": (("f", "theta"), kp.T, {"units": "rad m-1"}),
        },
        coords=coords,
        attrs={
            "k_cutoff": float(k[inds_k].max()),
            "total_slope_variance": float(slope_var),
            "total_elevation_variance": float(elev_var),
            "dk": float(dk),
            "df": float(df),
            "dtheta": float(dtheta),
            "heading_deg": float(heading_deg),
            "n_freq_avg": int(n_freq_avg),
            "dnu_base": float(dnu_base),
            "n_nu_log": 0 if n_nu_log is None else int(n_nu_log),
        },
    )
    return ds
