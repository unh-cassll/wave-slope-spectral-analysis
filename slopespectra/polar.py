"""Cartesian-to-polar conversion of directional wavenumber spectra."""

import numpy as np
from scipy.interpolate import PchipInterpolator

__all__ = ["polar_from_cartesian", "PolarBinner"]


class PolarBinner:
    """Reusable Cartesian (kx, ky) -> polar (k, theta) projection.

    Bin geometry depends only on the wavenumber grid, so a single instance
    is reused across many spectral slices (e.g. every frequency plane of a
    wavenumber-frequency spectrum).

    Direction convention: theta is the compass-style going-to direction
    (0 = +ky, 90 deg = +kx), offset by `heading_deg`, in radians on [0, 2*pi).

    Each ring's energy is set by its azimuthal-mean density times the
    analytic ring area 2*pi*k*dk, not by the discrete cell count, so the
    bin-to-bin sawtooth from the fluctuating cells-per-annulus (the Gauss
    circle problem) does not appear in the radial spectrum.
    """

    def __init__(self, kx, ky, dtheta=np.deg2rad(5.0), heading_deg=0.0):
        kx = np.asarray(kx, dtype=float)
        ky = np.asarray(ky, dtype=float)
        Nx, Ny = kx.shape
        self.num_bins = min(Nx, Ny) // 2
        self.kmin = float(np.mean(np.abs(np.diff(kx[0, :]))))
        self.dtheta = float(dtheta)
        kmax = self.num_bins * self.kmin
        self.k_vec = np.linspace(self.kmin, kmax, self.num_bins)

        # Compass angle of each (kx, ky) cell, plus camera heading
        theta_mat = np.mod(np.arctan2(kx, ky) + np.deg2rad(heading_deg),
                           2 * np.pi)

        dtheta_deg = np.rad2deg(self.dtheta)
        self.theta_vec = np.deg2rad(np.arange(0.0, 360.0, dtheta_deg))

        # Per-bin flat indices: full set (variance) and theta-unique sorted
        # set (interpolation support)
        kmat = np.sqrt(kx**2 + ky**2).ravel()
        theta_flat = theta_mat.ravel()
        self._bin_inds = []
        self._bin_unique_inds = []
        self._bin_unique_theta = []
        for n in range(1, self.num_bins + 1):
            k_val = n * self.kmin
            inds = np.flatnonzero((kmat > k_val - self.kmin * 0.5)
                                  & (kmat <= k_val + self.kmin * 0.5))
            order = np.argsort(theta_flat[inds], kind="stable")
            sorted_inds = inds[order]
            theta_sorted = theta_flat[sorted_inds]
            _, first = np.unique(theta_sorted, return_index=True)
            self._bin_inds.append(inds)
            self._bin_unique_inds.append(sorted_inds[first])
            self._bin_unique_theta.append(theta_sorted[first])

    def __call__(self, Skxy):
        """Project one Cartesian spectral slice onto the (k, theta) grid."""
        Skxy = np.asarray(Skxy, dtype=float)
        flat = Skxy.ravel()
        kmin, dtheta = self.kmin, self.dtheta
        out = np.full((self.num_bins, len(self.theta_vec)), np.nan)

        for n in range(self.num_bins):
            k_val = (n + 1) * kmin
            vals_full = flat[self._bin_inds[n]]
            finite = np.isfinite(vals_full)
            cnt = int(finite.sum())
            if cnt == 0:
                continue
            # Ring total variance = azimuthal-mean density x analytic ring
            # area (2*pi*k*dk). Using the analytic area instead of the
            # fluctuating cell count removes the bin-to-bin sawtooth (the
            # number of grid cells per annulus is noisy; the Gauss circle
            # problem).
            mean_density = vals_full[finite].sum() / cnt
            Svar = mean_density * 2 * np.pi * k_val * kmin
            theta_u = self._bin_unique_theta[n]
            vals_u = flat[self._bin_unique_inds[n]]
            good = np.isfinite(vals_u)
            theta_u = theta_u[good]
            vals_u = vals_u[good]
            if len(theta_u) < 2:
                continue
            # Periodic pchip support: triple the circle, evaluate the center copy
            theta_triple = np.concatenate((theta_u, theta_u + 2 * np.pi,
                                           theta_u + 4 * np.pi))
            vals_triple = np.concatenate((vals_u, vals_u, vals_u))
            interp = PchipInterpolator(theta_triple, vals_triple)(
                self.theta_vec + 2 * np.pi)
            # Conserve variance within the wavenumber bin
            denom = k_val * interp.sum() * kmin * dtheta
            if denom != 0:
                out[n, :] = interp * (Svar / denom)

        # Enforce global variance conservation
        var_cart = np.nansum(flat) * kmin * kmin
        var_polar = np.nansum(self.k_vec[:, None] * out) * kmin * dtheta
        if var_polar != 0 and np.isfinite(var_polar):
            out *= var_cart / var_polar
        return out


def polar_from_cartesian(Skxy, kx, ky, dtheta=np.deg2rad(5.0), heading_deg=0.0):
    """Convert a Cartesian directional wavenumber spectrum S(kx, ky) to
    polar coordinates S(k, theta), conserving variance per bin and globally.

    Args:
        Skxy        : (Ny, Nx) spectral density, NaN beyond the Nyquist circle
        kx, ky      : (Ny, Nx) wavenumber label arrays (kx along columns,
                      ky descending along rows)
        dtheta      : direction bin width [rad]
        heading_deg : camera heading added to the direction labels [deg]

    Returns:
        (Sktheta, k_vec, theta_vec)
    """
    binner = PolarBinner(kx, ky, dtheta=dtheta, heading_deg=heading_deg)
    return binner(Skxy), binner.k_vec, binner.theta_vec
