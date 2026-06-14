"""Figure styling and standard plots for wave slope spectral products.

Styling follows the E-PSS paper figure_style (seaborn ticks theme, Fira
Sans, fixed color cycle).
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

__all__ = ["figure_style", "plot_sub_spectra", "plot_kf_slice",
           "plot_omnidirectional"]

GRAV = 9.81          # m/s^2
SIGMA_W = 0.072      # N/m, clean-water surface tension
RHO_W = 1000.0       # kg/m^3


def figure_style(title_fontsize=10, label_fontsize=10, tick_fontsize=10):
    """E-PSS paper plot styling; returns (color_list, fullwidth, fullheight,
    fsize) with full-page dimensions for letter paper, 0.5 inch margins."""
    fsize = 10
    lw = 1.0

    sns.set_theme(style="ticks", palette="deep", font="Fira Sans")

    color_list = ['#4C2882', '#367588', '#A52A2A', '#C39953', '#2A52BE',
                  '#006611']
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)

    plt.rcParams.update({
        'axes.grid': True,
        'font.size': fsize,
        'axes.titlesize': title_fontsize,
        'axes.labelsize': label_fontsize,
        'xtick.labelsize': tick_fontsize,
        'ytick.labelsize': tick_fontsize,
        'legend.fontsize': label_fontsize,
        'grid.linewidth': lw,
        'xtick.major.width': lw,
        'ytick.major.width': lw,
    })

    fullwidth = 7.5
    fullheight = 10

    return color_list, fullwidth, fullheight, fsize


def _gravity_capillary_f(k):
    """Gravity-capillary linear dispersion frequency [Hz]."""
    return np.sqrt(GRAV * k + (SIGMA_W / RHO_W) * k**3) / (2 * np.pi)


def plot_sub_spectra(ds, fig=None):
    """Overview figure of the compute_sub_spectra products.

    Panels: (a) omnidirectional wavenumber spectra, (b) omnidirectional
    frequency spectra, (c) directional wavenumber slope spectrum,
    (d) directional frequency slope spectrum, (e) inverse phase speed
    spectra, (f) dispersion-shell phase speed difference.

    Args:
        ds : Dataset from compute_sub_spectra

    Returns:
        matplotlib Figure.
    """
    colors, fullwidth, fullheight, _ = figure_style()
    if fig is None:
        fig = plt.figure(figsize=(fullwidth, fullheight * 0.85))
    axs = fig.subplots(3, 2)

    k = ds["k"].values
    f = ds["f"].values
    nu = ds["nu"].values
    theta_deg = np.rad2deg(ds["theta"].values)

    ax = axs[0, 0]
    ax.loglog(k, ds["S_k"], label=r"$S(k)$ [slope$^2$/(rad m$^{-1}$)]")
    ax.loglog(k, ds["F_k"], label=r"$F(k)$ [m$^2$/(rad m$^{-1}$)]")
    ax.loglog(k, ds["B_k"], label=r"$B(k)$ [rad m$^{-1}$]")
    ax.axvline(ds.attrs.get("k_cutoff", np.nan), color="k", linestyle=":",
               linewidth=1)
    ax.set_xlabel(r"$k$ [rad m$^{-1}$]")
    ax.set_ylabel("spectral density")
    ax.legend(loc="lower left")
    ax.set_title("(a) wavenumber spectra")

    ax = axs[0, 1]
    ax.loglog(f, ds["S_f"], label=r"$S(f)$ [slope$^2$ Hz$^{-1}$]")
    ax.loglog(f, ds["F_f"], label=r"$F(f)$ [m$^2$ Hz$^{-1}$]")
    ax.set_xlabel(r"$f$ [Hz]")
    ax.set_ylabel("spectral density")
    ax.legend(loc="lower left")
    ax.set_title("(b) frequency spectra")

    ax = axs[1, 0]
    Sk_theta = ds["S_k_theta"].values
    pc = ax.pcolormesh(theta_deg, k,
                       np.log10(np.where(Sk_theta > 0, Sk_theta, np.nan)),
                       shading="nearest", cmap="viridis")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\theta$ [$^\circ$]")
    ax.set_ylabel(r"$k$ [rad m$^{-1}$]")
    fig.colorbar(pc, ax=ax, label=r"log$_{10} S(k,\theta)$")
    ax.set_title("(c) directional wavenumber spectrum")

    ax = axs[1, 1]
    Sf_theta = ds["S_f_theta"].values.T   # (f, theta)
    pc = ax.pcolormesh(theta_deg, f,
                       np.log10(np.where(Sf_theta > 0, Sf_theta, np.nan)),
                       shading="nearest", cmap="viridis")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\theta$ [$^\circ$]")
    ax.set_ylabel(r"$f$ [Hz]")
    fig.colorbar(pc, ax=ax, label=r"log$_{10} S(f,\theta)$")
    ax.set_title("(d) directional frequency spectrum")

    ax = axs[2, 0]
    ax.loglog(nu, ds["Qs"], label=r"$Q_s(\nu)$ [s m$^{-1}$]")
    ax.loglog(nu, ds["Qeta"], label=r"$Q_\eta(\nu)$ [m$^3$ s$^{-1}$]")
    ax.set_xlabel(r"$\nu$ [s m$^{-1}$]")
    ax.set_ylabel("spectral density")
    ax.legend(loc="lower left")
    ax.set_title("(e) inverse phase speed spectra")

    ax = axs[2, 1]
    dc = ds["delta_c_k_theta"].values
    pc = ax.pcolormesh(theta_deg, k, dc, shading="nearest", cmap="RdBu_r",
                       vmin=-np.nanmax(np.abs(dc)) if np.isfinite(dc).any()
                       else -1,
                       vmax=np.nanmax(np.abs(dc)) if np.isfinite(dc).any()
                       else 1)
    ax.set_yscale("log")
    ax.set_xlabel(r"$\theta$ [$^\circ$]")
    ax.set_ylabel(r"$k$ [rad m$^{-1}$]")
    fig.colorbar(pc, ax=ax, label=r"$\Delta c$ [m s$^{-1}$]")
    ax.set_title("(f) phase speed difference")

    fig.tight_layout()
    return fig


def plot_kf_slice(dirspect, axis="kx", klim=None, ax=None):
    """Wavenumber-frequency slice through the directional spectrum, with
    the gravity-capillary dispersion relation overlaid.

    Args:
        dirspect : DirectionalSpectrum from compute_slope_spectrum
        axis     : 'kx' for the cross-slope slice (along the kx axis) or
                   'ky' for the along-slope slice
        klim     : optional symmetric wavenumber axis limit [rad/m]

    Returns:
        matplotlib Axes.
    """
    figure_style()
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4.5))

    Skf = dirspect.Skf
    m = Skf.shape[0]
    center = m // 2
    dk = dirspect.dk
    kvec = dk * (np.arange(m) - center)
    f = dirspect.f

    if axis == "kx":
        spect_slice = Skf[center, :, :].T
    elif axis == "ky":
        spect_slice = Skf[:, center, :].T
        kvec = -kvec    # rows carry descending ky labels
    else:
        raise ValueError("axis must be 'kx' or 'ky'")

    pc = ax.pcolormesh(kvec, f,
                       np.log10(np.where(spect_slice > 0, spect_slice,
                                         np.nan)),
                       shading="nearest", cmap="viridis")
    kpos = dk * np.arange(1, m // 2 + 1)
    fd = _gravity_capillary_f(kpos)
    ax.plot(kpos, fd, "-r", linewidth=2)
    ax.plot(-kpos, fd, "-r", linewidth=2)
    if klim is not None:
        ax.set_xlim(-klim, klim)
    ax.set_xlabel(rf"$k_{axis[1]}$ [rad m$^{{-1}}$]")
    ax.set_ylabel(r"$f$ [Hz]")
    plt.colorbar(pc, ax=ax, label=r"log$_{10} S(k,f)$")
    return ax


def plot_omnidirectional(ds, ax=None, compare_dispersion=True, min_k_bins=0,
                         mark_fov_scale=True):
    """Omnidirectional frequency elevation spectrum F(f) compared against
    F(f) mapped from F(k) through deep-water linear dispersion.

    Args:
        ds          : Dataset from compute_sub_spectra
        min_k_bins  : lowest wavenumber bins excluded from the mapped
                      curve (waves longer than the FOV pile their energy
                      into these bins, inflating the mapped spectrum there)
        mark_fov_scale : mark f corresponding to the lowest resolved k

    Returns:
        matplotlib Axes.
    """
    figure_style()
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4.5))

    f = ds["f"].values
    ax.loglog(f, ds["F_f"], "-", linewidth=2, label=r"$F(f)$ direct")
    if compare_dispersion:
        # F(f) = F(k) dk/df with dk/df = 2*pi/cg (deep water)
        k = ds["k"].values[min_k_bins:]
        f_disp = np.sqrt(GRAV * k) / (2 * np.pi)
        cg = 0.5 * np.sqrt(GRAV / k)
        F_f_disp = ds["F_k"].values[min_k_bins:] * 2 * np.pi / cg
        ax.loglog(f_disp, F_f_disp, "-", linewidth=2,
                  label=r"$F(f)$ from $F(k)$, linear dispersion")
        if mark_fov_scale:
            ax.axvline(f_disp[0], color="k", linestyle=":", linewidth=1,
                       label=r"$f(k_{min})$, FOV scale")
    ax.set_xlabel(r"$f$ [Hz]")
    ax.set_ylabel(r"$F(f)$ [m$^2$ Hz$^{-1}$]")
    ax.legend(loc="lower left")
    return ax
