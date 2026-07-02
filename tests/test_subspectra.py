import numpy as np
import pytest

from slopespectra import compute_slope_spectrum, compute_sub_spectra


# Dominant component of the broadband fixture, in wavenumber bins
N_BIN_DOM = 8.4


@pytest.fixture(scope="module")
def subspectra(request):
    # Module-scoped broadband case: build the stack once
    from conftest import broadband_wave_stack
    sx, sy, meta = broadband_wave_stack()
    spec = compute_slope_spectrum(sx, sy, meta["dx"], meta["fs"])
    ds = compute_sub_spectra(spec.Skf, spec.dk, spec.df)
    return ds, spec, (sx, sy, meta)


def test_total_variance_matches_field(subspectra):
    ds, spec, (sx, sy, meta) = subspectra
    assert np.isclose(ds.attrs["total_slope_variance"],
                      sx.var() + sy.var(), rtol=0.02)


def test_directional_peak(subspectra):
    ds, _, _ = subspectra
    D = (ds["S_k_theta"] * ds["k"]).sum("k")
    az = np.rad2deg(float(ds["theta"][int(np.argmax(D.values))]))
    # Energy-weighted construction is centered near 90 deg (east)
    assert abs((az - 90 + 180) % 360 - 180) <= 15


def test_omni_peaks(subspectra):
    ds, spec, (_, _, meta) = subspectra
    k_pk = float(ds["k"][int(np.argmax(ds["S_k"].values))])
    # Dominant component is the 8th bin
    assert np.isclose(k_pk, N_BIN_DOM * spec.dk, atol=2 * spec.dk)
    f_pk = float(ds["f"][int(np.argmax(ds["S_f"].values))])
    f0_dom = np.sqrt(9.81 * N_BIN_DOM * spec.dk) / (2 * np.pi)
    assert np.isclose(f_pk, f0_dom, atol=2 * spec.df)


def test_elevation_saturation_relations(subspectra):
    ds, _, _ = subspectra
    k = ds["k"].values[:, None]
    assert np.allclose(ds["F_k_theta"].values,
                       k**-2 * ds["S_k_theta"].values, equal_nan=True)
    assert np.allclose(ds["B_k_theta"].values,
                       k**2 * ds["S_k_theta"].values, equal_nan=True)


def test_dispersion_shell_frequency(subspectra):
    ds, spec, _ = subspectra
    # Near the dominant wavenumber/direction, the shell frequency should
    # sit near linear dispersion (the shell is sparse for a synthetic
    # multi-component spectrum, so search a small neighborhood)
    k = ds["k"].values
    theta_deg = np.rad2deg(ds["theta"].values)
    i_k = np.argmin(np.abs(k - N_BIN_DOM * spec.dk))
    sel_k = slice(max(i_k - 2, 0), i_k + 3)
    sel_th = (theta_deg >= 75) & (theta_deg <= 105)
    f_shell = ds["f_down_k_theta"].values[sel_k][:, sel_th]
    f_lin = np.sqrt(9.81 * N_BIN_DOM * spec.dk) / (2 * np.pi)
    assert np.isfinite(f_shell).any()
    assert abs(np.nanmedian(f_shell) - f_lin) < 3 * spec.df


# Opposing directionally-spread wave groups at a common wavenumber: the
# phase-speed-difference cutoff requires energy on both sides of the k plane
OPPOSED_COMPONENTS = [(75.0, 0.04, 0.2), (90.0, 0.05, 1.1),
                      (105.0, 0.04, 2.3), (255.0, 0.03, 3.1),
                      (270.0, 0.035, 4.2), (285.0, 0.03, 5.3)]


def test_inverse_phase_speed_spectrum():
    from conftest import synthetic_wave_stack
    sx = sy = None
    for az, amp, ph in OPPOSED_COMPONENTS:
        sx_i, sy_i, meta = synthetic_wave_stack(n_bin=N_BIN_DOM,
                                                direction_deg=az,
                                                amplitude=amp, phase=ph)
        sx = sx_i if sx is None else sx + sx_i
        sy = sy_i if sy is None else sy + sy_i
    spec = compute_slope_spectrum(sx, sy, meta["dx"], meta["fs"])
    ds = compute_sub_spectra(spec.Skf, spec.dk, spec.df)
    # Dominant wave: c = 2*pi*f0/k0 -> nu = 1/c
    k0 = N_BIN_DOM * spec.dk
    f0 = np.sqrt(9.81 * k0) / (2 * np.pi)
    nu0 = k0 / (2 * np.pi * f0)
    assert ds.attrs["k_cutoff"] >= k0
    Qs = ds["Qs"].values
    nu_pk = float(ds["nu"][int(np.nanargmax(Qs))])
    assert abs(nu_pk - nu0) < 0.1
    assert np.nansum(Qs) > 0


def test_k_cutoff_within_grid(subspectra):
    ds, _, _ = subspectra
    assert 0 < ds.attrs["k_cutoff"] <= float(ds["k"].max())

def test_frequency_band_averaging(subspectra):
    ds, spec, (sx, sy, _) = subspectra
    ds5 = compute_sub_spectra(spec.Skf, spec.dk, spec.df, n_freq_avg=5)
    assert ds5.sizes["f"] == spec.Skf.shape[2] // 5
    assert ds5.attrs["df"] == pytest.approx(5 * spec.df)
    # Band centers stay within the original frequency range
    assert float(ds5["f"].min()) > 0
    assert float(ds5["f"].max()) <= spec.Skf.shape[2] * spec.df
    # Total variance is preserved up to the trimmed remainder bins
    assert ds5.attrs["total_slope_variance"] == pytest.approx(
        ds.attrs["total_slope_variance"], rel=1e-3)
    assert np.nansum(ds5["S_f"].values) * ds5.attrs["df"] == pytest.approx(
        np.nansum(ds["S_f"].values) * ds.attrs["df"], rel=1e-3)


def test_nu_grid_coarsening(subspectra):
    ds, spec, _ = subspectra
    ds2 = compute_sub_spectra(spec.Skf, spec.dk, spec.df, dnu_base=0.02)
    assert ds2.sizes["nu"] < ds.sizes["nu"]
    # Integrated inverse-phase-speed energy is grid-independent
    dnu = np.concatenate([[0.02], np.diff(ds2["nu"].values)])
    dnu0 = np.concatenate([[0.01], np.diff(ds["nu"].values)])
    tot2 = np.nansum(ds2["Qs"].values * dnu)
    tot0 = np.nansum(ds["Qs"].values * dnu0)
    assert tot2 == pytest.approx(tot0, rel=0.05)


def test_log_nu_grid(subspectra):
    ds, spec, _ = subspectra
    ds_log = compute_sub_spectra(spec.Skf, spec.dk, spec.df, n_nu_log=48)
    assert ds_log.sizes["nu"] == 48
    nu = ds_log["nu"].values
    # Contiguous log-spaced bins
    assert np.allclose(np.diff(np.log(nu)), np.diff(np.log(nu))[0])
    # Integrated inverse-phase-speed energy matches the default grid
    edges = np.geomspace(0.05, 4.0, 49)
    tot_log = np.nansum(ds_log["Qs"].values * np.diff(edges))
    dnu0 = np.concatenate([[0.01], np.diff(ds["nu"].values)])
    tot0 = np.nansum(ds["Qs"].values * dnu0)
    assert tot_log == pytest.approx(tot0, rel=0.05)
