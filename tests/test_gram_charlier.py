import numpy as np
import pytest

from slopespectra import gram_charlier_pdf, fit_gram_charlier_slope_pdf


def test_fit_recovers_known_coefficients():
    mss_u, mss_c = 0.02, 0.012
    truth = {"c21": -0.05, "c03": -0.2, "c40": 0.3, "c04": 0.4, "c22": 0.12}
    centers = np.linspace(-0.6, 0.6, 121)
    xi, zeta = np.meshgrid(centers / np.sqrt(mss_c), centers / np.sqrt(mss_u),
                           indexing="ij")
    P = gram_charlier_pdf(xi, zeta, mss_u, mss_c, **truth)
    out = fit_gram_charlier_slope_pdf(centers, P, mss_u, mss_c)
    for name, val in truth.items():
        assert out[name] == pytest.approx(val, abs=0.02), name
    assert out["R_squared"] > 0.999


def test_zero_coefficients_give_gaussian():
    mss_u = mss_c = 0.01
    centers = np.linspace(-0.4, 0.4, 81)
    xi, zeta = np.meshgrid(centers / np.sqrt(mss_c), centers / np.sqrt(mss_u),
                           indexing="ij")
    P = gram_charlier_pdf(xi, zeta, mss_u, mss_c)
    gauss = np.exp(-(xi**2 + zeta**2) / 2) / (2 * np.pi * mss_u**0.5
                                              * mss_c**0.5)
    assert np.allclose(P, gauss)
