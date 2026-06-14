"""Gram-Charlier expansion fits to wave slope joint PDFs (Cox & Munk [1954])."""

import numpy as np
from scipy.optimize import minimize

__all__ = ["gram_charlier_pdf", "fit_gram_charlier_slope_pdf"]


def gram_charlier_pdf(xi, zeta, mss_up, mss_cross,
                      c21=0.0, c03=0.0, c40=0.0, c04=0.0, c22=0.0):
    """Gram-Charlier expansion about a 2D Gaussian slope distribution.

    xi, zeta are crosswind and upwind slopes normalized by the respective
    slope standard deviations.
    """
    coeff = 1.0 / (2 * np.pi * np.sqrt(mss_up) * np.sqrt(mss_cross))
    return coeff * np.exp(-(xi**2 + zeta**2) / 2) * (
        1
        - 0.5 * c21 * (xi**2 - 1) * zeta
        - (1 / 6) * c03 * (zeta**3 - 3 * zeta)
        + (1 / 24) * c40 * (xi**4 - 6 * xi**2 + 3)
        + (1 / 24) * c04 * (zeta**4 - 6 * zeta**2 + 3)
        + (1 / 4) * c22 * (xi**2 - 1) * (zeta**2 - 1)
    )


def fit_gram_charlier_slope_pdf(slope_centers, P_slope_c_u, mss_u, mss_c):
    """Least-squares Gram-Charlier fit to a slope joint PDF.

    Args:
        slope_centers : 1D slope bin centers shared by both axes
        P_slope_c_u   : joint PDF, crosswind on axis 0, upwind on axis 1
        mss_u, mss_c  : upwind and crosswind mean square slope

    Returns:
        dict with the fitted PDF (P_fit), skewness coefficients (c21, c03),
        kurtosis coefficients (c40, c04, c22), R_squared and RMSE.
    """
    slope_centers = np.asarray(slope_centers, dtype=float)
    P_slope_c_u = np.asarray(P_slope_c_u, dtype=float)
    xi, zeta = np.meshgrid(slope_centers / np.sqrt(mss_c),
                           slope_centers / np.sqrt(mss_u), indexing="ij")

    def model(b):
        return gram_charlier_pdf(xi, zeta, mss_u, mss_c, *b)

    def cost(b):
        return np.sum((model(b) - P_slope_c_u) ** 2)

    result = minimize(cost, np.zeros(5), method="Nelder-Mead",
                      options={"xatol": 1e-8, "fatol": 1e-12, "maxiter": 5000})

    P_fit = model(result.x)
    residuals = P_fit - P_slope_c_u
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((P_slope_c_u - P_slope_c_u.mean())**2)

    return {
        "P_fit": P_fit,
        "c21": result.x[0],
        "c03": result.x[1],
        "c40": result.x[2],
        "c04": result.x[3],
        "c22": result.x[4],
        "R_squared": 1 - ss_res / ss_tot,
        "RMSE": np.sqrt(ss_res / P_slope_c_u.size),
    }
