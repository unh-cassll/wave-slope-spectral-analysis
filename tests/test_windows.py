import numpy as np

from slopespectra import blackman_harris, circular_tukey, taper_2d


def test_blackman_harris_endpoints_and_peak():
    w = blackman_harris(65)
    assert np.isclose(w[32], 1.0, atol=1e-3)
    assert abs(w[0]) < 1e-4
    assert abs(w[-1]) < 1e-4


def test_circular_tukey_geometry():
    arr = np.ones((64, 64))
    w = circular_tukey(arr, 0.2)
    # Corners are outside the inscribed circle
    assert w[0, 0] == 0.0
    assert w[-1, -1] == 0.0
    # Flat (scaled) interior
    interior = w[28:36, 28:36]
    assert np.allclose(interior, interior[0, 0])
    assert interior[0, 0] > 1.0  # spectral-norm rescaling exceeds unity


def test_circular_tukey_3d_repeats_window():
    arr = np.ones((32, 32, 5))
    w = circular_tukey(arr, 0.3)
    for i in range(1, 5):
        assert np.allclose(w[:, :, i], w[:, :, 0])


def test_taper_2d_cosine_zeros_borders():
    X = np.ones((40, 50))
    A, window = taper_2d(X, kind="cosine", w=0.25)
    assert A.shape == X.shape
    assert np.allclose(A[0, :], 0.0)
    assert np.allclose(A[:, 0], 0.0)
    assert np.isclose(window[20, 25], 1.0)


def test_taper_2d_3d_input():
    X = np.ones((30, 30, 3))
    A, _ = taper_2d(X, kind="cosine", w=0.2)
    assert A.shape == X.shape
    assert np.allclose(A[:, :, 0], A[:, :, 2])


def test_circular_tukey_power_normalization_unbiased():
    rng = np.random.default_rng(2)
    field = rng.standard_normal((128, 128, 8))
    w = circular_tukey(field, 0.2, normalization="power")
    assert np.isclose(w.var(), field.var(), rtol=0.05)


def test_circular_tukey_power_window_unit_power():
    w = circular_tukey(np.ones((96, 96, 12)), 0.25, normalization="power",
                       temporal_alpha=0.2)
    assert np.isclose(np.mean(w**2), 1.0, rtol=1e-10)


def test_circular_tukey_spectral_parity():
    # 'spectral' reproduces the MATLAB-library scaling
    arr = np.ones((64, 64))
    w = circular_tukey(arr, 0.2, normalization="spectral")
    w2d = circular_tukey(arr, 0.2, normalization="none")
    C = np.sqrt(64 * 64) / np.linalg.norm(w2d, 2)
    assert np.allclose(w, C * w2d)


def test_circular_tukey_temporal_taper():
    arr = np.ones((32, 32, 40))
    w = circular_tukey(arr, 0.2, normalization="none", temporal_alpha=0.5)
    # Temporal ends tapered, midpoint untapered
    assert w[16, 16, 0] == 0.0
    assert w[16, 16, 20] == 1.0
