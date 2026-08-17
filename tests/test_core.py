import numpy as np
import pytest

from pcca.core import cca, loo_projection, orthonormalize, pca, pca_from_gram


def test_pca_returns_orthonormal_columns():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(50, 200))
    u, _ = pca(X, dim=5)
    assert u.shape == (50, 5)
    assert np.allclose(u.T @ u, np.eye(5), atol=1e-10)


def test_pca_recovers_a_planted_low_rank_structure():
    """Two well-separated clusters should be visible in the first component."""
    rng = np.random.default_rng(1)
    n = 60
    signal = np.concatenate([np.ones(n // 2), -np.ones(n // 2)])
    loadings = rng.normal(size=100)
    X = np.outer(signal, loadings) + 0.01 * rng.normal(size=(n, 100))
    u, _ = pca(X, dim=1)
    assert abs(np.corrcoef(u[:, 0], signal)[0, 1]) > 0.99


def test_pca_from_gram_matches_pca():
    """The Gram-matrix route is exact, which is what makes leave-one-out cheap."""
    rng = np.random.default_rng(2)
    X = rng.normal(size=(40, 300))
    u_direct, s_direct = pca(X, dim=6)
    u_gram, s_gram = pca_from_gram(X @ X.T, dim=6)
    assert np.allclose(np.abs(u_direct), np.abs(u_gram), atol=1e-8)
    assert np.allclose(s_direct, s_gram, atol=1e-8)


def test_pca_from_gram_on_a_submatrix_equals_leaving_that_sample_out():
    """This equivalence is the entire justification for the LOO shortcut."""
    rng = np.random.default_rng(3)
    X = rng.normal(size=(30, 200))
    gram = X @ X.T
    hold_out = 7
    idx = np.array([i for i in range(30) if i != hold_out])
    u_direct, _ = pca(X[idx], dim=4)
    u_sub, _ = pca_from_gram(gram[np.ix_(idx, idx)], dim=4)
    assert np.allclose(np.abs(u_direct), np.abs(u_sub), atol=1e-8)


def test_orthonormalize_fixes_non_unit_columns():
    """plink2 .eigenvec columns are orthogonal but not unit-norm."""
    rng = np.random.default_rng(4)
    q, _ = np.linalg.qr(rng.normal(size=(20, 3)))
    scaled = q * np.array([5.0, 0.2, 13.0])
    assert not np.allclose(scaled.T @ scaled, np.eye(3))
    out = orthonormalize(scaled)
    assert np.allclose(out.T @ out, np.eye(3), atol=1e-10)


def test_cca_rejects_non_orthonormal_input():
    """The guard is deliberate: silent wrong answers are the failure mode."""
    rng = np.random.default_rng(11)
    q, _ = np.linalg.qr(rng.normal(size=(20, 3)))
    with pytest.raises(ValueError, match="orthonormal"):
        cca(q, q * 5.0, dim=2)


def test_cca_of_identical_subspaces_gives_unit_correlations():
    rng = np.random.default_rng(5)
    u, _ = np.linalg.qr(rng.normal(size=(40, 5)))
    _, _, rho = cca(u, u, dim=3)
    assert np.allclose(rho, 1.0, atol=1e-10)


def test_cca_of_orthogonal_subspaces_gives_zero_correlations():
    rng = np.random.default_rng(6)
    q, _ = np.linalg.qr(rng.normal(size=(40, 10)))
    _, _, rho = cca(q[:, :5], q[:, 5:], dim=3)
    assert np.allclose(rho, 0.0, atol=1e-10)


def test_cca_recovers_a_shared_latent_variable():
    """The point of PCCA: a latent shared by both views appears in both projections."""
    rng = np.random.default_rng(7)
    n = 80
    latent = rng.normal(size=n)
    x_view = np.column_stack([latent + 0.05 * rng.normal(size=n), rng.normal(size=n)])
    y_view = np.column_stack([latent + 0.05 * rng.normal(size=n), rng.normal(size=n)])
    u_x, _ = pca(x_view, dim=2)
    u_y, _ = pca(y_view, dim=2)
    l_x, l_y, rho = cca(u_x, u_y, dim=1)
    assert rho[0] > 0.95
    assert abs(np.corrcoef(l_x[:, 0], latent)[0, 1]) > 0.9


def test_cca_correlations_are_descending_and_bounded():
    rng = np.random.default_rng(8)
    u_x, _ = np.linalg.qr(rng.normal(size=(50, 6)))
    u_y, _ = np.linalg.qr(rng.normal(size=(50, 4)))
    _, _, rho = cca(u_x, u_y, dim=4)
    assert np.all(rho <= 1.0 + 1e-12) and np.all(rho >= -1e-12)
    assert np.all(np.diff(rho) <= 1e-12)


def test_loo_projection_shape_and_finiteness():
    rng = np.random.default_rng(9)
    n = 25
    latent = rng.normal(size=n)
    X = np.outer(latent, rng.normal(size=40)) + 0.1 * rng.normal(size=(n, 40))
    Y = np.outer(latent, rng.normal(size=60)) + 0.1 * rng.normal(size=(n, 60))
    proj = loo_projection(X, Y @ Y.T, dim_x=5, dim_y=3, dim_z=2)
    assert proj.shape == (n, 2)
    assert np.all(np.isfinite(proj))


def test_loo_projection_recovers_the_shared_latent_for_held_out_samples():
    """The point of Figure 2B: a held-out individual still lands by its group.

    This is the test that catches an unwhitened projection — dropping the 1/s**2
    term leaves the output dominated by the leading component and destroys the
    correlation with the latent.
    """
    rng = np.random.default_rng(21)
    n = 60
    group = np.concatenate([np.ones(n // 2), -np.ones(n // 2)])
    x_load, y_load = rng.normal(size=120), rng.normal(size=90)
    X = np.outer(group, x_load) + 1.5 * rng.normal(size=(n, 120))
    Y = np.outer(group, y_load) + 0.3 * rng.normal(size=(n, 90))
    X -= X.mean(axis=0)
    Y -= Y.mean(axis=0)

    proj = loo_projection(X, Y @ Y.T, dim_x=10, dim_y=3, dim_z=2)
    assert abs(np.corrcoef(proj[:, 0], group)[0, 1]) > 0.7


def test_loo_projection_signs_are_aligned_across_holdouts():
    """Unaligned holdouts scatter a group to both extremes of an axis.

    Each holdout's CCA returns singular vectors up to sign, so without alignment
    the held-out coordinates flip polarity at random and the group means collapse
    toward zero despite each individual projection being correct.
    """
    rng = np.random.default_rng(31)
    n = 50
    group = np.concatenate([np.ones(n // 2), -np.ones(n // 2)])
    X = np.outer(group, rng.normal(size=100)) + 1.2 * rng.normal(size=(n, 100))
    Y = np.outer(group, rng.normal(size=80)) + 0.3 * rng.normal(size=(n, 80))
    X -= X.mean(axis=0)
    Y -= Y.mean(axis=0)

    proj = loo_projection(X, Y @ Y.T, dim_x=8, dim_y=3, dim_z=2)
    lo = proj[group == 1, 0].mean()
    hi = proj[group == -1, 0].mean()
    # With aligned signs the two group means sit on opposite sides of the axis.
    assert abs(lo - hi) > 0.5 * np.abs(proj[:, 0]).mean()
