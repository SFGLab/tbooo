"""PCA and CCA primitives, ported from pachterlab/PCCA `pcca/util.py`.

PCCA is PCA on each view followed by CCA between the two sets of principal
components. The paper uses 30 expression PCs, 5 genotype PCs, and 2 canonical
components.
"""

from __future__ import annotations

import numpy as np


def pca(arr: np.ndarray, dim: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Left singular vectors of a samples x features matrix.

    Computed via the sample-by-sample Gram matrix, which is cheap when features
    vastly outnumber samples (14k genes or 6M SNPs against a few hundred people).
    Equivalent to PCA on a centered and scaled matrix, which is what the callers
    supply.
    """
    if dim is None:
        dim = arr.shape[0]
    return pca_from_gram(arr @ arr.T, dim)


def pca_from_gram(gram: np.ndarray, dim: int) -> tuple[np.ndarray, np.ndarray]:
    """PCA directly from a precomputed Gram matrix.

    Leave-one-out PCA is then just the same call on a submatrix, which turns the
    cross-validation from hours of 6M-SNP decompositions into seconds.
    """
    u, eigenvalues, _ = np.linalg.svd(gram)
    return u[:, :dim], np.sqrt(np.maximum(eigenvalues[:dim], 0.0))


def orthonormalize(mat: np.ndarray) -> np.ndarray:
    """Return a matrix with the same column space but orthonormal columns.

    Required for anything read from plink2 `.eigenvec`: those columns are
    orthogonal but scaled by the singular values. Feeding them to `cca` unfixed
    silently yields wrong canonical correlations rather than an error.
    """
    q, _ = np.linalg.qr(mat)
    return q


def cca(
    u_x: np.ndarray, u_y: np.ndarray, dim: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Canonical correlation analysis between two orthonormal bases.

    With orthonormal inputs, the canonical correlations are exactly the singular
    values of u_x.T @ u_y.

    Returns (l_x, l_y, rho): the two projections and the canonical correlations.
    """
    for name, mat in (("u_x", u_x), ("u_y", u_y)):
        identity = np.eye(mat.shape[1])
        if not np.allclose(mat.T @ mat, identity, atol=1e-6):
            raise ValueError(
                f"{name} does not have orthonormal columns; pass it through "
                "orthonormalize() first (plink2 .eigenvec columns are not unit-norm)."
            )
    u_c, rho, v_c = np.linalg.svd(u_x.T @ u_y, full_matrices=False)
    l_x = u_x @ u_c[:, :dim]
    l_y = u_y @ v_c.T[:, :dim]
    return l_x, l_y, rho[:dim]


def loo_projection(
    X: np.ndarray, gram_y: np.ndarray, dim_x: int, dim_y: int, dim_z: int
) -> np.ndarray:
    """Leave-one-out PCCA coordinates for every sample.

    For each individual: drop them, recompute both views' principal components on
    the remainder, re-run CCA, and project the held-out individual onto the
    resulting expression-side canonical basis. If individuals land near their own
    population despite never contributing to the fit, the structure is real rather
    than an artefact of overfitting.

    Args:
        X: samples x genes expression residuals.
        gram_y: samples x samples genotype Gram matrix.
        dim_x, dim_y, dim_z: expression PCs, genotype PCs, canonical components.

    Returns:
        (n, dim_z) array of held-out projections, row order matching X.
    """
    n = X.shape[0]
    gram_x = X @ X.T
    out = np.empty((n, dim_z))

    # Reference solution on the full data. Singular vectors are only defined up
    # to sign, so each holdout's canonical directions come back at arbitrary
    # polarity; without alignment a population scatters to both extremes of an
    # axis and the plot shows no structure even when the projection is correct.
    u_x_full, _ = pca_from_gram(gram_x, dim_x)
    u_y_full, _ = pca_from_gram(gram_y, dim_y)
    l_x_full, _, _ = cca(u_x_full, u_y_full, dim_z)

    for i in range(n):
        idx = np.delete(np.arange(n), i)
        u_x, s_x = pca_from_gram(gram_x[np.ix_(idx, idx)], dim_x)
        u_y, _ = pca_from_gram(gram_y[np.ix_(idx, idx)], dim_y)
        l_x, _, _ = cca(u_x, u_y, dim_z)

        signs = np.sign(np.sum(l_x * l_x_full[idx], axis=0))
        signs[signs == 0] = 1.0
        l_x = l_x * signs

        # Port of project_i_cca. Writing the training matrix as X = U S Vᵀ, the
        # held-out sample's score on component j is (x_i · v_j) / s_j. Since
        # v_j = Xᵀ u_j / s_j, that equals (X x_i · u_j) / s_j² — computable from
        # the Gram-side quantities alone, with no need to form V.
        #
        # The 1/s² whitening is the part that matters: without it the projection
        # is dominated by the leading component and carries no population signal.
        cross = X[idx] @ X[i]
        scores = (cross @ u_x) / np.where(s_x**2 == 0, np.inf, s_x**2)
        # l_x = u_x @ u_c, so u_xᵀ l_x recovers the canonical rotation u_c.
        out[i] = scores @ (u_x.T @ l_x)
    return out
