"""Vectorized cis-eQTL scan.

The comparison the paper draws is between two covariate strategies, so the scan
runs twice over identical genotypes and phenotypes with only the covariates
swapped. A naive per-(gene, SNP) loop is far too slow; instead the genotypes are
residualized against the covariates once per chromosome, after which each test is
a dot product.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

CIS_WINDOW = 1_000_000


def residualize(mat: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Residuals of `mat` (samples x k) after regressing on covariates C."""
    beta, *_ = np.linalg.lstsq(C, mat, rcond=None)
    return mat - C @ beta


def association_pvalues(
    y: np.ndarray, G: np.ndarray, dof: int
) -> tuple[np.ndarray, np.ndarray]:
    """Two-sided p-values for the correlation of y with every column of G.

    Both inputs must already be covariate-residualized (and therefore centered).
    For simple linear regression the t-statistic depends only on the correlation,
    so this is exact rather than an approximation.
    """
    y_norm = np.linalg.norm(y)
    g_norms = np.linalg.norm(G, axis=0)
    g_norms = np.where(g_norms == 0, np.inf, g_norms)
    r = (G.T @ y) / (g_norms * y_norm)
    r = np.clip(r, -1 + 1e-12, 1 - 1e-12)
    t = r * np.sqrt(dof / (1 - r**2))
    p = 2 * stats.t.sf(np.abs(t), dof)
    return r, p


def cis_scan(
    ann: pd.DataFrame,
    X: np.ndarray,
    chrom_dosages: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
    covariates: np.ndarray,
    window: int = CIS_WINDOW,
) -> pd.DataFrame:
    """Best cis association per gene.

    Args:
        ann: gene annotation with TargetID, Chr, Coord (GRCh37).
        X: samples x genes expression, columns aligned to `ann` rows.
        chrom_dosages: chromosome -> (variant_ids, positions, samples x variants).
        covariates: samples x k design matrix, must include an intercept.
        window: cis distance either side of the gene's coordinate.

    Returns:
        DataFrame with gene, chrom, best_p, n_snps — one row per gene that had at
        least one cis SNP.
    """
    n_samples = X.shape[0]
    # Rank, not column count: the PCCA covariate set mixes canonical components
    # with expression PCs drawn from the same 30-PC subspace, so the design can be
    # rank-deficient. Using the column count there would overstate the degrees of
    # freedom and quietly inflate every p-value in the comparison.
    dof = n_samples - int(np.linalg.matrix_rank(covariates)) - 1
    X_resid = residualize(X, covariates)

    records = []
    for chrom, (ids, positions, G) in chrom_dosages.items():
        gene_rows = np.where(ann["Chr"].to_numpy() == chrom)[0]
        if len(gene_rows) == 0:
            continue
        G_resid = residualize(G.astype(float), covariates)

        for row in gene_rows:
            coord = int(ann["Coord"].iloc[row])
            in_window = np.abs(positions - coord) <= window
            n_snps = int(in_window.sum())
            if n_snps == 0:
                continue
            _, p = association_pvalues(X_resid[:, row], G_resid[:, in_window], dof)
            records.append(
                {
                    "gene": ann["TargetID"].iloc[row],
                    "chrom": chrom,
                    "best_p": float(np.min(p)),
                    "n_snps": n_snps,
                }
            )
    return pd.DataFrame.from_records(records).reset_index(drop=True)


def count_egenes(scan: pd.DataFrame, alpha: float) -> int:
    """Genes with at least one cis SNP below alpha."""
    return int((scan["best_p"] < alpha).sum())
