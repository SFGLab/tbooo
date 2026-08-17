import numpy as np
import pandas as pd
import pytest

from pcca import paths
from pcca.expression import (
    build_expression_matrix,
    filter_genes,
    keep_gene,
    load_gd660,
    quantile_normalize,
    regress_out,
    standardize_genes,
)


def test_keep_gene_accepts_autosomes():
    assert keep_gene("1", 1_000_000)
    assert keep_gene("22", 500_000)


def test_keep_gene_rejects_sex_and_mito():
    assert not keep_gene("X", 1_000_000)
    assert not keep_gene("Y", 1_000_000)
    assert not keep_gene("MT", 1_000)


def test_keep_gene_rejects_mhc_window_on_chr6_only():
    assert not keep_gene("6", 30_000_000)  # inside MHC
    assert keep_gene("6", 28_999_999)  # just before
    assert keep_gene("6", 33_000_001)  # just after
    assert keep_gene("5", 30_000_000)  # same coord, different chromosome


def test_filter_genes_drops_low_expression_and_bad_locations():
    ann = pd.DataFrame(
        {
            "TargetID": ["keep", "low", "mhc", "chrX"],
            "Gene_Symbol": ["keep", "low", "mhc", "chrX"],
            "Chr": ["1", "1", "6", "X"],
            "Coord": [1_000_000, 2_000_000, 30_000_000, 1_000_000],
        }
    )
    expr = pd.DataFrame(
        np.array([[5.0, 5.0], [0.01, 0.01], [5.0, 5.0], [5.0, 5.0]]),
        columns=["s1", "s2"],
    )
    ann_f, expr_f = filter_genes(ann, expr, min_mean_rpkm=0.1)
    assert list(ann_f["TargetID"]) == ["keep"]
    assert expr_f.shape == (1, 2)


def test_quantile_normalize_gives_columns_identical_distributions():
    mat = np.array([[5.0, 4.0], [2.0, 1.0], [3.0, 9.0]])
    out = quantile_normalize(mat)
    assert np.allclose(np.sort(out[:, 0]), np.sort(out[:, 1]))


def test_quantile_normalize_preserves_within_column_ranks():
    mat = np.array([[5.0, 4.0], [2.0, 1.0], [3.0, 9.0]])
    out = quantile_normalize(mat)
    for j in range(mat.shape[1]):
        assert list(np.argsort(mat[:, j])) == list(np.argsort(out[:, j]))


def test_quantile_normalize_preserves_shape():
    rng = np.random.default_rng(0)
    mat = rng.random((50, 7))
    assert quantile_normalize(mat).shape == (50, 7)


def test_standardize_genes_zero_mean_unit_variance_per_row():
    rng = np.random.default_rng(1)
    mat = rng.normal(loc=10.0, scale=3.0, size=(20, 100))
    out = standardize_genes(mat)
    assert np.allclose(out.mean(axis=1), 0.0, atol=1e-12)
    assert np.allclose(out.std(axis=1), 1.0, atol=1e-12)


def test_standardize_genes_handles_constant_rows():
    mat = np.array([[3.0, 3.0, 3.0], [1.0, 2.0, 3.0]])
    out = standardize_genes(mat)
    assert np.all(np.isfinite(out))
    assert np.allclose(out[0], 0.0)


def test_regress_out_removes_the_confounder_signal():
    """A gene that is a pure function of the confounder becomes ~zero."""
    rng = np.random.default_rng(2)
    n = 100
    group = rng.integers(0, 2, size=n).astype(float)
    C = np.column_stack([np.ones(n), group])
    gene_driven_by_group = 5.0 * group
    gene_independent = rng.normal(size=n)
    X = np.column_stack([gene_driven_by_group, gene_independent])
    resid = regress_out(X, C)
    assert np.allclose(resid[:, 0], 0.0, atol=1e-10)
    assert np.corrcoef(resid[:, 1], gene_independent)[0, 1] > 0.9


def test_regress_out_leaves_residuals_orthogonal_to_confounders():
    rng = np.random.default_rng(3)
    n = 60
    C = np.column_stack([np.ones(n), rng.normal(size=n)])
    X = rng.normal(size=(n, 10))
    resid = regress_out(X, C)
    assert np.allclose(C.T @ resid, 0.0, atol=1e-9)


@pytest.mark.skipif(
    not paths.GD660_GENE_RPKM.exists(), reason="GD660 not downloaded"
)
def test_load_gd660_real_file():
    ann, expr = load_gd660(paths.GD660_GENE_RPKM)
    assert expr.shape[1] == 660
    assert len(ann) == len(expr)
    assert list(ann.columns) == ["TargetID", "Gene_Symbol", "Chr", "Coord"]
    assert expr.columns[0] == "HG00096.1.M_111124_6"


@pytest.mark.skipif(
    not (paths.GD660_GENE_RPKM.exists() and paths.GENCODE_GTF.exists()),
    reason="GD660 or GENCODE GTF not available",
)
def test_protein_coding_filter_is_what_brings_the_gene_count_into_range():
    """GD660 is GENCODE v12 across all biotypes; the paper used protein-coding only.

    Without the biotype filter we keep ~7k lncRNAs and pseudogenes the paper never
    quantified, so the filter is what makes the gene count comparable at all.
    """
    ann, expr = load_gd660(paths.GD660_GENE_RPKM)
    all_biotypes, _ = filter_genes(ann, expr)
    coding_only, _ = filter_genes(ann, expr, gtf_path=paths.GENCODE_GTF)

    assert len(all_biotypes) > 18_000, "expected the unfiltered count to overshoot"
    # 12,736 observed against the paper's 14,070. The residual gap is genes
    # retired between GENCODE v12 and v47, plus transcript- vs gene-level
    # expression thresholding.
    assert 11_000 <= len(coding_only) <= 16_000, (
        f"got {len(coding_only)} protein-coding genes, paper reports 14,070"
    )


@pytest.mark.skipif(
    not (paths.GD660_GENE_RPKM.exists() and paths.SDRF.exists()),
    reason="GEUVADIS inputs not downloaded",
)
def test_build_expression_matrix_replication_cohort():
    from pcca.metadata import build_confounders, load_sdrf, replication_cohort

    cohort = replication_cohort(load_sdrf(paths.SDRF))
    ann, expr = load_gd660(paths.GD660_GENE_RPKM)
    ann_f, X = build_expression_matrix(ann, expr, cohort)

    assert X.shape[0] == 343, "one row per individual in the paper's cohort"
    assert X.shape[1] == len(ann_f)
    assert np.all(np.isfinite(X))
    # Residuals are orthogonal to the confounder design by construction.
    assert np.allclose(build_confounders(cohort).T @ X, 0.0, atol=1e-8)
