import numpy as np
import pandas as pd
import pytest
from scipy import stats

from pcca.eqtl import association_pvalues, cis_scan, count_egenes, residualize


def test_association_pvalues_match_scipy_ols_without_covariates():
    """Ground truth: with only an intercept, this reduces to Pearson correlation."""
    rng = np.random.default_rng(0)
    n = 200
    G = rng.integers(0, 3, size=(n, 5)).astype(float)
    y = G[:, 0] * 0.4 + rng.normal(size=n)

    y_c = y - y.mean()
    G_c = G - G.mean(axis=0)
    _, p = association_pvalues(y_c, G_c, dof=n - 2)

    for j in range(5):
        expected = stats.pearsonr(y, G[:, j])[1]
        assert p[j] == pytest.approx(expected, rel=1e-6)


def test_association_pvalues_detect_a_planted_effect():
    rng = np.random.default_rng(1)
    n = 300
    G = rng.integers(0, 3, size=(n, 10)).astype(float)
    y = 1.5 * G[:, 3] + rng.normal(size=n)
    y_c = y - y.mean()
    G_c = G - G.mean(axis=0)
    _, p = association_pvalues(y_c, G_c, dof=n - 2)
    assert np.argmin(p) == 3
    assert p[3] < 1e-20


def test_residualize_leaves_output_orthogonal_to_covariates():
    rng = np.random.default_rng(2)
    n = 50
    C = np.column_stack([np.ones(n), rng.normal(size=n)])
    M = rng.normal(size=(n, 8))
    resid = residualize(M, C)
    assert np.allclose(C.T @ resid, 0.0, atol=1e-9)


def test_cis_scan_only_tests_snps_inside_the_window():
    """A SNP 15 Mb away must not be tested even if perfectly correlated."""
    rng = np.random.default_rng(3)
    n = 100
    ann = pd.DataFrame({"TargetID": ["G1"], "Chr": ["1"], "Coord": [5_000_000]})
    driver = rng.normal(size=n)
    X = (driver + 0.01 * rng.normal(size=n)).reshape(n, 1)
    # Column 0 is in-window noise; column 1 is far away but identical to the driver.
    G = np.column_stack([rng.normal(size=n), driver])
    positions = np.array([5_100_000, 20_000_000])
    ids = np.array(["1:5100000:A:G", "1:20000000:A:G"])

    C = np.ones((n, 1))
    result = cis_scan(ann, X, {"1": (ids, positions, G)}, C, window=1_000_000)
    assert len(result) == 1
    assert result.loc[0, "n_snps"] == 1


def test_cis_scan_reports_the_best_p_per_gene():
    rng = np.random.default_rng(4)
    n = 200
    ann = pd.DataFrame({"TargetID": ["G1"], "Chr": ["1"], "Coord": [1_000_000]})
    G = rng.integers(0, 3, size=(n, 4)).astype(float)
    X = (2.0 * G[:, 2] + rng.normal(size=n)).reshape(n, 1)
    ids = np.array([f"1:{1_000_000 + i}:A:G" for i in range(4)])
    positions = np.array([1_000_000 + i for i in range(4)])

    C = np.ones((n, 1))
    result = cis_scan(ann, X, {"1": (ids, positions, G)}, C, window=1_000_000)
    assert result.loc[0, "n_snps"] == 4
    assert result.loc[0, "best_p"] < 1e-10


def test_cis_scan_uses_rank_not_column_count_for_dof():
    """A rank-deficient covariate design must not inflate the degrees of freedom.

    Duplicating a covariate column changes the column count but not the rank, so
    the p-values must be identical.
    """
    rng = np.random.default_rng(5)
    n = 150
    ann = pd.DataFrame({"TargetID": ["G1"], "Chr": ["1"], "Coord": [1_000_000]})
    G = rng.integers(0, 3, size=(n, 6)).astype(float)
    X = (1.2 * G[:, 1] + rng.normal(size=n)).reshape(n, 1)
    ids = np.array([f"1:{1_000_000 + i}:A:G" for i in range(6)])
    positions = np.array([1_000_000 + i for i in range(6)])
    payload = {"1": (ids, positions, G)}

    covar = rng.normal(size=n)
    full_rank = np.column_stack([np.ones(n), covar])
    duplicated = np.column_stack([np.ones(n), covar, covar])

    a = cis_scan(ann, X, payload, full_rank)
    b = cis_scan(ann, X, payload, duplicated)
    assert a.loc[0, "best_p"] == pytest.approx(b.loc[0, "best_p"], rel=1e-12)


def test_count_egenes_counts_genes_below_alpha():
    scan = pd.DataFrame({"best_p": [1e-9, 1e-3, 1e-7]})
    assert count_egenes(scan, 1e-6) == 2
