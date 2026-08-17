import numpy as np
import pytest

from pcca import paths
from pcca.metadata import (
    build_confounders,
    load_sdrf,
    production_cohort,
    replication_cohort,
)

pytestmark = pytest.mark.skipif(
    not paths.SDRF.exists(), reason="SDRF not downloaded; see Task 1 Step 1"
)


@pytest.fixture(scope="module")
def sdrf():
    return load_sdrf(paths.SDRF)


def test_one_row_per_individual(sdrf):
    assert len(sdrf) == 462
    assert sdrf["sample_id"].is_unique


def test_seven_labs(sdrf):
    assert sorted(sdrf["lab"].unique()) == ["1", "2", "3", "4", "5", "6", "7"]


def test_population_codes(sdrf):
    counts = sdrf["pop"].value_counts().to_dict()
    assert counts == {"FIN": 95, "GBR": 94, "TSI": 93, "CEU": 91, "YRI": 89}


def test_replication_cohort_is_the_papers_343(sdrf):
    """Non-CEU individuals with Phase 1 genotypes reproduce the paper's cohort exactly."""
    cohort = replication_cohort(sdrf)
    assert len(cohort) == 343
    assert cohort["pop"].value_counts().to_dict() == {
        "FIN": 89,
        "TSI": 92,
        "GBR": 85,
        "YRI": 77,
    }
    assert "CEU" not in set(cohort["pop"])


def test_production_cohort_is_all_462(sdrf):
    assert len(production_cohort(sdrf)) == 462


def test_confounders_shape_and_rank(sdrf):
    """Intercept + 6 lab dummies (7 labs, first dropped) + sex = 8 columns."""
    cohort = replication_cohort(sdrf)
    conf = build_confounders(cohort)
    assert conf.shape == (343, 8)
    assert np.linalg.matrix_rank(conf) == 8
    assert np.allclose(conf[:, 0], 1.0)
