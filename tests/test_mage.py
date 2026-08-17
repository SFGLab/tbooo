import numpy as np
import pandas as pd
import pytest

from pcca import paths
from pcca.mage import (
    SWAPPED_SAMPLES,
    build_confounders,
    load_expression_bed,
    load_genotype_pcs,
    load_sample_metadata,
    mage_cohort,
)

needs_data = pytest.mark.skipif(
    not paths.MAGE_SAMPLE_METADATA.exists(), reason="MAGE not downloaded"
)


def test_mage_cohort_drops_the_swapped_pair():
    """SRR19762653 is labelled NA11919 but carries HG00237's material.

    HG00237 also has its own library, so relabelling would duplicate it. Both are
    excluded until MAGE publishes corrected data.
    """
    meta = pd.DataFrame(
        {"sample_id": ["HG00096", "HG00237", "NA11919", "NA12878"],
         "pop": ["GBR", "GBR", "CEU", "CEU"]}
    )
    cohort = mage_cohort(meta)
    assert list(cohort["sample_id"]) == ["HG00096", "NA12878"]
    assert not set(cohort["sample_id"]) & set(SWAPPED_SAMPLES)
    assert list(cohort.index) == [0, 1], "index must reset for row alignment"


def test_mage_cohort_can_keep_the_swapped_pair():
    meta = pd.DataFrame({"sample_id": ["HG00237", "NA11919"], "pop": ["GBR", "CEU"]})
    assert len(mage_cohort(meta, drop_swapped=False)) == 2


def test_build_confounders_is_full_rank_and_has_an_intercept():
    cohort = pd.DataFrame(
        {"batch": ["1", "1", "2", "2", "3", "3"],
         "sex": ["male", "female", "male", "female", "male", "female"]}
    )
    conf = build_confounders(cohort)
    # intercept + (3 batches - 1) dummies + sex
    assert conf.shape == (6, 4)
    assert np.allclose(conf[:, 0], 1.0)
    assert np.linalg.matrix_rank(conf) == 4


def test_build_confounders_handles_a_single_batch():
    """A one-batch cohort must not produce a rank-deficient design."""
    cohort = pd.DataFrame({"batch": ["1", "1"], "sex": ["male", "female"]})
    conf = build_confounders(cohort)
    assert conf.shape == (2, 2)
    assert np.linalg.matrix_rank(conf) == 2


@needs_data
def test_metadata_has_all_731_across_five_continental_groups():
    meta = load_sample_metadata()
    assert len(meta) == 731
    assert meta["sample_id"].is_unique
    assert meta["continental_group"].value_counts().to_dict() == {
        "AFR": 196, "EUR": 142, "EAS": 141, "SAS": 139, "AMR": 113,
    }
    assert meta["pop"].nunique() == 26


@needs_data
def test_cohort_is_729_after_the_swap_exclusion():
    assert len(mage_cohort(load_sample_metadata())) == 729


@needs_data
def test_sex_is_normalised_to_the_shared_vocabulary():
    """MAGE encodes XX/XY; build_confounders elsewhere keys on 'male'."""
    assert set(load_sample_metadata()["sex"]) == {"male", "female"}


@needs_data
def test_expression_annotation_is_in_the_shared_format():
    ann, expr = load_expression_bed()
    assert list(ann.columns) == ["TargetID", "Gene_Symbol", "Chr", "Coord"]
    assert ann["TargetID"].str.startswith("ENSG").all()
    assert "." not in ann["TargetID"].iloc[0], "version suffix must be stripped"
    assert not ann["Chr"].str.startswith("chr").any()
    assert expr.shape[1] == 731


@needs_data
def test_every_cohort_member_has_expression_and_genotype_pcs():
    cohort = mage_cohort(load_sample_metadata())
    _, expr = load_expression_bed()
    assert not set(cohort["sample_id"]) - set(expr.columns)
    pcs = load_genotype_pcs(cohort)
    assert pcs.shape == (729, 5)
    assert np.all(np.isfinite(pcs))


@needs_data
def test_genotype_pcs_follow_cohort_row_order():
    """Misordered PCs would couple mismatched individuals in the CCA."""
    cohort = mage_cohort(load_sample_metadata())
    straight = load_genotype_pcs(cohort)
    shuffled_cohort = cohort.iloc[::-1].reset_index(drop=True)
    reversed_pcs = load_genotype_pcs(shuffled_cohort)
    assert np.allclose(straight, reversed_pcs[::-1])


@needs_data
def test_missing_sample_in_covariates_raises():
    cohort = pd.DataFrame({"sample_id": ["HG00096", "NOT_A_SAMPLE"]})
    with pytest.raises(KeyError, match="NOT_A_SAMPLE"):
        load_genotype_pcs(cohort)
