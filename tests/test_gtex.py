import numpy as np
import pandas as pd
import pytest

from pcca import paths
from pcca.gtex import (
    build_cohort,
    build_confounders,
    donor_id,
    load_genotype_pcs,
    load_sample_attributes,
    load_subject_phenotypes,
    load_tpm,
)

needs_data = pytest.mark.skipif(
    not paths.GTEX_COVARIATES_TAR.exists(), reason="GTEx not downloaded"
)


def test_donor_id_strips_the_sample_suffix():
    """Expression columns are sample-level; genotype PCs are donor-level."""
    assert donor_id("GTEX-1117F-0005-SM-HL9SH") == "GTEX-1117F"
    assert donor_id("GTEX-111CU-0006-SM-5NQBE") == "GTEX-111CU"


def test_donor_id_is_idempotent_on_an_already_bare_donor():
    assert donor_id("GTEX-1117F") == "GTEX-1117F"


def test_build_cohort_keeps_only_donors_with_both_modalities():
    expr = pd.DataFrame({"GTEX-A": [1.0], "GTEX-B": [2.0], "GTEX-C": [3.0]})
    pcs = pd.DataFrame({"PC1": [0.1, 0.2]}, index=["GTEX-A", "GTEX-C"])
    attrs = pd.DataFrame({"sample_id": ["GTEX-A-1", "GTEX-C-1"],
                          "donor": ["GTEX-A", "GTEX-C"],
                          "tissue": ["Lung", "Lung"],
                          "centre": ["B1", "B1"], "batch": ["X", "Y"]})
    phen = pd.DataFrame({"donor": ["GTEX-A", "GTEX-C"],
                         "sex": ["male", "female"], "age": ["50-59", "60-69"]})
    cohort = build_cohort("lung", expr, pcs, attrs, phen)
    assert list(cohort["sample_id"]) == ["GTEX-A", "GTEX-C"]
    assert list(cohort["sex"]) == ["male", "female"]


def test_build_cohort_preserves_expression_column_order():
    """Row order must follow expression, or the CCA couples the wrong donors."""
    expr = pd.DataFrame({"GTEX-C": [1.0], "GTEX-A": [2.0]})
    pcs = pd.DataFrame({"PC1": [0.1, 0.2]}, index=["GTEX-A", "GTEX-C"])
    empty_attrs = pd.DataFrame(columns=["sample_id", "donor", "tissue", "centre", "batch"])
    empty_phen = pd.DataFrame(columns=["donor", "sex", "age"])
    cohort = build_cohort("lung", expr, pcs, empty_attrs, empty_phen)
    assert list(cohort["sample_id"]) == ["GTEX-C", "GTEX-A"]


def test_build_cohort_raises_when_nothing_overlaps():
    expr = pd.DataFrame({"GTEX-Z": [1.0]})
    pcs = pd.DataFrame({"PC1": [0.1]}, index=["GTEX-A"])
    empty_attrs = pd.DataFrame(columns=["sample_id", "donor", "tissue", "centre", "batch"])
    empty_phen = pd.DataFrame(columns=["donor", "sex", "age"])
    with pytest.raises(ValueError, match="no donors shared"):
        build_cohort("lung", expr, pcs, empty_attrs, empty_phen)


def test_build_confounders_collapses_high_cardinality_batches():
    """GTEx has hundreds of batches; an uncollapsed design would swamp n."""
    cohort = pd.DataFrame({
        "centre": ["B1"] * 60,
        "batch": [f"BP-{i}" for i in range(60)],
        "sex": ["male", "female"] * 30,
    })
    conf = build_confounders(cohort)
    assert conf.shape[1] < 30, "batch levels must be capped"
    assert np.allclose(conf[:, 0], 1.0)


@needs_data
def test_genotype_pcs_are_identical_across_tissues():
    """They are global per-donor PCs, which is what makes them a valid anchor."""
    blood = load_genotype_pcs("whole_blood")
    muscle = load_genotype_pcs("muscle_skeletal")
    shared = blood.index.intersection(muscle.index)
    assert len(shared) > 500
    assert np.allclose(blood.loc[shared].to_numpy(), muscle.loc[shared].to_numpy())


@needs_data
def test_genotype_pcs_have_five_components_per_donor():
    pcs = load_genotype_pcs("whole_blood")
    assert list(pcs.columns) == ["PC1", "PC2", "PC3", "PC4", "PC5"]
    assert pcs.index.is_unique
    assert pcs.index.str.startswith("GTEX-").all()


@needs_data
def test_unknown_tissue_raises():
    with pytest.raises(KeyError, match="not in"):
        load_genotype_pcs("not_a_tissue")


@needs_data
def test_tpm_columns_are_unique_donors():
    _, expr = load_tpm("whole_blood")
    assert expr.columns.is_unique, "donors must be collapsed to one column each"
    assert expr.columns.str.startswith("GTEX-").all()


@needs_data
def test_coordinates_come_from_gencode_and_filter_to_protein_coding():
    """GTEx TPM carries no coordinates; filter_genes needs them."""
    from pcca.gtex import attach_coordinates

    ann, expr = load_tpm("whole_blood")
    assert "Chr" not in ann.columns
    ann2, expr2 = attach_coordinates(ann, expr)
    assert list(ann2.columns) == ["TargetID", "Gene_Symbol", "Chr", "Coord"]
    assert len(ann2) == len(expr2)
    assert len(ann2) < len(ann), "non-coding genes should drop out"
    assert 15_000 <= len(ann2) <= 22_000
    assert ann2["Coord"].dtype.kind == "i"


@needs_data
def test_end_to_end_cohort_for_whole_blood():
    from pcca.gtex import attach_coordinates

    pcs = load_genotype_pcs("whole_blood")
    ann, expr = attach_coordinates(*load_tpm("whole_blood"))
    cohort = build_cohort("whole_blood", expr, pcs,
                          load_sample_attributes(), load_subject_phenotypes())
    assert len(cohort) > 700
    assert cohort["sample_id"].is_unique
    conf = build_confounders(cohort)
    assert conf.shape[0] == len(cohort)
    assert np.linalg.matrix_rank(conf) == conf.shape[1], "design must be full rank"
