import numpy as np
import pandas as pd
import pytest

from pcca import paths
from pcca.joint import (
    concordance,
    dataset_membership,
    infer_continental_group,
    load_joint_pcs,
    shared_individuals,
)

needs_data = pytest.mark.skipif(
    not paths.JOINT_PCS.exists(), reason="MAGE joint PCs not downloaded"
)


def _toy_joint():
    """Two tight labelled clusters plus unlabelled points sitting inside each."""
    rows = []
    for i in range(20):
        rows.append({"sampleID": f"AFR{i}", "continentalGroup": "AFR",
                     "in_MAGE": 1, "in_Geuvadis": 0, "in_GTEx": 0, "in_AFGR": 0,
                     "PC1": 1.0 + i * 0.001, "PC2": 0.0, "PC3": 0.0,
                     "PC4": 0.0, "PC5": 0.0})
        rows.append({"sampleID": f"EUR{i}", "continentalGroup": "EUR",
                     "in_MAGE": 1, "in_Geuvadis": 0, "in_GTEx": 0, "in_AFGR": 0,
                     "PC1": -1.0 - i * 0.001, "PC2": 0.0, "PC3": 0.0,
                     "PC4": 0.0, "PC5": 0.0})
    rows.append({"sampleID": "GTEX-1", "continentalGroup": np.nan,
                 "in_MAGE": 0, "in_Geuvadis": 0, "in_GTEx": 1, "in_AFGR": 0,
                 "PC1": 0.98, "PC2": 0.0, "PC3": 0.0, "PC4": 0.0, "PC5": 0.0})
    rows.append({"sampleID": "GTEX-2", "continentalGroup": np.nan,
                 "in_MAGE": 0, "in_Geuvadis": 0, "in_GTEx": 1, "in_AFGR": 0,
                 "PC1": -0.98, "PC2": 0.0, "PC3": 0.0, "PC4": 0.0, "PC5": 0.0})
    return pd.DataFrame(rows)


def test_infer_assigns_unlabelled_samples_to_their_nearest_cluster():
    result = infer_continental_group(_toy_joint(), k=5)
    groups = dict(zip(result["sampleID"], result["inferred_group"]))
    assert groups == {"GTEX-1": "AFR", "GTEX-2": "EUR"}
    assert (result["confidence"] == 1.0).all()


def test_infer_only_targets_unlabelled_members_of_the_named_cohort():
    result = infer_continental_group(_toy_joint(), k=5)
    assert set(result["sampleID"]) == {"GTEX-1", "GTEX-2"}


def test_infer_raises_without_labelled_anchors():
    joint = _toy_joint()
    joint["continentalGroup"] = np.nan
    with pytest.raises(ValueError, match="no labelled samples"):
        infer_continental_group(joint, k=5)


def test_shared_individuals_finds_the_intersection():
    joint = pd.DataFrame({
        "sampleID": ["a", "b", "c"],
        "in_MAGE": [1, 1, 0], "in_Geuvadis": [0, 1, 1],
        "in_GTEx": [0, 0, 0], "in_AFGR": [0, 0, 0],
    })
    assert shared_individuals(joint, "MAGE", "Geuvadis") == ["b"]


def test_concordance_is_sign_invariant():
    """Embedding signs are arbitrary; reproducibility is the question."""
    left = pd.DataFrame({"sample_id": list("abcdef"),
                         "pcca_c1": [1.0, 2, 3, 4, 5, 6],
                         "pcca_c2": [1.0, 0, 1, 0, 1, 0]})
    right = left.copy()
    right["pcca_c1"] = -right["pcca_c1"]  # flipped axis
    result, merged = concordance(left, right)
    assert result["n_shared"] == 6
    assert result["pcca_c1"] == 1.0
    assert len(merged) == 6


def test_concordance_reports_low_agreement_for_unrelated_embeddings():
    rng = np.random.default_rng(0)
    left = pd.DataFrame({"sample_id": [f"s{i}" for i in range(50)],
                         "pcca_c1": rng.normal(size=50),
                         "pcca_c2": rng.normal(size=50)})
    right = pd.DataFrame({"sample_id": [f"s{i}" for i in range(50)],
                          "pcca_c1": rng.normal(size=50),
                          "pcca_c2": rng.normal(size=50)})
    result, _ = concordance(left, right)
    assert result["pcca_c1"] < 0.5


@needs_data
def test_real_joint_space_spans_four_cohorts():
    joint = load_joint_pcs()
    counts = dataset_membership(joint)
    assert counts["in_MAGE"] == 731
    assert counts["in_GTEx"] == 837
    assert counts["in_Geuvadis"] == 452
    assert len(joint) == 2373


@needs_data
def test_mage_and_geuvadis_share_93_individuals():
    """The same people from two independent experiments — a free control."""
    assert len(shared_individuals(load_joint_pcs(), "MAGE", "Geuvadis")) == 93


@needs_data
def test_gtex_donors_are_unlabelled_but_inferable():
    joint = load_joint_pcs()
    gtex_rows = joint[joint["in_GTEx"] == 1]
    assert gtex_rows["continentalGroup"].isna().all(), "GTEx labels are protected"
    inferred = infer_continental_group(joint)
    assert len(inferred) == len(gtex_rows)
    # GTEx is ~85% European, so EUR must dominate any sane inference.
    assert inferred["inferred_group"].value_counts().idxmax() == "EUR"
