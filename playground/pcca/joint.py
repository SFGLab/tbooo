"""Cross-cohort analysis in the shared genotype PC space.

MAGE distributes `MAGE_Geuvadis_GTEx_AFGR.top_PCs.txt.gz`: 2,373 samples from four
cohorts placed in one jointly-generated genotype PC space, with membership flags
and 20 PCs.

Two things fall out of it. GTEx donors carry no population labels — those are
dbGaP-protected — but they occupy the same coordinates as labelled 1KGP samples,
so their ancestry can be positioned against labelled anchors. And 93 individuals
appear in both MAGE and Geuvadis: the same people, sequenced by two independent
labs, which is a free concordance test for any embedding built on them.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from pcca import paths

DATASET_FLAGS = ["in_MAGE", "in_Geuvadis", "in_GTEx", "in_AFGR"]


def load_joint_pcs(path=None) -> pd.DataFrame:
    """The shared PC table, one row per sample across all four cohorts."""
    return pd.read_csv(path or paths.JOINT_PCS, sep="\t")


def dataset_membership(joint: pd.DataFrame) -> pd.Series:
    """Sample counts per cohort."""
    return pd.Series({flag: int(joint[flag].sum()) for flag in DATASET_FLAGS})


def shared_individuals(joint: pd.DataFrame, left: str, right: str) -> list[str]:
    """Sample IDs present in both named cohorts."""
    mask = (joint[f"in_{left}"] == 1) & (joint[f"in_{right}"] == 1)
    return joint.loc[mask, "sampleID"].tolist()


def infer_continental_group(joint: pd.DataFrame, target_flag: str = "in_GTEx",
                            n_pcs: int = 5, k: int = 15) -> pd.DataFrame:
    """Assign each unlabelled sample the majority group of its labelled neighbours.

    GTEx ancestry is not observed, but the joint space puts its donors alongside
    labelled 1KGP samples, so a k-nearest-neighbour vote over labelled anchors
    gives a defensible estimate. The result is an *inference*, and every consumer
    is expected to present it as such.
    """
    pc_cols = [f"PC{i}" for i in range(1, n_pcs + 1)]
    labelled = joint[joint["continentalGroup"].notna()]
    if labelled.empty:
        raise ValueError("no labelled samples to anchor the inference")

    anchors = labelled[pc_cols].to_numpy(dtype=float)
    anchor_groups = labelled["continentalGroup"].to_numpy()

    targets = joint[(joint[target_flag] == 1) & joint["continentalGroup"].isna()]
    if targets.empty:
        return pd.DataFrame(columns=["sampleID", "inferred_group", "confidence"])

    coords = targets[pc_cols].to_numpy(dtype=float)
    # Full pairwise distance: a few hundred by a few thousand is trivial here and
    # avoids a dependency on a neighbour index.
    distances = np.linalg.norm(coords[:, None, :] - anchors[None, :, :], axis=2)
    nearest = np.argsort(distances, axis=1)[:, :k]

    inferred, confidence = [], []
    for row in nearest:
        groups, counts = np.unique(anchor_groups[row], return_counts=True)
        winner = int(np.argmax(counts))
        inferred.append(groups[winner])
        confidence.append(counts[winner] / k)

    return pd.DataFrame({
        "sampleID": targets["sampleID"].to_numpy(),
        "inferred_group": inferred,
        "confidence": np.round(confidence, 3),
    })


def concordance(left_coords: pd.DataFrame, right_coords: pd.DataFrame,
                id_column: str = "sample_id",
                value_columns: tuple[str, str] = ("pcca_c1", "pcca_c2")):
    """Correlate embeddings of the same individuals from two cohorts.

    Signs are arbitrary in any SVD-derived embedding, so correlations are
    compared in absolute value; the question is whether the coordinate is
    reproducible, not which way it points.
    """
    merged = left_coords.merge(right_coords, on=id_column, suffixes=("_l", "_r"))
    result = {"n_shared": len(merged)}
    for column in value_columns:
        if len(merged) < 3:
            result[column] = float("nan")
            continue
        result[column] = round(
            abs(float(np.corrcoef(merged[f"{column}_l"], merged[f"{column}_r"])[0, 1])),
            3,
        )
    return result, merged
