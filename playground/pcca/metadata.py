"""Cohort definition from the GEUVADIS SDRF.

The SDRF is the only source of the 7-lab assignment: the GD660 column headers
look like `HG00096.1.M_111124_6`, but all 660 runs share the prefix `M`, which
denotes the sequencing platform rather than the lab.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

ANCESTRY_TO_POP = {
    "British": "GBR",
    "Finnish": "FIN",
    "Tuscan": "TSI",
    "Utah": "CEU",
    "Yoruba": "YRI",
}


def load_sdrf(path) -> pd.DataFrame:
    """One row per individual, with the covariates the replication needs.

    The raw SDRF has two rows per individual (one per paired-end FASTQ);
    they are identical in every column used here.
    """
    raw = pd.read_csv(path, sep="\t", dtype=str)
    df = pd.DataFrame(
        {
            "sample_id": raw["Source Name"],
            "assay_name": raw["Assay Name"],
            "lab": raw["Characteristics[laboratory]"],
            "sex": raw["Characteristics[sex]"],
            "ancestry": raw["Characteristics[ancestry category]"],
            "has_phase1": raw["Comment[1000g Phase1 Genotypes]"] == "1",
        }
    )
    df["pop"] = df["ancestry"].map(ANCESTRY_TO_POP)
    if df["pop"].isna().any():
        unknown = sorted(df.loc[df["pop"].isna(), "ancestry"].unique())
        raise ValueError(f"Unmapped ancestry categories: {unknown}")
    return (
        df.drop_duplicates(subset="sample_id")
        .sort_values("sample_id")
        .reset_index(drop=True)
    )


def replication_cohort(sdrf: pd.DataFrame) -> pd.DataFrame:
    """The paper's exact 343 individuals: non-CEU, with Phase 1 genotypes.

    CEU is excluded in the paper for biased expression patterns.
    """
    mask = (sdrf["pop"] != "CEU") & sdrf["has_phase1"]
    return sdrf[mask].reset_index(drop=True)


def production_cohort(sdrf: pd.DataFrame) -> pd.DataFrame:
    """All 462 individuals including CEU — the input to the later fusion work."""
    return sdrf.reset_index(drop=True)


def build_confounders(meta: pd.DataFrame) -> np.ndarray:
    """Design matrix of intercept + 6 lab dummies + sex, shape (n, 8).

    One lab dummy is dropped to keep the matrix full rank alongside the intercept.
    """
    labs = sorted(meta["lab"].unique())
    dummies = np.column_stack(
        [(meta["lab"] == lab).to_numpy(dtype=float) for lab in labs[1:]]
    )
    sex = (meta["sex"] == "male").to_numpy(dtype=float)
    intercept = np.ones(len(meta))
    return np.column_stack([intercept, dummies, sex])
