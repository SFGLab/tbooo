"""MAGE cohort assembly.

MAGE is 731 LCL RNA-seq samples from 1000 Genomes individuals across 26
populations and five continental groups — the ancestry-diverse counterpart to
GEUVADIS, whose four populations included only one non-European group.

Expression arrives as BED4+ with TSS coordinates, and covariates (sex, 5 genotype
PCs, 60 PEER factors) ship separately from the expression matrix. The released
matrices are therefore not covariate-corrected, so the GEUVADIS trap — a
distributed file that was silently PEER-residualized — cannot recur here.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from pcca import paths

# Unresolved sample swap, MAGE README dated 2026-03-27 and not fixed in the
# Zenodo release: SRR19762653 is labelled NA11919 but carries HG00237's material.
# HG00237 also has its own library, so it is effectively duplicated and true
# NA11919 is absent. Relabelling would duplicate HG00237; excluding both costs
# 0.3% of the cohort and cannot propagate a mislabel.
SWAPPED_SAMPLES = ("HG00237", "NA11919")

BED_ANNOTATION_COLUMNS = ["#Chr", "start", "end", "ID"]


def load_sample_metadata(path=None) -> pd.DataFrame:
    """One row per MAGE individual with the covariates the analysis needs."""
    path = path or paths.MAGE_SAMPLE_METADATA
    raw = pd.read_csv(path, sep="\t")
    return pd.DataFrame(
        {
            "sample_id": raw["sample_kgpID"],
            "coriell_id": raw["sample_coriellID"],
            "sra": raw["SRA_accession"],
            "pop": raw["population"],
            "continental_group": raw["continentalGroup"],
            "batch": raw["batch"].astype(str),
            # MAGE encodes sex as XX/XY; normalize to the male flag used by
            # build_confounders elsewhere in this package.
            "sex": raw["sex"].map({"XY": "male", "XX": "female"}),
        }
    ).sort_values("sample_id").reset_index(drop=True)


def mage_cohort(metadata: pd.DataFrame, drop_swapped: bool = True) -> pd.DataFrame:
    """The analysis cohort: 729 individuals after excluding the swapped pair."""
    if drop_swapped:
        metadata = metadata[~metadata["sample_id"].isin(SWAPPED_SAMPLES)]
    return metadata.reset_index(drop=True)


def load_expression_bed(path=None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (annotation, expression) from the MAGE BED4+ matrix.

    Annotation columns are renamed to the GD660 convention (TargetID/Chr/Coord)
    so `pcca.expression.filter_genes` works unchanged across cohorts.
    """
    path = path or paths.MAGE_EXPRESSION
    df = pd.read_csv(path, sep="\t")
    ann = pd.DataFrame(
        {
            # Versioned Ensembl IDs; strip the suffix so the protein-coding
            # whitelist (keyed on bare accessions) matches.
            "TargetID": df["ID"].str.split(".").str[0],
            "Gene_Symbol": df["ID"],
            "Chr": df["#Chr"].astype(str).str.replace("^chr", "", regex=True),
            # BED is half-open and zero-based; `end` is the TSS in MAGE's layout.
            "Coord": df["end"].astype(int),
        }
    )
    expr = df.drop(columns=BED_ANNOTATION_COLUMNS)
    return ann.reset_index(drop=True), expr.reset_index(drop=True)


def load_genotype_pcs(cohort: pd.DataFrame, path=None, n_pcs: int = 5) -> np.ndarray:
    """MAGE's published genotype PCs, ordered to match `cohort` rows.

    The covariates file is covariates x samples. Row order must follow the
    expression matrix or the CCA silently couples mismatched individuals.
    """
    path = path or paths.MAGE_COVARIATES
    cov = pd.read_csv(path, sep="\t", index_col=0)
    wanted = [f"PC{i}" for i in range(1, n_pcs + 1)]
    missing_pcs = [pc for pc in wanted if pc not in cov.index]
    if missing_pcs:
        raise KeyError(f"genotype PCs absent from {path}: {missing_pcs}")

    pcs = cov.loc[wanted].T
    missing = [s for s in cohort["sample_id"] if s not in pcs.index]
    if missing:
        raise KeyError(
            f"{len(missing)} cohort samples absent from covariates: {missing[:5]}"
        )
    return pcs.loc[cohort["sample_id"]].to_numpy(dtype=float)


def build_confounders(cohort: pd.DataFrame) -> np.ndarray:
    """Intercept + batch dummies (first dropped) + sex.

    Mirrors the GEUVADIS design, where 7 laboratory indicators played the role
    MAGE's sequencing batch plays here.
    """
    batches = sorted(cohort["batch"].unique())
    dummies = np.column_stack(
        [(cohort["batch"] == b).to_numpy(dtype=float) for b in batches[1:]]
    ) if len(batches) > 1 else np.empty((len(cohort), 0))
    sex = (cohort["sex"] == "male").to_numpy(dtype=float)
    return np.column_stack([np.ones(len(cohort)), dummies, sex])
