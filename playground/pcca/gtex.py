"""GTEx cohort assembly, per tissue.

GTEx genotypes are dbGaP-protected, but `GTEx_Analysis_v10_eQTL_covariates.tar` is
open access and carries PC1-PC5 genotype PCs. Those PCs are identical across all
tissues for a given donor — verified — so they are global per-donor quantities and
serve as the CCA anchor without any dbGaP access.

Population labels are also protected, so GTEx ancestry is not observed directly.
The joint PC space in `pcca.joint` places these donors against labelled 1KGP
anchors instead.
"""

from __future__ import annotations

import gzip
import tarfile

import numpy as np
import pandas as pd

from pcca import paths

GCT_HEADER_LINES = 2  # "#1.2" then "<n_genes>\t<n_samples>"
_SEX_CODES = {1: "male", 2: "female"}


def donor_id(sample_id: str) -> str:
    """GTEX-1117F-0005-SM-HL9SH -> GTEX-1117F.

    Expression columns are sample-level; covariates and genotype PCs are
    donor-level. Everything is keyed on the donor.
    """
    return "-".join(sample_id.split("-")[:2])


def load_genotype_pcs(tissue: str, path=None, n_pcs: int = 5) -> pd.DataFrame:
    """Genotype PCs for one tissue's donors, indexed by donor ID."""
    path = path or paths.GTEX_COVARIATES_TAR
    member = f"{_covariate_member(tissue)}"
    with tarfile.open(path) as tar:
        try:
            handle = tar.extractfile(member)
        except KeyError:
            handle = None
        if handle is None:
            available = [m for m in tar.getnames() if m.endswith(".txt")][:3]
            raise KeyError(
                f"{member} not in {path.name}; examples present: {available}"
            )
        cov = pd.read_csv(handle, sep="\t", index_col=0)

    wanted = [f"PC{i}" for i in range(1, n_pcs + 1)]
    missing = [pc for pc in wanted if pc not in cov.index]
    if missing:
        raise KeyError(f"genotype PCs absent for {tissue}: {missing}")
    return cov.loc[wanted].T


def _covariate_member(tissue: str) -> str:
    """whole_blood -> Whole_Blood.v10.covariates.txt."""
    return "_".join(part.capitalize() for part in tissue.split("_")) + \
        ".v10.covariates.txt"


def load_tpm(tissue: str, path=None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (annotation, expression) from a GTEx per-tissue GCT.

    Columns are collapsed from sample-level to donor-level; a donor contributing
    more than one sample to a tissue is averaged, which keeps one row per person
    as every downstream step assumes.

    The GCT carries no coordinates, so `Chr`/`Coord` are left to be joined from
    GENCODE by `attach_coordinates`.
    """
    path = path or paths.gtex_tpm(tissue)
    with gzip.open(path, "rt") as handle:
        for _ in range(GCT_HEADER_LINES):
            handle.readline()
        df = pd.read_csv(handle, sep="\t")

    ann = pd.DataFrame(
        {
            "TargetID": df["Name"].str.split(".").str[0],
            "Gene_Symbol": df["Description"],
        }
    )
    expr = df.drop(columns=["Name", "Description"])
    donors = [donor_id(c) for c in expr.columns]
    if len(set(donors)) != len(donors):
        expr = expr.T.groupby(donors).mean().T
    else:
        expr.columns = donors
    return ann.reset_index(drop=True), expr.reset_index(drop=True)


def attach_coordinates(ann: pd.DataFrame, expr: pd.DataFrame, gtf_path=None):
    """Join GENCODE coordinates onto GTEx genes, dropping those without.

    GTEx per-tissue TPM gives only gene ID and symbol, but the autosomal and MHC
    filters need position. Genes absent from the protein-coding annotation are
    dropped here rather than later, which also applies the biotype filter.
    """
    from pcca.expression import gene_coordinates

    coords = gene_coordinates(str(gtf_path or paths.GENCODE_GTF))
    merged = ann.merge(coords, on="TargetID", how="left")
    keep = merged["Chr"].notna().to_numpy()
    merged = merged[keep].reset_index(drop=True)
    merged["Coord"] = merged["Coord"].astype(int)
    return merged, expr[keep].reset_index(drop=True)


def load_sample_attributes(path=None) -> pd.DataFrame:
    """Per-sample technical covariates: tissue, collection centre, batch."""
    path = path or paths.GTEX_SAMPLE_ATTRIBUTES
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    return pd.DataFrame(
        {
            "sample_id": raw["SAMPID"],
            "donor": raw["SAMPID"].map(donor_id),
            "tissue": raw["SMTSD"],
            "centre": raw.get("SMCENTER", pd.Series(index=raw.index, dtype=object)),
            "batch": raw.get("SMNABTCH", pd.Series(index=raw.index, dtype=object)),
        }
    )


def load_subject_phenotypes(path=None) -> pd.DataFrame:
    """Donor-level sex and age (the only openly released phenotypes)."""
    path = path or paths.GTEX_SUBJECT_PHENOTYPES
    raw = pd.read_csv(path, sep="\t")
    return pd.DataFrame(
        {
            "donor": raw["SUBJID"],
            "sex": raw["SEX"].map(_SEX_CODES),
            "age": raw["AGE"],
        }
    )


def build_cohort(tissue: str, expr: pd.DataFrame, genotype_pcs: pd.DataFrame,
                 attributes: pd.DataFrame, phenotypes: pd.DataFrame) -> pd.DataFrame:
    """Donors present in both expression and genotype PCs, with covariates.

    Returned rows are the intersection in expression-column order, so callers can
    index the expression matrix and the PC matrix with the same ordering. Every
    cohort here joins expression to genotype through a different identifier, and
    misalignment yields plausible-but-meaningless CCA rather than an error.
    """
    donors = [d for d in expr.columns if d in set(genotype_pcs.index)]
    if not donors:
        raise ValueError(f"no donors shared between expression and PCs for {tissue}")

    tissue_attrs = (
        attributes[attributes["donor"].isin(donors)]
        .drop_duplicates(subset="donor", keep="first")
        .set_index("donor")
    )
    pheno = phenotypes.drop_duplicates(subset="donor").set_index("donor")

    cohort = pd.DataFrame({"sample_id": donors})
    cohort["centre"] = cohort["sample_id"].map(tissue_attrs["centre"]).fillna("unknown")
    cohort["batch"] = cohort["sample_id"].map(tissue_attrs["batch"]).fillna("unknown")
    cohort["sex"] = cohort["sample_id"].map(pheno["sex"]).fillna("unknown")
    cohort["tissue"] = tissue
    return cohort.reset_index(drop=True)


def build_confounders(cohort: pd.DataFrame) -> np.ndarray:
    """Intercept + centre dummies + batch dummies + sex.

    GTEx's analogue of the GEUVADIS laboratory indicators. Batch has many levels,
    so it is collapsed to the most frequent ones plus an "other" bucket to keep
    the design from approaching the sample count.
    """
    columns = [np.ones(len(cohort))]
    for field, max_levels in (("centre", 12), ("batch", 20)):
        values = cohort[field].astype(str)
        top = values.value_counts().index[:max_levels]
        collapsed = values.where(values.isin(top), "other")
        levels = sorted(collapsed.unique())[1:]  # drop first against the intercept
        for level in levels:
            columns.append((collapsed == level).to_numpy(dtype=float))
    columns.append((cohort["sex"] == "male").to_numpy(dtype=float))
    return np.column_stack(columns)
