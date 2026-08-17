"""GD660 expression matrix construction.

Pipeline: load raw RPKM -> filter genes -> quantile normalize -> standardize
-> regress out lab and sex. The result is the expression side of PCCA.
"""

from __future__ import annotations

import functools
import gzip
import re

import numpy as np
import pandas as pd

ANNOTATION_COLUMNS = ["TargetID", "Gene_Symbol", "Chr", "Coord"]
AUTOSOMES = {str(i) for i in range(1, 23)}
MHC_START, MHC_END = 29_000_000, 33_000_000
MIN_MEAN_RPKM = 0.1


@functools.lru_cache(maxsize=2)
def _gencode_protein_coding_genes(gtf_path: str) -> tuple[tuple[str, str, int], ...]:
    """Protein-coding genes from a GENCODE GTF as (gene_id, chrom, tss).

    One pass, cached, serving two needs: the protein-coding whitelist and — for
    cohorts whose expression matrices carry no coordinates, such as GTEx per-tissue
    TPM — the chromosome and transcription start site.

    TSS is strand-aware: the gene start for + strand, the end for - strand.
    Chromosome names are returned without the "chr" prefix to match the
    GEUVADIS/MAGE convention.
    """
    genes: list[tuple[str, str, int]] = []
    with gzip.open(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[2] != "gene" or 'gene_type "protein_coding"' not in fields[8]:
                continue
            match = re.search(r'gene_id "([^".]+)', fields[8])
            if not match:
                continue
            chrom = fields[0][3:] if fields[0].startswith("chr") else fields[0]
            start, end, strand = int(fields[3]), int(fields[4]), fields[6]
            genes.append((match.group(1), chrom, start if strand == "+" else end))
    return tuple(genes)


@functools.lru_cache(maxsize=2)
def protein_coding_gene_ids(gtf_path: str) -> frozenset[str]:
    """Unversioned Ensembl gene IDs with biotype protein_coding.

    The paper quantified GENCODE v27 *protein-coding transcripts only*, while
    GD660 is GENCODE v12 across all biotypes — so without this filter we retain
    ~6k extra lncRNAs and pseudogenes the paper never saw. Ensembl accessions are
    stable across releases, so the version suffix is stripped before matching;
    genes retired between v12 and the local GTF simply drop out.
    """
    return frozenset(gene_id for gene_id, _, _ in
                     _gencode_protein_coding_genes(gtf_path))


def gene_coordinates(gtf_path: str) -> pd.DataFrame:
    """Protein-coding gene coordinates: TargetID, Chr, Coord.

    GTEx per-tissue TPM files carry only gene ID and symbol, but `filter_genes`
    needs chromosome and position for the autosomal and MHC filters. Column names
    match the GD660 annotation frame so downstream code is cohort-agnostic.
    """
    genes = _gencode_protein_coding_genes(str(gtf_path))
    return pd.DataFrame(
        {
            "TargetID": [g for g, _, _ in genes],
            "Chr": [c for _, c, _ in genes],
            "Coord": [t for _, _, t in genes],
        }
    )


def load_gd660(path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (annotation, expression) with expression as genes x 660 runs."""
    df = pd.read_csv(path, sep="\t", compression="gzip")
    ann = df[ANNOTATION_COLUMNS].copy()
    ann["Chr"] = ann["Chr"].astype(str)
    ann["Coord"] = ann["Coord"].astype(int)
    expr = df.drop(columns=ANNOTATION_COLUMNS)
    return ann.reset_index(drop=True), expr.reset_index(drop=True)


def keep_gene(chrom: str, coord: int) -> bool:
    """Autosomal and outside the MHC, matching _keep_gene in pachterlab/PCCA.

    GD660 gives only a start coordinate, so MHC membership is tested on the start
    site rather than on interval overlap.
    """
    if chrom not in AUTOSOMES:
        return False
    if chrom == "6" and MHC_START <= coord <= MHC_END:
        return False
    return True


def filter_genes(
    ann: pd.DataFrame,
    expr: pd.DataFrame,
    min_mean_rpkm: float = MIN_MEAN_RPKM,
    gtf_path: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drop weakly expressed, non-autosomal, MHC, and non-protein-coding genes.

    `gtf_path` supplies the protein-coding whitelist; pass None to skip that
    filter (used by unit tests that build tiny synthetic annotations).
    """
    expressed = expr.mean(axis=1).to_numpy() >= min_mean_rpkm
    located = np.array(
        [keep_gene(c, p) for c, p in zip(ann["Chr"], ann["Coord"])], dtype=bool
    )
    keep = expressed & located
    if gtf_path is not None:
        coding = protein_coding_gene_ids(str(gtf_path))
        keep &= ann["TargetID"].str.split(".").str[0].isin(coding).to_numpy()
    return (
        ann[keep].reset_index(drop=True),
        expr[keep].reset_index(drop=True),
    )


def quantile_normalize(mat: np.ndarray) -> np.ndarray:
    """Force every column (sample) to share the same distribution.

    Equivalent to the reference implementation's pandas version, but done in
    numpy to avoid stack/unstack row-ordering surprises. Ties break by first
    occurrence, matching rank(method='first').
    """
    order = np.argsort(mat, axis=0, kind="stable")
    ranks = np.argsort(order, axis=0, kind="stable")
    reference = np.sort(mat, axis=0).mean(axis=1)
    return reference[ranks]


def standardize_genes(mat: np.ndarray) -> np.ndarray:
    """Center and scale each gene (row) to mean 0, variance 1.

    Constant genes would divide by zero; they are left at zero instead, which is
    the correct centered value and carries no signal either way.
    """
    mean = mat.mean(axis=1, keepdims=True)
    sd = mat.std(axis=1, keepdims=True)
    sd = np.where(sd == 0, 1.0, sd)
    return (mat - mean) / sd


def regress_out(X: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Residuals of samples x genes matrix X after regressing on design C."""
    beta, *_ = np.linalg.lstsq(C, X, rcond=None)
    return X - C @ beta


def build_expression_matrix(
    ann: pd.DataFrame, expr: pd.DataFrame, cohort: pd.DataFrame
) -> tuple[pd.DataFrame, np.ndarray]:
    """Full expression side of PCCA for one cohort.

    Selects the cohort's nonredundant runs (the SDRF lists exactly one run per
    individual, so the 198 technical replicates in GD660 are simply not selected),
    filters genes, quantile normalizes, standardizes, and regresses out lab and sex.

    Returns the surviving annotation and a samples x genes residual matrix whose
    row order matches `cohort`.
    """
    from pcca import paths
    from pcca.metadata import build_confounders

    missing = set(cohort["assay_name"]) - set(expr.columns)
    if missing:
        raise KeyError(
            f"{len(missing)} cohort assay names absent from GD660, "
            f"e.g. {sorted(missing)[:3]}"
        )

    selected = expr[list(cohort["assay_name"])]
    ann_f, expr_f = filter_genes(ann, selected, gtf_path=paths.GENCODE_GTF)

    mat = quantile_normalize(expr_f.to_numpy(dtype=float))
    mat = standardize_genes(mat)

    X = mat.T  # samples x genes
    X = regress_out(X, build_confounders(cohort))
    return ann_f, X
