"""Assign GTEx donors synthetic EIDs and merge their expression PCs.

Unlike GEUVADIS and MAGE, GTEx donors are not 1000 Genomes samples, so they have
no existing EID and are appended to participant.parquet as new rows. TBOOO already
partitions participants by EID range — 1KGP at 1,000,000+, SGDP at 2,000,000+ —
and GTEx takes 3,000,000+.

GTEx rows carry expression PCs but null genetic-availability fields, because GTEx
genotypes are dbGaP-protected. This mirrors UK Biobank, where different assays
cover different subsets of participants.

Whole blood is the tissue represented in participant.parquet: it is the largest
cohort and the closest analogue to UK Biobank's own blood-derived assays. The
per-tissue analysis outputs live in data/metadata/ alongside.

Outputs:
    data/metadata/eid_map_gtex.tsv          columns: eid, donor
    data/metadata/gtex_expression_pcs.tsv   columns: eid, gtex_pc1..N
    data/Showcase/participant.parquet       (updated in place, rows appended)
"""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from sklearn.decomposition import PCA

from tbooo.config import Config
from tbooo.download.gtex import SUBJECT_PHENOTYPES_FILE
from tbooo.integrity import table_ok
from tbooo.utils import ensure_dirs, log

N_PCS = 20
N_HVG = 2000
EID_BASE = 3_000_000
REPRESENTATIVE_TISSUE = "whole_blood"
_GCT_HEADER_LINES = 2


def _donor_id(sample_id: str) -> str:
    return "-".join(sample_id.split("-")[:2])


def build_expression_pcs(cfg: Config) -> None:
    tpm_path = cfg.gtex_raw_dir() / f"gene_tpm_v10_{REPRESENTATIVE_TISSUE}.gct.gz"
    if not tpm_path.exists():
        raise FileNotFoundError(
            f"GTEx expression not found: {tpm_path}\nRun `tbooo download gtex` first."
        )

    ensure_dirs(cfg.metadata_dir())
    out_pcs = cfg.metadata_dir() / "gtex_expression_pcs.tsv"
    out_map = cfg.metadata_dir() / "eid_map_gtex.tsv"
    if table_ok(out_pcs, min_lines=2) and table_ok(out_map, min_lines=2):
        # The saved rows still go into participant.parquet, which may have been
        # built or rebuilt since they were written.
        log(f"  skip (valid): {out_pcs.name} + {out_map.name}")
        saved = pd.read_csv(out_pcs, sep="\t")
        _patch_parquet(cfg, saved, [c for c in saved.columns if c != "eid"])
        return

    log(f"Loading GTEx {REPRESENTATIVE_TISSUE} TPM…")
    with gzip.open(tpm_path, "rt") as handle:
        for _ in range(_GCT_HEADER_LINES):
            handle.readline()
        df = pd.read_csv(handle, sep="\t")

    expr = df.drop(columns=["Name", "Description"])
    donors = [_donor_id(c) for c in expr.columns]
    if len(set(donors)) != len(donors):
        expr = expr.T.groupby(donors).mean().T
    else:
        expr.columns = donors
    log(f"  {len(df)} genes × {expr.shape[1]} donors")

    values = np.nan_to_num(expr.to_numpy(dtype=float).T, nan=0.0,
                           posinf=0.0, neginf=0.0)
    X = np.log2(np.clip(values, 0.0, None) + 0.1)
    X = X[:, np.argsort(X.var(axis=0))[-N_HVG:]]
    log(f"  {X.shape[1]} highly variable genes selected")

    log(f"Running PCA (n_components={N_PCS})…")
    pca = PCA(n_components=N_PCS)
    scores = pca.fit_transform(X)
    log(f"  variance explained by {N_PCS} PCs: {pca.explained_variance_ratio_.sum():.1%}")

    # Deterministic EIDs: sorted donor order, so re-running is stable.
    donor_ids = sorted(expr.columns)
    eid_map = pd.DataFrame({
        "eid": [EID_BASE + i for i in range(len(donor_ids))],
        "donor": donor_ids,
    })
    eid_map.to_csv(out_map, sep="\t", index=False)
    log(f"  wrote {out_map} ({len(eid_map)} donors, EIDs from {EID_BASE})")

    pc_cols = [f"gtex_pc{i}" for i in range(1, N_PCS + 1)]
    scores_df = pd.DataFrame(scores, index=list(expr.columns), columns=pc_cols)
    scores_df.index.name = "donor"
    with_eids = eid_map.merge(scores_df.reset_index(), on="donor", how="inner")
    with_eids[["eid"] + pc_cols].to_csv(out_pcs, sep="\t", index=False)
    log(f"  wrote {out_pcs}")

    _patch_parquet(cfg, with_eids[["eid"] + pc_cols], pc_cols)


def _add_public_phenotypes(cfg: Config, rows: pd.DataFrame) -> None:
    """Fill the fields GTEx does release openly, in place.

    Sex is in SubjectPhenotypes and is genuinely observed; GTEx codes it 1=male,
    2=female, matching the convention already used for 1KGP rows. Ethnic
    background is *inferred* from the joint genotype PC space (see
    playground/pcca/joint.py) and is only filled when that analysis has been run;
    it is recorded in the docs as an inference so no consumer mistakes it for a
    reported value.
    """
    # Reuse the pipeline's existing Data-Coding 1001 mapping rather than a
    # second, divergent one.
    from tbooo.pipeline.phenotypes import _ETHNIC_BG
    donors = pd.read_csv(cfg.metadata_dir() / "eid_map_gtex.tsv", sep="\t")
    eid_to_donor = dict(zip(donors["eid"], donors["donor"]))
    donor_series = rows["eid"].map(eid_to_donor)

    phenotypes = cfg.gtex_raw_dir() / SUBJECT_PHENOTYPES_FILE
    if phenotypes.exists():
        pheno = pd.read_csv(phenotypes, sep="\t").drop_duplicates("SUBJID")
        sex = donor_series.map(dict(zip(pheno["SUBJID"], pheno["SEX"])))
        rows["p31"] = sex.astype("Int64")
        log(f"  sex filled for {int(sex.notna().sum())} GTEx donors")

    inferred_path = cfg.metadata_dir() / "gtex_inferred_ancestry.tsv"
    if inferred_path.exists():
        inferred = pd.read_csv(inferred_path, sep="\t")
        group = donor_series.map(dict(zip(inferred["sampleID"],
                                          inferred["inferred_group"])))
        rows["p21000_i0"] = group.map(_ETHNIC_BG).astype("Int64")
        log(f"  inferred ethnic background filled for "
            f"{int(group.notna().sum())} GTEx donors")
    else:
        log("  ethnic background left null (run the joint PC analysis to infer it)")


def _patch_parquet(cfg: Config, pc_df: pd.DataFrame, pc_cols: list[str]) -> None:
    """Append GTEx donors as new participant rows.

    Existing GTEx rows are removed first so the stage is idempotent; without that,
    re-running would duplicate every donor.
    """
    parquet_path = cfg.showcase_dir() / "participant.parquet"
    if not parquet_path.exists():
        log("  participant.parquet not found — PC scores saved to metadata only.")
        return

    existing = pq.read_table(str(parquet_path)).to_pandas()
    existing = existing[existing["eid"] < EID_BASE].copy()
    for column in pc_cols:
        if column not in existing.columns:
            existing[column] = np.nan

    # GTEx rows leave most fields null, and a plain int64 column cannot hold NA.
    # Promoting to pandas' nullable integer keeps sex/batch codes as integers
    # instead of silently becoming floats when the new rows are appended.
    for column, dtype in existing.dtypes.items():
        if column != "eid" and dtype.kind in "iu":
            existing[column] = existing[column].astype("Int64")

    new_rows = pd.DataFrame({"eid": pc_df["eid"].to_numpy()})
    for column in pc_cols:
        new_rows[column] = pc_df[column].to_numpy()
    _add_public_phenotypes(cfg, new_rows)

    # Every other field — including the genetic-availability flags — is left
    # absent so concat fills it with null: GTEx genotypes are not public.
    merged = pd.concat([existing, new_rows], ignore_index=True)[existing.columns]
    pq.write_table(pa.Table.from_pandas(merged), str(parquet_path),
                   compression="snappy")
    log(f"  patched {parquet_path.name}: {len(existing)} existing "
        f"+ {len(new_rows)} GTEx rows = {len(merged)}")
