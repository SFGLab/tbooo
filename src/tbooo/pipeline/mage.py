"""Compute MAGE expression PCs and merge them into participant.parquet.

MAGE individuals are all 1000 Genomes samples, so they already carry TBOOO EIDs;
this stage adds columns to existing rows rather than introducing new participants.

Steps:
    1. Read the MAGE TMM expression matrix (genes x 731 samples).
    2. Drop the two individuals affected by the unresolved sample swap.
    3. Filter genes, log-transform, select highly variable genes, PCA.
    4. Join to eid_map_1kg on sample_id.
    5. Write mage_expression_pcs.tsv and merge mage_pc* into participant.parquet.

Outputs:
    data/metadata/mage_expression_pcs.tsv   columns: eid, mage_pc1..N
    data/Showcase/participant.parquet       (updated in place)
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from sklearn.decomposition import PCA

from tbooo.config import Config
from tbooo.integrity import table_ok
from tbooo.utils import ensure_dirs, log

N_PCS = 20
N_HVG = 2000
_EXPRESSION_FILE = "TMM.filtered.TSS.MAGE.v1.0.bed.gz"
_METADATA_FILE = "sample.metadata.MAGE.v1.0.txt"
_BED_ANNOTATION = ["#Chr", "start", "end", "ID"]

# See playground/pcca/mage.py: SRR19762653 is labelled NA11919 but carries
# HG00237's material, and HG00237 has its own library, so relabelling would
# duplicate it. Both are excluded until MAGE releases corrected data.
_SWAPPED = ("HG00237", "NA11919")


def build_expression_pcs(cfg: Config) -> None:
    raw_dir = cfg.mage_raw_dir()
    expression_path = raw_dir / _EXPRESSION_FILE
    if not expression_path.exists():
        raise FileNotFoundError(
            f"MAGE expression matrix not found: {expression_path}\n"
            "Run `tbooo download mage` first."
        )

    eid_map_path = cfg.metadata_dir() / "eid_map_1kg.tsv"
    if not eid_map_path.exists():
        raise FileNotFoundError(
            f"EID map not found: {eid_map_path}\nRun `tbooo map eids` first."
        )

    ensure_dirs(cfg.metadata_dir())
    out_pcs = cfg.metadata_dir() / "mage_expression_pcs.tsv"
    if table_ok(out_pcs, min_lines=2):
        log(f"  skip (valid): {out_pcs.name}")
        return

    log("Loading MAGE expression matrix…")
    df = pd.read_csv(expression_path, sep="\t")
    expr = df.drop(columns=_BED_ANNOTATION)
    log(f"  {len(df)} genes × {expr.shape[1]} samples")

    dropped = [s for s in _SWAPPED if s in expr.columns]
    expr = expr.drop(columns=dropped)
    log(f"  dropped {len(dropped)} sample(s) affected by the unresolved swap: {dropped}")

    values = np.nan_to_num(expr.to_numpy(dtype=float).T, nan=0.0,
                           posinf=0.0, neginf=0.0)
    X = np.log2(np.clip(values, 0.0, None) + 0.1)

    gene_vars = X.var(axis=0)
    X = X[:, np.argsort(gene_vars)[-N_HVG:]]
    log(f"  {X.shape[1]} highly variable genes selected")

    log(f"Running PCA (n_components={N_PCS})…")
    pca = PCA(n_components=N_PCS)
    scores = pca.fit_transform(X)

    pc_cols = [f"mage_pc{i}" for i in range(1, N_PCS + 1)]
    scores_df = pd.DataFrame(scores, index=list(expr.columns), columns=pc_cols)
    scores_df.index.name = "sample_id"

    eid_map = pd.read_csv(eid_map_path, sep="\t")[["eid", "sample_id"]]
    with_eids = eid_map.merge(scores_df.reset_index(), on="sample_id", how="inner")
    with_eids[["eid"] + pc_cols].to_csv(out_pcs, sep="\t", index=False)
    log(f"  wrote {out_pcs} ({len(with_eids)} samples matched to EIDs)")
    log(f"  variance explained by {N_PCS} PCs: {pca.explained_variance_ratio_.sum():.1%}")

    _patch_parquet(cfg, with_eids[["eid"] + pc_cols], pc_cols)


def _patch_parquet(cfg: Config, pc_df: pd.DataFrame, pc_cols: list[str]) -> None:
    parquet_path = cfg.showcase_dir() / "participant.parquet"
    if not parquet_path.exists():
        log("  participant.parquet not found — PC scores saved to metadata only.")
        return

    existing = pq.read_table(str(parquet_path)).to_pandas()
    existing = existing.drop(columns=[c for c in existing.columns
                                      if c.startswith("mage_pc")], errors="ignore")
    merged = existing.merge(pc_df, on="eid", how="left")
    pq.write_table(pa.Table.from_pandas(merged), str(parquet_path),
                   compression="snappy")
    matched = int(pc_df["eid"].isin(existing["eid"]).sum())
    log(f"  patched {parquet_path.name}: {matched}/{len(existing)} rows have MAGE PCs")
