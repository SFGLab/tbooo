"""Compute gene expression PCA from GEUVADIS and merge scores into participant.parquet.

Steps:
    1. Read GD462 RPKM matrix (genes × samples).
    2. Filter genes: median RPKM >= MIN_MEDIAN_RPKM across all samples.
    3. Log2-transform: log2(RPKM + 0.1).
    4. PCA(n_components=N_PCS) on the (samples × genes) matrix.
    5. Join scores to eid_map_1kg on sample_id.
    6. Write geuvadis_expression_pcs.tsv and geuvadis_pca_variance.tsv.
    7. Merge PC columns into participant.parquet (samples not in GEUVADIS get null).

Outputs:
    data/metadata/geuvadis_expression_pcs.tsv   columns: eid, geuvadis_pc1..N
    data/metadata/geuvadis_pca_variance.tsv      columns: pc, variance_explained, cumulative_variance
    data/Showcase/participant.parquet            (updated in-place with geuvadis_pc* columns)
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from sklearn.decomposition import PCA

from tbooo.config import Config
from tbooo.utils import ensure_dirs, log

N_PCS = 10
_MIN_MEDIAN_RPKM = 0.1
_RPKM_FILENAME = "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"
# 1KGP sample IDs: 2 uppercase letters + 5 digits (HG00096, NA12878, …)
_SAMPLE_ID_LEN = 7


def build_expression_pcs(cfg: Config) -> None:
    rpkm_path = cfg.geuvadis_raw_dir() / _RPKM_FILENAME
    if not rpkm_path.exists():
        raise FileNotFoundError(
            f"GEUVADIS RPKM matrix not found: {rpkm_path}\n"
            "Run `tbooo download geuvadis` first."
        )

    eid_map_path = cfg.metadata_dir() / "eid_map_1kg.tsv"
    if not eid_map_path.exists():
        raise FileNotFoundError(
            f"EID map not found: {eid_map_path}\n"
            "Run `tbooo map eids` first."
        )

    ensure_dirs(cfg.metadata_dir())

    log("Loading GEUVADIS RPKM matrix…")
    rpkm = _load_rpkm(rpkm_path)
    log(f"  {len(rpkm)} genes × {len(rpkm.columns)} samples before filtering")

    rpkm = rpkm.loc[rpkm.median(axis=1) >= _MIN_MEDIAN_RPKM]
    log(f"  {len(rpkm)} genes after median RPKM ≥ {_MIN_MEDIAN_RPKM} filter")

    vals = np.clip(np.nan_to_num(rpkm.values.T.astype(float), nan=0.0, posinf=0.0, neginf=0.0), 0.0, None)
    X = np.log2(vals + 0.1)  # (n_samples, n_genes)
    sample_ids = list(rpkm.columns)

    log(f"Running PCA (n_components={N_PCS})…")
    pca = PCA(n_components=N_PCS)
    scores = pca.fit_transform(X)

    pc_cols = [f"geuvadis_pc{i}" for i in range(1, N_PCS + 1)]
    scores_df = pd.DataFrame(scores, index=sample_ids, columns=pc_cols)
    scores_df.index.name = "sample_id"

    eid_map = pd.read_csv(eid_map_path, sep="\t")[["eid", "sample_id"]]
    pc_with_eids = eid_map.merge(scores_df.reset_index(), on="sample_id", how="inner")

    out_pcs = cfg.metadata_dir() / "geuvadis_expression_pcs.tsv"
    pc_with_eids[["eid"] + pc_cols].to_csv(out_pcs, sep="\t", index=False)
    log(f"  wrote {out_pcs} ({len(pc_with_eids)} samples matched to EIDs)")

    var_df = pd.DataFrame({
        "pc": pc_cols,
        "variance_explained": pca.explained_variance_ratio_,
        "cumulative_variance": np.cumsum(pca.explained_variance_ratio_),
    })
    out_var = cfg.metadata_dir() / "geuvadis_pca_variance.tsv"
    var_df.to_csv(out_var, sep="\t", index=False)
    log(f"  wrote {out_var}")
    log(f"  variance explained by {N_PCS} PCs: {pca.explained_variance_ratio_.sum():.1%}")

    _patch_parquet(cfg, pc_with_eids[["eid"] + pc_cols], pc_cols)


def _load_rpkm(path) -> pd.DataFrame:
    """Return DataFrame with genes as rows, 1KGP sample IDs as columns."""
    df = pd.read_csv(path, sep="\t", compression="gzip", index_col=0)

    # Drop metadata columns (Chr, Start, End, Strand, Length, …).
    # Keep only columns whose names look like 1KGP sample IDs.
    sample_cols = [c for c in df.columns if _is_sample_id(c)]
    if not sample_cols:
        raise ValueError(
            f"No 1KGP-style sample ID columns found in {path}. "
            "Expected column names like 'HG00096' or 'NA12878'."
        )
    return df[sample_cols]


def _is_sample_id(col: str) -> bool:
    return len(col) == _SAMPLE_ID_LEN and col[:2].isalpha() and col[:2].isupper() and col[2:].isdigit()


def _patch_parquet(cfg: Config, pc_df: pd.DataFrame, pc_cols: list[str]) -> None:
    """Merge PC columns into participant.parquet, replacing any existing geuvadis_pc* columns."""
    parquet_path = cfg.showcase_dir() / "participant.parquet"
    if not parquet_path.exists():
        log(
            "  participant.parquet not found — PC scores saved to metadata only.\n"
            "  Run `tbooo map phenotypes` then re-run `tbooo map geuvadis` to embed scores."
        )
        return

    existing = pq.read_table(str(parquet_path)).to_pandas()
    existing = existing.drop(columns=[c for c in existing.columns if c.startswith("geuvadis_pc")], errors="ignore")
    merged = existing.merge(pc_df, on="eid", how="left")

    pq.write_table(pa.Table.from_pandas(merged), str(parquet_path), compression="snappy")
    n_matched = pc_df["eid"].isin(existing["eid"]).sum()
    log(f"  patched {parquet_path.name}: {n_matched}/{len(existing)} rows have expression PCs")
