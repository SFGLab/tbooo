"""Build UKB-mirrored sample QC and relatedness files.

ukb_sqc_v2.txt:
    One row per sample; key flags: excess.relatives, used.in.pca.calculation,
    in.white.british.ancestry.subset, het.missing.outliers.

ukb_rel.txt:
    Pairwise kinship estimates for pairs with kinship > 0.0442 (3rd-degree or closer).
    Columns: EID1, EID2, HetHet, IBS0, Kinship.

Steps:
  1. Merge per-chromosome array PLINK files into a genome-wide set
  2. Run PLINK2 --het for heterozygosity outlier detection
  3. Run KING --kinship for pairwise kinship estimation
  4. Parse KING output, replace original IDs with EIDs, filter to kinship > 0.0442
  5. Write ukb_sqc_v2.txt and ukb_rel.txt
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pandas as pd

from tbooo.config import Config
from tbooo.integrity import plink_ok, table_ok
from tbooo.utils import ensure_dirs, log, run

_KINSHIP_THRESHOLD = 0.0442  # 3rd-degree relatives and closer


def build_qc_files(cfg: Config) -> None:
    ensure_dirs(cfg.metadata_dir(), cfg.tmp_dir)

    sqc_out = cfg.metadata_dir() / "ukb_sqc_v2.txt"
    rel_out = cfg.metadata_dir() / "ukb_rel.txt"
    # ukb_rel.txt is header-only when no related pairs are found (still valid).
    if table_ok(sqc_out, min_lines=2) and table_ok(rel_out, min_lines=1):
        log(f"  skip (valid): {sqc_out.name} + {rel_out.name}")
        return

    eid_map = _load_combined_eid_map(cfg)
    merged_plink = _merge_array_plink(cfg)
    het_stats = _compute_het(cfg, merged_plink)
    king_kin = _run_king(cfg, merged_plink)

    _write_sqc(cfg, eid_map, het_stats)
    _write_rel(cfg, eid_map, king_kin)

    log("QC file generation complete.")


# ── Steps ─────────────────────────────────────────────────────────────────────

def _merge_array_plink(cfg: Config) -> Path:
    """Merge per-chromosome array PLINK files into a single genome-wide set."""
    merged = cfg.tmp_dir / "array_merged_all_chrs"
    if plink_ok(merged):
        log("  skip (valid): merged PLINK set")
        return merged

    # Build merge list
    merge_list = cfg.tmp_dir / "plink_merge_list.txt"
    stems: list[Path] = []
    for chrom in cfg.autosomes:
        stem = cfg.array_stem(str(chrom))
        if plink_ok(stem):
            stems.append(stem)

    if not stems:
        raise FileNotFoundError(
            "No array PLINK files found. Run `tbooo map array` first."
        )

    if len(stems) == 1:
        return stems[0]

    # plink2 merge-list: first file is the target, rest go into merge-list
    with open(merge_list, "w") as f:
        for s in stems[1:]:
            f.write(f"{s}.bed {s}.bim {s}.fam\n")

    run([
        cfg.tools.plink2,
        "--bfile", str(stems[0]),
        "--pmerge-list", str(merge_list),
        "--make-bed",
        "--out", str(merged),
    ])
    log(f"  merged {len(stems)} chromosomes → {merged}.bed")
    return merged


def _compute_het(cfg: Config, plink_stem: Path) -> pd.DataFrame:
    """Compute per-sample heterozygosity rate using PLINK2."""
    het_prefix = cfg.tmp_dir / "het_stats"
    het_file = het_prefix.parent / (het_prefix.name + ".het")

    if not table_ok(het_file, min_lines=2):
        run([
            cfg.tools.plink2,
            "--bfile", str(plink_stem),
            "--het",
            "--out", str(het_prefix),
        ])

    if not het_file.exists():
        log("  WARNING: PLINK --het output not found; het outlier flag set to 0 for all")
        return pd.DataFrame(columns=["IID", "het_outlier"])

    het = pd.read_csv(het_file, sep=r"\s+")
    # Columns: #FID IID O(HOM) E(HOM) OBS_CT F
    # Outlier: |F| > mean ± 5 SD
    f_col = "F"
    mean_f = het[f_col].mean()
    std_f = het[f_col].std()
    het["het_outlier"] = (
        (het[f_col] < mean_f - 5 * std_f) | (het[f_col] > mean_f + 5 * std_f)
    ).astype(int)
    return het[["IID", "het_outlier"]]


def _run_king(cfg: Config, plink_stem: Path) -> pd.DataFrame:
    """Run KING kinship estimation; return pairs above 3rd-degree threshold."""
    king_prefix = cfg.tmp_dir / "king_output"
    kin0_file = Path(str(king_prefix) + ".kin0")

    if not table_ok(kin0_file, min_lines=1):
        try:
            run([
                cfg.tools.king,
                "--bfile", str(plink_stem),
                "--kinship",
                "--prefix", str(king_prefix),
            ])
        except subprocess.CalledProcessError:
            log("  WARNING: KING failed; relatedness file will be empty")
            return pd.DataFrame(columns=["ID1", "ID2", "HetHet", "IBS0", "Kinship"])

    if not kin0_file.exists():
        return pd.DataFrame(columns=["ID1", "ID2", "HetHet", "IBS0", "Kinship"])

    kin = pd.read_csv(kin0_file, sep=r"\s+")
    # KING .kin0 columns: FID1 ID1 FID2 ID2 N_SNP HetHet IBS0 Kinship
    required = {"ID1", "ID2", "HetHet", "IBS0", "Kinship"}
    if not required.issubset(set(kin.columns)):
        log(f"  WARNING: unexpected KING output columns: {list(kin.columns)}")
        return pd.DataFrame(columns=["ID1", "ID2", "HetHet", "IBS0", "Kinship"])

    return kin[kin["Kinship"] > _KINSHIP_THRESHOLD][["ID1", "ID2", "HetHet", "IBS0", "Kinship"]]


# ── Output writers ────────────────────────────────────────────────────────────

def _write_sqc(cfg: Config, eid_map: pd.DataFrame, het_stats: pd.DataFrame) -> None:
    out = cfg.metadata_dir() / "ukb_sqc_v2.txt"

    # Build lookup from sample_id → eid row
    sqc = eid_map.copy()

    # used.in.pca.calculation: 1 for 1KGP unrelated samples, 0 for SGDP
    sqc["used.in.pca.calculation"] = (sqc["source"] == "1kg").astype(int)

    # in.white.british.ancestry.subset: 1 for GBR EUR samples
    sqc["in.white.british.ancestry.subset"] = (
        (sqc.get("pop", "") == "GBR") & (sqc.get("super_pop", "") == "EUR")
    ).astype(int)

    # excess.relatives: flag samples with known relatives in the dataset
    relative_eids = _find_excess_relatives(cfg)
    sqc["excess.relatives"] = sqc["eid"].isin(relative_eids).astype(int)

    # het.missing.outliers
    if not het_stats.empty and "IID" in het_stats.columns:
        # IID in het_stats is the EID (after rewriting FAM)
        outlier_set = set(het_stats[het_stats["het_outlier"] == 1]["IID"].astype(int))
        sqc["het.missing.outliers"] = sqc["eid"].isin(outlier_set).astype(int)
    else:
        sqc["het.missing.outliers"] = 0

    # putative sex chromosome aneuploidy (not computed; set to 0)
    sqc["putative.sex.chromosome.aneuploidy"] = 0

    out_cols = [
        "eid", "used.in.pca.calculation", "in.white.british.ancestry.subset",
        "excess.relatives", "het.missing.outliers", "putative.sex.chromosome.aneuploidy",
    ]
    sqc[out_cols].to_csv(out, sep="\t", index=False)
    log(f"  wrote {out} ({len(sqc)} samples)")


def _write_rel(cfg: Config, eid_map: pd.DataFrame, king_kin: pd.DataFrame) -> None:
    out = cfg.metadata_dir() / "ukb_rel.txt"

    if king_kin.empty:
        out.write_text("EID1\tEID2\tHetHet\tIBS0\tKinship\n")
        log(f"  wrote {out} (empty — no related pairs found)")
        return

    # Replace original sample IDs with EIDs
    # KING IDs will already be EIDs if plink FAM was rewritten; this is a safety pass
    id_to_eid = dict(zip(eid_map["sample_id"].astype(str), eid_map["eid"].astype(str)))

    def maybe_remap(x):
        return id_to_eid.get(str(x), str(x))

    rel = king_kin.copy()
    rel["EID1"] = rel["ID1"].map(maybe_remap)
    rel["EID2"] = rel["ID2"].map(maybe_remap)
    rel[["EID1", "EID2", "HetHet", "IBS0", "Kinship"]].to_csv(out, sep="\t", index=False)
    log(f"  wrote {out} ({len(rel)} related pairs above kinship threshold {_KINSHIP_THRESHOLD})")


def _find_excess_relatives(cfg: Config) -> set[int]:
    """
    Return set of EIDs with more relatives than expected.
    Uses 1KGP pedigree file to identify samples with known first-degree relatives.
    Samples with >10 relatives are flagged (arbitrary threshold matching UKB convention).
    """
    ped_path = cfg.kg_raw_dir() / "20130606_g1k.ped"
    eid_map_path = cfg.metadata_dir() / "eid_map_1kg.tsv"

    if not ped_path.exists() or not eid_map_path.exists():
        return set()

    ped = pd.read_csv(ped_path, sep="\t")
    eid_map = pd.read_csv(eid_map_path, sep="\t")
    id_to_eid = dict(zip(eid_map["sample_id"], eid_map["eid"]))

    # Count how many relatives each sample has in the pedigree
    from collections import Counter
    relative_counts: Counter = Counter()
    for col in ("Paternal ID", "Maternal ID"):
        if col in ped.columns:
            for sid in ped[col].dropna():
                if sid != "0" and sid in id_to_eid:
                    relative_counts[id_to_eid[sid]] += 1

    return {eid for eid, count in relative_counts.items() if count > 10}


def _load_combined_eid_map(cfg: Config) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for fname in ("eid_map_1kg.tsv", "eid_map_sgdp.tsv"):
        p = cfg.metadata_dir() / fname
        if p.exists():
            df = pd.read_csv(p, sep="\t")
            if "sample_id" not in df.columns and "ena_accession" in df.columns:
                df = df.rename(columns={"ena_accession": "sample_id"})
            frames.append(df)
    if not frames:
        raise FileNotFoundError("No EID maps found. Run `tbooo map eids` first.")
    combined = pd.concat(frames, ignore_index=True)
    for col in ("pop", "super_pop", "region", "source"):
        if col not in combined.columns:
            combined[col] = ""
    return combined
