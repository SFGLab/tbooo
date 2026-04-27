"""Build UKB-mirrored WGS data (Fields 23149, 23370).

Individual CRAMs:
  - Rename/symlink 1KGP NYGC 30x and SGDP CRAMs to <EID>_23149_0_0.cram
  - Place in EID-prefix subfolders under data/Bulk/Whole genome sequences/

Cohort pVCF:
  - Reheader NYGC 30x per-chromosome VCFs, replacing sample IDs with EIDs
  - Output: data/Bulk/Whole genome sequences/ukb23370_c{chrom}_b0_v1.pvcf.gz
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from tbooo.config import Config
from tbooo.utils import ensure_dirs, eid_prefix_dir, log, run


def rename_crams(cfg: Config) -> None:
    """Symlink all source CRAMs into the UKB-mirrored individual-file structure."""
    ensure_dirs(cfg.wgs_dir())

    kg_map = _load_eid_map(cfg, "eid_map_1kg.tsv")
    sgdp_map = _load_eid_map(cfg, "eid_map_sgdp.tsv")

    _rename_kg_crams(cfg, kg_map)
    _rename_sgdp_crams(cfg, sgdp_map)

    log("CRAM renaming complete.")


def build_pvcf(cfg: Config, chroms: list[str]) -> None:
    """Reheader NYGC 30x VCFs with EIDs and write as cohort pVCF."""
    ensure_dirs(cfg.wgs_dir())
    rename_file = cfg.metadata_dir() / "vcf_sample_rename_1kg.txt"

    if not rename_file.exists():
        raise FileNotFoundError(
            f"Sample rename file not found: {rename_file}\nRun `tbooo map eids` first."
        )

    for chrom in chroms:
        log(f"[wgs pVCF] chromosome {chrom}")
        src = cfg.nygc_vcf(chrom)
        if not src.exists():
            log(f"  SKIP: NYGC VCF not found: {src}")
            continue

        out = cfg.wgs_pvcf(chrom)
        if out.exists():
            log(f"  skip (exists): {out.name}")
            continue

        run([
            cfg.tools.bcftools, "reheader",
            "--samples", str(rename_file),
            "--output", str(out),
            str(src),
        ])
        run([cfg.tools.bcftools, "index", "--tbi", str(out)])
        log(f"  done → {out.name}")

    log("pVCF build complete.")


# ── Internals ─────────────────────────────────────────────────────────────────

def _rename_kg_crams(cfg: Config, kg_map: pd.DataFrame) -> None:
    if kg_map.empty:
        return
    raw = cfg.kg_raw_dir()
    log(f"  linking 1KGP NYGC CRAMs ({len(kg_map)} samples)…")
    for _, row in kg_map.iterrows():
        eid = int(row["eid"])
        sample_id = row["sample_id"]
        # NYGC CRAM naming varies; search for any CRAM containing the sample ID
        crms = list(raw.glob(f"*{sample_id}*.cram"))
        if not crms:
            continue
        src_cram = crms[0]
        src_crai = Path(str(src_cram) + ".crai")
        _symlink_cram(cfg, eid, src_cram, src_crai)


def _rename_sgdp_crams(cfg: Config, sgdp_map: pd.DataFrame) -> None:
    if sgdp_map.empty:
        return
    raw = cfg.sgdp_raw_dir()
    log(f"  linking SGDP CRAMs ({len(sgdp_map)} samples)…")
    for _, row in sgdp_map.iterrows():
        eid = int(row["eid"])
        acc = row["ena_accession"]
        crms = list(raw.glob(f"{acc}*.cram"))
        if not crms:
            continue
        src_cram = crms[0]
        src_crai = Path(str(src_cram) + ".crai")
        _symlink_cram(cfg, eid, src_cram, src_crai)


def _symlink_cram(cfg: Config, eid: int, src_cram: Path, src_crai: Path) -> None:
    prefix = eid_prefix_dir(eid)
    dest_dir = cfg.wgs_dir() / prefix
    ensure_dirs(dest_dir)

    dest_cram = dest_dir / f"{eid}_23149_0_0.cram"
    dest_crai = dest_dir / f"{eid}_23149_0_0.cram.crai"

    if not dest_cram.exists() and src_cram.exists():
        dest_cram.symlink_to(src_cram.resolve())
    if not dest_crai.exists() and src_crai.exists():
        dest_crai.symlink_to(src_crai.resolve())


def _load_eid_map(cfg: Config, filename: str) -> pd.DataFrame:
    path = cfg.metadata_dir() / filename
    if not path.exists():
        log(f"  WARNING: {filename} not found; skipping")
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")
