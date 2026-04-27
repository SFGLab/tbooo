"""Build UKB-mirrored WGS data (Fields 23149, 23151, 23370) and SGDP per-chromosome pVCFs.

Individual CRAMs (Field 23149):
  - Rename/symlink 1KGP NYGC 30x CRAMs to <EID>_23149_0_0.cram
  - Place in EID-prefix subfolders under data/Bulk/Whole genome sequences/

Individual gVCFs (Field 23151):
  - Symlink SGDP per-sample VCFs to <EID>_23151_0_0.g.vcf.gz
  - Same EID-prefix subdirectory layout as CRAMs

1KGP cohort pVCF (Field 23370):
  - Reheader NYGC 30x per-chromosome VCFs, replacing sample IDs with EIDs
  - Output: data/Bulk/Whole genome sequences/ukb23370_c{chrom}_b0_v1.pvcf.gz

SGDP per-chromosome pVCF:
  - Merge per-sample SGDP VCFs (data/raw/sgdp/vcf/) by chromosome
  - Reheader with EIDs
  - Output: data/raw/sgdp/pvcf/sgdp_c{chrom}.pvcf.gz
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from tbooo.config import Config
from tbooo.utils import ensure_dirs, eid_prefix_dir, log, run


def rename_crams(cfg: Config) -> None:
    """Symlink 1KGP NYGC 30x CRAMs into the UKB-mirrored individual-file structure."""
    ensure_dirs(cfg.wgs_dir())
    kg_map = _load_eid_map(cfg, "eid_map_1kg.tsv")
    _rename_kg_crams(cfg, kg_map)
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


def symlink_sgdp_gvcfs(cfg: Config) -> None:
    """Symlink SGDP per-sample VCFs into the UKB-mirrored gVCF structure (Field 23151).

    Each downloaded <accession>.vcf.gz becomes:
        data/Bulk/Whole genome sequences/<EID_prefix>/<EID>_23151_0_0.g.vcf.gz
    """
    vcf_dir = cfg.sgdp_vcf_dir()
    if not vcf_dir.exists() or not any(vcf_dir.glob("*.vcf.gz")):
        log("  No SGDP VCF files found; skipping. Run `tbooo download sgdp` first.")
        return

    eid_map = _load_eid_map(cfg, "eid_map_sgdp.tsv")
    if eid_map.empty:
        log("  WARNING: SGDP EID map not found; skipping. Run `tbooo map eids` first.")
        return

    log(f"[sgdp gVCF] symlinking {len(eid_map)} per-sample VCFs…")
    linked = 0
    for _, row in eid_map.iterrows():
        eid = int(row["eid"])
        ena_acc = row["ena_accession"]

        matches = list(vcf_dir.glob(f"*{ena_acc}*.vcf.gz"))
        if not matches:
            continue
        src_vcf = matches[0].resolve()
        src_tbi = Path(str(src_vcf) + ".tbi")

        prefix = eid_prefix_dir(eid)
        dest_dir = cfg.wgs_dir() / prefix
        ensure_dirs(dest_dir)

        dest_vcf = dest_dir / f"{eid}_23151_0_0.g.vcf.gz"
        dest_tbi = dest_dir / f"{eid}_23151_0_0.g.vcf.gz.tbi"

        if not dest_vcf.exists() and src_vcf.exists():
            dest_vcf.symlink_to(src_vcf)
            linked += 1
        if not dest_tbi.exists() and src_tbi.exists():
            dest_tbi.symlink_to(src_tbi)

    log(f"  linked {linked} gVCF(s) → {cfg.wgs_dir()}/")
    log("SGDP gVCF symlinking complete.")


def build_sgdp_pvcf(cfg: Config, chroms: list[str]) -> None:
    """Merge per-sample SGDP VCFs into per-chromosome multi-sample pVCFs with EID renaming."""
    vcf_dir = cfg.sgdp_vcf_dir()
    sample_vcfs = sorted(vcf_dir.glob("*.vcf.gz")) if vcf_dir.exists() else []
    if not sample_vcfs:
        log("  No SGDP VCF files found; skipping. Run `tbooo download sgdp` first.")
        return

    rename_file = cfg.metadata_dir() / "vcf_sample_rename_sgdp.txt"
    if not rename_file.exists():
        raise FileNotFoundError(
            f"SGDP rename file not found: {rename_file}\nRun `tbooo map eids` first."
        )

    pvcf_dir = cfg.sgdp_raw_dir() / "pvcf"
    ensure_dirs(pvcf_dir)
    ensure_dirs(cfg.tmp_dir)
    log(f"[sgdp pVCF] {len(sample_vcfs)} sample VCF(s) found")

    for chrom in chroms:
        out = cfg.sgdp_pvcf(chrom)
        if out.exists():
            log(f"  skip (exists): {out.name}")
            continue

        log(f"  chromosome {chrom}")
        tmp_per_sample: list[Path] = []

        for vcf in sample_vcfs:
            tmp = cfg.tmp_dir / f"sgdp_{vcf.stem}_chr{chrom}.vcf.gz"
            run([
                cfg.tools.bcftools, "view",
                "--regions", chrom,
                "--output-type", "z",
                "--output", str(tmp),
                str(vcf),
            ])
            run([cfg.tools.bcftools, "index", "--tbi", str(tmp)])
            tmp_per_sample.append(tmp)

        if len(tmp_per_sample) == 1:
            tmp_merged = tmp_per_sample[0]
        else:
            tmp_merged = cfg.tmp_dir / f"sgdp_merged_chr{chrom}.vcf.gz"
            run([
                cfg.tools.bcftools, "merge",
                "--output-type", "z",
                "--output", str(tmp_merged),
            ] + [str(v) for v in tmp_per_sample])
            run([cfg.tools.bcftools, "index", "--tbi", str(tmp_merged)])

        run([
            cfg.tools.bcftools, "reheader",
            "--samples", str(rename_file),
            "--output", str(out),
            str(tmp_merged),
        ])
        run([cfg.tools.bcftools, "index", "--tbi", str(out)])

        for tmp in tmp_per_sample:
            tmp.unlink(missing_ok=True)
            Path(str(tmp) + ".tbi").unlink(missing_ok=True)
        if tmp_merged not in tmp_per_sample:
            tmp_merged.unlink(missing_ok=True)
            Path(str(tmp_merged) + ".tbi").unlink(missing_ok=True)

        log(f"  done → {out.name}")

    log("SGDP pVCF build complete.")


# ── Internals ─────────────────────────────────────────────────────────────────

def _rename_kg_crams(cfg: Config, kg_map: pd.DataFrame) -> None:
    if kg_map.empty:
        return
    raw = cfg.kg_raw_dir()
    log(f"  linking 1KGP NYGC CRAMs ({len(kg_map)} samples)…")
    for _, row in kg_map.iterrows():
        eid = int(row["eid"])
        sample_id = row["sample_id"]
        crms = list(raw.glob(f"*{sample_id}*.cram"))
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
