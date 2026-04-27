"""Build UKB-mirrored BGEN imputed genotype files (Field 22828).

For each chromosome:
  1. Normalize 1KGP Phase 3 VCF (decompose multi-allelic, left-align)
  2. Replace original sample IDs with synthetic EIDs (bcftools reheader)
  3. Convert VCF → BGEN v1.2 (8-bit, zlib, dosages) via qctool
  4. Index the BGEN with bgenix (.bgi)
  5. Write a .sample file in GEN2 format

Outputs:
    data/Bulk/Imputed/ukb22828_c{chrom}_b0_v3.{bgen,bgen.bgi,sample}
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from tbooo.config import Config
from tbooo.utils import ensure_dirs, log, run

_GRCH37_FASTA = "GRCh37/human_g1k_v37.fasta.gz"


def run_imputed_pipeline(cfg: Config, chroms: list[str]) -> None:
    ensure_dirs(cfg.imputed_dir())
    eid_map = _load_eid_map(cfg)
    rename_file = cfg.metadata_dir() / "vcf_sample_rename_1kg.txt"

    # Write sample file (shared across all chromosomes)
    sample_path = cfg.imputed_dir() / "ukb22828_samples.sample"
    _write_sample_file(eid_map, sample_path)

    for chrom in chroms:
        log(f"[imputed] chromosome {chrom}")
        vcf = cfg.phase3_vcf(chrom)
        if not vcf.exists():
            log(f"  SKIP: Phase 3 VCF not found: {vcf}")
            continue
        _process_chrom(cfg, chrom, vcf, rename_file, sample_path)

    log("Imputed pipeline complete.")


def _process_chrom(
    cfg: Config,
    chrom: str,
    vcf: Path,
    rename_file: Path,
    sample_path: Path,
) -> None:
    stem = cfg.imputed_stem(chrom)
    tmp_norm = cfg.tmp_dir / f"imputed_chr{chrom}_norm.vcf.gz"
    tmp_reheadered = cfg.tmp_dir / f"imputed_chr{chrom}_eid.vcf.gz"
    ensure_dirs(cfg.tmp_dir)

    # Step 1: normalize (split multi-allelic, left-align indels)
    ref = cfg.reference_dir / _GRCH37_FASTA
    if ref.exists():
        norm_cmd = [
            cfg.tools.bcftools, "norm",
            "--multiallelics", "-",
            "--fasta-ref", str(ref),
            "--output-type", "z",
            "--output", str(tmp_norm),
            str(vcf),
        ]
    else:
        log(f"  WARNING: GRCh37 reference not found; skipping left-align")
        norm_cmd = [
            cfg.tools.bcftools, "norm",
            "--multiallelics", "-",
            "--output-type", "z",
            "--output", str(tmp_norm),
            str(vcf),
        ]
    run(norm_cmd)
    run([cfg.tools.bcftools, "index", "--tbi", str(tmp_norm)])

    # Step 2: replace sample IDs with EIDs
    if rename_file.exists():
        run([
            cfg.tools.bcftools, "reheader",
            "--samples", str(rename_file),
            "--output", str(tmp_reheadered),
            str(tmp_norm),
        ])
        run([cfg.tools.bcftools, "index", "--tbi", str(tmp_reheadered)])
        source_vcf = tmp_reheadered
    else:
        log("  WARNING: rename file not found; sample IDs not replaced")
        source_vcf = tmp_norm

    # Step 3: VCF → BGEN v1.2 via qctool
    bgen_path = Path(str(stem) + ".bgen")
    run([
        cfg.tools.qctool,
        "-g", str(source_vcf),
        "-filetype", "vcf",
        "-og", str(bgen_path),
        "-ofiletype", "bgen_v1.2",
        "-bgen-bits", "8",
        "-os", str(sample_path),
    ])

    # Step 4: index BGEN
    run([cfg.tools.bgenix, "-g", str(bgen_path), "-index"])

    # Step 5: symlink per-chromosome sample file
    chrom_sample = Path(str(stem) + ".sample")
    if not chrom_sample.exists():
        chrom_sample.symlink_to(sample_path.name)

    tmp_norm.unlink(missing_ok=True)
    tmp_reheadered.unlink(missing_ok=True)
    Path(str(tmp_norm) + ".tbi").unlink(missing_ok=True)
    Path(str(tmp_reheadered) + ".tbi").unlink(missing_ok=True)

    log(f"  done → {stem}.bgen/bgi/sample")


def _write_sample_file(eid_map: pd.DataFrame, path: Path) -> None:
    """Write a GEN2-format .sample file for BGEN."""
    if path.exists():
        log(f"  skip (exists): {path.name}")
        return
    lines = ["ID_1 ID_2 missing sex\n", "0 0 0 D\n"]
    for _, row in eid_map.iterrows():
        eid = row["eid"]
        sex = int(row.get("sex", 0))
        lines.append(f"{eid} {eid} 0 {sex}\n")
    path.write_text("".join(lines))
    log(f"  wrote sample file: {path} ({len(eid_map)} samples)")


def _load_eid_map(cfg: Config) -> pd.DataFrame:
    path = cfg.metadata_dir() / "eid_map_1kg.tsv"
    if not path.exists():
        raise FileNotFoundError(f"EID map not found: {path}\nRun `tbooo map eids` first.")
    return pd.read_csv(path, sep="\t")
