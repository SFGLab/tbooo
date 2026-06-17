"""Build UKB-mirrored WES files (Field 23157) by intersecting WGS VCFs with exome BED.

For each chromosome:
  1. Intersect NYGC 30x VCF with IDT xGen exome capture BED (GRCh38)
  2. Replace sample IDs with synthetic EIDs
  3. Decompose multi-allelic, normalize
  4. Convert to PLINK BED/BIM/FAM
  5. Convert to BGEN v1.2

Outputs:
    data/Bulk/Exome sequences/.../ukb23157_c{chrom}_b0_v1.{bed,bim,fam}
    data/Bulk/Exome sequences/.../ukb23157_c{chrom}_b0_v1.{bgen,bgen.bgi,sample}
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from tbooo.config import Config
from tbooo.integrity import bgen_ok, ensure, plink_ok, remove, vcf_ok
from tbooo.utils import ensure_dirs, log, run

_GRCH38_FASTA = "GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa"


def run_wes_pipeline(cfg: Config, chroms: list[str]) -> None:
    exome_bed = Path(cfg.exome_bed)
    if not exome_bed.exists():
        raise FileNotFoundError(
            f"Exome capture BED not found: {exome_bed}\n"
            "Run `tbooo download reference` to fetch it."
        )
    ensure_dirs(cfg.wes_dir(), cfg.wes_bgen_dir())

    rename_file = cfg.metadata_dir() / "vcf_sample_rename_1kg.txt"
    eid_map = _load_eid_map(cfg)
    sample_path = cfg.wes_bgen_dir() / "ukb23157_samples.sample"
    _write_sample_file(eid_map, sample_path)

    for chrom in chroms:
        log(f"[wes] chromosome {chrom}")
        src = cfg.nygc_vcf(chrom)
        if not src.exists():
            log(f"  SKIP: NYGC VCF not found: {src}")
            continue
        _process_chrom(cfg, chrom, src, exome_bed, rename_file, sample_path)

    log("WES pipeline complete.")


def _process_chrom(
    cfg: Config,
    chrom: str,
    src_vcf: Path,
    exome_bed: Path,
    rename_file: Path,
    sample_path: Path,
) -> None:
    plink_stem = cfg.wes_stem(chrom)
    bgen_stem = cfg.wes_bgen_stem(chrom)

    def _purge() -> None:
        remove(Path(str(plink_stem) + ".bed"), Path(str(plink_stem) + ".bim"),
               Path(str(plink_stem) + ".fam"), Path(str(plink_stem) + ".log"),
               Path(str(bgen_stem) + ".bgen"), Path(str(bgen_stem) + ".bgen.bgi"))

    ensure(
        f"{plink_stem.name} (wes chr{chrom})",
        check=lambda: plink_ok(plink_stem) and bgen_ok(bgen_stem),
        purge=_purge,
        build=lambda: _build_chrom(cfg, chrom, src_vcf, exome_bed,
                                   rename_file, sample_path, plink_stem, bgen_stem),
    )


def _build_chrom(
    cfg: Config,
    chrom: str,
    src_vcf: Path,
    exome_bed: Path,
    rename_file: Path,
    sample_path: Path,
    plink_stem: Path,
    bgen_stem: Path,
) -> None:
    ensure_dirs(cfg.tmp_dir)
    tmp_intersected = cfg.tmp_dir / f"wes_chr{chrom}_intersected.vcf.gz"
    tmp_reheadered = cfg.tmp_dir / f"wes_chr{chrom}_eid.vcf.gz"

    # Step 1: intersect with exome BED
    _intersect_exome(cfg, chrom, src_vcf, exome_bed, tmp_intersected)

    # Step 2: replace sample IDs
    if rename_file.exists():
        run([
            cfg.tools.bcftools, "reheader",
            "--samples", str(rename_file),
            "--output", str(tmp_reheadered),
            str(tmp_intersected),
        ])
        run([cfg.tools.bcftools, "index", "--tbi", str(tmp_reheadered)])
        source = tmp_reheadered
    else:
        source = tmp_intersected

    # Step 3: normalize
    ref = cfg.reference_dir / _GRCH38_FASTA
    tmp_norm = cfg.tmp_dir / f"wes_chr{chrom}_norm.vcf.gz"
    if ref.exists():
        run([
            cfg.tools.bcftools, "norm",
            "--multiallelics", "-",
            "--fasta-ref", str(ref),
            "--output-type", "z",
            "--output", str(tmp_norm),
            str(source),
        ])
    else:
        run([
            cfg.tools.bcftools, "norm",
            "--multiallelics", "-",
            "--output-type", "z",
            "--output", str(tmp_norm),
            str(source),
        ])
    run([cfg.tools.bcftools, "index", "--tbi", str(tmp_norm)])

    # Step 4: PLINK
    run([
        cfg.tools.plink2,
        "--vcf", str(tmp_norm),
        "--max-alleles", "2",
        "--snps-only", "just-acgt",
        "--make-bed",
        "--out", str(plink_stem),
        "--chr", chrom,
    ])

    # Step 5: BGEN
    bgen_path = Path(str(bgen_stem) + ".bgen")
    run([
        cfg.tools.qctool,
        "-g", str(tmp_norm),
        "-filetype", "vcf",
        "-og", str(bgen_path),
        "-ofiletype", "bgen_v1.2",
        "-bgen-bits", "8",
        "-os", str(sample_path),
    ])
    run([cfg.tools.bgenix, "-g", str(bgen_path), "-index"])

    # Cleanup
    for f in (tmp_intersected, tmp_reheadered, tmp_norm,
              Path(str(tmp_intersected) + ".tbi"),
              Path(str(tmp_reheadered) + ".tbi"),
              Path(str(tmp_norm) + ".tbi")):
        f.unlink(missing_ok=True)

    log(f"  done → {plink_stem.name}.bed and {bgen_stem.name}.bgen")


def _intersect_exome(
    cfg: Config, chrom: str, vcf: Path, bed: Path, out: Path
) -> None:
    if vcf_ok(out):
        log(f"  reusing intersected VCF (chr{chrom})")
        return
    remove(out, Path(str(out) + ".tbi"))  # clear any partial intersect from a prior crash
    log(f"  intersecting with exome BED (chr{chrom})…")
    # bcftools uses 1-based regions; BED is 0-based — bcftools handles this correctly
    # when using --regions-file with a BED file (it auto-detects BED format by filename)
    run([
        cfg.tools.bcftools, "view",
        "--regions-file", str(bed),
        "--output-type", "z",
        "--output", str(out),
        str(vcf),
    ])
    run([cfg.tools.bcftools, "index", "--tbi", str(out)])


def _write_sample_file(eid_map: pd.DataFrame, path: Path) -> None:
    if path.exists():
        return
    lines = ["ID_1 ID_2 missing sex\n", "0 0 0 D\n"]
    for _, row in eid_map.iterrows():
        eid = row["eid"]
        sex = int(row.get("sex", 0))
        lines.append(f"{eid} {eid} 0 {sex}\n")
    path.write_text("".join(lines))


def _load_eid_map(cfg: Config) -> pd.DataFrame:
    path = cfg.metadata_dir() / "eid_map_1kg.tsv"
    if not path.exists():
        raise FileNotFoundError(f"EID map not found: {path}\nRun `tbooo map eids` first.")
    return pd.read_csv(path, sep="\t")
