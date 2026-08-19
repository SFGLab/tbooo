"""Builds synthetic source trees. No real cohort data is ever needed to test."""

from __future__ import annotations

from pathlib import Path

from tbooo.cohorts.registry import Cohort


def make_flat_bam_cohort(tmp_path: Path, name: str = "cases", n: int = 3) -> Cohort:
    root = tmp_path / f"{name}_raw"
    root.mkdir(parents=True, exist_ok=True)
    # Deliberately awkward IDs: discovery must not assume a clean naming scheme.
    ids = ["s001_x_20", "s-002-TAG", "s___003_20"][:n]
    for sid in ids:
        (root / f"{sid}.bam").write_bytes(b"BAM\x01")
        (root / f"{sid}.bam.bai").write_bytes(b"BAI\x01")
    return Cohort(
        name=name,
        eid_base=6_000_000,
        raw_root=root,
        discovery="flat_bam",
        modalities=("bam",),
    )


def make_sample_dirs_cohort(
    tmp_path: Path, name: str = "controls", with_vcf: bool = True
) -> Cohort:
    root = tmp_path / f"{name}_raw"
    root.mkdir(parents=True, exist_ok=True)
    vcf_root = tmp_path / f"{name}_vcf"
    vcf_root.mkdir(parents=True, exist_ok=True)

    # p01: single run, no sample-ID prefix on the files (one naming convention)
    d = root / "p01"
    d.mkdir()
    (d / "AA1234_L01_7_1.fq.gz").write_bytes(b"\x1f\x8b")
    (d / "AA1234_L01_7_2.fq.gz").write_bytes(b"\x1f\x8b")

    # p02: multi-lane, files prefixed with the sample ID (the other convention)
    d = root / "p02"
    d.mkdir()
    for fc, lane in (("BB5678", "L1"), ("BB9999", "L2")):
        (d / f"p02_{fc}_{lane}_1.fq.gz").write_bytes(b"\x1f\x8b")
        (d / f"p02_{fc}_{lane}_2.fq.gz").write_bytes(b"\x1f\x8b")
    (d / "MD5.txt").write_text("ignored\n")

    # Non-sample entries that must be excluded
    (root / "manifest.csv").write_text("x\n")
    (root / "checksums.txt").write_text("x\n")

    if with_vcf:
        for sid in ("p01", "p02"):
            (vcf_root / f"calls__{sid}.vcf").write_text(
                "##fileformat=VCFv4.2\n"
                f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sid}\n"
                "chr1\t100\tsv_1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=200\tGT\t0/1\n"
            )

    return Cohort(
        name=name,
        eid_base=7_000_000,
        raw_root=root,
        discovery="sample_dirs",
        modalities=("fastq", "sv_vcf") if with_vcf else ("fastq",),
        vcf_root=vcf_root if with_vcf else None,
        vcf_pattern="calls__{sample}.vcf",
        exclude=frozenset({"manifest.csv", "checksums.txt"}),
        icd10=None,
        platform_patterns=(("^AA", "platform_one"), ("^BB", "platform_two")),
        flowcell_regex=r"^(?P<flowcell>[A-Z0-9]+)_L\d+",
        pedigree_regex=None,
    )
