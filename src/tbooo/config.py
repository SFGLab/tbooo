from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class Tools:
    bcftools: str = "bcftools"
    plink2: str = "plink2"
    qctool: str = "qctool"
    bgenix: str = "bgenix"
    samtools: str = "samtools"
    king: str = "king"
    wget: str = "wget"


@dataclass
class Config:
    # directories
    data_dir: Path
    reference_dir: Path
    tmp_dir: Path

    # chromosomes
    autosomes: list[int]
    sex_chromosomes: list[str]

    # 1KGP Phase 3
    kg_phase3_release_date: str
    kg_phase3_vcf_version: str
    kg_phase3_base_url: str

    # 1KGP NYGC 30x
    kg_nygc_date: str
    kg_nygc_base_url: str

    # SGDP (VCF + metadata; no CRAM downloads)
    sgdp_populations: list[str]

    # EID ranges
    kg_eid_start: int
    sgdp_eid_start: int

    # array simulation
    array_manifest: str
    array_proxy_maf: float

    # WES simulation
    exome_bed: str

    # tools
    tools: Tools

    # parallel downloads
    download_workers: int = 4

    # parallel workers for `bgzip -t` deep integrity checks (opt-in via --deep-check)
    deep_check_workers: int = 8

    @classmethod
    def load(cls, path: Path | str = "config.yaml") -> "Config":
        with open(path) as f:
            raw: dict[str, Any] = yaml.safe_load(f)
        raw["data_dir"] = Path(raw["data_dir"])
        raw["reference_dir"] = Path(raw["reference_dir"])
        raw["tmp_dir"] = Path(raw["tmp_dir"])
        raw["tools"] = Tools(**raw.get("tools", {}))
        return cls(**raw)

    # ── Derived path helpers ──────────────────────────────────────────────────

    @property
    def chromosomes(self) -> list[str]:
        return [str(c) for c in self.autosomes] + list(self.sex_chromosomes)

    def bulk_dir(self, *parts: str) -> Path:
        return self.data_dir / "Bulk" / Path(*parts)

    def array_dir(self) -> Path:
        return self.bulk_dir("Genotype Results", "Genotype calls")

    def imputed_dir(self) -> Path:
        return self.bulk_dir("Imputed")

    def wgs_dir(self) -> Path:
        return self.bulk_dir("Whole genome sequences")

    def wes_dir(self) -> Path:
        return self.bulk_dir(
            "Exome sequences",
            "Population level exome OQFE variants, PLINK format - 500k release",
        )

    def wes_bgen_dir(self) -> Path:
        return self.bulk_dir(
            "Exome sequences",
            "Population level exome OQFE variants, BGEN format - final release",
        )

    def showcase_dir(self) -> Path:
        return self.data_dir / "Showcase"

    def metadata_dir(self) -> Path:
        return self.data_dir / "metadata"

    def kg_raw_dir(self) -> Path:
        return self.data_dir / "raw" / "1kg"

    def sgdp_raw_dir(self) -> Path:
        return self.data_dir / "raw" / "sgdp"

    def geuvadis_raw_dir(self) -> Path:
        return self.data_dir / "raw" / "geuvadis"

    def sgdp_vcf_dir(self) -> Path:
        return self.sgdp_raw_dir() / "vcf"

    def sgdp_pvcf(self, chrom: str) -> Path:
        return self.sgdp_raw_dir() / "pvcf" / f"sgdp_c{chrom}.pvcf.gz"

    def nygc_pvcf(self, chrom: str) -> Path:
        """Intermediate per-chromosome NYGC pVCF (EID-reheadered, before merging with SGDP)."""
        return self.kg_raw_dir() / "pvcf" / f"nygc_c{chrom}.pvcf.gz"

    def phase3_vcf(self, chrom: str) -> Path:
        date = self.kg_phase3_release_date
        ver = self.kg_phase3_vcf_version
        if chrom == "X":
            return self.kg_raw_dir() / f"ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.{date}.genotypes.vcf.gz"
        if chrom == "Y":
            return self.kg_raw_dir() / f"ALL.chrY.phase3_integrated_v2b.{date}.genotypes.vcf.gz"
        return self.kg_raw_dir() / f"ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_{ver}.{date}.genotypes.vcf.gz"

    def nygc_vcf(self, chrom: str) -> Path:
        date = self.kg_nygc_date
        if chrom == "X":
            return self.kg_raw_dir() / f"CCDG_14151_B01_GRM_WGS_{date}_chrX.filtered.eagle2-phased.v2.vcf.gz"
        return self.kg_raw_dir() / f"CCDG_14151_B01_GRM_WGS_{date}_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz"

    def array_stem(self, chrom: str) -> Path:
        return self.array_dir() / f"ukb22418_c{chrom}_b0_v2"

    def imputed_stem(self, chrom: str) -> Path:
        return self.imputed_dir() / f"ukb22828_c{chrom}_b0_v3"

    def wes_stem(self, chrom: str) -> Path:
        return self.wes_dir() / f"ukb23157_c{chrom}_b0_v1"

    def wes_bgen_stem(self, chrom: str) -> Path:
        return self.wes_bgen_dir() / f"ukb23157_c{chrom}_b0_v1"

    def wgs_pvcf(self, chrom: str) -> Path:
        return self.wgs_dir() / f"ukb23370_c{chrom}_b0_v1.pvcf.gz"
