"""Download 1000 Genomes Phase 3 and NYGC 30x data from EBI FTP."""

from __future__ import annotations

from pathlib import Path

from tbooo.config import Config
from tbooo.utils import ensure_dirs, log, run, wget_download

# ── Phase 3 ──────────────────────────────────────────────────────────────────

_PHASE3_PANEL = "integrated_call_samples_v3.{date}.ALL.panel"
_PHASE3_PED = "20130606_g1k.ped"
_PHASE3_VCF = "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_{ver}.{date}.genotypes.vcf.gz"
_PHASE3_SV = "ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz"


def download_phase3(cfg: Config, chroms: list[str]) -> None:
    """Download Phase 3 per-chromosome VCFs, sample panel, and pedigree."""
    out = cfg.kg_raw_dir()
    ensure_dirs(out)
    base = cfg.kg_phase3_base_url
    date = cfg.kg_phase3_release_date
    ver = cfg.kg_phase3_vcf_version
    wget = cfg.tools.wget

    # Sample panel (needed for EID assignment and phenotype mapping)
    panel_name = _PHASE3_PANEL.format(date=date)
    _download_if_missing(f"{base}/{panel_name}", out / panel_name, wget)

    # Pedigree (needed for relatedness/family structure)
    _download_if_missing(
        f"https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/{_PHASE3_PED}",
        out / _PHASE3_PED,
        wget,
    )

    # Per-chromosome VCFs + tabix indices
    # Filename is derived from cfg.phase3_vcf() — sex chroms use different version strings.
    for chrom in chroms:
        vcf_path = cfg.phase3_vcf(chrom)
        vcf_url = f"{base}/{vcf_path.name}"
        _download_if_missing(vcf_url, vcf_path, wget)
        _download_if_missing(f"{vcf_url}.tbi", Path(str(vcf_path) + ".tbi"), wget)

    log(f"Phase 3 download complete → {out}")


# ── NYGC 30x ─────────────────────────────────────────────────────────────────

_NYGC_VCF = "CCDG_14151_B01_GRM_WGS_{date}_chr{chrom}.filtered.shapeit2-duohmm-phased.vcf.gz"
_NYGC_PANEL = "20130606_g1k_3202_samples_ped_population.txt"
_NYGC_PANEL_BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage"


def download_nygc(cfg: Config, chroms: list[str]) -> None:
    """Download NYGC 30x per-chromosome phased VCFs."""
    out = cfg.kg_raw_dir()
    ensure_dirs(out)
    base = cfg.kg_nygc_base_url
    date = cfg.kg_nygc_date
    wget = cfg.tools.wget

    # Extended sample panel (3,202 samples including related)
    _download_if_missing(
        f"{_NYGC_PANEL_BASE}/{_NYGC_PANEL}",
        out / _NYGC_PANEL,
        wget,
    )

    for chrom in chroms:
        vcf_name = _NYGC_VCF.format(date=date, chrom=chrom)
        vcf_url = f"{base}/{vcf_name}"
        vcf_path = out / vcf_name
        _download_if_missing(vcf_url, vcf_path, wget)
        _download_if_missing(f"{vcf_url}.tbi", Path(str(vcf_path) + ".tbi"), wget)

    log(f"NYGC 30x download complete → {out}")


# ── Genetic maps ─────────────────────────────────────────────────────────────

_GENETIC_MAP_TARBALL = "1000GP_Phase3_GRCh37_genetic_map.tar.gz"
_GENETIC_MAP_URL = (
    "https://github.com/joepickrell/1000-genomes-genetic-maps/archive/refs/heads/master.tar.gz"
)


def download_genetic_maps(cfg: Config) -> None:
    """Download 1KGP pedigree-based recombination maps from GitHub."""
    out = cfg.reference_dir / "genetic_maps"
    ensure_dirs(out)
    tarball = cfg.reference_dir / _GENETIC_MAP_TARBALL
    if not tarball.exists():
        log("Downloading genetic maps from GitHub…")
        wget_download(_GENETIC_MAP_URL, tarball, tool_wget=cfg.tools.wget)
        run(["tar", "-xzf", str(tarball), "-C", str(out), "--strip-components=1"])
    log(f"Genetic maps ready → {out}")


# ── Helpers ───────────────────────────────────────────────────────────────────

def _download_if_missing(url: str, dest: Path, wget_bin: str) -> None:
    if dest.exists():
        log(f"  skip (exists): {dest.name}")
        return
    log(f"  downloading: {dest.name}")
    wget_download(url, dest, tool_wget=wget_bin)
