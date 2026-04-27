"""Download reference files needed by the pipeline."""

from __future__ import annotations

from pathlib import Path

from tbooo.config import Config
from tbooo.utils import ensure_dirs, log, run, wget_download

# IDT xGen Exome Research Panel v1.0, GRCh38
# Distributed by Illumina as part of their design files.
_IDT_EXOME_URL = (
    "https://support.illumina.com/content/dam/illumina-support/documents/downloads/"
    "productfiles/truseq/truseq-exome-targeted-regions-manifest-v1-2.bed"
)
# Fallback: Broad Institute's exome interval list (same capture panel, publicly hosted)
_IDT_EXOME_FALLBACK = (
    "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/"
    "xgen-exome-research-panel-v2-targets-hg38.bed"
)

_GRCH37_FASTA_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
    "human_g1k_v37.fasta.gz"
)
_GRCH38_FASTA_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/"
    "GRCh38_full_analysis_set_plus_decoy_hla.fa"
)


def download_all(cfg: Config) -> None:
    ensure_dirs(cfg.reference_dir)

    _download_exome_bed(cfg)
    _download_genetic_maps(cfg)
    _download_grch37(cfg)
    _download_grch38(cfg)

    log("All reference files ready.")


def _download_exome_bed(cfg: Config) -> None:
    dest = Path(cfg.exome_bed)
    if dest.exists():
        log(f"  skip (exists): {dest.name}")
        return
    ensure_dirs(dest.parent)
    log("Downloading IDT xGen exome BED (GRCh38)…")
    try:
        wget_download(_IDT_EXOME_FALLBACK, dest, tool_wget=cfg.tools.wget)
    except Exception:
        log("  primary URL failed, trying fallback…")
        wget_download(_IDT_EXOME_FALLBACK, dest, tool_wget=cfg.tools.wget)


def _download_genetic_maps(cfg: Config) -> None:
    maps_dir = cfg.reference_dir / "genetic_maps"
    tarball = cfg.reference_dir / "1000gp_genetic_maps.tar.gz"
    if maps_dir.exists() and any(maps_dir.iterdir()):
        log(f"  skip (exists): genetic_maps/")
        return
    ensure_dirs(maps_dir)
    url = "https://github.com/joepickrell/1000-genomes-genetic-maps/archive/refs/heads/master.tar.gz"
    log("Downloading 1KGP pedigree-based genetic maps…")
    wget_download(url, tarball, tool_wget=cfg.tools.wget)
    run(["tar", "-xzf", str(tarball), "-C", str(maps_dir), "--strip-components=1"])
    tarball.unlink()


def _download_grch37(cfg: Config) -> None:
    dest = cfg.reference_dir / "GRCh37" / "human_g1k_v37.fasta.gz"
    if dest.exists():
        log(f"  skip (exists): {dest.name}")
        return
    ensure_dirs(dest.parent)
    log("Downloading GRCh37 reference FASTA (hs37d5-compatible)…")
    wget_download(_GRCH37_FASTA_URL, dest, tool_wget=cfg.tools.wget)
    # Build .fai index
    run([cfg.tools.samtools, "faidx", str(dest)])


def _download_grch38(cfg: Config) -> None:
    dest = cfg.reference_dir / "GRCh38" / "GRCh38_full_analysis_set_plus_decoy_hla.fa"
    if dest.exists():
        log(f"  skip (exists): {dest.name}")
        return
    ensure_dirs(dest.parent)
    log("Downloading GRCh38DH reference FASTA…")
    wget_download(_GRCH38_FASTA_URL, dest, tool_wget=cfg.tools.wget)
    run([cfg.tools.samtools, "faidx", str(dest)])
