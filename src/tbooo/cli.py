"""tbooo — top-level CLI.

Usage:
    tbooo download 1kg        Download 1000 Genomes Phase 3 + NYGC 30x data
    tbooo download sgdp       Download SGDP CRAM files
    tbooo download reference  Download reference files (exome BED, genetic maps)

    tbooo map eids            Assign synthetic EIDs to all samples
    tbooo map array           Build UKB-mirrored PLINK array files (Field 22418)
    tbooo map imputed         Build UKB-mirrored BGEN imputed files (Field 22828)
    tbooo map wgs             Build UKB-mirrored WGS CRAMs + pVCF (Field 23149/23370)
    tbooo map wes             Build UKB-mirrored WES PLINK/BGEN files (Field 23157)
    tbooo map phenotypes      Build synthetic phenotype Parquet table
    tbooo map qc              Build sample QC and relatedness files

    tbooo run                 Run the full pipeline via Snakemake
"""

from __future__ import annotations

from pathlib import Path

import click

CONFIG_OPTION = click.option(
    "-c", "--config",
    default="config.yaml",
    show_default=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to config.yaml",
)


def _cfg(config_path: str):
    from tbooo.config import Config
    return Config.load(config_path)


# ── Root ─────────────────────────────────────────────────────────────────────

@click.group()
def cli():
    """TBOOO: simulate UK Biobank structure from 1000 Genomes + SGDP data."""


# ── download ──────────────────────────────────────────────────────────────────

@cli.group()
def download():
    """Download source datasets and reference files."""


@download.command("1kg")
@CONFIG_OPTION
@click.option("--phase3/--no-phase3", default=True, help="Download Phase 3 VCFs")
@click.option("--nygc/--no-nygc", default=True, help="Download NYGC 30x VCFs")
@click.option("--chroms", default=None, help="Comma-separated chromosomes (default: all)")
def download_1kg(config, phase3, nygc, chroms):
    """Download 1000 Genomes Phase 3 and/or NYGC 30x VCFs."""
    from tbooo.download.kg import download_phase3, download_nygc
    cfg = _cfg(config)
    chrom_list = chroms.split(",") if chroms else cfg.chromosomes
    if phase3:
        download_phase3(cfg, chrom_list)
    if nygc:
        download_nygc(cfg, chrom_list)


@download.command("sgdp")
@CONFIG_OPTION
@click.option("--populations", default=None, help="Comma-separated population names to download")
@click.option("--workers", default=None, type=int, help="Parallel download workers")
def download_sgdp(config, populations, workers):
    """Download SGDP CRAM files from ENA."""
    from tbooo.download.sgdp import download_cramps
    cfg = _cfg(config)
    pop_list = populations.split(",") if populations else cfg.sgdp_populations
    n_workers = workers or cfg.sgdp_download_workers
    download_cramps(cfg, pop_list, n_workers)


@download.command("reference")
@CONFIG_OPTION
def download_reference(config):
    """Download reference files: exome BED, genetic maps, GRCh37/38 FASTAs."""
    from tbooo.download.reference import download_all
    download_all(_cfg(config))


# ── map ───────────────────────────────────────────────────────────────────────

@cli.group()
def map():
    """Transform source data into UKB-mirrored structure."""


@map.command("eids")
@CONFIG_OPTION
def map_eids(config):
    """Assign synthetic EIDs to all 1KGP and SGDP samples."""
    from tbooo.pipeline.eids import assign_eids
    assign_eids(_cfg(config))


@map.command("array")
@CONFIG_OPTION
@click.option("--chroms", default=None, help="Comma-separated chromosomes")
def map_array(config, chroms):
    """Build PLINK array files (ukb22418_c*_b0_v2.bed/bim/fam)."""
    from tbooo.pipeline.array import run_array_pipeline
    cfg = _cfg(config)
    chrom_list = chroms.split(",") if chroms else [str(c) for c in cfg.autosomes]
    run_array_pipeline(cfg, chrom_list)


@map.command("imputed")
@CONFIG_OPTION
@click.option("--chroms", default=None, help="Comma-separated chromosomes")
def map_imputed(config, chroms):
    """Build BGEN imputed files (ukb22828_c*_b0_v3.bgen/sample/bgi)."""
    from tbooo.pipeline.imputed import run_imputed_pipeline
    cfg = _cfg(config)
    chrom_list = chroms.split(",") if chroms else [str(c) for c in cfg.autosomes]
    run_imputed_pipeline(cfg, chrom_list)


@map.command("wgs")
@CONFIG_OPTION
@click.option("--chroms", default=None, help="Comma-separated chromosomes for pVCF")
@click.option("--croms/--no-croms", "do_croms", default=True, help="Rename individual CRAMs")
@click.option("--pvcf/--no-pvcf", default=True, help="Build cohort pVCF per chromosome")
def map_wgs(config, chroms, do_croms, pvcf):
    """Rename CRAMs and build cohort pVCF (Field 23149/23370)."""
    from tbooo.pipeline.wgs import rename_crams, build_pvcf
    cfg = _cfg(config)
    if do_croms:
        rename_crams(cfg)
    if pvcf:
        chrom_list = chroms.split(",") if chroms else cfg.chromosomes
        build_pvcf(cfg, chrom_list)


@map.command("wes")
@CONFIG_OPTION
@click.option("--chroms", default=None, help="Comma-separated chromosomes")
def map_wes(config, chroms):
    """Build WES PLINK + BGEN files by intersecting WGS with exome BED."""
    from tbooo.pipeline.wes import run_wes_pipeline
    cfg = _cfg(config)
    chrom_list = chroms.split(",") if chroms else [str(c) for c in cfg.autosomes]
    run_wes_pipeline(cfg, chrom_list)


@map.command("phenotypes")
@CONFIG_OPTION
def map_phenotypes(config):
    """Build synthetic phenotype Parquet table (Showcase/participant.parquet)."""
    from tbooo.pipeline.phenotypes import build_phenotype_table
    build_phenotype_table(_cfg(config))


@map.command("qc")
@CONFIG_OPTION
def map_qc(config):
    """Build sample QC file and relatedness file."""
    from tbooo.pipeline.qc import build_qc_files
    build_qc_files(_cfg(config))


# ── run ───────────────────────────────────────────────────────────────────────

@cli.command("run")
@CONFIG_OPTION
@click.option("-j", "--jobs", default=4, show_default=True, help="Snakemake parallel jobs")
@click.option("--dry-run", "-n", is_flag=True, help="Dry-run: print rules without executing")
@click.option("--target", default=None, help="Specific Snakemake target rule or file")
def run_pipeline(config, jobs, dry_run, target):
    """Run the full pipeline via Snakemake."""
    import subprocess, sys
    cmd = ["snakemake", "--configfile", config, "--cores", str(jobs)]
    if dry_run:
        cmd.append("--dry-run")
    if target:
        cmd.append(target)
    result = subprocess.run(cmd)
    sys.exit(result.returncode)
