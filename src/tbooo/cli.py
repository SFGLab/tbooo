"""tbooo — top-level CLI.

Usage:
    tbooo download 1kg        Download 1000 Genomes Phase 3 + NYGC 30x data
    tbooo download sgdp       Download SGDP CRAM files
    tbooo download geuvadis   Download GEUVADIS GD462 RPKM expression matrix
    tbooo download reference  Download reference files (exome BED, genetic maps)

    tbooo map eids            Assign synthetic EIDs to all samples
    tbooo map array           Build UKB-mirrored PLINK array files (Field 22418)
    tbooo map imputed         Build UKB-mirrored BGEN imputed files (Field 22828)
    tbooo map wgs             Build UKB-mirrored WGS CRAMs + pVCF (Field 23149/23370)
    tbooo map wes             Build UKB-mirrored WES PLINK/BGEN files (Field 23157)
    tbooo map geuvadis        Compute expression PCA and add geuvadis_pc* to participant.parquet
    tbooo map phenotypes      Build synthetic phenotype Parquet table
    tbooo map qc              Build sample QC and relatedness files
    tbooo map eda             EDA plots: sex/region demographics + RNA PCA

    tbooo run                 Run the full pipeline
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
@click.option("--vcf/--no-vcf", "do_vcf", default=True, help="Download per-sample VCF files")
def download_sgdp(config, do_vcf):
    """Download SGDP sample metadata and VCF files from ENA (no CRAMs)."""
    from tbooo.download.sgdp import download_metadata, download_vcfs
    cfg = _cfg(config)
    download_metadata(cfg)
    if do_vcf:
        download_vcfs(cfg)


@download.command("geuvadis")
@CONFIG_OPTION
def download_geuvadis(config):
    """Download GEUVADIS GD462 RPKM gene expression matrix."""
    from tbooo.download.geuvadis import download_expression
    download_expression(_cfg(config))


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
@click.option("--chroms", default=None, help="Comma-separated chromosomes")
@click.option("--croms/--no-croms", "do_croms", default=True, help="Rename individual CRAMs (Field 23149)")
@click.option("--gvcf/--no-gvcf", "do_gvcf", default=False,
              help="Build per-sample gVCFs from NYGC + SGDP (Field 23151)")
@click.option("--pvcf/--no-pvcf", "do_pvcf", default=True,
              help="Build cohort pVCF from NYGC + SGDP merged (Field 23370)")
def map_wgs(config, chroms, do_croms, do_gvcf, do_pvcf):
    """Build WGS outputs (Fields 23149/23151/23370) from NYGC and SGDP VCFs."""
    from tbooo.pipeline.wgs import rename_crams, build_gvcfs, build_pvcf
    cfg = _cfg(config)
    chrom_list = chroms.split(",") if chroms else cfg.chromosomes
    if do_croms:
        rename_crams(cfg)
    if do_gvcf:
        build_gvcfs(cfg, chrom_list)
    if do_pvcf:
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


@map.command("geuvadis")
@CONFIG_OPTION
def map_geuvadis(config):
    """Compute GEUVADIS expression PCA and merge scores into participant.parquet."""
    from tbooo.pipeline.geuvadis import build_expression_pcs
    build_expression_pcs(_cfg(config))


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


@map.command("eda")
@CONFIG_OPTION
def map_eda(config):
    """EDA plots: sex/region demographics + GEUVADIS RNA PCA."""
    from tbooo.pipeline.eda import run_eda
    run_eda(_cfg(config))


# ── run ───────────────────────────────────────────────────────────────────────

@cli.command("run")
@CONFIG_OPTION
@click.option("-j", "--jobs", default=4, show_default=True,
              help="Parallel jobs for per-chromosome pipeline steps")
@click.option("-n", "--dry-run", is_flag=True,
              help="Print steps without executing")
@click.option("--chroms", default=None,
              help="Comma-separated chromosomes to process (default: all)")
@click.option("--target", default=None,
              type=click.Choice(["array", "imputed", "wes", "wgs",
                                 "geuvadis", "phenotypes", "qc", "eda"]),
              help="Run only this stage (assumes prerequisites already exist)")
@click.option("--no-download", is_flag=True,
              help="Skip all download steps")
def run_pipeline(config, jobs, dry_run, chroms, target, no_download):
    """Run the full pipeline (or a single stage)."""
    from tbooo.pipeline.runner import run
    cfg = _cfg(config)
    chrom_list = chroms.split(",") if chroms else None
    run(cfg, jobs=jobs, dry_run=dry_run, chroms=chrom_list,
        target=target, skip_download=no_download)
