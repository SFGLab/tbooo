"""tbooo — top-level CLI.

Usage:
    tbooo download 1kg        Download 1000 Genomes Phase 3 + NYGC 30x data
    tbooo download sgdp       Download SGDP CRAM files
    tbooo download geuvadis   Download GEUVADIS GD462 RPKM expression matrix
    tbooo download mage       Download MAGE expression + genotype PCs (731 individuals)
    tbooo download gtex       Download GTEx per-tissue TPM + eQTL covariates
    tbooo download reference  Download reference files (exome BED, genetic maps)

    tbooo map eids            Assign synthetic EIDs to all samples
    tbooo map array           Build UKB-mirrored PLINK array files (Field 22418)
    tbooo map imputed         Build UKB-mirrored BGEN imputed files (Field 22828)
    tbooo map wgs             Build UKB-mirrored WGS gVCFs + pVCF (Field 24051/24310)
    tbooo map wes             Build UKB-mirrored WES PLINK/BGEN files (Field 23158/23159)
    tbooo map geuvadis        Compute expression PCA and add geuvadis_pc* to participant.parquet
    tbooo map mage            Compute MAGE expression PCA and add mage_pc* columns
    tbooo map gtex            Assign GTEx EIDs (3,000,000+) and append gtex_pc* rows
    tbooo map phenotypes      Build synthetic phenotype Parquet table
    tbooo map qc              Build sample QC and relatedness files
    tbooo map eda             EDA plots: sex/region demographics + RNA PCA

    tbooo restricted discover    Enumerate samples in each registered cohort
    tbooo restricted eids        Assign EIDs and write rename files
    tbooo restricted link        Build the symlink tree + provenance manifests
    tbooo restricted sv          Emit the SLURM script for the SV layer
    tbooo restricted phenotypes  Build participant.parquet (always all cohorts)

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
@click.option("--deep-check", is_flag=True,
              help="Run `bgzip -t` end-to-end integrity check on stale-index VCFs "
                   "(catches mid-stream corruption; uses deep_check_workers from config)")
def download_1kg(config, phase3, nygc, chroms, deep_check):
    """Download 1000 Genomes Phase 3 and/or NYGC 30x VCFs."""
    from tbooo.download.kg import download_phase3, download_nygc
    cfg = _cfg(config)
    chrom_list = chroms.split(",") if chroms else cfg.chromosomes
    if phase3:
        download_phase3(cfg, chrom_list, deep_check=deep_check)
    if nygc:
        download_nygc(cfg, chrom_list, deep_check=deep_check)


@download.command("sgdp")
@CONFIG_OPTION
@click.option("--vcf/--no-vcf", "do_vcf", default=True, help="Download per-sample VCF files")
@click.option("--deep-check", is_flag=True,
              help="Run `bgzip -t` end-to-end integrity check on every VCF "
                   "(catches mid-stream corruption; uses deep_check_workers from config)")
@click.option("--tolerate-broken", is_flag=True,
              help="If a VCF is still broken after re-download, drop it (and its sample) "
                   "from the dataset instead of crashing. Affected samples will be absent "
                   "from downstream outputs.")
def download_sgdp(config, do_vcf, deep_check, tolerate_broken):
    """Download SGDP sample metadata and VCF files from ENA (no CRAMs)."""
    from tbooo.download.sgdp import download_metadata, download_vcfs
    cfg = _cfg(config)
    download_metadata(cfg)
    if do_vcf:
        download_vcfs(cfg, deep_check=deep_check, tolerate_broken=tolerate_broken)


@download.command("geuvadis")
@CONFIG_OPTION
def download_geuvadis(config):
    """Download GEUVADIS GD462 RPKM gene expression matrix."""
    from tbooo.download.geuvadis import download_expression
    download_expression(_cfg(config))


@download.command("mage")
@CONFIG_OPTION
def download_mage_cmd(config):
    """Download MAGE expression, metadata and genotype PCs (selective, ~165 MB)."""
    from tbooo.download.mage import download_mage
    download_mage(_cfg(config))


@download.command("gtex")
@CONFIG_OPTION
def download_gtex_cmd(config):
    """Download GTEx per-tissue TPM, eQTL covariates and annotations."""
    from tbooo.download.gtex import download_gtex
    download_gtex(_cfg(config))


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
@click.option("--gvcf/--no-gvcf", "do_gvcf", default=True,
              help="Build per-sample gVCFs from NYGC + SGDP (Field 24051)")
@click.option("--pvcf/--no-pvcf", "do_pvcf", default=True,
              help="Build cohort pVCF from NYGC + SGDP merged (Field 24310)")
def map_wgs(config, chroms, do_gvcf, do_pvcf):
    """Build WGS outputs (Fields 24051/24310) from NYGC and SGDP VCFs."""
    from tbooo.pipeline.wgs import build_gvcfs, build_pvcf
    cfg = _cfg(config)
    chrom_list = chroms.split(",") if chroms else cfg.chromosomes
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


@map.command("mage")
@CONFIG_OPTION
def map_mage(config):
    """Compute MAGE expression PCA and merge mage_pc* into participant.parquet."""
    from tbooo.pipeline.mage import build_expression_pcs
    build_expression_pcs(_cfg(config))


@map.command("gtex")
@CONFIG_OPTION
def map_gtex(config):
    """Assign GTEx donors EIDs and append gtex_pc* rows to participant.parquet."""
    from tbooo.pipeline.gtex import build_expression_pcs
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


# ── restricted ────────────────────────────────────────────────────────────────

@cli.group()
def restricted():
    """Ingest externally-held cohorts defined in a private registry.

    Requires `cohorts_registry` in the config, pointing at a cohorts.yaml.
    See cohorts.example.yaml for the schema. The registry is gitignored.
    """


def _load_cohorts(cfg, only: str | None):
    """Discover every registered cohort (or just `only`) and assign EIDs."""
    import click as _click

    from tbooo.cohorts.discover import discover_samples
    from tbooo.cohorts.eids import assign_eids
    from tbooo.cohorts.phenotypes import CohortData
    from tbooo.cohorts.registry import RegistryError, load_registry

    if not cfg.cohorts_registry:
        raise _click.ClickException(
            "no cohorts_registry set in config. Copy cohorts.example.yaml to "
            "cohorts.yaml, fill it in, and point cohorts_registry at it."
        )
    try:
        registry = load_registry(cfg.cohorts_registry)
    except RegistryError as exc:
        raise _click.ClickException(str(exc)) from exc

    if only:
        if only not in registry:
            raise _click.ClickException(
                f"cohort {only!r} not in registry (have: {', '.join(sorted(registry))})"
            )
        registry = {only: registry[only]}

    data = []
    for name, cohort in sorted(registry.items()):
        try:
            samples = discover_samples(cohort)
        except RegistryError as exc:  # DiscoveryError included
            raise _click.ClickException(str(exc)) from exc
        data.append(CohortData(cohort, samples, assign_eids(cohort, samples)))
    return data


ONLY_OPTION = click.option(
    "--only", default=None, help="Restrict to a single cohort by name"
)


@restricted.command("discover")
@CONFIG_OPTION
@ONLY_OPTION
def restricted_discover(config, only):
    """Enumerate samples in each registered cohort (no writes)."""
    from tbooo.utils import log

    for item in _load_cohorts(_cfg(config), only):
        n_sv = sum(
            1 for s in item.samples if any(f.modality == "sv_vcf" for f in s.files)
        )
        log(f"{item.cohort.name}: {len(item.samples)} samples, {n_sv} with SV calls")


@restricted.command("eids")
@CONFIG_OPTION
@ONLY_OPTION
def restricted_eids(config, only):
    """Assign EIDs and write the eid map + reheader rename files."""
    from tbooo.cohorts.eids import write_eid_map, write_rename_file

    cfg = _cfg(config)
    for item in _load_cohorts(cfg, only):
        write_eid_map(cfg, item.cohort, item.samples, item.eids)
        write_rename_file(cfg, item.cohort, item.eids)


@restricted.command("link")
@CONFIG_OPTION
@ONLY_OPTION
def restricted_link(config, only):
    """Build the UKB-shaped symlink tree and provenance manifests."""
    from tbooo.cohorts.link import link_cohort

    cfg = _cfg(config)
    for item in _load_cohorts(cfg, only):
        try:
            link_cohort(cfg, item.cohort, item.samples, item.eids)
        except RuntimeError as exc:
            # Same surfacing as `restricted sv`: the containment and blind-prune
            # refusals are operator-facing errors, not bugs to show as tracebacks.
            raise click.ClickException(str(exc)) from exc


@restricted.command("sv")
@CONFIG_OPTION
@ONLY_OPTION
@click.option("--partition", default="short", show_default=True,
              help="SLURM partition for the generated job")
@click.option("--time-limit", default="01:00:00", show_default=True,
              help="SLURM time limit for the generated job")
def restricted_sv(config, only, partition, time_limit):
    """Emit the SLURM script that builds the SV layer. Runs nothing itself."""
    from tbooo.cohorts.sv import emit_sv_script

    cfg = _cfg(config)
    for item in _load_cohorts(cfg, only):
        try:
            emit_sv_script(
                cfg, item.cohort, item.samples, item.eids,
                partition=partition, time_limit=time_limit,
            )
        except RuntimeError as exc:
            raise click.ClickException(str(exc)) from exc


@restricted.command("phenotypes")
@CONFIG_OPTION
@ONLY_OPTION
def restricted_phenotypes(config, only):
    """Build participant.parquet, batch codes and putative pedigrees.

    Unlike every other restricted stage, this one always covers every cohort in
    the registry: participant.parquet and batch_codes.tsv are single tables
    spanning all of them and are rewritten in full on each run. `--only` is
    therefore rejected here rather than silently deleting the other cohorts.
    """
    from tbooo.cohorts.phenotypes import build_participants

    if only:
        raise click.ClickException(
            "`--only` is not supported by `restricted phenotypes`. "
            "participant.parquet and batch_codes.tsv are single tables covering "
            "every cohort, rewritten in full on each run, so building them from "
            "one cohort would delete every other cohort's rows. The other stages "
            "are per-cohort and accept --only; this one does not. Re-run without it."
        )

    cfg = _cfg(config)
    build_participants(cfg, _load_cohorts(cfg, None))


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
                                 "geuvadis", "mage", "gtex", "phenotypes", "qc", "eda"]),
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
