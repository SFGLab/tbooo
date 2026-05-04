"""Full pipeline runner — pure Python, no Snakemake.

Stage order and dependencies:
    downloads (parallel)
        └─ map eids
            └─ array / imputed / wes / wgs / geuvadis  (per-chrom parallel, stages sequential)
                └─ phenotypes
                    └─ qc
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Callable

from tbooo.config import Config
from tbooo.utils import log

# All recognised target names for --target
TARGETS = ("array", "imputed", "wes", "wgs", "geuvadis", "phenotypes", "qc")


def run(
    cfg: Config,
    *,
    jobs: int = 4,
    dry_run: bool = False,
    chroms: list[str] | None = None,
    target: str | None = None,
    skip_download: bool = False,
) -> None:
    chrom_list = chroms or cfg.chromosomes
    auto_list = [c for c in chrom_list if c in [str(a) for a in cfg.autosomes]]

    def step(label: str, fn: Callable, *args, **kwargs) -> None:
        if dry_run:
            log(f"  [dry-run] {label}")
            return
        log(f"── {label}")
        fn(*args, **kwargs)

    def parallel_chroms(label: str, fn: Callable, chroms: list[str]) -> None:
        """Run fn(cfg, [chrom]) for each chrom using up to `jobs` threads."""
        if not chroms:
            return
        if dry_run:
            for c in chroms:
                log(f"  [dry-run] {label} [chr{c}]")
            return
        log(f"── {label} ({len(chroms)} chromosomes, {jobs} parallel)")
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = {pool.submit(fn, cfg, [c]): c for c in chroms}
            for future in as_completed(futures):
                future.result()  # re-raise any exception immediately

    # ── Downloads ─────────────────────────────────────────────────────────────
    if not skip_download and target is None:
        from tbooo.download.kg import download_phase3, download_nygc
        from tbooo.download.sgdp import download_metadata as _sgdp_meta
        from tbooo.download.sgdp import download_vcfs as _sgdp_vcfs
        from tbooo.download.geuvadis import download_expression
        from tbooo.download.reference import download_all as download_reference

        step("download reference",      download_reference, cfg)
        step("download 1kg phase3",     download_phase3,    cfg, chrom_list)
        step("download 1kg nygc",       download_nygc,      cfg, chrom_list)
        step("download sgdp metadata",  _sgdp_meta,         cfg)
        step("download sgdp vcfs",      _sgdp_vcfs,         cfg)
        step("download geuvadis",       download_expression, cfg)

    # ── EID assignment ─────────────────────────────────────────────────────────
    if target is None:
        from tbooo.pipeline.eids import assign_eids
        step("map eids", assign_eids, cfg)

    # ── Array ──────────────────────────────────────────────────────────────────
    if target in (None, "array"):
        from tbooo.pipeline.array import run_array_pipeline
        parallel_chroms("map array", run_array_pipeline, auto_list)

    # ── Imputed ────────────────────────────────────────────────────────────────
    if target in (None, "imputed"):
        from tbooo.pipeline.imputed import run_imputed_pipeline
        parallel_chroms("map imputed", run_imputed_pipeline, auto_list)

    # ── WES ────────────────────────────────────────────────────────────────────
    if target in (None, "wes"):
        from tbooo.pipeline.wes import run_wes_pipeline
        parallel_chroms("map wes", run_wes_pipeline, auto_list)

    # ── WGS ────────────────────────────────────────────────────────────────────
    if target in (None, "wgs"):
        from tbooo.pipeline.wgs import rename_crams, build_pvcf
        step("map wgs crams", rename_crams, cfg)
        parallel_chroms("map wgs pvcf", build_pvcf, chrom_list)

    # ── GEUVADIS expression PCA ────────────────────────────────────────────────
    if target in (None, "geuvadis"):
        from tbooo.pipeline.geuvadis import build_expression_pcs
        step("map geuvadis", build_expression_pcs, cfg)

    # ── Phenotypes ─────────────────────────────────────────────────────────────
    if target in (None, "phenotypes"):
        from tbooo.pipeline.phenotypes import build_phenotype_table
        step("map phenotypes", build_phenotype_table, cfg)

    # ── QC ─────────────────────────────────────────────────────────────────────
    if target in (None, "qc"):
        from tbooo.pipeline.qc import build_qc_files
        step("map qc", build_qc_files, cfg)

    log("Pipeline complete." if not dry_run else "Dry-run complete.")
