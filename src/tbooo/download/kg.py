"""Download 1000 Genomes Phase 3 and NYGC 30x data from EBI FTP."""

from __future__ import annotations

from pathlib import Path

from tbooo.config import Config
from tbooo.utils import ensure_dirs, has_bgzf_eof, log, parallel_download, run, wget_download

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

    # Per-chromosome VCFs + tabix indices (parallel)
    # Filename is derived from cfg.phase3_vcf() — sex chroms use different version strings.
    vcf_tasks: list[tuple[str, Path, str]] = []
    for chrom in chroms:
        vcf_path = cfg.phase3_vcf(chrom)
        vcf_url = f"{base}/{vcf_path.name}"
        vcf_tasks.append((vcf_url, vcf_path, wget))
    tbi_tasks = [(url + ".tbi", Path(str(p) + ".tbi"), w) for url, p, w in vcf_tasks]
    parallel_download(vcf_tasks + tbi_tasks, cfg.download_workers)

    _repair_truncated_vcfs(cfg, vcf_tasks, label="Phase 3")
    _verify_indices(cfg, vcf_tasks, label="Phase 3")

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

    vcf_tasks: list[tuple[str, Path, str]] = []
    for chrom in chroms:
        vcf_path = cfg.nygc_vcf(chrom)
        vcf_url = f"{base}/{vcf_path.name}"
        vcf_tasks.append((vcf_url, vcf_path, wget))
    tbi_tasks = [(url + ".tbi", Path(str(p) + ".tbi"), w) for url, p, w in vcf_tasks]
    parallel_download(vcf_tasks + tbi_tasks, cfg.download_workers)

    _repair_truncated_vcfs(cfg, vcf_tasks, label="NYGC 30x")
    _verify_indices(cfg, vcf_tasks, label="NYGC 30x")

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


def _repair_truncated_vcfs(
    cfg: Config,
    vcf_tasks: list[tuple[str, Path, str]],
    *,
    label: str,
) -> None:
    """Detect VCFs missing the BGZF EOF marker (interrupted download) and re-fetch
    each broken VCF along with its .tbi index."""
    log(f"Validating {len(vcf_tasks)} {label} VCF(s) for completeness…")
    broken = [(url, path, wget_bin)
              for url, path, wget_bin in vcf_tasks
              if path.exists() and not has_bgzf_eof(path)]

    if not broken:
        log("  All VCFs intact.")
        return

    log(f"  Found {len(broken)} truncated {label} VCF(s) — removing and re-downloading…")
    redownload: list[tuple[str, Path, str]] = []
    for url, path, wget_bin in broken:
        log(f"    removing truncated: {path.name}")
        path.unlink(missing_ok=True)
        tbi = Path(str(path) + ".tbi")
        if tbi.exists():
            tbi.unlink()
        redownload.append((url, path, wget_bin))
        redownload.append((url + ".tbi", tbi, wget_bin))

    parallel_download(redownload, cfg.download_workers)

    still_broken = [p.name for _, p, _ in broken if not has_bgzf_eof(p)]
    if still_broken:
        sample = ", ".join(still_broken[:5])
        more = f" (+{len(still_broken) - 5} more)" if len(still_broken) > 5 else ""
        raise RuntimeError(
            f"Re-download failed: {len(still_broken)} {label} VCF(s) still truncated: {sample}{more}"
        )
    log(f"  Re-downloaded {len(broken)} {label} VCF(s) successfully.")


def _verify_indices(
    cfg: Config,
    vcf_tasks: list[tuple[str, Path, str]],
    *,
    label: str,
) -> None:
    """Ensure every VCF has a non-empty .tbi alongside it.
    Re-download missing/zero-byte indices; rebuild locally if .tbi is older than the VCF."""
    log(f"Checking {len(vcf_tasks)} {label} VCF index file(s)…")
    redownload: list[tuple[str, Path, str]] = []
    rebuilt = 0
    ok = 0
    for url, vcf, wget_bin in vcf_tasks:
        if not vcf.exists():
            continue
        tbi = Path(str(vcf) + ".tbi")
        if not tbi.exists() or tbi.stat().st_size == 0:
            log(f"  missing/empty index: {tbi.name}")
            if tbi.exists():
                tbi.unlink()
            redownload.append((url + ".tbi", tbi, wget_bin))
            continue
        if tbi.stat().st_mtime < vcf.stat().st_mtime:
            log(f"  rebuilding stale index: {tbi.name}")
            run([cfg.tools.bcftools, "index", "--tbi", "-f", str(vcf)])
            rebuilt += 1
            continue
        ok += 1

    if redownload:
        log(f"  Re-downloading {len(redownload)} index file(s)…")
        parallel_download(redownload, cfg.download_workers)

    log(f"  Indices: {ok} up-to-date, {rebuilt} rebuilt, {len(redownload)} re-downloaded.")
