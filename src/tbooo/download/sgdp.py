"""Download SGDP VCF files and sample metadata.

Metadata is fetched from the IGSR portal API (internationalgenome.org).
Per-sample phased VCF files are fetched from ENA analysis results for PRJEB9586.
SGDP CRAMs (~14 TB) are not downloaded.

Outputs:
    data/raw/sgdp/sgdp_samples.tsv
        columns: ena_accession, sample_alias, population, region, sex
    data/raw/sgdp/vcf/<accession>.vcf.gz  (one per sample)
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import pandas as pd

from tbooo.config import Config
from tbooo.utils import ensure_dirs, has_bgzf_eof, log, parallel_download, run, wget_download

# IGSR portal API — returns all samples across all collections
_IGSR_SAMPLE_API = "https://www.internationalgenome.org/api/beta/sample/_search/igsr_samples.tsv"
_IGSR_PAYLOAD = {
    "fields": [
        "name", "sex", "biosampleId",
        "populations.code", "populations.name",
        "populations.superpopulationCode", "populations.superpopulationName",
        "populations.elasticId",
        "dataCollections.title",
    ],
    "column_names": [
        "Sample name", "Sex", "Biosample ID",
        "Population code", "Population name",
        "Superpopulation code", "Superpopulation name",
        "Population elastic ID",
        "Data collections",
    ],
}

# ENA portal API — analysis-level records (processed files, including VCFs)
_ENA_ANALYSIS_API = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession=PRJEB9586"
    "&result=analysis"
    "&fields=analysis_accession,sample_accession,submitted_ftp"
    "&format=tsv"
    "&limit=0"
)


def download_metadata(cfg: Config) -> None:
    """Fetch SGDP sample metadata from the IGSR portal API."""
    out = cfg.sgdp_raw_dir() / "sgdp_samples.tsv"
    if out.exists():
        log(f"  skip (exists): {out.name}")
        return
    ensure_dirs(cfg.sgdp_raw_dir())

    log("Fetching SGDP sample metadata from IGSR portal…")
    try:
        body = json.dumps(_IGSR_PAYLOAD).encode()
        req = urllib.request.Request(
            _IGSR_SAMPLE_API,
            data=body,
            headers={"Content-Type": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=120) as resp:
            tsv = resp.read().decode("utf-8")
    except Exception as exc:
        log(f"  IGSR API failed: {exc}")
        log("  Writing empty metadata file — SGDP rows will be absent from phenotype table.")
        out.write_text("ena_accession\tsample_alias\tpopulation\tregion\tsex\n")
        return

    df = _parse_igsr_response(tsv, cfg.sgdp_populations)
    df.to_csv(out, sep="\t", index=False)
    log(f"  wrote {out} ({len(df)} samples)")


def download_vcfs(cfg: Config) -> None:
    """Download per-sample phased VCF files from ENA analysis results for PRJEB9586."""
    vcf_dir = cfg.sgdp_vcf_dir()
    ensure_dirs(vcf_dir)

    # Resolve population filter against the metadata table
    allowed: set[str] | None = None
    if cfg.sgdp_populations:
        samples_tsv = cfg.sgdp_raw_dir() / "sgdp_samples.tsv"
        if not samples_tsv.exists():
            download_metadata(cfg)
        if samples_tsv.exists() and samples_tsv.stat().st_size > 0:
            meta = pd.read_csv(samples_tsv, sep="\t")
            allowed = set(meta.loc[meta["population"].isin(cfg.sgdp_populations), "ena_accession"])
            log(f"  Population filter active: {len(allowed)} samples from {len(cfg.sgdp_populations)} populations")

    log("Querying ENA for SGDP VCF analysis files (PRJEB9586)…")
    try:
        url_map = _fetch_vcf_urls(allowed)
    except Exception as exc:
        log(f"  ENA analysis query failed: {exc}")
        log("  Cannot resolve SGDP VCF URLs automatically.")
        log("  Obtain files manually from:")
        log("    https://www.internationalgenome.org/data-portal/data-collection/SGDP/")
        return

    if not url_map:
        log("  No VCF analysis files found in ENA for PRJEB9586.")
        log("  Obtain files manually from:")
        log("    https://www.internationalgenome.org/data-portal/data-collection/SGDP/")
        return

    log(f"  Found {len(url_map)} VCF file(s) to download.")
    tasks = [(url, vcf_dir / filename, cfg.tools.wget) for filename, (_, url) in url_map.items()]
    parallel_download(tasks, cfg.download_workers)

    _repair_truncated_vcfs(cfg, url_map, vcf_dir)
    _index_vcfs(cfg, [vcf_dir / fn for fn in url_map])

    _link_vcfs_to_metadata(cfg, url_map)
    log("SGDP VCF download complete.")


def _repair_truncated_vcfs(cfg: Config, url_map: dict, vcf_dir: Path) -> None:
    """Detect VCFs missing the BGZF EOF marker (interrupted download) and re-fetch them."""
    log(f"Validating {len(url_map)} SGDP VCF(s) for completeness…")
    broken: list[tuple[str, Path]] = []
    for filename, (_, url) in url_map.items():
        path = vcf_dir / filename
        if path.exists() and not has_bgzf_eof(path):
            broken.append((url, path))

    if not broken:
        log("  All VCFs intact.")
        return

    log(f"  Found {len(broken)} truncated VCF(s) — removing and re-downloading…")
    for _, path in broken:
        log(f"    removing truncated: {path.name}")
        path.unlink(missing_ok=True)
        tbi = path.with_suffix(path.suffix + ".tbi")
        tbi.unlink(missing_ok=True)

    tasks = [(url, path, cfg.tools.wget) for url, path in broken]
    parallel_download(tasks, cfg.download_workers)

    still_broken = [p.name for _, p in broken if not has_bgzf_eof(p)]
    if still_broken:
        sample = ", ".join(still_broken[:5])
        more = f" (+{len(still_broken) - 5} more)" if len(still_broken) > 5 else ""
        raise RuntimeError(
            f"Re-download failed: {len(still_broken)} VCF(s) still truncated: {sample}{more}"
        )
    log(f"  Re-downloaded {len(broken)} VCF(s) successfully.")


def _index_vcfs(cfg: Config, vcfs: list[Path]) -> None:
    """Tabix-index each VCF. Skip if index is up-to-date; re-index with -f if stale or missing."""
    present = [v for v in vcfs if v.exists()]
    log(f"Indexing {len(present)} SGDP VCF(s)…")
    indexed = skipped = rebuilt = 0
    for i, vcf in enumerate(present, 1):
        tbi = vcf.with_suffix(vcf.suffix + ".tbi")
        if tbi.exists() and tbi.stat().st_mtime >= vcf.stat().st_mtime:
            skipped += 1
            continue
        force = tbi.exists()
        cmd = [cfg.tools.bcftools, "index", "--tbi"]
        if force:
            cmd.append("-f")
            rebuilt += 1
            log(f"  [{i}/{len(present)}] rebuilding stale index: {vcf.name}")
        else:
            indexed += 1
            log(f"  [{i}/{len(present)}] indexing: {vcf.name}")
        cmd.append(str(vcf))
        run(cmd)
    log(f"  Indexing done: {indexed} new, {rebuilt} rebuilt, {skipped} up-to-date.")


# ── Internals ─────────────────────────────────────────────────────────────────

def _fetch_url(url: str) -> str:
    with urllib.request.urlopen(url, timeout=60) as resp:
        return resp.read().decode("utf-8")


def _fetch_vcf_urls(allowed_accessions: set[str] | None) -> dict[str, tuple[str, str]]:
    """Return {filename: (sample_accession, https_url)} for all VCF analysis files in PRJEB9586."""
    tsv = _fetch_url(_ENA_ANALYSIS_API)
    lines = [l for l in tsv.splitlines() if l.strip()]
    if len(lines) < 2:
        return {}

    header = lines[0].split("\t")
    result: dict[str, tuple[str, str]] = {}
    for line in lines[1:]:
        parts = dict(zip(header, line.split("\t")))
        sample_acc = parts.get("sample_accession", "").strip()
        ftp_field = parts.get("submitted_ftp", "").strip()

        if allowed_accessions is not None and sample_acc not in allowed_accessions:
            continue

        for raw_url in ftp_field.split(";"):
            raw_url = raw_url.strip()
            if not raw_url.endswith(".vcf.gz"):
                continue
            filename = raw_url.rsplit("/", 1)[-1]
            https_url = raw_url.replace("ftp://ftp.sra.ebi.ac.uk", "https://ftp.sra.ebi.ac.uk")
            result[filename] = (sample_acc, https_url)

    return result


def _link_vcfs_to_metadata(cfg: Config, url_map: dict[str, tuple[str, str]]) -> None:
    """Write vcf_filename column into sgdp_samples.tsv, joining on ena_accession == sample_accession."""
    samples_tsv = cfg.sgdp_raw_dir() / "sgdp_samples.tsv"
    if not samples_tsv.exists():
        return

    # {sample_accession: filename}  — last file wins if multiple per sample
    acc_to_file: dict[str, str] = {sample_acc: fn for fn, (sample_acc, _) in url_map.items()}

    meta = pd.read_csv(samples_tsv, sep="\t")
    meta["vcf_filename"] = meta["ena_accession"].map(acc_to_file)
    meta.to_csv(samples_tsv, sep="\t", index=False)

    linked = meta["vcf_filename"].notna().sum()
    log(f"  Linked {linked}/{len(meta)} metadata rows to VCF files.")


def _parse_igsr_response(tsv: str, filter_populations: list[str]) -> pd.DataFrame:
    empty = pd.DataFrame(columns=["ena_accession", "sample_alias", "population", "region", "sex"])
    lines = [l for l in tsv.splitlines() if l.strip()]
    if len(lines) < 2:
        return empty

    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        parts = dict(zip(header, line.split("\t")))
        if "Simons Genome Diversity Project" not in parts.get("Data collections", ""):
            continue

        elastic_ids = parts.get("Population elastic ID", "").split(",")
        pop_names = parts.get("Population name", "").split(",")
        superpop_names = parts.get("Superpopulation name", "").split(",")

        # Pick the SGDP-specific entry (elastic ID ends with "SGDP")
        sgdp_idx = next(
            (i for i, eid in enumerate(elastic_ids) if eid.strip().endswith("SGDP")),
            None,
        )
        population = (
            pop_names[sgdp_idx].strip()
            if sgdp_idx is not None and sgdp_idx < len(pop_names)
            else (pop_names[0].strip() if pop_names else "Unknown")
        )
        region = (
            superpop_names[sgdp_idx].strip().replace(" (SGDP)", "")
            if sgdp_idx is not None and sgdp_idx < len(superpop_names)
            else (superpop_names[0].strip() if superpop_names else "Unknown")
        )

        rows.append({
            "ena_accession": parts.get("Biosample ID", "").strip(),
            "sample_alias": parts.get("Sample name", "").strip(),
            "population": population,
            "region": region,
            "sex": {"M": 1, "F": 2}.get(parts.get("Sex", "").strip(), 0),
        })

    if not rows:
        return empty

    df = pd.DataFrame(rows).drop_duplicates("ena_accession")
    if filter_populations:
        df = df[df["population"].isin(filter_populations)]
    return df
