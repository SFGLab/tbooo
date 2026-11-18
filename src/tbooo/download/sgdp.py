"""Download SGDP VCF files and sample metadata.

The per-sample VCF files are the project's own release, `<illumina id>.annotated.nh2.vcf.gz`,
fetched from the ENA analysis records of PRJEB9586. They hold one record for every base of the
genome, on hs37d5 (GRCh37), with unphased genotypes. `tbooo map sgdp-grch38` turns them into
phased GRCh38 variant calls. SGDP CRAMs (~14 TB) are not downloaded.

The sample table is built from two sources joined on the Illumina id, which is both the VCF
file's prefix and the sample name inside it: the ENA analysis report (which VCF belongs to which
sample) and the project's sample table (population, region, sex). The IGSR portal is not used:
its BioSample ids are not the ids ENA files the VCFs under, so nothing joined.

Outputs:
    data/raw/sgdp/sgdp_samples.tsv
        columns: sample_id, ena_accession, sample_alias, sgdp_id, population, region,
                 country, sex, latitude, longitude, vcf_filename
    data/raw/sgdp/vcf/<illumina id>.annotated.nh2.vcf.gz  (one per sample)
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import pandas as pd

from tbooo.config import Config
from tbooo.utils import (
    ensure_dirs,
    find_corrupt_bgzf,
    has_bgzf_eof,
    log,
    parallel_download,
    run,
    wget_download,
)

# ENA portal API — analysis-level records (processed files, including VCFs)
_ENA_ANALYSIS_API = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession=PRJEB9586"
    "&result=analysis"
    "&fields=analysis_accession,sample_accession,sample_alias,analysis_type,submitted_ftp"
    "&format=tsv"
    "&limit=0"
)

# The project's own sample table: one row per sequenced library, keyed by Illumina id.
_SGDP_TABLE_URL = (
    "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/"
    "SGDP_metadata.279public.21signedLetter.44Fan.samples.txt"
)

SAMPLE_COLUMNS = [
    "sample_id", "ena_accession", "sample_alias", "sgdp_id", "population", "region",
    "country", "sex", "latitude", "longitude", "vcf_filename",
]

# The project table writes regions without spaces; the phenotype stage keys on these names.
_REGION_NAMES = {
    "WestEurasia": "West Eurasia",
    "EastAsia": "East Asia",
    "SouthAsia": "South Asia",
    "CentralAsiaSiberia": "Central Asia / Siberia",
    "America": "Native Americas",
    "Oceania": "Oceania",
    "Africa": "Africa",
}


def download_metadata(cfg: Config) -> None:
    """Build the SGDP sample table from the ENA analysis report and the project's sample table."""
    out = cfg.sgdp_raw_dir() / "sgdp_samples.tsv"
    if out.exists():
        if "sample_id" in pd.read_csv(out, sep="\t", nrows=0).columns:
            log(f"  skip (exists): {out.name}")
            return
        # A table from the IGSR portal, keyed by an accession no VCF carries. Keep it aside.
        stale = out.with_name("sgdp_samples.igsr.tsv")
        out.replace(stale)
        log(f"  {out.name} predates the Illumina-id key; moved to {stale.name} and rebuilt")
    ensure_dirs(cfg.sgdp_raw_dir())

    log("Fetching SGDP VCF records from ENA and the project sample table…")
    analyses = _parse_ena_analyses(_fetch_url(_ENA_ANALYSIS_API))
    try:
        table = _parse_sgdp_table(_fetch_url(_SGDP_TABLE_URL, encoding="latin-1"))
    except Exception as exc:
        log(f"  project sample table unavailable ({exc}); populations will read 'Unknown'")
        table = {}

    df = _build_samples(analyses, table, cfg.sgdp_populations)
    df.to_csv(out, sep="\t", index=False)
    known = int((df["population"] != "Unknown").sum())
    log(f"  wrote {out} ({len(df)} samples, {known} with a population)")


def download_vcfs(
    cfg: Config,
    *,
    deep_check: bool = False,
    tolerate_broken: bool = False,
) -> None:
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

    _repair_truncated_vcfs(
        cfg, url_map, vcf_dir,
        deep_check=deep_check,
        tolerate_broken=tolerate_broken,
    )
    _index_vcfs(cfg, [vcf_dir / fn for fn in url_map])

    log("SGDP VCF download complete.")


def _repair_truncated_vcfs(
    cfg: Config,
    url_map: dict,
    vcf_dir: Path,
    *,
    deep_check: bool = False,
    tolerate_broken: bool = False,
) -> None:
    """Detect broken VCFs and re-fetch them.

    Always runs the cheap BGZF EOF check.
    With `deep_check=True`, additionally runs `bgzip -t` end-to-end on every VCF
    (parallelized across `cfg.deep_check_workers`) to catch mid-stream corruption.
    """
    # Build {path: url} for files that exist on disk
    existing: dict[Path, str] = {}
    for filename, (_, url) in url_map.items():
        path = vcf_dir / filename
        if path.exists():
            existing[path] = url

    log(f"Validating {len(existing)} SGDP VCF(s) — BGZF EOF check…")
    eof_broken = {p: u for p, u in existing.items() if not has_bgzf_eof(p)}
    if eof_broken:
        log(f"  {len(eof_broken)} missing EOF marker (truncated).")

    deep_broken: dict[Path, str] = {}
    if deep_check:
        intact_after_eof = [p for p in existing if p not in eof_broken]
        log(f"Deep-checking {len(intact_after_eof)} VCF(s) with `bgzip -t` "
            f"({cfg.deep_check_workers} workers)…")
        deep_broken_paths = find_corrupt_bgzf(intact_after_eof, cfg.deep_check_workers)
        deep_broken = {p: existing[p] for p in deep_broken_paths}

    broken: dict[Path, str] = {**eof_broken, **deep_broken}
    if not broken:
        log("  All VCFs intact.")
        return

    log(f"  Total {len(broken)} broken VCF(s) "
        f"({len(eof_broken)} truncated, {len(deep_broken)} corrupt) — re-downloading…")
    for path in broken:
        log(f"    removing: {path.name}")
        path.unlink(missing_ok=True)
        tbi = path.with_suffix(path.suffix + ".tbi")
        tbi.unlink(missing_ok=True)

    tasks = [(url, path, cfg.tools.wget) for path, url in broken.items()]
    parallel_download(tasks, cfg.download_workers)

    # Re-validate
    redownloaded = list(broken)
    still_broken_paths: list[Path] = [p for p in redownloaded if not has_bgzf_eof(p)]
    if deep_check:
        still_broken_paths.extend(find_corrupt_bgzf(
            [p for p in redownloaded if has_bgzf_eof(p)],
            cfg.deep_check_workers,
        ))

    if still_broken_paths:
        names = [p.name for p in still_broken_paths]
        sample = ", ".join(names[:5])
        more = f" (+{len(names) - 5} more)" if len(names) > 5 else ""
        msg = f"Re-download failed: {len(names)} VCF(s) still broken: {sample}{more}"
        if not tolerate_broken:
            raise RuntimeError(msg)
        log(f"  WARN: {msg}")
        log(f"  --tolerate-broken active: dropping broken VCF(s); "
            f"affected samples will be absent from downstream outputs.")
        for path in still_broken_paths:
            path.unlink(missing_ok=True)
            tbi = path.with_suffix(path.suffix + ".tbi")
            tbi.unlink(missing_ok=True)
            url_map.pop(path.name, None)

    survivors = len(broken) - len(still_broken_paths)
    log(f"  Re-downloaded {survivors}/{len(broken)} VCF(s) successfully"
        + (f"; {len(still_broken_paths)} dropped." if still_broken_paths else "."))


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

def _fetch_url(url: str, *, encoding: str = "utf-8") -> str:
    with urllib.request.urlopen(url, timeout=120) as resp:
        return resp.read().decode(encoding)


def _parse_ena_analyses(tsv: str) -> list[dict[str, str]]:
    """One row per VCF file of the ENA analysis report.

    Each row has `sample_id` (the Illumina id: the file name up to its first dot, which is
    also the sample name inside the VCF), `ena_accession`, `sample_alias`, `vcf_filename`
    and `url`.
    """
    lines = [l for l in tsv.splitlines() if l.strip()]
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    rows: list[dict[str, str]] = []
    for line in lines[1:]:
        parts = dict(zip(header, line.split("\t")))
        for raw_url in parts.get("submitted_ftp", "").split(";"):
            raw_url = raw_url.strip()
            if not raw_url.endswith(".vcf.gz"):
                continue
            filename = raw_url.rsplit("/", 1)[-1]
            if "://" not in raw_url:
                raw_url = "https://" + raw_url
            rows.append({
                "sample_id": filename.split(".", 1)[0],
                "ena_accession": parts.get("sample_accession", "").strip(),
                "sample_alias": parts.get("sample_alias", "").strip(),
                "vcf_filename": filename,
                "url": raw_url.replace("ftp://ftp.sra.ebi.ac.uk", "https://ftp.sra.ebi.ac.uk"),
            })
    return rows


def _parse_sgdp_table(text: str) -> dict[str, dict[str, str]]:
    """The project's sample table as {Illumina id: row}. Its header line starts with '#'."""
    lines = [l for l in text.splitlines() if l.strip()]
    if len(lines) < 2:
        return {}
    header = lines[0].lstrip("#").split("\t")
    table: dict[str, dict[str, str]] = {}
    for line in lines[1:]:
        row = dict(zip(header, line.split("\t")))
        key = row.get("Illumina_ID", "").strip()
        if key:
            table[key] = row
    return table


def _build_samples(
    analyses: list[dict[str, str]],
    table: dict[str, dict[str, str]],
    filter_populations: list[str],
) -> pd.DataFrame:
    """Join the VCF records to the project table on the Illumina id.

    A VCF whose id the table does not list keeps its ENA alias and reads 'Unknown' for
    population and region. Rows are sorted by sample id, so EIDs assigned from this table
    do not depend on the order ENA returns its report in.
    """
    rows = []
    for a in analyses:
        t = table.get(a["sample_id"], {})
        rows.append({
            "sample_id": a["sample_id"],
            "ena_accession": a["ena_accession"],
            "sample_alias": a["sample_alias"] or t.get("Sample_ID", "").strip(),
            "sgdp_id": t.get("SGDP_ID", "").strip(),
            "population": t.get("Population_ID", "").strip() or "Unknown",
            "region": _REGION_NAMES.get(t.get("Region", "").strip(), t.get("Region", "").strip() or "Unknown"),
            "country": t.get("Country", "").strip(),
            "sex": {"M": 1, "F": 2}.get(t.get("Gender", "").strip(), 0),
            "latitude": t.get("Latitude", "").strip(),
            "longitude": t.get("Longitude", "").strip(),
            "vcf_filename": a["vcf_filename"],
        })
    df = pd.DataFrame(rows, columns=SAMPLE_COLUMNS).drop_duplicates("sample_id")
    if filter_populations:
        df = df[df["population"].isin(filter_populations)]
    return df.sort_values("sample_id").reset_index(drop=True)


def _fetch_vcf_urls(allowed_accessions: set[str] | None) -> dict[str, tuple[str, str]]:
    """Return {filename: (sample_accession, https_url)} for all VCF analysis files in PRJEB9586."""
    result: dict[str, tuple[str, str]] = {}
    for row in _parse_ena_analyses(_fetch_url(_ENA_ANALYSIS_API)):
        if allowed_accessions is not None and row["ena_accession"] not in allowed_accessions:
            continue
        result[row["vcf_filename"]] = (row["ena_accession"], row["url"])
    return result
