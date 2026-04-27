"""Download SGDP VCF files and sample metadata from ENA.

SGDP raw CRAMs (~14 TB) are not downloaded. Instead we fetch:
  - Sample metadata (ENA sample API, ~KB)
  - Per-sample phased VCF files from ENA analysis results (~few hundred GB)

Outputs:
    data/raw/sgdp/sgdp_samples.tsv
        columns: ena_accession, sample_alias, population, region, sex
    data/raw/sgdp/vcf/<accession>.vcf.gz  (one per sample)
"""

from __future__ import annotations

import urllib.request
from pathlib import Path

import pandas as pd

from tbooo.config import Config
from tbooo.utils import ensure_dirs, log, parallel_download, run, wget_download

# ENA portal search API — sample-level metadata for SGDP study PRJEB9586
# The /search endpoint supports result=sample; /filereport is file-only and rejects it.
# `sex` is not a standard ENA sample field — we omit it and default to 0 (unknown).
_ENA_SAMPLE_API = (
    "https://www.ebi.ac.uk/ena/portal/api/search"
    "?query=study_accession%3DPRJEB9586"
    "&result=sample"
    "&fields=sample_accession,sample_alias,sample_title"
    "&format=tsv"
    "&limit=0"
)

# ENA portal API — analysis-level records (processed files, including VCFs)
_ENA_ANALYSIS_API = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession=PRJEB9586"
    "&result=analysis"
    "&fields=analysis_accession,sample_accession,submitted_ftp"
    "&format=tsv"
    "&limit=0"
)

# Mapping of known SGDP population name fragments → continental region
# (best-effort; full table is in the SGDP Nature 2016 supplementary)
_REGION_MAP: dict[str, str] = {
    # West Eurasia
    "Greek": "West Eurasia", "French": "West Eurasia", "Sardinian": "West Eurasia",
    "Spanish": "West Eurasia", "English": "West Eurasia", "Scottish": "West Eurasia",
    "Basque": "West Eurasia", "Italian": "West Eurasia", "Tuscan": "West Eurasia",
    "Finnish": "West Eurasia", "Norwegian": "West Eurasia", "Estonian": "West Eurasia",
    "Armenian": "West Eurasia", "Georgian": "West Eurasia", "Turkish": "West Eurasia",
    "Iranian": "West Eurasia", "Druze": "West Eurasia", "Palestinian": "West Eurasia",
    "Bedouin": "West Eurasia", "Maltese": "West Eurasia", "Cypriot": "West Eurasia",
    # Africa
    "Yoruba": "Africa", "Mandinka": "Africa", "Zulu": "Africa", "Ju_hoan": "Africa",
    "Dinka": "Africa", "Luo": "Africa", "Esan": "Africa", "Mende": "Africa",
    "Gambian": "Africa", "Somali": "Africa", "Ethiopian": "Africa", "Masai": "Africa",
    "Hadza": "Africa", "Sandawe": "Africa", "BantuKenya": "Africa", "BantuSA": "Africa",
    # East Asia
    "Han": "East Asia", "Japanese": "East Asia", "Korean": "East Asia",
    "Dai": "East Asia", "Vietnamese": "East Asia", "Cambodian": "East Asia",
    "Mongolian": "East Asia", "She": "East Asia", "Miao": "East Asia",
    # South Asia
    "Bengali": "South Asia", "Punjabi": "South Asia", "Tamil": "South Asia",
    "Telugu": "South Asia", "Sindhi": "South Asia", "Brahui": "South Asia",
    "Burusho": "South Asia", "Hazara": "South Asia", "Kalash": "South Asia",
    "Pathan": "South Asia", "Balochi": "South Asia",
    # Central Asia / Siberia
    "Buryat": "Central Asia / Siberia", "Yakut": "Central Asia / Siberia",
    "Nganasan": "Central Asia / Siberia", "Selkup": "Central Asia / Siberia",
    "Kazakh": "Central Asia / Siberia", "Kyrgyz": "Central Asia / Siberia",
    "Tuvinian": "Central Asia / Siberia", "Eskimo": "Central Asia / Siberia",
    # Oceania
    "Papuan": "Oceania", "Australian": "Oceania", "Bougainville": "Oceania",
    # Native Americas
    "Maya": "Native Americas", "Quechua": "Native Americas", "Aymara": "Native Americas",
    "Mixtec": "Native Americas", "Zapotec": "Native Americas", "Piapoco": "Native Americas",
    "Surui": "Native Americas", "Karitiana": "Native Americas",
}


def download_metadata(cfg: Config) -> None:
    """Fetch SGDP sample metadata from the ENA portal API."""
    out = cfg.sgdp_raw_dir() / "sgdp_samples.tsv"
    if out.exists():
        log(f"  skip (exists): {out.name}")
        return
    ensure_dirs(cfg.sgdp_raw_dir())

    log("Fetching SGDP sample metadata from ENA…")
    try:
        tsv = _fetch_url(_ENA_SAMPLE_API)
    except Exception as exc:
        log(f"  ENA API failed: {exc}")
        log("  Writing empty metadata file — SGDP rows will be absent from phenotype table.")
        out.write_text("ena_accession\tsample_alias\tpopulation\tregion\tsex\n")
        return

    df = _parse_sample_response(tsv, cfg.sgdp_populations)
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
    tasks = [(url, vcf_dir / filename, cfg.tools.wget) for filename, url in url_map.items()]
    parallel_download(tasks, cfg.download_workers)

    for filename in url_map:
        dest = vcf_dir / filename
        if dest.exists():
            run([cfg.tools.bcftools, "index", "--tbi", str(dest)])

    log("SGDP VCF download complete.")


# ── Internals ─────────────────────────────────────────────────────────────────

def _fetch_url(url: str) -> str:
    with urllib.request.urlopen(url, timeout=60) as resp:
        return resp.read().decode("utf-8")


def _fetch_vcf_urls(allowed_accessions: set[str] | None) -> dict[str, str]:
    """Return {filename: https_url} for all VCF analysis files in PRJEB9586."""
    tsv = _fetch_url(_ENA_ANALYSIS_API)
    lines = [l for l in tsv.splitlines() if l.strip()]
    if len(lines) < 2:
        return {}

    header = lines[0].split("\t")
    result: dict[str, str] = {}
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
            # Prefer HTTPS over FTP for wget compatibility
            https_url = raw_url.replace("ftp://ftp.sra.ebi.ac.uk", "https://ftp.sra.ebi.ac.uk")
            result[filename] = https_url

    return result


def _parse_sample_response(tsv: str, filter_populations: list[str]) -> pd.DataFrame:
    lines = [l for l in tsv.splitlines() if l.strip()]
    if not lines:
        return pd.DataFrame(columns=["ena_accession", "sample_alias", "population", "region", "sex"])

    header = lines[0].split("\t")
    rows = [dict(zip(header, l.split("\t"))) for l in lines[1:]]
    df = pd.DataFrame(rows)

    df = df.rename(columns={
        "sample_accession": "ena_accession",
        "sample_title": "population_raw",
    })

    # Population: SGDP aliases are formatted as "<Population>_<ID>" (e.g. "French_B_French-1")
    if "sample_alias" in df.columns:
        df["population"] = df["sample_alias"].str.extract(r"^([A-Za-z_]+)", expand=False)
    elif "population_raw" in df.columns:
        df["population"] = df["population_raw"]
    else:
        df["population"] = "Unknown"

    df["region"] = df["population"].apply(_infer_region)
    df["sex"] = 0  # ENA search does not expose sex as a standard field

    if filter_populations:
        df = df[df["population"].isin(filter_populations)]

    return df[["ena_accession", "sample_alias", "population", "region", "sex"]].drop_duplicates("ena_accession")


def _infer_region(population: str) -> str:
    for fragment, region in _REGION_MAP.items():
        if fragment.lower() in population.lower():
            return region
    return "Unknown"
