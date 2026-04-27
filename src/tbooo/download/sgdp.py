"""Download SGDP CRAM files from ENA using the FTP pointers file."""

from __future__ import annotations

import concurrent.futures
import urllib.request
from pathlib import Path

from tbooo.config import Config
from tbooo.utils import ensure_dirs, log, run, wget_download


def download_cramps(cfg: Config, populations: list[str], workers: int) -> None:
    """Download SGDP CRAM + index files listed in the ENA pointers file."""
    out = cfg.sgdp_raw_dir()
    ensure_dirs(out)

    pointers_path = out / "ena.ftp.pointers.txt"
    _fetch_pointers(cfg, pointers_path)

    entries = _parse_pointers(pointers_path, populations)
    log(f"SGDP: {len(entries)} samples to download (workers={workers})")

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(_download_sample, entry, out, cfg.tools.wget): entry
            for entry in entries
        }
        for fut in concurrent.futures.as_completed(futures):
            entry = futures[fut]
            try:
                fut.result()
            except Exception as exc:
                log(f"  ERROR downloading {entry['sample']}: {exc}")

    log(f"SGDP download complete → {out}")


# ── Internals ─────────────────────────────────────────────────────────────────

def _fetch_pointers(cfg: Config, dest: Path) -> None:
    if dest.exists():
        log(f"  skip (exists): {dest.name}")
        return
    log(f"Fetching SGDP ENA pointers file…")
    wget_download(cfg.sgdp_ena_pointers_url, dest, tool_wget=cfg.tools.wget)


def _parse_pointers(path: Path, populations: list[str]) -> list[dict]:
    """Return list of {sample, population, cram_url, crai_url} dicts."""
    entries: list[dict] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Expected columns: sample_name  population  ftp_url  md5  ...
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            sample, population, url = parts[0], parts[1], parts[2]
            if populations and population not in populations:
                continue
            if not url.endswith(".cram"):
                continue
            entries.append({
                "sample": sample,
                "population": population,
                "cram_url": url,
                "crai_url": url + ".crai",
            })
    return entries


def _download_sample(entry: dict, out: Path, wget_bin: str) -> None:
    sample = entry["sample"]
    cram_url = entry["cram_url"]
    crai_url = entry["crai_url"]

    # preserve original filename from URL
    cram_name = cram_url.rsplit("/", 1)[-1]
    crai_name = crai_url.rsplit("/", 1)[-1]

    cram_path = out / cram_name
    crai_path = out / crai_name

    if not cram_path.exists():
        log(f"  [{sample}] downloading CRAM…")
        wget_download(cram_url, cram_path, tool_wget=wget_bin)
    if not crai_path.exists():
        log(f"  [{sample}] downloading CRAI…")
        wget_download(crai_url, crai_path, tool_wget=wget_bin)
