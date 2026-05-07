from __future__ import annotations

import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Sequence


def run(
    cmd: str | Sequence[str],
    *,
    check: bool = True,
    capture: bool = False,
    cwd: Path | None = None,
) -> subprocess.CompletedProcess:
    if isinstance(cmd, str):
        cmd = cmd.split()
    kwargs: dict = {"cwd": cwd}
    if capture:
        kwargs["stdout"] = subprocess.PIPE
        kwargs["stderr"] = subprocess.PIPE
    return subprocess.run(list(cmd), check=check, **kwargs)


def check_tools(tools_dict: dict[str, str]) -> None:
    missing: list[str] = []
    for name, binary in tools_dict.items():
        if shutil.which(binary) is None:
            missing.append(f"{name} ({binary})")
    if missing:
        print("ERROR: missing required tools:", file=sys.stderr)
        for m in missing:
            print(f"  {m}", file=sys.stderr)
        sys.exit(1)


def ensure_dirs(*paths: Path) -> None:
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)


def bgzip_tabix(vcf_gz: Path, *, tool_bcftools: str = "bcftools") -> None:
    run([tool_bcftools, "index", "--tbi", str(vcf_gz)])


def wget_download(url: str, dest: Path, *, tool_wget: str = "wget", resume: bool = True) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        tool_wget,
        "--continue" if resume else "--no-continue",
        "--quiet", "--show-progress",
        "--read-timeout=10",   # treat 10s of silence as a stall
        "--tries=1000",          # retry up to 1000 times on stall or transient error
        "--waitretry=5",       # wait 5s before each retry (doubles up to 30s)
        "-O", str(dest), url,
    ]
    run(cmd)


def parallel_download(
    tasks: list[tuple[str, Path, str]],
    workers: int,
) -> None:
    """Download a list of (url, dest, wget_bin) tuples in parallel.

    Skips files that already exist. Re-raises the first exception encountered.
    """
    def _one(url: str, dest: Path, wget_bin: str) -> None:
        if dest.exists():
            log(f"  skip (exists): {dest.name}")
            return
        log(f"  downloading: {dest.name}")
        wget_download(url, dest, tool_wget=wget_bin)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_one, url, dest, wget): dest.name
                   for url, dest, wget in tasks}
        for future in as_completed(futures):
            future.result()


def eid_prefix_dir(eid: int) -> str:
    return str(eid)[:2]


def log(msg: str) -> None:
    print(f"[tbooo] {msg}", flush=True)
