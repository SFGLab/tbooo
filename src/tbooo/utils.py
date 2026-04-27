from __future__ import annotations

import shutil
import subprocess
import sys
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
    cmd = [tool_wget, "--continue" if resume else "--no-continue",
           "--quiet", "--show-progress", "-O", str(dest), url]
    run(cmd)


def eid_prefix_dir(eid: int) -> str:
    return str(eid)[:2]


def log(msg: str) -> None:
    print(f"[tbooo] {msg}", flush=True)
