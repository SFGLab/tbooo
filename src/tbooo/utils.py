from __future__ import annotations

import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Sequence


class CommandError(subprocess.CalledProcessError):
    """CalledProcessError that surfaces the captured stderr tail in its message.

    Subclasses CalledProcessError so existing `except subprocess.CalledProcessError`
    handlers keep working, but the traceback now shows *why* the command failed
    instead of a bare exit status.
    """

    def __str__(self) -> str:
        base = super().__str__()
        if self.stderr:
            tail = "".join(str(self.stderr).splitlines(keepends=True)[-30:])
            return f"{base}\n--- stderr (last 30 lines) ---\n{tail}"
        return base


def run(
    cmd: str | Sequence[str],
    *,
    check: bool = True,
    capture: bool = False,
    cwd: Path | None = None,
) -> subprocess.CompletedProcess:
    if isinstance(cmd, str):
        cmd = cmd.split()
    cmd = list(cmd)

    if capture:
        # Caller wants the output back, not on the console — full capture, no echo.
        return subprocess.run(
            cmd, check=check, cwd=cwd,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )

    # Default: stream output live (stdout inherits the terminal; stderr is teed so
    # progress bars using '\r' still render) AND capture stderr, so that on failure
    # the raised error carries the real diagnostic instead of just "exit status N".
    proc = subprocess.Popen(cmd, cwd=cwd, stderr=subprocess.PIPE)
    err = bytearray()
    assert proc.stderr is not None
    err_fd = proc.stderr.fileno()
    while True:
        chunk = os.read(err_fd, 65536)
        if not chunk:
            break
        sys.stderr.buffer.write(chunk)
        sys.stderr.buffer.flush()
        err.extend(chunk)
    proc.stderr.close()
    returncode = proc.wait()
    stderr_text = err.decode("utf-8", errors="replace")
    if check and returncode != 0:
        raise CommandError(returncode, cmd, output=None, stderr=stderr_text)
    return subprocess.CompletedProcess(cmd, returncode, None, stderr_text)


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
        "--tries=100000",          # retry up to 100000 times on stall or transient error
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


# Empty BGZF block — every valid bgzip-compressed file ends with these 28 bytes.
# Absence signals the download was truncated mid-stream.
_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def has_bgzf_eof(path: Path) -> bool:
    """True iff `path` exists and ends with the BGZF EOF marker (intact bgzip stream)."""
    try:
        if path.stat().st_size < len(_BGZF_EOF):
            return False
        with path.open("rb") as f:
            f.seek(-len(_BGZF_EOF), 2)
            return f.read() == _BGZF_EOF
    except OSError:
        return False


def bgzip_test(path: Path) -> bool | None:
    """End-to-end BGZF integrity check via `bgzip -t`.

    Returns True if intact, False if corrupt, None if `bgzip` is not on PATH.
    Catches mid-stream corruption that `has_bgzf_eof` misses.
    """
    try:
        proc = subprocess.run(
            ["bgzip", "-t", str(path)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    except FileNotFoundError:
        return None
    return proc.returncode == 0


def find_corrupt_bgzf(paths: list[Path], workers: int) -> list[Path]:
    """Run `bgzip -t` on every path in parallel; return those that fail (mid-stream corruption).

    Logs progress as `[i/N]`. If `bgzip` is missing, logs a warning once and returns [].
    """
    if not paths:
        return []
    corrupt: list[Path] = []
    bgzip_missing = False
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(bgzip_test, p): p for p in paths}
        for future in as_completed(futures):
            p = futures[future]
            result = future.result()
            done += 1
            if result is None:
                bgzip_missing = True
                continue
            if not result:
                corrupt.append(p)
                log(f"  [{done}/{len(paths)}] CORRUPT: {p.name}")
            elif done % 25 == 0 or done == len(paths):
                log(f"  [{done}/{len(paths)}] deep-checked")
    if bgzip_missing:
        log("  WARN: bgzip not on PATH — deep integrity check skipped")
    return corrupt


def log(msg: str) -> None:
    print(f"[tbooo] {msg}", flush=True)
