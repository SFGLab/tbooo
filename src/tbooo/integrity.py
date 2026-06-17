"""Cheap structural validation for pipeline artifacts → sanity-checked idempotency.

Re-running `tbooo run` must skip work that is already *done and intact*, and rebuild
only what is missing or broken — never the whole dataset from scratch.

The validators here answer "is this finished and intact?" in milliseconds, without
decompressing or parsing the full file. They catch the failure modes a killed/crashed
run leaves behind: truncated files (no BGZF/PLINK/parquet end-marker), and missing or
stale indexes (a build writes the data file first, then indexes it, so an index older
than its data means the prior run died between the two steps).

Deep end-to-end verification (`bgzip -t`) stays opt-in behind the `--deep-check` flag.

The `ensure_*` helpers wrap a build step so it only runs when the target is invalid,
purge any partial remains first, and re-validate afterwards (a build that produced
garbage fails loudly here instead of being trusted downstream).
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterable

from tbooo.utils import has_bgzf_eof, log

# PLINK 1 .bed magic: 0x6c 0x1b followed by mode byte 0x01 (variant-major).
_PLINK_BED_MAGIC = b"\x6c\x1b\x01"
# Apache Parquet files begin and end with the 4-byte magic "PAR1".
_PARQUET_MAGIC = b"PAR1"


# ── path helpers ────────────────────────────────────────────────────────────────

def _suf(path: Path, suffix: str) -> Path:
    """Append a literal suffix (handles multi-dot names like '.bgen.bgi')."""
    return Path(str(path) + suffix)


def _nonempty(p: Path) -> bool:
    try:
        return p.is_file() and p.stat().st_size > 0
    except OSError:
        return False


def _index_fresh(path: Path) -> bool:
    """True if a tabix/CSI index exists and is at least as new as its data file."""
    try:
        data_mtime = path.stat().st_mtime
    except OSError:
        return False
    for ext in (".tbi", ".csi"):
        idx = _suf(path, ext)
        try:
            if idx.exists() and idx.stat().st_mtime >= data_mtime:
                return True
        except OSError:
            continue
    return False


# ── validators ──────────────────────────────────────────────────────────────────

def vcf_ok(path: Path) -> bool:
    """bgzipped VCF: present, non-empty, intact BGZF stream, fresh index."""
    return _nonempty(path) and has_bgzf_eof(path) and _index_fresh(path)


def plink_ok(stem: Path) -> bool:
    """PLINK1 set: .bed/.bim/.fam all present & non-empty, .bed has its magic bytes."""
    bed, bim, fam = _suf(stem, ".bed"), _suf(stem, ".bim"), _suf(stem, ".fam")
    if not all(_nonempty(p) for p in (bed, bim, fam)):
        return False
    try:
        with bed.open("rb") as f:
            return f.read(3) == _PLINK_BED_MAGIC
    except OSError:
        return False


def bgen_ok(stem: Path) -> bool:
    """BGEN: .bgen present & non-empty with a fresh .bgen.bgi index.

    The .sample file is shared across chromosomes and validated separately.
    """
    bgen, bgi = _suf(stem, ".bgen"), _suf(stem, ".bgen.bgi")
    if not _nonempty(bgen):
        return False
    try:
        return bgi.exists() and bgi.stat().st_mtime >= bgen.stat().st_mtime
    except OSError:
        return False


def parquet_ok(path: Path) -> bool:
    """Parquet: present, non-empty, PAR1 magic at both head and tail."""
    if not _nonempty(path) or path.stat().st_size < 8:
        return False
    try:
        with path.open("rb") as f:
            head = f.read(4)
            f.seek(-4, 2)
            tail = f.read(4)
    except OSError:
        return False
    return head == _PARQUET_MAGIC and tail == _PARQUET_MAGIC


def table_ok(path: Path, *, min_lines: int = 1) -> bool:
    """Text table (TSV/TXT): present, non-empty, at least `min_lines` lines.

    `min_lines=1` accepts a header-only file (e.g. an empty-but-valid ukb_rel.txt).
    """
    if not _nonempty(path):
        return False
    n = 0
    try:
        with path.open("rb") as f:
            for _ in f:
                n += 1
                if n >= min_lines:
                    return True
    except OSError:
        return False
    return n >= min_lines


def symlink_ok(path: Path, target_check: Callable[[Path], bool] | None = None) -> bool:
    """Symlink: is a link, resolves to an existing target, target passes `target_check`."""
    if not path.is_symlink():
        return False
    target = path.resolve()
    if not target.exists():
        return False
    return target_check(target) if target_check is not None else True


# ── purge + ensure ───────────────────────────────────────────────────────────────

def remove(*paths: Path) -> None:
    """Delete files/symlinks if present (best effort; ignores dirs and missing files)."""
    for p in paths:
        try:
            if p.is_symlink() or p.is_file():
                p.unlink()
        except OSError:
            pass


def ensure(
    label: str,
    *,
    check: Callable[[], bool],
    build: Callable[[], None],
    purge: Callable[[], None],
) -> None:
    """Validated-idempotent build step.

    Skip if the artifact is already valid; otherwise purge any partial remains,
    rebuild, and re-validate — raising if the rebuild still fails the check.
    """
    if check():
        log(f"  skip (valid): {label}")
        return
    purge()
    build()
    if not check():
        raise RuntimeError(f"{label}: rebuilt but still failed validation")


# ── typed convenience wrappers ───────────────────────────────────────────────────

def ensure_vcf(label: str, path: Path, build: Callable[[], None],
               *, tmp: Iterable[Path] = ()) -> None:
    extra = list(tmp)
    ensure(
        label,
        check=lambda: vcf_ok(path),
        build=build,
        purge=lambda: remove(path, _suf(path, ".tbi"), _suf(path, ".csi"), *extra),
    )


def ensure_plink(label: str, stem: Path, build: Callable[[], None]) -> None:
    ensure(
        label,
        check=lambda: plink_ok(stem),
        build=build,
        purge=lambda: remove(_suf(stem, ".bed"), _suf(stem, ".bim"),
                             _suf(stem, ".fam"), _suf(stem, ".log")),
    )


def ensure_bgen(label: str, stem: Path, build: Callable[[], None],
                *, extra: Iterable[Path] = ()) -> None:
    extras = list(extra)
    ensure(
        label,
        check=lambda: bgen_ok(stem),
        build=build,
        purge=lambda: remove(_suf(stem, ".bgen"), _suf(stem, ".bgen.bgi"), *extras),
    )


def ensure_parquet(label: str, path: Path, build: Callable[[], None]) -> None:
    ensure(label, check=lambda: parquet_ok(path), build=build,
           purge=lambda: remove(path))


def ensure_table(label: str, path: Path, build: Callable[[], None],
                 *, min_lines: int = 1) -> None:
    ensure(label, check=lambda: table_ok(path, min_lines=min_lines), build=build,
           purge=lambda: remove(path))
