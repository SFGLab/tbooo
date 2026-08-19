"""Tests for src/tbooo/integrity.py using real bytes on disk (tmp_path).

No mocks: every failure mode here is constructed by actually writing/truncating/
touching files, because integrity.py exists to catch exactly these on-disk
corruption patterns for multi-TB pipeline artifacts.
"""

from __future__ import annotations

import os
import time

import pytest

from tbooo.integrity import (
    _index_fresh,
    bgen_ok,
    ensure,
    ensure_bgen,
    ensure_parquet,
    ensure_plink,
    ensure_table,
    ensure_vcf,
    parquet_ok,
    plink_ok,
    remove,
    symlink_ok,
    table_ok,
    vcf_ok,
)

# A minimal, complete, empty BGZF stream (the 28-byte EOF marker on its own is a
# valid, empty bgzip file). See tbooo.utils._BGZF_EOF.
_BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)


def _touch_older(path, seconds: float = 5) -> None:
    """Set both atime and mtime of `path` to `seconds` in the past."""
    now = time.time()
    os.utime(path, (now - seconds, now - seconds))


# ── vcf_ok ────────────────────────────────────────────────────────────────────

def test_vcf_ok_true_for_intact_bgzip_with_fresh_index(tmp_path):
    vcf = tmp_path / "x.vcf.gz"
    vcf.write_bytes(b"not really compressed but ends right" + _BGZF_EOF)
    tbi = tmp_path / "x.vcf.gz.tbi"
    tbi.write_bytes(b"tbi-index-bytes")
    # Index must be at least as new as the data.
    assert vcf_ok(vcf) is True


def test_vcf_ok_false_missing_file(tmp_path):
    vcf = tmp_path / "missing.vcf.gz"
    assert vcf_ok(vcf) is False


def test_vcf_ok_false_empty_file(tmp_path):
    vcf = tmp_path / "empty.vcf.gz"
    vcf.write_bytes(b"")
    tbi = tmp_path / "empty.vcf.gz.tbi"
    tbi.write_bytes(b"idx")
    assert vcf_ok(vcf) is False


def test_vcf_ok_false_truncated_bgzf_eof(tmp_path):
    vcf = tmp_path / "trunc.vcf.gz"
    # Write an intact-looking file then truncate the last 28 bytes (the EOF marker).
    vcf.write_bytes(b"payload" + _BGZF_EOF)
    data = vcf.read_bytes()
    vcf.write_bytes(data[:-28])
    tbi = tmp_path / "trunc.vcf.gz.tbi"
    tbi.write_bytes(b"idx")
    assert vcf_ok(vcf) is False


def test_vcf_ok_false_missing_index(tmp_path):
    vcf = tmp_path / "noidx.vcf.gz"
    vcf.write_bytes(b"payload" + _BGZF_EOF)
    assert vcf_ok(vcf) is False


def test_vcf_ok_false_stale_index(tmp_path):
    """The important case: index mtime older than data mtime means a run died
    between writing the data file and indexing it."""
    vcf = tmp_path / "stale.vcf.gz"
    tbi = tmp_path / "stale.vcf.gz.tbi"
    tbi.write_bytes(b"old-idx")
    _touch_older(tbi, 60)
    vcf.write_bytes(b"payload" + _BGZF_EOF)  # written (and thus mtime'd) after tbi
    assert vcf_ok(vcf) is False


def test_vcf_ok_true_with_csi_index_instead_of_tbi(tmp_path):
    vcf = tmp_path / "csi.vcf.gz"
    vcf.write_bytes(b"payload" + _BGZF_EOF)
    csi = tmp_path / "csi.vcf.gz.csi"
    csi.write_bytes(b"idx")
    assert vcf_ok(vcf) is True


def test_vcf_ok_true_when_index_mtime_exactly_equals_data_mtime(tmp_path):
    """Boundary case for the '>=' in _index_fresh: on a filesystem with coarse
    mtime resolution, a data file and its index written in the same tick can
    land on IDENTICAL mtimes. That must still count as fresh -- '>' instead of
    '>=' would wrongly call a just-built, perfectly intact artifact stale and
    trigger a needless TB-scale rebuild. Set both mtimes explicitly to one
    shared timestamp with os.utime; do not rely on write ordering, which would
    not exercise the exact-equality boundary at all.
    """
    vcf = tmp_path / "equal.vcf.gz"
    vcf.write_bytes(b"payload" + _BGZF_EOF)
    tbi = tmp_path / "equal.vcf.gz.tbi"
    tbi.write_bytes(b"idx")

    ts = time.time()
    os.utime(vcf, (ts, ts))
    os.utime(tbi, (ts, ts))
    assert vcf.stat().st_mtime == tbi.stat().st_mtime  # sanity: truly equal, not just close

    assert _index_fresh(vcf) is True
    assert vcf_ok(vcf) is True


# ── plink_ok ──────────────────────────────────────────────────────────────────

_PLINK_MAGIC = b"\x6c\x1b\x01"


def _write_plink_set(stem, *, magic=_PLINK_MAGIC, bed_body=b"rest-of-bed-data"):
    (stem.parent / (stem.name + ".bed")).write_bytes(magic + bed_body)
    (stem.parent / (stem.name + ".bim")).write_bytes(b"1\trs1\t0\t100\tA\tG\n")
    (stem.parent / (stem.name + ".fam")).write_bytes(b"fam1 ind1 0 0 1 -9\n")


def test_plink_ok_true_when_all_present_and_valid(tmp_path):
    stem = tmp_path / "geno"
    _write_plink_set(stem)
    assert plink_ok(stem) is True


@pytest.mark.parametrize("ext", [".bed", ".bim", ".fam"])
def test_plink_ok_false_when_one_file_missing(tmp_path, ext):
    stem = tmp_path / "geno"
    _write_plink_set(stem)
    (stem.parent / (stem.name + ext)).unlink()
    assert plink_ok(stem) is False


def test_plink_ok_false_wrong_magic_bytes(tmp_path):
    stem = tmp_path / "geno"
    _write_plink_set(stem, magic=b"\x00\x00\x00")
    assert plink_ok(stem) is False


def test_plink_ok_false_empty_bed(tmp_path):
    stem = tmp_path / "geno"
    _write_plink_set(stem)
    (stem.parent / (stem.name + ".bed")).write_bytes(b"")
    assert plink_ok(stem) is False


# ── bgen_ok ───────────────────────────────────────────────────────────────────

def test_bgen_ok_true_valid(tmp_path):
    stem = tmp_path / "chr1"
    (tmp_path / "chr1.bgen").write_bytes(b"bgen-payload")
    (tmp_path / "chr1.bgen.bgi").write_bytes(b"bgi-index")
    assert bgen_ok(stem) is True


def test_bgen_ok_false_missing_bgi(tmp_path):
    stem = tmp_path / "chr1"
    (tmp_path / "chr1.bgen").write_bytes(b"bgen-payload")
    assert bgen_ok(stem) is False


def test_bgen_ok_false_stale_bgi(tmp_path):
    stem = tmp_path / "chr1"
    bgi = tmp_path / "chr1.bgen.bgi"
    bgi.write_bytes(b"old-index")
    _touch_older(bgi, 60)
    (tmp_path / "chr1.bgen").write_bytes(b"bgen-payload")  # newer than bgi
    assert bgen_ok(stem) is False


def test_bgen_ok_false_missing_bgen(tmp_path):
    stem = tmp_path / "chr1"
    (tmp_path / "chr1.bgen.bgi").write_bytes(b"bgi-index")
    assert bgen_ok(stem) is False


def test_bgen_ok_true_when_bgi_mtime_exactly_equals_bgen_mtime(tmp_path):
    """Same '>=' boundary as vcf_ok/_index_fresh, but for bgen_ok's own inline
    comparison (bgi.stat().st_mtime >= bgen.stat().st_mtime). An index written
    in the same tick as its data file must still be treated as fresh.
    """
    stem = tmp_path / "chr1"
    bgen = tmp_path / "chr1.bgen"
    bgi = tmp_path / "chr1.bgen.bgi"
    bgen.write_bytes(b"bgen-payload")
    bgi.write_bytes(b"bgi-index")

    ts = time.time()
    os.utime(bgen, (ts, ts))
    os.utime(bgi, (ts, ts))
    assert bgen.stat().st_mtime == bgi.stat().st_mtime  # sanity: truly equal

    assert bgen_ok(stem) is True


# ── parquet_ok ────────────────────────────────────────────────────────────────

_PAR1 = b"PAR1"


def test_parquet_ok_true_valid(tmp_path):
    p = tmp_path / "x.parquet"
    p.write_bytes(_PAR1 + b"middle-payload-bytes" + _PAR1)
    assert parquet_ok(p) is True


def test_parquet_ok_false_head_only(tmp_path):
    p = tmp_path / "x.parquet"
    p.write_bytes(_PAR1 + b"middle-payload-bytes-no-tail-magic")
    assert parquet_ok(p) is False


def test_parquet_ok_false_tail_only(tmp_path):
    p = tmp_path / "x.parquet"
    p.write_bytes(b"no-head-magic-middle" + _PAR1)
    assert parquet_ok(p) is False


def test_parquet_ok_false_too_short(tmp_path):
    p = tmp_path / "x.parquet"
    p.write_bytes(b"PAR1PAR")  # 7 bytes, under the 8-byte minimum
    assert parquet_ok(p) is False


def test_parquet_ok_false_missing(tmp_path):
    p = tmp_path / "nope.parquet"
    assert parquet_ok(p) is False


# ── table_ok ──────────────────────────────────────────────────────────────────

def test_table_ok_header_only_passes_min_lines_1(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("col_a\tcol_b\n")
    assert table_ok(p, min_lines=1) is True


def test_table_ok_header_only_fails_min_lines_2(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("col_a\tcol_b\n")
    assert table_ok(p, min_lines=2) is False


def test_table_ok_false_missing(tmp_path):
    p = tmp_path / "nope.tsv"
    assert table_ok(p, min_lines=1) is False


def test_table_ok_false_empty(tmp_path):
    p = tmp_path / "empty.tsv"
    p.write_bytes(b"")
    assert table_ok(p, min_lines=1) is False


def test_table_ok_true_with_enough_lines(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("h1\th2\nr1\tr2\nr3\tr4\n")
    assert table_ok(p, min_lines=3) is True


# ── symlink_ok ────────────────────────────────────────────────────────────────

def test_symlink_ok_true_for_valid_link(tmp_path):
    target = tmp_path / "target.txt"
    target.write_text("hi")
    link = tmp_path / "link.txt"
    link.symlink_to(target)
    assert symlink_ok(link) is True


def test_symlink_ok_false_broken_link(tmp_path):
    target = tmp_path / "gone.txt"
    link = tmp_path / "link.txt"
    link.symlink_to(target)  # target never created -> broken link
    assert symlink_ok(link) is False


def test_symlink_ok_false_regular_file_not_a_symlink(tmp_path):
    p = tmp_path / "plain.txt"
    p.write_text("hi")
    assert symlink_ok(p) is False


def test_symlink_ok_applies_target_check_to_resolved_target(tmp_path):
    target = tmp_path / "target.txt"
    target.write_text("payload-content")
    link = tmp_path / "link.txt"
    link.symlink_to(target)

    seen = []

    def check(resolved):
        seen.append(resolved)
        return resolved.read_text() == "payload-content"

    assert symlink_ok(link, target_check=check) is True
    assert seen == [target.resolve()]

    def fails(resolved):
        return False

    assert symlink_ok(link, target_check=fails) is False


# ── remove ────────────────────────────────────────────────────────────────────

def test_remove_deletes_files(tmp_path):
    a = tmp_path / "a.txt"
    b = tmp_path / "b.txt"
    a.write_text("x")
    b.write_text("y")
    remove(a, b)
    assert not a.exists()
    assert not b.exists()


def test_remove_deletes_symlinks_not_their_targets(tmp_path):
    target = tmp_path / "target.txt"
    target.write_text("keep-me")
    link = tmp_path / "link.txt"
    link.symlink_to(target)
    remove(link)
    assert not link.exists() and not link.is_symlink()
    assert target.exists()  # remove() must not follow the link to delete the target


def test_remove_is_silent_on_missing_paths(tmp_path):
    missing = tmp_path / "never-existed.txt"
    remove(missing)  # must not raise


def test_remove_does_not_delete_directories(tmp_path):
    d = tmp_path / "a_dir"
    d.mkdir()
    remove(d)
    assert d.exists()


# ── ensure ────────────────────────────────────────────────────────────────────

class _Recorder:
    def __init__(self, check_results):
        self.check_results = list(check_results)
        self.check_calls = 0
        self.build_calls = 0
        self.purge_calls = 0

    def check(self):
        self.check_calls += 1
        return self.check_results[min(self.check_calls - 1, len(self.check_results) - 1)]

    def build(self):
        self.build_calls += 1

    def purge(self):
        self.purge_calls += 1


def test_ensure_skips_when_already_valid():
    rec = _Recorder([True])
    ensure("label", check=rec.check, build=rec.build, purge=rec.purge)
    assert rec.build_calls == 0
    assert rec.purge_calls == 0
    assert rec.check_calls == 1


def test_ensure_purges_then_builds_when_invalid():
    # check() is False first (invalid), True after build (rebuilt successfully).
    rec = _Recorder([False, True])
    ensure("label", check=rec.check, build=rec.build, purge=rec.purge)
    assert rec.purge_calls == 1
    assert rec.build_calls == 1
    assert rec.check_calls == 2


def test_ensure_raises_when_still_invalid_after_build():
    rec = _Recorder([False, False])
    with pytest.raises(RuntimeError, match="rebuilt but still failed validation"):
        ensure("mylabel", check=rec.check, build=rec.build, purge=rec.purge)
    assert rec.purge_calls == 1
    assert rec.build_calls == 1


def test_ensure_call_order_is_purge_then_build_then_recheck():
    order = []
    ensure(
        "label",
        check=lambda: (order.append("check"), len(order) > 3)[1],
        build=lambda: order.append("build"),
        purge=lambda: order.append("purge"),
    )
    # first check -> invalid, then purge, then build, then re-check -> valid
    assert order == ["check", "purge", "build", "check"]


# ── typed wrappers ────────────────────────────────────────────────────────────

def test_ensure_vcf_purges_tbi_csi_and_tmp_on_rebuild(tmp_path):
    vcf = tmp_path / "out.vcf.gz"
    tbi = tmp_path / "out.vcf.gz.tbi"
    csi = tmp_path / "out.vcf.gz.csi"
    tmp_extra = tmp_path / "scratch.partial"
    for p in (vcf, tbi, csi, tmp_extra):
        p.write_bytes(b"stale")

    built = []

    def build():
        built.append(True)
        vcf.write_bytes(b"payload" + _BGZF_EOF)
        new_tbi = tmp_path / "out.vcf.gz.tbi"
        new_tbi.write_bytes(b"fresh-idx")

    ensure_vcf("vcf-label", vcf, build, tmp=[tmp_extra])
    assert built == [True]
    # tmp scratch file was purged and never recreated
    assert not tmp_extra.exists()
    # the stale csi was purged (build didn't recreate it)
    assert not csi.exists()
    assert vcf_ok(vcf) is True


def test_ensure_plink_purges_bed_bim_fam_log_on_rebuild(tmp_path):
    stem = tmp_path / "geno"
    # Start invalid (wrong magic) so ensure() actually purges + rebuilds.
    _write_plink_set(stem, magic=b"\x00\x00\x00")
    log_file = tmp_path / "geno.log"
    log_file.write_text("stale log")

    def build():
        _write_plink_set(stem)

    ensure_plink("plink-label", stem, build)
    assert not log_file.exists()
    assert plink_ok(stem) is True


def test_ensure_bgen_purges_bgen_bgi_and_extra_on_rebuild(tmp_path):
    stem = tmp_path / "chr1"
    # Start invalid (no .bgi yet) so ensure() actually purges + rebuilds.
    (tmp_path / "chr1.bgen").write_bytes(b"stale")
    sample = tmp_path / "chr1.sample"
    sample.write_text("stale-sample")

    def build():
        (tmp_path / "chr1.bgen").write_bytes(b"fresh")
        (tmp_path / "chr1.bgen.bgi").write_bytes(b"fresh-idx")

    ensure_bgen("bgen-label", stem, build, extra=[sample])
    assert not sample.exists()
    assert bgen_ok(stem) is True


def test_ensure_parquet_purges_path_on_rebuild(tmp_path):
    p = tmp_path / "x.parquet"
    p.write_bytes(b"garbage-not-parquet")

    def build():
        p.write_bytes(_PAR1 + b"data" + _PAR1)

    ensure_parquet("parquet-label", p, build)
    assert parquet_ok(p) is True


def test_ensure_table_purges_path_on_rebuild(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("")  # empty -> invalid

    def build():
        p.write_text("h1\th2\n")

    ensure_table("table-label", p, build, min_lines=1)
    assert table_ok(p, min_lines=1) is True


# ── ensure: mutation-verified skip contract ─────────────────────────────────

def test_ensure_never_calls_build_when_check_already_true():
    """This is the load-bearing behaviour: skipping an already-valid multi-TB
    artifact must never trigger a rebuild."""

    def check():
        return True

    def build():
        raise AssertionError("build must not be called when check() is already True")

    def purge():
        raise AssertionError("purge must not be called when check() is already True")

    ensure("label", check=check, build=build, purge=purge)
