"""Tests for the subprocess-free functions in src/tbooo/utils.py."""

from __future__ import annotations

import sys

import pytest

from tbooo.utils import check_tools, ensure_dirs, eid_prefix_dir, has_bgzf_eof

_BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)


# ── eid_prefix_dir ────────────────────────────────────────────────────────────

@pytest.mark.parametrize(
    "eid, expected",
    [
        (1000000, "10"),
        (2000000, "20"),
        (6000000, "60"),
        (9999999, "99"),
    ],
)
def test_eid_prefix_dir_returns_first_two_digits(eid, expected):
    assert eid_prefix_dir(eid) == expected


# ── has_bgzf_eof ──────────────────────────────────────────────────────────────

def test_has_bgzf_eof_true_for_intact_file(tmp_path):
    p = tmp_path / "x.vcf.gz"
    p.write_bytes(b"some payload" + _BGZF_EOF)
    assert has_bgzf_eof(p) is True


def test_has_bgzf_eof_false_when_truncated(tmp_path):
    p = tmp_path / "x.vcf.gz"
    full = b"some payload" + _BGZF_EOF
    p.write_bytes(full[:-1])  # drop the last byte of the marker
    assert has_bgzf_eof(p) is False


def test_has_bgzf_eof_false_when_shorter_than_marker(tmp_path):
    p = tmp_path / "tiny.vcf.gz"
    p.write_bytes(b"short")
    assert has_bgzf_eof(p) is False


def test_has_bgzf_eof_false_when_missing(tmp_path):
    p = tmp_path / "nope.vcf.gz"
    assert has_bgzf_eof(p) is False


# ── ensure_dirs ───────────────────────────────────────────────────────────────

def test_ensure_dirs_creates_nested_paths(tmp_path):
    target = tmp_path / "a" / "b" / "c"
    ensure_dirs(target)
    assert target.is_dir()


def test_ensure_dirs_is_idempotent(tmp_path):
    target = tmp_path / "a" / "b"
    ensure_dirs(target)
    ensure_dirs(target)  # must not raise on an existing dir
    assert target.is_dir()


def test_ensure_dirs_accepts_multiple_paths(tmp_path):
    p1 = tmp_path / "one"
    p2 = tmp_path / "two" / "nested"
    ensure_dirs(p1, p2)
    assert p1.is_dir()
    assert p2.is_dir()


# ── check_tools ───────────────────────────────────────────────────────────────

def test_check_tools_exits_when_binary_missing(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(SystemExit) as exc_info:
        check_tools({"bcftools": "bcftools"})
    assert exc_info.value.code == 1


def test_check_tools_returns_normally_when_all_resolve(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: f"/usr/bin/{name}")
    check_tools({"bcftools": "bcftools", "plink2": "plink2"})  # must not raise/exit


def test_check_tools_reports_only_missing_ones(monkeypatch, capsys):
    def fake_which(name):
        return None if name == "missingtool" else f"/usr/bin/{name}"

    monkeypatch.setattr("shutil.which", fake_which)
    with pytest.raises(SystemExit):
        check_tools({"present": "presenttool", "absent": "missingtool"})
    err = capsys.readouterr().err
    assert "absent" in err and "missingtool" in err
    assert "present" not in err.split("missing required tools:")[1]
