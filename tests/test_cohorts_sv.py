from __future__ import annotations

import subprocess

import pytest

from conftest_cohorts import make_sample_dirs_cohort
from tbooo.cohorts.discover import discover_samples
from tbooo.cohorts.eids import assign_eids, write_rename_file
from tbooo.cohorts.sv import (
    emit_sv_script,
    merged_output_path,
    sv_output_path,
    sv_script_path,
)
from tbooo.utils import eid_prefix_dir


@pytest.fixture
def cfg(tmp_path):
    from tbooo.config import Config, Tools

    return Config(
        data_dir=tmp_path / "instance",
        reference_dir=tmp_path / "instance" / "reference",
        tmp_dir=tmp_path / "instance" / "tmp",
        autosomes=[1], sex_chromosomes=[],
        kg_phase3_release_date="", kg_phase3_vcf_version="", kg_phase3_base_url="",
        kg_nygc_date="", kg_nygc_base_url="", sgdp_populations=[],
        kg_eid_start=1_000_000, sgdp_eid_start=2_000_000,
        array_manifest="", array_proxy_maf=0.05, exome_bed="", tools=Tools(),
    )


@pytest.fixture(autouse=True)
def forbid_subprocess(monkeypatch):
    """Hard guarantee: this stage must never execute anything."""
    def boom(*args, **kwargs):
        raise AssertionError("sv stage must not spawn a subprocess on the login node")

    monkeypatch.setattr(subprocess, "run", boom)
    monkeypatch.setattr(subprocess, "Popen", boom)
    monkeypatch.setattr(subprocess, "check_output", boom)


def _emit(cfg, tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    samples = discover_samples(cohort)
    eids = assign_eids(cohort, samples)
    write_rename_file(cfg, cohort, eids)
    return cohort, samples, eids, emit_sv_script(cfg, cohort, samples, eids)


def test_emits_a_script_without_running_anything(cfg, tmp_path):
    cohort, _, _, script = _emit(cfg, tmp_path)
    assert script == sv_script_path(cfg, cohort)
    assert script.exists()


def test_script_has_sbatch_directives(cfg, tmp_path):
    _, _, _, script = _emit(cfg, tmp_path)
    text = script.read_text()
    assert text.startswith("#!/bin/bash")
    assert "#SBATCH --partition=short" in text
    assert "#SBATCH --job-name=" in text
    assert "set -euo pipefail" in text


def test_script_reheaders_bgzips_and_indexes_each_sample(cfg, tmp_path):
    cohort, samples, eids, script = _emit(cfg, tmp_path)
    text = script.read_text()
    for sid in eids:
        eid = eids[sid]
        assert f"{eid}_sv_0_0.vcf.gz" in text
    assert "bcftools reheader" in text
    assert "bgzip" in text
    assert "tabix" in text


def test_script_merges_into_one_cohort_pvcf(cfg, tmp_path):
    cohort, _, _, script = _emit(cfg, tmp_path)
    text = script.read_text()
    assert "bcftools merge" in text
    assert merged_output_path(cfg, cohort).name in text


def test_script_never_writes_to_the_source_tree(cfg, tmp_path):
    cohort, _, _, script = _emit(cfg, tmp_path)
    for line in script.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("#") or ">" not in stripped:
            continue
        redirect_target = stripped.split(">")[-1].strip().strip('"')
        assert str(cohort.raw_root) not in redirect_target
        assert str(cohort.vcf_root) not in redirect_target


def test_cohort_without_sv_modality_emits_nothing(cfg, tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path, name="nosv", with_vcf=False)
    samples = discover_samples(cohort)
    eids = assign_eids(cohort, samples)
    assert emit_sv_script(cfg, cohort, samples, eids) is None


def test_partition_and_time_limit_are_configurable(cfg, tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    samples = discover_samples(cohort)
    eids = assign_eids(cohort, samples)
    write_rename_file(cfg, cohort, eids)  # sv runs after eids
    script = emit_sv_script(
        cfg, cohort, samples, eids, partition="long", time_limit="12:00:00"
    )
    text = script.read_text()
    assert "#SBATCH --partition=long" in text
    assert "#SBATCH --time=12:00:00" in text


# ── ordering: `eids` must have run, and its output must still be current ──────

def test_missing_rename_file_fails_with_an_actionable_error(cfg, tmp_path):
    """Regression: the script interpolated the rename path without checking it
    existed, so the job died on the compute node under `set -e` instead of here."""
    cohort = make_sample_dirs_cohort(tmp_path)
    samples = discover_samples(cohort)
    eids = assign_eids(cohort, samples)
    with pytest.raises(RuntimeError, match="tbooo restricted eids"):
        emit_sv_script(cfg, cohort, samples, eids)
    assert not sv_script_path(cfg, cohort).exists()


def test_stale_rename_file_is_rejected(cfg, tmp_path):
    """Regression: if the sample set changed since `eids` last ran, the on-disk
    rename map carries OLD EIDs while the output paths use NEW ones — producing
    <newEID>_sv_0_0.vcf.gz whose internal sample column says a different EID."""
    cohort = make_sample_dirs_cohort(tmp_path)
    samples = discover_samples(cohort)
    eids = assign_eids(cohort, samples)
    write_rename_file(cfg, cohort, eids)
    shifted = {sid: eid + 1 for sid, eid in eids.items()}  # a sample was added
    with pytest.raises(RuntimeError, match="tbooo restricted eids"):
        emit_sv_script(cfg, cohort, samples, shifted)
    assert not sv_script_path(cfg, cohort).exists()


def test_rename_file_missing_a_sample_is_rejected(cfg, tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    samples = discover_samples(cohort)
    eids = assign_eids(cohort, samples)
    write_rename_file(cfg, cohort, {k: v for k, v in list(eids.items())[:1]})
    with pytest.raises(RuntimeError, match="tbooo restricted eids"):
        emit_sv_script(cfg, cohort, samples, eids)


# ── orphaned SV outputs: real files, so report but never delete ───────────────

def _stale_sv(cfg, eid):
    path = sv_output_path(cfg, eid)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\x1f\x8b")
    return path


def test_retired_eid_sv_outputs_are_reported_never_deleted(cfg, tmp_path, capsys):
    """Regression: _prune_orphaned_links only scans the link layers and only
    removes symlinks, so a stale SV VCF survived under a retired EID — the exact
    misattribution the prune exists to prevent, in the one layer holding real data."""
    cohort, samples, eids, _ = _emit(cfg, tmp_path)
    same_prefix = _stale_sv(cfg, max(eids.values()) + 1)
    other_prefix = _stale_sv(cfg, cohort.eid_base + 100_000)
    tbi = other_prefix.parent / (other_prefix.name + ".tbi")
    tbi.write_bytes(b"\x1f\x8b")
    foreign = _stale_sv(cfg, 6_000_123)  # another cohort's block — never ours

    capsys.readouterr()
    script = emit_sv_script(cfg, cohort, samples, eids)
    out = capsys.readouterr().out

    assert "WARN" in out
    for path in (same_prefix, other_prefix, tbi):
        assert str(path) in out
        assert path.exists(), "real SV data must never be auto-deleted"
    assert str(foreign) not in out
    assert foreign.exists()

    text = script.read_text()
    for path in (same_prefix, other_prefix, tbi):
        assert f'#   rm -f "{path}"' in text
    # the cleanup must be commented out, never executable
    for line in text.splitlines():
        if "rm -f" in line:
            assert line.lstrip().startswith("#")


def test_no_stale_warning_when_every_sv_output_is_current(cfg, tmp_path, capsys):
    cohort, samples, eids, _ = _emit(cfg, tmp_path)
    for eid in eids.values():
        _stale_sv(cfg, eid)  # these are the CURRENT EIDs, so not stale at all
    capsys.readouterr()
    script = emit_sv_script(cfg, cohort, samples, eids)
    assert "WARN" not in capsys.readouterr().out
    assert "rm -f" not in script.read_text()


# ── the stale scan must never leave the instance tree, and never inject shell ──

def test_symlinked_prefix_dir_never_reports_files_outside_the_instance_tree(
    cfg, tmp_path, capsys
):
    """Regression: the scan used iterdir() + is_dir(), both of which follow
    symlinks, so a prefix directory symlinked elsewhere made real files OUTSIDE
    the instance tree show up as `#   rm -f "..."` lines an operator might paste."""
    from tbooo.cohorts.sv import stale_sv_outputs

    cohort, samples, eids, _ = _emit(cfg, tmp_path)

    outside = tmp_path / "elsewhere"
    outside.mkdir()
    escaped = outside / "7000099_sv_0_0.vcf.gz"
    escaped.write_bytes(b"\x1f\x8b")

    link = cfg.sv_dir() / "70"
    assert not link.exists()
    link.symlink_to(outside, target_is_directory=True)

    # A genuinely stale file inside the tree, so the fix cannot just report nothing.
    genuine = _stale_sv(cfg, cohort.eid_base + 100_000)

    assert stale_sv_outputs(cfg, cohort, eids) == [genuine]

    capsys.readouterr()
    script = emit_sv_script(cfg, cohort, samples, eids)
    out = capsys.readouterr().out
    text = script.read_text()

    assert str(escaped) not in out
    assert str(escaped) not in text
    assert escaped.name not in text
    assert escaped.exists()
    assert str(genuine) in text


def test_symlinked_file_pointing_outside_the_tree_is_not_reported(cfg, tmp_path):
    """Same escape, one level down: a link inside a real prefix directory."""
    from tbooo.cohorts.sv import stale_sv_outputs

    cohort, _, eids, _ = _emit(cfg, tmp_path)
    outside = tmp_path / "elsewhere2"
    outside.mkdir()
    escaped = outside / "7000098_sv_0_0.vcf.gz"
    escaped.write_bytes(b"\x1f\x8b")

    prefix = cfg.sv_dir() / "70"
    prefix.mkdir(parents=True, exist_ok=True)
    (prefix / "7000098_sv_0_0.vcf.gz").symlink_to(escaped)

    assert stale_sv_outputs(cfg, cohort, eids) == []


def test_a_newline_in_a_directory_name_cannot_inject_a_live_rm(cfg, tmp_path):
    """Regression: a prefix directory named "71\\nrm -f PWNED" put the line
    `rm -f PWNED/7000098_sv_0_0.vcf.gz"` — uncommented — into a script the
    operator is told to `sbatch`. The path is disk-derived, so no test input
    guards it; it must be refused before it is interpolated anywhere."""
    cohort, samples, eids, _ = _emit(cfg, tmp_path)
    script_path = sv_script_path(cfg, cohort)
    script_path.unlink()

    evil = cfg.sv_dir() / "71\nrm -f PWNED"
    evil.mkdir(parents=True)
    (evil / "7000098_sv_0_0.vcf.gz").write_bytes(b"\x1f\x8b")

    with pytest.raises(RuntimeError, match="newline|carriage return|quote|unsafe"):
        emit_sv_script(cfg, cohort, samples, eids)

    assert not script_path.exists(), "no script may be written from an unsafe path"


def test_a_quote_in_a_directory_name_cannot_break_out_of_the_quoted_path(cfg, tmp_path):
    cohort, samples, eids, _ = _emit(cfg, tmp_path)
    script_path = sv_script_path(cfg, cohort)
    script_path.unlink()

    # No '/' in the name: it must be one directory, not a nested path.
    evil = cfg.sv_dir() / '72" ; rm -f PWNED ; echo "'
    evil.mkdir(parents=True)
    (evil / "7000097_sv_0_0.vcf.gz").write_bytes(b"\x1f\x8b")

    with pytest.raises(RuntimeError, match="newline|carriage return|quote|unsafe"):
        emit_sv_script(cfg, cohort, samples, eids)
    assert not script_path.exists()


def test_no_generated_line_can_ever_be_a_live_rm(cfg, tmp_path):
    """The general invariant, asserted over whatever the emitter does produce."""
    cohort, samples, eids, script = _emit(cfg, tmp_path)
    _stale_sv(cfg, max(eids.values()) + 1)
    script = emit_sv_script(cfg, cohort, samples, eids)
    for line in script.read_text().splitlines():
        assert not line.lstrip().startswith("rm ")
        if "rm -f" in line:
            assert line.lstrip().startswith("#")
