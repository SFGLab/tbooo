from __future__ import annotations

import dataclasses
import shutil
from pathlib import Path

import pytest

from conftest_cohorts import make_flat_bam_cohort, make_sample_dirs_cohort
from tbooo.cohorts.discover import discover_samples
from tbooo.cohorts.eids import assign_eids
from tbooo.cohorts.link import MANIFEST_COLUMNS, link_cohort, link_name, manifest_path
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


def _build(cfg, cohort):
    samples = discover_samples(cohort)
    eids = assign_eids(cohort, samples)
    return samples, eids, link_cohort(cfg, cohort, samples, eids)


def test_bam_links_use_ukb_field_24048_naming(cfg, tmp_path):
    cohort = make_flat_bam_cohort(tmp_path)
    _, eids, _ = _build(cfg, cohort)
    eid = eids[sorted(eids)[0]]
    link = cfg.wgs_dir() / eid_prefix_dir(eid) / f"{eid}_24048_0_0.bam"
    assert link.is_symlink() and link.resolve().exists()
    assert (cfg.wgs_dir() / eid_prefix_dir(eid) / f"{eid}_24048_0_0.bam.bai").is_symlink()


def test_fastq_links_carry_run_and_mate(cfg, tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    _, eids, _ = _build(cfg, cohort)
    eid = eids["p01"]
    d = cfg.raw_reads_dir() / eid_prefix_dir(eid)
    assert (d / f"{eid}_raw_AA1234_L01_7_1.fq.gz").is_symlink()
    assert (d / f"{eid}_raw_AA1234_L01_7_2.fq.gz").is_symlink()


def test_multi_lane_sample_gets_one_link_per_file(cfg, tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    _, eids, _ = _build(cfg, cohort)
    eid = eids["p02"]
    d = cfg.raw_reads_dir() / eid_prefix_dir(eid)
    assert len(list(d.glob(f"{eid}_raw_*"))) == 4


def test_sv_vcfs_are_not_symlinked(cfg, tmp_path):
    """SV VCFs are rewritten by the sv stage, never linked (they carry real IDs)."""
    cohort = make_sample_dirs_cohort(tmp_path)
    _build(cfg, cohort)
    assert not cfg.sv_dir().exists() or not list(cfg.sv_dir().rglob("*.vcf*"))


def test_manifest_records_provenance(cfg, tmp_path):
    cohort = make_flat_bam_cohort(tmp_path)
    _, eids, mpath = _build(cfg, cohort)
    assert mpath == manifest_path(cfg, cohort)
    lines = mpath.read_text().splitlines()
    assert lines[0].split("\t") == list(MANIFEST_COLUMNS)
    row = dict(zip(MANIFEST_COLUMNS, lines[1].split("\t")))
    assert row["source_path"].startswith(str(cohort.raw_root))
    assert int(row["source_size"]) > 0
    assert row["modality"] in {"bam", "bai"}
    assert int(row["eid"]) in eids.values()


def test_sources_are_never_modified(cfg, tmp_path):
    cohort = make_flat_bam_cohort(tmp_path)
    before = {
        p: (p.stat().st_size, p.stat().st_mtime_ns)
        for p in cohort.raw_root.iterdir()
    }
    _build(cfg, cohort)
    after = {
        p: (p.stat().st_size, p.stat().st_mtime_ns)
        for p in cohort.raw_root.iterdir()
    }
    assert before == after


def test_rerun_is_idempotent(cfg, tmp_path):
    cohort = make_flat_bam_cohort(tmp_path)
    samples, eids, mpath = _build(cfg, cohort)
    first = mpath.read_text()
    link_cohort(cfg, cohort, samples, eids)
    assert mpath.read_text() == first


def test_broken_link_is_rebuilt(cfg, tmp_path):
    cohort = make_flat_bam_cohort(tmp_path)
    samples, eids, _ = _build(cfg, cohort)
    eid = eids[sorted(eids)[0]]
    link = cfg.wgs_dir() / eid_prefix_dir(eid) / f"{eid}_24048_0_0.bam"
    link.unlink()
    link.symlink_to(tmp_path / "gone.bam")
    link_cohort(cfg, cohort, samples, eids)
    assert link.resolve().exists()


def _manifest_rows(mpath):
    return {
        row[6]: row  # source_path -> full row
        for row in (
            line.split("\t") for line in mpath.read_text().splitlines()[1:]
        )
    }


def test_changed_source_size_is_detected(cfg, tmp_path):
    """The manifest must be REBUILT on every run, not carried over: a source that
    grew between two runs has to show its new size, not the size recorded first."""
    cohort = make_flat_bam_cohort(tmp_path)
    samples, eids, mpath = _build(cfg, cohort)
    target = sorted(cohort.raw_root.glob("*.bam"))[0]
    before = _manifest_rows(mpath)[str(target)]
    old_size = int(before[4])

    target.write_bytes(b"BAM\x01" + b"EXTRA" * 40)
    new_size = target.stat().st_size
    assert new_size != old_size, "fixture failed to change the source size"

    link_cohort(cfg, cohort, samples, eids)
    after = _manifest_rows(mpath)[str(target)]
    assert int(after[4]) == new_size
    assert int(after[4]) != old_size


def test_prune_is_refused_when_discovery_yields_no_files(cfg, tmp_path):
    """Regression: samples present but zero manifest rows means discovery is
    broken, not that the cohort emptied. Pruning there wiped every link in the
    cohort's EID block and rewrote the manifest empty, with exit code 0."""
    cohort = make_flat_bam_cohort(tmp_path, n=3)
    samples, eids, mpath = _build(cfg, cohort)
    links = sorted(cfg.wgs_dir().rglob("*_24048_0_0.bam*"))
    assert links
    manifest_before = mpath.read_text()

    blind = [dataclasses.replace(s, files=()) for s in samples]
    with pytest.raises(RuntimeError, match="discovery"):
        link_cohort(cfg, cohort, blind, eids)

    assert all(p.is_symlink() for p in links)
    assert mpath.read_text() == manifest_before


def test_orphaned_eid_link_is_pruned_when_sample_count_shrinks(cfg, tmp_path):
    """Reviewer repro: assign_eids renumbers sequentially, so deleting the
    lowest-sorted sample shifts every later EID down by one. The link that used
    to sit at the old highest EID must not survive — otherwise it silently
    carries a duplicate genome under an EID that appears in no manifest.
    """
    cohort = make_flat_bam_cohort(tmp_path, n=3)
    samples1 = discover_samples(cohort)
    eids1 = assign_eids(cohort, samples1)
    link_cohort(cfg, cohort, samples1, eids1)

    highest_eid_run1 = max(eids1.values())
    orphan_link = (
        cfg.wgs_dir() / eid_prefix_dir(highest_eid_run1) / f"{highest_eid_run1}_24048_0_0.bam"
    )
    orphan_bai = Path(str(orphan_link) + ".bai")
    assert orphan_link.is_symlink()

    # Delete the sample that held the lowest EID — every later EID renumbers down.
    first_id = sorted(s.sample_id for s in samples1)[0]
    (cohort.raw_root / f"{first_id}.bam").unlink()
    (cohort.raw_root / f"{first_id}.bam.bai").unlink()

    samples2 = discover_samples(cohort)
    eids2 = assign_eids(cohort, samples2)
    assert highest_eid_run1 not in eids2.values()  # confirms the repro set up correctly
    mpath = link_cohort(cfg, cohort, samples2, eids2)

    assert not orphan_link.exists() and not orphan_link.is_symlink()
    assert not orphan_bai.exists() and not orphan_bai.is_symlink()

    manifest_eids = {line.split("\t")[0] for line in mpath.read_text().splitlines()[1:]}
    assert str(highest_eid_run1) not in manifest_eids


def test_duplicate_destination_within_one_run_raises(cfg, tmp_path):
    """Two SourceFiles resolving to the same link path must never silently
    last-write-wins — that would mix two samples' genomes under one EID.
    Unreachable via discovery today (it guards FASTQ collisions upstream), but
    this is exactly the failure mode the module must never allow through.
    """
    cohort = make_flat_bam_cohort(tmp_path, n=2)
    samples = discover_samples(cohort)
    # Force the collision directly: map both samples onto the same EID.
    eid = cohort.eid_base
    eids = {s.sample_id: eid for s in samples}

    with pytest.raises(RuntimeError, match="duplicate destination"):
        link_cohort(cfg, cohort, samples, eids)


def test_duplicate_destination_error_names_both_sources(cfg, tmp_path):
    cohort = make_flat_bam_cohort(tmp_path, n=2)
    samples = discover_samples(cohort)
    eid = cohort.eid_base
    eids = {s.sample_id: eid for s in samples}
    bam_paths = [
        str(f.path) for s in samples for f in s.files if f.modality == "bam"
    ]

    with pytest.raises(RuntimeError) as excinfo:
        link_cohort(cfg, cohort, samples, eids)

    message = str(excinfo.value)
    assert all(p in message for p in bam_paths)


def test_real_file_at_destination_is_never_clobbered(cfg, tmp_path):
    """The data-safety guard: if a real (non-symlink) file already occupies a
    link's destination, link_cohort must raise, not overwrite it.
    """
    cohort = make_flat_bam_cohort(tmp_path, n=1)
    samples = discover_samples(cohort)
    eids = assign_eids(cohort, samples)
    eid = next(iter(eids.values()))

    dest_dir = cfg.wgs_dir() / eid_prefix_dir(eid)
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{eid}_24048_0_0.bam"
    original = b"NOT A SYMLINK -- pre-existing real data, must not be touched"
    dest.write_bytes(original)

    with pytest.raises(RuntimeError):
        link_cohort(cfg, cohort, samples, eids)

    assert not dest.is_symlink()
    assert dest.read_bytes() == original


def test_data_dir_overlapping_raw_root_is_rejected(tmp_path):
    from tbooo.config import Config, Tools

    cohort = make_flat_bam_cohort(tmp_path)
    samples = discover_samples(cohort)
    eids = assign_eids(cohort, samples)

    # data_dir nested inside the source tree.
    bad_cfg = Config(
        data_dir=cohort.raw_root / "instance",
        reference_dir=cohort.raw_root / "instance" / "reference",
        tmp_dir=cohort.raw_root / "instance" / "tmp",
        autosomes=[1], sex_chromosomes=[],
        kg_phase3_release_date="", kg_phase3_vcf_version="", kg_phase3_base_url="",
        kg_nygc_date="", kg_nygc_base_url="", sgdp_populations=[],
        kg_eid_start=1_000_000, sgdp_eid_start=2_000_000,
        array_manifest="", array_proxy_maf=0.05, exome_bed="", tools=Tools(),
    )
    with pytest.raises(RuntimeError, match="overlap"):
        link_cohort(bad_cfg, cohort, samples, eids)

    # The reverse direction: raw_root nested inside data_dir.
    other_cfg = Config(
        data_dir=tmp_path / "instance2",
        reference_dir=tmp_path / "instance2" / "reference",
        tmp_dir=tmp_path / "instance2" / "tmp",
        autosomes=[1], sex_chromosomes=[],
        kg_phase3_release_date="", kg_phase3_vcf_version="", kg_phase3_base_url="",
        kg_nygc_date="", kg_nygc_base_url="", sgdp_populations=[],
        kg_eid_start=1_000_000, sgdp_eid_start=2_000_000,
        array_manifest="", array_proxy_maf=0.05, exome_bed="", tools=Tools(),
    )
    nested_cohort = make_flat_bam_cohort(tmp_path / "instance2" / "nested", name="nested")
    nested_samples = discover_samples(nested_cohort)
    nested_eids = assign_eids(nested_cohort, nested_samples)
    with pytest.raises(RuntimeError, match="overlap"):
        link_cohort(other_cfg, nested_cohort, nested_samples, nested_eids)


def test_link_name_bam_and_bai():
    assert link_name("bam", 6_000_000) == "6000000_24048_0_0.bam"
    assert link_name("bai", 6_000_000) == "6000000_24048_0_0.bam.bai"


def test_link_name_fastq_includes_run_mate_and_suffix():
    assert (
        link_name("fastq_r1", 6_000_000, run="AA1_L01", suffix=".fq.gz")
        == "6000000_raw_AA1_L01_1.fq.gz"
    )
    assert (
        link_name("fastq_r2", 6_000_000, run="AA1_L01", suffix=".fastq.gz")
        == "6000000_raw_AA1_L01_2.fastq.gz"
    )


def test_link_name_rejects_non_linkable_modality():
    with pytest.raises(ValueError, match="not linkable"):
        link_name("sv_vcf", 6_000_000)


def test_pruning_does_not_touch_another_cohorts_links_in_same_prefix_bucket(cfg, tmp_path):
    """Reviewer repro: cohorts based at 6,000,000 and 6,050,000 both map to
    EID-prefix directory "60" (a 2-digit prefix is a 100,000-wide bucket; a
    cohort's EID block is 1,000,000 wide, so distinct cohorts' blocks can share
    a prefix directory). Re-linking cohort B — with no sample removed or
    renumbered on either side — must never delete cohort A's live links.
    """
    cohort_a = make_flat_bam_cohort(tmp_path, name="cohort_a", n=2)
    cohort_b = dataclasses.replace(
        make_flat_bam_cohort(tmp_path, name="cohort_b", n=2), eid_base=6_050_000
    )
    assert eid_prefix_dir(cohort_a.eid_base) == eid_prefix_dir(cohort_b.eid_base) == "60"

    samples_a, eids_a, _ = _build(cfg, cohort_a)
    samples_b, eids_b, _ = _build(cfg, cohort_b)

    a_link_paths = []
    for eid in eids_a.values():
        d = cfg.wgs_dir() / eid_prefix_dir(eid)
        a_link_paths.append(d / f"{eid}_24048_0_0.bam")
        a_link_paths.append(d / f"{eid}_24048_0_0.bam.bai")
    assert all(p.is_symlink() for p in a_link_paths)
    a_targets_before = {p: p.resolve() for p in a_link_paths}

    # Re-link cohort B only. Nothing about cohort A changed.
    link_cohort(cfg, cohort_b, samples_b, eids_b)

    for p in a_link_paths:
        assert p.is_symlink(), f"cohort A link {p} was deleted while re-linking cohort B"
        assert p.resolve() == a_targets_before[p]


def test_prune_never_removes_a_link_outside_the_cohort_eid_block(cfg, tmp_path):
    """A same-prefix symlink whose EID falls outside [eid_base, eid_base +
    BLOCK_SIZE) must be left alone even though directory-scanning finds it —
    it can only belong to a different cohort's block.
    """
    from tbooo.cohorts.registry import BLOCK_SIZE

    # eid_base 6,050,000 is not bucket-aligned: prefix "60" spans 6,000,000-
    # 6,099,999, but this cohort's own block only starts at 6,050,000, so
    # 6,000,000-6,049,999 share the prefix while sitting outside the block.
    cohort = dataclasses.replace(make_flat_bam_cohort(tmp_path, n=1), eid_base=6_050_000)
    samples = discover_samples(cohort)
    eids = assign_eids(cohort, samples)
    link_cohort(cfg, cohort, samples, eids)

    prefix_dir = cfg.wgs_dir() / eid_prefix_dir(cohort.eid_base)
    assert prefix_dir.is_dir()
    foreign_eid = 6_000_005
    assert eid_prefix_dir(foreign_eid) == eid_prefix_dir(cohort.eid_base)
    assert not (cohort.eid_base <= foreign_eid < cohort.eid_base + BLOCK_SIZE)

    foreign_target = tmp_path / "someone_elses.bam"
    foreign_target.write_bytes(b"BAM\x01")
    foreign_link = prefix_dir / f"{foreign_eid}_24048_0_0.bam"
    foreign_link.symlink_to(foreign_target)

    link_cohort(cfg, cohort, samples, eids)  # re-run must not touch the foreign link

    assert foreign_link.is_symlink()
    assert foreign_link.resolve() == foreign_target.resolve()


def test_prune_never_removes_a_non_eid_named_symlink(cfg, tmp_path):
    """A stray symlink in an EID-prefix directory whose name doesn't parse as
    '<digits>_...' must never be treated as a pruning candidate.
    """
    cohort = make_flat_bam_cohort(tmp_path, n=1)
    samples = discover_samples(cohort)
    eids = assign_eids(cohort, samples)
    link_cohort(cfg, cohort, samples, eids)

    prefix_dir = cfg.wgs_dir() / eid_prefix_dir(cohort.eid_base)
    stray_target = tmp_path / "stray.bam"
    stray_target.write_bytes(b"BAM\x01")
    stray_link = prefix_dir / "README_not_an_eid.bam"
    stray_link.symlink_to(stray_target)

    link_cohort(cfg, cohort, samples, eids)

    assert stray_link.is_symlink()
    assert stray_link.resolve() == stray_target.resolve()


def test_prune_handles_the_fastq_side_too(cfg, tmp_path):
    """The fastq-side prune branch (cfg.raw_reads_dir()) needs its own
    exercise: same orphan-on-renumber scenario as the BAM regression test,
    but with a multi-lane FASTQ sample so an orphaned EID leaves 4 links.
    """
    cohort = make_sample_dirs_cohort(tmp_path, with_vcf=False)
    samples1 = discover_samples(cohort)
    eids1 = assign_eids(cohort, samples1)
    link_cohort(cfg, cohort, samples1, eids1)

    highest_eid_run1 = max(eids1.values())  # p02: 4 fastq links (2 lanes x 2 mates)
    orphan_dir = cfg.raw_reads_dir() / eid_prefix_dir(highest_eid_run1)
    orphans_before = list(orphan_dir.glob(f"{highest_eid_run1}_raw_*"))
    assert len(orphans_before) == 4
    assert all(p.is_symlink() for p in orphans_before)

    first_id = sorted(s.sample_id for s in samples1)[0]
    shutil.rmtree(cohort.raw_root / first_id)

    samples2 = discover_samples(cohort)
    eids2 = assign_eids(cohort, samples2)
    assert highest_eid_run1 not in eids2.values()  # confirms renumbering happened
    link_cohort(cfg, cohort, samples2, eids2)

    for p in orphans_before:
        assert not p.exists() and not p.is_symlink()
