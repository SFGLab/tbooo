from __future__ import annotations

import dataclasses

import pytest

from conftest_cohorts import make_flat_bam_cohort, make_sample_dirs_cohort
from tbooo.cohorts.discover import DiscoveryError, discover_samples, run_token


def test_flat_bam_discovers_every_bam_with_its_index(tmp_path):
    cohort = make_flat_bam_cohort(tmp_path)
    samples = discover_samples(cohort)
    assert [s.sample_id for s in samples] == sorted(
        ["s001_x_20", "s-002-TAG", "s___003_20"]
    )
    first = samples[0]
    assert {f.modality for f in first.files} == {"bam", "bai"}


def test_flat_bam_tolerates_missing_index(tmp_path):
    cohort = make_flat_bam_cohort(tmp_path)
    (cohort.raw_root / "s001_x_20.bam.bai").unlink()
    sample = next(s for s in discover_samples(cohort) if s.sample_id == "s001_x_20")
    assert {f.modality for f in sample.files} == {"bam"}


def test_sample_dirs_excludes_non_sample_entries(tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    assert [s.sample_id for s in discover_samples(cohort)] == ["p01", "p02"]


def test_fastq_pairs_get_matching_run_tokens(tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    p01 = discover_samples(cohort)[0]
    fq = sorted(
        (f for f in p01.files if f.modality.startswith("fastq")),
        key=lambda f: f.modality,
    )
    assert [f.modality for f in fq] == ["fastq_r1", "fastq_r2"]
    assert fq[0].run == fq[1].run == "AA1234_L01_7"


def test_sample_id_prefix_is_stripped_from_run_token(tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    p02 = discover_samples(cohort)[1]
    runs = sorted({f.run for f in p02.files if f.modality.startswith("fastq")})
    assert runs == ["BB5678_L1", "BB9999_L2"]


def test_multi_lane_sample_keeps_every_lane(tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    p02 = discover_samples(cohort)[1]
    assert len([f for f in p02.files if f.modality.startswith("fastq")]) == 4


def test_sv_vcf_attached_via_pattern(tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    for s in discover_samples(cohort):
        assert [f for f in s.files if f.modality == "sv_vcf"]


def test_missing_sv_vcf_is_tolerated(tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    (cohort.vcf_root / "calls__p01.vcf").unlink()
    p01 = discover_samples(cohort)[0]
    assert not [f for f in p01.files if f.modality == "sv_vcf"]


def test_run_token_collision_fails_loudly(tmp_path):
    """Two files collapsing to the same run token would silently clobber a symlink."""
    cohort = make_sample_dirs_cohort(tmp_path)
    d = cohort.raw_root / "p01"
    # Same run token as the existing pair, different extension spelling.
    (d / "AA1234_L01_7_1.fastq.gz").write_bytes(b"\x1f\x8b")
    with pytest.raises(DiscoveryError, match="collision"):
        discover_samples(cohort)


def test_unpaired_fastq_fails_loudly(tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    (cohort.raw_root / "p01" / "AA1234_L01_7_2.fq.gz").unlink()
    with pytest.raises(DiscoveryError, match="unpaired"):
        discover_samples(cohort)


def test_empty_cohort_fails_loudly(tmp_path):
    cohort = make_flat_bam_cohort(tmp_path)
    for f in cohort.raw_root.iterdir():
        f.unlink()
    with pytest.raises(DiscoveryError, match="no samples"):
        discover_samples(cohort)


@pytest.mark.parametrize(
    "filename,sample_id,expected_run,expected_mate",
    [
        ("AA1_L1_1.fq.gz", "p01", "AA1_L1", "1"),
        ("AA1_L1_2.fq.gz", "p01", "AA1_L1", "2"),
        ("p02_BB_L2_1.fq.gz", "p02", "BB_L2", "1"),
        ("p02_BB_L2_1.fastq.gz", "p02", "BB_L2", "1"),
        ("notes.txt", "p01", "", ""),
    ],
)
def test_run_token_parsing(filename, sample_id, expected_run, expected_mate):
    assert run_token(filename, sample_id) == (expected_run, expected_mate)


def _rename_to_illumina(sample_dir):
    """Re-spell a sample's FASTQ in Illumina bcl2fastq style."""
    for path in sorted(sample_dir.glob("*.fq.gz")):
        mate = "1" if path.name.endswith("_1.fq.gz") else "2"
        path.rename(sample_dir / f"{sample_dir.name}_S1_L001_R{mate}_001.fastq.gz")


def test_declared_fastq_modality_that_matches_nothing_fails_loudly(tmp_path):
    """Regression: a naming convention the default token parser does not know
    used to yield a Sample with zero files and no error at all, which made the
    link stage prune the cohort's entire EID block and exit 0."""
    cohort = make_sample_dirs_cohort(tmp_path)
    _rename_to_illumina(cohort.raw_root / "p01")
    with pytest.raises(DiscoveryError, match="p01"):
        discover_samples(cohort)


def test_fastq_pattern_from_registry_handles_illumina_names(tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    for sid in ("p01", "p02"):
        d = cohort.raw_root / sid
        for path in sorted(d.glob("*.fq.gz")):
            mate = "1" if path.name.endswith("_1.fq.gz") else "2"
            path.rename(d / f"{sid}_S1_L001_R{mate}_001.fastq.gz")
    cohort = dataclasses.replace(
        cohort,
        fastq_pattern=r"^(?P<run>.+?)_R(?P<mate>[12])_001\.fastq\.gz$",
    )
    samples = discover_samples(cohort)
    p01 = samples[0]
    fq = sorted(
        (f for f in p01.files if f.modality.startswith("fastq")),
        key=lambda f: f.modality,
    )
    assert [f.modality for f in fq] == ["fastq_r1", "fastq_r2"]
    # the sample-ID prefix is stripped from the run token, as with the default
    assert fq[0].run == fq[1].run == "S1_L001"


def test_default_fastq_naming_still_works_without_a_pattern(tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    assert cohort.fastq_pattern is None
    assert [s.sample_id for s in discover_samples(cohort)] == ["p01", "p02"]


def test_run_token_uses_registry_pattern_when_given():
    pattern = r"^(?P<run>.+?)_R(?P<mate>[12])_001\.fastq\.gz$"
    assert run_token("p02_S3_L004_R2_001.fastq.gz", "p02", pattern) == ("S3_L004", "2")
    assert run_token("p02_S3_L004_I1_001.fastq.gz", "p02", pattern) == ("", "")


# ── discovery decides what a read file is; the registry regex only names runs ──

def test_a_careless_fastq_pattern_cannot_link_non_read_files(tmp_path):
    """Regression: registry validation checks only that fastq_pattern is valid
    regex with 'run'/'mate' groups. With the plausible unanchored pattern
    '^(?P<run>.+)_R(?P<mate>[12])', qc_R1.pdf/qc_R2.pdf were linked into the
    raw-reads tree as a participant's reads, manifest modality fastq_r1/r2,
    exit 0."""
    cohort = make_sample_dirs_cohort(tmp_path, name="ill", with_vcf=False)
    for sid in ("p01", "p02"):
        d = cohort.raw_root / sid
        for path in sorted(d.glob("*.fq.gz")):
            mate = "1" if path.name.endswith("_1.fq.gz") else "2"
            path.rename(d / f"{path.name[:-len('_?.fq.gz')]}_R{mate}.fq.gz")
    d = cohort.raw_root / "p01"
    (d / "qc_R1.pdf").write_bytes(b"%PDF-1.4")
    (d / "qc_R2.pdf").write_bytes(b"%PDF-1.4")

    cohort = dataclasses.replace(
        cohort, fastq_pattern=r"^(?P<run>.+)_R(?P<mate>[12])"
    )
    p01 = discover_samples(cohort)[0]
    fastq = [f for f in p01.files if f.modality.startswith("fastq")]
    assert fastq, "the real reads must still be discovered"
    for f in fastq:
        assert f.path.suffix != ".pdf", f"non-read file linked as reads: {f.path}"
        assert f.path.name.endswith((".fq.gz", ".fastq.gz", ".fq", ".fastq"))


def test_a_fastq_pattern_cannot_match_a_subdirectory(tmp_path):
    """Regression: a pattern ending '_(?P<mate>[12])$' matched DIRECTORIES, so
    run_1/ and run_2/ were accepted as a participant's paired reads."""
    from tbooo.cohorts.registry import Cohort

    root = tmp_path / "dirs_raw"
    (root / "p01").mkdir(parents=True)
    for name in ("runA_1", "runA_2"):
        (root / "p01" / name).mkdir()

    cohort = Cohort(
        name="dirsonly",
        eid_base=9_000_000,
        raw_root=root,
        discovery="sample_dirs",
        modalities=("fastq",),
        fastq_pattern=r"^(?P<run>.+)_(?P<mate>[12])$",
    )
    with pytest.raises(DiscoveryError, match="fastq"):
        discover_samples(cohort)


# ── default path: an empty run token means "not a recognised read file" ───────

@pytest.mark.parametrize(
    "filename,sample_id",
    [
        ("_1.fq.gz", "p01"),
        ("_2.fastq.gz", "p01"),
        ("s1__1.fq.gz", "s1"),
    ],
)
def test_default_convention_skips_names_with_an_empty_run_token(filename, sample_id):
    """Regression: these started returning ("", "1"), so a stray unpaired
    _1.fq.gz hard-failed discovery where it had always been skipped."""
    assert run_token(filename, sample_id) == ("", "")


def test_a_stray_unpaired_underscore_one_file_is_ignored(tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path, with_vcf=False)
    (cohort.raw_root / "p01" / "_1.fq.gz").write_bytes(b"\x1f\x8b")
    p01 = discover_samples(cohort)[0]
    runs = {f.run for f in p01.files if f.modality.startswith("fastq")}
    assert runs == {"AA1234_L01_7"}
