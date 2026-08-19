from __future__ import annotations

import pytest

from conftest_cohorts import make_flat_bam_cohort, make_sample_dirs_cohort
from tbooo.cohorts.discover import discover_samples
from tbooo.cohorts.eids import (
    assign_eids,
    eid_map_path,
    rename_path,
    write_eid_map,
    write_rename_file,
)


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


def test_eids_start_at_base_and_are_contiguous(tmp_path):
    cohort = make_flat_bam_cohort(tmp_path)
    samples = discover_samples(cohort)
    eids = assign_eids(cohort, samples)
    assert sorted(eids.values()) == [
        cohort.eid_base + i for i in range(len(samples))
    ]


def test_assignment_follows_sorted_sample_id(tmp_path):
    cohort = make_flat_bam_cohort(tmp_path)
    samples = discover_samples(cohort)
    eids = assign_eids(cohort, samples)
    ordered = [sid for sid, _ in sorted(eids.items(), key=lambda kv: kv[1])]
    assert ordered == sorted(s.sample_id for s in samples)


def test_assignment_is_deterministic_across_runs(tmp_path):
    cohort = make_flat_bam_cohort(tmp_path)
    first = assign_eids(cohort, discover_samples(cohort))
    second = assign_eids(cohort, discover_samples(cohort))
    assert first == second


def test_no_eid_collision_between_cohorts(tmp_path):
    a = make_flat_bam_cohort(tmp_path, name="cases")
    b = make_sample_dirs_cohort(tmp_path, name="controls")
    ea = assign_eids(a, discover_samples(a))
    eb = assign_eids(b, discover_samples(b))
    assert not (set(ea.values()) & set(eb.values()))


def test_eid_map_has_expected_columns(cfg, tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    samples = discover_samples(cohort)
    out = write_eid_map(cfg, cohort, samples, assign_eids(cohort, samples))
    header = out.read_text().splitlines()[0].split("\t")
    assert header == ["eid", "sample_id", "cohort", "n_files", "has_sv"]


def test_rename_file_maps_sample_id_to_eid(cfg, tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    samples = discover_samples(cohort)
    eids = assign_eids(cohort, samples)
    out = write_rename_file(cfg, cohort, eids)
    lines = dict(l.split("\t") for l in out.read_text().splitlines())
    assert lines["p01"] == str(eids["p01"])


def test_paths_are_namespaced_per_cohort(cfg, tmp_path):
    cohort = make_flat_bam_cohort(tmp_path)
    assert cohort.name in eid_map_path(cfg, cohort).name
    assert cohort.name in rename_path(cfg, cohort).name
