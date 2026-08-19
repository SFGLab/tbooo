from __future__ import annotations

from dataclasses import replace

import pandas as pd
import pytest

from conftest_cohorts import make_flat_bam_cohort, make_sample_dirs_cohort
from tbooo.cohorts.discover import discover_samples
from tbooo.cohorts.eids import assign_eids
from tbooo.cohorts.phenotypes import (
    CohortData,
    build_participants,
    clinical_path,
    flowcell_of,
    participant_path,
    platform_of,
    putative_pedigree,
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


def _data(cohort) -> CohortData:
    samples = discover_samples(cohort)
    return CohortData(cohort, samples, assign_eids(cohort, samples))


def _cases(tmp_path):
    """A flat-BAM cohort carrying a diagnosis code."""
    item = _data(replace(make_flat_bam_cohort(tmp_path), icd10="A00.0"))
    return item


def test_one_row_per_sample_across_all_cohorts(cfg, tmp_path):
    data = [_data(make_flat_bam_cohort(tmp_path)), _data(make_sample_dirs_cohort(tmp_path))]
    df = pd.read_parquet(build_participants(cfg, data))
    assert len(df) == sum(len(d.samples) for d in data)
    assert df["eid"].is_unique


def test_expected_columns_present(cfg, tmp_path):
    df = pd.read_parquet(build_participants(cfg, [_data(make_flat_bam_cohort(tmp_path))]))
    for col in ("eid", "p31", "p22000", "p41270", "tbooo_cohort", "tbooo_has_sv"):
        assert col in df.columns


def test_icd10_comes_from_registry(cfg, tmp_path):
    df = pd.read_parquet(build_participants(cfg, [_cases(tmp_path)]))
    assert set(df["p41270"]) == {"A00.0"}


def test_null_icd10_cohort_is_the_control_group(cfg, tmp_path):
    df = pd.read_parquet(build_participants(cfg, [_data(make_sample_dirs_cohort(tmp_path))]))
    assert df["p41270"].isna().all()


def test_sex_is_null_pending_variant_calling(cfg, tmp_path):
    df = pd.read_parquet(build_participants(cfg, [_data(make_flat_bam_cohort(tmp_path))]))
    assert df["p31"].isna().all()


def test_has_sv_flag_reflects_available_calls(cfg, tmp_path):
    data = [_data(make_flat_bam_cohort(tmp_path)), _data(make_sample_dirs_cohort(tmp_path))]
    df = pd.read_parquet(build_participants(cfg, data)).set_index("tbooo_cohort")
    assert not df.loc["cases", "tbooo_has_sv"].any()
    assert df.loc["controls", "tbooo_has_sv"].all()


def test_batch_codes_are_written_and_interpretable(cfg, tmp_path):
    build_participants(cfg, [_data(make_sample_dirs_cohort(tmp_path))])
    codes = pd.read_csv(cfg.metadata_dir() / "batch_codes.tsv", sep="\t")
    assert list(codes.columns) == ["batch", "cohort", "platform", "flowcell"]
    assert codes["batch"].min() == 1
    assert set(codes["platform"]) == {"platform_one", "platform_two"}


def test_batch_code_distinguishes_flowcells(cfg, tmp_path):
    df = pd.read_parquet(build_participants(cfg, [_data(make_sample_dirs_cohort(tmp_path))]))
    # p02 has two flowcells; a sample maps to the batch of its first run.
    assert df["p22000"].notna().all()
    assert df["p22000"].nunique() >= 2


def test_cohort_without_flowcell_info_gets_one_batch(cfg, tmp_path):
    df = pd.read_parquet(build_participants(cfg, [_data(make_flat_bam_cohort(tmp_path))]))
    assert df["p22000"].nunique() == 1


def test_flowcell_and_platform_extraction(tmp_path):
    cohort = make_sample_dirs_cohort(tmp_path)
    assert flowcell_of(cohort, "AA1234_L01_7") == "AA1234"
    assert platform_of(cohort, "AA1234") == "platform_one"
    assert platform_of(cohort, "BB5678") == "platform_two"
    assert platform_of(cohort, "ZZ0000") == "unknown"


def test_clinical_file_is_joined_with_prefix(cfg, tmp_path):
    item = _data(make_sample_dirs_cohort(tmp_path))
    path = clinical_path(cfg, item.cohort)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("sample_id\tage\tseverity\np01\t41\tmild\np02\t55\tsevere\n")
    df = pd.read_parquet(build_participants(cfg, [item])).set_index("eid")
    assert df.loc[item.eids["p01"], "clin_age"] == 41
    assert df.loc[item.eids["p02"], "clin_severity"] == "severe"


def test_missing_clinical_file_is_a_no_op(cfg, tmp_path):
    df = pd.read_parquet(build_participants(cfg, [_data(make_sample_dirs_cohort(tmp_path))]))
    assert not [c for c in df.columns if c.startswith("clin_")]


def test_clinical_columns_are_null_for_other_cohorts(cfg, tmp_path):
    controls = _data(make_sample_dirs_cohort(tmp_path))
    path = clinical_path(cfg, controls.cohort)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("sample_id\tage\np01\t41\np02\t55\n")
    cases = _data(make_flat_bam_cohort(tmp_path))
    df = pd.read_parquet(build_participants(cfg, [cases, controls]))
    assert df[df["tbooo_cohort"] == "cases"]["clin_age"].isna().all()


def test_putative_pedigree_written_only_when_regex_given(cfg, tmp_path):
    item = _data(make_sample_dirs_cohort(tmp_path))
    assert putative_pedigree(cfg, item) is None

    cohort = replace(item.cohort, pedigree_regex=r"^(?P<stem>p0)(?P<member>\d)$")
    out = putative_pedigree(cfg, CohortData(cohort, item.samples, item.eids))
    assert out is not None and "putative" in out.name
    ped = pd.read_csv(out, sep="\t")
    assert list(ped.columns) == ["eid", "sample_id", "family_stem", "member_suffix"]
    assert set(ped["family_stem"]) == {"p0"}


def test_no_relatedness_file_is_produced(cfg, tmp_path):
    """Real relatedness needs called genotypes; a guessed one would be worse."""
    build_participants(cfg, [_data(make_sample_dirs_cohort(tmp_path))])
    assert not (cfg.metadata_dir() / "ukb_rel.txt").exists()


def test_shared_clinical_column_across_two_cohorts_does_not_collide(cfg, tmp_path):
    """Two cohorts each carry a clinical TSV with a shared column name ('age') plus
    a column unique to that cohort. build_participants must concatenate the
    per-cohort clinical frames first and merge once -- merging per-cohort in a
    loop instead produces pandas' default suffixes, clin_age_x / clin_age_y,
    instead of a single clin_age column with the union of both cohorts' rows.
    """
    cases = _data(make_flat_bam_cohort(tmp_path))
    controls = _data(make_sample_dirs_cohort(tmp_path))

    cases_path = clinical_path(cfg, cases.cohort)
    cases_path.parent.mkdir(parents=True, exist_ok=True)
    cases_ids = sorted(s.sample_id for s in cases.samples)
    cases_path.write_text(
        "sample_id\tage\tonset_year\n"
        + "".join(f"{sid}\t{30 + i}\t{2000 + i}\n" for i, sid in enumerate(cases_ids))
    )

    controls_path = clinical_path(cfg, controls.cohort)
    controls_path.parent.mkdir(parents=True, exist_ok=True)
    controls_ids = sorted(s.sample_id for s in controls.samples)
    controls_path.write_text(
        "sample_id\tage\tseverity\n"
        + "".join(
            f"{sid}\t{50 + i}\t{'mild' if i % 2 == 0 else 'severe'}\n"
            for i, sid in enumerate(controls_ids)
        )
    )

    df = pd.read_parquet(build_participants(cfg, [cases, controls])).set_index("eid")

    clin_cols = {c for c in df.columns if c.startswith("clin_")}
    assert clin_cols == {"clin_age", "clin_onset_year", "clin_severity"}
    assert "clin_age_x" not in df.columns
    assert "clin_age_y" not in df.columns

    for i, sid in enumerate(cases_ids):
        eid = cases.eids[sid]
        assert df.loc[eid, "clin_age"] == 30 + i
        assert df.loc[eid, "clin_onset_year"] == 2000 + i
        assert pd.isna(df.loc[eid, "clin_severity"])

    for i, sid in enumerate(controls_ids):
        eid = controls.eids[sid]
        assert df.loc[eid, "clin_age"] == 50 + i
        assert pd.isna(df.loc[eid, "clin_onset_year"])
        assert df.loc[eid, "clin_severity"] == ("mild" if i % 2 == 0 else "severe")
