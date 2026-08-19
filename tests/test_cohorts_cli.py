from __future__ import annotations

import yaml
from click.testing import CliRunner

from tbooo.cli import cli


def _instance(tmp_path):
    """A restricted config + registry pointing at a synthetic source tree."""
    from conftest_cohorts import make_sample_dirs_cohort

    cohort = make_sample_dirs_cohort(tmp_path)
    registry = tmp_path / "cohorts.yaml"
    registry.write_text(
        yaml.safe_dump(
            {
                "cohorts": {
                    "controls": {
                        "eid_base": 7000000,
                        "raw_root": str(cohort.raw_root),
                        "vcf_root": str(cohort.vcf_root),
                        "discovery": "sample_dirs",
                        "modalities": ["fastq", "sv_vcf"],
                        "vcf_pattern": "calls__{sample}.vcf",
                        "exclude": ["manifest.csv", "checksums.txt"],
                    }
                }
            }
        )
    )
    config = tmp_path / "config.restricted.yaml"
    config.write_text(
        yaml.safe_dump(
            {
                "data_dir": str(tmp_path / "instance"),
                "reference_dir": str(tmp_path / "instance" / "reference"),
                "tmp_dir": str(tmp_path / "instance" / "tmp"),
                "autosomes": [1], "sex_chromosomes": [],
                "kg_phase3_release_date": "", "kg_phase3_vcf_version": "",
                "kg_phase3_base_url": "", "kg_nygc_date": "", "kg_nygc_base_url": "",
                "sgdp_populations": [], "kg_eid_start": 1000000,
                "sgdp_eid_start": 2000000, "array_manifest": "",
                "array_proxy_maf": 0.05, "exome_bed": "",
                "cohorts_registry": str(registry),
            }
        )
    )
    return config


def test_restricted_group_is_registered():
    out = CliRunner().invoke(cli, ["restricted", "--help"])
    assert out.exit_code == 0
    for sub in ("discover", "eids", "link", "sv", "phenotypes"):
        assert sub in out.output


def test_discover_reports_sample_counts(tmp_path):
    config = _instance(tmp_path)
    out = CliRunner().invoke(cli, ["restricted", "discover", "-c", str(config)])
    assert out.exit_code == 0, out.output
    assert "controls" in out.output and "2" in out.output


def test_link_creates_the_tree(tmp_path):
    config = _instance(tmp_path)
    runner = CliRunner()
    assert runner.invoke(cli, ["restricted", "eids", "-c", str(config)]).exit_code == 0
    res = runner.invoke(cli, ["restricted", "link", "-c", str(config)])
    assert res.exit_code == 0, res.output
    assert (tmp_path / "instance" / "metadata" / "source_manifest_controls.tsv").exists()


def test_sv_prints_submit_command_and_runs_nothing(tmp_path):
    config = _instance(tmp_path)
    runner = CliRunner()
    runner.invoke(cli, ["restricted", "eids", "-c", str(config)])
    res = runner.invoke(cli, ["restricted", "sv", "-c", str(config)])
    assert res.exit_code == 0, res.output
    script = tmp_path / "instance" / "slurm" / "sv_controls.sbatch"
    assert script.exists()
    # Pin the actual submit instruction, not just any line ending in .sbatch --
    # the "wrote ... .sbatch" log line also matches a bare "sbatch" substring
    # check and would keep passing even if the submit-command line were deleted.
    assert f"submit with:  sbatch {script}" in res.output


def test_phenotypes_writes_parquet(tmp_path):
    config = _instance(tmp_path)
    runner = CliRunner()
    runner.invoke(cli, ["restricted", "eids", "-c", str(config)])
    res = runner.invoke(cli, ["restricted", "phenotypes", "-c", str(config)])
    assert res.exit_code == 0, res.output
    assert (tmp_path / "instance" / "Showcase" / "participant.parquet").exists()


def test_missing_registry_gives_an_actionable_error(tmp_path):
    config = _instance(tmp_path)
    (tmp_path / "cohorts.yaml").unlink()
    res = CliRunner().invoke(cli, ["restricted", "discover", "-c", str(config)])
    assert res.exit_code != 0
    assert "cohorts.example.yaml" in res.output


def test_only_flag_selects_one_cohort(tmp_path):
    config = _instance(tmp_path)
    res = CliRunner().invoke(
        cli, ["restricted", "discover", "-c", str(config), "--only", "nosuch"]
    )
    assert res.exit_code != 0
    assert "nosuch" in res.output


def _instance_two(tmp_path):
    """A registry with TWO cohorts, so a per-cohort run can be seen to clobber."""
    from conftest_cohorts import make_flat_bam_cohort, make_sample_dirs_cohort

    controls = make_sample_dirs_cohort(tmp_path)
    cases = make_flat_bam_cohort(tmp_path)
    registry = tmp_path / "cohorts.yaml"
    registry.write_text(
        yaml.safe_dump(
            {
                "cohorts": {
                    "controls": {
                        "eid_base": 7000000,
                        "raw_root": str(controls.raw_root),
                        "vcf_root": str(controls.vcf_root),
                        "discovery": "sample_dirs",
                        "modalities": ["fastq", "sv_vcf"],
                        "vcf_pattern": "calls__{sample}.vcf",
                        "exclude": ["manifest.csv", "checksums.txt"],
                    },
                    "cases": {
                        "eid_base": 6000000,
                        "raw_root": str(cases.raw_root),
                        "discovery": "flat_bam",
                        "modalities": ["bam"],
                    },
                }
            }
        )
    )
    config = tmp_path / "config.restricted.yaml"
    config.write_text(
        yaml.safe_dump(
            {
                "data_dir": str(tmp_path / "instance"),
                "reference_dir": str(tmp_path / "instance" / "reference"),
                "tmp_dir": str(tmp_path / "instance" / "tmp"),
                "autosomes": [1], "sex_chromosomes": [],
                "kg_phase3_release_date": "", "kg_phase3_vcf_version": "",
                "kg_phase3_base_url": "", "kg_nygc_date": "", "kg_nygc_base_url": "",
                "sgdp_populations": [], "kg_eid_start": 1000000,
                "sgdp_eid_start": 2000000, "array_manifest": "",
                "array_proxy_maf": 0.05, "exome_bed": "",
                "cohorts_registry": str(registry),
            }
        )
    )
    return config


def test_phenotypes_refuses_only_and_leaves_the_table_intact(tmp_path):
    """Regression: participant.parquet is ONE table across every cohort, so
    building it from a subset silently deleted the other cohorts' rows."""
    import pandas as pd

    config = _instance_two(tmp_path)
    runner = CliRunner()
    assert runner.invoke(cli, ["restricted", "eids", "-c", str(config)]).exit_code == 0
    res = runner.invoke(cli, ["restricted", "phenotypes", "-c", str(config)])
    assert res.exit_code == 0, res.output

    parquet = tmp_path / "instance" / "Showcase" / "participant.parquet"
    batches = tmp_path / "instance" / "metadata" / "batch_codes.tsv"
    before = pd.read_parquet(parquet)
    batches_before = batches.read_text()
    assert set(before["tbooo_cohort"]) == {"cases", "controls"}

    res = runner.invoke(
        cli, ["restricted", "phenotypes", "-c", str(config), "--only", "cases"]
    )
    assert res.exit_code != 0, res.output
    assert "--only" in res.output

    after = pd.read_parquet(parquet)
    assert len(after) == len(before)
    assert set(after["tbooo_cohort"]) == {"cases", "controls"}
    assert batches.read_text() == batches_before


def test_phenotypes_help_explains_it_covers_every_cohort(tmp_path):
    res = CliRunner().invoke(cli, ["restricted", "phenotypes", "--help"])
    assert res.exit_code == 0
    assert "every cohort" in res.output or "all cohorts" in res.output


def test_link_surfaces_refusals_as_a_clean_error_not_a_traceback(tmp_path):
    """`restricted sv` wraps RuntimeError in ClickException; `restricted link`
    did not, so _check_containment/_refuse_blind_prune surfaced as a raw
    traceback instead of an `Error: ...` line."""
    from conftest_cohorts import make_sample_dirs_cohort

    cohort = make_sample_dirs_cohort(tmp_path)
    registry = tmp_path / "cohorts.yaml"
    registry.write_text(
        yaml.safe_dump(
            {
                "cohorts": {
                    "controls": {
                        "eid_base": 7000000,
                        "raw_root": str(cohort.raw_root),
                        "vcf_root": str(cohort.vcf_root),
                        "discovery": "sample_dirs",
                        "modalities": ["fastq", "sv_vcf"],
                        "vcf_pattern": "calls__{sample}.vcf",
                        "exclude": ["manifest.csv", "checksums.txt"],
                    }
                }
            }
        )
    )
    config = tmp_path / "config.overlap.yaml"
    config.write_text(
        yaml.safe_dump(
            {
                # data_dir INSIDE raw_root -- _check_containment must refuse
                "data_dir": str(cohort.raw_root / "instance"),
                "reference_dir": str(tmp_path / "instance" / "reference"),
                "tmp_dir": str(tmp_path / "instance" / "tmp"),
                "autosomes": [1], "sex_chromosomes": [],
                "kg_phase3_release_date": "", "kg_phase3_vcf_version": "",
                "kg_phase3_base_url": "", "kg_nygc_date": "", "kg_nygc_base_url": "",
                "sgdp_populations": [], "kg_eid_start": 1000000,
                "sgdp_eid_start": 2000000, "array_manifest": "",
                "array_proxy_maf": 0.05, "exome_bed": "",
                "cohorts_registry": str(registry),
            }
        )
    )
    res = CliRunner().invoke(cli, ["restricted", "link", "-c", str(config)])
    assert res.exit_code != 0
    assert "refusing to link" in res.output
    assert not isinstance(res.exception, RuntimeError)
