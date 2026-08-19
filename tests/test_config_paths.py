"""Tests for src/tbooo/config.py: YAML loading and derived path helpers."""

from __future__ import annotations

from pathlib import Path

import yaml

from tbooo.config import Config, Tools

_MINIMAL_YAML: dict = {
    "data_dir": "/data/tbooo",
    "reference_dir": "/data/tbooo/reference",
    "tmp_dir": "/data/tbooo/tmp",
    "autosomes": list(range(1, 23)),
    "sex_chromosomes": ["X", "Y"],
    "kg_phase3_release_date": "20130502",
    "kg_phase3_vcf_version": "v5b",
    "kg_phase3_base_url": "https://example.org/phase3",
    "kg_nygc_date": "20190425",
    "kg_nygc_base_url": "https://example.org/nygc",
    "sgdp_populations": ["French", "Han"],
    "kg_eid_start": 1_000_000,
    "sgdp_eid_start": 2_000_000,
    "array_manifest": "manifest.bpm",
    "array_proxy_maf": 0.01,
    "exome_bed": "exome.bed",
}


def _write_config(tmp_path: Path, overrides: dict | None = None) -> Path:
    data = dict(_MINIMAL_YAML)
    if overrides:
        data.update(overrides)
    p = tmp_path / "config.yaml"
    p.write_text(yaml.safe_dump(data))
    return p


# ── Config.load ───────────────────────────────────────────────────────────────

def test_load_converts_dir_fields_to_path(tmp_path):
    cfg = Config.load(_write_config(tmp_path))
    assert isinstance(cfg.data_dir, Path)
    assert isinstance(cfg.reference_dir, Path)
    assert isinstance(cfg.tmp_dir, Path)
    assert cfg.data_dir == Path("/data/tbooo")
    assert cfg.reference_dir == Path("/data/tbooo/reference")
    assert cfg.tmp_dir == Path("/data/tbooo/tmp")


def test_load_builds_tools_instance(tmp_path):
    cfg = Config.load(_write_config(tmp_path, {"tools": {"bcftools": "/opt/bcftools"}}))
    assert isinstance(cfg.tools, Tools)
    assert cfg.tools.bcftools == "/opt/bcftools"
    assert cfg.tools.plink2 == "plink2"  # untouched default


def test_load_tools_defaults_when_key_omitted(tmp_path):
    cfg = Config.load(_write_config(tmp_path))
    assert isinstance(cfg.tools, Tools)
    assert cfg.tools.bcftools == "bcftools"
    assert cfg.tools.wget == "wget"


def test_load_omitted_optional_keys_take_their_defaults(tmp_path):
    cfg = Config.load(_write_config(tmp_path))
    assert cfg.download_workers == 4
    assert cfg.deep_check_workers == 8
    assert cfg.wgs_nygc_threads == 4
    assert cfg.wgs_sgdp_merge_threads == 8
    assert cfg.wgs_final_merge_threads == 4
    assert cfg.cohorts_registry == ""


def test_load_explicit_optional_keys_override_defaults(tmp_path):
    cfg = Config.load(_write_config(tmp_path, {"download_workers": 16, "cohorts_registry": "cohorts.yaml"}))
    assert cfg.download_workers == 16
    assert cfg.cohorts_registry == "cohorts.yaml"


# ── chromosomes ───────────────────────────────────────────────────────────────

def test_chromosomes_concatenates_autosomes_then_sex_as_strings(tmp_path):
    cfg = Config.load(_write_config(tmp_path))
    assert cfg.chromosomes == [str(c) for c in range(1, 23)] + ["X", "Y"]
    assert all(isinstance(c, str) for c in cfg.chromosomes)


def test_chromosomes_empty_sex_list(tmp_path):
    cfg = Config.load(_write_config(tmp_path, {"sex_chromosomes": []}))
    assert cfg.chromosomes == [str(c) for c in range(1, 23)]


# ── *_dir() helpers ───────────────────────────────────────────────────────────

def test_dir_helpers_derive_from_data_dir(tmp_path):
    cfg = Config.load(_write_config(tmp_path))
    assert cfg.showcase_dir() == cfg.data_dir / "Showcase"
    assert cfg.metadata_dir() == cfg.data_dir / "metadata"
    assert cfg.kg_raw_dir() == cfg.data_dir / "raw" / "1kg"
    assert cfg.array_dir() == cfg.data_dir / "Bulk" / "Genotype Results" / "Genotype calls"
    assert cfg.imputed_dir() == cfg.data_dir / "Bulk" / "Imputed"
    assert cfg.wgs_dir() == cfg.data_dir / "Bulk" / "Whole genome sequences"


def test_new_dir_helpers_added_by_this_plan(tmp_path):
    cfg = Config.load(_write_config(tmp_path))
    assert cfg.raw_reads_dir() == cfg.data_dir / "Bulk" / "Raw sequencing reads"
    assert cfg.sv_dir() == cfg.data_dir / "Bulk" / "Structural variants (ConsensuSV)"
    assert cfg.slurm_dir() == cfg.data_dir / "slurm"


# ── versioned stem/path helpers ───────────────────────────────────────────────

def test_array_stem_for_autosome(tmp_path):
    cfg = Config.load(_write_config(tmp_path))
    assert cfg.array_stem("7") == cfg.array_dir() / "ukb22418_c7_b0_v2"


def test_imputed_stem_for_autosome(tmp_path):
    cfg = Config.load(_write_config(tmp_path))
    assert cfg.imputed_stem("7") == cfg.imputed_dir() / "ukb22828_c7_b0_v3"


def test_wes_stem_for_autosome(tmp_path):
    cfg = Config.load(_write_config(tmp_path))
    assert cfg.wes_stem("7") == cfg.wes_dir() / "ukb23158_c7_b0_v1"


def test_wgs_pvcf_for_autosome(tmp_path):
    cfg = Config.load(_write_config(tmp_path))
    assert cfg.wgs_pvcf("7") == cfg.wgs_dir() / "ukb24310_c7_b0_v1.pvcf.gz"


def test_wgs_pvcf_for_x_and_y_use_same_pattern_as_autosomes(tmp_path):
    # wgs_pvcf has no special-cased X/Y branch, unlike phase3_vcf/nygc_vcf.
    cfg = Config.load(_write_config(tmp_path))
    assert cfg.wgs_pvcf("X") == cfg.wgs_dir() / "ukb24310_cX_b0_v1.pvcf.gz"
    assert cfg.wgs_pvcf("Y") == cfg.wgs_dir() / "ukb24310_cY_b0_v1.pvcf.gz"


def test_phase3_vcf_autosome(tmp_path):
    cfg = Config.load(_write_config(tmp_path))
    expected = (
        cfg.kg_raw_dir()
        / "ALL.chr7.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    )
    assert cfg.phase3_vcf("7") == expected


def test_phase3_vcf_chrx_uses_its_own_naming(tmp_path):
    cfg = Config.load(_write_config(tmp_path))
    expected = (
        cfg.kg_raw_dir()
        / "ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"
    )
    assert cfg.phase3_vcf("X") == expected
    # chrX branch hardcodes v1c, ignoring the configured vcf_version.
    assert "v1c" in str(cfg.phase3_vcf("X"))
    assert "v5b" not in str(cfg.phase3_vcf("X"))


def test_phase3_vcf_chry_uses_its_own_naming(tmp_path):
    cfg = Config.load(_write_config(tmp_path))
    expected = (
        cfg.kg_raw_dir()
        / "ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz"
    )
    assert cfg.phase3_vcf("Y") == expected
    # chrY naming drops "shapeit2_mvncall" entirely and hardcodes v2b.
    assert "shapeit2" not in str(cfg.phase3_vcf("Y"))
