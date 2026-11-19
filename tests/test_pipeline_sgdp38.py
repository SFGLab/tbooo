"""The parts of the SGDP GRCh38 stage that run without bcftools."""

import gzip
import subprocess

import pytest

from tbooo.pipeline import sgdp38


@pytest.fixture
def cfg(tmp_path):
    from tbooo.config import Config, Tools

    return Config(
        data_dir=tmp_path / "instance",
        reference_dir=tmp_path / "instance" / "reference",
        tmp_dir=tmp_path / "instance" / "tmp",
        autosomes=[1, 22], sex_chromosomes=["X"],
        kg_phase3_release_date="", kg_phase3_vcf_version="", kg_phase3_base_url="",
        kg_nygc_date="", kg_nygc_base_url="", sgdp_populations=[],
        kg_eid_start=1_000_000, sgdp_eid_start=2_000_000,
        array_manifest="", array_proxy_maf=0.05, exome_bed="", tools=Tools(),
    )


def test_outputs_are_named_by_chromosome(cfg):
    assert cfg.sgdp38_phased("22").name == "sgdp_c22.phased.vcf.gz"
    assert cfg.sgdp38_unphased("22").name == "sgdp_c22.unphased.vcf.gz"
    assert cfg.sgdp38_phased("22").parent == cfg.sgdp_raw_dir() / "grch38"


def test_added_markers_give_every_panel_sample_the_reference_allele(tmp_path):
    header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG1\tHG2\tHG3\n"
    out = tmp_path / "added.vcf.gz"
    sites = [["chr22", "100", ".", "A", "G"], ["chr22", "250", "rs1", "AT", "A"]]
    sgdp38.write_added_markers(header, sites, 3, out, bgzip="gzip")
    lines = gzip.open(out, "rt").read().splitlines()
    assert lines[:2] == header.splitlines()
    assert lines[2] == "chr22\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|0\t0|0"
    assert lines[3].split("\t")[:5] == ["chr22", "250", ".", "AT", "A"]
    assert len(lines) == 4


def test_a_pipeline_passes_data_through_and_returns_the_last_output():
    out = sgdp38._pipe(["printf", "b\\na\\nb\\n"], ["sort"], ["uniq"], capture=True)
    assert out == b"a\nb\n"


def test_a_pipeline_raises_where_a_middle_command_fails():
    with pytest.raises(subprocess.CalledProcessError) as err:
        sgdp38._pipe(["printf", "x\\n"], ["false"], ["cat"], capture=True)
    assert err.value.cmd == ["false"]


def test_the_x_chromosome_is_skipped_without_running_anything(cfg, monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError("nothing may run for a chromosome the stage does not build")
    monkeypatch.setattr(sgdp38, "_build_chrom", boom)
    sgdp38.build_sgdp_grch38(cfg, ["X"])


def test_a_chromosome_without_any_source_is_refused(cfg):
    with pytest.raises(FileNotFoundError, match="tbooo download sgdp"):
        sgdp38.build_sgdp_grch38(cfg, ["22"])
