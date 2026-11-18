"""Availability flags follow the files: Phase 3 has 2,504 of the 30x panel's 3,202 samples."""

import pandas as pd

from tbooo.pipeline import runner
from tbooo.pipeline.phenotypes import _build_1kg_rows, _build_sgdp_rows
from tbooo.pipeline.qc import _write_sqc
from tbooo.pipeline.samples import phase3_sample_ids


def _kg_map() -> pd.DataFrame:
    return pd.DataFrame({
        "eid": [1000000, 1000001, 1000002],
        "sample_id": ["HG00096", "HG00097", "HG00405"],
        "pop": ["GBR", "GBR", "CHS"], "super_pop": ["EUR", "EUR", "EAS"],
        "sex": [1, 2, 2], "source": ["1kg"] * 3, "batch": [1, 1, 3],
    })


def test_a_sample_outside_phase3_has_no_array_or_imputed_data():
    rows = _build_1kg_rows(_kg_map(), {"HG00096", "HG00097"})
    assert rows["p22418"].tolist() == [1, 1, 0]
    assert rows["p22828"].tolist() == [1, 1, 0]
    assert rows["p22020"].tolist() == [1, 1, 0]
    assert rows["p24051"].tolist() == [1, 1, 1]


def test_without_the_panel_every_sample_is_marked():
    rows = _build_1kg_rows(_kg_map(), None)
    assert rows["p22418"].tolist() == [1, 1, 1]


def test_phase3_sample_ids_reads_the_panel(tmp_path):
    class Cfg:
        kg_phase3_release_date = "20130502"
        def kg_raw_dir(self):
            return tmp_path
    assert phase3_sample_ids(Cfg()) is None
    # The release's file ends each line with two empty fields.
    (tmp_path / "integrated_call_samples_v3.20130502.ALL.panel").write_text(
        "sample\tpop\tsuper_pop\tgender\t\t\nHG00096\tGBR\tEUR\tmale\nHG00097\tGBR\tEUR\tfemale\n")
    assert phase3_sample_ids(Cfg()) == {"HG00096", "HG00097"}


def test_an_sgdp_sample_without_a_vcf_has_no_wgs_flag():
    df = pd.DataFrame({"eid": [2000000, 2000001], "sex": [1, None], "region": ["Africa", "Oceania"],
                       "vcf_filename": ["LP1.annotated.nh2.vcf.gz", None]})
    assert _build_sgdp_rows(df)["p24051"].tolist() == [1, 0]


def test_pca_flag_follows_the_array_files(tmp_path):
    class Cfg:
        def metadata_dir(self):
            return tmp_path
        def kg_raw_dir(self):
            return tmp_path
    eid_map = _kg_map()
    eid_map.loc[len(eid_map)] = [2000000, "LP1", "", "", 1, "sgdp", 0]
    _write_sqc(Cfg(), eid_map, pd.DataFrame(columns=["IID", "het_outlier"]), {"1000000", "1000001", "2000000"})
    sqc = pd.read_csv(tmp_path / "ukb_sqc_v2.txt", sep="\t")
    assert sqc["used.in.pca.calculation"].tolist() == [1, 1, 0, 0]


def test_every_target_of_the_command_is_one_the_runner_runs():
    from tbooo.cli import run_pipeline
    choices = next(p for p in run_pipeline.params if p.name == "target").type.choices
    assert set(choices) == set(runner.TARGETS)
