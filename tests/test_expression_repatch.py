"""Saved expression scores go into a participant table that was built after them."""

import gzip

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from tbooo.pipeline import geuvadis, gtex, mage


@pytest.fixture
def cfg(tmp_path):
    from tbooo.config import Config, Tools

    c = Config(
        data_dir=tmp_path / "instance",
        reference_dir=tmp_path / "instance" / "reference",
        tmp_dir=tmp_path / "instance" / "tmp",
        autosomes=[1], sex_chromosomes=[],
        kg_phase3_release_date="", kg_phase3_vcf_version="", kg_phase3_base_url="",
        kg_nygc_date="", kg_nygc_base_url="", sgdp_populations=[],
        kg_eid_start=1_000_000, sgdp_eid_start=2_000_000,
        array_manifest="", array_proxy_maf=0.05, exome_bed="", tools=Tools(),
    )
    c.metadata_dir().mkdir(parents=True)
    c.showcase_dir().mkdir(parents=True)
    pd.DataFrame({"eid": [1000000, 1000001], "sample_id": ["HG00096", "HG00097"]}).to_csv(
        c.metadata_dir() / "eid_map_1kg.tsv", sep="\t", index=False)
    return c


def _table(cfg) -> pd.DataFrame:
    return pq.read_table(str(cfg.showcase_dir() / "participant.parquet")).to_pandas()


def _fresh_table(cfg) -> None:
    df = pd.DataFrame({"eid": [1000000, 1000001], "p31": [1, 2]})
    pq.write_table(pa.Table.from_pandas(df), str(cfg.showcase_dir() / "participant.parquet"))


def test_geuvadis_scores_reach_a_table_built_after_them(cfg):
    cfg.geuvadis_raw_dir().mkdir(parents=True)
    (cfg.geuvadis_raw_dir() / geuvadis._RPKM_FILENAME).write_bytes(b"")
    pd.DataFrame({"eid": [1000000], "geuvadis_pc1": [0.5]}).to_csv(
        cfg.metadata_dir() / "geuvadis_expression_pcs.tsv", sep="\t", index=False)
    pd.DataFrame({"pc": ["geuvadis_pc1"], "variance_explained": [1.0], "cumulative_variance": [1.0]}).to_csv(
        cfg.metadata_dir() / "geuvadis_pca_variance.tsv", sep="\t", index=False)
    _fresh_table(cfg)
    geuvadis.build_expression_pcs(cfg)
    assert _table(cfg)["geuvadis_pc1"].tolist()[0] == 0.5
    assert pd.isna(_table(cfg)["geuvadis_pc1"].tolist()[1])


def test_mage_scores_reach_a_table_built_after_them(cfg):
    cfg.mage_raw_dir().mkdir(parents=True)
    (cfg.mage_raw_dir() / mage._EXPRESSION_FILE).write_bytes(b"")
    pd.DataFrame({"eid": [1000001], "mage_pc1": [1.5]}).to_csv(
        cfg.metadata_dir() / "mage_expression_pcs.tsv", sep="\t", index=False)
    _fresh_table(cfg)
    mage.build_expression_pcs(cfg)
    assert _table(cfg)["mage_pc1"].tolist()[1] == 1.5


def test_gtex_rows_reach_a_table_built_after_them(cfg):
    cfg.gtex_raw_dir().mkdir(parents=True)
    with gzip.open(cfg.gtex_raw_dir() / f"gene_tpm_v10_{gtex.REPRESENTATIVE_TISSUE}.gct.gz", "wt") as fh:
        fh.write("")
    pd.DataFrame({"eid": [gtex.EID_BASE], "gtex_pc1": [2.5]}).to_csv(
        cfg.metadata_dir() / "gtex_expression_pcs.tsv", sep="\t", index=False)
    pd.DataFrame({"eid": [gtex.EID_BASE], "donor": ["GTEX-1"]}).to_csv(
        cfg.metadata_dir() / "eid_map_gtex.tsv", sep="\t", index=False)
    _fresh_table(cfg)
    gtex.build_expression_pcs(cfg)
    gtex.build_expression_pcs(cfg)
    table = _table(cfg)
    assert table["eid"].tolist() == [1000000, 1000001, gtex.EID_BASE]
    assert table["gtex_pc1"].tolist()[2] == 2.5
