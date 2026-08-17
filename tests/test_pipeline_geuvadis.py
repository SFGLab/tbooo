import gzip

import numpy as np

from tbooo.pipeline.geuvadis import _is_sample_id, _load_rpkm


def test_is_sample_id_accepts_1kgp_ids():
    assert _is_sample_id("HG00096")
    assert _is_sample_id("NA12878")


def test_is_sample_id_rejects_annotation_columns():
    assert not _is_sample_id("TargetID")
    assert not _is_sample_id("Chr")
    assert not _is_sample_id("Coord")


def test_load_rpkm_collapses_gd660_run_columns_to_sample_ids(tmp_path):
    """GD660 headers are run IDs like HG00096.1.M_111124_6, not bare sample IDs."""
    path = tmp_path / "gd660.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(
            "TargetID\tGene_Symbol\tChr\tCoord\t"
            "HG00096.1.M_111124_6\tNA12878.2.M_120209_7\n"
        )
        handle.write("ENSG1\tENSG1\t1\t1000\t5.0\t7.0\n")
    frame = _load_rpkm(path)
    assert sorted(frame.columns) == ["HG00096", "NA12878"]
    assert frame.shape == (1, 2)


def test_load_rpkm_averages_technical_replicates(tmp_path):
    """660 runs cover 462 individuals; the surplus are replicates to average."""
    path = tmp_path / "gd660.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(
            "TargetID\tGene_Symbol\tChr\tCoord\t"
            "HG00096.1.M_111124_6\tHG00096.3.M_120223_6\tNA12878.2.M_120209_7\n"
        )
        handle.write("ENSG1\tENSG1\t1\t1000\t4.0\t8.0\t7.0\n")
    frame = _load_rpkm(path)
    assert sorted(frame.columns) == ["HG00096", "NA12878"]
    assert np.isclose(frame.loc["ENSG1", "HG00096"], 6.0)
    assert np.isclose(frame.loc["ENSG1", "NA12878"], 7.0)
