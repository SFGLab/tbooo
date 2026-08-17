import numpy as np
import pandas as pd
import pytest

from pcca.genotype import (
    fetch_phase3_urls,
    load_dosages,
    load_eigenvec,
    write_keep_file,
)


def test_write_keep_file_uses_an_iid_header(tmp_path):
    """The #IID header lets plink2 read a single column regardless of psam layout."""
    cohort = pd.DataFrame({"sample_id": ["HG00096", "NA12878"]})
    out = write_keep_file(cohort, tmp_path / "keep.txt")
    assert out.read_text().strip().split("\n") == ["#IID", "HG00096", "NA12878"]


def test_load_eigenvec_reorders_rows_to_match_cohort(tmp_path):
    """plink2 emits samples in .psam order; the CCA needs expression order."""
    path = tmp_path / "geno.eigenvec"
    path.write_text(
        "#FID\tIID\tPC1\tPC2\n"
        "NA12878\tNA12878\t0.5\t0.6\n"
        "HG00096\tHG00096\t0.1\t0.2\n"
    )
    cohort = pd.DataFrame({"sample_id": ["HG00096", "NA12878"]})
    pcs = load_eigenvec(path, cohort)
    assert pcs.shape == (2, 2)
    assert np.allclose(pcs[0], [0.1, 0.2])
    assert np.allclose(pcs[1], [0.5, 0.6])


def test_load_eigenvec_handles_iid_only_header(tmp_path):
    path = tmp_path / "geno.eigenvec"
    path.write_text("#IID\tPC1\tPC2\nHG00096\t0.1\t0.2\nNA12878\t0.5\t0.6\n")
    cohort = pd.DataFrame({"sample_id": ["NA12878", "HG00096"]})
    pcs = load_eigenvec(path, cohort)
    assert np.allclose(pcs[0], [0.5, 0.6])


def test_load_eigenvec_rejects_missing_samples(tmp_path):
    path = tmp_path / "geno.eigenvec"
    path.write_text("#FID\tIID\tPC1\nHG00096\tHG00096\t0.1\n")
    cohort = pd.DataFrame({"sample_id": ["HG00096", "MISSING"]})
    with pytest.raises(KeyError, match="MISSING"):
        load_eigenvec(path, cohort)


def test_load_dosages_parses_traw_and_mean_imputes(tmp_path):
    """plink2 --export A-transpose writes variants as rows, samples as columns."""
    path = tmp_path / "chr1.traw"
    path.write_text(
        "CHR\tSNP\t(C)M\tPOS\tCOUNTED\tALT\t0_HG00096\t0_HG00097\n"
        "1\t1:100:A:G\t0\t100\tG\tA\t0\t2\n"
        "1\t1:200:C:T\t0\t200\tT\tC\t1\tNA\n"
    )
    ids, positions, dosages, samples = load_dosages(path)
    assert list(ids) == ["1:100:A:G", "1:200:C:T"]
    assert list(positions) == [100, 200]
    assert samples == ["HG00096", "HG00097"]
    # samples x variants, so two rows (samples) and two columns (variants)
    assert dosages.shape == (2, 2)
    assert np.allclose(dosages[:, 0], [0.0, 2.0])
    # The NA is replaced by the column mean, which is 1.0 from the single value.
    assert np.allclose(dosages[:, 1], [1.0, 1.0])


@pytest.mark.network
def test_fetch_phase3_urls_finds_every_autosome():
    urls = fetch_phase3_urls()
    for chrom in range(1, 23):
        assert f"chr{chrom}_phase3.pgen.zst" in urls
        assert f"chr{chrom}_phase3.pvar.zst" in urls
    assert "phase3_orig.psam" in urls
    assert urls["chr1_phase3.pgen.zst"].startswith("https://")


def test_restrict_to_panel_drops_absent_samples(tmp_path):
    """Phase 3's build-37 release omits 15 of the paper's 343 individuals."""
    from pcca.genotype import restrict_to_panel

    psam = tmp_path / "phase3.psam"
    psam.write_text("#IID\tSEX\nHG00096\t1\nHG00097\t2\n")
    cohort = pd.DataFrame(
        {"sample_id": ["HG00096", "HG00097", "HG99999"], "pop": ["GBR", "GBR", "FIN"]}
    )
    out = restrict_to_panel(cohort, psam)
    assert list(out["sample_id"]) == ["HG00096", "HG00097"]
    assert list(out.index) == [0, 1], "index must be reset for row alignment"
