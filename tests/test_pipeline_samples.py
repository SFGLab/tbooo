"""Tests for the BGEN sample file: it follows the genotype file, not the EID map."""

from __future__ import annotations

import pandas as pd

from tbooo.pipeline.samples import sample_file_text

_EID_MAP = pd.DataFrame({"eid": [1000000, 1000001, 1000002, 1000003], "sex": [1, 2, 0, 2]})


def test_rows_follow_the_vcf_order_not_the_map_order():
    text = sample_file_text(["1000003", "1000000"], _EID_MAP)
    assert text.splitlines() == ["ID_1 ID_2 missing sex", "0 0 0 D", "1000003 1000003 0 2", "1000000 1000000 0 1"]


def test_a_map_with_more_samples_than_the_vcf_adds_no_rows():
    text = sample_file_text(["1000001"], _EID_MAP)
    assert len(text.splitlines()) == 3


def test_a_name_the_map_does_not_have_gets_sex_zero():
    assert sample_file_text(["HG00096"], _EID_MAP).splitlines()[-1] == "HG00096 HG00096 0 0"
