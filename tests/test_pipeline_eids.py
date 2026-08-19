"""Tests for the panel parsers in src/tbooo/pipeline/eids.py.

These are pure pandas transforms over small hand-written TSV/PED-style files,
covering the branching logic (PED-style detection, case-insensitive column
matching, pop/super_pop fallback) with real files on tmp_path.
"""

from __future__ import annotations

from pathlib import Path

from tbooo.pipeline.eids import _read_nygc_panel, _read_phase3_panel


def _write(tmp_path: Path, name: str, text: str) -> Path:
    p = tmp_path / name
    p.write_text(text)
    return p


# ── _read_phase3_panel ────────────────────────────────────────────────────────

def test_read_phase3_panel_maps_gender_and_returns_expected_columns(tmp_path):
    path = _write(
        tmp_path,
        "phase3.panel",
        "sample\tpop\tsuper_pop\tgender\n"
        "HG00096\tGBR\tEUR\tmale\n"
        "HG00097\tGBR\tEUR\tfemale\n"
        "HG00099\tGBR\tEUR\tunknown\n",
    )
    df = _read_phase3_panel(path)
    assert list(df.columns) == ["sample_id", "pop", "super_pop", "sex"]
    by_id = df.set_index("sample_id")
    assert by_id.loc["HG00096", "sex"] == 1
    assert by_id.loc["HG00097", "sex"] == 2
    assert by_id.loc["HG00099", "sex"] == 0


def test_read_phase3_panel_unrecognised_gender_falls_back_to_zero(tmp_path):
    path = _write(
        tmp_path,
        "phase3.panel",
        "sample\tpop\tsuper_pop\tgender\nHG00100\tGBR\tEUR\tweird\n",
    )
    df = _read_phase3_panel(path)
    assert df.set_index("sample_id").loc["HG00100", "sex"] == 0


# ── _read_nygc_panel: PED-style detection ────────────────────────────────────

def test_read_nygc_panel_detects_ped_style_and_maps_sex_numerically(tmp_path):
    path = _write(
        tmp_path,
        "nygc.ped",
        "FamilyID\tSampleID\tFatherID\tMotherID\tSex\tPhenotype\tPopulation\tSuperpopulation\n"
        "FAM1\tNA12878\t0\t0\t2\t-9\tCEU\tEUR\n"
        "FAM2\tNA12891\t0\t0\t1\t-9\tCEU\tEUR\n",
    )
    df = _read_nygc_panel(path)
    assert list(df.columns) == ["sample_id", "pop", "super_pop", "sex"]
    by_id = df.set_index("sample_id")
    assert by_id.loc["NA12878", "sex"] == 2
    assert by_id.loc["NA12891", "sex"] == 1
    assert by_id.loc["NA12878", "pop"] == "CEU"
    assert by_id.loc["NA12878", "super_pop"] == "EUR"


def test_read_nygc_panel_ped_style_is_case_insensitive_on_column_names(tmp_path):
    # Same PED-style shape but with differently-cased/spelled headers.
    path = _write(
        tmp_path,
        "nygc_mixed_case.ped",
        "familyid\tindividualid\tfatherid\tmotherid\tsex\tphenotype\tpop\tsuperpop\n"
        "FAM1\tNA12878\t0\t0\t2\t-9\tCEU\tEUR\n",
    )
    df = _read_nygc_panel(path)
    by_id = df.set_index("sample_id")
    assert by_id.loc["NA12878", "sex"] == 2
    assert by_id.loc["NA12878", "pop"] == "CEU"
    assert by_id.loc["NA12878", "super_pop"] == "EUR"


def test_read_nygc_panel_ped_style_missing_pop_columns_fill_unk(tmp_path):
    path = _write(
        tmp_path,
        "nygc_nopop.ped",
        "FamilyID\tSampleID\tFatherID\tMotherID\tSex\tPhenotype\n"
        "FAM1\tNA12878\t0\t0\t2\t-9\n",
    )
    df = _read_nygc_panel(path)
    by_id = df.set_index("sample_id")
    assert by_id.loc["NA12878", "pop"] == "UNK"
    assert by_id.loc["NA12878", "super_pop"] == "UNK"


# ── _read_nygc_panel: generic labelled TSV fallback ──────────────────────────

def test_read_nygc_panel_generic_tsv_fallback_maps_labelled_columns(tmp_path):
    path = _write(
        tmp_path,
        "generic.tsv",
        "SampleID\tPopulation\tSuperpopulation\tGender\n"
        "HG00096\tGBR\tEUR\tmale\n"
        "HG00097\tGBR\tEUR\tfemale\n",
    )
    df = _read_nygc_panel(path)
    assert list(df.columns) == ["sample_id", "pop", "super_pop", "sex"]
    by_id = df.set_index("sample_id")
    assert by_id.loc["HG00096", "sex"] == 1
    assert by_id.loc["HG00097", "sex"] == 2
    assert by_id.loc["HG00096", "pop"] == "GBR"
    assert by_id.loc["HG00096", "super_pop"] == "EUR"


def test_read_nygc_panel_generic_tsv_missing_sex_column_fills_zero(tmp_path):
    path = _write(
        tmp_path,
        "generic_nosex.tsv",
        "SampleID\tPopulation\tSuperpopulation\nHG00096\tGBR\tEUR\n",
    )
    df = _read_nygc_panel(path)
    assert df.set_index("sample_id").loc["HG00096", "sex"] == 0


def test_read_nygc_panel_generic_tsv_missing_pop_columns_fill_unk(tmp_path):
    path = _write(
        tmp_path,
        "generic_nopop.tsv",
        "SampleID\tGender\nHG00096\tfemale\n",
    )
    df = _read_nygc_panel(path)
    by_id = df.set_index("sample_id")
    assert by_id.loc["HG00096", "pop"] == "UNK"
    assert by_id.loc["HG00096", "super_pop"] == "UNK"
    assert by_id.loc["HG00096", "sex"] == 2
