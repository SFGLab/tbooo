"""Tests for the SGDP sample table: the ENA report joined to the project's own table."""

from __future__ import annotations

from tbooo.download import sgdp

_ENA = (
    "analysis_accession\tsample_accession\tsample_alias\tanalysis_type\tsubmitted_ftp\n"
    "ERZ125633\tSAMEA3302684\tHGDP01168\tREFERENCE_ALIGNMENT\t"
    "ftp.sra.ebi.ac.uk/vol1/analysis/ERZ125/ERZ125633/srt.aln.bam\n"
    "ERZ324400\tSAMEA3302684\tHGDP01168\tSEQUENCE_VARIATION\t"
    "ftp.sra.ebi.ac.uk/vol1/ERZ324/ERZ324400/LP6005441-DNA_F04.annotated.nh2.vcf.gz;"
    "ftp.sra.ebi.ac.uk/vol1/ERZ324/ERZ324400/LP6005441-DNA_F04.annotated.nh2.vcf.gz.tbi\n"
    "ERZ324401\tSAMEA3302700\tNA11201\tSEQUENCE_VARIATION\t"
    "ftp://ftp.sra.ebi.ac.uk/vol1/ERZ324/ERZ324401/LP6005677-DNA_F01.annotated.nh2.vcf.gz\n"
    "ERZ324402\tSAMEA3302701\tMystery1\tSEQUENCE_VARIATION\t"
    "ftp.sra.ebi.ac.uk/vol1/ERZ324/ERZ324402/LP6005441-DNA_A09.annotated.nh2.vcf.gz\n"
)

_TABLE = (
    "#Sequencing_Panel\tIllumina_ID\tSample_ID\tSample_ID(Aliases)\tSGDP_ID\tPopulation_ID\t"
    "Region\tCountry\tTown\tContributor\tGender\tLatitude\tLongitude\tDNA_Source\tEmbargo\n"
    "C\tLP6005441-DNA_F04\tHGDP01215\tHGDP01215\tS_Daur-2\tDaur\tEastAsia\tChina\tx\ty\tF\t48.5\t124\tz\tFullyPublic\n"
    "C\tLP6005677-DNA_F01\tNA11201\tNA11201\tS_Quechua-1\tQuechua\tAmerica\tPeru\tx\ty\tM\t-13.5\t-72\tz\tFullyPublic\n"
    "C\tLP0000000-DNA_Z99\tOther\tOther\tS_Other-1\tOther\tCentralAsiaSiberia\tNowhere\tx\ty\tU\t0\t0\tz\tFullyPublic\n"
)


def test_ena_report_yields_one_row_per_vcf_keyed_by_the_illumina_id():
    rows = sgdp._parse_ena_analyses(_ENA)
    assert [r["sample_id"] for r in rows] == ["LP6005441-DNA_F04", "LP6005677-DNA_F01", "LP6005441-DNA_A09"]
    assert rows[0]["ena_accession"] == "SAMEA3302684"
    assert rows[0]["vcf_filename"] == "LP6005441-DNA_F04.annotated.nh2.vcf.gz"
    # A bare host and an ftp:// URL both become https.
    assert rows[0]["url"] == "https://ftp.sra.ebi.ac.uk/vol1/ERZ324/ERZ324400/LP6005441-DNA_F04.annotated.nh2.vcf.gz"
    assert rows[1]["url"].startswith("https://ftp.sra.ebi.ac.uk/")


def test_alignment_records_and_index_files_are_not_samples():
    rows = sgdp._parse_ena_analyses(_ENA)
    assert all(r["vcf_filename"].endswith(".vcf.gz") for r in rows)
    assert len(rows) == 3


def test_project_table_is_keyed_by_illumina_id():
    table = sgdp._parse_sgdp_table(_TABLE)
    assert set(table) == {"LP6005441-DNA_F04", "LP6005677-DNA_F01", "LP0000000-DNA_Z99"}
    assert table["LP6005441-DNA_F04"]["Population_ID"] == "Daur"


def test_samples_join_on_the_name_inside_the_vcf():
    df = sgdp._build_samples(sgdp._parse_ena_analyses(_ENA), sgdp._parse_sgdp_table(_TABLE), [])
    assert list(df.columns) == sgdp.SAMPLE_COLUMNS
    assert list(df["sample_id"]) == sorted(df["sample_id"])
    daur = df.set_index("sample_id").loc["LP6005441-DNA_F04"]
    assert (daur["population"], daur["region"], daur["sex"], daur["sgdp_id"]) == ("Daur", "East Asia", 2, "S_Daur-2")
    quechua = df.set_index("sample_id").loc["LP6005677-DNA_F01"]
    assert (quechua["region"], quechua["sex"]) == ("Native Americas", 1)


def test_a_vcf_the_table_does_not_list_keeps_its_alias():
    df = sgdp._build_samples(sgdp._parse_ena_analyses(_ENA), sgdp._parse_sgdp_table(_TABLE), [])
    row = df.set_index("sample_id").loc["LP6005441-DNA_A09"]
    assert (row["sample_alias"], row["population"], row["region"], row["sex"]) == ("Mystery1", "Unknown", "Unknown", 0)


def test_a_table_row_without_a_vcf_is_not_a_sample():
    df = sgdp._build_samples(sgdp._parse_ena_analyses(_ENA), sgdp._parse_sgdp_table(_TABLE), [])
    assert "LP0000000-DNA_Z99" not in set(df["sample_id"])


def test_population_filter():
    df = sgdp._build_samples(sgdp._parse_ena_analyses(_ENA), sgdp._parse_sgdp_table(_TABLE), ["Daur"])
    assert list(df["sample_id"]) == ["LP6005441-DNA_F04"]


def test_regions_use_the_names_the_phenotype_stage_maps():
    from tbooo.pipeline.phenotypes import _ETHNIC_BG

    for name in sgdp._REGION_NAMES.values():
        assert name in _ETHNIC_BG
