from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from tbooo.cohorts.registry import (
    BLOCK_SIZE,
    Cohort,
    RegistryError,
    RESERVED_EID_BASES,
    load_registry,
)


def _write(tmp_path: Path, data: dict) -> Path:
    p = tmp_path / "cohorts.yaml"
    p.write_text(yaml.safe_dump(data))
    return p


def _minimal(tmp_path: Path) -> dict:
    root = tmp_path / "src_alpha"
    root.mkdir(exist_ok=True)
    return {
        "cohorts": {
            "alpha": {
                "eid_base": 6000000,
                "raw_root": str(root),
                "discovery": "flat_bam",
                "modalities": ["bam"],
            }
        }
    }


def test_loads_minimal_cohort(tmp_path):
    reg = load_registry(_write(tmp_path, _minimal(tmp_path)))
    assert set(reg) == {"alpha"}
    c = reg["alpha"]
    assert isinstance(c, Cohort)
    assert c.name == "alpha"
    assert c.eid_base == 6000000
    assert c.discovery == "flat_bam"
    assert c.modalities == ("bam",)
    assert c.icd10 is None
    assert c.exclude == frozenset()


def test_rejects_unknown_discovery_rule(tmp_path):
    data = _minimal(tmp_path)
    data["cohorts"]["alpha"]["discovery"] = "telepathy"
    with pytest.raises(RegistryError, match="discovery"):
        load_registry(_write(tmp_path, data))


def test_rejects_unknown_modality(tmp_path):
    data = _minimal(tmp_path)
    data["cohorts"]["alpha"]["modalities"] = ["bam", "proteomics"]
    with pytest.raises(RegistryError, match="modality"):
        load_registry(_write(tmp_path, data))


def test_rejects_eid_base_colliding_with_public_instance(tmp_path):
    data = _minimal(tmp_path)
    data["cohorts"]["alpha"]["eid_base"] = 3000000  # GTEx
    with pytest.raises(RegistryError, match="reserved"):
        load_registry(_write(tmp_path, data))
    assert 3000000 in RESERVED_EID_BASES


def test_rejects_duplicate_eid_base_between_cohorts(tmp_path):
    data = _minimal(tmp_path)
    other = tmp_path / "src_beta"
    other.mkdir()
    data["cohorts"]["beta"] = {
        "eid_base": 6000000,
        "raw_root": str(other),
        "discovery": "flat_bam",
        "modalities": ["bam"],
    }
    with pytest.raises(RegistryError, match="eid_base"):
        load_registry(_write(tmp_path, data))


def test_rejects_sv_modality_without_vcf_root(tmp_path):
    data = _minimal(tmp_path)
    data["cohorts"]["alpha"]["modalities"] = ["bam", "sv_vcf"]
    with pytest.raises(RegistryError, match="vcf_root"):
        load_registry(_write(tmp_path, data))


def test_rejects_missing_raw_root(tmp_path):
    data = _minimal(tmp_path)
    data["cohorts"]["alpha"]["raw_root"] = str(tmp_path / "nope")
    with pytest.raises(RegistryError, match="raw_root"):
        load_registry(_write(tmp_path, data))


def test_exclude_and_optional_fields_round_trip(tmp_path):
    data = _minimal(tmp_path)
    vcf = tmp_path / "vcf_alpha"
    vcf.mkdir()
    data["cohorts"]["alpha"].update(
        {
            "discovery": "sample_dirs",
            "modalities": ["fastq", "sv_vcf"],
            "vcf_root": str(vcf),
            "vcf_pattern": "calls__{sample}.vcf",
            "exclude": ["notes.txt", "index.csv"],
            "icd10": "Z00.0",
            "platform_patterns": [["^AA", "platform_one"], ["^BB", "platform_two"]],
            "flowcell_regex": r"^(?P<flowcell>[A-Z0-9]+)_L\d+$",
            "pedigree_regex": r"^(?P<stem>.+?)(?P<member>[a-z])$",
        }
    )
    c = load_registry(_write(tmp_path, data))["alpha"]
    assert c.exclude == frozenset({"notes.txt", "index.csv"})
    assert c.icd10 == "Z00.0"
    assert c.vcf_pattern == "calls__{sample}.vcf"
    assert c.platform_patterns == (("^AA", "platform_one"), ("^BB", "platform_two"))
    assert c.pedigree_regex is not None


def test_example_registry_ships_and_parses():
    """The committed example must load — it is the template users copy."""
    example = Path(__file__).resolve().parent.parent / "cohorts.example.yaml"
    assert example.exists(), "cohorts.example.yaml must be committed"
    raw = yaml.safe_load(example.read_text())
    assert "cohorts" in raw and raw["cohorts"], "example must define at least one cohort"


def test_rejects_overlapping_eid_blocks_between_cohorts(tmp_path):
    """Two cohorts don't need identical eid_base to collide -- their
    1,000,000-wide blocks just need to overlap. This is what let two cohorts
    end up sharing an EID-prefix directory and triggered link.py deleting a
    live cohort's symlinks while pruning another cohort's orphans.
    """
    data = _minimal(tmp_path)
    other = tmp_path / "src_beta"
    other.mkdir()
    data["cohorts"]["beta"] = {
        "eid_base": 6050000,  # overlaps alpha's [6000000, 7000000) block
        "raw_root": str(other),
        "discovery": "flat_bam",
        "modalities": ["bam"],
    }
    with pytest.raises(RegistryError, match="eid_base"):
        load_registry(_write(tmp_path, data))


def test_non_overlapping_eid_blocks_are_accepted(tmp_path):
    """Sanity check for the overlap logic: blocks that are actually disjoint
    (a full BLOCK_SIZE apart) must not be rejected.
    """
    data = _minimal(tmp_path)
    other = tmp_path / "src_beta"
    other.mkdir()
    data["cohorts"]["beta"] = {
        "eid_base": 6000000 + BLOCK_SIZE,  # exactly adjacent, not overlapping
        "raw_root": str(other),
        "discovery": "flat_bam",
        "modalities": ["bam"],
    }
    reg = load_registry(_write(tmp_path, data))
    assert set(reg) == {"alpha", "beta"}


def test_rejects_eid_base_overlapping_a_reserved_public_block(tmp_path):
    """An eid_base doesn't have to equal a reserved base to collide with it --
    e.g. 3,500,000 sits inside GTEx's reserved [3,000,000, 4,000,000) block.
    """
    data = _minimal(tmp_path)
    data["cohorts"]["alpha"]["eid_base"] = 3500000
    with pytest.raises(RegistryError, match="reserved"):
        load_registry(_write(tmp_path, data))


def test_rejects_flowcell_regex_without_flowcell_group(tmp_path):
    data = _minimal(tmp_path)
    data["cohorts"]["alpha"]["flowcell_regex"] = r"^(?P<lane>[A-Z0-9]+)_L\d+$"
    with pytest.raises(RegistryError, match="flowcell_regex"):
        load_registry(_write(tmp_path, data))


def test_rejects_pedigree_regex_missing_stem_group(tmp_path):
    data = _minimal(tmp_path)
    data["cohorts"]["alpha"]["pedigree_regex"] = r"^(?P<member>[a-z])$"
    with pytest.raises(RegistryError, match="pedigree_regex"):
        load_registry(_write(tmp_path, data))


def test_rejects_pedigree_regex_missing_member_group(tmp_path):
    data = _minimal(tmp_path)
    data["cohorts"]["alpha"]["pedigree_regex"] = r"^(?P<stem>.+?)$"
    with pytest.raises(RegistryError, match="pedigree_regex"):
        load_registry(_write(tmp_path, data))


def test_rejects_pedigree_regex_missing_both_named_groups(tmp_path):
    data = _minimal(tmp_path)
    data["cohorts"]["alpha"]["pedigree_regex"] = r"^(.+?)([a-z])$"
    with pytest.raises(RegistryError, match="pedigree_regex"):
        load_registry(_write(tmp_path, data))


def test_accepts_flowcell_and_pedigree_regex_with_required_groups(tmp_path):
    data = _minimal(tmp_path)
    data["cohorts"]["alpha"]["flowcell_regex"] = r"^(?P<flowcell>[A-Z0-9]+)_L\d+$"
    data["cohorts"]["alpha"]["pedigree_regex"] = r"^(?P<stem>.+?)(?P<member>[a-z])$"
    reg = load_registry(_write(tmp_path, data))
    assert reg["alpha"].flowcell_regex == r"^(?P<flowcell>[A-Z0-9]+)_L\d+$"


# ── fastq_pattern ────────────────────────────────────────────────────────────

def test_fastq_pattern_defaults_to_none(tmp_path):
    reg = load_registry(_write(tmp_path, _minimal(tmp_path)))
    assert reg["alpha"].fastq_pattern is None


def test_rejects_fastq_pattern_without_run_group(tmp_path):
    data = _minimal(tmp_path)
    data["cohorts"]["alpha"]["fastq_pattern"] = r"^(?P<mate>[12])\.fq\.gz$"
    with pytest.raises(RegistryError, match="fastq_pattern"):
        load_registry(_write(tmp_path, data))


def test_rejects_fastq_pattern_without_mate_group(tmp_path):
    data = _minimal(tmp_path)
    data["cohorts"]["alpha"]["fastq_pattern"] = r"^(?P<run>.+)\.fq\.gz$"
    with pytest.raises(RegistryError, match="fastq_pattern"):
        load_registry(_write(tmp_path, data))


def test_rejects_fastq_pattern_that_is_not_valid_regex(tmp_path):
    data = _minimal(tmp_path)
    data["cohorts"]["alpha"]["fastq_pattern"] = r"^(?P<run>.+_R(?P<mate>[12])$"
    with pytest.raises(RegistryError, match="fastq_pattern"):
        load_registry(_write(tmp_path, data))


def test_accepts_fastq_pattern_with_required_groups(tmp_path):
    data = _minimal(tmp_path)
    pattern = r"^(?P<run>.+?)_R(?P<mate>[12])_001\.fastq\.gz$"
    data["cohorts"]["alpha"]["fastq_pattern"] = pattern
    assert load_registry(_write(tmp_path, data))["alpha"].fastq_pattern == pattern


# ── the committed example must stay free of real deployment detail ───────────

PLACEHOLDER_PREFIX = "/path/to"


def _strings(node):
    """Every string in a parsed YAML document, keys included."""
    if isinstance(node, dict):
        for key, value in node.items():
            yield from _strings(key)
            yield from _strings(value)
    elif isinstance(node, list):
        for value in node:
            yield from _strings(value)
    elif isinstance(node, str):
        yield node


def test_example_registry_uses_only_placeholder_paths():
    """Guard: the committed example may carry no real filesystem path.

    Asserted positively — every absolute path, in values and in comments alike,
    must start with the placeholder prefix. Stating it this way keeps the guard's
    protective value without the guard itself naming anything real.
    """
    example = Path(__file__).resolve().parent.parent / "cohorts.example.yaml"
    text = example.read_text()

    placeholders = 0
    for value in _strings(yaml.safe_load(text)):
        if value.startswith("/"):
            assert value.startswith(PLACEHOLDER_PREFIX + "/"), (
                f"example registry contains a non-placeholder absolute path in a "
                f"value; every path must start with {PLACEHOLDER_PREFIX}/"
            )
            placeholders += 1
    assert placeholders, "example registry declares no paths at all — guard is vacuous"

    # Comments are not parsed, so scan the raw text for path-shaped tokens too.
    for token in text.split():
        token = token.strip("\"',;()[]")
        if token.startswith("/") and "/" in token[1:]:
            assert token.startswith(PLACEHOLDER_PREFIX + "/"), (
                f"example registry contains a non-placeholder absolute path; every "
                f"path must start with {PLACEHOLDER_PREFIX}/"
            )
