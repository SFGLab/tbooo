"""Load and validate the cohort registry.

The registry is the single place cohort-specific detail is allowed to live. It is
gitignored; `cohorts.example.yaml` documents the schema with placeholder values.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

# EID blocks already claimed by the public instance. A restricted cohort that
# reused one would produce EIDs meaning two different people.
RESERVED_EID_BASES: dict[int, str] = {
    1_000_000: "1KGP (pipeline/eids.py)",
    2_000_000: "SGDP (pipeline/eids.py)",
    3_000_000: "GTEx (pipeline/gtex.py:EID_BASE)",
}

# Every cohort (and every reserved public range above) owns a contiguous,
# half-open block of this many EIDs: [eid_base, eid_base + BLOCK_SIZE). Two
# cohorts whose blocks overlap would produce ambiguous EIDs, so registry
# validation and the symlink-pruning safety check both key off this constant.
BLOCK_SIZE = 1_000_000

DISCOVERY_RULES = ("flat_bam", "sample_dirs")
MODALITIES = ("bam", "fastq", "sv_vcf")


def _blocks_overlap(base_a: int, base_b: int) -> bool:
    """True if the half-open [base, base + BLOCK_SIZE) EID blocks intersect."""
    return abs(base_a - base_b) < BLOCK_SIZE


class RegistryError(Exception):
    """Raised when cohorts.yaml is malformed or internally inconsistent."""


@dataclass(frozen=True)
class SourceFile:
    modality: str          # bam | bai | fastq_r1 | fastq_r2 | sv_vcf
    path: Path
    run: str = ""          # run/lane token, FASTQ only


@dataclass(frozen=True)
class Sample:
    sample_id: str
    files: tuple[SourceFile, ...] = ()


@dataclass(frozen=True)
class Cohort:
    name: str
    eid_base: int
    raw_root: Path
    discovery: str
    modalities: tuple[str, ...]
    vcf_root: Path | None = None
    vcf_pattern: str = "{sample}.vcf"
    exclude: frozenset[str] = field(default_factory=frozenset)
    icd10: str | None = None
    platform_patterns: tuple[tuple[str, str], ...] = ()
    flowcell_regex: str | None = None
    pedigree_regex: str | None = None
    # Overrides the default FASTQ filename convention (see discover.run_token).
    # None means "use the default": <run>_1.fq.gz / <run>_2.fq.gz (or .fastq.gz).
    fastq_pattern: str | None = None

    def has(self, modality: str) -> bool:
        return modality in self.modalities


def load_registry(path: Path | str) -> dict[str, Cohort]:
    """Parse and validate cohorts.yaml, returning name -> Cohort."""
    path = Path(path)
    if not path.exists():
        raise RegistryError(
            f"cohort registry not found: {path}\n"
            f"Copy cohorts.example.yaml to {path} and fill in real values. "
            f"The real registry is gitignored by design."
        )
    raw: dict[str, Any] = yaml.safe_load(path.read_text()) or {}
    entries = raw.get("cohorts") or {}
    if not entries:
        raise RegistryError(f"{path}: no cohorts defined")

    cohorts: dict[str, Cohort] = {}
    seen_bases: dict[str, int] = {}
    for name, spec in entries.items():
        cohort = _build(name, spec)
        for reserved_base, label in RESERVED_EID_BASES.items():
            if _blocks_overlap(cohort.eid_base, reserved_base):
                raise RegistryError(
                    f"cohort {name!r}: eid_base {cohort.eid_base} (block "
                    f"[{cohort.eid_base}, {cohort.eid_base + BLOCK_SIZE})) overlaps "
                    f"the reserved block starting at {reserved_base} ({label})"
                )
        for other_name, other_base in seen_bases.items():
            if _blocks_overlap(cohort.eid_base, other_base):
                raise RegistryError(
                    f"cohort {name!r}: eid_base {cohort.eid_base} (block "
                    f"[{cohort.eid_base}, {cohort.eid_base + BLOCK_SIZE})) overlaps "
                    f"cohort {other_name!r}'s block starting at {other_base}"
                )
        seen_bases[name] = cohort.eid_base
        cohorts[name] = cohort
    return cohorts


def _build(name: str, spec: dict[str, Any]) -> Cohort:
    def req(key: str) -> Any:
        if key not in spec:
            raise RegistryError(f"cohort {name!r}: missing required key {key!r}")
        return spec[key]

    discovery = req("discovery")
    if discovery not in DISCOVERY_RULES:
        raise RegistryError(
            f"cohort {name!r}: unknown discovery rule {discovery!r} "
            f"(expected one of {', '.join(DISCOVERY_RULES)})"
        )

    modalities = tuple(req("modalities"))
    for m in modalities:
        if m not in MODALITIES:
            raise RegistryError(
                f"cohort {name!r}: unknown modality {m!r} "
                f"(expected one of {', '.join(MODALITIES)})"
            )

    raw_root = Path(req("raw_root"))
    if not raw_root.is_dir():
        raise RegistryError(f"cohort {name!r}: raw_root does not exist: {raw_root}")

    vcf_root = Path(spec["vcf_root"]) if spec.get("vcf_root") else None
    if "sv_vcf" in modalities:
        if vcf_root is None:
            raise RegistryError(
                f"cohort {name!r}: modality 'sv_vcf' declared but no vcf_root given"
            )
        if not vcf_root.is_dir():
            raise RegistryError(f"cohort {name!r}: vcf_root does not exist: {vcf_root}")

    # Named groups required by the code that later reads each regex's matches
    # (tbooo.cohorts.phenotypes: flowcell_of() reads 'flowcell'; putative_pedigree()
    # reads 'stem' and 'member'; tbooo.cohorts.discover: run_token() reads 'run' and
    # 'mate'). Checked here so a malformed registry fails loudly at load time
    # instead of silently degrading downstream — and, for fastq_pattern, instead of
    # discovering zero FASTQ files and taking that to mean the cohort is empty.
    _REQUIRED_GROUPS = {
        "flowcell_regex": ("flowcell",),
        "pedigree_regex": ("stem", "member"),
        "fastq_pattern": ("run", "mate"),
    }
    for key in ("flowcell_regex", "pedigree_regex", "fastq_pattern"):
        if spec.get(key):
            try:
                compiled = re.compile(spec[key])
            except re.error as exc:
                raise RegistryError(f"cohort {name!r}: {key} is not valid regex: {exc}") from exc
            missing = [g for g in _REQUIRED_GROUPS[key] if g not in compiled.groupindex]
            if missing:
                raise RegistryError(
                    f"cohort {name!r}: {key} is missing required named group(s) "
                    f"{', '.join(repr(g) for g in missing)}: {spec[key]!r}"
                )

    platform_patterns = tuple(
        (str(pat), str(label)) for pat, label in (spec.get("platform_patterns") or [])
    )
    for pat, _ in platform_patterns:
        try:
            re.compile(pat)
        except re.error as exc:
            raise RegistryError(
                f"cohort {name!r}: platform_patterns entry {pat!r} is not valid regex: {exc}"
            ) from exc

    return Cohort(
        name=name,
        eid_base=int(req("eid_base")),
        raw_root=raw_root,
        discovery=discovery,
        modalities=modalities,
        vcf_root=vcf_root,
        vcf_pattern=spec.get("vcf_pattern", "{sample}.vcf"),
        exclude=frozenset(spec.get("exclude") or []),
        icd10=spec.get("icd10"),
        platform_patterns=platform_patterns,
        flowcell_regex=spec.get("flowcell_regex"),
        pedigree_regex=spec.get("pedigree_regex"),
        fastq_pattern=spec.get("fastq_pattern"),
    )
