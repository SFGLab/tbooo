"""Enumerate a cohort's samples and their source files.

Discovery is the only stage that differs between cohorts, and even here the
variation is captured by two rules plus registry-supplied patterns. Nothing in
this module knows any real cohort.
"""

from __future__ import annotations

import re
from pathlib import Path

from tbooo.cohorts.registry import Cohort, RegistryError, Sample, SourceFile

_FASTQ_SUFFIXES = (".fq.gz", ".fastq.gz")

# What discovery will accept as a sequencing read file, whatever a registry
# `fastq_pattern` happens to match. This is deliberately a property of discovery
# rather than of the user's regex: registry validation can only check that the
# regex compiles and has 'run'/'mate' groups, so a careless but plausible
# pattern (say the unanchored '^(?P<run>.+)_R(?P<mate>[12])') otherwise linked
# qc_R1.pdf/qc_R2.pdf into the raw-reads tree as a participant's reads, and a
# pattern ending '_(?P<mate>[12])$' even matched SUBDIRECTORIES.
_FASTQ_EXTENSIONS = (".fq.gz", ".fastq.gz", ".fq", ".fastq")

# Which SourceFile modalities satisfy each registry-declared modality.
_MODALITY_FILE_KINDS: dict[str, tuple[str, ...]] = {
    "bam": ("bam",),
    "fastq": ("fastq_r1", "fastq_r2"),
    "sv_vcf": ("sv_vcf",),
}

# A sample may legitimately have no SV calls (the caller failed, the sample was
# dropped from the callset), so sv_vcf is advisory. Sequence modalities are not:
# a declared one that matched nothing means the registry and the tree disagree.
_OPTIONAL_MODALITIES = frozenset({"sv_vcf"})

# Mate tokens a `fastq_pattern` may capture, normalised to "1"/"2".
_MATE_ALIASES = {"1": "1", "2": "2", "R1": "1", "R2": "2", "r1": "1", "r2": "2"}


class DiscoveryError(RegistryError):
    """Raised when a source tree does not match what the registry declares."""


def run_token(
    filename: str, sample_id: str, fastq_pattern: str | None = None
) -> tuple[str, str]:
    """Split a FASTQ filename into (run token, mate).

    The run token is the filename with the mate suffix and any leading sample-ID
    prefix removed, which makes it unique within a sample across lanes and
    flowcells without inventing an ordering. Returns ("", "") for non-FASTQ names.

    With no `fastq_pattern` the default convention applies: `<run>_1.fq.gz` /
    `<run>_2.fq.gz` (`.fastq.gz` equally). A registry that sets `fastq_pattern`
    supplies a regex with named groups `run` and `mate`, matched against the whole
    filename, which is how any other convention — Illumina bcl2fastq's
    `<sample>_S1_L001_R1_001.fastq.gz`, say — is expressed. The sample-ID prefix
    is stripped from the captured run token either way.
    """
    if fastq_pattern:
        return _run_token_from_pattern(filename, sample_id, fastq_pattern)
    for suffix in _FASTQ_SUFFIXES:
        for mate in ("1", "2"):
            tail = f"_{mate}{suffix}"
            if filename.endswith(tail):
                run = _strip_sample_prefix(filename[: -len(tail)], sample_id)
                # An empty run token means the name is nothing but the mate
                # suffix ("_1.fq.gz"), or the sample ID plus it ("s1__1.fq.gz").
                # Under the default convention that is not a recognised read
                # file, so it is ignored rather than hard-failing discovery as
                # an unpaired read.
                return (run, mate) if run else ("", "")
    return "", ""


def _strip_sample_prefix(stem: str, sample_id: str) -> str:
    prefix = f"{sample_id}_"
    return stem[len(prefix) :] if stem.startswith(prefix) else stem


def _run_token_from_pattern(
    filename: str, sample_id: str, fastq_pattern: str
) -> tuple[str, str]:
    match = re.match(fastq_pattern, filename)
    if not match:
        return "", ""
    raw_mate = match.group("mate")
    mate = _MATE_ALIASES.get(raw_mate)
    if mate is None:
        raise DiscoveryError(
            f"fastq_pattern captured mate {raw_mate!r} from {filename!r}; the 'mate' "
            f"group must capture one of {', '.join(sorted(_MATE_ALIASES))}"
        )
    return _strip_sample_prefix(match.group("run"), sample_id), mate


def _check_declared_modalities(
    cohort: Cohort, sample_id: str, location: Path, files: list[SourceFile]
) -> None:
    """Fail loudly when a sample yields none of a modality the registry declares.

    Without this a naming convention the discovery rule does not recognise
    produced a Sample with zero files and no error at all: the link stage then
    wrote an empty manifest, pruned every link in the cohort's EID block as
    orphaned, and exited 0.
    """
    present = {f.modality for f in files}
    for modality in cohort.modalities:
        if modality in _OPTIONAL_MODALITIES:
            continue
        if present.intersection(_MODALITY_FILE_KINDS[modality]):
            continue
        hint = ""
        if modality == "fastq" and cohort.fastq_pattern:
            hint = (
                "\nThis cohort's registry fastq_pattern matched no filename here: "
                f"{cohort.fastq_pattern!r}"
            )
        elif modality == "fastq":
            hint = (
                "\nExpected the default convention '<run>_1.fq.gz'/'<run>_2.fq.gz' "
                "(or .fastq.gz). If this cohort names its FASTQ differently, set "
                "'fastq_pattern' for it in the registry — a regex over the whole "
                "filename with named groups 'run' and 'mate'."
            )
        raise DiscoveryError(
            f"cohort {cohort.name!r} sample {sample_id!r}: the registry declares "
            f"modality {modality!r} but no {modality} file was found in "
            f"{location}.{hint}"
        )


def discover_samples(cohort: Cohort) -> list[Sample]:
    """Enumerate samples for `cohort`, sorted by sample_id."""
    if cohort.discovery == "flat_bam":
        samples = _discover_flat_bam(cohort)
    elif cohort.discovery == "sample_dirs":
        samples = _discover_sample_dirs(cohort)
    else:  # pragma: no cover — registry validation rejects this earlier
        raise DiscoveryError(f"cohort {cohort.name!r}: unknown discovery rule")

    if not samples:
        raise DiscoveryError(
            f"cohort {cohort.name!r}: no samples found under {cohort.raw_root}"
        )
    return sorted(samples, key=lambda s: s.sample_id)


def _discover_flat_bam(cohort: Cohort) -> list[Sample]:
    samples: list[Sample] = []
    for bam in cohort.raw_root.glob("*.bam"):
        if bam.name in cohort.exclude:
            continue
        sample_id = bam.name[: -len(".bam")]
        files = [SourceFile("bam", bam)]
        index = Path(str(bam) + ".bai")
        if index.exists():
            files.append(SourceFile("bai", index))
        _check_declared_modalities(cohort, sample_id, cohort.raw_root, files)
        samples.append(Sample(sample_id, tuple(files)))
    return samples


def _discover_sample_dirs(cohort: Cohort) -> list[Sample]:
    samples: list[Sample] = []
    for entry in cohort.raw_root.iterdir():
        if entry.name in cohort.exclude or not entry.is_dir():
            continue
        files: list[SourceFile] = []
        if cohort.has("fastq"):
            files.extend(_fastq_files(cohort, entry))
        if cohort.has("sv_vcf"):
            files.extend(_sv_file(cohort, entry.name))
        _check_declared_modalities(cohort, entry.name, entry, files)
        samples.append(Sample(entry.name, tuple(files)))
    return samples


def _is_fastq(path: Path) -> bool:
    """A real file whose extension says it holds reads (see _FASTQ_EXTENSIONS)."""
    return path.is_file() and path.name.endswith(_FASTQ_EXTENSIONS)


def _fastq_files(cohort: Cohort, sample_dir: Path) -> list[SourceFile]:
    sample_id = sample_dir.name
    by_run: dict[str, dict[str, Path]] = {}
    for path in sample_dir.iterdir():
        if not _is_fastq(path):
            continue
        run, mate = run_token(path.name, sample_id, cohort.fastq_pattern)
        if not mate:
            continue
        mates = by_run.setdefault(run, {})
        if mate in mates:
            raise DiscoveryError(
                f"cohort {cohort.name!r} sample {sample_id!r}: run token collision "
                f"on {run!r} — {mates[mate].name} and {path.name} would produce the "
                f"same link name"
            )
        mates[mate] = path

    files: list[SourceFile] = []
    for run, mates in sorted(by_run.items()):
        if set(mates) != {"1", "2"}:
            raise DiscoveryError(
                f"cohort {cohort.name!r} sample {sample_id!r}: unpaired FASTQ for "
                f"run {run!r} (found mates: {sorted(mates) or 'none'})"
            )
        files.append(SourceFile("fastq_r1", mates["1"], run))
        files.append(SourceFile("fastq_r2", mates["2"], run))
    return files


def _sv_file(cohort: Cohort, sample_id: str) -> list[SourceFile]:
    assert cohort.vcf_root is not None  # guaranteed by registry validation
    path = cohort.vcf_root / cohort.vcf_pattern.format(sample=sample_id)
    return [SourceFile("sv_vcf", path)] if path.exists() else []
