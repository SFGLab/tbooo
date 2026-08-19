"""Build the UKB-shaped symlink tree and its provenance manifest.

Bulk sequence data is symlinked, never copied: sources are large and must not be
moved or modified. SV VCFs are deliberately NOT linked — they embed the original
sample ID, so the sv stage rewrites them under EIDs instead.
"""

from __future__ import annotations

import re
from datetime import datetime, timezone
from pathlib import Path

from tbooo.cohorts.registry import BLOCK_SIZE, Cohort, Sample
from tbooo.config import Config
from tbooo.utils import ensure_dirs, eid_prefix_dir, log

MANIFEST_COLUMNS = (
    "eid",
    "sample_id",
    "modality",
    "link_path",
    "source_size",
    "source_mtime",
    "source_path",
)

# UKB Field 24048 is "Whole genome CRAM files". These sources are BAM, so the
# field ID states what the data is and the extension states what it is.
_BAM_FIELD = "24048"


def manifest_path(cfg: Config, cohort: Cohort) -> Path:
    return cfg.metadata_dir() / f"source_manifest_{cohort.name}.tsv"


def link_name(modality: str, eid: int, run: str = "", suffix: str = "") -> str:
    """UKB-shaped filename for a linked source file."""
    if modality == "bam":
        return f"{eid}_{_BAM_FIELD}_0_0.bam"
    if modality == "bai":
        return f"{eid}_{_BAM_FIELD}_0_0.bam.bai"
    if modality in ("fastq_r1", "fastq_r2"):
        mate = "1" if modality.endswith("r1") else "2"
        return f"{eid}_raw_{run}_{mate}{suffix}"
    raise ValueError(f"modality {modality!r} is not linkable")


def _dest_dir(cfg: Config, modality: str, eid: int) -> Path:
    base = cfg.wgs_dir() if modality in ("bam", "bai") else cfg.raw_reads_dir()
    return base / eid_prefix_dir(eid)


# Longest first: ".fastq.gz" must win over ".fq.gz"/".gz".
_FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def _fastq_suffix(path: Path) -> str:
    """The source's own FASTQ extension, preserved on the link name."""
    for suffix in _FASTQ_SUFFIXES:
        if path.name.endswith(suffix):
            return suffix
    # Only reachable via a registry `fastq_pattern` naming some other extension;
    # mirror whatever the source uses rather than inventing one.
    return "".join(path.suffixes[-2:]) or path.suffix


def _check_containment(cfg: Config, cohort: Cohort) -> None:
    """Refuse to link if the instance root and the source root overlap.

    A misconfigured `data_dir` that sits inside (or wraps) `raw_root` would
    otherwise have `link_cohort` create `Bulk/`/`metadata/` directories inside
    the read-only source tree. Cheap insurance in the riskiest module here.
    """
    data_dir = cfg.data_dir.resolve()
    raw_root = cohort.raw_root.resolve()
    if (
        data_dir == raw_root
        or data_dir.is_relative_to(raw_root)
        or raw_root.is_relative_to(data_dir)
    ):
        raise RuntimeError(
            f"cohort {cohort.name!r}: instance data_dir ({data_dir}) and raw_root "
            f"({raw_root}) overlap; refusing to link into or under the source tree"
        )


def link_cohort(
    cfg: Config, cohort: Cohort, samples: list[Sample], eids: dict[str, int]
) -> Path:
    """Create/repair every symlink for `cohort`, prune orphaned ones, and rewrite its manifest."""
    _check_containment(cfg, cohort)
    ensure_dirs(cfg.metadata_dir())
    rows: list[tuple] = []
    created = repaired = 0
    seen: dict[Path, Path] = {}

    for sample in sorted(samples, key=lambda s: s.sample_id):
        eid = eids[sample.sample_id]
        for src in sample.files:
            if src.modality == "sv_vcf":
                continue  # rewritten by the sv stage, not linked
            suffix = _fastq_suffix(src.path) if src.modality.startswith("fastq") else ""
            dest_dir = _dest_dir(cfg, src.modality, eid)
            ensure_dirs(dest_dir)
            dest = dest_dir / link_name(src.modality, eid, src.run, suffix)

            if dest in seen:
                raise RuntimeError(
                    f"duplicate destination {dest}: both {seen[dest]} and "
                    f"{src.path} would link there — refusing to silently pick one"
                )
            seen[dest] = src.path

            state = _ensure_link(dest, src.path)
            created += state == "created"
            repaired += state == "repaired"

            stat = src.path.stat()
            rows.append(
                (
                    eid,
                    sample.sample_id,
                    src.modality,
                    str(dest),
                    stat.st_size,
                    datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
                    str(src.path),
                )
            )

    _refuse_blind_prune(cohort, samples, rows)
    pruned = _prune_orphaned_links(cfg, cohort, eids, {Path(row[3]) for row in rows})

    out = manifest_path(cfg, cohort)
    with out.open("w") as f:
        f.write("\t".join(MANIFEST_COLUMNS) + "\n")
        for row in rows:
            f.write("\t".join(str(v) for v in row) + "\n")

    log(
        f"  {cohort.name}: {len(rows)} links "
        f"({created} created, {repaired} repaired, {pruned} pruned) -> {out.name}"
    )
    return out


# Modalities that produce symlinks. `sv_vcf` deliberately does not (the sv stage
# rewrites those under EIDs), so a cohort declaring only sv_vcf legitimately
# produces no manifest rows.
_LINKABLE_MODALITIES = ("bam", "fastq")


def _refuse_blind_prune(cohort: Cohort, samples: list[Sample], rows: list[tuple]) -> None:
    """Refuse to prune when samples exist but none of them yielded a file.

    Samples-but-no-rows cannot mean "the cohort emptied" — an empty cohort has no
    samples at all, and `discover_samples` raises in that case. It means discovery
    matched no files, and pruning on that basis would delete every link in the
    cohort's EID block and rewrite the manifest empty, with exit code 0.
    """
    if rows or not samples:
        return
    if not any(cohort.has(m) for m in _LINKABLE_MODALITIES):
        return
    raise RuntimeError(
        f"cohort {cohort.name!r}: {len(samples)} samples discovered but not one "
        f"linkable file among them, while the registry declares "
        f"{', '.join(m for m in cohort.modalities if m in _LINKABLE_MODALITIES)}. "
        f"That means discovery is broken, not that the cohort emptied — refusing "
        f"to prune {cohort.name!r}'s links or rewrite its manifest. Check "
        f"raw_root and the file-naming keys in the registry, then re-run "
        f"`tbooo restricted discover`."
    )


# A link's filename always starts with its EID followed by an underscore
# (see link_name): "<eid>_24048_0_0.bam", "<eid>_raw_<run>_<mate>...". Anything
# that doesn't match this shape is not one of ours and must be left alone.
_EID_PREFIX_RE = re.compile(r"^(\d+)_")


def _parse_eid(name: str) -> int | None:
    m = _EID_PREFIX_RE.match(name)
    return int(m.group(1)) if m else None


def _prune_orphaned_links(
    cfg: Config, cohort: Cohort, eids: dict[str, int], valid: set[Path]
) -> int:
    """Remove symlinks orphaned by EID renumbering.

    `assign_eids` numbers sequentially over sorted sample_ids, so adding or
    removing a sample shifts every later EID and leaves the previous highest
    EID's link behind: a real, resolvable symlink that appears in no manifest
    and silently duplicates another sample's genome under a stale identity.

    A 2-digit EID-prefix directory is a 100,000-wide bucket, but a cohort's
    EID block (BLOCK_SIZE) is 1,000,000 wide — ten buckets. Two cohorts with
    different `eid_base`s can therefore share a prefix directory. Directory
    membership alone is NOT enough to prove a link belongs to this cohort:
    every candidate's EID is parsed out of its filename and must fall inside
    THIS cohort's own [eid_base, eid_base + BLOCK_SIZE) block before it is
    even considered for removal. Anything outside that block, or whose name
    doesn't parse as an EID, is left strictly alone — it may belong to a
    different cohort sharing this prefix directory. Only symlinks are ever
    removed (never a real file), and the target is never followed/resolved
    before deciding. No scanning happens at all when `eids` is empty.
    """
    if not eids:
        return 0

    lo = cohort.eid_base
    hi = cohort.eid_base + BLOCK_SIZE
    prefixes = {eid_prefix_dir(eid) for eid in eids.values()}
    bases = []
    if cohort.has("bam"):
        bases.append(cfg.wgs_dir())
    if cohort.has("fastq"):
        bases.append(cfg.raw_reads_dir())

    pruned = 0
    for base in bases:
        for prefix in sorted(prefixes):
            d = base / prefix
            if not d.is_dir():
                continue
            for item in sorted(d.iterdir()):
                if not item.is_symlink():
                    continue
                item_eid = _parse_eid(item.name)
                if item_eid is None or not (lo <= item_eid < hi):
                    continue  # not this cohort's EID block — never touch
                if item in valid:
                    continue
                item.unlink()
                pruned += 1
                log(f"  pruned orphaned link: {item}")
    return pruned


def _ensure_link(dest: Path, target: Path) -> str:
    """Create the symlink, repairing it if it is missing, stale or broken."""
    target = target.resolve()
    if dest.is_symlink():
        try:
            if dest.readlink() == target and dest.resolve().exists():
                return "ok"
        except OSError:
            pass
        dest.unlink()
        dest.symlink_to(target)
        return "repaired"
    if dest.exists():
        # A real file where a link belongs — refuse rather than clobber data.
        raise RuntimeError(f"{dest} exists and is not a symlink; refusing to replace")
    dest.symlink_to(target)
    return "created"
