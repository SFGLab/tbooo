"""Deterministic EID assignment for registry-defined cohorts.

Mirrors pipeline/eids.py (1KGP/SGDP) but is driven entirely by the registry.
Assignment is by sorted sample_id so that re-running produces identical EIDs.
"""

from __future__ import annotations

from pathlib import Path

from tbooo.cohorts.registry import Cohort, Sample
from tbooo.config import Config
from tbooo.utils import ensure_dirs, log


def eid_map_path(cfg: Config, cohort: Cohort) -> Path:
    return cfg.metadata_dir() / f"eid_map_{cohort.name}.tsv"


def rename_path(cfg: Config, cohort: Cohort) -> Path:
    return cfg.metadata_dir() / f"vcf_sample_rename_{cohort.name}.txt"


def assign_eids(cohort: Cohort, samples: list[Sample]) -> dict[str, int]:
    """Map sample_id -> EID, assigned sequentially from the cohort's base."""
    ordered = sorted(s.sample_id for s in samples)
    return {sid: cohort.eid_base + i for i, sid in enumerate(ordered)}


def write_eid_map(
    cfg: Config, cohort: Cohort, samples: list[Sample], eids: dict[str, int]
) -> Path:
    ensure_dirs(cfg.metadata_dir())
    out = eid_map_path(cfg, cohort)
    by_id = {s.sample_id: s for s in samples}
    with out.open("w") as f:
        f.write("eid\tsample_id\tcohort\tn_files\thas_sv\n")
        for sid, eid in sorted(eids.items(), key=lambda kv: kv[1]):
            sample = by_id[sid]
            has_sv = any(x.modality == "sv_vcf" for x in sample.files)
            f.write(f"{eid}\t{sid}\t{cohort.name}\t{len(sample.files)}\t{int(has_sv)}\n")
    log(f"  wrote {out} ({len(eids)} samples)")
    return out


def write_rename_file(cfg: Config, cohort: Cohort, eids: dict[str, int]) -> Path:
    """bcftools reheader --samples map: '<original_id>\\t<eid>' per line."""
    ensure_dirs(cfg.metadata_dir())
    out = rename_path(cfg, cohort)
    with out.open("w") as f:
        for sid, eid in sorted(eids.items(), key=lambda kv: kv[1]):
            f.write(f"{sid}\t{eid}\n")
    log(f"  wrote {out}")
    return out
