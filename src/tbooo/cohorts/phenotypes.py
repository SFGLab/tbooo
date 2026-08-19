"""Build the restricted instance's participant table.

The only phenotype available from the source trees is cohort membership, which is
mapped onto UKB Field 41270 (ICD-10 diagnoses) so the case/control axis lands in a
field real analysis code already reads. Anything richer arrives through the
clinical ingest hook.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from tbooo.cohorts.registry import Cohort, Sample
from tbooo.config import Config
from tbooo.utils import ensure_dirs, log

_UNKNOWN_PLATFORM = "unknown"


@dataclass(frozen=True)
class CohortData:
    cohort: Cohort
    samples: list[Sample]
    eids: dict[str, int]


def participant_path(cfg: Config) -> Path:
    return cfg.showcase_dir() / "participant.parquet"


def clinical_path(cfg: Config, cohort: Cohort) -> Path:
    return cfg.metadata_dir() / "clinical" / f"{cohort.name}.tsv"


def flowcell_of(cohort: Cohort, run: str) -> str:
    """Extract the flowcell token from a run token using the registry regex."""
    if not cohort.flowcell_regex:
        return run
    m = re.match(cohort.flowcell_regex, run)
    # The 'flowcell' group is guaranteed present by registry validation.
    return m.group("flowcell") if m else run


def platform_of(cohort: Cohort, flowcell: str) -> str:
    """First matching platform pattern wins; unmatched flowcells are 'unknown'."""
    for pattern, label in cohort.platform_patterns:
        if re.search(pattern, flowcell):
            return label
    return _UNKNOWN_PLATFORM


def _sample_batch_key(cohort: Cohort, sample: Sample) -> tuple[str, str, str]:
    """(cohort, platform, flowcell) for a sample, from its first FASTQ run."""
    runs = sorted({f.run for f in sample.files if f.run})
    if not runs:
        return (cohort.name, _UNKNOWN_PLATFORM, "")
    flowcell = flowcell_of(cohort, runs[0])
    return (cohort.name, platform_of(cohort, flowcell), flowcell)


def assign_batch_codes(
    cfg: Config, data: list[CohortData]
) -> dict[tuple[str, str, str], int]:
    """Assign a sequential integer per distinct (cohort, platform, flowcell)."""
    codes: dict[tuple[str, str, str], int] = {}
    for item in data:
        for sample in sorted(item.samples, key=lambda s: s.sample_id):
            key = _sample_batch_key(item.cohort, sample)
            if key not in codes:
                codes[key] = len(codes) + 1

    ensure_dirs(cfg.metadata_dir())
    out = cfg.metadata_dir() / "batch_codes.tsv"
    with out.open("w") as f:
        f.write("batch\tcohort\tplatform\tflowcell\n")
        for (cohort_name, platform, flowcell), code in sorted(
            codes.items(), key=lambda kv: kv[1]
        ):
            f.write(f"{code}\t{cohort_name}\t{platform}\t{flowcell}\n")
    log(f"  wrote {out} ({len(codes)} batches)")
    return codes


def putative_pedigree(cfg: Config, item: CohortData) -> Path | None:
    """Derive a PUTATIVE family grouping from sample-ID suffixes.

    This is filename inference, not kinship. It is written to a file whose name
    says 'putative', and it is never promoted to a relatedness file — that
    requires called genotypes, which the restricted instance does not have.
    """
    cohort = item.cohort
    if not cohort.pedigree_regex:
        return None

    rows: list[tuple[int, str, str, str]] = []
    for sample in sorted(item.samples, key=lambda s: s.sample_id):
        m = re.match(cohort.pedigree_regex, sample.sample_id)
        if not m:
            continue
        rows.append(
            (
                item.eids[sample.sample_id],
                sample.sample_id,
                m.group("stem"),
                m.group("member"),
            )
        )
    if not rows:
        return None

    ensure_dirs(cfg.metadata_dir())
    out = cfg.metadata_dir() / f"pedigree_{cohort.name}_putative.tsv"
    with out.open("w") as f:
        f.write("eid\tsample_id\tfamily_stem\tmember_suffix\n")
        for row in rows:
            f.write("\t".join(str(v) for v in row) + "\n")
    log(f"  wrote {out} ({len(rows)} samples) — PUTATIVE, inferred from sample IDs")
    return out


def _clinical_frame(cfg: Config, item: CohortData) -> pd.DataFrame | None:
    path = clinical_path(cfg, item.cohort)
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    key = df.columns[0]
    if key not in ("sample_id", "eid"):
        raise ValueError(
            f"{path}: first column must be 'sample_id' or 'eid', got {key!r}"
        )
    if key == "sample_id":
        df["eid"] = df["sample_id"].map(item.eids)
        missing = df["eid"].isna().sum()
        if missing:
            log(f"  WARN: {path.name}: {missing} rows have no matching sample")
        df = df.dropna(subset=["eid"]).drop(columns=["sample_id"])
    df["eid"] = df["eid"].astype("int64")
    df = df.rename(columns={c: f"clin_{c}" for c in df.columns if c != "eid"})
    log(f"  joined {path.name} ({len(df)} rows, {len(df.columns) - 1} columns)")
    return df


def build_participants(cfg: Config, data: list[CohortData]) -> Path:
    """Write participant.parquet for every sample across every cohort.

    This REPLACES participant.parquet and batch_codes.tsv wholesale, so `data`
    must be the full registry — passing a subset drops the omitted cohorts. The
    CLI enforces that by rejecting `--only` on this stage.
    """
    ensure_dirs(cfg.showcase_dir(), cfg.metadata_dir())
    codes = assign_batch_codes(cfg, data)

    rows: list[dict] = []
    for item in data:
        cohort = item.cohort
        for sample in sorted(item.samples, key=lambda s: s.sample_id):
            rows.append(
                {
                    "eid": item.eids[sample.sample_id],
                    "p31": pd.NA,  # sex — requires variant calling
                    "p22000": codes[_sample_batch_key(cohort, sample)],
                    "p41270": cohort.icd10,
                    "tbooo_cohort": cohort.name,
                    "tbooo_has_sv": any(f.modality == "sv_vcf" for f in sample.files),
                }
            )
        putative_pedigree(cfg, item)

    df = pd.DataFrame(rows).sort_values("eid").reset_index(drop=True)
    df["p31"] = df["p31"].astype("Int64")
    df["p41270"] = df["p41270"].astype("string")

    # Concatenate first, then merge once: two cohorts sharing a clinical column
    # name would otherwise collide into clin_age_x / clin_age_y.
    frames = [f for f in (_clinical_frame(cfg, item) for item in data) if f is not None]
    if frames:
        df = df.merge(pd.concat(frames, ignore_index=True), on="eid", how="left")

    out = participant_path(cfg)
    df.to_parquet(out, index=False)
    log(f"  wrote {out} ({len(df)} participants, {len(df.columns)} columns)")
    return out
