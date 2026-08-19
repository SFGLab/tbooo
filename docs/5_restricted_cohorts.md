# Restricted Cohort Ingestion

TBOOO's public instance is built from open datasets. Some deployments also need to
expose cohorts that cannot be redistributed — data held on an institutional cluster
under access control. This module ingests such cohorts into a *separate* TBOOO
instance without copying, moving, or modifying the source data.

## Two instances

A restricted instance is an ordinary TBOOO instance with `data_dir` pointing
somewhere else. It has its own `Bulk/`, `Showcase/` and `metadata/`. The public
instance is untouched and stays shareable; access control applies to the restricted
root as a single unit.

Keeping them separate is deliberate: merging restricted records into the public
`participant.parquet` would make the entire table restricted, including its public
rows.

## The registry

All cohort-specific knowledge lives in `cohorts.yaml`, which is **gitignored**. The
repository ships `cohorts.example.yaml` documenting the schema with placeholder
values. Code in `src/tbooo/cohorts/` knows only how to interpret a registry.

Each entry declares: an EID base, source roots, a discovery rule, which modalities
exist, and optional regexes for extracting flowcell, platform, FASTQ run/mate
tokens and putative family structure from filenames.

### Discovery rules

| Rule | Layout | `sample_id` |
|---|---|---|
| `flat_bam` | one directory of `<id>.bam` (+ optional `.bam.bai`) | BAM basename |
| `sample_dirs` | one subdirectory per sample containing paired FASTQ | directory name |

### Modalities

`bam`, `fastq`, `sv_vcf`. A cohort declares whichever it has. A declared modality
that matches no file for some sample is an error, not an empty result: discovery
fails and names the sample and the directory. (`sv_vcf` is the exception — a
sample may legitimately have no variant calls.)

### FASTQ naming

By default a FASTQ is recognised as `<run>_1.fq.gz` / `<run>_2.fq.gz`, and
`.fastq.gz` equally. Cohorts that name their reads differently set the optional
`fastq_pattern` key: a regex matched against the whole filename, with named groups
`run` and `mate`. That covers Illumina bcl2fastq output
(`<sample>_S1_L001_R1_001.fastq.gz`) and anything else. The registry rejects a
`fastq_pattern` missing either group. Omitting the key keeps the default.

## Stages

The stages are **ordered**, and each depends on the ones before it:

```
tbooo restricted discover      # 1. enumerate samples (read-only)
tbooo restricted eids          # 2. assign EIDs, write rename files
tbooo restricted link          # 3. symlink tree + provenance manifests
tbooo restricted sv            # 4. EMIT a SLURM script (runs nothing)
tbooo restricted phenotypes    # 5. participant.parquet
```

`discover` is a read-only preview and can be run at any time. `eids` must precede
`link`, `sv` and `phenotypes`, because all three name their outputs by EID. `sv`
additionally verifies that the rename file `eids` wrote still matches the EIDs it
was handed, and refuses to emit a script otherwise — a stale rename file would
produce a VCF whose filename says one EID and whose sample column says another.

`discover`, `eids`, `link` and `phenotypes` are metadata-only. `sv` writes an
sbatch script and prints the submit command — it never executes a subprocess,
because login nodes are typically for job submission only.

`--only <cohort>` restricts a stage to one cohort. It is accepted by `discover`,
`eids`, `link` and `sv`, all of which write strictly per-cohort outputs. It is
**rejected** by `phenotypes`: `participant.parquet` and `batch_codes.tsv` are
single tables covering every cohort and are rewritten in full on each run, so
building them from a subset would delete the other cohorts' rows.

## Re-running after the sample set changes

This is the most surprising behaviour in the system, and the one to understand
before re-running anything.

EIDs are assigned **by position**: samples are sorted by `sample_id` and numbered
sequentially from the cohort's `eid_base`. Assignment is therefore deterministic
and reproducible, but it is *not* stable across a change to the sample set. Adding
or removing a single sample shifts the EID of every sample that sorts after it.

Consequences of re-running the stages after the source tree gained or lost a
sample:

- **`eids` renumbers.** The same person can come out with a different EID than
  last time. Any EID recorded outside this tree — in an analysis script, a results
  file, a note — is silently invalidated.
- **`link` prunes.** Links left behind at EIDs no longer assigned are deleted, so
  that no retired EID keeps resolving to real sequence data. Only symlinks inside
  the cohort's own EID block are ever removed, and only when the run produced a
  manifest; a run that discovers samples but no files refuses to prune at all,
  because that means discovery is broken rather than that the cohort emptied.
- **`sv` reports but never deletes.** Its per-sample outputs are real files, not
  links, so a retired EID's VCF is listed in a warning and in commented-out `rm`
  lines in the generated sbatch script, for an operator to review and run. Nothing
  removes real data automatically.
- **`phenotypes` rebuilds** `participant.parquet` from scratch, over every cohort.

If stable identifiers matter for a cohort, freeze its source tree before the first
`eids` run, and treat the written `eid_map_<cohort>.tsv` as the record of record.

## Representation

Bulk data is **symlinked** under UKB-shaped names; sources are never copied or
written to. Variant calls are the exception: they embed the original sample ID, so
they are rewritten under EIDs to preserve pseudonymization.

Every link is recorded in `metadata/source_manifest_<cohort>.tsv` with the source
size and mtime, making provenance auditable and link rot detectable.

## Deviations from UKB conventions

| Deviation | Reason |
|---|---|
| Field 24048 files may carry a `.bam` extension | UKB 24048 is CRAM; a source may be BAM |
| `Bulk/Raw sequencing reads/` | UKB does not release FASTQ; no field exists |
| SV layer uses a custom `sv` label | UKB has no per-participant SV field |

## Phenotypes

Cohort membership maps to UKB Field 41270 (ICD-10) via the registry's `icd10` key,
so the case/control axis lands in a field analysis code already reads. A cohort
with `icd10: null` is a control group.

Richer phenotypes arrive through the clinical hook: drop a TSV at
`metadata/clinical/<cohort>.tsv` keyed on `sample_id` or `eid`, and every other
column joins into the parquet with a `clin_` prefix.

## What this module does not do

It performs no alignment and no variant calling, so SNV-based UKB layers (array,
imputed, WES, WGS pVCF) are not populated for restricted cohorts. It writes no
relatedness file — kinship requires called genotypes. Where sample IDs suggest
family structure, a `pedigree_<cohort>_putative.tsv` is emitted and clearly labelled
as inference.
