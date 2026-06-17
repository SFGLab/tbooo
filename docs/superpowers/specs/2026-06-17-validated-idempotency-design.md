# Validated Idempotency — Design

Date: 2026-06-17
Status: Approved (pending written-spec review)

## Problem

`tbooo run` crashed in `_merge_pvcfs` (chr5) with `bcftools merge` exit status 255:

```
subprocess.CalledProcessError: Command '['bcftools', 'merge', ...,
  'data/Bulk/Whole genome sequences/ukb23370_c5_b0_v1.pvcf.gz',
  'data/raw/1kg/pvcf/nygc_c5.pvcf.gz', 'data/raw/sgdp/pvcf/sgdp_c5.pvcf.gz']'
  returned non-zero exit status 255.
```

Two intertwined root causes:

1. **Failures are undiagnosable.** `utils.run()` runs subprocesses with `capture=False`, so the
   raised `CalledProcessError` carries no stderr. The actual `bcftools` message (which explains exit
   255) is lost in interleaved parallel output. We cannot confirm the precise trigger for chr5 from
   the traceback alone — exit 255 from `bcftools merge` is almost always a missing/stale index,
   incompatible contigs (`chr5` vs `5`), or a truncated input.

2. **Idempotency is `path.exists()`-based, which trusts partial/corrupt files.** Every skip guard in
   `wgs.py` is `if out.exists(): skip`. A run killed mid-write or mid-index (likely on a multi-week
   build) leaves a half-written `.pvcf.gz` or a missing/stale `.tbi`. The next run *skips* rebuilding
   it, then feeds the broken intermediate to `bcftools merge` → exit 255. Because one chromosome
   logged "pVCF build complete" before chr5 threw, cross-cohort contigs are likely fine — which points
   at a stale/partial chr5 intermediate or index being trusted.

   Conversely, `array.py` / `imputed.py` / `wes.py` have **no** per-chromosome skip guard at all —
   they rebuild every chromosome on every run (wasteful, but never trust stale data).

The requirement: `tbooo run` must be idempotent — re-running checks the sanity of completed work and
rebuilds only what is missing or broken, never the whole dataset from scratch.

## Goals / Non-goals

- **Goal:** Re-running skips artifacts that are present *and* structurally valid & complete.
- **Goal:** Detect partial/corrupt/unindexed artifacts and rebuild *only those*, then re-validate.
- **Goal:** Make subprocess failures diagnosable (surface stderr in the exception).
- **Non-goal:** Deep content validation by default. Full `bgzip -t` decompression stays opt-in behind
  the existing `--deep-check` flag.
- **Non-goal:** A standalone `tbooo check` command (explicitly declined — validation is inline only).
- **Non-goal:** Reworking download logic (already has resume + integrity handling).

## Design

### 1. New module `tbooo/integrity.py` — cheap structural validators

Each validator answers "is this finished and intact?" in milliseconds. No full decompression.

| Validator | Artifacts | Check |
|---|---|---|
| `vcf_ok(path)` | array/imputed/wes tmp VCFs; all wgs pvcfs (nygc, sgdp, final) | exists + non-empty + BGZF EOF marker (`utils.has_bgzf_eof`) + index (`.tbi` or `.csi`) exists **and** index mtime ≥ data mtime |
| `plink_ok(stem)` | array `ukb22418_c*`, wes `ukb23157_c*` | `.bed`, `.bim`, `.fam` all exist & non-empty + `.bed` begins with PLINK magic bytes `6c 1b 01` |
| `bgen_ok(stem)` | imputed `ukb22828_c*`, wes bgen `ukb23157_c*` | `.bgen` non-empty + `.bgen.bgi` exists & mtime ≥ `.bgen` + `.sample` non-empty |
| `parquet_ok(path)` | phenotypes / geuvadis `participant.parquet` | exists + non-empty + `PAR1` magic at head and tail |
| `table_ok(path)` | eid maps, rename `.txt`, qc `.txt`, geuvadis `.tsv` | exists + non-empty + at least one data line |
| `symlink_ok(path, target_check)` | wgs crams/gvcfs, single-source merge symlinks | is a symlink + target resolves (not dangling) + target passes `target_check` |

Notes:
- **Index freshness via mtime** is the key signal for VCF completeness: a build writes the data file,
  then indexes it, so a fresh index implies a complete data file. A missing or older index means the
  prior run died between the two steps.
- For symlinked VCF outputs (single-cohort merge), mtime comparison is unreliable; validate by
  resolving the link and running `vcf_ok` on the *target*.

### 2. One orchestration helper `ensure(...)`

```python
def ensure(label, *, check, purge, build):
    """Validated-idempotent build step.

    check():  -> bool  True if the existing artifact is valid & complete.
    purge():  -> None  Delete partial artifact + sidecars (.tbi/.csi/...) + stale tmp.
    build():  -> None  Produce the artifact.
    """
    if check():
        log(f"  skip (valid): {label}")
        return
    purge()              # remove anything partial so the rebuild starts clean
    build()
    if not check():      # build ran but output is still bad — fail loudly, do not continue
        raise RuntimeError(f"{label}: rebuilt but still failed validation")
```

This single helper:
- replaces every `if out.exists(): skip` in `wgs.py`,
- adds the missing skip-guards to `array.py` / `imputed.py` / `wes.py`,
- guarantees we never silently propagate a broken artifact downstream.

`purge()` for each artifact removes the primary file(s) plus all sidecars (`.tbi`, `.csi`, `.bgi`,
PLINK triple, and the matching `tmp_dir` intermediates) so the rebuild cannot reuse stale state.

### 3. `utils.run()` upgrade — surface stderr on failure

On non-zero exit, the raised error must include captured stderr. Approach: when `check=True` and not
explicitly capturing, capture stderr to a buffer while still echoing it to the console (tee), and on
`CalledProcessError` append the captured stderr tail to the exception message. Long-running commands
keep their live console output; the difference is the *exception* now contains the diagnostic text.

### 4. Stage-by-stage application

- **eids:** wrap each `eid_map_*.tsv` / `vcf_sample_rename_*.txt` with `ensure` + `table_ok`.
- **array:** wrap per-chrom PLINK output with `ensure` + `plink_ok`; wrap the filtered tmp VCF with
  `vcf_ok` (it already does a 0-variant diagnosis — keep it).
- **imputed:** wrap per-chrom BGEN output with `ensure` + `bgen_ok`; shared `.sample` with `table_ok`.
- **wes:** wrap PLINK output (`plink_ok`) and BGEN output (`bgen_ok`) per chrom.
- **wgs:** wrap `nygc_pvcf`, `sgdp_pvcf`, final merge (`vcf_ok` / `symlink_ok`); crams & gvcfs with
  `symlink_ok`.
- **geuvadis:** wrap `geuvadis_expression_pcs.tsv` / `geuvadis_pca_variance.tsv` (`table_ok`); the
  parquet patch validates with `parquet_ok`.
- **phenotypes:** wrap `participant.parquet` with `parquet_ok`.
- **qc:** wrap `ukb_sqc_v2.txt` / `ukb_rel.txt` (`table_ok`); the merged-PLINK / het / king tmp
  artifacts already guard on existence — upgrade to `plink_ok` / `table_ok`.

### 5. Effect on the reported crash

With these changes, on re-run chr5's stale/partial intermediate (or its missing index) fails
`vcf_ok`, is purged, and rebuilt for that chromosome only; the final merge then succeeds. If the true
cause is instead a genuine contig mismatch, the rebuild won't mask it — the now-captured `bcftools`
stderr surfaces the real error immediately instead of a bare exit 255.

## Verification

The project is a one-shot data-bank build with thin business logic and no test suite; a formal
test/TDD pass is overkill here. Correctness is verified by implication and by a runtime self-check
rather than unit tests:

- **Built-in re-validation.** `ensure()` re-runs `check()` after `build()` and raises `RuntimeError`
  if the artifact is still invalid. This is the runtime safety net: a build that produces garbage
  fails loudly at the point of failure instead of being trusted downstream. It is the behavioral
  guarantee unit tests would otherwise assert, enforced live on real data.
- **Validators are pure and self-evident.** Each `*_ok` predicate is a few lines of stat/magic-byte
  checks with no hidden state, reviewable by reading.
- **Smoke verification after implementation:** re-run an already-completed stage and confirm it logs
  `skip (valid)` for intact artifacts; corrupt one artifact by hand (truncate its bytes or delete its
  index) and confirm that exact artifact is purged and rebuilt while the rest are skipped. Confirm a
  forced non-zero subprocess surfaces stderr in the raised error.

## Risks

- mtime-based index freshness assumes builds write data-then-index (true throughout the codebase) and
  that the filesystem preserves mtime ordering. Acceptable for this single-host pipeline.
- Cheap checks won't catch mid-stream corruption that still ends with a valid BGZF EOF block;
  `--deep-check` remains the tool for that and is unchanged.
