# GEUVADIS transcriptional landscape via PCCA — design

**Date:** 2026-07-29
**Status:** approved, ready for implementation planning
**Deadline:** UKB meeting, week of 2026-08-03

## Goal

Replicate the PCCA (Principal Component Canonical Correlation Analysis) result of Brown, Bray &
Pachter, *Expression reflects population structure* (PLOS Genetics, 19 Dec 2018), on the
GEUVADIS dataset, and summarize the methodology of three related papers for presentation.

This concludes the GEUVADIS transcriptional arm of the project. The genetic arm (SNP-based
clustering) is being done separately by Antoni. Once both exist, the two vectors can be fused
into a single latent representation whose clusters should reproduce population structure as
well as the raw modalities do.

## Background: why the current approach is wrong

`src/tbooo/pipeline/geuvadis.py` currently runs plain PCA on
`GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz`. Two independent problems:

1. That file is **PEER-residualized** (`resk10` = residuals after removing 10 PEER factors) and
   quantile-normalized. PEER factors absorb population structure. The analysis removes the
   signal it is looking for.
2. The code applies `log2(x + 0.1)` to values that are already residuals centered near zero
   (observed range ≈ ±0.1). The transform is meaningless on that input.

Beyond the file choice, the paper's central negative result is that **plain PCA on expression
does not cluster by population at all**, even done correctly. PCA on genotype does. PCCA is the
method that recovers population structure from the expression side, by coupling expression PCs
to genotype PCs through CCA.

**Consequence:** PCCA is not a purely transcriptional method. It requires genotype PCs as an
anchor. The *deliverable* is transcriptional (a projection of the expression matrix), but
genotypes are a required input, independent of Antoni's track.

## The three papers

| Paper | Role here |
|---|---|
| Lappalainen et al., *Nature* 2013 (nature12531) | The GEUVADIS dataset itself: 462 LCL RNA-seq samples across CEU/FIN/GBR/TSI/YRI, sequenced at 7 labs. Source of the quantification matrices and the eQTL baseline. |
| Brown, Bray & Pachter, *PLOS Genetics* 2018 (pgen.1007841) | The method being replicated. PCA→CCA coupling of expression and genotype. |
| Fachrul et al., *Communications Biology* 2023 (s42003-023-05171-9) | RGStraP: infers genetic PCs *directly from RNA-seq reads* via STAR two-pass + GATK variant calling. Validated on a 376-sample Nepal cohort and on GEUVADIS. Summarized only — see "Out of scope". |

### Reference numbers from the PCCA paper

- 343 individuals, populations **FIN, GBR, TSI, YRI**; CEU excluded for biased expression patterns
- kallisto against GENCODE v27 protein-coding transcripts; transcript TPMs summed to genes
- transcripts with mean TPM < 0.1 dropped (58,012 of 95,659 retained); MHC and non-autosomal
  genes removed → **14,070 genes**
- genotypes: 1000 Genomes **Phase 1**, indels/non-biallelic/MAF<5% removed → 6,785,201 SNPs, **no
  LD pruning**
- confounder regression: **7 lab indicators + sex**
- PCCA dimensions: **30 expression PCs × 5 genotype PCs → 2 canonical components**
- **first two canonical correlations: 0.963 and 0.766**
- eQTL at α = 1e-6: **2,818 eGenes** (PCA covariates) vs **2,732 eGenes** (PCCA covariates);
  GTEx replication 49.92% vs 51.02%

## Scope decisions

| Decision | Choice |
|---|---|
| Fidelity | Published GEUVADIS RPKM matrices + 1KGP Phase 3 genotype PCs. No FASTQ download, no kallisto re-quantification. |
| Genotype source | plink2's prebuilt `all_phase3` resource (build 37, 2504 samples). |
| Figures | All four panels: PCCA projection, QQ plot, eGene-count curve, Venn. |
| Code placement | Analysis in a `playground/` notebook; `tbooo map geuvadis` fixed in place. |
| Paper summary | `docs/4_pcca_methods.md` plus a rendered PDF. |

## Cohort configurations

Two configurations sharing one code path:

- **Replication config** (primary): the paper's **exact 343 individuals**. Selecting SDRF rows
  with `Characteristics[ancestry category] != 'Utah'` and
  `Comment[1000g Phase1 Genotypes] == '1'` yields exactly 343 individuals — British 85,
  Finnish 89, Tuscan 92, Yoruba 77 — reproducing the paper's cohort rather than approximating
  it. All 7 labs are represented.
- **Production config**: all 462 individuals, 5 populations including CEU. This is what feeds
  the later genetic+transcriptional fusion, where dropping a population would be wrong.

**Primary success criterion:** in the replication config, the first two canonical correlations
land near 0.963 and 0.766, and the projection separates YRI from the European populations with
visible FIN separation. Divergence means preprocessing is wrong, detectable before the meeting.

## Inputs

| Input | Source | Size |
|---|---|---|
| `GD660.GeneQuantRPKM.txt.gz` | `ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1/analysis_results/` | ~100 MB |
| `E-GEUV-1.sdrf.txt` | same host, experiment root | small |
| `all_phase3_chr{1..22}.{pgen,pvar}.zst`, `phase3_corrected.psam` | cog-genomics plink2 resources | 2.25 GiB total compressed |
| `integrated_call_samples_v3.20130502.ALL.panel` | already in `data/raw/1kg/` | present |

`E-GEUV-1.sdrf.txt` is the authoritative sample table and supplies **all** per-individual
covariates: `Characteristics[laboratory]` (values 1–7 — exactly the paper's 7 labs),
`Characteristics[sex]`, `Characteristics[ancestry category]`, `Comment[1000g Phase1 Genotypes]`,
and `Assay Name`, which joins directly to GD660 column headers.

**Lab is not recoverable from the GD660 header.** Column names look like `HG00096.1.M_111124_6`,
but all 660 runs share the prefix `M` (it denotes the sequencing platform, not the lab), so the
SDRF join is mandatory rather than a cross-check.

GD660 has 660 columns: the 462-individual nonredundant set plus 198 technical replicates. The
SDRF lists **one run per individual** for all 462, so the analysis selects those 462 `Assay Name`
columns directly. No replicate averaging is required.

GD660's `Chr` and `Coord` columns are GENCODE v12 on GRCh37, consistent with both genotype
sources.

## Environment prerequisites

- Python: `/opt/anaconda3/bin/python` (numpy 1.26, pandas 1.5.3, scikit-learn 1.8, scipy 1.17,
  matplotlib, pyarrow, jupyter all present). The system `python3` lacks pandas — do not use it.
- **plink2 is not installed.** Install the macOS ARM64 build and record the path in
  `config.yaml` under `tools.plink2`.
- **Disk headroom is ~22 GiB.** The full decompressed `all_phase3.pgen` plus its pvar would
  consume most of that. Genotypes must therefore be processed **per chromosome**: decompress one
  chromosome, immediately subset to the GEUVADIS samples with MAF and biallelic-SNP filters,
  delete the decompressed intermediate, then move to the next. Peak temporary usage stays near
  1 GB. The final filtered 462-sample callset is a few hundred MB.

## Pipeline

### Step 1 — Fix `tbooo map geuvadis`

Change the stage's input from the `50FN.resk10` file to `GD660.GeneQuantRPKM.txt.gz`:

- add GD660 to `src/tbooo/download/geuvadis.py`; the `50FN.resk10` file stays downloaded but is
  no longer read by any stage (it remains the correct input for eQTL work that wants PEER
  correction, so it is not removed)
- select the 462 nonredundant `Assay Name` columns via the SDRF, dropping the 198 technical
  replicates
- filter genes to mean RPKM ≥ 0.1, autosomal, non-MHC
- `log2(RPKM + 0.1)` — now applied to genuine RPKM, where it is meaningful
- PCA as before, writing `geuvadis_pc1..20` and the variance table

Sample identity stays keyed on the 7-character 1KGP sample ID (the first dot-field of the assay
name), so the existing join to `eid_map_1kg.tsv` is unaffected.

This makes the existing stage's output the honest "plain expression PCA" baseline — the paper's
Figure 1A negative result — instead of PCA on PEER residuals. Existing outputs
(`geuvadis_expression_pcs.tsv`, `geuvadis_pca_variance.tsv`, the `geuvadis_pc*` columns in
`participant.parquet`) keep their schema, so nothing downstream breaks. The stage's existing
idempotency guard via `table_ok` is preserved; stale outputs from the old input must be deleted
on first run.

### Step 2 — Expression matrix (notebook)

1. select the cohort's `Assay Name` columns (343 for the replication config, 462 for production)
2. filter genes: mean RPKM ≥ 0.1, autosomal, non-MHC (chr6:29,000,000–33,000,000, matching
   `_keep_gene` in the reference implementation) → expect ≈14k genes
3. quantile-normalize across samples
4. center and scale each gene to mean 0, variance 1
5. regress out 7 lab indicators + sex, joined from the SDRF; keep residuals → N × G matrix

Because the nonredundant set has exactly one run per individual, the regression is
unambiguously per-individual and no replicate collapsing is involved.

### Step 3 — Genotype PCs

Per chromosome, under the disk constraint above:

```
plink2 --zst-decompress all_phase3_chr${i}.pgen.zst all_phase3_chr${i}.pgen
plink2 --pfile all_phase3_chr${i} \
       --keep geuvadis_samples.txt \
       --snps-only --max-alleles 2 --min-alleles 2 \
       --maf 0.05 --autosome \
       --make-pgen --out geuv_chr${i}
rm all_phase3_chr${i}.pgen
```

Then merge the per-chromosome sets and run `plink2 --pca 5`. No LD pruning — the paper does not
prune, and pruning would change the PC structure being replicated.

The filtered callset is also the input to the eQTL scan in Step 6, so it is kept.

#### Why not the paper's own Phase 1 genotypes

GEUVADIS hosts the exact Phase 1 genotypes the paper used, already subset to the 465
RNA-sequenced individuals, at
`E-GEUV-1/genotypes/GEUVADIS.chr{1..22}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz`.
Using them would eliminate the genotype-release deviation entirely.

They total **~50 GB** (chr1 alone is 3.8 GB), which this machine cannot accommodate. Phase 3 is
therefore the genotype source, and the Phase 1 vs Phase 3 difference is a **permanent, untestable
deviation** rather than a temporary one. It must be reported as such: if the canonical
correlations miss their targets, genotype release remains a candidate explanation that this
replication cannot rule out.

Phase 3 costs ~3.5 GB of downloading across all 22 chromosomes (chr1 272 MB, chr22 52 MB), each
deleted after filtering, with ~1.5 GB peak transient usage and ~500 MB retained.

### Step 4 — PCCA

Direct port of `pachterlab/PCCA` `pcca/util.py`:

```python
u_x, _ = pca(expression_residuals, dim=30)   # SVD of XXᵀ → left singular vectors
u_y     = orthonormalized genotype PCs (5)
l_x, l_y, rho = cca(u_x, u_y, dim=2)          # SVD of u_xᵀ u_y
```

`pca()` there is SVD of the sample-by-sample Gram matrix, which equals PCA on a centered and
scaled data matrix — satisfied because Step 2 centers and scales each gene.

**Implementation detail that will silently corrupt results if missed:** CCA via
`svd(u_xᵀ u_y)` requires both inputs to have **orthonormal** columns. `pca()` returns
orthonormal columns, but plink2 `.eigenvec` columns are orthogonal and *not* unit-norm. The
genotype PCs must be normalized (or re-orthonormalized via QR) before the CCA step.

Report `rho` and check against 0.963 / 0.766.

### Step 5 — Leave-one-out cross-validation

The naive form recomputes a 6M-SNP PCA 343 times. Instead:

- compute the genotype Gram/GRM (N×N) and the expression Gram (N×N) once
- leave-one-out PCA is the eigendecomposition of the (N−1)×(N−1) submatrix — exact, milliseconds
- re-run CCA on each holdout, project the held-out individual with `project_i_cca`

Reduces the experiment from hours to seconds while remaining exact.

### Step 6 — cis-eQTL scan, twice

- phenotype: the Step 2 expression matrix (batch-corrected but **structure retained** — removing
  structure here would defeat the comparison)
- cis window: ±1 Mb around each gene's `Coord`, GRCh37 throughout, consistent with the plink2
  build 37 callset
- **PCA strategy** covariates: 5 genotype PCs
- **PCCA strategy** covariates: first 5 components of *both* the gene and genotype projections
  (10 covariates). This asymmetry in covariate count is the paper's own; replicated as written
  and flagged in the methods doc.
- both strategies additionally include 10 expression PCs as a **PEER substitute** — see
  Deviations
- computation: residualize the genotype block against covariates once per chromosome, after
  which each gene's test reduces to a dot product; t-statistics from correlations
- eGene = gene with ≥1 cis SNP below α

Outputs: QQ plot of p-values for both strategies (panel B), eGene count vs α curve (panel C),
Venn of eGene sets at α = 1e-6 (panel D).

## Deviations from the published method

Every item here is recorded in `docs/4_pcca_methods.md` as an explicit replicated-vs-deviated
table.

| Paper | This replication | Rationale |
|---|---|---|
| kallisto TPM, GENCODE v27 protein-coding | GEUVADIS' own GEM/Flux gene RPKM (GD660), GENCODE v12 | Avoids a ~1 TB FASTQ download and multi-day requantification. |
| 1000 Genomes **Phase 1** | 1000 Genomes **Phase 3** | The paper's own Phase 1 callset is ~50 GB, which this machine cannot accommodate. Phase 3 costs ~3.5 GB. This deviation is permanent and cannot be tested away. |
| 10 PEER factors | Top 10 expression PCs | PEER is an unmaintained R package that is difficult to build; expression PCs are the standard documented surrogate. |

**Not a deviation:** the cohort is the paper's exact 343 individuals, recovered from the SDRF.

Because of the quantification and genotype-release differences, exact eGene counts will not
reproduce. The canonical correlations are the meaningful replication check, since they measure
the coupling the method is built on rather than the quantifier's idiosyncrasies.

## Risks

- **Canonical correlations come out materially below 0.963 / 0.766.** With the cohort now exact,
  the two candidate causes are quantification and genotype release — but only the first is
  testable here. For quantification, `GD660.TrQuantRPKM.txt.gz` on the same FTP allows summing
  transcript RPKMs to genes, reproducing the paper's transcript→gene step more closely; the
  expression PC count can also be swept around 30, since the paper's elbow was found on a
  different quantifier. Genotype release **cannot** be tested, because the Phase 1 callset is too
  large for this machine. If the numbers stay low after the expression-side mitigations, the
  honest report is the observed values plus genotype release named as an unresolved candidate —
  not a claimed match.
- **Disk exhaustion during the genotype step.** Mitigated by the per-chromosome design; the
  implementation must delete each decompressed `.pgen` before starting the next chromosome.
- **eQTL scan runtime.** ~14k genes × ~4k cis SNPs at n=343 is tractable once genotypes are
  residualized in chromosome-sized blocks, but a naive per-(gene, SNP) loop is not. The
  vectorized formulation is required, not optional.

## Deliverables

1. `playground/geuvadis_pcca.ipynb` — the analysis, run against `/opt/anaconda3/bin/python`
2. figures in `data/metadata/eda/pcca_*.png` — the four panels plus the Figure 1 baseline
3. `src/tbooo/pipeline/geuvadis.py` — corrected to use GD660 raw RPKM
4. `docs/4_pcca_methods.md` and a rendered PDF — three-paper methodology comparison plus the
   deviations table

## Out of scope

RGStraP (Communications Biology 2023) is summarized in the methods doc but **not executed**. It
calls variants from RNA-seq reads via STAR two-pass alignment and GATK, requiring the ~1 TB
E-GEUV-1 FASTQ set. It is the natural third arm once those reads are available, and would let us
compare RNA-derived genotype PCs against array-derived ones on the same individuals.

Fusion of the genetic and transcriptional vectors into a single latent representation is the
follow-on project, dependent on this work and on Antoni's SNP-based clustering.
