# Population structure from transcriptomes: three methodologies

This document summarizes the three papers that define the GEUVADIS transcriptional arm of this
project, and states exactly what our replication reproduces and where it departs.

The through-line is one question: **can you recover human population structure from gene
expression?** The three papers answer it differently. Lappalainen et al. built the dataset and
treated population structure largely as something to correct for. Brown, Bray & Pachter showed
that the structure is present in expression but invisible to naive PCA, and gave a method that
exposes it. Fachrul et al. sidestepped expression entirely and called genotypes from the RNA-seq
reads themselves.

---

## 1. Lappalainen et al., *Nature* 2013 — the GEUVADIS resource

*Transcriptome and genome sequencing uncovers functional variation in humans.* Nature 501:506–511.

### Design

Lymphoblastoid cell lines (LCLs) from 1000 Genomes individuals, mRNA and small RNA sequenced
across a seven-laboratory consortium.

| Property | Value |
|---|---|
| Individuals (mRNA, nonredundant) | 462 |
| QC-passed runs including technical replicates | 660 |
| Populations | CEU 91, FIN 95, GBR 94, TSI 93, YRI 89 |
| Sequencing laboratories | 7 |
| Annotation | GENCODE v12 |
| Genotypes | 1000 Genomes Phase 1 |
| Reference build | GRCh37 |

Population counts above are as recorded in `E-GEUV-1.sdrf.txt`, the authoritative sample table.

### Quantification and the normalization tiers

Expression was quantified with GEM/Flux Capacitor against GENCODE v12, and released at four
normalization levels. Understanding these levels matters, because choosing the wrong one silently
destroys the signal this project needs:

1. **Raw read counts** — no normalization.
2. **Library-depth scaled** — counts scaled by total mapped reads, adjusted to the sample-set
   median (45 M for mRNA).
3. **RPKM** — library depth and transcript length.
4. **Library depth + expressed + PEER** — as above, plus removal of units with zero counts in
   more than 50% of samples, plus removal of technical variation by PEER.

### Reading the filenames

The distributed filenames encode which tier a file belongs to:

| Fragment | Meaning |
|---|---|
| `GD660` | All 660 QC-passed runs, technical replicates included |
| `GD462` | The 462-individual nonredundant set used in most analyses |
| `50FN` | Filtered to units expressed in >50% of samples, then normalized |
| `resk10` | **Residuals after removing 10 PEER factors** |

> **This is the single most consequential fact in this document.** `GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz`
> is PEER-residualized. PEER factors absorb broad sources of sample-to-sample covariation, and
> population structure is exactly such a source. Running a population-structure analysis on that
> file removes the signal being looked for. Our pipeline previously did this; it now reads
> `GD660.GeneQuantRPKM.txt.gz` (tier 3, raw RPKM) instead.

### eQTL mapping

Cis-eQTLs were mapped separately in Europeans (`EUR373`) and Yoruba (`YRI89`), reported at FDR 5%,
with association strength given as a Spearman correlation derived from the linear regression
slope. These per-population files remain a useful external reference for any eQTL numbers we
generate.

### Relevance here

GEUVADIS is our input dataset. The paper treats population differences primarily as a covariate
to control, not as a target of inference — which is precisely the gap the next paper addresses.

---

## 2. Brown, Bray & Pachter, *PLOS Genetics* 2018 — PCCA

*Expression reflects population structure.* PLOS Genetics, 19 December 2018.
Code: [github.com/pachterlab/PCCA](https://github.com/pachterlab/PCCA)

### The central claim

Principal component analysis of a gene expression matrix **does not** cluster individuals by
population — the leading components are driven by batch and other technical variation. PCA of the
genotype matrix for the same individuals **does**. The structure is present in expression, but
naive PCA cannot see it.

PCCA recovers it by coupling the two views: take principal components of each, then find the
directions in which they are maximally correlated.

### Cohort

343 individuals: **FIN, GBR, TSI and YRI**. CEU is excluded for biased expression patterns.

This cohort is exactly reproducible from the SDRF — non-CEU individuals carrying the
`1000g Phase1 Genotypes` flag number precisely 343 (GBR 85, FIN 89, TSI 92, YRI 77).

### Expression preprocessing

1. kallisto against GENCODE v27 protein-coding transcripts
2. transcripts with mean TPM < 0.1 dropped — 58,012 of 95,659 retained
3. transcript TPMs summed to gene level
4. MHC-region and non-autosomal genes removed — **14,070 genes** remain
5. quantile normalization, then each gene centered and scaled to mean 0, variance 1
6. regression on a confounder matrix of **7 laboratory indicators + sex**, keeping residuals

### Genotype preprocessing

1000 Genomes Phase 1 in PLINK format; indels, non-biallelic variants and MAF < 5% removed,
leaving **6,785,201 SNPs**. Notably, **no LD pruning** is performed.

### The algorithm

From `pcca/util.py`, the method is compact:

```
u_x, _ = pca(expression, dim=30)      # SVD of XXᵀ → left singular vectors
u_y, _ = pca(genotype,   dim=5)
l_x, l_y, rho = cca(u_x, u_y, dim=2)  # SVD of u_xᵀ u_y
```

Two details carry the whole method:

- `pca` decomposes the **sample-by-sample Gram matrix** rather than the feature covariance. With
  14,070 genes or 6.8 M SNPs against a few hundred individuals, this is the only tractable
  direction — and it is what makes leave-one-out cross-validation cheap.
- Because `u_x` and `u_y` have orthonormal columns, the canonical correlations are exactly the
  singular values of `u_xᵀ u_y`. Feeding in non-orthonormal columns (for example, raw plink2
  `.eigenvec` values, which are orthogonal but not unit-norm) produces wrong correlations with no
  error raised.

Dimensions 30 / 5 / 2 were selected by the elbow method.

### Headline results

| Result | Value |
|---|---|
| First canonical correlation | **0.963** |
| Second canonical correlation | **0.766** |
| eGenes at α = 1e-6, genotype-PC covariates | 2,818 |
| eGenes at α = 1e-6, PCCA covariates | 2,732 |
| GTEx replication rate, genotype-PC strategy | 49.92% (1,407 / 2,818) |
| GTEx replication rate, PCCA strategy | 51.02% (1,394 / 2,732) |

Leave-one-out cross-validation removes each of the 343 individuals in turn, recomputes both sets
of principal components and the canonical variables on the remainder, and projects the held-out
individual. Individuals land near their own population despite never contributing to the fit,
which is what distinguishes real structure from overfitting.

The eQTL comparison uses 10 PEER factors in both strategies, plus either the 5 genotype PCs or
the first 5 components of **both** the gene and genotype PCCA projections. The covariate counts
are therefore asymmetric (5 vs 10); this is the paper's own design.

The conclusion is deliberately modest: PCCA covariates find slightly *fewer* eGenes, but those
eGenes replicate in GTEx at a slightly *higher* rate.

### Relevance here

This is the method we replicate. It is also the reason the transcriptional arm cannot be run in
isolation: **PCCA requires genotype PCs as an anchor.** The deliverable is a transcriptional
embedding, but genotypes are a required input.

---

## 3. Fachrul et al., *Communications Biology* 2023 — RGStraP

*Direct inference and control of genetic population structure from RNA sequencing data.*
Communications Biology, 2 August 2023.

### The central claim

If you have RNA-seq reads, you already have genotype information — the reads cover transcribed
variants. RGStraP calls those variants and computes genetic PCs from them, so population
structure can be controlled for without any array or WGS data.

### Pipeline

1. **Alignment:** STAR two-pass to GRCh38; merge with Picard; filter MAPQ < 20 with SAMtools
2. **Variant calling:** GATK4, following the recommended RNA-seq discovery workflow
3. **SNP filtering (PLINK):** remove duplicated and palindromic variants; autosomal SNPs only;
   drop individuals with >20% missingness; remove the HLA region; MAF > 0.05; LD r² < 0.05
4. **Refinement:** intersect with HapMap3 to retain well-established SNPs
5. **PCA:** FlashPCA on the filtered set

### Cohorts and results

| | Nepal (STRATAA) | GEUVADIS (validation) |
|---|---|---|
| Samples | 376 whole blood, 362 post-QC | 465 LCL, 463 post-filtering |
| Paired array genotypes | 299, 280 post-QC | 1000 Genomes Phase 3 |
| Variants called | 4,921,472 | — |
| After MAF and LD filters | 36,440 | — |
| After HapMap3 intersection | 4,887 | — |

Genotype concordance against arrays averaged **0.925**, with 232 of 280 samples (82.8%) above
0.90. Canonical correlation between RNA-derived and array-derived PCs reached Rc1 = 0.946 in
Nepal (95.6% shared variance over three canonical variates) and Rc1 = 0.994, Rc2 = 0.942,
Rc3 = 0.752 in GEUVADIS (99.8% cumulative).

Using RNA-derived PCs as covariates in differential expression reduced test-statistic inflation
from m = 0.935 to levels matching array PCs (m = 0.92), with 272 genes at FDR < 0.05 against 144
using array PCs and 325 using none.

### Why we summarize but do not run it

RGStraP operates on **reads**, not quantification matrices. Replicating it requires the E-GEUV-1
FASTQ set, roughly 1 TB, plus STAR and GATK runtime. That is out of scope here. It is the natural
third arm once those reads are available, and would let us compare RNA-derived genotype PCs
against array-derived ones on identical individuals.

---

## 4. Side-by-side

| | Lappalainen 2013 | Brown/Bray/Pachter 2018 | Fachrul 2023 |
|---|---|---|---|
| Input | RNA-seq + genotypes | Quantifications + genotypes | RNA-seq reads only |
| Genotypes required | Yes | **Yes** (as CCA anchor) | No — derived from reads |
| Population structure treated as | Covariate to control | Target of inference | Covariate to control |
| Recovered from | Genotypes | Expression, via coupling | Reads, via variant calling |
| Dimensionality reduction | PCA / PEER | PCA → CCA | PCA (FlashPCA) |
| Cost | Dataset generation | Minutes on a laptop | Alignment + variant calling |
| Replicated here | Input data | **Yes** | No — needs ~1 TB of reads |

---

## 5. What this replication does

**Cohort:** the paper's exact 343 individuals for replication; all 462 for the production
embedding that feeds the later genetic + transcriptional fusion.

**Success criterion:** the first two canonical correlations should approach 0.963 and 0.766. This
is the meaningful check — it measures the coupling the method is built on, rather than any
particular quantifier's idiosyncrasies.

### Replicated faithfully

- The exact 343-individual cohort, CEU excluded
- Gene filtering: mean expression ≥ 0.1, autosomal, non-MHC
- Quantile normalization, per-gene centering and scaling
- Confounder regression on 7 laboratory indicators + sex
- Biallelic SNPs at MAF ≥ 5%, no LD pruning
- PCCA dimensions 30 / 5 / 2, algorithm ported directly from `pcca/util.py`
- Leave-one-out cross-validation
- Both eQTL covariate strategies, including the paper's asymmetric covariate counts

### Deviations

| Paper | Here | Why |
|---|---|---|
| kallisto TPM, GENCODE v27 protein-coding | GEUVADIS GEM/Flux gene RPKM (GD660), GENCODE v12 | Avoids a ~1 TB FASTQ download and multi-day requantification |
| 1000 Genomes **Phase 1** | 1000 Genomes **Phase 3** | The paper's own Phase 1 callset is ~50 GB; this machine cannot accommodate it. **Permanent and untestable** — see below |
| 10 PEER factors | Top 10 expression PCs | PEER is an unmaintained R package that is difficult to build; expression PCs are the standard documented surrogate |

The genotype-release deviation deserves emphasis. GEUVADIS hosts the paper's actual Phase 1
genotypes, pre-subset to the 465 sequenced individuals, which would have removed this deviation
entirely — but at ~50 GB they exceed available disk. Consequently, **if the canonical
correlations fall short of 0.963 / 0.766, genotype release remains a candidate explanation this
replication cannot rule out.** Expression-side causes can be tested (transcript-level
quantification, confounder design, expression PC count); the genotype-side one cannot.

### Results

**The replication succeeds.**

| Quantity | Paper | This replication |
|---|---|---|
| First canonical correlation | 0.963 | **0.947** |
| Second canonical correlation | 0.766 | **0.782** |
| Individuals | 343 | 328 |
| Genes after filtering | 14,070 | 12,675 |
| SNPs after MAF ≥ 5%, biallelic | 6,785,201 | 6,695,351 |

The central claim reproduces cleanly. Measuring how much of the variance in each leading
coordinate is explained by population label:

| | Coordinate 1 | Coordinate 2 |
|---|---|---|
| PCA of expression | 0.048 | 0.071 |
| PCA of genotype | 1.000 | 0.975 |
| **PCCA projection of expression** | **0.897** | **0.605** |

Expression PCA carries essentially no population signal — 5% and 7%, consistent with the paper's
Figure 1A showing populations intermixed. The PCCA projection of the *same matrix* recovers 90%
and 61%. In the projection, YRI separates from the European populations along coordinate 1 and
FIN separates from GBR/TSI along coordinate 2, matching the paper's Figure 2A.

Leave-one-out cross-validation reproduces Figure 2B: individuals who contributed nothing to the
fit still land with their own population, so the structure is not an artefact of overfitting.

#### The eQTL arm does not reproduce the paper's direction

The cis-eQTL comparison ran genome-wide on 12,675 genes, ±1 Mb windows, under both covariate
strategies. The counts are the same order as the paper's, but **the direction of the difference
is reversed**:

| Cutoff | PCA covariates | PCCA covariates | Difference |
|---|---|---|---|
| 1e-4 | 5,434 | 5,441 | +7 |
| 1e-5 | 3,191 | 3,231 | +40 |
| **1e-6** | **2,277** | **2,288** | **+11** |
| 1e-7 | 1,764 | 1,794 | +30 |
| 1e-9 | 1,186 | 1,204 | +18 |

At α = 1e-6 the two sets overlap heavily: 2,157 shared, 120 unique to PCA, 131 unique to PCCA.

The paper reports PCCA finding slightly *fewer* eGenes (2,732 against 2,818) with *reduced*
inflation. Here PCCA finds slightly *more* at every cutoff, and genomic inflation is marginally
*higher*: λ = 1.119 for PCA covariates against λ = 1.127 for PCCA, over ~193,000 sampled cis
tests on chromosomes 1–2. The QQ curves are visually indistinguishable.

The differences are small in both the paper (3%) and here (0.5%), so this is a weak effect either
way rather than a contradiction of the structural result. The most likely cause of the reversal
is the **PEER substitution**: both strategies carry 10 expression PCs as the PEER stand-in, and
those PCs are drawn from the same 30-PC subspace the PCCA components come from. That shared
correction plausibly absorbs most of the confounding the two strategies are competing to remove,
leaving the residual difference near noise. Testing this would mean building PEER properly.

**What this means for the meeting:** the structural claim — that expression reflects population
structure and PCCA exposes it — replicates convincingly. The downstream claim that PCCA
covariates improve eQTL calibration does not, on this configuration.

**Additional cohort deviation discovered during execution.** The paper's 343 individuals are
recovered exactly from the SDRF, but 15 of them are absent from the 1000 Genomes Phase 3 build-37
release, which carries only the 2504 unrelated samples. The analysis therefore runs on **328**
individuals. The dropped samples span populations (GBR 8, FIN 3, YRI 2, TSI 2) rather than
concentrating in one, so this is unlikely to bias the structure.

#### Two implementation traps worth recording

Both produced plausible-looking output rather than errors, and both were caught by inspecting
figures rather than by tests passing.

1. **The leave-one-out projection needs `1/s²` whitening.** `project_i_cca` reduces to: take the
   held-out sample's PC scores in the *training* basis, whiten by the singular values, then apply
   the canonical rotation. A kernel-style weighted average of training coordinates looks
   reasonable and produces finite, correctly shaped output — with no population structure at all.
2. **Canonical vectors are only defined up to sign.** Each holdout refits independently, so
   projections come back at arbitrary polarity and a population scatters to *both* extremes of an
   axis. Each holdout must be sign-aligned against the full-data solution.

A third trap did not fire but is guarded: plink2 `.eigenvec` columns are orthogonal but **not**
unit-norm, and CCA via `svd(u_xᵀ u_y)` silently returns wrong correlations if fed non-orthonormal
input. `pcca.core.cca` raises rather than accepting them.

#### Figures

| File | Content |
|---|---|
| `pcca_fig1a_expression_pca.png` | Expression PCA — populations intermixed |
| `pcca_fig1b_genotype_pca.png` | Genotype PCA — four clean clusters |
| `pcca_fig2a_projection.png` | PCCA projection — the replication |
| `pcca_fig2b_loo.png` | Leave-one-out validation |
| `pcca_panelB_qq.png` | eQTL QQ, both covariate strategies |
| `pcca_panelC_curve.png` | eGenes vs significance cutoff |
| `pcca_panelD_venn.png` | eGene overlap at α = 1e-6 |

All in `data/metadata/eda/`. Analysis notebook: `playground/geuvadis_pcca.ipynb`.

---

## 5b. Extension to MAGE and GTEx

The GEUVADIS replication rested on a narrow cohort — four populations, one of them
non-European. Two further cohorts test whether the result generalizes along independent axes:
**ancestry breadth** (MAGE) and **tissue** (GTEx). Design:
`docs/superpowers/specs/2026-08-17-gtex-mage-pcca-design.md`.

| | GEUVADIS | MAGE | GTEx |
|---|---|---|---|
| Individuals | 328 | 729 | 601-816 per tissue |
| Tissue | LCL | LCL | 4 solid/blood tissues |
| Populations | 4 | **26** | unlabelled |
| Continental groups | 2 | **5, balanced** | EUR ~87%, AFR ~12% |
| First canonical correlation | 0.947 | **0.974** | 0.355-0.422 |
| Second canonical correlation | 0.782 | **0.795** | 0.247-0.277 |

### MAGE strengthens the result

MAGE spans 26 populations with no continental group dominating (AFR 196, EUR 142, EAS 141,
SAS 139, AMR 113). Canonical correlations are **0.974 and 0.795** — higher than GEUVADIS and
than the original paper's 0.963 / 0.766.

Variance in each leading coordinate explained by continental group:

| | Coordinate 1 | Coordinate 2 |
|---|---|---|
| PCA of expression | **0.035** | **0.004** |
| PCA of genotype | 0.968 | 0.945 |
| PCCA projection | 0.921 | 0.603 |

**This is the most informative single number in the extension.** Going in, the expectation was
that MAGE's far greater ancestry range might let plain expression PCA partly separate
populations, which would have meant the GEUVADIS negative result was an artefact of that
cohort's narrow composition. It does not: expression PCA explains 3.5% and 0.4%, even lower
than GEUVADIS's 4.8% and 7.1%. The negative result is not contingent on cohort composition.
Ancestry is present in expression at 731 individuals across five continental groups, and
ordinary PCA still cannot see it.

### GTEx: consistent direction, much weaker coupling

GTEx genotypes are dbGaP-protected, but the open-access eQTL covariates carry PC1-PC5 genotype
PCs, so PCCA is possible without authorization. Population labels are also protected, so the
contrast is measured against genotype PC1 as a continuous ancestry proxy rather than labels.

| Tissue | n | genes | rho1 | rho2 | expression PCA ~ gPC1 | PCCA ~ gPC1 |
|---|---|---|---|---|---|---|
| Whole blood | 800 | 13,857 | 0.399 | 0.251 | 0.132 | **0.391** |
| Muscle skeletal | 816 | 14,246 | 0.355 | 0.247 | 0.012 | **0.344** |
| Lung | 601 | 15,971 | 0.408 | 0.277 | 0.089 | **0.400** |
| Thyroid | 681 | 15,892 | 0.422 | 0.262 | 0.071 | **0.392** |

The direction is consistent in all four tissues: PCCA coordinates track ancestry 3-30x better
than the leading expression PCs. But the absolute coupling is roughly half that of the LCL
cohorts.

Two plausible reasons, neither tested here. GTEx is ~87% European with only ~12% African
ancestry, so there is far less ancestry variance available to couple to than in MAGE's balanced
design. And GTEx samples are post-mortem solid tissues, where cell-type composition, ischaemic
time, age and disease contribute expression variance that LCLs — a homogeneous cultured cell
type — do not have. Distinguishing these would need a cohort that varies one without the other.

The consistency across four tissues is itself the more useful finding: whatever drives the
weaker coupling is a property of the cohort rather than of any one tissue.

### Cross-cohort: one genotype PC space, four datasets

MAGE distributes genotype PCs generated jointly across MAGE, Geuvadis, GTEx and AFGR — 2,373
samples in a single coordinate system (GTEx 837, MAGE 731, AFGR 597, Geuvadis 452). Two things
follow.

**GTEx ancestry can be placed without dbGaP.** GTEx donors carry no population labels, but they
sit among labelled 1KGP samples in the shared space. A k-nearest-labelled-neighbour vote assigns
each donor a continental group at median confidence 1.00:

| Inferred group | Donors |
|---|---|
| EUR | 701 |
| AFR | 103 |
| AMR | 20 |
| EAS | 10 |
| SAS | 3 |

That is ~84% European and ~12% African, matching GTEx's documented composition — a reassurance
that the inference is behaving, not a validation of any individual assignment. These values
populate `p21000_i0` for GTEx rows and are flagged as inferred wherever they appear.

**The same people, two experiments.** 93 individuals are in both MAGE and Geuvadis; 91 survive
into both coordinate tables. Their PCCA coordinates can be compared directly, which is as close
to a reproducibility test as this design allows.

| | Overall | Within AFR (n=9) | Within EUR (n=82) |
|---|---|---|---|
| Coordinate 1 | **0.811** | 0.032 | 0.028 |
| Coordinate 2 | 0.273 | 0.434 | 0.329 |

The headline 0.811 is real but must not be over-read. Decomposing it shows the agreement is
**entirely** the between-group contrast: the two experiments agree on which continental group
each person belongs to, and within a group they agree essentially not at all (|r| ≈ 0.03).

So the embedding is reproducible at the level of ancestry group, not individual placement. For
the fused-representation goal that is the relevant guarantee — clusters should recover population
structure — but it does bound what a single individual's coordinate can be claimed to mean. Any
downstream use that treats per-individual PCCA position as a stable personal quantity is not
supported by this evidence.

## 6. What this enables

### The production embedding

The replication config deliberately drops CEU and the individuals without Phase 1 genotypes,
because that is what the paper did. For the fusion work we also ran the identical pipeline over
**all 462 GEUVADIS individuals** — 445 after intersecting with the Phase 3 panel — including CEU.

Canonical correlations there are **0.930 and 0.762**, marginally below the four-population
configuration. That is the expected direction: CEU is the population the paper excluded for
biased expression patterns, and adding it back dilutes the coupling slightly. The projection
still separates YRI cleanly on coordinate 1 and FIN on coordinate 2, with CEU overlapping GBR as
ancestry would predict.

Output: `data/metadata/geuvadis_pcca_coords.tsv` — 445 rows of

```
sample_id  pop  pcca_c1  pcca_c2  geno_pc1 … geno_pc5
```

This is the transcriptional vector the fusion work consumes, carrying its matched genotype PCs so
both modalities arrive on identical rows.

### Toward a fused representation

With the transcriptional vector defined, and Antoni's SNP-based clustering supplying the genetic
vector, the two modalities can be fused into a single latent representation. The test of that
representation is whether its clusters reproduce population structure as faithfully as the raw
modalities do individually.

One practical note for that work: the genotype side here is a filtered SNP table
(`data/raw/geuvadis/genotype/phase3/`, 6.7M biallelic autosomal SNPs at MAF ≥ 5%, no LD pruning).
If the genetic-clustering track consumes the same table, the two arms share one genotype input
rather than each building their own, and the fused rows are guaranteed to align.

PCCA is a natural precursor to that fusion, since it is already a two-view method — it finds the
shared subspace between expression and genotype. The fusion work generalizes this from a shared
subspace to a joint latent representation.

---

## References

1. Lappalainen T, Sammeth M, Friedländer MR, et al. Transcriptome and genome sequencing uncovers
   functional variation in humans. *Nature* 501:506–511 (2013).
   [doi:10.1038/nature12531](https://www.nature.com/articles/nature12531)
2. Brown BC, Bray NL, Pachter L. Expression reflects population structure. *PLOS Genetics*
   14(12):e1007841 (2018).
   [doi:10.1371/journal.pgen.1007841](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007841)
3. Fachrul M, Karkey A, Shakya M, et al. Direct inference and control of genetic population
   structure from RNA sequencing data. *Communications Biology* 6:804 (2023).
   [doi:10.1038/s42003-023-05171-9](https://www.nature.com/articles/s42003-023-05171-9)

Data: [GEUVADIS E-GEUV-1](https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1/) ·
Code: [pachterlab/PCCA](https://github.com/pachterlab/PCCA) ·
Design: [`docs/superpowers/specs/2026-07-29-geuvadis-pcca-design.md`](superpowers/specs/2026-07-29-geuvadis-pcca-design.md)
