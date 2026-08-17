# GTEx and MAGE integration with PCA/PCCA — design

**Date:** 2026-08-17
**Status:** implemented
**Builds on:** `docs/superpowers/specs/2026-07-29-geuvadis-pcca-design.md`

## Goal

Integrate two further transcriptomic cohorts into TBOOO — **MAGE** (731 LCL RNA-seq samples from
1000 Genomes individuals) and **GTEx** (bulk RNA-seq across tissues) — and run the PCA/PCCA
analysis already validated on GEUVADIS against both.

The GEUVADIS replication succeeded (canonical correlations 0.947 / 0.782 against the paper's
0.963 / 0.766) but rested on a narrow cohort: four populations, only one of them non-European.
MAGE and GTEx test whether the result generalizes along two independent axes — **ancestry
breadth** (MAGE spans 26 populations across five continental groups) and **tissue** (GTEx is not
LCLs).

## Why these two datasets

| | GEUVADIS | MAGE | GTEx |
|---|---|---|---|
| Individuals | 462 | 731 | 943 donors |
| Tissue | LCL | LCL | 50 tissues |
| Populations | 4 (+CEU) | **26** | unlabelled, ~2 groups |
| Continental groups | 2 | **5, balanced** | EUR ≫ AFR |
| Genotypes | 1KGP Phase 3 | 1KGP NYGC | dbGaP protected |
| Genotype PCs | computed by us | **published** | **published** |

MAGE's continental balance is AFR 196, EUR 142, EAS 141, SAS 139, AMR 113 — the first cohort here
where no group dominates.

## Key enabling discoveries

Three findings from exploration make this tractable; each was verified against the live data
rather than assumed.

### 1. GTEx genotype PCs are public

Individual GTEx genotypes require dbGaP authorization, which would normally rule out PCCA
entirely — the method needs genotype PCs as its CCA anchor.

But `GTEx_Analysis_v10_eQTL_covariates.tar` (20.8 MB, open access) contains **PC1–PC5 genotype
PCs** for every tissue. Verified: the PCs are **identical across all 50 tissues** for a shared
donor, confirming they are per-donor global genotype PCs rather than per-tissue quantities. The
union across tissues is **943 donors**.

PC1 is strongly skewed — median −0.0135, 95th percentile 0.088 — with roughly **119 donors** in
the far tail against **824** in the main mode. This matches GTEx's known composition of ~85%
European and ~12% African-American donors, so GTEx is effectively a two-group problem.

### 2. MAGE's 59.5 GB archive can be read selectively

MAGE ships as a single `MAGE.v1.0.data.zip` of 59,477,326,151 bytes on Zenodo. Zenodo serves HTTP
range requests (verified: `206 Partial Content`), and the archive is **Zip64**.

Parsing the Zip64 end-of-central-directory locator and central directory over range requests
enumerates all **111 members** without downloading the payload. The members we need total
**~165 MB** — roughly 0.3% of the archive. This is already proven working: the metadata and joint
PC files were extracted during exploration.

### 3. A shared genotype PC space spans all four cohorts

`MAGE.v1.0.data/dataset_comparison/MAGE_Geuvadis_GTEx_AFGR.top_PCs.txt.gz` (200 KB) holds **2,373
samples in one jointly-generated PC space**, with 20 PCs and membership flags:

| Dataset | Samples in the joint space |
|---|---|
| GTEx | 837 |
| MAGE | 731 |
| AFGR | 597 |
| Geuvadis | 452 |

The 837 GTEx donors carry no population labels (protected), but they occupy the *same* coordinate
system as labelled 1KGP samples. GTEx ancestry can therefore be positioned against labelled
anchors without dbGaP access.

**93 individuals appear in both MAGE and Geuvadis** — the same people, sequenced by two
independent labs. This is a free concordance test: PCCA coordinates for one person derived from
two separate experiments should agree.

## Data acquisition

### MAGE — selective extraction from the remote zip

Read the Zip64 central directory over range requests, then fetch only these members:

| Member | Size |
|---|---|
| `QTL_results/eQTL_results/expression_quants/TMM.filtered.TSS.MAGE.v1.0.bed.gz` | 107 MB |
| `QTL_results/eQTL_results/expression_quants/eQTL_covariates.tab.gz` | small |
| `sample_library_info/sample.metadata.MAGE.v1.0.txt` | 57 KB |
| `dataset_comparison/MAGE_Geuvadis_GTEx_AFGR.top_PCs.txt.gz` | 200 KB |

The expression matrix is BED4+: three TSS position columns, then Ensembl gene ID, then one column
per sample. **Gene coordinates are included**, so no external annotation is needed for MAGE.

MAGE's covariates are supplied separately from the expression matrix (sex, 5 genotype PCs, 60
PEER factors), so the released matrices are **not** covariate-corrected. The GEUVADIS trap —
where the distributed file was PEER-residualized and silently destroyed the signal — cannot recur
here.

### GTEx — per-tissue TPM plus covariates

| Object | Size |
|---|---|
| `bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_{tissue}.gct.gz` × 4 | 323 MB |
| `bulk-qtl/v10/single-tissue-cis-qtl/GTEx_Analysis_v10_eQTL_covariates.tar` | 20.8 MB |
| `annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt` | 36 MB |
| `annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt` | small |

**Per-tissue raw TPM is used in preference to `GTEx_Analysis_v10_eQTL_expression_matrices.tar`**
(1.4 GB). The eQTL matrices are TMM-normalized then inverse-normal transformed; the per-tissue
files are raw TPM. Raw is the right choice for three reasons: it matches the GEUVADIS treatment
(raw RPKM plus our own normalization), so all three cohorts traverse an identical code path and
cross-dataset comparisons are not confounded by differing upstream pipelines; the inverse-normal
transform equalizes every gene's distribution before we see it, duplicating a step our pipeline
already performs on its own terms; and it is a quarter of the download.

GCT format, verified: line 0 `#1.2`, line 1 `<genes>\t<samples>`, line 2 the header. Whole Blood
is 59,033 genes × 803 samples. Sample IDs are `GTEX-1117F-0005-SM-HL9SH`; the **donor** is the
first two dash-separated fields, `GTEX-1117F`.

## Cohorts

### MAGE: 729 individuals

The MAGE README records an unresolved sample swap, identified 2026-03-27 and **not yet fixed in
the Zenodo release**:

- `SRR19762530`, labelled HG00237, is actually NA11919 — not used in QTL mapping
- `SRR19762653`, labelled NA11919, is actually HG00237 — used in QTL mapping

Checking the released metadata: HG00237 appears once via `SRR19762247`, and NA11919 appears once
via `SRR19762653`. So the column labelled NA11919 carries HG00237's material, HG00237 is
effectively represented twice, and true NA11919 is absent from the released matrices.

**Both individuals are excluded**, leaving **729**. Relabelling would duplicate HG00237; dropping
costs 0.3% of the cohort and cannot propagate a mislabel. Revisit if MAGE publishes corrected
data.

All 731 MAGE individuals already carry TBOOO EIDs (verified against `eid_map_1kg.tsv`, 3,202
rows), so no new identifier mapping is required.

### GTEx: 4 tissues

**Whole Blood, Muscle Skeletal, Lung, Thyroid** — chosen for donor count and to test whether the
result holds across tissues rather than being an artefact of one.

Expression PCA across *all* tissues at once would be dominated by tissue identity, which says
nothing about ancestry; population structure is therefore only examined **within** a tissue.

## Expression preprocessing

Reuses `pcca.expression` unchanged. Only the confounder design differs:

| Cohort | Confounders regressed out |
|---|---|
| GEUVADIS | 7 laboratory indicators + sex |
| MAGE | batch (17 levels) + sex |
| GTEx | collection centre (`SMCENTER`) + nucleic-acid batch (`SMNABTCH`) + sex |

Pipeline per cohort: select cohort columns → filter genes (mean expression ≥ 0.1, autosomal,
non-MHC, protein-coding) → quantile-normalize across samples → standardize each gene → regress
out confounders.

### Gene coordinates for GTEx

GTEx per-tissue TPM files carry only `Name` (versioned Ensembl ID) and `Description` (symbol) —
**no chromosome or coordinate**, which `filter_genes` requires for its autosomal and MHC filters.

`pcca.expression.protein_coding_gene_ids` already parses `gencode.v47.annotation.gtf.gz` for the
protein-coding whitelist. That parse is extended to return a coordinate table — gene ID,
chromosome, TSS — reusing the same single pass and cache. No new dependency, no second GTF read.

MAGE needs none of this: its BED4+ matrix carries TSS coordinates directly.

## Genotype PCs

Published 5 PCs for both cohorts — the only option for GTEx, and consistency argues for it in
MAGE.

For MAGE only, genotype PCs are **additionally recomputed** from the 1KGP Phase 3 callset with
the existing plink2 code and correlated against the published values. MAGE is the one cohort
where independent verification is possible, and it costs ~10 minutes. Disagreement would indicate
a sample-alignment error, which is the failure mode that produces plausible-but-wrong CCA output.

Note the build difference: MAGE's published PCs derive from NYGC high-coverage calls on GRCh38,
while our recomputation uses Phase 3 on GRCh37. The two need not match exactly; a high
correlation on the leading PCs is the expected outcome, not identity.

## Analyses

### Per cohort

PCA, then PCCA at **30 expression PCs × 5 genotype PCs → 2 canonical components**, matching the
paper for comparability with the GEUVADIS arm, plus the leave-one-out projection.

For each cohort, the reported quantities are the canonical correlations and the variance in each
leading coordinate explained by population label — the same three-way contrast that made the
GEUVADIS result legible:

| | Expression PCA | Genotype PCA | PCCA projection |
|---|---|---|---|
| GEUVADIS (measured) | 0.048 / 0.071 | 1.000 / 0.975 | 0.897 / 0.605 |

GTEx has no public population labels, so its equivalent uses the inferred continental group from
the joint PC space, reported as inferred.

### Cross-dataset

1. **Joint PC space.** Place all four cohorts in the shared genotype PC coordinates and show
   where unlabelled GTEx donors fall relative to labelled 1KGP anchors. Assign each GTEx donor a
   continental group by nearest labelled neighbours, explicitly marked as inferred.
2. **MAGE ∩ GEUVADIS concordance.** For the 93 shared individuals, compare PCCA coordinates
   derived independently from each cohort. Agreement is evidence the embedding reflects the
   person rather than the experiment.

## TBOOO pipeline integration

New stages following the existing `geuvadis` pattern, including `table_ok` idempotency guards:

```
tbooo download mage     tbooo map mage
tbooo download gtex     tbooo map gtex
```

### participant.parquet

**MAGE** adds `mage_pc1..20` columns to existing 1KGP rows — all 729 individuals already have
EIDs, so this is a column addition with no new rows.

**GTEx** donors are not 1KGP samples and need a new synthetic EID range. The representative
tissue is whole blood, so the rows added are its **803** donors rather than the 943-donor union
across all 50 tissues. TBOOO already partitions
by range (1KGP at 1,000,000+, SGDP at 2,000,000+); GTEx takes **3,000,000+**, written to
`data/metadata/eid_map_gtex.tsv`. The 803 whole-blood donors append as new rows, taking the table from
**3,467 to 4,270**.

GTEx rows carry `gtex_pc*` expression PCs and **null genetic-availability fields** (`p22418`,
`p22828`, `p24051`), because GTEx genotypes are not public. This is faithful rather than
awkward: in UK Biobank different assays cover different participant subsets, and a participant
with expression but no genotype is a real pattern.

`p21000` (ethnic background) for GTEx donors is populated from the **inferred** continental group
derived in the cross-dataset analysis, and the inference is recorded in the docs so no downstream
consumer mistakes it for a reported value.

## Expected outcomes

**MAGE may show partial population structure in plain PCA**, unlike GEUVADIS. The MAGE paper
itself reports continental group explaining a measurable share of expression variance. If plain
PCA partly separates AFR from non-AFR, that is a finding rather than a failure: it would indicate
the GEUVADIS negative result is partly a consequence of that cohort's narrow ancestry range, and
that "expression PCA cannot see population structure" is contingent on cohort composition rather
than universal. PCCA should still sharpen the separation substantially.

**GTEx within a tissue is a two-group contrast** with ~119 African-ancestry donors against ~824
European. Expect a weaker but directionally clear result, and expect consistency across the four
tissues to be the more informative signal.

## Risks

- **Sample/donor alignment.** Every cohort here joins expression to genotype PCs through a
  different identifier: MAGE via `sample_kgpID`, GTEx via a donor prefix extracted from a
  sample ID. Misalignment produces plausible-looking but meaningless CCA output rather than an
  error. Mitigation: explicit intersection helpers that reset the index and raise on missing
  members, mirroring `restrict_to_panel` from the GEUVADIS arm, with tests.
- **GTEx tissue sample-to-donor multiplicity.** Whole Blood shows 803 sample columns; donors must
  be deduplicated before analysis. Mitigation: collapse to one sample per donor and assert
  uniqueness.
- **MAGE recomputed PCs disagreeing with published.** Expected to correlate highly but not
  match exactly given the GRCh37/GRCh38 and callset differences. A low correlation is a
  red flag for alignment, not a reason to abandon the published PCs.
- **The eQTL-arm caveat carries over.** The GEUVADIS work found the paper's eQTL covariate claim
  did not reproduce. No cis-eQTL scan is in scope here; this work is structural only.

## Deliverables

1. `src/tbooo/download/{mage,gtex}.py` and `src/tbooo/pipeline/{mage,gtex}.py`
2. `playground/pcca/{mage,gtex,joint}.py` plus tests
3. `playground/gtex_mage_pcca.ipynb`
4. Figures in `data/metadata/eda/` — per-cohort PCA/PCCA/LOO, joint PC space, concordance
5. `data/metadata/eid_map_gtex.tsv`, `geuvadis_pcca_coords.tsv` equivalents per cohort
6. Extension of `docs/4_pcca_methods.md` with the multi-cohort results

## Out of scope

cis-eQTL scans for either cohort; joint *expression* embedding across cohorts (quantification and
tissue differences would make batch effects dominate); GTEx tissues beyond the four selected;
AFGR, which appears in the joint PC space but has no expression data ingested here.
