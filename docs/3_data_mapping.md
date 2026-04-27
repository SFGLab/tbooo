# Data Mapping: 1000 Genomes & Simons → UK Biobank Structure

This document specifies exactly how data from the 1000 Genomes Project (1KGP) and Simons Genome Diversity Project (SGDP) is transformed and reorganized to mirror the UK Biobank (UKB) data structure on DNAnexus RAP. It is the authoritative reference for TBOOO's data pipeline decisions.

Cross-references: [1_ukb_structure.md](1_ukb_structure.md) for UKB format details, [2_data_sources.md](2_data_sources.md) for source data details.

---

## 1. Overall Mapping Strategy

TBOOO produces a directory tree and file set that structurally mirrors what a UKB-approved researcher would find on DNAnexus RAP, using public data as a stand-in for the restricted UKB cohort.

### Source → UKB data type assignments

| UKB Data Type | UKB Field | Source Dataset | Source Data |
|---------------|-----------|---------------|-------------|
| Genotyping array (PLINK) | 22418 | 1KGP Phase 3 | Phase 3 VCFs filtered to UKB array SNP positions → PLINK |
| Imputed genotypes (BGEN) | 22828 | 1KGP Phase 3 | Phase 3 VCFs → BGEN conversion (already population-level calls) |
| WGS individual CRAMs | 23149 | 1KGP NYGC 30x + SGDP | Individual CRAM files, renamed to UKB conventions |
| WGS cohort pVCF | 23370–23384 | 1KGP NYGC 30x | Per-chromosome multi-sample VCFs reformatted as pVCF blocks |
| WES cohort PLINK/BGEN | 23157 | 1KGP NYGC 30x | WGS VCFs intersected with UKB exome capture BED → PLINK/BGEN |
| Phenotypic data (Parquet) | various | 1KGP + SGDP metadata | Sample panel files → synthetic UKB field columns |
| Sample QC file | — | 1KGP + SGDP | Computed from KING + PCA on merged dataset |
| Relatedness file | — | 1KGP pedigree | `.ped` family structure + KING kinship estimates |

### Key design decisions

- **Array + imputation** use GRCh37 coordinates, matching 1KGP Phase 3 natively (no liftover needed).
- **WGS + WES** use GRCh38 coordinates, matching 1KGP NYGC 30x and SGDP natively.
- **SGDP samples** supplement 1KGP in the WGS layer only, adding geographic diversity beyond 1KGP's 26 populations.
- All samples are assigned synthetic 7-digit EIDs in the range 1,000,000–6,999,999, mirroring UKB pseudonymization.
- Only metadata-derivable phenotype fields are simulated; clinical phenotypes (disease diagnoses, hospital records) are left empty or omitted.

---

## 2. Participant ID (EID) Assignment

UKB pseudonymizes all participants with randomized 7-digit EIDs unique to each access application. TBOOO generates a deterministic synthetic EID for each source sample so that files and phenotype rows are consistently linked.

### EID generation scheme

```
EID = 1_000_000 + sequential_index
```

- 1KGP samples are assigned EIDs 1,000,000–1,003,201 (one per NYGC 3,202-sample set, ordered by the sample panel file)
- SGDP samples are assigned EIDs 2,000,000–2,000,278 (one per public sample, ordered by ENA accession)
- Index is zero-based within each source dataset

The mapping is stored in two TSV files:

```
data/metadata/eid_map_1kg.tsv     # columns: eid, sample_id, pop, super_pop, sex, source
data/metadata/eid_map_sgdp.tsv    # columns: eid, ena_accession, population, region, sex, source
```

### EID prefix subfolders

Individual-level bulk files are organized by the first two digits of the EID, exactly as in UKB:

| EID range | Subfolder |
|-----------|-----------|
| 1,000,000 – 1,099,999 | `10/` |
| 1,100,000 – 1,199,999 | `11/` |
| … | … |
| 2,000,000 – 2,099,999 | `20/` |

---

## 3. Output Directory Structure

The TBOOO output mirrors the `/Bulk/` tree and adds a `Showcase/` directory for phenotype data:

```
data/
├── Bulk/
│   ├── Genotype Results/
│   │   └── Genotype calls/
│   │       ├── ukb22418_c1_b0_v2.bed
│   │       ├── ukb22418_c1_b0_v2.bim
│   │       ├── ukb22418_c1_b0_v2.fam
│   │       ├── ukb22418_c2_b0_v2.bed
│   │       │   ...
│   │       └── ukb22418_c22_b0_v2.fam
│   ├── Imputed/
│   │   ├── ukb22828_c1_b0_v3.bgen
│   │   ├── ukb22828_c1_b0_v3.bgen.bgi
│   │   ├── ukb22828_c1_b0_v3.sample
│   │   │   ...
│   │   └── ukb22828_c22_b0_v3.sample
│   ├── Exome sequences/
│   │   └── Population level exome OQFE variants, PLINK format - 500k release/
│   │       ├── ukb23157_c1_b0_v1.bed
│   │       ├── ukb23157_c1_b0_v1.bim
│   │       ├── ukb23157_c1_b0_v1.fam
│   │       │   ...
│   │       └── ukb23157_c22_b0_v1.fam
│   └── Whole genome sequences/
│       ├── 10/
│       │   ├── 1000001_23149_0_0.cram
│       │   ├── 1000001_23149_0_0.cram.crai
│       │   ├── 1000002_23149_0_0.cram
│       │   │   ...
│       └── 20/
│           ├── 2000001_23149_0_0.cram
│           │   ...
├── Showcase/
│   └── participant.parquet        # synthetic phenotype table
└── metadata/
    ├── eid_map_1kg.tsv
    ├── eid_map_sgdp.tsv
    ├── ukb_sqc_v2.txt             # synthetic sample QC file
    └── ukb_rel.txt                # synthetic relatedness file
```

---

## 4. Genotyping Array Data (UKB Field 22418)

### Strategy

UKB array data covers 805,426 SNPs typed on Affymetrix arrays. The public data equivalent is constructed by subsetting the 1KGP Phase 3 VCFs to the positions present on the UKB Biobank Axiom array manifest, then converting to PLINK format.

### Source

| Attribute | Value |
|-----------|-------|
| Source dataset | 1KGP Phase 3 |
| Source files | `ALL.chr{1-22}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` |
| Source reference | GRCh37 — matches UKB array coordinates natively |
| Samples | All 2,504 unrelated Phase 3 individuals |

### Position filtering

The UKB Biobank Axiom array manifest (available from Thermo Fisher) provides GRCh37 positions for all 825,927 probes (805,426 after QC). Only variants present in both the 1KGP VCF and the array manifest are retained.

Array manifest file (used as filter):
```
data/reference/ukb_array_manifest_GRCh37.txt   # columns: chrom, pos, rsid, allele_A, allele_B
```

Filtering command (bcftools):
```bash
bcftools view \
  --regions-file data/reference/ukb_array_manifest_GRCh37.txt \
  --output-type z \
  ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
  | bcftools norm --multiallelics - \
  > tmp_chr1_array_sites.vcf.gz
```

### Format conversion: VCF → PLINK

```bash
plink2 \
  --vcf tmp_chr1_array_sites.vcf.gz \
  --make-bed \
  --out data/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c1_b0_v2
```

### Output files

| File | Contents |
|------|----------|
| `ukb22418_c<chr>_b0_v2.bed` | Binary genotype matrix; all 2,504 samples × array-overlapping variants on that chromosome |
| `ukb22418_c<chr>_b0_v2.bim` | SNP manifest; columns: chr, rsID (from 1KGP annotation), cM (from 1KGP genetic map), bp (GRCh37), A1, A2 |
| `ukb22418_c<chr>_b0_v2.fam` | Sample manifest; column 6 encodes a synthetic batch number (see below) |

### FAM file construction

The UKB FAM column 6 encodes batch. For TBOOO, batches are assigned by population group to preserve the population-stratified batch structure seen in real UKB data:

| Batch code | Samples included | Approximate size |
|-----------|-----------------|-----------------|
| 1 | EUR samples (CEU, TSI, FIN, GBR, IBS) | ~503 |
| 2 | AFR samples | ~661 |
| 3 | EAS samples | ~504 |
| 4 | SAS samples | ~489 |
| 5 | AMR samples | ~347 |

FAM column layout:
```
FID    IID       PAT  MAT  SEX  BATCH
0      1000001   0    0    1    1
0      1000002   0    0    2    1
...
```

- `FID`: always 0 (no family structure in Phase 3 unrelated set)
- `IID`: synthetic EID
- `PAT`/`MAT`: 0 (unknown / unrelated)
- `SEX`: from 1KGP panel file (`gender` column; 1=male, 2=female)
- `BATCH`: population-group batch code above

### BIM file notes

- cM positions populated from the 1KGP pedigree-based recombination maps (`joepickrell/1000-genomes-genetic-maps`); set to 0 if not available for a position
- rsIDs taken from the 1KGP VCF `ID` field; positions without rsIDs use the format `chr:pos:A1:A2`
- Alleles are harmonized to the array manifest strand orientation

---

## 5. Imputed Genotype Data (UKB Field 22828)

### Strategy

UKB imputed data contains ~92.7M variants across 487,409 participants in BGEN v1.2 format. The simulation uses 1KGP Phase 3 genotypes directly as a proxy for imputed calls: since 1KGP Phase 3 was itself one of the imputation reference panels used for UKB, its variant density (81.4M SNVs) is broadly representative of imputed data. Hard genotype calls from the VCF are encoded as genotype probabilities (0/1 rather than fractional).

For a more realistic simulation of dosage uncertainty, an optional imputation step can be performed: subset the 2,504 Phase 3 samples to array-overlapping sites only, then re-impute using Beagle or Minimac against the full Phase 3 reference panel.

### Source

| Attribute | Value |
|-----------|-------|
| Source dataset | 1KGP Phase 3 |
| Source files | `ALL.chr{1-22}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` |
| Source reference | GRCh37 — matches UKB imputed genotype coordinates natively |
| Samples | All 2,504 unrelated Phase 3 individuals |

### Format conversion: VCF → BGEN

```bash
# Decompose multi-allelic sites, then convert to BGEN v1.2
bcftools norm --multiallelics - --output-type z \
  ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
  | qctool \
    -g /dev/stdin -filetype vcf \
    -og data/Bulk/Imputed/ukb22828_c1_b0_v3.bgen \
    -ofiletype bgen_v1.2 \
    -bgen-bits 8
```

Index generation:
```bash
bgenix -g data/Bulk/Imputed/ukb22828_c1_b0_v3.bgen -index
# produces ukb22828_c1_b0_v3.bgen.bgi
```

### Sample file format

The `.sample` file accompanies each BGEN file and lists samples in the same order as the BGEN:

```
ID_1 ID_2 missing sex
0    0    0       D
1000001 1000001 0 1
1000002 1000002 0 2
...
```

- `ID_1`: synthetic EID (repeated in `ID_2` — UKB uses the same value in both columns)
- `missing`: 0 for all samples (no missingness indicator)
- `sex`: 1=male, 2=female, from 1KGP panel file

### Output files

| File | Naming pattern |
|------|---------------|
| BGEN | `ukb22828_c<chr>_b0_v3.bgen` |
| BGEN index | `ukb22828_c<chr>_b0_v3.bgen.bgi` |
| Sample file | `ukb22828_c<chr>_b0_v3.sample` (one shared sample file per chromosome) |

### Variant filtering

Apply MAF-dependent info-score filtering to match UKB practice. For the direct-conversion approach (hard calls → BGEN), info score = 1.0 for all variants (no uncertainty), so all variants pass any info threshold. If the optional imputation step is used, apply the UKB thresholds:

| MAF Range | Info Score Threshold |
|-----------|---------------------|
| > 3% | > 0.3 |
| 1% – 3% | > 0.6 |
| 0.5% – 1% | > 0.8 |
| 0.1% – 0.5% | > 0.9 |

MAF is computed from the 2,504-sample genotype matrix using PLINK2 `--freq`.

---

## 6. Whole Genome Sequencing — Individual CRAMs (UKB Field 23149)

### Strategy

UKB provides individual CRAM files for each participant at ~30x depth. The public equivalents are the existing per-sample CRAM files from 1KGP NYGC 30x and SGDP. They are renamed to the UKB individual-level file naming convention and placed in EID-prefix subfolders.

### Sources

| Source | Samples | Reference | Coverage |
|--------|---------|-----------|---------|
| 1KGP NYGC 30x | 3,202 | GRCh38 | 30x |
| SGDP Section A | 279 | GRCh38DH | ≥ 30x |

### File renaming

UKB individual CRAM naming:
```
<EID>_<FIELD-ID>_<INSTANCE-ID>_<ARRAY-ID>.cram
```

For the WGS CRAM field (23149), instance 0, array 0:
```
<EID>_23149_0_0.cram
<EID>_23149_0_0.cram.crai
```

**1KGP NYGC 30x rename example:**
```
# Source:
NA12878.alt_bwamem_GRCh38DH.20150826.CEU.simons.cram   (hypothetical NYGC naming)

# Target (EID 1000001 assigned to NA12878):
data/Bulk/Whole genome sequences/10/1000001_23149_0_0.cram
data/Bulk/Whole genome sequences/10/1000001_23149_0_0.cram.crai
```

**SGDP rename example:**
```
# Source:
SAMEA3302732.alt_bwamem_GRCh38DH.20200922.Greek.simons.cram

# Target (EID 2000001 assigned to SAMEA3302732):
data/Bulk/Whole genome sequences/20/2000001_23149_0_0.cram
data/Bulk/Whole genome sequences/20/2000001_23149_0_0.cram.crai
```

### Reference genome note

UKB WGS CRAMs use GRCh38. Both 1KGP NYGC 30x and SGDP CRAMs are also aligned to GRCh38 (or GRCh38DH, which is compatible). No re-alignment is needed.

---

## 7. Whole Genome Sequencing — Cohort pVCF (UKB Fields 23370–23384)

### Strategy

UKB WGS cohort data is released as multi-sample pVCF files split into 151,561 genomic blocks. TBOOO produces per-chromosome multi-sample VCFs from 1KGP NYGC 30x data (which are already per-chromosome phased VCFs) as a coarser-grained equivalent.

### Source

| Attribute | Value |
|-----------|-------|
| Source | 1KGP NYGC 30x, per-chromosome phased VCFs |
| Source files | `CCDG_14151_B01_GRM_WGS_2020-08-05_chr{N}.filtered.shapeit2-duohmm-phased.vcf.gz` |
| Reference | GRCh38 |
| Samples | 3,202 individuals |

### Sample ID replacement

The VCF sample headers use original 1KGP sample IDs (e.g., `NA12878`). These must be replaced with synthetic EIDs to match UKB conventions:

```bash
# Generate a sample rename file: original_id → EID
bcftools reheader \
  --samples data/metadata/vcf_sample_rename_1kg.txt \
  CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz \
  -o data/Bulk/Whole\ genome\ sequences/ukb23370_c1_b0_v1.pvcf.gz
```

The rename file format (one line per sample):
```
NA12878 1000001
NA12891 1000002
...
```

### Output naming

UKB uses field 23370–23384 for the 500k WGS pVCF release. TBOOO maps to a single field:

```
ukb23370_c<chr>_b0_v1.pvcf.gz
ukb23370_c<chr>_b0_v1.pvcf.gz.tbi
```

One file per chromosome (chr1–22, chrX). Block-level splitting (matching UKB's 151,561 blocks) is optional and only needed if downstream tools expect it.

---

## 8. Whole Exome Sequencing Simulation (UKB Field 23157)

### Strategy

Neither 1KGP nor SGDP has true exome-capture sequencing — both are WGS. To simulate UKB WES data, WGS variant calls are intersected with the UKB exome capture target BED file (IDT xGen Exome Research Panel v1.0, ~39 Mbp). This produces a variant set restricted to exonic regions, which is functionally equivalent for downstream analysis.

### Source

| Attribute | Value |
|-----------|-------|
| Source | 1KGP NYGC 30x per-chromosome VCFs |
| Reference | GRCh38 (matches UKB WES reference) |
| Capture BED | IDT xGen Exome Research Panel v1.0, GRCh38 coordinates |

### BED file

```
data/reference/idt_xgen_exome_v1_GRCh38.bed
```

Columns: `chrom`, `start` (0-based), `end`, `target_name`

Source: available from IDT website; also distributed with the UKB WES helper files on RAP.

### Variant intersection

```bash
bcftools view \
  --regions-file data/reference/idt_xgen_exome_v1_GRCh38.bed \
  CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz \
  | bcftools norm --multiallelics - \
  | bcftools reheader --samples data/metadata/vcf_sample_rename_1kg.txt \
  -o tmp_chr1_exome.vcf.gz
```

### PLINK conversion

```bash
plink2 \
  --vcf tmp_chr1_exome.vcf.gz \
  --make-bed \
  --out data/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ PLINK\ format\ -\ 500k\ release/ukb23157_c1_b0_v1
```

### BGEN conversion

```bash
qctool \
  -g tmp_chr1_exome.vcf.gz -filetype vcf \
  -og data/Bulk/Exome\ sequences/ukb23157_c1_b0_v1.bgen \
  -ofiletype bgen_v1.2 \
  -bgen-bits 8
bgenix -g data/Bulk/Exome\ sequences/ukb23157_c1_b0_v1.bgen -index
```

### Limitation

The intersection approach produces genuine WGS-quality variants in exonic regions, not array-capture variants. SNV density, read depth distribution, and strand bias patterns will differ from true exome capture. This is adequate for testing analysis pipelines but not for benchmarking exome-specific QC.

---

## 9. Phenotypic Data Mapping

### Strategy

UKB phenotypic data is a rich clinical and self-reported dataset that cannot be reproduced from 1KGP/SGDP. TBOOO generates a synthetic Parquet table containing only the fields derivable from sample metadata, plus placeholder columns (null/NA) for commonly-used clinical fields.

### Output

```
data/Showcase/participant.parquet
```

Schema follows the UKB column naming convention: `p<FIELD-ID>_i<INSTANCE-ID>_a<ARRAY-ID>`.

### Field mapping table

The following UKB fields can be directly or approximately populated:

| UKB Column | Field ID | Description | Source in 1KGP/SGDP | Notes |
|-----------|----------|-------------|---------------------|-------|
| `eid` | — | Participant identifier | Synthetic EID | Primary key |
| `p31` | 31 | Sex | `gender` column in panel file | 1=male, 2=female |
| `p21000_i0` | 21000 | Ethnic background | Superpopulation code (see §9.1) | Instance 0; approximate |
| `p22006` | 22006 | Genetic ethnic grouping — White British | 1 if super_pop=EUR and pop=GBR | Binary flag |
| `p22020` | 22020 | Used in PCA calculation | 1 for all unrelated samples | Used for QC flagging |
| `p22000` | 22000 | Genotyping array batch | Synthetic batch code (§4) | See FAM file batch scheme |
| `p22418` | 22418 | Genotype calls available | 1 for all samples | Binary flag |
| `p22828` | 22828 | Imputed genotypes available | 1 for all samples | Binary flag |
| `p23149` | 23149 | WGS CRAM available | 1 for all samples with CRAM | Binary flag |
| `p54_i0` | 54 | UK Biobank assessment centre | Population code (see §9.2) | Approximate; instance 0 |

### 9.1 Population → Ethnic Background (Field 21000)

UKB Field 21000 uses Data-Coding 1001 for self-reported ethnic background. The mapping from 1KGP/SGDP superpopulation to UKB ethnic background category is approximate (these are genetic ancestry groups, not self-reported ethnicity):

| Source | Superpop / Region | UKB Code | UKB Label |
|--------|------------------|----------|-----------|
| 1KGP | EUR | 1 | White |
| 1KGP | AFR | 4 | Black or Black British |
| 1KGP | EAS | 3 | Asian or Asian British |
| 1KGP | SAS | 3 | Asian or Asian British |
| 1KGP | AMR | 2 | Mixed |
| SGDP | West Eurasia | 1 | White |
| SGDP | Africa | 4 | Black or Black British |
| SGDP | East Asia | 3 | Asian or Asian British |
| SGDP | South Asia | 3 | Asian or Asian British |
| SGDP | Central Asia / Siberia | 6 | Other ethnic group |
| SGDP | Oceania | 6 | Other ethnic group |
| SGDP | Native Americas | 2 | Mixed |

> Note: This mapping is a structural approximation. Real UKB ethnic background is self-reported and does not directly correspond to genetic superpopulation labels.

### 9.2 Population → Assessment Centre (Field 54)

UKB Field 54 records which UK Biobank assessment centre the participant visited. There is no meaningful equivalent in 1KGP/SGDP. TBOOO assigns a synthetic centre code based on population:

| Population group | Synthetic centre code | Rationale |
|-----------------|----------------------|-----------|
| EUR / GBR | 11010 (Leeds) | UK-derived population |
| EUR / other | 11020 (Nottingham) | Arbitrary EUR centre |
| AFR | 11021 (Bristol) | Arbitrary non-EUR centre |
| EAS | 11022 (Hounslow) | Arbitrary non-EUR centre |
| SAS | 11023 (Croydon) | Arbitrary non-EUR centre |
| AMR / Other | 11024 (Birmingham) | Arbitrary non-EUR centre |

### 9.3 Fields left as null

The following commonly-used UKB fields have no equivalent in 1KGP/SGDP metadata and are included as null columns in the Parquet schema so downstream scripts referencing them do not break:

| Field | Description |
|-------|-------------|
| `p21022` | Age at recruitment |
| `p53_i0` | Date of attending assessment centre |
| `p21001_i0` | Body mass index |
| `p4079_i0` | Diastolic blood pressure |
| `p4080_i0` | Systolic blood pressure |
| `p41270` | Diagnoses – ICD10 |
| `p20002_i0_a0` | Non-cancer illness codes |
| `p20003_i0_a0` | Treatment/medication codes |
| `p40001_i0` | Cause of death – ICD10 |

---

## 10. Sample QC File (ukb_sqc_v2.txt)

### Format

The UKB sample QC file has one row per participant. The TBOOO equivalent is generated from computed statistics on the merged 1KGP + SGDP genotype data.

Output: `data/metadata/ukb_sqc_v2.txt`

### Columns populated

| Column | Source | Method |
|--------|--------|--------|
| `FID` | — | 0 for all |
| `IID` | EID map | Synthetic EID |
| `used.in.pca.calculation` | — | 1 for all unrelated Phase 3 samples; 0 for related (NYGC-added) and SGDP samples |
| `in.white.british.ancestry.subset` | Population | 1 if super_pop=EUR and pop=GBR |
| `excess.relatives` | 1KGP pedigree | 1 if sample has >10 identified relatives in the dataset (based on `.ped` family structure) |
| `putative.sex.chromosome.aneuploidy` | — | 0 for all (not computed from public data) |
| `het.missing.outliers` | PLINK | 1 if sample is >5 SD from mean heterozygosity rate on autosomes |

### Computation commands

Heterozygosity outliers (run on merged array-site PLINK files):
```bash
plink2 \
  --bfile data/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_merged_all_chrs \
  --het \
  --out data/metadata/het_stats
# Then flag samples >5 SD from mean as outliers
```

---

## 11. Relatedness File (ukb_rel.txt)

### Format

UKB provides pairwise kinship estimates for all pairs above a kinship threshold of 0.0442 (3rd-degree relatives and closer). TBOOO generates this using KING software on the merged array-site PLINK data.

Output: `data/metadata/ukb_rel.txt`

5-column tab-separated file:
```
EID1     EID2     HetHet   IBS0    Kinship
1000001  1000003  0.2310   0.0021  0.2498
...
```

### Source of known relatives

The 1KGP NYGC `.ped` file (updated 2025) documents 602 complete parent–child trios and other family relationships. These known first-degree pairs should have kinship ~0.25 in the output.

SGDP samples are unrelated to each other and to 1KGP samples by design; no SGDP pairs are expected above the kinship threshold.

### Computation

```bash
# Step 1: compute kinship on array-site PLINK data
king \
  --bfile data/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_merged_all_chrs \
  --kinship \
  --prefix data/metadata/king_output

# Step 2: filter to pairs above 3rd-degree threshold (kinship > 0.0442)
awk '$7 > 0.0442' data/metadata/king_output.kin0 > data/metadata/ukb_rel_king.txt

# Step 3: replace original sample IDs with EIDs
python src/data/replace_ids_in_rel.py \
  --rel data/metadata/ukb_rel_king.txt \
  --map data/metadata/eid_map_1kg.tsv \
  --out data/metadata/ukb_rel.txt
```

---

## 12. Reference Genome Handling

All TBOOO data is organized into two coordinate systems, matching UKB:

| Data layer | TBOOO reference | UKB reference |
|-----------|----------------|--------------|
| Array PLINK | GRCh37 | GRCh37 |
| Imputed BGEN | GRCh37 | GRCh37 |
| WES PLINK/BGEN | GRCh38 | GRCh38 |
| WGS pVCF | GRCh38 | GRCh38 |
| WGS individual CRAMs | GRCh38 | GRCh38 |

No coordinate liftover is required because:
- 1KGP Phase 3 VCFs are on GRCh37 → used for array and imputation layers
- 1KGP NYGC 30x VCFs and SGDP CRAMs are on GRCh38 → used for WGS and WES layers

Reference FASTA files needed:
```
data/reference/GRCh37/hs37d5.fa          # for decoding GRCh37 CRAMs
data/reference/GRCh38/GRCh38DH.fa        # for decoding GRCh38 CRAMs
```

---

## 13. Tool Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| bcftools | ≥ 1.17 | VCF filtering, normalization, reheadering |
| PLINK 2 | ≥ 2.00a3 | VCF → PLINK conversion, allele frequency, heterozygosity |
| qctool | ≥ 2.0.8 | VCF → BGEN v1.2 conversion |
| bgenix | ≥ 1.1.7 | BGEN indexing (.bgi) |
| KING | ≥ 2.3 | Kinship estimation for relatedness file |
| samtools | ≥ 1.17 | CRAM reheadering, index generation |
| Python | ≥ 3.10 | Metadata processing, ID remapping scripts |
| pandas / pyarrow | — | Parquet phenotype table generation |

---

## 14. What Cannot Be Simulated

The following UKB data types have no equivalent in 1KGP or SGDP and are out of scope for TBOOO:

| UKB Data Type | Reason not simulable |
|---------------|---------------------|
| Self-reported phenotypes (diagnoses, medications, lifestyle) | No equivalent questionnaire data in public genomics datasets |
| Hospital Episode Statistics (hesin*) | UK NHS linkage; no equivalent |
| GP records | UK primary care linkage; no equivalent |
| Death records | UK mortality register; no equivalent |
| Imaging (brain MRI, cardiac MRI, DEXA) | Not collected in 1KGP/SGDP |
| Proteomics (Olink) | Not collected |
| Metabolomics | Not collected |
| Accelerometry / wearable data | Not collected |
| COVID-19 test results | Not collected |
| Repeat assessment instances (i1, i2, i3) | All 1KGP/SGDP samples have only one timepoint |

---

## 15. Limitations and Caveats

| Limitation | Impact |
|-----------|--------|
| Array SNP overlap is partial | Not all 805,426 UKB array SNPs exist in 1KGP Phase 3 VCFs; only overlapping sites are included |
| No true imputation uncertainty | BGEN files from direct VCF conversion have dosage probabilities of 0 or 1; real UKB imputed data has fractional probabilities |
| Sample size | TBOOO has at most 3,481 samples (3,202 + 279) vs. UKB's 487,000–500,000; statistical power for rare-variant tests will be much lower |
| Population composition | UKB is ~94% European; TBOOO is globally diverse; allele frequency distributions will differ substantially |
| No WES capture simulation artifacts | Intersection of WGS variants with exome BED does not reproduce capture efficiency gradients, off-target reads, or GC-bias artifacts |
| Phenotype depth | Only ~10 phenotype fields are populated; all clinical fields are null |
| No longitudinal data | All samples have a single timepoint; no repeat assessments |
