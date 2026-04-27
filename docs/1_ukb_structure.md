# UK Biobank Data Structure on DNAnexus RAP

The UK Biobank (UKB) makes its data available to approved researchers through the **Research Analysis Platform (RAP)**, a cloud computing environment hosted on DNAnexus. This document describes how data is organized, named, and accessed on the RAP — the structure that TBOOO aims to simulate using public datasets.

---

## 1. Platform Overview

The RAP hosts two fundamentally different categories of data:

| Category | Storage Format | Access Method |
|----------|---------------|---------------|
| **Bulk data** | Files on DNAnexus filesystem | Direct file access via `/Bulk/` folder tree |
| **Phenotypic data** | Parquet tables (Spark SQL) | Cohort Browser GUI, `dx extract_dataset` CLI, JupyterLab Spark SQL |

Every approved access application gets its own **project** on DNAnexus. Within that project, participant identifiers (EIDs) are pseudonymized: each application receives a unique set of randomized 7-digit EIDs (typically in the range 1,000,000–6,000,000). The same biological participant will have a different EID in different applications.

---

## 2. File Naming Conventions

### 2.1 Cohort-Level Bulk Files

Files containing genotype data for the entire cohort (array calls, imputation, WES aggregates, WGS pVCFs) follow this pattern:

```
ukb<FIELD-ID>_c<CHROM>_b<BLOCK>_v<VERSION>.<SUFFIX>
```

| Token | Meaning | Example |
|-------|---------|---------|
| `FIELD-ID` | UK Biobank data field number | `22418` (array calls), `22828` (imputed) |
| `CHROM` | Chromosome | `1`–`22`, `X`, `Y`, `M` |
| `BLOCK` | Genomic block index, 0-based (for datasets split into segments) | `0` |
| `VERSION` | Dataset version assigned by UKB | `2`, `3` |
| `SUFFIX` | File extension | `.bed`, `.bgen`, `.pvcf.gz`, `.vcf.gz` |

Companion index files (`.tbi`, `.bgi`, `.crai`) share the same prefix as their parent file.

**Examples:**
```
ukb22418_c22_b0_v2.bed          # Array calls, chr22, block 0, v2
ukb22418_c22_b0_v2.bim
ukb22418_c22_b0_v2.fam
ukb22828_c1_b0_v3.bgen          # Imputed genotypes, chr1, v3
ukb22828_c1_b0_v3.bgen.bgi      # BGEN index
ukb22828_c1_b0_v3.sample        # Sample file
```

### 2.2 Individual-Level Bulk Files

Files containing data for a single participant (e.g., individual CRAM files) follow:

```
<EID>_<FIELD-ID>_<INSTANCE-ID>_<ARRAY-ID>.<SUFFIX>
```

These files are grouped into subfolders named by the first two digits of the EID (`10/` through `60/`) to avoid directory size limits.

**Example path:**
```
/Bulk/Whole genome sequences/10/1012345_23149_0_0.cram
/Bulk/Whole genome sequences/10/1012345_23149_0_0.cram.crai
```

---

## 3. Genotyping Array Data (Field 22418)

### Location
```
/Bulk/Genotype Results/Genotype calls/
```

### Format
PLINK binary format — three files per chromosome block:

| File | Contents |
|------|----------|
| `.bed` | Binary genotype data (2 bits per genotype, all samples × all SNPs) |
| `.bim` | SNP manifest: chromosome, SNP ID (rsID or Affymetrix ID), cM position, bp position, Allele1, Allele2 |
| `.fam` | Sample manifest: Family ID, Individual ID, Paternal ID, Maternal ID, Sex, Phenotype/Batch |

**Reference genome:** GRCh37 (hg19), 1-based coordinates.

### Coverage
- **805,426 markers** total
- **UK BiLEVE Axiom array** (~50,000 participants)
- **UK Biobank Axiom array** (~450,000 participants)

### Batches
- **106 genotyping batches**, ~4,700 individuals per batch
- Batch information is encoded in column 6 (Phenotype) of the `.fam` file
- **Field 22000** in the phenotype database stores the array type and batch designation per participant (Data-Coding 22000)

### Quality Control
- Markers are flagged as passing/failing QC; SNP lists for QC-passing markers are available (e.g., `final_array_snps_GRCh38_qc_pass.snplist`)
- The `.bim` SNP_ID is the rsID where available; otherwise the Affymetrix SNP ID is used

---

## 4. Imputed Genotype Data (Fields 22828–22829)

### Format
**BGEN v1.2** — the standard format for imputed genotype probabilities in UKB.

| File | Contents |
|------|----------|
| `.bgen` | Binary genotype probabilities, 8-bit encoded, zlib-compressed |
| `.sample` | Sample order and metadata |
| `.bgi` | BGEN index file for fast variant/sample lookup |

**Reference genome:** GRCh37, 1-based coordinates. Forward (+) strand orientation.

### Scale
- **487,409 individuals**
- **92,693,895 variants** across autosomes
- Total storage: ~2.1 TB across all chromosomes

### Versions
- **v2:** Initial release — known bugs; do not use
- **v3:** Corrected release — recommended for all analyses

### Imputation Panels
| Panel | Usage |
|-------|-------|
| Haplotype Reference Consortium (HRC) | Primary panel; largest set of variants |
| UK10K + 1000 Genomes Phase 3 | Fallback for SNPs absent from HRC |
| Genomics England (78,195 individuals) | Newer releases; improves accuracy especially for rarer variants |

### Phased Haplotypes
- Stored in a companion field (format: `ukb_hap_chr<chr>_v2.bgen`)
- Phased using **SHAPEIT3** statistical phasing
- Hard-called probabilities (values are 0 or 1, not fractional)

### MAF-Dependent Info-Score Filtering

Variants are released with MAF-conditional information score thresholds:

| MAF Range | Info Score Threshold |
|-----------|---------------------|
| > 3% | > 0.3 |
| 1% – 3% | > 0.6 |
| 0.5% – 1% | > 0.8 |
| 0.1% – 0.5% | > 0.9 |

SNP lists with info scores and MAF values are available as Resources 1965, 1671, and 1967. QCTOOL v2 was used to compute MAF and info statistics.

### Multi-Allelic Handling
All multi-allelic variants are split into bi-allelic records before storage.

---

## 5. Whole Exome Sequencing (WES) Data

### Field IDs

| Field | Contents |
|-------|----------|
| `23157` | Latest WES release (v11+) |
| `23148` | v7 and later releases |
| `23146` | Previous versions |
| `23141`–`23145` | WES metadata and auxiliary data |

### Location
```
/Bulk/Exome sequences/
/Bulk/Exome sequences_Previous exome releases/
/Bulk/Exome sequences_Alternative exome processing/
```

Example sub-paths:
```
/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - 500k release/
/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release/
/Bulk/Exome sequences/Population level exome OQFE variants, VCF format/
/Bulk/Exome sequences_Alternative exome processing/Exome variant call files (gnomAD) (VCFs)/helper_files/
```

### Scale
- **470,000 participants** (released in waves: 50k → 200k → 300k → 450k → 470k final)
- **Reference genome:** GRCh38

### Sequencing and Processing Pipeline (OQFE Protocol)
1. **Panel:** IDT xGen Exome Research Panel v1.0 + supplemental probes; targets ~39 Mbp of genome
2. **Alignment:** BWA-MEM with alt-aware mapping to GRCh38
3. **Variant calling:** NVIDIA Clara Parabricks DeepVariant (per-sample)
4. **Aggregation:** GLnexus multi-sample joint calling
5. **Output:** pVCF (project VCF), then derived PLINK and BGEN formats

### Available Formats

| Format | Description | Recommended Use |
|--------|-------------|-----------------|
| **pVCF** | Multi-sample VCF, split by chromosome and genomic segments | Direct variant inspection |
| **PLINK** | Bi-allelic; multi-allelic decomposed and normalized | General association testing |
| **BGEN** | Derived from pVCF with variant normalization | GWAS with Regenie |

### Annotation and Helper Files
- `ukb23158_helper_files.pdf` — documentation of annotation labels
- `ukb23158_500k_OQFE.masks` — masks combining annotation labels for burden tests
- SNP lists per chromosome
- Variant annotation: snpEff GRCh38.92, dbSNP build 154, gnomAD r2.1.1

---

## 6. Whole Genome Sequencing (WGS) Data

### Field IDs

| Field Range | Contents |
|-------------|----------|
| `23149`–`23151` | Individual CRAM files and initial WGS VCF data |
| `23370`–`23384` | 500k WGS expanded releases |

### Locations
```
/Bulk/Whole genome sequences/          # Individual CRAM files, organized by EID prefix
/Bulk/DRAGEN WGS/                      # DRAGEN-processed WGS data
/Bulk/DRAGEN WGS/ML-corrected DRAGEN population level WGS variants, pVCF format [500k release]/chr1/
/Bulk/DRAGEN WGS/ML-corrected DRAGEN population level WGS variants, pVCF format [500k release]/chr2/
...
/Bulk/DRAGEN WGS/ML-corrected DRAGEN population level WGS variants, pVCF format [500k release]/chrX/
```

### Scale
- **500,000+ participants**
- **Reference genome:** GRCh38

### Processing Pipelines

| Pipeline | Participants | Tools |
|----------|-------------|-------|
| Pipeline 1 (legacy) | 200k initial release | BWA-MEM + GATK; now merged into 500k enduring release |
| Pipeline 2 (current) | 500k enduring release | Illumina DRAGEN v3.7.8; ML-corrected variant calls |

### Available Formats

| Format | Description | Scale |
|--------|-------------|-------|
| **CRAM** | Individual-level compressed alignments + `.crai` index | Per participant |
| **pVCF** | Cohort multi-sample VCF, 151,561 genomic blocks | 1,473.85 TiB total; avg 9.96 GiB/file |
| **PLINK BED** | Cohort-level binary genotypes | Per chromosome |
| **BGEN** | Cohort-level genotype probabilities | Per chromosome |
| **GDS/aGDS** | Genomic Data Structure; 1,336× more compact than pVCF | ~1.10 TiB total |

### pVCF Notes
- Sample IDs in pVCF headers are pseudonymized EIDs specific to the access application
- Withdrawn participants appear as `W000001`, `W000002`, etc. in headers; their genotype records cannot be identified or removed from the data body

---

## 7. Phenotypic Data

### Storage on RAP
Phenotypic data on the RAP is stored as **Parquet-format tables** in a Spark SQL database — not as the legacy `.enc_ukb` encrypted files used in the classic download model.

The database is named: `app<APPLICATION-ID>_<CREATION-TIME>`

### Column Naming Convention

```
p<FIELD-ID>_i<INSTANCE-ID>_a<ARRAY-ID>
```

| Segment | Meaning | Notes |
|---------|---------|-------|
| `p<FIELD-ID>` | Data field number from UKB Showcase | Always present |
| `_i<INSTANCE-ID>` | Assessment visit (0=initial, 1–3=repeat visits) | Omitted for non-instanced fields |
| `_a<ARRAY-ID>` | Array index for multi-value fields | Omitted for non-arrayed fields |

**Examples:**

| Column | Field | Description |
|--------|-------|-------------|
| `p21022` | 21022 | Age at recruitment (single value) |
| `p53_i0`, `p53_i1` | 53 | Date of attending assessment centre (instanced) |
| `p41270` | 41270 | Diagnoses – ICD10 (array, non-instanced) |
| `p20003_i0_a0` | 20003 | Treatment/medication codes – first visit, first entry |

### Assessment Instances

| Instance | Meaning |
|----------|---------|
| 0 | Initial assessment centre visit |
| 1 | First repeat assessment (invitation-based) |
| 2 | Second repeat assessment |
| 3 | Third repeat assessment |

### Database Tables

| Table Pattern | Contents |
|---------------|----------|
| `participant_0001`–`participant_9999` | Main participant data (horizontally split across tables) |
| `hesin*` | Hospital Episode Statistics (inpatient records) |
| GP records | Primary care data |
| Death records | Death register linkage |
| COVID-19 | SARS-CoV-2 test results |
| Olink proteins | Proteomics data |

### Access Methods

| Method | Tool / Interface |
|--------|-----------------|
| GUI | Cohort Browser — visual field selection and cohort building |
| Export | Table Exporter app — outputs CSV/TSV |
| CLI | `dx extract_dataset` — programmatic field extraction |
| SQL | Spark SQL in JupyterLab — direct queries on Parquet tables |

---

## 8. Sample QC and Relatedness Files

### Sample QC File
**`ukb_sqc_v2.txt`** — one row per participant, key flags include:
- `excess.relatives` — participant has more relatives in dataset than expected
- `used.in.pca.calculation` — participant was included in the PCA computation

### Relatedness File
Pairwise kinship estimates for all related pairs, 5 columns:

| Column | Description |
|--------|-------------|
| EID1 | First participant ID |
| EID2 | Second participant ID |
| HetHet | Fraction of markers where both individuals are heterozygous (KING) |
| IBS0 | Fraction of markers sharing zero alleles (KING) |
| Kinship | Kinship coefficient estimate |

### Kinship Coefficient Thresholds

| Kinship Value | Relationship |
|---------------|-------------|
| > 0.354 | Duplicate or monozygotic twin |
| 0.177 – 0.354 | 1st-degree relative (parent–child, full sibling) |
| 0.0884 – 0.177 | 2nd-degree relative |
| 0.0442 – 0.0884 | 3rd-degree relative |

### Population Stratification Fields

| Field | Contents |
|-------|----------|
| `22006` | Genetic ethnic grouping — flags participants of White British ancestry |
| `22020` | Inclusion in PCA — flags participants used in the PCA computation |

PCA-derived principal components are available as separate fields in the phenotype database.

---

## 9. Key Field ID Reference

### Genotyping Array

| Field | Description |
|-------|-------------|
| `22418` | Genotype calls (PLINK format) |
| `22000` | Array type and batch designation per participant |
| `22006` | Genetic ancestry — White British flag |
| `22020` | Used in PCA calculation flag |

### Imputed Genotypes

| Field | Description |
|-------|-------------|
| `22828` | Imputed genotypes v3 (BGEN) — recommended |
| `22829` | Related imputation field |

### WES

| Field | Description |
|-------|-------------|
| `23157` | WES data — latest release (v11+) |
| `23148` | WES data — v7 and later |
| `23146` | WES data — previous versions |
| `23141`–`23145` | WES metadata |

### WGS

| Field | Description |
|-------|-------------|
| `23149`–`23151` | WGS CRAM files and initial VCF data |
| `23370`–`23384` | 500k WGS expanded releases |

### Data Category Codes

| Category | Contents |
|----------|----------|
| `100315` | Genotyping array data indicators (calls, confidences, intensities, CNV) |
| `100319` | Imputation data indicators |
| `170` | Exome sequencing data |

---

## 10. Reference Genome Summary

> This is critical when combining data types or performing coordinate-based operations.

| Data Type | Reference Genome |
|-----------|-----------------|
| Genotyping array | **GRCh37** (hg19) |
| Imputed genotypes | **GRCh37** (hg19) |
| WES | **GRCh38** (hg38) |
| WGS | **GRCh38** (hg38) |

Liftover between builds is required when merging array/imputed data with WES/WGS data. The DNAnexus app `liftover_plink_beds` (dnanexus-rnd) can convert PLINK files between builds.

---

## 11. Access Tools

| Tool | Purpose |
|------|---------|
| `gfetch` | Original genetic data download tool; chromosome-specific extraction |
| `ukbgene [type] -c[chrom]` | Extract specific data type by chromosome (types: `cal`, `imp`, `hap`) |
| `dx extract_dataset` | RAP CLI tool for phenotypic data extraction |
| DNAnexus SDK (`import dxpy`) | Programmatic access to files, jobs, and datasets |
| Cohort Browser | Web GUI for phenotype filtering and cohort construction |
| `dx get_data_dictionary` | Download field metadata dictionary from RAP |

Maximum simultaneous downloads per application: **10**.

All downloaded data must be encrypted (AES-256) when stored locally; all access is logged with IP addresses.

---

## 12. Participant Withdrawal

- Researchers are notified by email when a participant withdraws consent
- Researchers must remove withdrawn participant records from all local copies
- In pVCF files, withdrawn participants' sample IDs are replaced with `W000001`, `W000002`, etc. in the VCF header
- It is not possible to identify which genotype records within the data body correspond to withdrawn participants
- EID link files are dynamically updated when withdrawals occur

---

## 13. Data Release Versions (as of April 2026)

| Version | Date | Notable Changes |
|---------|------|-----------------|
| v20.5 | November 2025 | 26 new fields + expanded participant coverage |
| v19.1 | March 2025 | 12 new fields + genomic coordinate indices |
| v18.1 | November 2023 | 200k WGS merged into 500k enduring release |
| v17.1 | November 2023 | 29 new fields |
| v16.1 | August 2023 | 15 new fields + expanded imaging |
