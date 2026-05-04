# Public Data Sources: 1000 Genomes Project, Simons Genome Diversity Project & GEUVADIS

This document details the available data in the three public datasets used by TBOOO: the **1000 Genomes Project (1KGP)**, the **Simons Genome Diversity Project (SGDP)**, and the **GEUVADIS RNA-seq project**.

---

## 1. 1000 Genomes Project (1KGP)

### 1.1 Phases Overview

| Phase | Year | Samples | Depth | Key Characteristics |
|-------|------|---------|-------|---------------------|
| Phase 1 | ~2012 | ~1,000 | Low (mixed) | Mixed read lengths (36bp–100bp+); multiple platforms (Illumina, ABI SOLiD, 454); proof-of-concept |
| Phase 3 | 2013–2015 | 2,504 | Low (~7x) | Illumina-only; statistical variant calling (SHAPEIT2 phasing, MVNcall); reference standard for population genetics |
| NYGC 30x | 2022 | 3,202 | 30x | High-coverage resequencing; includes 698 additional related samples + 602 complete trios; GRCh38 primary |

### 1.2 Sample Counts and Population Distribution

**Phase 3 — 2,504 unrelated individuals, 26 populations:**

| Superpopulation | Code | Populations included | Sample count |
|-----------------|------|---------------------|-------------|
| African | AFR | YRI, LWK, GWD, MSL, ESN, ASW, ACB | ~661 |
| American | AMR | MXL, PUR, CLM, PEL | ~347 |
| East Asian | EAS | CHB, JPT, CHS, CDX, KHV | ~504 |
| European | EUR | CEU, TSI, FIN, GBR, IBS | ~503 |
| South Asian | SAS | GIH, PJL, BEB, STU, ITU | ~489 |

**NYGC 30x — 3,202 individuals (2,504 unrelated + 698 related):**

| Superpopulation | Sample count |
|-----------------|-------------|
| AFR | 893 |
| EUR | 633 |
| EAS | 601 |
| SAS | 585 |
| AMR | 490 |

The 698 additional samples include nearly all parent–child trios; **602 complete trios** are present in the collection.

### 1.3 Population Codes — 26 Populations

| Code | Population | Superpopulation |
|------|-----------|-----------------|
| CHB | Han Chinese in Beijing, China | EAS |
| JPT | Japanese in Tokyo, Japan | EAS |
| CHS | Southern Han Chinese | EAS |
| CDX | Chinese Dai in Xishuangbanna, China | EAS |
| KHV | Kinh in Ho Chi Minh City, Vietnam | EAS |
| CEU | Utah Residents (CEPH) with Northern and Western European Ancestry | EUR |
| TSI | Toscani in Italia | EUR |
| FIN | Finnish in Finland | EUR |
| GBR | British in England and Scotland | EUR |
| IBS | Iberian Population in Spain | EUR |
| YRI | Yoruba in Ibadan, Nigeria | AFR |
| LWK | Luhya in Webuye, Kenya | AFR |
| GWD | Gambian in Western Divisions in the Gambia | AFR |
| MSL | Mende in Sierra Leone | AFR |
| ESN | Esan in Nigeria | AFR |
| ASW | Americans of African Ancestry in SW USA | AFR |
| ACB | African Caribbeans in Barbados | AFR |
| MXL | Mexican Ancestry from Los Angeles USA | AMR |
| PUR | Puerto Ricans from Puerto Rico | AMR |
| CLM | Colombians from Medellín, Colombia | AMR |
| PEL | Peruvians from Lima, Peru | AMR |
| GIH | Gujarati Indian from Houston, Texas | SAS |
| PJL | Punjabi from Lahore, Pakistan | SAS |
| BEB | Bengali from Bangladesh | SAS |
| STU | Sri Lankan Tamil from the UK | SAS |
| ITU | Indian Telugu from the UK | SAS |

### 1.4 Variant Statistics

**Phase 3 integrated variant set (2,504 samples):**

| Variant Type | Count | Notes |
|-------------|-------|-------|
| Total variants | 88 million+ | All phased onto haplotypes |
| SNVs | 81.4 million | 99.6% biallelic |
| Short indels | 3.6 million | |
| Structural variants | ~60,000 | 5 SV types; DGVa accession nstd145 |

**NYGC 30x high-coverage set (3,202 samples):**

| Variant Type | Count | vs. Phase 3 |
|-------------|-------|-------------|
| SNVs | 96,950,998 | 1.24× more |
| INDELs | 13,132,415 | 4.05× more |
| SVs | Comprehensive ML-integrated set | Multi-technology; improved sensitivity |

### 1.5 VCF File Naming and Paths

**Phase 3 — per-chromosome multi-sample VCFs (GRCh37):**

```
ALL.chr{1-22}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
ALL.chrX.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz
ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
```

Corresponding tabix index files: `<filename>.vcf.gz.tbi`

Sites-only VCF (no genotypes):
```
ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz
```

Full EBI FTP path:
```
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
```

**NYGC 30x — per-chromosome phased VCFs (GRCh38):**

```
CCDG_14151_B01_GRM_WGS_2020-08-05_chr{N}.filtered.shapeit2-duohmm-phased.vcf.gz
```

Example full path (EBI):
```
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz
```

### 1.6 CRAM/BAM Alignment Files

Individual-level alignment files are available for each sample:

| File type | Extension | Description |
|-----------|-----------|-------------|
| CRAM | `.cram` | Compressed reference-based alignment (30–60% smaller than BAM) |
| CRAM index | `.crai` | Per-sample random-access index |
| BAM | `.bam` | Uncompressed alignment (legacy) |
| BAM index | `.bai` | Per-sample random-access index |
| Statistics | `.cram.bas` or `.bam.bas` | Tab-separated; one row per read group; alignment statistics |

CRAM files require the reference genome for decoding. Always download the matching `.crai` index alongside `.cram` files.

### 1.7 Structural Variant Data

| Attribute | Details |
|-----------|---------|
| SV types | DEL (deletions), DUP (duplications), INV (inversions), CNV (copy number variants), INS (insertions) |
| Database submission | DGVa, accession **nstd145** |
| Reference coordinates | GRCh37 primary; GRCh38 liftover in `supporting/GRCh38_positions/` subdirectory |
| Technologies integrated | Illumina WGS, Complete Genomics WGS, PacBio, MinION, aCGH, microarray platforms |
| NYGC SVs | Comprehensive ML-integrated set with improved sensitivity over Phase 3 |

SV VCF path (EBI FTP):
```
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/
```

### 1.8 Access Paths

**EBI FTP (primary):**
```
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/         # Phase 3 VCFs
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/ # Phase 3 SVs
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/ # NYGC 30x
```

**EBI HTTPS mirror:**
```
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
```

**AWS S3 (no authentication required):**
```
s3://1000genomes    # Region: us-east-1
```

Download example:
```bash
aws s3 cp s3://1000genomes/phase3/data/NA21144/sequence_read/ERR047877.filt.fastq.gz ./ --no-sign-request
```

No AWS account required for read access. Users pay only for compute resources.

**IGSR Data Portal:**
```
https://www.internationalgenome.org/data-portal/
https://www.internationalgenome.org/data-portal/data-collection/30x-grch38   # NYGC 30x
```

The portal supports filtering by individual, population, data type, and sequencing technology, and provides direct download links or export to download managers.

**Recommended transfer method for large volumes:** Globus (GridFTP protocol; fault-tolerant bulk transfers).

### 1.9 Imputation Reference Panels

1KGP Phase 3 is a standard imputation backbone for population genetics and GWAS:

| Format | Description |
|--------|-------------|
| HAP/LEGEND | `1000GP_Phase3_chr<N>_GRCh37.hap.gz` + `.legend.gz`; compatible with SHAPEIT4 |
| VCF | SHAPEIT4-compatible VCF reference panel |

Compatible imputation tools: **Beagle**, **IMPUTE2**, **Mach/Minimac**, **SHAPEIT2/SHAPEIT4**.

For UK Biobank imputation, 1KGP Phase 3 was used as a fallback panel (for variants absent from the HRC).

### 1.10 Recombination Maps

High-resolution pedigree-based genetic maps constructed from 1KGP data:

| Statistic | Value |
|-----------|-------|
| Crossovers analyzed | 3.38 million |
| Putative hotspots | 88,841 |
| Putative coldspots | 80,129 |
| Population specificity | Population-specific recombination rates reflecting ancestry |

**GitHub repository:** `joepickrell/1000-genomes-genetic-maps`

### 1.11 Sample Metadata Files

| File | Contents |
|------|----------|
| `integrated_call_samples_v3.20130502.ALL.panel` (54 KB) | All 2,504 samples: individual ID, population, superpopulation |
| `integrated_call_male_samples_v3.20130502.ALL.panel` (25 KB) | Male samples only |
| `.ped` files (2020, 2025 versions) | Full pedigree/family structure for 3,202-sample set |

Columns in `.panel` files: `sample`, `pop`, `super_pop`, `gender`

### 1.12 Reference Genomes

| Build | Usage in 1KGP | Notes |
|-------|--------------|-------|
| GRCh37 (hg19) | Phase 3 original coordinates | 1-based; all Phase 3 VCFs |
| GRCh38 (hg38) | NYGC 30x primary | Preferred for new analyses |
| hs37d5 | Some alignment files | GRCh37 + decoy sequences; improves mapping of reads from repetitive regions |

---

## 2. Simons Genome Diversity Project (SGDP)

### 2.1 Overview

The SGDP was designed to maximize anthropological, linguistic, and cultural diversity — filling geographic gaps underrepresented in datasets like 1KGP.

| Attribute | Value |
|-----------|-------|
| Total genomes | 300 |
| Publicly available | 279 (263 C-panel + 16 B-panel) |
| Restricted access | 21 genomes |
| Additional (Fan et al.) | 44 samples |
| Populations | 142 distinct populations |
| Continents | 6 |
| Sequencing depth | ≥ 30x |
| Sequencing platform | Illumina |
| Primary reference | hs37d5 (hg19 + decoy) |
| Secondary reference | GRCh38DH |

### 2.2 Population Distribution by Continent

| Continent / Region | Individuals |
|--------------------|------------|
| West Eurasia | 75 |
| East Asia | 47 |
| South Asia | 39 |
| Africa | 44 |
| Central Asia / Siberia | 27 |
| Oceania | 25 |
| Native Americas | 22 |
| **Total** | **279 (public)** |

### 2.3 Sequencing Specifications

- **Coverage:** All 300 genomes sequenced to minimum **30x** depth
- **Platform:** Illumina
- **Alignment strategy:** Customized BWA-MEM procedure optimized to eliminate reference bias; genotyping performed on a per-sample basis to avoid preferential variant calling in over-represented populations
- **Primary reference:** `hs37d5` — GRCh37 with decoy sequences (same reference as 1KGP Phase 3 alignments)
- **Secondary reference:** `GRCh38DH` — GRCh38 with decoy + HLA sequences; used for newer alignments available via UPPMAX

### 2.4 Data Sections and File Formats

The SGDP data is organized into labeled sections:

#### Section A — Raw Alignments

| File type | Extension | Description |
|-----------|-----------|-------------|
| CRAM | `.cram` | Compressed alignments against hs37d5 or GRCh38DH |
| CRAM index | `.crai` | Per-sample random-access index |
| BAM statistics | `.bas` | Tab-separated read group statistics |

**CRAM file naming pattern (GRCh38DH alignments):**
```
SAMEA[ENA-ID].alt_bwamem_GRCh38DH.20200922.[Population].simons.cram
```

Example:
```
SAMEA3302732.alt_bwamem_GRCh38DH.20200922.Greek.simons.cram
SAMEA3302732.alt_bwamem_GRCh38DH.20200922.Greek.simons.cram.crai
```

**Total CRAM data volume: ~14 TB** (all 279 public samples + indices)

#### Section B — Compact Formats

| Dataset | Description |
|---------|-------------|
| SGDP-lite v3 | ~140 GB compressed; hetfa format for rapid regional queries |
| Ctools | Population genetics software; operates on hetfa; available via GitHub |

#### Section E — Phased Genotypes (2021)

Imputed and phased genotype calls generated using **Glimpse** against the 1000 Genomes Phase 3 reference panel.

| VCF field | Contents |
|-----------|----------|
| `GT` | Phased genotypes |
| `DS` | Dosage scores |
| `PP` | Posterior genotype probabilities |

#### Section I — SNP Data

PLINK-format genotype matrix for all samples:

| File | Contents |
|------|----------|
| `.bim` | Variant manifest (chromosome, rsID, position, alleles) |
| `.bed` | Binary genotype matrix |
| `.fam` | Sample manifest |

Filtering: **MAF ≥ 0.001** applied before inclusion.

#### Mitochondrial and Y-Chromosome Data

| Dataset | Size | Notes |
|---------|------|-------|
| Y-chromosome BAMs | 129 GB tarball | Male samples only |
| Mitochondrial BAMs | 45 GB | 5,344× mean coverage per sample |

Short tandem repeats (STRs) are available through dbVar (accession **nstd128**).

### 2.5 Access Methods

#### ENA — Public Samples (279 genomes)

| Attribute | Value |
|-----------|-------|
| Study accession | `PRJEB9586` / `ERP010710` |
| Access | Open; no registration required |
| FTP pointers | `ena.ftp.pointers.txt` — lists direct BAM/CRAM FTP URLs for each sample |

#### EGA — Restricted Samples (21 genomes)

| Attribute | Value |
|-----------|-------|
| Study accession | `EGAS00001001959` |
| Access | Requires signed usage agreement (see §2.7) |

#### Harvard Reich Lab (institutional repository)

```
https://reichdata.hms.harvard.edu/pub/datasets/sgdp/
```
Alternative mirror:
```
https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/
```

Direct HTTP/FTP access to all public data files. Organized by section (A, B, E, I) and by sample.

#### Seven Bridges Cancer Genomics Cloud (CGC)

```
https://cgc-accounts.sbgenomics.com/auth/login
```

- **279 public genomes** from 130 diverse populations available
- Free account registration; new users receive $300 cloud credits
- Data copied to personal projects does not count against storage quotas
- Platform includes 900+ pre-built bioinformatic tools and workflows
- No special authorization required for public data

#### IGSR Data Portal

```
https://www.internationalgenome.org/data-portal/data-collection/SGDP/
```

SGDP is indexed as a data collection in the International Genome Sample Resource (IGSR), alongside 1KGP and HGDP. VCFs organized by chromosome; downloadable via Globus, FTP, or Aspera.

#### Globus (recommended for large transfers, 2024+)

Primary mechanism for both public and restricted samples. Managed endpoint operated by SGDP maintainers. Provides GridFTP for fault-tolerant bulk transfers and HTTPS for individual file access.

### 2.6 Access Requirements

**Public samples (279 genomes):**
- Fully open; no usage agreement required beyond standard open data licenses
- Download via ENA, Harvard repository, CGC, IGSR, or Globus

**Restricted samples (21 genomes):**
Require a dated, signed letter affirming:
- No redistribution to third parties
- Non-commercial use (or specified commercial terms)
- No linkage to personal identifiers

Contact:
```
shop@genetics.med.harvard.edu    # Shop Mallick (David Reich Lab)
reich@genetics.med.harvard.edu   # David Reich
```

Requests processed fastest through the Globus mechanism.

### 2.7 Companion Software

| Tool | Purpose | Source |
|------|---------|--------|
| **Ctools** | Population genetics analysis; designed for SGDP hetfa format | GitHub (C + Python) |
| **Glimpse** | Imputation and phasing from sparse data vs. reference panel | GitHub |

---

## 3. GEUVADIS RNA-seq

### 3.1 Overview

GEUVADIS (Genetic European Variation in Disease) is a collaborative RNA-seq study of lymphoblastoid cell lines (LCLs) derived from 1KGP individuals. Because every GEUVADIS sample is also a 1KGP sample, expression measurements can be directly linked to the matching WGS genotype data — making GEUVADIS the only publicly available dataset that pairs high-coverage WGS with bulk RNA-seq at population scale.

| Attribute | Value |
|-----------|-------|
| Publication | Lappalainen et al., *Nature* 2013 |
| Total samples | 462 |
| Tissue | LCL (EBV-transformed lymphoblastoid cell lines) |
| Sequencing | Illumina HiSeq 2000, strand-specific, ~20M 75 bp reads per sample |
| Reference genome | GRCh37 |
| Gene models | GENCODE v12 |
| All samples are | 1KGP individuals (direct overlap with Phase 3 genotype data) |

### 3.2 Population Distribution

| Population | Code | Superpop | Sample count |
|-----------|------|----------|-------------|
| Utah Residents (CEPH) with Northern/Western European Ancestry | CEU | EUR | 91 |
| Finnish in Finland | FIN | EUR | 95 |
| British in England and Scotland | GBR | EUR | 89 |
| Toscani in Italia | TSI | EUR | 98 |
| Yoruba in Ibadan, Nigeria | YRI | AFR | 89 |
| **Total** | | | **462** |

All 462 samples overlap with the 1KGP Phase 3 2,504-sample unrelated set. No non-Phase-3 individuals are included.

### 3.3 Data Types and Processed Files

#### Raw data (ENA)

| File type | Contents |
|-----------|----------|
| FASTQ | Paired-end raw reads, ~75 bp, strand-specific |
| BAM | Alignments to GRCh37 + GENCODE v12 |

#### Processed quantification (GEUVADIS FTP / ArrayExpress)

| File | Contents |
|------|----------|
| `GD462.GeneQuantRPKM.50FN.samplename.recast.txt.gz` | 462 samples × 23,722 genes, RPKM normalized gene-level expression |
| `GD462.GeneQuantCount.50FN.samplename.recast.txt.gz` | Same matrix, raw read counts |
| `GD660.GeneQuantRPKM.50FN.samplename.recast.txt.gz` | Expanded 660-sample release (includes technical replicates) |
| `EUR373.gene.cis.FDR5.all.rs137.txt.gz` | cis-eQTL results, European samples, FDR < 5% |
| `YRI89.gene.cis.FDR5.all.rs137.txt.gz` | cis-eQTL results, YRI samples, FDR < 5% |

TBOOO uses the **GD462** RPKM matrix as the source for expression PCA. The 660-sample release contains technical replicates and is not used.

### 3.4 Access Methods

**ArrayExpress (primary processed data):**
```
https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/   # EUR populations
https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-3/   # YRI replication set
```

**ENA (raw reads + alignments):**
```
Study accession: PRJEB3366
https://www.ebi.ac.uk/ena/browser/view/PRJEB3366
```

**GEUVADIS FTP (processed matrices, eQTL results):**
```
https://www.ebi.ac.uk/Tools/geuvadis-das/
ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1/
```

All data is fully open; no access agreement is required.

### 3.5 Sample Overlap with 1KGP

All 462 GEUVADIS sample IDs (e.g. `NA12878`, `HG00096`) match 1KGP sample IDs in the Phase 3 panel file. The `eid_map_1kg.tsv` generated by `tbooo map eids` is the join key: GEUVADIS sample ID → EID.

The remaining ~2,042 Phase 3 samples (non-CEU/FIN/GBR/TSI/YRI populations) and all SGDP samples have no GEUVADIS data and receive null values for expression PCA columns.

---

## 4. Comparative Summary

| Feature | 1KGP Phase 3 | 1KGP NYGC 30x | SGDP | GEUVADIS |
|---------|-------------|---------------|------|----------|
| **Total samples** | 2,504 | 3,202 | 300 | 462 |
| **Unrelated samples** | 2,504 | 2,504 | 279 (public) | 462 (all are Phase 3 unrelated) |
| **Family structure** | Limited trios | 602 complete trios | Not primary focus | None |
| **Populations** | 26 | 26 | 142 | 5 (CEU, FIN, GBR, TSI, YRI) |
| **Geographic scope** | 5 continental groups | 5 continental groups | 6 continents; maximally diverse | EUR + YRI only |
| **Data type** | WGS genotypes | WGS genotypes + alignments | WGS genotypes + alignments | RNA-seq (LCL) |
| **Sequencing depth** | Low (~7x, mixed) | 30x uniform | 30x uniform | ~20M reads/sample |
| **SNVs** | 81.4 million | 96.95 million | Adds diversity; imputed vs. 1KGP | — |
| **Short indels** | 3.6 million | 13.1 million | Included in variant data | — |
| **Structural variants** | ~60,000 (5 types) | Comprehensive ML set | Limited; STRs via dbVar nstd128 | — |
| **Primary reference** | GRCh37 | GRCh38 | hs37d5 (GRCh37 + decoy) | GRCh37 |
| **Secondary reference** | GRCh38 (liftover) | — | GRCh38DH | — |
| **Individual alignments** | CRAM/BAM per sample | CRAM/BAM per sample | CRAM per sample (~14 TB) | BAM per sample (ENA) |
| **Processed matrices** | — | — | — | RPKM gene expression (GD462) |
| **Cloud: AWS S3** | Yes (`s3://1000genomes`, free) | Yes | No | No |
| **Cloud: CGC** | Partial | Partial | Yes (279 public, free) | No |
| **FTP** | EBI FTP (primary) | EBI FTP | ENA, Harvard Reich Lab | EBI ArrayExpress FTP |
| **Globus** | Yes | Yes | Yes (2024+, primary) | No |
| **Access restrictions** | Open (all public) | Open (all public) | 279 open; 21 restricted (signed agreement) | Open (all public) |
| **mtDNA / Y data** | Yes | Yes | Yes | — |
| **TBOOO role** | Primary genotype + population backbone | WGS/WES layer | Population diversity beyond 1KGP | Expression PCA → phenotype-like signal |
