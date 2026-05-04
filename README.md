# TBOOO - The Biobank Of Our Own

A toolkit that reproduces the UK Biobank (UKB) data structure on DNAnexus RAP using freely available public datasets - **1000 Genomes Project** and the **Simons Genome Diversity Project** - so analysis pipelines can be developed, tested, and validated without requiring UKB access.

---

## What it produces

Running the full pipeline yields a directory tree that mirrors what a UKB researcher sees on DNAnexus RAP:

```
data/
├── Bulk/
│   ├── Genotype Results/Genotype calls/
│   │   └── ukb22418_c{1-22}_b0_v2.{bed,bim,fam}           # Field 22418 - array PLINK
│   ├── Imputed/
│   │   └── ukb22828_c{1-22}_b0_v3.{bgen,bgen.bgi,sample}  # Field 22828 - imputed BGEN
│   ├── Exome sequences/
│   │   ├── Population level exome OQFE variants, PLINK format - 500k release/
│   │   │   └── ukb23157_c{1-22}_b0_v1.{bed,bim,fam}        # Field 23157 - WES PLINK
│   │   └── Population level exome OQFE variants, BGEN format - final release/
│   │       └── ukb23157_c{1-22}_b0_v1.{bgen,bgen.bgi}      # Field 23157 - WES BGEN
│   └── Whole genome sequences/
│       ├── 10/1000001_23149_0_0.{cram,cram.crai}            # Field 23149 - individual CRAMs (1KGP)
│       ├── 20/2000001_23151_0_0.{g.vcf.gz,g.vcf.gz.tbi}    # Field 23151 - individual gVCFs (SGDP)
│       └── ukb23370_c{1-22,X}_b0_v1.{pvcf.gz,pvcf.gz.tbi} # Field 23370 - 1KGP cohort pVCF
├── Showcase/
│   └── participant.parquet                                   # synthetic phenotype table
├── raw/
│   └── sgdp/
│       ├── sgdp_samples.tsv                                  # SGDP sample metadata
│       ├── vcf/<accession>.vcf.gz                            # per-sample phased VCFs (from ENA)
│       └── pvcf/sgdp_c{1-22,X}.{pvcf.gz,pvcf.gz.tbi}       # merged SGDP pVCF per chromosome
└── metadata/
    ├── eid_map_1kg.tsv                                       # 1KGP sample → synthetic EID
    ├── eid_map_sgdp.tsv                                      # SGDP sample → synthetic EID
    ├── vcf_sample_rename_1kg.txt                             # bcftools reheader map (1KGP)
    ├── vcf_sample_rename_sgdp.txt                            # bcftools reheader map (SGDP)
    ├── ukb_sqc_v2.txt                                        # sample QC flags
    └── ukb_rel.txt                                           # pairwise kinship (KING)
```

---

## Data sources

| Dataset | Role in TBOOO | Samples | Reference |
|---------|--------------|---------|-----------|
| **1KGP Phase 3** | Array + imputed genotypes | 2,504 unrelated | GRCh37 |
| **1KGP NYGC 30x** | WGS CRAMs, cohort pVCF, WES simulation | 3,202 (incl. trios) | GRCh38 |
| **SGDP** | Per-sample VCFs + phenotype diversity | 300 (142 populations) | GRCh38 |
| **GEUVADIS** | Gene expression PCA → custom phenotype signal | 462 (5 populations, all in 1KGP) | GRCh37 |

SGDP CRAMs (~14 TB) are not downloaded. Per-sample phased VCF files are fetched from ENA analysis results and merged into per-chromosome pVCFs with EID renaming. GEUVADIS data is downloaded as the pre-computed gene-level RPKM matrix (no raw reads); PCA scores are stored as `geuvadis_pc1`–`geuvadis_pc10` columns in `participant.parquet`. All data is publicly available without application. See [docs/2_data_sources.md](docs/2_data_sources.md) for full access details.

---

## Requirements

### Python dependencies

```bash
pip install -e .
```

Requires Python ≥ 3.10. Installs: `click`, `pyyaml`, `pandas`, `pyarrow`, `tqdm`, `scikit-learn`, `numpy`.

### External bioinformatics tools

| Tool | Min version | Purpose |
|------|-------------|---------|
| `bcftools` | 1.17 | VCF filtering, normalization, reheadering, tabix indexing |
| `samtools` | 1.17 | FASTA indexing, CRAM operations |
| `plink2` | 2.00a3 | VCF → PLINK BED/BIM/FAM, heterozygosity stats, PLINK merge |
| `qctool` | 2.0.8 | VCF → BGEN v1.2 conversion |
| `bgenix` | 1.1.7 | BGEN index generation (.bgi files) |
| `king` | 2.3 | Pairwise kinship estimation |
| `wget` | any | File downloads (resumable) |
| `tar`, `gzip` | any | Archive extraction — pre-installed on both platforms |

All tool paths can be overridden in `config.yaml` under the `tools:` key.

#### Ubuntu

```bash
# ── 1. System basics ──────────────────────────────────────────────────────────
sudo apt-get update
sudo apt-get install -y wget curl unzip tar gzip build-essential zlib1g-dev

# ── 2. bcftools + samtools ────────────────────────────────────────────────────
# Ubuntu 24.04 ships bcftools 1.18 / samtools 1.19 — meets requirements:
sudo apt-get install -y bcftools samtools
# Ubuntu 22.04 ships 1.13 (too old). Use conda instead:
#   conda install -c bioconda bcftools samtools

# ── 3. plink2 ─────────────────────────────────────────────────────────────────
# Download the latest AVX2 build from https://www.cog-genomics.org/plink/2.0/
# Replace the filename with the current release date shown on that page.
wget https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_avx2_20260504.zip
unzip plink2_linux_avx2_20260504
sudo install plink2 /usr/local/bin/

# ── 4. qctool v2 ─────────────────────────────────────────────────────────────
# Prebuilt CentOS/Linux binary from https://www.well.ox.ac.uk/~gav/qctool_v2/
# (CentOS binaries run on Ubuntu via glibc compatibility)
wget "https://www.well.ox.ac.uk/~gav/resources/qctool_v2.2.0-CentOS_Linux7.8.2003-x86_64.tgz"
tar -xzf qctool_v2.2.0-CentOS_Linux7.8.2003-x86_64.tgz
sudo install qctool_v2.2.0-CentOS\ Linux7.8.2003-x86_64/qctool /usr/local/bin/

# ── 5. bgenix ─────────────────────────────────────────────────────────────────
# Part of the BGEN reference implementation — build from source.
# Project page: https://enkre.net/cgi-bin/code/bgen
# Only requires a C++11 compiler (already in build-essential); deps are bundled.
wget http://code.enkre.net/bgen/tarball/release/bgen.tgz
tar -xzf bgen.tgz
mv bgen.tgz bgen
cd bgen
# The release tarball uses std::ios::streampos which GCC 11+ rejects; patch before building.
sed -i 's/std::ios::streampos/std::streampos/g' src/View.cpp
./waf configure
./waf
sudo install build/apps/bgenix /usr/local/bin/
cd ..

# ── 6. KING ───────────────────────────────────────────────────────────────────
# Prebuilt Linux binary from https://www.kingrelatedness.com/
# Replace URL with the latest version shown on the downloads page.
wget https://www.kingrelatedness.com/executables/Linux-king231.tar.gz
tar -xzf Linux-king231.tar.gz
sudo install king /usr/local/bin/
```

---

## Configuration

All settings live in `config.yaml`:

```yaml
data_dir: data               # output root
reference_dir: data/reference
autosomes: [1, 2, ..., 22]  # chromosomes to process
sex_chromosomes: ["X"]

# optional: path to Thermo Fisher Axiom UKB array manifest
# if blank, common biallelic SNPs (MAF >= 0.05, rsID present) are used as proxy
array_manifest: ""
array_proxy_maf: 0.05

# optional: restrict SGDP metadata to specific populations
# empty list = all 300 samples across 142 populations
sgdp_populations: []
```

---

## Usage

### Step 1 - Download reference files

```bash
tbooo download reference
```

Downloads: IDT xGen exome capture BED (GRCh38), 1KGP pedigree-based genetic maps, GRCh37 and GRCh38 reference FASTAs.

### Step 2 - Download source data

```bash
# 1000 Genomes Phase 3 VCFs (GRCh37) + NYGC 30x VCFs (GRCh38)
tbooo download 1kg

# SGDP metadata + per-sample VCFs from ENA (no CRAMs)
tbooo download sgdp

# Skip VCF download (metadata only)
tbooo download sgdp --no-vcf

# Download only specific chromosomes
tbooo download 1kg --chroms 1,2,22
```

### Step 3 - Build the mirrored structure

Each step can be run individually or all at once via `tbooo run`.

```bash
# Must run first - assigns synthetic EIDs to all samples
tbooo map eids

# Build each data layer
tbooo map array       # Phase 3 VCF → PLINK (Field 22418, GRCh37)
tbooo map imputed     # Phase 3 VCF → BGEN  (Field 22828, GRCh37)
tbooo map wgs               # rename 1KGP CRAMs + build merged pVCF (NYGC+SGDP → Field 23370)
tbooo map wgs --gvcf        # also extract per-sample gVCFs from NYGC + SGDP (Field 23151)
tbooo map wgs --no-pvcf     # CRAMs only (skip pVCF build)
tbooo map wes         # NYGC VCF ∩ exome BED → PLINK + BGEN (Field 23157, GRCh38)
tbooo map geuvadis    # GEUVADIS RPKM → expression PCA → geuvadis_pc* in participant.parquet
tbooo map phenotypes  # EID maps → Parquet with UKB column naming (p<FIELD>_i<INST>_a<ARR>)
tbooo map qc          # PLINK --het + KING → ukb_sqc_v2.txt + ukb_rel.txt

# Process a subset of chromosomes
tbooo map array --chroms 1,2,22
tbooo map wes --chroms 22
```

### Run the full pipeline

```bash
# Full pipeline: downloads + all stages, 16 parallel chromosome jobs
tbooo run --jobs 16

# Dry run - print every step without executing
tbooo run --jobs 16 --dry-run

# Skip downloads (data already present), run all build stages
tbooo run --jobs 16 --no-download

# Run a single stage (prerequisites must already exist)
tbooo run --target array
tbooo run --target imputed
tbooo run --target wes
tbooo run --target wgs
tbooo run --target geuvadis
tbooo run --target phenotypes
tbooo run --target qc

# Limit to specific chromosomes
tbooo run --jobs 8 --chroms 1,2,3 --no-download
```

---

## Phenotype table

The synthetic phenotype table (`data/Showcase/participant.parquet`) follows UKB column naming (`p<FIELD-ID>_i<INSTANCE>_a<ARRAY>`). Populated fields:

| Column | UKB Field | Source |
|--------|-----------|--------|
| `eid` | - | Synthetic 7-digit ID |
| `p31` | 31 | Sex (from sample panel) |
| `p21000_i0` | 21000 | Ethnic background (mapped from superpopulation / SGDP region) |
| `p22006` | 22006 | White British ancestry flag |
| `p22020` | 22020 | Used in PCA calculation |
| `p22000` | 22000 | Genotyping batch code |
| `p22418` | 22418 | Array data available |
| `p22828` | 22828 | Imputed data available |
| `p23149` | 23149 | WGS CRAM available (1KGP samples only) |
| `p23151` | 23151 | Individual gVCF available (SGDP samples only) |
| `p54_i0` | 54 | Assessment centre (synthetic) |
| `geuvadis_pc1`–`geuvadis_pc10` | custom | Gene expression PC scores (462 1KGP samples only; null for all others) |

Clinical fields (diagnoses, medications, hospital records, imaging) are present as null columns so downstream scripts that reference them do not break.

---

## Documentation

| File | Contents |
|------|----------|
| [docs/1_ukb_structure.md](docs/1_ukb_structure.md) | UK Biobank data structure on DNAnexus RAP - file naming, field IDs, formats, reference genomes |
| [docs/2_data_sources.md](docs/2_data_sources.md) | 1KGP and SGDP data - phases, sample counts, VCF paths, access methods |
| [docs/3_data_mapping.md](docs/3_data_mapping.md) | Mapping specification - how each source file becomes each UKB-mirrored output |

---

## Limitations

- **Sample size**: ≤ 3,481 samples vs. UKB's ~500,000. Rare-variant power is substantially lower.
- **Array coverage**: Only variants present in both the 1KGP Phase 3 VCFs and the array manifest are included. Affymetrix-specific probes with no 1KGP equivalent are absent.
- **Imputation dosages**: BGEN files are built from hard genotype calls; dosage probabilities are 0 or 1. Real UKB imputed data has fractional uncertainty.
- **Population composition**: 1KGP/SGDP are globally diverse; UKB is ~94% European. Allele frequency distributions differ.
- **SGDP CRAMs not downloaded**: Only per-sample phased VCFs are fetched from ENA (much smaller than ~14 TB of CRAMs). SGDP genotype data feeds into `data/raw/sgdp/pvcf/` - it is not merged into the UKB-mirrored Bulk/ tree, which remains 1KGP-only.
- **No clinical phenotypes**: Only metadata-derivable fields are populated (sex, ancestry, batch). Disease diagnoses, hospital records, imaging, and longitudinal data are not available.
- **WES is simulated**: Intersection of WGS variants with an exome capture BED is used, not true capture sequencing. Capture efficiency gradients and strand bias differ from real exome data.
