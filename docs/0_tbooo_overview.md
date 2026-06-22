# TBOOO — Overview Report

**The Biobank Of Our Own** reproduces the UK Biobank (UKB) directory layout on DNAnexus RAP using only freely available public datasets, so analysis pipelines can be developed and tested without UKB access.

This report covers: what sources are used, what the user gets on disk, and how every populated field is derived.

---

## 1. Data sources

| Dataset | Role in TBOOO | Samples | Reference | Access |
|---|---|---|---|---|
| **1KGP Phase 3** | Array + imputed genotypes | 2,504 unrelated | GRCh37 | EBI FTP — `ftp.1000genomes.ebi.ac.uk` (`ALL.chr*.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`) |
| **1KGP NYGC 30x** | WGS CRAMs, cohort pVCF, simulated WES | 3,202 (incl. 602 trios) | GRCh38 | EBI FTP — `1000G_2504_high_coverage/working/20201028_3202_phased` |
| **SGDP** | WGS diversity (per-sample VCFs only) | ~300 across 142 populations | GRCh38 | ENA project **PRJEB9586** — per-sample phased VCFs only; the ~14 TB of CRAMs is *not* downloaded |
| **GEUVADIS** | Gene-expression PCA → continuous phenotype | 462 (subset of 1KGP, 5 populations) | n/a | EBI ArrayExpress — `E-GEUV-1/GD462.GeneQuantRPKM.50FN.samplename.recast.txt.gz` |

All four are publicly available with no application required. Reference assets fetched alongside: GRCh37 (`hs37d5`) and GRCh38 (`GRCh38DH`) FASTAs, IDT xGen exome v1 BED, and the 1KGP pedigree-based recombination maps.

---

## 2. What the user gets on disk

After `tbooo run --jobs N`, the output under `data/` mirrors what a UKB-approved researcher sees on DNAnexus RAP. The four top-level folders and what lives in each:

| Folder | File pattern | UKB field | Contents |
|---|---|---|---|
| `Bulk/Genotype Results/Genotype calls/` | `ukb22418_c{1-22}_b0_v2.{bed,bim,fam}` | 22418 | Array PLINK (GRCh37) |
| `Bulk/Imputed/` | `ukb22828_c{1-22}_b0_v3.{bgen,bgen.bgi,sample}` | 22828 | Imputed BGEN (GRCh37) |
| `Bulk/Exome sequences/…PLINK format - Final exome release/` | `ukb23158_c{1-22}_b0_v1.{bed,bim,fam}` | 23158 | WES PLINK (GRCh38) |
| `Bulk/Exome sequences/…BGEN format - Final exome release/` | `ukb23159_c{1-22}_b0_v1.{bgen,bgen.bgi}` | 23159 | WES BGEN (GRCh38) |
| `Bulk/Whole genome sequences/{10-29}/` | `<EID>_24051_0_0.{g.vcf.gz,g.vcf.gz.tbi}` | 24051 | Individual gVCFs (1KGP + SGDP), EID-prefix subfolders |
| `Bulk/Whole genome sequences/` | `ukb24310_c{1-22,X}_b0_v1.{pvcf.gz,pvcf.gz.tbi}` | 24310 | Cohort pVCF (NYGC + SGDP) |
| `Showcase/` | `participant.parquet` | various | Synthetic phenotype table, UKB column naming |
| `metadata/` | `eid_map_1kg.tsv`, `eid_map_sgdp.tsv` | — | Sample-ID → synthetic EID maps |
| `metadata/` | `vcf_sample_rename_1kg.txt`, `vcf_sample_rename_sgdp.txt` | — | `bcftools reheader` maps |
| `metadata/` | `ukb_sqc_v2.txt` | — | Sample QC flags |
| `metadata/` | `ukb_rel.txt` | — | Pairwise kinship (KING) |
| `raw/sgdp/` | `sgdp_samples.tsv`, `vcf/<accession>.vcf.gz`, `pvcf/sgdp_c{1-22,X}.pvcf.gz` | — | SGDP staging — **not** part of the UKB-mirrored tree |

Individual-level bulk files (CRAMs, gVCFs) are bucketed into `<first-two-digits-of-EID>/` subfolders (`10/`, `11/`, …, `20/`) — same convention UKB uses on RAP.

### Synthetic participant IDs (EIDs)

Every source sample is assigned a deterministic 7-digit EID, mirroring UKB pseudonymisation:

| Source | EID range | Ordered by |
|---|---|---|
| 1KGP (NYGC 3,202) | 1,000,000 – 1,003,201 | NYGC sample panel order |
| SGDP | 2,000,000 – 2,000,278 | ENA accession |

Individual-level bulk files are routed into `<first-two-digits>/` subfolders (`10/`, `11/`, `20/`, …) — same convention UKB uses.

### Coordinate systems

No liftover is performed; sources are already on the right build for each layer:

| Layer | Reference | Why |
|---|---|---|
| Array, Imputed | GRCh37 | 1KGP Phase 3 is native GRCh37 |
| WES, WGS (CRAM + pVCF) | GRCh38 | 1KGP NYGC and SGDP are native GRCh38 |

---

## 3. How each UKB field is derived

### 3.1 Genotype layers

| Field | What it represents in UKB | Source | Derivation |
|---|---|---|---|
| **22418** — Array PLINK | Affymetrix Axiom array calls, ~805k SNPs, GRCh37 | 1KGP Phase 3 VCFs | `bcftools view --regions-file <Axiom manifest positions>` → `bcftools norm -m-` → `plink2 --vcf --make-bed`. If no manifest is configured, falls back to common biallelic SNPs (MAF ≥ `array_proxy_maf`, rsID present) as a proxy. FAM column 6 (batch) is set to a population-group code (EUR=1, AFR=2, EAS=3, SAS=4, AMR=5) to preserve population-stratified batch structure. cM column populated from 1KGP pedigree-based recombination maps. |
| **22828** — Imputed BGEN v1.2 | ~92.7M variants, fractional dosages, GRCh37 | 1KGP Phase 3 VCFs | `bcftools norm -m-` → `qctool -ofiletype bgen_v1.2 -bgen-bits 8` → `bgenix -index`. Hard genotype calls used directly (info-score = 1.0 for all variants — no dosage uncertainty). `.sample` file written with EIDs in both `ID_1` and `ID_2` and `sex` from the 1KGP panel. |
| **24051** — Per-sample gVCFs | Individual gVCFs, GRCh38 | NYGC + SGDP per-sample VCFs | Built by default. 1KGP samples are extracted from the NYGC cohort VCFs with `bcftools +split` (one pass per chromosome) and concatenated; SGDP per-sample VCFs are symlinked. Output `<EID>_24051_0_0.g.vcf.gz`, tabix-indexed. (Field 24048 CRAMs are not mirrored — the 30x CRAMs are not downloaded.) |
| **24310** — Cohort pVCF | Multi-sample WGS pVCF, GRCh38, one file per chromosome | NYGC per-chromosome phased VCFs + merged SGDP pVCFs | `bcftools reheader --samples vcf_sample_rename_*.txt` swaps original sample IDs for EIDs, then merges NYGC and SGDP per chromosome → `ukb24310_c<chr>_b0_v1.pvcf.gz` + `.tbi`. UKB normally splits this into many blocks; TBOOO keeps one file per chromosome. |
| **23158 / 23159** — WES PLINK + BGEN | Population-level exome calls, GRCh38 | NYGC per-chr VCFs + IDT xGen exome BED | `bcftools view --regions-file <IDT exome BED>` → `bcftools norm -m-` → `bcftools reheader` (EIDs), then both `plink2 --make-bed` (PLINK output, `ukb23158`) and `qctool` + `bgenix` (BGEN output, `ukb23159`). Capture is *simulated by intersection*; true capture efficiency artefacts are not reproduced. |

### 3.2 Phenotype table — `data/Showcase/participant.parquet`

Schema follows UKB naming `p<FIELD>_i<INSTANCE>_a<ARRAY>`. Populated columns:

| Column | UKB Field | Derivation |
|---|---|---|
| `eid` | — | Synthetic 7-digit EID assigned in §2 |
| `p31` | 31 — Sex | `gender` from the 1KGP panel file (or SGDP metadata); encoded 1=male, 2=female |
| `p21000_i0` | 21000 — Ethnic background | Superpopulation → UKB Data-Coding 1001: EUR→1 (White), AFR→4 (Black/Black British), EAS/SAS→3 (Asian/Asian British), AMR→2 (Mixed). SGDP regions mapped analogously, with Central Asia/Siberia and Oceania → 6 (Other) |
| `p22006` | 22006 — Genetic ethnic grouping | 1 if `super_pop=EUR and pop=GBR` (used here as the White-British-ancestry indicator), else 0 |
| `p22020` | 22020 — Used in PCA calculation | 1 for unrelated Phase 3 samples; 0 for NYGC-added relatives and SGDP |
| `p22000` | 22000 — Genotyping batch | Population-group batch code (same scheme as the array FAM column) |
| `p22418` | 22418 — Array data available | 1 for all samples |
| `p22828` | 22828 — Imputed data available | 1 for all samples |
| `p24051` | 24051 — WGS gVCF available | 1 for all samples (per-sample gVCFs built for both cohorts) |
| `p54_i0` | 54 — Assessment centre | Synthetic centre code by population: EUR/GBR→11010 (Leeds), EUR/other→11020, AFR→11021, EAS→11022, SAS→11023, AMR/other→11024 |
| `geuvadis_pc1`…`geuvadis_pc10` | *custom* | See §3.3 |

**Null-but-present columns** (kept so downstream UKB scripts don't break): `p21022` (age), `p53_i0` (visit date), `p21001_i0` (BMI), `p4079_i0`/`p4080_i0` (blood pressure), `p41270` (ICD-10), `p20002_i0_a0`, `p20003_i0_a0`, `p40001_i0`.

### 3.3 GEUVADIS expression PCA (custom columns)

Not a UKB field — it's an additional continuous signal enabled by the GEUVADIS ∩ 1KGP overlap. Pipeline:

1. Download GD462 RPKM matrix (23,722 genes × 462 samples).
2. Filter genes with median RPKM < 0.1.
3. log₂(RPKM + 0.1).
4. Restrict samples to those in `eid_map_1kg.tsv`.
5. `sklearn.PCA(n_components=10)` on the (samples × genes) matrix, centred but not scaled.
6. Join the 10 PC scores onto `participant.parquet` by `eid` as `geuvadis_pc1`…`geuvadis_pc10`.
7. Variance-explained per PC written to `data/metadata/geuvadis_pca_variance.tsv`.

Samples outside the 5 GEUVADIS populations (CEU, FIN, GBR, TSI, YRI) and all SGDP samples receive null. PCs capture a mix of population ancestry, LCL biology, and gene-regulatory variation.

### 3.4 QC and relatedness files

| File | Derivation |
|---|---|
| `ukb_sqc_v2.txt` | Built by `qc.py`. Columns: `FID=0`, `IID=EID`, `used.in.pca.calculation` (1 for unrelated Phase 3, 0 for related NYGC + SGDP), `in.white.british.ancestry.subset` (1 if EUR+GBR), `excess.relatives` (1 if >10 relatives in the 1KGP `.ped`), `putative.sex.chromosome.aneuploidy=0` (not computed), `het.missing.outliers` (>5 SD from mean heterozygosity rate computed by `plink2 --het` on the merged array PLINK) |
| `ukb_rel.txt` | `king --bfile <merged array PLINK> --kinship`, filter pairs to kinship > 0.0442 (3rd-degree threshold, matching UKB), swap original IDs for EIDs. 1KGP NYGC's 602 trios appear here at kinship ≈ 0.25; SGDP samples are unrelated by design |

### 3.5 Identity-mapping helper files (`data/metadata/`)

| File | Contents |
|---|---|
| `eid_map_1kg.tsv` | `eid, sample_id, pop, super_pop, sex, source` |
| `eid_map_sgdp.tsv` | `eid, ena_accession, population, region, sex, source` |
| `vcf_sample_rename_1kg.txt` | `<original_id> <EID>` — used by `bcftools reheader` |
| `vcf_sample_rename_sgdp.txt` | same, for SGDP |

---

## 4. What cannot be derived

| Not produced | Why |
|---|---|
| Hospital Episode Statistics, GP records, death registry | UK NHS linkage; no public equivalent |
| Self-reported disease / medication / lifestyle | No questionnaire data in 1KGP/SGDP |
| Imaging (brain/cardiac MRI, DEXA) | Not collected |
| Proteomics, metabolomics, accelerometry | Not collected |
| Repeat-assessment instances (i1, i2, i3) | All public samples are single-timepoint |
| True imputation dosage uncertainty | Direct VCF→BGEN conversion gives 0/1 probabilities only |
| True exome capture artefacts | WES is simulated by intersecting WGS with the exome BED |

---

## 5. Limitations at a glance

- **Cohort size**: max 3,481 samples vs. UKB's ~500,000 → much lower rare-variant power.
- **Population composition**: globally diverse vs. UKB's ~94% European → allele-frequency distributions diverge.
- **Array overlap is partial**: only positions present in both 1KGP Phase 3 and the Axiom manifest survive.
- **SGDP is WGS-only in this build**: it does not feed the array, imputed, or WES layers.
- **GEUVADIS coverage is partial**: only 462 of ~3,481 samples have non-null `geuvadis_pc*` scores.

See [3_data_mapping.md](3_data_mapping.md) for the full per-command spec and [2_data_sources.md](2_data_sources.md) for source-level detail.
