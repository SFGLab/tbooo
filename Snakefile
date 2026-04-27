"""
TBOOO Snakemake workflow
========================
Produces a UK Biobank-mirrored data structure from 1000 Genomes + SGDP.

Usage:
    snakemake --configfile config.yaml --cores 8
    snakemake --configfile config.yaml --cores 8 --dry-run
    snakemake --configfile config.yaml --cores 8 all_array
    snakemake --configfile config.yaml --cores 8 all_imputed

Full pipeline DAG:
    download_phase3_vcfs
    download_nygc_vcfs    ─┐
    download_sgdp_meta     │   (metadata only — no CRAMs)
    download_reference     │
         │                 │
    assign_eids            │
         │        ─────────┘
    ┌────┴────┬───────┬──────────┬────┐
  array   imputed  wgs_crams  wgs_pvcf wes
    └────┬────┴───────┴──────────────────┘
      phenotypes
         │
        qc
"""

from pathlib import Path
import yaml

# ── Load config ───────────────────────────────────────────────────────────────

with open("config.yaml") as _f:
    _cfg = yaml.safe_load(_f)

DATA         = _cfg["data_dir"]
REF_DIR      = _cfg["reference_dir"]
TMP_DIR      = _cfg["tmp_dir"]
AUTOSOMES    = [str(c) for c in _cfg["autosomes"]]
SEX_CHROMS   = [str(c) for c in _cfg.get("sex_chromosomes", [])]
CHROMS       = AUTOSOMES + SEX_CHROMS

KG_DATE      = _cfg["kg_phase3_release_date"]
KG_VER       = _cfg["kg_phase3_vcf_version"]
KG_NYGC_DATE = _cfg["kg_nygc_date"]

KG_RAW       = f"{DATA}/raw/1kg"
SGDP_RAW     = f"{DATA}/raw/sgdp"
META         = f"{DATA}/metadata"
SHOWCASE     = f"{DATA}/Showcase"

ARRAY_DIR    = f"{DATA}/Bulk/Genotype Results/Genotype calls"
IMPUTED_DIR  = f"{DATA}/Bulk/Imputed"
WGS_DIR      = f"{DATA}/Bulk/Whole genome sequences"
WES_DIR      = f"{DATA}/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - 500k release"
WES_BGEN_DIR = f"{DATA}/Bulk/Exome sequences/Population level exome OQFE variants, BGEN format - final release"

# ── Target rules ──────────────────────────────────────────────────────────────

rule all:
    input:
        # array
        expand(f"{ARRAY_DIR}/ukb22418_c{{chrom}}_b0_v2.bed", chrom=AUTOSOMES),
        # imputed
        expand(f"{IMPUTED_DIR}/ukb22828_c{{chrom}}_b0_v3.bgen", chrom=AUTOSOMES),
        # WES
        expand(f"{WES_DIR}/ukb23157_c{{chrom}}_b0_v1.bed", chrom=AUTOSOMES),
        # WGS pVCF
        expand(f"{WGS_DIR}/ukb23370_c{{chrom}}_b0_v1.pvcf.gz", chrom=CHROMS),
        # phenotypes
        f"{SHOWCASE}/participant.parquet",
        # QC
        f"{META}/ukb_sqc_v2.txt",
        f"{META}/ukb_rel.txt",

rule all_array:
    input:
        expand(f"{ARRAY_DIR}/ukb22418_c{{chrom}}_b0_v2.bed", chrom=AUTOSOMES),

rule all_imputed:
    input:
        expand(f"{IMPUTED_DIR}/ukb22828_c{{chrom}}_b0_v3.bgen", chrom=AUTOSOMES),

rule all_wes:
    input:
        expand(f"{WES_DIR}/ukb23157_c{{chrom}}_b0_v1.bed", chrom=AUTOSOMES),

rule all_wgs:
    input:
        expand(f"{WGS_DIR}/ukb23370_c{{chrom}}_b0_v1.pvcf.gz", chrom=CHROMS),

rule all_sgdp_pvcf:
    input:
        expand(f"{SGDP_RAW}/pvcf/sgdp_c{{chrom}}.pvcf.gz", chrom=CHROMS),

# ── Download rules ────────────────────────────────────────────────────────────

rule download_phase3_vcf:
    # Autosome-only rule — sex chroms have different filename patterns (see rules below).
    wildcard_constraints:
        chrom = "|".join(str(c) for c in _cfg["autosomes"]),
    output:
        vcf   = f"{KG_RAW}/ALL.chr{{chrom}}.phase3_shapeit2_mvncall_integrated_{KG_VER}.{KG_DATE}.genotypes.vcf.gz",
        index = f"{KG_RAW}/ALL.chr{{chrom}}.phase3_shapeit2_mvncall_integrated_{KG_VER}.{KG_DATE}.genotypes.vcf.gz.tbi",
    params:
        base = _cfg["kg_phase3_base_url"],
    shell:
        """
        mkdir -p {KG_RAW}
        wget --continue -q --show-progress \
             -O {output.vcf} \
             {params.base}/ALL.chr{wildcards.chrom}.phase3_shapeit2_mvncall_integrated_{KG_VER}.{KG_DATE}.genotypes.vcf.gz
        wget --continue -q --show-progress \
             -O {output.index} \
             {params.base}/ALL.chr{wildcards.chrom}.phase3_shapeit2_mvncall_integrated_{KG_VER}.{KG_DATE}.genotypes.vcf.gz.tbi
        """

rule download_phase3_chrX_vcf:
    output:
        vcf   = f"{KG_RAW}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.{KG_DATE}.genotypes.vcf.gz",
        index = f"{KG_RAW}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.{KG_DATE}.genotypes.vcf.gz.tbi",
    params:
        base = _cfg["kg_phase3_base_url"],
    shell:
        """
        mkdir -p {KG_RAW}
        wget --continue -q --show-progress \
             -O {output.vcf} \
             {params.base}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.{KG_DATE}.genotypes.vcf.gz
        wget --continue -q --show-progress \
             -O {output.index} \
             {params.base}/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.{KG_DATE}.genotypes.vcf.gz.tbi
        """

rule download_phase3_chrY_vcf:
    output:
        vcf   = f"{KG_RAW}/ALL.chrY.phase3_integrated_v2b.{KG_DATE}.genotypes.vcf.gz",
        index = f"{KG_RAW}/ALL.chrY.phase3_integrated_v2b.{KG_DATE}.genotypes.vcf.gz.tbi",
    params:
        base = _cfg["kg_phase3_base_url"],
    shell:
        """
        mkdir -p {KG_RAW}
        wget --continue -q --show-progress \
             -O {output.vcf} \
             {params.base}/ALL.chrY.phase3_integrated_v2b.{KG_DATE}.genotypes.vcf.gz
        wget --continue -q --show-progress \
             -O {output.index} \
             {params.base}/ALL.chrY.phase3_integrated_v2b.{KG_DATE}.genotypes.vcf.gz.tbi
        """

rule download_phase3_panel:
    output:
        panel = f"{KG_RAW}/integrated_call_samples_v3.{KG_DATE}.ALL.panel",
        ped   = f"{KG_RAW}/20130606_g1k.ped",
    params:
        base = _cfg["kg_phase3_base_url"],
    shell:
        """
        mkdir -p {KG_RAW}
        wget --continue -q -O {output.panel} \
            {params.base}/integrated_call_samples_v3.{KG_DATE}.ALL.panel
        wget --continue -q -O {output.ped} \
            https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
        """

rule download_nygc_vcf:
    output:
        vcf   = f"{KG_RAW}/CCDG_14151_B01_GRM_WGS_{KG_NYGC_DATE}_chr{{chrom}}.filtered.shapeit2-duohmm-phased.vcf.gz",
        index = f"{KG_RAW}/CCDG_14151_B01_GRM_WGS_{KG_NYGC_DATE}_chr{{chrom}}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi",
    params:
        base = _cfg["kg_nygc_base_url"],
    shell:
        """
        mkdir -p {KG_RAW}
        wget --continue -q --show-progress \
             -O {output.vcf} \
             {params.base}/CCDG_14151_B01_GRM_WGS_{KG_NYGC_DATE}_chr{wildcards.chrom}.filtered.shapeit2-duohmm-phased.vcf.gz
        wget --continue -q --show-progress \
             -O {output.index} \
             {params.base}/CCDG_14151_B01_GRM_WGS_{KG_NYGC_DATE}_chr{wildcards.chrom}.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
        """

rule download_nygc_panel:
    output:
        f"{KG_RAW}/20130606_g1k_3202_samples_ped_population.txt",
    shell:
        """
        mkdir -p {KG_RAW}
        wget --continue -q -O {output} \
            https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt
        """

rule download_reference:
    output:
        exome_bed   = _cfg["exome_bed"],
        grch37_fai  = f"{REF_DIR}/GRCh37/human_g1k_v37.fasta.fai",
        grch38_fai  = f"{REF_DIR}/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai",
    shell:
        "tbooo --config config.yaml download reference"

# ── EID assignment ────────────────────────────────────────────────────────────

rule assign_eids:
    input:
        panel    = f"{KG_RAW}/20130606_g1k_3202_samples_ped_population.txt",
        sgdp_meta = f"{SGDP_RAW}/sgdp_samples.tsv",
    output:
        kg_map      = f"{META}/eid_map_1kg.tsv",
        kg_rename   = f"{META}/vcf_sample_rename_1kg.txt",
        sgdp_map    = f"{META}/eid_map_sgdp.tsv",
        sgdp_rename = f"{META}/vcf_sample_rename_sgdp.txt",
    shell:
        "tbooo --config config.yaml map eids"

# ── Array pipeline ────────────────────────────────────────────────────────────

rule array_chrom:
    input:
        vcf    = f"{KG_RAW}/ALL.chr{{chrom}}.phase3_shapeit2_mvncall_integrated_{KG_VER}.{KG_DATE}.genotypes.vcf.gz",
        index  = f"{KG_RAW}/ALL.chr{{chrom}}.phase3_shapeit2_mvncall_integrated_{KG_VER}.{KG_DATE}.genotypes.vcf.gz.tbi",
        eid_map = f"{META}/eid_map_1kg.tsv",
    output:
        bed = f"{ARRAY_DIR}/ukb22418_c{{chrom}}_b0_v2.bed",
        bim = f"{ARRAY_DIR}/ukb22418_c{{chrom}}_b0_v2.bim",
        fam = f"{ARRAY_DIR}/ukb22418_c{{chrom}}_b0_v2.fam",
    shell:
        "tbooo --config config.yaml map array --chroms {wildcards.chrom}"

# ── Imputed pipeline ──────────────────────────────────────────────────────────

rule imputed_chrom:
    input:
        vcf     = f"{KG_RAW}/ALL.chr{{chrom}}.phase3_shapeit2_mvncall_integrated_{KG_VER}.{KG_DATE}.genotypes.vcf.gz",
        index   = f"{KG_RAW}/ALL.chr{{chrom}}.phase3_shapeit2_mvncall_integrated_{KG_VER}.{KG_DATE}.genotypes.vcf.gz.tbi",
        rename  = f"{META}/vcf_sample_rename_1kg.txt",
    output:
        bgen   = f"{IMPUTED_DIR}/ukb22828_c{{chrom}}_b0_v3.bgen",
        bgi    = f"{IMPUTED_DIR}/ukb22828_c{{chrom}}_b0_v3.bgen.bgi",
    shell:
        "tbooo --config config.yaml map imputed --chroms {wildcards.chrom}"

# ── WGS pVCF ─────────────────────────────────────────────────────────────────

rule wgs_pvcf_chrom:
    input:
        vcf    = f"{KG_RAW}/CCDG_14151_B01_GRM_WGS_{KG_NYGC_DATE}_chr{{chrom}}.filtered.shapeit2-duohmm-phased.vcf.gz",
        rename = f"{META}/vcf_sample_rename_1kg.txt",
    output:
        pvcf  = f"{WGS_DIR}/ukb23370_c{{chrom}}_b0_v1.pvcf.gz",
        index = f"{WGS_DIR}/ukb23370_c{{chrom}}_b0_v1.pvcf.gz.tbi",
    shell:
        "tbooo --config config.yaml map wgs --no-croms --chroms {wildcards.chrom}"

rule wgs_crams:
    input:
        kg_map   = f"{META}/eid_map_1kg.tsv",
    output:
        flag = f"{META}/.crams_linked",
    shell:
        """
        tbooo --config config.yaml map wgs --no-pvcf
        touch {output.flag}
        """

rule download_sgdp:
    output:
        metadata = f"{SGDP_RAW}/sgdp_samples.tsv",
        vcf_flag = f"{SGDP_RAW}/vcf/.downloaded",
    shell:
        """
        tbooo --config config.yaml download sgdp
        touch {output.vcf_flag}
        """

rule sgdp_gvcfs:
    input:
        vcf_flag = f"{SGDP_RAW}/vcf/.downloaded",
        eid_map  = f"{META}/eid_map_sgdp.tsv",
    output:
        flag = f"{META}/.sgdp_gvcfs_linked",
    shell:
        """
        tbooo --config config.yaml map wgs --no-croms --no-pvcf --sgdp-gvcf
        touch {output.flag}
        """

rule sgdp_pvcf_chrom:
    input:
        vcf_flag  = f"{SGDP_RAW}/vcf/.downloaded",
        rename    = f"{META}/vcf_sample_rename_sgdp.txt",
    output:
        pvcf  = f"{SGDP_RAW}/pvcf/sgdp_c{{chrom}}.pvcf.gz",
        index = f"{SGDP_RAW}/pvcf/sgdp_c{{chrom}}.pvcf.gz.tbi",
    shell:
        "tbooo --config config.yaml map wgs --no-croms --no-pvcf --sgdp-pvcf --chroms {wildcards.chrom}"

# ── WES pipeline ─────────────────────────────────────────────────────────────

rule wes_chrom:
    input:
        vcf      = f"{KG_RAW}/CCDG_14151_B01_GRM_WGS_{KG_NYGC_DATE}_chr{{chrom}}.filtered.shapeit2-duohmm-phased.vcf.gz",
        exome_bed = _cfg["exome_bed"],
        rename   = f"{META}/vcf_sample_rename_1kg.txt",
    output:
        bed = f"{WES_DIR}/ukb23157_c{{chrom}}_b0_v1.bed",
        bim = f"{WES_DIR}/ukb23157_c{{chrom}}_b0_v1.bim",
        fam = f"{WES_DIR}/ukb23157_c{{chrom}}_b0_v1.fam",
    shell:
        "tbooo --config config.yaml map wes --chroms {wildcards.chrom}"

# ── Phenotypes ────────────────────────────────────────────────────────────────

rule phenotypes:
    input:
        kg_map = f"{META}/eid_map_1kg.tsv",
    output:
        parquet = f"{SHOWCASE}/participant.parquet",
    shell:
        "tbooo --config config.yaml map phenotypes"

# ── QC files ─────────────────────────────────────────────────────────────────

rule qc:
    input:
        expand(f"{ARRAY_DIR}/ukb22418_c{{chrom}}_b0_v2.bed", chrom=AUTOSOMES),
        kg_map = f"{META}/eid_map_1kg.tsv",
    output:
        sqc = f"{META}/ukb_sqc_v2.txt",
        rel = f"{META}/ukb_rel.txt",
    shell:
        "tbooo --config config.yaml map qc"
