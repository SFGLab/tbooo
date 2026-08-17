"""Canonical data locations and source URLs for the PCCA replication."""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

DATA = REPO_ROOT / "data"
GEUVADIS_RAW = DATA / "raw" / "geuvadis"
METADATA = DATA / "metadata"
FIGURES = METADATA / "eda"
GENOTYPE = GEUVADIS_RAW / "genotype"

MAGE_RAW = DATA / "raw" / "mage"
MAGE_EXPRESSION = MAGE_RAW / "TMM.filtered.TSS.MAGE.v1.0.bed.gz"
MAGE_COVARIATES = MAGE_RAW / "eQTL_covariates.tab.gz"
MAGE_SAMPLE_METADATA = MAGE_RAW / "sample.metadata.MAGE.v1.0.txt"
JOINT_PCS = MAGE_RAW / "MAGE_Geuvadis_GTEx_AFGR.top_PCs.txt.gz"

GTEX_RAW = DATA / "raw" / "gtex"
GTEX_COVARIATES_TAR = GTEX_RAW / "GTEx_Analysis_v10_eQTL_covariates.tar"
GTEX_SAMPLE_ATTRIBUTES = GTEX_RAW / "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
GTEX_SUBJECT_PHENOTYPES = GTEX_RAW / "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
GTEX_TISSUES = ["whole_blood", "muscle_skeletal", "lung", "thyroid"]


def gtex_tpm(tissue: str):
    return GTEX_RAW / f"gene_tpm_v10_{tissue}.gct.gz"

SDRF = GEUVADIS_RAW / "E-GEUV-1.sdrf.txt"
GD660_GENE_RPKM = GEUVADIS_RAW / "GD660.GeneQuantRPKM.txt.gz"
GD660_TR_RPKM = GEUVADIS_RAW / "GD660.TrQuantRPKM.txt.gz"
GENCODE_GTF = DATA / "reference" / "gencode.v47.annotation.gtf.gz"

_EBI = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1"
SDRF_URL = f"{_EBI}/E-GEUV-1.sdrf.txt"
GD660_GENE_RPKM_URL = f"{_EBI}/analysis_results/GD660.GeneQuantRPKM.txt.gz"
GD660_TR_RPKM_URL = f"{_EBI}/analysis_results/GD660.TrQuantRPKM.txt.gz"

# GEUVADIS also hosts the paper's own Phase 1 genotypes under E-GEUV-1/genotypes/,
# but they total ~50 GB and are deliberately not used here. See the Deviations
# section of the spec: Phase 3 is the genotype source for this replication.

# plink2 prebuilt 1000 Genomes Phase 3 (build 37, 2504 samples). The per-chromosome
# files are hosted on Dropbox behind opaque hashes that change between releases, so
# the URLs are scraped from the resources page at runtime rather than hardcoded.
# Files are named chr{N}_phase3.{pgen,pvar}.zst plus phase3_orig.psam.
PLINK2_RESOURCES_URL = "https://www.cog-genomics.org/plink/2.0/resources"
PLINK2_MAC_ARM64_URL = (
    "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_mac_arm64_20260504.zip"
)
