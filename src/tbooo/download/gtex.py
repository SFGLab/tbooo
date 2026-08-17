"""Download GTEx per-tissue expression, covariates and annotations.

Individual GTEx genotypes are dbGaP-protected, which would normally rule out PCCA
because the method needs genotype PCs as its CCA anchor. However
`GTEx_Analysis_v10_eQTL_covariates.tar` is open access and carries PC1-PC5 genotype
PCs, identical across all tissues for a given donor. Those PCs are the anchor.

Per-tissue raw TPM is used in preference to the eQTL expression matrices, which are
TMM-normalized then inverse-normal transformed. Raw TPM matches the GEUVADIS
treatment (raw RPKM plus our own normalization), so all cohorts traverse one code
path, and it is a quarter of the download.

Outputs (data/raw/gtex/):
    gene_tpm_v10_<tissue>.gct.gz             raw TPM, one file per tissue
    GTEx_Analysis_v10_eQTL_covariates.tar    genotype PCs (PC1-PC5) per tissue
    GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt   centre, batch, tissue
    GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt  sex, age
"""

from __future__ import annotations

from tbooo.config import Config
from tbooo.utils import ensure_dirs, log, wget_download

_BUCKET = "https://storage.googleapis.com/adult-gtex"

# Chosen for donor count and to test whether the result holds across tissues
# rather than being an artefact of one. Expression PCA across all tissues at once
# would be dominated by tissue identity, not ancestry.
TISSUES = ["whole_blood", "muscle_skeletal", "lung", "thyroid"]

COVARIATES_FILE = "GTEx_Analysis_v10_eQTL_covariates.tar"
SAMPLE_ATTRIBUTES_FILE = "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
SUBJECT_PHENOTYPES_FILE = "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"


def tissue_tpm_filename(tissue: str) -> str:
    return f"gene_tpm_v10_{tissue}.gct.gz"


def _urls() -> dict[str, str]:
    urls = {
        COVARIATES_FILE: (
            f"{_BUCKET}/bulk-qtl/v10/single-tissue-cis-qtl/{COVARIATES_FILE}"
        ),
        SAMPLE_ATTRIBUTES_FILE: (
            f"{_BUCKET}/annotations/v10/metadata-files/{SAMPLE_ATTRIBUTES_FILE}"
        ),
        SUBJECT_PHENOTYPES_FILE: (
            f"{_BUCKET}/annotations/v10/metadata-files/{SUBJECT_PHENOTYPES_FILE}"
        ),
    }
    for tissue in TISSUES:
        name = tissue_tpm_filename(tissue)
        urls[name] = f"{_BUCKET}/bulk-gex/v10/rna-seq/tpms-by-tissue/{name}"
    return urls


def download_gtex(cfg: Config) -> None:
    """Fetch GTEx expression for the selected tissues plus annotations."""
    out_dir = cfg.gtex_raw_dir()
    ensure_dirs(out_dir)

    for name, url in _urls().items():
        dest = out_dir / name
        if dest.exists():
            log(f"  skip (exists): {name}")
            continue
        log(f"Downloading {name}…")
        wget_download(url, dest, tool_wget=cfg.tools.wget)
        log(f"  wrote {dest}")
