"""Download MAGE expression, metadata and genotype PCs.

MAGE (Taylor et al., Nature 2024) is RNA-seq from lymphoblastoid cell lines of 731
1000 Genomes individuals spanning 26 populations and five continental groups. It
is the ancestry-diverse counterpart to GEUVADIS, whose four populations included
only one non-European group.

The processed data ships as a single 59.5 GB Zip64 archive on Zenodo. The members
needed here total ~165 MB, so they are pulled selectively over HTTP range requests
rather than downloading the archive. See tbooo.download.remote_zip.

Outputs (data/raw/mage/):
    TMM.filtered.TSS.MAGE.v1.0.bed.gz        gene expression, BED4+ with TSS coords
    eQTL_covariates.tab.gz                   sex, 5 genotype PCs, 60 PEER factors
    sample.metadata.MAGE.v1.0.txt            population, continental group, batch, sex
    MAGE_Geuvadis_GTEx_AFGR.top_PCs.txt.gz   joint genotype PCs across four cohorts
"""

from __future__ import annotations

from tbooo.config import Config
from tbooo.download.remote_zip import RemoteZip
from tbooo.utils import ensure_dirs, log

ZENODO_URL = (
    "https://zenodo.org/records/10535719/files/MAGE.v1.0.data.zip?download=1"
)

_ROOT = "MAGE.v1.0.data"

# Archive member -> local filename.
MEMBERS = {
    f"{_ROOT}/QTL_results/eQTL_results/expression_quants/"
    "TMM.filtered.TSS.MAGE.v1.0.bed.gz": "TMM.filtered.TSS.MAGE.v1.0.bed.gz",
    f"{_ROOT}/QTL_results/eQTL_results/expression_quants/"
    "eQTL_covariates.tab.gz": "eQTL_covariates.tab.gz",
    f"{_ROOT}/sample_library_info/"
    "sample.metadata.MAGE.v1.0.txt": "sample.metadata.MAGE.v1.0.txt",
    f"{_ROOT}/dataset_comparison/"
    "MAGE_Geuvadis_GTEx_AFGR.top_PCs.txt.gz": "MAGE_Geuvadis_GTEx_AFGR.top_PCs.txt.gz",
}


def download_mage(cfg: Config) -> None:
    """Fetch the MAGE members needed for population-structure analysis."""
    out_dir = cfg.mage_raw_dir()
    ensure_dirs(out_dir)

    missing = {m: name for m, name in MEMBERS.items() if not (out_dir / name).exists()}
    if not missing:
        log("  skip (exists): all MAGE files present")
        return

    log(f"Opening remote MAGE archive ({len(missing)} member(s) to fetch)…")
    archive = RemoteZip(ZENODO_URL)
    available = archive.entries
    log(f"  archive has {len(available)} members; reading selectively")

    for member, name in missing.items():
        if member not in available:
            raise KeyError(
                f"MAGE archive member not found: {member}\n"
                "The Zenodo release layout may have changed."
            )
        size = available[member].uncompressed_size
        log(f"  fetching {name} ({size / 1024 / 1024:.1f} MB)…")
        archive.extract(member, out_dir / name)

    log(f"  wrote {out_dir}")
