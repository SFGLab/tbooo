"""Download GEUVADIS gene expression data.

Downloads the GD462 gene-level RPKM matrix (462 samples × ~23 k genes)
from the GEUVADIS FTP. All samples are 1KGP individuals, enabling direct
join with eid_map_1kg.tsv downstream.

Outputs:
    data/raw/geuvadis/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz
"""

from __future__ import annotations

from tbooo.config import Config
from tbooo.utils import ensure_dirs, log, wget_download

_RPKM_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/"
    "GEUV/E-GEUV-1/analysis_results/"
    "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"
)
_RPKM_FILENAME = "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"


def download_expression(cfg: Config) -> None:
    """Download the GD462 RPKM matrix from the GEUVADIS FTP."""
    out_dir = cfg.geuvadis_raw_dir()
    ensure_dirs(out_dir)
    dest = out_dir / _RPKM_FILENAME
    if dest.exists():
        log(f"  skip (exists): {dest.name}")
        return
    log("Downloading GEUVADIS GD462 RPKM matrix…")
    wget_download(_RPKM_URL, dest, tool_wget=cfg.tools.wget)
    log(f"  wrote {dest}")
