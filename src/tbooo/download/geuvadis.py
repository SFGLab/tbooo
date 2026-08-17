"""Download GEUVADIS gene expression data and sample metadata.

Downloads:
  - GD660 raw gene-level RPKM (660 runs, incl. technical replicates). This is the
    input to all population-structure work: it is unnormalized, so it retains the
    structure that the PEER-residualized GD462 file removes.
  - GD462 PEER-residualized RPKM. Retained because it is the correct input for
    eQTL work that wants PEER correction; no pipeline stage reads it.
  - E-GEUV-1 SDRF: the only source of the 7-lab assignment, plus sex, ancestry
    and the Phase 1 genotype availability flag.

Outputs:
    data/raw/geuvadis/GD660.GeneQuantRPKM.txt.gz
    data/raw/geuvadis/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz
    data/raw/geuvadis/E-GEUV-1.sdrf.txt
"""

from __future__ import annotations

from tbooo.config import Config
from tbooo.utils import ensure_dirs, log, wget_download

_EBI = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1"

_FILES = {
    "GD660.GeneQuantRPKM.txt.gz": f"{_EBI}/analysis_results/GD660.GeneQuantRPKM.txt.gz",
    "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz": (
        f"{_EBI}/analysis_results/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"
    ),
    "E-GEUV-1.sdrf.txt": f"{_EBI}/E-GEUV-1.sdrf.txt",
}


def download_expression(cfg: Config) -> None:
    """Download the GEUVADIS quantification matrices and sample metadata."""
    out_dir = cfg.geuvadis_raw_dir()
    ensure_dirs(out_dir)
    for name, url in _FILES.items():
        dest = out_dir / name
        if dest.exists():
            log(f"  skip (exists): {name}")
            continue
        log(f"Downloading {name}…")
        wget_download(url, dest, tool_wget=cfg.tools.wget)
        log(f"  wrote {dest}")
