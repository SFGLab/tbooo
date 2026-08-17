"""plink2 driver for the genotype side of PCCA.

Everything works one chromosome at a time: the machine has ~22 GiB free on a
98%-full volume, and a whole-cohort decompressed callset does not fit. The .pvar
is left zstd-compressed and read via plink2's `vzs` modifier, which saves several
GB more.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from pcca import paths

PLINK2 = str(Path("~/bin/plink2").expanduser())
AUTOSOMES = list(range(1, 23))


def fetch_phase3_urls() -> dict[str, str]:
    """Scrape the plink2 resources page for build-37 Phase 3 download links.

    The files live on Dropbox behind opaque hashes that change between releases,
    so hardcoding them would rot silently.
    """
    with urllib.request.urlopen(paths.PLINK2_RESOURCES_URL, timeout=180) as response:
        html = response.read().decode("utf-8", errors="replace")
    urls: dict[str, str] = {}
    for match in re.finditer(r"https://www\.dropbox\.com/s/[^\"'>]+", html):
        url = match.group(0)
        filename = url.split("/")[-1].split("?")[0]
        urls.setdefault(filename, url)
    return urls


def write_keep_file(cohort: pd.DataFrame, path: Path) -> Path:
    """Write a plink2 --keep file listing the cohort's 1KGP sample IDs.

    The `#IID` header makes plink2 read a single-column file as sample IDs,
    which avoids assuming anything about whether the .psam carries an FID.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#IID\n" + "\n".join(cohort["sample_id"]) + "\n")
    return path


def restrict_to_panel(cohort: pd.DataFrame, psam_path: Path) -> pd.DataFrame:
    """Drop cohort members absent from the genotype panel.

    The paper used 1000 Genomes Phase 1. Phase 3's build-37 release carries only
    the 2504 unrelated samples, and 15 of the paper's 343 individuals are not
    among them. Silently letting plink2 drop them would leave the expression
    matrix with 343 rows and the genotype matrix with 328, misaligning every
    downstream row. This makes the intersection explicit and auditable.
    """
    psam = pd.read_csv(psam_path, sep="\t")
    psam.columns = [c.lstrip("#") for c in psam.columns]
    panel = set(psam["IID"].astype(str))
    return cohort[cohort["sample_id"].isin(panel)].reset_index(drop=True)


def _run(args: list) -> None:
    subprocess.run([str(a) for a in args], check=True)


def download_phase3_chromosome(chrom: int, urls: dict[str, str], out_dir: Path) -> Path:
    """Fetch one chromosome's pgen/pvar and the shared psam.

    Only the .pgen is decompressed; the .pvar stays zstd-compressed and is read
    with plink2's `vzs` modifier.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    psam = out_dir / "phase3_orig.psam"
    if not psam.exists():
        urllib.request.urlretrieve(urls["phase3_orig.psam"], psam)

    prefix = out_dir / f"chr{chrom}_phase3"

    pvar_zst = prefix.with_suffix(".pvar.zst")
    if not pvar_zst.exists():
        urllib.request.urlretrieve(urls[f"chr{chrom}_phase3.pvar.zst"], pvar_zst)

    pgen = prefix.with_suffix(".pgen")
    if not pgen.exists():
        pgen_zst = prefix.with_suffix(".pgen.zst")
        if not pgen_zst.exists():
            urllib.request.urlretrieve(urls[f"chr{chrom}_phase3.pgen.zst"], pgen_zst)
        _run([PLINK2, "--zst-decompress", pgen_zst, pgen])
        pgen_zst.unlink()

    shutil.copy(psam, prefix.with_suffix(".psam"))
    return prefix


def filter_chromosome(
    chrom: int, keep_file: Path, work_dir: Path, out_dir: Path, urls: dict[str, str]
) -> Path:
    """Subset one chromosome to the cohort and to common biallelic autosomal SNPs.

    Deletes the downloaded source immediately afterwards so the next chromosome
    has room. No LD pruning: the paper does not prune, and pruning would change
    the PC structure being replicated.
    """
    out_prefix = out_dir / f"geuv_chr{chrom}"
    if out_prefix.with_suffix(".pgen").exists():
        return out_prefix

    src = download_phase3_chromosome(chrom, urls, work_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _run(
        [
            PLINK2, "--pfile", src, "vzs",
            "--keep", keep_file,
            "--snps-only", "--max-alleles", "2", "--min-alleles", "2",
            "--maf", "0.05", "--autosome",
            "--set-all-var-ids", "@:#:$r:$a",
            "--new-id-max-allele-len", "100", "truncate",
            "--make-pgen", "--out", out_prefix,
        ]
    )
    for suffix in (".pgen", ".pvar.zst", ".psam"):
        Path(str(src) + suffix).unlink(missing_ok=True)
    return out_prefix


def merge_and_pca(
    out_dir: Path, prefix_glob: str = "geuv_chr*", n_pcs: int = 5
) -> Path:
    """Merge per-chromosome sets and compute genotype PCs."""
    prefixes = sorted(
        (p for p in out_dir.glob(f"{prefix_glob}.pgen")),
        key=lambda p: int(re.search(r"chr(\d+)", p.name).group(1)),
    )
    merge_list = out_dir / "merge_list.txt"
    merge_list.write_text(
        "\n".join(str(p)[: -len(".pgen")] for p in prefixes) + "\n"
    )

    merged = out_dir / "geuv_merged"
    # The trailing `pfile` tells plink2 the list holds pgen/pvar/psam prefixes.
    _run([PLINK2, "--pmerge-list", merge_list, "pfile", "--make-pgen", "--out", merged])
    _run([PLINK2, "--pfile", merged, "--pca", str(n_pcs), "--out", merged])
    return Path(str(merged) + ".eigenvec")


def load_eigenvec(path: Path, cohort: pd.DataFrame) -> np.ndarray:
    """Read plink2 PCs, reordered to match the cohort's row order.

    plink2 emits samples in .psam order, which is not the expression matrix order.
    Misaligned rows would produce a plausible-looking but meaningless CCA.
    """
    df = pd.read_csv(path, sep="\t")
    # plink2 writes a leading '#' on the header and may emit either "#IID" or
    # "#FID<TAB>IID" depending on the input .psam, so locate IID by name.
    df.columns = [c.lstrip("#") for c in df.columns]
    id_col = "IID" if "IID" in df.columns else df.columns[0]
    df[id_col] = df[id_col].astype(str)
    indexed = df.set_index(id_col)
    missing = [s for s in cohort["sample_id"] if s not in indexed.index]
    if missing:
        raise KeyError(
            f"{len(missing)} cohort samples absent from {path.name}: {missing[:5]}"
        )
    pc_cols = [c for c in indexed.columns if c.upper().startswith("PC")]
    return indexed.loc[cohort["sample_id"], pc_cols].to_numpy(dtype=float)


TRAW_META_COLUMNS = ["CHR", "SNP", "(C)M", "POS", "COUNTED", "ALT"]


def export_dosages(prefix: Path, out_prefix: Path) -> Path:
    """Export additive dosages for one chromosome as a plink2 .traw file.

    `A-transpose` puts variants in rows and samples in columns. The row-major
    `A` format would give a table ~500,000 columns wide, which pandas parses
    pathologically slowly — minutes per chromosome versus seconds here.
    """
    _run([PLINK2, "--pfile", prefix, "--export", "A-transpose", "--out", out_prefix])
    return Path(str(out_prefix) + ".traw")


def load_dosages(traw_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Parse a plink2 .traw into (variant_ids, positions, dosages, sample_ids).

    Dosages are returned samples x variants to match the expression matrix
    orientation. Missing calls are mean-imputed, which is standard for eQTL scans
    and keeps the regression well-defined.
    """
    df = pd.read_csv(traw_path, sep="\t")
    sample_cols = [c for c in df.columns if c not in TRAW_META_COLUMNS]
    ids = df["SNP"].to_numpy(dtype=object)
    positions = df["POS"].to_numpy(dtype=np.int64)

    mat = df[sample_cols].to_numpy(dtype=np.float32).T  # samples x variants
    if np.isnan(mat).any():
        col_means = np.nanmean(mat, axis=0)
        nan_rows, nan_cols = np.where(np.isnan(mat))
        mat[nan_rows, nan_cols] = np.take(col_means, nan_cols)

    # plink2 writes sample columns as "FID_IID"; the cohort keys on IID alone.
    sample_ids = [c.split("_", 1)[-1] for c in sample_cols]
    return ids, positions, mat, sample_ids
