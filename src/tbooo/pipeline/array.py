"""Build UKB-mirrored PLINK genotyping array files (Field 22418).

For each chromosome:
  1. Filter 1KGP Phase 3 VCF to array positions
     - If array_manifest is provided: use it as a sites filter
     - Otherwise: keep biallelic SNPs with MAF >= array_proxy_maf and an rsID
  2. Decompose multi-allelic sites, normalize
  3. Convert VCF → PLINK BED/BIM/FAM via plink2
  4. Rewrite FAM: replace original sample IDs with EIDs, add batch codes
  5. Add cM positions to BIM from genetic maps

Outputs:
    data/Bulk/Genotype Results/Genotype calls/ukb22418_c{chrom}_b0_v2.{bed,bim,fam}
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from tbooo.config import Config
from tbooo.integrity import ensure_plink, remove, vcf_ok
from tbooo.utils import bgzip_test, ensure_dirs, log, run

_CM_MAP_GLOB = "genetic_map_GRCh37_chr{chrom}.txt"


def run_array_pipeline(cfg: Config, chroms: list[str]) -> None:
    ensure_dirs(cfg.array_dir())
    eid_map = _load_eid_map(cfg)

    for chrom in chroms:
        log(f"[array] chromosome {chrom}")
        vcf = cfg.phase3_vcf(chrom)
        if not vcf.exists():
            log(f"  SKIP: Phase 3 VCF not found: {vcf}")
            continue
        _process_chrom(cfg, chrom, vcf, eid_map)

    log("Array pipeline complete.")


def _process_chrom(cfg: Config, chrom: str, vcf: Path, eid_map: pd.DataFrame) -> None:
    stem = cfg.array_stem(chrom)

    def _build() -> None:
        tmp = cfg.tmp_dir / f"array_chr{chrom}_filtered.vcf.gz"
        ensure_dirs(cfg.tmp_dir)

        # Step 1: filter to array sites
        _filter_to_array_sites(cfg, chrom, vcf, tmp)

        # Step 2: VCF → PLINK (plink2 handles MAF + biallelic decomposition)
        plink_cmd = [
            cfg.tools.plink2,
            "--vcf", str(tmp),
            "--max-alleles", "2",
            "--snps-only", "just-acgt",
            "--make-bed",
            "--out", str(stem),
            "--chr", chrom,
            "--no-psam-pheno",
        ]
        # Apply MAF cutoff at the plink2 stage (only when running the proxy filter —
        # an explicit manifest already restricts to chosen sites, no MAF gate there).
        if not (cfg.array_manifest and Path(cfg.array_manifest).exists()):
            plink_cmd += ["--maf", str(cfg.array_proxy_maf)]
        run(plink_cmd)

        # Step 3: rewrite FAM with EIDs and batch codes
        _rewrite_fam(stem.with_suffix(".fam"), eid_map)

        # Step 4: inject cM positions into BIM
        _inject_cm(cfg, chrom, stem.with_suffix(".bim"))

        tmp.unlink(missing_ok=True)
        log(f"  done → {stem}.bed/bim/fam")

    ensure_plink(f"{stem.name} (array chr{chrom})", stem, _build)


def _filter_to_array_sites(cfg: Config, chrom: str, vcf: Path, out: Path) -> None:
    manifest = cfg.array_manifest
    bcftools = cfg.tools.bcftools
    maf = cfg.array_proxy_maf

    if vcf_ok(out):
        log(f"  reusing filtered VCF (chr{chrom})")
        return
    remove(out, Path(str(out) + ".tbi"))  # clear any partial filter from a prior crash

    if manifest and Path(manifest).exists():
        # Use the manifest as a regions filter
        log(f"  filtering to array manifest sites (chr{chrom})…")
        regions_file = _manifest_to_regions(Path(manifest), chrom, cfg.tmp_dir)
        run([
            bcftools, "view",
            "--regions-file", str(regions_file),
            "--output-type", "z",
            "--output", str(out),
            str(vcf),
        ])
    else:
        # Proxy filter: biallelic SNPs only. MAF is applied later by plink2 (its
        # native allele-frequency machinery is more reliable than bcftools'
        # MAF/--min-af expressions on Phase 3 INFO, which silently yield 0 matches).
        log(f"  no manifest — proxy filter: biallelic SNPs (chr{chrom}); MAF>={maf} applied at plink2 stage…")
        run([
            bcftools, "view",
            "--min-alleles", "2", "--max-alleles", "2",
            "--types", "snps",
            "--output-type", "z",
            "--output", str(out),
            str(vcf),
        ])
    run([bcftools, "index", "--tbi", str(out)])

    # Sanity-check: a 0-variant filtered VCF crashes plink2 with a cryptic error.
    # Detect it now and diagnose whether the input VCF is corrupt.
    n = _count_variants(bcftools, out)
    if n == 0:
        _diagnose_empty_filter(cfg, chrom, vcf)
    log(f"  filtered chr{chrom}: {n:,} variants")


def _count_variants(bcftools: str, vcf: Path) -> int:
    proc = run([bcftools, "index", "--stats", str(vcf)], capture=True)
    total = 0
    for line in proc.stdout.decode().splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[-1].isdigit():
            total += int(parts[-1])
    return total


def _diagnose_empty_filter(cfg: Config, chrom: str, vcf: Path) -> None:
    """Raise with actionable guidance when the filtered VCF has 0 variants."""
    check = bgzip_test(vcf)
    if check is False:
        raise RuntimeError(
            f"chr{chrom}: input 1KGP VCF appears corrupt (bgzip integrity check failed).\n"
            f"  Path: {vcf}\n"
            f"  Fix:  rm '{vcf}' '{vcf}.tbi' && "
            f"tbooo download 1kg --chroms {chrom} --deep-check"
        )
    if check is None:
        hint = "Install `bgzip` to enable input integrity verification."
    else:
        hint = (
            f"Input VCF is intact — filter yielded 0 matches. "
            f"Check `array_manifest` or `array_proxy_maf={cfg.array_proxy_maf}` threshold."
        )
    raise RuntimeError(
        f"chr{chrom}: filter produced 0 variants from {vcf.name}. {hint}"
    )


def _manifest_to_regions(manifest: Path, chrom: str, tmp_dir: Path) -> Path:
    """Convert array manifest CSV to a bcftools regions BED file."""
    regions_path = tmp_dir / f"array_regions_chr{chrom}.bed"
    if regions_path.exists():
        return regions_path

    # Axiom manifest has columns: Chromosome, Position, ...
    df = pd.read_csv(manifest, comment="#", low_memory=False)
    col_chr = next((c for c in df.columns if "chrom" in c.lower()), None)
    col_pos = next((c for c in df.columns if "position" in c.lower() or c.lower() == "pos"), None)
    if not col_chr or not col_pos:
        raise ValueError(f"Cannot find Chromosome/Position columns in {manifest}")

    sub = df[df[col_chr].astype(str) == chrom][[col_chr, col_pos]].dropna()
    sub[col_pos] = sub[col_pos].astype(int)
    bed = sub.rename(columns={col_chr: "#CHROM", col_pos: "POS"})
    bed["START"] = bed["POS"] - 1
    bed = bed[["#CHROM", "START", "POS"]]
    bed.to_csv(regions_path, sep="\t", index=False)
    return regions_path


def _rewrite_fam(fam_path: Path, eid_map: pd.DataFrame) -> None:
    """Replace sample IDs with EIDs and set batch codes in FAM column 6."""
    fam = pd.read_csv(fam_path, sep=r"\s+", header=None,
                      names=["fid", "iid", "pat", "mat", "sex", "pheno"])
    lookup = eid_map.set_index("sample_id")[["eid", "sex", "batch"]]
    rows = []
    for _, r in fam.iterrows():
        if r["iid"] in lookup.index:
            info = lookup.loc[r["iid"]]
            eid = info["eid"]
            sex = info["sex"] if info["sex"] != 0 else r["sex"]
            batch = info["batch"]
        else:
            eid = r["iid"]
            sex = r["sex"]
            batch = 0
        rows.append({"fid": 0, "iid": eid, "pat": 0, "mat": 0, "sex": sex, "pheno": batch})
    pd.DataFrame(rows).to_csv(fam_path, sep="\t", index=False, header=False)


def _inject_cm(cfg: Config, chrom: str, bim_path: Path) -> None:
    """Fill BIM column 3 (cM) from 1KGP genetic maps."""
    map_glob = _CM_MAP_GLOB.format(chrom=chrom)
    maps_dir = cfg.reference_dir / "genetic_maps"
    map_files = list(maps_dir.glob(f"**/*chr{chrom}*.txt"))
    if not map_files:
        log(f"  WARNING: no genetic map found for chr{chrom}; cM column left as 0")
        return

    gmap = pd.read_csv(map_files[0], sep=r"\s+", header=0,
                       names=["pos", "rate", "cm"])
    gmap = gmap.sort_values("pos")
    positions = gmap["pos"].to_numpy()
    cms = gmap["cm"].to_numpy()

    bim = pd.read_csv(bim_path, sep="\t", header=None,
                      names=["chr", "snp", "cm", "bp", "a1", "a2"])
    import numpy as np
    bim["cm"] = np.interp(bim["bp"].to_numpy(), positions, cms)
    bim.to_csv(bim_path, sep="\t", index=False, header=False)


def _load_eid_map(cfg: Config) -> pd.DataFrame:
    path = cfg.metadata_dir() / "eid_map_1kg.tsv"
    if not path.exists():
        raise FileNotFoundError(f"EID map not found: {path}\nRun `tbooo map eids` first.")
    df = pd.read_csv(path, sep="\t")
    if "batch" not in df.columns:
        from tbooo.pipeline.eids import _BATCH_MAP
        df["batch"] = df["super_pop"].map(_BATCH_MAP).fillna(0).astype(int)
    return df
