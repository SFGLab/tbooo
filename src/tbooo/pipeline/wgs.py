"""Build UKB-mirrored WGS data (Fields 23149, 23151, 23370).

Individual CRAMs (Field 23149):
  - Rename/symlink 1KGP NYGC 30x CRAMs → <EID>_23149_0_0.cram

Per-sample gVCFs (Field 23151) — both cohorts:
  - 1KGP: extract per-sample VCF from NYGC 30x cohort VCFs → <EID>_23151_0_0.g.vcf.gz
  - SGDP: symlink downloaded per-sample VCFs        → <EID>_23151_0_0.g.vcf.gz

Cohort pVCF (Field 23370) — both cohorts merged:
  - 1KGP: reheader NYGC per-chromosome VCF with EIDs → raw/1kg/pvcf/nygc_c{chrom}.pvcf.gz
  - SGDP: merge per-sample VCFs by chromosome        → raw/sgdp/pvcf/sgdp_c{chrom}.pvcf.gz
  - Combined: bcftools merge nygc + sgdp             → Bulk/Whole genome sequences/ukb23370_c{chrom}_b0_v1.pvcf.gz
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from tbooo.config import Config
from tbooo.utils import ensure_dirs, eid_prefix_dir, log, run


# ── Public API ────────────────────────────────────────────────────────────────

def rename_crams(cfg: Config) -> None:
    """Symlink 1KGP NYGC 30x CRAMs into the UKB-mirrored individual-file structure."""
    ensure_dirs(cfg.wgs_dir())
    kg_map = _load_eid_map(cfg, "eid_map_1kg.tsv")
    _rename_kg_crams(cfg, kg_map)
    log("CRAM renaming complete.")


def build_gvcfs(cfg: Config, chroms: list[str]) -> None:
    """Build per-sample gVCFs (Field 23151) from both NYGC and SGDP.

    1KGP: extracts each sample from NYGC 30x cohort VCFs and concatenates chromosomes.
    SGDP: symlinks downloaded per-sample VCFs directly.
    """
    _build_nygc_gvcfs(cfg, chroms)
    _symlink_sgdp_gvcfs(cfg)
    log("gVCF build complete.")


def build_pvcf(cfg: Config, chroms: list[str]) -> None:
    """Build cohort pVCF (Field 23370) by merging NYGC and SGDP per-chromosome VCFs.

    Produces two intermediates (nygc_pvcf, sgdp_pvcf) then merges into the
    final UKB-mirrored output.  Either intermediate is used alone if the other
    is absent (e.g. SGDP VCFs not yet downloaded).
    """
    _build_nygc_pvcf(cfg, chroms)
    _build_sgdp_pvcf(cfg, chroms)
    _merge_pvcfs(cfg, chroms)
    log("pVCF build complete.")


# ── gVCF internals ────────────────────────────────────────────────────────────

def _build_nygc_gvcfs(cfg: Config, chroms: list[str]) -> None:
    kg_map = _load_eid_map(cfg, "eid_map_1kg.tsv")
    if kg_map.empty:
        log("  WARNING: 1KGP EID map not found; skipping NYGC gVCF extraction.")
        return

    avail = [(c, cfg.nygc_vcf(c)) for c in chroms if cfg.nygc_vcf(c).exists()]
    if not avail:
        log("  No NYGC VCFs found; skipping per-sample gVCF extraction.")
        return

    ensure_dirs(cfg.tmp_dir)
    log(f"[nygc gVCF] extracting {len(kg_map)} samples from {len(avail)} chromosome(s)…")

    for _, row in kg_map.iterrows():
        eid = int(row["eid"])
        sample_id = row["sample_id"]

        dest_dir = cfg.wgs_dir() / eid_prefix_dir(eid)
        ensure_dirs(dest_dir)
        dest_vcf = dest_dir / f"{eid}_23151_0_0.g.vcf.gz"

        if dest_vcf.exists():
            continue

        # Per-chromosome extraction
        tmp_chroms: list[Path] = []
        for chrom, src in avail:
            tmp = cfg.tmp_dir / f"nygc_{eid}_chr{chrom}.vcf.gz"
            run([cfg.tools.bcftools, "view", "--samples", sample_id,
                 "--output-type", "z", "--output", str(tmp), str(src)])
            run([cfg.tools.bcftools, "index", "--tbi", str(tmp)])
            tmp_chroms.append(tmp)

        # Concatenate chromosomes
        if len(tmp_chroms) == 1:
            tmp_all = tmp_chroms[0]
        else:
            tmp_all = cfg.tmp_dir / f"nygc_{eid}_concat.vcf.gz"
            run([cfg.tools.bcftools, "concat", "--output-type", "z",
                 "--output", str(tmp_all)] + [str(t) for t in tmp_chroms])
            run([cfg.tools.bcftools, "index", "--tbi", str(tmp_all)])

        # Reheader: original sample ID → EID
        rename_tmp = cfg.tmp_dir / f"rename_{eid}.txt"
        rename_tmp.write_text(f"{sample_id}\t{eid}\n")
        run([cfg.tools.bcftools, "reheader", "--samples", str(rename_tmp),
             "--output", str(dest_vcf), str(tmp_all)])
        run([cfg.tools.bcftools, "index", "--tbi", str(dest_vcf)])

        for t in tmp_chroms:
            t.unlink(missing_ok=True)
            Path(str(t) + ".tbi").unlink(missing_ok=True)
        if tmp_all not in tmp_chroms:
            tmp_all.unlink(missing_ok=True)
            Path(str(tmp_all) + ".tbi").unlink(missing_ok=True)
        rename_tmp.unlink(missing_ok=True)

        log(f"  → {dest_vcf.name}")

    log("NYGC gVCF extraction complete.")


def _symlink_sgdp_gvcfs(cfg: Config) -> None:
    vcf_dir = cfg.sgdp_vcf_dir()
    if not vcf_dir.exists() or not any(vcf_dir.glob("*.vcf.gz")):
        log("  No SGDP VCF files found; skipping. Run `tbooo download sgdp` first.")
        return

    eid_map = _load_eid_map(cfg, "eid_map_sgdp.tsv")
    if eid_map.empty:
        log("  WARNING: SGDP EID map not found; skipping. Run `tbooo map eids` first.")
        return

    log(f"[sgdp gVCF] symlinking {len(eid_map)} per-sample VCFs…")
    linked = 0
    for _, row in eid_map.iterrows():
        eid = int(row["eid"])
        ena_acc = row["ena_accession"]

        matches = list(vcf_dir.glob(f"*{ena_acc}*.vcf.gz"))
        if not matches:
            continue
        src_vcf = matches[0].resolve()
        src_tbi = Path(str(src_vcf) + ".tbi")

        dest_dir = cfg.wgs_dir() / eid_prefix_dir(eid)
        ensure_dirs(dest_dir)
        dest_vcf = dest_dir / f"{eid}_23151_0_0.g.vcf.gz"
        dest_tbi = dest_dir / f"{eid}_23151_0_0.g.vcf.gz.tbi"

        if not dest_vcf.exists() and src_vcf.exists():
            dest_vcf.symlink_to(src_vcf)
            linked += 1
        if not dest_tbi.exists() and src_tbi.exists():
            dest_tbi.symlink_to(src_tbi)

    log(f"  linked {linked} SGDP gVCF(s) → {cfg.wgs_dir()}/")


# ── pVCF internals ────────────────────────────────────────────────────────────

def _build_nygc_pvcf(cfg: Config, chroms: list[str]) -> None:
    rename_file = cfg.metadata_dir() / "vcf_sample_rename_1kg.txt"
    if not rename_file.exists():
        log("  WARNING: 1KGP rename file not found; skipping NYGC pVCF. Run `tbooo map eids` first.")
        return

    pvcf_dir = cfg.kg_raw_dir() / "pvcf"
    ensure_dirs(pvcf_dir)

    for chrom in chroms:
        out = cfg.nygc_pvcf(chrom)
        if out.exists():
            log(f"  skip (exists): {out.name}")
            continue

        src = cfg.nygc_vcf(chrom)
        if not src.exists():
            log(f"  SKIP: NYGC VCF not found: {src.name}")
            continue

        threads = max(1, cfg.wgs_nygc_threads)
        log(f"  [nygc pVCF] chr{chrom} (threads={threads})")
        run([cfg.tools.bcftools, "reheader",
             "--samples", str(rename_file),
             "--threads", str(threads),
             "--output", str(out), str(src)])
        run([cfg.tools.bcftools, "index", "--tbi",
             "--threads", str(threads), str(out)])
        log(f"    → {out.name}")


def _build_sgdp_pvcf(cfg: Config, chroms: list[str]) -> None:
    vcf_dir = cfg.sgdp_vcf_dir()
    sample_vcfs = sorted(vcf_dir.glob("*.vcf.gz")) if vcf_dir.exists() else []
    if not sample_vcfs:
        log("  No SGDP VCF files found; skipping SGDP pVCF. Run `tbooo download sgdp` first.")
        return

    rename_file = cfg.metadata_dir() / "vcf_sample_rename_sgdp.txt"
    if not rename_file.exists():
        log("  WARNING: SGDP rename file not found; skipping. Run `tbooo map eids` first.")
        return

    pvcf_dir = cfg.sgdp_raw_dir() / "pvcf"
    ensure_dirs(pvcf_dir, cfg.tmp_dir)
    log(f"[sgdp pVCF] {len(sample_vcfs)} sample VCF(s) found")

    # Detect chromosome naming once (chrN vs N) from the first file's tabix index.
    chrom_prefix = _detect_chrom_prefix(cfg, sample_vcfs[0])

    # Write a stable file-of-filenames so the merge command stays short.
    list_file = cfg.tmp_dir / "sgdp_vcf_list.txt"
    list_file.write_text("\n".join(str(v) for v in sample_vcfs) + "\n")

    threads = max(1, cfg.wgs_sgdp_merge_threads)

    for chrom in chroms:
        out = cfg.sgdp_pvcf(chrom)
        if out.exists():
            log(f"  skip (exists): {out.name}")
            continue

        region = f"{chrom_prefix}{chrom}"
        tmp_merged = cfg.tmp_dir / f"sgdp_merged_chr{chrom}.vcf.gz"

        log(f"  [sgdp pVCF] chr{chrom} — merging {len(sample_vcfs)} samples in one pass "
            f"(region={region}, threads={threads})")
        # One `bcftools merge` reads each sample's pre-built tabix index for the
        # requested region only — no per-sample view/index step required.
        run([cfg.tools.bcftools, "merge",
             "--file-list", str(list_file),
             "--regions", region,
             "--threads", str(threads),
             "--output-type", "z",
             "--output", str(tmp_merged)])

        # Reheader (rename samples → EIDs) writes the final output.
        run([cfg.tools.bcftools, "reheader",
             "--samples", str(rename_file),
             "--threads", str(threads),
             "--output", str(out),
             str(tmp_merged)])
        run([cfg.tools.bcftools, "index", "--tbi",
             "--threads", str(threads), str(out)])

        tmp_merged.unlink(missing_ok=True)
        log(f"    → {out.name}")


def _detect_chrom_prefix(cfg: Config, vcf: Path) -> str:
    """Return 'chr' if the VCF's contigs use a chr-prefix, else ''."""
    proc = run([cfg.tools.bcftools, "index", "--stats", str(vcf)], capture=True)
    for line in proc.stdout.decode().splitlines():
        parts = line.split()
        if parts and parts[0] not in ("", "*"):
            return "chr" if parts[0].startswith("chr") else ""
    return ""


def _merge_pvcfs(cfg: Config, chroms: list[str]) -> None:
    """Merge NYGC and SGDP intermediate pVCFs into the final Field 23370 output."""
    ensure_dirs(cfg.wgs_dir())

    for chrom in chroms:
        out = cfg.wgs_pvcf(chrom)
        if out.exists():
            log(f"  skip (exists): {out.name}")
            continue

        nygc = cfg.nygc_pvcf(chrom)
        sgdp = cfg.sgdp_pvcf(chrom)
        sources = [p for p in (nygc, sgdp) if p.exists()]

        if not sources:
            log(f"  SKIP chr{chrom}: no intermediate pVCFs available")
            continue

        log(f"  [merge pVCF] chr{chrom} ({'+'.join(p.stem for p in sources)})")

        if len(sources) == 1:
            # Only one cohort available — symlink rather than copy
            out.symlink_to(sources[0].resolve())
            tbi_src = Path(str(sources[0]) + ".tbi")
            tbi_dst = Path(str(out) + ".tbi")
            if tbi_src.exists() and not tbi_dst.exists():
                tbi_dst.symlink_to(tbi_src.resolve())
        else:
            threads = max(1, cfg.wgs_final_merge_threads)
            run([cfg.tools.bcftools, "merge",
                 "--threads", str(threads),
                 "--output-type", "z",
                 "--output", str(out)] + [str(s) for s in sources])
            run([cfg.tools.bcftools, "index", "--tbi",
                 "--threads", str(threads), str(out)])

        log(f"    → {out.name}")


# ── CRAM internals ────────────────────────────────────────────────────────────

def _rename_kg_crams(cfg: Config, kg_map: pd.DataFrame) -> None:
    if kg_map.empty:
        return
    raw = cfg.kg_raw_dir()
    log(f"  linking 1KGP NYGC CRAMs ({len(kg_map)} samples)…")
    for _, row in kg_map.iterrows():
        eid = int(row["eid"])
        sample_id = row["sample_id"]
        crms = list(raw.glob(f"*{sample_id}*.cram"))
        if not crms:
            continue
        src_cram = crms[0]
        src_crai = Path(str(src_cram) + ".crai")
        _symlink_cram(cfg, eid, src_cram, src_crai)


def _symlink_cram(cfg: Config, eid: int, src_cram: Path, src_crai: Path) -> None:
    dest_dir = cfg.wgs_dir() / eid_prefix_dir(eid)
    ensure_dirs(dest_dir)
    dest_cram = dest_dir / f"{eid}_23149_0_0.cram"
    dest_crai = dest_dir / f"{eid}_23149_0_0.cram.crai"
    if not dest_cram.exists() and src_cram.exists():
        dest_cram.symlink_to(src_cram.resolve())
    if not dest_crai.exists() and src_crai.exists():
        dest_crai.symlink_to(src_crai.resolve())


def _load_eid_map(cfg: Config, filename: str) -> pd.DataFrame:
    path = cfg.metadata_dir() / filename
    if not path.exists():
        log(f"  WARNING: {filename} not found; skipping")
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")
