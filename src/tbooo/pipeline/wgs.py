"""Build UKB-mirrored WGS data (Fields 24051, 24310).

Per-sample gVCFs (Field 24051, "Whole genome variant call files (GVCFs) (DRAGEN)") — both cohorts:
  - 1KGP: extract per-sample VCF from NYGC 30x cohort VCFs → <EID>_24051_0_0.g.vcf.gz
  - SGDP: symlink downloaded per-sample VCFs        → <EID>_24051_0_0.g.vcf.gz
  The two differ in build and content. A 1KGP file is on GRCh38 and holds the records its
  sample carries. An SGDP file is the project's own release on GRCh37, with every base.

Cohort pVCF (Field 24310, "DRAGEN population level WGS variants, pVCF format") — both cohorts merged:
  - 1KGP: reheader NYGC per-chromosome VCF with EIDs → raw/1kg/pvcf/nygc_c{chrom}.pvcf.gz
  - SGDP: the phased GRCh38 file of `tbooo map sgdp38` → raw/sgdp/grch38/sgdp_c{chrom}.phased.vcf.gz
  - Combined: bcftools merge nygc + sgdp             → Bulk/Whole genome sequences/ukb24310_c{chrom}_b0_v1.pvcf.gz
  Both inputs hold variant records alone. In the merged file a cohort's samples are missing at
  a record only the other cohort has: the cohort either carries the reference allele there or
  has no call, and the files do not say which.

Individual CRAMs (Field 24048) are not mirrored: the 30x CRAMs are not downloaded
(tens of GB per sample), so there is no source to link.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from tbooo.config import Config
from tbooo.integrity import ensure, ensure_vcf, remove, symlink_ok, vcf_ok
from tbooo.utils import ensure_dirs, eid_prefix_dir, has_bgzf_eof, log, run


# ── Public API ────────────────────────────────────────────────────────────────

def build_gvcfs(cfg: Config, chroms: list[str]) -> None:
    """Build per-sample gVCFs (Field 24051) from both NYGC and SGDP.

    1KGP: extracts each sample from NYGC 30x cohort VCFs and concatenates chromosomes.
    SGDP: symlinks downloaded per-sample VCFs directly.
    """
    _build_nygc_gvcfs(cfg, chroms)
    _symlink_sgdp_gvcfs(cfg)
    log("gVCF build complete.")


def build_pvcf(cfg: Config, chroms: list[str]) -> None:
    """Build cohort pVCF (Field 24310) by merging NYGC and SGDP per-chromosome VCFs.

    The NYGC side is built here. The SGDP side is the output of `tbooo map sgdp38`, on the
    same build as NYGC. Where it is absent the NYGC file is used alone.
    """
    _build_nygc_pvcf(cfg, chroms)
    _merge_pvcfs(cfg, chroms)
    log("pVCF build complete.")


# ── gVCF internals ────────────────────────────────────────────────────────────

# +split applies --include to each output file on its own, so every sample's file gets the
# records where that sample has a non-reference allele. --keep-tags drops INFO, whose allele
# counts describe the cohort and are wrong for a single sample, and keeps the genotype.
_SPLIT_CARRIED_ONLY = ["--include", 'GT="alt"', "--keep-tags", "FMT/GT"]


def _build_nygc_gvcfs(cfg: Config, chroms: list[str]) -> None:
    """Extract per-sample 1KGP gVCFs from the NYGC 30x cohort VCFs.

    Uses `bcftools +split` to fan a multi-sample chromosome VCF into per-sample
    files in ONE pass per chromosome (≈ n_chroms splits total), then concatenates
    each sample's per-chrom pieces. The previous per-sample `view` approach scanned
    every multi-GB chromosome VCF once per sample (n_samples × n_chroms full scans),
    which is intractable for 3,202 samples.

    Each file holds the records its sample carries, and the genotype alone. Without
    that filter +split writes every site of the cohort into every sample's file with
    the cohort's INFO on each line: a sample carries about 6% of them, so chromosome 1
    alone came to 2.5 TB and the whole genome to tens of terabytes, twice over once the
    pieces were concatenated. A phased panel has no reference blocks to give, so these
    are variant-only files under the gVCF field's name.
    """
    kg_map = _load_eid_map(cfg, "eid_map_1kg.tsv")
    if kg_map.empty:
        log("  WARNING: 1KGP EID map not found; skipping NYGC gVCF extraction.")
        return

    avail = [(c, cfg.nygc_vcf(c)) for c in chroms if cfg.nygc_vcf(c).exists()]
    if not avail:
        log("  No NYGC VCFs found; skipping per-sample gVCF extraction.")
        return

    samples = [(int(r["eid"]), str(r["sample_id"])) for _, r in kg_map.iterrows()]
    dest_of = {
        eid: cfg.wgs_dir() / eid_prefix_dir(eid) / f"{eid}_24051_0_0.g.vcf.gz"
        for eid, _ in samples
    }
    if all(vcf_ok(d) for d in dest_of.values()):
        log(f"[nygc gVCF] all {len(samples)} per-sample gVCFs present & valid; skipping")
        return

    ensure_dirs(cfg.tmp_dir)
    threads = max(1, cfg.wgs_nygc_threads)

    # `bcftools +split` holds one output file open per sample at once, so a
    # 3,202-sample cohort exceeds the open-file limit. Raise the soft FD limit
    # toward the hard cap, then split in sample batches sized to stay under it.
    sample_ids = [sid for _, sid in samples]
    batch_size = _fd_batch_size(len(sample_ids))

    # Step 1: `+split` each chromosome → tmp/nygc_split_chr<c>/<sample_id>.vcf.gz.
    # A `.split-complete` marker makes a killed run resumable without re-splitting
    # a chromosome that already finished. (+split is a plugin with no --threads.)
    split_dirs: list[Path] = []
    for chrom, src in avail:
        split_dir = cfg.tmp_dir / f"nygc_variants_chr{chrom}"
        marker = split_dir / ".split-complete"
        stale = cfg.tmp_dir / f"nygc_split_chr{chrom}"
        if stale.exists():
            log(f"  NOTE: {stale} holds unfiltered per-sample files from an earlier "
                "version (every cohort site in every file). It is not used and can be deleted.")
        if not marker.exists():
            ensure_dirs(split_dir)
            n_batches = (len(sample_ids) + batch_size - 1) // batch_size
            log(f"  [nygc gVCF] splitting chr{chrom} into per-sample VCFs "
                f"({len(sample_ids)} samples, {n_batches} batch(es))…")
            for i in range(0, len(sample_ids), batch_size):
                cmd = [cfg.tools.bcftools, "+split", str(src),
                       *_SPLIT_CARRIED_ONLY,
                       "--output-type", "z", "--output", str(split_dir)]
                if n_batches > 1:
                    sfile = cfg.tmp_dir / f"split_samples_chr{chrom}_{i}.txt"
                    sfile.write_text("\n".join(sample_ids[i:i + batch_size]) + "\n")
                    run(cmd + ["--samples-file", str(sfile)])
                    sfile.unlink(missing_ok=True)
                else:
                    run(cmd)
            marker.write_text("")
        split_dirs.append(split_dir)

    # Step 2: assemble each sample's gVCF by concatenating its per-chrom pieces,
    # then reheader the sample name → EID.
    log(f"[nygc gVCF] assembling {len(samples)} per-sample gVCFs from {len(avail)} chrom(s)…")
    for eid, sample_id in samples:
        dest_vcf = dest_of[eid]
        if vcf_ok(dest_vcf):
            continue
        ensure_dirs(dest_vcf.parent)
        remove(dest_vcf, Path(str(dest_vcf) + ".tbi"))

        pieces = [sd / f"{sample_id}.vcf.gz" for sd in split_dirs]
        pieces = [p for p in pieces if p.exists()]
        if not pieces:
            continue

        tmp_all = cfg.tmp_dir / f"nygc_{eid}_concat.vcf.gz"
        run([cfg.tools.bcftools, "concat", "--output-type", "z",
             "--threads", str(threads), "--output", str(tmp_all)]
            + [str(p) for p in pieces])

        rename_tmp = cfg.tmp_dir / f"rename_{eid}.txt"
        rename_tmp.write_text(f"{sample_id}\t{eid}\n")
        run([cfg.tools.bcftools, "reheader", "--samples", str(rename_tmp),
             "--output", str(dest_vcf), str(tmp_all)])
        run([cfg.tools.bcftools, "index", "--tbi",
             "--threads", str(threads), str(dest_vcf)])

        remove(tmp_all, rename_tmp)

    # Step 3: reclaim disk only once every sample is assembled — otherwise keep the
    # split dirs (with their markers) so a resumed run reuses them instead of
    # re-splitting every chromosome.
    if all(vcf_ok(d) for d in dest_of.values()):
        for sd in split_dirs:
            _rmtree(sd)
        log("NYGC gVCF extraction complete.")
    else:
        log("NYGC gVCF extraction incomplete (some samples missing source data); "
            "split dirs kept for resume.")


def _fd_batch_size(n_samples: int) -> int:
    """Largest sample batch that keeps `bcftools +split` under the open-file limit.

    +split opens one output file per sample simultaneously. Raise the soft
    RLIMIT_NOFILE toward the hard cap, then leave a margin for stdio/htslib FDs.
    Returns a batch size in [1, n_samples]; n_samples means a single pass suffices.
    """
    try:
        import resource
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        need = n_samples + 256
        target = need if hard == resource.RLIM_INFINITY else min(hard, need)
        if soft != resource.RLIM_INFINITY and soft < target:
            try:
                resource.setrlimit(resource.RLIMIT_NOFILE, (target, hard))
                soft = target
            except (ValueError, OSError):
                pass
        if soft == resource.RLIM_INFINITY:
            return n_samples
        return max(1, min(n_samples, soft - 128))
    except Exception:
        return max(1, min(n_samples, 256))


def _rmtree(d: Path) -> None:
    """Remove a split directory and its contents (best effort)."""
    if not d.exists():
        return
    for f in d.iterdir():
        try:
            f.unlink()
        except OSError:
            pass
    try:
        d.rmdir()
    except OSError:
        pass


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
        name = row.get("vcf_filename")
        if not isinstance(name, str) or not (vcf_dir / name).exists():
            continue
        src_vcf = (vcf_dir / name).resolve()
        src_tbi = Path(str(src_vcf) + ".tbi")

        dest_dir = cfg.wgs_dir() / eid_prefix_dir(eid)
        ensure_dirs(dest_dir)
        dest_vcf = dest_dir / f"{eid}_24051_0_0.g.vcf.gz"
        dest_tbi = dest_dir / f"{eid}_24051_0_0.g.vcf.gz.tbi"

        if _relink(dest_vcf, src_vcf):
            linked += 1
        _relink(dest_tbi, src_tbi)

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
        src = cfg.nygc_vcf(chrom)
        if not src.exists():
            log(f"  SKIP: NYGC VCF not found: {src.name}")
            continue

        def _build(out=out, src=src) -> None:
            threads = max(1, cfg.wgs_nygc_threads)
            log(f"  [nygc pVCF] chr{chrom} (threads={threads})")
            run([cfg.tools.bcftools, "reheader",
                 "--samples", str(rename_file),
                 "--threads", str(threads),
                 "--output", str(out), str(src)])
            run([cfg.tools.bcftools, "index", "--tbi",
                 "--threads", str(threads), str(out)])
            log(f"    → {out.name}")

        ensure_vcf(out.name, out, _build)


def _merge_pvcfs(cfg: Config, chroms: list[str]) -> None:
    """Merge NYGC and SGDP intermediate pVCFs into the final Field 24310 output."""
    ensure_dirs(cfg.wgs_dir())

    for chrom in chroms:
        out = cfg.wgs_pvcf(chrom)
        out_tbi = Path(str(out) + ".tbi")
        nygc = cfg.nygc_pvcf(chrom)
        sgdp = cfg.sgdp38_phased(chrom)
        sources = [p for p in (nygc, sgdp) if p.exists()]
        if not sgdp.exists():
            log(f"  chr{chrom}: no phased GRCh38 SGDP file; run `tbooo map sgdp38` to add SGDP")

        if not sources:
            log(f"  SKIP chr{chrom}: no intermediate pVCFs available")
            continue

        if len(sources) == 1:
            src = sources[0]

            def _check(out=out, out_tbi=out_tbi) -> bool:
                # Symlinked single-cohort output: link resolves to an intact BGZF
                # file and its .tbi symlink is not dangling.
                return symlink_ok(out, has_bgzf_eof) and symlink_ok(out_tbi)

            def _build(out=out, out_tbi=out_tbi, src=src) -> None:
                log(f"  [merge pVCF] chr{chrom} ({src.stem}, single cohort → symlink)")
                out.symlink_to(src.resolve())
                tbi_src = Path(str(src) + ".tbi")
                if tbi_src.exists():
                    out_tbi.symlink_to(tbi_src.resolve())
                log(f"    → {out.name}")
        else:
            def _check(out=out) -> bool:
                return vcf_ok(out)

            def _build(out=out, sources=tuple(sources)) -> None:
                threads = max(1, cfg.wgs_final_merge_threads)
                log(f"  [merge pVCF] chr{chrom} ({'+'.join(p.stem for p in sources)})")
                run([cfg.tools.bcftools, "merge",
                     "--threads", str(threads),
                     "--output-type", "z",
                     "--output", str(out)] + [str(s) for s in sources])
                run([cfg.tools.bcftools, "index", "--tbi",
                     "--threads", str(threads), str(out)])
                log(f"    → {out.name}")

        ensure(
            out.name,
            check=_check,
            build=_build,
            purge=lambda out=out, out_tbi=out_tbi: remove(
                out, out_tbi, Path(str(out) + ".csi")),
        )


# ── symlink helper ─────────────────────────────────────────────────────────────

def _relink(dest: Path, src: Path) -> bool:
    """Point `dest` at `src` via symlink; recreate it if missing or dangling.

    Returns True if a (re)link was performed. A valid existing symlink is left alone.
    """
    if not src.exists():
        return False
    if symlink_ok(dest):
        return False
    remove(dest)
    dest.symlink_to(src.resolve())
    return True


def _load_eid_map(cfg: Config, filename: str) -> pd.DataFrame:
    path = cfg.metadata_dir() / filename
    if not path.exists():
        log(f"  WARNING: {filename} not found; skipping")
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")
