"""Download reference files needed by the pipeline."""

from __future__ import annotations

import gzip
from pathlib import Path

from tbooo.config import Config
from tbooo.utils import ensure_dirs, log, run, wget_download

# AstraZeneca NGS reference_data repo — plain BED files, no login required
_AZ_REF_BASE = "https://raw.githubusercontent.com/AstraZeneca-NGS/reference_data/master/hg38/bed"

# IDT xGen Exome Research Panel v1.0 (hg38, matches UKB WES capture panel)
_IDT_V1_URL = f"{_AZ_REF_BASE}/Exome-IDT_V1.bed"

# GENCODE GRCh38 comprehensive annotation — fallback to derive exome intervals
# from protein-coding exon coordinates when IDT BED is unavailable
_GENCODE_GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_47/gencode.v47.annotation.gtf.gz"
)

_GRCH37_FASTA_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/"
    "human_g1k_v37.fasta.gz"
)
_GRCH38_FASTA_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/"
    "GRCh38_full_analysis_set_plus_decoy_hla.fa"
)


def download_all(cfg: Config) -> None:
    ensure_dirs(cfg.reference_dir)
    _download_exome_bed(cfg)
    _download_genetic_maps(cfg)
    _download_grch37(cfg)
    _download_grch38(cfg)
    log("All reference files ready.")


# ── Exome BED ─────────────────────────────────────────────────────────────────

def _download_exome_bed(cfg: Config) -> None:
    dest = Path(cfg.exome_bed)
    if dest.exists():
        log(f"  skip (exists): {dest.name}")
        return
    ensure_dirs(dest.parent)

    # Primary: IDT V1 from AstraZeneca-NGS/reference_data (plain BED, direct download)
    try:
        log(f"  downloading IDT xGen V1 exome BED from AstraZeneca-NGS/reference_data…")
        wget_download(_IDT_V1_URL, dest, tool_wget=cfg.tools.wget)
        log(f"  exome BED ready: {dest.name} ({dest.stat().st_size // 1024} KB)")
        return
    except Exception as exc:
        log(f"  IDT V1 download failed: {exc}")
        dest.unlink(missing_ok=True)

    # Fallback: generate from GENCODE v47 protein-coding exons (50 bp padding)
    log("  falling back to GENCODE v47 protein-coding exon intervals…")
    _generate_from_gencode(cfg, dest)


def _generate_from_gencode(cfg: Config, dest: Path) -> None:
    gtf_gz = cfg.reference_dir / "gencode.v47.annotation.gtf.gz"
    if not gtf_gz.exists():
        log("  downloading GENCODE v47 GTF…")
        wget_download(_GENCODE_GTF_URL, gtf_gz, tool_wget=cfg.tools.wget)

    log("  extracting protein-coding exon intervals (50 bp padding)…")
    pad = 50
    intervals: list[tuple[str, int, int]] = []
    with gzip.open(gtf_gz, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            if 'gene_type "protein_coding"' not in fields[8]:
                continue
            chrom = fields[0]
            start = max(0, int(fields[3]) - 1 - pad)
            end = int(fields[4]) + pad
            intervals.append((chrom, start, end))

    intervals.sort()
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in intervals:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(merged[-1][2], end))
        else:
            merged.append((chrom, start, end))

    with open(dest, "w") as out:
        for chrom, start, end in merged:
            out.write(f"{chrom}\t{start}\t{end}\n")

    log(f"  generated {len(merged)} intervals from GENCODE → {dest.name}")


# ── Genetic maps ──────────────────────────────────────────────────────────────

def _download_genetic_maps(cfg: Config) -> None:
    maps_dir = cfg.reference_dir / "genetic_maps"
    tarball = cfg.reference_dir / "1000gp_genetic_maps.tar.gz"
    if maps_dir.exists() and any(maps_dir.iterdir()):
        log("  skip (exists): genetic_maps/")
        return
    ensure_dirs(maps_dir)
    url = "https://github.com/joepickrell/1000-genomes-genetic-maps/archive/refs/heads/master.tar.gz"
    log("Downloading 1KGP pedigree-based genetic maps…")
    wget_download(url, tarball, tool_wget=cfg.tools.wget)
    run(["tar", "-xzf", str(tarball), "-C", str(maps_dir), "--strip-components=1"])
    tarball.unlink()
    log(f"  genetic maps ready → {maps_dir}")


# ── Reference FASTAs ──────────────────────────────────────────────────────────

def _download_grch37(cfg: Config) -> None:
    # The 1KGP file is in legacy RAZF format — samtools cannot index it directly.
    # We decompress to a plain FASTA and index that instead.
    fasta = cfg.reference_dir / "GRCh37" / "human_g1k_v37.fasta"
    fai = Path(str(fasta) + ".fai")
    if fai.exists():
        log(f"  skip (exists): {fasta.name}")
        return
    ensure_dirs(fasta.parent)

    fasta_gz = fasta.parent / "human_g1k_v37.fasta.gz"
    if not fasta_gz.exists():
        log("Downloading GRCh37 reference FASTA (RAZF-compressed)…")
        wget_download(_GRCH37_FASTA_URL, fasta_gz, tool_wget=cfg.tools.wget)

    # RAZF format: gzip data followed by an index trailer that gunzip doesn't
    # understand. gunzip exits 2 ("trailing garbage ignored") but the output is
    # valid, so we allow exit codes 0 and 2 and verify the file afterwards.
    log("  decompressing RAZF FASTA (trailing garbage warning from gunzip is expected)…")
    result = run(["gunzip", "--force", "--keep", str(fasta_gz)], check=False)
    if result.returncode not in (0, 2):
        raise RuntimeError(
            f"gunzip failed with exit code {result.returncode} on {fasta_gz}"
        )
    if not fasta.exists() or fasta.stat().st_size == 0:
        raise RuntimeError(f"gunzip produced no output for {fasta_gz}")

    run([cfg.tools.samtools, "faidx", str(fasta)])
    log(f"  GRCh37 FASTA ready → {fasta}")


def _download_grch38(cfg: Config) -> None:
    dest = cfg.reference_dir / "GRCh38" / "GRCh38_full_analysis_set_plus_decoy_hla.fa"
    fai = Path(str(dest) + ".fai")
    if fai.exists():
        log(f"  skip (exists): {dest.name}")
        return
    ensure_dirs(dest.parent)
    if not dest.exists():
        log("Downloading GRCh38DH reference FASTA…")
        wget_download(_GRCH38_FASTA_URL, dest, tool_wget=cfg.tools.wget)

    # The 1KGP GRCh38 FASTA has inconsistent sequence line lengths in some
    # contigs, which causes "Different line length" errors in samtools faidx.
    # Reformat to a uniform 60-char line width before indexing.
    log("  normalizing FASTA line widths to 60 chars (required for samtools faidx)…")
    _normalize_fasta(dest)

    run([cfg.tools.samtools, "faidx", str(dest)])
    log(f"  GRCh38 FASTA ready → {dest}")


def _normalize_fasta(fasta: Path, width: int = 60) -> None:
    """Rewrite a FASTA in-place with uniform sequence line length.

    Handles mixed line endings and guarantees every sequence line is exactly
    `width` characters wide (except the final line of each record).
    Uses a temp file so the original is never partially overwritten.
    """
    tmp = fasta.parent / (fasta.name + ".normalizing")
    with open(fasta, "r", newline="") as inp, open(tmp, "w") as out:
        seq_parts: list[str] = []

        def flush_seq():
            seq = "".join(seq_parts)
            for i in range(0, len(seq), width):
                out.write(seq[i : i + width] + "\n")
            seq_parts.clear()

        for raw_line in inp:
            line = raw_line.rstrip("\r\n")
            if line.startswith(">"):
                flush_seq()
                out.write(line + "\n")
            else:
                seq_parts.append(line)
        flush_seq()

    tmp.replace(fasta)
