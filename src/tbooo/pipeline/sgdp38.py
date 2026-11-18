"""SGDP on GRCh38, phased against the 1000 Genomes 30x panel.

The SGDP release is one file per sample on GRCh37 (hs37d5), with a record for every base of
the genome and unphased genotypes. For each autosome this stage writes one phased GRCh38 file
of the variant records of every sample:

  1. Variant records. A genotype at the project's filter level N or 0 is not a call and is set
     to missing. A record no sample then carries is dropped, and so are INFO and every FORMAT
     field but GT.
  2. Multiallelic records are split and indels left-aligned on GRCh37. The chromosome gets the
     chr prefix of GRCh38.
  3. Picard LiftoverVcf with the UCSC hg19ToHg38 chain. A record that does not lift, or lifts
     to another chromosome, is dropped. Exact duplicates are dropped.
  4. Beagle phases the samples against the panel. Beagle writes the reference's markers alone,
     so the SGDP records the panel lacks are first added to the reference as markers where
     every panel sample has the reference allele. They are then kept, and phased from the
     other SGDP samples that carry them. A variant one sample alone carries has no information
     for its phase, here or in any method without reads. Missing genotypes of a kept record
     are filled by Beagle. Markers of the panel that SGDP has no record of are not imputed.
  5. The samples are renamed to EIDs.

The X chromosome is not built: male genotypes there are haploid and need their own handling.

Outputs, per chromosome:
    data/raw/sgdp/grch38/sgdp_c{chrom}.unphased.vcf.gz   step 3, sample names as in the source
    data/raw/sgdp/grch38/sgdp_c{chrom}.phased.vcf.gz     step 5
    data/raw/sgdp/grch38/sgdp_c{chrom}.counts.tsv        records after each step
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from tbooo.config import Config
from tbooo.integrity import ensure_vcf, remove, vcf_ok
from tbooo.utils import ensure_dirs, log, run

_GRCH37_FASTA = "GRCh37/human_g1k_v37.fasta"
_GRCH38_FASTA = "GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa"
CHAIN = "liftover/hg19ToHg38.over.chain.gz"
BEAGLE_MAP = "genetic_maps/beagle_GRCh38/chr_in_chrom_field/plink.chrchr{chrom}.GRCh38.map"

# FL is the SGDP filter level of a genotype, 0 to 9 or N. The project's README reads N and 0
# as no call.
_NOT_A_CALL = 'FMT/FL="N" | FMT/FL="0"'


def build_sgdp_grch38(cfg: Config, chroms: list[str]) -> None:
    """Build the phased GRCh38 file of each autosome in `chroms`."""
    autosomes = {str(a) for a in cfg.autosomes}
    for chrom in chroms:
        if chrom not in autosomes:
            log(f"  SKIP chr{chrom}: the SGDP GRCh38 stage builds autosomes only")
            continue
        _build_chrom(cfg, chrom)
    log("SGDP GRCh38 build complete.")


def _build_chrom(cfg: Config, chrom: str) -> None:
    out = cfg.sgdp38_phased(chrom)
    unphased = cfg.sgdp38_unphased(chrom)
    work = cfg.tmp_dir / f"sgdp38_chr{chrom}"
    counts: dict[str, int] = {}

    def _lifted() -> None:
        ensure_dirs(work, unphased.parent)
        variants = work / "1.variants.bcf"
        normed = work / "2.norm.vcf.gz"
        _variant_records(cfg, chrom, variants)
        counts["variant records on GRCh37"] = _count(cfg, variants)
        _normalise(cfg, chrom, variants, normed)
        counts["after the split of multiallelic records"] = _count(cfg, normed)
        _lift(cfg, chrom, normed, unphased, work, counts)

    def _phased() -> None:
        ensure_dirs(work)
        raw = work / "5.beagle.vcf.gz"
        _phase(cfg, chrom, unphased, raw, work, counts)
        _rename_to_eids(cfg, raw, out)
        counts["phased records"] = _count(cfg, out)

    ensure_vcf(unphased.name, unphased, _lifted)
    ensure_vcf(out.name, out, _phased)
    if counts:
        table = cfg.sgdp38_counts(chrom)
        table.write_text("".join(f"{k}\t{v}\n" for k, v in counts.items()))
        log(f"  wrote {table.name}")


# ── Steps ─────────────────────────────────────────────────────────────────────

def _variant_records(cfg: Config, chrom: str, out: Path) -> None:
    """Step 1, from the merged all-site file where one exists, else from the sample files."""
    bcftools, threads = cfg.tools.bcftools, str(max(1, cfg.wgs_sgdp_merge_threads))
    merged = cfg.sgdp_pvcf(chrom)
    if vcf_ok(merged):
        log(f"  [sgdp38] chr{chrom}: variant records of {merged.name}")
        source = [bcftools, "view", "--threads", threads, "--min-ac", "1", "-Ou", str(merged)]
    else:
        files = sorted(cfg.sgdp_vcf_dir().glob("*.vcf.gz")) if cfg.sgdp_vcf_dir().exists() else []
        if not files:
            raise FileNotFoundError(
                f"No SGDP VCFs in {cfg.sgdp_vcf_dir()}. Run `tbooo download sgdp` first.")
        listing = out.parent / "sample_files.txt"
        listing.write_text("".join(f"{f}\n" for f in files))
        log(f"  [sgdp38] chr{chrom}: merging {len(files)} sample files, every base of the "
            "chromosome. This is the slow route; nothing of it is written to disk.")
        source = [bcftools, "merge", "--file-list", str(listing), "--regions", chrom,
                  "--threads", threads, "-Ou"]
    _pipe(
        source,
        [bcftools, "+setGT", "-Ou", "--", "-t", "q", "-n", ".", "-i", _NOT_A_CALL],
        [bcftools, "annotate", "--remove", "INFO,^FORMAT/GT", "-Ou"],
        [bcftools, "view", "--threads", threads, "--min-ac", "1", "-Ob", "-o", str(out)],
    )
    run([bcftools, "index", "--force", "--threads", threads, str(out)])


def _normalise(cfg: Config, chrom: str, src: Path, out: Path) -> None:
    """Step 2."""
    bcftools, threads = cfg.tools.bcftools, str(max(1, cfg.wgs_sgdp_merge_threads))
    fasta = cfg.reference_dir / _GRCH37_FASTA
    if not fasta.exists():
        raise FileNotFoundError(f"{fasta} not found. Run `tbooo download reference` first.")
    names = out.parent / "chr.map"
    names.write_text(f"{chrom} chr{chrom}\n")
    _pipe(
        [bcftools, "norm", "--threads", threads, "--multiallelics", "-any",
         "--fasta-ref", str(fasta), "--check-ref", "w", "-Ou", str(src)],
        [bcftools, "annotate", "--rename-chrs", str(names), "-Oz", "-o", str(out)],
    )
    run([bcftools, "index", "--force", "--tbi", str(out)])


def _lift(cfg: Config, chrom: str, src: Path, out: Path, work: Path, counts: dict[str, int]) -> None:
    """Step 3."""
    bcftools = cfg.tools.bcftools
    fasta = cfg.reference_dir / _GRCH38_FASTA
    chain = cfg.reference_dir / CHAIN
    for needed in (fasta, chain):
        if not needed.exists():
            raise FileNotFoundError(f"{needed} not found. Run `tbooo download reference` first.")
    heap = f"-Xmx{cfg.sgdp38_java_gb}g"
    dictionary = fasta.with_suffix(".dict")
    if not dictionary.exists():
        run([cfg.tools.picard, heap, "CreateSequenceDictionary", f"R={fasta}", f"O={dictionary}"])
    lifted, rejected = work / "3.lifted.vcf.gz", work / "3.rejected.vcf.gz"
    run([cfg.tools.picard, heap, "LiftoverVcf", f"I={src}", f"O={lifted}", f"CHAIN={chain}",
         f"REJECT={rejected}", f"R={fasta}", "RECOVER_SWAPPED_REF_ALT=true",
         "WARN_ON_MISSING_CONTIG=true", "MAX_RECORDS_IN_RAM=200000"])
    counts["not lifted"] = _count_lines(cfg, rejected)
    # Picard's own index holds no record count.
    run([bcftools, "index", "--force", "--tbi", str(lifted)])
    counts["lifted"] = _count(cfg, lifted)
    _pipe(
        [bcftools, "view", "--regions", f"chr{chrom}", "-Ou", str(lifted)],
        [bcftools, "sort", "--temp-dir", str(work / "sort"), "-Ou"],
        [bcftools, "norm", "--rm-dup", "exact", "-Oz", "-o", str(out)],
    )
    run([bcftools, "index", "--force", "--tbi", str(out)])
    counts["on this chromosome of GRCh38, duplicates removed"] = _count(cfg, out)


def _phase(cfg: Config, chrom: str, src: Path, out: Path, work: Path, counts: dict[str, int]) -> None:
    """Step 4."""
    bcftools, threads = cfg.tools.bcftools, str(max(1, cfg.wgs_sgdp_merge_threads))
    panel = cfg.nygc_vcf(chrom)
    genetic_map = cfg.reference_dir / BEAGLE_MAP.format(chrom=chrom)
    for needed in (panel, genetic_map):
        if not needed.exists():
            raise FileNotFoundError(f"{needed} not found. Run `tbooo download 1kg --nygc` and "
                                    "`tbooo download reference` first.")

    # One record per position on both sides: Beagle matches markers by position.
    target = work / "4.target.vcf.gz"
    _pipe([bcftools, "norm", "--rm-dup", "all", "--threads", threads, "-Oz", "-o", str(target), str(src)])
    run([bcftools, "index", "--force", "--tbi", str(target)])
    counts["one record per position"] = _count(cfg, target)

    # Records of the target the panel lacks, by exact alleles and by position. The second set
    # becomes new markers. The difference sits at a position the panel has with other alleles.
    lacking = _pipe([bcftools, "isec", "--complement", "--write", "1", "--collapse", "none",
                     str(target), str(panel)],
                    [bcftools, "query", "-f", "%POS\n"], capture=True).count(b"\n")
    added = [line.split("\t") for line in _pipe(
        [bcftools, "isec", "--complement", "--write", "1", "--collapse", "all", str(target), str(panel)],
        [bcftools, "query", "-f", "%CHROM\t%POS\t%ID\t%REF\t%ALT\n"], capture=True).decode().splitlines()]
    counts["records the panel lacks"] = lacking
    counts["of them at a position the panel has with other alleles, not phased"] = lacking - len(added)

    reference = work / "4.reference.vcf.gz"
    extra = work / "4.added.vcf.gz"
    header = run([bcftools, "view", "--header-only", str(panel)], capture=True).stdout.decode()
    n_panel = len(header.rstrip("\n").rsplit("\n", 1)[-1].split("\t")) - 9
    write_added_markers(header, added, n_panel, extra, cfg.tools.bgzip)
    run([bcftools, "index", "--force", "--tbi", str(extra)])
    run([bcftools, "concat", "--allow-overlaps", "--threads", threads, "-Oz", "-o", str(reference),
         str(panel), str(extra)])
    run([bcftools, "index", "--force", "--tbi", str(reference)])

    prefix = str(out)[: -len(".vcf.gz")]
    run([cfg.tools.beagle, f"-Xmx{cfg.sgdp38_java_gb}g", f"gt={target}", f"ref={reference}",
         f"map={genetic_map}", f"out={prefix}", f"nthreads={threads}", "impute=false"])
    run([bcftools, "index", "--force", "--tbi", str(out)])
    remove(reference, Path(str(reference) + ".tbi"), extra, Path(str(extra) + ".tbi"))


def write_added_markers(header: str, sites: list[list[str]], n_samples: int, out: Path, bgzip: str = "bgzip") -> None:
    """Write a VCF of `sites` where each of `n_samples` samples carries the reference allele.

    @param header the panel's VCF header, whose sample columns the body follows
    @param sites CHROM, POS, ID, REF and ALT of each record, sorted by position
    """
    genotypes = "\t0|0" * n_samples
    with subprocess.Popen([bgzip, "-c"], stdin=subprocess.PIPE, stdout=out.open("wb")) as proc:
        assert proc.stdin is not None
        proc.stdin.write(header.encode())
        for chrom, pos, _id, ref, alt in sites:
            proc.stdin.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT{genotypes}\n".encode())
        proc.stdin.close()
        if proc.wait() != 0:
            raise RuntimeError(f"bgzip failed writing {out}")


def _rename_to_eids(cfg: Config, src: Path, out: Path) -> None:
    """Step 5. Refuses a file whose samples are not all in the EID map."""
    rename = cfg.metadata_dir() / "vcf_sample_rename_sgdp.txt"
    if not rename.exists():
        raise FileNotFoundError(f"{rename} not found. Run `tbooo map eids` first.")
    pairs = dict(line.split()[:2] for line in rename.read_text().splitlines() if line.strip())
    eids = set(pairs.values())
    names = run([cfg.tools.bcftools, "query", "--list-samples", str(src)], capture=True).stdout.decode().split()
    unknown = [n for n in names if n not in pairs and n not in eids]
    if unknown:
        raise RuntimeError(f"{len(unknown)} sample(s) of {src.name} are not in {rename.name}, "
                           f"the first is {unknown[0]}")
    listing = out.parent / f"{out.name}.samples.txt"
    listing.write_text("".join(f"{pairs.get(n, n)}\n" for n in names))
    run([cfg.tools.bcftools, "reheader", "--samples", str(listing), "--output", str(out), str(src)])
    run([cfg.tools.bcftools, "index", "--force", "--tbi", str(out)])
    listing.unlink()


# ── helpers ───────────────────────────────────────────────────────────────────

def _pipe(*commands: list[str], capture: bool = False) -> bytes:
    """Run the commands as one pipeline and raise where any of them fails.

    @param capture return the last command's output where True, else leave it on the console
    """
    procs: list[subprocess.Popen] = []
    for i, command in enumerate(commands):
        last = i == len(commands) - 1
        stdin = procs[-1].stdout if procs else None
        stdout = subprocess.PIPE if (capture or not last) else None
        procs.append(subprocess.Popen(command, stdin=stdin, stdout=stdout))
        if stdin is not None:
            stdin.close()
    out = procs[-1].stdout.read() if capture and procs[-1].stdout is not None else b""
    codes = [p.wait() for p in procs]
    for command, code in zip(commands, codes):
        if code != 0:
            raise subprocess.CalledProcessError(code, command)
    return out


def _count(cfg: Config, vcf: Path) -> int:
    """Records of an indexed VCF or BCF."""
    out = run([cfg.tools.bcftools, "index", "--nrecords", str(vcf)], capture=True).stdout.decode()
    return int(out.strip())


def _count_lines(cfg: Config, vcf: Path) -> int:
    """Records of a VCF without an index."""
    out = run([cfg.tools.bcftools, "view", "--no-header", str(vcf)], capture=True).stdout
    return out.count(b"\n")
