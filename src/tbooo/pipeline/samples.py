"""BGEN sample files written in the order of the genotype file they describe."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from tbooo.config import Config
from tbooo.utils import log, run


def vcf_samples(cfg: Config, vcf: Path) -> list[str]:
    """The sample names of a VCF, in column order."""
    proc = run([cfg.tools.bcftools, "query", "--list-samples", str(vcf)], capture=True)
    return [line for line in proc.stdout.decode().splitlines() if line]


def sample_file_text(samples: list[str], eid_map: pd.DataFrame) -> str:
    """A GEN2 sample file with one row per name in `samples`, in that order.

    Sex comes from the EID map where the name is an EID in it, else 0. The file must follow
    the genotype file row for row: a BGEN carries no names of its own that a reader checks, so
    a sample file in any other order assigns every genotype to the wrong participant.
    """
    sex_of = {str(e): int(s) for e, s in zip(eid_map["eid"], eid_map.get("sex", pd.Series(0, index=eid_map.index)).fillna(0))}
    lines = ["ID_1 ID_2 missing sex\n", "0 0 0 D\n"]
    lines += [f"{name} {name} 0 {sex_of.get(str(name), 0)}\n" for name in samples]
    return "".join(lines)


def write_sample_file(cfg: Config, vcf: Path, eid_map: pd.DataFrame, path: Path) -> None:
    """Write the sample file of the BGEN made from `vcf`, replacing one that differs."""
    text = sample_file_text(vcf_samples(cfg, vcf), eid_map)
    if path.exists() and path.read_text() == text:
        return
    path.write_text(text)
    log(f"  wrote sample file: {path} ({text.count(chr(10)) - 2} samples, in VCF order)")
