"""Phase of the held out sites of one phasing arm against the panel's phase of the same people.

    python3 phase_eval.py <truth.tsv> <arm.tsv> <held.tsv> <counts.tsv> <pairs.tsv>

truth.tsv and arm.tsv are `bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n'` outputs with the
samples in the order of pairs.tsv, the panel's file for truth and the arm's for the estimate.
held.tsv lists the positions kept out of the reference. counts.tsv gives the count of the
alternate allele among the SGDP samples per record. pairs.tsv has the SGDP name and the panel
name of each person both projects sequenced.

A site counts for a person where both files have it heterozygous and phased. Two numbers per
arm. The switch rate over the scaffold sites, the ones the reference kept: the share of
neighbouring pairs whose relative phase differs from the panel's. And for each held out site the
agreement of its phase, relative to the scaffold site before it, with the panel's. Half is
chance. The second number is split by the allele count among the SGDP samples, since a variant
one sample carries has no other carrier to be phased from.
"""

from __future__ import annotations

import sys
from collections import defaultdict


def read(path: str) -> dict[tuple[int, str, str], list[str]]:
    out: dict[tuple[int, str, str], list[str]] = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[(int(f[0]), f[1], f[2])] = f[3:]
    return out


def side(gt: str) -> int:
    """0 for 0|1, 1 for 1|0, -1 for anything that is not a phased heterozygote."""
    return {"0|1": 0, "1|0": 1}.get(gt, -1)


def band(ac: int) -> str:
    return "1" if ac <= 1 else "2 to 5" if ac <= 5 else "6 to 20" if ac <= 20 else "over 20"


def main() -> None:
    truth, arm = read(sys.argv[1]), read(sys.argv[2])
    held = {int(line.split("\t")[1]) for line in open(sys.argv[3]) if line.strip()}
    counts: dict[tuple[int, str, str], int] = {}
    with open(sys.argv[4]) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            counts[(int(f[0]), f[1], f[2])] = int(f[3])
    people = sum(1 for line in open(sys.argv[5]) if line.strip())
    keys = sorted(k for k in truth if k in arm)
    in_arm_held = sum(1 for k in arm if k[0] in held)
    pairs = switches = 0
    agree: dict[str, int] = defaultdict(int)
    total: dict[str, int] = defaultdict(int)
    for i in range(people):
        last = None
        for k in keys:
            t, e = side(truth[k][i]), side(arm[k][i])
            if t < 0 or e < 0:
                continue
            rel = t ^ e
            if k[0] in held:
                if last is not None:
                    b = band(counts.get(k, 0))
                    total[b] += 1
                    agree[b] += rel == last
                continue
            if last is not None:
                pairs += 1
                switches += rel != last
            last = rel
    print(f"  records of the arm: {len(arm)}, of them at held out positions: {in_arm_held}")
    print(f"  scaffold: {pairs} neighbouring pairs, switch rate {switches / max(pairs, 1):.4f}")
    for b in ("1", "2 to 5", "6 to 20", "over 20"):
        if total[b]:
            print(f"  held out, allele count {b}: {total[b]} sites, {agree[b] / total[b]:.4f} agree")
    n = sum(total.values())
    if n:
        print(f"  held out, all: {n} sites, {sum(agree.values()) / n:.4f} agree")


if __name__ == "__main__":
    main()
