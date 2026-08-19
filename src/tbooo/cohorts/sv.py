"""Emit a SLURM batch script that builds the pseudonymized SV layer.

This module deliberately executes nothing. The login node is an SSH jumpbox for
scheduling jobs; all processing belongs on compute nodes. `emit_sv_script` writes
a script and returns its path — the caller prints the submit command.

The source VCFs are read-only inputs. Every write goes into the instance tree.
"""

from __future__ import annotations

import re
from pathlib import Path

from tbooo.cohorts.eids import rename_path
from tbooo.cohorts.registry import BLOCK_SIZE, Cohort, Sample
from tbooo.config import Config
from tbooo.utils import ensure_dirs, eid_prefix_dir, log

# Per-sample SV outputs, and the index files that accompany them.
_SV_OUTPUT_RE = re.compile(r"^(?P<eid>\d+)_sv_0_0\.vcf\.gz(?:\.tbi|\.csi)?$")

# Characters that make a path impossible to interpolate safely into the
# double-quoted shell words this script is built from. A newline ends the line —
# and therefore any `#` comment it sits in — so the remainder of a directory name
# becomes a LIVE command in a script the operator is told to `sbatch`; a double
# quote closes the quoting early to the same effect. These names arrive from
# `iterdir()`, so no amount of care in the caller's input guards them. A path
# like that inside this tree is a red flag in its own right, so it is refused
# loudly rather than quietly mangled.
_UNSAFE_PATH_CHARS = {"\n": "a newline", "\r": "a carriage return", '"': "a double-quote"}


def _shell_safe(path: Path) -> str:
    """Return `path` as a string, refusing anything unquotable. See above."""
    text = str(path)
    for char, name in _UNSAFE_PATH_CHARS.items():
        if char in text:
            raise RuntimeError(
                f"refusing to name a path containing {name} in a generated script "
                f"or log line: {text!r}. The sv stage interpolates paths into a "
                f"shell script an operator is told to submit, where such a "
                f"character breaks out of the quoted path and can turn a "
                f"commented-out cleanup line into a live command. Rename the "
                f"offending directory or file, then re-run."
            )
    return text


def sv_script_path(cfg: Config, cohort: Cohort) -> Path:
    return cfg.slurm_dir() / f"sv_{cohort.name}.sbatch"


def sv_output_path(cfg: Config, eid: int) -> Path:
    return cfg.sv_dir() / eid_prefix_dir(eid) / f"{eid}_sv_0_0.vcf.gz"


def _read_rename_map(path: Path) -> dict[str, str]:
    """Parse a bcftools reheader --samples file into {original_id: eid}."""
    mapping: dict[str, str] = {}
    for lineno, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            raise RuntimeError(
                f"{path}:{lineno}: expected '<sample_id>\\t<eid>', got {line!r}. "
                f"Re-run `tbooo restricted eids` to regenerate it."
            )
        mapping[parts[0].strip()] = parts[1].strip()
    return mapping


def _describe_drift(on_disk: dict[str, str], expected: dict[str, str]) -> str:
    """A short, specific summary of how the two sample -> EID maps differ."""
    missing = sorted(set(expected) - set(on_disk))
    extra = sorted(set(on_disk) - set(expected))
    changed = sorted(s for s in set(on_disk) & set(expected) if on_disk[s] != expected[s])
    parts = []
    if missing:
        parts.append(f"{len(missing)} sample(s) absent from the file (e.g. {missing[0]!r})")
    if extra:
        parts.append(f"{len(extra)} sample(s) in the file that no longer exist (e.g. {extra[0]!r})")
    if changed:
        first = changed[0]
        parts.append(
            f"{len(changed)} sample(s) mapped to a different EID (e.g. {first!r}: "
            f"file says {on_disk[first]}, current assignment is {expected[first]})"
        )
    return "; ".join(parts) or "the two maps differ"


def check_rename_file(cfg: Config, cohort: Cohort, eids: dict[str, int]) -> Path:
    """Verify the on-disk rename map exists and matches the EIDs we were handed.

    The sv stage runs AFTER eids. If it does not, the script interpolates a path
    that is not there and the job dies on the compute node under `set -e`. Worse,
    if the sample set changed since eids last ran, the rename map still holds the
    OLD EIDs while the output paths use the NEW ones — producing an
    `<newEID>_sv_0_0.vcf.gz` whose internal sample column names a different EID.
    That is a pseudonymization correctness bug, so it is checked, not assumed.
    """
    rename = rename_path(cfg, cohort)
    if not rename.exists():
        raise RuntimeError(
            f"cohort {cohort.name!r}: the sample rename map {rename} does not exist. "
            f"The sv stage must run after the eids stage — run "
            f"`tbooo restricted eids --only {cohort.name}` first."
        )
    on_disk = _read_rename_map(rename)
    expected = {sid: str(eid) for sid, eid in eids.items()}
    if on_disk != expected:
        raise RuntimeError(
            f"cohort {cohort.name!r}: {rename} disagrees with the EIDs currently "
            f"assigned — {_describe_drift(on_disk, expected)}. The sample set has "
            f"changed since the eids stage last ran; emitting now would write SV "
            f"files named for one EID whose sample column says another. Re-run "
            f"`tbooo restricted eids --only {cohort.name}` first."
        )
    return rename


def stale_sv_outputs(cfg: Config, cohort: Cohort, eids: dict[str, int]) -> list[Path]:
    """SV outputs inside this cohort's EID block whose EID is no longer assigned.

    `link._prune_orphaned_links` cannot cover these: it scans only the link layers
    and only ever removes symlinks, while these are REAL files. Renumbering
    therefore strands a VCF under a retired EID — the very misattribution the
    prune exists to prevent, in the one layer holding real data. They are reported
    rather than deleted; see `emit_sv_script`.

    As with the link prune, only EIDs inside [eid_base, eid_base + BLOCK_SIZE) are
    ever considered, because a 2-digit prefix directory is a 100,000-wide bucket
    and can be shared with another cohort. Every prefix directory is scanned, not
    just those of the current EIDs: a retired EID can sit in a bucket that the
    current assignment no longer reaches.

    The scan never leaves the instance tree. `is_dir()` follows symlinks, so a
    prefix directory symlinked elsewhere would otherwise list real files OUTSIDE
    `sv_dir` and report them for removal — irreplaceable data an operator might
    delete by pasting the suggested command. Symlinked prefix directories are
    skipped outright, and every reported item's RESOLVED path must still sit
    inside the resolved `sv_dir`. Anything failing either check is skipped, never
    reported.
    """
    sv_dir = cfg.sv_dir()
    if not sv_dir.is_dir():
        return []
    root = sv_dir.resolve()
    lo, hi = cohort.eid_base, cohort.eid_base + BLOCK_SIZE
    current = set(eids.values())
    stale: list[Path] = []
    for prefix_dir in sorted(
        p for p in sv_dir.iterdir() if p.is_dir() and not p.is_symlink()
    ):
        for item in sorted(prefix_dir.iterdir()):
            match = _SV_OUTPUT_RE.match(item.name)
            if not match:
                continue
            eid = int(match.group("eid"))
            if not (lo <= eid < hi and eid not in current):
                continue
            try:
                resolved = item.resolve()
            except OSError:  # pragma: no cover — unresolvable link, treat as outside
                continue
            if not resolved.is_relative_to(root):
                continue
            stale.append(item)
    return stale


def _warn_stale(cohort: Cohort, stale: list[Path]) -> None:
    # Validated before the first log line, so an unquotable path is refused
    # whole rather than half-printed.
    names = [_shell_safe(path) for path in stale]
    log(
        f"  WARN: {cohort.name}: {len(stale)} SV output file(s) carry an EID this "
        f"cohort no longer assigns (the sample set changed since they were written):"
    )
    for name in names:
        log(f"  WARN:   {name}")
    log(
        "  WARN: these are real files, not symlinks, so nothing deletes them "
        "automatically. Review them; the generated script lists the removal "
        "commands as comments."
    )


def merged_output_path(cfg: Config, cohort: Cohort) -> Path:
    return cfg.sv_dir() / f"consensusv_{cohort.name}.vcf.gz"


def emit_sv_script(
    cfg: Config,
    cohort: Cohort,
    samples: list[Sample],
    eids: dict[str, int],
    *,
    partition: str = "short",
    time_limit: str = "01:00:00",
) -> Path | None:
    """Write the sbatch script for `cohort`. Returns None if it has no SV calls."""
    with_sv = [
        (s, next((f for f in s.files if f.modality == "sv_vcf"), None))
        for s in sorted(samples, key=lambda s: s.sample_id)
    ]
    with_sv = [(s, f) for s, f in with_sv if f is not None]
    stale = stale_sv_outputs(cfg, cohort, eids)
    if stale:
        _warn_stale(cohort, stale)

    if not with_sv:
        log(f"  {cohort.name}: no SV calls — nothing to emit")
        return None

    rename = check_rename_file(cfg, cohort, eids)
    ensure_dirs(cfg.slurm_dir(), cfg.sv_dir())
    bcftools = cfg.tools.bcftools
    merged = merged_output_path(cfg, cohort)
    file_list = cfg.slurm_dir() / f"sv_{cohort.name}_files.txt"

    lines: list[str] = [
        "#!/bin/bash",
        f"#SBATCH --job-name=tbooo_sv_{cohort.name}",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --time={time_limit}",
        "#SBATCH --cpus-per-task=4",
        "#SBATCH --mem=8G",
        f"#SBATCH --output={cfg.slurm_dir()}/sv_{cohort.name}_%j.out",
        "",
        "set -euo pipefail",
        "",
        "# Generated by `tbooo restricted sv`. Do not edit by hand.",
        "# Reads source VCFs read-only; every write goes into the instance tree.",
        "",
    ]

    if stale:
        lines += [
            "# WARNING: the following SV output files carry an EID this cohort no",
            "# longer assigns, left behind when the sample set changed. They are real",
            "# files holding real data, so nothing here removes them. Review them and,",
            "# if they are genuinely obsolete, run these by hand:",
            *[f'#   rm -f "{_shell_safe(path)}"' for path in stale],
            "",
        ]

    lines += [
        f': > "{_shell_safe(file_list)}"',
        "",
    ]

    for sample, src in with_sv:
        eid = eids[sample.sample_id]
        out = sv_output_path(cfg, eid)
        lines += [
            f'mkdir -p "{_shell_safe(out.parent)}"',
            f'{bcftools} reheader --samples "{_shell_safe(rename)}" '
            f'"{_shell_safe(src.path)}" | bgzip -c > "{_shell_safe(out)}"',
            f'tabix -p vcf "{_shell_safe(out)}"',
            f'echo "{_shell_safe(out)}" >> "{_shell_safe(file_list)}"',
            "",
        ]

    lines += [
        "# Merge into one multi-sample SV pVCF for the cohort.",
        "# NOTE: this is a convenience view, not a joint callset — breakpoint",
        "# imprecision leaves near-duplicate records across samples.",
        f"{bcftools} merge --output-type z "
        f'--file-list "{_shell_safe(file_list)}" --output "{_shell_safe(merged)}"',
        f'tabix -p vcf "{_shell_safe(merged)}"',
        "",
        f'echo "done: {len(with_sv)} samples -> {_shell_safe(merged)}"',
        "",
    ]

    script = sv_script_path(cfg, cohort)
    # The submit line below is meant to be pasted into a shell, so it gets the
    # same treatment as the script's own contents.
    safe_script = _shell_safe(script)
    script.write_text("\n".join(lines))
    script.chmod(0o755)
    log(f"  {cohort.name}: wrote {safe_script} ({len(with_sv)} samples)")
    log(f"  submit with:  sbatch {safe_script}")
    return script
