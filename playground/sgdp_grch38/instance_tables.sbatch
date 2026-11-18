#!/bin/bash
#SBATCH --job-name=tbooo_tables
#SBATCH --partition=short
#SBATCH --account=sfglab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/mnt/evafs/groups/sfglab/data/tbooo/data/slurm/sgdp_grch38/instance_tables.%j.log
# The two tables of the instance that carry SGDP ids, built again after the id map changed:
# metadata/ukb_sqc_v2.txt with ukb_rel.txt from `tbooo map qc` and Showcase/participant.parquet
# from `tbooo map phenotypes`.
#
# Both stages skip an output that is there, so the old files are moved aside first under the
# suffix .before-illumina-key. A file already moved is left where it is. plink2 and king are
# in an environment of their own, tools/env-qc, so the stage's environment is not changed.
set -euo pipefail
T=/mnt/evafs/groups/sfglab/data/tbooo
D=$T/data
SRC=${SRC:?the checkout of the sv-panel branch to run}
N=${SLURM_CPUS_PER_TASK:-8}
export PATH=$T/tools/env/bin:$T/tools/env-qc/bin:$PATH
for tool in plink2 king bcftools; do command -v $tool > /dev/null || { echo "$tool is not on the path"; exit 1; }; done
S=$D/tmp/instance_tables
mkdir -p $S
sed -e "s|^data_dir:.*|data_dir: $D|" -e "s|^reference_dir:.*|reference_dir: $D/reference|" -e "s|^tmp_dir:.*|tmp_dir: $S|" $SRC/config.yaml > $S/config.yaml
for f in "$D/metadata/ukb_sqc_v2.txt" "$D/metadata/ukb_rel.txt" "$D/Showcase/participant.parquet"; do
  if [ -e "$f" ] && [ ! -e "$f.before-illumina-key" ]; then mv -n "$f" "$f.before-illumina-key"; fi
done
ls -la $D/metadata/*.before-illumina-key $D/Showcase/*.before-illumina-key
cd $SRC
run() { PYTHONPATH=src $T/.venv/bin/python -c "from tbooo.cli import cli; cli()" "$@" -c $S/config.yaml 2>&1 | grep -vE "^\s*$" | tail -40; }
date +"=== %H:%M:%S qc"
run map qc
date +"=== %H:%M:%S phenotypes"
run map phenotypes
date +"=== %H:%M:%S the new tables against the old"
$T/.venv/bin/python - "$D" <<'PY'
import sys

import pandas as pd

d = sys.argv[1]
ids = set(pd.read_csv(f"{d}/metadata/eid_map_sgdp.tsv", sep="\t").eid)
for name in ("ukb_sqc_v2.txt", "ukb_rel.txt"):
    for tag, path in (("new", f"{d}/metadata/{name}"), ("old", f"{d}/metadata/{name}.before-illumina-key")):
        t = pd.read_csv(path, sep=r"\s+", header=None, low_memory=False)
        first = pd.to_numeric(t[0], errors="coerce")
        print(f"{name} {tag}: {len(t)} rows, {t.shape[1]} columns, {int(first.between(2_000_000, 2_999_999).sum())} with an SGDP id, {int(first.isin(ids).sum())} in the id map")
for tag, path in (("new", f"{d}/Showcase/participant.parquet"), ("old", f"{d}/Showcase/participant.parquet.before-illumina-key")):
    p = pd.read_parquet(path)
    sg = p[p.eid.between(2_000_000, 2_999_999)]
    print(f"participant {tag}: {len(p)} rows, {len(p.columns)} columns, {len(sg)} SGDP rows, {int(sg.eid.isin(ids).sum())} in the id map")
new = pd.read_parquet(f"{d}/Showcase/participant.parquet")
old = pd.read_parquet(f"{d}/Showcase/participant.parquet.before-illumina-key")
print("columns only in the old table:", sorted(set(old.columns) - set(new.columns))[:12])
print("columns only in the new table:", sorted(set(new.columns) - set(old.columns))[:12])
kg = new[new.eid < 2_000_000].merge(old[old.eid < 2_000_000], on="eid", suffixes=("_new", "_old"))
print(f"1000 Genomes rows in both: {len(kg)} of {int((new.eid < 2_000_000).sum())} new and {int((old.eid < 2_000_000).sum())} old")
PY
date +"=== %H:%M:%S done"
