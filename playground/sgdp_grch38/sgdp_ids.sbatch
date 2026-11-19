#!/bin/bash
#SBATCH --job-name=sgdp_ids
#SBATCH --partition=short
#SBATCH --account=sfglab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/mnt/evafs/groups/sfglab/data/tbooo/data/slurm/sgdp_grch38/sgdp_ids.%j.log
# The SGDP sample table and id map keyed on the name inside the VCF, the Illumina id.
#
#     MODE=trial sbatch sgdp_ids.sbatch    builds both in data/tmp/sgdp_ids_trial and compares
#     MODE=apply sbatch sgdp_ids.sbatch    builds both in the instance, after copies of the old
#
# The instance's table and map of June are keyed on the ENA accession, which no VCF carries, so
# its rename file renames nothing. `tbooo download sgdp --no-vcf` moves the old table aside as
# sgdp_samples.igsr.tsv and `tbooo map eids` removes the old map and rename file before it writes
# the new ones. The apply mode therefore copies all three first, under the suffix
# .before-illumina-key, and stops if a copy is missing. Only the SGDP half of the id stage is
# called, so the 1000 Genomes map is not read or written.
#
# The login node cannot import numpy, so this is a job.
set -euo pipefail
T=/mnt/evafs/groups/sfglab/data/tbooo
D=$T/data
SRC=${SRC:-$T/tools/sv-panel-src}
MODE=${MODE:-trial}
export PATH=$T/tools/env/bin:$PATH
PY=$T/.venv/bin/python

if [ "$MODE" = trial ]; then
  S=$D/tmp/sgdp_ids_trial
  mkdir -p $S/data/raw/sgdp $S/data/metadata $S/tmp
  DATA=$S/data
elif [ "$MODE" = apply ]; then
  S=$D/tmp/sgdp_ids_apply
  mkdir -p $S/tmp
  DATA=$D
  for f in $D/metadata/eid_map_sgdp.tsv $D/metadata/vcf_sample_rename_sgdp.txt $D/raw/sgdp/sgdp_samples.tsv; do
    [ -e "$f" ] && cp -n -p "$f" "$f.before-illumina-key"
    [ ! -e "$f" ] || [ -s "$f.before-illumina-key" ] || { echo "no copy of $f, stopping"; exit 1; }
  done
  ls -la $D/metadata/*.before-illumina-key $D/raw/sgdp/*.before-illumina-key
else
  echo "MODE is trial or apply"; exit 1
fi
sed -e "s|^data_dir:.*|data_dir: $DATA|" -e "s|^reference_dir:.*|reference_dir: $D/reference|" -e "s|^tmp_dir:.*|tmp_dir: $S/tmp|" \
    $SRC/config.yaml > $S/config.yaml

cd $SRC
date +"=== %H:%M:%S $MODE, code $(md5sum src/tbooo/pipeline/eids.py | cut -c1-8)"
PYTHONPATH=src $PY - "$S/config.yaml" <<'EOF'
import sys
from tbooo.cli import _cfg
from tbooo.download.sgdp import download_metadata
from tbooo.pipeline.eids import _assign_sgdp
cfg = _cfg(sys.argv[1])
download_metadata(cfg)
_assign_sgdp(cfg)
EOF

date +"=== %H:%M:%S the new map against the old one and the names in the VCF"
# The old map is the copy in apply mode and the instance's own file in trial mode.
OLD=$D/metadata/eid_map_sgdp.tsv
[ "$MODE" = apply ] && OLD=$OLD.before-illumina-key
$PY - "$DATA" "$OLD" "$D/raw/sgdp/pvcf/sgdp_c22.pvcf.gz" <<'EOF'
import subprocess
import sys

import pandas as pd

data, old_path, pvcf = sys.argv[1:4]
new = pd.read_csv(f"{data}/metadata/eid_map_sgdp.tsv", sep="\t", dtype=str, keep_default_na=False)
old = pd.read_csv(old_path, sep="\t", dtype=str, keep_default_na=False)
rename = dict(line.split()[:2] for line in open(f"{data}/metadata/vcf_sample_rename_sgdp.txt") if line.strip())
vcf = subprocess.run(["bcftools", "query", "-l", pvcf], capture_output=True, text=True, check=True).stdout.split()
print(f"samples: new map {len(new)}, old map {len(old)}, names in the pVCF {len(vcf)}, rename lines {len(rename)}")
print(f"eids: {new.eid.iloc[0]} to {new.eid.iloc[-1]}, unique {new.eid.is_unique}, sample ids unique {new.sample_id.is_unique}")
missing = sorted(set(vcf) - set(rename))
extra = sorted(set(new.sample_id) - set(vcf))
print(f"pVCF names the rename file lacks: {len(missing)} {missing[:5]}")
print(f"map ids the pVCF lacks: {len(extra)} {extra[:5]}")
print(f"with a population: {int((new.population != 'Unknown').sum())}, with a sex: {int((new.sex != '0').sum())}, with an ENA accession: {int((new.ena_accession != '').sum())}")
print(f"ENA accessions on more than one row: {int(new.ena_accession[new.ena_accession != ''].duplicated().sum())}")
m = new.merge(old, on="ena_accession", how="outer", suffixes=("_new", "_old"), indicator=True)
print("by ENA accession:", m["_merge"].value_counts().to_dict())
both = m[m["_merge"] == "both"]
print(f"people of both maps whose eid stays: {int((both.eid_new == both.eid_old).sum())} of {len(both)}")
print(f"people of both maps whose population differs: {int((both.population_new != both.population_old).sum())}")
print(f"people of both maps whose sex differs: {int((both.sex_new != both.sex_old).sum())}")
print(new.head(3).to_string())
EOF
ls -la $DATA/metadata/eid_map_sgdp.tsv $DATA/metadata/vcf_sample_rename_sgdp.txt $DATA/raw/sgdp/sgdp_samples*.tsv*
date +"=== %H:%M:%S done"
