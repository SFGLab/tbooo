#!/bin/bash
#SBATCH --job-name=sgdp38
#SBATCH --partition=short
#SBATCH --account=sfglab
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --array=1-22%8
#SBATCH --output=/mnt/evafs/groups/sfglab/data/tbooo/data/slurm/sgdp_grch38/sgdp38_genome.%A_%a.log
# `tbooo map sgdp38` on the instance, one autosome a task, into data/raw/sgdp/grch38.
#
#     sbatch --dependency=afterok:<panel job> sgdp38_genome.sbatch
#
# Phasing reads the chromosome's file of the 30x panel, so a task first reads the panel job's
# log and stops unless every file passed its md5. The work of a chromosome is in its own
# directory under data/tmp and its outputs carry its number, so tasks do not meet. The stage
# skips an output that is whole, so a task sent again goes on.
set -euo pipefail
T=/mnt/evafs/groups/sfglab/data/tbooo
D=$T/data
SRC=${SRC:?the checkout of the sv-panel branch to run}
PANEL_LOG=${PANEL_LOG:?the log of the panel download}
C=$SLURM_ARRAY_TASK_ID
N=${SLURM_CPUS_PER_TASK:-8}
export PATH=$T/tools/env/bin:$PATH
grep -q "^files: 46, mismatches: 0$" "$PANEL_LOG" || { echo "the panel is not whole: $(grep "^files:\|MISMATCH\|FAILED" "$PANEL_LOG" | head -5)"; exit 1; }
S=$D/tmp/sgdp38_genome
mkdir -p $S
sed -e "s|^data_dir:.*|data_dir: $D|" -e "s|^reference_dir:.*|reference_dir: $D/reference|" -e "s|^tmp_dir:.*|tmp_dir: $D/tmp|" \
    -e "s|^wgs_sgdp_merge_threads:.*|wgs_sgdp_merge_threads: $N|" $SRC/config.yaml > $S/config.$C.yaml
cd $SRC
date +"=== %H:%M:%S chr$C, code $(md5sum src/tbooo/pipeline/sgdp38.py | cut -c1-8)"
PYTHONPATH=src $T/.venv/bin/python -c "from tbooo.cli import cli; cli()" map sgdp38 -c $S/config.$C.yaml --chroms $C 2>&1 \
    | grep -vE "^\s*$|bcf_hdr_check_sanity|INFO\s+20|Progress|hts_idx_load3" | tail -60
date +"=== %H:%M:%S checks"
O=$D/raw/sgdp/grch38
P=$O/sgdp_c$C.phased.vcf.gz
cat $O/sgdp_c$C.counts.tsv
echo "samples: $(bcftools query -l $P | wc -l), first $(bcftools query -l $P | head -2 | tr '\n' ' ')"
echo "records: $(bcftools index -n $P)   with an unphased or missing genotype: $(bcftools view -H -i 'GT~"/" | GT~"\."' $P | wc -l)"
echo "chromosomes: $(bcftools index -s $P | cut -f1,3 | tr '\n' ' ')"
ls -la $O/sgdp_c$C.*
date +"=== %H:%M:%S done"
