#!/bin/bash
#SBATCH --job-name=sgdp38_stage
#SBATCH --partition=short
#SBATCH --account=sfglab
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH --time=10:00:00
#SBATCH --output=/mnt/evafs/groups/sfglab/data/tbooo/data/slurm/sgdp_grch38/stage_smoke_chr22.%j.log
# The first real run of `tbooo map sgdp38`, on chromosome 22, from a copy of the sv-panel branch
# in tools/sv-panel-src. It runs on a data directory of its own under data/tmp/sgdp38_smoke that
# links the raw files, so nothing of the instance is written. The sample ids are stand-ins made
# from the names in the merged file, since the instance's rename file predates the fix of its key.
set -euo pipefail
T=/mnt/evafs/groups/sfglab/data/tbooo
export PATH=$T/tools/env/bin:$PATH
D=$T/data; S=$D/tmp/sgdp38_smoke; N=${SLURM_CPUS_PER_TASK:-8}
mkdir -p $S/data/raw/sgdp $S/data/metadata $S/tmp
ln -sfn $D/raw/sgdp/pvcf $S/data/raw/sgdp/pvcf
ln -sfn $D/raw/sgdp/vcf $S/data/raw/sgdp/vcf
ln -sfn $D/raw/1kg $S/data/raw/1kg
bcftools query -l $D/raw/sgdp/pvcf/sgdp_c22.pvcf.gz | awk '{ print $1 "\t" 2000000 + NR - 1 }' > $S/data/metadata/vcf_sample_rename_sgdp.txt
sed -e "s|^data_dir:.*|data_dir: $S/data|" -e "s|^reference_dir:.*|reference_dir: $D/reference|" -e "s|^tmp_dir:.*|tmp_dir: $S/tmp|" \
    -e "s|^wgs_sgdp_merge_threads:.*|wgs_sgdp_merge_threads: $N|" $T/tools/sv-panel-src/config.yaml > $S/config.yaml
cd $T/tools/sv-panel-src
date +"=== %H:%M:%S start"
PYTHONPATH=src $T/.venv/bin/python -c "from tbooo.cli import cli; cli()" map sgdp38 -c $S/config.yaml --chroms 22,X 2>&1 | grep -vE "^\s*$|bcf_hdr_check_sanity|INFO\s+20|Progress|hts_idx_load3" | tail -60
date +"=== %H:%M:%S checks"
O=$S/data/raw/sgdp/grch38
ls -la $O; cat $O/sgdp_c22.counts.tsv
P=$O/sgdp_c22.phased.vcf.gz
echo "samples: $(bcftools query -l $P | wc -l), first $(bcftools query -l $P | head -2 | tr '\n' ' ')"
echo "records: $(bcftools index -n $P)   with an unphased or missing genotype: $(bcftools view -H -i 'GT~"/" | GT~"\."' $P | wc -l)"
echo "chromosomes: $(bcftools index -s $P | cut -f1,3 | tr '\n' ' ')"
echo "the trial's arm B wrote $(bcftools index -n $D/tmp/sgdp_grch38_phase/b/armB.vcf.gz) records with the shared people's relatives left out of the reference"
date +"=== %H:%M:%S a second call builds nothing"
PYTHONPATH=src $T/.venv/bin/python -c "from tbooo.cli import cli; cli()" map sgdp38 -c $S/config.yaml --chroms 22 2>&1 | tail -4
du -sh $S
