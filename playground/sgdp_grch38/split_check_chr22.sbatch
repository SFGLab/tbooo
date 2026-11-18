#!/bin/bash
#SBATCH --job-name=split_check
#SBATCH --partition=short
#SBATCH --account=sfglab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/mnt/evafs/groups/sfglab/data/tbooo/data/slurm/sgdp_grch38/split_check_chr22.%j.log
# A check of the per sample split of the 1KGP panel on real data, the options of wgs.py's
# _SPLIT_CARRIED_ONLY. Twenty samples of chromosome 22: each file must hold the records its
# sample carries and no other, and the size of all 3,202 files follows from the twenty.
set -euo pipefail
T=/mnt/evafs/groups/sfglab/data/tbooo
export PATH=$T/tools/env/bin:$PATH
D=$T/data; W=$D/tmp/split_check_chr22
PANEL=$D/raw/1kg/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
mkdir -p $W/carried $W/all; cd $W
bcftools query -l $PANEL | awk 'NR % 160 == 1' | head -20 > samples.txt; wc -l < samples.txt
echo "=== $(date +%H:%M:%S) split with the filter"
bcftools +split $PANEL --include 'GT="alt"' --keep-tags FMT/GT --output-type z --output $W/carried --samples-file samples.txt
echo "=== $(date +%H:%M:%S) split without it, two samples"
head -2 samples.txt > two.txt
bcftools +split $PANEL --output-type z --output $W/all --samples-file two.txt
echo "=== $(date +%H:%M:%S) per sample"
total=$(bcftools index -n $PANEL); echo "panel records: $total"
bad=0
while read -r s; do
  f=carried/$s.vcf.gz
  n=$(bcftools view -H $f | wc -l)
  want=$(bcftools view -H -s $s -i 'GT="alt"' $PANEL | wc -l)
  ref=$(bcftools view -H -e 'GT="alt"' $f | wc -l)
  echo "$s records $n, carried in the panel $want, without the alternate allele $ref, bytes $(stat -c %s $f)"
  [ "$n" = "$want" ] && [ "$ref" = 0 ] || bad=$((bad+1))
done < samples.txt
echo "samples that differ: $bad"
bcftools view -h carried/$(head -1 samples.txt).vcf.gz | grep -c "^##INFO" || true
bcftools view -H carried/$(head -1 samples.txt).vcf.gz | head -2 | cut -c1-160
a=$(du -sb carried | cut -f1); b=$(du -sb all | cut -f1)
echo "20 filtered files: $a bytes, so 3,202 files of chromosome 22: $(( a / 20 * 3202 / 1000000000 )) GB"
echo "2 unfiltered files: $b bytes, so 3,202 files of chromosome 22: $(( b / 2 * 3202 / 1000000000 )) GB"
echo "chromosome 22 holds $(stat -c %s $PANEL) of the panel's bytes"
