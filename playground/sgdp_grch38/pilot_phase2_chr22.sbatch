#!/bin/bash
#SBATCH --job-name=sgdp38_phase2
#SBATCH --partition=short
#SBATCH --account=sfglab
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH --time=10:00:00
#SBATCH --output=/mnt/evafs/groups/sfglab/data/tbooo/data/slurm/sgdp_grch38/pilot_phase2_chr22.%j.log
# The second part of the phasing trial. In the first part the added markers stood in for held
# out shared sites, most of them common, and the switch rate over the shared sites doubled. The
# stage adds the SGDP records the panel lacks, four in five of them carried by one sample. This
# part measures that case: every shared site stays in the reference, and the switch rate over
# them is compared with and without the added markers.
#   R0  the panel without the shared people and their relatives
#   B   the same with the SGDP records the panel lacks added
# Reads the outputs of pilot_phase_chr22.sbatch and writes under its directory, in b/.
set -euo pipefail
T=/mnt/evafs/groups/sfglab/data/tbooo
export PATH=$T/tools/env/bin:$PATH
D=$T/data; W=$D/tmp/sgdp_grch38_phase; C=22; N=${SLURM_CPUS_PER_TASK:-8}
cd $W; mkdir -p b
PANEL=$D/raw/1kg/1kGP_high_coverage_Illumina.chr$C.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
MAP=$D/reference/genetic_maps/beagle_GRCh38/chr_in_chrom_field/plink.chrchr$C.GRCh38.map
PAIRS=$D/tmp/sgdp_grch38_pilot/shared.samples
E=$T/playground/sgdp_grch38/phase_eval.py
step() { echo; echo "=== $(date +%H:%M:%S) $*"; }
count() { bcftools index -n "$1"; }

step "references"
bcftools view -S ref.samples --threads $N -Oz -o b/ref0.vcf.gz $PANEL; bcftools index -t b/ref0.vcf.gz
awk 'NR==FNR { h[$2] = 1; next } !($2 in h)' held.tsv added.sites > b/only.sites
( bcftools view -h b/ref0.vcf.gz
  awk -v n=$(wc -l < ref.samples) 'BEGIN { for (i = 0; i < n; i++) g = g "\t0|0" } { print $1 "\t" $2 "\t.\t" $3 "\t" $4 "\t.\tPASS\t.\tGT" g }' b/only.sites ) | bgzip -@ $N > b/only.vcf.gz
bcftools index -t b/only.vcf.gz
bcftools concat -a --threads $N -Oz -o b/refB.vcf.gz b/ref0.vcf.gz b/only.vcf.gz; bcftools index -t b/refB.vcf.gz
echo "reference R0 records: $(count b/ref0.vcf.gz)   added: $(wc -l < b/only.sites)   reference B records: $(count b/refB.vcf.gz)"
: > b/none.tsv

run_arm() {
  local a=$1; shift
  step "arm $a"
  beagle -Xmx56g gt=g.vcf.gz map=$MAP out=b/arm$a nthreads=$N impute=false "$@" 2>&1 | grep -E "Study     markers|time" | tail -4 || true
  bcftools index -f -t b/arm$a.vcf.gz
  bcftools view -S shared.sgdp -Ob -o b/arm$a.shared.bcf b/arm$a.vcf.gz
  bcftools query -l b/arm$a.shared.bcf | diff -q - shared.sgdp
  bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' b/arm$a.shared.bcf > b/arm$a.tsv
  python3 $E truth.tsv b/arm$a.tsv b/none.tsv counts.tsv $PAIRS
}
run_arm R0 ref=b/ref0.vcf.gz
run_arm B ref=b/refB.vcf.gz
step "done"
