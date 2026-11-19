#!/bin/bash
#SBATCH --job-name=sgdp38_phase
#SBATCH --partition=short
#SBATCH --account=sfglab
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH --time=10:00:00
#SBATCH --output=/mnt/evafs/groups/sfglab/data/tbooo/data/slurm/sgdp_grch38/pilot_phase_chr22.%j.log
# A trial of the phasing step of the SGDP to GRCh38 stage on chromosome 22. Beagle with the panel
# as its reference writes the panel's sites alone, so a third of the SGDP records are lost. Three
# arms are compared on the people both projects sequenced, whose phase the panel knows:
#   R  the panel as the reference, the first trial's way
#   A  the panel with the SGDP records it lacks added as sites where every panel sample has the
#      reference allele, so that Beagle keeps them and phases them from the other SGDP carriers
#   N  no reference, the 277 samples phased among themselves
# A fifth of the shared positions is kept out of the reference and stands in for SGDP records the
# panel lacks, since the panel knows their phase. The shared people, their parents and their
# children are left out of the reference. Reads the first trial's outputs in
# data/tmp/sgdp_grch38_pilot and writes under data/tmp/sgdp_grch38_phase.
set -euo pipefail
T=/mnt/evafs/groups/sfglab/data/tbooo
export PATH=$T/tools/env/bin:$PATH
D=$T/data; P=$D/tmp/sgdp_grch38_pilot; W=$D/tmp/sgdp_grch38_phase; C=22; N=${SLURM_CPUS_PER_TASK:-8}
mkdir -p $W; cd $W
G=$P/4.grch38.vcf.gz
PANEL=$D/raw/1kg/1kGP_high_coverage_Illumina.chr$C.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
MAP=$D/reference/genetic_maps/beagle_GRCh38/chr_in_chrom_field/plink.chrchr$C.GRCh38.map
PED=$D/raw/1kg/20130606_g1k_3202_samples_ped_population.txt
PAIRS=$P/shared.samples
E=$T/playground/sgdp_grch38/phase_eval.py
step() { echo; echo "=== $(date +%H:%M:%S) $*"; }
count() { bcftools index -n "$1"; }

step "reference samples"
cut -f1 $PAIRS > shared.sgdp; cut -f2 $PAIRS > shared.panel
awk 'NR==FNR { s[$1] = 1; next } FNR > 1 { if (($2 in s) || ($3 in s) || ($4 in s)) x[$2] = 1; if ($2 in s) { x[$3] = 1; x[$4] = 1 } }
     END { for (k in x) if (k != "0") print k }' shared.panel $PED | sort > excluded.samples
bcftools query -l $PANEL | grep -v -x -F -f excluded.samples > ref.samples
echo "left out: $(wc -l < excluded.samples)   reference samples: $(wc -l < ref.samples)"

step "the SGDP input with one record per position"
bcftools norm -d all --threads $N -Oz -o g.vcf.gz $G; bcftools index -t g.vcf.gz
echo "records: $(count $G)   one per position: $(count g.vcf.gz)"
bcftools +fill-tags g.vcf.gz -Ou -- -t AC | bcftools query -f '%POS\t%REF\t%ALT\t%INFO/AC\n' > counts.tsv

step "held out positions and the two references"
bcftools isec -n=2 -w1 -c none -Oz -o shared.vcf.gz g.vcf.gz $PANEL --threads $N; bcftools index -t shared.vcf.gz
bcftools query -f '%CHROM\t%POS\n' shared.vcf.gz | awk '$2 % 5 == 0' > held.tsv
echo "shared records: $(count shared.vcf.gz)   held out: $(wc -l < held.tsv)"
bcftools view -S ref.samples -T ^held.tsv --threads $N -Oz -o refR.vcf.gz $PANEL; bcftools index -t refR.vcf.gz
echo "reference R records: $(count refR.vcf.gz)"
bcftools query -f '%POS\n' refR.vcf.gz | sort -u > refR.pos
bcftools isec -C -w1 -c none g.vcf.gz refR.vcf.gz --threads $N \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
  | awk 'NR==FNR { p[$1] = 1; next } { if ($2 in p) clash++; else print } END { print clash + 0 > "clash.count" }' refR.pos - > added.sites
echo "SGDP records the reference lacks: $(wc -l < added.sites)   at a position the reference has with other alleles, not added: $(cat clash.count)"
awk 'NR==FNR { h[$2] = 1; next } { if ($2 in h) a++; else b++ } END { print "of the added: " a + 0 " held out, " b + 0 " absent from the panel" }' held.tsv added.sites
awk -F'\t' 'NR==FNR { h[$2] = 1; next } FILENAME == "added.sites" { if (!($2 in h)) only[$2" "$3" "$4] = 1; next }
     (($1" "$2" "$3) in only) { n++; if ($4 <= 1) s++; else if ($4 <= 5) f++; else m++ }
     END { printf "records absent from the panel by allele count among SGDP: %d with 1, %d with 2 to 5, %d with more, of %d\n", s, f, m, n }' held.tsv added.sites counts.tsv
( bcftools view -h refR.vcf.gz
  awk -v n=$(wc -l < ref.samples) 'BEGIN { for (i = 0; i < n; i++) g = g "\t0|0" } { print $1 "\t" $2 "\t.\t" $3 "\t" $4 "\t.\tPASS\t.\tGT" g }' added.sites ) | bgzip -@ $N > added.vcf.gz
bcftools index -t added.vcf.gz
bcftools concat -a --threads $N -Oz -o refA.vcf.gz refR.vcf.gz added.vcf.gz; bcftools index -t refA.vcf.gz
echo "reference A records: $(count refA.vcf.gz)"

step "truth: the panel's phase of the shared people"
bcftools view -S shared.panel --regions-overlap pos -R <(bcftools query -f '%CHROM\t%POS\n' shared.vcf.gz) -Ob -o truth.bcf $PANEL
bcftools query -l truth.bcf | diff -q - shared.panel
bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' truth.bcf > truth.tsv

run_arm() {  # name, then the beagle arguments
  local a=$1; shift
  step "arm $a"
  beagle -Xmx56g gt=g.vcf.gz map=$MAP out=arm$a nthreads=$N "$@" 2>&1 | grep -E "markers|samples|time|ERROR|Exception|rror" | tail -12 || true
  if [ -s arm$a.vcf.gz ]; then
    bcftools index -f -t arm$a.vcf.gz
    bcftools view -S shared.sgdp -Ob -o arm$a.shared.bcf arm$a.vcf.gz
    bcftools query -l arm$a.shared.bcf | diff -q - shared.sgdp
    bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' arm$a.shared.bcf > arm$a.tsv
    python3 $E truth.tsv arm$a.tsv held.tsv counts.tsv $PAIRS
  else
    echo "  arm $a wrote no output"
  fi
}
run_arm R ref=refR.vcf.gz impute=false
run_arm A ref=refA.vcf.gz impute=false
run_arm N

step "does Beagle take several records at one position"
beagle -Xmx16g gt=$G map=$MAP out=dup chrom=chr$C:20000000-24000000 nthreads=$N 2>&1 | grep -E "markers|ERROR|Exception|rror|uplicate" | tail -5 || true
[ -s dup.vcf.gz ] && { bcftools index -f -t dup.vcf.gz; echo "input records in the window: $(bcftools view -H -r chr$C:20000000-24000000 $G | wc -l)   output: $(count dup.vcf.gz)"; }

step "done"; du -sh $W
