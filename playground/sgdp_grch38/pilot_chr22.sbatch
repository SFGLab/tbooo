#!/bin/bash
#SBATCH --job-name=sgdp38_pilot
#SBATCH --partition=short
#SBATCH --account=sfglab
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH --time=10:00:00
#SBATCH --output=/mnt/evafs/groups/sfglab/data/tbooo/data/slurm/sgdp_grch38/pilot_chr22.%j.log
# A trial of the SGDP to GRCh38 stage on chromosome 22. It reads the merged all-site file and the
# references, and writes under data/tmp/sgdp_grch38_pilot and three new reference downloads.
set -euo pipefail
T=/mnt/evafs/groups/sfglab/data/tbooo
export PATH=$T/tools/env/bin:$PATH
D=$T/data; W=$D/tmp/sgdp_grch38_pilot; C=22; N=${SLURM_CPUS_PER_TASK:-8}
mkdir -p $W $D/reference/liftover $D/reference/genetic_maps/beagle_GRCh38
cd $W
step() { echo; echo "=== $(date +%H:%M:%S) $*"; }
count() { bcftools index -n "$1"; }

step "downloads"
PANEL=$D/raw/1kg/1kGP_high_coverage_Illumina.chr$C.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
B=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
[ -s $PANEL ] || wget -q -O $PANEL.part $B/$(basename $PANEL) && [ -s $PANEL ] || mv $PANEL.part $PANEL
[ -s $PANEL.tbi ] || wget -q -O $PANEL.tbi $B/$(basename $PANEL).tbi
CHAIN=$D/reference/liftover/hg19ToHg38.over.chain.gz
[ -s $CHAIN ] || wget -q -O $CHAIN https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
MAPDIR=$D/reference/genetic_maps/beagle_GRCh38
[ -s $MAPDIR/chr_in_chrom_field/plink.chrchr$C.GRCh38.map ] || { wget -q -O $MAPDIR/plink.GRCh38.map.zip https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip && python3 -c "import zipfile; zipfile.ZipFile('$MAPDIR/plink.GRCh38.map.zip').extractall('$MAPDIR')"; }
ls -la $PANEL $CHAIN $MAPDIR | head -8; MAP=$MAPDIR/chr_in_chrom_field/plink.chrchr$C.GRCh38.map; head -2 $MAP

step "the all-site file"
SRC=$D/raw/sgdp/pvcf/sgdp_c$C.pvcf.gz
bcftools view -h $SRC | grep -E "^##(FORMAT|contig=<ID=$C,|reference)" | head -8
bcftools query -l $SRC | wc -l

step "variant sites, low filter levels set to missing"
# FL is the project's filter level, 0 to 9 or N. Genotypes at N or 0 are not calls.
bcftools view --threads $N -c 1 -Ou $SRC \
  | bcftools +setGT -Ou -- -t q -n . -i 'FMT/FL="N" | FMT/FL="0"' \
  | bcftools annotate -x INFO,^FORMAT/GT -Ou \
  | bcftools view -c 1 -Ob -o 1.variants.bcf --threads $N
bcftools index --threads $N 1.variants.bcf; echo "variant records: $(count 1.variants.bcf)"

step "normalise on GRCh37, rename the chromosome"
echo "$C chr$C" > chr.map
bcftools norm --threads $N -m -any -f $D/reference/GRCh37/human_g1k_v37.fasta --check-ref w -Ou 1.variants.bcf 2> norm.log \
  | bcftools annotate --rename-chrs chr.map -Oz -o 2.norm.vcf.gz
bcftools index -t 2.norm.vcf.gz; tail -2 norm.log; echo "after norm: $(count 2.norm.vcf.gz)"

step "lift to GRCh38"
REF38=$D/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
[ -s ${REF38%.fa}.dict ] || picard -Xmx8g CreateSequenceDictionary R=$REF38 O=${REF38%.fa}.dict 2>&1 | tail -1
picard -Xmx48g LiftoverVcf I=2.norm.vcf.gz O=3.lifted.vcf.gz CHAIN=$CHAIN REJECT=3.rejected.vcf.gz R=$REF38 \
  RECOVER_SWAPPED_REF_ALT=true WARN_ON_MISSING_CONTIG=true MAX_RECORDS_IN_RAM=200000 2>&1 | grep -E "INFO.*(Processed|lifted|failed|rejected)|Exception|ERROR" | tail -8
echo "lifted: $(count 3.lifted.vcf.gz)   rejected: $(bcftools view -H 3.rejected.vcf.gz | wc -l)"
bcftools view -H 3.rejected.vcf.gz | cut -f7 | sort | uniq -c | sort -rn | head -6
bcftools view -H 3.lifted.vcf.gz | cut -f1 | sort | uniq -c | sort -rn | head -5

step "this chromosome only, sorted, duplicates removed"
bcftools view -r chr$C -Ou 3.lifted.vcf.gz | bcftools sort -T $W/sort -Ou | bcftools norm -d exact -Oz -o 4.grch38.vcf.gz
bcftools index -t 4.grch38.vcf.gz; echo "GRCh38 chr$C records: $(count 4.grch38.vcf.gz)"

step "sites shared with the panel"
bcftools isec -n=2 -w1 -c none -Oz -o 5.shared.vcf.gz 4.grch38.vcf.gz $PANEL --threads $N; bcftools index -t 5.shared.vcf.gz
echo "panel records: $(count $PANEL)   shared: $(count 5.shared.vcf.gz)"
echo "panel symbolic records: $(bcftools view -H -i 'ALT~"<"' $PANEL | wc -l)"

step "samples the panel also has"
bcftools query -l $PANEL > panel.samples
awk 'NR==FNR {p[$1]=1; next} ($2 in p) {print $1"\t"$2}' panel.samples $T/playground/sgdp_grch38/pilot_aliases.tsv > shared.samples; wc -l < shared.samples; head -3 shared.samples
if [ -s shared.samples ]; then
  L=$(head -1 shared.samples | cut -f1); K=$(head -1 shared.samples | cut -f2)
  bcftools query -s $L -f '%POS\t%REF\t%ALT[\t%GT]\n' 5.shared.vcf.gz > a.gt
  bcftools query -s $K -R <(bcftools query -f '%CHROM\t%POS\n' 5.shared.vcf.gz) -f '%POS\t%REF\t%ALT[\t%GT]\n' $PANEL > b.gt
  awk -F'\t' 'function dose(g,  n,a) { if (g ~ /\./) return -1; n = split(g, a, /[\/|]/); return (a[1] != "0") + (a[2] != "0") }
    NR==FNR { k[$1" "$2" "$3] = dose($4); next } (($1" "$2" "$3) in k) { x = k[$1" "$2" "$3]; y = dose($4); if (x >= 0 && y >= 0) { n++; if (x == y) same++; else if (x > 0 || y > 0) { } ; if (x > 0 || y > 0) { nv++; if (x == y) sv++ } } }
    END { printf "%s against %s: %d sites compared, %.4f agree; where either carries the variant %d sites, %.4f agree\n", "'$L'", "'$K'", n, same / n, nv, sv / nv }' a.gt b.gt
fi

step "phase and impute against the panel"
beagle -Xmx56g gt=4.grch38.vcf.gz ref=$PANEL map=$MAP out=6.beagle nthreads=$N impute=true 2>&1 | grep -vE "^$" | tail -25 || echo "BEAGLE FAILED"
if [ -s 6.beagle.vcf.gz ]; then
  bcftools index -t 6.beagle.vcf.gz
  echo "beagle output records: $(count 6.beagle.vcf.gz)"
  echo "of them in the SGDP input: $(bcftools isec -n=2 -w1 -c none 6.beagle.vcf.gz 4.grch38.vcf.gz | grep -vc '^#')"
  echo "symbolic records: $(bcftools view -H -i 'ALT~"<"' 6.beagle.vcf.gz | wc -l)"
  echo "DR2 of imputed symbolic records, deciles:"; bcftools query -i 'ALT~"<"' -f '%INFO/DR2\n' 6.beagle.vcf.gz | sort -g | awk '{a[NR]=$1} END {for (q = 1; q <= 9; q++) printf "%s ", a[int(NR*q/10)+1]; print ""}'
  echo "symbolic records carried by at least one SGDP sample: $(bcftools view -H -c 1 -i 'ALT~"<"' 6.beagle.vcf.gz | wc -l)"
  bcftools view -h 6.beagle.vcf.gz | grep -E "^##(INFO|FORMAT)" | cut -c1-110
fi
step "done"; du -sh $W
