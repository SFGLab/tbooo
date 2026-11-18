#!/bin/bash
#SBATCH --job-name=kg_panel_2022
#SBATCH --partition=short
#SBATCH --account=sfglab
#SBATCH --cpus-per-task=6
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
#SBATCH --output=/mnt/evafs/groups/sfglab/data/tbooo/data/slurm/sgdp_grch38/download_panel.%j.log
# The 1000 Genomes 30x phased panel of 2022-04-22, with structural variants, under the names
# `tbooo download 1kg` gives it. Resumes a partial file and checks every file against the
# release's own manifest.
set -uo pipefail
T=/mnt/evafs/groups/sfglab/data/tbooo
OUT=$T/data/raw/1kg
B=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
cd $OUT
wget -q -O 20220804_manifest.txt $B/20220804_manifest.txt
names() { for c in $(seq 1 22); do echo 1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz; done; echo 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz; }
{ names; names | sed 's/$/.tbi/'; } > panel_files.txt
# The manifest has two columns, the file name and its md5.
md5_of() { awk -v f="$1" '$1 == f {print $2}' 20220804_manifest.txt | head -1; }
get() { f=$1; want=$(md5_of "$f"); if [ -s "$f" ] && [ -n "$want" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$want" ]; then echo "have $f"; return; fi; wget -q -c -O "$f" "$2/$f" && echo "got $f $(stat -c %s "$f")" || echo "FAILED $f"; }
export -f md5_of
export -f get
xargs -P 6 -I{} bash -c 'get {} '"$B" < panel_files.txt
echo "=== check against the manifest"
bad=0
while read -r f; do
  want=$(md5_of "$f")
  if [ -z "$want" ]; then echo "not in manifest: $f"; bad=$((bad+1)); continue; fi
  [ "$(md5sum "$f" 2>/dev/null | cut -d' ' -f1)" = "$want" ] || { echo "MD5 MISMATCH $f"; bad=$((bad+1)); }
done < panel_files.txt
echo "files: $(wc -l < panel_files.txt), mismatches: $bad"; du -ch $(cat panel_files.txt) | tail -1
