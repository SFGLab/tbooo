#!/bin/bash
#SBATCH --job-name=kg_chr7
#SBATCH --partition=short
#SBATCH --account=sfglab
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=08:00:00
#SBATCH --output=/mnt/evafs/groups/sfglab/data/tbooo/data/slurm/sgdp_grch38/panel_chr7_ranges.%j.log
# The last file of the 30x panel, fetched in eight byte ranges side by side.
#
# The panel job 1945120 had every file but chromosome 7, whose one connection fell to 46 KB/s.
# Here each range has its own connection, a connection under 300 KB/s for 15 seconds is dropped
# and the range goes on from the byte it reached. The parts are joined under a new name and
# checked against the release's manifest. Only when the md5 agrees is the slow job stopped, its
# partial file moved aside as .partial-<job> and the whole file given its name. The check of
# every panel file against the manifest follows, with the line the later jobs read.
set -uo pipefail
T=/mnt/evafs/groups/sfglab/data/tbooo
OUT=$T/data/raw/1kg
B=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV
F=1kGP_high_coverage_Illumina.chr7.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
SLOW=${SLOW_JOB:-1945120}
N=8
W=$OUT/ranges_chr7
mkdir -p $W
cd $W
md5_of() { awk -v f="$1" '$1 == f {print $2}' $OUT/20220804_manifest.txt | head -1; }
size=$(curl -sI --connect-timeout 30 $B/$F | awk 'tolower($1) == "content-length:" {print $2}' | tr -d '\r')
want=$(md5_of $F)
[ -n "$size" ] && [ -n "$want" ] || { echo "no size or no md5 for $F"; exit 1; }
part=$(( (size + N - 1) / N ))
date +"=== %H:%M:%S $F, $size bytes in $N ranges"
fetch() {
  local i=$1 a b need have
  a=$(( i * part )); b=$(( a + part - 1 )); [ $b -ge $size ] && b=$(( size - 1 )); need=$(( b - a + 1 ))
  for _ in $(seq 1 600); do
    have=$(stat -c %s part.$i 2>/dev/null || echo 0)
    [ $have -ge $need ] && return 0
    curl -s --fail --connect-timeout 30 --speed-limit 300000 --speed-time 15 -r $(( a + have ))-$b $B/$F >> part.$i
  done
  return 1
}
for i in $(seq 0 $(( N - 1 ))); do fetch $i & done
while sleep 60; do
  jobs -r | grep -q . || break
  echo "$(date +%H:%M:%S) $(cat part.* 2>/dev/null | wc -c) of $size bytes"
done
wait
total=0
for i in $(seq 0 $(( N - 1 ))); do total=$(( total + $(stat -c %s part.$i 2>/dev/null || echo 0) )); done
[ $total -eq $size ] || { echo "the parts hold $total of $size bytes"; exit 1; }
cat $(for i in $(seq 0 $(( N - 1 ))); do echo part.$i; done) > $OUT/$F.ranges
got=$(md5sum $OUT/$F.ranges | cut -d' ' -f1)
date +"=== %H:%M:%S joined, md5 $got, the manifest has $want"
[ "$got" = "$want" ] || { echo "MD5 MISMATCH of the joined file, nothing is changed"; exit 1; }
# The slow job is still writing the file of this name. It is stopped before the swap.
scancel $SLOW
for _ in $(seq 1 60); do squeue -h -j $SLOW 2>/dev/null | grep -q . || break; sleep 5; done
squeue -h -j $SLOW 2>/dev/null | grep -q . && { echo "job $SLOW is still there, nothing is swapped"; exit 1; }
[ -e $OUT/$F ] && mv -n $OUT/$F $OUT/$F.partial-$SLOW
mv -n $OUT/$F.ranges $OUT/$F
ls -la $OUT/$F $OUT/$F.partial-$SLOW
cd $OUT
echo "=== check against the manifest"
bad=0
while read -r f; do
  w=$(md5_of "$f")
  if [ -z "$w" ]; then echo "not in manifest: $f"; bad=$((bad+1)); continue; fi
  [ "$(md5sum "$f" 2>/dev/null | cut -d' ' -f1)" = "$w" ] || { echo "MD5 MISMATCH $f"; bad=$((bad+1)); }
done < panel_files.txt
echo "files: $(wc -l < panel_files.txt), mismatches: $bad"
du -ch $(cat panel_files.txt) | tail -1
date +"=== %H:%M:%S done"
[ $bad -eq 0 ]
