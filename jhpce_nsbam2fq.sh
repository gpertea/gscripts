#!/bin/bash

## tasks.tfa was generated for this script like this:
# ls /dcl02/lieber/ajaffe/PublicData/CMC/BAM/BAM_sorted/*_mapped.namesorted.bam | \
#  perl -pe 's{.+/([^/]+)$}{$1};s/(.+)_mapped\.namesorted\.bam$/>$. $1/' > tasks.tfa
# cdbfasta tasks.tfa

## This script should be placed in the same working directory with tasks.tfa
##    (where ./fastq/ is a subdir already created for output).
## Submit like this:
# qsub -cwd -l mem_free=4G,h_vmem=4G -pe local 4 -N bam2fqz -t 1-613 ./jhpce_nsbam2fq.sh tasks.tfa

indir=/dcl02/lieber/ajaffe/PublicData/CMC/BAM/BAM_sorted
outdir=./fastq #set the output directory 
fdb="$1"
if [[ -z "$fdb" ]]; then
  echo "Error: no cdb file given"; exit 1
fi
if [[ ! -f "$fdb" ]]; then
  echo "Error: $fdb file not found!" ; exit 1
  fdb="$fdb.cidx"
  if [[ ! -f "$fdb" ]]; then
    echo "Error: $fdb index file not found!"; exit 1
  fi
fi

id=$SGE_TASK_ID
if [[ -z "$id" ]]; then
 id="$2"
 if [[ -z "$id" ]]; then
  echo "Error: no task ID given!"
  exit 1
 fi 
fi
# get file prefix to use to get mapped/unmapped bam files
r=($(cdbyank -a $id $fdb))
fpre=${r[1]}
flog="$outdir/$fpre.log"
cmd="samtools merge -nu - $indir/${fpre}_mapped.namesorted.bam $indir/${fpre}_unmapped.namesorted.bam |\
 samtools fastq -0 $outdir/${fpre}_r0.fastq.gz -s $outdir/${fpre}_s.fastq.gz \
 -1 $outdir/${fpre}_r1.fastq.gz -2 $outdir/${fpre}_r2.fastq.gz -c 8 -@ 4 -"
echo "Running: $cmd" > $flog
eval "$cmd" >> $flog 2>&1

# singleton and r0 gz files are 20 bytes when empty, remove those
for f in $outdir/${fpre}_r0.fastq.gz $outdir/${fpre}_s.fastq.gz; do
 if [[ $(stat -c %s "$f") -le 20 ]]; then 
   unlink $f
 fi
done

