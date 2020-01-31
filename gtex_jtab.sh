#!/bin/bash
#guides=genome/ref_ann.gff
ref=genome/ref_genome.fa
outdir=stringtie
if [[ -z "$1" || $1 = '-h' || $1 = '--help' ]]; then
 echo "Usage: "
 echo "    gtex_jtab.sh path/to/tissue/SRAxxxxxx.bam|cram"
 echo "Must be run in the base directory, using $outdir/<tissue> as output directory!"
 exit 1
fi
if [[ ! -d $outdir ]]; then
  echo "Error: folder $outdir not found!"
  exit 1
fi

#if [[ ! -f $guides ]]; then
#  echo "Error: reference annotation file $guides not found!"
#  exit 1
#fi

fp="$1"
fn="${fp##*/}"
if [ ! -s $fp ]; then
   echo "Error: alignment file empty or not found: $fp"
   ls -l $fp
   exit 1
fi
tissue="${fp%/*}"
tissue="${tissue##*/}"
ext="${fn##*.}"
sra="${fn%.*}"
#echo "Tissue is <$tissue>, SRA accession is <$sra>"
if [[ ! -d "$outdir/$tissue" ]]; then
  mkdir -p "$outdir/$tissue"
fi

# run stringtie-j now
if [[ $ext == "bam" ]]; then
 stringtie-j $fp > $outdir/$tissue/$sra.jtab
else #cram assumed
 samtools view -uh -T $ref $fp | stringtie-j --bam - > $outdir/$tissue/$sra.jtab
fi

