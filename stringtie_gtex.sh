#!/bin/bash
guides=genome/ref_ann.gff
ref=genome/ref_genome.fa
outdir=stringtie
if [[ -z "$1" || $1 = '-h' || $1 = '--help' ]]; then
 echo "Usage: "
 echo "    stringtie_gtex.sh path/to/tissue/SRAxxxxxx.bam|cram"
 echo "Must be run in the base directory, using $outdir/<tissue> as output directory!"
 echo "Using guides: $guides, reference genome $ref"
 exit 1
fi
if [[ ! -d $outdir ]]; then
  echo "Error: folder $outdir not found!"
  exit 1
fi

if [[ ! -f $guides ]]; then
  echo "Error: reference annotation file $guides not found!"
  exit 1
fi


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

# run stringtie now
if [[ $ext == "bam" ]]; then
(stringtie -p2 -o $outdir/$tissue/$sra.gtf -G $guides $fp) > $outdir/$tissue/$sra.stringie_err.log 2>&1
else #cram assumed
 (samtools view -uh -T $ref $fp | stringtie -p2 -o $outdir/$tissue/$sra.gtf -G $guides --bam -) > $outdir/$tissue/$sra.stringie_err.log 2>&1
fi

