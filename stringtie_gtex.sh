#!/bin/bash
if [[ -z "$1" || $1 = '-h' || $1 = '--help' ]]; then
 echo "Usage: "
 echo "    stringtie_gtex.sh /full/path/to/tissue/SRAxxxxxx.bam"
 echo "Warning: must be run in the base output directory!"
 exit 1
fi
fp="$1"
fn="${fp##*/}"
if [ ! -s $fp ]; then
   echo "Error: file empty or not found: $fp"
   ls -l $fp
   exit 1
fi
tissue="${fp%/*}"
tissue="${tissue##*/}"
sra="${fn%.*}"
echo "Tissue is <$tissue> , SRA accession is <$sra>"
if [[ ! -d "$tissue" ]]; then
  mkdir -p "$tissue"
fi

# run stringtie now
cd "$tissue"
stringtie -p2 -o $sra.gtf $fp > $sra.stringie_err.log 2>&1

