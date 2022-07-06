#!/bin/bash
## can run with parallel:
# find . -name '*.cram' | \
#  parallel --jobs 6 rnum_cram_rnaseqc.sh '{}'
ref=/ccb/salz8-2/chess-brain/ref/hg38mod_noPARs.fa
qcgtf=$HOME/work/ref/gencode.v39_collapsed_genes.gtf
if [[ -z "$1" || $1 = '-h' || $1 = '--help' ]]; then
 echo "Usage: "
 echo "    run_rnaseqc.sh path/to/R\d+/R\d+_xxxxxx.bam|cram"
 echo "Using ref: $qcgtf"
 exit 1
fi

## format MUST be:   ...... /R\d+/R*.cram
fp="$1"
rnum="${fp%/*}" # remove file name
rnum="${rnum##*/}" # keep only the last directory
#echo "RNum is: <$rnum>"

if [[ ! -d "$rnum" ]]; then
  mkdir -p "$rnum"
fi

# run rnaseqc
cd $rnum
(samtools view --input-fmt-option filter='rname=~"^chr[0-9MXY]+$"' -u -T $ref $fp | \
 rnaseqc $qcgtf - rqc --mapping-quality=60 --coverage -s $rnum) > rnaseqc.log 2>&1
 
#cd rnaseq_out
#rename 's/^\-/r/' *.*

