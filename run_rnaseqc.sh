#!/bin/bash
ref=$HOME/work/cbrain/ref/hg38mod_noPARs.fa
qcgtf=$HOME/work/cbrain/ref/gencode.v39_collapsed_genes.gtf
if [[ -z "$1" || $1 = '-h' || $1 = '--help' ]]; then
 echo "Usage: "
 echo "    run_rnaseqc.sh path/to/something.bam|cram"
 echo "Using ref: $qcgtf"
 exit 1
fi

## format MUST be:   ...... /R\d+/R*.cram
fp="$1"
outdir="${fp%/*}" # remove file name (or last dir)
fbam="${fp##*/}" # remove the directories 
fo="${fbam%.*}" # remove extension (bam/cram)
fo="${fo%_sorted*}" # remove "_sorted" suffix

#output will be in the same directory with input BAM, rqc_* subdirs
echo "Output will be in $outdir/rqc/$fo.*" 

# run rnaseqc
cd $outdir
# use samtools view with -T for CRAM files
#(samtools view --input-fmt-option filter='rname=~"^chr[0-9MXY]+$"' -u -T $ref $fp | \
rnaseqc $qcgtf $fbam rqc --mapping-quality=60 --coverage -s $fo >& rqc_$fo.log
echo "done." >> rqc_$fo.log
# when using samtools view| files have bad names ( '-' file name)
#cd rnaseq_out
#rename 's/^\-/r/' *.*

