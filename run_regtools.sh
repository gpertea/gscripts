#!/bin/bash
## script to be run with parallel like this:
# find . -name '*.bed' | parallel --jobs 6 run_regtools.sh '{}'
# or: ls `pwd -P`/*.bed | head -2 | parallel --jobs 1 run_regtools.sh '{}'

function err_exit {
 echo -e "Error: $1"
 exit 1
}

ref_fa=/faststore/ref/hg38c.fa
ref_ann=/faststore/ref/gencode25.exonCDS.gtf

#must be absolute path:
outdir=$(pwd -P)/jann

fbed="$1" #full path to the input bed file, need to parse R\d+_flowcell from it
## /data/gpertea/work/R/dbwrk/BSP1_junctions/R2809_D2A01ACXX_junctions_primaryOnly_regtools.bed
## (or it might lack _flowcell)
fname="${fbed##*/}"
fext="${fname##*.}"
if [[ $fext != "bed" ]]; then
 err_exit "File $fname does not have .bed extension!"
fi
## -- if format is RNum_flowcell then this is needed:
#re="(R[0-9]+_[0-9A-Z]+)*"
## but for BSP1 reprocessing, it's simply RNum:
re="(R[0-9]+)*"
rid=""
if [[ $fname =~ $re ]]; then 
  rid="${BASH_REMATCH[1]}"
else 
 err_exit "could not figure out a RNum ID from $fname"
fi
if [[ -z "$rid" ]]; then
 err_exit "could not parse a RNum from $fname"
fi
#output file name:
fjann=$rid.jann.tab
rlog=$outdir/$rid.run.log
##if bed has extra non-chr contigs/scaffolds:
#cmd="grep '^chr' $fbed  | regtools junctions annotate -o $outdir/$fjann - $ref_fa $ref_ann"
cmd="regtools junctions annotate -o $outdir/$fjann $fbed $ref_fa $ref_ann"
echo -e "Running command:\n$cmd" > $rlog
eval "$cmd" >> $rlog 2>&1
if [[ $? -ne 0 ]]; then
  echo "Error code detected. ">> $rlog
fi
