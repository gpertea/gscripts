#!/usr/bin/env bash

if [[ -z $1 || "$1" == "-h" || "$1" == "--help" ]]; then
   echo -e "Usage:\n gt_annotate.sh input_gt.bcf [annDataFile]"
   exit 1
fi

function err_exit {
 echo -e "Error: $1"
 exit 1
}
host=$(hostname -s)

fann=$HOME/genotyping/Annotation/GCF_000001405.38_chrUpdate.bcf.gz
if [[ -n $2 ]]; then
  fann=$2
fi

if [[ ! -f $fann ]]; then
    err_exit "cannot find $fann"
fi
f=$1

if [[ ! -f $f ]]; then
    err_exit "cannot find input file $f"
fi
ext="${f##*.}"
fb="${f%.*}" #remove extension
fout=$fb.rsann.bcf
ot="-Ob"
if [[ ext == "gz" ]]; then
  ##output vcf.gz as well
  fb="${f%.*}"
  fout=$fb.rsann.vcf.gz
  ot="-Oz"
fi

bcftools annotate -a $fann --threads 4 -c INFO/RS $ot -o $fout $f



