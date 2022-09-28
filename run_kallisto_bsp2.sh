#!/bin/bash
## MUST be run in the output base directory where the output sampleID 
##   sub-directories (will) reside
########
### Can run this with parallel like this, with a list of sample_ids
## whose directories must be within cramdir, e.g.
##  ls -d $HOME/work/cbrain/aln/bsp2_hippo/*/*.cram | \
##     grep -oP '/\K[^\/]+$' | cut -f1 -d_ > sampleIDs
#### then:
## cat sampleIDs | parallel --delay .01 -j 10 ./run_kallisto.sh '{}'
##
function err_exit {
 echo -e "Error: $1"
 exit 1
}

if [[ -z "$1" || $1 = '-h' || $1 = '--help' ]]; then
 echo "Usage: "
 echo "    $0 sampleID"
 exit 1
fi

sid=$1 ## sampleID passed as parameter

gfa=$HOME/work/cbrain/ref/hg38mod_noPARs.fa
g41refdir=/ccb/salz8-3/gpertea/cbrain/gencode41/ref
qcgtf=$g41refdir/gencode41.nri.collapsed.gtf # for rnaseqc
saf=$g41refdir/gencode41.nri.flattened.saf # for gene counts
xgtf=$g41refdir/gencode41.nri.exonsOnly.gtf
ann=$g41refdir/gencode41.nri.gtf
ktx=$g41refdir/gencode41.nri.ktx

## REQUIRED input directories
fqdir=$HOME/work/cbrain/fastq/libd_bsp2/fastq_hippo
#cramdir=$HOME/work/cbrain/aln/bsp2_hippo

if [[ ! -f $ktx ]]; then err_exit "kallisto index $ktx not found!"; fi

fqfiles=$(ls $fqdir/$sid*.f*q.gz)

if [[ -z "$fqfiles" ]]; then
 err_exit "no fastq files found as $fqdir/$sid*.f*q.gz"
fi

if [[ ! -d $sid/kallisto ]]; then
 mkdir -p $sid/kallisto
fi
cd $sid/kallisto || err_exit "could not cd $sid/kallisto"

af="abundance.tsv"
if [[ -f $af ]]; then
  echo "$af already exists, skipping."
  exit 0
fi

ksflag="--rf-stranded" ## these are all ribo-zero samples

kallisto quant -t 4 -i $ktx $ksflag -o . $fqfiles >& k.log 

echo "Done." >> k.log


