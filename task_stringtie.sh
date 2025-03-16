#!/usr/bin/env bash
##x mem=24G
##x cpus=4

#### Script for running stringtie to assemble RNA-seq samples

## requires a merged manifest file (use fq_merge_manifest.pl on trimmed_samples.manifest)
## sampleID folders with cram files are expected in the current directory

## run with this command on JHPCE:
# arx sub -m24G -c4 --qos=shared-200-2 -t 6:00:00 -a1- -j80 --cfg ../strg_g41.cfg task_stringtie.sh ../merged.manifest
#

gref=${GENOME_FA:-$refdir/fa/assembly_hg38_gencode_v25_main.fa}
gann=${GENOME_ANN:-$refdir/gtf/gencode25.main.gtf}

function err_exit {
 echo -e "Error: $1"
 exit 1
}

jobid=$SLURM_JOBID
taskid="$2" # the task# could be given directly (e.g. by parallel)
if [[ -z $jobid ]]; then
 jobid="$$"
fi
if [[ -z $taskid ]]; then
   taskid=$SLURM_ARRAY_TASK_ID
   if [[ -z "$taskid" ]]; then
     err_exit "no task index number given or found in $SLURM_ARRAY_TASK_ID!"
   fi
fi

rlog=strg_a_t${taskid}.tlog

host=${HOSTNAME%%.*}
echo '['$(date '+%m/%d %H:%M')"] task ${jobid}.${taskid} on $host:${PWD}"

# requires a task file as the 1st arg
fdb="$1"
if [[ -z "$fdb" ]]; then
  err_exit "no task file given!"
fi

for f in $gref $gann ; do
 if [[ ! -f $f ]]; then
    err_exit "cannot find $f"
 fi
done

##path having libd_bsp1, cmc1 etc. dirs 
# only needed when relative paths are given in taskdb
pwd=$(pwd -P) # current directory, absolute path

line=$(linix $fdb $taskid | cut -f1,3,5)
#expected format:     read_1.fq[,read_1b.fq.gz,...] read_2.fq[,read_2b.fq.gz,...] sampleID
t=( $line )
fn1=${t[0]}
fn2=${t[1]}
sid=${t[2]}
if [[ -z $sid ]]; then
 err_exit "could not parse base sampleID!"
fi

mcram=$sid.cram # merged CRAM across flowcells/lanes

if [[ ! -f $sid/$mcram ]]; then
  err_exit "$sid/$mcram does not exist!"
fi

outpath=$sid
cd $outpath || err_exit "failed at: cd $outpath"
#crams=( $(ls ${sid}*.cram) ) # could be one per flowcell
## need actual files, not symlinks

if [[ -f "$rlog" ]]; then
  bk=1
  while [[ -f "$rlog.$bk" ]]; do
    ((bk++))
  done
  mv $rlog "$rlog.$bk"
fi

echo "["$(date '+%m/%d %H:%M')"] task ${jobid}.${taskid} starting on $host:${PWD}" | tee -a $rlog


run=''

odir=strg
if [[ ! -d $odir ]]; then
 mkdir -p $odir
fi

cd $odir || err_exit "failed at cd $odir"
agtf=$sid.gtf
egtf=$sid.eB.gtf

if [[ ! -s $agtf ]]; then
  cmd="stringtie --cram-ref $gref -o $agtf -G $gann ../$mcram"
  echo -e "running stringtie-asm:\n$cmd" | tee -a $rlog
  run="${run}a"
  eval "$cmd" |& tee -a $rlog & 
fi
if [[ ! -s $egtf ]]; then
  cmd="stringtie --cram-ref $gref -eB -o $egtf -A $sid.gabund.tab -G $gann ../$mcram"
  echo -e "running stringtie-e:\n$cmd" | tee -a $rlog
  run="${run}e"
  eval "$cmd" |& tee -a $rlog &
fi

## wait for all background tasks to finish (rnaseqc, regtools, featureCounts exon+gene
if [[ $run ]]; then
 wait
fi

cmd="gffcompare -r $gann -j gcmp_novel_jx.tab -o gcmp $agtf" 
eval "$cmd" |& tee -a $rlog

echo 1 > "$sid.strg.done"
echo -e "["$(date '+%m/%d %H:%M')"]\tstrg task done [$taskid]." | tee -a $rlog

