#!/usr/bin/env bash
##x mem=18G
##x cpus=8

#### Script for transcript quantification running salmon 
###  requires: salmon in the path

## requires a merged manifest file (fq_merge_manifest.pl)
## sampleID folders with cram files are expected in the current directory

## run with:
# arx sub -m16G -c6 -a1- --cfg ../counts.cfg task_salmon.sh ../merged.manifest
refbase="/dcs04/lieber/lcolladotor/dbDev_LIBD001/ref"
salmidx=${SALMON_IDX:-$refbase/salmon_idx/gencode25.main}
ncpus=${SALMON_CPUS:-6}

function err_exit {
 echo -e "Error: $1"
 exit 1
}

SW='/dcs04/lieber/lcolladotor/dbDev_LIBD001/sw'
export PATH=$SW/gscripts:$SW/bin:$PATH
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

rlog=salmon_t${taskid}.tlog

pwd=$(pwd -P) # current directory, absolute path

host=${HOSTNAME%%.*}
#if [[ -n $SLURM_ARRAY_TASK_ID ]]; then
# echo ">task $SLURM_ARRAY_TASK_ID running on $host"
#fi
echo '['$(date '+%m/%d %H:%M')"] task ${jobid}.${taskid} on $host:${pwd}"

# requires a task file as the 1st arg
fdb="$1"
if [[ -z "$fdb" ]]; then
  err_exit "no task file given!"
fi

for f in $salmidx/pos.bin ; do
 if [[ ! -f $f ]]; then
    err_exit "cannot find $f"
 fi
done

## on JHPCE use $MYSCRATCH
if [[ $host == transfer-* || $host == compute-* ]]; then
 tmpdir=$MYSCRATCH/${jobid}_$taskid
else
 if [[ $host == srv05 ]]; then
   tmpdir=/dev/shm/${USER}-${jobid}_$taskid
 else
   tmpdir=$pwd/tmp/${jobid}_$taskid
 fi
fi

mkdir -p $tmpdir || err_exit "failed to create $tmpdir"

line=$(linix $fdb $taskid | cut -f1,3,5)
#expected format:     read_1.fq[,read_1b.fq.gz,...] read_2.fq[,read_2b.fq.gz,...] sampleID
t=( $line )
fn1=${t[0]} # path(s) to read_1 fastq.gz files
fn2=${t[1]}  # sampleID_1.fastq.gz | RNum_*_R1_*.fastq.gz
oid=${t[2]}  # output directory = base sampleID
if [[ -z $oid ]]; then
 err_exit "could not parse base sampleID!"
fi

fqs1=( ${fn1//,/ } )
fqs2=( ${fn2//,/ } )

if [[ ! -f "${fqs1[0]}" ]]; then
   err_exit "FASTQ file ${fqs1[0]} cannot be found!"
fi

outpath=$oid
mkdir -p $outpath
cd $outpath || err_exit "failed at: cd $outpath"
sid=$oid #for these counts, unify sid (merge flowcells)

if [[ -f "$rlog" ]]; then
  bk=1
  while [[ -f "$rlog.$bk" ]]; do
    ((bk++))
  done
  mv $rlog "$rlog.$bk"
fi

echo "["$(date '+%m/%d %H:%M')"] task ${jobid}.${taskid} starting on $host:${pwd}" | tee -a $rlog

## start Salmon
fsalm="salmon/quant.sf"
if [[ ! -s $fsalm ]]; then
  cmd="salmon quant -p $ncpus -lA -1 ${fqs1[@]} -2 ${fqs2[@]} \
   -i $salmidx --numBootstraps 100 --validateMappings -d -o salmon >& salmon.log"
  echo -e "running salmon:\n$cmd" | tee -a $rlog
  run="${run}s"
  eval "$cmd" |& tee -a $rlog &
fi

if [[ $run ]]; then
 wait
fi

if [[ ! -s $fsalm ]]; then
   err_exit "$fsalm have zero size! Check $rlog" |& tee -a $rlog
fi

echo -e "["$(date '+%m/%d %H:%M')"]\tsalmon task done [$taskid]." | tee -a $rlog

