#!/usr/bin/env bash
##x mem=24G
##x cpus=4

#### Script for running stringtie3 
## requires a list of sample IDs (RNums) which are expected to be sub-directories
## i.e. sampleID folders with cram files expected in $ALN_DIR variable

## run with:
# arx sub -m24G -c4 -a1- -t 02:00:00 -j30 --cfg ../strg_g41.cfg task_strg_asm.sh ../merged.manifest
## or for local SMP:
# arx sub -P -a1- -j24 --cfg strg_gencode47_refseq.cfg ./task_strg3nasc.sh RNums.lst

## these will be overriden by --cfg file
refdir=$HOME/work/cbrain/ref
alndir=${ALN_DIR:-hisat2}
gref=${GENOME_FA:-$refdir/hg38mod_noPARs.fa}
gann=${GENOME_ANN:-$refdir/gencode49/gencode49.nri.gff}
PWD=$(pwd -P)
outbase=$PWD/strg3 ## it will create strg3/sampleID dirs

alndir=$(realpath -- $alndir)

function err_exit {
 echo -e "Error: $1"
 exit 1
}

if [[ ! -d $alndir ]]; then
  err_exit "$alndir not found!"
fi

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

rlog=strg3_${taskid}.tlog

host=${HOSTNAME%%.*}
echo '['$(date '+%m/%d %H:%M')"] task ${jobid}.${taskid} on $host:${PWD}"

# requires a task file as the 1st arg
fdb="$1"
if [[ -z "$fdb" ]]; then
  err_exit "no task file given!"
fi
gref=$(realpath -- $gref)
gann=$(realpath -- $gann)
for f in $gref $gann ; do
 if [[ ! -f $f ]]; then
    err_exit "cannot find $f"
 fi
done

### special case: expects the RNum (subdir name) as the first field
sid=$(linix $fdb $taskid | awk '{print $1}')
#rnum=$sid

if [[ -z $sid ]]; then
 err_exit "could not parse base sampleID!"
fi

#mcram=$alndir/$sid/${sid}*.cram
## grab the first cram file matching that file mask
shopt -s nullglob
matches=("$alndir/$sid/${id}"*.cram)
if (( ${#matches[@]} == 0 )); then
    err_exit "No files match: $alndir/$sid/${id}*.cram"
fi
mcram="${matches[0]}"
if [[ ! -f $mcram ]]; then
  err_exit "$mcram does not exist!"
fi

if [[ -f "$rlog" ]]; then
  bk=1
  while [[ -f "$rlog.$bk" ]]; do
    ((bk++))
  done
  mv $rlog "$rlog.$bk"
fi

odir=$outbase/$sid
if [[ ! -d $odir ]]; then
 mkdir -p $odir
fi

cd $odir || err_exit "failed at cd $odir"
echo "["$(date '+%m/%d %H:%M')"] task ${jobid}.${taskid} starting on $host:${PWD}" | tee -a $rlog
run=''

agtf=$sid.gtf
egtf=$sid.eB.gtf
maingtf=$sid.mainchr.gtf
## run only the assembly
if [[ ! -s $agtf ]]; then
  cmd="stringtie3 -N --cram-ref $gref -o $agtf -G $gann $mcram"
  echo -e "running stringtie3-asm:\n$cmd" | tee -a $rlog
  run="${run}a"
  eval "$cmd" |& tee -a $rlog & 
fi
#if [[ ! -s $egtf ]]; then
#  cmd="stringtie --cram-ref $gref -eB -o $egtf -A $sid.gabund.tab -G $gann ../$mcram"
#  echo -e "running stringtie-e:\n$cmd" | tee -a $rlog
#  run="${run}e"
#  eval "$cmd" |& tee -a $rlog &
#fi

## wait for all background tasks to finish (rnaseqc, regtools, featureCounts exon+gene
if [[ $run ]]; then
 wait
fi
## extract the assemblies on main chromosomes only
grep -P '^chr[0-9XYM]+\t' $agtf > $maingtf

cmd="gffcompare -r $gann -j gcmp_novel_jx.tab -o gcmp $maingtf" 
eval "$cmd" |& tee -a $rlog

echo 1 > "$sid.strg3N.done"
echo -e "["$(date '+%m/%d %H:%M')"]\tstrg task done [$taskid]." | tee -a $rlog

