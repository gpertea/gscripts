#!/usr/bin/env bash
function err_exit {
 echo -e "Error: $1"
 exit 1
}
fdb="$1"
if [[ -z "$fdb" || ! -f $fdb ]]; then
  err_exit "no line file given!"
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
rlog=degstrg_t${taskid}.tlog

host=${HOSTNAME%%.*}
echo '['$(date '+%m/%d %H:%M')"] task ${jobid}.${taskid} on $host:${PWD}"

gann=/dcs04/lieber/lcolladotor/chessBrain_LIBD4085/gencode41/ref/gencode41.nri.gff

fbam=$(linix $fdb $taskid)
if [[ -z "$fbam" ]]; then
  err_exit "cannot get bam file for task $taskid from $fdb"
fi

for f in $gann $fbam ; do
 if [[ ! -f $f ]]; then
    err_exit "cannot find $f"
 fi
done

sid=$(echo $fbam | perl -lane '($m)=(m/(Br\d+_[ABCDEF]\d_[^_]+)/);print $m')
if [[ ! $sid ]]; then
  err_exit "cannot parse sample id from $fbam"
fi
mkdir -p $sid || err_exit "failed at mkdir -p $sid"
cd $sid || err_exit "failed at: cd $sid"

fstrg=$sid.strg.gtf
## run stringtie
rlog=run.log
(samtools view --input-fmt-option filter='rname=~"^chr[0-9MXY]+$"' -u $fbam | stringtie -G $gann -o $fstrg --bam - ) |& tee $rlog
## run gffcompare
gffcompare -r $gann -o gfcmp $fstrg |& tee -a $rlog
echo '['$(date '+%m/%d %H:%M')"] task #${taskid} done." | tee -a $rlog

