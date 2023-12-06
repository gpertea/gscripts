#!/usr/bin/env bash
#$ -cwd
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -e logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_fsize=200G
#$ -l mem_free=2G,h_vmem=4G
function err_exit {
 echo -e "Error: $1"
 exit 1
}
#host=$(hostname -s)
f="$1"
if [[ -z "$f" ]]; then
  err_exit "no cmdfile given!"
fi
if [[ -z "$2" ]]; then
  taskid=$SLURM_ARRAY_TASK_ID
else
  taskid=$2
fi
if [[ -z $taskid ]]; then
 err_exit "no SLURM_ARRAY_TASK_ID or task parameter found!"
fi
cmd=$(sed "${taskid}q;d" $f)
if [[ -n $cmd ]]; then
 eval "$cmd"
fi
