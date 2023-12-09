#!/usr/bin/env bash
##x mem=10G
##x cpus=8

#### Script for adapter trimming with adapter autodetection 
###  requires: fastp in the path

## requires the original samples.manifest file
## output is written in a <sampleID> folder for each pair of fastq files
## as <SeqSampleID>_1.trimmed.fq.gz
##
##   where <SeqSampleID> may be <sampleID>_<flowcell>
## run with:
# arx sub -m8G -c6 -a1- task_fqtrim.sh ../samples.manifest

ncpus=8

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
     echo "Usage example (local execution with parallel):"
     echo "arx sub -P -a1- -j 6 task_fqtrim.sh ../samples.manifest"
     err_exit "no task index number given or found in \$SLURM_ARRAY_TASK_ID!"
   fi
fi

rlog=fqtrim_t${taskid}.tlog

pwd=$(pwd -P) # current directory, absolute path

host=${HOSTNAME%%.*}
echo '['$(date '+%m/%d %H:%M')"] task ${jobid}.${taskid} on $host:${pwd}"

# requires a task file as the 1st arg
fdb="$1"
if [[ -z "$fdb" ]]; then
  err_exit "no task file given!"
fi
fastp=$(which fastp)
if [[ -z "$fastp" || ! -f $fastp ]]; then
   err_exit "cannot find fastp program!"
fi

line=$(linix $fdb $taskid | cut -f1,3,5)
#expected format:     read_1.fq read_2.fq sampleID
t=( $line )
fq1=${t[0]} # path to first fastq.gz file
fq2=${t[1]}  # sampleID_1.fastq.gz | RNum_*_R1_*.fastq.gz
sid=${t[2]}  # output directory = base sampleID
if [[ -z $sid ]]; then
 err_exit "could not parse base sampleID!"
fi

if [[ ! -f "$fq1" || ! -f "$fq2" ]]; then
   err_exit "Files cannot be found!"
fi

fn=${fq1##*/} # remove path, keep only filename
## if it's RNum_flowcellXX_*_R1_*.fastq.gz
## R4225_D1AAPACXX_ATTCCT_L005_R1_001.fastq.gz

if [[ $fn == *.f*q.[gb]z* ]]; then
  fbase=${fn/.f*q.[gb]z*/}
else
  fbase=${fn/.f*q/}
fi
fbase=${fbase/_R1_/_}
fbase=${fbase/_1[._]/_} # just in case it was _1.token.fastq.gz
[[ $fbase == *_1 ]] && fbase="${fbase%_*}"

ofn=$fbase # keep flowcell and lane etc. if they were there

## output directory (sid) is the merged/common prefix as given in the manifest
## but ofn / output file base names may have a flowcell and lane suffix etc.
## NOTE: FASTQ files are not merged across flowcells and lanes at this stage, they are trimmed and QC-ed separately
## if there are multiple flowcells/lanes, they will be aligned independently
## so $ofn may include a flowcell and lane suffixes !

outpath=$sid ## directory might be shared by multiple lines in samples.manifest
mkdir -p $outpath
cd $outpath || err_exit "failed at: cd $outpath"

if [[ -f "$rlog" ]]; then
  bk=1
  while [[ -f "$rlog.$bk" ]]; do
    ((bk++))
  done
  mv $rlog "$rlog.$bk"
fi

echo "["$(date '+%m/%d %H:%M')"] task ${jobid}.${taskid} starting on $host:${pwd}" | tee -a $rlog

## start fastp quality and adapter trimming
ofq1=${ofn}_1.trimmed.fq.gz
ofq2=${ofn}_2.trimmed.fq.gz
if [[ -s $ofq1 ]]; then
  err_exit "output file $outpath/$ofq1 already exists. Skipping."
fi

cmd="fastp -i $fq1 -I $fq2 -l 30 -Q -w $ncpus --dont_eval_duplication  -j ${ofn}_fastp.json -h ${ofn}_fastp.html \
    --overlap_diff_percent_limit 10 --detect_adapter_for_pe -z 7 -o $ofq1 -O $ofq2 >& ${ofn}_fastp.log"
echo -e "running fastp:\n$cmd" | tee -a $rlog
eval "$cmd" |& tee -a $rlog 
if [[ ! -s $ofq1 ]]; then
   err_exit "$outpath/$ofq1 has zero size! Check $rlog" |& tee -a $rlog
fi

echo -e "["$(date '+%m/%d %H:%M')"]\tfastp trimming task done [$taskid]." | tee -a $rlog

