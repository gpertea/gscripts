#!/usr/bin/env bash
##x mem=32G
##x cpus=6
#### change refdir below accordingly

## run with arx:
#  arx sub -a 1- -J bsp45 -m32G -c6 -M geo.pertea@libd.org task_fq2cram.sh samples.manifest

## add -P option to run with parallel

## or with sbatch: 
##   sbatch -a 1-210 --mem=32G -c 6 task_fq2cram.sh samples.manifest


kvalue=${HISAT_K:-40} # -k option of HISAT2
## TODO: these should be pulled from a config file passed to this script
refdir=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/SPEAQeasy/Annotation/reference/hg38
gref_base=gencode_v25_main
gref=${GENOME_FA:-$refdir/assembly/fa/assembly_hg38_${gref_base}.fa}
hsidx=${HISAT_IDX:-$refdir/assembly/index/hisat2_assembly_hg38_${gref_base}}
hscpus=${HISAT_CPUS:-6}
#gref_tx="$refdir/transcripts/kallisto/kallisto_index_hg38_gencode_v25"
#salm_tidx="$refdir/transcripts/salmon/salmon_index_hg38_gencode_v32"

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
rlog=fq2cram_t${taskid}.tlog
#JMDIR="$HOME/_jobs/$jobid"
#mkdir -p $JMDIR
#
#rlog=$JMDIR/$flog
## ---- start main task execution
#/bin/rm -f $flog
#ln -s $rlog

host=${HOSTNAME%%.*}
#if [[ -n $SLURM_ARRAY_TASK_ID ]]; then
# echo ">task $SLURM_ARRAY_TASK_ID running on $host"
#fi
echo '['$(date '+%m/%d %H:%M')"] task ${jobid}.${taskid} on $host:${PWD}"

# requires a task file as the 1st arg
fdb="$1"
if [[ -z "$fdb" ]]; then
  err_exit "no task file given!"
fi

for f in $fdb $gref $hsidx.1.ht2 ; do
 if [[ ! -f $f ]]; then
    err_exit "cannot find $f"
 fi
done

##path having libd_bsp1, cmc1 etc. dirs 
# only needed when relative paths are given in taskdb
pwd=$(pwd -P) # current directory, absolute path

## on JHPCE use $MYSCRATCH
if [[ $host == transfer-* || $host == compute-* ]]; then
 tmpdir=$MYSCRATCH/fq2cram/${jobid}_$taskid
else
 if [[ $host == srv05 ]]; then
   tmpdir=/dev/shm/${USER}-${jobid}_$taskid
 else
   tmpdir=$pwd/tmp/${jobid}_$taskid
 fi
fi

mkdir -p $tmpdir || err_exit "failed to create $tmpdir"

line=$(linix $fdb $taskid | cut -f1,3,5)
#expected format:     read_1.fq read_2.fq sampleID
t=( $line )
fn1=${t[0]} # path to first fastq.gz file
fn2=${t[1]}  # sampleID_1.fastq.gz | RNum_*_R1_*.fastq.gz
oid=${t[2]}  # output directory = base sampleID
if [[ -z $oid ]]; then
 err_exit "could not parse base sampleID!"
fi
#oid=$sid ## but we keep flowcells in separate output files for the same RNum 
if [[ ! -f "$fn1" || ! -f "$fn2" ]]; then
   err_exit "Files cannot be found!"
fi

sid=$oid # for now
fn=${fn1##*/} # remove path
## if it's RNum_flowcellXX_*_R1_*.fastq.gz
if [[ $fn == *_*_*.f*q.gz ]]; then
  a=${fn%%_*}       #remove any _ tokens after first
  rest=${fn#*_}     #remove first _ token (rnum)
  b=${rest%%_*} 
  sid=${a}_$b
fi

ofn=$sid # output file base name (may include flowcell)
## NOTE: hisat2 bam, cram, log outputs will NOT be merged at this stage
## if there are multiple flowcells/lanes, they will be aligned independently
## sid and ofn may include a flowcell suffix !

outpath=$oid
mkdir -p "$outpath" || err_exit "failed at mkdir -p $outpath"

cd $outpath || err_exit "failed at: cd $outpath"

cram=$ofn.cram
bam=$ofn.bam
if [[ -f "$rlog" ]]; then
  bk=1
  while [[ -f "$rlog.$bk" ]]; do
    ((bk++))
  done
  mv $rlog "$rlog.$bk"
fi
echo "["$(date '+%m/%d %H:%M')"] task ${jobid}.${taskid} starting on $host:${PWD}" | tee -a $rlog

if [[ -f "$cram.done" && -s "$cram.crai" ]]; then
 echo "  $ofn cram already done, skipping this task ($taskid)"
 exit 0
else
 /bin/rm -f ${ofn}*cram*
fi

echo "processing: $sid $ofn" | tee -a $rlog

if [[ ! -f $bam.done ]]; then
  params="--mm -x $hsidx -1 $fn1 -2 $fn2 -k $kvalue 2>${ofn}.align_summary.txt"
  ## - full search, but if you want to use strandness:
  #if [[ $dataset != bsp1* ]]; then
  #  params="--rna-strandness RF $params"
  #  #ribo-zero samples are always reverse-forward
  #fi
  #echo "$line" | tee -a $rlog
  echo "["$(date '+%m/%d %H:%M')"] starting hisat2 into $bam:" | tee -a $rlog
  cmd="hisat2 -p $hscpus --phred33 --min-intronlen 20 $params |\
   samtools view -b -o $bam -"
  echo -e $cmd | tee -a $rlog
  eval "$cmd" |& tee -a $rlog
  if [[ $? -ne 0  || $(stat -c %s $bam 2>/dev/null || echo 0) -lt 100000 ]]; then
    echo "error exit detected (or BAM file too small) aborting" | tee -a $rlog
    exit 1
  fi
  echo 1 > $bam.done
fi

## ---- sort and convert to CRAM
echo '['$(date '+%m/%d %H:%M')"] start sorting+conversion to CRAM" | tee -a $rlog
tmpsrt=$tmpdir/$ofn.bam_srt_tmp
/bin/rm -f ${tmpsrt}*
cmd="samtools sort -O cram,version=3.1 --reference=$gref -T $tmpsrt -o $cram --write-index -m 7G --no-PG -@ 4 $bam"
echo -e "$cmd" | tee -a $rlog

eval "$cmd" |& tee -a $rlog

if [[ $? -ne 0 || $(stat -c %s $cram 2>/dev/null || echo 0) -lt 100000 ]]; then
   echo '['$(date '+%m/%d %H:%M')"] error exit or cram file too small" | tee -a $rlog  
   /bin/rm -rf $tmpdir
   exit 1
fi
/bin/rm -rf $tmpdir
echo 1 > $cram.done

echo '['$(date '+%m/%d %H:%M')"] task #${taskid} done." | tee -a $rlog
