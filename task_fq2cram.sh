#!/usr/bin/env bash
##x mem=32G
##x cpus=6
#### change refdir below accordingly
## run with: 
##   sbatch -a 1-210 --mem=32G -c 6 task_fq2cram.sh samples.manifest

## OR run with parallel:
# parallel --delay .01 -j 4 task_fq2cram.sh samples.manifest {1} ::: {1..210}

## OR with arx:
#  arx sub -a 1- -J bsp45 -m32G -c6 -M geo.pertea@libd.org task_fq2cram.sh samples.manifest

kvalue=40 # -k option of HISAT2
## TODO: these should be pulled from a config file passed to this script
refdir='/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/SPEAQeasy/Annotation/reference/hg38'
gref_base='gencode_v25_main'
gref="$refdir/assembly/fa/assembly_hg38_${gref_base}.fa"
hsref="$refdir/assembly/index/hisat2_assembly_hg38_${gref_base}"
gref_tx="$refdir/transcripts/kallisto/kallisto_index_hg38_gencode_v25"
#salm_tidx="$refdir/transcripts/salmon/salmon_index_hg38_gencode_v32"

function err_exit {
 echo -e "Error: $1"
 exit 1
}

#host=$(hostname -s)
#jobid=$JOB_ID 
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
echo "task ${jobid}.${taskid} on $host:${PWD}["$(date '+%m/%d %H:%M')"]" 

if [[ ! -f $gref ]]; then err_exit "not found: $gref"; fi
if [[ ! -f $hsref ]]; then err_exit "not found: $hsref"; fi

##path having libd_bsp1, cmc1 etc. dirs 
# only needed when relative paths are given in taskdb
pwd=$(pwd -P) # current directory, absolute path

## final .cram files will be written under outdir,
## one subdir per sampleID
outdir=$pwd/aln

# requires a task db pseudo-fasta file as the 1st arg
fdb="$1"
if [[ -z "$fdb" ]]; then
  err_exit "no cdb task db given!"
fi
if [[ ! -f "$fdb" ]]; then
  err_exit "$fdb file not found!"
fi
fdbix="$fdb.cidx"
if [[ ! -f "$fdbix" ]]; then
    # cdbfasta $fdb -- no, this would mess up grid/parallel runs
    err_exit "$fdbix file must exist!"
fi
fdb=$fdbix

## on JHPCE use $MYSCRATCH
if [[ $host == transfer-* || $host == compute-* ]]; then
 tmpdir=$MYSCRATCH/fq2cram/$jobid/$taskid
else
 tmpdir=$pwd/tmp/$jobid/$taskid
fi

mkdir -p $tmpdir || err_exit "failed to create $tmpdir"

line=$(linix $fdb $taskid | cut -f1,3,5)
#expected format:     read_1.fq read_2.fq sampleID
#           >6 dataset_dir sampleID_1.fastq.gz M Schizo AA 52.02
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

echo "task ${jobid}.${taskid} starting on $host:${PWD}["$(date '+%m/%d %H:%M')"]" | tee -a $rlog

if [[ -f "$cram.done" ]]; then
 echo "  $ofn cram already done, skipping this task ($taskid)"
 exit 0
fi

#fqarr=($(ls $fdir/${sid}[_.]*f*q.gz))
#nfq=${#fqarr[@]}
#if ((nfq==0)); then
#  err_exit "could not find $fdir/${sid}[_.]*f*q.gz"
#fi
#if ((nfq<2)); then
#  err_exit "not matching paired fastq: ls $fdir/${sid}[_.]*f*q.gz"
#fi
#fqs1=${fqarr[0]}
#fqs2=${fqarr[1]}
## there could be multiple lanes per sample, merge them 
#n=2
#if ((nfq>2)); then
#  until ((n>=nfq)); do
#    fqs1=$fqs1,${fqarr[$n]}
#    ((n=n+1))
#    fqs2=$fqs2,${fqarr[$n]}
#    ((n=n+1))
#  done
#fi
echo "processing: $sid $ofn" | tee -a $rlog

if [[ ! -f $bam.done ]]; then
  params="--mm -x $hsref -1 $fn1 -2 $fn2 -k $kvalue 2>${ofn}.align_summary.txt"
  ## - full search, but if you want to use strandness:
  #if [[ $dataset != bsp1* ]]; then
  #  params="--rna-strandness RF $params"
  #  #ribo-zero samples are always reverse-forward
  #fi
  #echo "$line" | tee -a $rlog
  echo "["$(date '+%m/%d %H:%M')"] starting hisat2 into $bam:" | tee -a $rlog
  cmd="hisat2 -p 6 --phred33 --min-intronlen 20 $params |\
   samtools view -b -o $bam -"
  echo -e $cmd | tee -a $rlog
  tmpsrt=$tmpdir/$fn.bam_srt_tmp
  /bin/rm -f ${tmpsrt}*
  eval "$cmd" |& tee -a $rlog
  if [[ $? -ne 0  || $(stat -c %s $bam 2>/dev/null || echo 0) -lt 100000 ]]; then
    echo "error exit detected (or BAM file too small) aborting" | tee -a $rlog
    exit 1
  fi
  echo 1 > $bam.done
fi

## ---- sort and convert to CRAM
echo '['$(date '+%m/%d %H:%M')"] start sorting+conversion to CRAM" | tee -a $rlog
#cmd="samtools sort -O cram,use_lzma=1,use_tok=1,use_fqz=1,seqs_per_slice=50000 \
# --reference=$gref -T $tmpsrt -o $cram -m 7G --no-PG -@ 4 $bam"
cmd="samtools sort -O cram,version=3.1 --reference=$gref -T $tmpsrt -o $cram -m 7G --no-PG -@ 4 $bam"
## "scramble -P -B -I bam -O cram -8 -r $gref -X small -t 4 -! - $cram"
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
