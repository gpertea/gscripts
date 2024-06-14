#!/usr/bin/env bash
##x mem=36G
##x cpus=8
##x time=96:00:00 (or time=4-00:00)
#### change refdir below accordingly

## run with arx on SLURM:
#    arx sub -a 1- -j 20 -t 4-00:00 -m36G -c8 --cfg fq2cram.cfg task_fq2cram.sh samples.manifest
## add -P option to run with GNU parallel

kvalue=${HISAT_K:-40} # -k option of HISAT2
refdir=${GREF_DIR:-/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/SPEAQeasy/Annotation/reference/hg38}
gref_base=${GREF_BASE:-gencode_v25_main}
gref=${GENOME_FA:-$refdir/assembly/fa/assembly_hg38_${gref_base}.fa}
hsidx=${HISAT_IDX:-$refdir/assembly/index/hisat2_assembly_hg38_${gref_base}}
hscpus=${HISAT_CPUS:-6}

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
 elif [[ $host == srv16 ]]; then
   tmpdir=/scratch/tmp/${USER}-${jobid}_$taskid
 else
   tmpdir=$pwd/tmp/${jobid}_$taskid
 fi
fi

mkdir -p $tmpdir || err_exit "failed to create $tmpdir"

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

outpath=$sid
mkdir -p "$outpath" || err_exit "failed at mkdir -p $outpath"

cd $outpath || err_exit "failed at: cd $outpath"

cram=$ofn.cram
bam=$ofn.bam ## only primary alignments, unsorted (for featureCounts)
ucram=$ofn.unsorted.cram
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
  params="--mm -x $hsidx -1 $fq1 -2 $fq2 -k $kvalue 2>${ofn}.align_summary.txt"
  ## - full search, but if you want to use strandness:
  #if [[ $dataset != bsp1* ]]; then
  #  params="--rna-strandness RF $params"
  #  #ribo-zero samples are always reverse-forward
  #fi
  #echo "$line" | tee -a $rlog
  echo "["$(date '+%m/%d %H:%M')"] starting hisat2 into $bam, $ucram:" | tee -a $rlog
  
  cmd="hisat2 -p $hscpus --phred33 --min-intronlen 20 $params |\
   samtools view -u - | tee >(samtools view -F 260 -b -o $bam -) |\
   samtools view -O cram,version=3.1 --threads 2 --reference=$gref -o $ucram -"
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
tmpsrt=$tmpdir/$ofn.cram_srt_tmp
/bin/rm -f ${tmpsrt}*
cmd="samtools sort -O cram,version=3.1 --reference=$gref -T $tmpsrt -o $cram --write-index -m 7G -@ 4 $ucram"
echo -e "$cmd" | tee -a $rlog

eval "$cmd" |& tee -a $rlog

if [[ $? -ne 0 || $(stat -c %s $cram 2>/dev/null || echo 0) -lt 100000 ]]; then
   echo '['$(date '+%m/%d %H:%M')"] error exit or cram file too small" | tee -a $rlog  
   /bin/rm -rf $tmpdir
   exit 1
fi
/bin/rm -rf $tmpdir
# should remove the unsorted. cram file as well
# let's just rename it for now to get it out of the way
mv $ucram $ofn.cram.unsorted
echo 1 > $cram.done

echo '['$(date '+%m/%d %H:%M')"] task #${taskid} done." | tee -a $rlog
