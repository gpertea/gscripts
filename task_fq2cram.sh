#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -e logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_fsize=100G
#$ -l mem_free=16G,h_vmem=18G
#$ -pe local 4

kvalue=80 # this is the -k option of HISAT2

## run with: 
#   qsub -t 1-56 ../../../code/task_hisat2.sh taskdb.tfa
## taskdb expected format (4th field, output SAMPLE_ID, is optional)
# >3 AN00000904/Br8667_Mid_Nuc AN00000904_Br8667_Mid_Nuc_1.fastq.gz F       AN00000904_Br8667_Mid_Nuc
#  0               1 (dir)           2 (mate 1 file mask)           3 (sex) 4 (output SAMPLE_ID base file name)
## for multi-flowcell:
# >3 libd_bsp4and5 R11338_*_R1_*.fastq.gz M R11338_HCCFWBBXX

## OR run with parallel:
# parallel --delay .01 -j 8 task_fq2cram.sh taskdb.tfa {1} ::: {2..210}

## run with arx from base dataset directory with ./fastq
# code=$HOME/work/cbrain/code/01_prep 
#  $code/prep_taskdb.pl libd_bsp4and5 brainseq_phase4and5: > taskdb.cfa
#  cdbfasta taskdb.cfa
#  $code/../arx sub -t 1-527 -N bsp45 -m geo.pertea@libd.org $code/task_fq2cram.sh taskdb.cfa

function err_exit {
 echo -e "Error: $1"
 exit 1
}

#host=$(hostname -s)
jobid=$JOB_ID 
taskid="$2" # the task# could be given directly (e.g. by parallel)
if [[ -z $jobid ]]; then
 jobid="$$"
fi
if [[ -z $taskid ]]; then
   taskid=$SGE_TASK_ID
   if [[ -z "$taskid" ]]; then
     err_exit "no task index number given or found in SGE_TASK_ID!"
   fi
fi
flog=t_${taskid}.tlog
JMDIR="$HOME/_jobs/$jobid"
mkdir -p $JMDIR

rlog=$JMDIR/$flog
## ---- start main task execution
/bin/rm -f $flog
#ln -s $rlog

echo "task ${jobid}.${taskid} starting on $host:${PWD}["$(date '+%m/%d %H:%M')"]" | tee $rlog

host=${HOSTNAME%%.*}
if [[ -n $SGE_TASK_ID ]]; then
 echo ">task $SGE_TASK_ID running on $host"
 module load conda_R/4.1.x
fi

##dataset="deconvo"

#base dir for ref data
basedir=$HOME/work/cbrain
refbase=$basedir/ref
hsrefM=$refbase/hisat2_hg38mod_noPARs/hg38mod_noPARs
hsrefF=$refbase/hisat2_hg38mod_noY/hg38mod_noY
gref=$refbase/hg38mod_noPARs.fa
if [[ ! -f $gref ]]; then err_exit "not found: $gref"; fi
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

line=$(cdbyank -a $taskid $fdb)
#expected format:     1            2           3   4     5   6
#           >6 dataset_dir sampleID_1.fastq.gz M Schizo AA 52.02
if [[ $line != '>'* ]]; then
 err_exit "invalid line pulled for $taskid:\n$line"
fi
t=( $line )
fdir=${t[1]} # path to fastq.gz file
fn1=${t[2]}  # sampleID_1.fastq.gz | RNum_*_R1_*.fastq.gz
sx=${t[3]}
oid=${t[4]}  # output sample id, only used if present and matches the sample id from fn1

rnum='' ## deconvo files, CMC data lack RNum !
sid='' ## we could keep flowcells separate for the same RNum 
if [[ $fn1 == *_1.f*q.gz ]]; then
  fn2=${fn1/_1./_2.}
  if [[ "$fn1" == "$fn2" ]]; then 
   err_exit "Paired files identical? ($fn1 vs $fn2)"
  fi
  sid=${fn1%%.*}  # remove all extensions
  sid=${sid/%_1/} # remove _1 at the end
else ## assume it's RNum_flowcellXX_*_R1_*.fastq.gz
  if [[ $fn1 == R*_*_R1_*.f*q.gz ]]; then
    rnum=${fn1%%_*} #remove any _ tokens after first
    rest=${fn1#*_}  #remove first _ token (rnum)
    #fcell=${rest%%_*}
    ##sid=${rnum}_$fcell ## do this only if RNum_flowcell should be kept separate
    sid=$rnum # merge outputs from multiple flowcells and lanes by $rnum
    fn2=${fn1/_R1_/_R2_}
  fi
fi

ofn=$sid #output file base name 
if [[ -n $oid ]]; then
  if [[ "$oid" == "$sid"* ]]; then
    ofn=$oid;
  fi
fi

if [[ -z $sid ]]; then
 err_exit "could not parse sampleID from $fn1"
fi


## -- get the new aln path pattern:
#dsdir=${fdir//libd_/}
#dsdir=${dsdir//\/fastq/}/$sid
outpath=$outdir/$sid
if [[ -n "$rnum" ]]; then
   outpath=$outdir/$rnum
fi

mkdir -p "$outpath" || err_exit "failed at mkdir -p $outpath"

#proto="ribo"
#if [[ $dataset = *"bsp1"* ]]; then 
# proto="polyA"
#fi

hsref=$hsrefM
if [[ $sx == "F" ]]; then
 hsref=$hsrefF
fi
# indir -> pwd in this case
#fqarr=($(ls  $indir/$sub/${sid}*.f*q.gz))
if [[ $fdir != '/'* ]]; then
  if [[ -d $fdir ]]; then #relative path
     fdir=$pwd/$fdir # convert to absolute path
  else #MUST be relative to $basedir/fastq/
   if [[ -d $basedir/fastq/$fdir ]]; then
     fdir=$basedir/fastq/$fdir # convert to absolute path
   else
     err_exit "$fdir not located in ./ or $basedir/fastq"
   fi
  fi
fi

## -- assume absolute path was given

fqarr=($(ls $fdir/${sid}[_.]*f*q.gz))
nfq=${#fqarr[@]}
if ((nfq==0)); then
  err_exit "could not find $fdir/${sid}[_.]*f*q.gz"
fi
if ((nfq<2)); then
  err_exit "not matching paired fastq: ls $fdir/${sid}[_.]*f*q.gz"
fi
fqs1=${fqarr[0]}
fqs2=${fqarr[1]}
## there could be multiple lanes per sample, merge them 
n=2
if ((nfq>2)); then
  until ((n>=nfq)); do
    fqs1=$fqs1,${fqarr[$n]}
    ((n=n+1))
    fqs2=$fqs2,${fqarr[$n]}
    ((n=n+1))
  done
fi

#fn=$sid
#fn=$ofn

params="--mm -x $hsref -1 $fqs1 -2 $fqs2 -k $kvalue 2>${ofn}.align_summary.txt"

## - full search, but if you want to use strandness:
#if [[ $dataset != bsp1* ]]; then
#  params="--rna-strandness RF $params"
#  #ribo-zero samples are always reverse-forward
#fi

cd $outpath || err_exit "failed at: cd $outpath"
ln -s $rlog

mkdir -p fastq
cd fastq #just to keep symlinks to all source FASTQ files
ln -s $fdir/${sid}*.f*q.gz . > /dev/null 2>&1
cd ..

cram=$ofn.cram
bam=$ofn.bam

if [[ -f $cram ]]; then
 if [[ $(stat -c %s $bam 2>/dev/null || echo 0) -gt 100000 ]]; then
  err_exit " $cram already exists in $PWD!"
 else
  echo " warning: removing existing .cram file"
  unlink $cram
 fi
fi
if [[ -f $bam ]]; then
 echo " warning: removing existing .bam file"
 /bin/rm -f $bam ## WARNING! existing BAM will be lost
fi

#rlog=$sid.log
echo "$line" | tee -a $rlog
echo "processing sample: $sid" | tee -a $rlog
echo "["$(date '+%m/%d %H:%M')"] starting:" | tee -a $rlog
cmd="hisat2 -p 4 --phred33 --min-intronlen 20 $params |\
 samtools view -b -o $bam -"
echo -e $cmd | tee -a $rlog

tmpsrt=$tmpdir/$fn.srt_tmp
/bin/rm -f ${tmpsrt}*

### DEBUG ONLY
#echo "task done." | tee -a $rlog
#exit 1
###

eval "$cmd" |& tee -a $rlog

if [[ $? -ne 0  || $(stat -c %s $bam 2>/dev/null || echo 0) -lt 100000 ]]; then
  echo "error exit detected (or BAM file too small) aborting" | tee -a $rlog
  exit 1
fi
## ---- sort and convert to CRAM
echo '['$(date '+%m/%d %H:%M')"] start sorting/conversion to CRAM" | tee -a $rlog
cmd="samtools sort -T $tmpsrt -l 6 -m 8G --no-PG -@ 4 $bam | \
 scramble -P -B -I bam -O cram -8 -r $gref -X small -t 4 -! - $cram"
echo -e "$cmd" | tee -a $rlog

eval "$cmd" |& tee -a $rlog
if [[ $? -eq 0 ]]; then
  if [[ $(stat -c %s $cram 2>/dev/null || echo 0) -gt 100000 ]]; then
   #/bin/rm -f $bam
   echo "keeping unsorted bam file: $bam"
  else
   echo "error: cram file too small" | tee -a $rlog
  fi
else
   echo "sort|scramble error exit detected!" | tee -a $rlog
   exit 1
fi

echo '['$(date '+%m/%d %H:%M')"] task #${taskid} done." | tee -a $rlog
