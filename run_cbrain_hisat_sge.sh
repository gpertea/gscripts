#!//usr/bin/env bash
#$ -cwd
#$ -N cmc1_aln
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -e logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_fsize=100G
#$ -l mem_free=8G,h_vmem=10G
#$ -pe local 4
module load conda_R/4.1.x

## IMPORTANT: change this, output will be in: $basedir/aln/$dataset/$sid/
dataset="cmc1" 
## run with: 
## qsub -t 613 -m ea -M geo.pertea@gmail.com scripts/run_cbrain_hisat_sge.sh aln_taskdb.cfa

function err_exit {
 echo -e "Error: $1"
 exit 1
}
host=$(hostname -s)

#base dir for the directory structure of this project (fastq, aln and strg folders required)
basedir=$HOME/work/cbrain
if [[ ! -d $basedir/aln || ! -d $basedir/fastq || ! -d $basedir/strg ]]; then
 err_exit "fastq/aln/strg folders must exist under basedir $basedir"
fi
## final .cram files will be written in subdirectories here,
## one subdir per sampleID
outdir=$basedir/aln/$dataset

refbase=$basedir/ref
hsrefM=$refbase/hisat2_hg38mod/hg38mod_noPARs
hsrefF=$refbase/hisat2_hg38mod_noY/hg38mod_noY

gref=$refbase/hg38mod_noPARs.fa
if [[ ! -f $gref ]]; then err_exit "not found: $gref"; fi
##path having libd_bsp1, cmc1 etc. dirs 
# only needed when relative paths are given in taskdb
indir=$basedir/fastq/$dataset

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
taskid="$2" # the task# could be given directly (e.g. by parallel)
jobid="$$"
if [[ -z $taskid ]]; then
  #taskid=$SLURM_ARRAY_TASK_ID
  #jobid=$SLURM_JOB_ID
  #if [[ -z "$taskid" ]]; then
   taskid=$SGE_TASK_ID
   jobid=$JOB_ID
   if [[ -z "$taskid" ]]; then
     err_exit "no task index number given or found!"
  # fi
  fi
fi

## could be on /dev/shm/$USER/$SLURM_JOB_ID/ on MARCC
## on JHPCE use $MYSCRATCH 
tmpdir=$MYSCRATCH/cbrain_aln/$dataset/$JOB_ID/$SGE_TASK_ID

mkdir -p $tmpdir || err_exit "failed to create $tmpdir"

line=$(cdbyank -a $taskid $fdb)
#assumes absolute path is given
#expected format:       1                 2   3   4      5    6   7
#           >6 /path/to/libd_bsp1/fastq R2828 M DLPFC Schizo AA 52.02
if [[ $line != '>'* ]]; then
 err_exit "invalid line pulled for $taskid:\n$line"
fi
t=( $line )
fdir=${t[1]}
sid=${t[2]} # sample ID (RNum for LIBD data)

sx=${t[3]}
reg=${t[4]}
#dx=${t[5]}
#race=${t[6]}
#age=${t[7]}

## -- get the new aln path pattern:
#dsdir=${fdir//libd_/}
#dsdir=${dsdir//\/fastq/}/$sid
outpath=$outdir/$sid
mkdir -p "$outpath" || err_exit "failed at mkdir -p $outpath"

proto="ribo"
if [[ $dataset = *"bsp1"* ]]; then 
 proto="polyA"
fi

hsref=$hsrefM
if [[ $sx == "F" ]]; then
 hsref=$hsrefF
fi

#fqarr=($(ls  $indir/$sub/${sid}*.f*q.gz))
## -- assume absolute path was given
fqarr=($(ls  $fdir/${sid}*.f*q.gz))
if [[ ${#fqarr[@]} -eq 0 ]]; then
  err_exit "could not find $fdir/${sid}*.f*q.gz"
fi

fn=$sid

params="-x $hsref -1 ${fqarr[0]} -2 ${fqarr[1]} -k 30 2>${fn}.align_summary.txt"

## - search both strands for now --
#if [[ $dataset != bsp1* ]]; then
#  params="--rna-strandness RF $params"
#  #ribo-zero samples
#fi
cd $outpath || err_exit "failed at: cd $outpath"

cram=$fn.cram
bam=$fn.bam 

if [[ -f $cram ]]; then
  err_exit " $cram already exists in $PWD!"
fi

/bin/rm -f $bam ## WARNING! existing BAM will be lost

rlog=$sid.log

## ---- start main task execution
echo "$line" | tee $rlog

echo ">Task ${jobid}.${taskid} on $host:$outpath" | tee -a $rlog
echo "["$(date '+%m/%d %H:%M')"] starting:" | tee -a $rlog
cmd="hisat2 -p 4 --phred33 --min-intronlen 20 $params |\
 samtools view -b -o $bam -"
echo -e $cmd | tee -a $rlog

tmpsrt=$tmpdir/$fn.srt_tmp
/bin/rm -f ${tmpsrt}*

eval "$cmd" |& tee -a $rlog

## ---- sort and convert to BAM
echo '['$(date '+%m/%d %H:%M')"] start sorting/conversion to CRAM" | tee -a $rlog
cmd="samtools sort -T $tmpsrt -l 6 -m 8G --no-PG -@ 4 $bam | \
 scramble -B -I bam -O cram -8 -r $gref -X small -t 4 -n -! - $cram"
echo -e "$cmd" | tee -a $rlog

eval "$cmd" |& tee -a $rlog

if [[ $? -eq 0 ]]; then
 rm -f $bam
else
 echo "sort|scramble error exit detected!" | tee -a $rlog
 exit 1
fi
echo '['$(date '+%m/%d %H:%M')"] done." | tee -a $rlog
