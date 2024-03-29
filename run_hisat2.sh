#!//usr/bin/env bash
#
refbase=/data/gpertea/work/ref/hg38mod
hsrefM=$refbase/hisat2_hg38mod/hg38mod_noPARs
hsrefF=$refbase/hisat2_hg38mod_noY/hg38mod_noY
gref=$refbase/hg38mod_noPARs.fa

#path having libd_bsp1 subdirs with fastq dirs 
indir=/data2/cbrain/data

## final .cram files will be written in subdirectories here:
outdir=/data/gpertea/work/cbrain-alns/new

# could be on /dev/shm/$USER/$SLURM_JOB_ID/
#tmpdir=$HOME/scratch/$SLURM_JOB_ID/$SLURM_ARRAY_TASK_ID
tmpdir=/data2/cbrain/data/tmp

function err_exit {
 echo -e "Error: $1"
 exit 1
}

mkdir -p $tmpdir || err_exit "failed to create $tmpdir"

fdb="$1"
if [[ -z "$fdb" ]]; then
  err_exit "no cdb task database file given!"
fi

if [[ ! -f "$fdb" ]]; then
  err_exit "$fdb file not found!"
fi

fdb="$fdb.cidx"

if [[ ! -f "$fdb" ]]; then
    err_exit "$fdb.cidx file not found (with cdbfasta)!"
fi

#id=$SGE_TASK_ID
id=$SLURM_ARRAY_TASK_ID
if [[ -z "$id" ]]; then
 id="$2"
 if [[ -z "$id" ]]; then
   err_exit "no task ID given or found in the environment!"
 fi
fi

rnum=""
sub=""
line=$(cdbyank -a $id $fdb)
#format is: >342 relpath rnum
if [[ $line != '>'* ]]; then
 err_exit "invalid line pulled for $id:\n$line"
fi
t=( $line )
sub=${t[1]}
rnum=${t[2]}
hsref=$hsrefM
if [[ ${t[3]} == "F" ]]; then
 hsref=$hsrefF
fi
osub="${sub/#libd_/}"
osub="${osub//\/fastq/}"
outpath="$outdir/$osub/$rnum"
mkdir -p "$outpath" || err_exit "failed at mkdir -p $outpath"

#cd "$indir/$sub" || err_exit "failed at cd $indir/$sub"
arr=($(ls  $indir/$sub/${rnum}_*_R1*.f*q*)) 
# if that fails:  arr=($(ls  ${rnum}_*_1.f*q*))
if [[ ${#arr[@]} -eq 0 ]]; then
  err_exit "could not find $indir/$sub/${rnum}_*_R1*.f*q* in $indir/$sub"
fi
rid=${arr[0]}
rid="${rid##*/}" # get the file name
re="(R[0-9]+_[0-9A-Z]+)*"
if [[ $rid =~ $re ]]; then 
  rid="${BASH_REMATCH[1]}"
else 
 err_exit "could not figure out a RNum ID from $rid"
fi

printf -v jn '%s,' "${arr[@]}"
m1="${jn%,}"

arr2=($(ls  $indir/$sub/${rnum}_*_R2*.f*q*))
if [[ ${#arr2[@]} -eq 0 ]]; then
  err_exit "could not find $indir/$sub/${rnum}_*_R2*.f*q* in $indir/$sub"
fi

printf -v jn '%s,' "${arr2[@]}"
m2="${jn%,}"

if [[ ! ${#arr[@]} -eq ${#arr2[@]} ]]; then
 err_exit " arrays are not equal length:\n${arr[@]}\n${arr2[@]}"
fi

fn=$rid

params="-x $hsref -1 $m1 -2 $m2 2>${fn}.align_summary.txt"

if [[ $osub != bsp1* ]]; then
  params="--rna-strandness RF $params"
  #ribo-Zero sample
fi
#echo -e "will running hisat2 with params: \n$params"
#echo "Output dir: $outpath/ with file name base: $fno"
# -- file info collected, ready to run the script ---

cd $outpath || err_exit "failed at cd $outpath"

# echo '['$(date '+%m/%d %H:%M')"] starting.."
cram=$fn.cram
#bam could be in a different/temp path 
#  e.g. /dev/shm/$USER/$SLURM_JOB_ID/$SLURM_ARRAY_TASK_ID
bam=$fn.bam 

if [[ -f $cram ]]; then
  err_exit " $cram already exists in $PWD!"
fi
rlog=$fn.rlog
echo '['$(date '+%m/%d %H:%M')"]:" > $rlog
echo -e "hisat2 --mm -p 4 --phred33 --min-intronlen 20 $params | \
samtools view -b -o $fn.bam -" >> $rlog
tmpsrt=$tmpdir/$fn.srt_tmp

hisat2 --mm -p 4 --phred33 --min-intronlen 20 $params | \
 samtools view -b -o $bam
## now sort it and convert to BAM
 echo '['$(date '+%m/%d %H:%M')"] running:" >> $rlog
 echo "samtools sort -T $fn.srttmp -l 4 -m 7G -@ 4 $bam | \
scramble -B -I bam -O cram -8 -r $gref -X small -t 4 -! - $cram" >> $rlog
 samtools sort -T $tmpsrt -l 4 -m 7G --no-PG -@ 4 $bam | \
  scramble -B -I bam -O cram -8 -r $gref -X small -t 4 -! - $cram
if [[ $? -eq 0 ]]; then
 rm -f $bam
else
 echo "sort|scramble error exit detected!" >> $rlog
fi
echo '['$(date '+%m/%d %H:%M')"] done." >> $rlog

#echo "$cmd" >> $ofn.info
#eval "$cmd" >> $ofn.log 2>&1
