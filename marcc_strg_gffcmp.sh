#!//usr/bin/env bash
#SBATCH -p shared
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=24G

## how to submit this on MARCC
# sbatch --array=1-400 --mail-type=FAIL --mail-user=geo.pertea@gmail.com ./run_strg_gffcmp.sh uniq_rnum_info.cfa
#

refbase=/work-zfs/ssalzbe1/chess-brain/ref
gref=$refbase/hg38mod_noPARs.fa
ann=$refbase/hg38mod_noPARs.gencode37xrefseq_nopseudo.gff
gencode=$HOME/work/ref/gencode.v40.annotation.gtf

## note: both alns and stringtie directories have libd_ prefix and the intermediate
## folder '/fastq' REMOVED, and RNum is created instead for each sample data
## so the mappings are like such under ~/cbrain/:
#
#data/libd_bsp1/fastq/R*.fastq.gz          aln/bsp1/R*/*.cram
#data/libd_bsp2/fastq_dlpfc/R*.fastq.gz    aln/bsp2_dlpfc/R*/*.cram
#data/libd_bsp2/fastq_hippo/R*.fastq.gz    aln/bsp2_hippo/R*/*.cram
#data/libd_bsp3/fastq/R*.fastq.gz          aln/bsp3/R*/*.cram

###- for stringtie output, they will be:
# strg/bsp1/R*/*.gtf
# strg/bsp2_dlpfc/R*/*.gtf
# strg/bsp2_hippo/R*/*.gtf
# strg/bsp3/R*/*.gtf

#base dir for alignment data
mydir=/scratch/groups/ssalzbe1/gpertea

alndir=$mydir/cbrain-aln
## output files are written in dataset_path/subdirectories here:
outdir=$mydir/cbrain-stringtie

function err_exit {
 echo -e "Error: $1"
 exit 1
}
# requires a task pseudo-fasta file

fdb="$1"
if [[ -z "$fdb" ]]; then
  err_exit "no cdb task db given!"
fi
if [[ ! -f "$fdb" ]]; then
  err_exit "$fdb file not found!"
fi
fdb="$fdb.cidx"
if [[ ! -f "$fdb" ]]; then
    err_exit "$fdb file not found!"
fi

id=$SLURM_ARRAY_TASK_ID # or $SGE_TASK_ID on SGE (JHPCE)
if [[ -z "$id" ]]; then
 id="$2"
 if [[ -z "$id" ]]; then
   err_exit "no task ID found or given!"
 fi
fi

line=$(cdbyank -a $id $fdb)
#format is:       1              2   3   4      5    6   7
#           >6 libd_bsp1/fastq R2828 M DLPFC Schizo AA 52.02

if [[ $line != '>'* ]]; then
 err_exit "invalid line pulled for $id:\n$line"
fi
t=( $line )
fdir=${t[1]}
sid=${t[2]} # RNum for LIBD data
#rnum=${t[3]}
sx=${t[3]}
reg=${t[4]}
#dx=${t[5]}
#race=${t[6]}
#age=${t[7]}

## -- get the new aln and stringtie path pattern:
dsdir=${fdir//libd_/}
dsdir=${dsdir//\/fastq/}/$sid

proto="ribo"
if [[ $dsdir = *"bsp1"* ]]; then 
 proto="polyA"
fi

## for strintie we assume the CRAM file exists already
cram=$(ls $alndir/$dsdir/${sid}*.cram)
if [[ ! -f $cram ]]; then
  err_exit "could not find CRAM file: $cram"
fi

outpath="$outdir/$dsdir"
mkdir -p "$outpath" || err_exit "failed at: mkdir -p $outpath"
cd "$outpath" || err_exit "failed at: cd $outpath"

params="--cram-ref $gref -G $ann -o $sid.gtf $cram"

#/bin/rm -f $rnum.gtf *.ctab

rlog=$sid.log

## ---- start main task execution
echo "$line" | tee $rlog
echo "task ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID} on $HOSTNAME $outpath" | tee -a $rlog
echo "["$(date '+%m/%d %H:%M')"] starting command:" | tee -a $rlog
cmd="stringtie $params"
echo -e $cmd | tee -a $rlog

eval "$cmd" |& tee -a $rlog

if [[ $? -eq 0 ]]; then
 # rm -f $bam
 echo '['$(date '+%m/%d %H:%M')"] initial stringtie run done." | tee -a $rlog
else
 echo '['$(date '+%m/%d %H:%M')"] error exit detected!" | tee -a $rlog
 exit 1
fi

#annotate with gffcmp vs Gencode 
cmd1="gffcompare -r $gencode -D -o cmp $sid.gtf"
cmd2="trmap $gencode $sid.gtf -J -o $sid.jtab"
eval "$cmd1" &
pids+="$! "
eval "$cmd2" &
pids+="$! "

# also run with -eB to generate stringtie .ctab files on its own output (junction counts)
mkdir -p 'eB'
eparams="--cram-ref $gref -G $ann -o eB/$sid.gtf -e -b eB $cram"
echo "["$(date '+%m/%d %H:%M')"] starting command:" | tee -a $rlog
cmd="stringtie $eparams"
echo -e $cmd | tee -a $rlog
eval "$cmd" |& tee -a $rlog &
pids+="$! "

for pid in $pids; do
  wait $pid
  if [[ $? -eq 0 ]]; then
    # rm -f $bam
    echo "cmd finished."
  else
    echo '['$(date '+%m/%d %H:%M')"] error exit detected for pid $pid!" | tee -a $rlog
    exit 1
  fi
done

echo '['$(date '+%m/%d %H:%M')"] done." | tee -a $rlog

#if [[ $? -eq 0 ]]; then
# # rm -f $bam
# echo '['$(date '+%m/%d %H:%M')"] done." | tee -a $rlog
#else
# echo '['$(date '+%m/%d %H:%M')"] error exit detected!" | tee -a $rlog
# exit 1
#fi
