#!/usr/bin/env bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G

## submit with :
##  sbatch --job-name=kallistoHTT --array=1-487%25 ./run_kallisto.sh

export PATH=/dcs04/lieber/lcolladotor/dbDev_LIBD001/sw/bin:$PATH
function err_exit {
 echo -e "Error: $1"
 exit 1
}

fqdir=/dcs04/lieber/lcolladotor/chessBrain_LIBD4085/raw-data/libd_bsp3
ktx=gencode25_xHTT.kix
id=$SLURM_ARRAY_TASK_ID
if [[ -z "$id" ]]; then
 id="$1"
 if [[ -z "$id" ]]; then
   err_exit "no task ID given or found in the environment!"
 fi
fi
# get RNum based on line number
sid=$(sed "${id}q;d" samples.lst)
if [[ -z $sid ]]; then
  err_exit  "no RNum found!"
fi
fqfiles=( $(ls $fqdir/$sid*.f*q.gz) )
if [[ -z "${fqfiles[1]}" ]]; then
 err_exit "no paired fastq files found as $fqdir/$sid*.f*q.gz"
fi
#echo "input files: ${fqfiles[@]}"
odir=$sid
rlog=$odir/run.log
mkdir -p $odir || err_exit "failed to create $odir"
psbam=$odir/pseudoalignments.bam
cmd="kallisto quant -t 8 -i $ktx --rf-stranded -o $odir --pseudobam ${fqfiles[@]}"
echo -e ">task_$id : running \n$cmd" > $rlog
eval "$cmd" >> $rlog 2>&1
samtools view -h $psbam | grep -P '[:\t]HTT-' | samtools sort | samtools view -b -o $odir/htt_alns.bam
samtools index $odir/htt_alns.bam
echo "done." >> $rlog
