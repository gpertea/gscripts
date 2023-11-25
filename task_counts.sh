#!/usr/bin/env bash
##x mem=18G
##x cpus=8

#### Script for generating QC with rnaseqc and expression data (counts) 
#### using featureCounts (gene, exon), regtools (junctions), salmon (tx)
### requires: featureCounts v2.0.3, regtools v0.5.33, salmon

## requires a merged manifest file (fq_merge_manifest.pl)
## sampleID folders with cram files are expected in the current directory

## run with:
# arx sub -m18G -c8 -a1- --cfg ../counts.cfg task_counts.sh ../merged.manifest
refbase="/dcs04/lieber/lcolladotor/dbDev_LIBD001/ref"
gfa=${GENOME_FA:-$refbase/fa/assembly_hg38_gencode_v25_main.fa}
gann=${GENOME_ANN:-$refbase/gtf/gencode25.main.gtf}
fcexon=${FCOUNTS_EXON:-$refbase/gtf/gencode25.main.exons.gtf}
fcgene=${FCOUNTS_GENE:-$refbase/gtf/gencode25.main.flattened.saf}
refqc=${RQC_REF:-$refbase/gtf/gencode25.main.collapsed.gtf}
salmidx=${SALMON_IDX:-$refbase/salmon_idx/gencode25.main}

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

rlog=counts_t${taskid}.tlog

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

for f in $gfa $gann $fcexon $fcgene $refqc $salmidx/pos.bin ; do
 if [[ ! -f $f ]]; then
    err_exit "cannot find $f"
 fi
done

##path having libd_bsp1, cmc1 etc. dirs 
# only needed when relative paths are given in taskdb
pwd=$(pwd -P) # current directory, absolute path

## on JHPCE use $MYSCRATCH
if [[ $host == transfer-* || $host == compute-* ]]; then
 tmpdir=$MYSCRATCH/fq2cram/$jobid/$taskid
else
 if [[ $host == srv05 ]]; then
   tmpdir=/dev/shm/${USER}-tmp/$jobid/$taskid
 else
   tmpdir=$pwd/tmp/$jobid/$taskid
 fi
fi

mkdir -p $tmpdir || err_exit "failed to create $tmpdir"

line=$(linix $fdb $taskid | cut -f1,3,5)
#expected format:     read_1.fq[,read_1b.fq.gz,...] read_2.fq[,read_2b.fq.gz,...] sampleID
t=( $line )
fn1=${t[0]} # path(s) to read_1 fastq.gz files
fn2=${t[1]}  # sampleID_1.fastq.gz | RNum_*_R1_*.fastq.gz
oid=${t[2]}  # output directory = base sampleID
if [[ -z $oid ]]; then
 err_exit "could not parse base sampleID!"
fi

fqs1=( ${fn1/,/ } )
fqs2=( ${fn2/,/ } )

if [[ ! -f "${fqs1[0]}" ]]; then
   err_exit "FASTQ file ${fqs1[0]} cannot be found!"
fi

outpath=$oid
cd $outpath || err_exit "failed at: cd $outpath"
crams=( $(ls ${oid}*.cram) ) # could be one for each flowcell
ncram=${#crams[@]}
if ((ncram==0)); then
  err_exit "could not find cram files in $outpath"
fi

sid=$oid #for these counts, unify sid (merge flowcells)

echo "["$(date '+%m/%d %H:%M')"] task ${jobid}.${taskid} starting on $host:${PWD}" | tee -a $rlog

fpri=$oid.pri.cram

run=''

## Round 1: get primary alignments, run Salmon to get transcript counts
#sampipe="samtools view --threads 2 --input-fmt-option filter='rname=~\"^chr[0-9MXY]+$\"' -u $cref ${faln[@]} |"
sampipe="samtools view --threads 2 -T $gfa -u ${crams[0]} |"
sampri="samtools view --threads 2 -F 260 -T $gfa -o $fpri --write-index -O cram,version=3.1 ${crams[0]}"
if [[ $ncram -gt 1 ]]; then 
 sampipe="samtools merge -T $gfa --threads 2 -u -o - ${crams[@]} |"
 ## | samtools view --input-fmt-option filter='rname=~\"^chr[0-9MXY]+$\"' -u -|"
 sampri="$sampipe samtools view --threads 2 -F 260 -T $gfa -o $fpri --write-index -O cram,version=3.1 -"
fi

## write only primary alignments
if [[ ! -f $fpri ]]; then
    echo "building primary alignments file: $sampri"
    run="${run}p"
    eval "$sampri" |& tee -a $rlog &
fi


## start Salmon
fsalm="salmon/quant.sf"
if [[ ! -f $fsalm ]]; then
  cmd="salmon quant -p 6 -lA -1 <(gunzip -c ${fqs1[@]}) -2 <(gunzip -c ${fqs2[@]}) \
   -i $salmidx --gcBias -q --numGibbsSamples 20 --thinningFactor 40 -d -o salmon >& salmon.log"
  echo -e "running salmon:\n$cmd" | tee -a $rlog
  run="${run}s"
  eval "$cmd" |& tee -a $rlog &
fi

if [[ $run ]]; then
 wait
fi

######################
## Round 2:  rnaseqc, regtools, featureCounts
run=''

rqc="rqc"; # rnaseqc output subdir
fmet="$rqc/$oid.metrics.tsv"
if [[ ! -f $fmet || $(stat -c%s $fmet) -lt 200 ]]; then
  cmd="$sampipe rnaseqc $refqc - rqc --mapping-quality=60 --coverage -s $oid"
  echo -e "running rnaseqc:\n$cmd" | tee -a $rlog
  eval "$cmd" |& tee -a $rlog &
  run="${run}q"
fi

## let regtools and featureCounts run on the primary alignments only
sampipe="samtools view --threads 2 -T $gfa -u $fpri |"

fjtab=$oid.regtools.ctab
## generate regtools counts
if [[ ! -f $fjtab || $(stat -c%s $fjtab) -lt 1200 ]]; then
  cmd="$sampipe regtools junctions extract -s 0 -c $fjtab -o $oid.regtools.bed -"
  echo -e "running regtools extract:\n$cmd" | tee -a $rlog
  run="${run}j"
  eval "$cmd" |& tee -a $rlog &  
fi

## run featureCounts for exon and gene counts, on primary alignments only
## -s $sflag omitted (default: 0=unstranded)
fexsum=$sid.exon.fcounts.summary
if [[ ! -f $fexsum || $(stat -c%s $fexsum) -lt 200 ]]; then
 ## --- Getting exon counts:
 cmd="$sampipe featureCounts --tmpDir $tmpdir -p -T 2 -O -f \
 -a $fcexon -o $sid.exon.fcounts"
 echo -e " running featureCounts/exon:\n$cmd" | tee -a $rlog
 run="${run}e"
 eval "$cmd" |& tee -a $rlog &
fi

fgsum=$sid.gene.fcounts.summary
if [[ ! -f $fgsum || $(stat -c%s $fgsum) -lt 200 ]]; then
## --- Getting gene counts:
 cmd="$sampipe featureCounts --tmpDir $tmpdir -p -T 2 -F SAF -a $fcgene \
 -o $sid.gene.fcounts"
 echo -e "running featureCounts/gene:\n$cmd" | tee -a $rlog
 run="${run}g"
 eval "$cmd" |& tee -a $rlog &
fi

## wait for background tasks to finish
if [[ $run ]]; then
 wait
fi


## adding Mito_mapped metrics if not found
if ! fgrep -q Mito_mapped $fmet; then
  echo "adding Mito metrics" | tee -a $rlog
  mtcount=0
  samtools idxstats $fpri > $oid.pri.chrstats
  mtcount=$(grep '^chrM\b' $oid.pri.chrstats | cut -f3)
  echo -e "Mito_mapped\t$mtcount" >> $fmet
fi

## regtools annotate - not needed
#fjann=$sid.regtools.ann
#if [[ ! -f $fjann ]]; then
#  cmd="regtools junctions annotate -o $fjann -S $oid.regtools.bed $gfa $gann"
#  echo -e "running regtool jx annotation:\n$cmd" | tee -a $rlog
#  run="${run}J"
#  eval "$cmd" |& tee -a $rlog &
#fi

fstrand='strandness.tab'
## adding strandess info if not found
if [[ ! -f $fstrand ]]; then
 fgrep 'Sense Rate' $fmet | \
   perl -ne '($i,$v)=m/([12]) Sense Rate\s+([\d\.]+)$/;$e[$i]=$v; 
    END {$d=$e[1]-$e[2]; print abs($d)<0.2 ? "unstranded" : $d<0 ? "reverse" : "forward";
    print "\t$e[1]\t$e[2]\n"}' > $fstrand
fi


echo -e "["$(date '+%m/%d %H:%M')"]\tcounts task done [$taskid]." | tee -a $rlog

