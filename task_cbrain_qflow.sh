#!/usr/bin/env bash
#$ -cwd
#$ -o logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -e logs/$JOB_NAME.$JOB_ID.$TASK_ID.log
#$ -l h_fsize=900G
#$ -l mem_free=8G,h_vmem=12G
#$ -pe local 4

#### Script for generating QC and expression data (counts) from CHESS-Brain 
#### CRAM alignments and FASTQ files (the latter only needed for transcript 
#### abundance estimation with kallisto)

# Run with: 
#   qsub -N qfdlpfc -t 1-56 task_cbrain_qflow.sh taskdb.tfa
#

# taskdb.tfa expected line format:
# >3 sampleID aln_dir
# e.g. 
# >3 R12408 aln/bsp2_dlpfc
# if aln_dir is a relative path, it must be relative to the project
# directory base path (~/work/cbrain)
## Note: this script only uses the sampleID to locate the cram file as
##        aln_dir/sampleID/sampleID*.cram
## the last directory part in aln_dir is taken as the dataset name (ds)
#  and it is used as the output sub-directory
# Output will be in a dataset sub-directory in the working directory, 
# one subdir per sample, e.g. 
##    bsp2_dlpfc/R12408/....
#
### Requires: featureCounts v2.0.3, regtools v0.5.33, kallisto 

### runing with parallel:
## parallel --delay .01 -j 8 task_cbrain_fcounts.sh bsp2dlpfc_sel.tfa {1} ::: {1..299}

function err_exit {
 echo -e "Error: $1"
 exit 1
}
host=$(hostname -s)
onjh=''
if [[ $host == compute-* ]]; then
 onjh=1
fi 
basedir="$HOME/work/cbrain"; #base project directory
gfa="$basedir/ref/hg38mod_noPARs.fa"
## on salz:
g41dir="$HOME/work/cb3/gencode41"
if [[ $onjh ]]; then
  g41dir="$HOME/work/cbrain/gencode41"
fi

qcgtf="$g41dir/ref/gencode41.nri.collapsed.gtf" # for rnaseqc
saf="$g41dir/ref/gencode41.nri.flattened.saf" # for gene counts
xgtf="$g41dir/ref/gencode41.nri.exonsOnly.gtf"
ann="$g41dir/ref/gencode41.nri.gtf"
ktx="$g41dir/ref/gencode41.nri.ktx"
for f in $gfa $qcgtf $saf $xgtf $ann $ktx; do
 if [[ ! -f $f ]]; then
    err_exit "cannot find $f"
 fi
done

# requires a task db pseudo-fasta file as the 1st arg
fdb="$1"
if [[ -z "$fdb" ]]; then
  err_exit "no cdb task db given!"
fi
if [[ ! -f "$fdb" ]]; then
  err_exit "$fdb file not found!"
fi
if [[ ! $fdbidx =~ \.cidx$ ]]; then
  fdb="$fdb.cidx"
  if [[ ! -f "$fdb" ]]; then
    err_exit "$fdb file must exist!"
  fi
fi

taskid="$2" # the task# could be given directly (e.g. by parallel)
jobid="$$"
if [[ -z $taskid ]]; then
   taskid=$SGE_TASK_ID
   jobid=$JOB_ID
   if [[ -z "$taskid" ]]; then
     err_exit "no task index number given or found in SGE_TASK_ID!"
   fi
fi

tmpdir="/dev/shm/qflow/$jobid/$taskid"
if [[ $onjh ]]; then
  tmpdir=$MYSCRATCH/qflow/$jobid/$taskid
fi

pwd=$(pwd -P) # current directory, absolute path

if [[ ! -d $tmpdir ]]; then
  mkdir -p $tmpdir || err_exit "failed to create $tmpdir"
fi

line=$(cdbyank -a $taskid $fdb)
#expected format:       1        2
#                  >6 R11234 aln/bsp2_dlpfc
if [[ $line != '>'* ]]; then
 err_exit "got invalid line for $taskid:\n$line"
fi

t=( $line )
sid=${t[1]} # sampleID
alndir=${t[2]}  # base directory for aln files
ds="${alndir##*/}" #last directory in the path
fqdir="" # fastq path for dataset, special deconvo case
if [[ $ds == 'aln' ]]; then 
  #polyA_vs_RiboZero base directory, grandparent dir is the dataset
  ud="${alndir%/*}"
  ds="${ud##*/}"
  if [[ $ds == 'ribo'* || $ds == 'poly'* ]]; then
    ud="${ud%/*}"
    ds="${ud##*/}"
  fi
fi
echo "Dataset is: $ds"

if [[ $ds == 'deconvo' ]]; then
  #deconvo dataset, fastq path is tricky
  fqbase="${alndir%/*}" #parent directory of aln dir
  fqpre=$(echo "$sid" | cut -f1 -d_)
  fqrest=$(echo "$sid" | cut -f2- -d_)
  fqdir="$fqbase/$fqpre/$fqrest"
fi

if [[ $alndir != '/'* ]]; then
 alndir="$basedir/$alndir" # convert into absolute path
fi

fcram=$(ls $alndir/$sid/$sid*.cram | head -1);
if [[ ! -f $fcram || $(stat -c%s $fcram) -lt 100000 ]]; then
  err_exit "aln file not found ($fcram)"
fi
if [[ -z $fqdir ]]; then
  if [[ $ds == 'bsp1' ]]; then
     fqdir="$basedir/fastq/libd_bsp1/fastq"
     if [[ $onjh ]]; then
       fqdir="$basedir/fastq/libd_bsp1/$sid"
     fi
   elif [[ $ds == 'bsp2_'* ]]; then
     r=$(echo $ds | cut -f2 -d_)
     fqdir="$basedir/fastq/libd_bsp2/fastq_$r"
     if [[ $onjh ]]; then
       fqdir="$basedir/fastq/libd_bsp2_$r"
     fi
   elif [[ $ds == 'bsp3' ]]; then
     fqdir="$basedir/fastq/libd_bsp3/fastq"
     if [[ $onjh ]]; then
       fqdir="$basedir/fastq/libd_bsp3"
     fi
   else
    fqdir="$basedir/fastq/$ds"
  fi
fi
fqfiles=( $(ls $fqdir/$sid*.f*q.gz) )
if [[ -z "${fqfiles[1]}" ]]; then
 err_exit "no paired fastq files found as $fqdir/$sid*.f*q.gz"
fi

odir="$ds/$sid"
if [[ ! -d $odir ]]; then
  mkdir -p $odir || err_exit "failed to create $odir"
fi

cd $odir || err_exit "failed at: cd $odir"

rlog=$sid.fcounts.log

echo ">Task ${jobid}.${taskid} starting in $host:$pwd/$odir ["$(date '+%m/%d %H:%M')"]" | tee $rlog
if [[ ! -f $fcram.'.crai' ]]; then
  samtools index $fcram
fi

sampipe="samtools view --threads 2 --input-fmt-option filter='rname=~\"^chr[0-9MXY]+$\"' -u -T $gfa $fcram |"
cmd=''
rqc="rqc"; # rnaseqc output subdir
fmet="$rqc/$sid.metrics.tsv"

if [[ ! -f $fmet || $(stat -c%s $fmet) -lt 200 ]]; then
  cmd="$sampipe rnaseqc $qcgtf - rqc --mapping-quality=60 --coverage -s $sid"
  echo -e "running rnaseqc:\n$cmd" | tee -a $rlog
  eval "$cmd" |& tee -a $rlog
fi

fstrand='strandness.tab'
## adding strandess info if not found
if [[ ! -f $fstrand ]]; then
 fgrep 'Sense Rate' $fmet | \
   perl -ne '($i,$v)=m/([12]) Sense Rate\s+([\d\.]+)$/;$e[$i]=$v; 
    END {$d=$e[1]-$e[2]; print abs($d)<0.2 ? "unstranded" : $d<0 ? "reverse" : "forward";
    print "\t$e[1]\t$e[2]\n"}' > $fstrand
fi

## adding Mito_mapped metrics if not found 
if ! fgrep -q Mito_mapped $fmet; then
  echo "adding Mito metrics" | tee -a $rlog
  mtcount=$(samtools view -T $gfa -u $indir/$fcram chrM | samtools stats - | \
  grep -oP 'reads mapped:\s+\K\d+')
  echo -e "Mito_mapped\t$mtcount" >> $fmet
fi

strand=$(cut -f1 $fstrand) #assume we have this after running rnaseqc
sflag=0
if [[ $strand == 'forward' ]]; then
   sflag=1
 elif [[ $strand == 'reverse' ]]; then
   sflag=2
 else
   sflag=0
fi
fexsum=$sid.exon.fcounts.summary
if [[ ! -f $fexsum || $(stat -c%s $fexsum) -lt 200 ]]; then
 ## --- Getting exon counts:
 cmd="$sampipe featureCounts --tmpDir $tmpdir -s $sflag -p -T 2 -O -f \
 -a $xgtf -o $sid.exon.fcounts"
 echo -e " running featureCounts for exons:\n$cmd" | tee -a $rlog
 eval "$cmd" |& tee -a $rlog &
fi

fgsum=$sid.gene.fcounts.summary
if [[ ! -f $fgsum || $(stat -c%s $fgsum) -lt 200 ]]; then
## --- Getting gene counts:
 cmd="$sampipe featureCounts --tmpDir $tmpdir -s $sflag -p -T 2 -F SAF -a $saf \
 -o $sid.gene.fcounts"
 echo -e "running featureCounts for genes:\n$cmd" | tee -a $rlog
 eval "$cmd" |& tee -a $rlog &
fi

frbed=$sid.regtools.bed

## generate regtools bed
if [[ ! -f $frbed || $(stat -c%s $frbed) -lt 1200 ]]; then
  cmd="$sampipe regtools junctions extract -s XS -c $sid.regtools.ctab -o $frbed -"
  echo -e "running regtools extract:\n$cmd" | tee -a $rlog
  eval "$cmd" |& tee -a $rlog &  
fi

## also start stringtie -eB
edir=strge
sout=$edir/$sid.strge.gtf

if [[ ! -f $sout ]]; then
  mkdir -p $edir
  cmd="stringtie --cram-ref $gfa -o $sout -e -b$edir -A $edir/$sid.genes.tab -G $ann $fcram"
  echo -e "running stringtie -eB:\n$cmd" | tee -a $rlog
  eval "$cmd" |& tee -a $rlog &  
fi

wait
## -- remove tmpdir files -- should no longer be needed
#find $tmpdir -type f -delete
if [[ $tmpdir == '/dev/shm'* || $tmpdir == *myscratch* ]]; then
 /usr/bin/rm -f $tmpdir/*
fi

fjann=$sid.regtools.ann
if [[ ! -f $fjann ]]; then
  cmd="regtools junctions annotate -o $fjann -S $frbed $gfa $ann"
  echo -e "running regtool jx annotation:\n$cmd" | tee -a $rlog
  eval "$cmd" |& tee -a $rlog &
fi

## also run kallisto if needed
fk=kallisto/abundance.tsv
if [[ ! -f $fk || $(stat -c%s $fk) -lt 4200 ]]; then
 mkdir -p kallisto
 strand=$(cut -f1 $fstrand)
 if [[ $strand == 'forward' ]]; then
    ksflag='--fr-stranded'
  elif [[ $strand == 'reverse' ]]; then
    ksflag='--rf-stranded'
  else
    ksflag=''
 fi
 fqfiles=( $(ls $fqdir/$sid*.f*q.gz) ) #should have been checked before
 #using 3 CPUs (-t 3)
 cmd="kallisto quant -t 3 -i $ktx $ksflag -o kallisto ${fqfiles[@]}"
 echo -e "running kallisto:\n$cmd" | tee -a $rlog
 eval "$cmd" |& tee -a $rlog &
fi

wait

if [[ $tmpdir == '/dev/shm'* || $tmpdir == *myscratch* ]]; then
 /usr/bin/rm -rf $tmpdir
fi
echo -e "["$(date '+%m/%d %H:%M')"]\ttask done [$taskid]." | tee -a $rlog
