#!//usr/bin/env bash
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
#   qsub -N qfdlpfc -t 1-56 redo_regtools.sh taskdb.tfa
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
### Requires: regtools v0.5.33

### runing with parallel:
## parallel --delay .01 -j 8 redo_regtools.sh bsp2dlpfc_sel.tfa {1} ::: {1..299}

function err_exit {
 echo -e "Error: $1"
 exit 1
}
host=$(hostname -s)
onjh=''
onsalz=''
if [[ $host == compute-* ]]; then
 onjh=1
fi
if [[ $host == salz* ]]; then
 onsalz=1
fi 

basedir="$HOME/work/cbrain"; #base project directory
gfa="$basedir/ref/hg38mod_noPARs.fa"
g41dir="$HOME/work/cbrain/gencode41"
## on salz:
if [[ $onsalz ]]; then
  g41dir="$HOME/work/cb3/gencode41"
fi

ann="$g41dir/ref/gencode41.nri.gtf"
for f in $gfa $ann; do
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

pwd=$(pwd -P) # current directory, absolute path

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

if [[ $alndir != '/'* ]]; then
 alndir="$basedir/$alndir" # convert into absolute path
fi

fcram=$(ls $alndir/$sid/$sid*.cram | head -1);
if [[ ! -f $fcram || $(stat -c%s $fcram) -lt 100000 ]]; then
  err_exit "aln file not found ($fcram)"
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

frbed=$sid.regtools.bed
## redo - delete existing BED !
/bin/rm -f $frbed
## -s 0 means ALWAYS use the XS tag to get the strand
## generate regtools bed
#if [[ ! -f $frbed || $(stat -c%s $frbed) -lt 1200 ]]; then
  cmd="$sampipe regtools junctions extract -s 0 -c $sid.regtools.ctab -o $frbed -"
  echo -e "running regtools extract:\n$cmd" | tee -a $rlog
  eval "$cmd" |& tee -a $rlog
#fi

fjann=$sid.regtools.ann
## always overwrite regtools.ann
#if [[ ! -f $fjann ]]; then
  cmd="regtools junctions annotate -o $fjann -S $frbed $gfa $ann"
  echo -e "running regtool jx annotation:\n$cmd" | tee -a $rlog
  eval "$cmd" |& tee -a $rlog
#fi

echo -e "["$(date '+%m/%d %H:%M')"]\ttask done [$taskid]." | tee -a $rlog
