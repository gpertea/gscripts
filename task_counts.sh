#!/usr/bin/env bash
##x mem=18G
##x cpus=6

#### Script for generating QC with rnaseqc and expression data (counts) 
#### using featureCounts (gene, exon), regtools (junctions), salmon (tx)
### requires: featureCounts v2.0.3, regtools v0.5.33, salmon

## requires a merged manifest file (use fq_merge_manifest.pl on trimmed_samples.manifest)
## sampleID folders with bam and cram files are expected in the current directory

## run with:
# arx sub -m18G -c8 -a1- -j 6 --cfg ../counts.cfg task_counts.sh ../merged.manifest

## if multiple CRAM/BAM files are present for a sample, they will be merged into
## a single CRAM, with RG groups assigned accordingly

refdir=${GREF_DIR:-/dcs04/lieber/lcolladotor/dbDev_LIBD001/ref}
gref=${GENOME_FA:-$refdir/fa/assembly_hg38_gencode_v25_main.fa}
snvs=${SNV:-$refdir/Genotyping/common_missense_SNVs_hg38.bed}
gann=${GENOME_ANN:-$refdir/gtf/gencode25.main.gtf}
fcexon=${FCOUNTS_EXON:-$refdir/gtf/gencode25.main.exons.gtf}
fcgene=${FCOUNTS_GENE:-$refdir/gtf/gencode25.main.flattened.saf}
refqc=${RQC_REF:-$refdir/gtf/gencode25.main.collapsed.gtf}
salmidx=${SALMON_IDX:-$refdir/salmon_idx/gencode25.main}
ncpus=${SALMON_CPUS:-6}
## note 3 more cpus will be used 
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
echo '['$(date '+%m/%d %H:%M')"] task ${jobid}.${taskid} on $host:${PWD}"

# requires a task file as the 1st arg
fdb="$1"
if [[ -z "$fdb" ]]; then
  err_exit "no task file given!"
fi

for f in $gref $snvs $gann $fcexon $fcgene $refqc $salmidx/pos.bin ; do
 if [[ ! -f $f ]]; then
    err_exit "cannot find $f"
 fi
done

##path having libd_bsp1, cmc1 etc. dirs 
# only needed when relative paths are given in taskdb
pwd=$(pwd -P) # current directory, absolute path

## on JHPCE use $MYSCRATCH
if [[ $host == transfer-* || $host == compute-* ]]; then
 tmpdir=$MYSCRATCH/${jobid}_$taskid
else
 if [[ $host == srv05 ]]; then
   tmpdir=/dev/shm/${USER}-${jobid}_$taskid
 elif [[ $host == srv16 ]]; then
   tmpdir=/tmp/scratch/${USER}-${jobid}_$taskid
 else
   tmpdir=$pwd/tmp/${jobid}_$taskid
 fi
fi

mkdir -p $tmpdir || err_exit "failed to create $tmpdir"

line=$(linix $fdb $taskid | cut -f1,3,5)
#expected format:     read_1.fq[,read_1b.fq.gz,...] read_2.fq[,read_2b.fq.gz,...] sampleID
t=( $line )
fn1=${t[0]} # path(s) to read_1 fastq.gz files
fn2=${t[1]}  # sampleID_1.fastq.gz | RNum_*_R1_*.fastq.gz
sid=${t[2]}  # output directory = base sampleID that should merge all the parts
if [[ -z $sid ]]; then
 err_exit "could not parse base sampleID!"
fi

#fqs1=( ${fn1//,/ } )
IFS=',' read -ra fqs1 <<< "$fn1"
#fqs2=( ${fn2//,/ } )
IFS=',' read -ra fqs2 <<< "$fn2"

if [[ ! -f "${fqs1[0]}" ]]; then
   err_exit "FASTQ file ${fqs1[0]} cannot be found!"
fi

outpath=$sid
cd $outpath || err_exit "failed at: cd $outpath"
#crams=( $(ls ${sid}*.cram) ) # could be one per flowcell
## need actual files, not symlinks
crams=( $(find . -maxdepth 1 -type f -name "${sid}*.cram" | sed 's/.\///') )

bams=( $(ls ${sid}*.bam) ) # the unsorted primary alignments to use for featureCounts

ncram=${#crams[@]}
if ((ncram==0)); then
  err_exit "could not find ${sid}*.cram files in $outpath"
fi

if ((ncram>1)); then
  ## if there is more than one, the merged or unsorted one could be among them
  #exclude the unsorted and those likely generated by a previous run
  new_crams=()  # Create an empty array to store filtered results
  for file in "${crams[@]}"; do 
    if [[ "$file" != "${sid}.cram" && "$file" != "${sid}.unsorted.cram" ]]; then
       new_crams+=( "$file" )  # Add valid files to the new array
    fi
  done
  crams=( "${new_crams[@]}" )  # replace the original array
fi

if [[ -f "$rlog" ]]; then
  bk=1
  while [[ -f "$rlog.$bk" ]]; do
    ((bk++))
  done
  mv $rlog "$rlog.$bk"
fi

echo "["$(date '+%m/%d %H:%M')"] task ${jobid}.${taskid} starting on $host:${PWD}" | tee -a $rlog

mcram=$sid.cram # merged CRAM across flowcells/lanes

run=''

# function to check if a file exists and is not a symbolic link
check_cram_nolink() {
  local fcram="$1"  # Get the file path as an argument
  if [[ -e "$fcram" ]]; then
    if [[ -L "$fcram" ]]; then
      unlink "$fcram" 
      unlink "$fcram.crai" 
    else
      err_exit "$fcram exists and is not a symbolic link."
    fi
  fi
}

## salmon can be quite quick, so we start regtools and rnaseqc first in the background
## then we start salmon followed by featureCounts in the main script

## Round 1: rnaseqc running on the CRAM alignments; regtools will also write the merged CRAM and get junction counts, 
#           Salmon to get transcript counts
#sampipe="samtools view --threads 2 --input-fmt-option filter='rname=~\"^chr[0-9MXY]+$\"' -u $cref ${faln[@]} |"
#sampri="samtools view --threads 2 -F 260 -T $gref -o $fpri --write-index -O cram,version=3.1 ${crams[0]}"
teemerge=""
if ((ncram==1)); then
  fn=${crams[0]/.cram/} ## remove .cram
  if [[ "$fn" != "$sid" ]]; then
     check_cram_nolink "$sid.cram"
     ln -s ${crams[0]} "$sid.cram"
     ln -s ${crams[0]}.crai "$sid.cram.crai"
     crams[0]="$sid.cram"
  fi
  crampipe="samtools view -T $gref --threads 2 -u ${crams[0]} |" # for regtools, which does not take newer CRAM files
  bampri="samtools view -F 0x900 --threads 2 -u ${bams[0]} |" # for featureCounts only use primary alignments
else #if [[ $ncram -gt 1 ]]; then 
  rgf=rg.txt
  /bin/rm $rgf
  rcrams=()
  for fc in "${crams[@]}"; do
    rg=${fc/.cram/} # remove cram extension - this will be the RG id
    rg=${rg/_trimmed/}
    rg=${rg/${sid}_/}    
    check_cram_nolink "$rg.cram"
    ln -s "$fc" "$rg.cram"
    ln -s "$fc.crai" "$rg.cram.crai"
    rcrams+=("$rg.cram")
    echo -e "@RG\tID:$rg" >> $rgf
  done
  crams=("${rcrams[@]}")
  crampipe="samtools merge -rh $rgf --threads 2 --reference=$gref -u - ${crams[@]} | "
  if [[ ! -s $mcram ]]; then
    teemerge="tee >(samtools view -O cram,version=3.1 --threads 2 --write-index --reference=$gref -o $mcram -) |"
  fi
  ## | samtools view --input-fmt-option filter='rname=~\"^chr[0-9MXY]+$\"' -u -|"
  bampri="samtools merge --threads 2 -u - ${bams[@]} | samtools view -F 0x900 -u -|"
fi

rqc="rqc"; # rnaseqc output subdir
#crampipe="samtools view --threads 2 -T $gref -u $mcram |"
fmet="$rqc/$sid.metrics.tsv"
if [[ ! -f $fmet || $(stat -c%s $fmet) -lt 200 ]]; then
## -->>> start rnaseqc &&&
  cmd="$crampipe rnaseqc $refqc - rqc --mapping-quality=60 --coverage -s $sid"
  echo -e "running rnaseqc:\n$cmd" | tee -a $rlog
  eval "$cmd" |& tee -a $rlog &
  run="${run}q"
fi

# run=''

fjtab=$sid.regtools.ctab
if [[ ! -f $fjtab || $(stat -c%s $fjtab) -lt 1200 ]]; then
## -->>> start regtools &&&
## this also writes the merged $sid.cram file if that's not the input
  cmd="$crampipe $teemerge regtools junctions extract -s 0 -c $fjtab -o $sid.regtools.bed -"
  echo -e "running regtools extract:\n$cmd" | tee -a $rlog
  run="${run}j"
  eval "$cmd" |& tee -a $rlog &  
fi

fsalm="salmon/quant.sf"
if [[ ! -s $fsalm ]]; then
## -->>> run Salmon (and wait for it to finish)
  cmd="salmon quant -p $ncpus -lA -1 <(gunzip -c ${fqs1[@]}) -2 <(gunzip -c ${fqs2[@]}) \
   -i $salmidx --gcBias -q --numGibbsSamples 30 --thinningFactor 40 -d -o salmon >& salmon.log"
  echo -e "running salmon:\n$cmd" | tee -a $rlog
  #run="${run}s" -- no need, we wait for salmon to finish
  eval "$cmd" |& tee -a $rlog 
fi

## - salmon finished, but regtools and rnaseqc might still be running

## -s $sflag omitted (default: 0=unstranded)
fexsum=$sid.exon.fcounts.summary
if [[ ! -f $fexsum || $(stat -c%s $fexsum) -lt 200 ]]; then
 ## -->>> start featureCounts exon level &&& 
 /bin/rm -f $fexsum 
 cmd="$bampri featureCounts --tmpDir $tmpdir -p -T 3 -O -f \
 -a $fcexon -o $sid.exon.fcounts"
 echo -e " running featureCounts/exon:\n$cmd" | tee -a $rlog
 run="${run}e"
 eval "$cmd" |& tee -a $rlog &
fi

fgsum=$sid.gene.fcounts.summary
if [[ ! -f $fgsum || $(stat -c%s $fgsum) -lt 200 ]]; then
 ## -->>> start featureCounts gene level &&&
 /bin/rm -f $fgsum
 cmd="$bampri featureCounts --tmpDir $tmpdir -p -T 3 -F SAF -a $fcgene \
 -o $sid.gene.fcounts"
 echo -e "running featureCounts/gene:\n$cmd" | tee -a $rlog
 run="${run}g"
 eval "$cmd" |& tee -a $rlog &
fi

## wait for all background tasks to finish (rnaseqc, regtools, featureCounts exon+gene
if [[ $run ]]; then
 wait
fi
##---------------  all processed finished
## merged cram file $mcram must be created
/bin/rm -rf $tmpdir

if [[ ! -s $fmet ]]; then
  echo "rnaseqc output $fmet has zero size! Check $rlog" | tee -a $rlog
  exit 1
fi
if [[ ! -s $fjtab ]]; then
  echo "regtools output $fjtab has zero size! Check $rlog" | tee -a $rlog
  exit 1
fi

if [[ ! -s $fgsum || ! -s $fexsum ]]; then
  echo "one of featureCounts outputs has zero size! Check $rlog" | tee -a $rlog
  exit 1
fi

function get_fsize() {
    local file_path="$1"
    # Follow symbolic links
    if [[ -L "$file_path" ]]; then
        file_path=$(readlink -f "$file_path")
    fi
    # Get file size, providing 0 if an error occurs
    stat -c %s "$file_path" 2>/dev/null || echo 0
}
mcram_crai="$mcram.crai"

## crampipe="samtools view --threads 2 -T $gref -u $mcram |"
#if [[ $(stat -c %s $mcram 2>/dev/null || echo 0) -lt 100000 || $(stat -c %s $mcram.crai 2>/dev/null || echo 0) -lt 1200 ]]; then
#   echo '['$(date '+%m/%d %H:%M')"] merged cram file too small" | tee -a $rlog  
#   exit 1
#fi

if [[ $(get_fsize "$mcram") -lt 100000 || $(get_fsize "$mcram_crai") -lt 1200 ]]; then
  echo '['$(date '+%m/%d %H:%M')"] Warning: merged cram or cram index suspiciously small" | tee -a $rlog  
  exit 1
fi

## adding Mito_mapped metrics if not found
if ! fgrep -q Mito_mapped $fmet; then
  echo "adding Mito metrics" | tee -a $rlog
  mtcount=0
  mtcount=$(samtools view -T $gref -F 260 $mcram chrM | wc -l)
  if [[ -z "$mtcount" ]]; then 
     echo "Mito counts could not be retrieved. Check $rlog. " |& tee -a $rlog
     exit 1
  fi
  #samtools idxstats $fpri > $sid.pri.chrstats
  #mtcount=$(grep '^chrM\b' $sid.pri.chrstats | cut -f3)
  echo -e "Mito_mapped\t$mtcount" >> $fmet
fi

fstrand='strandness.tab'
## adding strandess info if not found
if [[ ! -s $fstrand ]]; then
 fgrep 'Sense Rate' $fmet | \
   perl -ne '($i,$v)=m/([12]) Sense Rate\s+([\d\.]+)$/;$e[$i]=$v; 
    END {$d=$e[1]-$e[2]; print abs($d)<0.2 ? "unstranded" : $d<0 ? "reverse" : "forward";
    print "\t$e[1]\t$e[2]\n"}' > $fstrand
fi

## finally, run bcftools to call the SNPs

fvcf="$sid.common_SNVs.vcf.gz"
if [[ ! -s $fvcf.csi ]]; then
## as seen in SPEAQeasy:
 bcf_pileup="-q 0 -Q 13 -d 1000000 -AB"
 bcf_call="-mvOz"
 cmd="bcftools mpileup -R $snvs $bcf_pileup --threads 3 -Ou -f $gref $mcram | bcftools call $bcf_call --threads 3 -o $fvcf - "
 echo -e ">> calling variants:\n$cmd" | tee -a $rlog
 eval "$cmd" |& tee -a $rlog 
 if [[ -s $fvcf ]]; then
    cmd="bcftools index --csi $fvcf"
    echo -e "  indexing vcf: $cmd" | tee -a $rlog
    eval $cmd |& tee -a $rlog
 else 
    echo "Error: $fvcf file not valid" |& tee -a $rlog
    exit 1
 fi
fi

echo 1 > "$sid.counts.done"
echo -e "["$(date '+%m/%d %H:%M')"]\tcounts task done [$taskid]." | tee -a $rlog

