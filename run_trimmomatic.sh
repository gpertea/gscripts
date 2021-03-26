#!/bin/bash
prog="$HOME/sw/Trimmomatic-0.39/trimmomatic-0.39.jar"
adir="${prog%/*}/adapters"
if [[ ! -f $prog ]]; then
  echo "Error: path $prog is invalid!"
  exit 1
fi

read -r -d '' USAGE << EOM
Usage:
 run_trimmomatic.sh [-S] [-l NN] ..files to trim..
 
Options:
  -S      single end reads (using TruSeq3-SE.fa adapters)
          default: paired end reads (using TruSeq3-PE-2.fa)

  -l  NN  specify the minimum length for a read (default:75)
EOM
tmode="PE"
adfile="TruSeq3-PE-2.fa"
args="2:30:10:1:TRUE"
qargs="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75"
while getopts "hSl:" OPT; do
  case $OPT in
  h)
    echo "$USAGE"
    exit 1
    ;;
  S)
    tmode="SE"
    adfile="TruSeq3-SE.fa"
    args="2:30:10"
    ;;
  l)
    qargs="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$OPTARG"
    ;;
  *)
    echo "$USAGE"
    echo "Error: incorrect options provided"
    exit 1
    ;;
   esac
done

## on the cluster:
##  java -Xmx512M -jar ... PE -threads 4 -phred33 *.f*q* \
##  ${fpre}_trimmed_paired_1.fastq ${fpre}_unpaired_1.fastq ${fpre}_trimmed_paired_2.fastq ${fpre}_unpaired_2.fastq \
##  ILLUMINACLIP:$adapter_fa:2:30:10:1:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:45
##

shift $((OPTIND-1)) #shift $@
f1="$1"
f2="$2"
fouts=""
if [[ "$tmode" == "SE" ]]; then
  #single argument expected
  if [[ $# -ne 1 ]]; then
    echo "Error: a single file argument is expected!"
    exit 2
  fi
  if [[ ! -f $f1 ]]; then
     echo "Error: $f1 not found!"
     exit 2
  fi
  fouts="${f1%.f*q*}_trimmed.fastq.gz"
else # PE
  #two arguments expected
  if [[ $# -ne 2 ]]; then
    echo "Error: two arguments expected!"
    exit 2
  fi
  if [[ ! -f $f1 || ! -f $f2 ]]; then
     echo "Error: $f1 or $f2 not found!"
     exit 2
  fi
  f1o="${f1%.f*q*}_trimmed_paired_1.fastq.gz"
  f2o="${f2%.f*q*}_trimmed_paired_2.fastq.gz"
  f1un="${f2%.f*q*}_trimmed_unpaired1.fastq.gz"
  f2un="${f2%.f*q*}_trimmed_unpaired2.fastq.gz"
  fouts="$f1o $f1un $f2o $f2un"
fi

java -Xmx512M \
 -jar $prog $tmode \
 -threads 4 -quiet \
 -phred33 \
 -summary "${f1%.f*q*}.summary.txt" \
 $@ \
 $fouts \
 "ILLUMINACLIP:$adir/$adfile:$args" \
 $qargs
 
