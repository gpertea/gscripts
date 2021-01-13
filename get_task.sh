#!/bin/bash
#
ann=${HOME}/work/ref/gencode.v32.gff
fdb="$1"
if [[ -z "$fdb" ]]; then
  echo "Error: no cdb file given"
  exit 1
fi

if [[ ! -f "$fdb" ]]; then
  echo "Error: $fdb file not found!"
  exit 1
  fdb="$fdb.cidx"
  if [[ ! -f "$fdb" ]]; then
    echo "Error: $fdb index file not found!"
    exit 1
  fi
fi

id=$SGE_TASK_ID
if [[ -z "$id" ]]; then
 id="$2"
 if [[ -z "$id" ]]; then
  echo "Error: no task ID given!"
  exit 1
 fi 
fi

bams="" # original bam paths
bls="" # bam symlinks
multi=""
rnum=""
prot=""
while read -r line; do
 if [[ $line == '>'* ]]; then
   t=( $line )
   rnum=${t[1]}
   prot=${t[3]}
   mkdir -p "${t[2]}/$rnum"
   #echo "created ${t[2]}/$rnum (protocol $prot)"
   cd "${t[2]}/$rnum"
   continue
 fi
 t=( $line )
 bl="${t[2]}.bam"
 ln -s "${t[3]}" $bl
 if [[ -z $bls ]]; then
   bams=${t[3]}
   bls=$bl
 else
   bams+=" ${t[3]}"
   bls+=" $bl"
   multi=1
 fi
done < <(cdbyank -a $id $fdb)

ofn=$rnum
#cho "bls=<$bls>"
if [[ -f $ofn.info ]]; then
  echo "Error: $ofn.info already exists in $PWD!"
  exit 1
fi
echo "$rnum $prot $bams" > $ofn.info

if [[ -z "$multi" ]]; then
 cmd="stringtie -o $ofn.gtf -A $ofn.gabund.tab -G $ann $bls"
else
 cmd="samtools merge -O SAM - $bls | stringtie -o $ofn.gtf -A $ofn.ga -G $ann -"
fi

echo "Running: $cmd" > $ofn.log
echo "$cmd" >> $ofn.info
eval "$cmd" >> $ofn.log 2>&1
