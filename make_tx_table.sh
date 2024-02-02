#!/bin/bash
fgff=$1
if [[ -z $fgff ]]; then
 echo "Error: input GTF/GFF required!"
 exit 1
fi
if [[ ! -f $fgff ]]; then
  echo "$fgff file not found!\n"
  exit 1
fi
tabcols='@id,@chr,@start,@end,@strand,@numexons,@covlen,@cdslen,@track,@geneid,gene_name,gene_type,transcript_name,transcript_type,level,tag'
ftab=${fgff/.gtf/.tab}
if [[ $ftab == $fgff ]]; then
 ftab=${fgff/.gff/.tab}
fi
echo "$tabcols" |  tr -d '@'  | tr "," "\t" > $ftab
#zcat $fgff | gffread --table $tabcols >> $ftab
cat $fgff | gffread --table $tabcols >> $ftab
