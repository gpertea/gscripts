#!/bin/bash
for r in R*; do
 cd $r
samtools view htt_alns.bam | fgrep 'AS:i:101' | cut -f1,3,10 | grep -P '(CAG){17,}' |\
 grep -P '[\tT]CCAGCAG' | fgrep 'CAGCAA' | perl -nle \
 'BEGIN{use File::Basename} print basename($ENV{PWD}),"\t",(length($1)/3) if m/((CAG){13,})/' \
 | sort | uniq -c | perl -nle 'chomp;@t=split;print join("\t", @t[1,2], $t[0])' > CAG_counts.tab
nmax=17
if [[ -s CAG_counts.tab ]]; then
  nmax=$(sort -k2,2n CAG_counts.tab | tail -1 | cut -f2)
fi
samtools view htt_alns.bam | fgrep 'AS:i:101' | cut -f1,3,10 | grep -P "(CAG){$nmax,}" |\
 perl -nle 'BEGIN{use File::Basename} print basename($ENV{PWD}),"\t",(length($1)/3) if m/((CAG){13,})/'|\
 sort | uniq -c | perl -nle 'chomp;@t=split;print join("\t", @t[1,2], $t[0])' > CAG_counts_unbounded.tab
 cd ..
done

#samtools view htt_alns.bam | fgrep 'AS:i:10' | cut -f1,3,10 | grep -P '(CAG){17,}' | grep -P '[\tT]CCAGCAG' | fgrep 'CAGCAA' |  perl -nle 'print ((length($1)/3),"\t",$_) if m/((CAG){13,})/'
#samtools view htt_alns.bam | fgrep 'AS:i:101' | cut -f1,3,10 | grep -P '(CAG){27,}' | perl -nle 'print ((length($1)/3),"\t",$_) if m/((CAG){27,})/'
