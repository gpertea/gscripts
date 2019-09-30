#!/bin/bash

tmp='/ccb/salz4-4/gpertea/rice/tmp'
rd=rapdb
asmfa=cg_pseudochr.fa
reffa=$rd/IRGSP-1.0_genome.fa
#for c in 02 03 04 05 06 07 08 09 10 11 12; do
for c in {01..12}; do
 chr=chr$c
 echo "processing $chr .."
 acc=$(grep $chr refseq/chr2acc.txt | cut -f2 | cut -f1 -d '.')
 echo " acc=$acc"
 mdir=ann_$chr
 if [[ ! -d $mdir ]]; then
   mkdir $mdir
 fi
 if [[ ! -d $tmp/$chr ]]; then
   mkdir -p $tmp/$chr
 fi
 set -e
 cfa=cg.$chr.fa
 if [[ ! -s $mdir/$cfa ]]; then
   cdbyank -l $asmfa.cidx | fgrep "$acc" | cdbyank $asmfa.cidx > $mdir/$cfa
 fi
 if [[ ! -s $mdir/ref.$chr.fa ]]; then
   rm $mdir/ref.$chr.fa.fai
   cdbyank -a "$chr" $reffa.cidx > $mdir/ref.$chr.fa
   samtools faidx $mdir/ref.$chr.fa
 fi
 if [[ ! -s $mdir/gffr.$chr.prot.fa ]]; then
   grep '^'$chr'\b' $rd/IRGSP_1.0.40.chr.simplified.gff3 | gffread -g $mdir/ref.$chr.fa -y $mdir/gffr.$chr.prot.fa -w $mdir/gffr.$chr.rna.fa
   cdbfasta $mdir/gffr.$chr.prot.fa
   cdbfasta $mdir/gffr.$chr.rna.fa
 fi
 if [[ ! -s $mdir/ref.$chr.rna.gff3 ]]; then
   grep '^'$chr'\b' $rd/IRGSP_1.0.40.chr.simplified.gff3 > $mdir/ref.$chr.rna.gff3
 fi
 if [[ ! -s $mdir/ref.$chr.rna.fa ]]; then
   cdbyank -l $mdir/gffr.$chr.prot.fa.cidx | \
     cdbyank $rd/Oryza_sativa.IRGSP-1.0.pep.all.fa.cidx | seqmanip -l70 > $mdir/ref.$chr.pep.fa
   cdbyank -l $mdir/gffr.$chr.rna.fa.cidx | \
     cdbyank $rd/Oryza_sativa.IRGSP-1.0.cdna.all.fa.cidx | seqmanip -l70 > $mdir/ref.$chr.rna.fa
   cdbyank -l $mdir/gffr.$chr.rna.fa.cidx | \
     cdbyank $rd/Oryza_sativa.IRGSP-1.0.ncrna.fa.cidx | seqmanip -l70 >> $mdir/ref.$chr.rna.fa
 fi
 echo ".. $chr done."
done
