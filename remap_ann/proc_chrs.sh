#!/bin/bash
#if [[ -z $1 ]]; then
# echo -e "Usage:\n $0 nn[ nn ...]"
# echo "Where nn are 2-digit chromosome numbers, e.g. 01 02 03 etc."
# exit 1
#fi
rd=rapdb

#for c in $@; do
for c in {01..12}; do
 chr=chr$c
 #mdir=ann_$chr
 echo ">$chr begin processing."
 cd ann_$chr
  #--------------
  #../gmap_run.sh $chr &
  #minimap2 -cx asm5 -t4 --cs ref.$chr.fa cg.$chr.fa > cg2ref_asm5.paf
  #minimap2 -cx asm10 -t4 --cs ref.$chr.fa cg.$chr.fa > cg2ref_asm10.paf
  #
  #wait
  #--------------
  #gff_liftover.pl -C ref2cg_asm5.paf ref_${chr}.gffread.gff3 > cg.asm5.liftover.gff3
  if [[ ! -f ref.${chr}.pep.fa.cidx ]]; then
    cdbfasta ref.${chr}.pep.fa
  fi
  if [[ ! -f ref.$chr.fa.fai ]]; then
    samtools faidx ref.$chr.fa
  fi
  if [[ ! -f cg.$chr.fa.fai ]]; then
    samtools faidx cg.$chr.fa
  fi
  set -e
  grep '^'$chr'\b' ../$rd/IRGSP_1.0.40.chr.simplified.gff3 | gffread -F -o ref_${chr}.gffread.gff3
  grep '^'$chr'\b' ../$rd/IRGSP_1.0.40.chr.simplified.gff3 | \
   gffread -g ref.$chr.fa -F -P --tlf -o ref.ann.tlf
  cdbyank -l ref.${chr}.pep.fa.cidx | sort -u > ref.rna.protein_coding.lst
  gff_liftover.pl -C ref2cg_asm10.paf ref_${chr}.gffread.gff3 > cg.asm10.liftover.gff3
  paftools.js splice2bed minimap2.rna2cg_${chr}.paf > minimap2.rna2cg.paf.bed
  #take care of duplicates (multiple mappings)
  perl -ne '@t=split(/\t/);$v=(++$h{$t[3]});$t[3].=".mrna$v";print join("\t",@t)' \
   minimap2.rna2cg.paf.bed > minimap2.rna2cg.dedup.bed
  gffread -o minimap2.rna2cg.paf.gff minimap2.rna2cg.dedup.bed
  gffcount -l minimap2.rna2cg.paf.gff | cut -f1 -d '.' | sort -u > minimap2_rna2cg.lst
  gffcount -l gmap.rna2cg_$chr.gff | cut -f1 -d '.' | sort -u > gmap_rna2cg.lst
  
  #gffcount -l cg.asm10.liftover.gff3 | sort -u > liftover_cds2cg.lst
  comm -23 minimap2_rna2cg.lst gmap_rna2cg.lst > gmap_missed.lst
  #comm -23 liftover_cds2cg.lst gmap_rna2cg.lst >> gmap_missed.lst
  grep '^'$chr'\b' ../$rd/IRGSP_1.0.40.chr.simplified.gff3 | \
  grep -E '(description|biotype)=' | grep -v -P 'gene\t' | perl -ne \
  'chomp;($id)=(m/ID=([^;]+)/);($b)=(m/biotype=([^;]+)/);
  ($p)=(m/Parent=([^;]+)/);print join("\t",$id, $b, $p)."\n"' \
  | sort > ref_transcripts.tab
  perl -ne 'chomp;($id)=(m/ID=([^;]+)/);($b)=(m/biotype=([^;]+)/);
  ($g)=(m/geneID=([^;]+)/);($d)=(m/description=([^;]+)/);
  print join("\t",$id, $b, $g, $d)."\n"' < ref.ann.tlf > ref_transcripts_descr.tab
  #
  gffread --tlf -o rna2cg.mappings.tlf gmap.rna2cg_$chr.gff
  fgrep -w -f gmap_missed.lst minimap2.rna2cg.paf.gff | gffread --tlf -o- >> rna2cg.mappings.tlf
  fgrep -w -f ref.rna.protein_coding.lst rna2cg.mappings.tlf > rna2cg.protein_coding.tlf
  gffread -g cg.$chr.fa -W -w rna2cg.protein_coding.fa rna2cg.protein_coding.tlf
  cds2rna_blat.pl rna2cg.protein_coding.fa ref.${chr}.pep.fa.cidx > rna2cg.CDS_reassign.tlf
  gffread -V -F -P -W -y rna2cg.CDS_reassign.pep.fa -g cg.$chr.fa rna2cg.CDS_reassign.tlf \
   -o rna2cg.coding.validated.gff
   
  #filter ncrna mappings
  gffread -F --tlf -o gmap.mappings.tlf gmap.rna2cg_$chr.gff
  fgrep -v 'protein_coding' ref_transcripts_descr.tab | cut -f1 > ref.rna.noncoding.lst
  fgrep -w -f ref.rna.noncoding.lst gmap.mappings.tlf | perl -ne \
  'if (s/\.mrna1;/;/) { s/(geneID|CDS|Name|gene_name|CDSphase)=[^;]+;?//g; print $_ }' \
   > gmap.noncoding.top.tlf
  ##--this should be further filtered by coverage & identity to avoid bad mappings!
  perl -ne '($c)=(m/coverage=([\d\.]+)/);
  ($i)=(m/identity=([\d\.]+)/);print $_ if $i>90 && $c>90' gmap.noncoding.top.tlf > gmap.noncoding.filtered.tlf
  ##fgrep -w -f ref.rna.noncoding.lst gmap_missed.lst | fgrep -w -f - minimap2.rna2cg.paf.gff
  ## careful above, there might be 0 noncoding transcripts missed by gmap!
  ##
  ## noncoding transcripts shall be gathered in: rna2cg.noncoding.gff
  ## then we should transfer gene and descr from ref_transcripts_descr.tab
  ## to both rna2cg.coding.validated.gff and rna2cg.noncoding.gff
  echo -e "\t$chr done."
 cd ..
done
echo "All Done."
