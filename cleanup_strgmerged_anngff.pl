#!/usr/bin/perl
## this cleans up a file obtained by merging multiple stringtie --merge outputs
## which was also annotated with gffcompare
## Steps
##>1) two large lists of gtfs were merged with stringtie after filtering the trmap -J output
####    to discard all the class_code=r and the transcripts linking multiple genes
####  The list of gtfs was split in half to avoid the "too many files open" limit 
##$  stringtie-mflt --mflt=../dlpfc.n10_noR_noMrg.gtf --merge -G ~/ref/gencode.v32.annotation.gtf merge1.lst -o dlpfc_merge1.gtf &
##$  stringtie-mflt --mflt=../dlpfc.n10_noR_noMrg.gtf --merge -G ~/ref/gencode.v32.annotation.gtf merge2.lst -o dlpfc_merge2.gtf &

##>2) collapse the two merge results:
##$  gffread -MQ dlpfc_merge2.gtf dlpfc_merge1.gtf > dlpfc_merge.gff

##>3) annotate the merged result:
##$ gffcompare -r ~/ref/gencode.v32.annotation.gtf dlpfc_merge.gff
## -> this created gffcmp.annotated.gtf

##>4) convert to gff 
##$ gffread -F gffcmp.annotated.gtf -o dlpfc.merged.ann.gff
## ===> the input to this script should be dlpfc.merged.ann.gff
use strict;
while (<>) {
 if (s/;class_code==//) { 
   my ($t)=(m/ID=(EN[^;]+)/);
   s/;cmp_ref=\Q$t// if $t
 }
 if (s/;ref_gene_id=(EN[^;]+)/;gene_id=$1/) {
   my $g=$1;
   s/;geneID=[^;]+/;geneID=$g/;
 }
 print $_;
}
