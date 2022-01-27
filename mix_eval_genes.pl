#!/usr/bin/perl
## evaluates accuracy (matched reference transcripts) for each ref gene name
## for 3 assembly results: --mix, long only, short only

## to spot entries where mixed did much better in retrieving transcript structures:
##  gawk '$2<$4 && $3<$4' mix_eval_table.tab
## -- or even better, when ONLY mix retrieved any tx in a locus:
##  gawk '$2==0 && $3==0' mix_eval_table.tab

use strict;

my $usage=q{
 Usage:
  mix_eval_genes.pl short.ann.tbl long.ann.tbl mix.ann.tbl
 Generates a table with number of matching assembled transcripts and total assembled
 transcripts for each gene locus and each of the short/long/mix assembly results
 
 The ann.gtf files are generated by first running gffcompare with -MNT option on each of the
 stringtie gtf outputs (short/long/mix):

 gffcompare -MNT -r ref_ann.gff strg_(short/long/mix).gtf -o cmpMN_(short/long/mix)
 
 ..then running gffread on the annotated.gtf files like this:
 gffread --table transcript_id,gene_id,@chr,@start,@end,@strand,@numexons,gene_name,cmp_ref,class_code \
        cmpMN_(short/long/mix).annotated.gtf > (short/long/mix).ann.tbl
};
die $usage unless @ARGV==3;


my (%mgh, %lgh, %sgh, %allg);
# *gh: gene => [tref_match_count, tcount, max_exons_per_tx]
# %allg : gene => tref_match_count (across ALL tbl)

my (%mg2strg, %sg2strg, %lg2strg);# gene => [STRG.N, .. ] #mapping gene names to STRG loci

## first run gffcompare on each stringtie gtf output (--mix , -L, default w/ short reads only)
# gffcompare -MN -r chr21.gff -T chr21_{short,mix,long}_wguides.gtf -o cmpMN_{short,mix,long}_wguides
### this script operates on tables generated by this command:
###                     0            1     2     3     4     5        6         7        8       9        10
## gffread --table transcript_id,gene_id,@chr,@start,@end,@strand,@numexons,gene_name,cmp_ref,class_code,ref_gene_id \
##         cmpMN_{short,mix,long}_wguides.annotated.gtf 
my $exc=4; #minimum number of exons to consider for matching transcripts in a gene
loadGenes($ARGV[0], \%sgh, \%sg2strg);
loadGenes($ARGV[1], \%lgh, \%lg2strg);
loadGenes($ARGV[2], \%mgh, \%mg2strg);
#             gawk:   1     2         3       4    5       6        7        8          9        10       11          12         13
print join("\t", qw(gene short_nm long_nm mix_nm short_nt long_nt mix_nt short_nloc long_nloc mix_nloc short_numex long_numex mix_numex))."\n";
# for each gene and assembly type it reports:
#    * the number of matching transcripts (*_nm)
#    * the number of total transcripts (*_nt)
while (my($g, $c) = each(%allg)) { 
  my $d;
  my ($nsm, $nst, $nsx) = ($d=$sgh{$g}) ? @$d : (0,0,0); # (matching, total, max_exons)
  my ($nlm, $nlt, $nlx) = ($d=$lgh{$g}) ? @$d : (0,0,0);
  my ($nmm, $nmt, $nmx) = ($d=$mgh{$g}) ? @$d : (0,0,0); 
  next unless $nmx>=$exc; #only keep results with at least $exc exons in mix
  my $nmstrg = ($d=$mg2strg{$g}) ? scalar(@$d) : 0;
  my $nsstrg = ($d=$sg2strg{$g}) ? scalar(@$d) : 0;
  my $nlstrg = ($d=$lg2strg{$g}) ? scalar(@$d) : 0;
  print join("\t", $g, $nsm, $nlm, $nmm, $nst, $nlt, $nmt, $nsstrg, $nlstrg, $nmstrg, $nsx, $nlx, $nmx)."\n";
}

###---------------
sub loadGenes {
  my ($fn, $href, $hs)=@_;
  open(FH, $fn) || die("Error opening $fn\n");
  while (<FH>) {
    chomp;
    my @t=split(/\t/); 
    # 0:tid   1:gid(STRG.N)  2:chr  3:start  4:end  5:strand  6:numexons  7:gname  8:tref  9:class 10:ref_gene_id
    next if ($t[7] eq '.');
    my $gref=$t[10].'|'.$t[7];
    my $d=$$href{$gref};
    unless($d) {
      $d=[0,0, 0]; # num_tx_matching, num_tx, max_exon_count_for_matching
      $$href{$gref}=$d;
    }
    $$d[1]++;
    if ($t[9] eq '=') {
       $$d[0]++;
       $allg{$gref}++;
       $$d[2]=$t[6] if ($$d[2] < $t[6]); #update max exon count for matching tx
    }
    my $sd=$$hs{$gref};
    unless($sd) {
     $sd=[];
     $$hs{$gref}=$sd; # list of distinct STRG.N (subloci) for this gene (STRG.N)
    }
    # only add STRG.N if not already there
    my @si=grep { $$sd[$_] eq $t[1] } (0 .. $#$sd);
    push(@$sd, $t[1]) if (@si==0); 
    
  }
  close(FH)
}

