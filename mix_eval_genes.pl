#!/usr/bin/perl
## evaluates accuracy (matched reference transcripts) for each ref gene name
## for 3 assembly results: --mix, long only, short only

## problem genes: SLC37A1
## spot them with: gawk '$2<$4 || $2<$3' mix_eval_table.tab
use strict;
my (%mgh, %lgh, %sgh, %allg);
# *gh: gene => [tref_match_count, tcount]
# %allg : gene => tref_match_count (across ALL tbl)

my (%mg2strg, %sg2strg, %lg2strg);# gene => [STRG.N, .. ] #mapping gene names to STRG loci

## first run gffcompare on each stringtie gtf output (--mix , -L, default w/ short reads only)
# gffcompare -MN -r chr21.gff -T chr21_{short,mix,long}_wguides.gtf -o cmpMN_{short,mix,long}_wguides
### this script operates on tables generated by this command:
###                     0            1     2     3     4     5        6         7        8       9
## gffread --table transcript_id,gene_id,@chr,@start,@end,@strand,@numexons,gene_name,cmp_ref,class_code \
##         cmpMN_{short,mix,long}_wguides.annotated.gtf 
loadGenes('cmpMN_mix_wG.tbl', \%mgh, \%mg2strg);
loadGenes('cmpMN_long_wG.tbl', \%lgh, \%lg2strg);
loadGenes('cmpMN_short_wG.tbl', \%sgh, \%sg2strg);
#TODO: we should really gather STRG.N loci overlapping each gene(s)
#TODO:    and use a gene_name:loc gene identifier, gathering/merging locus coordinates etc.

print join("\t", qw(gene mix_nm short_nm long_nm mix_nt short_nt long_nt mix_nloc short_nloc long_nloc))."\n";
while (my($g, $c) = each(%allg)) { 
  my $d;
  my ($nmm, $nmt) = ($d=$mgh{$g}) ? @$d : (0,0); # (matching, total)
  my ($nsm, $nst) = ($d=$sgh{$g}) ? @$d : (0,0);
  my ($nlm, $nlt) = ($d=$lgh{$g}) ? @$d : (0,0);
  my $nmstrg = ($d=$mg2strg{$g}) ? scalar(@$d) : 0;
  my $nsstrg = ($d=$sg2strg{$g}) ? scalar(@$d) : 0;
  my $nlstrg = ($d=$lg2strg{$g}) ? scalar(@$d) : 0;
  print join("\t", $g, $nmm, $nsm, $nlm, $nmt, $nst, $nlt, $nmstrg, $nsstrg, $nlstrg)."\n";
}

###---------------
sub loadGenes {
  my ($fn, $href, $hs)=@_;
  open(FH, $fn) || die("Error opening $fn\n");
  while (<FH>) {
    chomp;
    my @t=split(/\t/); 
    # 0:tid   1:gid  2:chr  3:start  4:end  5:strand  6:numexons  7:gref  8:tref  9:class
    next if ($t[7] eq '.');
    
    my $d=$$href{$t[7]};
    unless($d) {
      $d=[0,0];
      $$href{$t[7]}=$d;
    }
    $$d[1]++;
    if ($t[9] eq '=') {
       $$d[0]++;
       $allg{$t[7]}++;
    }
    my $sd=$$hs{$t[7]};
    unless($sd) {
     $sd=[];
     $$hs{$t[7]}=$sd;
    }
    # only add STRG.N if not already there
    my @si=grep { $$sd[$_] eq $t[1] } (0 .. $#$sd);
    push(@$sd, $t[1]) if (@si==0); 
    
  }
  close(FH)
}

