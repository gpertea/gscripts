#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q{Usage:
  update_gene_by_besTmatch.pl trmap_TB.ttab trmap_J.jtab tofix.gtf/gff

  Change gene_id and gene_name to match those of best overlapping
  entries in the given trmap -TB and -J output files
  
  To be used to bring refseq added entried into Gencode annotation and gene namespace
  or to just assign a gene_id/gene_name to an arbitrary transfrag file
};
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
die($usage."\n") if @ARGV==0 || $ARGV[0]=~m/^\-+h/;
my $tmap=shift(@ARGV);
die("${usage} Error: input tmap file expected!\n") unless $tmap && -f $tmap;
my $jtab=shift(@ARGV);
die("${usage} Error: input jtab file expected!\n") unless $jtab && -f $jtab;

my %tr; ## tid=>[ gene_name, gene_id ]
my %tj; ## tid=> novel_junction_info
open(T, $tmap) || die("Error opening file $tmap\n");
while (<T>) {
  my @t=split(/\t/);
  my ($t, $tg)=split(/\|/, $t[0]);
  my $c=$t[1];
  my ($r, $rg, $rgid)=split(/\|/, $t[2]);
  if ($rg eq $tg || $c=~m/[joec=qkmn]/) {
   $tr{$t}=[$rg, $rgid];
  }
}
close(T);
open(J, $jtab) || die("Error opening file $jtab\n");
while (<J>) {        #  0    1        2                     3          4
  chomp;
  my @t=split(/\t/); # tid exons codes|refTx|ref_gnames, ovl_genes,  njx_info,
  $tj{$t[0]} = ($t[4] eq '.') ? '' : $t[4];
}
close(J);

my $gtf=0;

# read/transform gff/gtf
while (<>) {
  if (m/^#/) { print $_;next }
  chomp;
  if (m/\bID=([^;]+)/) {
  # GFF
    my $t=$1;
    if (my $td=$tr{$t}) {
      my ($gid, $gn)=@$td;
      my $repl=0;
      $repl |= s/gene=([^;]+)/gene=$gn;ogene=$1/;
      $repl |= s/gene_name=([^;])/gene_name=$gn;ogene_name=$1/;
      s/;?$/;gene=$gn/ unless $repl;
      $repl=0;
      $repl |= s/Parent=([^;]+)/Parent=$gid;oParent=$1/;
      $repl |= s/geneID=([^;]+)/geneID=$gid;ogeneID=$1/; 
      s/;?$/;geneID=$gid/ unless $repl;
    }
    if (my $jd=$tj{$t}) {
      s/;?$/;njx=$jd/;
    }
  } elsif (m/transcript_id\s+"([^"]+)/) { 
  # GTF
    my $t=$1;
    if (my $td=$tr{$t}) {
      my ($gid, $gn)=@$td;      
      unless (s/gene_id\s+"([^"]+)/gene_id "$gid"; ogene_id "$1"/) {
         s/;?$/;gene_id "$gid"/;
      }
      unless( s/gene_name\s+"([^"]+)/gene_name "$gn; ogene_name "$1"/ ) {
         s/;?$/;gene_name "$gn"/;
      }
    }
    if (m/\ttranscript\t/ && (my $jd=$tj{$t})) {
      s/;?$/;njx "$jd"/;
    }
  }

  print "$_\n";
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

