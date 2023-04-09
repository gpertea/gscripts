#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q{Usage:
  update_gene_by_besTmatch.pl trmap_TB.tab file.gtf/gff

  Change gene_id and gene_name to match those of best overlapping
  entries in the given trmap -TB output file.
  
  To be used to bring refseq added entried into Gencode annotation and gene namespace
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
die("${usage} Error: input file expected!\n") unless $tmap && $tmap;

my %tr; ## tid=>[ gene_name, gene_id ]

open(T, $tmap) || die("Error opening file $tmap\n");
while (<T>) {
  my @t=split(/\t/);
  my ($t, $tg)=split(/\|/, $t[0]);
  my $c=$t[1];
  my ($r, $rg, $rgid)=split(/\|/, $t[2]);
  if ($rg eq $tg || $c=~m/[joec=kmn]/) {
   $tr{$t}=[$rg, $rgid];
  }
}
close(T);
my $gtf=0;

# read/transform gff/gtf
while (<>) {
  if (m/^#/) { print $_;next }
  chomp;
  if (m/ID=([^;]+)/) {
  # GFF
    my $t=$1;
    if (my $td=$tr{$t}) {
      my ($gid, $gn)=@$td;
      unless( s/gene=[^;]+/gene=$gn/ || 
              s/gene_name=[^;]/gene_name=$gn/) {
         s/;?$/;gene=$gn/;
      }
      unless( s/Parent=([^;]+)/Parent=$gid;oParent=$1/ || 
            s/geneID=[^;]+/geneID=$gid/) {
         s/;?$/;geneID=$gid/;
      }
    }
  } elsif (m/transcript_id\s+"([^"]+)/) { 
  # GTF
    my $t=$1;
    if (my $td=$tr{$t}) {
      my ($gid, $gn)=@$td;
      unless( s/gene_id\s+"[^"]+/gene_id "$gid/ ) {
         s/;?$/;gene_id "$gid"/;
      }
      unless( s/gene_name\s+"[^"]+/gene_name "$gn/ ) {
         s/;?$/;gene_name "$gn"/;
      }
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

