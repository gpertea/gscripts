#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 gb2chrnames.pl GRCh_assembly_report.txt hgXX.fa.fai > gb2chr.txt
 
 Attempts to match the chromosome names from the fasta file
 with the ones from the NCBI assembly report, creating a translation
 table which can be used by gffread's -m option.
/;
umask 0002;

getopts('o:') || die($usage."\n");
die("$usage\n") unless @ARGV==2 && -f $ARGV[0] && -f $ARGV[1];

my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --

my %gb2chr; # GB accession => fasta name (multiple GB accs may point to the same chr name)
my %chr2gb; # fasta chr name => [GBacc1, GBacc2] (multiple GB accs possible)
my @fachr;
open(CHRNAMES, $ARGV[1]) || die("Error opening file $ARGV[1]\n");
open(GBTBL, $ARGV[0]) || die("Error opening assembly report file $ARGV[0]");
while(<GBTBL>) {
  next if m/^#/;
  my @t=split(/\t/);
  if ($t[1] eq 'assembled-molecule') {
    if ($t[0] eq 'MT') { #mouse file has Mt chromosome in there
      $t[0]='M';
    }
    $chr2gb{'chr'.$t[0]}=$t[6];
    next;
  }
  #anything else should be a scaffold:
  if ($t[1]=~m/\-scaffold/) {
    my $vname=$t[4];
    $vname=~tr/./v/;
    $chr2gb{$vname}=$t[6];
  }
}
close(GBTBL);

#now go through fasta names and match them in %chr2gb
while (<CHRNAMES>) {
  next if m/^#/;
  my @t=split(/\t/); # $t[0] = UCSC chr name
  if ($t[0]=~m/chr[\dXYMUn]{1,2}_([A-Z][A-Z]\d+v\d)/) {
    my $vacc=$1;
    my $refacc=$chr2gb{$vacc};
    if ($refacc) {
      print "$refacc\t$t[0]\n";
    }
    else {
      print STDERR "WARNING: no GenBank accession match found for $t[0] ($vacc)!\n";
    }
  } else { # assembled molecule
    my $refacc = $chr2gb{$t[0]};
    if ($refacc) {
      print "$refacc\t$t[0]\n";
    }
    else {
      print STDERR "WARNING: unrecognized chromosome/scaffold name: $t[0] !\n";
    }
  }
}
close(CHRNAMES);

#foreach my $acc (sort keys(%gb2chr)) {
# my $chr=$gb2chr{$acc};
# print "$acc\t$chr\n";
#}


# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

