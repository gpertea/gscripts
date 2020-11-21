#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
  gff2tab.pl [-o output] <miRge.gff> ...
/;
umask 0002;
getopts('o:') || die($usage."\n");
die($usage."\n") unless @ARGV>0;
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --            0          1             2              3        4      5         6       7
my @colnames=("miRNA", "database", "isoform_or_can", "start", "stop", "dot",  "strand", "dot2", 
  "Read", "UID", "Name", "Parent", "Variant", "Cigar", "Filter");
my @rnums; # Rxxxx in order they are parsed, $rnums[0]=first file (R1) etc.
my %mirs; # mirID => [ metadata, [ count_in_R1, count_in_R2,  ... ] ]
foreach my $fn (@ARGV) {
  open(F, $fn) || die("Error opening file $fn\n");
  my $rnum;
  while (<F>) {
    if (m/^##/) {
       $rnum=$1 if (m/COLDATA:\s+(R\d+)/);
       next;
    }
    chomp;
    my @t=split(/\t/);
    map { push(@t, (split(/ /, $_))[1]) } (split(/; /,pop(@t)));
    my $c=splice(@t, 14,1); #count
    my $md=join("\t",@t);
    my $d=$mirs{$t[0]};
    if ($d) { #validate metadata
       die("Error: metadata mismatch for file $fn:\n$md\n$$d[0]\n") unless $md eq $$d[0];
    }
    else {
       $mirs{$t[0]}=[$md, []];
       $d=$mirs{$t[0]};
    }
    $$d[1][scalar(@rnums)]=$c;
  }
  die("Error: $rnum not parsed!\n") unless $rnum;
  push(@rnums, $rnum);
  close(F);
}

print join("\t", @colnames, @rnums)."\n"; # print header
foreach my $k (keys %mirs) {
  my $d=$mirs{$k};
  print join("\t", $k, $$d[0],
     ( map { int($_) } @{$$d[1]} ) )."\n";
}

# --
if ($outfile) {
  select(STDOUT);
  close(OUTF);
}

#************ Subroutines **************

