#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 hh [-n<numlines> [-c<numcols>] [-q] file.tab
 Options:
  -n : number of lines to show
  -c : number of columns to show
  -q : do not show the info line
/;
umask 0002;
getopts('qto:n:c:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $nl=$Getopt::Std::opt_n || 6;
my $nc=$Getopt::Std::opt_c || 6;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
# --
$nc--;
my $tc;
while(<>) {
 chomp;
 my @t=split(/\t/);
 unless ($tc) {
  $tc=scalar(@t);
  print "(found $tc columns.)\n" unless $Getopt::Std::opt_q;
 }
 print join("\t", @t[0..$nc])."\n";
 last if $.==$nl;
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

