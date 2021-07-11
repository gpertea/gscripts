#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 hh [-n<numlines> [-c<numcols>] [-q] [-t] file.tab
 Options:
  -t : <numlines> are shown from the end of the file
  -q : do not show the info line at the end
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

