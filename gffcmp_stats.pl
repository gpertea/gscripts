#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 gffcmp_stats.pl <gffcmp.stats>
 
 Collect the relevant stats from the gffcompare *.stats output
 in a tab-delimited table:
 qfname rtnum rlocnum qtnum qlocnum matchICnum matchTnum iSn iPr icSn icPr tSn tPr
 Use -H option to print the header like the above
/;
umask 0002;
getopts('Ho:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
if ($Getopt::Std::opt_H) {
 print join("\t", 
  split(/\s+/, q/queryFileName rtnum rlocnum qtnum qlocnum matchIC matchT iSn iPr icSn icPr tSn tPr/))."\n";
}
my ($qf, $rtnum, $rlocnum, $qtnum, $qlocnum, $icMatch, $tMatch, $iSn, $iPr, $icSn, $icPr, 
  $tSn, $tPr);
while (<>) {
  if (m/^#= Summary for dataset:\s+(\S+)/) {
    my $f=$1;
    if ($qf) {
      print join("\t", ($qf, $rtnum, $rlocnum, $qtnum, $qlocnum, $icMatch, $tMatch, $iSn, $iPr, 
                       $icSn, $icPr, $tSn, $tPr))."\n";
    }
    $tMatch='';
    $qf=$f;
    $qf=~s/\.(\w+)$//;
    next;
  }
  if (m/^#\s+Query mRNAs\s*:\s*(\d+)\s+in\s+(\d+)\s+loci/) {
    ($qtnum, $qlocnum)=($1, $2);
    next;
  }
  if (m/^#\s+Reference mRNAs\s*:\s*(\d+)\s+in\s+(\d+)\s+loci/) {
    ($rtnum, $rlocnum)=($1, $2);
    next;
  }
  if (m/^\s+Matching intron chains\s*:\s*(\d+)/) {
    $icMatch=$1;
    next;
  }
  if (m/^\s+Matching transcripts\s*:\s*(\d+)/) {
    $tMatch=$1;
    next;
  }
  if (m/^\s*Intron level\s*:\s*([\d\.]+)[\s\|]+([\d\.]+)/) {
    ($iSn, $iPr)=($1, $2);
    next;
  }
  if (m/^\s*Intron chain level\s*:\s*([\d\.]+)[\s\|]+([\d\.]+)/) {
    ($icSn, $icPr)=($1, $2);
    next;
  }
  if (m/^\s*Transcript level\s*:\s*([\d\.]+)[\s\|]+([\d\.]+)/) {
    ($tSn, $tPr)=($1, $2);
    next;
  }

} #while <>

if ($qf && $tMatch) {
  print join("\t", ($qf, $rtnum, $rlocnum, $qtnum, $qlocnum, $icMatch, $tMatch, $iSn, $iPr, $icSn, $icPr, 
     $tSn, $tPr))."\n";
}


# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

