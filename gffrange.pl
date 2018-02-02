#!/usr/bin/perl
use strict;
use Getopt::Std;
my $usage=q{
 A basic GFF/GTF filter by genomic location.

 Usage:
   ... | gffrange [-C] seqname:start-end[strand]
   
 The program expects a GFF or GTF formatted stream at stdin
 and outputs only the features found to overlap the provided interval.
 The strand character (+ or -) is optional.
 
 Example:
   zcat genomic.gff.gz | gffrange.pl chr8:21736219-21749705+ > c8r.gff
 
 Use the -C option to restrict the output to only features fully "contained"
 in the given region (excluding those just partially intersecting the region).
 };

getopts('C') || die "$usage\n";
my $strict=$Getopt::Std::opt_C;
my $loc=shift(@ARGV);
die "$usage\n" unless $loc;
my $strand;
if ($loc!~m/\d$/) {
  $strand=chop($loc);
  die("$usage\nError: strand must be '-' or '+', nothing else)\n") 
     if ($strand ne '-' && $strand ne '+')
}
my ($chr, $start, $end)=($loc=~m/(.+)\:(\d+)[\-\.]+(\d+)$/);
die("$usage\nError: cannot parse the provided range ($chr:$start-$end)\n") 
     unless $chr && $start>0 && $end>0;
($start, $end)=($end, $start) if $end<$start;
while (<>) {
 next if m/^#/;
 my @t=split(/\t/);
 next unless @t>8;
 next if $strand && $strand ne $t[6];
 next if $t[0] ne $chr;
 if ($strict) {
   next unless ($t[3]>=$start && $t[4]<=$end);
 }
 else {
   next unless ($t[3]<=$end && $start <= $t[4]);
 }
 print $_;
}
