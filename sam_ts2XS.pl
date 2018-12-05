#!/usr/bin/perl
use strict;
my $usage = q/Usage:
 sam_ts2XS.pl minimap2.sam > converted.sam
 
 Converts the "ts" tags from spliced alignments produced by minimap2
 into "XS" tags for Cufflinks\/StringTie.
 It also discards unmapped reads, if any.
/;

while (<>) {
  if (m/^@(HD|SQ)\t/) {
    print $_;
    next;
  }
  my $samline=$_;
  chomp;
  my ($qname, $flags, $gseq, $pos, $mapq, $cigar, $rnext, $pnext, 
      $tlen, $seq, $qual, @extra)= split(/\t/);
  next if ($flags & 0x04)!=0;
  my $strand= (($flags & 0x10)==0) ? '+' : '-';
  #my @cigdata=($cigar=~m/(\d+[A-Z,=])/g);
  my $changed;
  foreach my $tag (@extra) {
    if ($tag=~m/ts:A:([\+\-])/) {
      my $xs=$1;
      if (($flags & 0x10)!=0) {
        $xs = ($xs eq '-') ? '+' : '-';
      }
      $tag = "XS:A:$xs";
      $changed=1;
      last;
    }
  }
  print $changed ? join("\t", $qname, $flags, $gseq, $pos, $mapq, 
     $cigar, $rnext, $pnext, $tlen, $seq, $qual, @extra)."\n" : $samline;
}
#--
