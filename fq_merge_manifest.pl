#!/usr/bin/perl
use strict;
use Getopt::Std;
use strict;
my $usage = q/Usage:
  fq_merge_manifest.pl samples.manifest > samples.manifest.mrg
 Merges the lines belonging to the same sample ID (column 5) from the given
 manifest file, creating comma delimited lists for read 1 (column 1) and
 read 2 (column 3)
/;
umask 0002;
#getopts('o:') || die($usage."\n");
#my $outfile=$Getopt::Std::opt_o;
my @smp; #list of unique sample IDs in the order they are encountered
my %smpfq1; # sampleID => [ read1_fastqs... ]
my %smpfq2; # sampleID => [ read2_fastqs... ]
while (<>) {
 chomp;
 my @t=split(/\s+/);
 die("No sampleID parsed for line:\n$_\n") unless $t[4];
 my ($f1, $f2, $s)=@t[0,2,4];
 if (!exists($smpfq1{$s})) {
   push(@smp,$s);
   $smpfq1{$s}=[$f1];
   $smpfq2{$s}= $f2 ? [ $f2 ] : [];
 } else {
   push(@{$smpfq1{$s}}, $f1);
   push(@{$smpfq2{$s}}, $f2) if $f2;
 }
}

foreach my $s (@smp) {
  my $f1=$smpfq1{$s} || die("Error: no read 1 fastq files for sample $s\n");
  my $f2=$smpfq2{$s};
  my $lst2 = ($f2 && @$f2>0) ? join(',', @$f2) : '.' ;
  print join("\t", join(',',@$f1), 0, $lst2, 0, $s)."\n";
}

