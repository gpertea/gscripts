#!/usr/bin/perl
use strict;

my $usage = q/Usage:
 ann_fasta_from_gff.pl <transcripts.fa> <annotation.gff3>

 Outputs a fasta file based on the input <transcripts.fa> 
 but with protein product info added to the defline, as found
 in the provided annotation file.
 (By default attributes protSim, protSimAcc and protSimData are used)
/;

my ($fa, $gff)=@ARGV;
die($usage) unless $gff;
die($usage."Error: file $fa not found!\n") unless -f $fa;
die($usage."Error: file $gff not found!\n") unless -f $gff;

my %ann; #transcript_id => defline

open(GFF, $gff) || die("Error opening $gff ($!)\n");
open(FA, $fa) || die("Error opening $fa ($!)\n");

while(<GFF>) {
  next unless (m/protSim=([^;]+)/);
  my $product=$1;
  ##-- remove stray CDS=coords
  my ($tid)=(m/\bID=([^;]+)/);
  my ($acc)=(m/protSimAcc=([^;]+)/);
  my ($simdata)=(m/protSimData=([^;]+)/);
  my ($ev, $cov, $sim);
  if ($simdata) {
    ($ev)=($simdata=~m/e\-val:([^\,]+)/);
    ($cov)=($simdata=~m/cov:([^\,]+)/);
    ($sim)=($simdata=~m/sim:([^\,]+)/);
  }
  $ann{$tid}="similar to $acc $product (E-value=$ev, psim=$sim, coverage=$cov)";
}
close(GFF);


while (<FA>) {
 my $seqid;
 if (($seqid)=(m/^>(\S+)/)) {
   s/ CDS=\d+\-\d+//;
   my $pdescr=$ann{$seqid};
   if ($pdescr) {
     chomp;
     $_.=" $pdescr\n";
   }
 }
 print $_;
}
close(FA);
