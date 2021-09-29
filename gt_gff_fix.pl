#!/usr/bin/perl
## restoring the IDs from genometools' gt gtf_to_gff3 conversion
##
use strict;
while (<>) {
  if (m/^#/) { print $_; next }
  chomp;
  my ($chr, $track, $ftype, $start, 
      $end, $sc, $strand, $phase, $attrs)=split(/\t/);
  if ($ftype eq 'gene') {
    s/;gene_id=([^;]+)//;
    my $id=$1;
    s/ID=([^;]+)/ID=$id/;
    print "$_\n";
    next;
  }
  if ($ftype eq 'mRNA') {
    s/;transcript_id=([^;]+)//;
    my $id=$1;
    s/;gene_id=([^;]+)//;
    my $p=$1;
    s/\bID=([^;]+)/ID=$id/;
    s/\bParent=([^;]+)/Parent=$p/;
    print "$_\n";
    next;
  }
  ## only exon and CDS from here on
  s/;transcript_id=([^;]+)//;
  my $p=$1;
  s/;gene_id=([^;]+)//;
  #my $p=$1;
  s/\bParent=([^;]+)/Parent=$p/;
  print "$_\n";
}
