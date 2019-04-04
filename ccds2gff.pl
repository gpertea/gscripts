#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:

  ccds2gff.pl [-o output.gff][-B][-L] CCDS.txt CCDS2Sequence.txt

  Options:
    -B :   output BED format instead of GFF3
    -L :   output TLF (Transcript Line Format) instead of GFF3
/;
umask 0002;
getopts('BLo:') || die($usage."\n");
my $bed=$Getopt::Std::opt_B; #bed output instead of gff
my $tlf=$Getopt::Std::opt_L; #tlf output 
die($usage."\n") unless -f $ARGV[0] && -f $ARGV[1];
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
# --
my %c2s; # ccds_id => [ [ncbi_src,...] , [ebi_src] ] #only current members
open(CSEQ, $ARGV[1]) || die("Error opening file $ARGV[1]\n");
while (<CSEQ>) {
  #ccds, original_member, current_member, source, nucleotide_ID, protein_ID, status_in_CCDS, sequence_status
  next if m/^#/;
  chomp;
  my @t=split(/\t/);
  next if $t[2] eq '0'; #exclude the Updated or Removed 
                        #  entries which are no longer current
  my $d=$c2s{$t[0]};
  if (!$d) {
    $d=[ [] , [] ];
    $c2s{$t[0]}=$d;
  }
  my $si=($t[3]=~m/NCBI/)? 0 : 1;
  push(@{$d->[$si]}, $t[4]);
}
close(CSEQ);
my $numkeys=keys(%c2s);
print STDERR "Loaded $numkeys CCDSs from source data file $ARGV[1]\n";
#                                0      1      2     3    4      5
#my %cdata;  # chr|CCDS#.# => [ chr, strand, start, end, segs, gene ]
            #  where segs = [ [start, end], ... ]
#my @ccds; # list of chr|CCDS#.#
#                           0      1     2      3     4      5      6
#passing to writeCCDS: [ccds_id, chr, strand, start, end, [segs], gene ]
open(CCDS, $ARGV[0]) || die("Error opening file $ARGV[0]\n");
my $count=0;
while (<CCDS>) {
  #chr, chr_acc, gene, gene_id#, ccds_id, ccds_status, cds_strand, cds_from, cds_to, cds_segs, match_type
  # 0       1      2       3         4        5           6           7         8        9         10
  next if m/^#/;
  chomp;
  my @t=split(/\t/);
  next if $t[5]=~m/Withdrawn/ || $t[9] eq '-' || $t[7]==0;
  my $chr='chr'.$t[0];
  my $id=$t[4];
  #my $chrId=$chr.'|'.$id;
  #push(@ccds, $chrID);
  #parse CDS segments:
  $t[9]=~tr/[] //d;
  my @segs=map { [ split(/\-/) ] } (split(/\,/, $t[9]));
  print STDERR "Warning: start coord not matching lowest segment coord for $id!\n"
    unless $t[7]==$segs[0]->[0];
  print STDERR "Warning: end coord not matching highest segment coord for $id!\n"
    unless $t[8]==$segs[-1]->[1];
  $count++;
  writeCCDS([$id, $chr, $t[6], $t[7], $t[8], \@segs, $t[2]]);
}
close(CCDS);
print STDERR "$count CCDS records written.\n";
# --
if ($outfile) {
  select(STDOUT);
  close(OUTF);
}

#************ Subroutines **************

sub writeCCDS {
 my $rdata=$_[0];
 my ($id, $chr, $strand, $start, $end, $segs, $gene)=@$rdata;
 # if ($bed)
 # if ($tlf)
 my $src=$c2s{$id};
 #-- print GFF
 my $chrId=$chr.'|'.$id;
 print join("\t", $chr, 'CCDS', 'mRNA', $start, $end, ".\t$strand\t.", "ID=$chrId");
 print ";gene=$gene" if $gene;
 if ($src) {
  print ';ncbi_src='.join(',', @{$src->[0]}) if @{$src->[0]}>0;
  print ';ebi_src='.join(',', @{$src->[1]}) if @{$src->[1]}>0;
 }
 else {
  print STDERR "Warning: no sequence source info found for $id ?\n";
 }
 print "\n";
 #print CDS segments:
 foreach my $c (@$segs) {
   print join("\t", $chr, 'CCDS', 'CDS', $$c[0], $$c[1], ".\t$strand\t.", "Parent=$chrId")."\n";
 }
}
