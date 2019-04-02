#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 gb2chrnames.pl GCF_..._assembly_report.txt genome.fa.fai > gb2chr.txt

 Matches the chromosome names from GenBank assembly report to the names found
 in the given fasta file index (which must be of UCSC or Gencode origin).
 This will create a 2-column translation table which can be used by 
 gffread's -m option.
/;
umask 0002;

getopts('o:') || die($usage."\n");
die("$usage\n") unless @ARGV==2 && -f $ARGV[0] && -f $ARGV[1];

my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --

my %chr2gb; # chrNum or GB accession => RefSeq accession
my %chr2gblen; # chrNum or GB accession => RefSeq length
my @fachr;
open(GBTBL, $ARGV[0]) || die("Error opening assembly report file $ARGV[0]");
open(CHRNAMES, $ARGV[1]) || die("Error opening file faidx file $ARGV[1]\n");

while(<GBTBL>) {
  next if m/^#/;
  my @t=split(/\t/);
  if ($t[1] eq 'assembled-molecule') {
    $t[0]='M' if lc($t[0]) eq 'mt';
    print STDERR "Warning: suspicious chromosome name: chr$t[0] !\n"
       unless ($t[0]=~m/^\d\d/ || $t[0]=~m/[\dXYM]/);
    $chr2gb{'chr'.$t[0]}=$t[6];
    $chr2gblen{'chr'.$t[0]}=$t[8];
    next;
  }
  # --anything else should be a one of: 
  #  alt-scaffold, unlocalized-scaffold, unplaced-scaffold, novel-patch, fix-patch
  my $gbacc=$t[4];
  $chr2gb{$gbacc}=$t[6];
  $chr2gblen{$gbacc}=$t[8];
}
close(GBTBL);

#now go through fasta names and match them in %chr2gb
while (<CHRNAMES>) {
  next if m/^#/;
  my @t=split(/\t/); # $t[0] = UCSC/Gencode chr name
  my $gbid; #ID to find in %chr2gb
  if ($t[0]=~m/^chr[\dXYM][\dTt]?$/) {
    # assembled molecule (chromosome)
    $gbid=$t[0];
    $gbid='chrM' if lc($t[0]) eq 'chrmt';
  }
  elsif ($t[0]=~m/^[A-Z][A-Z]\d+\.\d+$/) {
   # Gencode: plain accession.version
   $gbid=$t[0];
  }
  elsif ($t[0]=~m/chr[\dXYMTUn]{1,2}_([A-Z][A-Z]\d+v\d)/) {
    $gbid=$1;
    $gbid=~s/v/./;
  }
  else {
    print STDERR "WARNING: unexpected sequence id format: $t[0]\n";
    next;
  }
  my $racc = $chr2gb{$gbid};
  if ($racc eq 'na') {
   print STDERR "Warning: $gbid has 'na' mapping in GenBank sequence data, converting to self-mapping\n";
   $racc=$gbid;
  }
  if ($racc) {
     print "$racc\t$t[0]\n";
     print STDERR "WARNING: $racc / $t[0] sequence length mismatch !\n" 
         if $t[1]!=$chr2gblen{$gbid};
  } else {
     print STDERR "WARNING: $gbid chromosome/scaffold name not found!\n";
  }
}
close(CHRNAMES);
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
}
