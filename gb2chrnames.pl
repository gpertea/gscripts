#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 gb2chrnames.pl GCF_..._assembly_report.txt [gencode.fa.fai] > ids2ucsc.chrmap

 Provides a mapping of the accessions from the GenBank assembly report to 
 UCSC genome sequence names
 -- OR, if a fasta index file is also given as a 2nd parameter, attempts 
 to provide a mapping table between those sequences names and 
 UCSC sequence names.

 The output is a 2-column translation table which can be used with 
 gffread's -m option to rename GFF annotation data.

/;
umask 0002;

getopts('o:') || die($usage."\n");
die("$usage\n") unless (@ARGV==1 || @ARGV==2);

my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my @gaccs;  # genome assembly accessions from the assembly report (7th column)
my %gacc2ucsc; # genome assembly accessions to UCSC sequence names
#my %chr2gacc; # chrName or GB accession => Refseq genomic assembly accession
my %chr2glen; # chrName or GB accession => seq length
my %gb2ucsc;  # chrName or GB accession => UCSC name 
my @fachr;
open(GBTBL, $ARGV[0]) || die("Error opening assembly report file $ARGV[0]");

while(<GBTBL>) {
  next if m/^#/;
  s/[\n\r]+$//; #stupid DOS file
  my @t=split(/\t/);
  if ($t[1] eq 'assembled-molecule') {
    $t[0]='M' if lc($t[0]) eq 'mt';
    print STDERR "Warning: suspicious chromosome name: chr$t[0] !\n"
       unless ($t[0]=~m/^\d\d/ || $t[0]=~m/[\dXYM]/);
    my $k='chr'.$t[0];
    #$chr2gacc{$k}=$t[6];
    $chr2glen{$k}=$t[8];
    push(@gaccs, $t[6]);
    $gacc2ucsc{$t[6]}=$k;
    $gb2ucsc{$k}=$k;
    $gb2ucsc{$t[4]}=$k;
    if ($t[9] ne 'na') {
      print STDERR "Warning: UCSC style name ($t[9]) not matching expected name $k\n"
       if $t[9] ne $k;
    }
    next;
  }
  # --anything else should be a one of: 
  #  alt-scaffold, unlocalized-scaffold, unplaced-scaffold, novel-patch, fix-patch
  my $gbacc=$t[4];
  my $uid=$t[9];
  #$chr2gacc{$gbacc}=$t[6];
  if (length($uid)<3) {
    #$chr2gacc{$t[9]}=$t[6];
    $uid='';
    if ($t[1]=~m/(novel|fix)\-patch/) {
      my $suf=($1 eq 'novel') ? '_alt' : '_fix';
      $t[2]='M' if $t[2] eq 'MT';
      my $pre=($t[2] eq 'na') ? 'chrUn_' : 'chr'.$t[2].'_';
      my $accv=$gbacc;
      $accv=~tr/./v/;
      $uid=$pre.$accv.$suf;
    }
  }
  print STDERR "Warning: no UCSC id can be guessed from the assembly report for $t[1] $t[4]\n" if !$uid;
  my $gacc=$t[6];
  $gacc=$t[4] if ($gacc eq 'na');
  $gb2ucsc{$gbacc}=$uid;
  push(@gaccs, $gacc);
  $gacc2ucsc{$gacc}=$uid;
  $chr2glen{$gbacc}=$t[8];
}
close(GBTBL);
if (!$ARGV[1]) {
  foreach my $a (@gaccs) {
    my $ucsc=$gacc2ucsc{$a};
    die("Error: no UCSC name parsed for genomic accession $a\n") unless $ucsc;
    print "$a\t$ucsc\n";
  }
  exit(0);  
}
# second file given, try to recognize and map those sequence names to ucsc names
open(CHRNAMES, $ARGV[1]) || die("Error opening file faidx file $ARGV[1]\n");
#now go through fasta names and match them in %gb2ucsc
while (<CHRNAMES>) {
  next if m/^#/;
  my @t=split(/\t/); # $t[0] = UCSC/Gencode chr name
  my $gbid; #ID to find in %gb2ucsc
  if ($t[0]=~m/^[\dXYM][\dTt]?$/) {
    $gbid='chr'.$t[0];
    $gbid='chrM' if lc($t[0]) eq 'chrmt';
  }
  elsif ($t[0]=~m/^chr[\dXYM][\dTt]?$/) {
    # assembled molecule (chromosome)
    $gbid=$t[0];
    $gbid='chrM' if lc($t[0]) eq 'chrmt';
  }
  elsif ($t[0]=~m/^[A-Z][A-Z]_?\d+\.\d+$/) {
   # plain accession.version entry
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
  my $uid = $gb2ucsc{$gbid};
  if (!$uid) {
    $uid = $gacc2ucsc{$gbid};
  }
  if ($uid eq 'na') {
   print STDERR "Warning: $gbid has 'na' mapping in GenBank sequence data, converting to self-mapping\n";
   $uid=$gbid;
  }
  if ($uid) {
     print "$t[0]\t$uid\n";
     print STDERR "WARNING: $uid / $t[0] sequence length mismatch !\n" 
         if $t[1]!=$chr2glen{$gbid};
  } else {
     print STDERR "WARNING: $gbid ($t[0]) sequence name could not be mapped!\n";
  }
}
close(CHRNAMES);
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
}
