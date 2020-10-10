#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q{Usage:
  ls /full/path/*/fastq.gz | ls2manifest.pl > samples.manifest
  
  Quickly build a samples.manifest file for SPEAQeasy (for paired reads)
};
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
# --
my $lastpair;
while (<>) {
 chomp;
 my ($fn)=(m/\/([^\/]+)$/);
 die("Error: cannot parse file name!($_)\n") unless $fn;
 $fn=~s/\.gz$//;
 $fn=~s/[._](fq|fastq)$//;
 my @s=split(/_/,$fn);
 my $pair=$fn;
 my $mate;
 if ($pair=~s/_([12])$//) {
   $mate=$1;
 } elsif ($pair=~s/_R([12])(_[^_]+)$/$2/) {
   $mate=$1;
 } else { die("Error parsing the mate # from $fn \n"); }
 if ($mate==2) {
   die("Error matching mates in $fn ($pair vs $lastpair)\n")
       if $pair ne $lastpair;
   my $sID=(@s>2)?$s[0].'_'.$s[1] : $s[0];
   print "\t$_\t0\t$sID\n";
 } else { print "$_\t0" }
 $lastpair=$pair;
}



# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

