#!/usr/bin/perl
use strict;
use IO::Zlib;  #requires Compress::Zlib, of course
use Getopt::Std;
# this might need command :
#     ulimit -n 10000
my $usage = q/Usage:
 genovcf_split.pl [-m id2br] [-p dir_prefix] [-n numfiles] genotypes.vcf[.gz]
 Options:
   -m : file mapping sample ID to BrNum
   -p : prefix for the output sub-directories
   -n : number of files (samples) per directory
   
/;
getopts('p:n:m:') || die($usage."\n");
my $pre=$Getopt::Std::opt_p || 'vcfsplit';
my $n=$Getopt::Std::opt_n || 500;
my $brf=$Getopt::Std::opt_m;
my $fin=$ARGV[0];
die($usage) if ($fin =~ m/^\-+h/) || ! -f $fin;
my %mbr; # maps sampleID => brnum
if ($brf && -f $brf) {
 open(M, $brf) || die("Error opening $brf\n");
 while (<M>) {
   chomp;
   my ($sid, $brnum)=split;
   $mbr{$sid}=$brnum;
 }
 close(M)
} else { $brf=0 }

my $fh = IO::Zlib->new($fin) || die("Error opening $fin\n");
my $cs=0; #count samples written
my $hdr; #keep common header lines here
my $fdata;
my @samples; 
while (<$fh>) {
  if (@samples) {
   ## reading data
   exit(0); #for now
   next
  }
  ## still in header lines
  if (m/^#CHROM\S+/) { #column headers, we get sample IDs here
    chomp;
    my @t=split(/\t/);
    @samples=@t[9..$#t];
    print STDERR "Found $#samples samples.\n";
    ##TODO: write header here!
    next
  }
  next if m/^##bcftools/;
  $hdr.=$_;
}
