#!/usr/bin/perl
# -- Example: prepare a partial path/SAMPLE_ID list in ~/cbrain/ like this:
# ls data/libd_bsp*/fastq*/R*.gz | perl -pe 's{^data/(\S+/R\d+).+}{\1}' | \
#   sort -u > bsp_fastq_path_rnums.lsid
# -- a line would look like this: libd_bsp2/fastq_hippo/R11140
# -- then run this script with the resulting lsid file as input
# -- the output file will have records like this:
# >523 libd_bsp2/fastq_hippo R11140 Sex subjID Dx Region Race Age

use Getopt::Std;
use strict;
my $usage = q/Build a grid array job task db for RNAseq samples based on their path and a sample ID prefix for each
Usage:
  build_taskdb.pl [-o taskdb.cfa] [-p pheno_data.tab] sample_basepaths.lsid
/;
umask 0002;
getopts('o:p:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}

## could pass another phenodata sample here
my $host=`uname -n`;
my $basepath=$ENV{HOME}.'/work/cbrain';
## -- expect this same folder structure in all work environments
my $pdf=$Getopt::Std::opt_p || "$basepath/phenodata/samples_phenodata.tab";
die("Error: could not find phenodata file $pdf!\n") unless -f $pdf;

my %pd; # SAMPLE_ID => [subj, region, dx, sex, race, age]
                      #   0       1    2    3    4     5
open(PD, $pdf) || die("Error opening $pdf\n");
while (<PD>) {
  next if $.==1 && m/(dx|age)/i; # skip header
  chomp;
  my @f=split(/[\s\,]/); #expected: sample_id, subj_id, region, dx, sex, race, age
  $pd{$f[0]}=[@f[1..6]];
}
close(PD);
my $n=0;
while (<>) {
  chomp; #expects a line like this: libd_bsp2/fastq_hippo/R11140
  my ($path, $sid)=(m{(.+)/([^/]+)$});
  die("Error parsing sample_id from: $_\n") unless $sid;
  my $sd=$pd{$sid} || die("Error: could not get pheno data for sample $sid\n");
  $n++;
  print join(' ', ">$n", $path, $sid, $$sd[3], @$sd[1,2,4,5])."\n";
}

# -- the end 
if ($outfile) {
 select(STDOUT);
 close(OUTF);
}
