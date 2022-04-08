#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 zcat infinium-omni2-5-8v1-5-a1-manifest-file-csv.zip | \
    gt_prep.pl map [-z] > infinium-omni2-5-8v1-5-a1.map
 or 
    gt_prep.pl rep [-l][-s] gtReport.txt > gtrep.csv

 Commands:
    map: generate plink map file from SNP array manifest csv
         Options:
          -z  discard zero-chromosome entries
          -i  just show header info and exit
    rep: parse gtReport into tab delimited output where 1st column has
         variant IDs, the other columns have the per-sample genotype calls
         Options:
          -s  just list samples and exit
          -i  just show header info and exit
/;
umask 0002;
my %cmds=( 'map'=>1, 'rep'=>1 );
my $cmd=shift(@ARGV) || die($usage, "Error: no command given!\n");
die($usage, "Error: command '$cmd' not recognized!\n") if !$cmds{$cmd};
getopts('zislo:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my $mapCmd=$cmd eq 'map';
my $repCmd=$cmd eq 'rep';
my $info=$Getopt::Std::opt_i;
if ($mapCmd) {
  my $noz=$Getopt::Std::opt_z;
  my $h;
  while (<>) {
  }
  exit;
}

if ($repCmd) {
  my $li=$Getopt::Std::opt_l;
  my $ls=$Getopt::Std::opt_s;
  my ($h, $d);
  die("$usage Error: only one of -l, -s is allowed!\n") if ($ls && $li);
  while (<>) {
    if ($d) {
      s/[\r\n]+$//;
      print $_."\n";
      next;
    }
    die("Error: could not find data!\n") if $.>12;
    die("Error: could not find header!\n") if (!$h && $.>3);
    if (m/\[Header/) { $h=1;next }
    if (m/\[Data/) { 
       exit if $info;
       $d=1;
       $_=<>; # read sample IDs
       s/[\r\n]+$//;
       s/^\s+//;
       my @s=split(/\t/);
       if ($ls) { print join("\n", @s)."\n"; exit }
       print join("\t", 'id', @s)."\n"; #header
    }
  }
  exit;
}


# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

