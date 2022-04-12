#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 zcat infinium-omni2-5-8v1-5-a1-manifest-file-csv.zip | \
    gt_prep.pl map [-z] > infinium-omni2-5-8v1-5-a1.map
 or 
    gt_prep.pl rep [-l][-v][-i] gtReport.txt > gtrep.csv

 Commands:
    map: generate plink map file from the SNP array manifest (CSV)
         Options:
          -z  discard zero-chromosome entries
          -i  just show header info and exit
          -v  list variant IDs

    rep: parse gtReport into tab delimited output where 1st column has
         variant IDs and the other columns have the per-sample genotype calls
         Options:
          -l  just list samples and exit
          -i  just show header info and exit
          -v  list variant IDs

/;
umask 0002;
my %cmds=( 'map'=>1, 'rep'=>1 );
my $cmd=shift(@ARGV) || die($usage, "Error: no command given!\n");
die($usage, "Error: command '$cmd' not recognized!\n") if !$cmds{$cmd};
getopts('zivlo:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my $mapCmd=$cmd eq 'map';
my $repCmd=$cmd eq 'rep';
my $info=$Getopt::Std::opt_i;
my $lv=$Getopt::Std::opt_v;
## genome build is assumed to be 37 !
if ($mapCmd) {
  my %mk; @mk{'Name', 'Chr', 'MapInfo'}=();
  my $noz=$Getopt::Std::opt_z;
  my ($h, $d);
  my @hix; # indexes for columns of interest (mk)
  while (<>) {
    if ($d) {
      last if (m/^\s*\[Control/);
      my @t=split(/\,/);
      next if $noz && $t[$hix[1]] eq '0'; # assumes Chr is 2nd index!
      if ($lv) {
        print $t[$hix[0]]."\n"; #assume Name is the first index
      } else { # print map lines:
        print join("\t",@t[$hix[1],$hix[0]],0,$t[$hix[2]])."\n";
      }
      next;
    } # data 
    # before data
    die("Error: could not find data!\n") if $.>12;
    die("Error: could not find header!\n") if (!$h && $.>3);
    if (m/\[Head/) { $h=1;next }
    if (m/\[Assay/) { 
       last if $info;
       $d=1;
       $_=<>; # read header
       s/[\r\n]+$//;
       my @hdr=split(/\,/);
       my $i=-1;
       @hix=map { $i++; exists($mk{$_}) ? ($i) : () } @hdr;
       die("Error parsing assay header!\n") if @hix!=3;
    }
    if ($h && $info) {
      s/[\r\n]+$//;
      s/\s*\,\s*/:\t/;
      print $_."\n";
    }
  }  
}

if ($repCmd) {
  my $ls=$Getopt::Std::opt_l;
  my ($h, $d); # in_header, in_data
  die("$usage Error: only one of -l, -v, -i is allowed!\n") if ($ls && $lv) || ($info && ($ls || $lv));
  while (<>) {
    if ($d) {
      if ($lv) {
        ($_)=(m/^(\S+)/)
      } else {
        s/[\r\n]+$//;
      }
      print $_."\n";
      next;
    }
    # before data:
    die("Error: could not find data!\n") if $.>12;
    die("Error: could not find header!\n") if (!$h && $.>3);
    if (m/\[Head/) { $h=1;next }
    if (m/\[Data/) { 
       last if $info;
       $d=1;
       $_=<>; # read sample IDs
       s/[\r\n]+$//;
       s/^\s+//;
       my @s=split(/\t/);
       if ($ls) { print join("\n", @s)."\n"; last }
       print join("\t", 'id', @s)."\n" unless $lv; #header
       next;
    }
    if ($h && $info) {
      s/[\r\n]+$//;
      print $_."\n";
    }
  }
}


# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

