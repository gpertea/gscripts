#!/usr/bin/perl
use strict;
use Getopt::Std;
#use IO::Handle;
#use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
  tabix_merge.pl [-a col2_append ] qry_to_merge tabix_db.gz > merged_tabix_db
  
  Merges a given set of intervals into reference intervals from a tabix 
  database. 
  
  tabix_db.gz is the tabix indexed database.

  qry_to_merge is a BED format set of intervals to be merged into 
  overlapping intervals of tabix_db.gz (if such overlaps exist); if 
  this parameter is '-' then input is expected from stdin
  
  If -a is provided, the provided string is appended to the 2nd column 
  of the merged_tabix_db output.
/;
umask 0002;
getopts('Do:a:') || die($usage."\n");
my $debug=$Getopt::Std::opt_D;
my $outfile=$Getopt::Std::opt_o;
my $mappend=$Getopt::Std::opt_a;
die($usage) unless @ARGV==2;
my ($fqry, $tabixdb) = @ARGV;
die($usage."Error: could not find tabix file $tabixdb\n") unless -f $tabixdb;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
# --
my %mrepl; #hash keeping track of replaced tabix db lines

my $pid=open(TABIXRD, "-|");
defined($pid) || die("Error: cannot fork!\n");
if ($pid) { #parent process
  #stdout from child process comes into TABIXRD
  my ($qctg, $qstart, $qend);
  my ($pmctg, $pmstart, $pmend, $pmline); #previously merged region and line entry to replace
  while(<TABIXRD>) { 
    if (m/^##Qry\|(.+):(\d+)\-(\d+)$/) {
      ($qctg, $qstart, $qend)=($1,$2,$3);
    }
    else { #db interval
      print STDERR "# [$qctg $qstart..$qend] : $_" if $debug;
      chomp;
      my $mline=$_;
      my @t=split(/\t/);
      $t[1].=','.$mappend if $mappend;
      my ($mctg, $mstart, $mend)=($t[0],$t[3],$t[4]);
      #assume overlap! no need to test for it
      $mstart=$qstart if ($mstart>$qstart);
      $mend=$qend if ($mend<$qend);
      ($t[3],$t[4])=($mstart, $mend);
      if ($mctg eq $pmctg) {
        #check for overlap and merge the regions as needed
        if ($mstart<=$pmend && $pmstart <= $mend) {
          #overlap with previous merge
          $mstart=$pmstart if ($mstart>$pmstart);
          $mend=$pmend if ($mend<$pmend);
          ($t[3],$t[4])=($mstart, $mend);
          $mrepl{$pmline}=join("\t",@t);
        }
      }
      $mrepl{$mline}=join("\t",@t);
      ($pmctg, $pmstart, $pmend, $pmline)=($mctg, $mstart, $mend, $mline);
      print STDERR join("\t",@t)."\n" if $debug;
    }
  }
  close(TABIXRD);
  #now output the modified tabix db with replaced/merged entries
  open(TABIXDB, "gzip -cd $tabixdb |") || die("Error opening $tabixdb pipe\n");
  my $last;
  while (<TABIXDB>) {
    my $l=$_;
    chomp($l);
    my $rl=$mrepl{$l};
    if ($rl) {
      print $rl."\n" if ($rl ne $last);
    }
    else { print $l."\n"; }
  }
  close(TABIXDB);
}
else { #child process
 #STDOUT from here goes to TABIXRD above
 my $qfh;
 my $tabixcmd="tabix -Q -R - $tabixdb";
 if ($fqry eq '-') {
   $qfh=*STDIN;
 }
 else {
  open($qfh, $fqry) || die("Error: could not open file $fqry for reading!\n");
 }
 open(TABIXRUN, "| $tabixcmd") || die("Error creating pipe to $tabixcmd\n");
 while (<$qfh>) {
   print TABIXRUN $_;
 }
 close(TABIXRUN);
 close($qfh) if $fqry ne '-';
 exit(0);
}



# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
}

#************ Subroutines **************

