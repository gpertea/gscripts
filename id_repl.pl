#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 id_repl.pl <ID2replID.tab>  ...
 Replace IDs with new IDs based on a 2 column correspondence table
 mapping the IDs to its replacement (tab, space or comma delimited).
 Outputs the result of such replacement applied to the input data.
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --

my $rtab=shift(@ARGV);
die("${usage}Error: no ID mapping table provided!\n") unless $rtab && -f $rtab;

my %map; # oldID => newID
open(TAB, $rtab) || die("Error opening file $rtab!\n$!\n");
while (<TAB>) {
  chomp;
  next unless length>2;
  my ($id, $nid)=split(/[ \t\,]/);
  die("Error: invalid ID mapping on line:\n$_\n") unless length($id)>0;
  die("Error: duplicate entry for ID $id \n") if exists($map{$id});
  $map{$id}=$nid;
}
close(TAB);
my @ids=keys(%map);
while (<>) {
 foreach my $id (@ids) {
  my $nid=$map{$id};
  s/\b$id\b/$nid/g;
 }
 print $_;
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
}

#************ Subroutines **************

