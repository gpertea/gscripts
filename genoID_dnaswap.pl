#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
  genoID_dnaswap.pl <dnaswap_cols> <brnum2genoID.ctab>
  
  dnaswap_cols must have tab delimited 5+ cols format:
  
  brnum1 genoID1 brnum_fix genoID_fix ... (drop | swap resolved)
  
  .. where brnum_fix can be "drop"
  
  The script removes entries having "drop" as brnum_fix and replaces "brnum1 genoID1" pairs
  with their "resolved" pairs
/;
umask 0002;
getopts('o:') || die($usage."\n");
my %drop; # "brnum1 genoID1" => 1
my %res;  # "brnum1 genoID1" => "brnum_fix genoID_fix"
die($usage."\n") unless @ARGV==2;
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}

# --
my $fs=shift(@ARGV);
open(F, $fs) || die("Error opening DNA swap table $fs\n");
while (<F>) {
 chomp;
 my @t=split(/\t/);
 if (m/\tdrop/) {
   $drop{$t[0]."\t".$t[1]}=1;
 } elsif (m/swap resolved/) {
   $res{$t[0]."\t".$t[1]}=$t[2]; # new BrNum assigned
         # $t[2]."\t".$t[3];
 }
}
close(F);
#- now parse the input stream
while (<>) {
 chomp;
 my @t=split(/\t/);
 my $brid=$t[0]."\t".$t[1];
 next if ($drop{$brid}); #keep track of dropped?
 my $newbr=$res{$brid};
 if ($newbr) {
   my $oldbr=$t[0];
   $t[0]=$newbr;
   #splice(@t, 1, 1);
   $oldbr=~s/^Br/Brx/;
   push(@t, 'dnaswap:'.$oldbr);
 }
 print join("\t", @t)."\n";
}



# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

