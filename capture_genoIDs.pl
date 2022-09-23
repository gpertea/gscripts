#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
   capture_genoIDs.pl genoIDs.lst *.tab > genoIDs2Brnum.tab 
/;
umask 0002;
getopts('o:') || die($usage."\n");
die($usage."\n") unless @ARGV>1;
my $fids=shift(@ARGV) || die($usage."\n");
my %gids;
open(F, $fids) || die("Error opening genoIDs file: $fids\n");
while(<F>) {
 chomp;
 next unless length>1;
 my @t=split;
 $gids{$t[0]}=1;
}
close(F);

my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
foreach my $ft (@ARGV) {
  open(T, $ft) || die("Error opening $ft\n");
  my $nobrnum;
  while(<T>) {
     chomp;
     tr/"//d; #"
     tr/,/\t/;
     my ($brnum)=(m/\b(Br\d+)/);
     my ($dnum)=(m/\b(D\d+)\b/);
     my ($s)=(m/\b([MFU])\b/);
     my @t=split(/\t/);
     # first look for the genoID just in case it's there already
     my ($gid, @aft);
     foreach my $w (@t) {
       if ($gid) {
          push(@aft, $w);
          next;
       }
       $gid=$w if $gids{$w};
     }
     if (!$gid) {
       my $n=scalar(@t)-1;
       for(my $i=0;$i<$n;$i++) {
         if ($gid) {
           push(@aft, $t[$i]);
           next;
         }
         my $w=$t[$i].'_'.$t[$i+1];
         if ($gids{$w}) {
             $gid=$w;
             $n++;
         }
       }
     }
     next unless ($gid);
     my $extra='';
     unless ($brnum) {
        $brnum='nobr';
        $extra="\t".join("\t", @aft);
     }
     $s='.' unless $s;
     print join("\t",$brnum, $gid, $s, $ft).$extra."\n";
     
  } #while<T>
  close(T);
}
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

