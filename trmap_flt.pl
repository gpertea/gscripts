#!/usr/bin/perl
use strict;
my $usage=q/Usage: 
  trmap_flt.pl 'codes' in.ovltab 
Filters trmap's output by the overlap codes given
in the 'codes' string. Example:
  trmap_flt.pl '=c' < trmap.ovltab > trmap.match_or_contained.ovltab
This example will output only the overlaps with class codes '=' (match)
or 'c' (contained)./;

my $codes=shift(@ARGV);
die($usage."\n") if (!$codes || $codes eq '-h' || 
 $codes eq '--help' || -f $codes);
my %c= map { $_=>1 } (split(//,$codes));
my ($h, $hpr);
while (<>) {
 my $l=$_;
 if (m/^>/) {
   #query transcript header
   ($h, $hpr) = ($_, 1);
 }
 else {
   #overlapped reference info
   my @d=split(/\t/);
   if (exists($c{$d[0]})) {
     if ($hpr) {
       print $h;
       $hpr=0;
     }
     print $l;
   }
 }
}
