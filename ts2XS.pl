#!/usr/bin/perl
use strict;
my %refs;
while (<>) {
 my $line=$_;
 #check for header:
 if (m/^@[A-Z][A-Z]\t/) {
   print $_;
   #keep refseq length
   if (m/^\@SQ\tSN:(\S+)\tLN:(\d+)/) {
     $refs{$1}=$2;
   }
   next;
 }
 chomp;
 my ($rname, $flags, $refname, $pos, $mapq, $cigarstr, 
     $rnext, $pnext, $tlen, $seq, $quals, $tags)=split(/\t/, $_, 12);
 $flags=int($flags);
 my $isrev = ($flags & 0x10) != 0;
 my $sflag = $flags & 0xc0;
 my $frag_idx=0;
 if ($sflag == 0x40) {
    $rname.='/1';
    $frag_idx=1;
 } elsif ($sflag == 0x80) {
    $rname.='/2';
    $frag_idx=2;
 }
 my ($ts)=($tags=~m/\bts:A:(.)/);
 my $xs;
 if ($ts eq '+') {
   $xs = $isrev ? '-' : '+';
 } elsif ($ts eq '-') {
   $xs = $isrev ? '+' : '-';
 }
 if ($xs) {
  $tags.="\tXS:A:$xs";
 }
 print join("\t", $rname, $flags, $refname, $pos, $mapq, $cigarstr, 
     $rnext, $pnext, $tlen, $seq, $quals, $tags)."\n";
}
