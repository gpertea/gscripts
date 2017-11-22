#!/usr/bin/perl
use strict;

my $usage = q/Usage:
 wgs_ann_update.pl <ann_update.tab> <annotation.gff3>
 The update table <ann_update.tab> has the following format
 
 seqid trim_l
 
 Where trim_l is the amount to trim at the beginning of seqid
 If trim_l is missing or 0 then the annotation for seqid will be 
 discarded from the input gff3
/;
umask 0002;
die($usage."\n") unless @ARGV==2;
my $utab=shift(@ARGV);

my %useq;
open(UTAB, $utab) || die("Error opening the update table file $utab!\n");
while (<UTAB>) {
 chomp;
 my @s=split;
 my $v=int($s[1]) if @s>1;
 $useq{$s[0]}=$v;
}
close(UTAB);
shift (@ARGV) if ($ARGV[1] eq '-' or $ARGV[1] eq 'stdin');

while (<>) {
  my @t=split(/\t/);
  if (exists($useq{$t[0]})) {
    my $v=$useq{$t[0]};
    next unless $v;
    $t[3]-=$v;
    $t[4]-=$v;
    print join("\t",@t);
  }
  else { print $_; }
}
