#!/usr/bin/perl
use strict;

my $usage = q/Usage:
 ann_remap.pl <ann_update.tab> <annotation.gff3>
 The update table <ann_update.tab> has the following format
 
 seqid trim_l seqlen [new_seqid]
 
 Where trim_l is the amount to trim at the beginning of seqid, 
 while seqlen is the new sequence length (any annotation past
 that coordinate will be discarded instead of being shifted by trim_l)

 trim_l can be a negative value which means that new sequence has 
 been inserted at the beginning (or remapped).
 
 Optionally the reference sequence can be renamed in the output if a 4th
 column value is specified as the new ID.
 
 Special cases:
  if trim_l is 'd' or 'del': any annotation on seqid is discarded 
                             (seqlen is not required)
  if trim_l is 'r' or 'rc': annotation is reverse-complemented
                            (seqlen is required)
 
/;
umask 0002;
die($usage."\n") unless @ARGV==2;
my $utab=shift(@ARGV);

my %useq;
open(UTAB, $utab) || die("Error opening the update table file $utab!\n");
while (<UTAB>) {
 chomp;
 my @s=split;
 #my $v=int($s[1]) if @s>1;
 my $id=shift(@s);
 if (@s>0) {
   if ($s[0] eq '0') {
     $id.=':0';
   }
   $s[0]='d' if lc($s[0]) eq 'del';
   $useq{$id}=[@s];
 }
}
close(UTAB);
shift (@ARGV) if ($ARGV[1] eq '-' or $ARGV[1] eq 'stdin');

while (<>) {
  my @t=split(/\t/);
  my $d=$useq{$t[0]};
  next if $d && $$d[0] eq 'd';
  my $d0=$useq{$t[0].':0'};
  if ($d) {
    my ($v, $rmax, $nid)=@$d;
    if ($v eq 'r' || $v eq 'rc') {
      die("Error: cannot reverse complement without sequence length!\n") 
        unless $rmax && $rmax>=$t[4];
     my ($ol, $or)=($t[3], $t[4]);
     $t[6] = ($t[6] eq '-') ? '+' : '-';
     $t[4] = $rmax-$ol+1;
     $t[3] = $rmax-$or+1;
    } else {
     #regular coordinate shift
     if ($d0) {
       my ($v0, $rmax0, $nid0)=@$d0;
       die("Error: wrong entry for $t[0]:0 (v=$v0)!\n") unless $v0 eq '0' && $rmax0>0;
       if ($t[3]<$rmax0 && $t[4]>$rmax0) {
         print STDERR "Warning: $t[2] ($t[3]-$t[4]) annotation discarded".
            "due to crossing split boundary $t[8]\n";
         next;
       }
       if ($t[4]<=$rmax0) {
          $nid=$nid0;
          $v=$v0; #0
       }
     }
     next if $rmax && $t[4]>$rmax;
     next if $t[3]-$v<1;
     $t[3]-=$v;
     $t[4]-=$v;
    }
    if ($nid) { 
       my $oid=$t[0];
       $t[0]=$nid; 
       $t[8]=~s/$oid/$nid/g;
    }
    print join("\t",@t);
  }
  else { print $_; }
}
