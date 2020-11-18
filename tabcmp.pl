#!/usr/bin/perl
use strict;
use IO::File;
use IO::Zlib qw(:gzip_external 0);
use Getopt::Std;
use POSIX qw(locale_h);

my $usage = q/Compare accuracy of sam alignments.
Input is compressed samread -T alignment data, sorted 
alphanumerically (with LC_ALL=C sort)
Usage:
  tabcmp.pl [-o summary.txt] [-c col1,col2..] ref.tab.gz test.tab.gz
/;
umask 0002;
getopts('o:c:') || die($usage."\n");

my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
# --
die("$usage\n") unless @ARGV==2 && -f $ARGV[0] && -f $ARGV[1];

my ($fr, $fq);
if ($ARGV[0]=~m/\.gz$/) {
  $fr=new IO::Zlib;
  $fr->open($ARGV[0], "rb") || die("Error: could not open $ARGV[0]!\n");
} else {
  $fr=IO::File->new();
  $fr->open("< $ARGV[0]") || die("Error: could not open $ARGV[0]!\n");
}
if ($ARGV[1]=~m/\.gz$/) {
  $fq=new IO::Zlib;
  $fq->open($ARGV[1], "rb") || die("Error: could not open $ARGV[1]!\n");
} else {
  $fq=IO::File->new();
  $fq->open("< $ARGV[1]") || die("Error: could not open $ARGV[1]!\n");
}
my @cols;
my $cols=$Getopt::Std::opt_c;
@cols=split(/\,/,$cols);

my ($rl, $prev_rl, @rd, $rni, $ql, $prev_ql, @qd, $qni);
my ($r_total, $rj_total, $q_total, $qj_total);
my $t_match=0; #total alignment matches (query TP)
my $tj_match=0; #total spliced alignment matches
my $t_qbs=0; #query badly spliced (should be unspliced)
my $t_qus=0; #query unspliced (should be spliced)
my $t_ronly=0;
my $t_qonly=0;
my $qeof=0;
while (1) {
  $prev_rl=$rl;
  $rl=<$fr>;
  last unless $rl;
  $r_total++;
  #exit(0) if ($r_total==3);
  #print STDERR ">$r_total\n";
  if ($qeof) {
    $t_ronly++;
    next;
  }
  chomp($rl);
  @rd=split(/\t/,$rl);#readID, start, strand, cigar, exons
  $rl=join("\t",@rd[0..4]);
  $rni=($rd[5]=~tr/,//);
  $rj_total++ if $rni;
  die("Error: ref data out of order:\n$prev_rl\n$rl\n")
        if $prev_rl && ($prev_rl cmp $rl)>0;
  #                       0      1      2       3      4
  #print STDERR "rl=$rl\n";
NEXTQ:
  unless ($ql) { #read next qry line
     $ql=<$fq>;
     unless ($ql) {
       $qeof=1;
       $t_ronly++;
       next;
     }
     $q_total++;
     chomp($ql);
     @qd=split(/\t/,$ql);
     if ($qd[4]=~m/^(\d+)S\d+M/) {
      $qd[2]-=$1;
      $qd[4]='101M';
     }
     $ql=join("\t",@qd[0..4]);
     $qni=($qd[5]=~tr/,//);
     $qj_total++ if $qni>0;
     #print STDERR "ql=$ql\n";
     #order check:
     die("Error: query data out of order:\n$prev_ql\n$ql\n")
        if $prev_ql && ($prev_ql cmp $ql)>0;
  }
  my $cmp = ($rl cmp $ql);
  if ($cmp == 0) {#same alignment
    $t_match++;
    $tj_match++ if $rni>0;
    $prev_ql=$ql;
    undef($ql);@qd=();undef($qni);
    next;
  }
  elsif ($qd[0] eq $rd[0]) {
    #same read, check how it does not match
    $t_qus++ if $rni>$qni; #missed spliced alignment
    $t_qbs++ if $qni>$rni; #incorrectly spliced alignment
    if ($cmp>0) { #advance to the next query alignment
       $t_qonly++; #incorrect (misplaced) alignment
       $prev_ql=$ql;
       undef($ql);@qd=();undef($qni);
       goto NEXTQ;
    }
    #$cmp < 0 
    $t_ronly++; #missed alignment
    next;
  }
  else { #different reads
    if ($cmp>0) {
      # r>q, read next q
      # extra q line found; should only happen if there are 
      #   multiple q alignments for a read
      $t_qonly++;
      $prev_ql=$ql;
      undef($ql);@qd=();undef($qni);
      goto NEXTQ;
    } else {
      # r<q, read next r 
      # q alignment missing for rl, check same ql vs next rl
      $t_ronly++;
      next;
    }
  }
}
$fr->close;
undef $fr;
$fq->close;
undef $fq;

print join("\t", @cols, $r_total, $rj_total, $q_total, $qj_total,  
                $t_match, $tj_match,  $t_ronly, $t_qonly, $t_qbs, $t_qus)."\n";
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }


#************ Subroutines **************

