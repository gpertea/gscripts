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
  tabcmp.pl [-o summary.txt] ref.tab.gz test.tab.gz
/;
umask 0002;
getopts('o:') || die($usage."\n");

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

my ($rl, $prev_rl, @rd, $ql, $prev_ql, @qd);
my $r_total;
my $q_total;
my $t_match=0; #total matches (query TP)
my $t_qbs=0; #query badly spliced (should align as unspliced)
my $t_qus=0; #query unspliced (should align as spliced)
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
  $rl=join("\t",@rd[0..5]);
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
     $ql=join("\t",@qd[0..5]);
     #print STDERR "ql=$ql\n";
     #order check:
     die("Error: query data out of order:\n$prev_ql\n$ql\n")
        if $prev_ql && ($prev_ql cmp $ql)>0;
  }
  my $cmp = ($rl cmp $ql);
  if ($cmp == 0) {#same alignment
    $t_match++;
    
    $prev_ql=$ql;
    undef($ql);@qd=();
    next;
  }
  elsif ($qd[0] eq $rd[0]) {
    #same read, check how it does not match
    #TODO: check if it's a qbs or qus?
    #also try disregard the soft clipping
    if ($cmp>0) { #advance to the next query alignment
       $t_qonly++;
       $prev_ql=$ql;
       undef($ql);@qd=();
       goto NEXTQ;
    }
    #$cmp < 0 
    $t_ronly++; 
    next;
  }
  else { #different reads
    if ($cmp>0) {
      # r>q, read next q
      # extra q line found; should only happen if there are 
      #   multiple q alignments for a read
      $t_qonly++;

      $prev_ql=$ql;
      undef($ql);@qd=();
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

print join("\t", $t_ronly, $t_qonly, $t_match)."\n";
print "qtotal=$q_total, rtotal=$r_total\n";



# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }


#************ Subroutines **************

