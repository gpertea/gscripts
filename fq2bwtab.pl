#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q/Usage:
 fq2bwtab.pl [-I] [-C] [-o <outfile>] reads_1.fastq[.gz] [ reads_2.fastq[.gz] ]
 
 Takes FASTQ input reads (single, paired or interleaved) and generates
 tab delimited bowtie raw input (for the --12 input option)
 Use -I to specify that the input consists of interleaved paired reads.
 Use -C to disable mate name consistency check (expecting paired to have the same name).
 The input can be 'stdin' or '-' to specify streaming at stdin.
/;
umask 0002;
getopts('ICo:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $interleaved=$Getopt::Std::opt_I;
my $nameCheck = !$Getopt::Std::opt_C;
$outfile='' if $outfile eq '-';
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }

my @f=@ARGV;
my $paired=((@f==2) || $interleaved);
my @fh=(undef, undef);
#determine compression, if any
my @fz=('','');
for (my $i=0;$i<@f;$i++) {
  last unless $f[$i];
  $f[$i]='-' if $f[$i] eq 'stdin';
  if ($f[$i] eq '-') {
    $fh[$i] = *STDIN;
    next;
  }
  if ($f[$i]=~m/\.bzi?p?2$/) {
   $fz[$i]='bzip2';
   }
   elsif ($f[$i]=~m/.g?zi?p?$/) {
   $fz[$i]='gzip';
   }
  if ($fz[$i]) {
    open($fh[$i], $fz[$i]." -cd '".$f[$i]."'|") || 
        die("Error creating decompression pipe: $fz[0] -cd '$f[$i]' !\n");
    }
   else {
    open($fh[$i], $f[$i]) || die("Error opening file $f[$i] ! \n");
    }
}

my $endin;
my $skipget;
my @rbuf;
while (!$endin) {
  my ($rname, $rseq, $rquals)=$skipget ? @rbuf : getFastq($fh[0]);
  $skipget=0;
  if (!$rname) {
    $endin=1;
    last;
  }
  if (!$paired) {
   printBWTab(\$rname, \$rseq, \$rquals);
   next;
  }
  # paired from here on
  my $fhin = $interleaved ? $fh[0] : $fh[1];
  my ($mname, $mseq, $mquals) = getFastq($fhin);
  die("Error: no mate found for read $rname!\n") unless $mname;
  if ($nameCheck && ($rname ne $mname)) { 
     #if ($interleaved) {
     #  #treat this as unpaired
     #  $skipget=1;
     #  @rbuf=($mname, $mseq, $mquals);
     # ($mname, $mseq, $mquals)=(undef, undef, undef);
     #}
     $rname=substr($rname, 0, -2); 
     die("Error: mate names do not match ($rname vs $mname)\n")
        if ($rname ne substr($mname, 0, -2));
      #   #expect mates to be named similarly
  }
  
  #print paired reads
  printBWTab(\$rname, \$rseq, \$rquals, \$mseq, \$mquals);
}

for (my $i=0;$i<@f;$i++) {
  last unless $f[$i];
  close($fh[$i]) unless $f[$i] eq '-';
}

if ($outfile) {
 select(STDOUT);
 close(OUTF);
}

sub getFastq {
 my $fh=$_[0]; # file handle
 #parses next FASTQ record
 #returns ($readname, $readseq, $readquals)
 my ($rname, $rseq, $rquals);
 while (<$fh>) {
   ($rname)=(m/^@(\S+)/);
   last if $rname;
   }
 if ($_) {
   while (<$fh>) {
     last if m/^\+/;
     chomp;
     $rseq.=$_;
     }
   if ($_) {
     while (<$fh>) {
       chomp;
       $rquals.=$_;
       last if (length($rseq)<=length($rquals));
       }
   }
  }
 return ($rname, $rseq, $rquals);
}

sub printBWTab {
 my ($rn, $rs, $rq, $ms, $mq)=@_;
 if ($ms && $$ms) {
   #print pair
   print join("\t", $$rn, $$rs, $$rq, $$ms, $$mq)."\n";
 }
 else {
   #print single
   print join("\t", $$rn, $$rs, $$rq)."\n";
 }
}

sub printSAM {
 my ($rn, $rs, $rq, $ms, $mq)=@_;
 if ($ms && $$ms) {
   #print pair
   print join("\t", $$rn, 77,  '*','0','0','*','*','0','0',$$rs, $$rq)."\n";
   print join("\t", $$rn, 141, '*','0','0','*','*','0','0',$$ms, $$mq)."\n";
 }
 else {
   #print single
   print join("\t", $$rn, 4,  '*','0','0','*','*','0','0',$$rs, $$rq)."\n";
 }
}
