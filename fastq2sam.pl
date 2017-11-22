#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q/Usage:
 fastq2sam.pl [-o <out.sam>] reads_1.fastq[.gz] [ reads_2.fastq[.gz] ]
 
 Takes FASTQ input reads (single or paired) and generates
 SAM output (which can be piped into other conversion programs)
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
$outfile='' if $outfile eq '-';
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }

my @f=@ARGV;
my $paired=(@f==2);
my @fh=(undef, undef);
#determine compression, if any
my @fz=('','');
for (my $i=0;$i<@f;$i++) {
  last unless $f[$i];
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

my ($ended1, $ended2);
print "\@HD\tVN:1.4\tSO:unsorted\n";
print "\@PG\tID:fastq2sam.pl\tPN:fastq2sam\n";

while (1) {
  my ($rname, $rseq, $rquals);
  unless ($ended1) {
   ($rname, $rseq, $rquals)=getFastq($fh[0]);
   $ended1=1 unless $rname;
   }
  if (!$paired) {
   last if $ended1;
   printSAM(\$rname, \$rseq, \$rquals);
   next;
  }
  my ($mname, $mseq, $mquals);
  unless ($ended2) {
  ($mname, $mseq, $mquals)=getFastq($fh[1]);
   $ended2=1 unless $mname;
  }
  last if ($ended1 && $ended2);
  if ($rname && $mname && $rname ne $mname) {
     $rname=substr($rname, 0, -2); 
     if ($rname ne substr($mname, 0, -2)) {
        #expect mates to be named similarly
        die("Error: mate names do not match ($rname vs $mname)\n");
     }
  }
  if (!$rname) { #print single read
    printSAM(\$mname, \$mseq, \$mquals);
    
  }
  else { #print paired reads
  printSAM(\$rname, \$rseq, \$rquals, \$mseq, \$mquals);
  }
}

for (my $i=0;$i<@f;$i++) {
  last unless $f[$i];
  close($fh[$i]);
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
