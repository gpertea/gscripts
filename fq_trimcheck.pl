#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q{Usage:
 fq_trimcheck.pl [-o report.txt] reads_1.trimmed.fastq[.gz] reads_1.fastq[.gz] 
};

umask 0002;
getopts('o:Nm:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
#my $nocheck=$Getopt::Std::opt_N;
#my $forcemate=$Getopt::Std::opt_m;
#die("${usage}Error: invalid -m option!\n") if ($forcemate && $forcemate!=1 && $forcemate!=2);
$outfile='' if $outfile eq '-';
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}

my @f=@ARGV;
#my $paired=(@f==2);
die("$usage") unless @f>1;
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

# >Illumina TruSeq Adapter Read 1
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
# >Illumina TruSeq Adapter Read 2
# AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

my ($ended1, $ended2);
## assumed trimmed output is in the same order
my ($no, $nt); #number original, trimmed
while (1) {
  my ($rname, $rseq, $rquals, $rl);
  unless ($ended1) { # trimmed file, might be missing
   ($rname, $rseq, $rquals, $rl)=getFastq($fh[0]);
   if ($rname) { $nt++ } else { $ended1=1 }
  }
  last if $ended1;
  # find original read in fastq
  my ($mname, $mseq, $mquals, $ml);
  do { # catch up with trimmed read
    ($mname, $mseq, $mquals, $ml)=getFastq($fh[1]);
    if ($mname) { $no++ } else { $ended2=1 }
  } until ($ended2 || $mname eq $rname);
  if ($ended2 && $rname) {
    die("Error: could not find matching original read for ($nt) $rname\n");
  }
  last if $ended2;
  if ($rl!=$ml) {
   my ($d, $tseq)=($ml-$rl, '');
   #if ($d<8) { # show what was trimmed
     my $p=index($mseq, $rseq);
     if ($p==0) { #trimmed at the end
       $tseq='E: '.substr($mseq, -$d);
     } else {
       $tseq='S: '.substr($mseq, 0, $p);
     }
   #}
   print join("\t", $rname, $ml, $rl, $d, $tseq)."\n";
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
 my $slen;
 if ($_) {
   while (<$fh>) {
     last if m/^\+/;
     chomp;
     $rseq.=$_;
     }
   $slen=length($rseq);
   if ($_) {
     while (<$fh>) {
       chomp;
       $rquals.=$_;
       last if ($slen<=length($rquals));
       }
   }
 }
 return ($rname, $rseq, $rquals, $slen);
}
