#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 fastq_deinterleave.pl [-o <outprefix> ] [-F] reads.fq[.{gz|bz2}]
 Converts paired reads given as a single fastq file with interleaved pairs
 into two fastq\/fasta files.
 Use the -F option to output FASTA format instead of FASTQ (i.e. discard
 quality values)
 Default output files names are the input file name with _1 and _2 suffixes
 added (unless option -o is provided).
/;
umask 0002;
getopts('Fo:') || die($usage."\n");
my $outprefix=$Getopt::Std::opt_o;
die($usage."Error: one input file required!\n") unless @ARGV==1;
my $fasta=$Getopt::Std::opt_F;
my $
if ($outprefix) {
  open(OUTF, '>'.$outfile1) || die("Error creating output files $outfiles\n");
  select(OUTF);
  }

my @f=@ARGV;
#my $paired=(@f==2);

my @fh=(undef, undef);
#determine compression, if any
my @fz=('','');
my @foh=(undef, undef);
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
while (1) {
 my ($rname, $rseq, $rquals);
 unless ($ended1) {
  ($rname, $rseq, $rquals)=getFastq($fh[0]);
  $ended1=1 unless $rname;
  }
 my ($mname, $mseq, $mquals);
 unless ($ended2) {
 ($mname, $mseq, $mquals)=getFastq($fh[1]);
  $ended2=1 unless $mname;
 }
 last if ($ended1 && $ended2);
 if ($rname eq $mname) {
   $rname.='/1';
   $mname.='/2';
 }
 if ($fasta) {
  print '>'.$rname."\n$rseq\n" if $rname;
  print '>'.$mname."\n$mseq\n" if $mname;
 }
 else {
  print '@'.$rname."\n$rseq\n+\n$rquals\n" if $rname;
  print '@'.$mname."\n$mseq\n+\n$mquals\n" if $mname;
  }
}

for (my $i=0;$i<@f;$i++) {
  last unless $f[$i];
  close($fh[$i]);
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }


#************ Subroutines **************

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
