#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q{Usage:
 fastq_deinterleave.pl [-o <outprefix> ] [-F] [-I] reads.fq[.(gz|bz2)]
 Converts paired reads given as a single fastq file with interleaved pairs
 into two fastq\/fasta files.
 
 Use -I option to ignore read name suffixes (i.e. don't look for /1 or /2)
 and simply write alternatively the input reads into each output file.
 
 Use the -F option to output FASTA format instead of FASTQ (i.e. discard
 quality values)
  
 Default output files names are the input file name with _1 and _2 suffixes
 added (unless option -o is provided).
 
};
umask 0002;
getopts('FWIo:') || die($usage."\n");
my $outprefix=$Getopt::Std::opt_o;
my $ignoreSuffix=$Getopt::Std::opt_I;
my $noWarn=$Getopt::Std::opt_W;
die($usage."Error: one input file required!\n") unless @ARGV==1;
my $fasta=$Getopt::Std::opt_F;
my ($ih, $oh1, $oh2);
my ($outfile1, $outfile2);
my $fn=shift(@ARGV);
my $fnpre;
if ($fn eq '-' || $fn eq 'stdin' ) {
    $ih=*STDIN;
    $fnpre='stdin';
}
else {
  my $fz;
  if ($fn=~m/\.bzi?p?2$/i) {
   $fz='bzip2';
  }
   elsif ($fn=~m/.g?zi?p?$/i) {
   $fz='gzip';
  }
  $fnpre=$fn;
  if ($fz) {
    open($ih, $fz." -cd '".$fn."'|") || 
        die("Error creating decompression pipe: $fz -cd '$fn' !\n");
    $fnpre=~s/\.[bgzip2]+$//i;
  }
  else {
    open($ih, $fn) || die("Error opening file $fn !\n");
  }
  $fnpre=~s/\.f[ast]*q$//i;
}
if ($outprefix) {
   $outfile1=$outprefix.'_1.fq';
   $outfile2=$outprefix.'_2.fq';
}
else {
  $outfile1=$fnpre.'_1.fq';
  $outfile2=$fnpre.'_2.fq';
}

open($oh1, '>'.$outfile1) || die("Error creating output file $outfile1\n");
open($oh2, '>'.$outfile2) || die("Error creating output file $outfile2\n");

my $counter=0;
while (1) {
 my ($rname, $rseq, $rquals)=getFastq($ih);
 last unless $rname;
 $counter++;
 my $oh;
 if ($ignoreSuffix) {
   $oh = (($counter % 2) == 0) ? $oh2 : $oh1;
 }
 else {
   my ($suf)=($rname=~m/[\/\.\:_]([12])$/);
   if ($suf) {
       $oh=($suf eq '1') ? $oh1 : $oh2;
       if (!$noWarn) {
         if (($suf eq '2' && ($counter % 2)!=0) ||
              ($suf eq '1' && ($counter % 2)!=1)) {
           print STDERR "Warning: suffix for $rname (read $counter) does not match expected interleaved order.\n";
         }
       }
   }
   else {
    die("Error: no /1 or /2 suffix found for input reads, use -I to ignore this check.\n");
   }
 }
 
 if ($fasta) {
  print $oh '>'.$rname."\n$rseq\n" if $rname;
 }
 else {
  print $oh '@'.$rname."\n$rseq\n+\n$rquals\n" if $rname;
 }
}

close($oh1);
close($oh2);

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
