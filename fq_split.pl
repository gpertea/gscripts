#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 fq_split.pl [-o <outprefix> ] <num_reads_per_file> reads_1.fq[.gz] [ reads_2.fq[.gz] ]
 Extract [paired] reads into numbered files <outprefix>_0N_1.fq, etc. 
 with at most <num_reads_per_file>
/;

umask 0002;
getopts('o:') || die($usage."\n");
my $outprefix=$Getopt::Std::opt_o || 'fq_split';
die($usage."\n") unless @ARGV>0 && @ARGV<4;
my $rpf=shift(@ARGV); #reads per file
$rpf=~tr/,//d;
die($usage."\n") unless $rpf>1;
my $partno=0; #output part#
my $ircount=0; #total input reads counter
my $orcount=0; #total output reads counter

my $rcount=0; #read count in current output file
my @fh=(undef, undef);
#determine compression, if any
my @fz=('','');
my @fno;
#=($outprefix.'_1.fq', $outprefix.'_2.fq');

my @ifn=@ARGV;
my $paired=(@ifn==2);
my @foh=(undef, undef);

initOutPart($partno);
$partno++;
for (my $i=0;$i<@ifn;$i++) {
  last unless $ifn[$i];
  if ($ifn[$i]=~m/\.bzi?p?2$/) {
   $fz[$i]='bzip2';
   }
   elsif ($ifn[$i]=~m/.g?zi?p?$/) {
   $fz[$i]='gzip';
   }
  if ($fz[$i]) {
    open($fh[$i], $fz[$i]." -cd '".$ifn[$i]."'|") || 
        die("Error creating decompression pipe: $fz[0] -cd '$ifn[$i]' !\n");
    }
   else {
    open($fh[$i], $ifn[$i]) || die("Error opening file $ifn[$i] ! \n");
    }
}

while (1) {
 my ($rname, $rseq, $rquals)=getFastq($fh[0]);
 last unless $rname;
 $ircount++;
 my $iname=$rname;
 $iname=~s/\/[12]$//;
 my ($mname, $mseq, $mquals);
 if ($ifn[1]) {
     ($mname, $mseq, $mquals)=getFastq($fh[1]);
     my $check=$mname;
     $check=~s/\/[12]$//;
     die("Error: cannot find mate for $rname (found $mname instead)!\n") 
          unless $check eq $iname;
     }
 #next unless exists($h{$iname});
 print { $foh[0] } '@'.$rname."\n$rseq\n+\n$rquals\n";
 $orcount++;
 if ($mname) {
     print { $foh[1] } '@'.$mname."\n$mseq\n+\n$mquals\n";
     }
 $rcount++;
 if ($rcount>=$rpf) {
   initOutPart($partno);
   $partno++;
   $rcount=0;
 }
} # for each input read/pair

print STDERR "Closing files with input read count $ircount (out count: $orcount).\n";
for (my $i=0;$i<@ifn;$i++) {
   last unless $ifn[$i];
   close($foh[$i]);
   close($fh[$i]);
}




#************ Subroutines **************

sub closeFo {
  #print STDERR "closing.\n";
  close($foh[0]);
  close($foh[1]) if $paired;
}

sub initOutPart {
  my $p=$_[0]+1;
  my $part=sprintf('p%02d', $p);
  print STDERR "..writing part $p (${outprefix}_$part*.fq) \[i:$ircount, o:$orcount\]\n";
  closeFo() if $p>1;
  if ($paired) {
     @fno=($outprefix."_${part}_1.fq", $outprefix."_${part}_2.fq");
  }
  else {
   $fno[0]=$outprefix."_$part.fq";
  }
  for (my $i=0;$i<@ifn;$i++) {
    last unless $ifn[$i]; 
    #print STDERR "opening file $i ($fno[$i])\n";
    open($foh[$i], '>'.$fno[$i]) || die("Error opening file $fno[$i] !\n");
  }
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
