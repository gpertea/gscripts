#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q{Usage:
  ls `pwd -P`/*/*.fastq.gz | ls2manifest.pl > samples.manifest
  
  Quickly build a samples.manifest file for SPEAQeasy (for paired reads)
};
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
## -- WARNING: this script assumes that input lines are grouped by directories and then
##      by filename (alphanumerically sorted within each directory)
##    AND that the first 2 tokens before the first '_' are a sample identifier! (not a brain identifier!)
## -- filename format can also be like this, 
##    with a different suffix AFTER the _R[12]_ part:
# R2890_C00L0ACXX_TTAGGC_L002_R1_001.fastq.gz
# R2890_C00L0ACXX_TTAGGC_L002_R1_002.fastq.gz
# R2890_C00L0ACXX_TTAGGC_L002_R1_003.fastq.gz
# R2890_C00L0ACXX_TTAGGC_L002_R2_001.fastq.gz
# R2890_C00L0ACXX_TTAGGC_L002_R2_002.fastq.gz
# R2890_C00L0ACXX_TTAGGC_L002_R2_003.fastq.gz
## best way to deal with this is to first collect all the entries 
## with the same prefix before _R[12]_ and sort the group out 
my (@r1, @r2); #list of file names for R1 and R2 respectively
my (@sid, @fc1, @fc2); #sample_id, common pattern between pairs
my (@fn1, @fn2); #sample_id, common pattern between pairs
my $ppre; #previous prefix
while (<>) {
 chomp;
 my ($fn)=(m{/([^/]+)$}); #file name
 die("Error: cannot parse file name!($_)\n") unless $fn;
 my $fp=$_;
 $fp=~s{/[^/]+$}{}; #file path
 my ($d)=($fp=~m{/([^/]+)$}); #last directory
 $fn=~s/\.[gx]z$//i;
 $fn=~s/\.bz2$//i;
 $fn=~s/[._](fq|fastq)$//i;
 my $fc=$fn; #common pattern 
 my $mate;
 if ($fc=~s/_([12])$//) { 
   $mate=$1;
 } elsif ($fc=~s/_R([12])(_[^_]+)$/$2/) {
   $mate=$1;
 }# else { die("Error parsing the mate # from $fn \n"); }

 my @s=split(/_/,$fc);
 my $si=(@s>1)?$s[0].'_'.$s[1] : $s[0];
 if (@s>2 && $s[0]=~m/^\d+$/) { 
   ## for situations like this: 
   #   10_Br5460_SLC17A7pos_HGW2VBBXY_S43..
   #   10_Br5460_SLC17A7pos_HGW2VBBXY_S43..
   #   11_Br5460_SNAP25pos_HGW2VBBXY_S44..
   #   11_Br5460_SNAP25pos_HGW2VBBXY_S44..
   $si.='_'.$s[2];
 }
 my $pre=$d.'/'.$si; # last dir + rnum + flow cell
 #print STDERR "processing: $pre $fn ($fc)\n";
 if ($pre ne $ppre) {
   if (@r1>0) {
     flushData();
     @r1=();@r2=();@sid=();
     @fc1=();@fc2=();
     @fn1=();@fn2=();
   }
   $ppre=$pre;
 }

 if ($mate) {
   if ($mate==1) {
     push(@r1, $_);
     push(@fc1, $fc);
     push(@fn1, $fn);
     push (@sid, $si);
   } else {
     push(@r2, $_);
     # mate 2 always follows mate 1, but it may be after other mate 1s
     my $m=scalar(@fc2);
     die("Error: pair mismatch between $fc ($fn) and\n".
         "                             $fc1[$m] ($fn1[$m]\n")
       unless $fc eq $fc1[$m];
     push(@fc2, $fc);
     push(@fn2, $fn);
   }
   #print "$_\t0\t$_\t0\t$sID\n";
 } else {
   push(@r1, $_);
   push (@sid, $si);
   #print "" 
 }
} #while (<>)

if (@r1>0) {
  flushData();
}
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
}

#************ Subroutines **************
sub flushData {
  next unless @r1>0;
  die("Error: mismatch between (".join(', ', @r1).") and (".join(', ', @sid).")!\n") 
    unless scalar(@r1)==scalar(@sid);
  if (@r2>0) { #paired data, validate pairs
    my $msg='Error: mate mismatch? @r1=('.join(', ',@r1).")\n".
            '                   vs @r2=('.join(', ',@r2).")\n";
    die($msg) unless scalar(@r1)==scalar(@r2);
    #validate and print pairs:
    for (my $i=0;$i<@r1;$i++) {
      my ($m1, $m2, $si)=($r1[$i], $r2[$i], $sid[$i]);
      #validate naming substring prefix ?
      print join("\t",$m1,'0',$m2,'0',$si)."\n";
    }
  }
  else { #assume single-end reads
   for (my $i=0;$i<@r1;$i++) {
     print join("\t",$r1[$i],'0',$sid[$i])."\n";
   }
  }
}
