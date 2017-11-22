#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 mgbl_flt.pl [-p<minpid>] [-l<minovl>[%]] [-Q] [-v<maxovh>] [-T] < mgblast_hits

 Use -Q with -l minovl% to enforce minimum overlap percentage of the 
 query length (otherwise the percentage of the shorter length between query
  and subject is used)
 Use -C to inject a 13th column in the output with coverage percentage
 of the query
 
 Use -T to only keep the first ("top") hit for each query sequence (assumes
 hits are already sorted by query_id and a desired score metric !)
/;
umask 0002;
getopts('TCQp:l:v:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $minpid=$Getopt::Std::opt_p || 50.00;
my $minovl=$Getopt::Std::opt_l;
my $Qmin=$Getopt::Std::opt_Q;
my $qtop=$Getopt::Std::opt_T;
my $insQcov=$Getopt::Std::opt_C;
my $ovlperc = ($minovl=~s/\%$//);
my $maxovh=$Getopt::Std::opt_v;

if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my $lastq;
while (<>) {
 s{\t1/1\t}{\t+\t};
 s{\t1/\-1\t}{\t-\t};
 my $line=$_;
 my @t=split(/\t/);
 #my ($q, $qlen, $qstart, $qend, $s, $slen, $sstart, $send,
 #    $pid, $bitscore, $pvalue, $strand, $qgaps, $sgaps)=split(/\t/);

 next if $t[8]<$minpid;
 my ($ql, $qr)=$t[2]<$t[3] ? ($t[2], $t[3]) : ($t[3], $t[2]);
 my $qlen=$t[1];
 my $hlen=$t[5];
 if ($minovl) {
   if ($ovlperc) {
     my $minlen =  $qlen;
     if (!$Qmin) {
        $minlen= ( $qlen<$hlen ? $qlen : $hlen );
     }
     my $movl=int(($minlen*$minovl)/100);
     next if $qr-$ql+1<$movl;
   } else {
      next if $qr-$ql+1<$minovl;
   }
 }
 if ($maxovh) {
   my ($hl, $hr) = $t[6]<$t[7] ? ($t[6], $t[7]) : ($t[7], $t[6]);
   my $l_ovh = $ql>$hl ? $hl-1 : $ql-1;
   my $r_ovh = $qlen-$qr > $hlen-$hr ? $hlen-$hr : $qlen-$qr;
   next if $l_ovh>$maxovh || $r_ovh>$maxovh;
 }

 if ($insQcov) {
  my $qcov=sprintf("%.2f", (100.0*($qr-$ql+1))/$qlen);
  $line=join("\t",@t[0..11],$qcov,@t[12..$#t]);
 }
 next if ($qtop && $t[0] eq $lastq);
 $lastq=$t[0];
 print $line;
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

