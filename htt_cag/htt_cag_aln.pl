#!/usr/bin/perl
use strict;
my $st=0; #min ref position
#my @bam=glob('R*.cag_reg.bam');
#die("Error: could not find cag_reg BAM file") unless @bam>0;
#my ($sid)=($bam[0]=~m/(R\d+)/);

my $cmd="samtools view htt_alns.bam | grep -P '(CAG){17,}' | fgrep 'AS:i:101' | cut -f 1-3,4,10 |";
my %qs; #keep track of read names

open(S, $cmd) || die("Cannot open stream from cmd:\n$cmd\n");
my $last='';
my %freq; # freq of >=9 x CAGs
while (<S>) {
 chomp;
 #     1   2    3    4    5       6      7       8       9     10     11
 #my ($q, $f, $ref, $p, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual)=split(/\t/);
 my ($q, $f, $ref, $p, $seq)=split(/\t/);
 my $sflag = $f & 0xc0; # mate #
 if ($sflag == 0x40) {
    $q.='/1';
 } elsif ($sflag == 0x80) {
    $q.='/2';
 }
 next if exists($qs{$q}); #already seen this read
 $qs{$q}=1;
 my ($cag)=($seq=~m/((CAG){17,})/);
 my $c=0;
 if ($cag) {
    $c=int(length($cag)/3);
    $freq{$c}++;
 } else { next }
 if ($last ne $ref) {
   $st=$p; #if $st==0;
   $last=$ref;
   my $e=$st+110;
   my $s=`seqmanip -LX -r $st-$e ../../$ref.tfa`;
   chomp($s);
   print ">$ref $st-$e\n$s\n";
 }
 my $ns=$p-$st;
 print ' ' x $ns;
 print $seq."\n";
}
#print "sample\tct\tfreq\n";
#showFreq();

#sub showFreq {
#  my @counts=sort {$b <=> $a} (keys(%freq));
#  foreach my $c (@counts) {
#     print "$sid\t$c\t".$freq{$c}."\n";
#  }
#}
