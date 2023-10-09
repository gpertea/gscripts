#!/usr/bin/perl
use strict;
#my $st=0; #min ref position
my @bam=glob('R*.cag_reg.bam');
die("Error: could not find cag_reg BAM file") unless @bam>0;

my ($sid)=($bam[0]=~m/(R\d+)/);

my $cmd="samtools view $bam[0] 'ENST00000355072.9|ENSG00000197386.10|OTTHUMG00000159916.5|OTTHUMT00000358234.3|HTT-001|HTT|13474|protein_coding|:192-253'".
  " 'ENST00000999999.9|ENSG00000197386.10|OTTHUMG00000159916.5|OTTHUMT00000358234.3|HTT-099|HTT|13474|protein_coding|:197-307' |";
my %qs; #keep track of read names
open(S, $cmd) || die("Cannot open stream from cmd:\n$cmd\n");
my $last='';
my %freq; # freq of >=9 x CAGs
while (<S>) {
 chomp;
 my ($q, $f, $ref, $p, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual)=split(/\t/);
 #if ($last && $last ne $ref) {
 #}
 my $sflag = $f & 0xc0; # mate #
 if ($sflag == 0x40) {
    $q.='/1';
 } elsif ($sflag == 0x80) {
    $q.='/2';
 }
 next if exists($qs{$q}); #already seen this read
 $qs{$q}=1;
 my ($cag)=($seq=~m/((CAG){13,})/);
 if ($cag) {
    my $c=int(length($cag)/3);
    $freq{$c}++;
 }
 #if ($last ne $ref) {
 #  $st=$p if $st==0;
 #  $last=$ref;
 #  my $e=$st+240;
 #  my $s=`samtools faidx ../../htt_2tx.tfa $ref:$st-$e|tail -1`;
 #  chomp($s);
 #  print ">$ref\n$s\n";
 #}
 #my $ns=$p-$st;
 #print ' ' x $ns;
 #print $seq."\n";
}
print "sample\tct\tfreq\n";
showFreq();

sub showFreq {
  my @counts=sort {$b <=> $a} (keys(%freq));
  foreach my $c (@counts) {
     print "$sid\t$c\t".$freq{$c}."\n";
  }
}
