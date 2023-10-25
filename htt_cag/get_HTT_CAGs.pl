#!/usr/bin/perl
use strict;
use Cwd qw();
my ($sid) = (Cwd::cwd()=~m/([^\/]+)$/);

my $cmd="samtools view htt_alns.bam | grep -P '(CAG){15,}' | fgrep 'AS:i:10' | cut -f 1-3,4,10,12- |";
my %qs; #keep track of read names, discard duplicate entries (keep the one with higher score)

my $showaln = $ARGV[0]=~m/\-?[aA]/;
## first pass - keep only exact match entries across the CAG STR (flanked or not)
#                          0       1     2        3     4      5     6   7
my @alns; ## array of [HTT-nnn, qname, qpos, flanked, xCAG, flags, $as, seq] 
## flankmatch can be 0, 1 or 2 (no flanks matched, one flank matched, both flanks matched
open(S, $cmd) || die("Cannot open stream from cmd:\n$cmd\n");
my %refspan; ## HTT-nn => [min-start, max-end] keep track of region with alignments in each ref
           
my %freq; # frequency of flanked cag counts (>= 15x)
my %ufreq; # unbounded cag count frequences only
my $last='';
my $maxflanked=0;
my $lastp;
my ($q, $f, $ref, $p, $seq, $cag);
while (<S>) {
 chomp;
 my ($as)=(m/AS:i:(\d+)/); #alignment score
 #     1   2    3    4    5       6      7       8       9     10     11
 #my ($q, $f, $ref, $p, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual)=split(/\t/);
 ($q, $f, $ref, $p, $seq)=split(/\t/);
 # reject non-exact matches
 $_=$seq;
 next if ($as==100 && !(m/^N/ || m/N$/));
 my $sflag = $f & 0xc0; # mate #
 if ($sflag == 0x40) {
    $q.='/1';
 } elsif ($sflag == 0x80) {
    $q.='/2';
 }
 next if exists($qs{$q}); #already seen this read with a full match
 $qs{$q}=1;
 ($cag)=(m/((CAG){15,})/);
 my $cagcount=0;
 if ($cag) {
    $cagcount=int(length($cag)/3);
 } else { next }
 my $lflanked=m/T?CCAGCAGCAGCAGCAG/;
 my $rflanked=m/CAGCAGCAGCAA/;
 my $flanked=$rflanked+$lflanked;
 if ($flanked == 2) {
   $maxflanked=$cagcount if $cagcount>$maxflanked;
   $freq{$cagcount}++;
 } else {  $ufreq{$cagcount}++ }
 if ($last ne $ref) {
   $refspan{$ref}=[$p, 0];
   $refspan{$last}->[1]=$lastp+100 if $last;
   $last=$ref;
 }
 $lastp=$p;
 push(@alns, [$ref, $q, $p, $flanked, $cagcount, $f, $as, $seq]);
}
$refspan{$last}->[1]=$lastp+100 if $last;

## second pass: show full alignments, or show CAG count support
#@alns = sort { return $$a[0] eq $$b[0] ? $$b[3] <=> $$a[3] : $$a[0] cmp $$b[0] } @alns
$last='';
my ($st, $e);
if ($showaln) {
  foreach my $aln (@alns) {
    my ($ref, $q, $p, $flanked, $cagcount, $f, $as, $seq)=@$aln;
   if ($last ne $ref) { #reference block start
     my $span=$refspan{$ref} || die("Error getting refspan for $ref!\n");
     $last=$ref;
     ($st, $e)=@$span;
     my $s=`seqmanip -LX -r $st-$e ../../$ref.tfa`;
     chomp($s);
     print ">$ref $st-$e\n$s\n";
   }
   my $ns=$p-$st;
   print ' ' x $ns;
   my $es=' ' x ($e-$p-100);
   my $x=($flanked==2) ? "[$cagcount]" : "$cagcount+";
   print "$seq $es $x\n";
  }
} else { # just show frequencies:
  my @counts=sort {$a <=> $b} (keys(%freq));
  foreach my $c (@counts) {
     print join("\t", $sid,$c,$freq{$c}, 'F')."\n";
  }
  ## unbounded
  @counts=sort {$a <=> $b} (keys(%ufreq));
  foreach my $c (@counts) {
     print join("\t", $sid,$c,$ufreq{$c}, 'u')."\n" if $c>=$maxflanked;
  }
}

