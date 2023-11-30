#!/usr/bin/perl
use strict;
use Cwd qw();
my ($sid) = (Cwd::cwd()=~m/([^\/]+)$/);
my $trefpath='/dcs05/lieber/gpergola/Madhur_Projects/htt/tref';
die("Cannot find tx refs in $trefpath/\n") unless -f "$trefpath/HTT-017.tfa";

my $cmd="samtools view htt_alns.bam | grep -P '(CAG){15,}' | cut -f 1-4,6,10,12- |";

my $showaln = $ARGV[0]=~m/^\-?[aA]$/;
## first pass - keep only exact match entries across the CAG STR (flanked or not)
#                          0       1     2        3     4      5     6   7
my @alns; ## array of [HTT-nnn, qname, qpos, flanked, xCAG, flags, AS, seq] 
## flankmatch can be 0, 1 or 2 (no flanks matched, one flank matched, both flanks matched
my %qs; #keep track of read names, discard duplicate entries (keep the one with higher score)
##  qs{read_name} => [ AS, @alns-idx ]

my %refspan; ## HTT-nn => [min-start, max-end] keep track of region with alignments in each ref

my $rlen;
open(S, $cmd) || die("Cannot open stream from cmd:\n$cmd\n");
my ($q, $f, $ref, $p, $cig, $seq, $tags); # parsing SAM 1-4,6,10,12-
my $cag;
## the BAM file was sorted so alignments are grouped by tx and position
## however, reads can be reassigned to another ref out of turn if a better alignment score is found
while (<S>) {
 chomp;
 my ($as)=(m/AS:i:(\d+)/); #alignment score
 #     1   2    3    4    5       6      7       8       9     10     11
 #my ($q, $f, $ref, $p, $mapq, $cigar, $rnext, $pnext, $tlen, $seq, $qual)=split(/\t/);
 ($q, $f, $ref, $p, $cig, $seq, $tags)=split(/\t/, $_, 7);
 $_=$seq;
 $rlen=length;
 my $minscore=$rlen-1; # allow one mismatch!
 $minscore-- if (m/^N/ || m/N$/);
 next if $as<$minscore; #ignore if more than 1 mismatch
 my $sflag = $f & 0xc0; # mate#
 if ($sflag == 0x40) {
    $q.='/1';
 } elsif ($sflag == 0x80) {
    $q.='/2';
 }
 my $raln=$qs{$q}; # did this read have an alignment already?
 next if $raln && $$raln[0]>=$as; #already seen this read with a better/same score

 ($cag)=(m/((CAG){15,})/);
 next unless $cag;
 my $cagcount=int(length($cag)/3);

 my $lflanked=m/T?CCAGCAGCAGCAGCAG/;
 my $rflanked=m/CAGCAGCAGCAGCAGCAA/;
 my $flanked=$rflanked+$lflanked;
 if ($raln) { #replace previous alignment (better score)
   $alns[$$raln[1]]=[$ref, $q, $p, $flanked, $cagcount, $f, $as, $seq];
 } else { # it'll be pushed into @alns
   $qs{$q}=[$as, scalar(@alns)];
   push(@alns, [$ref, $q, $p, $flanked, $cagcount, $f, $as, $seq]);
 }
}

## show full alignments or show CAG count support
if ($showaln) {
  ## sort by ref, position first
  @alns = sort { return ($$a[0] eq $$b[0]) ? $$a[2] <=> $$b[2] : $$a[0] cmp $$b[0] } @alns;
  ## first pass: get ref spans
  my $last='';
  my $lastp;
  foreach my $aln (@alns) {
    my ($ref, $q, $p, @rest)=@$aln;
    if ($last ne $ref) {
      $refspan{$ref}=[$p, 0];
      $refspan{$last}->[1]=$lastp+$rlen-1 if $last;
      $last=$ref;
    }
    $lastp=$p;
  }
  $refspan{$last}->[1]=$lastp+$rlen-1 if $last;
  my ($st, $e);
  foreach my $aln (@alns) {
    my ($ref, $q, $p, $flanked, $cagcount, $f, $as, $seq)=@$aln;
    if ($last ne $ref) { #reference block start
       my $span=$refspan{$ref} || die("Error getting refspan for $ref!\n");
       $last=$ref;
       ($st, $e)=@$span;
       my $s=`seqmanip -LX -r $st-$e $trefpath/$ref.tfa` || die("Error getting range $st-$e from trefpath/$ref.tfa\n");
       chomp($s);
       print ">$ref $st-$e\n$s\n";
    }
    my $ns=$p-$st;
    print ' ' x $ns;
    my $es=' ' x ($e-$p-$rlen+2);
    my $x=($flanked==2) ? "[$cagcount]" : " $cagcount+";
    #print "$seq${es}$x $as\n"; -- print alignment score?
    print "$seq${es}$x\n";
  }
  exit;
}
# just show cagcounts and their read support frequencies:
my %freq; # frequency of flanked cag counts (>= 17x)
my %ufreq; # unbounded cag count frequences only
my $maxfl=0;
## sort by ref, decreasing flanked status
## @alns = sort { return ($$a[0] eq $$b[0]) ? $$b[3] <=> $$a[3] : $$a[0] cmp $$b[0] } @alns;
foreach my $aln (@alns) {
    my ($ref, $q, $p, $flanked, $cagcount, $f, $as, $seq)=@$aln;
    if ($flanked == 2) {
      $maxfl=$cagcount if $cagcount>$maxfl;
      $freq{$cagcount}++;
    } else { $ufreq{$cagcount}++ }
}
my @counts=sort {$a <=> $b} (keys(%freq));
foreach my $c (@counts) {
     print join("\t", $sid,$c,$freq{$c}, 'F')."\n";
}
## unbounded
@counts=sort {$a <=> $b} (keys(%ufreq));
foreach my $c (@counts) {
     print join("\t", $sid,$c,$ufreq{$c}, 'u')."\n" if $c>=$maxfl;
}
