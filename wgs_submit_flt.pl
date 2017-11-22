#!/usr/bin/perl
use strict;
#Process WGS submission screening 
#Expected input format (tab delimited): 
# <scaffold> <scaff_len> <match_start> <match_end> <match_target>
# which should be obtained from NCBIs trim report like this:
# perl -ne 'split(/\t/,$_,4);foreach(split(/\,/,$_[2]))'\
# '{s/\.\./\t/;print join("\t",@_[0..1],$_,$_[3])}' trim_report.tab | \
#  sort -k1,1 -k3,3n > to_trim.srt.tab
my $usage=q{
 Usage:
  wgs_submit_flt.pl <cdbfasta_index> [<sorted_match_spans_file>]
  Output is a multi-FASTA with the trim/masked/split records
  (some of the sequences will be discarded in the process).
  
  The input can be given like this:
   perl -pe 's/\.\./\t/g' < trim_report.tab | \
    sort -k1,1 -k3,3n | wgs_submit_flt.pl genome.fasta.cidx \
    > trim_split.out.fasta
};
die $usage unless @ARGV;
my $cdb=shift(@ARGV);
die("$usage Error: cannot find $cdb cdbfasta index.\n") unless -f $cdb;
my $gfa=$cdb;
$gfa=~s/\.cidx$//;
open(LOG, ">$gfa.flt.log") || 
  die("Error creating log file wgs_submit_flt.log\n");
my $sid; #current contig/scaffold being analyzed
my $slen; #length of current contig/scaffold;
my (@matches, @ranges);
while (<>) {
 chomp;
 next if m/^#/;
 my ($id, $len, $m_start, $m_end, $smatch)=split(/\t/);
 die("Error: input tab file seems invalid!\n") unless $len && $smatch;
 if ($sid ne $id) {
    processMatches();
    ($sid, $slen)=($id, $len);
 }
 #my @matches=split(/\,/, $smatches);
 push (@matches, $smatch);
 #my @ranges=map { $_=[ split(/[\-\.]+/) ] } ( split(/\,/, $sranges) );
 push (@ranges, [$m_start, $m_end]);
} #while match lines

processMatches();
close(LOG);

sub processMatches {
 return unless $sid;
 my @vmatches=grep { m/adaptor/ || m/vector/ } @matches;
 my $vecmatch = ( scalar(@vmatches)==scalar(@matches) );
 #print $_."\n" if @vmatches>1;
 #cleanup and merge overlapping or close ranges
 my $numr=@ranges;
 @ranges=mergeIntervals(\@ranges);
 if ($numr>@ranges) {
   print LOG "#####info: $sid : $numr ranges merged into ".scalar(@ranges)."\n";
 }
 
 my $mlen = 0;
 map { $mlen+=$$_[1]-$$_[0]+1 } @ranges;
 if ($mlen/$slen > 0.42 || $slen-$mlen<300) {
    print LOG ">$sid discarded (len $slen, match len=$mlen)\n";
    @matches=(); @ranges=(); 
    return;
 }
 #get the sequence to look for N spans
 my $farec=`cdbyank -a '$sid' $cdb`;
 die("Error getting FASTA record for $sid\n") 
           unless length($farec)>length($sid)+$slen;
 my ($faid, $seq)=($farec=~m/^>(\S+)[^\n]*\n(.+)/s);
 $seq=~tr/\n \t//d;
 die("Error: non matching id $faid for FASTA record for $sid\n") 
     unless $faid eq $sid;
 my $faseqlen=length($seq);
 die("Error: non matching sequence length $slen != $faseqlen for FASTA record $sid\n") 
    unless ($slen == $faseqlen);
 #print STDOUT "pulled $sid ($slen) sequence\n";
 #if (@nr>0) {
 #  my @nrs= map { $_="$$_[0]-$$_[1]" } @nr;
 #  print STDOUT "$sid N-ranges: ".join(',',@nrs)."\n";
 #}
 my ($ltrim, $rtrim)=(0,0);
 if ($ranges[0]->[0] < 300) {
      $ltrim=$ranges[0]->[1];
      #shift(@ranges); 
      print LOG "## $sid match trim left $ltrim\n";
 }
 if ($ranges[-1]->[1]+300>$slen) {
      $rtrim=$slen-$ranges[-1]->[0];
      print LOG "## $sid match trim right $rtrim\n";
      #pop(@ranges);
 }
 if (@ranges==1 && $ltrim && $rtrim) {
   @ranges=();
 }
 else {
   if ($ltrim) {
     shift(@ranges);
     #should adjust remaining ranges by ltrim ?
   }
   pop(@ranges) if $rtrim;
 }
 if (@ranges==0) {
   my $newlen=$slen-$rtrim-$ltrim;
   if ($newlen<300) {
      print LOG ">$sid (new len=$newlen) discarded after trimming only ($ltrim, $rtrim)\n";
      @matches=(); @ranges=();
      return;
   }
 }
 my @nr;
 getNpos(\@nr, \$seq, $slen);
 my @rngsplit; #only ranges to split
 #check if any match ranges are near a N-gap, 
 # and if so extend the gap into the range
 # otherwise, push it into @rngsplit
 if ($ltrim) {
   unshift(@ranges, [1,$ltrim]);
 }
 foreach my $rng (@ranges) {
   my ($ml, $mr)=@$rng;
   my $nExt=0; #rng stuck to a gap
   my $prevn; # previous range to delete if overlapping
   for my $n (@nr) {
     if ($$n[1]<$ml && $ml-$$n[1]<100) {
       my $xlen=$mr-$$n[1]+1;
       substr($seq, $$n[1]-1, $xlen) = 'N' x $xlen;
       print LOG "## $sid gap $$n[0]-$$n[1] extended r to $$n[0]-$mr\n";
       $$n[1]=$mr;
       $nExt=1;
     }
     if ($$n[0]>$mr && $$n[0]-$mr<100) {
       my $xlen=$$n[0]-$ml+1;
       substr($seq, $ml-1, $xlen) = 'N' x $xlen;
       print LOG "## $sid gap $$n[0]-$$n[1] extended l to $ml-$$n[1]\n";
       $$n[0]=$ml;
       #check if a previous range is being overlapped
       $nExt=1;
     }
     if ($prevn && $$n[0]< $$prevn[1]) {
        $$n[0]=$$prevn[0];
        $$prevn[0]=0; #delete previous N-range
       }
    $prevn=$n; 
   }
   @nr = grep { $_->[0]>0 } @nr;
   #split around the match range if not merged into a gap
   push (@rngsplit, [$ml, $mr]) if $nExt==0;
 }
 if (@rngsplit>0) {
    my @prng; #ranges to print
    #check first part for a ltrim too close to a gap
    if (@nr>0 && $nr[0]->[0]-$ltrim<250) {
     $ltrim=$nr[0]->[1];
     print LOG "## $sid first part gap-extend ltrim to $ltrim\n"
    }
    push(@prng, [$ltrim+1, $rngsplit[0]->[0]]) if $rngsplit[0]->[0]+250>$ltrim;
    for (my $r=0;$r<@rngsplit;$r++) {
      if ($r+1<@rngsplit) {
          push(@prng, [$rngsplit[$r]->[1]+1, $rngsplit[$r+1]->[0]-1] );
      }
      else { #last part
          if (@nr>0 && $nr[-1]->[1]<250) {
            $rtrim=$slen-$nr[-1]->[0]+1;
            print LOG "## $sid last part gap-extend rtrim to $rtrim\n"
          }
          push(@prng, [$rngsplit[$r]->[1]+1, $slen-$rtrim]);
      }
    }
    my $p=1;
    print LOG ">$sid split in ".scalar(@prng)." parts:\n";
    foreach my $pr (@prng) {
      my $plen=$$pr[1]-$$pr[0]+1;
      my $toosmall=($plen<300);
      print LOG "## $sid.p$p = $sid\:$$pr[0]-$$pr[1]".
          ($toosmall?" (too small, discarded)":"")."\n";
      unless ($toosmall) {
        print ">$sid.p$p\n";
        printFasta(substr($seq, $$pr[0]-1, $plen));
      }
      $p++;
    }
 }
 else { #print whole sequence
      my $newlen=$slen-$rtrim-$ltrim;
      if ($newlen>=300 && @nr>0) {
        while (@nr>0 && $nr[0]->[0]-$ltrim<250) {
           $ltrim=$nr[0]->[1];
           print LOG "## $sid gap-extend ltrim to $ltrim\n";
           shift(@nr);
        }
        while (@nr>0 && $slen-$rtrim-$nr[-1]->[1]<250) {
           $rtrim=$slen-$nr[-1]->[0]+1;
           print LOG "## $sid gap-extend rtrim to $rtrim\n";
           pop(@nr);
        }
        $newlen=$slen-$rtrim-$ltrim;
      }
      if ($newlen>=300) {
        print LOG ">$sid (new len = $newlen) trimmed (ltrim=$ltrim, rtrim=$rtrim)\n";
        print ">$sid\n";
        printFasta(substr($seq, $ltrim, $newlen));
      }
      else {
        print LOG ">$sid (new len = $newlen) discarded after trimming (ltrim=$ltrim, rtrim=$rtrim)\n";
      }
 }
 @matches=(); @ranges=();
}

sub mergeIntervals { #merge [almost] overlapping intervals
 my $aref=shift(@_);
 my @r=@$aref;
 my $merged=1;
 #intervals are given in increasing begin coordinate
 while ($merged) {
   $merged=0;
   my $i=0;
   while ($i<$#r) {
    if ($r[$i+1]->[0]<=$r[$i]->[1]+100) {
      $r[$i]->[1]=$r[$i+1]->[1];
      splice(@r,$i+1,1);
      $merged=1;
      last;
    }
    ++$i;
   }
 }
 return @r;
}

sub getNpos {
  my ($ret, $sref, $slen)=@_;
  #print "[DBG]:>SEQ>$seqname\t".length($seq)."\n";
  #my @nranges; #array of [$p_start, $p_end] for each base in @PChars (A,C,G,T,N)
  my $cseq=$$sref;
  my $p=0;
  while (length($cseq)>=10) {
     if ($cseq=~m/(N{10,})/i) {
      my ($p_start, $p_end)=($-[0], $+[0]);
      $p+=$p_start;
      my $pl=$p_end-$p_start;
      push(@$ret, [$p+1,  $p+$pl]);
      $cseq=substr($cseq, $p_end);
      $p+=$pl;
     }
     else { last; }
  } #while matching a poly
}

sub printFasta {
 my $slen=length($_[0]);
 my @lines=unpack("A100" x (int(($slen-1)/100)+1),$_[0]);
 print join("\n",@lines)."\n";
}
