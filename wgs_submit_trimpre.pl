#!/usr/bin/perl
use strict;
#prepare the trimming data  WGS submission screening 
#Expected input format (tab delimited): 
# <scaffold> <scaff_len> <match_start> <match_end> <match_target>
# which should be obtained from NCBIs trim report like this:
# perl -ne 'split(/\t/,$_,4);foreach(split(/\,/,$_[2]))'\
# '{s/\.\./\t/;print join("\t",@_[0..1],$_,$_[3])}' trim_report.tab | \
#  sort -k1,1 -k3,3n > to_trim.srt.tab
my $usage=q{
 Usage:
  wgs_submit_trimpre.pl genome.fa.cidx to_trim.tab | sort -k1,1 -k3,3n > to_trim.sorted.tab
  Takes the NCBI screening results for trimming and produces
  a tab delimited format with coordinates mapped on scaffolds properly:
    scaffold <sequence_length> trim_start trim_end vector/adaptor

};
die $usage unless @ARGV ;
die $usage if $ARGV[0] eq '-h';
my $faidx=shift(@ARGV);
die("$usage\nError finding genome index $faidx\n") unless -f $faidx;
#die $usage unless -f $ARGV[0];

while (<>) {
  s/\.\./\t/g;
  my @t=split(/\t/);
  if ($t[0]=~m/\~(\d+)$/) {
    my $cstart=$1;
    $t[0]=~s/\~\d+$//;
    my $cend=$t[1];
    splice(@t,1,1);
    $t[1]=getSeqLen($t[0]);
    $t[2]=$cstart+$t[2]-1;
    $t[3]=$cstart+$t[3]-1;
    print join("\t",@t);
  }
  else {
    print join("\t",@t);
  }
}

sub getSeqLen {
 my $id=$_[0];
 open(CDB, "cdbyank -a '$id' $faidx |") || die("Error at cdbyank -a '$id' $faidx \n");
 my $r=0;
 while (<CDB>) {
   next if m/^>/;
   chomp;
   $r+=length($_);
 }
 close(CDB);
 return $r;
}
