#!/usr/bin/perl
### input: arg1 : FASTQ base directory - relative to current directory 
##                                      OR to ~/work/cbrain/fastq/
##         arg2 : file with a list of RNum to use
##
##         arg3 (optional) - path to phenodata file, default:
##                           ~/work/cbrain/phenodata/libd_samples_phenodata.tab
##              if '<dataset>:' prefix is found, it is expected to match
##                             the first column in the phenodata file !

### outputs a CHESS-brain taskdb.tfa with the following format:
# >3 dataset_fastq_dir mate1_filemask sex sample_ID
### which means [1]=FASTQ directory, [2]=mate1 file mask, [3]= M or F   [4]=sample_ID
## sample_ID should be used as the proposed output base file name for the alignment files etc.
### Expects:
#### dataset_fastq_dir : relative to ~/work/cbrain/fastq/ or current directory
#### a file with a list of RNums to use
########
#### mate1_filemask :  RNum_flowcellXX_*_R1_*.fastq.gz
####                 

### requires: ~/work/cbrain/phenodata/libd_samples_phenodata.tab
##            (unless arg2 is provided)
use strict;
my $usage=q/
 prep_rnums_taskdb.pl <FASTQ_basepath> <rnum_subset.lst> [[dataset:]<phenodata_file>]
/;
## pheno data expected as:
##  0       1      2        3       4      5      6    7   8    9    10   11
##dataset rnum  sample_id region protocol rin subj_id dx  sex  race age  pmi
my $HOME=$ENV{'HOME'};
my $basedir="$HOME/work/cbrain"; #base project directory
my ($fdir, $frnums, $pdf)=@ARGV;
die("$usage\n") unless $fdir && $frnums && -d $fdir && -f $frnums;
my $dataset;
if ($pdf) {
  my @fds=split(/:/, $pdf);
  ($dataset, $pdf)=@fds if ($fds[0] ne $pdf); # dataset given
}
if ($dataset) {
  print STDERR "Dataset requested: $dataset\n";
  $dataset=lc($dataset);
}

$pdf="$basedir/phenodata/libd_samples_phenodata.tab" unless $pdf;
open(FR, $frnums) || die("Error opening $frnums");
my %rnums; # keep track of RNums to keep
while (<FR>) {
 chomp;
 my ($rnum)=(m/^(\S+)/);
 $rnums{$rnum}=1 if length($rnum)>0;
}
close(FR);
my $nrnums=keys(%rnums);
print STDERR "$nrnums RNums requested.\n";

my %pd; # RNum => [subj, region, dx, sex, race, age, sampleID]
                #   0       1    2    3    4     5      6
my $libd_pd=($pdf=~m/libd_samples_phenodata/);

die("$usage\n") unless -f $pdf;
my $fpath=$fdir;
my $cbpath;
if (! -d $fpath) { #not a relative path to current dir?
  $fpath="$basedir/fastq/$fdir";
  $cbpath=1; # bool flag for absolute fqs paths
  die("Error: path $fpath not found!\n") unless -d $fpath;
}
my @fqs=glob("$fpath/*.f*q.gz");
if (@fqs<2) {
  @fqs=glob("$fpath/*/*.f*q.gz");
  die("Error: cannot find FASTQ gz files under $fpath!\n") unless @fqs>=2;
}

open(PD, $pdf) || die("Error opening $pdf\n");
my $pdrnums=0;
while (<PD>) {
  next if $.==1 && m/(dx|age)/i; # skip header
  chomp;  my @f=split(/[\s\,]/); 
  if ($libd_pd) {
    next if ($dataset && (lc($f[0]) ne $dataset));
    next if !exists($rnums{$f[1]});
    $pdrnums++;
    ### libd phenodata fields:
    ## dataset rnum sample_id region protocol rin p_id dx sex race age pmi
    ## 0        1       2       3        4     5   6    7  8   9    10 11
    $pd{$f[1]}=[@f[6,3,8,7, 9,10,2 ]];  # subjID, region, sex, dx, race, age, sample_ID
  }# else {
   # #for anything else but default LIBD phenodata: sample_id, subj_id, region, dx, sex, race, age
   # #                                                 0        1         2     3    4    5     6
   # $pd{$f[0]}=[@f[1,2,4,3,5,6]]; # subj_id, region, sex, dx, race, age
  #}
}
close(PD);

die("Error: multiple sample/pheno data found for the same RNum\n".
"A dataset identifier must be used. ($pdrnums vs $nrnums)\n") if ($pdrnums>$nrnums);

my $rn=0;
my %sproc; ## processed rnums : rnum => 1
foreach my $fq (@fqs) {
  if ($cbpath) { #absolute path
    substr($fq, 0, length("$basedir/fastq/"), ''); #remove base path
  }
  my ($fn)=($fq=~m{([^/]+)$});
  my $fqdir=$fq;

  $fqdir=~s{/[^/]+$}{}; # relative path to fastq file to write
  #$fqdir=~s/[_\.\-]f[ast]*q$//;
  my ($sid, $rnum, $fc, $spd, $fqmask);
  if ($fn=~m/^(R\d+)_.+_R[12]_.+/) {
    ($rnum, $fc)=($1,$2);
    next if !exists($rnums{$rnum});
    #$sid="${rnum}_$fc";
    $spd=$pd{$rnum}  || die("Error retrieving phenodata for $rnum ($fn)\n");
    $sid=$$spd[6] || die("Error getting sample_id for $rnum ($fn)\n");
    $fqmask="${rnum}_*_R1_*.f*q.gz";
  } elsif ($fn=~m/(.+)_[12]\.f[ast]*q\.gz$/) {
    $sid=$1;
    $spd=$pd{$sid} || die("Error retrieving phenodata for $sid ($fn)\n");
    $fqmask="${sid}_1.f*q.gz";
  }
  unless ($sid) {
    print STDERR "Warning: cannot parse sampleID from $fn!\n";
    next;
  }

  if (!exists($sproc{$rnum})) {
     $sproc{$rnum}=1;
     my ($subj, $region, $sex, $dx, $race, $age, $sampleID)=@$spd;
     $rn++;
     print ">$rn $fqdir $fqmask $sex $sampleID $dx $race $age $subj $region\n";
    
  }
}

if ($nrnums != $rn) {
 print STDERR "Warning: $rn records found for $nrnums RNums given\n";
} else {
 print STDERR "$rn records written.\n";
}
