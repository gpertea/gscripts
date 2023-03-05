#!/usr/bin/perl
### input: arg1 : FASTQ base directory - relative to current directory 
##                                      OR to ~/work/cbrain/fastq/
##         arg2 (optional) - path to phenodata file, default:
##                           ~/work/cbrain/phenodata/libd_samples_phenodata.tab
##              if '<dataset>:' prefix is found, it is expected to match
##                             the first column in the phenodata file !

### outputs a CHESS-brain taskdb.tfa with the following format:
# >3 dataset_fastq_dir mate1_filemask sex
### which means [1]=FASTQ directory, [2]=mate1 file mask, [3]= M or F
### Expects:
#### dataset_fastq_dir : relative to ~/work/cbrain/fastq/ or current directory
#### mate1_filemask :   sampleID_1.fastq.gz
####                 or RNum_flowcell_*_R1_*.fastq.gz

### requires: ~/work/cbrain/phenodata/libd_samples_phenodata.tab
##            (unless arg2 is provided)
use strict;
my $usage=q/
 prep_taskdb.pl <FASTQ_basepath> [[dataset:]<phenodata_file>]
/;

my $HOME=$ENV{'HOME'};
my $basedir="$HOME/work/cbrain"; #base project directory
my ($fdir, $pdf)=@ARGV;
die("$usage\n") unless $fdir;
my $dataset;
if ($pdf) {
  my @fds=split(/:/, $pdf);
  ($dataset, $pdf)=@fds if ($fds[0] ne $pdf); # dataset given
  $dataset=lc($dataset);
}
$pdf="$basedir/phenodata/libd_samples_phenodata.tab" unless $pdf;

my %pd; # SAMPLE_ID => [subj, region, dx, sex, race, age]
                      #   0       1    2    3    4     5
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
while (<PD>) {
  next if $.==1 && m/(dx|age)/i; # skip header
  chomp;  my @f=split(/[\s\,]/); 
  if ($libd_pd) {
    next if ($dataset && lc($f[0]) ne $dataset);
    ## default libd phenodata fields
    ## dataset rnum sample_id region protocol rin p_id dx sex race age pmi
    ## 0        1       2       3        4     5   6    7  8   9    10 11
    $pd{$f[1]}=[@f[6,3,8,7, 9,10]];
  } else {
    #for anything else but default LIBD phenodata: sample_id, subj_id, region, dx, sex, race, age
    #                                                 0        1         2     3    4    5     6
    $pd{$f[0]}=[@f[1,2,4,3,5,6]]; # p_id, region, sex, dx, race, age
  }
}
close(PD);
my $rn=1;
my %sproc; ## processed sample IDs : sampleID => 1
foreach my $fq (@fqs) {
  if ($cbpath) { #absolute path
    substr($fq, 0, length("$basedir/fastq/"), ''); #remove base path
  }
  my ($fn)=($fq=~m{([^/]+)$});
  my $fqdir=$fq;
  $fqdir=~s{/[^/]+$}{}; # relative path to fastq file to write
  my ($sid, $rnum, $fc, $spd, $fqmask);
  if ($fn=~m/^(R\d+)_([A-Z,0-9]+XX)_.+_R[12]_.+/) {
    ($rnum, $fc)=($1,$2);
    $sid="${rnum}_$fc";
    $spd=$pd{$rnum}  || die("Error retrieving phenodata for $rnum ($fn)\n");
    $fqmask="${sid}_*_R1_*.f*q.gz";
  } elsif ($fn=~m/(.+)_[12]\.f[ast]*q\.gz$/) {
    $sid=$1;
    $spd=$pd{$sid} || die("Error retrieving phenodata for $sid ($fn)\n");
    $fqmask="${sid}_1.f*q.gz";
  }
  unless ($sid) {
    print STDERR "Warning: cannot parse sampleID from $fn!\n";
    next;
  }

  if (!exists($sproc{$sid})) {
     $sproc{$sid}=1;
     my ($subj, $region, $sex, $dx, $race, $age)=@$spd;
     print ">$rn $fqdir $fqmask $sex $dx $race $age\n";
     $rn++;
  }
}
