#!/usr/bin/perl
use strict;
use Getopt::Std;
##
my $usage = q/Usage:
  prep_taskdb.pl [-p pheno_data.tab] fastq_basedir
  
  Generate a pseudo-fasta records of the format:

>234 sampleID F
path sampleID_r1.fq.gz sampleID_r2.fq.gz

The sex (M or F) will only follow the record number and 
sampleID if a pheno_data.tab file was provided with -p
/;
umask 0002;
getopts('p:o:') || die($usage."\n");
die($usage."\nError: fastq directory missing ($ARGV[0])!\n") 
  if (@ARGV==0 || ! -d $ARGV[0]); 
my $outfile=$Getopt::Std::opt_o;
my $pdfile=$Getopt::Std::opt_p;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
my $fdir=$ARGV[0];
my %pd; # sample_id => [subj_id, dx, sex, age, race]
my $pdat=0;
if ($pdfile && -f $pdfile) {
  open(PD, $pdfile) || die("Error opening $pdfile\n");
  #header expected:
  $_=<PD>; chomp;
  my @h=split(/\t/);
  #find sex column
  my ($sidx) = grep { lc($h[$_]) eq 'sex' } (0 .. $#h);
  die("Error: Sex column not found!\n") unless defined $sidx;
  my ($dxidx) = grep { lc($h[$_]) eq 'dx' } (0 .. $#h);
  die("Error: Dx column not found!\n") unless defined $dxidx;
  my ($ageidx) = grep { lc($h[$_]) eq 'age' } (0 .. $#h);
  die("Error: Age column not found!\n") unless defined $ageidx;

  my ($raceidx) = grep { lc($h[$_]) eq 'race' } (0 .. $#h);
  die("Error: Race column not found!\n") unless defined $raceidx;
  ## assume first column is sample_id, 2nd is subj_id
  while(<PD>) {
    chomp;
    my @t=split(/\t/);
    $pd{$t[0]}=[ @t[1, $dxidx, $sidx, $ageidx, $raceidx] ];
  }
  close(PD);
  $pdat=scalar(%pd);
}

# --
my $cmd="find $fdir -name '*.gz'";
#print STDERR " running: $cmd\n";
my @ffns=sort(split(/\n/, `$cmd`));
#print join("\n", @ffns)."\n";
my ($prevd, $prevf)=('','');
my $c=1;
foreach my $fp (@ffns) {
  chomp($fp);
  my ($d, $f)=($fp=~m{(.+)/([^/]+)$});
  utf8::downgrade($f);
  if ($d eq $prevd && length($f)==length($prevf)) {
    my $xor = ($f ^ $prevf);
    my @r;
    if ((@r=($xor=~m/[^\0]/g))==1) {
       my $p=$-[0];
       if (substr($prevf, $p, 1) eq '1' && substr($f, $p, 1) eq '2') {
         my ($r1, $r2)=(substr($f, 0, $p+1), substr($f, 0, $p+1));
         $r1=~s/[_\.]r?[12]$//i; 
         $r2=~s/[_\.]r?[12]$//i;
         die("Error: could not get sample ID from $r1 vs $r2\n") 
             unless $r1 eq $r2;
         my $sid=$r1;
         if ($pdat) {
           #my ($r1, $r2)=($prevf, $f);
           #$r1=~s/^([^\.]+).+/$1/;
           #$r2=~s/^([^\.]+).+/$1/;
           #if ($r1 ne $r2) {
             # remove _r?[12] at the end to get the sample ID
           #}
           my $sd=$pd{$sid};
           die("Error: could not get pheno data for $r1\n") unless $sd;
           my $s=$$sd[2];
           if ($s=~m/^([MF])/i) {
             $s=uc($1);
           } elsif ($s=~m/^[XY]+$/) { 
             $s = ($s =~ m/Y/) ? 'M' : 'F';
           } else { die("Error: cannot recognize sex : $s\n"); }
           print ">$c $sid $s\n".join(" ", $d, $prevf, $f)."\n";
         } else {
           print ">$c $sid\n".join(" ", $d, $prevf, $f)."\n";
         }
        $c++;
      }
    }
    
  }
  $prevf=$f;
  $prevd=$d;
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

