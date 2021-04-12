#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
  sampSel2taskdb.pl [-r region] [-d dataset] [-p protocol]\
       samples.manifest rnaseq_phenodata.tab
  Outputs a pseudo-fasta file with tasks info based on the 
  given selection criteria
/;
umask 0002;
getopts('o:r:d:p:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
my $selreg=$Getopt::Std::opt_r;
my $seldset=$Getopt::Std::opt_d;
my $selproto=$Getopt::Std::opt_p;
my %spd; # rnum_flowcell => [dataset, protocol, dx, sex, race, age, brain_region]
my %rnum2sn; # rnum => [ list of rnum_flowcell, .. ]
my ($fman, $fpd)=@ARGV;
die("${usage}Error: no phenodata file provided\n") unless $fpd;
die("${usage}Error: cannot locate samples manifest file ($fman)!\n") unless -f $fman;
die("${usage}Error: cannot locate phenodata file ($fpd)!\n") unless -f $fpd;

# -- load phenodata
open(PD, $fpd) || die("Error opening $fpd !\n");
my $n=0;
while (<PD>) {
  next if m/^dataset\t/;
  chomp;
  #tr/\n\r//d; #if chomp fails for DOS file :
  my ($ds, $rnum, $sid, $reg, $proto, $rin, $pid, $dx, $sex, $race, $age, $pmi, $bf)=split(/\t/);
  my @b=split(/\,/, $bf);
  #next if ($b[0]=~m/_accepted_hits[\.\w]+bam$/);
  #next unless @b>1;
  my ($b1)=($b[0]=~m{([^/]+)$}); #remove path
  $b[0]=$b1;
  $b1=~s/_accepted_hits[\.\w]+bam$//;
  #choose only the first 2 _-delimited fields if there
  my @f=split(/_/, $b1);
  my $bsid= $f[0];$bsid.='_'.$f[1] if @f>1;
  #next if $bsid eq $sid;
  #print "b1='$b1', sid='$sid' : ".join(" | ", @b)."\n";
  #-- debug only:
  #last if $n++>5;
  my $sd=[ $ds, $proto, $dx, $sex, $race, $age, $reg ];
  foreach my $bn (@b) {
     $bn=~s/_accepted_hits[\.\w]+bam$//;
     my @p=split(/_/, $bn);
     my $bbn= $p[0];$bbn.='_'.$p[1] if @p>1;
     $spd{$bbn}=$sd;
  }
}
close(PD); 

## now load and check the samples.manifest
my $drnum=''; # full_path/rnum 
my @sids; #set of sample_ids (rnum_flowCell)
my (@r1, @r2); #full paths to all r1 and r2 files for $drnum
open(SM, $fman) || die("Error opening $fman !\n");
my $tc=0; #task counter

sub flushDb {
  #print STDERR "entering flushDb() with scalar(\@r1)=".scalar(@r1)."\n";
  return if @r1==0;
  $tc++;
  if (scalar(@r1)!=scalar(@sids)) {
   die("Error: number of sample_ids not matchin number of @r1 ($drnum)!\n");
  }
  my $s_id=$sids[0];
  my $sd=$spd{$s_id};
  unless ($sd) {
    print STDERR "WARNING: could not retrieve data for sample id $s_id ($drnum)!\n";
    return;
  }
  my ($rnum)=($s_id=~m/^(R\d+)/);
  my ($ds, $proto, $dx, $sex, $race, $age, $reg)=@$sd;
  return if ($selreg && $reg ne $selreg);
  return if ($seldset && $ds ne $seldset);
  return if ($selproto && $proto ne $selproto);
  my $fp=$drnum;
  $fp=~s{/[^/]+$}{};
  print ">$tc $fp $rnum $proto $reg $dx $sex $race $age\n";
  for (my $i=0;$i<@r1;$i++) {
   print ((@r2>0) ? "$r1[$i] $r2[$i]\n" : "$r1[$i]\n");
  }
}

while(<SM>) {
  next if m/^#/;
  chomp;
  my @t=split(/\t/);
  my $fp=$t[0];
  $fp=~s{/[^/]+$}{};
  my ($fn)=($t[0]=~m{([^/]+)$});
  my ($rnum)=($fn=~m/^(R\d+)/);
  if ($drnum ne $fp.'/'.$rnum) {
     if ($drnum && @r1>0) {
        #print "flushing ($drnum)\n";
        flushDb();
        @r1=();@r2=();@sids=();
     }
     $drnum=$fp.'/'.$rnum;
  }
  
  my ($sid, $f2)=(@t>3)? ($t[4], $t[2]) : ($t[2], '');
  #print "drnum=$drnum : rnum=$rnum, sid=$sid\n";
  die("Error: could not parse RNum from $fn :\n$_\n") unless $rnum;
  if (index($fn, $sid)!=0) { die("Error: sample id ($sid) not matching file prefix:\n$_\n")}
  push(@r1, $fn);
  push(@sids, $sid);
  if ($f2) {
     my ($fn2)=($f2=~m{([^/]+)$});
     if (index($fn2, $sid)!=0) { die("Error: sample id ($sid) not matching 2nd read file prefix:\n$_\n")}
     push(@r2, $fn2);
  }
  #last if $n++>5;
}
if ($drnum && @r1>0) {
        flushDb();
}
close(SM);

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

