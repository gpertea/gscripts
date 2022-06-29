#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

## assume rnum data file has the following columns:
###  sample_id,rnum,dataset,protocol,region,brnum,dx,sex,race,age,dropped
##       0       1     2       3        4     5    6  7    8   9   10

my $usage = q/Usage:
  Show subjects, Dx, Region summaries for a list of RNums 
  Input can be a list of RNums or file names with a R\\d+_ format 
  given at stdin or with the -f option.  
  Options:
      -l   just list the unique rnums found
      -m   output the metadata table for the unique rnums given
      -f   file with RNums to consider (default: stdin)
      -d   rnaseq_samples.tab file (default: ~\/rnaseq_samples.tab
/;
umask 0002;
getopts('o:f:d:lm') || die($usage."\n");
my $list=$Getopt::Std::opt_l;
my $meta=$Getopt::Std::opt_m;
my $fdb=$Getopt::Std::opt_d || $ENV{HOME}.'/rnaseq_samples.tab';
die("Error locating data file $fdb\nCheck -d option?\n") unless -f $fdb;
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my $rnfile=$Getopt::Std::opt_f;
die("Error locating input file $rnfile\n") if $rnfile && ! -f $rnfile;
die($usage) if (! $rnfile && -t STDIN ); #no input given 

my %urn; # rnums are keys
my @rnums; # list of rnums in the order given

my $in; # read either from $rnfile or STDIN
if ($rnfile) {
  open($in, '<', $rnfile) || die("Error opening $rnfile for reading!\n");
}
else {
  $in=*STDIN;
}
my $tr=0;
while($_=<$in>) {
 my ($r)=(m/\b(R\d\d+)/);
 if (!$r) {
    print STDERR "Warning: no RNum could be parsed from line:\n$r\n";
    next;
 }
 $tr++;
 if (!exists($urn{$r})) {
   $urn{$r}=1;
   push(@rnums,$r);
 }
}
close($in) if $rnfile;
if ($list) {
 print join("\n", @rnums)."\n";
 exit(0);
}
print STDERR scalar(@rnums)." unique RNums given ($tr parsed).\n";


my %rdb; # rnum => [ 
#         [ sample_id, dataset,protocol,region,brnum,dx,sex,race,age,dropped ],
#         [ sample_id, dataset,protocol,region,brnum,dx,sex,race,age,dropped ],
#        ]

# read the database file
open(FDB, $fdb) || die("Error opening database file $fdb\n");
if ($meta) {
 print join("\t", qw(sample_id rnum dataset protocol region brnum dx sex race age dropped))."\n"
}

while(<FDB>) {
 next if m/rnum/i; #skip header
 chomp;
 my @t=split(/[\t\,]/); #sample_id,rnum,dataset,protocol,region,brnum,dx,sex,race,age,dropped
 if ($meta && exists($urn{$t[1]})) {
   print "$_\n";
   next;
 }
 my $sid=shift(@t);
 my $rnum=$t[0];
 $t[0]=$sid; 
 my $rd=$rdb{$rnum};
 if (!$rd) { #no entry for this RNum
  $rdb{$rnum}=[ [@t] ];
  next;
 }
 push(@$rd, [ @t ])
}
close(FDB);
exit(0) if $meta;
my (%tds, %tdx, %treg, $numdrop, @nf); # counts per dataset, Dx, region, and list of missing (not found)
# sample_id, dataset, protocol, region, brnum, dx, sex, race, age, dropped
#     0        1        2         3       4     5   6    7     8     9
foreach my $rn (@rnums) {
 my $rds=$rdb{$rn};
 unless($rds) { push(@nf, $rn); next }
 foreach my $rd (@$rds) {
    $tds{$$rd[1]}++;
    $treg{$$rd[3]}++;
    $tdx{$$rd[5]}++;
    $numdrop++ if $$rd[9];
 }
}

## show the stats, list the missing to stdout
my @ds=map { "$_ (".$tds{$_}.')' } keys(%tds);
my @dx=map { "$_ (".$tdx{$_}.')' } keys(%tdx);
my @reg=map { "$_ (".$treg{$_}.')' } keys(%treg);

print " Datasets: ".join(', ',@ds)."\n";
print "       Dx: ".join(', ',@dx)."\n";
print "  Regions: ".join(', ',@reg)."\n";

print "  Dropped: $numdrop\n";
print "Not found: ",scalar(@nf), "  ",(@nf>0 ? join(',',@nf)."\n" : "\n");

# -- 
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

