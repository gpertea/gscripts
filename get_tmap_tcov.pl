#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Gather coverage for matching reference transcripts 
in given tmap files. Appends this info as a new column to the tinfo.tab file 
Usage:
 get_tcov.pl [-t tinfo.tab] gffcmp.gtf.tmap ..
/;
umask 0002;
getopts('t:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
my $tbl=$Getopt::Std::opt_t;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
# --
my @td; #list of [$fn, \%tcov]
my %allt; #hash with all ref transcript IDs in the tmaps parsed
foreach my $fn (@ARGV) {
  open(F, $fn) || die("Error opening file $fn!\n");
  my %tc;
  while (<F>) {
    next if m/^ref_gene_id/;
    my @t=split(/\t/);
    next unless $t[2] eq '=';
    $tc{$t[1]}=$t[8];
    $allt{$t[1]}=1;
  }
  $fn=~s/\.gtf\.tmap$//;
  $fn=~s/[^\/]+\///;
  $fn=~s/^cmp_//;
  $fn=~s/test/t/;
  $fn=~s/t.wsplnoise/tws/;
  $fn=~tr{./}{__};
  push(@td, [$fn, \%tc]);
  close(F);
}

if ($tbl) {
  open (T, $tbl) || die("Error opening $tbl\n");
  print join("\t", ('tid', 'num_exons', 'len', 'sim_tpm', 'read_count', 'sim_cov') )."\t";
  print join("\t", (map { $_->[0] } @td))."\n";
  while (<T>) {
    my @t=split(/\t/);
    chomp;
    print $_."\t";
    print join("\t", (map { sprintf('%.2f', ${$_->[1]}{$t[0]}) } @td))."\n";
  }
  close(T);
}
else {
  print "tid\t".join("tid\t", (map { $_->[0] } @td))."\n";
  foreach my $tid (keys %allt) {
     print "$tid\t".join("\t", (map { sprintf('%.2f', ${$_->[1]}{$tid}) } @td))."\n";
  }
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
}

#************ Subroutines **************

