#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage=q/Usage: 
  trmap_flt.pl [-t] [-s] [-m] 'codes' in.ovltab 
Filters trmap's output by the overlap codes given in the 'codes' string and
other criteria.

Example:
  trmap_flt.pl '=c' < trmap.ovltab > trmap.match_or_contained.ovltab

This example will output only the overlaps with class codes '=' (match)
or 'c' (contained).

Options:

 -s     only show overlaps with single-exon query transcripts
 -m     only show overlaps with multi-exon query transcripts
 -S     only show overlaps with single-exon reference transcripts
 -M     only show overlaps with multi-exon reference transcripts
 -t     output a simple 3 column table format (tab delimited): 
        query_id, ovl_code, reference_id
/;

die($usage."\n") if (@ARGV==0 || $ARGV[0] eq '-h' || 
 $ARGV[0] eq '--help');

getopts('smSMt') || die "$usage\n";

my $tblout=$Getopt::Std::opt_t;
my $rset=$Getopt::Std::opt_S;
my $rmet=$Getopt::Std::opt_M;
my $qset=$Getopt::Std::opt_s;
my $qmet=$Getopt::Std::opt_m;

my $codes=shift(@ARGV);

my %c= map { $_=>1 } (split(//,$codes));
my ($h, $hpr, $qid);
my $skip;
while (<>) {
 my $l=$_;
 if (m/^>/) {
   # header : query transcript data
   $skip=0;
   ($h, $hpr) = ($_, 1);
   if ($tblout) {
     ($qid)=(m/^>(\S+)/);
   }
   if ($qset || $qmet) {
     chomp;
     my ($x)=(m/([\d\-\,]+)$/);
     my $qc=scalar(split(/\,/,$x));
     if ($qset) {
       $skip=($qc>1);
     } else { #multi-exon only
       $skip=($qc==1);
     }
   }
 }
 else {
   next if $skip;
   #overlapped reference info
   my @d=split(/\t/);
   if (exists($c{$d[0]})) { #code match
     if ($rset || $rmet) {
       chomp;
       my ($x)=(m/([\d\-\,]+)$/);
       my $rc=scalar(split(/\,/,$x));
       if ($rset) {
         next if $rc>1;
       } else { #multi-exon only
         next if $rc==1;
       }
     }
     if ($tblout) {
        print join("\t",$qid, $d[0], $d[5])."\n";
        next;
     }
     if ($hpr) {
       print $h;
       $hpr=0;
     }
     print $l;
   }
 }
}
