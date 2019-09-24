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
          qry_id  qryExonCount ovl_code ref_id refExonCount [CDSmatch]
        (the CDSmatch column is only present if trmap was run with --show-cds
        and it will be T if CDS coordinates match or F otherwise
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
my ($h, $hpr, $qid, $qxn, $qcd, $qcstart, $qcstop);
my $skip;
while (<>) {
 my $l=$_;
 if (m/^>/) {
   # header line
   # query transcript data: tid, loc, strand, exons, [CDS:start:end]
   $skip=0;
   ($h, $hpr) = ($_, 1);
   chomp;
   my ($qloc, $qs, $qx);
   ($qid, $qloc, $qs, $qx, $qcd)=split(/\s+/);
   $qid=substr($qid,1);
   if ($qcd && $tblout) {
     ($qcstart)=($qcd=~m/^[CDS:]*(\d+)/);
     ($qcstop)=($qcd=~m/(\d+)$/);
     ($qcstart, $qcstop)=($qcstop, $qcstart) if ($qs eq '-');
   }
   $qxn=scalar(split(/\,/,$qx));
   if ($qset || $qmet) {
     chomp;
     if ($qset) {
       $skip=($qxn>1);
     } else { #multi-exon only
       $skip=($qxn==1);
     }
   }
 }
 else { #ref hit line:
   #  ovlCode, chr, strand, start, end, rID, exons, [CDS:start:end]
   next if $skip;
   #overlapped reference info
   chomp;
   my ($oc, $chr, $rs, $r0, $r1, $rid, $rx, $rcd)=split(/\t/);
   if (exists($c{$oc})) { #ovl code match
   my $rxn=scalar(split(/\,/,$rx));
   my ($rcstart, $rcstop);
   if ($rset || $rmet) {
      if ($rset) {
        next if $rxn>1;
      } else { #multi-exon only
        next if $rxn==1;
      }
   }
   if ($tblout) {
      my $cm='.';
      if ($rcd && $qcd) {
        $cm='same_CDS' if $rcd eq $qcd;
        if ($cm eq '.') {
           ($rcstart)=($rcd=~m/[CDS:]*(\d+)/);
           ($rcstop)=($rcd=~m/(\d+)$/);
           ($rcstart, $rcstop)=($rcstop, $rcstart) if ($rs eq '-');
        }
        $cm='same_START' if ($rcstart && $rcstart==$qcstart);
        if ($rcstop && $rcstop==$qcstop) {
           $cm= (length($cm)>1) ? $cm.',same_STOP' : 'same_STOP';
        }
      }
      print join("\t",$qid, $qxn, $oc, $rid, $rxn, $cm)."\n";
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
