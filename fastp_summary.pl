#!/usr/bin/perl
## quick parsing of basic adapter trimming summary from _fastp.json and _fastp.log output
## provide only the prefix, _fastp.json and _fastp.log are automatically added
use strict;
my $fpre=shift(@ARGV);
die("Error: file name prefix expected!\n")
  if !defined($fpre) || ($fpre=~m/^\-\-?h/) || ! (-f $fpre.'_fastp.json' && -f $fpre.'_fastp.log');
my ($bef, $aft, $fr, $ins, $adcut, # sections
    $bef_r1r, $bef_r1b, $bef_r2r,  $bef_r2b,  $bef_r1_avglen, $bef_r2_avglen, 
    $aft_r1r, $aft_r1b, $aft_r2r, $aft_r2b,  $aft_r1_avglen, $aft_r2_avglen, 
    $fr_passed, $fr_short, $fr_ins_peak, $fr_trim_r, $fr_trim_b);
open(L, $fpre.'_fastp.log') || die(" Error opening log file\n");
open(J, $fpre.'_fastp.json') || die(" Error opening json file\n");
## most info is obtained from the log file, we only need JSON for read1/2_mean_length
my %lv=( 'before.r1.reads'=>\$bef_r1r, 'before.r1.bases'=>\$bef_r1b,
         'before.r2.reads'=>\$bef_r2r, 'before.r2.bases'=>\$bef_r2b,
         'after.r1.reads'=>\$aft_r1r, 'after.r1.bases'=>\$aft_r1b,
         'after.r2.reads'=>\$aft_r2r, 'after.r2.bases'=>\$aft_r2b,
          );
while (<L>) {
  if ($fr) {
    if (m/Insert size peak/) {
      ($fr_ins_peak)=(m/(\d+)$/);
      last;
    }
    if (m/reads passed filter:/) {
      ($fr_passed)=(m/(\d+)$/);
      next
    }
    if (m/reads failed due to too short:/) {
      ($fr_short)=(m/(\d+)$/);
      next
    }
    if (m/reads with adapter trimmed:/) {
      ($fr_trim_r)=(m/(\d+)$/);
      next
    }
    if (m/bases trimmed due to adapters:/) {
      ($fr_trim_b)=(m/(\d+)$/);
      next
    }
    next;
  } elsif (m/^Filtering result:/) {
    $fr=1;
    next;
  }
  ## not in filtering results yet:
  for my $ba (('before', 'after')) {
    for my $rn ((1,2)) {
      if (m/Read$rn $ba filtering:/) {
        for my $rb (('reads', 'bases')) {
          $_=<L>;
          if (m/otal $rb:\s*(\d+)/) {
            my $v=$1;
            my $k="$ba.r$rn.$rb";
            ${$lv{$k}}=$v;
            my $kv=${$lv{$k}};
          } else { expErr("r$rn: total $rb"); }
        } # reads, bases
      }
    }
  }
}
close(L);
($bef, $aft)=(0,0);
while (<J>) {
   if ($aft) {
     if (m/"read1_mean_length":\s*(\d+)/) {
       $aft_r1_avglen=$1;
       next;
     }
     if (m/"read2_mean_length":\s*(\d+)/) { 
       $aft_r2_avglen=$1;
       last;  ## quit parsing
     }
     next;
   } elsif (m/"after_filtering":/) {
     $aft=1;
     next;
   }
   if ($bef) {
     if (m/"read1_mean_length":\s*(\d+)/) {
       $bef_r1_avglen=$1;
       next;
     }
     if (m/"read2_mean_length":\s*(\d+)/) { 
       $bef_r2_avglen=$1;
       next;
     }
   } elsif (m/"before_filtering":/) {
     $aft=0;
     $bef=1;
     next;
   }
}
close(J);

print join("\t", qw(r1_reads r1_bases r2_reads r2_bases));
print "\t".join("\t", qw(trim_r1_reads trim_r1_bases trim_r2_reads trim_r2_bases));
print "\t".join("\t", qw(trim_passed trim_too_short trimmed_reads trimmed_bases fastp_insert_peak));
print "\t".join("\t", qw(r1_avg_len r2_avg_len trim_r1_avg_len trim_r2_avg_len))."\n";
print join("\t", $bef_r1r, $bef_r1b, $bef_r2r, $bef_r2b);
print "\t".join("\t", $aft_r1r, $aft_r1b, $aft_r2r, $aft_r2b);
print "\t".join("\t", $fr_passed, $fr_short, $fr_trim_r, $fr_trim_b, $fr_ins_peak);
print "\t".join("\t", $bef_r1_avglen, $bef_r2_avglen, $aft_r1_avglen, $aft_r2_avglen)."\n";

sub expErr {
 die("Error: not getting expected log line with $_[0]\n");
}
