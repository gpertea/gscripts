#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q{Usage:
  cbrain_check_jobs.pl [-w tasks_to_rerun.cfa] dataset taskdb.cfa
  
  dataset is the base directory name where output was directed.
  It assumes the log files are alongside output files as:
    in ~/work/cbrain/[aln|strg]/dataset/sampleID/sampleID.log
  ..and that sampleID is the 3rd token in the taskdb.cfa header, e.g.:
      >5 /path/to/fastq_or_aln_dir MSSM_RNA_BP_PFC_10 ... 
      
  Option -w allows creation of a copy and renumbered list of tasks entries
  that were found having issues.
};
umask 0002;
getopts('o:w:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my $ds=shift(@ARGV) || die($usage." No dataset name provided\n");
my $tfa=shift(@ARGV) || die($usage." No taskdb file provided\n");
my $tdbout=$Getopt::Std::opt_w;
if ($tdbout) {
  die("Error: -w option cannot be the same with input taskdb!\n")
    if $tfa eq $tdbout;
 open(WT, '>'.$tdbout) || die("Error creating $tdbout\n");
}
open(TF, $tfa) || die(" Error opening taskdb file $tfa\n" );
my $wtid=0;
my $wrest=0;
while(<TF>) {
 my $line=$_;
 if (m/^>(\d+)/) {
   my $tid=$1; #task ID   
   $wrest=0;
   chomp;
   my @t=split();
   my $sid=$t[2];
   my $odir=$ENV{HOME}."/work/cbrain/aln/$ds/$sid";
   die("Error: no sample directory: $odir\n") unless -d $odir;
   my $fbam="$odir/$sid.bam";
   my $fcram="$odir/$sid.cram";
   my $nobam= (! -f $fbam);
   my $nocram= (! -f $fcram);
   my $bamsize= -s $fbam;
   my $cramsize= -s $fcram;
   my $flog="$odir/$sid.log";
   my $badcram=($cramsize<300000000); # should be >300M
   my $badbam=($bamsize>0 && $bamsize<300000000);
   my $cramckfail=0;
   my $bamckfail=0;
   if (!$badcram) { # check cram integrity
        $cramckfail=system("samtools quickcheck $fcram");
        $badcram=1 if ($cramckfail);
   }
   if ($bamsize>0 && !$badbam) {
        $bamckfail=system("samtools quickcheck $fbam");
        $badbam=1 if ($bamckfail);
   }
   my $errmsg='.';
   # also check logs for anything suspicious
   my $errs=`egrep -i 'erro|warn|fail|couldn|invalid|unable|cann|not found|broken|dump|fault' $flog`;
   if (length($errs)>2) {
     $errs=~s/[\n\r]+$//;
     $errmsg=join(";", (split(/[\n\r]+/, $errs) ));
   }
   if ($badcram || $badbam || length($errmsg)>1) {
     my $bst=$nobam ? 'nobam' : ($badbam ? 'badbam' : 'bamOK');
     my $cst=$nocram ? 'nocram' : ($badcram ? 'badcram' : 'cramOK');
     print join("\t",$tid, $cst, $bst, $errmsg)."\n";
     if ($tdbout) {
        $wtid++;
        $line=~s/^>\d+/>$wtid/;
        print WT $line;
        $wrest=1; #to print the other lines for this record, if any
     }
   }
 } elsif ($wrest) {
    print $line;
 }
} 
close(WT) if ($tdbout);
close(TF);

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

