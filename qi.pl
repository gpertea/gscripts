#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
  qi.pl [-u <user>] [-F]
  
  Show SGE jobs termination data from qacct -j output
  
  Options:
    -u <user>  specify the username to show jobs info 
               (default is the current username)
    -F        shows also failed jobs (default: omitted)
/;
umask 0002;
getopts('Fhu:o:') || die($usage."\n");
die($usage."\n") if $Getopt::Std::opt_h;
my $outfile=$Getopt::Std::opt_o;
my $incFailed=$Getopt::Std::opt_F;

my $user=$Getopt::Std::opt_u || $ENV{USER};

my ($owner, $jobnam, $jobnum, $taskid, $slots, $failed, $exit, 
    $start, $end, $wallclock, $maxrss, $maxvmem, $categ);
open(P, 'qacct -j |') || die("Error opening pipe!\n");
while (<P>) {
 chomp;
 if (m/^=====/) {
  if ($owner && index($owner, $user)==0) {
   $exit=~tr/ / /s;
   printf('>%s job:%d task:%d slots:%d %s - %s'."\n".
      ' CMD: %s |%s exit_status: %s'."\n".' Req: %s'."\n".
      ' Res: wall:.%s maxrss:%s maxvmem:%s'."\n",
      $owner, $jobnum,  $taskid, $slots, $start, $end,
      $jobnam, $failed>0 ?'FAILED':'', $exit, $categ, 
      $wallclock, $maxrss, $maxvmem ) if ($incFailed || $failed>0)
  }
  ($owner, $failed, $exit)=('',0,0);
  next;
 }
 if (m/^owner\s+(\S+)/) {
   $owner=$1; next; }
 if (m/^jobname\s+(\S.+)/) {
   $jobnam=$1; next; }
 if (m/^jobnumber\s+(\S.+)/) {
   $jobnum=$1; next; }
 if (m/^taskid\s+(\S.+)/) {
   $taskid=$1; next; }
 if (m/^slots\s+(\S.+)/) {
   $slots=$1; next; }
 if (m/^failed\s+(\S.+)/) {
   $failed=$1; next; }
 if (m/^exit_status\s+(\S.+)/) {
   $exit=$1; next; }
 if (m/^start_time\s+(\S.+)/) {
   $start=$1; next; }
 if (m/^end_time\s+(\S.+)/) {
   $end=$1; next; }
 if (m/^ru_wallclock\s+(\S.+)/) {
   $wallclock=$1; next; }
 if (m/^ru_maxrss\s+(\S.+)/) {
   $maxrss=$1; next; }
 if (m/^maxvmem\s+(\S.+)/) {
   $maxvmem=$1; next; }
 if (m/^category\s+(\S.+)/) {
   $categ=$1; next; }
}

close(P);

if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --



# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

