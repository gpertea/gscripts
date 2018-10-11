#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 gff_liftover.pl [-C] ref2query.paf input.gff
 If input.gff is '-', the gff is expected at stdin
 Options: 
   -C  only CDS segments are converted (anything else is discarded)
/;
umask 0002;
getopts('Co:') || die($usage."\n");
die("${usage}Error: not enough parameters given!\n") unless @ARGV==2;
my $outfile=$Getopt::Std::opt_o;
my $paf=shift(@ARGV);
die("${usage}Error: cannot find PAF file $paf !\n") unless -f $paf;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --

my $cdsonly=$Getopt::Std::opt_C;
my $fbed='gff_liftover.bed';
open(BED, ">$fbed") || die("Error creating file gff_liftover.bed\n");
my @exons;
my ($lstrand, $ltid, $lchr); #last (previous) values
while (<>) {
 next if m/^\s*#/;
 my $line=$_;
 chomp;
 next if m/^\s*$/;
 my ($chr, $track, $f, $fstart, $fend, $fscore, 
     $strand, $frame, $lnum)=split(/\t/);
 my ($ID)=($lnum=~m/\bID=([^;]+)/);
 my ($Parent)=($lnum=~m/\bParent=([^;]+)/);
 die("Error: no ID or Parent found for:\n$line\n") 
    unless ($ID || $Parent);
 #next unless ($ID || $Parent);
 my $is_transcript=($f =~ m/RNA$/i || $f=~m/transcript$/i);
 my $is_exon=($f =~ m/exon$/i || $f =~ m/utr$/i || $f eq 'CDS');
 my $is_gene=($f =~ m/gene$/i);
 next unless $is_exon;
 next if $cdsonly && $f ne 'CDS';
 next if $f eq 'CDS' && !$cdsonly; #convert only exons
 if ($Parent ne $ltid) {
   flushBED(); #uses globals
   ($lstrand, $ltid, $lchr)=($strand, $Parent, $chr);
   @exons=();
 }
 push(@exons, [$fstart-1, $fend]);
}
flushBED();
close(BED);
my $pcmd="paftools.js liftover $paf $fbed";
open(MAP, "$pcmd|") || die("Error opening pipe $pcmd|\n");
open(BED, $fbed) || die("Error opening $fbed for reading!\n");
my $bedline=<BED>;
my @exd; #list of [$start, $end] for exons of $ltid
my ($lnumexons, $lstatus);
while(<MAP>) {
  chomp;
  my ($mchr, $mstart, $mend, $oloc, $anum, $ostrand)=split(/\t/);
  my $oseg=$oloc;
  while ($oseg=~s/_t[53]$//) {}; #cleanup the mapped segment info
  my ($bchr, $bstart, $bend, $t_info, $bnumexons, $bstrand);
  my $found=0;
  while ($bedline) {
    chomp($bedline);
    ($bchr, $bstart, $bend, $t_info, $bnumexons, $bstrand)=split(/\t/, $bedline);
    if ($oseg eq $bchr.'_'.$bstart.'_'.$bend) {
      $found++;
      my ($tid, $exno)=split(/\~/, $t_info);
      if ($tid ne $ltid) { #change of $t
       writeMTranscript($ltid, $lnumexons, $lchr, $lstrand, \@exd);
       ($ltid, $lnumexons, $lchr, $lstrand)=($tid, $bnumexons, $mchr, $bstrand);
       @exd=();
      }
      if ($found==1) {
        my $status=($oseg ne $oloc)?'a':''; # a = ambiguous mapping
        push(@exd, [$mstart+1, $mend, $exno.$status]);
      } elsif ($found==2) {
        $exd[-1]->[2].='m';  #m=multiple mappings
      }
    }
    else { #different source segment
      last if $found;
    }
    $bedline=<BED>;
  } #while src BED lines to check for the current mapping
} #while MAPpings
writeMTranscript($ltid, $lnumexons, $lchr, $lstrand, \@exd);

close(MAP);
close(BED);

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************
sub flushBED{
 my $n=scalar(@exons);
 return if $n==0;
 my $en=0;
 foreach my $ed (@exons) {
   $en++;
   print BED join("\t", $lchr, @$ed, $ltid.'~'.$en, $n, $lstrand)."\n";
 }
}

sub writeMTranscript {
  my ($tid, $numexons, $chr, $strand, $exr)=@_;
  return if @$exr==0;
  my @ex=sort { $main::a->[0] <=> $main::b->[0] } @$exr;
  my ($tstart, $tend)=($ex[0]->[0], $ex[-1]->[1]);
  my @aex=grep { $_->[2]=~m/[am]/ } @ex;
  my $attrs="ID=$tid;exoncount=$numexons";
  $attrs.=';ambiguous=1' if @aex>0;
  $attrs.=';partial=1' if @ex!=$numexons;
  $attrs.=';ok=1' if @aex==0 && @ex==$numexons;
  print join("\t", $chr, 'liftover', 'transcript', $tstart, $tend, '.', $strand, '.', $attrs)."\n";
  my $f=$cdsonly?'CDS':'exon';
  foreach my $e (@ex) {
   $attrs="Parent=$tid;exoninfo=$$e[2]\n";
   print join("\t",  $chr, 'liftover', $f, $$e[0], $$e[1], '.', $strand, '.', $attrs);
  }
  @$exr=();
}
