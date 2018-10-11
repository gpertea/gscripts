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
my %bed; # tid~exon# => [chr, exon_start, exon_end, strand, num_exons]
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
#                                        0       1       2          3
my %mt; #remapped transcripts: $ltid->[ chr, strand,  exoncount, \@exd]
#where @exd = list of [$chr, $start, $end, $orig_eno, $flags] = remapped exons of $tid
#                        0       1     2        3        4
my %idc; # interval remapping count
while(<MAP>) {
  chomp;
  my ($mchr, $mstart, $mend, $idloc, $zero, $ostrand)=split(/\t/);
  $idloc=~s/~(\d+)\|/~$1\t/;
  my ($eid, $oloc)=split(/\t/, $idloc);
  my $oseg=$oloc;
  my $tm=0; #truncated mapping for this exon?
  while ($oseg=~s/_t[53]$//) { ++$tm; }; #cleanup the mapped segment info
  my $rbed=$bed{$eid};
  die("Error: could not find interval ID $eid in BED file!\n") unless $rbed;
  my ($bchr, $bstart, $bend, $bstrand, $bnumexons)=@$rbed;
  my ($tid, $eno)=split(/~/,$eid);
  my $md=$mt{$tid};
  my $found=($idc{$eid}++);
  my $eflags=0;
  $eflags|=2 if $tm;
  if ($found>1) { #this exon has multiple mappings
    die("Error: exon $eid found before but transcript not stored?!\n") unless $md;
    $eflags|=1;
    #update exising exon
    my $fex=0; #found exon?
    for my $ed (@{$md->[3]}) {
       if ($ed->[3]==$eno) {
         $fex=1;
         $ed->[4] |= $eflags;
         last;
       }
    }
    die("Error: exon $eid was not found for multiple mapping ($found)!\n") unless $fex;
  }
  else { #first time seeing this interval
    if ($md) {
      #update exons
      #check if exon mapping already stored -> should not be!
      for my $ed (@{$md->[3]}) {
        die("Error: exon $eid already found for first mapping!\n")
         if ($ed->[3]==$eno);
      }
      push(@{$md->[3]}, [$mchr, $mstart, $mend, $eno, $eflags]);
    } else {
      #create transcript entry in %mt
      $mt{$tid}=['', $bstrand, $bnumexons, [[$mchr, $mstart, $mend, $eno, $eflags]] ];
    }
  }
} #while MAPpings
while (my ($tid,$td)=each(%mt)) {
 printRemapped($tid, $td);
}

close(MAP);

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
   my $iid=$ltid.'~'.$en;
   $bed{$iid}=[$lchr, $$ed[0], $$ed[1], $lstrand, $n];
   print BED join("\t", $lchr, @$ed, $iid, $n, $lstrand)."\n";
 }
}

sub printRemapped {
  my ($tid, $td)=@_;
  my ($tchr, $tstrand, $numexons, $exr)=@$td;
  my @ex=sort { $main::a->[1] <=> $main::b->[1] } @$exr;
  my ($tstart, $tend)=($ex[0]->[1], $ex[-1]->[2]);
  my @tex=grep { $_->[4] & 2 } @ex; #any truncated exon mappings?
  my @mex=grep { $_->[4] & 1 } @ex; #any multiple exon mappings?
  $tchr=$ex[0]->[0];
  my @xex=grep { $_->[0] ne $tchr } @ex; #multi-chromosome mappings?
  my $attrs="ID=$tid;exoncount=$numexons;m_exoncount=".scalar(@ex);
  my $al=length($attrs);
  $attrs.=';truncated=1' if @tex>0;
  $attrs.=';ambiguous=1' if @mex>0;
  $attrs.=';partial=1' if @ex!=$numexons;
  $attrs.=';multi_chr=1' if @xex>0;
  $attrs.=';ok=1' if $al==length($attrs);
  print join("\t", $tchr, 'liftover', 'transcript', $tstart, $tend, '.', $tstrand, '.', $attrs)."\n";
  my $f=$cdsonly?'CDS':'exon';
  foreach my $e (@ex) {
   $attrs="Parent=$tid;src_exon_number=".$e->[3];
   $attrs.=';truncated=1' if ($e->[4] & 2);
   $attrs.=';multimap=1' if ($e->[4] & 1);
   $attrs.=';diff_chr=1' if ($e->[0] ne $tchr);
   print join("\t",  $e->[0], 'liftover', $f, $$e[1], $$e[2], '.', $tstrand, '.', $attrs)."\n";
  }
  @$exr=();
}
