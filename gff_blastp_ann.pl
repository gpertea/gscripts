#!/usr/bin/perl
use strict;
use Getopt::Std;
my $usage=q{Usage:
gff_blastp_ann.pl [-o <outpuf.gff3>] <models.gff> <prot_hits.ptab.fa.cidx>
 
Create a gff file with the same entries like <models.gff> but with
an annotation attribute ("product=...") added to the last field, based
on the best protein similarity hit found in <prot_hits.tab>

The <prot_hits.ptab.fa> should be FASTA-ized from blastp tab output 
sorted by transcripts ID and bit-score (best hits first, for each transcripts ID).

blastp -outfmt expected:
"qseqid qlen qstart qend sseqid slen sstart send ppos bitscore evalue frames salltitles"
Sorted like this:
 sort -k1,1 -k10,10nr -k11,11g 
then FASTA-ized with a command like this:
 perl -pe '@t=split(/\t/);print ">$t[0]\n" if $t[0] ne $p; $p=$t[0]'
..and finally use cdbfasta to create the .cidx file for this output.
};
umask 0002;
getopts('D:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
#my $ptabfile=$Getopt::Std::opt_p;
my $gff=shift(@ARGV) || die($usage."Error: missing argumemts!\n");
my $ptabf=shift(@ARGV) || die($usage."Error: no protein hits file given!\n");
die($usage."Error: cannot open protein hits file $ptabf!\n") unless -f $ptabf;
#my $pfh;
#if ($ptabf eq '-') {
#  $pfh = \*STDIN;
#}
#else { open($pfh, $ptabf) || die("Error opening file $ptabf !\n"); }

#-------- local array holdin uniformative annotation patterns:
my @uninformative=(
'\bunknown\b',
'\bhypothetical\b',
'[A-Z]+\d+ cDNA sequence',
'cDNA sequence [A-Z]+\d+',
'uncharacterized protein',
'unnamed protein product',
'open reading frame',
'\borf\b',
#'\bputative\b',
#'\bhomologue\b',
#'\bsimilar to',
'^expressed sequence \S+$',
'\bHA\d{4}\b',
'\bDKFZP\S+\b',
'PROTEIN FOR MGC:\d+',
'PROTEIN FOR IMAGE:\d+',
'\bR\d{5}\_\d\b',
'\bPRO\d{4}\b',
# 'KIAA\d+ GENE PRODUCT',
#'KIAA\d+ PROTEIN',
#'\bKIAA\d+\b',
'\bHSPC\d+\b',
# HSPC\d+ PROTEIN
#'\bC\d+ORF\d+\b',
'FLJ\d+ PROTEIN',
'\bDJ\d+[A-Z]\d+(\.\d+)*',
'NOVEL PROTEIN',
'CG\d+ PROTEIN',
'CG\d+ GENE PRODUCT',
'^\s*CG\d+\s*$',
'CGI\-\d+ PROTEIN',
'CGI\-\d+',
'CDNA:? FLJ\d+ FIS, CLONE \w+',
'BA\d+[A-Z]\d+[A-Z]?\.\d(\.\d)?',
#'\bRIKEN CDNA .{10} GENE\b',
'\bRIKEN.+?CDNA\b',
'MRNA, COMPLETE CDS, CLONE:\d+(\+\d[A-Z])?\-\d+',
'MRNA, COMPLETE CDS, CLONE:SMAP\d+\-\w+',
'BRAIN CDNA, CLONE MNCB-\d+',
'.{10}RIK PROTEIN',
'^MY\d{3}\s*$',
'MY\d{3} PROTEIN^',
'^probable\b',
'BRAIN MY\d{3}$',
'NPD\d{3} PROTEIN',
'[A-Z]\d{2}[A-Z0-9]+\.\d+ PROTEIN',
'WUGSC:H_\w+\.\w+ PROTEIN',
#'DNA SEGMENT, CHR [0-9XY]+, WAYNE STATE UNIVERSITY \d+, EXPRESSED',
#'DNA SEGMENT, CHR [0-9XY]+, KL MOHLKE \d+',
#'DNA SEGMENT, CHR [0-9XY]+, BAYLOR \d+',
'\bDNA SEGMENT\b',
'PROTEIN HSPC\d+',
#'HYPOTHETICAL [\.\d]+\s*KDA PROTEIN \S+ IN CHROMOSOME \S+',
'EG:[0-9A-Z\.]+ PROTEIN',
'GENOMIC DNA, CHROMOSOME \d+, P1 CLONE:\S+',
'[^,]+, RIKEN FULL-LENGTH ENRICHED LIBRARY, CLONE:.{10}, FULL INSERT SEQUENCE',
'ZK\d+\.\d+ PROTEIN',
'\bEST \w+',
'B2 ELEMENT'
);

open(GFF, $gff) || die ("Error: cannot open $gff\n");

#expects the models to be given in the proper order
# 'CDS' exons are the only ones used;
# the first 'mRNA' entry is going to be annotated
my @linebuf; #current gff gene/mRNA line buffer
#my ($gchr, $gstrand, @m_cds, $m_cdslen) ;  #current predicted mRNA info
my $curmodel; #current transcript ID
#                        0     1    2       3      4        5         6
my @genes;  # list of [$gID, $chr, $start, $end, $strand, $gffline, [@tIDs] ]
my %hgenes; # $gID => ^^^
#                       0     1      2       3      4        5         6
my @ts;     # list of [$tID, $chr, $start, $end, $strand, $gffline, [@sub-gfflines] ]
my %hts;    # $tID => ^^^
if ($outfile) {
  open(OUTF, '>'.$outfile) || die ("Error creating $outfile\n");
  select OUTF;
}

my $gffline=1;
while ($gffline) {
 $gffline=<GFF>;
 next if (defined($gffline) && ($gffline=~m/^\s*#/ || $gffline=~m/^\s+$/ || length($gffline)<4));
 chomp($gffline);
 my ($chr, $v, $f, $fstart, $fend, $fscore, $strand, $frame, $info)=split(/\t/, $gffline);
 next unless $info;
 my ($id)=($info=~m/ID=([^;]+)/);
 my ($parent)=($info=~m/Parent=([^;]+)/);
 if ($f eq 'gene') {
     if (!$id) {
       print STDERR "Warning: gene feature has no ID?\n$gffline\n";
       next;
     }
     my $gd=$hgenes{$id};
     print STDERR "Warning: gene $id already encountered?!\n" if $gd;
     $gd=[$id, $chr, $fstart, $fend, $strand, $gffline, []];
     $hgenes{$id}=$gd;
     push(@genes, $gd);
 } elsif ($f eq 'mRNA') { #new transcript starting
    if (!$id) {
      print STDERR "Warning: mRNA feature has no ID?\n$gffline\n";
      next;
    }
    my $td=$hts{$id};
    die("Error: transcript $id already encountered?!\n") if $td;
    if ($parent) {
      my $gd=$hgenes{$parent};
      die("Error: no gene parent defined for transcript $id!\n") unless $gd;
      push(@{$$gd[6]}, $id);
    }
    else { print STDERR "Warning: no parent gene defined for transcript $id ?!\n" }
    my ($acc, $descr, $org, $eval, $psim, $cov)=annModel($id);
    if ($cov>=25) {
       $gffline=addAnnotation($gffline, $acc, $descr, $eval, $psim, $cov);
    } else {
       if ($descr)
            { print STDERR "Warning: $id does not have any >=25% qcov protein matches (best match: $cov $acc $descr)\n" }
       else { print STDERR "Warning: $id lacks informative protein matches!\n" }
    }
    $td=[$id, $chr, $fstart, $fend, $strand, $gffline, []];
    $hts{$id}=$td;
    push(@ts, $td);
 }
 else { #must be a subfeature
   unless ($parent) {
      print STDERR "Warning: no Parent found for sub-feature $f ?\n$gffline\n";
      next;
   }
   my $td=$hts{$parent};
   die("Error: no transcript stored for $parent!\n") unless $td;
   push(@{$$td[6]}, $gffline);
 }
} #while <GFF>

close(GFF);
foreach my $gd (@genes) {
  print $$gd[5]."\n";
  foreach my $tid (@{$$gd[6]}) {
     my $td=$hts{$tid};
     die("Error: could not retrieve transcript data for $tid !\n") unless $td;
     print $$td[5]."\n";
     print join("\n", @{$$td[6]})."\n";
  }
}

if ($outfile) {
 select STDOUT;
 close(OUTF);
}

#----------------------------------

sub addAnnotation {
  my ($gffline, $acc, $descr, $eval, $psim, $cov)=@_;
  my $add="protSim=$descr;protSimAcc=$acc;protSimData=e-val:$eval,psim:$psim\%,qcov:$cov;";
  $gffline=~s/(ID=[^;]+;)/$1$add/ unless ($gffline=~s/(Parent=[^;]+;)/$1$add/);
  return $gffline;
}

sub annModel {
 my ($tid)=@_;
 open(PGET, "cdbyank -a '$tid' $ptabf  |") ||
  die("Error opening pipe from: cdbyank -a '$tid' $ptabf\n");
 my ($ann, $acov, $besthit, $bestcov);
 my @anndata; # list of [$acc, $descr, $org, $cov]
 while (<PGET>) {
  if (m/^>(\S+)\s*/) {
      die("Error: cdbyank did not retrieve the record for $tid ?\n") unless $tid eq $1;
      next;
  }
  chomp;
  my ($qid, $qlen, $qstart, $qend, $pid, $plen, $pstart, $pend, $psim, $bitscore, $ev, $frames, $allt)=split(/\t/);
  unless ($allt) {
     print STDERR "Warning: no protein description found for match:\n$_\n";
     next;
  }
  $allt=~s/(\s+)\>(\w+)/$1>~$2/g;
  $allt=">$pid ".$allt;
  ##split descriptions if needed:
  my @pdesc=split(/\s\>~/, $allt);
  if (@pdesc>1) {
     for (my $i=0;$i<@pdesc;$i++) {
       $pdesc[$i]='>'.$pdesc[$i] if $i>0;
       #$pdesc[$i].=']' if $i<$#pdesc;
     }
  }
  ($qstart, $qend)=($qend, $qstart) if $qend<$qstart;
  my $cov=sprintf('%d%', (($qend-$qstart+1)*100)/$qlen);
  foreach my $d (@pdesc) {
     #reduce imbricated square brackets
     #$d=~s/(\[[^\[]*)\[/$1{/;$d=~s/\]([^\]]*\])/$1/;
     #$d=~s/(\[[^\[]*)\[/$1{/;$d=~s/\]([^\]]*\])/$1/;
     #my ($acc, $descr, $org)=($d=~m/^>(\S+)\s([^\[]+)\[([^\]]+)\]$/);
     ## TIL: perl cannot do regex parsing from right to left -- forget about anchors..
     ## so doing a reverse trick here to get the organism name within imbricated patterns
     $d=reverse($d);
     #get the species name out of the way after parsing it
     $d=~s/(\]([^\[\]]|(?R))*\[)//;
     my $org=reverse($1);
     $d=reverse($d);
     if ($org) {
        $org=~s/^\[\s*//;$org=~s/\s*\]$//;
     } #else {
       # print STDERR "Warning:can't parse organism properly: $d\n$allt\n" unless $org;
     #}
     my ($acc, $descr)=($d=~m/^>(\S+)\s+(.+)/);
     $acc=~s/\|+$//;
     $descr=~s/\s+$//;
     if ($descr=~m/\sFull=\s*([^;]+)/) {
        $descr=$1;
     } else {
        $descr=$1 if ($descr=~m/^([^;]+);/);
     }
     die("Error:can't parse description properly (org=$org): $d\n$allt\n") unless $descr;
     push(@anndata, [$acc, $descr, $org, $ev, $psim, $cov]) if isInformative($descr);
  }
 } #while <PGET>
 close(PGET);
 #return getBestHits($hits, $mrna);
 return (@{$anndata[0]}) if @anndata>0;
 return ('','','','', '', '', 0);
}

sub ovlTest {
 my ($a1, $a2, $b1, $b2)=@_;
 #return ($a1<$b1) ? $b1<$a2 : $a1<$b2;
 return ($a1 <= $b2 && $b1 <= $a2);
}

#===============================================
# bool isInformative($description) 
# expects only the descripts - not the accession
#===============================================
sub isInformative {
 local $_=$_[0];
 s/^\s+//g;s/\s+$//g;
 return 0 if length($_)<2;
 foreach my $pat (@uninformative) {
   if (m/$pat/i) {
     #&flog("uninformative by /$pat/i : '$_'") if ($debug);
     return 0;
     }
   }
return 1;
}
