#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 gff_get.pl [-T|-I|-C|-E] [-m <minlen>] [-K] <gff_data_stream..>
 Extract feature coordinates from given GFF data stream.
 Output this tabulated fomat:
 
 <chr> <feature> <start> <end> <strand> <transcript_info>
 
 By default the program extracts only exon features from all transcripts
 found, unless one of the following options are used:
  -T : extract transcript start,end coordinates
  -I : extract intron coordinates
  -C : extract CDS coordinates
  -E : extract exon features (default)
 
 Options:
  -X try to exclude microRNAs and pseudogenes (NCBI annotation  only)
  -K only extract features from transcripts on "core" chromosomes (chr[\dXY]+)
  -m minium feature length to extract (default 4)
  -U only output unique genomic intervals (e.g. no same exon will be shown
     twice)

/;
umask 0002;
getopts('UXETICKm:o:') || die($usage."\n");
my $minlen=$Getopt::Std::opt_m || 4;
my $coreOnly=$Getopt::Std::opt_K; #core chromosomes only
my $getIntrons=$Getopt::Std::opt_I;
my $getCDS=$Getopt::Std::opt_C;
my $getExons=$Getopt::Std::opt_E;
my $uniqSegs=$Getopt::Std::opt_U;
my $getTranscripts=$Getopt::Std::opt_T;
$getExons=1 unless $getTranscripts || $getIntrons || $getCDS;
my $no_miRNA=$Getopt::Std::opt_X;
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# ---------------------  0       1       2           3      4      5       6      7         8
my %banned; #hash with parent IDs to discard
my %recs; # recID =>  [ chr,  strand, feat_type,  gname, tdescr, fstart, fend, [@exons], [@cds] ]

my %xgenes; 
my @tlst;
my @exd; #exons for current model
while (<>) {
   next if m/^\s*#/;
   chomp;
   my ($chr, $track, $f, $fstart, $fend, $fscore, $strand, $frame, $lnum)=split(/\t/);
   next unless $fstart>1 && $lnum;
   #next if $f =~m/gene$/i; # skipping any 'gene' features, unconditionally?
   next if $coreOnly && $chr!~m/^chr[\dXY]+$/;
   my ($gname, $tdescr);
   my $gff3_ID;
   my $gff3_Parent;
   ($fstart, $fend)=($fend, $fstart) if $fend<$fstart;
   ($gff3_ID)=($lnum=~m/\bID=([^;]+)/);
   ($gff3_Parent)=($lnum=~m/\bParent=([^;]+)/);
   my $is_gff3=($gff3_ID || $gff3_Parent);
   my $is_transcript=($f =~ m/RNA$/i || $f=~m/transcript$/i);
   my $is_exon=($f =~ m/exon$/i || $f =~ m/utr$/i || $f eq 'CDS');
   my $is_gene=($f =~ m/gene$/i);
   next unless $is_exon || $is_transcript || $is_gene;
   my $ID;
   if ($is_gff3) {
      $gff3_ID=~tr/"//d; #"
      $gff3_Parent=~tr/"//d; #"
      if ($gff3_Parent && !$is_exon ) {
        $gff3_Parent=''; #don't care about parents of transcripts or genes
      }
      $ID = $is_exon ? $gff3_Parent : $gff3_ID;
      next unless $ID; #should not happen
   }
   else { # GTF
    ($ID)=($lnum=~m/transcript_id\s+\"([^"]+)/);
    next unless $ID; #should not happen
   }
   my $class='';
   $class='class=lncRNA;' if $lnum=~m/class=lncRNA/i;
   my $pseudo=($f=~m/pseudo/i || $lnum=~m/pseudo=true/i || $lnum=~m/pseudogene/i);
   $class.='class=pseudogene;' if $pseudo;
   my $miRNA=($lnum=~m/~\bmiRNA/ || $lnum=~m/\bmicroRNA/);
   $class.='class=miRNA;' if $miRNA;
   if ($no_miRNA && ($miRNA || $pseudo)) {
       $banned{"$chr|$ID"}=1 if !$is_exon;
       next;
   }
   next if ($no_miRNA && exists($banned{"$chr|$ID"}));
   if ($is_gene) {
      if ($lnum=~m/\bgene_name[\s=]+"?([^;"]+)/i) {
         $gname=$1;
      } elsif ($lnum=~m/\bgene\w*[\s=]+"?([^;"]+)/i) {
         $gname=$1;
      }
      $xgenes{"$chr|$ID"}=[$class, $gname];
      #print STDERR "storing xgenes{$chr|$ID}=[$class, $gname]\n";
      next;
   }
   if ($is_transcript) {
         die("Error: duplicate feature $ID on $chr\n") if (exists($recs{"$chr|$ID"}));
         # try to parse the description, if any
         $tdescr='';
         my $tdn='';
         $gname='';
         if ($lnum=~m/\b(descr|description|tophit|info|product)\s*=\s*"?([^;"]+)/i) {
           $tdn=lc($1);
           $tdescr=$2;
         }
         elsif ($lnum=~m/(\w*name)\s*=\s*"?([^;"]+)/i) {
           $tdn=lc($1);
           $tdescr=$2;
         }
         if ($lnum=~m/\bgene_name[\s=]+"?([^;"]+)/i) {
           $gname=$1;
         } elsif ($lnum=~m/\bgene\w*[\s=]+"?([^;"]+)/i) {
           $gname=$1;
         }
         # elsif ($lnum=~m/Name\s*=\s*"?([^;"]+)/i) {
         #  $gname=$1;
         #}
         $tdescr='' if ($tdescr eq $gname);
         $gname='' if $gname eq $ID;
         $tdescr=$tdn.'='.$tdescr if $tdescr;
         if ($class) {
          $tdescr = $tdescr ? $class.$tdescr : $class;
         }
         my $recID="$chr|$ID";
         push(@tlst, $recID);
         $recs{$recID} = [$chr, $strand, $f, $gname, $tdescr, $fstart, $fend, [], [] ];
         next;
   }
   # -------------- exon/CDS line here:
   my $recID=$ID;
   if (!$gname && $lnum=~m/gene_name[= ]+(['"\:\w\.\|\-]+)/) {
      $gname=$1;
      $gname=~tr/"//d; #"
      }
   $tdescr='' if index($recID, $tdescr)>=0;
   $gname='' if index($recID, $gname)>=0;
   #$curtag=$chr.$strand;
   $recID=$chr.'|'.$recID;
   my $ld = $recs{$recID};
   if ($ld) { #existing entry
     my $i=($f eq 'CDS') ? 8 : 7;
     my ($lstart, $lend)=($$ld[5], $$ld[6]);
     $$ld[5]=$fstart if $fstart<$lstart;
     $$ld[6]=$fend if $fend>$lend;
     push(@{$$ld[$i]}, [$fstart, $fend, $fscore]);
     }
    else { # first time storing this transcript/gene
     my $gd=$xgenes{$recID};
     if ($gd) {
        ($class, $gname)=@$gd;
        if ($class) {
          $tdescr = $tdescr ? $class.$tdescr : $class;
        }
     }
     push(@tlst, $recID);
     $recs{$recID} = ($f eq 'CDS') ? 
           [$chr, $strand, $f, $gname, $tdescr, $fstart, $fend, [], [[$fstart, $fend, $fscore]] ] :
           [$chr, $strand, $f, $gname, $tdescr, $fstart, $fend, [[$fstart, $fend, $fscore]], [] ] ;
     }
} #while <>

writeModels();

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************
sub writeModels {
 my $num=keys(%recs);
 my %uexons;
 my %uintrons;
 my %utrans;
 return if $num==0;
 if ($num!=scalar(@tlst)) {
   die("Error: number of hash keys not matching list size!\n");
 }
 #while ( my ($l, $ld)=each(%recs) ) {
 foreach my $l (@tlst) {
   my $ld=$recs{$l} || die ("Error: locus $l found in list but not in hash!\n");
   my ($chr, $strand, $ftype, $gname, $descr, $lstart, $lend, $er, $cr) = @$ld;
   #my ($tstart,$tend)=($lstart, $lend);
   my @ex;
   if (@$er<1 && @$cr>0) { 
     @ex = sort { $main::a->[0] <=> $main::b->[0] } @$cr;
     }
    else {
      if (@$er>0) {
        @ex = sort { $main::a->[0] <=> $main::b->[0] } @$er ;
      }
      else {
        @ex = ([$lstart, $lend]);
      }
     }
   my ($tstart, $tend) = ($ex[0]->[0], $ex[-1]->[1]);
   my $ann=$descr if $descr;
   $ann.=';' if $ann && $gname;
   $ann.="gene=$gname" if $gname;
   my $t_id=substr($l, length($chr)+1);
   if ($getTranscripts) {
    #use min-max exon coordinates
     print join("\t", $chr, 'transcript', $tstart, $tend, $strand, $t_id, $ann)."\n"
      if !$uniqSegs || (++$utrans{$chr.'|'.$tstart.'-'.$tend}) == 1;
   }
   if ($getExons) {
    foreach my $x (@ex) {
       print join("\t", $chr, 'exon', $x->[0], $x->[1], $strand, $t_id, $ann)."\n"
          if !$uniqSegs || (++$uexons{$chr.'|'.$x->[0].'-'.$x->[1]}) == 1;
    }
   }
   if ($getCDS && @$cr>0) {
    foreach my $x (@$cr) {
       print join("\t", $chr, 'CDS', $x->[0], $x->[1], $strand, $t_id, $ann)."\n"
               if !$uniqSegs || (++$uexons{$chr.'|'.$x->[0].'-'.$x->[1]}) == 1;
    }
   }
   if ($getIntrons) {
     next if @ex<2;
     my @in; #introns to be printed stored here
     #my $icount=0; #intron output counter
     for (my $i=1;$i<@ex;$i++) {
       my ($istart, $iend)=($ex[$i-1]->[1]+1,$ex[$i]->[0]-1);
       next unless ($iend-$istart+1>=$minlen);
       #push(@in, [$istart, $iend]);
       #$icount++;
       print join("\t", $chr, 'intron', $istart, $iend, $strand, $t_id, $ann)."\n"
         if !$uniqSegs || (++$uintrons{$chr.'|'.$istart.'-'.$iend}) == 1;
       }
    }
  } #for each stored transcript
}

