#!/usr/bin/perl
use strict;
use Getopt::Std;
use IPC::Open2;

my $usage = q/Usage:
  cds2rna_blat.pl [-M][-w outblat.psl] [-b inblat.psl] \
                  tmappings.fa refprot.fa.cidx > tmappings_wCDS.tlf
  
  tmappings.fa must be a transcripts sequence file 
                 generated with gffread -W -w
  refprot.fa.cidx must be the cdb indexed file with reference proteins
                  corresponding to the transcripts in tmappings.fa
  Output will be a CDS annotated compact gff3 (tlf)
/;
umask 0002;
getopts('w:b:MDo:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
my $debug=$Getopt::Std::opt_D;
my $wblat=$Getopt::Std::opt_w;
my $inblat=$Getopt::Std::opt_b;
my $mito_codons=$Getopt::Std::opt_M;
# --
die($usage."\n") unless @ARGV==2;
my ($tfa, $pepdb)=@ARGV;
die("${usage}File $tfa not found!\n") unless -f $tfa;
die("${usage}File $pepdb not found!\n") unless -f $pepdb;
pop(@ARGV);
my ($tid, $defline, $tseq);
my %tdescr; # tid => [prot_len, description]
my %tm; #a transcript can have multiple mappings (e.g. gmap)
        #we may have multiple blat alignments
        # tid=>[ [tseqid, chr, strand, start, end, [exons], bscore, CDstart, CDend, bstatus], ...]
        #          0       1      2     3      4       5       6       7      8       9
my @stop_codons=qw(TAA TAG TAR TGA TRA);
my @start_codons=qw(ATG);
my %stopCodons;
my %startCodons;
@stopCodons{@stop_codons}=('.') x @stop_codons;
@startCodons{@start_codons}=('M')  x @start_codons;
if ($mito_codons) {
  delete($stopCodons{'TGA'}); # TGA is W (Trp) in mitochondria
  @stopCodons{'AGA', 'AGG'}=('.', '.');
  $startCodons{'ATA'}='M';
}

my $totalprot=0;

mkdir('/dev/shm/gpertea');
my $ftmp="/dev/shm/gpertea/tpep.$$.fa";
my $counter=0;
my @Bnextaln; # used for getNextAln() from blat

if ($inblat) {
 die("Error: input blat file must be different from output blat file!\n") if $wblat eq $inblat;
 open(BLAT, $inblat) || die("Error opening psl file $inblat !\n");
}

if ($wblat) {
 open(WBLAT, ">$wblat") || die("Error: cannot create file $wblat !\n");
}
{
 local $/="\n>";
 while (<>) {
    s/^>//;
    chomp;
    $counter++;
    my ($tseqid, $ann, $trawseq)=(m/^(\S+)[ \t\x01]*([^\n]*)\n?(.*)/s);
    #my @nr=split(/\x01/, $ann, 2);
    # example: Os01t0100100-01.mrna1 CDS=382-2490 loc:NC_029256.1.R|6946-14778|+ \
    # exons:6946-7231,7317-7579,8320-8418,9420-9523,11099-11907,11991-12113,12195-12283,12371-12571,13173-13580,14067-14150,14237-14393,14467-14778 \
    # segs:1-286,287-549,550-648,649-752,753-1561,1562-1684,1685-1773,1774-1974,1975-2382,2383-2466,2467-2623,2624-2935
    my ($tloc)=($ann=~m/loc:(\S+)/);
    die("Error: no location could be parsed from $ann\n") unless $tloc;
    my @l=split(/\|/, $tloc);
    my ($tstrand, $tstart_end)=(pop(@l), pop(@l));
    my $chr=join('|', @l);
    my ($tstart, $tend)=split(/\-/,$tstart_end);
    my ($exons)=($ann=~m/exons:([\d\-\,]+)/);
    my @ex=map { [split(/\-/)] } (split(/\,/, $exons));
    my ($CDstart, $CDend);
    if ($ann=~m/\bCDS=(\d+)\-(\d+)/) {
      ($CDstart, $CDend)=($1, $2);
    }
    my $tseq = uc($trawseq);
    $tseq=~tr/\t \n\r//d;
    my $tlen=length($tseq);
    my $tid=$tseqid;
    $tid=~s/\.[a-z]+\d+$//;
    my $protcmd="cdbyank -a '$tid' $pepdb > $ftmp";
    system($protcmd) && die("Error running command: $protcmd\n");
    die("Error: protein file $ftmp has 0 size\n") unless (-s $ftmp);
    #--parse $ftmp to get the description, protein sequence etc.
    my ($pepid, $pdescr, $pseq, $plen);
    {
     local $/="\n";
     open(PEP, $ftmp) || die("Error opening $ftmp!\n");
     local $_=<PEP>;
     if (m/^>(\S+)/) {
       chomp;
       $pepid=$1;
       die("Error: protein id $pepid differs from query tid $tid!\n") if $pepid ne $tid;
       if (m/description:([^;]+)/) {
         $pdescr=$1;
         $pdescr=~s/\.? \(\S+\)$//;
       }
       else { $pdescr=''; }
     }
     else {
       die("Error getting defline from $ftmp\n");
     }
     while (<PEP>) {
       chomp;
       $pseq.=$_;
     }
     $plen=length($pseq);
     $totalprot+=$plen; # reference protein space
     $tdescr{$tid}=[$plen, $pdescr];
     close(PEP);
    } #reading protein
    #--
    my $prcid;
    unless ($inblat) {
      my $blatcmd="blat -t=dnax -q=prot -out=psl -noHead stdin $ftmp stdout";
      print STDERR "Running blat cmd:\n$blatcmd\n" if $debug;
      $prcid=open2(\*BLAT, \*DNA_IN, $blatcmd) || die("Error: open2() failed $!\n");
      print DNA_IN ">$tseqid\n$trawseq\n";
      close(DNA_IN); #flush input to blat command
    }
    #now we can read the output
    {
     local $/="\n";
     # blat PSL format
     # #match #mism #repm #Ns #qgaps #qgapb #tgaps #tgapb strand q_id q_len q_start q_end t_name t_len t_start t_end #blocks blkSizes qStarts tStarts
     #  0       1     2    3    4      5      6       7      8     9    10    11     12     13    14      15    16     17       18       19      21
     my @psl;
     while (getNextAln(\@psl, $tseqid)) {
       #last if ($inblat && m/^#.#/);
       print WBLAT join("\t",@psl)."\n" if $wblat;
       #chomp;
       #my @psl=split(/\t/);
       die("Error: expect PSL alignment for $tseqid but got $psl[13] instead !\n") if $tseqid ne $psl[13];
       #die("Error: unexpected blat strand specification: -+ for $tseqid\n") if ($psl[8] eq '-+');
       my $revmap=($psl[8] eq '+-');
       my ($cds_start, $cds_end); #0-based CDS start, end offsets on the transcript sequence
       my $bscore = (3 * ($psl[0] + $psl[2])) - (3 * $psl[1]) - $psl[4] - $psl[6];
       my $bstatus='gmap';
       my ($bCDstart, $bCDend) = $revmap ? ($psl[14]-$psl[16], $psl[14]-$psl[15]) : ($psl[15], $psl[16]); #0-based PSL start coordinate for the CDS
       #my ($bCDstart, $bCDend) = ($psl[15], $psl[16]);
       if ($debug) {
          print STDERR ">analyzing $tseqid mapping: $chr:$tstart-$tend|$tstrand CDS:$CDstart-$CDend: blat\[$psl[8]\] $bCDstart-$bCDend\n";
       }
       if ($CDstart && $bCDstart+1==$CDstart) {
         #agreeing with GMAP here, use its coordinates
         $cds_start=$bCDstart;
         $cds_end=$CDend;
         $bstatus='gmap';
       } else { #take blat's t_start coordinate and find the longest reading frame from there
         $cds_start=$bCDstart;
         $cds_end=$bCDend;
         $bstatus='blat';
       }
       my ($seq, $ofs)=$revmap ? (reverseComplement($tseq), $bCDstart) : ($tseq, $bCDstart);
       if ($debug) {
         print STDERR "psl $psl[8] alignment of qry($psl[11]-$psl[12]) to target($psl[15]-$psl[16])  ";
         print STDERR 'bCDstart='.($bCDstart+1),", bCDend=$bCDend\n";
         print STDERR "sequence at offset $ofs:\n".substr($seq, $ofs, 72)."\n";
       }
       my $partial='';
       my @cods=unpack('(A3)*', substr($seq, $ofs));
       my $ci=0;
       $partial.='5' if !($psl[11]==0 && exists($startCodons{$cods[0]}));
       foreach my $c (@cods) {
         last if exists($stopCodons{$c});
         $ci+=length($c);
       }
       my $bstop_adj=0;
       my $p_end=$cds_end;
       if ($ci+$ofs==$tlen) {
          $partial.='3';
          print STDERR "\tno stop codon found, cds_end set to tlen $tlen\n" if $debug;
          $cds_end=$tlen;
       } else { #stop codon found, include it in CDS
          $bstop_adj=3;
          $cds_end=$ofs+$ci+3;
       }
       if ($p_end!=$cds_end-$bstop_adj) {
         if ($debug) {
           print STDERR "\t$tseqid mapping adjusted CDS end from $p_end to $cds_end\n";
         }
         $bscore-=3*($p_end-$cds_end);
         $bstatus.='_adj';
       }
       if ($cds_end-$cds_start<10 && $plen>10) {
          #print STDERR "Warning: $tseqid mapping has a very short ORF ($cds_start-$cds_end) for protein of length $plen\n";
          $bstatus.='.shrt';
       }
       $bstatus.='|partial'.$partial if $partial;
       if ($debug) {
          print STDERR "\tCDS settled to $tseqid:".($cds_start+1)."-$cds_end".($revmap?'-':'')." status=$bstatus\n";
       }
       my $td=$tm{$tid};
       if (!$td) {
         $td=[];
         $tm{$tid}=$td;
       }
       #adjust CDS offsets by original transcript strand:
       my $adjstrand=$tstrand;
       if ($revmap && $bstatus=~m/^blat/) {
         $adjstrand = $tstrand eq '-' ? '+' : '-';
       }
       ($cds_start, $cds_end)=($tlen-$cds_end, $tlen-$cds_start) if $adjstrand eq '-';
       push(@{$td}, [$tseqid, $chr, $adjstrand, $tstart, $tend, [@ex], $bscore, $cds_start, $cds_end, $bstatus]);
     }#for each BLAT mapping
    } #reading blat results
    #print WBLAT "#.#\n";
    if (!$inblat) {
       close(BLAT);
       waitpid($prcid,0);
    }
 } #while FASTA transcript dna records for the mappings
}

if ($wblat) {
  close(WBLAT);
}
if ($inblat) {
 close(BLAT);
} else {
 unlink($ftmp) unless $debug;
}

#now pick the best mapping and show the results:
my $totalcds=0;
while (my ($tid, $td) = each %tm) {
  my @srtd=sort { $b->[6] <=> $a->[6] } @$td;
  #first entry should be the largest score!
  my $pd=$tdescr{$tid};
  die("Error: could not retrieve protein data for $tid\n") unless $pd;
  my ($plen, $pdescr)=@$pd;
  #cstart and $cend are 0-based, already adjusted for strand
  my ($tseqid, $chr, $strand, $start, $end, $exs, $bscore, $cstart, $cend, $bstatus)=@{$srtd[0]};
  my $prr=(($cend-$cstart)*100.00)/(3*$plen);
  $totalcds+=$cend-$cstart;
  $totalcds-=3 if ($bstatus!~m/\|partial.?3$/);
  my $tacc=0; #accumulating length
  my @cds; #CDS segments to populate
  my ($CDstart, $CDend);
  $cend--; #adjust CDend calculation
  foreach my $e (@$exs) {
    my ($cs, $ce)=@$e;
    my $tinc=$ce-$cs+1;
    my $last=0;
    if ($cstart>=$tacc && $cstart-$tacc<$tinc) {
      $CDstart=$$e[0]+$cstart-$tacc;
      $cs=$CDstart;
      #push(@cds, [$CDstart, $$e[1]]);
    }
    if ($cend>=$tacc && $cend-$tacc<$tinc) {
       $CDend=$$e[0]+$cend-$tacc;
       $ce=$CDend;
       #push(@cds, [$$e[0], $CDend]);
       #last;
    }
    push(@cds, [$cs, $ce]) if $CDstart;
    last if $CDend;
    $tacc+=$tinc;
  }
 my $attrs="ID=$tid;mstatus=$bstatus";
 $attrs.=";exonCount=".scalar(@$exs);
 $attrs.=";exons=".join(',', (map { $$_[0].'-'.$$_[1] } @$exs));
 $attrs.=";CDS=$CDstart:$CDend";
 if ($prr<50.0) {
    $attrs.= $prr<30.0 ? ';prr=low' : ';prr=poor';
 }
 $attrs.=";descr=$pdescr" if $pdescr;
 print join("\t", $chr, 'remap', 'transcript', $start, $end, '.', $strand, '.', $attrs)."\n";
}
printf '## Protein Recovery Rate (PRR) : %.2f'."\n", (($totalcds*100.00)/(3*$totalprot));

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
}

#************ Subroutines **************
sub getNextAln {
  #returns 1 if indeed 
  my ($rpsl, $t_id)=@_;
  if (@Bnextaln>0) { #pushed earlier because it's a new target id
    if ($Bnextaln[13] eq $t_id) {
        @$rpsl=@Bnextaln;
        @Bnextaln=();
        return 1;
    }
    return 0;
  }
  do {
    $_=<BLAT>;
  } until (!m/^#/);
  if (!$_) { #eof
    @$rpsl=();
    return 0;
  }
  my @r=split(/\t/);
  if ($r[13] eq $t_id) {
     @$rpsl=@r;
     return 1;
  }
  #new target id, keep it for later
  @Bnextaln=@r;
  @$rpsl=();
  return 0;
}

sub reverseComplement {
  my $s=reverse($_[0]);
  $s =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBb/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVv/;
  return $s;
 }
