#!/usr/bin/perl
use strict;
use Getopt::Std;
use IPC::Open2;

my $usage = q/Usage:
  cds2rna_blat.pl tmappings.fa refprot.fa.cidx > tmappings_wCDS.tlf
  
  tmappings.fa must be a transcripts sequence file 
                 generated with gffread -W -w
  refprot.fa.cidx must be the cdb indexed file with reference proteins
                  corresponding to the transcripts in tmappings.fa
  Output will be a CDS annotated compact gff3 (tlf)
/;
umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
# --
die($usage."\n") unless @ARGV==2;
my ($tfa, $pepdb)=@ARGV;
die("${usage}File $tfa not found!\n") unless -f $tfa;
die("${usage}File $pepdb not found!\n") unless -f $pepdb;
pop(@ARGV);
my ($tid, $defline, $tseq);
my %tdescr; # tid => description
my %tm; #a transcript can have multiple mappings (e.g. gmap)
        #we may have multiple blat alignments
        # tid=>[ [tseqid, chr, strand, start, end, [exons], bscore, CDstart, CDend, bstatus], ...]
        #          0       1      2     3      4       5       6       7      8       9
my @stopCodons=qw(TAA TAG TAR TGA TRA);
my %stops;
@stops{@stopCodons}=('*') x @stopCodons;
mkdir('/dev/shm/gpertea');
my $ftmp="/dev/shm/gpertea/tpep.$$.fa";
my $counter=0;
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
       $tdescr{$tid}=$pdescr if $pdescr;
     }
     else {
       die("Error getting defline from $ftmp\n");
     }
     while (<PEP>) {
       chomp;
       $pseq.=$_;
     }
     $plen=length($pseq);
     close(PEP);
    } #reading protein
    #--
    my $blatcmd="blat -t=dnax -q=prot -out=psl -noHead stdin $ftmp stdout";
    my $prcid=open2(\*BL_OUT, \*DNA_IN, $blatcmd) || die("Error: open2() failed $!\n");
    print DNA_IN ">$tseqid\n$trawseq\n";
    close(DNA_IN); #flush input to blat command
    #now we can read the output
    {
     local $/="\n";
     # blat PSL format
     # #match #mism #repm #Ns #qgaps #qgapb #tgaps #tgapb strand q_id q_len q_start q_end t_name t_len t_start t_end #blocks blkSizes qStarts tStarts
     #  0       1     2    3    4      5      6       7      8     9    10    11     12     13    14      15    16     17       18       19      21
     while (<BL_OUT>) {
        chomp;
        my @psl=split(/\t/);
        my $revmap=($psl[8] eq '+-');
        my ($cds_start, $cds_end); #0-based CDS start, end offsets on the transcript sequence
        my $bscore = (3 * ($psl[0] + $psl[2])) - (3 * $psl[1]) - $psl[4] - $psl[6];
        my $bstatus='gmap';
        my ($bCDstart, $bCDend) = $revmap ? ($psl[14]-$psl[16], $psl[14]-$psl[15]) : ($psl[15], $psl[16]); #0-based PSL start coordinate for the CDS
        print STDERR ">analyzing $tseqid mapping: $chr:$tstart-$tend|$tstrand CDS:$CDstart-$CDend: blat: $bCDstart, $bCDend\n";
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
        my ($seq, $ofs)=$revmap ? (reverseComplement($tseq), $psl[15]) : ($tseq, $bCDstart);
        my @cods=unpack('(A3)*', substr($seq, $ofs));
        my $ci=0;
        foreach my $c (@cods) {
          $ci++;
          last if exists($stops{$c});
        }
        #stop found at $ofs+$ci*3
        $ci=scalar(@cods)-1 if $ci==@cods;
        my $adj_end=$ofs+$ci*3;
        $adj_end=$psl[14]-$adj_end if $revmap;
        if ($adj_end!=$cds_end) {
          print STDERR "\t$tseqid mapping adjusted ORF end from $cds_end to $adj_end\n";
          $cds_end=$adj_end;
          $bscore=3*($cds_end-$cds_start);
          $bstatus.='_adj';
        }
        if ($cds_end-$cds_start<10 && $plen>10) {
           print STDERR "Warning: $tseqid mapping has a very short ORF ($cds_start-$cds_end) for protein of length $plen\n";
           $bstatus='.tooshort';
        }
       my $td=$tm{$tid};
       if (!$td) {
         $td=[];
         $tm{$tid}=$td;
       }
       #adjust CDS offsets by original transcript strand:
       ($cds_start, $cds_end)=($tlen-$cds_end, $tlen-$cds_start) if $tstrand eq '-';
       push(@{$td}, [$tseqid, $chr, $tstrand, $tstart, $tend, [@ex], $bscore, $cds_start, $cds_end, $bstatus]);
     }#for each BLAT mapping
    } #reading blat results
    close(BL_OUT);
    waitpid($prcid,0);
 #-- for now --
 #last if ($counter>5); #for testing now
 } #while FASTA dna records
}
unlink($ftmp);

#now pick the best mapping and show the results:
while (my ($tid, $td) = each %tm) {
  my @srtd=sort { $b->[6] <=> $a->[6] } @$td;
  #first entry should be the largest score!
  my $descr=$tdescr{$tid};
  #cstart and $cend are 0-based, already adjusted for strand
  my ($tseqid, $chr, $strand, $start, $end, $exs, $bscore, $cstart, $cend, $bstatus)=@{$srtd[0]};
  my $tacc=0; #accumulating length
  my @cds; #CDS segments to populate
  my ($CDstart, $CDend);
  $cend--; #adjust CDend calculation
  foreach my $e (@$exs) {
    my $tinc=$$e[1]-$$e[0]+1;
    if ($cstart>=$tacc && $cstart-$tacc<$tinc) {
      $CDstart=$$e[0]+$cstart-$tacc;
      push(@cds, [$CDstart, $$e[1]]);
      $tacc+=$tinc;
      next;
    }
    if ($cend>=$tacc && $cend-$tacc<$tinc) {
       $CDend=$$e[0]+$cend-$tacc;
       push(@cds, [$$e[0], $CDend]);
       last;
    }
    push(@cds, [@$e]) if $CDstart;
    $tacc+=$tinc;
  }
 my $attrs="ID=$tid;mstatus=$bstatus";
 $attrs.=";exonCount=".scalar(@$exs);
 $attrs.=";exons=".join(',', (map { $$_[0].'-'.$$_[1] } @$exs));
 $attrs.=";CDS=$CDstart:$CDend";
 $attrs.=";descr=$descr" if $descr;
 print join("\t", $chr, 'remap', 'transcript', $start, $end, '.', $strand, '.', $attrs)."\n";
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
}

#************ Subroutines **************
sub reverseComplement {
  my $s=reverse($_[0]);
  $s =~ tr/AaCcTtGgUuMmRrWwSsYyKkVvHhDdBb/TtGgAaCcAaKkYyWwSsRrMmBbDdHhVv/;
  return $s;
 }

