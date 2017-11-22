#!/usr/bin/perl
use strict;
use Getopt::Std;
my $usage = q/gb_gff_filter.pl [options] <input_GFF_stream>

Processes genomic annotation GFF files from NCBI
filtering out anything but the following gene_biotype genes:
   lncRNA
   protein_coding

/;


getopts('ho:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
 open(OUTFILE, '>'.$outfile) || die ("Error creating file $outfile\n");
 select(OUTFILE);
 }
my @OKbiotypes=('lncRNA', 'protein_coding');
my %biotypes;
@biotypes{@OKbiotypes}= (1) x scalar(@OKbiotypes);
my %OKgenes; # chr|ID for "gene" features with OK biotype
my %OKtranscripts; # chr|ID for children of OK genes
my ($numgenes, $numtranscripts);
while (<>) {
 next if m/^\s*#/;
 my $line=$_;
 chomp;
 my ($chr, $track, $f, $fstart, $fend, $fscore, 
     $strand, $frame, $lnum)=split(/\t/);
 my ($ID)=($lnum=~m/\bID=([^;]+)/);
 my ($Parent)=($lnum=~m/\bParent=([^;]+)/);
 die("Error: no ID or Parent found for:\n$line\n") 
    unless ($ID || $Parent);
 next unless ($ID || $Parent);
 my ($g_id, $t_id);
 if ($f eq 'gene') {
   my ($btype)=($lnum=~m/\bgene_biotype=([^;]+)/);
   die("Error: biotype not found for gene line:\n$line\n") unless $btype;
   #print STDERR "biotype=<$btype>\n";
   next unless ($btype && exists($biotypes{$btype}));
   $g_id=$chr.'|'.$ID;
   $OKgenes{$g_id}=[] unless exists($OKgenes{$g_id});
   print STDOUT $line;
   next;
 }
 if ($f=~m/transcript$/ || $f=~m/RNA$/) {
   $t_id=$chr.'|'.$ID;
   $g_id=$chr.'|'.$Parent;
   next unless (exists($OKgenes{$g_id}));
   print STDOUT $line;
   $OKtranscripts{$t_id}=[] unless exists($OKtranscripts{$t_id});
   next;
 }
 if ($f=~m/exon$/ || $f=~m/CDS$/) {
   $t_id=$chr.'|'.$Parent;
   if (exists($OKgenes{$t_id})) {
     print STDERR "Gene-transcript exception (exon/CDS parented by gene $t_id)\n";
     $OKtranscripts{$t_id}=[];
   }
   next unless (exists($OKtranscripts{$t_id}));
   print STDOUT $line;
 }
}
