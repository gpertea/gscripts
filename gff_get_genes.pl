#!/usr/bin/perl
use strict;
use Getopt::Std;
my $usage = q/gff_get_genes.pl geneID1,geneID2,... <input_GFF_stream>
 Let pass all the genes, along with their transcripts, whose IDs are 
 in the given comma-delimited list
/;


getopts('ho:') || die($usage."\n");
die($usage."\n") unless @ARGV;
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
 open(OUTFILE, '>'.$outfile) || die ("Error creating file $outfile\n");
 select(OUTFILE);
}
my $geneIDs=shift(@ARGV);
my @gids=split(/\,\s*/, $geneIDs);
my %tofetch; # IDs for "gene" features with OK biotype
foreach my $gid (@gids) {
 $tofetch{$gid}=[];
}

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
 if (exists($tofetch{$Parent})) {
   print STDOUT $line;
   $tofetch{$ID}=$Parent;
   next;
 }
 print $line if exists($tofetch{$ID});
}
