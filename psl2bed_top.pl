#!/usr/bin/env perl

use strict;
use warnings;

=head1 NAME

This scripts converts a psl file into a bed file written by Dave Tang.

=head1 SYNOPSIS

Usage: psl2bed_top.pl <infile.psl> <outfile.bed>

=head1 DESCRIPTION

A psl file is the default output of blat. This script reads in a psl file,
keeps only the best scoring hit and outputs the entries in bed format. If
there are two or more equal best hits, they are all outputted. The score is
calculated as per http://genome.ucsc.edu/FAQ/FAQblat.html#blat4

psl_to_bed_best_score.pl created on Tue Nov  5 14:58:46 JST 2013.

=head1 AUTHOR

Please report bugs to <me@davetang.org>.

=cut

my $usage = "Usage: $0 <infile> <outfile>\n";
my $infile = shift or die $usage;

#store for psl results
my %psl_result = ();
#colour for bed file
my $red = '255,0,0';

#open psl file
open (IN, "<", $infile) || die "Could not open $infile: $!\n";
while (<IN>){
   chomp;
   #skip psl header lines and blanks
   next if (/^psL/ || /^match/ || /^\s+match/ || /^----/ || /^$/);
   my $current_line = $_;

   my ($matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts) = parse_psl($current_line);

   #Calculate score for this match
   #I know I shouldn't hardcode this in but I mostly work with DNA sequences
   #change sizeMul to 3 if you are working with amino acid sequences
   my $sizeMul = 1; #This is 1 for dna and 3 for protein
   my $score = ($sizeMul * ($matches + $repMatches)) - ($sizeMul * $misMatches) - $qNumInsert - $tNumInsert;

   #if query has been stored before
   if (exists $psl_result{$qName}){
      #go through previous score
      foreach my $prevScore (keys %{$psl_result{$qName}}){
         if ($score > $prevScore){
            delete $psl_result{$qName}{$prevScore};
            $psl_result{$qName}{$score}->['0'] = $current_line;
         }
         #if the score is equal to the previous score
         #store both entries
         elsif ($score == $prevScore){
            my $num = scalar @{$psl_result{$qName}{$score}};
            $psl_result{$qName}{$score}->[$num] = $current_line;
         }
      }
   }
   #if query is new
   else {
      $psl_result{$qName}{$score}->['0'] = $current_line;
   }

}
close(IN);

#report number of entries stored
my $entries = scalar(keys %psl_result);
warn("Stored $entries entries\n");

my $outfile = shift or die $usage;
open (OUT, ">", $outfile) || die "Could not open $outfile for writing: $!\n";

#go through the best scoring entries
foreach my $query (keys %psl_result){
   #should only have one score, the best score
   foreach my $score (keys %{$psl_result{$query}}){
      #go through entries with the best score
      for (my $i = 0; $i < scalar (@{$psl_result{$query}{$score}}); ++$i){
         my $result = $psl_result{$query}{$score}->[$i];
         my ($matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts) = parse_psl($result);
         my @splitBlockSizes = split(/,/,$blockSizes);
         my @split_target_block_start = split(/,/,$tStarts);

         #The tBlocks are given in the positive direction, regardless of mapping to the negative strand
         #tStarts + $blockSizes will be the exon
         #However, BED chromStarts offsets must be relative to chromStart, not absolute.
         #Subtract chromStart from each offset in chromStarts.
         my $relative_t_starts = '';
         foreach my $psl_t_starts (@split_target_block_start){
            my $sub = $psl_t_starts - $tStart;
            $relative_t_starts .= "$sub,";
         }
         $relative_t_starts =~ s/,$//;
         #BED  chrom,chromStart,chromEnd,name  ,score ,strand ,thickStart,thickEnd,itemRgb,blockCount ,blockSizes ,blockStarts
         print OUT join("\t", $tName, $tStart, $tEnd, $qName, $score, $strand, $tStart, $tEnd, $red, $blockCount, $blockSizes, $relative_t_starts),"\n";
      }
   }
}

close(OUT);

exit(0);

#I've included this subroutine just as a means to document the psl format
sub parse_psl {
   my ($line) = @_;
   my @result = split(/\t/,$line);
   my $matches     = $result[0];    # Number of bases that matches the query matched to the target
   my $misMatches  = $result[1];    # Number of bases that don't match
   my $repMatches  = $result[2];    # Number of bases that match but are part of repeats
   my $nCount      = $result[3];    # Number of 'N' bases
   my $qNumInsert  = $result[4];    # Number of inserts in query
   my $qBaseInsert = $result[5];    # Number of bases inserted in query
   my $tNumInsert  = $result[6];    # Number of inserts in target
   my $tBaseInsert = $result[7];    # Number of bases inserted in target
   my $strand      = $result[8];    # '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
   my $qName       = $result[9];    # Query sequence name
   my $qSize       = $result[10];   # Query sequence size
   my $qStart      = $result[11];   # Alignment start position in query
   my $qEnd        = $result[12];   # Alignment end position in query
   my $tName       = $result[13];   # Target sequence name
   $tName =~ s/\.fa$//;
   my $tSize       = $result[14];   # Target sequence size
   my $tStart      = $result[15];   # Alignment start position in target
   my $tEnd        = $result[16];   # Alignment end position in target
   my $blockCount  = $result[17];   # Number of blocks in the alignment (a block contains no gaps)
   my $blockSizes  = $result[18];   # Comma-separated list of sizes of each block
   my $qStarts     = $result[19];   # Comma-separated list of starting positions of each block in query
   my $tStarts     = $result[20];   # Comma-separated list of starting positions of each block in target
   return($matches,$misMatches,$repMatches,$nCount,$qNumInsert,$qBaseInsert,$tNumInsert,$tBaseInsert,$strand,$qName,$qSize,$qStart,$qEnd,$tName,$tSize,$tStart,$tEnd,$blockCount,$blockSizes,$qStarts,$tStarts);
}
