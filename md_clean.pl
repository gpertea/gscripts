#!/usr/bin/perl
use strict;
use Getopt::Std;
#use FindBin;use lib $FindBin::Bin;
use IPC::Open3;
use Symbol 'gensym';

my $usage = q/Usage:
  md_clean.pl [-C] [-W] [-P] [-h 2] input.md > output.txt
  Cleans up markdown documents as output by ChatGPT
  Options:
  -C     clean-up non-ascii characters
  -W     use pandoc to convert the output to mediawiki syntax
  -P     properly embed perplexity citations in the text
  -h <n> force top heading level to <n>
/;
## pandoc markdown.md -f markdown -t mediawiki -o output.wiki.txt

umask 0002;
getopts('oCWPh:') || die($usage."\n");
my $hmin=$Getopt::Std::opt_h || 2; # minimum heading level
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
die($usage) unless @ARGV>0;
# --
my $incb=0; ## in-code-block flag
my $incit=0; ## in the citation section of perplexity
my $doc; # full document, needed for multi-line fixes
my @pcit; # perplexity citations (from the end of the document)
my $addh=0; ## how many # to add/subtract to reach $hmmin
my $firstHeading=1;
while(<>) {
  if (s/^\s*```(\S*)/```$1/) {
    $incb=(length($1)>0);
  }
  if (m/^\s*Citations:\s*$/) {
    s/$/<br\/>/;
    $doc.=$_;
    $_=<>;
    if (m/^\s*\[\d+\]\s+http\S+$/) { $incit=1 }
  }
  if (!$incb) {
   if ($incit) {
     if (s/^\s*\[(\d+)\]\s+(http\S+)\s*$/$1 $2<br\/>\n/) {
       @pcit[$1]=$2;
     }
   }# incit
   else { # not in the end citation block
     s/\[(\d+)\]/ [$1](prplx~Cit~$1)/g; # to be updated later
     if (m/^(#+) /) {
       my $h=$1;
       if ($firstHeading) {
         $addh=$hmin-length($h);
         print STDERR "setting addh as $hmin-".length($h)." = $addh\n";
         $firstHeading=0;
       }
       s/\*+//g;
       if ($addh) {
         my $newh=($addh<0) ? substr($h, 0, $addh) : $h.=('#' x $addh);
         s/^#+/$newh/x;
       }
     }
   }
  } # not in a code block
  $doc.=$_;
}
## full document processing:
$doc=~s/\[\s+/[/sg;$doc=~s/\s+\]/]/sg;
$doc=~s/\(prplx~Cit~(\d+)\)/($pcit[$1])/xg;
if ($Getopt::Std::opt_W) {
  my $err = gensym;  # for capturing STDERR
  # Open a bidirectional pipe with pandoc.
  my $pid = open3(my $in, my $out, $err, "pandoc", "-", "-f", "markdown", "-t", "mediawiki");
  # Print the whole document to the process's stdin.
  print $in $doc;
  close $in; # Important: signal EOF to the pandoc command.

  # Read and process the output line by line.
  while (my $line = <$out>) {
      #chomp($line);
      # Process each line as needed.
      print "$line";
  }
  # Optionally, handle any errors.
  while (my $line = <$err>) {
      chomp($line);
      warn "Error encountered: $line\n";
  }
  # Wait for the process to finish.
  waitpid($pid, 0);
} else {
 print $doc;
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
}

#************ Subroutines **************



