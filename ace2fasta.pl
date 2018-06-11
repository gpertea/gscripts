#!/usr/bin/perl
use strict;
use Getopt::Std;

my $usage = q/Usage:
 ace2fasta.pl [-o <contigs_fasta>] [-p <prefix>] [-c <components_file>] <ACE>

 In the output FASTA contigs will be renamed with the prefix 
 <prefix>_ctg# if -p <prefix> option was provided.
/;
umask 0002;
getopts('o:c:p:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
$outfile='' if $outfile eq '-';
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}
my $prefix=$Getopt::Std::opt_p;
my $lnkfile=$Getopt::Std::opt_c;
if ($lnkfile) {
 open(COMP, '>'.$lnkfile) || die("Error creating file $lnkfile!\n");
 }

my $ctg; #current contig name
my $ctgseq; #current contig sequence
my $ctglen; #current contig seq length
my $numseqs;
my @ctgcomp; #list of component sequence names
my %seqs; # seq => [ seqlen, strand, seqL, seqR, asmL, asmR ]
my %ctgnames;
while (<>) {
 if (m/^CO\s+(\S+)\s+(\d+)\s+(\d+)/) {
   my ($c, $l, $n)=($1, $2, $3);
   #if ($c ne $ctg) {
   &writeCtg() if $ctg;
   undef %seqs;
   undef @ctgcomp;
   my $ctgdup=++$ctgnames{$c};
   if ($prefix) {
     my ($cn)=($c=~m/(\d+)$/);
     die("Error: could not parse contig # from contig name: $c\n") unless $cn>0;
     $c=$prefix.'_ctg'.sprintf('%04d',$cn);
   }
   #$c='ASM_'.$c;
   $c.='.'.($ctgdup-1) if ($ctgdup>1);
   ($ctg, $ctglen, $numseqs)=($c, $l, $n);
   #  }
   $ctgseq='';
   my $seqline;
   while (<>) {
     ($seqline)=(/^(\S+)/);
     last unless $seqline;
     $ctgseq.=$seqline;
     }
   } #contig start line
  elsif (m/^AF (\S+) ([UC]) ([\-\d]+)/) {
    my $strand=($2 eq 'U')?'+':'-';
    $seqs{$1}=[0, $strand, 0, 0, $3, 0]; #the untrimmed asmL position for now
    }
  elsif (m/^RD (\S+) (\d+) (\d+) (\d+)/) {
    my ($seqname, $seqlen)=($1, $2);
    my $seqd=$seqs{$seqname};
    push(@ctgcomp, $seqname);
    die("Error at ACE parsing: no sequence found for RD $seqname\n")
      unless $seqd;
    $seqd->[0]=$seqlen;
    my $rdseq='';
    my $seqline;
    while ( defined($_=<>) && (($seqline)=(/(\S+)/))) {
      $rdseq.=$seqline;
      }
    do { 
      $_=<>;
      } until m/^QA (\d+) (\d+)/ || !defined($_);
    die("Couldn't find the QA entry for RD $seqname!\n") unless defined($_);
    ($seqd->[2], $seqd->[3]) = ($1, $2);
    my ($trimL, $trimR)=($seqd->[2]-1, $seqlen-$seqd->[3]);
    $seqd->[5]=$seqd->[4]+$seqlen-$trimR-1;
    $seqd->[4]+=$trimL;
    }
}

writeCtg() if $ctg;

if ($outfile) {
 select(STDOUT);
 close(OUTF);
}

close(COMP) if $lnkfile;

sub writeCtg {
 #write FASTA formatted
 print ">$ctg $numseqs\n";
 $ctgseq=~tr/-*\n\r//d; #removing gaps, if any!
 print join("\n", (unpack('(A72)*',$ctgseq)))."\n";
 if ($lnkfile) {
   print COMP ">$ctg $numseqs $ctglen\n";
   foreach my $seqname (@ctgcomp) {
     my $sd=$seqs{$seqname} || 
        die("Error: no seqdata for component $seqname!\n");
     print COMP $seqname.' '.join(' ',@$sd)."\n";
     }
   }

}
