#!/usr/bin/perl
use strict;
use Getopt::Std;
my $usage=q{Usage:
gff_mod_attrs.pl [-o <outpuf.gff3>] [-d attr_list] [-f id_attrs.tab] \
  [ -r attr_list] <input.gff> 
 
Modify the GFF attributes for GFF records.
Options:
  -d delete attributes specified in comma-delimited list attr_list
  -f for GFF records with IDs specified in the 1st column of the
     tab delimited file id_attrs.tab, add the attributes as given in
     the 2nd column of the file
  -r reorder the attributes according to their order in comma-delimited
     list attr_list
};
umask 0002;
getopts('r:f:d:o:') || die($usage."\n");
die("${usage}Error: only one input GFF file expected. Use '-' to specify a GFF stream at stdin\n")
  if (@ARGV!=1);

my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
}

shift(@ARGV) if $ARGV[0] eq '-' || $ARGV[0] eq 'stdin';

my ($ftab, $adel, $aorder)=($Getopt::Std::opt_f, $Getopt::Std::opt_d, $Getopt::Std::opt_r);
my %fattrs;
if ($ftab) {
  open(IN, $ftab) || die ("Error opening $ftab. $!\n");
  while (<IN>) {
    next if m/^#/;
    chomp;
    my ($id, $attrs)=split(/\t/);
    my %ah;
    my @alst;
    my @astr=split(/\s*\;\s*/, $attrs);
    if (@astr>0) {
      $astr[0]=~s/^\s+//;
      $astr[-1]=~s/\s+$//;
      foreach my $avstr (@astr) {
          my ($an, $av)=split(/\s*=\s*/, $avstr);
          $an=~s/^\s+//;$av=~s/\s+$//;
          push(@alst, $an);
          $ah{$an}=$av;
      }
    }
    $fattrs{$id}=[[@alst], \%ah];
  }
  close(IN);
}
my %dlist;
if ($adel) {
  $adel=~tr/,;:/ /;
  $adel=~tr/ / /s;
  %dlist=map { $_=>1 } (split(/\s+/, $adel));
}
my @olist;
my %ordh;
if ($aorder) {
  $aorder=~tr/,;:/ /;
  $aorder=~tr/ / /s;
  @olist=split(/\s+/, $aorder);
  @ordh{@olist}=();
}
while (<>) {
  #just store all ids with their attrs
  my $srcline=$_;
  chomp;
  my ($chr, $track, $f, $fstart, $fend, $fscore, $strand, $phase, $ats)=split(/\t/);
  next if m/^#/ || !$ats;
  my ($id)=($ats=~m/\bID=([^;]+)/);
  my %ah;
  my @astr=split(/\s*\;\s*/, $ats);
  if (@astr>0) {
    $astr[0]=~s/^\s+//;
    $astr[-1]=~s/\s+$//;
  }
  my @alst; #list of attributes stored, in order they were found
  foreach my $avstr (@astr) {
    my ($an, $av)=split(/\s*=\s*/, $avstr);
    $an=~s/^\s+//;$av=~s/\s+$//;
    next if $adel && exists($dlist{$an});
    push(@alst, $an);
    $ah{$an}=$av;
  }
  my $fd;
  if ($id && $ftab && ($fd=$fattrs{$id})) {
      my @fas=@{$$fd[0]};
      foreach my $fa (@fas) {
        my $fv=${$$fd[1]}{$fa};
        push(@alst, $fa) unless exists($ah{$fa});
        $ah{$fa}=$fv;
      }
  }
  my @olst=@alst; #output list of attributes
  if ($aorder) {
    @olst=grep { !exists($ah{$_}) } @olist;
    @alst=grep { !exists($ordh{$_}) } @alst;
    push(@olst, @alst);
  }
  print join("\t", $chr, $track, $f, $fstart, $fend, $fscore, $strand, $phase)."\t";
  my @oa;
  foreach my $a (@olst) {
    next if $adel && exists($dlist{$a});
    my $v=$ah{$a};
    die("Error: value for attribute $a not found!\n$srcline\n") unless $v;
    push(@oa, "$a=$v") if $v;
  }
  print join(';',@oa)."\n";
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
}
