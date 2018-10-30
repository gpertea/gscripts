#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
  gff_cp_attrs.pl [-a 'attr_lst'] [ -r 'remove_lst'] attrs_src.gff target.gff
  Copies the attributes from attrs_src.gff to the transcripts with the 
  same IDs in target.gff.
  
  Output a GFF which is the same with target.gff but with attributes from
  attrs_src.gff copied over for the matching ID.
  
  Option -a takes a comma, period or colon delimited list of attribute names
  which should be copied over from attrs_src.gff (if -a is not given, all 
  attributes will be copied, except ID,exons,CDS,CDSphase)
  Options -r takes a similar list of attributes but these attributes will be
  discarded from the final output.
  Option -t <trackname> can be used to override the GFF track name
  in the output.
/;

umask 0002;
getopts('t:r:a:o:') || die($usage."\n");

my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my %xcludeAttr= ( 'ID' => 1, 'exonCount'=> 1, 'exons'=> 1, 'CDS'=>1, CDSphase=>1 );
my %removeAttr;
my $gfftrack=$Getopt::Std::opt_t;
my $attrlist=$Getopt::Std::opt_a;
my $removeattrs=$Getopt::Std::opt_r;
my %attr;
if ($attrlist) {
 my @l=split(/[:\.,]/, $attrlist);
 map { $attr{$_}=1 } @l;
}
if ($removeattrs) {
 my @l=split(/[:\.,]/, $removeattrs);
 map { $removeAttr{$_}=1 } @l;
}
die("${usage}") if @ARGV!=2;
my $srcgff=shift(@ARGV);
my %sd; #ID => { attr => value, ... }
open(SRCIN, $srcgff)||die("Error: cannot open $srcgff! ($!)\n");

while(<SRCIN>) {
  #just store all ids with their attrs
  chomp;
  my @t=split(/\t/);
  next if m/^#/ && @t<9;
  my ($tid)=($t[8]=~m/\bID=([^;]+)/);
  my %ah;
  #my %ah= map { @a=split(/\s*=\s*/);$a[0]=>$a[1] } @astr;
  my @astr=split(/\s*\;\s*/, $t[8]);
  next unless $tid; #ignore
  foreach my $avstr (@astr) {
    my ($an, $av)=split(/\s*=\s*/, $avstr);
    if ($an eq 'ID') {
      die("Error: parsed ID ($tid) not matching attr value for ID ($av)\n")
         if $tid ne $av;
      next;
    }
    if ($attrlist) { #specific attr list given
      next if !exists($attr{$an});
    } else {
     #check excluded attributes only if there's no specific attr list
     next if exists($xcludeAttr{$an});
    }
    $ah{$an}=$av;
  }
  #store feature type as well:
  $ah{'ftype'}=$t[2];
  $tid=~s/\.([a-z]+)\d+$//;
  $sd{$tid}={ %ah };
}
close(SRCIN);

shift(@ARGV) if $ARGV[0] eq '-';
#open(TGTIN, $tgtgff)||die("Error: cannot open $targetgff! ($!)\n");
#while (<TGTIN>) {
while (<>) {
  my $line=$_;
  chomp;
  my @t=split(/\t/);
  if (m/^#/ && @t<9) {
    print $line;
    next;
  }
  if ($gfftrack) {
    $t[1]=$gfftrack;
    $line=join("\t",@t)."\n";
  }
  my @astr=split(/\s*\;\s*/, $t[8]);
  my @a;
  my %ah= map { @a=split(/\s*=\s*/);$a[0]=>$a[1] } @astr;
  my $tid=$ah{'ID'};
  if (!$tid) { #exon line?
    print $line;
    next;
   }
 $tid=~s/\.([a-z]+)\d+$//;
 my $sa=$sd{$tid};
 if (!$sa) {
  print STDERR "Warning: ID $tid not found in source data!\n";
  print $line;
  next;
 }
 #important: source feature name too:
 $t[2]=delete($sa->{'ftype'}) || die("Error: could not find ftype attribute!\n");
 print join("\t",@t[0..7])."\tID=$tid";
 foreach my $attrval (@astr) {
   my ($an,$av)=split(/\s*=\s*/, $attrval);
   next if $an eq 'ID';
   ## no longer needed, we filtered earlier
   next if ($removeattrs && exists($removeAttr{$an}));
   my $sv=$sa->{$an};
   if ($sv) {
      print ";$an=".$sv;
      delete $sa->{$an}; #so we don't print it later
   }
   else {
      print ";$an=$av";
   }
 }
 #now print the remaining attributes from source
 foreach my $a (sort keys(%{$sa})) {
   next if ($removeattrs && exists($removeAttr{$a}));
   my $av=$sa->{$a};
   print ";$a=$av";
 }
 print "\n";
}
#close(TGTIN);
# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

