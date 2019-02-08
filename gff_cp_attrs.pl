#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 gff_cp_attrs.pl [-a 'attr_lst'] [ -r 'remove_lst'] attrs_src.gff target.gff
 Copies the additional transcript attributes from attrs_src.gff to the 
 transcripts with the same IDs from target.gff. Exon and CDS features are
 unchanged.
  
 The output will be a GFF3 which is the same with target.gff but with
 attributes from attrs_src.gff copied over for the matching transcript ID.
  
 Options:
  -a denotes a list of attribute names (comma, period or colon delimited)
     which should be copied over from attrs_src.gff (default is to transfer
     all new attributes)
  -r denotes a list of transcript attributes to be discarded from the output
  -t override the GFF track column in the output with the given value
  -C checks if there is matching of the boundary and exon definitions 
     of transcripts with the same ID in the two gff input files and 
     transfer the CDS features from attrs_src.gff to target.gff
/;

umask 0002;
getopts('Ct:r:a:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my %xcludeAttr= ( 'ID' => 1, 'exonCount'=> 1, 'exons'=> 1, 'CDS'=>1, CDSphase=>1 );
my %removeAttr;
my $CDStransfer=$Getopt::Std::opt_C;
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
my %satrs; #ID => { attr => value, ... }
my %se; #ID => [ [@exons], [@CDS] ] - just plain [start, end]
open(SRCIN, $srcgff)||die("Error: cannot open $srcgff! ($!)\n");

while(<SRCIN>) { #source file to copy attributes for each ID
  #just store all ids with their attrs
  my $srcline=$_;
  chomp;
  my ($chr, $track, $f, $fstart, $fend, $fscore, $strand, $phase, $ats)=split(/\t/);
  next if m/^#/ || !$ats;
  my ($id)=($ats=~m/\bID=([^;]+)/);
  my %ah;
  #my %ah= map { @a=split(/\s*=\s*/);$a[0]=>$a[1] } @astr;
  my @astr=split(/\s*\;\s*/, $ats);
  if (@astr>0) {
    $astr[0]=~s/^\s+//;
    $astr[-1]=~s/\s+$//;
  }
  #next unless $tid; #ignore lines without ID
  next unless ($id || $CDStransfer);
  my @alst; #list of attributes stored, in order they were found
  foreach my $avstr (@astr) {
    my ($an, $av)=split(/\s*=\s*/, $avstr);
    if ($an eq 'ID') {
      die("Error: parsed ID ($id) not matching attr value for ID ($av)\n")
         if $id ne $av;
      next;
    }
    
    if ($attrlist) { #specific attr list given
      next if !exists($attr{$an});
    } else {
      #check excluded attributes only if there's no specific attr list
      next if exists($xcludeAttr{$an});
    }
    push(@alst, $an);
    $ah{$an}=$av;
  }
  #special storing of feature type:
  $ah{'~f'}=$f;
  #special storing of attribute list in order
  $ah{';;'}=[@alst];
  #$id=~s/\.([a-z]+)\d+$//; #trim GMAP suffixes?
  $satrs{$id}={ %ah };
  #my ($parent)=($ats=~m/\bParent=([^;]+)/);
  my $parent=$ah{'Parent'};
  if ($parent && $CDStransfer) {
    my $ed=$se{$parent}; #exon data for this parent ID
    #print STDERR "[DBG]>> trying to add $f \[$fstart, $fend\] to \$se\{$parent\}\n";
    if ($f=~m/(exon|utr)/i) {
       unless ($ed) { $ed=[[],[]]; $se{$parent}=$ed; }
       push(@{$$ed[0]}, [$fstart, $fend]);
    }
    elsif ($f=~m/(cds|codon)/i) {
       unless ($ed) { $ed=[[],[]]; $se{$parent}=$ed; }
       push(@{$$ed[1]}, [$fstart, $fend, $srcline]);
       push(@{$$ed[0]}, [$fstart, $fend]); #merge CDS into exons (to be normalized, just in case)
    }
  }

}
close(SRCIN);
#normalize exon data:
if ($CDStransfer) {
 foreach my $ed (values (%se)) {
   normalizeExons($$ed[0]);
 }
}
# read the target file
shift(@ARGV) if $ARGV[0] eq '-';
my @ts; #transcripts or gene IDs, in the order they were encountered
my %tdata; # id=>[line, [@nexon/CDS_lines], [@exons], [@cds]];
############        0         1[]              2[]     3[]   
while (<>) {
  my $line=$_;
  chomp;
  my @t=split(/\t/);
  if (m/^#/) {
    push(@ts, $line);
    next;
  }
  die("Error parsing line in target GFF:\n$line") unless $t[8];
  #my ($id)=($t[8]=~m/\bID=([^;]+)/);
  if ($gfftrack) {
    $t[1]=$gfftrack;
    $line=join("\t",@t)."\n";
  }
  my @avs=split(/\s*\;\s*/, $t[8]);
  if (@avs>0) { #should be
    $avs[0]=~s/^\s+//;
    $avs[-1]=~s/\s+$//;
  }
  my %ah;
  my @alst;
  foreach my $av (@avs) {
    my ($a, $v)=split(/\s*=\s*/, $av);
    $a=~s/^\s+//;
    $v=~s/\s+$//;
    push(@alst, $a);
    #print STDERR "[DBG]>> found <$a>:<$v> pair\n";
    $ah{$a}=$v;
  }
  my $id=$ah{'ID'};
  my $pid=$ah{'Parent'};
  my $isExon=$pid && $t[2]=~m/(exon|utr|cds|codon)/i;
  my $td;
  if ($isExon) { #exon or CDS
     $td=$tdata{$pid};
     if (!$td) { $td=['', [], [], [] ]; $tdata{$pid}=$td; }
     if ($CDStransfer) {
       if ($t[2]=~m/(cds|codon)/i) {
         push(@{$$td[3]}, [$t[3], $t[4]]);
         push(@{$$td[2]}, [$t[3], $t[4]]);
       }
       else { #exon/UTR
         push(@{$$td[2]}, [$t[3], $t[4]]);
       }
     }
  } else { #non-exon (i.e. transcript or gene)
    die("Error parsing no-ID line:\n$line\n") unless $id;
    $td=$tdata{$id};
    die("Error: ID $id duplicated?\n") if $td;
    $td=['', [], [], [] ];
    $tdata{$id}=$td;
    push(@ts, $id);
  }
  #$id=~s/\.([a-z]+)\d+$//; #remove GMAP prefix?
  my $sa=$satrs{$id} if $id; #source attributes data for this ID (exons too, if they have ID)
  if (!$sa) {
    #ID data not found in source, just keep it as is
    if ($isExon) { push(@{$$td[1]}, $line) }
          else { $$td[0]=$line }
    next;
  }
  #important: source feature name should be transferred too, could be different
  $t[2]=delete($sa->{'~f'}) || die("Error: could not find feature type attribute for $id!\n");
  my $bl=join("\t",@t[0..7])."\tID=$id"; #building line -- with adjusted attributes
  foreach my $a (@alst) {
    if ($a eq 'ID' 
         || ($removeattrs && exists($removeAttr{$a})) ) {
      delete $sa->{$a};
      next;
    }
    my $v=$ah{$a};
    my $sv=$sa->{$a};
    if ($sv) {
      $bl.=";$a=$sv";
      delete $sa->{$a}; #so we don't print it later
    }
    else {
      $bl.=";$a=$v"; #print the original value
    }
  }
  #now print the remaining attributes from source
  foreach my $a (@{$sa->{';;'}}) {
    #next if ($removeattrs && exists($removeAttr{$a}));
    my $av=$sa->{$a};
    $bl.=";$a=$av" if $av;
  }
  if ($isExon) { push(@{$$td[1]}, "$bl\n") }
         else { $$td[0]="$bl\n" }
} # while <target>

foreach my $tid (@ts) {
   if ($tid=~m/^#[^\n]+\n$/) {
     print $tid;
     next;
   }
   my $td=$tdata{$tid};
   if ($CDStransfer && @{$$td[2]}>0) {
      #die("Warning: record ID $tid has no exons?!\n") if (@{$$td[2]}==0);
      normalizeExons($$td[2]);
      print STDERR "[DBG]>> ".scalar(@{$$td[2]})." exons normalized for $tid\n";
      ## -- now validate if exon structure is matching the source
      ##    and if so, transfer CDS
      my $sd=$se{$tid};
      if ($sd) { #only valid if source ID exists
        if (matchingExons($$td[2], $$sd[0])) {
            print STDERR "[DBG]>> exon structures matching for $tid\n";
            if (@{$$sd[1]}>0) {
              print STDERR "WARNING: existing target CDS override for $tid\n"
                if (@{$$td[3]}>0);
              #@{$$td[3]}=@{$$sd[1]};
              #now we MUST delete any existing CDS lines in $$td[1]
              my @nocds=grep { !m/^\S+\t[^\t]+\t(\w*CDS\w*|\w*codon)\t/i } @{$$td[1]};
              foreach my $cd (@{$$sd[1]}) {
                 if ($gfftrack) {
                   my @t=split(/\t/, $$cd[2]);
                   $t[1]=$gfftrack;
                   $$cd[2]=join("\t",@t);
                 }
                 push(@nocds, $$cd[2]); #print line as stored in the source
              }
              @{$$td[1]}=@nocds;
            }
        }
        else { print STDERR "WARNING: exons structure mismatch for $tid!\n" }
      } else {
        print STDERR "[DBG]>> Warning: no source data found for $tid\n";
      }#source ID to compare to
   } #-- if $CDStransfer
   print $$td[0];
   foreach my $el (@{$$td[1]}) {
     print $el;
   }
} #for each target ID

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
}

#************ Subroutines **************

sub normalizeExons {
 my ($rex)=@_;
 @$rex=sort { $a->[0] <=> $b->[0] } @$rex;
 my $i=0;
 while ($i+1<@$rex) {
 # $$rev[$i] vs $$rex[$i+1]
   my $exdist=$$rex[$i+1]->[0]-$$rex[$i]->[1]; #inter-exon distance
   if ($exdist<=1) { #overlapping/adjacent exons
       $$rex[$i]->[1]=$$rex[$i+1]->[1] if $$rex[$i+1]->[1]>$$rex[$i]->[1];
       splice(@$rex, $i+1, 1);
   }
   else { ++$i }
 }
}

sub matchingExons {
 my ($re, $rd)=@_;
 return 0 if scalar(@$re) != scalar(@$rd);
 for (my $i=0;$i<@$re;$i++) {
   return 0 if ($$re[$i]->[0]!=$$rd[$i]->[0] ||
         $$re[$i]->[1]!=$$rd[$i]->[1]);
 }
 return 1;
}
