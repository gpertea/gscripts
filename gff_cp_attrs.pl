#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
 gff_cp_attrs.pl [-a 'attr_lst'] [ -r 'remove_lst'] [-F] [-E | -C] [-V] \
   [-m ID.map] src.gff target.gff
 Transfer the additional transcript attributes from src.gff to the 
 transcripts with the same IDs from target.gff (or mapped IDs with -m option).
 Optionally, it can also transfer or replace all exons and\/or CDS in 
 target.gff with the exons and CDS from src.gff
  
 The output will be a GFF3 which is an "enriched" version of target.gff 
 with new attributes (or exons\/CDS) added from src.gff
  
 Options:
  -a a list of attribute names (comma, period or colon delimited)
     which should be transferred from src.gff (default is to transfer
     only attributes that are not already present in target.gff)
  -r a list of transcript attributes to be discarded from the GFF output
  -O override (replace) existing attribute values in target
  -m instead of looking for a direct match of transcript\/gene IDs, map the
     target.gff IDs to the src.gff IDs using a mapping file ID.map (with 
     the target IDs in the 1st column and corresponding src.gff IDs in the 2nd)
  -t override the GFF track column in the output with the given value
  -F overwrite original target feature name with the one in the source file
  -C checks if there is exon compatibility between the ID matched source and
     target transcripts, and if so, transfer the CDS features from src
  -E transfer all exon and CDS features from src.gff to target.gff, 
     completely replacing existing exon\/CDS features in target.gff
  -A for -E, do not replace exons for target transcripts with ASSEMBLED=yes
  -K for -E, never remove CDS (when the source has no CDS, but the target does)
  -V verbose processing (logging the applied changes)

Note: input is expected to be GFF3 "normalized", with full 'exon' and\/or
'CDS' features (no *codon* or *UTR* features)
/;

umask 0002;
my $cmdline=$0.' '.join(' ',@ARGV);
getopts('OVFEAKCt:r:a:m:o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
my %xcludeAttr= ( 'ID' => 1, 'exonCount'=> 1, 'exons'=> 1, 'CDS'=>1, CDSphase=>1 );
my %removeAttr;
my %idmap; # srcID => targetID
my $CDStransfer=$Getopt::Std::opt_C;
my $exonReplace=$Getopt::Std::opt_E; #completely replace exon/CDS lines
my $protectAssembled=$Getopt::Std::opt_A;
my $protectCDS=$Getopt::Std::opt_K;
my $attrReplace=$Getopt::Std::opt_O;
my $featSrc=$Getopt::Std::opt_F;
my $gfftrack=$Getopt::Std::opt_t;
my $attrlist=$Getopt::Std::opt_a;
my $mapfile=$Getopt::Std::opt_m;
my $removeattrs=$Getopt::Std::opt_r;
my $verbose=$Getopt::Std::opt_V;
print STDERR "Command line: $cmdline\n" if ($verbose);
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
####vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv####
my %satrs; ##source attributes: ID => { attr => value, ... }
my %se;    ## source data:      ID => [ [@exons], [@CDS], strand ]
            # (@exons and @CDS are lists of [start, end, gffline])
####^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^####
if ($mapfile) {
  open(MAP, $mapfile) || die("Error opening $mapfile ($!)!\n");
  while(<MAP>) {
    next if m/^#/;
    chomp;
    my @t=split;
    $idmap{$t[1]}=$t[0];
  }
  close(MAP);
}

open(SRCIN, $srcgff)||die("Error: cannot open $srcgff! ($!)\n");

while(<SRCIN>) { #source file to copy attributes for each ID
  #just store all ids with their attrs
  my $srcline=$_;
  chomp;
  my @t=split(/\t/);
  if ($gfftrack) {
    $t[1]=$gfftrack;
    $srcline=join("\t",@t)."\n";
  }
  next if m/^#/ || !$t[8];
  my ($id)=($t[8]=~m/\bID=([^;]+)/);
  my ($pid)=($t[8]=~m/\bParent=([^;]+)/);
  if ($mapfile) {
   if ($id) {
      my $tid=$idmap{$id};
      if ($tid) {
         $id=$tid;
         $t[8]=~s/\bID=([^;]+)/ID=$id/;
         $srcline=join("\t",@t)."\n";
      }
   }
   if ($pid) {
      my $p=$idmap{$pid};
      if ($p) {
         $pid=$p;
         $t[8]=~s/\bParent=([^;]+)/Parent=$pid/;
         $srcline=join("\t",@t)."\n";
      }
   }
  }
  #     0      1     2      3       4      5        6        7      8
  my ($chr, $track, $f, $fstart, $fend, $fscore, $strand, $phase, $ats)=@t;
  my %ah;
  #my %ah= map { @a=split(/\s*=\s*/);$a[0]=>$a[1] } @astr;
  my @astr=split(/\s*\;\s*/, $ats);
  if (@astr>0) {
    $astr[0]=~s/^\s+//;
    $astr[-1]=~s/\s+$//;
  }
  my $isExon=($f=~m/^(exon|CDS)$/i);
  my $isCDS=($f=~m/^CDS$/i);
  #next unless ($id || $CDStransfer || $exonReplace);
  my @alst; #list of attributes stored, in order they were found
  foreach my $avstr (@astr) {
    my ($an, $av)=split(/\s*=\s*/, $avstr);
    next if $an eq 'ID';
    $av=$pid if $an eq 'Parent';
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
  #special storing of attribute list in order they were found
  $ah{';;'}=[@alst];
  #$id=~s/\.([a-z]+)\d+$//; #trim GMAP suffixes?
  $satrs{$id}={ %ah };
  #my $parent=$ah{'Parent'};
  if ($pid && $isExon && ($CDStransfer || $exonReplace)) {
    my $ed=$se{$pid}; #exon data for this parent ID
    #print STDERR "[DBG]>> trying to add $f \[$fstart, $fend\] to \$se\{$parent\}\n";
    unless ($ed) { $ed=[[],[], $strand]; $se{$pid}=$ed; }
    push(@{$$ed[ $isCDS ? 1 : 0 ]}, [$fstart, $fend, $srcline]);
  }
}
close(SRCIN);

#normalize exon data ?
#if ($CDStransfer || $exonReplace) {
# foreach my $ed (values (%se)) {
#   normalizeExons($$ed[0]);
# }
#}

# read the target file
shift(@ARGV) if $ARGV[0] eq '-';
####vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv####
my @ts;    #transcripts or gene IDs, in the order they were found in the input stream

my %tdata; # id=>[transcript/gene line,  [@exons], [@cds], strand];
############              0                 1[]     2[]     3
####^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^####

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
  if (@avs>0) {
    $avs[0]=~s/^\s+//;
    $avs[-1]=~s/\s+$//;
  }
  my %ah; #hash of current target attributes
  my @alst; #list of current target attributes
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
  my $isExon=$pid && $t[2]=~m/^(exon|CDS)$/i;
  my $isCDS=$isExon && $t[2]=~m/^CDS$/i;
  my $td;
  if ($isExon) { #exon or CDS
     $td=$tdata{$pid};
     if (!$td) { #shouldn't happen, exons shouldn't be in %se
         my $sed=$se{$pid};
         if ($sed && $$sed[2] ne $t[6]) {
           $t[6]=$$sed[2];
           $line=join("\t",@t)."\n";
         }
         $td=['', [], [], $t[6] ]; 
         $tdata{$pid}=$td;
     } elsif ($td->[3] ne $t[6]) {
       $t[6]=$td->[3];
       $line=join("\t",@t)."\n";
     }
     if ($CDStransfer || $exonReplace) {
         push(@{$$td[ $isCDS ? 2 : 1 ]}, [$t[3], $t[4], $line]);
     }
     next; #exon/CDS lines are stored as they are, replaced later if requested
  }
  #else { # non-exon (i.e. transcript or gene)
  die("Error parsing no-ID line:\n$line\n") unless $id;
  $td=$tdata{$id};
  die("Error: ID $id duplicated?\n") if $td;
  my $sed=$se{$id};
  if ($sed && $$sed[2] ne $t[6]) {
     $t[6]=$$sed[2];
     $line=join("\t",@t)."\n";
  }
  $td=[$line, [], [], $t[6] ];
  $tdata{$id}=$td;
  push(@ts, $id);
  #}
  #$id=~s/\.([a-z]+)\d+$//; #remove GMAP prefix?
  my $sa=$satrs{$id} if $id; #source attributes data for this ID
  if (!$sa) {
    ##ID data not found in source, just keep it as is
    #if ($isExon) { push(@{$$td[1]}, $line) }
    #      else { $$td[0]=$line }
    next;
  }

  if ($featSrc) {
    ##NOTE: source feature name should be transferred too, could be different!
     # exception: if the source feature is gene but it has exons or CDS data, keep the target feature name
     # or make it a "transcript"
     my $ofname=$t[2]; #original feature name
     $t[2]=delete($sa->{'~f'}) || die("Error: could not find feature type attribute for $id!\n");
     if ($t[2] eq 'gene' && (@{$$sed[0]}>0 || @{$$sed[1]}>0)) {
       $t[2]= ($ofname ne 'gene') ? $ofname : 'transcript'; #if it has exons, it's a transcript
     }
  }
  else  {
     delete($sa->{'~f'});
     # still don't allow 'gene' as a feature name if it has exons/introns and a parent
     $t[2]='transcript' if ($pid && $t[2] eq 'gene' && (@{$$sed[0]}>0 || @{$$sed[1]}>0));
  }
  my $bl=join("\t",@t[0..7])."\tID=$id"; #building the GFF line -- with adjusted attributes
  my $tattrs=''; #gather target attributes here, to append at the end of the line
  foreach my $a (@alst) {
    if ($a eq 'ID' 
         || ($removeattrs && exists($removeAttr{$a})) ) {
      delete $sa->{$a};
      next;
    }
    my $v=$ah{$a};
    my $sv=$sa->{$a};
    if (length($sv) && $attrReplace) {
      $bl.=";$a=$sv";
    }
    else {
      $tattrs.=";$a=$v"; #print the original value
    }
    delete $sa->{$a} if length($sv); #so we don't print it later
  }
  #now print the remaining attributes from source
  foreach my $a (@{$sa->{';;'}}) {
    next if ($removeattrs && exists($removeAttr{$a}));
    my $av=$sa->{$a};
    $bl.=";$a=$av" if length($av);
  }
  $bl.=$tattrs;
  #if ($isExon) { push(@{$$td[1]}, "$bl\n") }
  # else { $$td[0]="$bl\n" }
  $$td[0]="$bl\n";
} # while <target>

#now for each target gene/transcript, perform any exon/CDS replacements, as requested
# and print it
foreach my $tid (@ts) {
   if ($tid=~m/^#[^\n]+\n$/) {
     print $tid;
     next;
   }
   my $td=$tdata{$tid};
   if ($CDStransfer && @{$$td[1]}>0) {
      #die("Warning: record ID $tid has no exons?!\n") if (@{$$td[1]}==0);
      #normalizeExons($$td[1]);
      ## -- for $CDStransfer, check that the exon structure is matching the source
      ##    and if so, transfer CDS
      my $sd=$se{$tid};
      if ($sd) { #only valid if source ID exists
        if (compatExons($$td[1], $$sd[0])) {
            if (@{$$sd[1]}>0) {
              print STDERR "Info: existing CDS override for $tid\n"
                if (@{$$td[2]}>0);
              @{$$td[2]}=@{$$sd[1]};
            }
        }
        else { print STDERR "Warning: exons structure mismatch for $tid, CDS unchanged!\n" }
       #} #else {
       #print STDERR "[DBG]>> Warning: no source data found for $tid\n";
      }#source ID to compare to
   } #-- if $CDStransfer
   elsif ($exonReplace) {
       my $sd=$se{$tid};
       if ($sd) {
         my $rex=1; #replace exons?
         my $rCDS=1; #replace CDS?
         if ($protectAssembled && $$td[0]=~m/ASSEMBLED=yes/) {
            if (compatExons($$td[1], $$sd[0])) {
               my ($rd, $re)=($$td[1], $$sd[0]);
               if ($$rd[0]->[0]!=$$re[0]->[0] ||
                         $$rd[-1]->[1]!=$$re[-1]->[1]) {
                 $rex=0;
                 print STDERR "Warning: ASSEMBLED $tid has imperfect exon match with source, exons unchanged\n";
               }
            }
            else {
               $rex=0;
               $rCDS=0;
               print STDERR "Warning: ASSEMBLED $tid has incompatible exon structure with source, unchanged\n";
            }
            if ($rCDS && @{$$sd[1]}==0 && @{$$td[2]}>0) {
              $rCDS=0;
              print STDERR "Warning: ASSEMBLED $tid has CDS in target but not source, prevent CDS removal\n";
            }
         }
         if ($rex) {
            print STDERR "Info: exons replaced for $tid\n";
            @{$$td[1]}=@{$$sd[0]};
         }
         if ($rCDS && $protectCDS && @{$$td[2]}>0 && @{$$sd[1]}==0) {
            print STDERR "Warning: prevent CDS removal for $tid\n";
            $rCDS=0;
         }
         if ($rCDS) {
           my $cop;
           my $msg='Info';
           if (@{$$sd[1]}>0) {
             $cop=(@{$$td[2]}>0) ? 'replaced' : 'added';
           } elsif (@{$$td[2]}>0) {
              $cop='removed';
              $msg='Warning';
           }
           if ($cop) {
             @{$$td[2]}=@{$$sd[1]};
             print STDERR "$msg: CDS $cop for $tid\n";
           }
         }
       }
   }
   # -- print target GFF record and its exons/CDS data
   print $$td[0];
   foreach my $el (@{$$td[1]}) {
      print $el->[2];
   }
   foreach my $cl (@{$$td[2]}) {
      print $cl->[2];
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

sub compatExons {
  my ($re, $rd)=@_;
  return 0 if $#$re != $#$rd;
  foreach my $i (0 .. $#$re) {
    my ($nl, $nr);
    $nl = ($i==0) ? ($$re[0]->[0] > $$rd[0]->[0]) : ($$re[$i]->[0] != $$rd[$i]->[0]);
    $nr = ($i==$#$re) ? ($$re[$i]->[1] < $$rd[$i]->[1]) : ($$re[$i]->[1] != $$rd[$i]->[1]); 
    return 0 if ($nl || $nr);
  }
  return 1;
}
