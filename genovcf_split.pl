#!/usr/bin/perl
use strict;
use IO::Zlib;  #requires Compress::Zlib, of course
use Getopt::Std;
# this might need command :
#     ulimit -n 10000
my $usage = q/Usage:
 genovcf_split.pl [-m id2br] [-p dir_prefix] [-n numfiles] genotypes.vcf[.gz]
 Options:
   -m : file mapping sample ID to BrNum
   -p : prefix for the output sub-directories
   -n : number of files (samples) per directory
   
/;
getopts('p:n:m:') || die($usage."\n");
my $pre=$Getopt::Std::opt_p || 'vcfsplit';
my $n=$Getopt::Std::opt_n || 500;
my $brf=$Getopt::Std::opt_m;
my $fin=$ARGV[0];
die($usage) if ($fin =~ m/^\-+h/) || ! -f $fin;
my %mbr; # maps sampleID => brnum
if ($brf && -f $brf) {
 open(M, $brf) || die("Error opening $brf\n");
 while (<M>) {
   chomp;
   my ($sid, $brnum)=split;
   $mbr{$sid}=$brnum;
 }
 close(M)
} else { $brf=0 }

my $fh = IO::Zlib->new($fin, 'rb') || die("Error opening $fin\n");
my $hdr=''; #keep common header lines here
my @fws; #array of writing file handles
my @samples; 
while (<$fh>) {
  if (@samples>0) {
   ## reading data
   #closeAll();exit(0); #for now
   chomp;
   my @t=split(/\t/);
   my $o=join("\t",@t[0..8])."\t";
   my $i=0;
   foreach my $fw (@fws) {
     print $fw $o.$t[9+$i]."\n";
     $i++;
   }
   next
  }
  ## still in header lines
  if (m/^#CHROM\s+/) { #column headers, we get sample IDs here
    chomp;
    my @t=split(/\t/);
    @samples=@t[9 .. $#t];
    print STDERR "Found ".scalar(@samples)." samples.\n";    
    ## write header here and sample list
    open(H, ">$pre.header.vcf") || die("Error writing header file\n");
    print H $hdr;
    close(H);
    open(SL, ">$pre.samples.lst") || die("Error writing samples lst\n");
    print SL join("\n", @samples)."\n";
    close(SL);
    # -- create the directories and open all files for write:
    my $cs=0; #count samples written
    my @dnrange=(1, $n);
    $dnrange[1]=scalar(@samples) if $n>@samples;
    my $dir=join("_", $pre, @dnrange);
    mkdir($dir);
    foreach my $sid (@samples) {
      $cs++;
      if ($cs>1 && ($cs % $n) == 1) { #next dir
        @dnrange=($cs, $cs+$n-1);
        $dnrange[1]=scalar(@samples) if $n>@samples;
        $dir=join("_", $pre, @dnrange);
        mkdir($dir);
      }
      my $sid=$samples[$cs-1];
      my $fn="$dir/$sid".'.vcf.gz';
      my $fw=IO::Zlib->new($fn, 'wb') || die "Error: cannot create $fn : $!";
      print $fw "##fileformat=VCFv4.1\n";
      print $fw join("\t", '#CHROM', 'POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',$sid)."\n";
      push(@fws, $fw);
    }
    next
  }
  next if m/^##bcftools/;
  $hdr.=$_;
}

closeAll();


##======================================================

sub closeAll {
  print STDERR "closing files...\n";
  $fh->close;
  foreach my $fw (@fws) {
     $fw->close;
  }
}
