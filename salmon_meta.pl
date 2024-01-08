#!/usr/bin/perl
use strict;
use Getopt::Std;
use File::Basename;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
  salmon_meta.pl [dir1, dir2, ...]
This script searches recursively for all meta_info.json files in the given 
paths (or the current directory if no paths to search in are given) and 
outputs a table with the relevant metrics reported per sample.
/;

umask 0002;
getopts('o:') || die($usage."\n");
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --
if (scalar(@ARGV)==0) {
  push(@ARGV, '.');
} else {
 foreach my $d (@ARGV) {
   die("Error: path $d not found!\n") unless -d $d;
 }
}

my $ffind='meta_info.json';
print join("\t", qw{sampleID frag_len_mean frag_len_sd num_reads num_mapped perc_mapped})."\n";
foreach my $dir_path (@ARGV) {
   my @meta_files;
   get_meta($dir_path, $ffind, \@meta_files);
   my $ns=scalar(@meta_files);
   next if $ns==0; #number of samples
   #print STDERR join("\n", @meta_files), "\n";
   my (@uniqIDs,  $idlevel);
   my @dsplit= map { [ split(/\/+/) ] } @meta_files;
   my @depths = map { scalar(@$_) } @dsplit;
   my $d=$depths[0];
   for (my $i=1;$i<=$#depths;$i++) {
     if ($depths[$i] ne $d) {
       die("Error: metadata files not found at the same level\n".
       $meta_files[0]."\n".$meta_files[$i]."\n");
     }
   }
   $d-=2;
   for(my $i=$d;$i>=0;$i--) {
     my %uk; #unique keys
     @uniqIDs=();
     $idlevel=$d-$i+1;
     my $r=0;
     foreach my $dspl (@dsplit) {
       my $id=$$dspl[$i];
       $uk{$id}=1;
       $r++;
       if (scalar(keys(%uk))<$r) {
          @uniqIDs=();
          undef $idlevel;
          last
       }
       push(@uniqIDs, $id);
     }
     last if defined($idlevel);
   }
   #if (defined($idlevel)) {
   #  print STDERR "Found IDs at level $idlevel\n",
   #     join(", ", @uniqIDs),"\n";
   #}
   die("Error: could not find unique IDs for the $ns meta_info files!\n") unless defined($idlevel);
   die("Error: unique IDs (".scalar(@uniqIDs)." and pathnames ($ns) do not match!\n")
     unless $ns==scalar(@uniqIDs);
  my ($frag_len_mean, $frag_len_sd, $num_reads, $num_mapped, $perc_mapped);
   for (my $i=0;$i<$ns;$i++) {
     open(F, '<'.$meta_files[$i]) || die("Error opening $meta_files[$i]\n");
     while (<F>) {
       if (m/"frag_length_mean":\s*([\d\.]+)/) {
         $frag_len_mean=floatfmt($1);
         next
       }
       if (m/"frag_length_sd":\s*([\d\.]+)/) {
         $frag_len_sd=floatfmt($1);
         next
       }
       if (m/"num_processed":\s*([\d\.]+)/) {
         $num_reads=$1;
         next
       }
       if (m/"num_mapped":\s*([\d\.]+)/) {
         $num_mapped=$1;
         next
       }
       if (m/"percent_mapped":\s*([\d\.]+)/) {
         $perc_mapped=floatfmt($1);
         #next
         last  #the last one we need
       }
     }
     close(F);
   print join("\t", $uniqIDs[$i], $frag_len_mean, $frag_len_sd, $num_reads, $num_mapped, $perc_mapped)."\n";
   }
}



# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************
sub get_meta {
    my ($dir_path, $ff, $meta_files) = @_;
    opendir(my $dh, $dir_path) || die "Can't open directory: $!";
    my @files = readdir($dh);
    closedir($dh);

    foreach my $fd (@files) {
        next if ($fd eq '.' || $fd eq '..');
        if (-d "$dir_path/$fd") {
            get_meta("$dir_path/$fd", $ff, $meta_files);
        } elsif ($fd eq $ff) {
            push @$meta_files, "$dir_path/$fd";
        }
    }
}

sub floatfmt {
 local $_=shift(@_);
 my $n=shift(@_);
 $n=3 unless defined($n);
 $_=sprintf("%.${n}f", $_);
 s/0+$//; s/\.$//;
 return $_;
}
