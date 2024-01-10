#!/usr/bin/perl
use strict;
use Getopt::Std;
use File::Find;

umask 0002;
## quick parsing of basic adapter trimming summary from _fastp.json and _fastp.log output
## provide only the prefix, _fastp.json and _fastp.log are automatically added
my $usage = q/Usage:
  fastp_meta.pl [dir1, dir2, ...]
This script searches recursively for all _fastp.json and _fastp.log files in the given 
paths (or the current directory if no paths to search in are given) and 
outputs a table with the relevant metrics reported per sample.
/;

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

my $fmask='.+_fastp\.json';
my @ffiles; # found _fastp json files
my $hdr=1;
foreach my $path (@ARGV) {
   my @ffiles;
   # Get a custom subroutine for this iteration
   my $wanted = make_wanted(\@ffiles);
   find($wanted, $path);
   my (@uniqIDs,  $idlevel);
   $idlevel=getUIDs(\@ffiles, \@uniqIDs); ## this populates \@uniqIDs
   ## now process @ffiles results
   my $i=-1;
   foreach my $fj (@ffiles) {
     my $flog=$fj;
     $flog=~s/\.json$/.log/;
     $i++;
     next unless -f $flog;
     #print STDERR "$flog - $fj\n";
     my $sid=$uniqIDs[$i];
     my ($bef, $aft, $fr, $ins, $adcut, # sections
         $bef_r1r, $bef_r1b, $bef_r2r,  $bef_r2b,  $bef_r1_avglen, $bef_r2_avglen, 
         $aft_r1r, $aft_r1b, $aft_r2r, $aft_r2b,  $aft_r1_avglen, $aft_r2_avglen, 
         $fr_passed, $fr_short, $fr_ins_peak, $fr_trim_r, $fr_trim_b);
     open(L, $flog) || die(" Error opening $flog\n");
     open(J, $fj) || die(" Error opening $fj\n");
     ## most info is obtained from the log file, we only need JSON for read1/2_mean_length
     my %lv=( 'before.r1.reads'=>\$bef_r1r, 'before.r1.bases'=>\$bef_r1b,
              'before.r2.reads'=>\$bef_r2r, 'before.r2.bases'=>\$bef_r2b,
              'after.r1.reads'=>\$aft_r1r, 'after.r1.bases'=>\$aft_r1b,
              'after.r2.reads'=>\$aft_r2r, 'after.r2.bases'=>\$aft_r2b,
               );
     while (<L>) {
       if ($fr) {
         if (m/Insert size peak/) {
           ($fr_ins_peak)=(m/(\d+)$/);
           last;
         }
         if (m/reads passed filter:/) {
           ($fr_passed)=(m/(\d+)$/);
           next
         }
         if (m/reads failed due to too short:/) {
           ($fr_short)=(m/(\d+)$/);
           next
         }
         if (m/reads with adapter trimmed:/) {
           ($fr_trim_r)=(m/(\d+)$/);
           next
         }
         if (m/bases trimmed due to adapters:/) {
           ($fr_trim_b)=(m/(\d+)$/);
           next
         }
         next;
       } elsif (m/^Filtering result:/) {
         $fr=1;
         next;
       }
       ## not in filtering results yet:
       for my $ba (('before', 'after')) {
         for my $rn ((1,2)) {
           if (m/Read$rn $ba filtering:/) {
             for my $rb (('reads', 'bases')) {
               $_=<L>;
               if (m/otal $rb:\s*(\d+)/) {
                 my $v=$1;
                 my $k="$ba.r$rn.$rb";
                 ${$lv{$k}}=$v;
                 my $kv=${$lv{$k}};
               } else { expErr("r$rn: total $rb"); }
             } # reads, bases
           }
         }
       }
     }
     close(L);
     ($bef, $aft)=(0,0);
     while (<J>) {
        if ($aft) {
          if (m/"read1_mean_length":\s*(\d+)/) {
            $aft_r1_avglen=$1;
            next;
          }
          if (m/"read2_mean_length":\s*(\d+)/) { 
            $aft_r2_avglen=$1;
            last;  ## quit parsing
          }
          next;
        } elsif (m/"after_filtering":/) {
          $aft=1;
          next;
        }
        if ($bef) {
          if (m/"read1_mean_length":\s*(\d+)/) {
            $bef_r1_avglen=$1;
            next;
          }
          if (m/"read2_mean_length":\s*(\d+)/) { 
            $bef_r2_avglen=$1;
            next;
          }
        } elsif (m/"before_filtering":/) {
          $aft=0;
          $bef=1;
          next;
        }
     }
     close(J);
     if ($hdr) {
       print join("\t", qw(sampleID r1_reads r1_bases r2_reads r2_bases));
       print "\t".join("\t", qw(trim_r1_reads trim_r1_bases trim_r2_reads trim_r2_bases));
       print "\t".join("\t", qw(trim_passed trim_too_short trimmed_reads trimmed_bases fastp_insert_peak));
       print "\t".join("\t", qw(r1_avg_len r2_avg_len trim_r1_avg_len trim_r2_avg_len))."\n";
       $hdr=0;
     }
     print join("\t", $sid, $bef_r1r, $bef_r1b, $bef_r2r, $bef_r2b);
     print "\t".join("\t", $aft_r1r, $aft_r1b, $aft_r2r, $aft_r2b);
     print "\t".join("\t", $fr_passed, $fr_short, $fr_trim_r, $fr_trim_b, $fr_ins_peak);
     print "\t".join("\t", $bef_r1_avglen, $bef_r2_avglen, $aft_r1_avglen, $aft_r2_avglen)."\n";
   }
}

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }

#************ Subroutines **************

## Subroutine to generate a custom wanted subroutine for File::Find find 
## passing a reference to an array
sub make_wanted {
    my $faref = shift;
    return sub {
        # Skip directories
        return if -d;
        # Check if file matches the pattern
        if ($_ =~ m/$fmask$/) {
            # Store the full path in the provided array
            push @$faref, $File::Find::name;
        }
    }
}

sub getUIDs {
   my ($mfiles, $uids)=@_;
   my $ns=scalar(@$mfiles);
   my $ilevel;
   my @dsplit = map { [ split(/\/+/) ] } @$mfiles;
   my @depths = map { scalar(@$_) } @dsplit;
   my $d=$depths[0];
   for (my $i=1;$i<=$#depths;$i++) {
     if ($depths[$i] ne $d) {
       die("Error: metadata files not found at the same level\n".
       $$mfiles[0]."\n".$$mfiles[$i]."\n");
     }
   }
   $d-=2;
   for(my $i=$d;$i>=0;$i--) {
     my %uk; #unique keys
     @$uids=();
     $ilevel=$d-$i+1;
     my $r=0;
     foreach my $dspl (@dsplit) {
       my $id=$$dspl[$i];
       $uk{$id}=1;
       $r++;
       if (scalar(keys(%uk))<$r) {
          @$uids=();
          undef $ilevel;
          last
       }
       push(@$uids, $id);
     }
     last if defined($ilevel);
   }
   die("Error: could not find unique IDs for the $ns metadata files!\n") unless defined($ilevel);
   die("Error: unique IDs (".scalar(@$uids)." and pathnames ($ns) do not match!\n") 
       unless $ns==scalar(@$uids);
   return $ilevel;
}



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

sub expErr {
 die("Error: not getting expected log line with $_[0]\n");
}
