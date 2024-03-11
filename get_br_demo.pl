#!/usr/bin/env perl
use strict;
use Getopt::Std;
use Sys::Hostname;
use FindBin;use lib $FindBin::Bin;
use dbPg;

my $usage = q/Usage:
  get_br_demo.pl [-o outfile.tab] brnums_list
  Takes a list of BrNum IDs (in a file or at stdin if '-') and outputs
  a tab delimited table with basic demographics as found in the database.
/;
umask 0002;
getopts('ho:') || die($usage."\n");
die($usage) if $Getopt::Std::opt_h || scalar(@ARGV)==0;
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
  open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
  select(OUTF);
  }
# --

my ($host)=split(/\./, hostname());
my $srv='srv16:5432';
dbLogin($srv, 'rse', 'ruser');
my $qry = q/SELECT brnum, dx, sex, race, age, pmi, mod
  FROM subjects s, dx d where brint=? and d.id=s.dx_id/;
dbPrep($qry);
my @notfound;
print join("\t", qw{ BrNum Dx Sex Race Age PMI MoD })."\n";
while (<>) {
 my @brnums=(m/(Br\d+)/g);
 next unless @brnums>0;
 #print STDERR " got brnums: ",join(", ", @brnums)."\n";
 foreach my $b (@brnums) {
    my ($brint)=($b=~m/(\d+)/);
    $brint=int($brint);
    my $res = dbExecPrep( ($brint) );
    my $found=0;
    while (my $row = dbFetch()) {
      $found=1;
      print join("\t", (map {
           ref ? '{'.join(',',@$_).'}' : $_
          } @$row) )."\n";
    } #while rows
    push(@notfound, $b) unless $found;
 } #for each BrNum
}
if (@notfound>0) {
  print STDERR " ~~~ Warning: ",scalar(@notfound)," BrNum IDs not found in database:\n",
     join(', ',@notfound)."\n";
}

dbLogout();

# --
if ($outfile) {
 select(STDOUT);
 close(OUTF);
 }


#dbPrint($qry);
#dbLogout();



