#!/usr/bin/env perl
use strict;
use Getopt::Std;
use Sys::Hostname;
use FindBin;use lib $FindBin::Bin;
use dbPg;

my $usage = q/Usage:
  get_br_demo.pl [-o outfile.tab] brnums_list|all
  Takes a list of BrNum IDs (in a file or at stdin if '-') and outputs
  a tab delimited table with basic demographics as found in the database.
  Use "all" instead of a BrNum list to export all subjects and all populated
  subject fields, including subjects marked as dropped.
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
my $all_subjects=scalar(@ARGV)==1 && $ARGV[0] eq 'all';

if ($all_subjects) {
  # discover columns so future populated fields are exported automatically
  my $schema_cols=dbQuery(q/SELECT column_name FROM information_schema.columns
    WHERE table_schema='public' AND table_name='subjects'
    ORDER BY ordinal_position/);
  my (@headers, @select_cols);
  foreach my $colrow (@$schema_cols) {
    my $column=$colrow->[0];
    die("Error: unsafe subjects column name: $column\n")
      unless $column=~m/^[A-Za-z_][A-Za-z0-9_]*$/;
    my $has_value=dbQuery(qq/SELECT EXISTS
      (SELECT 1 FROM subjects WHERE "$column" IS NOT NULL)/)->[0]->[0];
    next unless $has_value;

    push(@headers, $column);
    push(@select_cols, qq/s."$column"/);
    if ($column eq 'dx_id') {
      push(@headers, qw{dx dx_name});
      push(@select_cols, 'd.dx', 'd.name AS dx_name');
    }
  }
  my $all_qry='SELECT '.join(', ', @select_cols).
    ' FROM subjects s JOIN dx d ON d.id=s.dx_id ORDER BY s.brint';

  print join("\t", @headers)."\n";
  my ($sth)=dbExec($all_qry);
  while (my $row=dbFetch($sth)) {
    print join("\t", map {
      !defined($_) ? '' : ref($_) ? '{'.join(',', @$_).'}' : $_
    } @$row)."\n";
  }
  dbLogout();

  if ($outfile) {
    select(STDOUT);
    close(OUTF);
  }
  exit 0;
}

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

