#!/usr/bin/perl
use strict;
my $usage=q{
  db_filter_sql.pl tbl_list.txt schema.sql

  Filters schema.sql, which is expected to include CREATE TABLE
  ALTER TABLE or CREATE SEQUENCE etc. to only let through SQL statements
  relating to the tables in tbl_list.txt
};

foreach my $f (@ARGV) {
  die("Error: file $f not found!\n") unless -f $f;
}

my ($ftbl, $fsql)=@ARGV;

my %th; # table => 1

open(F, $ftbl) || die("Error opening $ftbl!\n");
while(<F>) {
 chomp;
 my @t=split();
 $th{lc($t[0])}=1
}
close(F);

open(SQL, $fsql) || die("Error opening $fsql!\n");
my $ml=0; #multiline flag
my $k=0; #keep flag
while(<SQL>) {
  # check for multi-line SQL to preserve
  if ($ml) { #multi-line entry
    if ($k) {
       print $_;
    }
    if (m/\);\s*$/) {
     $ml=0;
     $k=0;
    }
    next;
  }
  if (m/CREATE\s+TABLE\s+"public"\."([^"]+)/) {
    my $tbl=lc($1);
    $ml=1;
    if ($th{$tbl}) {
      print $_;
      $k=1;
    }
    next;
  }
  if (m/ALTER\s+TABLE\s+"public"\."([^"]+)/ || 
         m/OWNED BY "public"\."([^"]+)/ ) {
    my $tbl=lc($1);
    if ($th{$tbl}) {
      print $_;
    }
    next;
  }
}
close(SQL);
