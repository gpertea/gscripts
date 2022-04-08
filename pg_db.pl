#!/usr/bin/env perl
## a few basic routines to connect to the postgres server
use strict;
use Sys::Hostname;
use FindBin;use lib $FindBin::Bin;
use dbPg;

my ($host)=split(/\./, hostname());

my $qry = "SELECT * FROM regions";
#my $srv = $host eq 'linwks34' ? 'localhost' : 'gdebsrv';
my $srv='localhost:5432';
dbLogin($srv);
dbPrint($qry);
dbLogout();
exit;

dbExec($qry);
while (my $row = dbFetch()) {
  print join("\t", (map {
       ref ? '{'.join(',',@$_).'}' : $_
      } @$row) )."\n";
}

#dbPrint("select * from datasets");

dbLogout();


