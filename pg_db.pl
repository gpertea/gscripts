#!/usr/bin/perl
## a few basic routines to connect to the postgres server
use strict;
use DBI;
use DBD::Pg qw(:pg_types);
use Sys::Hostname;
my ($host)=split(/\./, hostname());
my $dbh; # global DB handle object

my $qry = "SELECT * FROM regions";
print STDERR '<'.$host.">\n";
my $srv = $host eq 'linwks34' ? 'localhost' : 'gdebsrv';
dbLogin($srv);
my ($sth, $r)=dbExec($qry);
while (my $row = dbFetch($sth)) {
  print join("\t", @$row)."\n";
}
$sth->finish(); # only needed if NOT all the rows were fetched

dbPrint("select * from datasets");

$dbh->disconnect();

#------------ subroutines -----
## see: https://docs.mojolicious.org/DBD/Pg

## methods common to all handles ($dbh, $sth):
## $h->err : current error code 
# 0  Empty query string
# 1  A command that returns no data successfully completed.
# 2  A command that returns data successfully completed.
# 3  A COPY OUT command is still in progress.
# 4  A COPY IN command is still in progress.
# 5  A bad response was received from the backend.
# 6  A nonfatal error occurred (a notice or warning message)
# 7  A fatal error was returned: the last query failed.

## $h->errstr : last error reported by Postgres

##  $h->state : 5 character state code from Postgres
#  00000 Successful completion
#  25P01 No active SQL transaction
#  25P02 In failed SQL transaction
#  S8006 Connection failure

sub dbErr {
 #my $dbLastError="@_\n";
 #print STDERR $_[0]."\n";
 #exit(1) unless defined($_[1]);
 die join("\n",@_)."\n";
}

#sub onErrExit {
# $dbExitSub=$_[0];
#}

# $dbh =  global database handler
sub dbLogin {
  my ($server, $db, $user)= @_;
  $server='localhost' unless $server;
  $db='rse' unless $db;
  $user='ruser' unless $user;
  open(PGPASS, "$ENV{HOME}/.pgpass") || die("Error opening $ENV{HOME}.pgpass\n");
  #hostname:port:database:username:password
  my ($pass);
  while (<PGPASS>) {
    chomp;
    next if m/^#/;
    my ($host,$port,$d,$u,$p)=split(/:/);
    if ($host eq $server && $d eq $db) {
       ($user,$pass)=($u, $p);
       last;
    }
  }
  close(PGPASS);
  die("Error: could not retrieve pass for user $user, db $db on $server\n") unless $pass;
  $dbh = DBI -> connect("dbi:Pg:dbname=$db;host=$server",  
                            $user, $pass,
                            {AutoCommit => 0, RaiseError => 1,
                             pg_server_prepare => 1 }
                         ) or die $DBI::errstr;
  # The AutoCommit attribute should always be explicitly set
}

sub dbQuery {
 #Execute a query and returns ALL the results as reference to an array 
 # of references to field value lists
 my ($query)=@_;
 my $aref=$dbh->selectall_arrayref($query) 
      || dbErr("Select all failed for:\n$query"); 
 
 return $aref;
}

sub dbPrep {
 my $sth = $dbh->prepare($_[0]);
 return $sth; 
}

sub dbExec { # execute non-query statement (update, insert)
  my ($req, @parm)=@_; # $req could be a query or a $sth
  my $sth;
  if (ref($req)) { # it's a sth
    $sth=$req;
  }
  else { #prep it first
    $sth=dbPrep($req);
  }
  ## execute should return the number of rows affected
  my $r=$sth->execute(@parm);
  dbErr(" *** execute failed! ".$DBI::errstr) if length($r)==0;
  return wantarray ? ($sth, $r) : $r;
}

sub dbDo {
 return dbExec(@_);
}

sub dbFetch {
  my ($sth)=@_;
   #given a sth for an executed statement
   #return an arrayref with row data, or undef if none left
  return $sth->fetchrow_arrayref();
}


sub dbFetchAll {
 # fetch all rows and return a ref to array of arrays
  my ($sth)=@_;
  return $sth->fetchall_arrayref();
}

sub dbRun {
  my $qry=$_[0];
  if ($qry=~m/\b(insert|update|delete|alter|create|drop)\b/i && 
      $qry!~m/\breturning\b/i) {
    return $dbh->do($qry); #do($q, \%attr, [@bind_values]);
  }
  return dbQuery($qry)
}


sub dbPrint {
 my ($q, $csep)=@_;
 $csep="\t" unless $csep;
 #if (ref($qry)) # pass a sth directly? nah..
 my ($s,$res)=dbExec($q);
 while (my $rd=dbFetch($s)) {
    print(join($csep, @$rd)."\n");
 }
}

sub dbDrop {
 return dbExec("drop table if exists $_[0]");
}
