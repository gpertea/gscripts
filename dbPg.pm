package dbPg;
## a few basic routines to connect to the postgres server
use strict;
use Exporter;
use DBI;
use DBD::Pg qw(:pg_types);

=head1 NAME

dbSession - a DBI wrapper with various db helpers and utils

=cut

=head1 SYNOPSIS

  use dbPg;

=head1 DESCRIPTION

A wrapper module for DBI, with various helper subroutines

=cut

our ($VERSION, @ISA, @EXPORT);
@ISA = qw(Exporter);
@EXPORT = qw( dbErr dbLogin dbQuery dbPrep dbFetch dbFetchAll 
              dbExec dbExecPrep dbFetchPrep dbDo dbRun dbPrint dbDrop dbLogout);

#------------ subroutines -----
## see: https://docs.mojolicious.org/DBD/Pg

## methods common to all handles:
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

our ($pg_dbh_, $pg_sth_); # last dbh, sth generated
our ($file_log, $stdout_log, $stderr_log);
local *__OLDERR;
local *__OLDSTD;

return 1;

sub dbErr {
 #my $dbLastError="@_\n";
 #print STDERR $_[0]."\n";
 #exit(1) unless defined($_[1]);
 die join("\n",@_)."\n";
}

#sub onErrExit {
# $dbExitSub=$_[0];
#}

sub parseDbSpec { #parses: [user@]server[:port]/database
  my $cstr=shift(@_);
  my ($server, $db, $user);
  if (index($cstr, '@')>=0) {
    ($user)=(m/^([^@]+)@/);
    $cstr=~s/^([^@]+)//;
  }
  if (index($cstr, '/')>=0) {
    ($server, $db)=split(/\//, $cstr);
  } else {  $db=$cstr }
  return ($server, $db, $user);
}

# $pg_dbh_ =  global database handler
sub dbLogin {  
  my ($server, $db, $user)= @_;
  if (index($server,'/')>=0 && length($db)==0) {
    ($server, $db, $user)=parseDbSpec($server);
  }
  die("Error database required!\n") unless $db;

  $server='localhost' unless $server;  
  $user='ruser' unless $user;
  my $port=5432;
  my @sp=split(/:/, $server);
  ($server, $port)=@sp[0,1] if (@sp>1);
  #hostname:port:database:username:password
  my ($pass);
  open(PGPASS, "$ENV{HOME}/.pgpass") || die("Error opening $ENV{HOME}.pgpass\n");
  while (<PGPASS>) {
    chomp;
    next if m/^#/;
    my ($host,$port,$d,$u,$p)=split(/:/);
    if ($host eq $server && $user eq $u) {
       $pass=$p;
       last;
    }
  }
  close(PGPASS);
  die("Error: could not retrieve pass for user $user on $server\n") unless $pass;
  $pg_dbh_ = DBI -> connect("dbi:Pg:dbname=$db;host=$server;port=$port",
                            $user, $pass,
                            {AutoCommit => 0, RaiseError => 1,
                             pg_server_prepare => 1 }
                         ) or die $DBI::errstr;
  # The AutoCommit attribute should always be explicitly set
  #return ($server.':'.$port, $db, $user);
  return $pg_dbh_;
}

sub dbQuery {
 #Execute a query and returns ALL the results as reference to an array 
 # of references to field value lists
 my ($query)=@_;
 my $aref=$pg_dbh_->selectall_arrayref($query) 
      || dbErr("Select all failed for:\n$query"); 
 return $aref;
}

sub dbPrep {
 if ($pg_sth_) {
  $pg_sth_->finish();
  $pg_sth_=0;
 }
 $pg_sth_ = $pg_dbh_->prepare($_[0]);
 return $pg_sth_; 
}

sub dbExec { # execute query or prepared statement
  my ($req, @parm)=@_; # $req could be a query or a $sth
  my $sth;
  if (ref($req)) { # it's a sth
    $sth=$req;
  }
  else { #prep it first
    dbPrep($req);
    $sth=$pg_sth_;
  }
  ## execute should return the number of rows affected
  my $r=$sth->execute(@parm);
  dbErr(" *** execute failed! ".$DBI::errstr) if length($r)==0;
  return wantarray ? ($sth, $r) : $r;
}

sub dbExecPrep { #execute prepared statement $pg_sth_
 my $r=$pg_sth_->execute(@_);
 dbErr(" *** execute failed! ".$DBI::errstr) if length($r)==0;
 return wantarray ? ($pg_sth_, $r) : $r; 
}

sub dbFetchPrep { #execute and fetch prepared statement in $pg_sth_
  my $r=$pg_sth_->execute(@_);
  dbErr(" *** execute failed! ".$DBI::errstr) if length($r)==0;
  ## fetch all rows and return a ref to array of arrays
  return $pg_sth_->fetchall_arrayref();
}

sub dbDo {
 return dbExec(@_);
}

sub dbFetch {
  my $sth=shift(@_) || $pg_sth_;
  #given a sth for an executed statement
  #return an arrayref with row data, or undef if none left
  return $sth->fetchrow_arrayref();
}

sub dbFetchAll {
  # fetch all rows and return a ref to array of arrays
  my $sth = shift(@_) || $pg_sth_;
  return $sth->fetchall_arrayref();
}

sub dbRun {
  my $qry=$_[0];
  if ($qry=~m/\b(insert|update|delete|alter|create|drop)\b/i && 
      $qry!~m/\breturning\b/i) {
    return $pg_dbh_->do($qry); #do($q, \%attr, [@bind_values]);
  }
  return dbQuery($qry)
}


sub dbPrint {
 my ($q, $csep)=@_;
 $csep="\t" unless $csep;
 #if (ref($qry)) # pass a sth directly? nah..
 my ($s,$res)=dbExec($q);
 while (my $rd=dbFetch($s)) {
    print(join($csep, map {ref ? '{'.join(',',@$_).'}' : $_} @$rd)."\n");
 }
}

sub dbDrop {
 return dbExec("drop table if exists $_[0]");
}

sub dbLogout {
 if ($pg_sth_) {
  $pg_sth_->finish();
  $pg_sth_=0;
 }
 $pg_dbh_->disconnect();
 $pg_dbh_=0;
 $pg_sth_=0;
}
