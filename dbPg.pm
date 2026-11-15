package dbPg;
## a few basic routines to connect to the postgres server
use strict;
use Exporter;
use DBI;
use DBD::Pg qw(:pg_types);

=head1 NAME

dbPg - a DBI wrapper with various db helpers and utils

=cut

=head1 SYNOPSIS

  use dbPg;
  dbLogin({host=>'/var/run/postgresql',dbname=>'rse',user=>'gpertea',
           autocommit=>1,read_only=>1});
  my $rows=dbQueryHash('SELECT id FROM subjects WHERE brnum=?','Br0001');

=head1 DESCRIPTION

A wrapper module for DBI, with various helper subroutines

=cut

our ($VERSION, @ISA, @EXPORT);
$VERSION='0.2';
@ISA = qw(Exporter);
@EXPORT = qw( dbErr dbLogin dbQuery dbPrep dbFetch dbFetchAll 
              dbQueryHash dbExec dbExecPrep dbFetchPrep dbDo dbRun dbPrint dbDrop dbLogout);

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

sub parseDbSpec { ## [user@]host[:port]/database; socket paths use hash options
  my ($spec)=@_;
  my ($user,$server,$db);
  $user=$1 if $spec=~s/^([^@\/]+)@//;
  if ($spec=~m{^([^/]*)/([^/]+)$}) { ($server,$db)=($1,$2) }
  else { $db=$spec }
  return ($server,$db,$user);
}

sub dbLogin {
  my ($server,$db,$user)=@_;
  my %opt;
  if (ref($server) eq 'HASH') {
    %opt=%$server;
    ($server,$db,$user)=@opt{qw(host dbname user)};
    $db=$opt{database} unless defined($db);
  } else {
    if (@_==1) { ($server,$db,$user)=parseDbSpec($server) }
    $server='localhost' unless defined($server) && length($server);
    $user='ruser' unless defined($user) && length($user);
    $opt{autocommit}=0;
    $opt{port}=5432;
  }
  die "Error: database required!\n" unless defined($db) && length($db);
  my $port=$opt{port};
  ## only TCP host:port is split; Unix socket directory names stay intact
  if (defined($server) && $server=~s/^([^:\/]+):(\d+)$/$1/) { $port=$2 }
  my %dsn=(dbname=>$db,host=>$server,port=>$port,service=>$opt{service});
  $dsn{options}='--default_transaction_read_only=on' if $opt{read_only};
  my @parts;
  for my $key (qw(dbname host port service options)) {
    next unless defined($dsn{$key}) && length($dsn{$key});
    die "Invalid connection parameter $key\n" if $dsn{$key}=~/[;\r\n\0]/;
    push @parts,"$key=$dsn{$key}";
  }
  ## libpq handles peer auth, PGPASSFILE, escaped entries and all matching fields
  $pg_dbh_=DBI->connect('dbi:Pg:'.join(';',@parts),$user,$opt{password},
      {AutoCommit=>($opt{autocommit}//0),RaiseError=>1,PrintError=>0,
       pg_server_prepare=>1,pg_enable_utf8=>1}) or die $DBI::errstr;
  return $pg_dbh_;
}

sub dbQuery {
  my ($query,@bind)=@_;
  return $pg_dbh_->selectall_arrayref($query,undef,@bind)
      || dbErr("Select all failed");
}

sub dbQueryHash {
  my ($query,@bind)=@_;
  return $pg_dbh_->selectall_arrayref($query,{Slice=>{}},@bind)
      || dbErr("Select hashes failed");
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
