use strict;
use warnings;
use Test::More;
use FindBin;
use lib "$FindBin::Bin/..";
use dbPg;

is_deeply([dbPg::parseDbSpec('alice@server:5444/rse')],['server:5444','rse','alice'],'user and port parsing');
is_deeply([dbPg::parseDbSpec('server/rse')],['server','rse',undef],'legacy host/database');
is_deeply([dbPg::parseDbSpec('rse')],[undef,'rse',undef],'database-only form');
{
    no warnings 'redefine';
    my @captured;
    local *DBI::connect=sub { @captured=@_; return bless {},'FakeDB' };
    dbLogin('host:5444','mydb','legacy');
    is($captured[1],'dbi:Pg:dbname=mydb;host=host;port=5444','legacy DSN');
    ok(!defined($captured[3]),'password selection delegated to libpq');
    is($captured[4]{AutoCommit},0,'legacy transaction default preserved');
    dbLogin({host=>'/var/run/postgresql',dbname=>'rse',user=>'gpertea',read_only=>1,autocommit=>1});
    like($captured[1],qr/host=\/var\/run\/postgresql/,'socket directory preserved');
    like($captured[1],qr/options=--default_transaction_read_only=on/,'read-only enforced at connection');
    is($captured[4]{AutoCommit},1,'explicit autocommit');
    dbLogin('alice@host/rse');
    is($captured[2],'alice','connection string user honored');
    like($captured[1],qr/port=5432/,'legacy default port preserved');
}
if ($ENV{DBPG_LIVE_TEST}) {
    dbLogin({host=>'/var/run/postgresql',dbname=>'rse',user=>'gpertea',read_only=>1,autocommit=>1});
    is(dbQuery('SELECT ?::integer',17)->[0][0],17,'live bound query');
    is(dbQueryHash('SELECT ?::text AS value',q{x' OR true --})->[0]{value},q{x' OR true --},'placeholder preserves literal SQL-like data');
    is(dbQuery(q{SHOW default_transaction_read_only})->[0][0],'on','live server enforces read-only');
    dbLogout();
}
done_testing();
