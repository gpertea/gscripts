# dbPg database helper

Source: `dbPg.pm`, version 0.2. Paths are relative to the gscripts repository.

Legacy calls keep their connection and transaction defaults:

```perl
use dbPg;
dbLogin('localhost', 'rse', 'ruser');
my $rows=dbQuery('SELECT id FROM subjects WHERE brnum=?', 'Br1536');
dbLogout();
```

Explicit options support local peer authentication and read-only utilities:

```perl
dbLogin({host=>'/var/run/postgresql', dbname=>'rse', user=>'gpertea',
         autocommit=>1, read_only=>1});
my $rows=dbQueryHash('SELECT id,brnum FROM subjects WHERE brnum=?', 'Br1536');
dbLogout();
```

`dbQuery` returns array rows; `dbQueryHash` returns hash rows. Both bind extra
arguments through DBI. Existing prepared-statement helpers remain available.
`parseDbSpec` handles `[user@]host[:port]/database`; Unix socket paths use the
explicit options form. Legacy AutoCommit remains off. `read_only` configures
PostgreSQL at connection startup, before any application statement.

libpq handles peer authentication, password files, services and environment
settings. The helper no longer manually selects a password by host/user while
ignoring database/port. A password can be supplied explicitly in the options
hash; credentials must not be printed or stored in project files.

```bash
prove tests/dbPg.t
DBPG_LIVE_TEST=1 prove tests/dbPg.t
```

The optional live checks only query rse, testing placeholders and the server's
read-only setting. Other tests use a mocked DBI connection and need no
credentials. Workflow usage is documented in the RNA genotype collection at
`/dbdata/cdb/genotyping/imputation_pipeline/rna-check/docs/brain-genotype.md`.
