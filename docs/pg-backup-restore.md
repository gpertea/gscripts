# PostgreSQL backup and restore

The authoritative scripts are `pg-backup.sh` and `pg-restore.sh` in this repository.
They target PostgreSQL 17 on srv16. Backups are logical database dumps, with role,
database, extension, tablespace and server-configuration metadata.

## Scheduled use

```bash
/home/gpertea/gscripts/pg-backup.sh --if-changed genotypes rse sib storage_minder
```

No cron entry is installed by these scripts. Cron must supply a PATH containing
PostgreSQL clients, zstd, rsync, Python 3 and standard Linux utilities. Capture
stdout/stderr and check the exit status. Exit 0 includes an unchanged/skip result;
a failed secondary copy returns nonzero while retaining the complete local backup.

Change inspection requires this extension in the cluster's `postgres` database:

```sql
CREATE EXTENSION IF NOT EXISTS pg_walinspect WITH SCHEMA public;
```

`--if-changed` compares against the latest completed backup of each database:

- Cluster identity, timeline, server start time and database OID must match.
- Sequence values, effective settings and configuration files must match.
- WAL since the start of the previous dump must contain only explicitly allowed
  checkpoint, hint-page, heap-page pruning and running-transaction bookkeeping records.
- The existing backup must still pass manifest, size and SHA-256 verification.

All other WAL records, including commits, cause a backup. This is deliberately
cluster-wide: changes to one database can cause backups of the others. It also
covers transactions whose row changes preceded the previous dump but committed
later. Ordinary read-only queries do not themselves require a backup.

Missing/recycled WAL, an unavailable inspection function, a 30-second inspection
timeout, a restart, identity mismatch or unlogged tables cause a backup. Sequence
values are checked separately because cached/prelogged sequence changes do not
always generate new WAL. Statistics counters are not used. Maintenance can cause
extra backups; uncertainty never permits skipping. No replication slot or extra
WAL retention is needed. The first run after upgrading old backups takes a new
backup because old sets have no change-detection state.

This skips dump creation, but reads existing backup files to verify their hashes.
It does not compare every table's contents or promise to suppress every redundant
backup. External files and foreign-table data are outside a normal pg_dump backup.

## Backup controls

```bash
pg-backup.sh --check genotypes rse sib storage_minder
pg-backup.sh rse                         ## force a fresh dump
PG_BACKUP_SECONDARY='' pg-backup.sh rse  ## local backup only
```

Defaults: user `gpertea`, database `rse`, root `/data/backups/postgres`, compression
`zstd:9`, 20 completed local sets per database, secondary
`linwks34:/nfs/gdata/postgres`. Override with `PG_BACKUP_USER`, `PG_BACKUP_DB`,
`PG_BACKUP_ROOT`, `PG_BACKUP_COMPRESS`, `PG_BACKUP_KEEP`, `PG_BACKUP_SECONDARY`,
`PGHOST` and `PGPORT`. The secondary accepts one existing absolute local path or
`host:/absolute/path`; an empty value disables it.

A lock prevents overlapping jobs sharing the backup root. Local and secondary
sets are published by renaming staging directories after completion. Failed
copies are retried on later unchanged runs. Retention only removes complete sets
of the same database, after replication succeeds. Secondary retention is manual.

`PG_BACKUP_MIN_FREE_GB` defaults to 5. The script checks available space before
pg_dump and once per second while it runs, cancelling the dump below that floor.
This reduces disk-full risk; it is not a filesystem quota or a guarantee against
other processes consuming the remaining space. Failed partial local dumps are
removed on normal error/signal exit; SIGKILL or power loss can leave hidden partial
directories, which are never counted as completed backups.

## Restore

```bash
pg-restore.sh --verify-only /path/to/backup-set candidate
pg-restore.sh --precheck-only /path/to/backup-set candidate
pg-restore.sh /path/to/backup-set candidate 4
pg-restore.sh --replace /path/to/backup-set rse 4
```

Verification checks every recorded file and rejects incomplete manifests,
unrecorded payloads, unsafe paths and symlinks. Precheck additionally reads the
whole archive and checks target extensions/default versions and tablespaces.
It does not perform a trial SQL restore; object/role compatibility is established
by the staging restore when an actual restore is requested.

Actual restoration creates a randomly named staging database. Owners, grants,
locale, database settings and database privileges are preserved for new backups.
Legacy sets remain readable, but lack database-level metadata and may require
manual adjustment of template0 defaults after restoration.

`--replace` keeps the existing database online while staging is restored. At the
final switch it disables original connections, terminates clients and renames
both databases in one transaction. The original is retained under the printed
`pgrollback_*` name with connections disabled. Failed staging restores leave the
original untouched; failed switches restore its original connection setting on
normal error/signal exit. SIGKILL or a lost server connection during the switch
may require an administrator to re-enable connections manually. Retained originals
consume storage until explicitly removed by an administrator.

By default required roles and tablespaces must exist. `--globals` transactionally
reconciles role definitions and memberships from the new roles-only artifact,
including existing roles. This is a cluster-wide change: use it only when intended.
It does not create tablespaces; create those manually at appropriate target paths.
Legacy globals files require manual review/application. Configuration is archived
as a compressed TSV of absolute server paths and base64 file contents, alongside
effective settings; it is never automatically installed on the target server.

## Validation

`bash tests/pg-backup-restore.sh` starts a disposable local PostgreSQL cluster and
checks restore fidelity, replacement failure/success, manifests, change detection,
replication recovery, retention and the free-space guard. It does not use the live
server or production backup paths.

References: [PostgreSQL 17 WAL inspection](https://www.postgresql.org/docs/17/pgwalinspect.html)
and [pg_restore](https://www.postgresql.org/docs/17/app-pgrestore.html).
