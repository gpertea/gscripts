# PostgreSQL backup and restore

Base directory: `~/gscripts/`. Sources: `pg-backup.sh`, `pg-backup-wal.py` and
`pg-restore.sh`, at the same Git revision as this document.
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

`--if-changed` considers each database separately against its latest completed backup:

- Cluster identity, timeline, server start time and database OID must match.
- Sequence values, effective settings, configuration files and shared globals must match.
- WAL block references identify the affected database by its OID. Ordinary writes
  to another database do not trigger a backup of this one.
- Parent/child transaction IDs connect commits, including prepared commits, to
  the databases whose data was written. Shared dependency records follow the
  database of their transaction, rather than marking every database changed.
- Known physical-maintenance records are excluded. A commit with no associated
  writes is ignored only when a saved or in-range XID allocation boundary proves
  that its writes could not predate the inspected WAL interval.
- The existing backup must still pass manifest, size and SHA-256 verification.

For example, inserts into `storage_minder` do not require another `genotypes`,
`rse` or `sib` dump. Each decision prints its reason. Existing September 13 state
archives remain usable; a fresh full baseline is not required solely for this
upgrade. New backups also reserve and save a transaction ID after the WAL marker;
this allocation boundary includes older active/prepared transactions correctly.
It does not modify database rows.

If a transaction wrote before the previous backup's WAL boundary and committed
later, its database may not be assignable from the available records. That case
requires a backup, with an explicit uncertainty reason; commits are never blindly
ignored. Shared role/tablespace or configuration changes can affect multiple
databases and still require backups. Other shared-catalog changes, unknown WAL
record formats, missing/recycled WAL, an inspection timeout, a restart, identity
mismatch or unlogged tables also require a backup. The decoder targets PostgreSQL
17 and fails conservatively on other major versions.

Sequence values are checked separately because cached/prelogged changes do not
always generate new WAL. Statistics counters are not used. No replication slot
or additional WAL retention is needed. Archives without any saved detection state
require a first full backup.

```bash
pg-backup.sh --check-changes genotypes rse sib storage_minder
```

`--check-changes` reports selection decisions without dumping, copying or pruning.
It reads saved metadata and WAL but does not hash all payload files; actual skip
operations perform that additional integrity check. It uses the normal backup-root
lock and removes its temporary state files on exit.

This skips dump creation, but reads existing backup files to verify their hashes.
It does not compare every table's contents or promise to suppress every redundant
backup. External files and foreign-table data are outside a normal pg_dump backup.

## Backup controls

```bash
pg-backup.sh --check genotypes rse sib storage_minder
pg-backup.sh rse                         ## force a fresh dump
PG_BACKUP_SECONDARY='' pg-backup.sh rse  ## local backup only
```

Keep `pg-backup-wal.py` alongside `pg-backup.sh` when deploying; change-detection
modes stop if the helper is missing.

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

`python3 tests/pg-backup-wal.py` checks WAL attribution and transaction edge cases.
`bash tests/pg-backup-restore.sh` starts a disposable local PostgreSQL cluster and
checks restore fidelity, replacement failure/success, manifests, change detection,
replication recovery, retention, the free-space guard, isolation between databases,
shared dependencies, subtransactions and commits across the backup boundary. It does not use the live
server or production backup paths.

References: [PostgreSQL 17 WAL inspection](https://www.postgresql.org/docs/17/pgwalinspect.html)
and [pg_restore](https://www.postgresql.org/docs/17/app-pgrestore.html).
