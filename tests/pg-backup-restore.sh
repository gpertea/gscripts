#!/usr/bin/env bash
set -Eeuo pipefail
repo=$(cd -- "$(dirname -- "$0")/.." && pwd -P)
fixture=$(mktemp -d)
PATH="$(pg_config --bindir):$PATH"
export PATH
export PGHOST=$fixture PGPORT=5432 PGUSER
PGUSER=$(id -un)
export PG_BACKUP_USER=$PGUSER PG_BACKUP_ROOT=$fixture/backups PG_BACKUP_SECONDARY=$fixture/secondary
export PG_BACKUP_COMPRESS=zstd:1
unset PGDATABASE PGSERVICE PGOPTIONS
cleanup() {
  pg_ctl -D "$fixture/data" -m immediate stop >/dev/null 2>&1 || true
  rm -rf -- "$fixture"
}
trap cleanup EXIT
fail() { printf 'FAIL: %s\n' "$*" >&2; exit 1; }
q() { psql -X -w -Atq -v ON_ERROR_STOP=1 "$@"; }
count() { find "$PG_BACKUP_ROOT" -mindepth 2 -maxdepth 2 -name COMPLETE | wc -l; }
latest() { find "$PG_BACKUP_ROOT" -mindepth 2 -maxdepth 2 -name COMPLETE -printf '%h\n' | sort | tail -n 1; }
backup() { sleep 1; "$repo/pg-backup.sh" "$@" > "$fixture/backup.log" 2>&1 || { cat "$fixture/backup.log"; return 1; }; }
initdb -D "$fixture/data" -A trust --no-locale -E UTF8 > "$fixture/init.log"
pg_ctl -D "$fixture/data" -l "$fixture/server.log" -o "-c listen_addresses='' -c unix_socket_directories='$fixture' -c autovacuum=off -c max_prepared_transactions=10" start >/dev/null
mkdir "$PG_BACKUP_SECONDARY"
q -d postgres -c 'CREATE EXTENSION pg_walinspect; CREATE ROLE backup_owner; CREATE ROLE backup_reader; CREATE ROLE role_from_backup;'
createdb sample
q -d sample <<'SQL'
CREATE TABLE items(id serial PRIMARY KEY, value text);
INSERT INTO items(value) VALUES ('one'),('two');
ALTER TABLE items OWNER TO backup_owner;
GRANT SELECT ON items TO backup_reader;
ALTER DATABASE sample OWNER TO backup_owner;
ALTER DATABASE sample SET work_mem='8MB';
REVOKE CONNECT ON DATABASE sample FROM PUBLIC;
GRANT CONNECT ON DATABASE sample TO backup_reader;
SQL
"$repo/pg-backup.sh" --check sample > "$fixture/check.log"
backup sample
first=$(latest)
[[ $(count) == 1 && -f $PG_BACKUP_SECONDARY/${first##*/}/COMPLETE ]] || fail publication
"$repo/pg-restore.sh" --precheck-only "$first" restored > "$fixture/precheck.log"
"$repo/pg-restore.sh" "$first" restored > "$fixture/restore.log"
[[ $(q -d restored -c 'select count(*) from items') == 2 ]] || fail data
[[ $(q -d restored -c "select tableowner from pg_tables where tablename='items'") == backup_owner ]] || fail ownership
[[ $(q -d restored -c "select has_table_privilege('backup_reader','items','SELECT')") == t ]] || fail grants
[[ $(q -d restored -c 'show work_mem') == 8MB ]] || fail database_settings
[[ $(q -d postgres -c "select has_database_privilege('backup_reader','restored','CONNECT')") == t ]] || fail database_grants
printf 'PASS: restore preserves data, owners, grants and database settings\n'

## restore itself wrote WAL, so establish a fresh baseline before testing a quiet cluster.
backup sample
before=$(count)
backup --if-changed sample
[[ $(count) == "$before" ]] || { cat "$fixture/backup.log"; fail unchanged; }
printf 'PASS: unchanged database skips without creating a backup\n'
## losing a replication receipt after publication must not create another dump.
rm -f -- "$PG_BACKUP_ROOT/.replicated-$(basename "$(latest)")"
backup --if-changed sample
[[ $(count) == "$before" ]] || fail lost_receipt
## statistics resets cannot hide changes; change detection does not trust these counters.
q -d sample -c 'select pg_stat_reset()' >/dev/null
backup --if-changed sample
[[ $(count) == "$before" ]] || fail stats_reset
checkpoint_lsn=$(zstd -dc "$(latest)/$(basename "$(latest)").state.tar.zst" | tar -xOf - lsn)
q -d postgres -c 'CHECKPOINT' >/dev/null
backup --if-changed sample
[[ $(count) == "$before" ]] || { cat "$fixture/backup.log"; q -d postgres -c "select resource_manager,record_type,count(*) from pg_get_wal_records_info('$checkpoint_lsn',pg_current_wal_lsn()) group by 1,2"; fail checkpoint_skip; }
printf 'PASS: lost receipt recovery, statistics resets and checkpoints do not cause duplicate backups\n'
q -d sample -c "select nextval('items_id_seq')" >/dev/null
backup --if-changed sample
[[ $(count) == $((before+1)) ]] || fail sequence_change
q -d sample -c "INSERT INTO items(value) VALUES ('three')" >/dev/null
backup --if-changed sample
[[ $(count) == $((before+2)) ]] || fail row_change
q -d sample -c 'ALTER TABLE items ADD COLUMN extra integer' >/dev/null
backup --if-changed sample
[[ $(count) == $((before+3)) ]] || fail schema_change
printf '\n## test configuration fingerprint\n' >> "$fixture/data/postgresql.conf"
backup --if-changed sample
[[ $(count) == $((before+4)) ]] || fail config_change
printf 'PASS: sequence, row, schema and configuration changes require backups\n'
q -d postgres -c 'DROP EXTENSION pg_walinspect' >/dev/null
backup --if-changed sample
[[ $(count) == $((before+5)) ]] || fail unavailable_wal_inspection
q -d postgres -c 'CREATE EXTENSION pg_walinspect' >/dev/null
printf 'PASS: unavailable WAL inspection forces a backup\n'

## an injected restore failure must leave the existing target intact and remove staging data.
q -d restored -c 'CREATE TABLE sentinel(id integer)' >/dev/null
mkdir "$fixture/bin"
export REAL_PG_RESTORE
REAL_PG_RESTORE=$(command -v pg_restore)
cat > "$fixture/bin/pg_restore" <<'SH'
#!/usr/bin/env bash
for arg in "$@"; do [[ $arg != --exit-on-error ]] || exit 42; done
exec "$REAL_PG_RESTORE" "$@"
SH
chmod +x "$fixture/bin/pg_restore"
if PATH="$fixture/bin:$PATH" "$repo/pg-restore.sh" --replace "$(latest)" restored > "$fixture/failure.log" 2>&1; then fail injected_failure; fi
[[ $(q -d restored -c "select to_regclass('sentinel')") == sentinel ]] || fail original_lost
[[ $(q -d postgres -c "select count(*) from pg_database where datname like 'pgrestore_%'") == 0 ]] || fail staging_leak
## prepared transactions prevent the final rename; connection access must be restored on failure.
q -d restored -c "BEGIN; INSERT INTO sentinel VALUES (1); PREPARE TRANSACTION 'restore_blocker'" >/dev/null
if "$repo/pg-restore.sh" --replace "$(latest)" restored > "$fixture/switch-failure.log" 2>&1; then fail blocked_switch; fi
[[ $(q -d postgres -c "select datallowconn from pg_database where datname='restored'") == t ]] || { cat "$fixture/switch-failure.log"; fail connection_recovery; }
[[ $(q -d restored -c "select to_regclass('sentinel')") == sentinel ]] || fail switch_original_lost
q -d restored -c "ROLLBACK PREPARED 'restore_blocker'" >/dev/null
"$repo/pg-restore.sh" --replace "$(latest)" restored > "$fixture/replace.log"
rollback=$(q -d postgres -c "select datname from pg_database where datname like 'pgrollback_%'")
[[ -n $rollback && $(q -d restored -c 'select count(*) from items') == 3 ]] || fail replacement
[[ $(q -d postgres -c "select datallowconn from pg_database where datname='$rollback'") == f ]] || fail rollback_access
q -d postgres -c "ALTER DATABASE \"$rollback\" ALLOW_CONNECTIONS true" >/dev/null
[[ $(q -d "$rollback" -c "select to_regclass('sentinel')") == sentinel ]] || fail rollback_data
printf 'PASS: failed replacement preserves original; successful replacement retains rollback database\n'

q -d postgres -c 'DROP ROLE role_from_backup' >/dev/null
"$repo/pg-restore.sh" --globals "$(latest)" globals_restore > "$fixture/globals.log"
[[ $(q -d postgres -c "select count(*) from pg_roles where rolname='role_from_backup'") == 1 ]] || fail globals
printf 'PASS: transactional role restore tolerates existing roles and restores missing ones\n'

mkdir "$fixture/corrupt"
cp -a "$(latest)" "$fixture/corrupt/"
bad=$fixture/corrupt/$(basename "$(latest)")
printf '\n' >> "$bad/${bad##*/}.dump"
if "$repo/pg-restore.sh" --replace "$bad" restored > "$fixture/corrupt.log" 2>&1; then fail corrupt_backup_accepted; fi
[[ $(q -d restored -c 'select count(*) from items') == 3 ]] || fail corrupt_touched_target
printf 'PASS: corrupt backup rejected before target changes\n'
cp -- "$(latest)/$(basename "$(latest)").dump" "$bad/${bad##*/}.dump"
printf 'kind\tpath\tsize\tsha256\n' > "$bad/${bad##*/}.manifest.tsv"
if "$repo/pg-restore.sh" --verify-only "$bad" restored > "$fixture/empty-manifest.log" 2>&1; then fail empty_manifest; fi
printf 'PASS: empty manifest cannot bypass verification\n'

export PG_BACKUP_SECONDARY=$fixture/retry
if backup sample; then fail missing_secondary_accepted; fi
last=$(latest); before=$(count)
mkdir "$PG_BACKUP_SECONDARY"
backup --if-changed sample
[[ $(count) == "$before" && -f $PG_BACKUP_SECONDARY/${last##*/}/COMPLETE ]] || fail secondary_retry
printf 'PASS: failed secondary copy retains local backup and retries without a new dump\n'

q -d sample -c 'CREATE UNLOGGED TABLE transient(id integer)' >/dev/null
backup sample
before=$(count)
backup --if-changed sample
[[ $(count) == $((before+1)) ]] || fail unlogged_skip
printf 'PASS: unlogged tables disable skipping\n'
q -d sample -c 'DROP TABLE transient' >/dev/null
backup sample
before=$(count)
pg_ctl -D "$fixture/data" restart -m fast -l "$fixture/server.log" >/dev/null
backup --if-changed sample
[[ $(count) == $((before+1)) ]] || fail restart_skip
printf 'PASS: server restart forces a new backup\n'
mkdir "$PG_BACKUP_ROOT/pgbackup_$(hostname -s)_p5432_sample_20990101_000000"
mkdir "$PG_BACKUP_ROOT/20990101_000000"
printf 'complete\t2099\n' > "$PG_BACKUP_ROOT/20990101_000000/COMPLETE"
PG_BACKUP_KEEP=2 backup sample
[[ $(count) == 3 && -f $PG_BACKUP_ROOT/20990101_000000/COMPLETE && -d $PG_BACKUP_ROOT/pgbackup_$(hostname -s)_p5432_sample_20990101_000000 ]] || fail retention
if PG_BACKUP_KEEP=0 backup sample; then fail invalid_retention; fi
printf 'PASS: retention ignores incomplete sets and rejects zero retention\n'
before=$(count)
if PG_BACKUP_MIN_FREE_GB=99999 backup sample; then fail free_space_guard; fi
[[ $(count) == "$before" ]] || fail guard_published_partial
printf 'PASS: free-space guard rejects backup without publishing a partial set\n'
printf 'All backup/restore integration checks passed.\n'
