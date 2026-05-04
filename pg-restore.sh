#!/usr/bin/env bash
set -Eeuo pipefail
umask 077

replace=0
restore_globals=0
args=()

die() { printf 'error: %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "missing command: $1"; }
sql_lit() { local s; s=$(printf '%s' "$1" | sed "s/'/''/g"); printf "'%s'" "$s"; }

usage() {
  printf 'usage: %s [--replace] [--globals] <backup-set-dir|dump-file> <target-db> [jobs]\n' "${0##*/}"
}

while (($#)); do
  case $1 in
    --replace) replace=1 ;;
    --globals) restore_globals=1 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; args+=( "$@" ); break ;;
    -*) die "unknown option: $1" ;;
    *) args+=( "$1" ) ;;
  esac
  shift
done

((${#args[@]} >= 2 && ${#args[@]} <= 3)) || { usage; exit 1; }
src=${args[0]}
target_db=${args[1]}
jobs=${args[2]:-4}
[[ $jobs =~ ^[0-9]+$ && $jobs -gt 0 ]] || die "jobs must be a positive integer"
export PGUSER=${PG_BACKUP_USER:-gpertea}

for cmd in psql pg_restore zstd sha256sum stat awk grep basename dirname createdb dropdb; do need "$cmd"; done

if [[ -d $src ]]; then
  ## a backup set is a directory named exactly like its contained dump.
  setdir=${src%/}
  base=$(basename "$setdir")
  dump=$setdir/$base.dump
elif [[ -f $src ]]; then
  ## direct dump paths are accepted, but still require nearby metadata.
  dump=$src
  setdir=$(dirname "$dump")
  base=${dump##*/}
  base=${base%.dump}
else
  die "backup set or dump not found: $src"
fi

[[ $base =~ ^pgbackup_[A-Za-z0-9_.-]+_p[0-9]+_[A-Za-z0-9_.-]+_[0-9]{8}_[0-9]{6}$ ]] || die "unrecognized backup name: $base"
## COMPLETE is written last by pg-backup.sh, after metadata and rsync input exist.
[[ -f $setdir/COMPLETE ]] || die "backup set is incomplete: missing COMPLETE"
[[ -f $dump ]] || die "dump file missing: $dump"
manifest=$setdir/$base.manifest.tsv
[[ -f $manifest ]] || die "manifest missing: $manifest"

## verify files recorded in the manifest before touching the target database.
while IFS=$'\t' read -r kind path size sum; do
  [[ $kind == file ]] || continue
  [[ -f $setdir/$path ]] || die "manifest file missing: $path"
  [[ $(stat -c '%s' "$setdir/$path") == "$size" ]] || die "size mismatch: $path"
  [[ $(sha256sum "$setdir/$path" | awk '{print $1}') == "$sum" ]] || die "hash mismatch: $path"
done < "$manifest"

pg_restore -l "$dump" >/dev/null || die "pg_restore cannot read dump: $dump"

exists=$(psql -d postgres -Atq -v ON_ERROR_STOP=1 -c "select 1 from pg_database where datname = $(sql_lit "$target_db")")
if [[ $exists == 1 ]]; then
  ## replacement is explicit because dropdb is irreversible.
  ((replace)) || die "target database exists; use --replace to drop it: $target_db"
  psql -d postgres -v ON_ERROR_STOP=1 -c "select pg_terminate_backend(pid) from pg_stat_activity where datname = $(sql_lit "$target_db")" >/dev/null
  dropdb "$target_db"
fi

if ((restore_globals)); then
  ## globals are optional because roles/tablespaces are often site-managed.
  globals=$setdir/$base.globals.sql.zst
  [[ -f $globals ]] || die "globals file missing: $globals"
  zstd -dc "$globals" | psql -d postgres -v ON_ERROR_STOP=1
fi

createdb -T template0 "$target_db"
pg_restore --exit-on-error --no-owner --no-privileges -j "$jobs" -d "$target_db" "$dump"

if [[ -d $setdir/config ]]; then
  printf 'config backups are in %s/config; review manually, do not auto-install them\n' "$setdir"
fi
printf 'restored %s into %s\n' "$dump" "$target_db"
