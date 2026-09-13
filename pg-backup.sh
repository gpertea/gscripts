#!/usr/bin/env bash
set -Eeuo pipefail
umask 077

backup_root=${PG_BACKUP_ROOT:-/data/backups/postgres}
secondary=${PG_BACKUP_SECONDARY-linwks34:/nfs/gdata/postgres}
keep=${PG_BACKUP_KEEP:-20}
compress=${PG_BACKUP_COMPRESS:-zstd:9}
min_free_gb=${PG_BACKUP_MIN_FREE_GB:-5}
dump_pid=''
check_only=0
if_changed=0
dbs=()
work=''
failed=0
export PGUSER=${PG_BACKUP_USER:-gpertea} PGCONNECT_TIMEOUT=${PGCONNECT_TIMEOUT:-10}

usage() {
  printf 'usage: %s [--check] [--if-changed] [db ...]\n' "${0##*/}"
  printf '%s\n' 'env: PG_BACKUP_USER, PGHOST, PGPORT, PG_BACKUP_ROOT, PG_BACKUP_DB, PG_BACKUP_KEEP,' \
    '     PG_BACKUP_COMPRESS, PG_BACKUP_MIN_FREE_GB (default 5),' \
    '     PG_BACKUP_SECONDARY (empty disables the secondary copy)' \
    '--if-changed skips only when WAL, sequence and configuration checks prove no changes.'
}
die() { printf 'error: %s\n' "$*" >&2; exit 1; }
warn() { printf 'warning: %s\n' "$*" >&2; }
psqlq() { psql -X -w -Atq -F $'\t' -v ON_ERROR_STOP=1 "$@"; }
trap 'if [[ -n $dump_pid ]]; then kill "$dump_pid" 2>/dev/null || true; wait "$dump_pid" 2>/dev/null || true; fi; [[ -z $work ]] || rm -rf -- "$work"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

while (($#)); do
  case $1 in
    --check) check_only=1 ;;
    --if-changed) if_changed=1 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; dbs+=( "$@" ); break ;;
    -*) die "unknown option: $1" ;;
    *) dbs+=( "$1" ) ;;
  esac
  shift
done
((${#dbs[@]})) || dbs=( "${PG_BACKUP_DB:-rse}" )
[[ $keep =~ ^[1-9][0-9]{0,5}$ ]] || die 'PG_BACKUP_KEEP must be between 1 and 999999'
[[ $min_free_gb =~ ^[1-9][0-9]{0,4}$ ]] || die 'PG_BACKUP_MIN_FREE_GB must be between 1 and 99999'
min_free=$((min_free_gb * 1024 * 1024 * 1024))
for db in "${dbs[@]}"; do
  [[ $db =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]{0,62}$ ]] || die "unsupported database name: $db"
done
for cmd in psql pg_dump pg_dumpall pg_restore zstd sha256sum stat find sort tail flock mktemp python3 tar cmp sync hostname date cut cat mv rm mkdir df awk sleep readlink dirname; do
  command -v "$cmd" >/dev/null || die "missing command: $cmd"
done

## destination syntax is deliberately restricted so remote shell arguments stay literal.
remote_host=''
remote_path=$secondary
if [[ -n $secondary ]]; then
  command -v rsync >/dev/null || die 'missing command: rsync'
  if [[ $secondary == *:* ]]; then
    remote_host=${secondary%%:*}; remote_path=${secondary#*:}
    [[ $remote_host =~ ^[A-Za-z0-9_][A-Za-z0-9_.@-]*$ ]] || die 'invalid secondary host'
    command -v ssh >/dev/null || die 'missing command: ssh'
  fi
  [[ $remote_path =~ ^/[A-Za-z0-9_./-]+$ && $remote_path != / ]] || die 'secondary path must be an absolute simple path'
fi
secondary_exists() {
  if [[ -n $remote_host ]]; then
    ssh -o BatchMode=yes -o ConnectTimeout=10 "$remote_host" "test -d '$remote_path'"
  else
    [[ -d $remote_path ]]
  fi
}

## identify the actual cluster, including restarts and timeline changes.
cluster=$(psqlq -d postgres -c "select system_identifier || E'\t' || timeline_id || E'\t' || pg_postmaster_start_time() from pg_control_system(), pg_control_checkpoint()")
server=${PGHOST:-}
if [[ -z $server || $server == /* || $server == localhost || $server == 127.* || $server == ::1 ]]; then
  server=$(hostname -s)
fi
[[ $server =~ ^[A-Za-z0-9_.-]+$ ]] || die 'unsupported host name'
port=$(psqlq -d postgres -c "show port")
[[ $port =~ ^[0-9]+$ ]] || die 'invalid server port'
mkdir -p -- "$backup_root"
backup_root=$(cd -- "$backup_root" && pwd -P)
exec 9>"$backup_root/.pg-backup.lock"
flock -n 9 || die 'another backup process holds the backup-root lock'

if ((check_only)); then
  for db in "${dbs[@]}"; do psqlq -d "$db" -c 'select current_database()'; done
  [[ -z $secondary ]] || secondary_exists || die 'secondary destination unavailable'
  printf 'preflight ok: server=%s port=%s root=%s secondary=%s\n' "$server" "$port" "$backup_root" "$secondary"
  exit 0
fi

free_space_ok() {
  local free
  free=$(df -PB1 -- "$backup_root" | awk 'NR==2 {print $4}')
  [[ $free =~ ^[0-9]+$ && $free -ge $min_free ]]
}

run_dump() {
  free_space_ok || die "less than $min_free_gb GiB free on the backup filesystem"
  pg_dump -w -d "$db" -Fc -Z "$compress" -f "$work/$base.dump" &
  dump_pid=$!
  ## a live guard protects the filesystem even when compression/size cannot be predicted.
  while kill -0 "$dump_pid" 2>/dev/null; do
    if ! free_space_ok; then
      kill "$dump_pid" 2>/dev/null || true
      wait "$dump_pid" 2>/dev/null || true
      dump_pid=''
      die "backup stopped: free space fell below $min_free_gb GiB"
    fi
    sleep 1
  done
  if wait "$dump_pid"; then dump_pid=''; else dump_pid=''; die "pg_dump failed for $db"; fi
}

replicate() {
  local setdir=$1 base receipt staging differences
  [[ -n $secondary ]] || return 0
  base=${setdir##*/}
  receipt=$backup_root/.replicated-$base
  ## failed copies are retried on the next run, even when the database is unchanged.
  if [[ -f $receipt && $(cat "$receipt") == "$secondary" ]]; then
    if [[ -n $remote_host ]]; then
      ssh -o BatchMode=yes -o ConnectTimeout=10 "$remote_host" "test -f '$remote_path/$base/COMPLETE'" && return 0
    elif [[ -f $remote_path/$base/COMPLETE ]]; then
      return 0
    fi
  fi
  secondary_exists || return 1
  ## a crash after remote publication but before the receipt must be recoverable without overwriting a set.
  if [[ -n $remote_host ]]; then
    if ssh -o BatchMode=yes -o ConnectTimeout=10 "$remote_host" "test -e '$remote_path/$base'"; then
      differences=$(rsync -rcn --delete --out-format='%i' -e 'ssh -o BatchMode=yes -o ConnectTimeout=10' \
        "$setdir/" "$remote_host:$remote_path/$base/") || return 1
      [[ -z $differences ]] || { warn "existing secondary set differs: $base"; return 1; }
      printf '%s\n' "$secondary" > "$receipt.tmp"
      mv -f -- "$receipt.tmp" "$receipt"
      return 0
    fi
  elif [[ -e $remote_path/$base ]]; then
    differences=$(rsync -rcn --delete --out-format='%i' "$setdir/" "$remote_path/$base/") || return 1
    [[ -z $differences ]] || { warn "existing secondary set differs: $base"; return 1; }
    printf '%s\n' "$secondary" > "$receipt.tmp"
    mv -f -- "$receipt.tmp" "$receipt"
    return 0
  fi
  staging=$remote_path/.partial-$base
  if [[ -n $remote_host ]]; then
    rsync -a --no-group --delete --partial -e 'ssh -o BatchMode=yes -o ConnectTimeout=10' \
      "$setdir/" "$remote_host:$staging/" || return 1
    ssh -o BatchMode=yes -o ConnectTimeout=10 "$remote_host" \
      "test ! -e '$remote_path/$base' && sync -f '$staging' && mv -T '$staging' '$remote_path/$base' && sync -f '$remote_path'" || return 1
  else
    [[ ! -e $remote_path/$base ]] || return 1
    rsync -a --no-group --delete --partial "$setdir/" "$staging/" || return 1
    sync -f "$staging" || return 1
    mv -T -- "$staging" "$remote_path/$base" || return 1
    sync -f "$remote_path" || return 1
  fi
  printf '%s\n' "$secondary" > "$receipt.tmp"
  mv -f -- "$receipt.tmp" "$receipt"
}

list_sets() {
  local prefix=$1 dir suffix
  while IFS= read -r dir; do
    [[ ${dir##*/} == "$prefix"* ]] || continue
    suffix=${dir##*/}; suffix=${suffix#"$prefix"}
    [[ $suffix =~ ^[0-9]{8}_[0-9]{6}$ ]] && printf '%s\n' "$dir"
  done < <(find "$backup_root" -mindepth 2 -maxdepth 2 -name COMPLETE -type f -printf '%h\n')
  return 0
}

prune() {
  local prefix=$1 old
  ## only published, complete sets count toward retention; partial directories never do.
  mapfile -t old < <(list_sets "$prefix" | sort -r | tail -n +"$((keep + 1))")
  for old in "${old[@]}"; do
    rm -rf -- "$old"
    rm -f -- "$backup_root/.replicated-${old##*/}"
  done
}

for db in "${dbs[@]}"; do
  work=$(mktemp -d "$backup_root/.partial-XXXXXXXX")
  prefix=pgbackup_${server}_p${port}_${db}_
  dbinfo=$(psqlq -d "$db" -c "select oid from pg_database where datname=current_database()")
  ## record the WAL boundary before any snapshot; concurrent commits remain dirty next time.
  lsn=$(psqlq -d postgres -c 'select pg_current_wal_insert_lsn()')
  printf '%s\t%s\n' "$cluster" "$dbinfo" > "$work/identity"
  printf '%s\n' "$lsn" > "$work/lsn"
  psqlq -d "$db" > "$work/sequences" <<'SQL'
select format('select %L, last_value, is_called from %I.%I;', n.nspname || '.' || c.relname, n.nspname, c.relname)
from pg_class c join pg_namespace n on n.oid=c.relnamespace where c.relkind='S' order by n.nspname,c.relname
\gexec
SQL
  psqlq -d postgres -c "select name,setting from pg_settings order by name" > "$work/settings"
  ## read server-side paths, including includes and ALTER SYSTEM, even over a remote connection.
  psqlq -d postgres > "$work/config" <<'SQL'
with paths as (
  select setting as path from pg_settings where name in ('config_file','hba_file','ident_file')
  union select sourcefile from pg_file_settings
  union select file_name from pg_hba_file_rules union select file_name from pg_ident_file_mappings
  union select current_setting('data_directory') || '/postgresql.auto.conf'
)
select path || E'\t' || coalesce(replace(encode(pg_read_binary_file(path,true),'base64'),E'\n',''),'MISSING')
from paths where path is not null order by path;
SQL
  unlogged=$(psqlq -d "$db" -c "select count(*) from pg_class where relpersistence='u' and relkind in ('r','m','S')")
  latest=$(list_sets "$prefix" | sort | tail -n 1)
  unchanged=0
  if ((if_changed)) && [[ -n $latest && $unlogged == 0 ]]; then
    previous=$latest/${latest##*/}.state.tar.zst
    if [[ -f $previous ]]; then
      mkdir "$work/previous"
      ## only our five fixed state members are read; no archive paths are extracted.
      state_ok=1
      for member in identity lsn sequences settings config; do
        if ! zstd -dc "$previous" | tar -xOf - "$member" > "$work/previous/$member"; then state_ok=0; break; fi
      done
      for member in identity sequences settings config; do
        cmp -s "$work/$member" "$work/previous/$member" || state_ok=0
      done
      if ((state_ok)); then
        old_lsn=$(cat "$work/previous/lsn")
        if [[ $old_lsn == "$lsn" ]]; then
          unchanged=1
        elif [[ $old_lsn =~ ^[0-9A-F]+/[0-9A-F]+$ ]] && [[ $(psqlq -d postgres -c "select pg_current_wal_flush_lsn() >= '$lsn'::pg_lsn") == t ]]; then
          ## ignore only known bookkeeping records; unknown records and missing WAL require a backup.
          if dirty=$(psqlq -d postgres -c "set statement_timeout='30s'; select exists(select 1 from public.pg_get_wal_records_info('$old_lsn','$lsn') where not ((resource_manager='XLOG' and record_type in ('CHECKPOINT_ONLINE','CHECKPOINT_SHUTDOWN','CHECKPOINT_REDO','FPI_FOR_HINT')) or (resource_manager='Standby' and record_type='RUNNING_XACTS') or (resource_manager='Heap2' and record_type='PRUNE_ON_ACCESS')))" 2>"$work/wal-error"); then
            [[ $dirty == f ]] && unchanged=1
          else
            warn "WAL history unavailable for $db; taking a backup"
          fi
        fi
      fi
    fi
  fi
  if ((unchanged)); then
    verifier=$(dirname -- "$(readlink -f -- "$0")")/pg-restore.sh
    if ! "$verifier" --verify-only "$latest" backup_integrity_check; then
      warn "previous backup failed integrity verification; taking a new backup"
      unchanged=0
    fi
  fi
  if ((unchanged)); then
    printf 'unchanged: %s (last backup %s)\n' "$db" "${latest##*/}"
    if replicate "$latest"; then prune "$prefix"; else warn "secondary copy failed: $latest"; failed=1; fi
    rm -rf -- "$work"; work=''
    continue
  fi

  base=$prefix$(date +%Y%m%d_%H%M%S)
  setdir=$backup_root/$base
  [[ ! -e $setdir ]] || die "backup name collision: $setdir"
  printf 'creating %s\n' "$setdir"
  run_dump
  pg_restore -l "$work/$base.dump" >/dev/null
  pg_dumpall -w -g | zstd -q -9 -c > "$work/$base.globals.sql.zst"
  pg_dumpall -w --roles-only | zstd -q -9 -c > "$work/$base.roles.sql.zst"
  psqlq -d "$db" -c "select extname,extversion from pg_extension order by extname" > "$work/$base.extensions.tsv"
  psqlq -d "$db" > "$work/$base.tablespaces.tsv" <<'SQL'
select spcname from pg_tablespace where oid in (
  select reltablespace from pg_class where not relisshared and reltablespace<>0
  union select dattablespace from pg_database where datname=current_database()
) order by spcname;
SQL
  ## psql's quoted variable allows the database metadata to follow a new restore name.
  psqlq -d "$db" > "$work/$base.database.sql" <<'SQL'
select format('CREATE DATABASE :"restore_db" TEMPLATE template0 OWNER %I ENCODING %L LC_COLLATE %L LC_CTYPE %L LOCALE_PROVIDER %s%s%s TABLESPACE %I;',
  pg_get_userbyid(datdba),pg_encoding_to_char(encoding),datcollate,datctype,
  case datlocprovider when 'c' then 'libc' when 'i' then 'icu' when 'b' then 'builtin' end,
  case when datlocprovider='i' then format(' ICU_LOCALE %L',to_jsonb(d)->>'datlocale')
    when datlocprovider='b' then format(' BUILTIN_LOCALE %L',to_jsonb(d)->>'datlocale') else '' end,
  case when daticurules is not null then format(' ICU_RULES %L',daticurules) else '' end,
  (select spcname from pg_tablespace where oid=dattablespace)) from pg_database d where datname=current_database();
select format('COMMENT ON DATABASE :"restore_db" IS %L;',shobj_description(oid,'pg_database'))
from pg_database where datname=current_database();
select format('ALTER DATABASE :"restore_db" CONNECTION LIMIT %s;',datconnlimit)
from pg_database where datname=current_database();
select 'REVOKE ALL ON DATABASE :"restore_db" FROM PUBLIC;' from pg_database
where datname=current_database() and datacl is not null;
select format('GRANT %s ON DATABASE :"restore_db" TO %s%s;',a.privilege_type,
  case when a.grantee=0 then 'PUBLIC' else quote_ident(pg_get_userbyid(a.grantee)) end,
  case when a.is_grantable then ' WITH GRANT OPTION' else '' end)
from pg_database d cross join lateral aclexplode(d.datacl) a where d.datname=current_database();
select case when setrole=0 then format('ALTER DATABASE :"restore_db" SET %I TO %L;',split_part(v,'=',1),substr(v,strpos(v,'=')+1))
  else format('ALTER ROLE %I IN DATABASE :"restore_db" SET %I TO %L;',pg_get_userbyid(setrole),split_part(v,'=',1),substr(v,strpos(v,'=')+1)) end
from pg_db_role_setting cross join lateral unnest(setconfig) v
where setdatabase=(select oid from pg_database where datname=current_database());
SQL
  tar -C "$work" -cf - identity lsn sequences settings config | zstd -q -9 -c > "$work/$base.state.tar.zst"
  zstd -q -9 -c "$work/config" > "$work/$base.config.tsv.zst"
  zstd -q -9 -c "$work/settings" > "$work/$base.settings.tsv.zst"
  rm -rf -- "$work/previous"
  rm -f -- "$work/identity" "$work/lsn" "$work/sequences" "$work/settings" "$work/config" "$work/wal-error"
  {
    printf 'kind\tpath\tsize\tsha256\n'
    printf 'meta\tcreated_at\t\t%s\n' "$(date -Is)"
    printf 'meta\tdatabase\t\t%s\n' "$db"
    while IFS= read -r file; do
      [[ $file == "$base.manifest.tsv" ]] && continue
      printf 'file\t%s\t%s\t%s\n' "$file" "$(stat -c %s "$work/$file")" "$(sha256sum "$work/$file" | cut -d ' ' -f1)"
    done < <(find "$work" -maxdepth 1 -type f -printf '%f\n' | sort)
  } > "$work/$base.manifest.tsv"
  printf 'complete\t%s\n' "$(date -Is)" > "$work/COMPLETE"
  ## publish only after payloads and manifest have reached stable storage.
  sync -f "$work"
  mv -T -- "$work" "$setdir"; work=''
  sync -f "$backup_root"
  if replicate "$setdir"; then prune "$prefix"; else warn "secondary copy failed; local backup retained: $setdir"; failed=1; fi
  printf 'done %s\n' "$setdir"
done
exit "$failed"
