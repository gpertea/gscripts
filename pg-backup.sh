#!/usr/bin/env bash
set -Eeuo pipefail
umask 077

## defaults are intentionally small and local to this site.
backup_root=${PG_BACKUP_ROOT:-/data/backups/postgres}
secondary_dests=( "linwks34:/nfs/gdata/postgres" )
default_db=${PG_BACKUP_DB:-rse}
keep_local=${PG_BACKUP_KEEP:-20}
compress=${PG_BACKUP_COMPRESS:-zstd:9}
check_only=0
dbs=()

die() { printf 'error: %s\n' "$*" >&2; exit 1; }
warn() { printf 'warning: %s\n' "$*" >&2; }
need() { command -v "$1" >/dev/null 2>&1 || die "missing command: $1"; }
safe_name() { printf '%s' "$1" | sed -E 's/[^A-Za-z0-9_.-]+/_/g; s/^_+//; s/_+$//'; }
sql_lit() { local s; s=$(printf '%s' "$1" | sed "s/'/''/g"); printf "'%s'" "$s"; }

usage() {
  printf 'usage: %s [--check] [db ...]\n' "${0##*/}"
  printf 'env: PG_BACKUP_USER, PGHOST, PGPORT, PG_BACKUP_ROOT, PG_BACKUP_DB, PG_BACKUP_KEEP, PG_BACKUP_COMPRESS\n'
  exit 0;
}

if (( $# == 0 )); then usage; fi

while (($#)); do
  case $1 in
    --check) check_only=1 ;;
    -h|--help) usage ;;
    --) shift; dbs+=( "$@" ); break ;;
    -*) die "unknown option: $1" ;;
    *) dbs+=( "$1" ) ;;
  esac
  shift
done


((${#dbs[@]})) || dbs=( "$default_db" )
export PGUSER=${PG_BACKUP_USER:-gpertea}

for cmd in psql pg_dump pg_dumpall pg_restore zstd rsync ssh sha256sum stat sed awk find hostname date; do need "$cmd"; done

## use the real host name for local connections, not localhost.
conn_host=${PGHOST:-}
if [[ -z $conn_host || $conn_host == localhost || $conn_host == 127.* || $conn_host == ::1 || $conn_host == /* ]]; then
  server=$(hostname -s)
else
  server=$conn_host
fi
server_safe=$(safe_name "$server")
port=$(psql -d postgres -Atq -v ON_ERROR_STOP=1 -c "select current_setting('port')" 2>/dev/null || printf '%s' "${PGPORT:-5432}")
port_safe=$(safe_name "$port")

check_dest() {
  local dest=$1 host path
  ## remote destinations must already exist; this avoids silent bad copies.
  if [[ $dest == *:* ]]; then
    host=${dest%%:*}; path=${dest#*:}
    ssh -o BatchMode=yes -o ConnectTimeout=5 "$host" "test -d '$path'" || return 1
  else
    [[ -d $dest ]] || return 1
  fi
}

copy_configs() {
  local cfgdir=$1 f b
  mkdir -p "$cfgdir"
  while IFS= read -r f; do
    [[ -n $f && -r $f ]] || { [[ -n $f ]] && warn "cannot read config: $f"; continue; }
    b=$(basename "$f")
    zstd -q -9 -c "$f" > "$cfgdir/$b.zst"
  done < <(psql -d postgres -Atq -v ON_ERROR_STOP=1 -c \
    "select current_setting('config_file') union all select current_setting('hba_file') union all select current_setting('ident_file')" 2>/dev/null)
  find "$cfgdir" -type f -print -quit | grep -q . || rmdir "$cfgdir"
}

write_manifest() {
  local setdir=$1 base=$2 db=$3 rel size sum
  local manifest=$setdir/$base.manifest.tsv
  {
    printf 'kind\tpath\tsize\tsha256\n'
    printf 'meta\tserver\t\t%s\n' "$server"
    printf 'meta\tport\t\t%s\n' "$port"
    printf 'meta\tdatabase\t\t%s\n' "$db"
    printf 'meta\tcreated_at\t\t%s\n' "$(date -Is)"
    printf 'meta\tpg_dump\t\t%s\n' "$(pg_dump --version)"
    printf 'meta\tpg_restore\t\t%s\n' "$(pg_restore --version)"
    printf 'meta\tcompression\t\t%s\n' "$compress"
  } > "$manifest"
  ## record every payload file so restore can catch partial or damaged sets.
  while IFS= read -r rel; do
    [[ ${rel##*/} == "$base.manifest.tsv" || ${rel##*/} == COMPLETE ]] && continue
    size=$(stat -c '%s' "$setdir/$rel")
    sum=$(sha256sum "$setdir/$rel" | awk '{print $1}')
    printf 'file\t%s\t%s\t%s\n' "$rel" "$size" "$sum" >> "$manifest"
  done < <(cd "$setdir" && find . -type f -printf '%P\n' | sort)
}

prune_local() {
  local db_safe=$1 pattern old
  pattern="pgbackup_${server_safe}_p${port_safe}_${db_safe}_*"
  ## retention is local only; remote retention can have its own policy.
  mapfile -t old < <(find "$backup_root" -maxdepth 1 -type d -name "$pattern" | sort -r | tail -n +"$((keep_local + 1))")
  ((${#old[@]})) && rm -rf -- "${old[@]}"
}

db_exists() {
  local db=$1
  psql -d postgres -Atq -v ON_ERROR_STOP=1 -c "select 1 from pg_database where datname = $(sql_lit "$db") and datallowconn" | grep -qx 1
}

mkdir -p "$backup_root" || die "cannot create $backup_root"
[[ -w $backup_root ]] || die "backup root is not writable: $backup_root"
check_dest "$backup_root" || die "bad local backup root: $backup_root"
for dest in "${secondary_dests[@]}"; do check_dest "$dest" || die "secondary destination unavailable: $dest"; done

if ((check_only)); then
  printf 'server=%s port=%s user=%s root=%s compression=%s\n' "$server" "$port" "$PGUSER" "$backup_root" "$compress"
  printf 'databases=%s\n' "${dbs[*]}"
  printf 'secondary=%s\n' "${secondary_dests[*]}"
  for db in "${dbs[@]}"; do db_exists "$db" || die "database not found or not connectable: $db"; done
  exit 0
fi

for db in "${dbs[@]}"; do
  db_exists "$db" || die "database not found or not connectable: $db"
  ts=$(date +%Y%m%d_%H%M%S)
  db_safe=$(safe_name "$db")
  base="pgbackup_${server_safe}_p${port_safe}_${db_safe}_${ts}"
  setdir="$backup_root/$base"
  mkdir "$setdir" || die "cannot create backup set: $setdir"

  printf 'creating %s\n' "$setdir"
  pg_dump -d "$db" -Fc -Z "$compress" -f "$setdir/$base.dump"
  pg_dumpall -g | zstd -q -9 -c > "$setdir/$base.globals.sql.zst"
  { psql -d postgres -Atq -v ON_ERROR_STOP=1 -c "select name || E'\t' || setting from pg_settings order by name"; psql -d postgres -Atq -v ON_ERROR_STOP=1 -c "show all"; } | zstd -q -9 -c > "$setdir/$base.settings.tsv.zst"
  copy_configs "$setdir/config"
  write_manifest "$setdir" "$base" "$db"
  printf 'complete\t%s\n' "$(date -Is)" > "$setdir/COMPLETE"

  for dest in "${secondary_dests[@]}"; do
    printf 'syncing %s to %s\n' "$base" "$dest"
    rsync -a --partial "$setdir/" "$dest/$base/" || die "rsync failed for $dest"
  done
  prune_local "$db_safe"
  printf 'done %s\n' "$setdir"
done
