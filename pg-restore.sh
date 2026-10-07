#!/usr/bin/env bash
set -Eeuo pipefail
umask 077
replace=0
restore_globals=0
precheck_only=0
verify_only=0
args=()
work=''
stage=''
stage_oid=''
blocked_oid=''
old_allow=''
export PGUSER=${PG_BACKUP_USER:-gpertea} PGCONNECT_TIMEOUT=${PGCONNECT_TIMEOUT:-10}

die() { printf 'error: %s\n' "$*" >&2; exit 1; }
psqlq() { psql -X -w -Atq -v ON_ERROR_STOP=1 "$@"; }
sql_lit() { local s=${1//\\/\\\\}; printf "E'%s'" "${s//\'/\'\'}"; }
usage() {
  printf 'usage: %s [--replace] [--globals] [--precheck-only|--verify-only] <backup-set-dir|dump-file> <target-db> [jobs]\n' "${0##*/}"
  printf '%s\n' '--replace validates a staging restore first and retains the original under a rollback name.' \
    '--globals applies role definitions/memberships cluster-wide, transactionally; tablespaces must already exist.' \
    'Ownership and grants are preserved. Configuration files require manual review.'
}
cleanup() {
  local status=$?
  trap - EXIT
  ## use OIDs so a failed rename never drops the requested database or somebody else's replacement.
  if [[ -n $blocked_oid ]]; then
    psqlq -d postgres <<SQL || printf 'warning: could not restore connection access for original database OID %s\n' "$blocked_oid" >&2
select format('ALTER DATABASE %I ALLOW_CONNECTIONS %s;',datname,'$old_allow') from pg_database
where oid=$blocked_oid and datname=$(sql_lit "$target_db")
\gexec
SQL
  fi
  if [[ -n $stage_oid ]]; then
    if [[ $(psqlq -d postgres -c "select oid from pg_database where datname=$(sql_lit "$stage")" 2>/dev/null) == "$stage_oid" ]]; then
      dropdb -w --force -- "$stage" || printf 'warning: failed staging database retained: %s\n' "$stage" >&2
    fi
  fi
  [[ -z $work ]] || rm -rf -- "$work"
  exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
while (($#)); do
  case $1 in
    --replace) replace=1 ;;
    --globals) restore_globals=1 ;;
    --precheck-only) precheck_only=1 ;;
    --verify-only) verify_only=1 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; args+=( "$@" ); break ;;
    -*) die "unknown option: $1" ;;
    *) args+=( "$1" ) ;;
  esac
  shift
done
((${#args[@]} >= 2 && ${#args[@]} <= 3)) || { usage; exit 1; }
src=${args[0]}; target_db=${args[1]}; jobs=${args[2]:-4}
[[ $jobs =~ ^[1-9][0-9]{0,2}$ ]] || die 'jobs must be between 1 and 999'
[[ $target_db =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]{0,62}$ ]] || die 'unsupported target database name'
[[ $target_db != postgres && $target_db != template* ]] || die 'refusing a maintenance/template target database'
for cmd in psql pg_restore zstd python3 mktemp createdb dropdb head tail cp awk date rm; do
  command -v "$cmd" >/dev/null || die "missing command: $cmd"
done
if [[ -d $src ]]; then
  setdir=$(cd -- "$src" && pwd -P)
  base=${setdir##*/}
elif [[ -f $src ]]; then
  setdir=$(cd -- "$(dirname -- "$src")" && pwd -P)
  base=${src##*/}; base=${base%.dump}
else
  die "backup not found: $src"
fi
[[ $base =~ ^pgbackup_[A-Za-z0-9_.-]+_p[0-9]+_[A-Za-z0-9_.-]+_[0-9]{8}_[0-9]{6}$ ]] || die 'unrecognized backup name'
dump=$setdir/$base.dump

## reject malformed/incomplete manifests, symlinks, traversal and unrecorded payloads.
python3 - "$setdir" "$base" <<'PY'
import csv, hashlib, pathlib, re, sys
root=pathlib.Path(sys.argv[1]); base=sys.argv[2]
def fail(message):
    sys.exit('error: '+message)
for name in ('COMPLETE',base+'.manifest.tsv'):
    p=root/name
    if not p.is_file() or p.is_symlink(): fail('missing or unsafe '+name)
if not (root/'COMPLETE').read_text().startswith('complete\t'): fail('invalid COMPLETE marker')
seen=set()
with (root/(base+'.manifest.tsv')).open(newline='') as f:
    rows=csv.reader(f,delimiter='\t')
    if next(rows,None)!=['kind','path','size','sha256']: fail('invalid manifest header')
    for row in rows:
        if len(row)!=4: fail('malformed manifest row')
        kind,name,size,digest=row
        if kind in ('meta','extension'): continue
        if kind!='file': fail('unknown manifest record')
        p=pathlib.PurePosixPath(name)
        if p.is_absolute() or '..' in p.parts or str(p)!=name or name in seen: fail('unsafe/duplicate manifest path')
        if not re.fullmatch(r'[0-9]+',size) or not re.fullmatch(r'[0-9a-f]{64}',digest): fail('invalid size/hash')
        full=root/p
        if any((root.joinpath(*p.parts[:i])).is_symlink() for i in range(1,len(p.parts)+1)):
            fail('symlink in manifest path')
        if not full.is_file() or full.stat().st_size!=int(size): fail('missing file or size mismatch: '+name)
        h=hashlib.sha256()
        with full.open('rb') as data:
            for block in iter(lambda:data.read(8*1024*1024),b''): h.update(block)
        if h.hexdigest()!=digest: fail('hash mismatch: '+name)
        seen.add(name)
if base+'.dump' not in seen: fail('dump absent from manifest')
actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() or p.is_symlink()}
if actual!=seen|{'COMPLETE',base+'.manifest.tsv'}: fail('unrecorded or missing payload')
PY
work=$(mktemp -d)
pg_restore -l "$dump" > "$work/toc"
((verify_only)) && { printf 'verified %s\n' "$setdir"; exit 0; }

## a complete data-stream read catches archive errors that the TOC alone cannot detect.
pg_restore --file=/dev/null "$dump"
if [[ -f $setdir/$base.extensions.tsv ]]; then
  cp -- "$setdir/$base.extensions.tsv" "$work/extensions"
else
  awk '$0 ~ / EXTENSION - / && $0 !~ /COMMENT - EXTENSION/ {print $NF "\t"}' "$work/toc" > "$work/extensions"
fi
while IFS=$'\t' read -r ext version; do
  [[ -n $ext ]] || continue
  available=$(psqlq -d postgres -c "select default_version from pg_available_extensions where name=$(sql_lit "$ext")")
  [[ -n $available ]] || die "missing target extension: $ext"
  ## pg_restore CREATE EXTENSION uses the target default, so availability alone is insufficient.
  [[ -z $version || $version == "$available" ]] || die "extension default version differs: $ext (backup $version, target $available)"
done < "$work/extensions"
if [[ -f $setdir/$base.tablespaces.tsv ]]; then
  while IFS= read -r tablespace; do
    [[ $(psqlq -d postgres -c "select 1 from pg_tablespace where spcname=$(sql_lit "$tablespace")") == 1 ]] \
      || die "missing target tablespace: $tablespace (create it manually before restore)"
  done < "$setdir/$base.tablespaces.tsv"
fi

exists=$(psqlq -d postgres -c "select oid from pg_database where datname=$(sql_lit "$target_db")")
[[ -z $exists || $replace == 1 ]] || die 'target database exists; use --replace to retain and replace it'
if [[ -n $exists ]]; then
  old_allow=$(psqlq -d postgres -c "select case when datallowconn then 'true' else 'false' end from pg_database where oid=$exists and not datistemplate")
  [[ -n $old_allow ]] || die 'cannot replace a template database'
fi
if ((restore_globals)); then
  roles=$setdir/$base.roles.sql.zst
  [[ -f $roles ]] || die 'legacy backup lacks a roles-only file; review/apply its globals SQL manually'
  zstd -dc "$roles" > "$work/roles.sql"
  ## existing bootstrap roles must not make a globals restore abort; other errors still roll back.
  python3 - "$work/roles.sql" > "$work/roles-merge.sql" <<'PY'
import pathlib,re,sys,uuid
for line in pathlib.Path(sys.argv[1]).read_text().splitlines():
    m=re.fullmatch(r'CREATE ROLE (.+);',line)
    if m:
        ident=m[1]
        name=ident[1:-1].replace('""','"') if ident.startswith('"') else ident
        literal="E'"+name.replace('\\','\\\\').replace("'","''")+"'"
        tag='$role_'+uuid.uuid4().hex+'$'
        line=f'DO {tag} BEGIN IF NOT EXISTS (SELECT FROM pg_roles WHERE rolname={literal}) THEN CREATE ROLE {ident}; END IF; END {tag};'
    print(line)
PY
fi
((precheck_only)) && { printf 'precheck ok: %s -> %s (no restore performed)\n' "$base" "$target_db"; exit 0; }

## --globals is explicit: it may change existing role attributes across the cluster.
if ((restore_globals)); then psqlq -d postgres --single-transaction -f "$work/roles-merge.sql"; fi
stage=pgrestore_$(date +%Y%m%d_%H%M%S)_${RANDOM}_$$
rollback=pgrollback_$(date +%Y%m%d_%H%M%S)_${RANDOM}_$$
[[ -z $(psqlq -d postgres -c "select oid from pg_database where datname=$(sql_lit "$stage")") ]] || die 'staging name collision'
if [[ -f $setdir/$base.database.sql ]]; then
  ## create separately so cleanup knows the new OID even if later metadata restoration fails.
  head -n 1 "$setdir/$base.database.sql" > "$work/create.sql"
  psqlq -d postgres -v restore_db="$stage" -f "$work/create.sql"
else
  printf 'warning: legacy backup lacks database-level metadata; using template0 defaults\n' >&2
  createdb -w -T template0 -- "$stage"
fi
stage_oid=$(psqlq -d postgres -c "select oid from pg_database where datname=$(sql_lit "$stage")")
[[ $stage_oid =~ ^[0-9]+$ ]] || die 'could not identify staging database'
pg_restore -w --exit-on-error -j "$jobs" -d "$stage" "$dump"
if [[ -f $setdir/$base.database.sql ]]; then
  tail -n +2 "$setdir/$base.database.sql" > "$work/properties.sql"
  psqlq -d postgres --single-transaction -v restore_db="$stage" -f "$work/properties.sql"
fi

## the original stays online throughout the expensive restore. Only the final switch interrupts clients.
if [[ -n $exists ]]; then
  blocked_oid=$exists
  psqlq -d postgres <<SQL
select pg_advisory_lock(hashtext('pg-restore'),hashtext($(sql_lit "$target_db")));
DO \$guard\$ BEGIN
  IF (select oid from pg_database where datname=$(sql_lit "$target_db")) IS DISTINCT FROM $exists::oid THEN
    RAISE EXCEPTION 'target changed during restore';
  END IF;
END \$guard\$;
ALTER DATABASE "$target_db" ALLOW_CONNECTIONS false;
select pg_terminate_backend(pid,5000) from pg_stat_activity where datid=$exists;
BEGIN;
ALTER DATABASE "$target_db" RENAME TO "$rollback";
ALTER DATABASE "$stage" RENAME TO "$target_db";
COMMIT;
SQL
  blocked_oid=''
  printf 'original retained with connections disabled: %s\n' "$rollback"
else
  psqlq -d postgres -c "ALTER DATABASE \"$stage\" RENAME TO \"$target_db\""
fi
stage_oid=''
printf 'restored %s into %s at %s\n' "$base" "$target_db" "$(date -Is)"
