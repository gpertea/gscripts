#!/usr/bin/env python3
"""Classify PostgreSQL 17 WAL for one database; unknown evidence requires a backup."""
import argparse
import json
import re
import subprocess
import sys

# these records maintain physical storage or transaction bookkeeping, not logical contents.
MAINTENANCE = {
    'XLOG': {'CHECKPOINT_ONLINE', 'CHECKPOINT_SHUTDOWN', 'CHECKPOINT_REDO',
             'FPI_FOR_HINT', 'NOOP', 'SWITCH', 'NEXTOID'},
    'Heap2': {'PRUNE_ON_ACCESS', 'PRUNE_VACUUM_SCAN', 'PRUNE_VACUUM_CLEANUP', 'VISIBLE'},
    'Standby': {'RUNNING_XACTS', 'LOCK', 'INVALIDATIONS'},
    'Btree': {'REUSE_PAGE'},
}
PHYSICAL = {'Heap', 'Heap2', 'Btree', 'Hash', 'Gin', 'Gist', 'SPGist', 'BRIN', 'Sequence', 'Generic', 'XLOG'}
BOOKKEEPING = {'CLOG', 'MultiXact', 'CommitTs'}


def path_database(path):
    """Decode PostgreSQL's documented relation-path forms, rejecting unfamiliar ones."""
    if re.fullmatch(r'global/\d+(?:_(?:fsm|vm|init))?', path):
        return 0
    match = re.fullmatch(r'(?:base|pg_tblspc/\d+/PG_[^/]+)/([0-9]+)/\d+(?:_(?:fsm|vm|init))?', path)
    return int(match[1]) if match else None


def classify(records, database, xid_floor=None):
    parents = {}
    locations = {}
    transactions = []
    unknown = []
    changed = []
    dependencies = []

    def root(xid):
        parents.setdefault(xid, xid)
        while parents[xid] != xid:
            parents[xid] = parents[parents[xid]]
            xid = parents[xid]
        return xid

    def join(a, b):
        parents[root(b)] = root(a)

    for record in records:
        manager, kind = record['manager'], record['kind']
        xid, description = int(record['xid']), record['description'] or ''
        databases = set(record['databases'])
        dependency_only = record.get('dependency_only', False)
        label = f'{manager}/{kind} at {record["lsn"]}'
        if manager == 'Standby' and kind == 'RUNNING_XACTS' and xid_floor is None:
            match = re.search(r'\bnextXid ([0-9]+)\b', description)
            if match:
                xid_floor = int(match[1])
        if manager == 'Transaction':
            # prepared commits carry the original XID in their description, not necessarily the header.
            if kind in ('COMMIT_PREPARED', 'ABORT_PREPARED'):
                match = re.match(r'([0-9]+): ', description)
                if not match:
                    unknown.append(label)
                    continue
                xid = int(match[1])
            elif kind == 'ASSIGNMENT':
                match = re.match(r'xtop ([0-9]+): ', description)
                if not match:
                    unknown.append(label)
                    continue
                xid = int(match[1])
            subxacts = re.search(r'(?:^|[;:]) subxacts:([0-9 ]+)(?:;|$)', description)
            if 'subxacts:' in description and not subxacts:
                unknown.append(label)
            if xid and subxacts:
                for child in subxacts[1].split():
                    join(xid, int(child))
            if kind in ('COMMIT', 'COMMIT_PREPARED', 'PREPARE'):
                # an XID allocated after the saved/in-range boundary cannot have earlier writes.
                known_new = xid_floor is not None and xid >= 3 and ((xid - xid_floor) & 0xffffffff) < 0x80000000
                transactions.append((xid, label, known_new))
            elif kind not in ('ABORT', 'ABORT_PREPARED', 'ASSIGNMENT', 'INVALIDATION'):
                unknown.append(label)
            continue

        if manager in BOOKKEEPING or kind in MAINTENANCE.get(manager, set()):
            continue
        if manager == 'Storage' and kind in ('CREATE', 'TRUNCATE'):
            # storage creation/truncation records have no block reference.
            match = re.fullmatch(r'(\S+)(?: to [0-9]+ blocks flags [0-9]+)?', description)
            db = path_database(match[1]) if match else None
            if db is not None:
                databases.add(db)
        if manager not in PHYSICAL | {'Storage'} or not databases:
            unknown.append(label)
            continue
        if 0 in databases and dependency_only:
            databases.remove(0)
            dependencies.append((xid, label, False))
        if xid:
            locations.setdefault(xid, set()).update(databases)
        if database in databases:
            changed.append(label)
        if 0 in databases:
            unknown.append('shared catalog change: ' + label)

    # shared dependency entries belong to their transaction's database, not every database.
    # a parent commit can refer to writes logged with child XIDs after SAVEPOINT.
    grouped = {}
    for xid, databases in locations.items():
        grouped.setdefault(root(xid), set()).update(databases)
    for xid, label, known_new in transactions + dependencies:
        databases = grouped.get(root(xid), set()) if xid else set()
        if not databases and known_new:
            continue
        if not databases:
            # writes may precede the saved boundary and commit after the dump snapshot.
            unknown.append('transaction cannot be assigned to a database: ' + label)
        elif database in databases:
            changed.append(label)
        elif 0 in databases:
            unknown.append('shared catalog transaction: ' + label)
    if changed:
        return 'changed', changed[0]
    if unknown:
        return 'uncertain', unknown[0]
    return 'unchanged', 'no changes for this database in the inspected WAL'


def inspect(database, start, end, xid_floor=None):
    if not all(re.fullmatch(r'[0-9A-F]+/[0-9A-F]+', value) for value in (start, end)):
        raise ValueError('invalid WAL boundary')
    # use typed database OIDs from pg_walinspect, not regexes over block descriptions.
    sql = f"""
SET statement_timeout='30s';
SELECT CASE WHEN current_setting('server_version_num')::int / 10000 = 17
  AND pg_current_wal_flush_lsn() >= '{end}'::pg_lsn THEN 'ready' ELSE 'unsupported or unflushed' END;
WITH dependency_files AS (
  SELECT pg_relation_filenode('pg_shdepend'::regclass) AS node
  UNION SELECT pg_relation_filenode(indexrelid) FROM pg_index WHERE indrelid='pg_shdepend'::regclass
), blocks AS MATERIALIZED (
  SELECT start_lsn,array_agg(DISTINCT reldatabase::bigint) AS databases,
    bool_and(reldatabase<>0 OR relfilenode IN (SELECT node FROM dependency_files)) AS dependency_only
  FROM public.pg_get_wal_block_info('{start}','{end}',false) GROUP BY start_lsn
)
SELECT json_build_object('lsn',r.start_lsn,'xid',r.xid::text,'manager',r.resource_manager,
  'kind',r.record_type,'description',r.description,'databases',coalesce(b.databases,ARRAY[]::bigint[]),
  'dependency_only',coalesce(b.dependency_only,false))
FROM public.pg_get_wal_records_info('{start}','{end}') r LEFT JOIN blocks b USING(start_lsn) ORDER BY r.start_lsn;
"""
    result = subprocess.run(['psql', '-X', '-w', '-Atq', '-v', 'ON_ERROR_STOP=1', '-d', 'postgres'],
                            input=sql, text=True, capture_output=True, timeout=65, check=True)
    lines = result.stdout.splitlines()
    if not lines or lines[0] != 'ready':
        raise ValueError('PostgreSQL 17 and flushed WAL are required')
    return classify((json.loads(line) for line in lines[1:]), database, xid_floor)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('database_oid', type=int)
    parser.add_argument('start_lsn')
    parser.add_argument('end_lsn')
    parser.add_argument('--xid-floor', type=int)
    args = parser.parse_args()
    try:
        status, reason = inspect(args.database_oid, args.start_lsn, args.end_lsn, args.xid_floor)
    except (ValueError, KeyError, subprocess.SubprocessError) as error:
        detail = ((error.stderr or 'psql failed').strip().splitlines()[-1]
                  if isinstance(error, subprocess.CalledProcessError) else str(error))
        status, reason = 'uncertain', 'WAL inspection failed: ' + detail
    print(f'{status}: {reason}')
    return 0 if status == 'unchanged' else 1


if __name__ == '__main__':
    sys.exit(main())
