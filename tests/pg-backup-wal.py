#!/usr/bin/env python3
"""Regression checks for WAL attribution and incomplete transaction evidence."""
import importlib.util
from pathlib import Path
import sys
import unittest

sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('backup_wal', Path(__file__).resolve().parents[1] / 'pg-backup-wal.py')
wal = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wal)


def record(manager, kind, xid=0, databases=(), description=''):
    return dict(manager=manager, kind=kind, xid=xid, databases=list(databases),
                description=description, lsn='0/1234')


class Attribution(unittest.TestCase):
    def status(self, records, database=10):
        return wal.classify(records, database)[0]

    def test_other_database_commit(self):
        records = [record('Heap', 'INSERT', 100, [20]), record('Transaction', 'COMMIT', 100)]
        self.assertEqual(self.status(records), 'unchanged')
        self.assertEqual(self.status(records, 20), 'changed')

    def test_late_commit_without_writes_in_window(self):
        self.assertEqual(self.status([record('Transaction', 'COMMIT', 100)]), 'uncertain')

    def test_subtransaction_commit(self):
        records = [record('Heap', 'INSERT', 101, [20]),
                   record('Transaction', 'COMMIT', 100, description='2026-09-17; subxacts: 101')]
        self.assertEqual(self.status(records), 'unchanged')
        self.assertEqual(self.status(records, 20), 'changed')

    def test_assignment_and_prepared_commit(self):
        records = [record('Heap', 'INSERT', 101, [20]),
                   record('Transaction', 'ASSIGNMENT', description='xtop 100: subxacts: 101'),
                   record('Transaction', 'PREPARE', 100),
                   record('Transaction', 'COMMIT_PREPARED', description='100: 2026-09-17')]
        self.assertEqual(self.status(records), 'unchanged')
        self.assertEqual(self.status(records, 20), 'changed')

    def test_shared_catalogs(self):
        self.assertEqual(self.status([record('Heap', 'UPDATE', 100, [0])]), 'uncertain')

    def test_unknown_record_fails_closed(self):
        self.assertEqual(self.status([record('custom128', 'UNKNOWN', 100, [20])]), 'uncertain')

    def test_storage_without_blocks(self):
        records = [record('Storage', 'CREATE', 100, description='base/20/123'),
                   record('Storage', 'TRUNCATE', 100, description='base/20/123 to 0 blocks flags 7'),
                   record('Transaction', 'COMMIT', 100)]
        self.assertEqual(self.status(records), 'unchanged')
        self.assertEqual(self.status(records, 20), 'changed')

    def test_dependency_records_follow_database(self):
        shared = record('Heap', 'INSERT', 100, [0])
        shared['dependency_only'] = True
        records = [shared, record('Heap', 'INSERT', 100, [20]), record('Transaction', 'COMMIT', 100)]
        self.assertEqual(self.status(records), 'unchanged')
        self.assertEqual(self.status(records, 20), 'changed')
        self.assertEqual(self.status([shared]), 'uncertain')

    def test_write_free_commit_allocated_after_boundary(self):
        records = [record('Standby', 'RUNNING_XACTS', description='nextXid 100 oldestRunningXid 90'),
                   record('Transaction', 'COMMIT', 101)]
        self.assertEqual(self.status(records), 'unchanged')
        self.assertEqual(wal.classify([record('Transaction', 'COMMIT', 101)], 10, 100)[0], 'unchanged')
        self.assertEqual(wal.classify([record('Transaction', 'COMMIT', 99)], 10, 100)[0], 'uncertain')

    def test_later_snapshot_does_not_clear_earlier_unknown_commit(self):
        records = [record('Transaction', 'COMMIT', 99),
                   record('Standby', 'RUNNING_XACTS', description='nextXid 100 oldestRunningXid 90')]
        self.assertEqual(self.status(records), 'uncertain')

    def test_abort_alone_does_not_commit_changes(self):
        self.assertEqual(self.status([record('Transaction', 'ABORT', 100)]), 'unchanged')

    def test_maintenance(self):
        self.assertEqual(self.status([record('Heap2', 'PRUNE_ON_ACCESS', databases=[10]),
                                      record('XLOG', 'CHECKPOINT_REDO')]), 'unchanged')


if __name__ == '__main__':
    unittest.main()
