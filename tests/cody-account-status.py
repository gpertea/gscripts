#!/usr/bin/env python3
"""Exercise live status context with fictional processes and private Mesh data."""

import json
import os
from pathlib import Path
import runpy
import shlex
import shutil
import sqlite3
import subprocess
import tempfile
import time
import unittest
from unittest import mock

MANAGER = Path(__file__).resolve().parents[1] / 'cody-account'
module = runpy.run_path(str(MANAGER))


@unittest.skipUnless(shutil.which('cc') and Path('/proc').is_dir(), 'requires Linux and cc')
class Status(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.build = tempfile.TemporaryDirectory()
        root = Path(cls.build.name)
        source = root / 'fixture.c'
        source.write_text('#include <unistd.h>\nint main(void) { for (;;) pause(); }\n')
        cls.binary = root / 'codex'
        subprocess.run(['cc', str(source), '-o', str(cls.binary)], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.build.cleanup()

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.home = self.root / '.codex'
        self.home.mkdir()
        self.work = self.root / 'fictional project' / 'src'
        self.work.mkdir(parents=True)
        self.mesh = self.root / 'mesh-state'
        self.mesh.mkdir()
        self.env = {'HOME': str(self.root), 'PATH': os.environ['PATH'],
                    'CODEX_HOME': str(self.home), 'MESH_STATE_DIR': str(self.mesh)}
        module['initialize'](self.home)
        self.process = subprocess.Popen([str(self.binary)], cwd=self.work, env=self.env)
        self.addCleanup(self.stop_process)
        self.wait_process(self.process.pid)

    def stop_process(self):
        self.process.terminate()
        self.process.wait(timeout=5)

    def wait_process(self, pid):
        for _ in range(100):
            if str(pid) in module['codex_processes'](self.home):
                return
            time.sleep(.01)
        self.fail('fictional Codex process was not detected')

    def status(self, *args, **environment):
        result = subprocess.run([str(MANAGER), *args], env={**self.env, **environment},
                                capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.stdout

    def catalog(self, records):
        path = self.mesh / 'catalog.sqlite'
        with sqlite3.connect(path) as db:
            db.executescript('PRAGMA user_version=2; CREATE TABLE meta(key TEXT PRIMARY KEY,value TEXT);'
                             'CREATE TABLE entities(kind TEXT,id TEXT,payload TEXT,updated REAL,PRIMARY KEY(kind,id));')
            db.execute('INSERT INTO meta VALUES (?,?)', ('node_id', json.dumps('fictional-local')))
            for index, record in enumerate(records):
                record = {'node_id': 'fictional-local', **record}
                db.execute('INSERT INTO entities VALUES (?,?,?,?)',
                           ('project', str(index), json.dumps(record), 0))
        return path

    def legacy(self):
        path = self.mesh / 'registry.sqlite'
        with sqlite3.connect(path) as db:
            db.executescript('PRAGMA user_version=1; CREATE TABLE projects(id TEXT,label TEXT);'
                             'CREATE TABLE checkouts(id TEXT,project_id TEXT,path TEXT);')
            db.execute('INSERT INTO projects VALUES (?,?)', ('fictional-project', 'Legacy label'))
            db.execute('INSERT INTO checkouts VALUES (?,?,?)',
                       ('fictional-checkout', 'fictional-project', str(self.work.parent)))
        return path

    def test_default_and_explicit_status_show_pid_and_cwd(self):
        for args in ((), ('--status',), ('--home', str(self.home), '--status')):
            with self.subTest(args=args):
                output = self.status(*args)
                self.assertIn(str(self.process.pid), output)
                self.assertIn(str(self.work), output)
                self.assertIn('WORKING DIRECTORY', output)
                self.assertIn('TMUX', output)
                self.assertNotIn('MESH PROJECT', output)

    def test_unmanaged_home_still_lists_process_context(self):
        shutil.rmtree(self.home / 'cody-accounts')
        output = self.status('--status')
        self.assertIn('Active account: unmanaged', output)
        self.assertIn(str(self.process.pid), output)
        self.assertIn(str(self.work), output)

    def test_catalog_uses_nearest_local_project_and_is_read_only(self):
        path = self.catalog([
            {'path': str(self.root), 'label': 'Outer project'},
            {'path': str(self.work.parent), 'label': 'Fictional project'},
            {'path': str(self.work), 'label': 'Remote wrong', 'node_id': 'fictional-remote'},
            {'path': str(self.work) + '-other', 'label': 'Prefix wrong'},
        ])
        before = path.read_bytes()
        output = self.status('--status')
        self.assertIn('MESH PROJECT (CWD)', output)
        self.assertIn('Fictional project', output)
        self.assertNotIn('Remote wrong', output)
        self.assertNotIn('Outer project', output)
        self.assertNotIn('Prefix wrong', output)
        self.assertIn('cached', output)
        self.assertEqual(path.read_bytes(), before)

    def test_legacy_registry_and_catalog_precedence(self):
        self.legacy()
        self.assertIn('Legacy label', self.status('--status'))
        self.catalog([{'path': str(self.work.parent), 'label': 'Catalog label'}])
        output = self.status('--status')
        self.assertIn('Catalog label', output)
        self.assertNotIn('Legacy label', output)

    def test_bad_catalog_falls_back_to_registry(self):
        (self.mesh / 'catalog.sqlite').write_text('fictional broken database')
        self.legacy()
        output = self.status('--status')
        self.assertIn('Legacy label', output)
        self.assertIn('unavailable', output)
        self.assertIn(str(self.work), output)

    def test_locked_catalog_preserves_process_output(self):
        path = self.catalog([{'path': str(self.work), 'label': 'Locked project'}])
        with sqlite3.connect(path) as db:
            db.execute('BEGIN EXCLUSIVE')
            output = self.status('--status')
        self.assertIn(str(self.work), output)
        self.assertIn('unavailable', output)

    @unittest.skipIf(os.getuid() == 0, 'root bypasses directory permission checks')
    def test_unreadable_mesh_directory_preserves_process_output(self):
        self.legacy()
        self.addCleanup(lambda: self.mesh.chmod(0o700))
        self.mesh.chmod(0)
        output = self.status('--status')
        self.assertIn(str(self.work), output)
        self.assertIn('unavailable', output)

    def test_control_characters_in_project_labels_are_escaped(self):
        self.catalog([{'path': str(self.work), 'label': 'Fictional\n\x1b[31mproject'}])
        output = self.status('--status')
        self.assertIn(r'Fictional\n\x1b[31mproject', output)
        self.assertNotIn('\x1b', output)

    def test_no_processes_and_other_homes(self):
        other = self.root / 'other-home'
        other.mkdir()
        self.stop_process()
        self.process = subprocess.Popen([str(self.binary)], cwd=self.work,
                                        env={**self.env, 'CODEX_HOME': str(other)})
        for _ in range(100):
            if str(self.process.pid) in module['codex_processes'](other):
                break
            time.sleep(.01)
        else:
            self.fail('other-home fixture did not start')
        output = self.status('--status')
        self.assertIn('this home: none', output)
        self.assertNotIn(str(self.process.pid), output)

    def test_unsupported_catalog_and_malformed_payload_do_not_hide_processes(self):
        path = self.catalog([{'path': str(self.work), 'label': 'Wrong version'}])
        with sqlite3.connect(path) as db:
            db.execute('PRAGMA user_version=999')
        output = self.status('--status')
        self.assertNotIn('Wrong version', output)
        self.assertIn(str(self.work), output)
        self.assertIn('unavailable', output)
        with sqlite3.connect(path) as db:
            db.execute('PRAGMA user_version=2')
            db.execute("UPDATE entities SET payload='[]'")
        self.assertIn(str(self.work), self.status('--status'))

    def test_state_discovery_xdg_default_and_missing(self):
        self.legacy()
        del self.env['MESH_STATE_DIR']
        xdg = self.root / 'xdg'
        xdg.mkdir()
        self.mesh.rename(xdg / 'codex-mesh')
        self.assertIn('Legacy label', self.status('--status', XDG_STATE_HOME=str(xdg)))
        default = self.root / '.local/state'
        default.mkdir(parents=True)
        (xdg / 'codex-mesh').rename(default / 'codex-mesh')
        self.assertIn('Legacy label', self.status('--status'))
        missing = self.root / 'absent'
        self.assertIn(str(self.work), self.status('--status', MESH_STATE_DIR=str(missing)))
        self.assertFalse(missing.exists())

    @unittest.skipUnless(shutil.which('tmux'), 'requires tmux')
    def test_real_tmux_private_socket_and_renamed_session(self):
        socket = self.root / 'tmux.sock'
        command = ['tmux', '-S', str(socket)]
        launch = ['env', 'CODEX_HOME=' + str(self.home), str(self.binary)]
        subprocess.run([*command, '-f', '/dev/null', 'new-session', '-d', '-s', 'fictional',
                        '-c', str(self.work), 'exec ' + shlex.join(launch)],
                       env=self.env, check=True)
        self.addCleanup(lambda: subprocess.run([*command, 'kill-server'], capture_output=True))
        pane_pid = int(subprocess.check_output([*command, 'display-message', '-p', '#{pane_pid}']))
        self.wait_process(pane_pid)
        output = self.status('--status')
        self.assertIn('fictional:0.0', output)
        self.assertIn('%0', output)
        # a valid but inherited-outside-tmux hint must not claim another process's pane.
        outsider = subprocess.Popen([str(self.binary)], cwd=self.work,
                                    env={**self.env, 'TMUX': str(socket) + ',1,0', 'TMUX_PANE': '%0'})
        try:
            self.wait_process(outsider.pid)
            output = self.status('--status')
            row = next(line for line in output.splitlines() if line.startswith(str(outsider.pid)))
            self.assertIn('unavailable', row)
            self.assertNotIn('fictional:0.0', row)
        finally:
            outsider.terminate()
            outsider.wait(timeout=5)
        subprocess.run([*command, 'rename-session', '-t', 'fictional', 'renamed'], check=True)
        self.assertIn('renamed:0.0', self.status('--status'))

        # descendants of a shell in another pane must resolve through their ancestry.
        subprocess.run([*command, 'new-window', '-t', 'renamed', '-c', str(self.work),
                        shlex.join(launch) + ' & wait'], env=self.env, check=True)
        for _ in range(100):
            output = self.status('--status')
            if 'renamed:1.0' in output:
                break
            time.sleep(.01)
        self.assertIn('renamed:1.0', output)

    def test_stale_tmux_environment_is_not_reported_as_live(self):
        self.stop_process()
        self.process = subprocess.Popen([str(self.binary)], cwd=self.work,
                                        env={**self.env, 'TMUX': str(self.root / 'gone') + ',1,0',
                                             'TMUX_PANE': '%0'})
        self.wait_process(self.process.pid)
        output = self.status('--status')
        self.assertIn(str(self.work), output)
        self.assertIn('unavailable', output)

    def test_tmux_failure_and_metadata_errors_do_not_change_pid_detection(self):
        # an unreadable or vanished cwd must never remove a conservatively detected PID.
        globals_ = module['show_processes'].__globals__
        with mock.patch.dict(globals_, process_context=mock.Mock(return_value={
                'pid': str(self.process.pid), 'cwd': None, 'tmux': None})), \
                mock.patch('builtins.print') as output:
            module['show_processes'](self.home)
        rendered = '\n'.join(' '.join(str(x) for x in call.args) for call in output.call_args_list)
        self.assertIn(str(self.process.pid), rendered)
        self.assertIn('unavailable', rendered)

    def test_tmux_missing_and_timeout_leave_context_unknown(self):
        for error in (FileNotFoundError(), subprocess.TimeoutExpired('tmux', 1)):
            process = {'pid': '123', 'cwd': '/fictional', 'tmux': None,
                       'socket': '/fictional/socket', 'pane': '%0', 'ancestors': {123: 1}}
            with self.subTest(error=type(error).__name__), \
                    mock.patch.object(module['subprocess'], 'run', side_effect=error):
                module['tmux_context']([process])
            self.assertIsNone(process['tmux'])
            self.assertEqual(process['cwd'], '/fictional')


if __name__ == '__main__':
    unittest.main()
