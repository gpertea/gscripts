#!/usr/bin/env python3
"""Exercise crash-safe credential moves and process coordination using fixtures."""

import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import tempfile
import time
import unittest
from unittest import mock

MANAGER = Path(__file__).resolve().parents[1] / 'cody-account'
module = runpy.run_path(str(MANAGER))


def secret(path, value):
    path.write_text(value)
    path.chmod(0o600)


class Recovery(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.home = Path(self.temp.name) / '.codex'
        self.home.mkdir()
        module['initialize'](self.home)
        self.root = module['private_dir'](self.home)
        self.auth = self.home / 'auth.json'
        self.old = self.root / 'default.json'
        self.new = self.root / 'g.json'

    def tearDown(self):
        self.temp.cleanup()

    def journal(self, source_id, target_id):
        module['write_json'](self.root / 'journal.json',
                             {'from': 'default', 'to': 'g', 'source_id': source_id, 'target_id': target_id})

    def test_recovery_at_every_file_move(self):
        for source_present in (False, True):
            for target_present in (False, True):
                for stage in range(6):
                    with self.subTest(source=source_present, target=target_present, stage=stage):
                        for path in self.root.iterdir():
                            path.unlink()
                        if self.auth.exists():
                            self.auth.unlink()
                        module['write_json'](self.root / 'state.json', {'active': 'default'})
                        if source_present:
                            secret(self.auth, 'fictional-work')
                        if target_present:
                            secret(self.new, 'fictional-g')
                        self.journal(module['identity'](self.auth), module['identity'](self.new))
                        if stage >= 1 and source_present:
                            os.link(self.auth, self.old)
                        if stage >= 2 and source_present:
                            self.auth.unlink()
                        if stage >= 3 and target_present:
                            os.link(self.new, self.auth)
                        if stage >= 4 and target_present:
                            self.new.unlink()
                        if stage == 5:
                            module['write_json'](self.root / 'state.json', {'active': 'g'})
                        module['recover'](self.home, self.root)
                        self.assertEqual(self.old.exists(), source_present)
                        self.assertEqual(self.auth.exists(), target_present)
                        if source_present:
                            self.assertEqual(self.old.read_text(), 'fictional-work')
                        if target_present:
                            self.assertEqual(self.auth.read_text(), 'fictional-g')
                        self.assertFalse(self.new.exists())
                        self.assertFalse((self.root / 'journal.json').exists())
                        self.assertEqual(module['state'](self.root)['active'], 'g')

    def test_logged_out_account_and_conflicts_are_preserved(self):
        self.journal(None, None)
        module['recover'](self.home, self.root)
        self.assertFalse(self.auth.exists())
        self.assertEqual(module['state'](self.root)['active'], 'g')
        module['write_json'](self.root / 'state.json', {'active': 'default'})
        secret(self.auth, 'fictional-unexpected')
        secret(self.old, 'fictional-work')
        secret(self.new, 'fictional-g')
        self.journal(module['identity'](self.old), module['identity'](self.new))
        with self.assertRaises(module['AccountError']):
            module['recover'](self.home, self.root)
        self.assertEqual(self.auth.read_text(), 'fictional-unexpected')
        self.assertEqual(self.old.read_text(), 'fictional-work')
        self.assertEqual(self.new.read_text(), 'fictional-g')

    def test_status_reports_label_without_credential_contents(self):
        secret(self.auth, 'fictional-secret-never-display')
        result = subprocess.run([str(MANAGER), '--home', str(self.home), '--status'],
                                capture_output=True, text=True, timeout=5)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('Active account: default', result.stdout)
        self.assertNotIn('fictional-secret-never-display', result.stdout + result.stderr)

    def test_administrator_keyring_setting_blocks_file_swaps(self):
        policy = self.root / 'requirements.toml'
        policy.write_text('cli_auth_credentials_store = "keyring"\n')
        with self.assertRaises(module['AccountError']):
            module['require_file_store']([policy])
        policy.write_text('cli_auth_credentials_store = "file"\n')
        module['require_file_store']([policy])

    def test_logged_out_and_saved_account_recovery(self):
        for source_present, target_present in ((False, True), (True, False), (False, False)):
            with self.subTest(source_present=source_present, target_present=target_present):
                for path in self.root.iterdir():
                    path.unlink()
                if self.auth.exists():
                    self.auth.unlink()
                module['write_json'](self.root / 'state.json', {'active': 'default'})
                if source_present:
                    secret(self.auth, 'fictional-work')
                if target_present:
                    secret(self.new, 'fictional-g')
                self.journal(module['identity'](self.auth), module['identity'](self.new))
                module['recover'](self.home, self.root)
                self.assertEqual(self.old.exists(), source_present)
                self.assertEqual(self.auth.exists(), target_present)
                if target_present:
                    self.assertEqual(self.auth.read_text(), 'fictional-g')

    def test_symlink_and_open_permissions_fail_closed(self):
        self.auth.symlink_to(self.root / 'missing')
        with self.assertRaises(module['AccountError']):
            module['swap'](self.home, self.root, 'default', 'g')
        self.auth.unlink()
        self.auth.write_text('fictional-work')
        self.auth.chmod(0o644)
        with self.assertRaises(module['AccountError']):
            module['swap'](self.home, self.root, 'default', 'g')


class RunningSessions(unittest.TestCase):
    def test_shared_lock_blocks_switch(self):
        with tempfile.TemporaryDirectory() as temporary:
            home = Path(temporary) / '.codex'
            home.mkdir()
            module['initialize'](home)
            root = module['private_dir'](home)
            secret(home / 'auth.json', 'fictional-work')
            first = subprocess.Popen([str(MANAGER), '--home', str(home), '--account', 'default',
                                      '--', 'sleep', '30'])
            try:
                for _ in range(100):
                    fd = module['lockfile'](root)
                    try:
                        import fcntl
                        try:
                            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        except BlockingIOError:
                            break
                    finally:
                        os.close(fd)
                    time.sleep(.02)
                else:
                    self.fail('first session did not hold the shared lock')
                changed = subprocess.run([str(MANAGER), '--home', str(home), '--account', 'g',
                                          '--', 'true'], capture_output=True, text=True, timeout=5)
                self.assertNotEqual(changed.returncode, 0)
                self.assertIn('close it before switching', changed.stderr)
                self.assertEqual((home / 'auth.json').read_text(), 'fictional-work')
            finally:
                first.terminate()
                first.wait(timeout=5)
            changed = subprocess.run([str(MANAGER), '--home', str(home), '--account', 'g',
                                      '--', 'true'], capture_output=True, text=True, timeout=5)
            self.assertEqual(changed.returncode, 0, changed.stderr)
            self.assertFalse((home / 'auth.json').exists())
            self.assertEqual((root / 'default.json').read_text(), 'fictional-work')

    @unittest.skipUnless(shutil.which('cc'), 'a C compiler is required for the npm fixture')
    def test_npm_launcher_uses_only_a_matching_native_binary(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / 'package'
            (root / 'bin').mkdir(parents=True)
            shim = root / 'bin/codex.js'
            shim.write_text('#!/usr/bin/env node\n')
            shim.chmod(0o700)
            (root / 'package.json').write_text('{"version":"1.2.3"}')
            vendor = root / 'node_modules/@openai/codex-linux-x64/vendor/x86_64-unknown-linux-musl'
            (vendor / 'bin').mkdir(parents=True)
            (vendor / 'codex-package.json').write_text('{"version":"1.2.3"}')
            source = root / 'fixture.c'
            source.write_text('#include <stdio.h>\nint main(void) { puts("codex-cli 1.2.3"); return 0; }\n')
            native = vendor / 'bin/codex'
            subprocess.run(['cc', str(source), '-o', str(native)], check=True)
            with mock.patch.dict(os.environ, {'CODEX_MANAGED_BY_BUN': '1'}):
                command = module['native_command']([str(shim), 'resume'])
                self.assertEqual(command, [str(native), 'resume'])
                self.assertEqual(os.environ['CODEX_MANAGED_BY_NPM'], '1')
                self.assertNotIn('CODEX_MANAGED_BY_BUN', os.environ)
                self.assertEqual(os.environ['CODEX_MANAGED_PACKAGE_ROOT'], str(root.resolve()))
                (root / 'package.json').write_text('{"version":"2.0.0"}')
                with self.assertRaises(module['AccountError']):
                    module['native_command']([str(shim), 'resume'])
                (root / 'package.json').write_text('{"version":"1.2.3"}')
                source.write_text('#include <stdio.h>\nint main(void) { puts("codex-cli 9.9.9"); return 0; }\n')
                subprocess.run(['cc', str(source), '-o', str(native)], check=True)
                with self.assertRaises(module['AccountError']):
                    module['native_command']([str(shim), 'resume'])

    @unittest.skipUnless(shutil.which('cc'), 'a C compiler is required for the signal fixture')
    def test_native_cli_receives_default_pipe_signal(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            home = root / '.codex'
            home.mkdir()
            module['initialize'](home)
            source = root / 'fixture.c'
            source.write_text('#include <signal.h>\nint main(void) { struct sigaction a; sigaction(SIGPIPE, 0, &a); return a.sa_handler == SIG_DFL ? 0 : 1; }\n')
            binary = root / 'codex'
            subprocess.run(['cc', str(source), '-o', str(binary)], check=True)
            result = subprocess.run([str(MANAGER), '--home', str(home), '--account', 'default',
                                     '--', str(binary)], capture_output=True, text=True, timeout=5)
            self.assertEqual(result.returncode, 0, result.stderr)

    @unittest.skipUnless(shutil.which('cc'), 'a C compiler is required for the process fixture')
    def test_unmanaged_codex_blocks_switch(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            home = root / '.codex'
            home.mkdir()
            module['initialize'](home)
            secret(home / 'auth.json', 'fictional-work')
            source = root / 'fixture.c'
            source.write_text('#include <unistd.h>\nint main(void) { sleep(30); return 0; }\n')
            binary = root / 'codex'
            subprocess.run(['cc', str(source), '-o', str(binary)], check=True)
            process = subprocess.Popen([str(binary)], env={'CODEX_HOME': str(home), 'HOME': str(root)})
            try:
                for _ in range(100):
                    if module['active_codex'](home):
                        break
                    time.sleep(.02)
                else:
                    self.fail('unmanaged process was not detected')
                result = subprocess.run([str(MANAGER), '--home', str(home), '--account', 'g',
                                         '--', 'true'], capture_output=True, text=True, timeout=5)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('another Codex session or server', result.stderr)
                self.assertEqual((home / 'auth.json').read_text(), 'fictional-work')
                login = subprocess.run([str(MANAGER), '--home', str(home), '--account', 'default',
                                        '--mutates-auth', '--', 'true'], capture_output=True, text=True, timeout=5)
                self.assertNotEqual(login.returncode, 0)
                self.assertEqual((home / 'auth.json').read_text(), 'fictional-work')
            finally:
                process.terminate()
                process.wait(timeout=5)


if __name__ == '__main__':
    unittest.main()
