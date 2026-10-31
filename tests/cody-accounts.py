#!/usr/bin/env python3
"""Account selection without real logins or credentials."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class Accounts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.bin = self.root / 'bin'
        self.bin.mkdir()
        self.log = self.root / 'record.json'
        self.wrapper = self.root / 'cody'
        self.wrapper.write_bytes((ROOT / 'cody').read_bytes())
        self.wrapper.chmod(0o700)
        for name in ('gcody', 'bcody'):
            (self.root / name).symlink_to('cody')
        for name in ('.codex-g', '.codex-b'):
            (self.root / name).mkdir(mode=0o700)
        fake = self.bin / 'codex'
        fake.write_text('''#!/usr/bin/env python3
import os,json,sys
with open(os.environ['TEST_RECORD'],'w') as out:
 json.dump({'args':sys.argv[1:],'home':os.getenv('CODEX_HOME'),
            'sqlite':os.getenv('CODEX_SQLITE_HOME')},out)
''')
        fake.chmod(0o700)
        self.env = {'PATH': str(self.bin) + os.pathsep + os.environ['PATH'],
                    'HOME': str(self.root), 'HOSTNAME': 'fixture', 'TEST_RECORD': str(self.log)}

    def call(self, entry='cody', *args, **overrides):
        env = {**self.env, **overrides}
        r = subprocess.run([str(self.root / entry), *args], env=env, capture_output=True, text=True)
        return r, json.loads(self.log.read_text()) if self.log.exists() else None

    def test_account_aliases_and_flags_select_separate_homes_for_login_and_resume(self):
        for entry, flags, suffix in [('cody', ('-g',), 'g'), ('gcody', (), 'g'),
                                      ('cody', ('-b',), 'b'), ('bcody', (), 'b')]:
            for command in (('login',), ('resume', '--last')):
                with self.subTest(entry=entry, command=command):
                    r, data = self.call(entry, *flags, '--noupdate', *command,
                                        CODEX_HOME='/fictional/other', CODEX_SQLITE_HOME='/fictional/other')
                    self.assertEqual(r.returncode, 0, r.stderr)
                    self.assertEqual(data['home'], str(self.root / ('.codex-' + suffix)))
                    self.assertEqual(data['sqlite'], data['home'])
                    self.assertIn('cli_auth_credentials_store="file"', data['args'])
                    self.assertIn('sqlite_home="' + data['home'] + '"', data['args'])
                    self.assertEqual('--no-daemon' in data['args'], command[0] == 'resume')
                    self.log.unlink()

    def test_unselected_home_and_environment_are_preserved(self):
        r, data = self.call('cody', '--noupdate', CODEX_HOME='/fictional/home',
                            CODEX_SQLITE_HOME='/fictional/db')
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual((data['home'], data['sqlite']), ('/fictional/home', '/fictional/db'))
        self.assertNotIn('cli_auth_credentials_store="file"', data['args'])

    def test_selectors_mix_with_wrapper_options_and_conflicts_are_rejected(self):
        r, data = self.call('cody', '--noupdate', '-g', '--docs')
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(data['home'], str(self.root / '.codex-g'))
        self.log.unlink()
        for entry, args in [('cody', ('-g', '-b')), ('gcody', ('-b',))]:
            r, data = self.call(entry, *args)
            self.assertEqual(r.returncode, 2)
            self.assertIsNone(data)

    def test_selected_account_rejects_overriding_its_storage(self):
        for args in [('-c', 'sqlite_home="/fictional/other"'),
                     ('--config', 'cli_auth_credentials_store="keyring"'),
                     ('--config=sqlite_home="/fictional/other"',),
                     ('--config=sqlite_home = "/fictional/other"',),
                     ('--config=  sqlite_home = "/fictional/other"',),
                     ('-c', ' "sqlite_home" = "/fictional/other"'),
                     ('-ccli_auth_credentials_store="keyring"',)]:
            with self.subTest(args=args):
                r, data = self.call('gcody', '--noupdate', *args)
                self.assertEqual(r.returncode, 2)
                self.assertIsNone(data)

    def test_missing_or_symlink_home_and_environment_credentials_fail_closed(self):
        for env in ({'CODEX_ACCESS_TOKEN': 'fictional'}, {'CODEX_API_KEY': 'fictional'},
                    {'OPENAI_IDENTITY_TOKEN_FILE': '/fictional/token'}):
            r, data = self.call('gcody', **env)
            self.assertEqual(r.returncode, 2)
            self.assertIsNone(data)
        (self.root / '.codex-g').rmdir()
        r, _ = self.call('gcody')
        self.assertNotEqual(r.returncode, 0)
        (self.root / '.codex-g').symlink_to(self.root / '.codex-b')
        r, _ = self.call('gcody')
        self.assertNotEqual(r.returncode, 0)


if __name__ == '__main__':
    unittest.main()
