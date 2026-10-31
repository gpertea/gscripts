#!/usr/bin/env python3
"""Exercise account selectors with fictional homes and an executable fixture."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

CODY = Path(__file__).resolve().parents[1] / 'cody'
MANAGER = CODY.with_name('cody-account')


@unittest.skipUnless(shutil.which('cc'), 'a C compiler is required for the executable fixture')
class Accounts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.home = self.root / 'home'
        self.home.mkdir()
        self.store = self.home / '.codex'
        self.store.mkdir()
        self.bin = self.root / 'bin'
        self.bin.mkdir()
        source = self.root / 'fixture.c'
        source.write_text('''#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char **argv) {
    for (int i = 1; i < argc; i++) if (!strcmp(argv[i], "--version")) {
        puts("codex-cli 1.2.3"); return 0;
    }
    const char *home = getenv("CODEX_HOME");
    printf("HOME=%s\\n", home ? home : "");
    const char *display = getenv("DISPLAY"), *xauth = getenv("XAUTHORITY");
    const char *browser = getenv("LIN_BROWSER_USE_AUTH_DISPLAY");
    printf("DISPLAY=%s\\nXAUTHORITY=%s\\nBROWSER=%s\\n", display ? display : "", xauth ? xauth : "", browser ? browser : "");
    for (int i = 1; i < argc; i++) printf("ARG=%s\\n", argv[i]);
    return 0;
}
''')
        subprocess.run(['cc', str(source), '-o', str(self.bin / 'codex')], check=True)
        self.env = {'HOME': str(self.home), 'HOSTNAME': 'fixture',
                    'PATH': str(self.bin) + ':' + os.environ['PATH']}
        subprocess.run([str(MANAGER), '--home', str(self.store), '--init'], env=self.env, check=True)
        (self.store / 'auth.json').write_text('fictional-work')
        (self.store / 'auth.json').chmod(0o600)

    def tearDown(self):
        self.tmp.cleanup()

    def run_cody(self, *args, **extra):
        return subprocess.run([str(CODY), '--noupdate', *args], env={**self.env, **extra},
                              capture_output=True, text=True, timeout=10)

    def test_both_accounts_use_shared_home(self):
        for selector in ('g', 'b'):
            result = self.run_cody('-' + selector, 'resume')
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('ARG=--no-daemon', result.stdout)
            self.assertIn('ARG=cli_auth_credentials_store="file"', result.stdout)
            self.assertFalse((self.home / ('.codex-' + selector)).exists())
            state = json.loads((self.store / 'cody-accounts/state.json').read_text())
            self.assertEqual(state['active'], selector)
        result = self.run_cody('resume')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.store / 'auth.json').read_text(), 'fictional-work')

    def test_aliases_select_accounts(self):
        for name, account in (('gcody', 'g'), ('bcody', 'b')):
            alias = CODY.with_name(name)
            result = subprocess.run([str(alias), '--noupdate', 'resume'], env=self.env,
                                    capture_output=True, text=True, timeout=10)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(json.loads((self.store / 'cody-accounts/state.json').read_text())['active'], account)

    def test_conflicting_accounts_are_rejected(self):
        self.assertEqual(self.run_cody('-g', '-b').returncode, 2)

    def test_external_credentials_are_rejected(self):
        for key in ('CODEX_ACCESS_TOKEN', 'CODEX_API_KEY', 'OPENAI_API_KEY', 'OPENAI_IDENTITY_TOKEN_FILE'):
            with self.subTest(key=key):
                self.assertEqual(self.run_cody('-g', **{key: 'fixture'}).returncode, 2)

    def test_nonprimary_home_is_rejected_for_selection(self):
        self.assertEqual(self.run_cody('-g', CODEX_HOME='/fictional/elsewhere').returncode, 2)

    def test_equivalent_primary_home_paths_remain_coordinated(self):
        alias = self.root / 'home-alias'
        alias.symlink_to(self.store, target_is_directory=True)
        result = self.run_cody('-g', 'resume', CODEX_HOME=str(alias))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads((self.store / 'cody-accounts/state.json').read_text())['active'], 'g')
        self.assertIn('HOME=' + str(alias), result.stdout)

    def test_storage_overrides_are_rejected(self):
        for arg in ('sqlite_home="/fictional"', 'cli_auth_credentials_store="keyring"'):
            with self.subTest(arg=arg):
                self.assertEqual(self.run_cody('-g', '-c', arg).returncode, 2)

    def test_managed_session_preserves_display_settings(self):
        result = self.run_cody('resume', DISPLAY=':12', XAUTHORITY='/fictional/auth',
                               LIN_BROWSER_USE_AUTH_DISPLAY=':88')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('DISPLAY=:12', result.stdout)
        self.assertIn('XAUTHORITY=/fictional/auth', result.stdout)
        self.assertIn('BROWSER=:88', result.stdout)

    def test_other_codex_home_is_unchanged(self):
        result = self.run_cody('resume', CODEX_HOME='/fictional/home')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('HOME=/fictional/home', result.stdout)
        self.assertEqual((self.store / 'auth.json').read_text(), 'fictional-work')

    def test_help_does_not_switch(self):
        result = self.run_cody('-g', '--help')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.store / 'auth.json').read_text(), 'fictional-work')
        self.assertEqual(json.loads((self.store / 'cody-accounts/state.json').read_text())['active'], 'default')


if __name__ == '__main__':
    unittest.main()
