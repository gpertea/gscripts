#!/usr/bin/env python3
"""Exercise embedded launches and update consent using isolated PTY doubles."""

import errno
import json
import os
from pathlib import Path
import pty
import select
import shutil
import subprocess
import sys
import tempfile
import time
import unittest


REPO = Path(__file__).resolve().parents[1]


class CodyUpdateTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="cody-update-test-")
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.bin = self.root / "bin"
        self.bin.mkdir()
        self.cody = self.root / "cody"
        self.cody.write_bytes((REPO / "cody").read_bytes())
        self.cody.chmod(0o755)
        self.log = self.root / "calls"
        self.version = self.root / "version"
        self.version.write_text("0.159.2")
        self.environment = dict(os.environ)
        self.environment.update(PATH=f"{self.bin}:/usr/bin:/bin", HOME=str(self.root),
            HOSTNAME="fictional-client", CODY_TEST_ROOT=str(self.root), CODY_TEST_LATEST="0.159.3")
        for key in ("DISPLAY", "CODEX_HOME", "LIN_BROWSER_USE_AUTH_DISPLAY", "XAUTHORITY",
                    "XDG_RUNTIME_DIR", "DBUS_SESSION_BUS_ADDRESS", "FORCE_COLOR", "NO_COLOR"):
            self.environment.pop(key, None)
        self.write_double("uname", "uname", "print(os.environ.get('CODY_TEST_OS', 'Linux'))")
        self.write_double("module", "module", "")
        for name in ("pgrep", "ps"):
            self.write_double(name, "FORBIDDEN-" + name, "raise SystemExit(97)")
        for name in ("cody-display", "cody-resume"):
            self.write_double(self.root / name, "FORBIDDEN-" + name, "raise SystemExit(97)")
        self.write_codex(self.bin / "codex", "codex")
        self.write_double("npm", "npm", '''
if args[0] == 'view':
    if os.environ.get('CODY_TEST_CHECK_FAIL') == '1': raise SystemExit(6)
    print(os.environ['CODY_TEST_LATEST'])
elif args[0] == 'install':
    if os.environ.get('CODY_TEST_UPDATE_FAIL') == '1': raise SystemExit(7)
    (root/'version').write_text(os.environ['CODY_TEST_LATEST'])
else:
    raise SystemExit(98)
''')

    def tearDown(self):
        self.assertFalse([row for row in self.calls() if row['name'].startswith('FORBIDDEN-')])
        for row in self.calls():
            if row['name'] in ('codex', 'standalone'):
                self.assertNotEqual(row['args'][:1], ['app-server'])

    def write_double(self, name, label, body):
        path = Path(name) if isinstance(name, Path) else self.bin / name
        path.parent.mkdir(parents=True, exist_ok=True)
        preamble = f'''#!{sys.executable}
import json,os,sys
from pathlib import Path
root=Path(os.environ['CODY_TEST_ROOT']);args=sys.argv[1:]
fields=('PATH','PSModulePath','DISPLAY','LIN_BROWSER_USE_AUTH_DISPLAY','XAUTHORITY',
        'XDG_RUNTIME_DIR','DBUS_SESSION_BUS_ADDRESS','FORCE_COLOR','NO_COLOR')
with (root/'calls').open('a') as stream:
    stream.write(json.dumps({{'name':{label!r},'args':args,'cwd':os.getcwd(),
        'environment':{{key:os.environ.get(key) for key in fields}}}})+'\\n')
'''
        path.write_text(preamble + body + "\n")
        path.chmod(0o755)

    def write_codex(self, path, label, version_file="version"):
        body = '''
if args == ['--version']:
    if (os.environ.get('CODY_TEST_POST_UPDATE_VERSION_FAIL') == '1' and
            (root/'VERSION_FILE').read_text() == os.environ['CODY_TEST_LATEST']):
        raise SystemExit(23)
    print('codex-cli '+(root/'VERSION_FILE').read_text())
elif args and args[0] == 'app-server':
    raise SystemExit(97)
else:
    print('LAUNCH '+' '.join(args))
    raise SystemExit(int(os.environ.get('CODY_TEST_LAUNCH_EXIT','0')))
'''
        self.write_double(path, label, body.replace('VERSION_FILE', version_file))

    def windows(self):
        self.windows_root = self.root / 'Windows Root'
        self.system32 = self.windows_root / 'System32'
        self.local_data = self.root / 'Local App Data'
        self.native_cli = self.local_data / 'Programs/OpenAI/Codex/bin/codex.exe'
        self.write_codex(self.native_cli, 'standalone')
        (self.root/'npm-version').write_text('9.9.9')
        self.write_codex(self.bin/'codex', 'codex', version_file='npm-version')
        self.powershell = self.system32 / 'WindowsPowerShell/v1.0/powershell.exe'
        self.environment.update(CODY_TEST_OS='MINGW64_NT-10.0', SystemRoot=r'C:\Windows',
            LOCALAPPDATA=r'C:\Fictional User\AppData\Local', PSModulePath='fictional-inherited-modules',
            CODY_TEST_WINDOWS_ROOT=str(self.windows_root), CODY_TEST_LOCAL_DATA=str(self.local_data),
            CODY_TEST_LOCAL_NATIVE=r'C:\Fictional User\AppData\Local')
        self.write_double('cygpath', 'cygpath', '''
assert args[0] == '-u'
if args[1] == os.environ['SystemRoot']: print(os.environ['CODY_TEST_WINDOWS_ROOT'])
elif args[1] == os.environ['CODY_TEST_LOCAL_NATIVE']: print(os.environ['CODY_TEST_LOCAL_DATA'])
else: raise SystemExit(96)
''')
        self.write_double(self.powershell, 'powershell', '''
command=args[args.index('-Command')+1]
if 'GetFolderPath' in command:
    print(os.environ['CODY_TEST_LOCAL_NATIVE'])
elif 'registry.npmjs.org' in command:
    if os.environ.get('CODY_TEST_CHECK_FAIL') == '1': raise SystemExit(6)
    print(os.environ['CODY_TEST_LATEST'])
elif 'install.ps1' in command:
    if os.environ.get('CODY_TEST_UPDATE_FAIL') == '1': raise SystemExit(7)
    (root/'version').write_text(os.environ['CODY_TEST_LATEST'])
else:
    raise SystemExit(95)
''')

    def run_cody(self, arguments=(), answers=b"", interactive=True):
        if not interactive:
            result = subprocess.run([str(self.cody), *arguments], env=self.environment,
                input=answers, capture_output=True, timeout=10, cwd=self.root)
            return result.returncode, (result.stdout + result.stderr).decode()

        ## use a real terminal while every executable and update remains fictional
        master, slave = pty.openpty()
        process = subprocess.Popen([str(self.cody), *arguments], env=self.environment,
            stdin=slave, stdout=slave, stderr=slave, cwd=self.root)
        os.close(slave)
        output = bytearray()
        try:
            if answers:
                os.write(master, answers)
            deadline = time.monotonic() + 10
            while time.monotonic() < deadline:
                ready, _, _ = select.select([master], [], [], 0.1)
                if ready:
                    try:
                        chunk = os.read(master, 65536)
                    except OSError as error:
                        if error.errno == errno.EIO:
                            break
                        raise
                    if not chunk:
                        break
                    output.extend(chunk)
            else:
                self.fail("cody did not complete within ten seconds")
            return process.wait(timeout=2), output.decode()
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
            os.close(master)

    def calls(self, name=None):
        rows = [json.loads(line) for line in self.log.read_text().splitlines()] if self.log.exists() else []
        return rows if name is None else [row for row in rows if row['name'] == name]

    def launches(self):
        return [row for row in self.calls() if row['name'] in ('codex', 'standalone') and row['args'] != ['--version']]

    def assert_embedded(self, launch):
        self.assertEqual(launch['args'][:2], ['--yolo', '--no-daemon'])
        self.assertEqual(launch['cwd'], str(self.root))
        self.assertNotIn('--remote', launch['args'])

    def test_accepted_update_requires_one_consent_then_embedded_launch(self):
        code, output = self.run_cody(answers=b"y\n")
        self.assertEqual(code, 0, output)
        self.assertEqual(self.version.read_text(), '0.159.3')
        self.assertEqual([row['args'][0] for row in self.calls('npm')], ['view', 'install'])
        self.assertEqual(self.calls('npm')[1]['args'], ['install', '--global', '@openai/codex@latest'])
        self.assertEqual(output.count('[y/N]'), 1)
        self.assertIn('Installed CLI: 0.159.3', output)
        self.assertNotIn('Running app-server', output)
        self.assertNotIn('Daemon package', output)
        self.assert_embedded(self.launches()[0])
        self.assertIn('check_for_update_on_startup=false', self.launches()[0]['args'])

    def test_update_declined(self):
        code, output = self.run_cody(answers=b'n\n')
        self.assertEqual(code, 0, output)
        self.assertEqual([row['args'][0] for row in self.calls('npm')], ['view'])
        self.assertEqual(self.version.read_text(), '0.159.2')
        self.assert_embedded(self.launches()[0])
        self.assertIn('check_for_update_on_startup=false', self.launches()[0]['args'])

    def test_noupdate_disables_both_checks_with_other_wrapper_options(self):
        code, output = self.run_cody(('--noupdate', '--docs', 'resume'))
        self.assertEqual(code, 0, output)
        self.assertFalse(self.calls('npm'))
        self.assert_embedded(self.launches()[0])
        self.assertIn('check_for_update_on_startup=false', self.launches()[0]['args'])
        self.assertIn('mcp_servers.openaiDeveloperDocs.enabled=true', self.launches()[0]['args'])
        self.assertEqual(self.launches()[0]['args'][-1], 'resume')

    def test_failed_cli_update_does_not_launch(self):
        self.environment['CODY_TEST_UPDATE_FAIL'] = '1'
        code, output = self.run_cody(answers=b'y\n')
        self.assertEqual(code, 7, output)
        self.assertFalse(self.launches())
        self.assertEqual(self.version.read_text(), '0.159.2')

    def test_offline_check_still_launches_embedded(self):
        self.environment['CODY_TEST_CHECK_FAIL'] = '1'
        code, output = self.run_cody()
        self.assertEqual(code, 0, output)
        self.assertIn('update check unavailable', output)
        self.assertEqual([row['args'][0] for row in self.calls('npm')], ['view'])
        self.assert_embedded(self.launches()[0])

    def test_current_and_newer_cli_are_not_updated(self):
        for version in ('0.159.3', '0.160.0'):
            with self.subTest(version=version):
                self.log.unlink(missing_ok=True)
                self.version.write_text(version)
                code, output = self.run_cody()
                self.assertEqual(code, 0, output)
                self.assertEqual([row['args'][0] for row in self.calls('npm')], ['view'])
                self.assertNotIn('update available', output)
                self.assertNotIn('Installed CLI:', output)
                self.assert_embedded(self.launches()[0])

    def test_noninteractive_launch_does_no_network_work(self):
        code, output = self.run_cody(interactive=False)
        self.assertEqual(code, 0, output)
        self.assertFalse(self.calls('npm'))
        self.assert_embedded(self.launches()[0])

    def test_explicit_update_installs_reports_and_exits_in_both_terminal_modes(self):
        for interactive in (False, True):
            with self.subTest(interactive=interactive):
                self.log.unlink(missing_ok=True)
                code, output = self.run_cody(('--update',), interactive=interactive)
                self.assertEqual(code, 0, output)
                self.assertEqual([row['args'][0] for row in self.calls('npm')], ['install'])
                self.assertIn('Installed CLI: 0.159.3', output)
                self.assertFalse(self.launches())
                self.assertNotIn('[y/N]', output)

    def test_explicit_update_failure_preserves_exit_code(self):
        self.environment['CODY_TEST_UPDATE_FAIL'] = '1'
        code, output = self.run_cody(('--update',), interactive=False)
        self.assertEqual(code, 7, output)
        self.assertFalse(self.launches())

    def test_post_install_version_failure_is_not_reported_as_success(self):
        self.environment['CODY_TEST_POST_UPDATE_VERSION_FAIL'] = '1'
        for explicit in (True, False):
            with self.subTest(explicit=explicit):
                self.log.unlink(missing_ok=True); self.version.write_text('0.159.2')
                args = ('--update',) if explicit else ()
                code, output = self.run_cody(args, answers=b'y\n' if not explicit else b'',
                                             interactive=not explicit)
                self.assertEqual(code, 23, output)
                self.assertEqual(self.version.read_text(), '0.159.3')
                self.assertFalse(self.launches())
                self.assertNotIn('Installed CLI:', output)
                self.assertTrue(any(row['args'][0] == 'install' for row in self.calls('npm')))

    def test_help_and_version_skip_network_and_embedded_flag(self):
        for argument in ('--help', '--version'):
            with self.subTest(argument=argument):
                self.log.unlink(missing_ok=True)
                code, output = self.run_cody((argument,))
                self.assertEqual(code, 0, output)
                self.assertFalse(self.calls('npm'))
                call = [row for row in self.calls('codex') if row['args'] != ['--version']][0]
                self.assertNotIn('--no-daemon', call['args'])

    def test_conflicting_options_make_no_update_or_launch(self):
        code, output = self.run_cody(('--noupdate', '--update'))
        self.assertEqual(code, 2, output)
        self.assertFalse(self.calls('npm'))
        self.assertFalse(self.launches())

    def test_native_administrative_commands_skip_updates_and_daemon_flag(self):
        commands = ('exec', 'e', 'review', 'login', 'logout', 'mcp', 'plugin', 'app-server',
            'remote-control', 'completion', 'update', 'doctor', 'sandbox', 'debug', 'apply',
            'a', 'cloud', 'exec-server', 'features', 'help', 'agents', 'queue', 'archive',
            'delete', 'unarchive', 'migrate-rollouts', 'migrate-nfs-db')
        for command in commands:
            with self.subTest(command=command):
                self.log.unlink(missing_ok=True)
                code, output = self.run_cody((command,), answers=b'n\n')
                self.assertEqual(code, 0, output)
                self.assertFalse(self.calls('npm'))
                self.assertEqual(self.launches()[0]['args'], ['--yolo', command])
                self.assertEqual(len(self.calls('codex')), 1)

    def test_retired_remote_options_are_rejected_before_update_work(self):
        for arguments in (('--attach', ':2'), ('--server-status',), ('--server-endpoint',),
                          ('--server-stop',), ('--remote', 'ws://fictional.invalid')):
            with self.subTest(arguments=arguments):
                self.log.unlink(missing_ok=True)
                code, output = self.run_cody(arguments)
                self.assertEqual(code, 2, output)
                self.assertFalse(self.calls('npm'))
                self.assertFalse(self.launches())

    def test_windows_standalone_update_uses_installer_environment_only(self):
        self.windows()
        code, output = self.run_cody(answers=b'y\n')
        self.assertEqual(code, 0, output)
        self.assertFalse(self.calls('npm'))
        self.assertFalse(self.calls('codex'))
        self.assertIn(str(self.native_cli), output)
        self.assert_embedded(self.launches()[0])
        self.assertEqual(self.launches()[0]['name'], 'standalone')
        ps = self.calls('powershell')
        check = next(row for row in ps if 'registry.npmjs.org' in row['args'][-1])
        installer = next(row for row in ps if 'install.ps1' in row['args'][-1])
        self.assertIn('-TimeoutSec 5', check['args'][-1])
        self.assertIn('CODEX_NON_INTERACTIVE', installer['args'][-1])
        self.assertEqual(check['environment']['PATH'], self.environment['PATH'])
        prefix = f'{self.system32}:{self.system32}/WindowsPowerShell/v1.0:'
        self.assertEqual(installer['environment']['PATH'], prefix + self.environment['PATH'])
        self.assertTrue(all(row['environment']['PSModulePath'] is None for row in ps))
        launch = self.launches()[0]
        self.assertEqual(launch['environment']['PATH'], self.environment['PATH'])
        self.assertEqual(launch['environment']['PSModulePath'], 'fictional-inherited-modules')
        self.assertIn('features.unified_exec=true', launch['args'])
        self.assertIn('allow_login_shell=false', launch['args'])

    def test_windows_current_standalone_is_silent_and_local_data_can_be_discovered(self):
        self.windows(); self.environment.pop('LOCALAPPDATA'); self.version.write_text('0.159.3')
        (self.root/'npm-version').write_text('0.100.0')
        code, output = self.run_cody()
        self.assertEqual(code, 0, output)
        self.assertNotIn('update available', output)
        self.assertNotIn('Installed CLI:', output)
        self.assertFalse(self.calls('codex'))
        self.assertTrue(any('GetFolderPath' in row['args'][-1] for row in self.calls('powershell')))
        self.assertEqual(self.launches()[0]['name'], 'standalone')
        self.assertTrue(all(row['environment']['PSModulePath'] is None for row in self.calls('powershell')))

    def test_windows_decline_offline_noupdate_and_installer_failure(self):
        self.windows()
        for mode in ('decline', 'offline', 'noupdate', 'failure'):
            with self.subTest(mode=mode):
                self.log.unlink(missing_ok=True); self.version.write_text('0.159.2')
                for key in ('CODY_TEST_CHECK_FAIL', 'CODY_TEST_UPDATE_FAIL'): self.environment.pop(key, None)
                if mode == 'offline': self.environment['CODY_TEST_CHECK_FAIL'] = '1'
                if mode == 'failure': self.environment['CODY_TEST_UPDATE_FAIL'] = '1'
                args = ('--noupdate',) if mode == 'noupdate' else ()
                answer = b'y\n' if mode == 'failure' else b'n\n'
                code, output = self.run_cody(args, answers=answer)
                self.assertEqual(code, 7 if mode == 'failure' else 0, output)
                self.assertFalse(self.calls('npm'))
                if mode == 'failure': self.assertFalse(self.launches())
                else:
                    self.assertEqual(self.launches()[0]['name'], 'standalone')
                    self.assert_embedded(self.launches()[0])
                    self.assertFalse([r for r in self.calls('powershell') if 'install.ps1' in r['args'][-1]])
                if mode == 'noupdate': self.assertFalse(self.calls('powershell'))

    def test_windows_explicit_update_and_path_fallback_preserve_direct_launch(self):
        self.windows()
        code, output = self.run_cody(('--update',), interactive=False)
        self.assertEqual(code, 0, output)
        self.assertFalse(self.calls('npm'))
        self.assertFalse(self.launches())
        self.assertIn('Installed CLI: 0.159.3', output)
        self.native_cli.unlink(); self.log.unlink()
        self.environment['DISPLAY'] = ':2'
        code, output = self.run_cody(('--noupdate',), interactive=False)
        self.assertEqual(code, 0, output)
        self.assertEqual(self.launches()[0]['name'], 'codex')
        self.assertIsNone(self.launches()[0]['environment']['LIN_BROWSER_USE_AUTH_DISPLAY'])

    def test_jhpce_displayless_update_prefers_user_cli(self):
        for hostname, display in (('transfer-01.cm.cluster', None), ('compute-001', '')):
            with self.subTest(hostname=hostname, display=display):
                self.log.unlink(missing_ok=True); self.version.write_text('0.159.2')
                self.environment['HOSTNAME'] = hostname
                if display is None: self.environment.pop('DISPLAY', None)
                else: self.environment['DISPLAY'] = display
                user_bin = self.root / '.local/bin'; user_bin.mkdir(parents=True, exist_ok=True)
                (user_bin/'codex').unlink(missing_ok=True); (user_bin/'codex').symlink_to(self.bin/'codex')
                code, output = self.run_cody(('resume',), answers=b'y\n')
                self.assertEqual(code, 0, output)
                rows = self.calls(); module = next(i for i,r in enumerate(rows) if r['name']=='module')
                install = next(i for i,r in enumerate(rows) if r['name']=='npm' and r['args'][0]=='install')
                self.assertLess(module, install)
                self.assertIn(f'({user_bin}/codex)', output)
                self.assert_embedded(self.launches()[0])
                self.assertEqual(self.launches()[0]['args'][-1], 'resume')

    def test_jhpce_recreates_node_local_helpers(self):
        local = Path(tempfile.mkdtemp(prefix='cody-codex-arg0-'))
        self.addCleanup(shutil.rmtree, local, ignore_errors=True); local.rmdir()
        test_uid = local.name.removeprefix('cody-codex-arg0-')
        self.write_double('id', 'id', f'print({test_uid!r})')
        home_tmp = self.root / '.codex/tmp'; home_tmp.mkdir(parents=True)
        (home_tmp/'arg0').symlink_to(local)
        self.environment['HOSTNAME'] = 'transfer-01.cm.cluster'
        code, output = self.run_cody(('--noupdate',))
        self.assertEqual(code, 0, output)
        self.assertEqual(local.stat().st_mode & 0o777, 0o700)
        self.assert_embedded(self.launches()[0])

    def test_jhpce_rejects_symlinked_local_helpers(self):
        local = Path(tempfile.mkdtemp(prefix='cody-codex-arg0-')); local.rmdir()
        local.symlink_to(self.root, target_is_directory=True); self.addCleanup(local.unlink)
        test_uid = local.name.removeprefix('cody-codex-arg0-')
        self.write_double('id', 'id', f'print({test_uid!r})')
        home_tmp = self.root / '.codex/tmp'; home_tmp.mkdir(parents=True)
        (home_tmp/'arg0').symlink_to(local)
        self.environment['HOSTNAME'] = 'compute-001'
        code, output = self.run_cody(('--noupdate',))
        self.assertEqual(code, 1, output)
        self.assertIn('expected a private helper directory', output)
        self.assertFalse(self.launches())


if __name__ == '__main__':
    unittest.main()
