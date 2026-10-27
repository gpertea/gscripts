#!/usr/bin/env python3
"""Exercise macOS routing with fake SSH and TurboVNC; no network or GUI."""

import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]
SSH = '''#!/usr/bin/env python3
import json, os, sys
args = sys.argv[1:]
with open(os.environ["VNC_TEST_LOG"], "a") as log:
    log.write(json.dumps(["ssh", *args]) + "\\n")
if "-G" in args:
    if os.environ.get("MOCK_CONFIG_FAIL"):
        sys.exit(7)
    host = args[-1].split("@")[-1]
    for key, default in [("hostname", host), ("port", "22"),
                         ("proxyjump", "none"), ("proxycommand", "none")]:
        print(key, os.environ.get("MOCK_" + key.upper(), default))
elif "-O" in args:
    sys.exit(0 if args[args.index("-O") + 1] == "exit" else
             (0 if os.environ.get("MOCK_REUSE") else 1))
elif "-fNT" in args:
    sys.exit(7 if os.environ.get("MOCK_FORWARD_FAIL") else 0)
elif "list" in args:
    if os.environ.get("MOCK_LIST_FAIL"):
        sys.exit(7)
    print("NAME\\tDISPLAY\\tSTATE\\tLISTEN")
    print(os.environ.get("MOCK_SESSIONS", ""), end="")
elif "ensure" in args or "display" in args:
    if "ensure" in args and os.environ.get("MOCK_ENSURE_FAIL"):
        sys.exit(7)
    selector = next((arg for arg in args if arg.startswith(":") and arg[1:].isdigit()), "1")
    print(os.environ.get("MOCK_DISPLAY", selector.lstrip(":")))
'''
VIEWER = '''#!/usr/bin/env python3
import json, os, sys
with open(os.environ["VNC_TEST_LOG"], "a") as log:
    log.write(json.dumps(["viewer", *sys.argv[1:]]) + "\\n")
'''


class MacTunnelTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="vnc-lan-test-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / ".vnc").mkdir()
        (self.root / ".vnc/passwd").touch()
        for name, source in [("ssh", SSH), ("viewer", VIEWER)]:
            path = self.root / name
            path.write_text(source)
            path.chmod(0o755)
        self.log = self.root / "calls"
        self.env = {key: value for key, value in os.environ.items()
                    if not key.startswith(("VNC_", "MOCK_"))}
        self.env.update(HOME=str(self.root), TMPDIR=str(self.root),
                        PATH=str(self.root) + os.pathsep + os.environ["PATH"],
                        VNC_VIEWER=str(self.root / "viewer"),
                        VNC_TEST_LOG=str(self.log), OSTYPE='darwin')

    def run_client(self, *args, program=None, **env):
        self.log.write_text("")
        result = subprocess.run(["bash", str(program or REPO / "vnc-lan"), *args],
                                env=self.env | env, capture_output=True,
                                text=True, timeout=5)
        self.calls = [json.loads(line) for line in self.log.read_text().splitlines()]
        self.forwards = [call for call in self.calls if "-fNT" in call]
        self.viewers = [call for call in self.calls if call[0] == "viewer"]
        return result

    def assert_connected(self, result, target, forward=None):
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(self.viewers), 1)
        self.assertEqual(self.viewers[0][-1], target)
        self.assertIn("-Scale", self.viewers[0])
        if forward is None:
            self.assertEqual(self.forwards, [])
        else:
            self.assertEqual(len(self.forwards), 1)
            self.assertEqual(self.forwards[0][self.forwards[0].index("-L") + 1], forward)
            self.assertIn("ExitOnForwardFailure=yes", self.forwards[0])

    def test_direct_lan_and_explicit_target(self):
        for args, target in [(('server',), 'server:1'),
                             (('user@server',), 'server:1'),
                             (('gateway', 'desktop', '1600x900'), 'desktop:1')]:
            with self.subTest(args=args):
                self.assert_connected(self.run_client(*args, VNC_TUNNEL='auto'), target)

    def test_auto_tunnels_for_ssh_routing(self):
        for env in [dict(MOCK_HOSTNAME="other-server"), dict(MOCK_PORT="2222"),
                    dict(MOCK_PROXYJUMP="gateway"),
                    dict(MOCK_PROXYCOMMAND="ssh gateway -W %h:%p")]:
            with self.subTest(env=env):
                self.assert_connected(self.run_client('server', VNC_TUNNEL='auto', **env),
                                      'localhost:11', '127.0.0.1:5911:127.0.0.1:5901')

    def test_override_modes_and_separate_target(self):
        self.assert_connected(self.run_client('server', VNC_TUNNEL='never', MOCK_PORT='2222'),
                              'server:1')
        self.assert_connected(self.run_client('gateway', 'desktop', VNC_TUNNEL='auto', MOCK_PROXYJUMP='jump'),
                              'desktop:1')
        self.assert_connected(self.run_client('gateway', 'desktop', VNC_TUNNEL='always'),
                              'localhost:11', '127.0.0.1:5911:desktop:5901')

    def test_named_sessions_keep_tunneling(self):
        for args in [('server', '--session', 'work'), ('--session', 'work', 'server')]:
            with self.subTest(args=args):
                self.assert_connected(self.run_client(*args, VNC_TUNNEL='never',
                                                      MOCK_DISPLAY='2', VNC_LOCAL_DISPLAY='21'),
                                      'localhost:21', '127.0.0.1:5921:127.0.0.1:5902')
                self.assertTrue(any('ensure' in call and 'work' in call for call in self.calls))

    def test_status_and_stop_do_not_require_viewer_or_password(self):
        (self.root / '.vnc/passwd').unlink()
        for args in [('server', '--status', 'work'), ('--status', 'work', 'server'),
                     ('server', '--stop', '--force', 'work'),
                     ('--stop', '--force', 'work', 'server')]:
            with self.subTest(args=args):
                result = self.run_client(*args, VNC_VIEWER='missing-viewer')
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertFalse(self.viewers)
                self.assertTrue(any('status' in call or 'stop' in call for call in self.calls))
                if '--stop' in args:
                    self.assertTrue(any('--force' in call for call in self.calls))
                    self.assertTrue(any('-O' in call and 'exit' in call for call in self.calls))

    def test_legacy_forward_identity_includes_target_and_local_port(self):
        paths = []
        for host, port in [('desktop', '11'), ('other-desktop', '11'), ('desktop', '12')]:
            result = self.run_client('gateway', host, VNC_TUNNEL='always', VNC_LOCAL_DISPLAY=port)
            self.assertEqual(result.returncode, 0, result.stderr)
            paths.append(self.forwards[0][self.forwards[0].index('-S') + 1])
        self.assertEqual(len(set(paths)), 3)

    def test_reuses_existing_tunnel(self):
        self.assert_connected(self.run_client('server', VNC_TUNNEL='always', MOCK_REUSE='1'),
                              'localhost:11')

    def test_control_socket_fits_macos_with_long_tmpdir(self):
        result = self.run_client('server', VNC_TUNNEL='always',
                                 TMPDIR=str(self.root / ('long-macos-tmpdir-' * 5)))
        self.assertEqual(result.returncode, 0, result.stderr)
        control_path = self.forwards[0][self.forwards[0].index('-S') + 1]
        # OpenSSH expands %C to 40 characters and adds a temporary suffix.
        self.assertLess(len(control_path.replace('%C', '0' * 40)) + 17, 104)

    def test_errors_do_not_launch_viewer(self):
        for env, message in [(dict(VNC_TUNNEL='invalid'), 'invalid VNC_TUNNEL'),
                             (dict(MOCK_CONFIG_FAIL='1', VNC_TUNNEL='auto'), 'could not read SSH configuration'),
                             (dict(VNC_TUNNEL='always', VNC_LOCAL_DISPLAY='0'), 'invalid VNC_LOCAL_DISPLAY'),
                             (dict(VNC_TUNNEL='always', MOCK_DISPLAY='60000'), 'VNC display is too large'),
                             (dict(VNC_TUNNEL='always', MOCK_FORWARD_FAIL='1'), ''),
                             (dict(MOCK_DISPLAY='no display'), 'invalid remote display'),
                             (dict(MOCK_ENSURE_FAIL='1'), '')]:
            with self.subTest(env=env):
                result = self.run_client('server', **env)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.viewers)
                self.assertIn(message, result.stderr)

    def alias(self, name):
        path = self.root / name
        path.symlink_to(REPO / 'vnc-lan')
        return path

    def test_host_alias_ensures_exact_display(self):
        alias = self.alias('vnc-srv16')
        for args, display in [((), '1'), ((':2',), '2'), (('3',), '3')]:
            with self.subTest(args=args):
                self.assert_connected(self.run_client(*args, program=alias),
                                      'localhost:' + str(int(display) + 10),
                                      f'127.0.0.1:{5910 + int(display)}:127.0.0.1:{5900 + int(display)}')
                self.assertTrue(any('ensure' in call and ':' + display in call
                                    and 'srv16' in call for call in self.calls))
                self.assertFalse(any('list' in call for call in self.calls))

    def test_user_alias_reuses_first_owned_running_desktop(self):
        alias = self.alias('vnc-srv16-umango')
        sessions = 'old\t:2\tstopped\tloopback\nlater\t:10\trunning\tloopback\n-\t:3\trunning\tunknown\n'
        self.assert_connected(self.run_client(program=alias, MOCK_SESSIONS=sessions),
                              'localhost:13', '127.0.0.1:5913:127.0.0.1:5903')
        self.assertTrue(any('list' in call and 'umango@srv16' in call for call in self.calls))
        self.assertTrue(any('ensure' in call and ':3' in call and 'umango@srv16' in call
                            for call in self.calls))

    def test_user_alias_allocates_only_when_no_owned_desktop(self):
        alias = self.alias('vnc-srv16-umango')
        self.assert_connected(self.run_client(program=alias, MOCK_DISPLAY='7'),
                              'localhost:17', '127.0.0.1:5917:127.0.0.1:5907')
        self.assertTrue(any('ensure' in call and 'self' in call for call in self.calls))
        self.assertFalse(any(':1' in call for call in self.calls))

    def test_user_alias_explicit_display_skips_discovery(self):
        self.assert_connected(self.run_client(':2', program=self.alias('vnc-srv16-umango')),
                              'localhost:12', '127.0.0.1:5912:127.0.0.1:5902')
        self.assertFalse(any('list' in call for call in self.calls))

    def test_direct_exact_display_and_restart_are_explicit(self):
        result = self.run_client('--direct', '--force-restart', 'server:2')
        self.assert_connected(result, 'server:2')
        self.assertTrue(any('ensure' in call and '--direct' in call and '--force-restart' in call
                            and ':2' in call for call in self.calls))
        result = self.run_client('--force-restart', 'server:2')
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.viewers)
        self.assertFalse(any('ensure' in call for call in self.calls))

    def test_list_needs_no_viewer_and_filters_stopped_desktops(self):
        (self.root / '.vnc/passwd').unlink()
        result = self.run_client('--list', 'server', VNC_VIEWER='missing-viewer',
                                 MOCK_SESSIONS='old\t:2\tstopped\tloopback\n-\t:3\trunning\tunknown\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn('stopped', result.stdout)
        self.assertIn(':3', result.stdout)
        self.assertFalse(self.viewers)

    def test_user_discovery_failure_does_not_create_desktop(self):
        result = self.run_client(program=self.alias('vnc-srv16-umango'), MOCK_LIST_FAIL='1')
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.viewers)
        self.assertFalse(any('ensure' in call for call in self.calls))


if __name__ == '__main__':
    unittest.main(verbosity=2)
