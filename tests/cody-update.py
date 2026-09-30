#!/usr/bin/env python3
"""Exercise update/restart consent in a PTY using isolated CLI/npm doubles."""

import errno
import os
from pathlib import Path
import pty
import select
import subprocess
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
        self.log = self.root / "calls"
        self.version = self.root / "version"
        self.version.write_text("0.159.2")
        (self.root / "server-status").write_text("running")
        (self.root / "server-version").write_text("0.155.1")
        self.environment = dict(os.environ)
        self.environment.update(
            PATH=f"{self.bin}:{os.environ['PATH']}",
            HOME=str(self.root),
            CODY_TEST_ROOT=str(self.root),
            CODY_TEST_LATEST="0.159.3",
        )
        self.environment.pop("DISPLAY", None)
        self.write_executable("pgrep", "#!/bin/bash\nexit 1\n")
        self.write_executable("codex", r'''#!/bin/bash
set -eu
printf 'codex %s\n' "$*" >> "$CODY_TEST_ROOT/calls"
if [[ "${1:-}" == --version ]]; then
  printf 'codex-cli %s\n' "$(cat "$CODY_TEST_ROOT/version")"
elif [[ "$*" == 'app-server daemon version' ]]; then
  printf '{"status":"%s","appServerVersion":"%s","managedCodexVersion":"%s"}\n' \
    "$(cat "$CODY_TEST_ROOT/server-status")" \
    "$(cat "$CODY_TEST_ROOT/server-version")" \
    "$(cat "$CODY_TEST_ROOT/server-version")"
elif [[ "$*" == 'app-server daemon stop' ]]; then
  printf stopped > "$CODY_TEST_ROOT/server-status"
elif [[ "$*" == 'app-server daemon update --from-cli --yes' ]]; then
  if [[ "${CODY_TEST_DAEMON_FAIL:-0}" == 1 ]]; then exit 8; fi
  cp "$CODY_TEST_ROOT/version" "$CODY_TEST_ROOT/server-version"
elif [[ "$*" == 'app-server daemon start' ]]; then
  printf running > "$CODY_TEST_ROOT/server-status"
else
  printf 'LAUNCH %s\n' "$*"
fi
''')
        self.write_executable("npm", r'''#!/bin/bash
set -eu
printf 'npm %s\n' "$*" >> "$CODY_TEST_ROOT/calls"
case "$1" in
  view)
    if [[ "${CODY_TEST_CHECK_FAIL:-0}" == 1 ]]; then exit 6; fi
    printf '%s\n' "$CODY_TEST_LATEST"
    ;;
  install)
    if [[ "${CODY_TEST_UPDATE_FAIL:-0}" == 1 ]]; then exit 7; fi
    printf '%s' "$CODY_TEST_LATEST" > "$CODY_TEST_ROOT/version"
    ;;
esac
''')

    def write_executable(self, name, source):
        path = self.bin / name
        path.write_text(source)
        path.chmod(0o755)

    def run_cody(self, arguments=(), answers=b"", interactive=True):
        if not interactive:
            result = subprocess.run(
                [str(REPO / "cody"), *arguments], env=self.environment,
                input=answers, capture_output=True, timeout=10,
            )
            return result.returncode, (result.stdout + result.stderr).decode()

        ## give Bash a real terminal without invoking installed Codex or npm
        master, slave = pty.openpty()
        process = subprocess.Popen(
            [str(REPO / "cody"), *arguments], env=self.environment,
            stdin=slave, stdout=slave, stderr=slave,
        )
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

    def calls(self):
        return self.log.read_text() if self.log.exists() else ""

    def test_accepted_update_and_restart(self):
        code, output = self.run_cody(answers=b"y\ny\n")
        self.assertEqual(code, 0, output)
        calls = self.calls()
        sequence = [
            "npm install --global @openai/codex@latest",
            "codex app-server daemon stop",
            "codex app-server daemon update --from-cli --yes",
            "codex app-server daemon start",
            "codex --yolo -c check_for_update_on_startup=false",
        ]
        positions = [calls.index(command) for command in sequence]
        self.assertEqual(positions, sorted(positions))
        self.assertIn("Installed CLI: 0.159.3", output)
        self.assertIn("Running app-server: 0.155.1", output)
        self.assertIn("Running app-server: 0.159.3", output)

    def test_restart_declined(self):
        code, output = self.run_cody(answers=b"y\nn\n")
        self.assertEqual(code, 0, output)
        self.assertEqual(self.version.read_text(), "0.159.3")
        self.assertNotIn("codex app-server daemon stop", self.calls())
        self.assertNotIn("codex app-server daemon update", self.calls())
        self.assertIn("App-server left running", output)

    def test_update_declined(self):
        code, output = self.run_cody(answers=b"n\n")
        self.assertEqual(code, 0, output)
        self.assertNotIn("npm install", self.calls())
        self.assertNotIn("app-server daemon", self.calls())
        self.assertIn("check_for_update_on_startup=false", self.calls())

    def test_noupdate_disables_both_checks(self):
        code, output = self.run_cody(("--noupdate", "--docs", "resume"))
        self.assertEqual(code, 0, output)
        self.assertNotIn("npm ", self.calls())
        self.assertNotIn("app-server daemon", self.calls())
        self.assertIn("check_for_update_on_startup=false", output)
        self.assertIn("mcp_servers.openaiDeveloperDocs.enabled=true", output)
        self.assertIn("resume", output)

    def test_no_server_is_started(self):
        (self.root / "server-status").write_text("stopped")
        code, output = self.run_cody(answers=b"y\n")
        self.assertEqual(code, 0, output)
        self.assertIn("No running app-server to restart", output)
        self.assertNotIn("codex app-server daemon start", self.calls())

    def test_failed_cli_update_leaves_server_alone(self):
        self.environment["CODY_TEST_UPDATE_FAIL"] = "1"
        code, output = self.run_cody(answers=b"y\n")
        self.assertEqual(code, 7, output)
        self.assertNotIn("app-server daemon", self.calls())
        self.assertNotIn("LAUNCH", output)

    def test_failed_daemon_update_does_not_start_it(self):
        self.environment["CODY_TEST_DAEMON_FAIL"] = "1"
        code, output = self.run_cody(answers=b"y\ny\n")
        self.assertEqual(code, 8, output)
        self.assertNotIn("codex app-server daemon start", self.calls())

    def test_offline_check_still_launches(self):
        self.environment["CODY_TEST_CHECK_FAIL"] = "1"
        code, output = self.run_cody()
        self.assertEqual(code, 0, output)
        self.assertIn("update check unavailable", output)
        self.assertIn("LAUNCH", output)
        self.assertNotIn("npm install", self.calls())

    def test_newer_local_version_is_not_downgraded(self):
        self.version.write_text("0.160.0")
        code, output = self.run_cody()
        self.assertEqual(code, 0, output)
        self.assertNotIn("npm install", self.calls())
        self.assertIn("LAUNCH", output)

    def test_noninteractive_launch_does_no_network_work(self):
        code, output = self.run_cody(interactive=False)
        self.assertEqual(code, 0, output)
        self.assertNotIn("npm ", self.calls())
        self.assertIn("LAUNCH --yolo", output)

    def test_noninteractive_explicit_update_reports_without_restart(self):
        code, output = self.run_cody(("--update",), interactive=False)
        self.assertEqual(code, 0, output)
        self.assertIn("Installed CLI: 0.159.3", output)
        self.assertNotIn("codex app-server daemon stop", self.calls())
        self.assertNotIn("LAUNCH", output)

    def test_conflicting_options_make_no_changes(self):
        code, output = self.run_cody(("--noupdate", "--update"))
        self.assertEqual(code, 2, output)
        self.assertNotIn("npm ", self.calls())


if __name__ == "__main__":
    unittest.main()
