#!/usr/bin/env python3
"""Check Linux viewer lifetime through real host aliases, without network/GUI."""

import json
import os
from pathlib import Path
import signal
import subprocess
import tempfile
import time
import unittest


REPO = Path(__file__).resolve().parents[1]
SSH = '''#!/usr/bin/env python3
import os, sys
args = sys.argv[1:]
if "-G" in args:
    print("hostname", args[-1])
    if "-S" in args:
        print("controlpath", args[args.index("-S") + 1].replace("%C", "fixture"))
elif "-O" in args:
    sys.exit(0 if args[args.index("-O") + 1] == "check" else 1)
elif "--help" in args:
    print("ensure [--direct|--tunnel]; stop [--force] NAME|:DISPLAY")
elif "ensure" in args:
    if os.environ.get("MOCK_ENSURE_FAIL"):
        print("fixture session failure", file=sys.stderr)
        sys.exit(7)
    print("1")
elif "list" in args:
    print("NAME DISPLAY STATE LISTEN")
    print("desktop :1 running loopback")
elif "display" in args:
    print("1")
elif "status" in args:
    print("State: running")
elif "stop" in args:
    print("Stopped")
else:
    sys.exit(2)
'''
VIEWER = '''#!/usr/bin/env python3
import json, os, sys, time
from pathlib import Path
record = Path(os.environ["VIEWER_RECORD"])
temporary = record.with_suffix(".tmp")
temporary.write_text(json.dumps({
    "pid": os.getpid(), "args": sys.argv[1:],
    "stdin": os.readlink("/proc/self/fd/0"),
    "stdout": os.readlink("/proc/self/fd/1"),
    "stderr": os.readlink("/proc/self/fd/2"),
}))
temporary.replace(record)
print("fixture viewer started", flush=True)
time.sleep(60)
'''


class ViewerBackgroundTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="vnc-background-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for directory in ("bin", "runtime", ".vnc"):
            (self.root / directory).mkdir()
        (self.root / ".vnc/passwd").touch()
        self.record = self.root / "viewer.json"
        self.addCleanup(self.stop_viewer)
        for name, source in (("ssh", SSH), ("viewer", VIEWER)):
            path = self.root / "bin" / name
            path.write_text(source)
            path.chmod(0o755)

        # fake only SSH PID discovery; exercise actual parsing and launch paths
        source = (REPO / "vnc-session-view").read_text()
        marker = 'if [[ "$action" == "list" ]]; then'
        self.assertEqual(source.count(marker), 1)
        source = source.replace(marker,
                                'remember_tunnel() { local_port=45901; }\n\n' + marker)
        for name, content in (("vnc-session-view", source),
                              ("vnc-host", (REPO / "vnc-host").read_text())):
            path = self.root / "bin" / name
            path.write_text(content)
            path.chmod(0o755)
        for host in ("srv16", "srv05", "futurehost"):
            (self.root / "bin" / ("vnc-" + host)).symlink_to("vnc-host")
        self.env = {key: value for key, value in os.environ.items()
                    if not key.startswith(("VNC_", "MOCK_"))}
        self.env.update(HOME=str(self.root), XDG_STATE_HOME=str(self.root / "state"),
                        XDG_RUNTIME_DIR=str(self.root / "runtime"),
                        PATH=str(self.root / "bin") + os.pathsep + os.environ["PATH"],
                        VNC_VIEWER=str(self.root / "bin/viewer"),
                        VIEWER_RECORD=str(self.record), OSTYPE="linux-gnu")

    def stop_viewer(self):
        if self.record.exists():
            pid = json.loads(self.record.read_text())["pid"]
            try:
                os.kill(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            self.record.unlink()

    def run_alias(self, host, *args, **env):
        process = subprocess.Popen([str(self.root / "bin" / ("vnc-" + host)), *args],
                                   env=self.env | env, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True,
                                   start_new_session=True)
        try:
            output, error = process.communicate(timeout=2)
        except subprocess.TimeoutExpired:
            self.stop_viewer()
            process.kill()
            process.communicate()
            self.fail("host command blocked while its viewer was open")
        return process.returncode, output, error

    def test_aliases_return_with_viewer_alive_for_both_transports(self):
        for host in ("srv16", "srv05", "futurehost"):
            for options in ((), ("--direct",)):
                with self.subTest(host=host, options=options):
                    code, output, error = self.run_alias(host, *options)
                    self.assertEqual(code, 0, error)
                    deadline = time.monotonic() + 2
                    while not self.record.exists() and time.monotonic() < deadline:
                        time.sleep(0.01)
                    self.assertTrue(self.record.exists(), output + error)
                    viewer = json.loads(self.record.read_text())
                    expected = host + ":1" if options else "localhost::45901"
                    self.assertEqual(viewer["args"],
                                     [f"passwd={self.root}/.vnc/passwd", expected])
                    self.assertEqual(viewer["stdin"], "/dev/null")
                    log = str(self.root / "state/vnc-session/viewer.log")
                    self.assertEqual(viewer["stdout"], log)
                    self.assertEqual(viewer["stderr"], log)
                    os.kill(viewer["pid"], signal.SIGHUP)
                    time.sleep(0.05)
                    os.kill(viewer["pid"], 0)
                    self.assertNotEqual(Path(f'/proc/{viewer["pid"]}/stat').read_text()
                                        .split(") ", 1)[1].split()[0], "Z")
                    self.stop_viewer()

    def test_setup_failure_remains_visible_and_launches_no_viewer(self):
        code, output, error = self.run_alias("srv05", MOCK_ENSURE_FAIL="1")
        self.assertNotEqual(code, 0)
        self.assertIn("fixture session failure", error)
        self.assertFalse(self.record.exists(), output)

    def test_management_actions_remain_foreground_and_launch_no_viewer(self):
        for options in (("--list",), ("--status",), ("--stop", "--force")):
            with self.subTest(options=options):
                code, output, error = self.run_alias("srv16", *options)
                self.assertEqual(code, 0, error)
                self.assertFalse(self.record.exists())
                if options == ("--list",):
                    self.assertIn("desktop :1 running loopback", output)
                elif options == ("--status",):
                    self.assertIn("State: running", output)


if __name__ == "__main__":
    unittest.main()
