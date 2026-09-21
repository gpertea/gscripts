#!/usr/bin/env bash
set -euo pipefail

repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
fixture="$(mktemp -d)"
server_pid=""
cleanup() {
  [[ -z "$server_pid" ]] || kill "$server_pid" 2>/dev/null || true
  rm -rf -- "$fixture"
}
trap cleanup EXIT

mkdir -p "$fixture/bin" "$fixture/home" "$fixture/runtime"
cat > "$fixture/bin/codex" <<'PY'
#!/usr/bin/env python3
import os
import signal
import socket
import sys
import time

log = os.environ["CODY_TEST_LOG"]
args = sys.argv[1:]
if args and args[0] == "app-server":
    endpoint = args[args.index("--listen") + 1]
    path = endpoint.removeprefix("unix://")
    with open(log, "a", encoding="ascii") as handle:
        handle.write(f"server DISPLAY={os.environ.get('DISPLAY')} PIN={os.environ.get('LIN_BROWSER_USE_AUTH_DISPLAY')}\n")
    listener = socket.socket(socket.AF_UNIX)
    listener.bind(path)
    listener.listen()
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(0))
    while True:
        time.sleep(1)

with open(log, "a", encoding="ascii") as handle:
    handle.write("client " + " ".join(args) + "\n")
PY
chmod +x "$fixture/bin/codex"

export PATH="$fixture/bin:$PATH"
export HOME="$fixture/home"
export XDG_RUNTIME_DIR="$fixture/runtime"
export CODY_TEST_LOG="$fixture/calls.log"

DISPLAY=:42.0 XAUTHORITY="$fixture/xauth" "$repo/cody" first
server_pid="$(cat "$fixture/runtime/cody/display-42/app-server.pid")"
grep -Fxq 'server DISPLAY=:42.0 PIN=:42' "$fixture/calls.log"
grep -Eq '^client --yolo --remote unix://.*/display-42/app-server.sock first$' "$fixture/calls.log"

DISPLAY=:42.0 XAUTHORITY="$fixture/xauth" "$repo/cody" second
[[ "$(grep -c '^server ' "$fixture/calls.log")" == 1 ]]
grep -Eq '^client --yolo --remote unix://.*/display-42/app-server.sock second$' "$fixture/calls.log"
printf '%s\n' 'PASS: numeric displays reuse one pinned app-server'

env -u DISPLAY "$repo/cody" terminal-only
grep -Fxq 'client --yolo terminal-only' "$fixture/calls.log"
printf '%s\n' 'PASS: launches without DISPLAY retain direct behavior'
