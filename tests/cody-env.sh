#!/usr/bin/env bash
set -euo pipefail

repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
fixture="$(mktemp -d)"
trap 'rm -rf -- "$fixture"' EXIT
mkdir -p "$fixture/bin" "$fixture/home" "$fixture/runtime"
cp "$repo/cody" "$fixture/cody"
## any invocation of the retired helper is a regression
cat > "$fixture/cody-display" <<'STUB'
#!/bin/bash
exit 99
STUB
chmod +x "$fixture/cody-display"
cat > "$fixture/bin/codex" <<'STUB'
#!/bin/bash
{
  printf 'args='; printf '<%s>' "$@"; printf '\n'
  printf 'cwd=%s\n' "$PWD"
  for key in DISPLAY XAUTHORITY LIN_BROWSER_USE_AUTH_DISPLAY XDG_RUNTIME_DIR DBUS_SESSION_BUS_ADDRESS TERM COLORTERM NO_COLOR FORCE_COLOR CLICOLOR_FORCE CODEX_HOME; do
    printf '%s=%s\n' "$key" "${!key-}"
  done
} > "$CODY_TEST_LOG"
exit "${CODY_TEST_EXIT:-0}"
STUB
chmod +x "$fixture/bin/codex"
export PATH="$fixture/bin:$PATH" HOME="$fixture/home" HOSTNAME=fixture
export XDG_RUNTIME_DIR="$fixture/runtime" CODY_TEST_LOG="$fixture/calls.log"

DISPLAY=:42.0 TERM=xterm-256color COLORTERM=rxvt-xpm NO_COLOR=1 \
 env -u LIN_BROWSER_USE_AUTH_DISPLAY -u XAUTHORITY -u DBUS_SESSION_BUS_ADDRESS -u FORCE_COLOR -u CLICOLOR_FORCE "$fixture/cody" first
grep -Fxq 'args=<--yolo><--no-daemon><first>' "$CODY_TEST_LOG"
grep -Fxq 'DISPLAY=:42.0' "$CODY_TEST_LOG"
grep -Fxq "XAUTHORITY=$fixture/home/.Xauthority" "$CODY_TEST_LOG"
grep -Fxq 'LIN_BROWSER_USE_AUTH_DISPLAY=:42' "$CODY_TEST_LOG"
grep -Fxq "DBUS_SESSION_BUS_ADDRESS=unix:path=$fixture/runtime/bus" "$CODY_TEST_LOG"
for field in TERM=xterm-256color COLORTERM=rxvt-xpm NO_COLOR=1 FORCE_COLOR= CLICOLOR_FORCE=; do
  grep -Fxq "$field" "$CODY_TEST_LOG"
done
## explicit browser/session environment survives rather than being replaced
DISPLAY=:12 LIN_BROWSER_USE_AUTH_DISPLAY=:88 XAUTHORITY=/fictional/auth DBUS_SESSION_BUS_ADDRESS=unix:path=/fictional/bus CODEX_HOME=/fictional/codex \
 "$fixture/cody" resume --last
for field in LIN_BROWSER_USE_AUTH_DISPLAY=:88 XAUTHORITY=/fictional/auth DBUS_SESSION_BUS_ADDRESS=unix:path=/fictional/bus CODEX_HOME=/fictional/codex; do
  grep -Fxq "$field" "$CODY_TEST_LOG"
done
for display in '' localhost:10.0 :9; do
 env -u LIN_BROWSER_USE_AUTH_DISPLAY DISPLAY="$display" "$fixture/cody" resume --all
 grep -Fxq 'args=<--yolo><--no-daemon><resume><--all>' "$CODY_TEST_LOG"
done
env -u DISPLAY -u LIN_BROWSER_USE_AUTH_DISPLAY -u XAUTHORITY -u XDG_RUNTIME_DIR -u DBUS_SESSION_BUS_ADDRESS \
 "$fixture/cody" terminal-only
for key in DISPLAY LIN_BROWSER_USE_AUTH_DISPLAY XAUTHORITY XDG_RUNTIME_DIR DBUS_SESSION_BUS_ADDRESS; do
  grep -Fxq "$key=" "$CODY_TEST_LOG"
done
## native cwd/resume scoping needs no injected --cd or permission proxy
DISPLAY=:42 "$fixture/cody" -m resume -C '/fictional/path with spaces' resume --last
grep -Fxq 'args=<--yolo><--no-daemon><-m><resume><-C></fictional/path with spaces><resume><--last>' "$CODY_TEST_LOG"
grep -Fxq "cwd=$PWD" "$CODY_TEST_LOG"
DISPLAY=:42 "$fixture/cody" fork
grep -Fxq 'args=<--yolo><--no-daemon><fork>' "$CODY_TEST_LOG"
"$fixture/cody" --no-daemon resume
grep -Fxq 'args=<--yolo><--no-daemon><resume>' "$CODY_TEST_LOG"
"$fixture/cody" -- '--remote'
grep -Fxq 'args=<--yolo><--no-daemon><--><--remote>' "$CODY_TEST_LOG"
for flag in --attach --server-status --server-stop --server-endpoint --remote; do
 if "$fixture/cody" resume "$flag" > "$fixture/output" 2>&1; then exit 1; fi
 grep -q 'remote/display-server options are retired' "$fixture/output"
done
"$fixture/cody" mcp list
grep -Fxq 'args=<--yolo><mcp><list>' "$CODY_TEST_LOG"
if CODY_TEST_EXIT=23 "$fixture/cody"; then exit 1; else [[ $? == 23 ]]; fi
printf '%s\n' 'PASS: direct embedded routing, display/environment preservation, cwd/arguments, retired flags, exit status'
