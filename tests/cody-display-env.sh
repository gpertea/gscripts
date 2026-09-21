#!/usr/bin/env bash
set -euo pipefail

repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
fixture="$(mktemp -d)"
trap 'rm -rf -- "$fixture"' EXIT

mkdir -p "$fixture/bin" "$fixture/home" "$fixture/runtime"
cat > "$fixture/bin/codex" <<'SH'
#!/usr/bin/env bash
if [[ "${1:-}" == "app-server" ]]; then
  printf 'daemon_arg=%s\n' "$@" >> "$CODY_TEST_LOG"
  exit 0
fi
{
  printf 'arg=%s\n' "$@"
  printf 'DISPLAY=%s\n' "${DISPLAY-}"
  printf 'XAUTHORITY=%s\n' "${XAUTHORITY-}"
  printf 'PIN=%s\n' "${LIN_BROWSER_USE_AUTH_DISPLAY-}"
  printf 'TERM=%s\n' "${TERM-}"
  printf 'COLORTERM=%s\n' "${COLORTERM-}"
  printf 'NO_COLOR=%s\n' "${NO_COLOR-}"
  printf 'FORCE_COLOR=%s\n' "${FORCE_COLOR-}"
  printf 'CLICOLOR_FORCE=%s\n' "${CLICOLOR_FORCE-}"
} >> "$CODY_TEST_LOG"
SH
chmod +x "$fixture/bin/codex"

cat > "$fixture/bin/bru" <<'SH'
#!/usr/bin/env bash
printf 'bru_arg=%s\n' "$@" >> "$CODY_TEST_LOG"
SH
chmod +x "$fixture/bin/bru"

export PATH="$fixture/bin:$PATH"
export HOME="$fixture/home"
export XDG_RUNTIME_DIR="$fixture/runtime"
export CODY_TEST_LOG="$fixture/calls.log"

: > "$CODY_TEST_LOG"
DISPLAY=:42.0 TERM=xterm-256color COLORTERM=rxvt-xpm NO_COLOR=1 \
  env -u XAUTHORITY -u FORCE_COLOR -u CLICOLOR_FORCE "$repo/cody" first
grep -Fxq 'daemon_arg=app-server' "$fixture/calls.log"
grep -Fxq 'daemon_arg=daemon' "$fixture/calls.log"
grep -Fxq 'daemon_arg=start' "$fixture/calls.log"
grep -Fxq 'arg=--yolo' "$fixture/calls.log"
grep -Fxq 'arg=--remote' "$fixture/calls.log"
grep -Fxq 'arg=unix://' "$fixture/calls.log"
grep -Fxq 'arg=-C' "$fixture/calls.log"
grep -Fxq "arg=$repo" "$fixture/calls.log"
grep -Fxq 'arg=shell_environment_policy.set.DISPLAY=":42.0"' "$fixture/calls.log"
grep -Fxq 'bru_arg=--select' "$fixture/calls.log"
grep -Fxq 'bru_arg=:42' "$fixture/calls.log"
grep -Fxq 'arg=first' "$fixture/calls.log"
grep -Fxq 'DISPLAY=:42.0' "$fixture/calls.log"
grep -Fxq "XAUTHORITY=$fixture/home/.Xauthority" "$fixture/calls.log"
grep -Fxq 'PIN=:42' "$fixture/calls.log"
grep -Fxq 'TERM=xterm-256color' "$fixture/calls.log"
grep -Fxq 'COLORTERM=rxvt-xpm' "$fixture/calls.log"
grep -Fxq 'NO_COLOR=1' "$fixture/calls.log"
grep -Fxq 'FORCE_COLOR=' "$fixture/calls.log"
grep -Fxq 'CLICOLOR_FORCE=' "$fixture/calls.log"
printf '%s\n' 'PASS: display pin and terminal colors are preserved'

: > "$CODY_TEST_LOG"
DISPLAY=:42.0 "$repo/cody" resume --last
grep -Fxq 'arg=--remote' "$fixture/calls.log"
grep -Fxq 'arg=resume' "$fixture/calls.log"
grep -Fxq 'arg=--last' "$fixture/calls.log"
if grep -Fxq 'arg=--yolo' "$fixture/calls.log"; then
  printf '%s\n' 'FAIL: remote resume included --yolo' >&2
  exit 1
fi
printf '%s\n' 'PASS: remote resume omits unsupported permission overrides'

: > "$CODY_TEST_LOG"
env -u DISPLAY -u LIN_BROWSER_USE_AUTH_DISPLAY "$repo/cody" terminal-only
grep -Fxq 'arg=--yolo' "$fixture/calls.log"
grep -Fxq 'arg=terminal-only' "$fixture/calls.log"
if grep -Fxq 'arg=--remote' "$fixture/calls.log"; then
  printf '%s\n' 'FAIL: display-less launch used app-server' >&2
  exit 1
fi
grep -Fxq 'PIN=' "$fixture/calls.log"
printf '%s\n' 'PASS: display-less launches remain local and unpinned'
