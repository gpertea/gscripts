#!/usr/bin/env bash
set -euo pipefail

repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
fixture="$(mktemp -d)"
trap 'rm -rf -- "$fixture"' EXIT

mkdir -p "$fixture/bin" "$fixture/home" "$fixture/runtime"
cat > "$fixture/bin/codex" <<'SH'
#!/usr/bin/env bash
{
  printf 'args='
  printf '<%s>' "$@"
  printf '\n'
  printf 'DISPLAY=%s\n' "${DISPLAY-}"
  printf 'XAUTHORITY=%s\n' "${XAUTHORITY-}"
  printf 'PIN=%s\n' "${LIN_BROWSER_USE_AUTH_DISPLAY-}"
  printf 'TERM=%s\n' "${TERM-}"
  printf 'COLORTERM=%s\n' "${COLORTERM-}"
  printf 'NO_COLOR=%s\n' "${NO_COLOR-}"
  printf 'FORCE_COLOR=%s\n' "${FORCE_COLOR-}"
  printf 'CLICOLOR_FORCE=%s\n' "${CLICOLOR_FORCE-}"
} > "$CODY_TEST_LOG"
SH
chmod +x "$fixture/bin/codex"

export PATH="$fixture/bin:$PATH"
export HOME="$fixture/home"
export XDG_RUNTIME_DIR="$fixture/runtime"
export CODY_TEST_LOG="$fixture/calls.log"

DISPLAY=:42.0 TERM=xterm-256color COLORTERM=rxvt-xpm NO_COLOR=1 \
  env -u XAUTHORITY -u FORCE_COLOR -u CLICOLOR_FORCE "$repo/cody" first
grep -Fxq 'args=<--yolo><first>' "$fixture/calls.log"
grep -Fxq 'DISPLAY=:42.0' "$fixture/calls.log"
grep -Fxq "XAUTHORITY=$fixture/home/.Xauthority" "$fixture/calls.log"
grep -Fxq 'PIN=:42' "$fixture/calls.log"
grep -Fxq 'TERM=xterm-256color' "$fixture/calls.log"
grep -Fxq 'COLORTERM=rxvt-xpm' "$fixture/calls.log"
grep -Fxq 'NO_COLOR=1' "$fixture/calls.log"
grep -Fxq 'FORCE_COLOR=' "$fixture/calls.log"
grep -Fxq 'CLICOLOR_FORCE=' "$fixture/calls.log"
printf '%s\n' 'PASS: display pin and terminal colors are preserved'

env -u DISPLAY -u LIN_BROWSER_USE_AUTH_DISPLAY "$repo/cody" terminal-only
grep -Fxq 'args=<--yolo><terminal-only>' "$fixture/calls.log"
grep -Fxq 'PIN=' "$fixture/calls.log"
printf '%s\n' 'PASS: display-less launches remain local and unpinned'
