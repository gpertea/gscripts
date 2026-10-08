#!/usr/bin/env bash
set -euo pipefail

repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
fixture="$(mktemp -d)"
trap 'rm -rf -- "$fixture"' EXIT

mkdir -p "$fixture/bin" "$fixture/home" "$fixture/runtime"
## isolate routing checks from real X sessions and app-servers
cp "$repo/cody" "$fixture/cody"
cat > "$fixture/cody-display" <<'SH'
#!/usr/bin/env bash
printf 'manager=' > "$CODY_TEST_MANAGER_LOG"
printf '<%s>' "$@" >> "$CODY_TEST_MANAGER_LOG"
printf '\n' >> "$CODY_TEST_MANAGER_LOG"
while [[ "$1" != -- ]]; do shift; done
shift
exec codex "$@"
SH
chmod +x "$fixture/cody-display"
cat > "$fixture/bin/codex" <<'SH'
#!/usr/bin/env bash
{
  printf 'args='
  printf '<%s>' "$@"
  printf '\n'
  printf 'DISPLAY=%s\n' "${DISPLAY-}"
  printf 'XAUTHORITY=%s\n' "${XAUTHORITY-}"
  printf 'PIN=%s\n' "${LIN_BROWSER_USE_AUTH_DISPLAY-}"
  printf 'XDG_RUNTIME_DIR=%s\n' "${XDG_RUNTIME_DIR-}"
  printf 'DBUS_SESSION_BUS_ADDRESS=%s\n' "${DBUS_SESSION_BUS_ADDRESS-}"
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
export CODY_TEST_MANAGER_LOG="$fixture/manager.log"

DISPLAY=:42.0 TERM=xterm-256color COLORTERM=rxvt-xpm NO_COLOR=1 \
  env -u XAUTHORITY -u FORCE_COLOR -u CLICOLOR_FORCE "$fixture/cody" first
grep -Fxq 'manager=<local><--><--yolo><first>' "$fixture/manager.log"
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

## startup resume filtering must not change plain launches or explicit choices
DISPLAY=:42 "$fixture/cody" resume
grep -Fxq "manager=<local><--restore-yolo><--><--cd><$PWD><resume>" "$fixture/manager.log"
DISPLAY=:42 "$fixture/cody" --attach :42 resume
grep -Fxq "manager=<attach><--display><:42><--restore-yolo><--><--cd><$PWD><resume>" "$fixture/manager.log"
for override in '--all' '--cd=/tmp/another-project' '-C/tmp/another-project'; do
  DISPLAY=:42 "$fixture/cody" resume "$override"
  grep -Fxq "args=<resume><$override>" "$fixture/calls.log"
done
DISPLAY=:42 "$fixture/cody" resume --cd /tmp/another-project
grep -Fxq 'args=<resume><--cd></tmp/another-project>' "$fixture/calls.log"
DISPLAY=:42 "$fixture/cody" -m resume --no-alt-screen resume
grep -Fxq "args=<--cd><$PWD><-m><resume><--no-alt-screen><resume>" "$fixture/calls.log"
DISPLAY=:42 "$fixture/cody" -C /tmp/another-project resume
grep -Fxq 'args=<-C></tmp/another-project><resume>' "$fixture/calls.log"
DISPLAY=:42 "$fixture/cody" fork
grep -Fxq 'args=<fork>' "$fixture/calls.log"
DISPLAY=:42 "$fixture/cody" -m resume new-prompt
grep -Fxq 'args=<--yolo><-m><resume><new-prompt>' "$fixture/calls.log"
printf '%s\n' 'PASS: resume scope defaults to cwd and preserves explicit overrides'

## the terminal-only wrapper must work even with no display manager installed
rm "$fixture/cody-display" "$fixture/manager.log"
for display_state in unset empty; do
  display_env=(-u DISPLAY)
  [[ "$display_state" != empty ]] || display_env=(DISPLAY=)
  env -u LIN_BROWSER_USE_AUTH_DISPLAY -u XAUTHORITY -u XDG_RUNTIME_DIR \
    -u DBUS_SESSION_BUS_ADDRESS "${display_env[@]}" "$fixture/cody" terminal-only
  grep -Fxq 'args=<--yolo><terminal-only>' "$fixture/calls.log"
  for key in DISPLAY PIN XAUTHORITY XDG_RUNTIME_DIR DBUS_SESSION_BUS_ADDRESS; do
    grep -Fxq "$key=" "$fixture/calls.log"
  done
done
[[ ! -e "$fixture/manager.log" ]]
printf '%s\n' 'PASS: unset and empty DISPLAY launch directly without desktop setup or manager'
