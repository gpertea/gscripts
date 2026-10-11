#!/usr/bin/env bash
set -euo pipefail

repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
fixture="$(mktemp -d)"
trap 'rm -rf -- "$fixture"' EXIT
mkdir -p "$fixture/bin"
export VNC_CLI_FIXTURE="$fixture"
export PATH="$fixture/bin:$PATH"
unset VNC_HOST_ALIAS VNC_HOST_DEFAULT_TRANSPORT

## mock only SSH so real option parsing and platform routing are exercised
cat > "$fixture/bin/ssh" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >> "$VNC_CLI_FIXTURE/calls"
[[ "$*" == '-n srv16 "$HOME/bin/vnc-session" list' ]] || exit 2
[[ ! -e "$VNC_CLI_FIXTURE/fail" ]] || exit 7
cat "$VNC_CLI_FIXTURE/sessions"
MOCK
chmod +x "$fixture/bin/ssh"
printf 'NAME\tDISPLAY\tSTATE\tLISTEN\nbrowser\t:2\trunning\tloopback\nold\t:3\tstopped\tlan\n-\t:4\trunning\tunknown\n' > "$fixture/sessions"
printf 'NAME\tDISPLAY\tSTATE\tLISTEN\nbrowser\t:2\trunning\tloopback\n-\t:4\trunning\tunknown\n' > "$fixture/expected"
ln -s "$repo/vnc-host" "$fixture/bin/vnc-srv16"

for program in vnc-host vnc-session-view vnc-win-lan; do
  "$repo/$program" -h > "$fixture/short-help"
  "$repo/$program" --help > "$fixture/long-help"
  cmp "$fixture/short-help" "$fixture/long-help"
done
printf '%s\n' 'PASS: -h matches --help for wrapper and both helpers'

for client_os in linux-gnu msys; do
  for option in -l --list; do
    : > "$fixture/calls"
    OSTYPE="$client_os" "$repo/vnc-host" "$option" srv16 > "$fixture/output"
    cmp "$fixture/expected" "$fixture/output"
    [[ "$(wc -l < "$fixture/calls")" == 1 ]]
    OSTYPE="$client_os" "$fixture/bin/vnc-srv16" "$option" > "$fixture/output"
    cmp "$fixture/expected" "$fixture/output"
  done
done
printf '%s\n' 'PASS: Linux and Windows list active managed/unmanaged sessions via host and alias'

printf 'NAME\tDISPLAY\tSTATE\tLISTEN\nold\t:3\tstopped\tlan\n' > "$fixture/sessions"
printf 'NAME\tDISPLAY\tSTATE\tLISTEN\n' > "$fixture/expected"
"$repo/vnc-host" --list srv16 > "$fixture/output"
cmp "$fixture/expected" "$fixture/output"
printf '%s\n' 'PASS: no active sessions yields only the header'

touch "$fixture/fail"
for client_os in linux-gnu msys; do
  if OSTYPE="$client_os" "$repo/vnc-host" --list srv16 > "$fixture/output"; then
    echo 'FAIL: SSH failure was hidden' >&2
    exit 1
  fi
done
printf '%s\n' 'PASS: remote listing failures propagate on both platforms'
