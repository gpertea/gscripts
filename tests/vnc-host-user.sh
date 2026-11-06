#!/usr/bin/env bash
set -euo pipefail

repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
fixture="$(mktemp -d)"
trap 'rm -rf -- "$fixture"' EXIT
export VNC_USER_FIXTURE="$fixture"
export PATH="$fixture:$PATH"
unset VNC_HOST_ALIAS VNC_HOST_DEFAULT_TRANSPORT

cat > "$fixture/ssh" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
if [[ "$1" == -G ]]; then
  printf 'hostname %s\nport 22\n' "${2#*@}"
  exit
fi
printf '%s\n' "$*" >> "$VNC_USER_FIXTURE/ssh-calls"
[[ "$*" == '-n '*@*' "$HOME/bin/vnc-session" list' ]] || exit 2
[[ ! -e "$VNC_USER_FIXTURE/fail" ]] || exit 7
cat "$VNC_USER_FIXTURE/sessions"
MOCK
chmod +x "$fixture/ssh"
printf '#!/bin/sh\nexit 1\n' > "$fixture/nc"
chmod +x "$fixture/nc"
printf 'NAME\tDISPLAY\tSTATE\tLISTEN\n' > "$fixture/sessions"

## use adjacent viewer probes to inspect every argument passed by the real wrapper
cp "${VNC_HOST_TEST_SCRIPT:-$repo/vnc-host}" "$fixture/vnc-host"
cp "$repo/vnc-mac-route" "$fixture/vnc-mac-route"
chmod +x "$fixture/vnc-host"
for helper in vnc-session-view vnc-win-lan vnc-lan; do
  cat > "$fixture/$helper" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "${0##*/}" "$@" > "$VNC_USER_FIXTURE/actual"
MOCK
  chmod +x "$fixture/$helper"
done
for alias in vnc-srv16 vnc-srv16-umango vnc-lab-server-umango vnc-srv16-; do
  ln -s vnc-host "$fixture/$alias"
done

expect_route() {
  local client_os="$1" expected="$2"
  shift 2
  OSTYPE="$client_os" "$@"
  printf '%s\n' "$expected" > "$fixture/expected"
  diff -u "$fixture/expected" "$fixture/actual"
}

for client_os in linux-gnu msys cygwin darwin; do
  helper=vnc-session-view
  host_transport=--tunnel
  user_transport=--tunnel
  if [[ "$client_os" == darwin ]]; then
    helper=vnc-lan
    host_transport=--auto
    user_transport=--auto
  elif [[ "$client_os" != linux-gnu ]]; then
    helper=vnc-win-lan
    host_transport=--direct
  fi
  expect_route "$client_os" "$helper"$'\n--status\numango@srv16:self' \
    "$fixture/vnc-srv16-umango" --status
  expect_route "$client_os" "$helper"$'\n'"$user_transport"$'\numango@srv16:self' \
    "$fixture/vnc-srv16-umango"
  expect_route "$client_os" "$helper"$'\n'"$user_transport"$'\numango@lab-server:self' \
    "$fixture/vnc-lab-server-umango"
  expect_route "$client_os" "$helper"$'\n'"$user_transport"$'\numango@srv16:3' \
    "$fixture/vnc-srv16-umango" :3
  expect_route "$client_os" "$helper"$'\n'"$user_transport"$'\numango@srv16:browser' \
    "$fixture/vnc-srv16-umango" browser
  expect_route "$client_os" "$helper"$'\n--list\numango@srv16:self' \
    "$fixture/vnc-srv16-umango" --list
  expect_route "$client_os" "$helper"$'\n--stop\n--force\numango@srv16:self' \
    "$fixture/vnc-srv16-umango" --stop --force
  expect_route "$client_os" "$helper"$'\n--direct\n--force-restart\numango@srv16:self' \
    "$fixture/vnc-srv16-umango" --direct --force-restart
  expect_route "$client_os" "$helper"$'\n--direct\numango@srv16:self' \
    env VNC_HOST_DEFAULT_TRANSPORT=--direct "$fixture/vnc-srv16-umango"
  expect_route "$client_os" "$helper"$'\n--tunnel\numango@srv16:self' \
    env VNC_HOST_DEFAULT_TRANSPORT=--direct "$fixture/vnc-srv16-umango" --tunnel
  expect_route "$client_os" "$helper"$'\n'"$user_transport"$'\numango@srv16:self' \
    "$fixture/vnc-host" umango@srv16
  expect_route "$client_os" "$helper"$'\n'"$user_transport"$'\numango@srv16:3' \
    "$fixture/vnc-host" umango@srv16:3
  expect_route "$client_os" "$helper"$'\n'"$user_transport"$'\numango@srv16:self' \
    env VNC_HOST_ALIAS=umango@srv16 "$fixture/vnc-host"
  expect_route "$client_os" "$helper"$'\n'"$host_transport"$'\nsrv16:1' \
    "$fixture/vnc-srv16"
  expect_route "$client_os" "$helper"$'\n'"$host_transport"$'\nlab-server' \
    "$fixture/vnc-host" lab-server
  expect_route "$client_os" "$helper"$'\n'"$host_transport"$'\nlab-server:1' \
    env VNC_HOST_ALIAS=lab-server "$fixture/vnc-srv16-umango"
  expect_route "$client_os" "$helper"$'\n'"$host_transport"$'\nvnc-host:1' \
    env VNC_HOST_ALIAS=vnc-host "$fixture/vnc-host"
done
printf '%s\n' 'PASS: user aliases, explicit targets, actions and transports on Linux/macOS/MSYS2/Cygwin'

if "$fixture/vnc-srv16-" > "$fixture/stdout" 2> "$fixture/stderr"; then
  echo 'FAIL: empty alias username accepted' >&2
  exit 1
fi
grep -q 'invalid host/user alias' "$fixture/stderr"
printf '%s\n' 'PASS: malformed user alias is rejected'

printf 'NAME\tDISPLAY\tSTATE\tLISTEN\nold\t:2\tstopped\tloopback\nlater\t:10\trunning\tloopback\n-\t:3\trunning\tunknown\n' > "$fixture/sessions"
for client_os in linux-gnu darwin msys; do
  case "$client_os" in
    linux-gnu) helper=vnc-session-view ;;
    darwin) helper=vnc-lan ;;
    msys) helper=vnc-win-lan ;;
  esac
  user_transport=--tunnel
  [[ "$client_os" != darwin ]] || user_transport=--auto
  expect_route "$client_os" "$helper"$'\n'"$user_transport"$'\numango@srv16:3' \
    "$fixture/vnc-srv16-umango"
  expect_route "$client_os" "$helper"$'\n'"$user_transport"$'\numango@srv16:self' \
    "$fixture/vnc-srv16-umango" self
done
printf '%s\n' 'PASS: user default chooses the lowest running owned display; explicit self stays named'

printf 'NAME\tDISPLAY\tSTATE\tLISTEN\nlater\t:10\tstopped\tloopback\nold\t:2\tstopped\tloopback\n' > "$fixture/sessions"
expect_route darwin $'vnc-lan\n--auto\numango@srv16:2' "$fixture/vnc-srv16-umango"
printf '%s\n' 'PASS: a retained stopped desktop is preferred over allocating another display'

touch "$fixture/fail"
rm -f "$fixture/actual"
if OSTYPE=darwin "$fixture/vnc-srv16-umango" > "$fixture/stdout" 2> "$fixture/stderr"; then
  echo 'FAIL: SSH listing failure was ignored' >&2
  exit 1
fi
[[ ! -e "$fixture/actual" ]]
printf '%s\n' 'PASS: listing failure does not create a replacement session or open a viewer'
