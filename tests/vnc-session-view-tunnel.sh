#!/usr/bin/env bash
set -euo pipefail
repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
fixture="$(mktemp -d)"
trap 'rm -rf -- "$fixture"' EXIT
source <(sed -n '/^init_tunnel_state()/,/^if \[\[ "\$action"/p' "$repo/vnc-session-view" | sed '$d')
die() { printf '%s\n' "$*" >&2; exit 1; }
control_dir="$fixture"
control_path="$fixture/control-%C"
remote_host=test-host
remote_port=5901

allocate_tunnel_port() { local_port=$((40000 + ++allocations)); }
tunnel_control() { [[ -f "$fixture/running" ]]; }
remember_tunnel() { local_port="$(cat "$fixture/running")"; }
stop_tunnel() { :; }
ssh() {
  printf '%s\n' "$*" >> "$fixture/calls"
  case "$scenario" in
    race)
      if [[ "$local_port" == 40001 ]]; then
        echo 'bind [127.0.0.1]:40001: Address already in use' >&2
        return 1
      fi ;;
    occupied)
      echo 'bind [127.0.0.1]:40001: Address already in use' >&2
      return 1 ;;
    auth)
      echo 'Permission denied (publickey).' >&2
      return 1 ;;
  esac
  printf '%s\n' "$local_port" > "$fixture/running"
}
reset_case() {
  allocations=0
  scenario="$1"
  requested_port=""
  rm -f "$fixture/running"
  : > "$fixture/calls"
}

reset_case race
ensure_tunnel
[[ "$local_port" == 40002 && "$allocations" == 2 ]]
printf '%s\n' 'PASS: bind race allocates another port'
ensure_tunnel
[[ "$local_port" == 40002 && "$allocations" == 2 ]]
[[ "$(wc -l < "$fixture/calls")" == 2 ]]
printf '%s\n' 'PASS: retained tunnel reuses its actual port'

requested_port=5911
if (ensure_tunnel) 2>"$fixture/error"; then exit 1; fi
grep -q 'existing tunnel uses port 40002' "$fixture/error"
printf '%s\n' 'PASS: conflicting fixed port preserves the active tunnel'

reset_case auth
if (ensure_tunnel) 2>"$fixture/error"; then exit 1; fi
[[ "$(wc -l < "$fixture/calls")" == 1 ]]
grep -q 'Permission denied (publickey)' "$fixture/error"
printf '%s\n' 'PASS: authentication failure is not retried'

reset_case occupied
requested_port=5911
if (ensure_tunnel) 2>"$fixture/error"; then exit 1; fi
[[ "$(wc -l < "$fixture/calls")" == 1 ]]
printf '%s\n' 'PASS: occupied fixed port is not silently changed'

reset_case occupied
if (ensure_tunnel) 2>"$fixture/error"; then exit 1; fi
[[ "$(wc -l < "$fixture/calls")" == 5 ]]
grep -q 'after five attempts' "$fixture/error"
printf '%s\n' 'PASS: allocation retries are bounded'

reset_case success
requested_port=5911
ensure_tunnel
[[ "$local_port" == 5911 && "$allocations" == 0 ]]
printf '%s\n' 'PASS: fixed local-display port is honored'
