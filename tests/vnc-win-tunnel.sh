#!/usr/bin/env bash
set -euo pipefail
repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
fixture="$(mktemp -d)"
trap 'rm -rf -- "$fixture"' EXIT
source <(sed -n '/^is_ssh_pid()/,/^if \[\[ "\$action"/p' "$repo/vnc-win-lan" | sed '$d')
die() { printf '%s\n' "$*" >&2; exit 1; }
control_dir="$fixture"
control_path="$fixture/control-%C"
remote_host=test-host
remote_port=5901

allocate_local_port() {
  local_port=$((40000 + ++allocations))
}
ssh() {
  printf '%s\n' "$*" >> "$fixture/calls"
  case "$scenario" in
    collision)
      if [[ "$*" == *127.0.0.1:40001:* ]]; then
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
}
reset_case() {
  allocations=0
  scenario="$1"
  unset VNC_LOCAL_PORT
  : > "$fixture/calls"
}

reset_case collision
start_control_master
[[ "$local_port" == 40002 && "$allocations" == 2 ]]
[[ "$(wc -l < "$fixture/calls")" == 2 ]]
printf '%s\n' 'PASS: allocation race retries with a different port'

reset_case auth
if (start_control_master) 2>"$fixture/error"; then exit 1; fi
[[ "$(wc -l < "$fixture/calls")" == 1 ]]
grep -q 'Permission denied (publickey)' "$fixture/error"
printf '%s\n' 'PASS: authentication failure is reported without retries'

reset_case occupied
VNC_LOCAL_PORT=40001
if (start_control_master) 2>"$fixture/error"; then exit 1; fi
[[ "$(wc -l < "$fixture/calls")" == 1 ]]
printf '%s\n' 'PASS: explicit occupied port fails without changing the request'

reset_case occupied
if (start_control_master) 2>"$fixture/error"; then exit 1; fi
[[ "$(wc -l < "$fixture/calls")" == 5 ]]
grep -q 'after five attempts' "$fixture/error"
printf '%s\n' 'PASS: repeated allocation races have a bounded retry count'

reset_case success
VNC_LOCAL_PORT=42000
start_control_master
[[ "$local_port" == 42000 && "$allocations" == 0 ]]
grep -q '127.0.0.1:42000:127.0.0.1:5901' "$fixture/calls"
printf '%s\n' 'PASS: explicit available port is honored'

bounded_control_command() { return 0; }
remember_control_master() { :; }
read_master_port() { local_port=43210; }
start_control_master() { echo 'unexpected new master' >&2; exit 1; }
reset_case success
ensure_control_master
[[ "$local_port" == 43210 && "$allocations" == 0 ]]
printf '%s\n' 'PASS: existing tunnel reuses its actual forwarding port'

VNC_LOCAL_PORT=43211
if (ensure_control_master) 2>"$fixture/error"; then exit 1; fi
grep -q 'existing tunnel uses port 43210' "$fixture/error"
printf '%s\n' 'PASS: conflicting override does not replace an active tunnel'

## Control identity and cleanup must never include a wildcard across hosts.
ssh() {
  printf 'controlpath %s/ssh-display-%s-endpoint-%s\n' \
    "$control_dir" "$remote_display" "$remote_host"
}
TMPDIR="$fixture"
remote_display=1
remote_host=host-a
init_control_state
first_socket="$control_socket"
flock -u 9
exec 9>&-
remote_host=host-b
init_control_state
[[ "$control_socket" != "$first_socket" && "$control_socket" != *%* ]]
flock -u 9
exec 9>&-
printf '%s\n' 'PASS: same display on different hosts has distinct control state'

perl -MSocket -e '
  for my $path (@ARGV) {
    socket(my $socket, PF_UNIX, SOCK_STREAM, 0) or die "socket: $!";
    bind($socket, sockaddr_un($path)) or die "bind: $!";
    close($socket);
  }
' "$first_socket" "$control_socket"
[[ -S "$first_socket" && -S "$control_socket" ]]
remove_control_socket
[[ -S "$first_socket" && ! -e "$control_socket" ]]
printf '%s\n' 'PASS: cleanup preserves another host control socket'
