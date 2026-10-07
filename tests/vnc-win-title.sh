#!/usr/bin/env bash
set -euo pipefail

repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
fixture="$(mktemp -d)"
trap 'rm -rf -- "$fixture"' EXIT
export VNC_TITLE_FIXTURE="$fixture"
mkdir -p "$fixture/bin"

## Exercise the exact remote program, remapping only the distro utility path.
sed -n "/^set -euo pipefail$/,/^REMOTE_VNC_TITLE$/p" "$repo/vnc-win-lan" |
  sed -n '/^display="\$1"$/,/^REMOTE_VNC_TITLE$/p' |
  sed '$d' | sed "s|/usr/bin/|$fixture/bin/|g" > "$fixture/title.sh"
[[ -s "$fixture/title.sh" ]]

cat > "$fixture/utility" <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
name="${0##*/}"
mode="$(cat "$VNC_TITLE_FIXTURE/$name.mode")"
printf '%s %s\n' "$name" "$*" >> "$VNC_TITLE_FIXTURE/calls"
[[ "$mode" != incompatible ]] || { echo 'No VNC extension' >&2; exit 1; }
case "$3" in
  -get)
    [[ "$mode" != empty ]] || { echo 'Getting param desktop failed' >&2; exit 0; }
    cat "$VNC_TITLE_FIXTURE/title"
    ;;
  -set)
    [[ "$mode" != silent_failure ]] || { echo 'Setting param failed' >&2; exit 0; }
    printf '%s\n' "${4#desktop=}" > "$VNC_TITLE_FIXTURE/title"
    ;;
  *) exit 2 ;;
esac
MOCK
chmod +x "$fixture/utility"
cp "$fixture/utility" "$fixture/bin/tigervncconfig"
ln -s tigervncconfig "$fixture/bin/vncconfig"
cp "$fixture/utility" "$fixture/bin/tigervncconfig.ubuntu-1.13"
export PATH="$fixture/bin:$PATH"

reset_case() {
  printf '%s\n' "$1" > "$fixture/vncconfig.mode"
  printf '%s\n' "$1" > "$fixture/tigervncconfig.mode"
  printf '%s\n' "$2" > "$fixture/tigervncconfig.ubuntu-1.13.mode"
  printf '%s\n' "${3:-old-title}" > "$fixture/title"
  : > "$fixture/calls"
}

run_title() {
  bash -euo pipefail "$fixture/title.sh" :1 srv16:1 > "$fixture/stdout" 2> "$fixture/stderr"
}

reset_case compatible incompatible
run_title
[[ "$(cat "$fixture/title")" == srv16:1 ]]
! grep -q ubuntu "$fixture/calls"
printf '%s\n' 'PASS: matching current utility sets and verifies title'

reset_case incompatible compatible
run_title
[[ "$(cat "$fixture/title")" == srv16:1 ]]
grep -q 'tigervncconfig.ubuntu-1.13.*-set' "$fixture/calls"
! grep -q '^tigervncconfig ' "$fixture/calls"
printf '%s\n' 'PASS: old-server fallback and alias deduplication'

reset_case compatible incompatible srv16:1
run_title
! grep -q -- '-set' "$fixture/calls"
printf '%s\n' 'PASS: already-correct title needs no write'

reset_case empty compatible
run_title
[[ "$(cat "$fixture/title")" == srv16:1 ]]
printf '%s\n' 'PASS: failed read with exit zero triggers fallback'

reset_case silent_failure incompatible
if run_title; then echo 'FAIL: unverified title accepted' >&2; exit 1; fi
[[ "$(cat "$fixture/title")" == old-title ]]
grep -q 'title update failed' "$fixture/stderr"
printf '%s\n' 'PASS: failed write with exit zero is rejected'

reset_case incompatible incompatible
if run_title; then echo 'FAIL: incompatible utilities accepted' >&2; exit 1; fi
grep -q 'No compatible vncconfig' "$fixture/stderr"
[[ ! -s "$fixture/stdout" ]]
printf '%s\n' 'PASS: incompatible utilities produce a diagnostic and failure'
