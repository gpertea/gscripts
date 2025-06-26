#!/usr/bin/env bash
# Usage:  rsync_by_prefix.sh  <DEST>  <SDIR>  [GO]
#   DEST   Brain‐region keyword *and* destination directory name (e.g. Caudate)
#   SDIR   Absolute path that holds the FASTQ files
#   GO     Optional literal string “GO” – actually copies instead of dry-run
#
# The script:
#   - parses ../../metadata_n151_rnaseq425.tab for lines containing DEST
#   - extracts the sample prefixes (column 2) and rsyncs matching files
#   - fails if SDIR is missing or if no matching files are found
#
# Exit codes: 0 success ▏ 1 usage ▏ 2 SDIR missing ▏ 3 no matches

set -euo pipefail

function err_exit {
  echo -e "Error: $1" >&2;  local err=${2:-1};  exit $err; }

### basic CLI check
[[ $# -ge 2 && $# -le 3 ]] || err_exit "\nUsage: $(basename "$0") DEST SDIR [GO]"

dest=$1                     # keyword + output dir
sdir=$2                     # absolute source dir
mode=${3:-DRY}              # default DRY unless 3rd arg is GO

## metadata file location hard coded here, adjust as neeeded:
meta=../metadata_n151_rnaseq425.tab   # adjust if necessary

[[ -d $sdir ]] || err_exit "source dir '$sdir' not found." 2

### Collect prefixes (single fgrep call)
mapfile -t prefixes < <( fgrep "$dest" "$meta" | cut -f2 )
(( ${#prefixes[@]} )) || err_exit "no '$dest' entries in $meta." 3

### Quick probe: does *any* file in SDIR match?
shopt -s nullglob         # empty globs disappear instead of staying literal
found=''
for p in "${prefixes[@]}"; do
  for _ in "$sdir"/"$p"*; do   # inner loop runs only if a match exists
    found=1
    break 2                    # jump out of both loops
  done
done
shopt -u nullglob

[[ $found ]] || err_exit "no files in '$sdir' match any of the ${#prefixes[@]} prefixes." 4

### Build rsync filter rules from the prefixes ──────────────────────
inc_rules=()
for p in "${prefixes[@]}"; do
  inc_rules+=( "+ /${p}*" )
done
inc_rules+=( "- *" )        # exclude everything else

### Run rsync (dry-run unless GO) ───────────────────────────────────
opts=(-av --prune-empty-dirs --filter="merge -")
[[ $mode == GO ]] || { opts+=(-n); echo "Dry-run (add 'GO' to copy)." >&2; }

mkdir -p "$dest"
printf '%s\n' "${inc_rules[@]}" | rsync "${opts[@]}" "$sdir"/ "$dest"/
