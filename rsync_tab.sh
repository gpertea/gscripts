#!/bin/bash
## expected input: file with 3 tab delimited columns
##    src_path groups destination_parent_dir
## Where "groups" is a comma delimited list of groups; the first entry is
## primary group that should get ownership; secondary groups are optional
## and they are currently not used in any way
###
## the script should be run with stdbuf -oL to enforce line-level buffering
##    stdbuf -oL rsync_tab.sh dir_info.tab > rsync_dir_info.log 

tbl=$1
if [[ -z $tbl || ! -f $tbl ]]; then
  echo "rsync table not provided!"
  exit 1
fi
## redirect stderr to stdout from now on
exec 2>&1

# read each input line into 3 variables (using tab as field delimiter)
while IFS=$'\t' read -r sourceDir groupNames destDir
do
  # Split the comma-separated group names into an array
  IFS=',' read -ra groups <<< "$groupNames"

  sourceDir=${sourceDir%/}
  lastDir=$(basename "$sourceDir")
  destDir=${destDir%/}
  if [[ ! -d $destDir ]]; then
     mkdir -p $destDir
  fi
  destPath="$destDir/$lastDir"
  destDir="$destDir/"
  # Get the last directory name from $sourceDir
  echo ">srcdir: $sourceDir"
  # Destination directory path
  echo "   dest: $destPath"
  echo "     primary group: ${groups[0]} | secondary groups: ${groups[@]:1}"
  ## copy directories while preserving executable attributes
  ##   but change permissions and ownership to the primary group
  cmd="rsync -a --chown :${groups[0]} --chmod=D770,F660 --executability '$sourceDir' '$destDir'"
  echo "  running: $cmd" 
  eval "stdbuf -oL $cmd"
  exitcode=$?
  # running the command here
  # check if rsync command succeeded
  if [ $exitcode -ne 0 ]; then
    echo " rsync to $destPath finished with non-zero exit code $exitcode"
  else 
    echo " rsync to $destPath finished. "
  fi
  #if needed, assign read permissions for all members of the secondary groups
  #for group in "${groups[@]:1}"; do
  #  ## use setfacl here
  #  # nfs4_setfacl -Rm g:"$group":rX "$destPath"
  #  echo "  2nd group: $group" >&2
  #done
done < $tbl
