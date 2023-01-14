#!/usr/bin/env bash
if [[ $# -ne 2 ]]; then
 echo "usage:"
 echo " sshfs_mount.sh remote_path local_mount_path"
 exit 1
fi
if [[ ! -d $2 ]]; then
 echo "Error: local mount directory '$2' does not exist!"
 exit 1
fi
sshfs $1 $2
