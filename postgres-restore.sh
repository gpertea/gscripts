#!/bin/bash
read -r -d '' USAGE << EOM
Usage:
 postgres-restore.sh server pgdump.file
 
 Restore a postgres backup made with the -Fc option
 into the rse database
 
EOM

if [[ $# -lt 2 || $1 == '-h' || $1 == '--help' ]]; then
  echo "$USAGE"
  exit 1
fi

PGSRV="$1"
pgdump="$2"

if [[ ! -f $pgdump ]]; then
 echo "Backup file not found: $pgdump"
 exit 1
fi

echo "Restoring database backup $pgdump to server $PGSRV"

pg_restore -h $PGSRV -c -U postgres -d rse -v $pgdump

