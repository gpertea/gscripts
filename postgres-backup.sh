#!/bin/bash
#set -x

read -r -d '' USAGE << EOM
Usage:
 postgres-backup.sh server [db]
 
 The default db is 'rse'
EOM

if [[ $# -lt 1 || $1 == '-h' || $1 == '--help' ]]; then
  echo "$USAGE"
  exit 1
fi
PGSRV="$1"
PGDB="$2"
PGDB=${PGDB:=rse}
echo "Backing up database $PGDB on server $PGSRV"

HOST=$(hostname)

##TODO: make it work here with msmtp-mta ?
err_mail () {
  subj="postgres-backup error on $HOST"
  #if [[ -f "$1" ]]; then
     #mail -s "$subj" 'geo.pertea@gmail.com' < "$1"
  #else
  #  #echo "Error: $1" | mail -s "$subj" 'geo.pertea@gmail.com'
  #fi
  echo "Error: $1"
  exit 1
}

## wake-up the external drives:
## only needed on my setup:
#mount /mnt/p2box_m
#mount /mnt/p2box_f
##
#local work/backup directory:
bd=postgres
LDIR=/data/backups/$bd
db=$PGDB
H=/home/gpertea
## for cron execution:
export PATH=$H/gscripts:/opt/geo/bin:$H/bin:/usr/bin:/bin:/sbin:/usr/sbin
fpre="pgbck.$PGSRV.$db."
bfn=${fpre}$(date +%y%m%d_%H-%M)
host=$(hostname -s)
lastb="pgbck.$PGSRV.$db.last"
##TODO: the backup and (remote) copy directories should be taken from a 
##    local config file ~/.postgres-backup or so
#rdirs=( "/data2/backups/$bd" "srv02:/mnt/Data2/gpertea/backups/$bd" )
#rdirs=( "/data/gdebsrv_data1/backups/$bd" )
rdirs=( "/data2/backups/$bd" )

ferr=${fpre}stderr
echo "Backing up $db into: $LDIR/${bfn}*"
cd $LDIR

#mkdir $bfn
#if [[ $? -ne 0 ]] ; then
#    err_mail "could not create dir $bfn in $LDIR"
#fi
## postgres password for server $PGSRV must be in ~/.pgpass
pg_dump -Fc -h $PGSRV -U postgres -f $bfn  2>$ferr

## copy to remote directories, if specified
for rdest in "${rdirs[@]}" ; do
  #scp $dbf $rdest/ 2>>$frerr
   rsync -v $bfn $rdest/ 2>>$ferr
   if [[ $? != 0 ]]; then 
      echo "Warning: rsync failed for: rsync -v $bfn $rdest/" >> $ferr
   fi
done


# clean up files older than 30 days unless those backup files are smaller
#dbfs=$(stat -c'%s' "$dbf")
#oldbf=$(ls -1t *.${wiki}_db.sql.gz 2>/dev/null | tail -n +2 | head -1)
#oldbfs=$(stat -c'%s' "$oldbf")

#if [[ $cleanup_old -gt 0 ]]; then
#  #echo "Cleaning up backups older than 30 days"
#  if [[ $prevff ]]; then 
#     /bin/rm -f $prevff
#     fi
#  nice -n19 find . -maxdepth 1 -name 'bck_*.'$wiki'_*.*' -mtime +30 -exec /bin/rm -rf {} \;
#fi

err=$(grep -E -i 'error|fail|cannot|couldn|\bfault|timeout|discon' $ferr 2>/dev/null)
if  grep -q -i -E 'error|fail|cannot|couldn|\bfault|timeout|discon' $ferr 2>/dev/null; then
 err_mail $ferr
fi

rm -rf ./$lastb
mv $bfn $lastb

