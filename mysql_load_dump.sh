#!/bin/sh
srv="$1"
db="$2"
user="$3"
dumpfile="$4"

if [[ -z "$dumpfile" ]] ; then
 echo -e "Usage:\n mysql_load_dump.sh <srv> <db> <user> <mysqldump_file>"
 echo " (<mysqldump_file> can be a compressed file with a .gz extension)"
 echo "Example:"
 echo "mysql_load_dump.sh localhost igmwiki igmusr \ "
 echo " /mylocal/geo/backups/igmwiki/bck_140826_03-15.ccbwiki_db.sql.gz"
 exit 0
fi
dext=${dumpfile##*.}
if [[ $dext == "gz" ]]; then
 #echo "Using: mysql -p -h $srv -D $db -u $user"
 gzip -cd $dumpfile | mysql -p -h $srv -D $db -u $user
else
 mysql -p -h $srv -D $db -u $user < $dumpfile
fi

