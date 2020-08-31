#!/bin/bash
#set -x
#exec > $HOME/backup_wiki.stderr.log 2>&1
H=/home/gpertea
#local work/backup directory:
wiki=gwiki
BDIR=/mylocal/geo/backups/wiki
export PATH=$H/gscripts:/opt/geo/bin:$H/bin:/usr/bin:/bin:/sbin:/usr/sbin
bn=bck_$(date +%y%m%d_%H-%M)
host=$(hostname -s)
wikidir=/mylocal/geo/httpd/html/mediawiki
rdirs=( "/nfs/p2box/m/_backup/$wiki" "/nfs/p2box/n/_backup/$wiki" "gdebsrv:$H/backups/$wiki" "maestro:/data/backups/$wiki" \
 "salz:$H/backups/$wiki" "salz:/ccb/salz8-3/gpertea/_backups/$wiki")

ferr=$BDIR/$bn.stderr
echo "Backing up into: $BDIR/${bn}*"
cd $wikidir

##- check if LocalSettings.php has changed
files_changed=''
if ! cmp -s ./LocalSettings.php $BDIR/$wiki/LocalSettings.php; then
 cp -p ./LocalSettings.php $BDIR/$wiki/
 files_changed=1
fi

##- make the whole wiki read-only:
#nice -n 19 perl -i -pe 's/^[#\s]+(\$wgReadOnly[ =])/$1/' LocalSettings.php

#clean up pChart4mw cache files older than 20 days
if [[ -d uploads/pChart4mw_cache ]]; then
 cd uploads/pChart4mw_cache
 cache_size=$(du -s | cut -f1)
 if [[ $cache_size -gt 20000000 ]]; then
   echo "Cleaning up cache.." > $ferr
   ( nice -n19 find . -maxdepth 1 -mtime +20 -exec /bin/rm -rf {} \; ) 2>>$ferr
 fi
fi
##### --- see if wiki files have been updated (uploads, extensions, configuration etc.)
cd $wikidir
ff=$bn.${wiki}_files.tar.bz2
fsync=$BDIR/$bn.rsync.info
(
nice -n19 rsync -Wai --delete --exclude=uploads/pChart4mw_cache/'*' \
 --exclude=LocalSettings.php ./ $BDIR/$wiki > $fsync
) 2>>$ferr
if grep -q -E '^(>f|c)' $fsync 2>/dev/null; then
  #we have changes, build the tar file
  files_changed=1
#else
# rsync -aicI ./LocalSettings.php  ~/backups/wiki/gwiki/LocalSettings.php > $fsync
# if grep -q '^>f' $fsync 2>/dev/null; then
#    files_changed=1
# fi
fi

dbf="$bn.${wiki}_db.sql.gz"
echo "Backing up database with mysqldump ($dbf).." >> $ferr
#using ~/.my.cnf instead, removed: -u "$dbuser" -p"$dbpass"
( nice -n 19 mysqldump -h localhost --default-character-set=binary \
   --single-transaction $wiki -c | gzip > $BDIR/$dbf
) 2>>$ferr

### -- xml dump:
## --disable for now
#xmlf="$bn.${wiki}_xml.bz2"
#echo "Backing up pages as xml ($xmlf).." >> $ferr
#( nice -n19 php -d error_reporting=E_ERROR maintenance/dumpBackup.php --current | \
#  nice -n19 bzip2 -9 > ~/backups/wiki/$xmlf 
#) 2>>$ferr

##- restore wiki to read-write:
#nice -n 19 perl -i -pe 's/^\s*(\$wgReadOnly[ =])/## $1/' LocalSettings.php

##- prepare and send the backups to remote backup locations
cd $BDIR
prevff=$(ls -1t *.${wiki}_files.tar.bz2 2>/dev/null | head -1)
if [[ -z $prevff ]]; then
 files_changed=1
fi
if [[ $files_changed ]]; then
  echo "Backing up files ($ff).." >> $ferr
  #echo "Files changed, building $ff" >> $ferr
  nice -n 19 tar cfj $ff $wiki 2>> $ferr
  c=0
  for fn in $( ls -1t *.${wiki}_files.tar.bz2 2>/dev/null ) ; do
    (( c++ ))
    if [[ $c -gt 5 ]]; then
      /bin/rm -f $fn
    fi
  done
else
  if [[ $prevff ]]; then
    touch $prevff #make it current
    ff=$prevff
    prevff=''
  fi
fi
# clean up files older than 30 days unless those backup files are smaller
dbfs=$(stat -c'%s' "$dbf")
oldbf=$(ls -1t *.${wiki}_db.sql.gz 2>/dev/null | tail -n +2 | head -1)
oldbfs=$(stat -c'%s' "$oldbf")

cleanup_old=1
(( oldbfs -= 50000 ))
#echo "Adjusted oldb file size: $oldbfs"

errtext=$(grep -E -i 'error|fail|cannot|space|couldn|\bfault' $ferr 2>/dev/null)
#echo "errtext='$errtext'"


if [[ $dbfs -lt 1024 || $dbfs -lt $oldbfs || $errtext ]]; then
  cleanup_old=0
  echo "Warning: issues encountered during backup ($dbf, size $dbfs, prev. $oldbf size $oldbfs), check $ferr" >> $ferr
     mail -s "${wiki} backup issue" -r 'backuper@'$host.jhu.edu 'geo.pertea@gmail.com' < $ferr
fi

if [[ $cleanup_old -gt 0 ]]; then
  #echo "Cleaning up backups older than 30 days"
  if [[ $prevff ]]; then 
     /bin/rm -f $prevff
     fi
  nice -n19 find . -maxdepth 1 -name 'bck_*.'$wiki'_*.*' -mtime +30 -exec /bin/rm -rf {} \;
fi
frerr=$bn.ssh.stderr

#ff=$bn'.${wiki}_files.tar.bz2'
frerr='ssh_test.log'
/bin/rm -f $frerr
scperr=''

for rdest in "${rdirs[@]}" ; do
 #scp $dbf $rdest/ 2>>$frerr
 rsync -v $dbf $rdest/ 2>>$frerr
 if [[ $? != 0 ]]; then 
    scperr=1
    echo "Warning: rsync failed for: rsync -v $dbf $rdest/" >> $frerr
 fi
 #scp -q $xmlf $rdest/ 2>>$frerr
 #if [[ $? != 0 ]]; then scperr=1; fi
 if [[ $files_changed ]]; then
   #scp $ff $rdest/ 2>>$frerr
   rsync -v $ff $rdest/ 2>>$frerr
   if [[ $? != 0 ]]; then 
      echo "Warning: rsync failed for: rsync -v $ff $rdest/" >> $frerr
      scperr=1
   fi
 fi
 oifs=$IFS
 IFS=':'
 arr=($rdest)
 IFS=$oifs
 if [[ ${#arr[@]} -gt 1 ]]; then
   rhost=${arr[0]}
   rdir=${arr[1]}
 else
   rhost=''
   rdir=$rdest
 fi

 rbasedir=${rdir%/*} #the remote "backups" directory
 #rtdir=${rdir##*/} #target subdirectory under the "backups" directory
 if [[ -z $rhost ]]; then
   echo "Notifying local target: $rbasedir/backup_received.sh '$rdir' '$bn' '$cleanup_old' '$ff'\"" | tee -a $frerr
   $rbasedir/backup_received.sh "$rdir" "$bn" "$cleanup_old" "$ff" 2>&1 | tee -a $frerr
   echo " -- target $rbasedir done --" | tee -a $frerr
 else
   echo "Notifying remote target: ssh $rhost \"$rbasedir/backup_received.sh '$rdir' '$bn' '$cleanup_old' '$ff'\"" | tee -a $frerr
   ssh $rhost "$rbasedir/backup_received.sh '$rdir' '$bn' '$cleanup_old' '$ff'" 2>&1 | tee -a $frerr
   echo " -- remote target $rdest done --" | tee -a $frerr
 fi
done

err=$(grep -E -i 'error|fail|cannot|couldn|\bfault|timeout|discon' $frerr 2>/dev/null)
if [[ $err || $scperr ]]; then
  echo "err = <$err>" >> $frerr
  echo "scperr = <$scperr>" >> $frerr
  mail -s "$wiki backup target issue"  -r 'backuper@'$host.jhu.edu 'geo.pertea@gmail.com' < $frerr
fi
