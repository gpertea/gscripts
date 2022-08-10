#!/bin/bash
#shopt -s globstar
#mkdir idats
#ls -d ~/genotyping_iDATs/*_idat/**/*.idat | \
#  fgrep -f <(cut -f2 rebecca_33samples.tab) | xargs -I {} cp {} idats/

#ls -d ~/genotyping_iDATs/PennCNV_input_*/** | \
# fgrep -f <(cut -f2 rebecca_33samples.tab) > sifiles.ls
cat sifiles.ls | perl -mPOSIX -pe \
 'chomp;@t=split(/\./);
 $t=POSIX::strftime("\%y\%m\%d", localtime((stat $_)[9]));
 $_=$t[-1]."\t$t\t$_\n"' | \
 sort -k1,1 -k2,2nr > sifiles.tab
if [[ ! -d sifiles ]]; then
 mkdir sifiles
fi

## now rename every file to add a BrNum prefix
lastid=''
while read -r line; do
 IFS=$'\2' read -ra t <<< "${line//$'\t'/$'\2'}" # split by tabs into array t
 IFS=' '$'\t'$'\n' # restore IFS immediately
 echo ">>> processing gid: ${t[0]}"
 if [[ ${t[0]} != $lastid ]]; then
  echo "keeping $line"
  gid=${t[0]}
  lastid=$gid
  br=$(fgrep "$gid" rebecca_33samples.tab | cut -f1)
  if [[ -z $br ]]; then
    echo "No BrNum found for $gid"
    exit 1
  fi
  f=${t[2]} # full path to file
  fn=${f##*/} # just the filename
  sed "1 s/$gid\.//g" $f | \
   perl -pe 'chomp;@t=split(/\t/);@t[4,5]=($t[5],$t[4]);$_=join("\t",@t)."\n"' \
  > sifiles/$br.$fn.txt
 fi
 #line in $(cut -f1,2 ../rebecca_33samples.tab); do
  
done < sifiles.tab
