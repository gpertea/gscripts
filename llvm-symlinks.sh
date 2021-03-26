#!/bin/bash
read -r -d '' USAGE << EOM
Usage:
 llvm-symlinks.sh <version>
 
 Symlink llvm and clang programs in /usr/bin/ for the given version.

Example:
 llvm-symlinks.sh 11
 
 It will symlink in /usr/bin/ the following programs:
 /usr/lib/llvm-11/bin/* (except scan-build and scan-view)
 /usr/lib/llvm-11/lib/clang/11.*/bin/hwasan_symbolize
 /usr/lib/llvm-11/share/clang/*.py
 /usr/share/clang/scan-build-11/bin/scan-build
 /usr/share/clang/scan-build-py-11/bin/scan-build-py
 /usr/share/clang/scan-view-11/bin/scan-view

EOM

if [ $# -ne 1 ]; then
  echo "$USAGE"
  exit 1
fi
ver="$1"
if [[ $ver == "-h" || $ver == "--help" ]]; then
  echo "$USAGE"
  exit 1
fi


cd /usr/bin

for fp in $(find /usr/lib/llvm-$ver/bin -maxdepth 1 -perm -u=x -a -not -type l -a -not -type d); do
  fn="${fp##*/}"
  if [[ $fn == "scan-view" || $fn == "scan-build" ]]; then
    continue 
  fi
  if [[ ! -f $fn || -L $fn ]]; then
    sudo /bin/rm -f $fn
    sudo ln -s $fp .
  fi
done

for fp in /usr/lib/llvm-$ver/lib/clang/$ver.*/bin/hwasan_symbolize \
  /usr/lib/llvm-$ver/share/clang/*.py /usr/share/clang/scan-build-$ver/bin/scan-build \
  /usr/share/clang/scan-build-py-$ver/bin/scan-build-py \
  /usr/share/clang/scan-view-$ver/bin/scan-view ; do
  fn="${fp##*/}"
  if [[ ! -f $fn || -L $fn ]]; then
    sudo /bin/rm -f $fn
    sudo ln -s $fp .
  fi
done

#special case:
fn=scan-build-py
if [[ ! -f $fn || -L $fn ]]; then
  sudo /bin/rm -f $fn
  sudo ln -s /usr/share/clang/scan-build-py-$ver/bin/scan-build $fn
fi

