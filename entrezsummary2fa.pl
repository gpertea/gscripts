#!/usr/bin/perl
use strict;
my ($acc, $gi, $descr, $plen);
while (<>) {
  if (m/^\s*$/) { 
     if ($acc && $gi && $descr && $plen) {
       print ">$acc\t$plen\tGI:$gi\t$descr\n";
       ($acc, $gi, $descr, $plen)=('', 0, '', 0);
     }
     next
  }
  if (m/^\d+\.\s+(.+)/) {
     $descr=$1;
     next;
  }
  if (m/^(\S+)\sGI:(\d+)$/) {
    ($acc, $gi)=($1, $2);
    next;
  }
  if (m/^(\d+)\s+/) {
    $plen=$1;
  }
}

print ">$acc\t$plen\tGI:$gi\t$descr\n"
   if ($acc && $gi && $descr && $plen);

