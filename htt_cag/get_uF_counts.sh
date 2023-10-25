#!/bin/bash
## run this in an interactive session or as a one-shot job
for r in R*; do
 cd $r
   ../get_HTT_CAGs.pl > CAG_counts_uF.tab
 cd ..
done
