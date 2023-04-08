#!/usr/bin/env bash
## CHESS-like basic filtering of a gffcompare combined gtf resulted
## from merging multiple StringTie outputs
function err_exit {
 echo -e "Error: $1"
 exit 1
}

ref=$HOME/work/ref/gencode43.nri.main.gtf
#xglst=$HOME/work/cbrain/chess3/chess3m_xgenes.glst
xglst=$HOME/work/cbrain/chess3/vs_gencode43/g43_xgenes.glst

read -r -d '' USAGE << EOT
 Filter large merged combined gffcompare gtf based on number of samples
 and (optional) classification codes, excluding transfrags spanning across 
 multiple genes.
 
 Usage:
  gcmp_comb_flt.sh min_num_samples gffcmp.combined.gtf
 
   Using reference: $ref
 Using xgenes list: $xglst
EOT

if [[ $# -ne 2 || ! -f $2  ]]; then
  echo "$USAGE"
  exit 1
fi

for f in $ref $xglst; do
 if [[ ! -f $f ]]; then
    err_exit "cannot find $f"
 fi
done

ns=$1
ingtf=$2

if ((ns<2)); then
  err_exit "invalid minimum number of samples $ns, must be >=2"
fi
fb=${ingtf/.annotated.gtf/}
fb="${fb/.combined.gtf/}" #just in case
echo "$ingtf : base file name=$fb"

# get all transcripts found with the same structure in at least 2 samples, include SET
# this discards repeat and contained transfrags!
ofb=$fb.flt_n${ns}
gtfout=$ofb.gtf
fjt=$ofb.jtable
if [[ -f $gtfout && $(stat -c%s $gtfout) -gt 100000 ]]; then
  echo -e "["$(date '+%m/%d %H:%M')"]\t$gtfout already exists, skipping"
else 
  echo -e "["$(date '+%m/%d %H:%M')"]\tfiltering by num_smp>=$ns to $gtfout"
  gffcmp_flt -M -n${ns} -S $ingtf > $gtfout
fi

# map this GTF to reference:
echo -e "["$(date '+%m/%d %H:%M')"]\tmapping $gtfout to reference annotation"
trmap -J $ref $gtfout > $fjt

# of these, keep only transfrags NOT spanning multiple genes (unless known to be overlapping)
echo -e "["$(date '+%m/%d %H:%M')"]\tgetting list of transfrags to keep"
jtable_flt $xglst $fjt $ofb.excl.jtable > $ofb.tokeep.tlst

## fish out those transcripts:
fkept=$ofb.kept.gtf
echo -e "["$(date '+%m/%d %H:%M')"]\twriting final tranfrags to $fkept"
fgrep -wf $ofb.tokeep.tlst $gtfout > $fkept
echo -e "["$(date '+%m/%d %H:%M')"]\t"$(gtfcount $fkept)" transfrags kept ($ofb.tokeep.lst)"

