#!/bin/bash
## checks the given bcf/vcf query against unimputed genotypes
## reports only the best matches for each query sample
## NOTE: this is better to be run in a new clean folder, e.g. mkdir gtcheck; cd gtcheck

gtpath=$HOME/genotyping/raw_data/unimputed/postmortem

function err_exit {
    echo -e "Error: $1"
    exit 1
}

# Check if no argument is given or if the argument is "-h" or "--help"
if [ $# -eq 0 ] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    err_exit "Usage: $0 <qryfile>"
fi

qryfile="$1"

# Check if $qryfile exists and if ${qryfile}.csi exists
if [[ ! -f "$qryfile" ]]; then
    err_exit "Query file $qryfile does not exist."
fi

if [[ ! -f "${qryfile}.csi" ]]; then
    err_exit "Index file ${qryfile}.csi does not exist."
fi

# Iterate through all files matching the pattern
for gtbcf in "$gtpath"/gt_*_n*/hg38_*.bcf; do
    # Extract the second token from the directory name
    gtn=$(basename $(dirname "$gtbcf") | cut -d'_' -f2)
    gtn=${gtn/.5Mv1./M}
    # Prepare the output file name
    gtcheck_out="gtcheck_${gtn}.tab"
    echo "output file: $gtcheck_out"
    # Run the bcftools gtcheck command
    {
      bcftools gtcheck -uGT,GT -E0 -g "$gtbcf" -Ot "$qryfile" |& grep -v -P '^(#|INFO|\[W::)' | \
       cut -f2- | sort -k1,1 -k6,6nr | awk '$3<$6 && $3<120' > "$gtcheck_out"
      if [[ ! -s $gtcheck_out ]]; then
         unlink "$gtcheck_out"
      fi 
    } &
    ### DEBUG/TEST ONLY
    ##exit
done
wait
gtcheck_out_hm="gtcheckBest.himatch.tab"
gtcheck_out_lmis="gtcheckBest.lomis.tab"
cat gtcheck_*.tab | sort -k1,1 -k6,6nr -k3,3n | awk -F'\t' '{
 if (!seen_hm[$1]++) {
     print $0
     first[$1] = $3"."$6
 } else if ($3"."$6 == first[$1]) {
     print $0
 }
}' > "$gtcheck_out_hm"
cat gtcheck_*.tab | sort -k1,1 -k3,3n -k6,6nr |  awk -F'\t' '{
 if (!seen_lmis[$1]++) {
     print $0
     first[$1] = $3"."$6
 } else if ($3"."$6 == first[$1]) {
     print $0
 }
}' > "$gtcheck_out_lmis"

sort -u <(cut -f1-3,5,6 "$gtcheck_out_hm") <(cut -f1-3,5,6 "$gtcheck_out_lmis") | sed 's/.cram//' > gtcheck_best_matches.tab

cut -f2 gtcheck_best_matches.tab | fgrep -wf- $HOME/work/genotyping/geno2brnum.tab | sort -k2,2 > gtcheck_geno2brnum.tab

### to further check the mappings with an existing phenodata:
#cut -f1 gtcheck_best_matches.tab | sort -u | fgrep -wf- ~/work/cbrain/processed_data/gencode41/bsp2_dlpfc.phenodata.tab | \
# cut -f2,3 | perl -pe 's/Br(\d\d\d)\b/Br0$1/' | sort -k2,2 > gtcheck_rnum2brnum.tab
## get BrNums that were not reported at all as best matches in gtcheck, but they are in the phenodata:
# comm -23 <(cut -f2 gtcheck_rnum2brnum.tab) <(cut -f2 gtcheck_geno2brnum.tab) > gtcheck_BrNums_not_found.lst
## what RNums they were supposed to match?
#fgrep -f gtcheck_BrNums_not_found.lst gtcheck_rnum2brnum.tab > gtcheck_rnum2brnum_missing.tab
## what genotypes were found to match instead?
# cut -f1 gtcheck_rnum2brnum_missing.tab | fgrep -wf- gtcheck_best_matches.tab > gtcheck_rnum_missing_bestmatch.tab
## see if these matching genotypes were reported before as mapped to a brain having multiple genotypes
#cut -f2 gtcheck_rnum_missing_bestmatch.tab | fgrep -f- ~/work/genotyping/brnum2genoID_postDNAswap_n2584.ctab
