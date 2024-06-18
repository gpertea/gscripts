#!/bin/bash
## checks the given bcf/vcf query against unimputed genotypes
## reports only the best matches for each query sample

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
       cut -f2- | sort -k1,1 -k6,6nr > "$gtcheck_out"
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
