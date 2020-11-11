#!/bin/bash
## running the pipeline locally on a multi-core machine

SPQZ="/opt/SPEAQeasy" #this should not be changed usually

OUTDIR="/media/Data1_/gpertea/RNAsp/2002UNHP-0326"
INDIR="$OUTDIR/input" # samples.manifest must be here
STRAND=reverse # or "unstranded" or "forward" if needed
#make sure output directories exist
mkdir -p "$OUTDIR/wrk"
mkdir -p "$OUTDIR/results"

SPLOG=$PWD/SPEAQeasy_output.log
export _JAVA_OPTIONS="-Xms5g -Xmx7g"
$SPQZ/Software/nextflow $SPQZ/main.nf \
    --sample "paired" \
    --reference hg38 \
    --strand $STRAND \
    --annotation "$SPQZ/Annotation" \
    --input  "$INDIR" \
         -w  "$OUTDIR/wrk" \
    --output "$OUTDIR/results" \
    --no_biomart \
    -with-report execution_reports/pipeline_report.html \
    -profile docker_local \
    2>&1 | tee -a $SPLOG
    
echo "Generating per-sample logs for debugging..."
bash $SPQZ/scripts/generate_logs.sh $SPLOG
