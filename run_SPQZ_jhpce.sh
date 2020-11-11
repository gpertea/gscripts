#!/bin/bash
#$ -l bluejay,mem_free=40G,h_vmem=40G,h_fsize=150G
#$ -o ./SPEAQeasy_output.log
#$ -e ./SPEAQeasy_output.log
#$ -cwd
### you can also add  #$ -M <youremail> and submit this script 
### with qsub -m ea to be notified by email when the job ends or aborts

SPQZ=/dcl02/lieber/ajaffe/gpertea/SPEAQeasy
STRAND=reverse # or "unstranded" or "forward" if needed
OUTDIR=/dcl02/lieber/ajaffe/gpertea/RNAsp_work/2002UNHP-0326
INDIR=/dcl02/lieber/ajaffe/gpertea/data/2002UNHP-0326
## $INDIR/samples.manifest must exist

module use /jhpce/shared/jhpce/modulefiles/libd
module load nextflow
export _JAVA_OPTIONS="-Xms8g -Xmx10g"

#make sure these output directories exist
mkdir -p "$OUTDIR/wrk"
mkdir -p "$OUTDIR/results"

## -- now start the pipeline:
nextflow $SPQZ/main.nf \
    --sample "paired" \
    --reference hg38 \
    --strand $STRAND --force_strand \
    --trim_mode force \
    --annotation "/dcl01/lieber/ajaffe/Nick/SPEAQeasy/Annotation" \
    -with-report execution_reports/JHPCE_run.html \
    --no_biomart \
    -profile jhpce \
    --input  "$INDIR" \
         -w  "$OUTDIR/wrk" \
    --output "$OUTDIR/results"

## if you want to resume an interrupted run, add the '-resume \' line
## anywhere above the --output line (and after the nextflow line)


#  Produces a report for each sample tracing the pipeline steps
#  performed (can be helpful for debugging).
#
#  Note that the reports are generated from the output log produced in the above
#  section, and so if you rename the log, you must also pass replace the filename
#  in the bash call below.
echo "Generating per-sample logs for debugging..."
bash $SPQZ/generate_logs.sh $PWD/SPEAQeasy_output.log
