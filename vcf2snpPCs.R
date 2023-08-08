#!/bin/env Rscript
## script to get MDS and snpPCs for a set of genotypes in a VCF.gz input data
## on JHPCE, the plink module must be loaded (plink must be available in PATH)!
library(data.table)

## change here the path to the vcf file to calculate MDS for:
vcf='./raw_n128.vcf.gz'

## write the bed files with this bedout prefix 
## change the output prefix as needed
bedout="plink_out"
#bedout=sub('.vcf.gz', '', vcf, fixed=T)

## Optional QC filter (add it for )
qcflt='--geno 0.1 --maf 0.05 --hwe 0.000001'
## should not be applied to an already extracted VCF that has all the samples of interest
## do NOT apply the maf QC filter for very small sets 
## for a file with 68 genotypes, using qcflt like this eliminates over half of the variants!


## convert vcf to plink BED:
## NOTE: remove qcflt if this is a small sub-sample, not a large cohort
cmd=paste("plink --make-bed --output-chr chrM --keep-allele-order --vcf",
          vcf, qcflt, '--out',   bedout)
## add '--keep samples_to_extract.txt' if only a subset of genotyping samples are needed
system(cmd)

indfile=paste0(bedout, '_indep')

### independent and cluster
## --indep <window size>['kb'] <step size (variant ct)> <VIF threshold>
## produce a pruned subset of markers that are in approximate linkage equilibrium with each other
## --indep requires three parameters:
##          a window size in variant count or kilobase (if the 'kb' modifier is present) units,
##          a variant count to shift the window at the end of each step,
##          a variance inflation factor (VIF) threshold
cmd = paste("plink --bfile ", bedout, "--indep 100 10 1.25 --out", bedout)
system(cmd)

## MDS components
## use "--mds-plot 15" below if you want 15 SNPs
cmd = paste0("plink --bfile ", bedout,
             " --cluster --mds-plot 10 --extract ", bedout, ".prune.in --out ", bedout)
system(cmd)

# ## A transpose
cmd = paste("plink --bfile", bedout, "--recode A-transpose --out", bedout)
system(cmd)

#### read in MDS
mds = read.table(paste0(bedout, ".mds"), header=TRUE,as.is=TRUE)
rmds = mds[,-(1:3)] #remove FID, IID and SOL columns 1-3
colnames(rmds) = paste0("snpPC",1:ncol(rmds))
rmds$SAMPLE_ID=paste0(mds$FID, '_', mds$IID)
setcolorder(rmds, 'SAMPLE_ID')

## write snpPCs file as csv
fwrite(rmds, file=paste0(bedout, ".snpPCs.csv"), sep='\t')


