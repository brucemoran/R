#! /usr/env/R

##used as a follow-on from aligning with STAR, featureCounts
##as per standard nf-core/rnaseq pipeline (v1.4 currently)
argsIn <- commandArgs(trailingOnly = TRUE)

source(argsIn[1])

##create output dir
OUTDR <- "results/DESeq2"
dir.create(OUTDR, showWarnings = FALSE)

print(paste0("Arguments: ", argsIn))

nf_core_rnaseq_DESeq2(TAG = argsIn[2], SAMPLE_COND = argsIn[3], COND_DESIGN = argsIn[4], REF_COND = argsIn[5], OUTDIR = argsIn[6], ENS_VERSION = argsIn[7], ORG_PREFIX = argsIn[8], DELIM_SAMP = argsIn[9])
