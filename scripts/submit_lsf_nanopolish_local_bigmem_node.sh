#!/usr/bin/env bash

MEMORY=1500000
LSF_CORES=110
QUEUE=bigmem
SNAKEMAKE_LOCAL_CORES=112
PROFILE="lsf"
LOG_DIR=logs/
JOB_NAME="snakemake_master_process."$(date --iso-8601='minutes')

mkdir -p $LOG_DIR

bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY] span[hosts=1]" \
    -M "$MEMORY" \
    -n "$LSF_CORES" \
    -P "$QUEUE" \
    -o "$LOG_DIR"/"$JOB_NAME".o \
    -e "$LOG_DIR"/"$JOB_NAME".e \
    -J "$JOB_NAME" \
      snakemake --snakefile Snakefile_nanopolish_local \
                --local-cores "$SNAKEMAKE_LOCAL_CORES" \
                --profile "$PROFILE" \
                --keep-going "$@"

exit 0
