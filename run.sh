#!/bin/bash
# USAGE: ./run.sh <local|slurm> <homo|triticum> <ntimes:int> <nthreads:int>

set -euo pipefail

if ! command -v snakemake &> /dev/null; then
    echo "snakemake could not be found. Please load the appropriate module."
    exit 1
fi

# Parse and validate arguments
COMPUTE=${1:-local}
if [[ ! "$COMPUTE" =~ ^(local|slurm)$ ]]; then
    echo "Error: Compute option must be 'local' or 'slurm'"
    exit 1
fi

FIRST_ARG=""
if [[ $COMPUTE == "slurm" ]]; then
    FIRST_ARG="--profile snk-slurm-exe"
elif [[ $COMPUTE == "local" ]]; then
    FIRST_ARG="--use-conda"
fi

ORGANISM=${2:-homo}
NTIMES=${3:-3}
NTHREADS=${4:-8}

# Directories
OUTPUT_DIR="output"
LOGS_DIR="logs"
SNAKEMAKE_DIR=".snakemake"
TRASH_DIR="trash"
CONDA_DIR="${SNAKEMAKE_DIR}/conda"

# Validate input files
echo "Checking input BAM files..."
if ! samtools quickcheck zenodo/*.bam; then
    echo "Error: BAM file check failed. Please verify your input files."
    exit 1
fi
sleep 3s

# Run the pipeline
echo "Running snakemake pipeline for organism: ${ORGANISM}"
snakemake $FIRST_ARG --benchmark-extended --config organism=$ORGANISM ntimes=$NTIMES nthreads=$NTHREADS

# Archive results
FINISHTIME=$(date +%Y%m%d)
echo "Cleaning up output directory..."
find "${OUTPUT_DIR}" -mindepth 1 -maxdepth 1 -type d | xargs rm -rf
find "${OUTPUT_DIR}" -type f \( -name "*_done_*.txt" -o -name "iter_count_*.txt" \) -delete

# Move log files
if [ -d "${SNAKEMAKE_DIR}/slurm_logs" ]; then
    mv "${SNAKEMAKE_DIR}/slurm_logs" "${LOGS_DIR}/"
fi
mv "${SNAKEMAKE_DIR}/log/"*.snakemake.log "${LOGS_DIR}/"
mv "${LOGS_DIR}" "${OUTPUT_DIR}/${LOGS_DIR}_bs99_${ORGANISM}"

# Rename output directory with timestamp
FINAL_DIR="run_t${NTHREADS}n${NTIMES}_${FINISHTIME}_${ORGANISM}"
mv "${OUTPUT_DIR}" "${FINAL_DIR}"
echo "Finished: ${FINAL_DIR}"
