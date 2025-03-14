#!/bin/bash
set -euo pipefail

# Define the default organism (centralized here as the single source of truth)
DEFAULT_ORGANISM="homo"

# Use the provided organism or fall back to the default
ORGANISM=${1:-$DEFAULT_ORGANISM}
ORGANISM_ARG="--config organism=$ORGANISM"
ORGANISM_TAG="_$ORGANISM"

# Rest of the file, modified to use the arguments
rm -rf trash && mv .snakemake trash
mkdir .snakemake && mv trash/conda .snakemake/
mkdir -p {output,logs} && mv output trash/ && mv logs trash/
rm -rf trash

# Check input files were transferred okay, then execute the pipeline
samtools quickcheck zenodo/*.bam && sleep 3s
snakemake --profile snk-slurm-exe --benchmark-extended $ORGANISM_ARG

# Save everything except for data output files and flags
find output -mindepth 1 -maxdepth 1 -type d | xargs rm -rf
find output -type f \( -name "*_done_*.txt" -o -name "bechmark_iteration_*.txt" \) -delete
mv .snakemake/log/*.snakemake.log logs/
mv .snakemake/slurm_logs logs/
mv logs output/
mv output run_$(date +%Y%m%d%H%M%S)${ORGANISM_TAG}