#!/bin/bash
set -euo pipefail

# Tidy-up, start fresh
rm -rf trash && mv .snakemake trash
mkdir .snakemake && mv trash/conda .snakemake/
mkdir -p {output,logs} && mv output trash/ && mv logs trash/
rm -rf trash

## This is how we would run the pipeline on a local machine
#git pull && \
#	snakemake --use-conda --forcerun --benchmark-extended \
#		--cores $(awk -F'=' '/Nthreads/ {print $2; exit}' Snakefile)
## On an HPC cluster, we need to submit the job to the scheduler (see next block)

# Check input files were transferred okay, then execute the pipeline
samtools quickcheck zenodo/*.bam && sleep 3s
snakemake --profile snk-slurm-exe --benchmark-extended

# Save everything except for data output files and flags
find output -mindepth 1 -maxdepth 1 -type d | xargs rm -rf
find output -type f \( -name "*_done_*.txt" -o -name "bechmark_iteration_*.txt" \) -delete
mv .snakemake/log/*.snakemake.log logs/
mv .snakemake/slurm_logs logs/
mv logs output/
mv output run_$(date +%Y%m%d%H%M%S)