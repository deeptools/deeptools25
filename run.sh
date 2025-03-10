#!/bin/bash
set -euo pipefail
mkdir -p .snakemake && find .snakemake -maxdepth 1 -mindepth 1 -not -name conda -type d -delete
rm -rf ./{output,logs} && mkdir -p ./{output,logs}
#git pull && \
#	snakemake --use-conda --forcerun --benchmark-extended \
#		--cores $(awk -F'=' '/Nthreads/ {print $2; exit}' Snakefile)
samtools quickcheck zenodo/*.bam
snakemake --profile snk-slurm-exe --benchmark-extended
