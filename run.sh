#!/bin/bash
set -euo pipefail
mkdir -p ./{output,logs} && rm -rf ./{output,logs}/*
#git pull &&
#	snakemake --use-conda --forcerun --benchmark-extended \
#		--cores $(awk -F'=' '/Nthreads/ {print $2; exit}' Snakefile)
snakemake --profile snk-slurm-exe --benchmark-extended
