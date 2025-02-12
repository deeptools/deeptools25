#!/bin/bash
set -euo pipefail
mkdir -p ./{output,logs} && rm -rf ./{output,logs}/*
git pull && snakemake --use-conda --latency-wait 30 \
    --forcerun --benchmark-extended \
    --cores $(nproc --all)
