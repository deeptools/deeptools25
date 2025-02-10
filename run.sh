#!/bin/bash
\rm -rf ./output/*
\rm -rf ./logs/*
git pull || :
snakemake \
    --cores $(nproc --all) \
    --latency-wait 120 \
    --forcerun \
    --use-conda \
    --benchmark-extended
