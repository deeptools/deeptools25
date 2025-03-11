#!/bin/bash
set -euo pipefail
rm -rf trash && mv .snakemake trash
mkdir .snakemake && mv trash/conda .snakemake/
mv output trash/ && mv logs trash/
rm -rf trash
#git pull && \
#	snakemake --use-conda --forcerun --benchmark-extended \
#		--cores $(awk -F'=' '/Nthreads/ {print $2; exit}' Snakefile)
samtools quickcheck zenodo/*.bam && sleep 3s
snakemake --profile snk-slurm-exe --benchmark-extended
