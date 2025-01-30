### Dependencies

Just use an IDE that supports devcontainers. Or...

#### Manually

1. Install `rustc` >= 1.73
1. Install `Clang` and `Snakemake`. For example: `conda create -n snkmk clang-19 snakemake`

### Download data

> WIP.

1. Get GTF files from Ensembl for human and/ or wheat. If you are not inclined into benchmarking on full transcriptome, downsample them: `grep 'transcript_id' homo.v91.full.gtf | shuf | head -n 1000 | bedtools sort -i - > homo.v91.sample.gtf`

## Run benchmark

```{bash}
conda activate snkmk
snakemake --forcerun --use-conda --benchmark-extended --cores 12
```
