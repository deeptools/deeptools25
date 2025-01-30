### Dependencies

Just use an IDE that supports devcontainers. Or...

#### Manually

1. Install `rustc` >= 1.73
1. Install `Clang` and `Snakemake`. For example: `conda create -n snkmk clang-19 snakemake`

### Download data

1. Download the data from [here](https://zenodo.org/record/14760356) (`10.5281/zenodo.14760356`), put these BAM files under `zenodo/` directory.
1. Get GTF files from Ensembl for human and/ or wheat, put them under `regions/` directory. If you are not inclined into benchmarking on full transcriptome, downsample: `grep 'transcript_id' homo.v91.full.gtf | shuf | head -n 1000 | bedtools sort -i - > homo.v91.sample.gtf`

## Run benchmark

```{bash}
conda activate snkmk
snakemake --forcerun --use-conda --benchmark-extended --cores 12
```
