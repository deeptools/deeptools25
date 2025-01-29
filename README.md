### Dependencies

Just use an IDE that supports devcontainers. Or...

#### Manually

1. Install `rustc` >= 1.73
1. Install `Clang` and `Snakemake`. For example: `conda create -n snkmk clang-19 snakemake`

### Download data

> WIP.

## Run benchmark

```{bash}
conda activate snkmk
snakemake --forcerun --use-conda --benchmark-extended --cores 12
```
