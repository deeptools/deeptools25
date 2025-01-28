### Preparation

1. Install `rustc` >= 1.73
1. Install `Clang` and `Snakemake`. For example:

```{bash}
conda create -n snakemake clang-19 snakemake
conda activate snakemake
```

#### Download data

> WIP.

## Run benchmark

`snakemake --forcerun --use-conda --benchmark-extended --cores 12`
