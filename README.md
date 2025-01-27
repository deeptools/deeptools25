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

> [!NOTE]
> First, adjust `./bench.cfg` as desired with a symlink to any of the CFG files bundled within,
> for example, to use our `olddata/`, execute this: `ln -sf old.cfg bench.cfg`

`snakemake --forcerun --use-conda --conda-frontend mamba --benchmark-extended --cores 12`


#### Clean

```{bash}
rm ./output/*.{png,txt}
cp -v ./data/output/*.{png,txt} ./output/
```
